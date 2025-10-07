#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module for graph distance calculation and analysis.

This module provides utilities for reading, converting, and analyzing graph
data using NetworkX and igraph libraries. It includes functions to cache
shortest distance matrices for efficiency, perform pairwise distance calculations,
and compute high-level distance metrics for graph analyses. The module also
supports multiprocessing for distance metric computations.

Dependencies
------------
- NetworkX
- igraph
- numpy
- pathlib
- concurrent.futures
- logging
- time
- os
- random

"""
from __future__ import annotations
from pathlib import Path
import logging
import time
import igraph as ig
from concurrent.futures import ThreadPoolExecutor, as_completed
import os, random, numpy as np
from typing import Sequence, Dict, Any
import networkx as nx

###############################################################################
# 로깅 설정
###############################################################################
logging.basicConfig(
    format="[%(levelname)s] %(message)s",
    level=logging.INFO,
)
log = logging.getLogger(__name__)

###############################################################################
# 상수 & 헬퍼
###############################################################################
CACHE_DIR = Path(".cache"); CACHE_DIR.mkdir(exist_ok=True)
SENTINEL   = np.uint16(0xFFFF)   # unreachable 거리 표시 (65 535)

###############################################################################
# 0. NetworkX → igraph, 거리행렬 캐시
###############################################################################

def nx_to_igraph(G: nx.Graph):
    g = ig.Graph.TupleList(G.edges(), directed=False)
    g.vs["name"] = list(G.nodes())
    name2idx = {n: i for i, n in enumerate(g.vs["name"])}
    return g, name2idx
SENT = 65535          # DistMatrix에서 쓰는 sentinel 값

def _shortest_mean(mat: np.ndarray, rows: np.ndarray, cols: np.ndarray) -> float:
    """
    A×B 부분행렬 평균. sentinel(SENT) 값은 제외.
    """
    sub = mat[rows][:, cols]          # (|A|, |B|) 뷰
    valid = sub != SENT               # True = 유한 거리
    if not valid.any():               # 완전히 단절된 경우
        return float("inf")
    return float(sub[valid].mean())

def build_matrix(g: ig.Graph, out_file: Path) -> np.ndarray:
    log.info("⚙️  building all‑pairs shortest‑path matrix … (one‑time)")
    t0 = time.perf_counter()
    dist = np.array(g.shortest_paths_dijkstra(), dtype=np.float64)
    dist[np.isinf(dist)] = SENTINEL               # ∞ → sentinel
    dist_u16 = dist.astype(np.uint16)
    np.save(out_file, dist_u16)
    log.info("✅ matrix saved → %s  [%.1f MB, %.1fs]", out_file, dist_u16.nbytes/1e6, time.perf_counter()-t0)
    return dist_u16


def load_graph_assets(
    edgelist: str | Path,
    *,
    matrix_path: str | Path | None = None,
    rebuild: bool = False,
):
    """(G_nx, DistMatrix) 반환. 캐시가 있으면 메모리‑맵, 없으면 생성."""
    edgelist = Path(edgelist)
    log.info("📂 reading edgelist: %s", edgelist)
    G = nx.read_edgelist(edgelist, data=False)
    log.info("   ↳ %d nodes / %d edges", G.number_of_nodes(), G.number_of_edges())

    if matrix_path is None:
        matrix_path = CACHE_DIR / f"{edgelist.stem}_dist.npy"
    else:
        matrix_path = Path(matrix_path)
        matrix_path.parent.mkdir(parents=True, exist_ok=True)

    if matrix_path.exists() and not rebuild:
        log.info("🗄  using cached matrix %s", matrix_path)
        mat = np.load(matrix_path, mmap_mode="r")
        # Keep matrix indices aligned with a cached node-order file
        names_path = Path(str(matrix_path).replace(".npy", ".names.npy"))
        if names_path.exists():
            cached_names = np.load(names_path, allow_pickle=True).tolist()
            name2idx = {str(n): i for i, n in enumerate(cached_names)}
            log.info("🧭 loaded cached node order: %s", names_path.name)
        else:
            # Fallback (less safe): assume current G.nodes() order matches the matrix
            name2idx = {n: i for i, n in enumerate(G.nodes())}
            log.warning("⚠️ node-order cache missing (%s). Using G.nodes() order; indices may misalign.", names_path.name)
    else:
        log.info("🚧 cache miss → calculating matrix")
        g_ig, name2idx = nx_to_igraph(G)
        mat = build_matrix(g_ig, matrix_path)
        # Persist node order next to the matrix for future exact alignment
        names_path = Path(str(matrix_path).replace(".npy", ".names.npy"))
        np.save(names_path, np.array(g_ig.vs["name"], dtype=object))
        log.info("🧭 saved node order → %s", names_path.name)

    return G, DistMatrix(mat, name2idx)

###############################################################################
# 1. 거리 조회 래퍼
###############################################################################
class DistMatrix:
    """NPY 행렬 기반 최단거리 look‑up (µs)"""

    def __init__(self, mat: np.ndarray, name2idx: Dict[str, int]):
        self.mat = mat
        self.idx = name2idx

    def get(self, u: Any, v: Any) -> float:
        try:
            i, j = self.idx[str(u)], self.idx[str(v)]
        except KeyError:
            return float("inf")
        d = self.mat[i, j]
        return float("inf") if d == SENTINEL else float(d)

###############################################################################
# 2. PairwiseLengths
###############################################################################
class PairwiseLengths:
    def __init__(self, dist: DistMatrix, A: Sequence[Any], B: Sequence[Any]):
        self.d = dist; self.A = list(A); self.B = list(B)

    def build(self):
        len_AB, len_BA, len_AA, len_BB = {}, {}, {}, {}
        for a in self.A:
            len_AB[a] = {b: self.d.get(a, b) for b in self.B}
            len_AA[a] = {x: self.d.get(a, x) for x in self.A if x != a}
        for b in self.B:
            len_BA[b] = {a: len_AB[a][b] for a in self.A}
            len_BB[b] = {y: self.d.get(b, y) for y in self.B if y != b}
        return len_AB, len_BA, len_AA, len_BB

###############################################################################
# 3. NetworkDistance (필요 지표만 유지)
###############################################################################
class NetworkDistance:
    def __init__(self, tables, A, B):
        self.len_AB, self.len_BA, self.len_AA, self.len_BB = tables
        self.A, self.B = A, B

    def _closest(self):
        dA = [min(self.len_AB[a].values()) for a in self.A]
        dB = [min(self.len_BA[b].values()) for b in self.B]
        return float(np.mean(dA + dB))

    def _separation_inner(self, group, table):
        if len(group) <= 1: return 0.0
        return float(np.mean([min(table[n].values()) for n in group]))

    def separation(self):
        return self._closest() - (
            self._separation_inner(self.A, self.len_AA) +
            self._separation_inner(self.B, self.len_BB)
        ) / 2

    def shortest(self):
        return float(np.mean([
            *[d for a in self.A for d in self.len_AB[a].values()],
            *[d for b in self.B for d in self.len_BA[b].values()],
        ]))

###############################################################################
# 4. high‑level API
###############################################################################

# Z-score 계산을 위한 CPU 사용───────────────────────────────────────────────────────────

def compute_network_distances_CPU_NPY(G: nx.Graph,
                                  dist: DistMatrix,
                                  A: Sequence[Any],
                                  B: Sequence[Any],
                                  *, random_time=100,
                                  seed=42,
                                  max_workers=None
                                  ):
    import os, random
    from concurrent.futures import ThreadPoolExecutor, as_completed

    dist_mat_np = dist.mat
    name2idx = dist.idx
    A = [n for n in A if n in G]
    B = [n for n in B if n in G]
    if not A or not B: raise ValueError("노드 없음")

    A_idx = np.asarray([name2idx[n] for n in A], dtype=np.int32)
    B_idx = np.asarray([name2idx[n] for n in B], dtype=np.int32)

    # ── 1) 원본 거리 ─────────────────────────────────────
    d_orig = _shortest_mean(dist_mat_np, A_idx, B_idx)

    # ── 2) degree bucket 준비 ───────────────────────────
    deg2 = {}
    for n in G:
        deg2.setdefault(G.degree[n], []).append(name2idx[n])

    rng = random.Random(seed)
    def one_sample(k):
        rA = np.fromiter((rng.choice(deg2[G.degree[a]]) for a in A), int, len(A))
        rB = np.fromiter((rng.choice(deg2[G.degree[b]]) for b in B), int, len(B))
        return _shortest_mean(dist_mat_np, rA, rB)

    # ── 3) 병렬 실행 (스레드) ───────────────────────────
    max_workers = max_workers or min(32, (os.cpu_count() or 1)*4)
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        rnd_vals = [f.result()
                    for f in as_completed(ex.submit(one_sample, i)
                                          for i in range(random_time))]

    rnd_vals = [v for v in rnd_vals if np.isfinite(v)]  # ← 추가: ∞ 샘플 제거
    if not rnd_vals:
        raise RuntimeError("모든 무작위 샘플이 연결되지 않았습니다.")

    mean, std = float(np.mean(rnd_vals)), float(np.std(rnd_vals))
    z = 0.0 if std == 0 else (d_orig - mean) / std
    p = sum(v <= d_orig for v in rnd_vals) / len(rnd_vals)  # 표본 수 변경
    # print(f"🔍 d_orig = {d_orig:.4f}")
    # print(f"🔍 mean(random) = {mean:.4f}")
    # print(f"🔍 std(random) = {std:.4f}")
    # print(f"🔍 z = {z:.4f}, p = {p:.4f}, #samples = {len(rnd_vals)}")

    return {"shortest": d_orig,
            "Z_score": {"d": d_orig,"z":z,"mean":mean,"std":std,"p":p}}

def compute_network_distances_CPU(G: nx.Graph,
                                  dist: DistMatrix | None,
                                  A: Sequence[Any],
                                  B: Sequence[Any],
                                  *, random_time=100,
                                  seed=42,
                                  max_workers=None
                                  ):
    """
    Compute A–B network distance & Z-score.

    Incremental caching added:
      - If `dist` is provided (DistMatrix), use the original fast matrix path (no cache).
      - If `dist` is None, compute on-demand using BFS *and* cache per-source rows
        under "./Data/Human_PPI" so subsequent calls get faster.

    Cache layout (when `dist is None`)
    ----------------------------------
      ./Data/Human_PPI/nodes_order.json    # stable node order for row vectors
      ./Data/Human_PPI/rows/<b64(name)>.npy  # uint16 distances for one source row
        - value 65535 denotes unreachable (∞)
    """
    import os, random, json, base64
    from pathlib import Path as _Path
    from concurrent.futures import ThreadPoolExecutor, as_completed

    # sanitize inputs
    A = [str(n) for n in A if n in G]
    B = [str(n) for n in B if n in G]
    if not A or not B:
        raise ValueError("노드 없음")

    # ------------------------------
    # Helper for matrix fast path
    # ------------------------------
    def _fast_path_matrix(dist: DistMatrix):
        dist_mat_np = dist.mat
        name2idx = dist.idx
        A_idx = np.asarray([name2idx[n] for n in A if n in name2idx], dtype=np.int32)
        B_idx = np.asarray([name2idx[n] for n in B if n in name2idx], dtype=np.int32)
        if A_idx.size == 0 or B_idx.size == 0:
            raise ValueError("거리행렬 인덱스에 없는 노드만 남았습니다.")
        # 1) observed distance
        d_orig = _shortest_mean(dist_mat_np, A_idx, B_idx)
        # 2) degree buckets (indices)
        deg2: dict[int, list[int]] = {}
        for n in G:
            if n in name2idx:
                deg2.setdefault(G.degree[n], []).append(name2idx[n])
        rng = random.Random(seed)
        def one_sample_matrix(_k):
            rA = np.fromiter((rng.choice(deg2[G.degree[a]]) for a in A if G.degree[a] in deg2), int)
            rB = np.fromiter((rng.choice(deg2[G.degree[b]]) for b in B if G.degree[b] in deg2), int)
            if rA.size == 0 or rB.size == 0:
                return float("inf")
            return _shortest_mean(dist_mat_np, rA, rB)
        max_workers_eff = max_workers or min(32, (os.cpu_count() or 1)*4)
        with ThreadPoolExecutor(max_workers=max_workers_eff) as ex:
            rnd_vals = [f.result() for f in as_completed(ex.submit(one_sample_matrix, i)
                                                         for i in range(random_time))]
        return d_orig, rnd_vals

    # ------------------------------
    # Helper for incremental row cache (single memmap file)
    # ------------------------------
    import threading, time  # local import for locks & fsync pacing

    def _ensure_cache_and_order():
        """
        Prepare cache directory and node order file. Returns:
          base_dir, mmap_path, nodes(list[str]), name2pos(dict[str,int])
        """
        base = _Path("./Data/Human_PPI")
        base.mkdir(parents=True, exist_ok=True)
        mmap_path = base / "rows_mm.npy"            # <— single file cache
        order_path = base / "nodes_order.json"

        nodes_now = [str(n) for n in G.nodes()]
        if order_path.exists():
            with open(order_path, "r", encoding="utf-8") as f:
                nodes_saved = json.load(f)
            # exact node-set match required for consistent indexing
            if len(nodes_saved) != len(nodes_now) or set(nodes_saved) != set(nodes_now):
                raise RuntimeError("nodes_order.json과 현재 그래프 노드셋이 일치하지 않습니다. 캐시를 비우거나 재생성하세요: ./Data/Human_PPI")
            nodes = nodes_saved
        else:
            nodes = nodes_now
            with open(order_path, "w", encoding="utf-8") as f:
                json.dump(nodes, f, ensure_ascii=False)

        name2pos = {n: i for i, n in enumerate(nodes)}
        return base, mmap_path, nodes, name2pos

    def _open_or_init_memmap(mmap_path: _Path, N: int) -> np.memmap:
        """
        Open existing memmap or initialize a new NxN uint16 matrix with SENTINEL(65535).
        """
        exists = mmap_path.exists()
        mode = "r+" if exists else "w+"
        mm = np.memmap(mmap_path, dtype=np.uint16, mode=mode, shape=(N, N))
        if not exists:
            mm[:] = 65535  # ∞ sentinel
            mm.flush()
            # small sleep to ensure filesystem metadata settles (macOS)
            time.sleep(0.01)
        return mm

    # per-row locks to prevent concurrent writes to the same row
    _row_locks: dict[int, threading.Lock] = {}
    _row_locks_guard = threading.Lock()

    def _row_lock(idx: int) -> threading.Lock:
        with _row_locks_guard:
            lock = _row_locks.get(idx)
            if lock is None:
                lock = threading.Lock()
                _row_locks[idx] = lock
            return lock

    def _ensure_row_built(src_idx: int, src_name: str,
                          mm: np.memmap,
                          name2pos: dict[str, int]) -> np.ndarray:
        """
        Ensure that row for src_idx is materialized in the memmap.
        Uses a lock + atomic write pattern to avoid torn rows.
        Returns the row view (np.ndarray of shape (N,)).
        """
        row = mm[src_idx]
        # fast path: any non‑sentinel present means row already built
        if (row != 65535).any():
            return row

        lock = _row_lock(src_idx)
        with lock:
            row = mm[src_idx]
            if (row != 65535).any():
                return row  # another thread built it meanwhile

            # Build via BFS
            lengths = nx.single_source_shortest_path_length(G, src_name)
            # write directly into the memmap row
            for dst, d in lengths.items():
                j = name2pos.get(str(dst))
                if j is not None and d < 65535:
                    row[j] = np.uint16(d)
            mm.flush()
            return row

    def _shortest_mean_rows_memmap(A_idx: np.ndarray, B_idx: np.ndarray,
                                   nodes: list[str], mm: np.memmap, name2pos: dict[str,int]) -> float:
        vals = []
        for i in A_idx:
            row = _ensure_row_built(i, nodes[i], mm, name2pos)
            sub = row[B_idx]
            m = sub.min(initial=np.uint16(65535))
            if m != 65535:
                vals.append(int(m))
        if not vals:
            return float("inf")
        return float(np.mean(vals))

    def _slow_path_cached():
        base, mmap_path, nodes, name2pos = _ensure_cache_and_order()
        N = len(nodes)
        mm = _open_or_init_memmap(mmap_path, N)

        # indices for A, B within the fixed node order
        A_idx = np.asarray([name2pos[a] for a in A if a in name2pos], dtype=np.int32)
        B_idx = np.asarray([name2pos[b] for b in B if b in name2pos], dtype=np.int32)
        if A_idx.size == 0 or B_idx.size == 0:
            return float("inf"), []

        # 1) observed distance using cached memmap rows
        d_orig = _shortest_mean_rows_memmap(A_idx, B_idx, nodes, mm, name2pos)

        # 2) degree buckets (by indices)
        deg2: dict[int, list[int]] = {}
        for n in G:
            j = name2pos.get(str(n))
            if j is not None:
                deg2.setdefault(G.degree[n], []).append(j)

        rng = random.Random(seed)

        def one_sample_cached(_k: int) -> float:
            # degree‑matched resampling by indices
            try:
                rA_idx = np.fromiter((rng.choice(deg2[G.degree[nodes[i]]]) for i in A_idx), dtype=np.int32, count=A_idx.size)
                rB_idx = np.fromiter((rng.choice(deg2[G.degree[nodes[j]]]) for j in B_idx), dtype=np.int32, count=B_idx.size)
            except KeyError:
                return float("inf")
            return _shortest_mean_rows_memmap(rA_idx, rB_idx, nodes, mm, name2pos)

        max_workers_eff = max_workers or min(32, (os.cpu_count() or 1)*4)
        with ThreadPoolExecutor(max_workers=max_workers_eff) as ex:
            rnd_vals = [f.result() for f in as_completed(ex.submit(one_sample_cached, i)
                                                         for i in range(random_time))]
        return d_orig, rnd_vals

    # ------------------------------
    # Choose path & compute
    # ------------------------------
    if dist is not None:
        d_orig, rnd_vals = _fast_path_matrix(dist)
    else:
        d_orig, rnd_vals = _slow_path_cached()

    # summarize
    rnd_vals = [v for v in rnd_vals if np.isfinite(v)]
    if not rnd_vals:
        raise RuntimeError("모든 무작위 샘플이 연결되지 않았습니다.")
    mean, std = float(np.mean(rnd_vals)), float(np.std(rnd_vals))
    z = 0.0 if std == 0 else (d_orig - mean) / std
    p = sum(v <= d_orig for v in rnd_vals) / len(rnd_vals)

    return {"shortest": d_orig,
            "Z_score": {"d": d_orig, "z": z, "mean": mean, "std": std, "p": p}}
