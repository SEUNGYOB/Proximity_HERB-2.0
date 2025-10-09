"""
Network proximity (Z-score) computation CLI runner
--------------------------------------------------

이 스크립트는 HERB 2.0 – STRING Human PPI 네트워크를 기반으로,
두 노드 집합(A, B) 간의 **Network Proximity (근접도)** 를 계산합니다.

주요 기능:
-----------
1. SQLite DB(`DB.db`)에서 Human_PPI 테이블을 불러와 NetworkX 그래프를 생성합니다.
2. `Utils.get_ensp_ids()`를 사용하여 HERB/HBIN/HBDIS 등의 ID를 ENSP 노드 집합으로 매핑합니다.
3. `proximity_util.compute_network_distances_CPU()`를 호출하여
   - 평균 최단거리 (shortest path mean distance)
   - Z-score (degree-matched 랜덤 네트워크 기반)
   를 계산합니다.
4. 결과를 콘솔(JSON 형식)로 출력하거나, `--out` 인자를 사용하여 JSON 파일로 저장할 수 있습니다.

인자 설명:
-----------
--db           : SQLite 데이터베이스 경로 (기본값: ./Data/DB.db)
--table        : PPI 테이블 이름 (기본값: Human_PPI)
--src-col      : PPI 테이블의 source 컬럼명 (기본값: protein1)
--tgt-col      : PPI 테이블의 target 컬럼명 (기본값: protein2)
--a-id         : 첫 번째 노드 집합 식별자 (예: HERB002168)
--b-id         : 두 번째 노드 집합 식별자 (예: HBDIS001345)
--random-time  : degree-matched resampling 반복 횟수 (기본값: 100)
--seed         : 랜덤 시드 (기본값: 42)
--workers      : 병렬 워커 수 (기본값: 16)
--out          : 결과 JSON 저장 경로 (선택사항)

참고:
------
- 거리 행렬 캐시는 증분(memmap) 방식으로 `./Data/Human_PPI/rows_mm.npy` 경로에 자동 생성됩니다.
- 풀 매트릭스 계산(`build_dist_matrix_from_graph`)은 포함되어 있지 않으며,
  대규모 실험 시 필요할 경우 별도 구현 가능.
"""

#!/usr/bin/env python3
# run_proximity.py
# Network proximity (Z-score) CLI runner — incremental cache only version

import argparse
import json
import os
import sqlite3
from pathlib import Path

import pandas as pd
import networkx as nx
import Utils
import proximity_util as pu


# =========================================
# 1. Graph Loading
# =========================================
def load_edges_from_db(db_path, table_name="Human_PPI", source_col="protein1", target_col="protein2"):
    """DB에서 PPI 엣지 로드 ('9606.' 접두사 제거 포함)"""
    if not Path(db_path).exists():
        raise FileNotFoundError(f"DB not found: {db_path}")

    conn = sqlite3.connect(db_path)
    query = f"SELECT {source_col}, {target_col} FROM {table_name}"
    df = pd.read_sql_query(query, conn)
    conn.close()

    def strip_prefix(x):
        if isinstance(x, str) and x.startswith("9606."):
            return x.split("9606.", 1)[-1]
        return x

    df[source_col] = df[source_col].map(strip_prefix)
    df[target_col] = df[target_col].map(strip_prefix)
    edges = list(df.itertuples(index=False, name=None))
    print(f"✅ Loaded {len(edges):,} edges from '{table_name}' (prefix '9606.' removed)")
    return edges


def create_network_from_edges(edges):
    """엣지 리스트로부터 networkx 그래프 생성"""
    G = nx.Graph()
    G.add_edges_from(edges)
    print(f"✅ Graph built: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")
    return G


# =========================================
# 2. Helper Functions
# =========================================
def ids_to_ensp_set(identifier: str):
    """Utils.get_ensp_ids() 래퍼 (빈 세트 경고 포함)"""
    ensp_ids = set(Utils.get_ensp_ids(identifier) or [])
    print(f"🔎 {identifier}: mapped to {len(ensp_ids)} ENSP nodes")
    if not ensp_ids:
        print(f"⚠️ Warning: {identifier} → empty ENSP set (check mapping or DB)")
    return ensp_ids


# =========================================
# 3. Main Logic
# =========================================
def run(args):
    # 1) 그래프 로드
    edges = load_edges_from_db(
        db_path=args.db,
        table_name=args.table,
        source_col=args.src_col,
        target_col=args.tgt_col,
    )
    G = create_network_from_edges(edges)

    # 2) ENSP 노드 집합 변환
    A = ids_to_ensp_set(args.a_id)
    B = ids_to_ensp_set(args.b_id)
    if not A or not B:
        raise SystemExit("❌ Error: Empty A/B set detected — check HERB/DISEASE ID validity")

    # 3) 근접도 계산 (memmap incremental cache 사용)
    print("ℹ️ Using incremental memmap cache (./Data/Human_PPI/rows_mm.npy)")
    res = pu.compute_network_distances_CPU(
        G=G,
        dist=None,  # 증분 캐시 모드
        A=A,
        B=B,
        random_time=args.random_time,
        seed=args.seed,
        max_workers=args.workers,
    )

    # 4) 결과 출력
    result = {
        "A_id": args.a_id,
        "B_id": args.b_id,
        "random_time": args.random_time,
        "seed": args.seed,
        "result": res,
    }

    print("\n==== Proximity Result ====")
    print(json.dumps(result, ensure_ascii=False, indent=2))

    # 5) 결과 저장
    if args.out:
        out_path = Path(args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(result, ensure_ascii=False, indent=2))
        print(f"💾 Saved result to: {out_path.resolve()}")


# =========================================
# 4. CLI Parser
# =========================================
def parse_args():
    p = argparse.ArgumentParser(description="Compute network proximity (Z-score) between two node sets.")
    p.add_argument("--db", default="./Data/DB.db", help="Path to SQLite DB (default: ./Data/DB.db)")
    p.add_argument("--table", default="Human_PPI", help="Edge table name (default: Human_PPI)")
    p.add_argument("--src-col", default="protein1", help="Source column (default: protein1)")
    p.add_argument("--tgt-col", default="protein2", help="Target column (default: protein2)")
    p.add_argument("--a-id", required=True, help="A set identifier (e.g., HERB002168)")
    p.add_argument("--b-id", required=True, help="B set identifier (e.g., HBDIS001345)")
    p.add_argument("--random-time", type=int, default=100, help="Degree-matched resampling count (default: 100)")
    p.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    p.add_argument("--workers", type=int, default=16, help="Max worker threads (default: 16)")
    p.add_argument("--out", help="Optional: save result JSON to this path")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args)
