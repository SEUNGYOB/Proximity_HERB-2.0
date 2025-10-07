#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A program for extracting, processing, and scraping data, including parallel
scraping, from HERB APIs and outputs in designated formats like CSV files.

The script includes functionalities such as data validation, ID extraction,
HTTP requests handling, file writing, and batch processing of entity data
(e.g., Herbs, Ingredients, and Diseases). It supports parallel execution to
speed up operations with large data sets.

Functions:
- resolve_out_dir: Resolves the absolute path for the output directory.
- _extract_ids_from_string: Extracts specific ID patterns from a string.
- _walk_and_collect_ids: Recursively walks through JSON-like objects to
  collect IDs.
- _post_detail_api: Sends POST requests to fetch entity details from HERB API.
- fetch_ids_only_raw: Fetches and optionally saves IDs from the HERB API.
- _write_edges_bulk: Writes or appends unique edges in CSV format.
- _read_codes_file: Reads code values from a text file.
- _build_code_list: Builds a list of codes from various input parameters.
- scrape_parallel: Scrapes data from HERB API in parallel for specified entities.
"""
from __future__ import annotations
from pathlib import Path
import argparse
import ast
import json
import re
import sys
from typing import Any, Dict, List, Tuple

import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

# ============================ Constants ============================
BASE = "http://47.92.70.12"
API  = f"{BASE}/chedi/api/?"
UA   = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
    "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126.0 Safari/537.36"
)

# 엔티티 메타데이터: label / code prefix / 엣지 파일명
ENTITY_INFO: Dict[str, Dict[str, str]] = {
    "Herb":       {"label": "Herb",       "prefix": "HERB", "edge": "HERB_edge.csv"},
    "Ingredient": {"label": "Ingredient", "prefix": "HBIN", "edge": "ING_edge.csv"},
    "Formula":    {"label": "Formula",    "prefix": "HBFO", "edge": "FORMULA_edge.csv"},
    "Target":     {"label": "Target",     "prefix": "HBTAR", "edge": "TARGET_edge.csv"},
    "Disease":    {"label": "Disease",    "prefix": "HBDIS", "edge": "DISEASE_edge.csv"},
}

SCRIPT_DIR = Path(__file__).resolve().parent
HB_ID_RE = re.compile(r"\bHB[A-Z0-9]+\b")

# Only keep IDs that match these 9 prefixes followed by digits
ALLOWED_PREFIXES = ("HBIN", "HBDIS", "HBTAR", "HBCT", "HBMA", "HERB", "HBFO", "HBREF", "HBEXP")
ALLOWED_RE = re.compile(r"^(?:" + "|".join(ALLOWED_PREFIXES) + r")\d+$")

_DEBUG_ON = False

# ============================ Paths ============================

def resolve_out_dir(out_dir: str | Path) -> Path:
    p = Path(out_dir)
    return p if p.is_absolute() else (SCRIPT_DIR / p)

# ============================ ID extraction ============================

def _extract_ids_from_string(s: str) -> set[str]:
    if not s:
        return set()
    return set(HB_ID_RE.findall(s))

def _walk_and_collect_ids(obj: Any) -> set[str]:
    """JSON-like 구조를 재귀 탐색하여 HB* ID 추출."""
    ids: set[str] = set()
    if obj is None:
        return ids
    if isinstance(obj, str):
        ids |= _extract_ids_from_string(obj)
        if "{" in obj and "}" in obj:
            try:
                maybe = ast.literal_eval(obj)
                if isinstance(maybe, (dict, list, tuple)):
                    ids |= _walk_and_collect_ids(maybe)
            except Exception:
                pass
        return ids
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(k, str):
                ids |= _extract_ids_from_string(k)
            ids |= _walk_and_collect_ids(v)
        return ids
    if isinstance(obj, (list, tuple)):
        for it in obj:
            ids |= _walk_and_collect_ids(it)
        return ids
    try:
        return _extract_ids_from_string(str(obj))
    except Exception:
        return ids

# ============================ HTTP ============================

def _post_detail_api(code: str, *, label: str, timeout: int = 60, cookie: str | None = None) -> Dict[str, Any] | List[Any]:
    payload = {"v": code, "label": label, "key_id": code, "func_name": "detail_api"}
    headers = {
        "User-Agent": UA,
        "Accept": "*/*",
        "Content-Type": "application/json",
        "Origin": BASE,
        "Referer": f"{BASE}/Detail/?v={code}&label={label}",
    }
    if cookie:
        headers["Cookie"] = cookie

    resp = requests.post(API, json=payload, headers=headers, timeout=timeout)
    resp.raise_for_status()
    try:
        return resp.json()
    except ValueError:
        return json.loads(resp.text)

def fetch_ids_only_raw(code: str, label: str, *, timeout: int = 60, cookie: str | None = None,
                       save_json: bool = False, out_dir: str = "./Data/HERB 2.0/api", debug: bool = False) -> List[str]:
    data = _post_detail_api(code, label=label, timeout=timeout, cookie=cookie)
    if save_json:
        base = resolve_out_dir(out_dir); base.mkdir(parents=True, exist_ok=True)
        (base / f"{code}__raw.json").write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    ids = sorted(_walk_and_collect_ids(data))
    if debug:
        print(f"[DEBUG] {label.upper()} {code}: {len(ids)} ids")
    return ids

# ============================ CSV write (single pass) ============================

def _canon_pair(a: str, b: str) -> Tuple[str, str]:
    a = a.strip()
    b = b.strip()
    return (a, b) if a <= b else (b, a)

def _write_edges_bulk(pairs: List[Tuple[str, str]], *, out_dir: str, total_name: str) -> Path:
    base = resolve_out_dir(out_dir)
    base.mkdir(parents=True, exist_ok=True)
    total_path = base / total_name

    if not pairs:
        return total_path

    df_new = pd.DataFrame(pairs, columns=["Node_1", "Node_2"], dtype="string")
    # Canonicalize undirected edges so (A,B) and (B,A) are the same
    df_new[["Node_1", "Node_2"]] = df_new.apply(lambda r: pd.Series(_canon_pair(r["Node_1"], r["Node_2"])), axis=1)
    if total_path.exists():
        try:
            df_old = pd.read_csv(total_path, dtype="string")
        except Exception:
            df_old = pd.DataFrame(columns=["Node_1", "Node_2"], dtype="string")
        # Canonicalize old edges as well before concatenation
        if not df_old.empty and set(["Node_1", "Node_2"]).issubset(df_old.columns):
            df_old[["Node_1", "Node_2"]] = df_old.apply(lambda r: pd.Series(_canon_pair(r["Node_1"], r["Node_2"])), axis=1)
        df_all = pd.concat([df_old, df_new], ignore_index=True)
        df_all.drop_duplicates(subset=["Node_1", "Node_2"], inplace=True)
    else:
        df_all = df_new

    tmp = total_path.with_suffix(".tmp.csv")
    df_all.to_csv(tmp, index=False)
    tmp.replace(total_path)
    return total_path

# ============================ Code list builder ============================

def _read_codes_file(path: Path) -> List[str]:
    codes: List[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        codes.append(s)
    return codes

def _build_code_list(*, prefix: str, code: str | None, codes: List[str] | None,
                     codes_file: str | None, range_str: str | None) -> List[str]:
    final_codes: List[str] = []
    if codes:
        final_codes.extend(codes)
    if codes_file:
        p = Path(codes_file)
        if not p.exists():
            raise FileNotFoundError(f"codes file not found: {p}")
        final_codes.extend(_read_codes_file(p))
    if range_str:
        try:
            s, e = range_str.split("-", 1)
            s_i, e_i = int(s), int(e)
        except Exception:
            raise ValueError(f"invalid range format: {range_str} (expected like 1-500)")
        lo, hi = min(s_i, e_i), max(s_i, e_i)
        final_codes.extend([f"{prefix}{n:06d}" for n in range(lo, hi + 1)])
    if code:
        final_codes.append(code)
    if not final_codes:
        raise ValueError("No codes provided. Use --code / --codes-file / --range.")
    return final_codes

# ============================ Unified parallel scraper ============================

def scrape_parallel(*,
                    entity: str,
                    code: str | None = None,
                    codes: List[str] | None = None,
                    codes_file: str | None = None,
                    range_str: str | None = None,
                    out_dir: str = "./Data/HERB 2.0/api",
                    timeout: int = 60,
                    cookie: str | None = None,
                    save_json: bool = False,
                    debug: bool = False,
                    max_workers: int = 8,
                    edge_name: str | None = None) -> Dict[str, int]:
    """범용 병렬 스크레이퍼 — Herb/Ingredient/Formula/Target/Disease"""
    global _DEBUG_ON
    _DEBUG_ON = bool(debug)

    if entity not in ENTITY_INFO:
        raise ValueError(f"Unknown entity '{entity}'. Choose one of: {', '.join(ENTITY_INFO)}")

    meta = ENTITY_INFO[entity]
    label = meta["label"]
    prefix = meta["prefix"]
    total_edge = edge_name or meta["edge"]

    final_codes = _build_code_list(prefix=prefix, code=code, codes=codes,
                                   codes_file=codes_file, range_str=range_str)

    summary: Dict[str, int] = {}
    pairs: List[Tuple[str, str]] = []

    def flush_partial():
        nonlocal pairs
        if pairs:
            _write_edges_bulk(pairs, out_dir=out_dir, total_name=total_edge)
            if debug:
                print(f"[FLUSH] wrote {len(pairs)} pairs (partial)")
            pairs = []

    def task(c: str) -> Tuple[str, List[str]]:
        ids = fetch_ids_only_raw(c, label, timeout=timeout, cookie=cookie,
                                 save_json=save_json, out_dir=out_dir, debug=debug)
        return c, ids

    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        future_map = {ex.submit(task, c): c for c in final_codes}
        for fut in as_completed(future_map):
            c = future_map[fut]
            try:
                code_ret, ids = fut.result()
                # keep only 9 kinds: prefix + digits
                ids_f = [x for x in ids if isinstance(x, str) and ALLOWED_RE.match(x.strip())]
                summary[code_ret] = len(ids_f)
                for x in ids_f:
                    pairs.append(_canon_pair(code_ret, x))
                # 100개 단위로 flush
                if len(pairs) >= 100:
                    flush_partial()
            except Exception as e:
                print(f"[ERROR] {c}: {e}")
                summary[c] = 0

    # 남은 데이터 마지막 flush
    flush_partial()

    print(f"[EDGE] finished → {resolve_out_dir(out_dir) / total_edge}")
    return summary



def scrape_all_entities(out_dir: str = "./Data/HERB 2.0/api",
                        max_workers: int = 12,
                        debug: bool = True) -> dict[str, dict[str, int]]:
    """
    HERB 2.0의 엔티티(Herb, Ingredient, Formula, Target, Disease)에 대해서 edge를 수집하는 함수
    max_workers 는 CPU 상황에 맞춰서 변경
    모든 엔티티(Herb, Ingredient, Formula, Target, Disease)를 순차적으로 실행.
    각 엔티티의 summary를 모아 반환.
    """
    results: dict[str, dict[str, int]] = {}

    # 엔티티별 최대 범위 사전
    ENTITY_RANGES = {
        "Herb": 7263,
        "Ingredient": 49258,
        "Formula": 6743,
        "Target": 15515,
        "Disease": 30170,
    }

    for entity, max_num in ENTITY_RANGES.items():
        print(f"\n[RUN] {entity} (1-{max_num}) 시작…")
        summary = scrape_parallel(
            entity=entity,
            range_str=f"1-{max_num}",
            out_dir=out_dir,
            max_workers=max_workers,
            debug=debug,
        )
        results[entity] = summary
        print(f"[DONE] {entity}: {len(summary)} codes processed")
    return results

def merge_edge_csvs(input_dir: str = "./Data/HERB 2.0/api",
                    output_file: str = "./Data/HERB 2.0/api/HERB_total_edges.csv") -> Path:
    """
    엔티티별 edge CSV 파일들을 하나로 합쳐 저장 (무방향 그래프). Node_1/Node_2 순서와 무관하게 중복 제거.

    Parameters
    ----------
    input_dir : str
        edge CSV 파일들이 들어 있는 디렉토리
    output_file : str
        최종 병합된 CSV 경로

    Returns
    -------
    Path
        최종 병합된 CSV 파일 경로
    """
    input_dir = resolve_out_dir(input_dir)
    csv_files = list(input_dir.glob("*_edge.csv"))

    if not csv_files:
        raise FileNotFoundError(f"No edge CSV files found in {input_dir}")

    dfs = []
    for f in csv_files:
        try:
            df = pd.read_csv(f, dtype="string")
            # Canonicalize undirected edges before appending
            if not df.empty and set(["Node_1", "Node_2"]).issubset(df.columns):
                df[["Node_1", "Node_2"]] = df.apply(lambda r: pd.Series(_canon_pair(r["Node_1"], r["Node_2"])), axis=1)
            dfs.append(df)
            print(f"[OK] Loaded {f} ({len(df)} rows)")
        except Exception as e:
            print(f"[WARN] Failed to load {f}: {e}")

    merged = pd.concat(dfs, ignore_index=True).drop_duplicates(subset=["Node_1", "Node_2"])

    # Ensure parent directory exists and normalize path relative to script dir
    out_path = Path(output_file)
    if not out_path.is_absolute():
        parent = resolve_out_dir(Path(output_file).parent)
        parent.mkdir(parents=True, exist_ok=True)
        out_path = parent / Path(output_file).name
    else:
        out_path.parent.mkdir(parents=True, exist_ok=True)

    merged.to_csv(out_path, index=False)
    print(f"[DONE] Merged {len(merged)} unique edges → {out_path}")

    return out_path

# ============================ CLI ============================

def main(argv: List[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Unified HERB API parallel scraper (Herb/Ingredient/Formula/Target/Disease)")
    p.add_argument("--entity", required=True, choices=list(ENTITY_INFO.keys()), help="Entity type")
    p.add_argument("--code", help="Single code (e.g., HERB000001 / HBIN000001 / HBFO000001 …)")
    p.add_argument("--codes-file", help="Text file with one code per line")
    p.add_argument("--codes", nargs="*", help="Explicit list of codes")
    p.add_argument("--range", dest="range_str", help="Inclusive numeric range like 1-500")
    p.add_argument("--out", default="./Data/HERB 2.0/api", help="Output directory")
    p.add_argument("--timeout", type=int, default=60, help="HTTP timeout seconds")
    p.add_argument("--cookie", help="Cookie header value if needed")
    p.add_argument("--save-json", action="store_true", help="Save raw JSON per code")
    p.add_argument("--debug", action="store_true", help="Verbose debug logs")
    p.add_argument("--max-workers", type=int, default=8, help="Thread parallelism")
    p.add_argument("--edge-name", help="Override edge CSV filename (optional)")

    args = p.parse_args(argv)
    print(f"[OUT] Files will be saved under: {resolve_out_dir(args.out)}")

    try:
        summary = scrape_parallel(
            entity=args.entity,
            code=args.code,
            codes=args.codes,
            codes_file=args.codes_file,
            range_str=args.range_str,
            out_dir=args.out,
            timeout=args.timeout,
            cookie=args.cookie,
            save_json=args.save_json,
            debug=args.debug,
            max_workers=args.max_workers,
            edge_name=args.edge_name,
        )
    except Exception as e:
        print(f"[ERROR] {e}")
        return 2

    any_saved = any(v > 0 for v in summary.values()) if summary else False
    return 0 if any_saved else 1

if __name__ == "__main__":

    #엔티티별 edge 정보 수집

    # If CLI args are provided, use CLI; otherwise, run all entities batch.
    if len(sys.argv) > 1:
        sys.exit(main())
    else:
        all_summary = scrape_all_entities(
            out_dir="./Data/HERB 2.0/api",
            max_workers=12,
            debug=True,
        )
        print("\n=== 전체 결과 ===")
        for ent, summ in all_summary.items():
            print(f"{ent}: {sum(summ.values())} IDs saved")

    #각 엔티티별로 edge를 다운 받은 다음 하나의 csv로 합침.
    merge_edge_csvs("./Data/HERB 2.0/api", "./Data/HERB 2.0/api/HERB 2.0_total_edges.csv")
