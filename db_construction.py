#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
HERB 2.0 × STRING Human PPI Database Builder
--------------------------------------------

Builds a unified SQLite database (`DB.db`) integrating
HERB 2.0 and STRING v12.0 Human PPI data for network-based
analysis (e.g., Network Proximity, GNN).

Main Steps
-----------
1. **Human_PPI Construction**
   - Load STRING file (`9606.protein.links.v12.0.txt`)
   - Filter by `combined_score ≥ 700`
   - Save as table `Human_PPI`

2. **HERB 2.0 Construction**
   - Load CSVs under `./Data/HERB 2.0/`
   - Infer table: Herb_info, Formula_info, Ingredient_info,
     Target_info, Disease_info
   - Load edge file (`HERB 2.0_total_edges.csv`) → `HERB_edges`

3. **Edge Cleaning**
   - Keep biologically meaningful pairs:
     Formula–Herb, Herb–Ingredient, Ingredient–Target, Target–Disease
   - Optionally apply ADME filters (default: off)
   - Save as `HERB_edges_clean`

Output
-------
SQLite file: `./Data/DB.db`
  ├── Human_PPI
  ├── Herb_info / Formula_info / Ingredient_info / Target_info / Disease_info
  ├── HERB_edges
  └── HERB_edges_clean

Author: Seungyob Yi
Last updated: 2025-10
"""

import pandas as pd
from pathlib import Path
import sqlite3


#1. Human_PPI_Construction

def build_ppi_db(filepath: str | Path | None = None, cutoff: int = 700) -> pd.DataFrame:
    """
    ... (docstring 동일) ...
    Returns
    -------
    pd.DataFrame
        필터링된 데이터프레임
    """
    if filepath is None:
        filepath = "./Data/Human_PPI/9606.protein.links.v12.0.txt"
    filepath = Path(filepath)

    # 공백 다중 구분자 안전
    df = pd.read_csv(filepath, sep=r"\s+", engine="python")
    df = df[df["combined_score"] >= cutoff].reset_index(drop=True)

    # DB 보장
    db_path = Path("./Data/DB.db")
    db_path.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(db_path) as conn:
        df[["protein1","protein2","combined_score"]].to_sql("Human_PPI", conn, if_exists="replace", index=False)
        # 인덱스 추천
        cur = conn.cursor()
        cur.execute("CREATE INDEX IF NOT EXISTS ix_hppi_p1 ON Human_PPI(protein1)")
        cur.execute("CREATE INDEX IF NOT EXISTS ix_hppi_p2 ON Human_PPI(protein2)")
        conn.commit()

    print(f"[DONE] Human_PPI written ({len(df):,} edges, cutoff={cutoff})")
    return df

#2-1. Herb 2.0 Construction - node
herb_info_path = "./Data/HERB 2.0/api/HERB 2.0_total_edges.csv"
db_path = "./Data/DB.db"

#2-2. Herb 2.0 Construction - node

# Mapping from file-kind to target table name
TABLE_NAME_MAP = {
    "herb": "Herb_info",
    "formula": "Formula_info",
    "target": "Target_info",
    "disease": "Disease_info",
    "ingredient": "Ingredient_info",
}

def _infer_kind_from_filename(p: Path) -> str | None:
    """
    Infer which entity a CSV belongs to based on filename.
    We match whole tokens between separators to avoid the 'HERB_' prefix being read as 'herb'.
    Priority order: disease > formula > ingredient > target > herb
    """
    name = p.stem.lower()
    # split on common separators
    import re
    tokens = re.split(r"[^\w]+", name)  # keep alnum/underscore as token chars
    tokens = [t for t in tokens if t]   # drop empties

    # priority list ensures we don't misclassify 'HERB_disease_info' as 'herb'
    priority = ["disease", "formula", "ingredient", "target", "herb"]
    for kind in priority:
        if kind in tokens:
            return kind

    # secondary: regex for *_{kind}_* pattern
    for kind in priority:
        if re.search(rf"(?:^|[_\-]){kind}(?:[_\-]|$)", name):
            return kind

    return None

def _infer_pk_column(df: pd.DataFrame) -> str | None:
    """
    Heuristically infer a primary key-like column (for indexing) from DataFrame columns.
    Preference: columns ending with '_id', else 'id' variants.
    """
    for c in df.columns:
        if str(c).lower().endswith("_id"):
            return c
    for c in ("id", "ID", "Id"):
        if c in df.columns:
            return c
    return None

def _load_csv_clean(path: Path) -> pd.DataFrame:
    """
    Load delimited text robustly:
      - tries automatic delimiter detection (engine='python', sep=None)
      - falls back to common delimiters: ',', '\\t', '|', ';', whitespace
      - uses UTF-8 with BOM handling
      - skips malformed lines as last resort
    Then:
      - cast to pandas 'string' dtype where applicable
      - strip whitespace on string columns
      - drop exact duplicate rows
    """
    # Ordered attempts: (kwargs, description)
    attempts = [
        (dict(sep=None, engine="python"), "auto-detect"),
        (dict(sep=","), "comma"),
        (dict(sep="\t"), "tab"),
        (dict(sep="|"), "pipe"),
        (dict(sep=";"), "semicolon"),
        (dict(sep=r"\s+", engine="python"), "whitespace"),
    ]

    last_err = None
    for kwargs, label in attempts:
        try:
            df = pd.read_csv(path, dtype="string", encoding="utf-8-sig", **kwargs)
            # basic sanity: at least 1 column and 1 row
            if df.shape[1] == 0:
                raise ValueError("Parsed 0 columns")
            print(f"[READ] {path.name}: parsed with {label}")
            break
        except Exception as e:
            last_err = e
            continue
    else:
        # final fallback: skip bad lines
        try:
            df = pd.read_csv(
                path,
                dtype="string",
                encoding="utf-8-sig",
                sep=None,
                engine="python",
                on_bad_lines="skip",
            )
            print(f"[READ] {path.name}: parsed with auto-detect + skip bad lines")
        except Exception as e:
            raise RuntimeError(f"Failed to parse {path} — last error: {last_err}") from e

    # Normalize string columns
    for col in df.columns:
        if pd.api.types.is_string_dtype(df[col]):
            df[col] = df[col].str.strip()
    df = df.drop_duplicates().reset_index(drop=True)
    return df

def build_herb_info_tables(info_dir: str | Path = None,
                           table_name_map: dict[str, str] = None,
                           include_unknown: bool = False) -> dict[str, int]:
    """
    Scan `herb_info_dir` for CSVs and load them into SQLite as:
      - Herb_info, Formula_info, Target_info, Disease_info, Ingredient_info
    using filename-based inference (e.g., any CSV whose name contains 'herb' -> Herb_info).

    Parameters
    ----------
    info_dir : str | Path | None
        Directory containing entity CSV files. Defaults to global `herb_info_dir`.
    table_name_map : dict[str,str] | None
        Optional override of kind->table mapping.
    include_unknown : bool
        If True, CSVs with unknown kind will be loaded into tables using their stem as table name.

    Returns
    -------
    dict[str,int] : {table_name: row_count_loaded}
    """
    herb_info_dir = "./Data/HERB 2.0"

    if info_dir is None:
        info_dir = herb_info_dir
    info_dir = Path(info_dir)
    if table_name_map is None:
        table_name_map = TABLE_NAME_MAP

    if not info_dir.exists():
        raise FileNotFoundError(f"Info directory not found: {info_dir}")

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    # 1) 수집 단계: 파일 → 대상 테이블명으로 그룹핑
    files_by_table: dict[str, list[Path]] = {}

    csv_files = sorted(list(info_dir.glob("*.csv")) + list(info_dir.glob("*.txt")))
    if not csv_files:
        print(f"[WARN] No CSV/TXT files found in {info_dir}")

    for csv_path in csv_files:
        kind = _infer_kind_from_filename(csv_path)
        if kind is None:
            if not include_unknown:
                print(f"[SKIP] Unknown kind for file: {csv_path.name}")
                continue
            table_name = csv_path.stem  # fallback: 파일 이름 그대로
        else:
            table_name = table_name_map.get(kind)
            if not table_name:
                print(f"[SKIP] No table mapping for kind='{kind}' ({csv_path.name})")
                continue

        files_by_table.setdefault(table_name, []).append(csv_path)

    # 2) 적재 단계: 테이블 단위로 모두 합쳐서 1번만 write (replace)
    results: dict[str, int] = {}
    for table_name, paths in files_by_table.items():
        dfs = []
        for p in paths:
            df_part = _load_csv_clean(p)
            dfs.append(df_part)
            print(f"[LOAD] {p.name} → {table_name}: {len(df_part):,} rows")
        if not dfs:
            continue
        df = pd.concat(dfs, ignore_index=True)
        # 전역 중복 제거 (완전 동일 행)
        df = df.drop_duplicates().reset_index(drop=True)

        # 쓰기 (테이블당 1회 replace)
        df.to_sql(table_name, conn, if_exists="replace", index=False)
        results[table_name] = len(df)
        print(f"[DONE] Wrote {len(df):,} rows to {table_name} (from {len(paths)} files)")

        # 3) 인덱스 생성: *_id 우선, 없으면 id
        pk = _infer_pk_column(df)
        if pk:
            try:
                cur.execute(f'CREATE INDEX IF NOT EXISTS idx_{table_name}_{pk} ON "{table_name}"("{pk}")')
            except Exception as e:
                print(f"[WARN] Failed to create index on {table_name}({pk}): {e}")

    conn.commit()
    conn.close()
    if not results:
        print("[WARN] No tables were created. Check filenames and directory path.")
    return results
# build_herb_info_tables()


#2-2. Herb 2.0 Construction - edge

def build_herb_edges_db(filepath: str | Path | None = None, table_name: str = "HERB_edges") -> None:
    """
    HERB 2.0 edge CSV 파일을 읽어 SQLite DB에 저장.

    Parameters
    ----------
    filepath : str | Path | None
        HERB 2.0 edge CSV 경로 (예: "HERB 2.0_total_edges.csv").
        None일 경우 기본 경로 사용.
    table_name : str
        저장할 SQLite 테이블 이름 (default="HERB_edges")

    Returns
    -------
    None
    """
    if filepath is None:
        filepath = herb_info_path
    filepath = Path(filepath)
    df = pd.read_csv(filepath, dtype="string")
    print(f"[INFO] Loaded {len(df)} edges from {filepath}")

    conn = sqlite3.connect(db_path)
    df.to_sql(table_name, conn, if_exists="replace", index=False)
    conn.close()
    print(f"[DONE] Saved {len(df)} edges into table '{table_name}' in {db_path}")

    return

# build_herb_edges_db("./Data/HERB 2.0/api/HERB 2.0_total_edges.csv")  # 기본 경로 사용



# 2-3. Cleaning HERB_edges


from pathlib import Path
import sqlite3

def build_edges_clean_adme(db_path_local: str | Path,
                           out_table: str = "HERB_edges_clean",
                           drop_original: bool = False,
                           *,
                           use_adme: bool = False,
                           adme_OB: float = 30.0,
                           adme_DL: float = 0.18,
                           adme_MW: float = 500.0,
                           adme_HA: int = 10,
                           adme_HD: int = 5,
                           adme_LogP: float = 5.5,
                           adme_RB: int = 10) -> None:
    """
    HERB_edges를 정리해 {out_table} 생성:
      - (옵션) Ingredient에 ADME 필터 적용
      - Formula–Herb / Herb–Ingredient / Ingredient–Target / Target–Disease만 유지
      - 무방향 정규화 + 중복 제거
    """
    if db_path_local is None:
        raise ValueError("db_path_local을 반드시 지정하세요.")
    db_path_local = str(db_path_local)

    conn = sqlite3.connect(db_path_local)
    cur = conn.cursor()
    try:
        # 0) 타입 라벨링 뷰
        cur.execute("DROP VIEW IF EXISTS EdgeTyped")
        cur.execute("""
        CREATE VIEW EdgeTyped AS
        WITH e AS (SELECT Node_1, Node_2 FROM HERB_edges)
        SELECT
          e.Node_1, e.Node_2,
          h1.Herb_id       AS h1,  i1.Ingredient_id AS i1,
          t1.Target_id     AS t1,  t1.Ensembl_id    AS g1,
          d1.Disease_id    AS d1,  f1.Formula_id    AS f1,
          h2.Herb_id       AS h2,  i2.Ingredient_id AS i2,
          t2.Target_id     AS t2,  t2.Ensembl_id    AS g2,
          d2.Disease_id    AS d2,  f2.Formula_id    AS f2
        FROM e
        LEFT JOIN Herb_info        h1 ON e.Node_1 = h1.Herb_id
        LEFT JOIN Ingredient_info  i1 ON e.Node_1 = i1.Ingredient_id
        LEFT JOIN Target_info      t1 ON e.Node_1 = t1.Target_id
        LEFT JOIN Disease_info     d1 ON e.Node_1 = d1.Disease_id
        LEFT JOIN Formula_info     f1 ON e.Node_1 = f1.Formula_id
        LEFT JOIN Herb_info        h2 ON e.Node_2 = h2.Herb_id
        LEFT JOIN Ingredient_info  i2 ON e.Node_2 = i2.Ingredient_id
        LEFT JOIN Target_info      t2 ON e.Node_2 = t2.Target_id
        LEFT JOIN Disease_info     d2 ON e.Node_2 = d2.Disease_id
        LEFT JOIN Formula_info     f2 ON e.Node_2 = f2.Formula_id;
        """)
        conn.commit()

        # 1) (옵션) ADME 필터 뷰
        if use_adme:
            cur.execute("DROP VIEW IF EXISTS I_FILTER")
            cur.execute(f"""
            CREATE VIEW I_FILTER AS
            SELECT Ingredient_id
            FROM Ingredient_info
            WHERE
              CAST(NULLIF(OB_score,'')              AS REAL) >= {adme_OB}
              AND CAST(NULLIF(Drug_likeness,'')     AS REAL) >= {adme_DL}
              AND CAST(NULLIF(MolWt,'')             AS REAL) <= {adme_MW}
              AND CAST(NULLIF(NumHAcceptors,'')     AS REAL) <= {adme_HA}
              AND CAST(NULLIF(NumHDonors,'')        AS REAL) <= {adme_HD}
              AND CAST(NULLIF(MolLogP,'')           AS REAL) <= {adme_LogP}
              AND CAST(NULLIF(NumRotatableBonds,'') AS REAL) <= {adme_RB}
            """)
            conn.commit()

        # 2) KEEP 세트 (무방향 정규화 포함)
        cur.execute(f'DROP TABLE IF EXISTS "{out_table}_tmp"')
        hi_join      = "JOIN I_FILTER F ON F.Ingredient_id = i2" if use_adme else ""
        hi_join_rev  = "JOIN I_FILTER F ON F.Ingredient_id = i1" if use_adme else ""
        it_join_left = "JOIN I_FILTER F ON F.Ingredient_id = i1" if use_adme else ""
        it_join_right= "JOIN I_FILTER F ON F.Ingredient_id = i2" if use_adme else ""

        cur.execute(f"""
        CREATE TABLE "{out_table}_tmp" AS
        WITH
        FH AS (
          SELECT DISTINCT f1 AS A, h2 AS B FROM EdgeTyped WHERE f1 IS NOT NULL AND h2 IS NOT NULL
          UNION
          SELECT DISTINCT f2 AS A, h1 AS B FROM EdgeTyped WHERE f2 IS NOT NULL AND h1 IS NOT NULL
        ),
        HI AS (
          SELECT DISTINCT h1 AS A, i2 AS B
          FROM EdgeTyped {hi_join}
          WHERE h1 IS NOT NULL AND i2 IS NOT NULL
          UNION
          SELECT DISTINCT h2 AS A, i1 AS B
          FROM EdgeTyped {hi_join_rev}
          WHERE h2 IS NOT NULL AND i1 IS NOT NULL
        ),
        IT AS (
          SELECT DISTINCT i1 AS A, COALESCE(t2, g2) AS B
          FROM EdgeTyped {it_join_left}
          WHERE i1 IS NOT NULL AND (t2 IS NOT NULL OR g2 IS NOT NULL)
          UNION
          SELECT DISTINCT i2 AS A, COALESCE(t1, g1) AS B
          FROM EdgeTyped {it_join_right}
          WHERE i2 IS NOT NULL AND (t1 IS NOT NULL OR g1 IS NOT NULL)
        ),
        TD AS (
          SELECT DISTINCT COALESCE(t1, g1) AS A, d2 AS B
          FROM EdgeTyped
          WHERE (t1 IS NOT NULL OR g1 IS NOT NULL) AND d2 IS NOT NULL
          UNION
          SELECT DISTINCT COALESCE(t2, g2) AS A, d1 AS B
          FROM EdgeTyped
          WHERE (t2 IS NOT NULL OR g2 IS NOT NULL) AND d1 IS NOT NULL
        ),
        ALL_KEEP AS (
          SELECT 'Formula–Herb' AS tp, A, B FROM FH
          UNION ALL
          SELECT 'Herb–Ingredient' AS tp, A, B FROM HI
          UNION ALL
          SELECT 'Ingredient–Target' AS tp, A, B FROM IT
          UNION ALL
          SELECT 'Target–Disease' AS tp, A, B FROM TD
        )
        SELECT tp,
               CASE WHEN A <= B THEN A ELSE B END AS Node_1,
               CASE WHEN A <= B THEN B ELSE A END AS Node_2
        FROM ALL_KEEP
        WHERE A IS NOT NULL AND B IS NOT NULL AND A <> B
        """)
        conn.commit()

        # 3) 중복 제거 + 인덱스
        cur.execute(f'DROP TABLE IF EXISTS "{out_table}"')
        cur.execute(f'''
            CREATE TABLE "{out_table}" AS
            SELECT tp, Node_1, Node_2
            FROM "{out_table}_tmp"
            GROUP BY tp, Node_1, Node_2
        ''')
        cur.execute(f'DROP TABLE IF EXISTS "{out_table}_tmp"')

        cur.execute(f'CREATE INDEX IF NOT EXISTS ix_{out_table}_tp  ON "{out_table}"(tp)')
        cur.execute(f'CREATE INDEX IF NOT EXISTS ix_{out_table}_n1  ON "{out_table}"(Node_1)')
        cur.execute(f'CREATE INDEX IF NOT EXISTS ix_{out_table}_n2  ON "{out_table}"(Node_2)')
        cur.execute(f'CREATE UNIQUE INDEX IF NOT EXISTS ux_{out_table}_pair ON "{out_table}"(tp, Node_1, Node_2)')
        conn.commit()

        if drop_original:
            cur.execute('DROP TABLE IF EXISTS HERB_edges')
            conn.commit()

        print(f"[DONE] Built {out_table} ({'with' if use_adme else 'without'} ADME filters)")
    finally:
        conn.close()
