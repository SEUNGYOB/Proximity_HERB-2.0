from __future__ import annotations

import sqlite3
from pathlib import Path
import pandas as pd
import re

# Default SQLite DB path for this project
db_path = "./Data/DB.db"
def _get_table_columns(conn: sqlite3.Connection, table: str) -> set[str]:
    try:
        rows = conn.execute(f'PRAGMA table_info("{table}")').fetchall()
        return {r[1] for r in rows}
    except Exception:
        return set()

def get_ensp_ids_by_ingredient(
    ing_query: str,
    db_path_local: str | Path = None,
    distinct: bool = True,
    with_meta: bool = False,
    exact: bool = True,
    profile: str = "balanced",
    debug: bool = False
) -> list[str] | pd.DataFrame | dict:
    """
    Resolve an Ingredient identifier/name to Ingredient_id. Check ADME eligibility
    under a chosen profile ("balanced" or "lenient"). If eligible, traverse edges
    using **HERB_edges** (raw) to collect targets (HBTAR/Ensembl), map to ENSG,
    then to ENSP. If not eligible, return a report of failed criteria.

    Parameters
    ----------
    ing_query : str
        Ingredient_id (HBIN...), or ingredient name/alias/SMILES/formula.
    db_path_local : str | Path | None
        SQLite DB path.
    distinct : bool
        De-duplicate outputs where applicable.
    with_meta : bool
        If True, return DataFrame with columns [ensg, ensp, hgnc_symbol, uniprot].
    exact : bool
        Exact vs partial matching for ingredient lookup.
    profile : str
        "balanced" or "lenient" – ADME screen profile.
    debug : bool
        If True, return diagnostics dict (always), instead of mapping results when helpful.

    Returns
    -------
    list[str] | pd.DataFrame | dict
        ENSP list/DF, or a dict report when not eligible / when debug=True.
    """
    if db_path_local is None:
        db_path_local = db_path

    # --- helpers -----------------------------------------------------------
    def _fetch_adme(conn: sqlite3.Connection, ing_id: str) -> dict:
        fields = [
            "OB_score", "Drug_likeness", "MolWt", "NumHAcceptors",
            "NumHDonors", "MolLogP", "NumRotatableBonds"
        ]
        cols_info = _get_table_columns(conn, "Ingredient_info")
        if {"Ingredient_id"}.issubset(cols_info):
            select_cols = [f'"{c}"' if c in cols_info else 'NULL AS "'+c+'"' for c in fields]
            sql = f"""
            SELECT {', '.join(select_cols)}
            FROM Ingredient_info WHERE Ingredient_id = ?
            """
            row = pd.read_sql_query(sql, conn, params=[ing_id])
            if not row.empty:
                return {k: (None if pd.isna(row.iloc[0][k]) else row.iloc[0][k]) for k in fields}
        if _table_exists(conn, "Ingredient_adme"):
            cols_adme = _get_table_columns(conn, "Ingredient_adme")
            select_cols = [f'"{c}"' if c in cols_adme else 'NULL AS "'+c+'"' for c in fields]
            sql = f"""
            SELECT {', '.join(select_cols)} FROM Ingredient_adme WHERE Ingredient_id = ?
            """
            row = pd.read_sql_query(sql, conn, params=[ing_id])
            if not row.empty:
                return {k: (None if pd.isna(row.iloc[0][k]) else row.iloc[0][k]) for k in fields}
        return {k: None for k in fields}

    def _adme_check(metrics: dict, profile: str) -> tuple[bool, list[str], dict]:
        p = (profile or "balanced").lower().strip()
        if p not in {"balanced", "lenient"}:
            p = "balanced"
        ob_min   = 30.0 if p == "balanced" else 20.0
        dl_min   = 0.18 if p == "balanced" else 0.10
        logp_max = 5.5  if p == "balanced" else 6.0
        max_viol = 1    if p == "balanced" else 2
        OB   = metrics.get("OB_score")
        DL   = metrics.get("Drug_likeness")
        MW   = metrics.get("MolWt")
        HBA  = metrics.get("NumHAcceptors")
        HBD  = metrics.get("NumHDonors")
        LogP = metrics.get("MolLogP")
        fails = []
        if OB is not None and float(OB) < ob_min:
            fails.append(f"OB_score<{ob_min}")
        if DL is not None and float(DL) < dl_min:
            fails.append(f"Drug_likeness<{dl_min}")
        vio = 0
        if MW  is not None and float(MW)  > 500: vio += 1
        if HBA is not None and float(HBA) > 10:  vio += 1
        if HBD is not None and float(HBD) > 5:   vio += 1
        if LogP is not None and float(LogP) > logp_max: vio += 1
        if vio > max_viol:
            fails.append(f"Lipinski_violations>{max_viol}")
        eligible = (len(fails) == 0)
        details = {
            "profile": p,
            "thresholds": {
                "OB_min": ob_min, "DL_min": dl_min,
                "LogP_max": logp_max, "Lipinski_max_viol": max_viol
            },
            "metrics": metrics,
            "violations_count": vio
        }
        return eligible, fails, details

    # --- 1) Resolve ingredient(s) to Ingredient_id ------------------------
    if not ing_query:
        return [] if not with_meta else pd.DataFrame(columns=["ensg","ensp","hgnc_symbol","uniprot"])

    conn = sqlite3.connect(db_path_local)
    diagnostics = {"query": ing_query, "profile": profile}
    try:
        is_hbin = bool(re.fullmatch(r"(?i)HBIN\d+", (ing_query or "").strip()))
        params = {"q_exact": ing_query, "q_lower": str(ing_query).lower()}
        name_cols = ["Ingredient_name", "Ingredient_alias_name", "Canonical_smiles", "Isomeric_smiles", "Molecular_formula"]
        if is_hbin:
            sql_ing_resolve = """
            SELECT DISTINCT Ingredient_id FROM Ingredient_info
            WHERE Ingredient_id = :q_exact
            """
        else:
            conds = ['Ingredient_id = :q_exact']
            if exact:
                for c in name_cols:
                    conds.append(f'LOWER(TRIM(COALESCE("{c}", ""))) = :q_lower')
                sql_ing_resolve = f"""
                WITH cand AS (
                  SELECT Ingredient_id FROM Ingredient_info
                  WHERE {' OR '.join(conds)}
                )
                SELECT DISTINCT Ingredient_id FROM cand
                """
            else:
                params["q_like"] = f"%{params['q_lower']}%"
                for c in name_cols:
                    conds.append(f'LOWER(TRIM(COALESCE("{c}", ""))) LIKE :q_like')
                sql_ing_resolve = f"""
                WITH cand AS (
                  SELECT Ingredient_id FROM Ingredient_info
                  WHERE {' OR '.join(conds)}
                )
                SELECT DISTINCT Ingredient_id FROM cand
                """
        ing_ids = pd.read_sql_query(sql_ing_resolve, conn, params=params)["Ingredient_id"].dropna().astype(str).tolist()
        diagnostics["resolved_ing_ids"] = ing_ids
        if not ing_ids:
            report = {"eligible": False, "reason": "ingredient_not_found", "diagnostics": diagnostics}
            return report if debug else ([] if not with_meta else pd.DataFrame(columns=["ensg","ensp","hgnc_symbol","uniprot"]))

        # --- 2) ADME eligibility per ingredient ---------------------------
        adme_reports = {}
        for iid in ing_ids:
            metrics = _fetch_adme(conn, iid)
            eligible, fails, details = _adme_check(metrics, profile)
            adme_reports[iid] = {"eligible": eligible, "fails": fails, **details}
        diagnostics["adme_reports"] = adme_reports
        all_pass = all(r["eligible"] for r in adme_reports.values())
        if not all_pass:
            fail_report = {"eligible": False, "failed": adme_reports}
            if debug:
                return {**fail_report, "adme_reports": adme_reports, "diagnostics": diagnostics}
            return fail_report

        # --- 3) Single-shot SQL: ING→TARGET→ENSG→ENSP ---------------------
        edge_cols = _get_table_columns(conn, "HERB_edges")
        has_tp = ("tp" in edge_cols)
        placeholders = ",".join(["?"]*len(ing_ids))
        tp_clause = ("AND LOWER(e.tp) LIKE '%ingredient%target%'" if has_tp else "")

        if is_hbin:
            # Fast path: exact Ingredient_id matching on edges
            seeds_cte = f"""
            WITH seeds AS (
              SELECT Ingredient_id FROM Ingredient_info
              WHERE Ingredient_id IN ({placeholders})
            ),
            it_raw AS (
              SELECT CASE
                       WHEN e.Node_1 = s.Ingredient_id THEN e.Node_2
                       WHEN e.Node_2 = s.Ingredient_id THEN e.Node_1
                     END AS tgt
              FROM HERB_edges e
              JOIN seeds s ON e.Node_1 = s.Ingredient_id OR e.Node_2 = s.Ingredient_id
              {tp_clause}
            ),
            tg AS (
              SELECT DISTINCT ti.Ensembl_id AS ensg
              FROM it_raw r
              LEFT JOIN Target_info ti ON ti.Target_id = r.tgt OR ti.Ensembl_id = r.tgt
              WHERE ti.Ensembl_id IS NOT NULL
            )
            """
        else:
            # Name/alias tolerant path (case/space-insensitive)
            seeds_cte = f"""
            WITH seeds AS (
              SELECT Ingredient_id, Ingredient_name, Ingredient_alias_name
              FROM Ingredient_info
              WHERE Ingredient_id IN ({placeholders})
            ),
            ing_nodes AS (
              SELECT DISTINCT TRIM(Ingredient_id)   AS node FROM seeds WHERE TRIM(Ingredient_id)   <> ''
              UNION
              SELECT DISTINCT TRIM(Ingredient_name) AS node FROM seeds WHERE TRIM(Ingredient_name) <> ''
              UNION
              SELECT DISTINCT TRIM(Ingredient_alias_name) AS node FROM seeds WHERE TRIM(Ingredient_alias_name) <> ''
            ),
            it_raw AS (
              SELECT CASE
                       WHEN LOWER(TRIM(e.Node_1)) IN (SELECT LOWER(node) FROM ing_nodes) THEN e.Node_2
                       WHEN LOWER(TRIM(e.Node_2)) IN (SELECT LOWER(node) FROM ing_nodes) THEN e.Node_1
                     END AS tgt
              FROM HERB_edges e
              WHERE (LOWER(TRIM(e.Node_1)) IN (SELECT LOWER(node) FROM ing_nodes)
                     OR LOWER(TRIM(e.Node_2)) IN (SELECT LOWER(node) FROM ing_nodes))
                {tp_clause}
            ),
            tg AS (
              SELECT DISTINCT ti.Ensembl_id AS ensg
              FROM it_raw r
              LEFT JOIN Target_info ti ON ti.Target_id = r.tgt OR ti.Ensembl_id = r.tgt
              WHERE ti.Ensembl_id IS NOT NULL
            )
            """

        if with_meta:
            sql_final = f"""
            {seeds_cte}
            SELECT DISTINCT tg.ensg, x.ensp, x.hgnc_symbol, x.uniprot
            FROM tg
            JOIN Ensembl_xref x ON x.ensg = tg.ensg
            """
            df = pd.read_sql_query(sql_final, conn, params=ing_ids)
            if debug:
                cnt_sql = f"""
                {seeds_cte}
                SELECT COUNT(DISTINCT ensg) AS n FROM tg
                """
                ensg_count = int(pd.read_sql_query(cnt_sql, conn, params=ing_ids).iloc[0,0])
                diagnostics["ensg_count"] = ensg_count
            if df.empty:
                return pd.DataFrame(columns=["ensg","ensp","hgnc_symbol","uniprot"]) if with_meta else []
            if distinct:
                df = df.drop_duplicates(subset=["ensg","ensp"])  # safety
            if debug:
                return {
                    "eligible": True,
                    "targets": None,
                    "ensg": df["ensg"].nunique(),
                    "ensp": df["ensp"].nunique(),
                    "adme_reports": diagnostics.get("adme_reports"),
                    "diagnostics": diagnostics,
                    "result": df,
                }
            return df
        else:
            sql_final = f"""
            {seeds_cte}
            SELECT DISTINCT x.ensp
            FROM tg
            JOIN Ensembl_xref x ON x.ensg = tg.ensg
            """
            ensp_list = pd.read_sql_query(sql_final, conn, params=ing_ids)["ensp"].dropna().astype(str).tolist()
            if debug:
                cnt_sql = f"""
                {seeds_cte}
                SELECT COUNT(DISTINCT ensg) AS n FROM tg
                """
                ensg_count = int(pd.read_sql_query(cnt_sql, conn, params=ing_ids).iloc[0,0])
                diagnostics["ensg_count"] = ensg_count
                return {
                    "eligible": True,
                    "targets": None,
                    "ensg": ensg_count,
                    "ensp": len(set(ensp_list)),
                    "adme_reports": diagnostics.get("adme_reports"),
                    "diagnostics": diagnostics,
                    "result": ensp_list,
                }
            return list(dict.fromkeys(ensp_list)) if distinct else ensp_list

    finally:
        conn.close()

def _table_exists(conn: sqlite3.Connection, name: str) -> bool:
    try:
        return pd.read_sql_query(
            "SELECT 1 FROM sqlite_master WHERE type IN ('table','view') AND name=?",
            conn, params=(name,)
        ).shape[0] > 0
    except Exception:
        return False


def _normalize_ensembl_id(x: str | None) -> str | None:
    if x is None:
        return None
    s = str(x).strip()
    if not s:
        return None
    # strip version suffix like ENSP00000312345.2 → ENSP00000312345
    return s.split(".", 1)[0]

# === ENSG → ENSP helpers ===
from typing import Iterable, List, Dict, Any, Sequence

def ensg_list_to_ensp(ensgs: Sequence[str],
                      db_path_local: str | Path = None,
                      distinct: bool = True,
                      with_meta: bool = False) -> List[str] | pd.DataFrame:
    """
    Batch version: map a list of ENSG IDs to ENSP IDs.
    - Normalizes IDs (drops version suffix)
    - Chunks queries to respect SQLite variable limits

    Parameters
    ----------
    ensgs : sequence of str
        Input gene IDs.
    db_path_local : str | Path | None
        SQLite DB path. Defaults to module-level db_path.
    distinct : bool
        Drop duplicate ENSP per gene when True.
    with_meta : bool
        If True, return a DataFrame with columns [ensg, ensp, hgnc_symbol, uniprot].
        If False, return a flat list of ENSP strings.
    """
    if db_path_local is None:
        db_path_local = db_path
    if not ensgs:
        return [] if not with_meta else pd.DataFrame(columns=["ensg","ensp","hgnc_symbol","uniprot"])

    # Normalize and unique the inputs
    ensg_norm = []
    seen = set()
    for g in ensgs:
        n = _normalize_ensembl_id(g)
        if n and n not in seen:
            seen.add(n)
            ensg_norm.append(n)
    if not ensg_norm:
        return [] if not with_meta else pd.DataFrame(columns=["ensg","ensp","hgnc_symbol","uniprot"])

    # Chunk to avoid SQLite parameter limits (default ~999)
    chunk_size = 800
    frames = []
    conn = sqlite3.connect(db_path_local)
    try:
        # If mapping table/view is missing, return ENSG passthrough (no crash)
        if not _table_exists(conn, "Ensembl_xref"):
            # Fallback: return input ENSG as-is (optionally with meta columns)
            if with_meta:
                return pd.DataFrame({
                    "ensg": ensg_norm,
                    "ensp": [],
                    "hgnc_symbol": [],
                    "uniprot": []
                })
            return ensg_norm

        for i in range(0, len(ensg_norm), chunk_size):
            chunk = ensg_norm[i:i+chunk_size]
            placeholders = ",".join(["?"] * len(chunk))
            sql = f'''
                SELECT ensg, ensp, hgnc_symbol, uniprot
                FROM Ensembl_xref
                WHERE ensg IN ({placeholders}) AND ensp IS NOT NULL
            '''
            frames.append(pd.read_sql_query(sql, conn, params=chunk))
    finally:
        conn.close()
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["ensg","ensp","hgnc_symbol","uniprot"])
    if df.empty:
        return [] if not with_meta else df
    if distinct:
        df = df.drop_duplicates(subset=["ensg","ensp"]).reset_index(drop=True)
    if with_meta:
        return df
    return df["ensp"].dropna().astype(str).unique().tolist()

def _guess_col(df: pd.DataFrame, candidates) -> str | None:
    """
    Find a column in df matching any candidate name (case-insensitive, substring ok).
    """
    cols = list(df.columns)
    low = [c.lower() for c in cols]
    for cand in candidates:
        c = cand.lower()
        # exact case-insensitive
        if c in low:
            return cols[low.index(c)]
    # substring fallback
    for cand in candidates:
        c = cand.lower()
        for i, name in enumerate(low):
            if c in name:
                return cols[i]
    return None


def import_ensembl_xref_from_csv(
    csv_path: str | Path,
    db_path_local: str | Path = None,
    overwrite: bool = True
) -> dict:
    """
    Import an Ensembl cross-reference CSV (ENSG→ENSP [+ optional HGNC/UniProt])
    into SQLite table `Ensembl_xref`.

    The CSV can be downloaded from Ensembl BioMart. Typical column names include:
      - Ensembl Gene ID (ENSG)
      - Ensembl Protein ID (ENSP)
      - HGNC symbol
      - UniProtKB/Swiss-Prot ID

    Parameters
    ----------
    csv_path : str | Path
        Path to the CSV exported from BioMart.
    db_path_local : str | Path | None
        SQLite DB path. Defaults to module-level db_path.
    overwrite : bool
        When True, replace the table. When False, append new rows.

    Returns
    -------
    dict with summary counts: {'n_rows', 'n_valid', 'n_unique', 'table': 'Ensembl_xref'}
    """
    if db_path_local is None:
        db_path_local = db_path

    # 1) Load CSV
    df_raw = pd.read_csv(csv_path)
    if df_raw.empty:
        return {"n_rows": 0, "n_valid": 0, "n_unique": 0, "table": "Ensembl_xref"}

    # 2) Guess columns
    ensg_col = _guess_col(df_raw, ["Ensembl Gene ID", "Gene stable ID", "Gene ID", "ensg", "Ensembl_ID", "Ensembl_id"])
    ensp_col = _guess_col(df_raw, ["Ensembl Protein ID", "Protein stable ID", "Protein ID", "ensp"])
    sym_col  = _guess_col(df_raw, ["HGNC symbol", "Gene name", "HGNC Symbol", "hgnc_symbol", "Gene", "Symbol"])
    uni_col  = _guess_col(df_raw, ["UniProtKB/Swiss-Prot ID", "UniProtKB ID", "Uniprot ID", "uniprot", "UniProtKB accession"])

    if ensg_col is None or ensp_col is None:
        raise ValueError(f"Could not find required columns in CSV. Detected -> ENSG: {ensg_col}, ENSP: {ensp_col}")

    # 3) Normalize & filter
    df = pd.DataFrame({
        "ensg": df_raw[ensg_col].astype(str).str.strip().str.upper().map(_normalize_ensembl_id),
        "ensp": df_raw[ensp_col].astype(str).str.strip().str.upper().map(_normalize_ensembl_id),
    })
    if sym_col is not None:
        df["hgnc_symbol"] = df_raw[sym_col].astype(str).str.strip()
    else:
        df["hgnc_symbol"] = None
    if uni_col is not None:
        df["uniprot"] = df_raw[uni_col].astype(str).str.strip()
    else:
        df["uniprot"] = None

    # keep only plausible ENSG/ENSP
    df = df[df["ensg"].str.startswith("ENSG", na=False) & df["ensp"].str.startswith("ENSP", na=False)].copy()
    n_valid = len(df)

    # 4) Deduplicate
    df = df.drop_duplicates(subset=["ensg", "ensp"]).reset_index(drop=True)
    n_unique = len(df)

    # 5) Write to SQLite
    if n_unique == 0:
        return {"n_rows": len(df_raw), "n_valid": 0, "n_unique": 0, "table": "Ensembl_xref"}

    with sqlite3.connect(db_path_local) as conn:
        if overwrite:
            if _table_exists(conn, "Ensembl_xref"):
                try:
                    conn.execute('DROP TABLE IF EXISTS Ensembl_xref')
                except Exception:
                    pass
        df.to_sql("Ensembl_xref", conn, if_exists="replace" if overwrite else "append", index=False)

        # Indexes
        try:
            conn.execute("CREATE UNIQUE INDEX IF NOT EXISTS ux_Ensembl_xref_pair ON Ensembl_xref(ensg, ensp)")
        except Exception:
            pass
        try:
            conn.execute("CREATE INDEX IF NOT EXISTS idx_Ensembl_xref_ensg ON Ensembl_xref(ensg)")
        except Exception:
            pass
        try:
            conn.execute("CREATE INDEX IF NOT EXISTS idx_Ensembl_xref_ensp ON Ensembl_xref(ensp)")
        except Exception:
            pass

    return {"n_rows": len(df_raw), "n_valid": int(n_valid), "n_unique": int(n_unique), "table": "Ensembl_xref"}

def get_ensp_ids_by_herb(
    herb_query: str,
    db_path_local: str | Path = None,
    distinct: bool = True,
    with_meta: bool = False
) -> list[str] | pd.DataFrame:
    """
    Resolve a herb identifier/name to Herb_id, traverse edges Herb–Ingredient → Ingredient–Target,
    collect ENSG IDs, and convert them to ENSP IDs via `ensg_list_to_ensp`.

    Parameters
    ----------
    herb_query : str
        Any of the following that identifies the herb:
        [Herb_id, Herb_pinyin_name, Herb_cn_name, Herb_alias_name, Herb_en_name, Herb_latin_name]
        Matching for non-ID names is case-insensitive and exact (no partial matching).
    db_path_local : str | Path | None
        SQLite DB path. Defaults to module-level db_path.
    distinct : bool
        Drop duplicate mappings when True (applies both to ENSG aggregation and ENSP results).
    with_meta : bool
        If True, return a DataFrame with columns from `ensg_list_to_ensp` (ensg, ensp, hgnc_symbol, uniprot).
        If False, return a flat list of ENSP strings.

    Returns
    -------
    list[str] | pd.DataFrame
        ENSP list (default) or metadata DataFrame from `ensg_list_to_ensp(with_meta=True)`.
    """
    if db_path_local is None:
        db_path_local = db_path

    # 1) Resolve herb(s) → Herb_id
    #    We allow multiple matches if herb_query matches across multiple columns.
    herb_ids: list[str] = []
    conn = sqlite3.connect(db_path_local)
    # prefer cleaned edges if present, otherwise fall back to raw
    edge_table = "HERB_edges_clean" if _table_exists(conn, "HERB_edges_clean") else "HERB_edges"
    edge_cols = _get_table_columns(conn, edge_table)
    has_tp = ("tp" in edge_cols)
    try:
        params = {
            "q_exact": herb_query,
            "q_lower": str(herb_query).lower()
        }
        sql_herb = """
        WITH cand AS (
          SELECT Herb_id
          FROM Herb_info
          WHERE Herb_id = :q_exact
             OR LOWER(Herb_pinyin_name) = :q_lower
             OR LOWER(Herb_cn_name)     = :q_lower
             OR LOWER(Herb_alias_name)  = :q_lower
             OR LOWER(Herb_en_name)     = :q_lower
             OR LOWER(Herb_latin_name)  = :q_lower
        )
        SELECT DISTINCT Herb_id FROM cand
        """
        herb_ids = pd.read_sql_query(sql_herb, conn, params=params)["Herb_id"].dropna().astype(str).tolist()
        if not herb_ids:
            # Nothing to do
            return [] if not with_meta else pd.DataFrame(columns=["ensg","ensp","hgnc_symbol","uniprot"])

        # 2) Traverse H–I (role-agnostic) → get Ingredient_ids
        sql_ing = f"""
        WITH herbs(Herb_id) AS (
          SELECT DISTINCT Herb_id FROM Herb_info
          WHERE Herb_id IN ({{placeholders}})
        ),
        hi_raw AS (
          SELECT
            CASE
              WHEN e.Node_1 = h.Herb_id THEN e.Node_2
              WHEN e.Node_2 = h.Herb_id THEN e.Node_1
              ELSE NULL
            END AS maybe_ing
          FROM {edge_table} e
          JOIN herbs h
            ON (e.Node_1 = h.Herb_id OR e.Node_2 = h.Herb_id)
           { "AND LOWER(e.tp) LIKE '%herb%ingredient%'" if has_tp else "" }
        ),
        ingredients AS (
          SELECT DISTINCT i.Ingredient_id
          FROM hi_raw r
          JOIN Ingredient_info i ON i.Ingredient_id = r.maybe_ing
          WHERE r.maybe_ing IS NOT NULL
        )
        SELECT Ingredient_id FROM ingredients
        """
        placeholders = ",".join(["?"] * len(herb_ids))
        sql_ing = sql_ing.format(placeholders=placeholders)
        ing_ids = pd.read_sql_query(sql_ing, conn, params=herb_ids)["Ingredient_id"].dropna().astype(str).tolist()
        if not ing_ids:
            return [] if not with_meta else pd.DataFrame(columns=["ensg","ensp","hgnc_symbol","uniprot"])

        # 3) Traverse I–T (role-agnostic) → get Target nodes and resolve to Target_info (ENSG)
        sql_tg = f"""
        WITH ingredients(Ingredient_id) AS (
          SELECT DISTINCT Ingredient_id FROM Ingredient_info
          WHERE Ingredient_id IN ({{placeholders}})
        ),
        it_raw AS (
          SELECT
            CASE
              WHEN e.Node_1 = ing.Ingredient_id THEN e.Node_2
              WHEN e.Node_2 = ing.Ingredient_id THEN e.Node_1
              ELSE NULL
            END AS maybe_t
          FROM {edge_table} e
          JOIN ingredients ing
            ON (e.Node_1 = ing.Ingredient_id OR e.Node_2 = ing.Ingredient_id)
           { "AND LOWER(e.tp) LIKE '%ingredient%target%'" if has_tp else "" }
        ),
        targets AS (
          SELECT DISTINCT
            COALESCE(t.Ensembl_id, t.Target_id) AS any_id,
            t.Ensembl_id AS ensg
          FROM it_raw r
          LEFT JOIN Target_info t
            ON t.Target_id = r.maybe_t OR t.Ensembl_id = r.maybe_t
          WHERE r.maybe_t IS NOT NULL
        )
        SELECT DISTINCT ensg
        FROM targets
        WHERE ensg IS NOT NULL
        """
        placeholders_it = ",".join(["?"] * len(ing_ids))
        sql_tg = sql_tg.format(placeholders=placeholders_it)
        ensgs = pd.read_sql_query(sql_tg, conn, params=ing_ids)["ensg"].dropna().astype(str).tolist()

    finally:
        conn.close()

    if not ensgs:
        return [] if not with_meta else pd.DataFrame(columns=["ensg","ensp","hgnc_symbol","uniprot"])

    # Normalize ENSG IDs and map to ENSP via existing helper
    ensgs_norm = []
    seen = set()
    for g in ensgs:
        n = _normalize_ensembl_id(g)
        if n and n not in seen:
            seen.add(n)
            ensgs_norm.append(n)

    return ensg_list_to_ensp(
        ensgs_norm,
        db_path_local=db_path_local,
        distinct=distinct,
        with_meta=with_meta
    )


def get_ensp_ids_by_disease(
    dis_query: str,
    db_path_local: str | Path = None,
    distinct: bool = True,
    with_meta: bool = False,
    exact: bool = True,
    debug: bool = False,
) -> list[str] | pd.DataFrame | dict:
    """
    Resolve a disease identifier to **Disease_id (HBDIS...)** from Disease_info using
    DisGeNET_id → Disease_id, then via HERB_edges infer connected targets and map to
    ENSP (or return meta).

    Flow
    ----
    1) Resolve input:
       - If input matches `HBDIS/d+`, treat as Disease_id.
       - Else, resolve `DisGeNET_id` → `Disease_id` (exact or LIKE per `exact`).
    2) SQL one-shot traversal (no Python round-trips):
       seeds(Disease_id) → HERB_edges(Node_1/2) → Target_info(Target_id/Ensembl_id)
       → tg(ENSG) → Ensembl_xref → ENSP

    Parameters
    ----------
    dis_query : str
        Disease_id (HBDIS...) or DisGeNET_id.
    db_path_local : str | Path | None
        SQLite DB path. Defaults to module-level `db_path`.
    distinct : bool
        De-duplicate outputs where applicable.
    with_meta : bool
        If True, return DataFrame with columns [ensg, ensp, hgnc_symbol, uniprot].
    exact : bool
        If False and `dis_query` is not HBDIS, perform LIKE on DisGeNET_id.
    debug : bool
        If True, return diagnostics dict (counts + resolved IDs + result).
    """
    if db_path_local is None:
        db_path_local = db_path

    q = (dis_query or "").strip()
    if not q:
        return [] if not with_meta else pd.DataFrame(columns=["ensg","ensp","hgnc_symbol","uniprot"])

    conn = sqlite3.connect(db_path_local)
    diagnostics = {"query": dis_query}
    try:
        # 1) Resolve to Disease_id (HBDIS...)
        is_hbdis = bool(re.fullmatch(r"(?i)HBDIS\d+", q))
        if is_hbdis:
            sql_resolve = """
            SELECT DISTINCT Disease_id FROM Disease_info WHERE Disease_id = :q
            """
            disease_ids = pd.read_sql_query(sql_resolve, conn, params={"q": q})["Disease_id"].dropna().astype(str).tolist()
        else:
            if exact:
                sql_resolve = """
                SELECT DISTINCT Disease_id FROM Disease_info
                WHERE DisGeNET_id = :q
                """
                disease_ids = pd.read_sql_query(sql_resolve, conn, params={"q": q})["Disease_id"].dropna().astype(str).tolist()
            else:
                sql_resolve = """
                SELECT DISTINCT Disease_id FROM Disease_info
                WHERE LOWER(TRIM(COALESCE(DisGeNET_id,''))) LIKE :q_like
                """
                disease_ids = pd.read_sql_query(sql_resolve, conn, params={"q_like": f"%{q.lower()}%"})["Disease_id"].dropna().astype(str).tolist()
        diagnostics["resolved_disease_ids"] = disease_ids
        if not disease_ids:
            rep = {"ok": False, "reason": "disease_not_found", "diagnostics": diagnostics}
            return rep if debug else ([] if not with_meta else pd.DataFrame(columns=["ensg","ensp","hgnc_symbol","uniprot"]))

        # 2) Single-shot SQL: DISEASE→TARGET→ENSG→ENSP
        placeholders = ",".join(["?"] * len(disease_ids))
        seeds_cte = f"""
        WITH seeds AS (
          SELECT Disease_id FROM Disease_info
          WHERE Disease_id IN ({placeholders})
        ),
        dt_raw AS (
          SELECT CASE
                   WHEN e.Node_1 = s.Disease_id THEN e.Node_2
                   WHEN e.Node_2 = s.Disease_id THEN e.Node_1
                 END AS tgt
          FROM HERB_edges e
          JOIN seeds s ON e.Node_1 = s.Disease_id OR e.Node_2 = s.Disease_id
        ),
        tg AS (
          SELECT DISTINCT ti.Ensembl_id AS ensg
          FROM dt_raw r
          LEFT JOIN Target_info ti ON ti.Target_id = r.tgt OR ti.Ensembl_id = r.tgt
          WHERE ti.Ensembl_id IS NOT NULL
        )
        """

        if with_meta:
            sql_final = f"""
            {seeds_cte}
            SELECT DISTINCT tg.ensg, x.ensp, x.hgnc_symbol, x.uniprot
            FROM tg
            JOIN Ensembl_xref x ON x.ensg = tg.ensg
            """
            df = pd.read_sql_query(sql_final, conn, params=disease_ids)
            if df.empty:
                return pd.DataFrame(columns=["ensg","ensp","hgnc_symbol","uniprot"]) if with_meta else []
            if distinct:
                df = df.drop_duplicates(subset=["ensg","ensp"])  # safety
            if debug:
                cnt = pd.read_sql_query(f"""
                    {seeds_cte}
                    SELECT COUNT(DISTINCT ensg) AS n FROM tg
                """, conn, params=disease_ids).iloc[0,0]
                diagnostics["ensg_count"] = int(cnt)
                return {"ok": True, "ensg": int(cnt), "ensp": df["ensp"].nunique(), "result": df, "diagnostics": diagnostics}
            return df
        else:
            sql_final = f"""
            {seeds_cte}
            SELECT DISTINCT x.ensp
            FROM tg
            JOIN Ensembl_xref x ON x.ensg = tg.ensg
            """
            ensp_list = pd.read_sql_query(sql_final, conn, params=disease_ids)["ensp"].dropna().astype(str).tolist()
            if debug:
                cnt = pd.read_sql_query(f"""
                    {seeds_cte}
                    SELECT COUNT(DISTINCT ensg) AS n FROM tg
                """, conn, params=disease_ids).iloc[0,0]
                diagnostics["ensg_count"] = int(cnt)
                return {"ok": True, "ensg": int(cnt), "ensp": len(set(ensp_list)), "result": ensp_list, "diagnostics": diagnostics}
            return list(dict.fromkeys(ensp_list)) if distinct else ensp_list

    finally:
        conn.close()

def get_ensp_ids(
    query: str,
    db_path_local: str | Path = None,
    distinct: bool = True,
    with_meta: bool = False,
    exact: bool = True,
    profile: str = "balanced",   # passed to ingredient path
    debug: bool = False           # forwarded to ingredient/disease paths
) -> list[str] | pd.DataFrame | dict:
    """
    Prefix router with early-returns.

    Branching rules
    ---------------
    - **HERB***  → `get_ensp_ids_by_herb`
    - **HBIN***  → `get_ensp_ids_by_ingredient` (passes `profile`, `exact`, `debug`)
    - **HBDIS*** → `get_ensp_ids_by_disease`    (passes `exact`, `debug`)
    - otherwise → empty result (list or DataFrame depending on `with_meta`)
    """
    if db_path_local is None:
        db_path_local = db_path

    q = (query or "").strip()
    if not q:
        return [] if not with_meta else pd.DataFrame(columns=["ensg","ensp","hgnc_symbol","uniprot"])

    q_upper = q.upper()

    # HERB → Herb path
    if q_upper.startswith("HERB"):
        return get_ensp_ids_by_herb(
            q,
            db_path_local=db_path_local,
            distinct=distinct,
            with_meta=with_meta,
        )

    # HBIN → Ingredient path (fast path + ADME profile supported)
    if q_upper.startswith("HBIN"):
        return get_ensp_ids_by_ingredient(
            q,
            db_path_local=db_path_local,
            distinct=distinct,
            with_meta=with_meta,
            exact=exact,
            profile=profile,
            debug=debug,
        )

    # HBDIS → Disease path (DisGeNET_id → Disease_id → targets)
    if q_upper.startswith("HBDIS"):
        return get_ensp_ids_by_disease(
            q,
            db_path_local=db_path_local,
            distinct=distinct,
            with_meta=with_meta,
            exact=exact,
            debug=debug,
        )

    # Fallback: unsupported prefix
    return [] if not with_meta else pd.DataFrame(columns=["ensg","ensp","hgnc_symbol","uniprot"])


