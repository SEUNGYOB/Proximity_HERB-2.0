# 🧬 Proximity-Final

Network-based **Herbal Medicine – Target – Disease** relationship analysis using **STRING Human PPI** and **HERB 2.0 (HIT)** database.  
Implements efficient **Network Proximity** (A–B set distance & Z-score) calculation with **incremental memmap caching**.

---

## 📁 1. Database Download

### 🔹 1.1 Human PPI (STRING)
- **Source:** [STRING v12.0](https://string-db.org/cgi/download?sessionId=bxgqh7YBMvBD&species_text=Homo+sapiens&settings_expanded=0&min_download_score=0&filter_redundant_pairs=0&delimiter_type=txt)
- **File:** `9606.protein.links.v12.0.txt.gz`
- **Processed with:** `build_ppi_db()`
- **Uploaded Table:** `Human_PPI`
  - Total edges: **473,860**
  - Node format: `9606.ENSPXXXX` → prefix `9606.` removed during DB upload

---

### 🔹 1.2 HERB 2.0 (HIT) Data
- **Source:** [HERB 2.0 Download Portal](http://47.92.70.12/Download/)
- **Categories:**
  | Type | Range | Example ID | Note |
  |------|--------|-------------|------|
  | Herb | `HERB000001` ~ `HERB007263` | — |  |
  | Ingredient | `HBIN000001` ~ `HBIN049258` | — | ADME filtering applied |
  | Formula | `HBFO000001` ~ `HBFO006743` | — |  |
  | Target | `HBTAR000001` ~ `HBTAR015515` | ENSG format |  |
  | Disease | `HBDIS000001` ~ `HBDIS030170` | — |  |
  | Meta-analysis | `HBMA000001` | — |  |
  | Clinical trial | `HBCT000001` | — |  |
  | Reference | `HBREF000001` | — |  |
  | Experiment | `HBEXP000002` | — |  |

---

## 🧩 2. Database Construction

**Tables**
- `Human_PPI`
- `Herb_2_0_edge`
- `Herb_2_0_node`

**Filtering**
- Retained edges:  
  `Formula–Herb–Ingredient–Target (F–H–I–T)` and `Disease–Target (D–T)`
- Only **ADME-passing ingredients** are kept.

### ⚗️ ADME Screening Criteria (Balanced Profile)

| Property | Threshold | Rule |
|-----------|------------|------|
| OB_score | ≥ 30 | None = neutral |
| Drug_likeness | ≥ 0.18 |  |
| MolWt | ≤ 500 |  |
| NumHAcceptors | ≤ 10 |  |
| NumHDonors | ≤ 5 |  |
| MolLogP | ≤ 5.5 |  |
| NumRotatableBonds | ≤ 10 |  |

---
## 📊 Database Information (DB.db)`

SQLite `DB.db`의 주요 테이블을 워크북 시트로 내보낸 파일입니다.  
각 시트의 출처·의미·주요 식별자를 아래에 정리했습니다. *(괄호 안은 현재 행 수 예시)*

| Table | Source | What it contains | Key IDs / Notes |
|---|---|---|---|
| **BIOMART_mapping_raw** (533,785) | Ensembl **BioMart** | BioMart에서 그대로 받은 원시 매핑 | `Gene stable ID(ENSG)`, `Protein stable ID(ENSP)`, `Ensembl Canonical`, `APPRIS`, `Gene name` |
| **Ensembl_xref** (23,869) | Derived (from BioMart) | **대표 1:1 매핑** (ENSG→ENSP, Canonical/APPRIS 우선 규칙 적용) | `ensg`, `ensp`, `hgnc_symbol`, `uniprot`, `picked_by` |
| **Ensembl_xref_all** (245,535) | Derived (from BioMart) | **전체 1:N 매핑** (모든 ENSG–ENSP 쌍) | `ensg`, `ensp`, `hgnc_symbol`, `uniprot` |
| **Human_PPI** (473,860) | **STRING v12.0** | 인간 단백질–단백질 상호작용 엣지 | `protein1`, `protein2` (형태는 `ENSP`; 업로드 시 `9606.` 접두사 제거) |
| **Disease_info** (30,170) | HERB 2.0 | 질병 메타데이터 | `Disease_id (HBDIS...)`, `Disease_name`, `DisGeNET_id`, 분류 체계/식별자 |
| **Formula_info** (6,743) | HERB 2.0 | 처방(방제) 메타데이터 | `Formula_id (HBFO...)`, 한/영/병음명, 범주/출전 |
| **Herb_info** (6,892) | HERB 2.0 | 본초(약재) 메타데이터 | `Herb_id (HERB...)`, 한/영/병음/이명 |
| **Ingredient_info** (44,595) | HERB 2.0 | 성분(화합물) 메타데이터 + **ADME 필드** | `Ingredient_id (HBIN...)`, `Canonical_smiles`, `Molecular_formula`, `OB_score`, `Drug_likeness`, `MolWt`, `NumHAcceptors/Donors`, `MolLogP`, `NumRotatableBonds` |
| **Target_info** (15,515) | HERB 2.0 | 타깃 메타데이터 | `Target_id (HBTAR...)`, `Ensembl_id (ENSG)` |
| **HERB_edges** (859,757) | **Scraped from HERB 2.0** | HERB 2.0에서 수집한 **원시 네트워크 엣지** | `Node_1`, `Node_2` (노드에는 `HERB/HBIN/HBTAR/HBDIS/HBFO` 등 코드가 혼재) |

> 🧭 **파이프라인 요약**  
> - **BioMart**로 ENSG→ENSP 매핑을 만들고(`Ensembl_xref`, `Ensembl_xref_all`),  
> - **HERB_edges**의 Target(주로 ENSG)을 **ENSP**로 연결해 **STRING Human_PPI** 네트워크에서 **Network Proximity**를 계산합니다.  
> - PPI 노드 표기는 `ENSP`; STRING 원본의 `9606.ENSPxxxxx`는 업로드 시 `9606.` 접두사를 제거했습니다.


## 🧬 3. BioMart (ENSG → ENSP Mapping)

- **Portal:** [Ensembl BioMart](https://www.ensembl.org/biomart/martview/8403dac70986fcef75d758c4bc648d6f)
- **Dataset:** *Ensembl Genes 115, Human (GRCh38.p14)*
- **Columns Used:**
  - Gene stable ID (ENSG)
  - Protein stable ID (ENSP)
  - APPRIS annotation
  - Ensembl Canonical
  - Gene name
- **Output:** `ENSG → ENSP` mapping table  
  Used to link HERB Target (ENSG) to STRING PPI (ENSP).

---

## 🔗 4. Network Proximity Calculation (Example)

아래 예시는 DB에서 STRING PPI를 불러와 그래프를 만들고,
HERB/DISEASE에 해당하는 **ENSP 노드 집합**을 구한 뒤
`compute_network_distances_CPU()`로 **근접도(Z-score)** 를 계산합니다.

> 캐싱 전략  
> - 빠르게 시작하려면 `dist=None` → **증분(memmap) 캐시**(`./Data/Human_PPI/rows_mm.npy`) 자동 사용  
> - 대규모 반복 실험이면 `build_dist_matrix_from_graph()`로 **풀 매트릭스**를 만들어 `dist=`에 넘기면 더 빠름

```python
# example: examples/proximity_demo.py
import sqlite3
import pandas as pd
import networkx as nx
import numpy as np

import Utils
import proximity_util as pu
from Execution_Proximity import load_edges_from_db, create_network_from_edges
# (필요시) from Execution_Proximity import build_dist_matrix_from_graph

DB_PATH = "./Data/DB.db"

# 1) Human PPI 로드 → Graph
edges = load_edges_from_db(DB_PATH, table_name="Human_PPI",
                           source_col="protein1", target_col="protein2")
G = create_network_from_edges(edges)

# 2) 쿼리 → ENSP 노드 집합
herb_ensp = Utils.get_ensp_ids("HERB002168")
ing_ensp  = Utils.get_ensp_ids("HBIN046526")
dis_ensp  = Utils.get_ensp_ids("HBDIS001345")

print(len(herb_ensp), len(ing_ensp), len(dis_ensp))  # sanity check

# 3A) 증분 캐싱 경로 (권장: 처음 사용할 때)
res = pu.compute_network_distances_CPU(
    G=G,
    dist=None,              # memmap 캐시(./Data/Human_PPI/rows_mm.npy)에 행을 점진적으로 쌓음
    A=herb_ensp,
    B=dis_ensp,
    random_time=100,        # degree-matched resampling 횟수
    seed=42,                # 재현성(다양성 원하면 None)
    max_workers=16
)
print(res)
# 출력 예:
# {
#   'shortest': 1.7488,
#   'Z_score': {'d': 1.7488, 'z': -7.79, 'mean': 2.0093, 'std': 0.0334, 'p': 0.0}
# }

# 3B) (선택) 풀 매트릭스 경로: 반복 실험/배치에 유리
# dist_full = build_dist_matrix_from_graph(G)
# res = pu.compute_network_distances_CPU(G=G, dist=dist_full, A=herb_ensp, B=dis_ensp, random_time=100, seed=42)
# print(res)



### 🧠 Output Format

```python
{
  "shortest": <observed mean distance>,
  "Z_score": {
      "d": <observed>,
      "z": <Z-score>,
      "mean": <random mean>,
      "std": <random std>,
      "p": <empirical p-value>
  }
}
