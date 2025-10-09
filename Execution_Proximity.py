import Utils
import proximity_util as pu
import sqlite3
import pandas as pd
import networkx as nx
import numpy as np
import os
import json

db_path = "./Data/DB.db"



def load_edges_from_db(db_path, table_name="Human_PPI", source_col="protein1", target_col="protein2"):
    """
    데이터베이스의 Human_PPI 테이블에서 엣지를 로딩합니다.
    - 노드 ID가 '9606.ENSP000...' 형식이면 '9606.' 접두사를 제거합니다.
    """
    conn = sqlite3.connect(db_path)
    query = f"SELECT {source_col}, {target_col} FROM {table_name}"
    df = pd.read_sql_query(query, conn)
    conn.close()

    # 접두사 제거
    def strip_prefix(x):
        if isinstance(x, str) and x.startswith("9606."):
            return x.split("9606.")[-1]
        return x

    df[source_col] = df[source_col].map(strip_prefix)
    df[target_col] = df[target_col].map(strip_prefix)

    edges = list(df.itertuples(index=False, name=None))
    print(f"✅ '{table_name}' 테이블에서 {len(edges)}개의 엣지를 가져왔습니다. (9606. 접두사 제거 완료)")
    return edges
def create_network_from_edges(edges):
    """
    로딩된 엣지를 기반으로 networkx Graph 객체를 생성합니다.
    """
    G = nx.Graph ()
    G.add_edges_from (edges)
    print (f"✅ 네트워크 생성 완료: {len (G.nodes ())}개의 노드와 {len (G.edges ())}개의 엣지 추가됨.")
    return G
def build_dist_matrix_from_graph(G):
    """
    그래프 G로부터 DistMatrix 객체를 생성합니다.
    - 모든 노드 쌍의 최단 거리 계산
    - name → index 매핑 저장
    - compute_network_distances_CPU() 에서 dist 인자로 사용 가능
    캐시가 존재하면 캐시에서 로드하고, 없으면 계산 후 저장합니다.
    """
    cache_dir = "./Data/Human_PPI"
    dist_matrix_path = os.path.join(cache_dir, "dist_matrix.npy")
    node_index_path = os.path.join(cache_dir, "node_index.json")

    if os.path.exists(dist_matrix_path) and os.path.exists(node_index_path):
        mat = np.load(dist_matrix_path)
        with open(node_index_path, "r") as f:
            name2idx = json.load(f)
        print("✅ 거리 행렬을 캐시에서 불러왔습니다.")
    else:
        nodes = list(G.nodes())
        name2idx = {n: i for i, n in enumerate(nodes)}
        N = len(nodes)
        mat = np.full((N, N), np.inf, dtype=float)
        np.fill_diagonal(mat, 0.0)

        for src, lengths in nx.all_pairs_shortest_path_length(G):
            i = name2idx[src]
            for dst, d in lengths.items():
                mat[i, name2idx[dst]] = float(d)

        os.makedirs(cache_dir, exist_ok=True)
        np.save(dist_matrix_path, mat)
        with open(node_index_path, "w") as f:
            json.dump(name2idx, f)
        print("✅ 거리 행렬을 새로 계산하여 캐시에 저장했습니다.")

    class DistMatrix:
        def __init__(self, mat, idx):
            self.mat = mat
            self.idx = idx

    return DistMatrix(mat, name2idx)



Whole_G = create_network_from_edges (load_edges_from_db (db_path))

Herb_node_lists = Utils.get_ensp_ids("HERB002168")
Ingredient_node_lists = Utils.get_ensp_ids("HBIN046526")
Disease_node_lists = Utils.get_ensp_ids("HBDIS001345")

dist = pu.compute_network_distances_CPU(
    G=Whole_G,
    dist=None,                  # 점진 캐싱 경로
    A=Herb_node_lists,
    B=Disease_node_lists,
    random_time=1,
    max_workers=16
)

