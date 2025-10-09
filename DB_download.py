# download_db_gdrive.py
from pathlib import Path
import gdown

FILE_ID = "1MceeVsg9utgvFRjC73QPgJWPh9o9J_US"  # <- 너의 파일 ID
OUT_NAME = "DB.db"

DATA_DIR = Path("Data")
DATA_DIR.mkdir(exist_ok=True)
out_path = DATA_DIR / OUT_NAME

url = f"https://drive.google.com/uc?id={FILE_ID}"
gdown.download(url, str(out_path), quiet=False, fuzzy=True)

# 간단 무결성 체크: SQLite 헤더 검사
with open(out_path, "rb") as f:
    head = f.read(16)
assert head.startswith(b"SQLite format 3\0"), "다운로드된 파일이 SQLite DB로 보이지 않습니다."
print(f"✅ Saved: {out_path.resolve()}")
