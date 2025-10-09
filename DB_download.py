#!/usr/bin/env python3
"""
download_db_gdrive.py
---------------------
Google Drive에서 SQLite DB 파일(DB.db)을 자동으로 다운로드합니다.
"""

from pathlib import Path
import gdown
import sys

FILE_ID = "1MceeVsg9utgvFRjC73QPgJWPh9o9J_US"  # ← 네 DB 파일 ID
OUT_NAME = "DB.db"

def main():
    data_dir = Path("Data")
    data_dir.mkdir(exist_ok=True)
    out_path = data_dir / OUT_NAME

    url = f"https://drive.google.com/uc?id={FILE_ID}"
    print(f"📥 Downloading from Google Drive...\nURL: {url}")

    try:
        gdown.download(url, str(out_path), quiet=False, fuzzy=True)
    except Exception as e:
        print(f"❌ 다운로드 실패: {e}")
        sys.exit(1)

    # 무결성 검사 (SQLite 헤더 확인)
    try:
        with open(out_path, "rb") as f:
            head = f.read(16)
        assert head.startswith(b"SQLite format 3\0"), "다운로드된 파일이 SQLite DB로 보이지 않습니다."
        print(f"✅ 성공적으로 저장되었습니다: {out_path.resolve()}")
    except Exception as e:
        print(f"⚠️ 파일 무결성 오류: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
