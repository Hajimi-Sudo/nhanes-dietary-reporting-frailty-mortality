#!/usr/bin/env python3
"""Download the declared NHANES source files and record SHA-256 hashes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import urllib.request
from datetime import date
from pathlib import Path


NHANES_BASE = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles"
MORTALITY_BASE = "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality"

FILES = [
    "DEMO_J.xpt", "DR1IFF_J.xpt", "DR1TOT_J.xpt", "BMX_J.xpt", "PFQ_J.xpt",
    "MCQ_J.xpt", "DPQ_J.xpt", "BPQ_J.xpt", "DIQ_J.xpt", "HUQ_J.xpt",
    "RXQ_RX_J.xpt", "KIQ_U_J.xpt", "GHB_J.xpt", "CBC_J.xpt", "PAQ_J.xpt",
    "SMQ_J.xpt", "ALQ_J.xpt",
]
MORTALITY_FILE = "NHANES_2017_2018_MORT_2019_PUBLIC.dat"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, target: Path) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": "survey-frailty-reproducibility/0.1"})
    with urllib.request.urlopen(request, timeout=120) as response, target.open("wb") as output:
        shutil.copyfileobj(response, output)
    if target.stat().st_size < 1000:
        raise RuntimeError(f"Downloaded file is unexpectedly small: {target}")
    prefix = target.read_bytes()[:256].lower()
    if b"<!doctype" in prefix or b"<html" in prefix:
        raise RuntimeError(f"Downloaded HTML instead of data: {target}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    raw_root = project_root / "data" / "raw" / "nhanes_2017_2018"
    log_root = project_root / "logs"
    raw_root.mkdir(parents=True, exist_ok=True)
    log_root.mkdir(parents=True, exist_ok=True)

    rows = []
    declared = [(name, f"{NHANES_BASE}/{name}") for name in FILES]
    declared.append((MORTALITY_FILE, f"{MORTALITY_BASE}/{MORTALITY_FILE}"))
    for name, url in declared:
        target = raw_root / name
        print(f"Downloading {name}")
        download(url, target)
        rows.append({
            "file": name,
            "source_url": url,
            "downloaded_at": date.today().isoformat(),
            "bytes": target.stat().st_size,
            "sha256": sha256(target),
        })

    manifest = log_root / "nhanes_download_manifest.csv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Downloaded {len(rows)} files")
    print(manifest)


if __name__ == "__main__":
    main()
