#!/usr/bin/env python
"""Download the T-REMD trialanine (Ala3) tutorial dataset from Google Drive.

Usage:
    python -m genepie.tests.download_tremd_data

This fetches ``remd_alat_tutorial.zip`` (~50 MB) and extracts it into the test
data directory as ``data/remd_alat_tutorial/``. The archive is a self-contained
subset of a Temperature-REMD simulation of trialanine (42-atom solute, 20
temperature states, 5000 frames each) plus the reference MBAR results, used by
the MBAR reweighting/resampling tutorial and regression test.

The data is NOT committed to the repository (it lives on Google Drive and is
downloaded on demand, the same way the chignolin integration data is).

Note: Requires the ``gdown`` package (pip install gdown) for reliable large
file downloads from Google Drive.
"""
# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))
    __package__ = "genepie.tests"
# --------------------------------------------

import zipfile

from .conftest import DATA_DIR, TREMD_DIR

# Google Drive file ID of remd_alat_tutorial.zip (see
# scripts/pack_remd_alat_tutorial.py for how the archive is built).
FILE_ID = "150tf8E7ONAhdNfckqwo9mA047izM03tw"
# The archive's single top-level directory. A README.md at its root signals a
# complete extraction.
TOP_LEVEL = "remd_alat_tutorial"


def is_available() -> bool:
    """Return True if the extracted tutorial data looks complete."""
    return (TREMD_DIR / "README.md").is_file()


def download() -> int:
    """Download and extract the T-REMD tutorial dataset."""
    if is_available():
        print(f"[SKIP] {TOP_LEVEL} (already extracted at {TREMD_DIR})")
        return 0

    try:
        import gdown
    except ImportError:
        print("Error: gdown package is required.")
        print("Install with: pip install gdown")
        return 1

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    zip_path = DATA_DIR / "remd_alat_tutorial.zip"

    url = f"https://drive.google.com/uc?id={FILE_ID}"
    print(f"[DOWNLOAD] {TOP_LEVEL}.zip from Google Drive...")
    try:
        gdown.download(url, str(zip_path), quiet=False)
    except Exception as e:
        print(f"  [ERROR] Failed to download: {e}")
        if zip_path.exists():
            zip_path.unlink()
        return 1

    # Guard against HTML error pages masquerading as the zip.
    if not zip_path.exists() or zip_path.stat().st_size < 10000:
        print("  [ERROR] Downloaded file is missing or too small; "
              "Google Drive may have returned an error page.")
        if zip_path.exists():
            zip_path.unlink()
        return 1

    print(f"[EXTRACT] -> {DATA_DIR}")
    try:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(DATA_DIR)
    except zipfile.BadZipFile:
        print("  [ERROR] Downloaded file is not a valid zip archive.")
        zip_path.unlink()
        return 1
    finally:
        if zip_path.exists():
            zip_path.unlink()

    if not is_available():
        print(f"  [ERROR] Extraction did not produce {TREMD_DIR / 'README.md'}.")
        return 1

    print()
    print(f"Download complete! Data available at: {TREMD_DIR}")
    return 0


def main():
    import sys
    sys.exit(download())


if __name__ == "__main__":
    main()
