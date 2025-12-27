#!/usr/bin/env python
"""Download test data for integration tests from Google Drive.

Usage:
    python -m genepie.tests.download_test_data

This script downloads the chignolin test data (PDB, PSF, DCD) from Google Drive.
These files are required for running the integration tests (test_integration.py).

Note: Requires gdown package (pip install gdown) for reliable large file downloads.
"""
# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))
    __package__ = "genepie.tests"
# --------------------------------------------

import pathlib

from .conftest import DATA_DIR

# Google Drive file IDs and filenames
FILES = [
    ("1WyFzvhuMjlwp2pNjga9B8RvTKoygBh-a", "chignolin.pdb"),
    ("1L1Y7YdSz46sTI1lQ7PoQJIqqbzM4F9Vh", "chignolin.psf"),
    ("1DZFUbCBVdCsfKzzrroIslre0eSctMaY-", "chignolin.dcd"),
]


def download():
    """Download chignolin test data from Google Drive."""
    try:
        import gdown
    except ImportError:
        print("Error: gdown package is required.")
        print("Install with: pip install gdown")
        return 1

    data_dir = DATA_DIR / "chignolin"
    data_dir.mkdir(parents=True, exist_ok=True)

    print(f"Downloading test data to: {data_dir}")
    print()

    for file_id, filename in FILES:
        dest = data_dir / filename
        if dest.exists():
            # Check if it's a valid file (not an HTML error page)
            if dest.stat().st_size > 10000 or not filename.endswith('.dcd'):
                print(f"[SKIP] {filename} (already exists)")
                continue
            else:
                # Remove corrupted file
                dest.unlink()

        url = f"https://drive.google.com/uc?id={file_id}"
        print(f"[DOWNLOAD] {filename}...")

        try:
            gdown.download(url, str(dest), quiet=False)
            print(f"  -> {dest}")
        except Exception as e:
            print(f"  [ERROR] Failed to download {filename}")
            print(f"          {e}")
            # Remove partial file if exists
            if dest.exists():
                dest.unlink()
            continue

    print()
    print("Download complete!")
    return 0


def main():
    import sys
    sys.exit(download())


if __name__ == "__main__":
    main()
