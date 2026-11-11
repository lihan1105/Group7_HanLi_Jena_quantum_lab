"""Verify original data integrity and compare reproduced fits with report values."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs")
args = parser.parse_args()
for line in (ROOT / "raw" / "SHA256SUMS.txt").read_text().splitlines():
    expected, name = line.split("  ", 1)
    actual = hashlib.sha256((ROOT / "raw" / name).read_bytes()).hexdigest()
    if actual != expected:
        raise ValueError(f"Raw file checksum mismatch: {name}")

reference = json.loads((ROOT / "reference" / "results.json").read_text())
actual = json.loads((args.output_dir / "results.json").read_text())
if len(reference) != len(actual):
    raise ValueError("Scan count mismatch")
for ref, result in zip(reference, actual):
    if ref.keys() != result.keys():
        raise ValueError("Result field mismatch")
    for key, expected in ref.items():
        value = result[key]
        if isinstance(expected, str):
            if value != expected:
                raise ValueError(f"Scan {ref['scan']}, {key}: mismatch")
        elif not np.allclose(value, expected, rtol=1e-5, atol=1e-10):
            raise ValueError(f"Scan {ref['scan']}, {key}: differs from report")

figures = [f"fit_{i:02}.png" for i in range(1, 11)] + [
    "glass.png", "waveplate.png", "visibility.png"
]
for name in figures:
    if not (args.output_dir / name).is_file():
        raise ValueError(f"Missing plot: {name}")
print("Verified: 11 original files, 10 fitted scans, and 13 plots.")
print("All numerical results agree with the report within the stated tolerances.")
