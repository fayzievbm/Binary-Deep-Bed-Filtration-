#!/usr/bin/env python3
"""Check the integrity of the recorded research code, figures, and numerical data."""
from pathlib import Path
import hashlib
import json
import sys

root = Path(__file__).resolve().parents[1]
manifest = json.loads((root / 'SHA256SUMS.json').read_text(encoding='utf-8'))
failed = []
for relative, expected in manifest.items():
    path = root / relative
    if not path.is_file():
        failed.append(f'MISSING: {relative}')
    elif hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        failed.append(f'CHANGED: {relative}')
if failed:
    print('\n'.join(failed))
    print(f'Integrity check failed for {len(failed)} of {len(manifest)} files.')
    sys.exit(1)
print(f'Verified {len(manifest)} files; every SHA-256 checksum matches.')
