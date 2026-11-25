"""Check saved hydraulic predictions after excluding overlapping symbols.

Reads fixed observation/prediction CSVs only. It does not call a solver, fit
parameters, change observations, or treat digitization bounds as experimental
measurement uncertainty. Output is alem/results/alem_nonoverlap_check.json.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / 'alem' / 'results'


def statistics(rows):
    if not rows:
        raise ValueError('Cannot score an empty observation subset.')
    residuals = [float(row['residual']) for row in rows]
    return {
        'n': len(rows),
        'point_ids': [int(row['point_id']) for row in rows],
        'PV': [float(row['pore_volumes']) for row in rows],
        'RMSE': math.sqrt(sum(value * value for value in residuals) / len(rows)),
        'MAE': sum(abs(value) for value in residuals) / len(rows),
        'bias': sum(residuals) / len(rows),
        'minimum_absolute_error': min(abs(value) for value in residuals),
        'maximum_absolute_error': max(abs(value) for value in residuals),
        'digitization_value_bounds': [
            float(row['digitization_value_bound']) for row in rows
        ],
    }


def main():
    output = {
        'scope': (
            'Sensitivity to exclusion of overlapping constant-flow Alem Figure 6a '
            'symbols, using saved N=240 predictions without refitting.'
        ),
        'observation': 'K_over_K0',
        'residual_definition': 'prediction minus observation',
        'digitization_bounds_are_not_experimental_uncertainties': True,
        'models': {},
    }
    for model in ('linear', 'multistage'):
        source = RESULTS / f'{model}_observation_predictions.csv'
        source_bytes = source.read_bytes()
        with source.open(newline='') as handle:
            rows = [
                row for row in csv.DictReader(handle)
                if row['mode'] == 'constant_flow' and row['quantity'] == 'K_over_K0'
            ]
        if len(rows) != 7:
            raise ValueError(f'{model}: expected seven constant-flow hydraulic symbols.')
        unambiguous = [row for row in rows if row['overlapping_series'].lower() == 'false']
        overlapping = [row for row in rows if row['overlapping_series'].lower() == 'true']
        if len(unambiguous) != 3 or len(overlapping) != 4:
            raise ValueError(f'{model}: expected three unambiguous and four overlapping symbols.')
        for row in rows:
            reconstructed = float(row['predicted']) - float(row['value'])
            if not math.isclose(reconstructed, float(row['residual']), rel_tol=1e-12, abs_tol=1e-14):
                raise ValueError(f'{model}: saved residual does not match prediction minus observation.')
        output['models'][model] = {
            'all_seven': statistics(rows),
            'unambiguous_three': statistics(unambiguous),
            'overlapping_four': statistics(overlapping),
            'source_csv': str(source.relative_to(ROOT)),
            'source_sha256': hashlib.sha256(source_bytes).hexdigest(),
        }
        if source.read_bytes() != source_bytes:
            raise RuntimeError(f'{model}: source CSV changed during the check.')
    destination = RESULTS / 'alem_nonoverlap_check.json'
    destination.write_text(json.dumps(output, indent=2) + '\n')
    print(destination)


if __name__ == '__main__':
    main()
