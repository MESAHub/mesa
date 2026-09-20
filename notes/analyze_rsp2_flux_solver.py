"""Summarize RSP2_report_flux_solver output without compiling or running MESA."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from statistics import median


def distribution(values):
    values = sorted(values)
    if not values:
        return {}
    return dict(count=len(values), minimum=values[0], median=median(values),
                p90=values[int(0.9*(len(values)-1))], maximum=values[-1])


def analyze(path):
    contents = path.read_bytes()
    rows, iterations, options = {}, [], []
    failed = []
    for line in contents.decode().splitlines():
        fields = line.split()
        if not fields:
            continue
        name = fields[0]
        if name in ('RSP2_flux_residual', 'RSP2_flux_Y', 'RSP2_flux_before',
                    'RSP2_flux_after', 'RSP2_flux_dR'):
            key = tuple(map(int, fields[1:7]))
            rows.setdefault(key, {})[name[10:]] = list(map(float, fields[7:]))
        elif name == 'RSP2_flux_iteration':
            iterations.append((*map(int, fields[1:4]), fields[4], *map(float, fields[5:])))
        elif name == 'RSP2_flux_options':
            options.append(fields[1:])
        elif name == 'RSP2_flux_trial_failed':
            failed.append(line)

    # Select the larger of the two reported flux residuals for each trial.
    maxima = {}
    for key, row in rows.items():
        if 'residual' not in row:
            continue
        trial = key[:4] + key[5:]
        if trial not in maxima or abs(row['residual'][3]) > abs(maxima[trial][1]['residual'][3]):
            maxima[trial] = (key, row)

    # Backtracked trials legitimately have a nonzero linear prediction.
    late = [row for key, row in rows.items() if key[2] >= 4 and key[4] >= 300
            and all(name in row for name in ('residual', 'Y', 'before', 'after', 'dR'))
            and row['residual'][0] == 1]
    report = {
        'input': str(path.resolve()),
        'sha256': hashlib.sha256(contents).hexdigest(),
        'calls': len(options),
        'model_range': [min(int(o[0]) for o in options), max(int(o[0]) for o in options)],
        'options_after_model_solver_call': dict(Counter(' '.join(o[2:]) for o in options)),
        'tolerance_pass_iterations': dict(sorted(Counter(i[2] for i in iterations if i[3] == 'T').items())),
        'failed_equation_trials': failed,
        'reported_retry_counts': sorted(set(key[5] for key in rows)),
        'dt_seconds': distribution(row['residual'][1] for row in rows.values()),
        'by_iteration': {},
        'full_step_inner_trials_iteration_4_onward': {
            'actual_residual': distribution(abs(r['residual'][3]) for r in late),
            'predicted_residual': distribution(abs(r['residual'][4]) for r in late),
            'raw_Newton_residual': distribution(abs(r['residual'][5]) for r in late),
            'Y_update_error_in_residual_units': distribution(abs((r['Y'][5]-r['Y'][6])*r['Y'][7]) for r in late),
            'relative_pressure_difference': distribution(abs(r['after'][4]-r['after'][5])/r['after'][5] for r in late),
            'relative_temperature_difference': distribution(abs(r['after'][6]-r['after'][7])/r['after'][7] for r in late),
            'all_Lc_Lt_zero': all(r['after'][2] == r['after'][3] == 0 for r in late),
        },
    }
    for iteration in sorted(set(key[2] for key in rows)):
        subset = [(key, row) for key, row in maxima.values() if key[2] == iteration]
        report['by_iteration'][iteration] = {
            'max_flux_residual': distribution(abs(row['residual'][3]) for _, row in subset),
            'worst_faces': Counter(key[4] for key, _ in subset).most_common(10),
        }
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trace', type=Path)
    arguments = parser.parse_args()
    print(json.dumps(analyze(arguments.trace), indent=2))
