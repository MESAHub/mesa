"""Read a saved solver trace/profile; never compile or execute MESA."""
from collections import Counter
from pathlib import Path
import argparse
import json
import math
import re
import statistics

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('trace', type=Path)
parser.add_argument('--profile', type=Path)
args = parser.parse_args()
pattern = re.compile(
    r'^\s*(\d+)\s+(\d+)\s+coeff\s+([\d.]+) avg resid\s+(\S+)\s+'
    r'max resid\s+(\S+)\s+(\d+)\s+(\S+) mix type (\S+) avg corr\s+(\S+)\s+'
    r'max corr\s+(\S+)\s+(\d+)\s+(\S+) mix type (\S+)\s+(.*)', re.M)
rows = []
for match in pattern.finditer(args.trace.read_text()):
    v = match.groups()
    rows.append(dict(model=int(v[0]), iteration=int(v[1]), coeff=float(v[2]),
                     avg_residual=float(v[3]), equation=v[4], zone=int(v[5]),
                     max_residual=float(v[6]), mixing=v[7], avg_correction=float(v[8]),
                     variable=v[9], correction_zone=int(v[10]),
                     max_correction=float(v[11]), status=v[13]))
if not rows:
    raise SystemExit('No solver-iteration rows found.')

summary = dict(models=[min(r['model'] for r in rows), max(r['model'] for r in rows)],
               iterations=len(rows),
               acceptance_iterations=dict(sorted(Counter(
                   r['iteration'] for r in rows if 'okay!' in r['status']).items())))
summary['by_iteration'] = []
for iteration in range(1, 11):
    group = [r for r in rows if r['iteration'] == iteration]
    if not group:
        continue
    summary['by_iteration'].append(dict(
        iteration=iteration, count=len(group),
        median_average_correction=statistics.median(r['avg_correction'] for r in group),
        median_maximum_residual=statistics.median(r['max_residual'] for r in group),
        zones=[min(r['zone'] for r in group), max(r['zone'] for r in group)],
        meets_printed_final_tolerances=sum(
            r['avg_correction'] < 3e-5 and r['avg_residual'] < 1e-8 and
            r['max_residual'] < 1e-5 for r in group)))

if args.profile:
    lines = args.profile.read_text().splitlines()
    header = dict(zip(lines[1].split(), lines[2].split()))
    columns = lines[5].split()
    cells = [dict(zip(columns, map(float, line.split()))) for line in lines[6:] if line.strip()]
    model = int(header['model_number'])
    summary['profile_model'] = model
    Lmax = max(abs(cell['luminosity']) for cell in cells)
    comparisons = []
    for row in rows:
        if row['model'] != model or row['iteration'] < 4 or row['equation'] != 'rsp2_fl':
            continue
        k = row['zone']
        if not 1 < k <= len(cells):
            continue
        cell, outer = cells[k-1], cells[k-2]
        alpha = outer['dq']/(outer['dq']+cell['dq'])
        beta = 1-alpha
        w = [10**(x['log_etrb']/2) if x['log_etrb'] > -98 else 0 for x in (cell, outer)]
        w_face = alpha*w[0]+beta*w[1]
        if w_face == 0 or cell['lum_conv'] == 0:
            continue
        # Use the final profile's L for a close estimate of the unavailable L_start.
        # Lc/w_face is directly measured from the saved state, avoiding an EOS model.
        scale = max(abs(cell['luminosity']), 1e-3*Lmax)
        response = abs(cell['lum_conv'])/(w_face*scale)
        predicted = [response*5e-5*weight for weight in (alpha, beta)]
        match = min(predicted, key=lambda p: abs(p-row['max_residual']))
        comparisons.append(dict(iteration=row['iteration'], zone=k,
                                observed=row['max_residual'],
                                predicted_reset_in_cell=predicted[0],
                                predicted_reset_in_outer_cell=predicted[1],
                                closest_relative_difference=abs(match/row['max_residual']-1)))
    summary['negative_w_repair_comparison'] = comparisons
    inner, outer = cells[-1], cells[-2]
    summary['inner_cell'] = dict(dq=inner['dq'],
                                relative_pressure_difference=abs(math.expm1(
                                    math.log(10)*(outer['logP']-inner['logP']))),
                                gradT=inner['gradT'], gradL=inner['gradL'])
print(json.dumps(summary, indent=2))
