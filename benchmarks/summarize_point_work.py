#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Aggregate a work audit without exporting cell IDs, indices or coordinates."""
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

FIELDS = ('calls', 'instructions', 'primitives', 'reject_instructions', 'reject_calls')

def summarize(path, ray_paths):
    groups = defaultdict(Counter)
    cells = {}
    total = Counter()
    seen = set()
    with path.open(newline='') as stream:
        for row in csv.DictReader(stream):
            index = int(row['cell_index'])
            key = (index, row['reason'])
            if key in seen:
                raise ValueError('duplicate cell/reason row')
            seen.add(key)
            values = Counter({field: int(row[field]) for field in FIELDS})
            if any(v < 0 or v >= 2**32-1 for v in values.values()):
                raise ValueError('negative or saturated counter')
            if values['reject_calls'] > values['calls'] or values['reject_instructions'] > values['instructions']:
                raise ValueError('inconsistent rejection counts')
            identity = tuple(row[k] for k in ('cell_id', 'material_id', 'program_instructions', 'kind', 'point_certificate', 'skip_ops', 'threadable_ops'))
            cell = cells.setdefault(index, {'identity': identity, 'cost': 0, 'skips': int(row['skip_ops']), 'threadable': int(row['threadable_ops'])})
            if cell['identity'] != identity:
                raise ValueError('inconsistent cell metadata')
            if not 0 <= cell['threadable'] <= cell['skips']:
                raise ValueError('invalid branch metadata')
            cell['cost'] += values['instructions']
            groups[(row['kind'], row['point_certificate'], row['reason'])].update(values)
            total.update(values)
    if not total['instructions']:
        raise ValueError('empty work profile')
    ranking = sorted(cells, key=lambda i: (-cells[i]['cost'], i))
    ranks = (1, 10, 100, 1000)
    result = {
        'totals': dict(total),
        'visited_cells': len(cells),
        'original_branch_opportunities': {
            'cells_with_chains': sum(c['threadable'] > 0 for c in cells.values()),
            'skip_ops': sum(c['skips'] for c in cells.values()),
            'threadable_ops': sum(c['threadable'] for c in cells.values()),
            'work_in_cells_with_chains_percent': 100*sum(c['cost'] for c in cells.values() if c['threadable'])/total['instructions'],
        },
        'groups': [dict(kind=k[0], certificate=k[1], reason=k[2], **v,
                        work_percent=100*v['instructions']/total['instructions'])
                   for k, v in sorted(groups.items(), key=lambda pair: -pair[1]['instructions'])],
        'top_cell_work_percent': {str(n): 100*sum(cells[i]['cost'] for i in ranking[:n])/total['instructions'] for n in ranks},
        'ray_overlap': {},
    }
    top = {n: set(ranking[:n]) for n in ranks}
    for label, ray_path in ray_paths:
        work = matched = chain_work = 0
        overlap = Counter()
        kinds = Counter()
        observed = set()
        with ray_path.open(newline='') as stream:
            for row in csv.DictReader(stream):
                i = int(row['cell_index'])
                if i in observed:
                    raise ValueError('duplicate ray cell')
                observed.add(i)
                cost = int(row['instructions'])
                if cost < 0 or cost >= 2**32-1:
                    raise ValueError('negative or saturated ray counter')
                work += cost
                if i not in cells:
                    continue
                if cells[i]['identity'][:3] != tuple(row[k] for k in ('cell_id', 'material_id', 'program_instructions')):
                    raise ValueError('point/ray metadata mismatch: incompatible snapshots')
                matched += cost
                if cells[i]['threadable']:
                    chain_work += cost
                kinds[cells[i]['identity'][3]] += cost
                for n in ranks:
                    if i in top[n]:
                        overlap[n] += cost
        result['ray_overlap'][label] = {
            'instructions': work, 'matched_instructions': matched,
            'matched_kind_instructions': dict(kinds),
            'work_in_chain_cells_percent': 100*chain_work/work if work else 0,
            'work_in_top_point_cells_percent': {str(n): 100*overlap[n]/work if work else 0 for n in ranks},
        }
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('points', type=Path)
    p.add_argument('--rays', action='append', default=[], metavar='LABEL=CSV')
    p.add_argument('--output', required=True, type=Path)
    args = p.parse_args()
    rays = []
    for item in args.rays:
        label, path = item.split('=', 1)
        # Labels become report keys. Keep private paths out of generated summaries.
        if not label or not all(c.isascii() and (c.isalnum() or c in '-_') for c in label):
            p.error('ray labels must be ASCII letters, digits, hyphens or underscores')
        if label in [entry[0] for entry in rays]:
            p.error('duplicate ray label')
        rays.append((label, Path(path)))
    args.output.write_text(json.dumps(summarize(args.points, rays), indent=2)+'\n')


if __name__ == '__main__':
    main()
