#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Reconcile ray profiles and emit aggregates without geometry identifiers."""
import argparse
import csv
import json
from pathlib import Path

FIELDS = ('cell_calls', 'instructions', 'primitive_intersections', 'interval_steps')
RAY_FIELDS = ('bvh_visits', 'hierarchy_steps', 'parent_checks', 'parent_cache_hits', 'rejected_ranges')
LIMIT = 2**32 - 1


def read_rows(path, required):
    with path.open(newline='') as source:
        reader = csv.DictReader(source)
        if not set(required).issubset(reader.fieldnames or []):
            raise ValueError('Missing profile columns')
        rows = list(reader)
    for row in rows:
        for key in required:
            value = int(row[key])
            if not 0 <= value < LIMIT:
                raise ValueError('Invalid or saturated profile counter')
            row[key] = value
    return rows


def summarize(prefix):
    rays = read_rows(Path(str(prefix) + '-rays.csv'), FIELDS + RAY_FIELDS)
    cells = read_rows(Path(str(prefix) + '-cells.csv'), FIELDS + ('is_container',))
    if not rays:
        raise ValueError('Empty ray profile')
    totals = {key: sum(row[key] for row in rays) for key in FIELDS + RAY_FIELDS}
    for key in FIELDS:
        if totals[key] != sum(row[key] for row in cells):
            raise ValueError('Ray/cell totals disagree for ' + key)
    units = {row.get('instruction_unit', 'unspecified') for row in rays + cells}
    modes = {row.get('csg_groups', 'unspecified') for row in rays + cells}
    if len(units) != 1 or len(modes) != 1:
        raise ValueError('Inconsistent execution metadata')
    if (next(iter(modes)), next(iter(units))) not in {('0', 'bytecode'), ('1', 'bytecode_and_group_ops'), ('unspecified', 'unspecified')}:
        raise ValueError('Invalid execution metadata')
    roles = {}
    for container, name in ((0, 'terminal'), (1, 'container')):
        selected = [row for row in cells if row['is_container'] == container]
        row = {key: sum(item[key] for item in selected) for key in FIELDS}
        row['active_cells'] = len(selected)
        row['intervals_per_call'] = row['interval_steps'] / row['cell_calls'] if row['cell_calls'] else 0
        row['primitive_share_percent'] = 100 * row['primitive_intersections'] / totals['primitive_intersections'] if totals['primitive_intersections'] else 0
        roles[name] = row
    if sum(role['active_cells'] for role in roles.values()) != len(cells):
        raise ValueError('Invalid cell role')
    ranked = sorted((row['primitive_intersections'] for row in cells), reverse=True)
    concentration = {str(k): 100 * sum(ranked[:k]) / totals['primitive_intersections'] if totals['primitive_intersections'] else 0 for k in (1, 10, 100)}
    return {
        'schema_version': 1, 'rays': len(rays), 'active_cells': len(cells),
        'csg_groups': next(iter(modes)), 'instruction_unit': next(iter(units)),
        'totals': totals, 'roles': roles,
        'intervals_per_call': totals['interval_steps'] / totals['cell_calls'] if totals['cell_calls'] else 0,
        'primitives_per_interval': totals['primitive_intersections'] / totals['interval_steps'] if totals['interval_steps'] else 0,
        'top_cells_primitive_share_percent': concentration,
        'parent_cache_hit_percent': 100 * totals['parent_cache_hits'] / totals['parent_checks'] if totals['parent_checks'] else 0,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prefix', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    try:
        result = summarize(args.prefix)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, str(error) + '\n')
    text = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
