#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Summarize aggregate region eligibility counters; emit no model identities."""
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

FIELDS = ('calls', 'instructions', 'primitives', 'reject_instructions', 'reject_calls')
GROUP_FIELDS = ('group_instructions', 'group_primitives', 'group_entries', 'group_skips', 'group_combines')
KINDS = ('general', 'single', 'intersection')
REASONS = ('specialized', 'no_program', 'uncertified', 'outside_domain', 'unchanged')


def summarize(path):
    total = Counter()
    groups = defaultdict(Counter)
    seen = set()
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream)
        group_fields = set(GROUP_FIELDS).intersection(reader.fieldnames or ())
        if group_fields and group_fields != set(GROUP_FIELDS):
            raise ValueError('incomplete group counters')
        fields = FIELDS + GROUP_FIELDS if group_fields else FIELDS
        for row in reader:
            key = tuple(row[k] for k in ('kind', 'depth', 'reason'))
            if key in seen or key[0] not in KINDS or key[1] not in ('root', 'child') or key[2] not in REASONS:
                raise ValueError('duplicate or unknown eligibility group')
            seen.add(key)
            values = Counter({k: int(row[k]) for k in fields})
            if any(v < 0 or v >= 2**64 for v in values.values()):
                raise ValueError('invalid aggregate counter')
            if values['reject_calls'] > values['calls'] or values['reject_instructions'] > values['instructions']:
                raise ValueError('inconsistent rejection counts')
            if group_fields:
                if (values['group_instructions'] != values['group_primitives'] + values['group_skips'] + values['group_combines']
                        or values['group_instructions'] > values['instructions']
                        or values['group_primitives'] > values['primitives']
                        or values['group_entries'] > values['group_primitives']):
                    raise ValueError('inconsistent group counters')
            total.update(values)
            groups[key[0]].update(values)
    if not total['instructions']:
        raise ValueError('empty audit')
    result = {'totals': dict(total), 'classes': {
        kind: dict(groups[kind], instruction_percent=100*groups[kind]['instructions']/total['instructions'],
                   primitive_percent=100*groups[kind]['primitives']/total['primitives'] if total['primitives'] else 0)
        for kind in KINDS}}

    if group_fields:
        result['group_coverage'] = {
            'instruction_percent': 100*total['group_instructions']/total['instructions'],
            'primitive_percent': 100*total['group_primitives']/total['primitives'] if total['primitives'] else 0,
            'mixed_instruction_percent_of_total': 100*groups['general']['group_instructions']/total['instructions'],
            'control_instruction_percent_of_total': 100*(total['group_skips']+total['group_combines'])/total['instructions'],
            'mean_primitives_per_entry': total['group_primitives']/total['group_entries'] if total['group_entries'] else 0,
        }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('csv', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(summarize(args.csv), indent=2) + '\n')


if __name__ == '__main__':
    main()
