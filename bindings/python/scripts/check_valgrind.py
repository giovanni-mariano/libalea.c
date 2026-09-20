# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Reject repeated definite leaks originating in the pyAlea binding."""

import re
import sys
from pathlib import Path


text = Path(sys.argv[1]).read_text(encoding="utf-8", errors="replace")
records = re.split(r"\n(?===\d+== \d[\d,]* bytes in \d[\d,]* blocks are definitely lost)", text)
binding_leaks = []

for record in records:
    match = re.search(
        r"bytes in ([\d,]+) blocks are definitely lost", record
    )
    if not match or int(match.group(1).replace(",", "")) <= 10:
        continue
    if re.search(r"PyAlea|grid_array_from_owned_data|"
                 r"ray_slice_raster_array_from_owned_data", record):
        binding_leaks.append(record.strip())

if binding_leaks:
    print("Repeated definite leaks found in pyAlea:\n", file=sys.stderr)
    print("\n\n".join(binding_leaks), file=sys.stderr)
    raise SystemExit(1)

print("No repeated definite leaks found in pyAlea")
