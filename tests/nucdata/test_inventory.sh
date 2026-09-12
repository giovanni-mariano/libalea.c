#!/bin/sh
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0

set -eu

inventory=../../bin/nuc_inventory

neutron=$($inventory njoy-test07/xsdir 92235.00c)
printf '%s\n' "$neutron" | awk -F '\t' '
    NR == 2 && $1 == "92235.00c" && $2 == "neutron" &&
    $3 == "yes" && $4 == "yes" && $5 == "yes" &&
    $6 ~ /fission/ && $6 ~ /delayed-neutron/ &&
    $6 ~ /photon-production/ && $7 == "ok" { found = 1 }
    END { exit !found }
'

law67=$($inventory njoy-test08/xsdir 28061.00c)
printf '%s\n' "$law67" | awk -F '\t' '
    NR == 2 && $1 == "28061.00c" && $2 == "neutron" &&
    $3 == "yes" && $4 == "yes" && $5 == "yes" &&
    $6 ~ /neutron-emission/ && $6 ~ /photon-production/ &&
    $7 == "ok" { found = 1 }
    END { exit !found }
'

photon=$($inventory njoy-test59/xsdir 92000.31p)
printf '%s\n' "$photon" | awk -F '\t' '
    NR == 2 && $1 == "92000.31p" && $2 == "photoatomic" &&
    $3 == "yes" && $4 == "yes" && $5 == "yes" &&
    $6 == "photon" && $7 == "ok" { found = 1 }
    END { exit !found }
'

thermal=$($inventory njoy-test25/xsdir lwtr.10t)
printf '%s\n' "$thermal" | awk -F '\t' '
    NR == 2 && $1 == "lwtr.10t" && $2 == "thermal" &&
    $3 == "yes" && $4 == "yes" && $5 == "yes" &&
    $6 == "thermal-sab" && $7 == "ok" { found = 1 }
    END { exit !found }
'

if $inventory njoy-test25/xsdir missing.10t >/dev/null 2>&1; then
    echo "missing ZAID unexpectedly succeeded" >&2
    exit 1
fi

echo "Capability inventory checks passed."
