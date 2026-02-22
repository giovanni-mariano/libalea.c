# Surface Deduplication Pipeline

This document describes the surface deduplication system: how primitives are
hashed, canonicalized, deduplicated, and how the correct cell sense is
preserved through export to MCNP and OpenMC formats.

## Overview

The deduplication pipeline has three stages:

1. **Primitive-level dedup** (parse time) — identical geometry gets the same
   `primitive_id` via a hash table with tolerance-based equality.
2. **Surface-level dedup** (export time) — multiple surface entries that
   reference the same `primitive_id` collapse to a single canonical surface ID.
3. **Sense correction** (export time) — the cell region expression adjusts
   each surface's sign to account for orientation differences introduced by
   canonicalization and merging.

```
  MCNP/OpenMC input
        │
        ▼
  ┌─────────────────────────┐
  │  Parse surface          │
  │  → canonicalize         │   alea_canonicalize_primitive()
  │  → hash                 │   alea_compute_primitive_hash()
  │  → dedup or create      │   alea_get_or_create_primitive()
  │  → store inverted flag  │
  └─────────┬───────────────┘
            │
            ▼
  ┌─────────────────────────┐
  │  Build canonical        │
  │  surface map            │   alea_build_canonical_surface_map()
  │  prim_id → surface_id   │
  │  prim_id → inverted     │
  └─────────┬───────────────┘
            │
            ▼
  ┌─────────────────────────┐
  │  Emit cell regions      │
  │  XOR sense correction   │   effective_sense = sense ⊕ (node_inv ≠ surf_inv)
  └─────────────────────────┘
```

---

## 1. Canonicalization

**File:** `src/core/alea_tolerance.c` — `alea_canonicalize_primitive()`

Before hashing or comparing, every primitive is put into canonical form. This
ensures that two surfaces describing the same geometry always produce the same
representation regardless of how the input specified them.

### Planes

A plane `ax + by + cz + d = 0` is canonicalized in two steps:

1. **Normalize** the normal vector to unit length:
   ```
   norm = sqrt(a² + b² + c²)
   a /= norm;  b /= norm;  c /= norm;  d /= norm;
   ```

2. **Ensure the first non-zero coefficient is positive:**
   ```
   if      (a < -ε)                          → negate (a,b,c,d), inverted=1
   else if (|a| < ε  and  b < -ε)            → negate (a,b,c,d), inverted=1
   else if (|a| < ε  and  |b| < ε  and c<-ε) → negate (a,b,c,d), inverted=1
   else                                       → inverted=0
   ```

The `inverted` output flag records whether the plane normal was flipped to
achieve canonical form. This flag is stored on every CSG node that references
the primitive and is used later to correct the cell sense during export.

### Other primitives

Spheres, cylinders, cones, tori, and general quadrics have a unique parametric
representation and require no canonicalization (`inverted` is always 0).

RPP (rectangular parallelepiped) swaps min/max bounds if inverted. All other
macrobodies (RCC, BOX, WED, etc.) are already canonical by construction.

---

## 2. Hash Function

**File:** `src/core/alea_primitive_dedup.c` — `alea_compute_primitive_hash()`

### Algorithm: FNV-1a with tolerance quantization

The hash serves as a fast filter before the full tolerance-based equality
check. Two primitives that are equal within tolerance must produce the same
hash, so floating-point values are quantized before hashing.

#### Double quantization (`hash_double`)

```c
static uint64_t hash_double(double value, const alea_config_t* tol) {
    if (fabs(value) < tol->zero_threshold)
        value = 0.0;                              // snap near-zero to 0
    else
        value = round(value / tol->abs_tol) * tol->abs_tol;  // grid-snap

    uint64_t bits;
    memcpy(&bits, &value, sizeof(double));         // bit-exact repr
    return bits;
}
```

Each double is snapped to a grid of spacing `abs_tol` (default 1e-14). Values
below `zero_threshold` (default 1e-14) collapse to 0. The resulting canonical
double is reinterpreted as a 64-bit integer for hashing.

#### FNV-1a core

```
hash = 14695981039346656037  (FNV offset basis)
for each byte b:
    hash ^= b
    hash *= 1099511628211    (FNV prime)
```

#### Per-type hashing

The hash starts with `fnv1a(type_enum)`, then XORs in the quantized doubles
for each type's parameters:

| Type       | Hashed fields                                    |
|------------|--------------------------------------------------|
| Plane      | a, b, c, d                                       |
| Sphere     | center_x/y/z, radius                             |
| Cylinder   | center_x/y, radius, z_min, z_max                 |
| Cone       | apex_x/y/z, tan²θ, axis_min/max, sheet_selection |
| Torus      | axis, center_x/y/z, R_major, r_minor, B          |
| Quadric    | 10 coefficients                                  |
| RPP        | min_x/y/z, max_x/y/z                             |
| RCC        | base_x/y/z, height_x/y/z, radius                |
| BOX        | corner_x/y/z, v1/v2/v3 (3×3)                    |
| TRC        | base_x/y/z, height_x/y/z, base_r, top_r         |
| ELL        | v1_x/y/z, v2_x/y/z, major_axis_len              |
| WED        | vertex_x/y/z, v1/v2/v3 (3×3)                    |
| RHP        | vertex_x/y/z, height_x/y/z, side1/2/3 (3×3)     |
| ARB        | num_corners, num_faces, all corners, face indices |

### Hash table structure

```
primitive_hash_table_t
├── buckets[]           — array of bucket chains
├── bucket_count        — initially 256
└── entry_count         — grows; resize at load factor > 0.75

primitive_hash_entry_t
├── primitive_id        — index into sys->primitives
├── hash                — cached 64-bit hash
└── *next               — collision chain (separate chaining)
```

**Lookup** (`primitive_hash_table_find`): compute bucket index as
`hash % bucket_count`, walk the chain, and for each entry with matching hash
call `alea_primitives_equal()` for the full tolerance check.

**Insert**: prepend to bucket chain; double bucket count when
`entry_count > bucket_count * 3/4`.

---

## 3. Tolerance-Based Equality

**File:** `src/core/alea_tolerance.c` — `alea_primitives_equal()`

After a hash match, the full equality check compares every parameter using a
three-tier tolerance:

```c
bool alea_numbers_equal(double a, double b, const alea_config_t* tol) {
    if (a == b) return true;                           // exact
    if (|a| < zero_threshold && |b| < zero_threshold)  // both ~0
        return true;
    if (|a - b| <= abs_tol) return true;               // absolute
    if (|a - b| <= max(|a|,|b|) * rel_tol) return true; // relative
    return false;
}
```

Default tolerances: `abs_tol = 1e-14`, `rel_tol = 1e-12`,
`zero_threshold = 1e-14`.

### Plane opposite-normal detection

For planes, the equality check also detects opposite normals:

```
dot = a1·a2 + b1·b2 + c1·c2
if |dot| ≈ 1:
    if dot > 0: same direction  → check d1 ≈  d2
    if dot < 0: opposite        → check d1 ≈ -d2, set match_inverted = 1
```

When `match_inverted = 1`, the caller (`alea_get_or_create_primitive`) toggles
the `inverted` flag on the node: `*inverted = !(*inverted)`. This encodes that
the matched primitive has an opposite normal to the one being registered.

---

## 4. Primitive Creation and Dedup

**File:** `src/core/alea_system.c` — `alea_get_or_create_primitive()`

This is the single entry point for all primitive creation. Every surface
parsed from MCNP or OpenMC goes through this function.

```
alea_get_or_create_primitive(sys, type, data, &inverted):
    1. canonicalize(type, data, &inverted)
    2. hash = compute_hash(type, data)
    3. existing = hash_table_find(hash, type, data, &match_inverted)
    4. if existing found:
         primitives[existing].ref_count++
         if match_inverted: inverted = !inverted    ← opposite normal
         return existing
    5. else:
         id = allocate new primitive(type, data)
         hash_table_insert(id, hash)
         return id
```

The returned `inverted` flag is stored on the CSG tree node as
`node->primitive.inverted`. It encodes the full chain of orientation
transformations: canonicalization flip (if any) combined with opposite-normal
dedup match (if any).

---

## 5. Surface Entries

**Defined in:** `include/alea_types.h`

```c
typedef struct {
    int mcnp_surface_id;              // unique ID in exported file
    alea_primitive_id_t primitive_id;  // index into sys->primitives
    alea_node_id_t pos_node;          // node with sense=+1 (outside)
    alea_node_id_t neg_node;          // node with sense=-1 (inside)
    int transform_id;                 // TRn (0 = none)
    bool transform_applied;
    alea_boundary_type_t boundary_type;
    int periodic_surface_id;
    alea_node_id_t expanded_pos_node; // macrobody expansion
    alea_node_id_t expanded_neg_node;
} alea_surface_entry_t;
```

The `pos_node` field is critical for deduplication: it identifies the CSG node
whose `inverted` flag represents the canonical orientation of this surface.
The canonical surface map reads `pos_node->primitive.inverted` and stores it
as `prim_to_surface_inverted[prim_id]`.

### Synthetic surface entries

Macrobody expansion decomposes a single MCNP surface (e.g., an RPP) into
multiple primitive surfaces (6 planes). These child primitives initially have
`mcnp_surface_id = 0` because they don't correspond to any input surface.

Both export paths (`alea_export.c` and `openmc_export.c`) walk all cell trees
and assign synthetic surface IDs to these primitives via
`assign_surface_ids_recursive`. When a new synthetic ID is assigned, a surface
entry is also pushed to `sys->surfaces`:

```c
alea_surface_entry_t* entry = alea_vec_push_uninit(&sys->surfaces, ...);
entry->mcnp_surface_id = new_id;
entry->primitive_id    = prim_id;
entry->pos_node        = node_id;   // ← records the canonical inverted flag
entry->neg_node        = ALEA_NODE_ID_INVALID;
```

Without this registration, `alea_build_canonical_surface_map` would not know
about the primitive, leaving `prim_to_surface_inverted[prim_id] = 0` and
causing incorrect sense in cell regions.

---

## 6. Canonical Surface Map

**File:** `src/core/alea_export.c` — `alea_build_canonical_surface_map()`

After all surface entries exist (including synthetic ones), this function
builds two lookup arrays indexed by `primitive_id`:

```
prim_to_surface[prim_id]          → canonical mcnp_surface_id
prim_to_surface_inverted[prim_id] → inverted flag from canonical pos_node
```

### Algorithm

```
for each surface entry in sys->surfaces:
    prim_id = entry.primitive_id
    if this is the first entry for prim_id, or entry.mcnp_surface_id is
    lower than the current canonical:
        prim_to_surface[prim_id] = entry.mcnp_surface_id
        prim_to_surface_inverted[prim_id] = nodes[entry.pos_node].inverted
```

The **lowest** `mcnp_surface_id` wins as the canonical representative for each
primitive. All other surfaces with the same `primitive_id` are dedup'd away —
they won't appear in the exported surface block.

### Example

```
Surface 10: prim_id=5, pos_node has inverted=0
Surface 20: prim_id=5, pos_node has inverted=1

→ prim_to_surface[5] = 10        (lowest ID wins)
→ prim_to_surface_inverted[5] = 0  (from surface 10's pos_node)
→ Surface 20 is dedup'd out
```

---

## 7. Cell Sense Correction (XOR Formula)

When multiple surfaces collapse to the same canonical surface, a node that
originally referenced a different orientation must have its sense corrected.

### The problem

Consider two nodes referencing the same primitive but with different `inverted`
flags:

```
Node A: prim_id=5, sense=+1, inverted=0   (original surface 10)
Node B: prim_id=5, sense=+1, inverted=1   (original surface 20)
```

Without dedup, Node A emits `+10` and Node B emits `+20`. Both describe
"outside their respective surface", which is correct because the surfaces have
opposite normals.

With dedup, both use surface 10. But Node B's "outside" was relative to the
*opposite* normal. Emitting `+10` for Node B would be wrong — it should be
`-10`.

### The formula

Both MCNP and OpenMC exporters use the same XOR correction:

```c
int effective_sense = node->primitive.sense;

int8_t surface_inverted = prim_to_surface_inverted[prim_id];
if (node->primitive.inverted != surface_inverted) {
    effective_sense = -effective_sense;  // flip
}
```

In the example: Node B has `inverted=1`, canonical surface has
`surface_inverted=0`. Since `1 ≠ 0`, the sense is flipped: `+1 → -1`,
producing `-10`. Correct.

### Truth table

```
node.inverted  surface_inverted  sense  →  effective_sense
     0              0              +1        +1  (no flip)
     0              0              -1        -1  (no flip)
     1              0              +1        -1  (flip)
     1              0              -1        +1  (flip)
     0              1              +1        -1  (flip)
     0              1              -1        +1  (flip)
     1              1              +1        +1  (no flip)
     1              1              -1        -1  (no flip)
```

The rule is: **flip sense when `node.inverted ⊕ surface_inverted = 1`** (XOR).

### MCNP vs OpenMC difference

The XOR formula is the same in both exporters. The difference is in how surface
*coefficients* are emitted:

- **MCNP exporter**: un-canonicalizes plane coefficients by flipping signs
  according to `pos_node->inverted`. The exported surface may have opposite
  coefficients from the canonical primitive. The XOR formula accounts for this.

- **OpenMC exporter**: emits canonical coefficients as-is (no
  un-canonicalization). The XOR formula still works because
  `prim_to_surface_inverted` records the same `pos_node->inverted` flag used
  by the canonical surface map.

---

## 8. End-to-End Example

### Input (MCNP)

```
c Surfaces
10  pz  5.0      $ z = 5 (normal pointing +z)
20  pz  5.0      $ same plane, duplicate
```

### Parse time

Surface 10:
```
plane (a=0, b=0, c=1, d=-5)
canonicalize → c>0, no flip → inverted=0
hash(0,0,1,-5) → H1
hash_table: empty → create primitive 0
surface_entry: {id=10, prim=0, pos_node=N1}
node N1: {prim_id=0, inverted=0, sense=+1}
```

Surface 20:
```
plane (a=0, b=0, c=1, d=-5)
canonicalize → c>0, no flip → inverted=0
hash(0,0,1,-5) → H1 (same!)
hash_table: found prim 0, match_inverted=0
→ reuse prim 0, ref_count=2
surface_entry: {id=20, prim=0, pos_node=N2}
node N2: {prim_id=0, inverted=0, sense=+1}
```

### Export time (with dedup)

Canonical map:
```
prim_to_surface[0] = 10         (lowest ID)
prim_to_surface_inverted[0] = 0 (from N1.inverted)
```

Cell using `+20`:
```
node.inverted=0, surface_inverted=0 → no flip → emit +10   ✓
```

### Opposite normal variant

If surface 20 had been defined with opposite normal (`pz -5.0` interpreted as
`c=-1, d=5`):

```
canonicalize → c<0, flip → (0,0,1,-5), inverted=1
hash(0,0,1,-5) → H1 (same!)
hash_table: found prim 0 with dot=-1 → match_inverted=1
inverted = !1 = 0... wait, let's trace carefully:

canonicalize sets inverted=1 (because we flipped the normal)
match_inverted=1 (opposite normal match)
inverted = !inverted = !1 = 0

Hmm — but this node's "outside" is on the -z side (original normal was -z).
The canonical primitive's +1 sense means +z side. So the XOR should flip.

Actually: the original surface had normal -z. Canonicalization flipped to +z
and set inverted=1. Then the hash match found it's the same primitive with
the SAME canonical normal (not opposite — the canonicalization already
aligned them). So match_inverted=0 in this case.

Final: inverted=1, surface_inverted=0 → flip → correct.
```

This shows why the `inverted` flag must flow through both canonicalization AND
the dedup match to encode the complete orientation relationship.
