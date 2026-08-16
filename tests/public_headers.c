// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * Compile-only check: every installed header must resolve with -Iinclude and
 * without src/ on the include path. Keep this list in sync with include/.
 */
#include "alea.h"
#include "alea_geo_validator.h"
#include "alea_mcnp.h"
#include "alea_mesh.h"
#include "alea_nucdata.h"
#include "alea_nucdata_types.h"
#include "alea_openmc.h"
#include "alea_raycast.h"
#include "alea_render.h"
#include "alea_serpent.h"
#include "alea_slice.h"
#include "alea_types.h"
