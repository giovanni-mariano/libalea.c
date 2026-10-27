// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_EXPAT_CONFIG_H
#define ALEA_EXPAT_CONFIG_H

#include "expat_prefix.h"

/* Project-owned configuration; upstream sources under vendor/expat are intact. */
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
#define BYTEORDER 4321
#else
#define BYTEORDER 1234
#endif
#define XML_GE 0
#define XML_NS 1
#define XML_CONTEXT_BYTES 1024
#define XML_STATIC 1
#define HAVE_STDINT_H 1
#define HAVE_STDLIB_H 1
#define HAVE_STRING_H 1
#if defined(__EMSCRIPTEN__)
#define HAVE_GETENTROPY 1
#elif !defined(_WIN32)
#define XML_DEV_URANDOM 1
#endif

#endif
