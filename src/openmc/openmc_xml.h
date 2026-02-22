// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file openmc_xml.h
 * @brief XML string builder tailored for OpenMC input files
 */

#ifndef OPENMC_XML_H
#define OPENMC_XML_H

#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>
#include "util/str_builder.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  str_builder_t sb;       /* Underlying buffer management */
  int     col;            /* Current column (1-based) */
  int     max_col;        /* Max column before wrapping (default 100) */
  int     indent_width;   /* Spaces per indent level (default 2) */
  int     indent_level;   /* Current element nesting depth */
  bool    pretty;         /* Pretty-print mode (newlines + indentation) */
  bool    no_wrap;        /* Suppress attribute wrapping (for cell elements) */
  bool    in_open_tag;    /* We wrote "<name" but not yet '>' or "/>" */
} openmc_xml_t;

/* ---- lifecycle ---- */
void   openmc_xml_init(openmc_xml_t* x, arena_t* arena, size_t initial_capacity);
void   openmc_xml_init_ex(openmc_xml_t* x, arena_t* arena, size_t initial_capacity,
                          int max_col, int indent_width, bool pretty);

/* ---- document ---- */
bool   openmc_xml_start_document(openmc_xml_t* x, const char* version, const char* encoding);

/* ---- elements & attributes ---- */
bool   openmc_xml_start_element(openmc_xml_t* x, const char* name);
bool   openmc_xml_attribute(openmc_xml_t* x, const char* name, const char* value);
bool   openmc_xml_end_start_tag(openmc_xml_t* x, bool self_close);
bool   openmc_xml_end_element(openmc_xml_t* x, const char* name);

/* ---- text ---- */
bool   openmc_xml_text(openmc_xml_t* x, const char* text);

/* ---- numeric helpers ---- */
bool   openmc_xml_attribute_i(openmc_xml_t* x, const char* name, int value);
bool   openmc_xml_attribute_f(openmc_xml_t* x, const char* name, double value,
                              int precision, bool scientific);
/* Writes a whitespace-separated list of numbers as element text.
   If max_per_line > 0, will insert a wrap after that many values regardless of width. */
bool   openmc_xml_numbers(openmc_xml_t* x, const double* vals, size_t n,
                          int precision, bool scientific, size_t max_per_line);

/* ---- finalization ---- */
size_t openmc_xml_finish(openmc_xml_t* x);
bool   openmc_xml_error(const openmc_xml_t* x);
const char* openmc_xml_get(openmc_xml_t* x);
bool   openmc_xml_write(openmc_xml_t* x, FILE* out);

#ifdef __cplusplus
}
#endif
#endif /* OPENMC_XML_H */
