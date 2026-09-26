// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_xml/xml_writer.h
 * @brief XML string builder for native ALEA model files
 */

#ifndef ALEA_XML_WRITER_H
#define ALEA_XML_WRITER_H

#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>
#include "util/str_builder.h"


typedef struct {
  str_builder_t sb;       /* Underlying buffer management */
  int     col;            /* Current column (1-based) */
  int     max_col;        /* Max column before wrapping (default 100) */
  int     indent_width;   /* Spaces per indent level (default 2) */
  int     indent_level;   /* Current element nesting depth */
  bool    pretty;         /* Pretty-print mode (newlines + indentation) */
  bool    no_wrap;        /* Suppress attribute wrapping (for cell elements) */
  bool    in_open_tag;    /* We wrote "<name" but not yet '>' or "/>" */
} alea_xml_writer_t;

/* ---- lifecycle ---- */
void   alea_xml_writer_init(alea_xml_writer_t* x, arena_t* arena, size_t initial_capacity);
void   alea_xml_writer_init_ex(alea_xml_writer_t* x, arena_t* arena, size_t initial_capacity,
                          int max_col, int indent_width, bool pretty);

/**
 * @brief Initialize in streaming mode.
 *
 * Uses a fixed @p buf_size malloc'd buffer.  Content is flushed directly to
 * @p stream whenever the buffer fills, keeping peak memory bounded.
 * Call alea_xml_writer_write() at the end to flush any remaining bytes.
 */
void   alea_xml_writer_init_stream_ex(alea_xml_writer_t* x, FILE* stream, size_t buf_size,
                                 int max_col, int indent_width, bool pretty);

/* ---- document ---- */
bool   alea_xml_writer_start_document(alea_xml_writer_t* x, const char* version, const char* encoding);

/* ---- elements & attributes ---- */
bool   alea_xml_writer_start_element(alea_xml_writer_t* x, const char* name);
bool   alea_xml_writer_attribute(alea_xml_writer_t* x, const char* name, const char* value);
bool   alea_xml_writer_end_start_tag(alea_xml_writer_t* x, bool self_close);
bool   alea_xml_writer_end_element(alea_xml_writer_t* x, const char* name);

/* ---- text ---- */
bool   alea_xml_writer_text(alea_xml_writer_t* x, const char* text);

/* ---- numeric helpers ---- */
bool   alea_xml_writer_attribute_i(alea_xml_writer_t* x, const char* name, int value);
bool   alea_xml_writer_attribute_f(alea_xml_writer_t* x, const char* name, double value,
                              int precision, bool scientific);
/* Writes a whitespace-separated list of numbers as element text.
   If max_per_line > 0, will insert a wrap after that many values regardless of width. */
bool   alea_xml_writer_numbers(alea_xml_writer_t* x, const double* vals, size_t n,
                          int precision, bool scientific, size_t max_per_line);

/* ---- finalization ---- */
size_t alea_xml_writer_finish(alea_xml_writer_t* x);
bool   alea_xml_writer_error(const alea_xml_writer_t* x);
const char* alea_xml_writer_get(alea_xml_writer_t* x);
bool   alea_xml_writer_write(alea_xml_writer_t* x, FILE* out);

#endif /* ALEA_XML_WRITER_H */
