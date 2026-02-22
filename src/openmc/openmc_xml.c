// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file openmc_xml.c
 * @brief Implementation of XML string builder for OpenMC inputs
 *
 * Option A:
 *  - puts_raw(x, s, n) rejects any NUL byte within the first n bytes (sticky error).
 *  - Helpers:
 *      PUTS_LIT(x, "literal")  -> copies sizeof(literal)-1 (no NUL).
 *      puts_cstr(x, s)         -> copies strlen(s).
 */

#include "openmc_xml.h"
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <limits.h>
#include <ctype.h>
#include <assert.h>

/* ---- helpers for safe writing of strings -------------------------------- */
static bool puts_raw(openmc_xml_t* x, const char* s, size_t n);
/* Literal helper: never copies the terminating NUL of a string literal */
#ifndef PUTS_LIT
#  define PUTS_LIT(x, lit) puts_raw((x), (lit), sizeof(lit) - 1)
#endif

/* C-string helper: uses strlen; safe for zero-terminated runtime strings */
static inline bool puts_cstr(openmc_xml_t* x, const char* s)
{
  return (s && *s) ? puts_raw(x, s, strlen(s)) : true;
}

/* ---- internal buffer management ----------------------------------------- */

static bool putc_raw(openmc_xml_t* x, char c)
{
  if (x->sb.error) return false;
  if (!str_builder_putc(&x->sb, c)) return false;
  if (c == '\n') x->col = 1;
  else x->col++;
  return true;
}

/* Option A: strict span writer (rejects embedded NUL in the first n bytes). */
static bool puts_raw(openmc_xml_t* x, const char* s, size_t n)
{
  if (x->sb.error) return false;
  if (!s || n == 0) return true;

  /* Reject any NUL byte inside the span. This catches sizeof(literal) mistakes. */
  const void* nulpos = memchr(s, '\0', n);
  if (nulpos != NULL) {
    assert(!"puts_raw: span contains NUL byte (probable sizeof(literal) used)");
    x->sb.error = true;
    return false;
  }

  if (!str_builder_write(&x->sb, s, n)) return false;

  /* update column accounting */
  for (size_t i = 0; i < n; ++i) {
    if (s[i] == '\n') x->col = 1;
    else x->col++;
  }
  return true;
}

/* ---- formatting helpers -------------------------------------------------- */

static bool newline_and_indent(openmc_xml_t* x)
{
  if (!putc_raw(x, '\n')) return false;
  if (!x->pretty) return true;
  return true;
}

static bool maybe_wrap_for_attr(openmc_xml_t* x, size_t incoming_len)
{
  if (!x->pretty) return true;
  if (x->no_wrap) return true;
  if (x->col + (int)incoming_len <= x->max_col) return true;

  if (!putc_raw(x, '\n')) return false;
  int base_spaces = x->indent_width * x->indent_level;
  int cont_spaces = base_spaces + x->indent_width;
  for (int i = 0; i < cont_spaces; ++i)
    if (!putc_raw(x, ' ')) return false;
  return true;
}

static size_t escaped_len(const char* s, bool for_attr)
{
  size_t n = 0;
  for (; *s; ++s) {
    switch (*s) {
      case '&': n += 5; break;      /* &amp;  */
      case '<': n += 4; break;      /* &lt;   */
      case '>': n += 4; break;      /* &gt;   */
      case '"': n += for_attr ? 6 : 1; break; /* &quot; in attrs */
      case '\'': n += for_attr ? 6 : 1; break;/* &apos; in attrs */
      default: n += 1; break;
    }
  }
  return n;
}

static bool put_escaped(openmc_xml_t* x, const char* s, bool for_attr)
{
  for (; *s; ++s) {
    switch (*s) {
      case '&':  if (!PUTS_LIT(x, "&amp;")) return false; break;
      case '<':  if (!PUTS_LIT(x, "&lt;"))  return false; break;
      case '>':  if (!PUTS_LIT(x, "&gt;"))  return false; break;
      case '"':  if (for_attr) { if (!PUTS_LIT(x, "&quot;")) return false; }
                 else { if (!putc_raw(x, '"')) return false; }
                 break;
      case '\'': if (for_attr) { if (!PUTS_LIT(x, "&apos;")) return false; }
                 else { if (!putc_raw(x, '\'')) return false; }
                 break;
      default:   if (!putc_raw(x, *s)) return false; break;
    }
  }
  return true;
}

/* Optional soft-wrap for text nodes on spaces when pretty-printing */
static bool put_text_wrapped(openmc_xml_t* x, const char* s)
{
  if (!x->pretty || x->max_col == INT_MAX) {
    return put_escaped(x, s, false);
  }

  const char* p = s;
  while (*p) {
    /* handle whitespace as wrap points (normalize spaces; keep newlines as breaks) */
    if (isspace((unsigned char)*p)) {
      if (*p == '\n') { if (!newline_and_indent(x)) return false; }
      else {
        if (x->col >= x->max_col) { if (!newline_and_indent(x)) return false; }
        else { if (!putc_raw(x, ' ')) return false; }
      }
      ++p;
      continue;
    }

    /* gather a word (non-space run) */
    const char* w = p;
    while (*p && !isspace((unsigned char)*p)) ++p;
    size_t wl = (size_t)(p - w);

    if (x->pretty) {
      int base = x->indent_width * x->indent_level + 1;
      if (x->col + (int)wl > x->max_col && x->col > base) {
        if (!newline_and_indent(x)) return false;
      }
    }

    /* emit escaped word (only &, <, > need escaping in text) */
    for (size_t i = 0; i < wl; ++i) {
      char c = w[i];
      switch (c) {
        case '&': if (!PUTS_LIT(x, "&amp;")) return false; break;
        case '<': if (!PUTS_LIT(x, "&lt;"))  return false; break;
        case '>': if (!PUTS_LIT(x, "&gt;"))  return false; break;
        default:  if (!putc_raw(x, c)) return false; break;
      }
    }
  }
  return true;
}

/* Format a float with requested style. Returns length or -1 on failure. */
static int format_float(char* tmp, size_t tmpsz, double value, int precision, bool scientific)
{
  if (precision < 0) precision = 6;
  if (scientific) return snprintf(tmp, tmpsz, "%.*e", precision, value);
  return snprintf(tmp, tmpsz, "%.*g", precision, value);
}

/* ---- public API --------------------------------------------------------- */

void openmc_xml_init(openmc_xml_t* x, arena_t* arena, size_t initial_capacity)
{
  openmc_xml_init_ex(x, arena, initial_capacity, 100, 2, true);
}

void openmc_xml_init_ex(openmc_xml_t* x, arena_t* arena, size_t initial_capacity,
                        int max_col, int indent_width, bool pretty)
{
  str_builder_init(&x->sb, arena, initial_capacity);
  x->col = 1;
  x->max_col = max_col;
  x->indent_width = indent_width;
  x->indent_level = 0;
  x->pretty = pretty;
  x->no_wrap = false;
  x->in_open_tag = false;
}

bool openmc_xml_start_document(openmc_xml_t* x, const char* version, const char* encoding)
{
  if (x->sb.error) return false;
  if (!PUTS_LIT(x, "<?xml version=\"")) return false;
  if (!puts_cstr(x, version ? version : "1.0")) return false;
  if (!PUTS_LIT(x, "\" encoding=\"")) return false;
  if (!puts_cstr(x, encoding ? encoding : "utf-8")) return false;
  if (!PUTS_LIT(x, "\"?>")) return false;
  if (x->pretty) return newline_and_indent(x);
  return true;
}

bool openmc_xml_start_element(openmc_xml_t* x, const char* name)
{
  if (x->sb.error || !name) return false;

  /* If a previous start tag is still open, close it and step into its content */
  if (x->in_open_tag) {
    if (!PUTS_LIT(x, ">")) return false;
    if (x->pretty) { if (!newline_and_indent(x)) return false; }
    x->in_open_tag = false;
    x->indent_level++; /* inside the previous element's content */
  }

  if (x->pretty) {
    int spaces = x->indent_width * x->indent_level;
    for (int i = 0; i < spaces; ++i)
      if (!putc_raw(x, ' ')) return false;
  }

  if (!PUTS_LIT(x, "<")) return false;
  if (!puts_cstr(x, name)) return false;
  x->in_open_tag = true;
  return true;
}

bool openmc_xml_attribute(openmc_xml_t* x, const char* name, const char* value)
{
  if (x->sb.error || !x->in_open_tag || !name || !value) return false;

  size_t val_len = escaped_len(value, true);
  size_t projected = 1 + strlen(name) + 2 + val_len + 1; /* ' ' + name + '="' + val + '"' */
  if (!maybe_wrap_for_attr(x, projected)) return false;

  if (!putc_raw(x, ' ')) return false;
  if (!puts_cstr(x, name)) return false;
  if (!PUTS_LIT(x, "=\"")) return false;
  if (!put_escaped(x, value, true)) return false;
  if (!PUTS_LIT(x, "\"")) return false;
  return true;
}
 
bool openmc_xml_end_start_tag(openmc_xml_t* x, bool self_close)
{
  if (x->sb.error || !x->in_open_tag) return false;

  if (self_close) {
    if (!PUTS_LIT(x, " />")) return false;
    x->in_open_tag = false;
    if (x->pretty) { if (!newline_and_indent(x)) return false; }
    /* No indent level change for self-closing tags */
  } else {
    if (!PUTS_LIT(x, ">")) return false;
    x->in_open_tag = false;
    if (x->pretty) { if (!newline_and_indent(x)) return false; }
    x->indent_level++; /* entering content of this element */
  }
  return true;
}

bool openmc_xml_end_element(openmc_xml_t* x, const char* name)
{
  if (x->sb.error || !name) return false;

  /* If the element had no content and tag is still open, convert to self-closing */
  if (x->in_open_tag) {
    if (!PUTS_LIT(x, " />")) return false;
    x->in_open_tag = false;
    if (x->pretty) { if (!newline_and_indent(x)) return false; }
    return true;
  }

  if (x->indent_level > 0) x->indent_level--;

  if (x->pretty) {
    if (x->col != 1) { 
      // if (!newline_and_indent(x)) return false;  // do we really need it?
      }
    else {
      int spaces = x->indent_width * x->indent_level;
      for (int i = 0; i < spaces; ++i)
        if (!putc_raw(x, ' ')) return false;
    }
  }

  if (!PUTS_LIT(x, "</")) return false;
  if (!puts_cstr(x, name)) return false;
  if (!PUTS_LIT(x, ">")) return false;
  if (x->pretty) { if (!newline_and_indent(x)) return false; }
  return true;
}

bool openmc_xml_text(openmc_xml_t* x, const char* text)
{
  if (x->sb.error || !text) return false;

  /* Close a pending start tag to enter content */
  if (x->in_open_tag) {
    if (!PUTS_LIT(x, ">")) return false;
    x->in_open_tag = false;
    if (x->pretty && x->col > x->indent_width * x->indent_level + 1) {
      /* continue on same line; wrapping may add newlines below */
    }
    x->indent_level++;
  }

  return put_text_wrapped(x, text);
}

/* ---- numeric helpers ---------------------------------------------------- */

bool openmc_xml_attribute_i(openmc_xml_t* x, const char* name, int value)
{
  if (x->sb.error || !name || !x->in_open_tag) return false;
  char buf[32];
  int n = snprintf(buf, sizeof(buf), "%d", value);
  if (n <= 0) return false;
  return openmc_xml_attribute(x, name, buf);
}

bool openmc_xml_attribute_f(openmc_xml_t* x, const char* name, double value,
                            int precision, bool scientific)
{
  if (x->sb.error || !name || !x->in_open_tag) return false;
  char buf[64];
  int n = format_float(buf, sizeof(buf), value, precision, scientific);
  if (n <= 0) return false;
  return openmc_xml_attribute(x, name, buf);
}

bool openmc_xml_numbers(openmc_xml_t* x, const double* vals, size_t n,
                        int precision, bool scientific, size_t max_per_line)
{
  if (x->sb.error || !vals) return false;

  /* Ensure we are in content context */
  if (x->in_open_tag) {
    if (!PUTS_LIT(x, ">")) return false;
    x->in_open_tag = false;
    x->indent_level++;
  }

  size_t on_line = 0;
  for (size_t i = 0; i < n; ++i) {
    char tmp[64];
    int  len = format_float(tmp, sizeof(tmp), vals[i], precision, scientific);
    if (len <= 0) return false;

    bool need_space = (i > 0);
    int token_len = len + (need_space ? 1 : 0);

    if (x->pretty) {
      int base = x->indent_width * x->indent_level + 1;
      bool width_wrap = (x->col + token_len > x->max_col && x->col > base);
      bool count_wrap = (max_per_line > 0 && on_line >= max_per_line);
      if (width_wrap || count_wrap) {
        if (!newline_and_indent(x)) return false;
        on_line = 0;
        need_space = false; /* first token on new line has no leading space */
      }
    }

    if (need_space) { if (!putc_raw(x, ' ')) return false; }
    if (!puts_raw(x, tmp, (size_t)len)) return false;
    on_line++;
  }
  return true;
}

/* ---- finalization ------------------------------------------------------- */

size_t openmc_xml_finish(openmc_xml_t* x)
{
  return str_builder_finish(&x->sb);
}

bool openmc_xml_error(const openmc_xml_t* x)
{
  return x->sb.error;
}

const char* openmc_xml_get(openmc_xml_t* x)
{
  openmc_xml_finish(x);
  return x->sb.buf;
}

bool openmc_xml_write(openmc_xml_t* x, FILE* out)
{
  if (x->sb.error || !out) return false;
  openmc_xml_finish(x);

  size_t wrote = fwrite(x->sb.buf, 1, x->sb.len, out);
  if (wrote != x->sb.len) return false;
  fputc('\n', out);
  return true;
}
