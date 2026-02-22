// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file openmc_parse.c
 * @brief Minimal XML parser implementation for OpenMC geometry files
 *
 * Hand-rolled XML parser. OpenMC XML is simple - no need for libxml2.
 * Handles: elements, attributes, text content, basic entity escapes.
 * Ignores: comments, processing instructions, CDATA, namespaces, DTD.
 */

#include "openmc_parse.h"
#include "util/str_builder.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* ============================================================================
 * LEXER
 * ============================================================================ */

typedef enum {
    XML_TOK_EOF,
    XML_TOK_ERROR,
    XML_TOK_OPEN_TAG,       /* < */
    XML_TOK_CLOSE_TAG,      /* </ */
    XML_TOK_TAG_END,        /* > */
    XML_TOK_SELF_CLOSE,     /* /> */
    XML_TOK_NAME,           /* tag/attr name */
    XML_TOK_EQUALS,         /* = */
    XML_TOK_STRING,         /* "value" or 'value' */
    XML_TOK_TEXT,           /* text content */
    XML_TOK_COMMENT,        /* <!-- ... --> */
    XML_TOK_PI              /* <?...?> */
} xml_token_type_t;

typedef struct {
    const char* src;
    size_t len;
    size_t pos;
    int line;
    int col;

    /* Current token */
    xml_token_type_t token;
    str_builder_t token_sb;

    /* Arena for allocations */
    arena_t* arena;

    /* Error state */
    char error_msg[256];
    bool has_error;
} xml_lexer_t;

static void lexer_init(xml_lexer_t* lex, const char* src, size_t len, arena_t* arena) {
    lex->src = src;
    lex->len = len;
    lex->pos = 0;
    lex->line = 1;
    lex->col = 1;
    lex->token = XML_TOK_EOF;
    str_builder_init(&lex->token_sb, arena, 64);
    lex->arena = arena;
    lex->error_msg[0] = '\0';
    lex->has_error = false;
}

static char peek(xml_lexer_t* lex) {
    if (lex->pos >= lex->len) return '\0';
    return lex->src[lex->pos];
}

static char peek_n(xml_lexer_t* lex, size_t n) {
    if (lex->pos + n >= lex->len) return '\0';
    return lex->src[lex->pos + n];
}

static char advance(xml_lexer_t* lex) {
    if (lex->pos >= lex->len) return '\0';
    char c = lex->src[lex->pos++];
    if (c == '\n') {
        lex->line++;
        lex->col = 1;
    } else {
        lex->col++;
    }
    return c;
}

static void skip_whitespace(xml_lexer_t* lex) {
    while (isspace((unsigned char)peek(lex))) {
        advance(lex);
    }
}

static void set_error(xml_lexer_t* lex, const char* msg) {
    snprintf(lex->error_msg, sizeof(lex->error_msg),
             "Line %d, col %d: %s", lex->line, lex->col, msg);
    lex->has_error = true;
    lex->token = XML_TOK_ERROR;
}

static void token_append(xml_lexer_t* lex, char c) {
    if (!str_builder_putc(&lex->token_sb, c)) {
        set_error(lex, "Out of memory");
    }
}

static void token_clear(xml_lexer_t* lex) {
    str_builder_reset(&lex->token_sb);
}

/* Decode basic XML entity references */
static bool decode_entity(xml_lexer_t* lex) {
    /* Already consumed '&' */
    char buf[16];
    int i = 0;

    while (i < 15 && peek(lex) != '\0' && peek(lex) != ';') {
        buf[i++] = advance(lex);
    }
    buf[i] = '\0';

    if (peek(lex) != ';') {
        set_error(lex, "Unterminated entity reference");
        return false;
    }
    advance(lex); /* consume ';' */

    if (strcmp(buf, "lt") == 0) token_append(lex, '<');
    else if (strcmp(buf, "gt") == 0) token_append(lex, '>');
    else if (strcmp(buf, "amp") == 0) token_append(lex, '&');
    else if (strcmp(buf, "quot") == 0) token_append(lex, '"');
    else if (strcmp(buf, "apos") == 0) token_append(lex, '\'');
    else if (buf[0] == '#') {
        /* Numeric character reference */
        int code;
        if (buf[1] == 'x') {
            code = (int)strtol(buf + 2, NULL, 16);
        } else {
            code = (int)strtol(buf + 1, NULL, 10);
        }
        if (code > 0 && code < 128) {
            token_append(lex, (char)code);
        }
    } else {
        /* Unknown entity - just copy literally */
        token_append(lex, '&');
        for (int j = 0; j < i; j++) token_append(lex, buf[j]);
        token_append(lex, ';');
    }

    return true;
}

static bool is_name_start_char(char c) {
    return isalpha((unsigned char)c) || c == '_' || c == ':';
}

static bool is_name_char(char c) {
    return isalnum((unsigned char)c) || c == '_' || c == ':' || c == '-' || c == '.';
}

static void lex_name(xml_lexer_t* lex) {
    token_clear(lex);
    while (is_name_char(peek(lex))) {
        token_append(lex, advance(lex));
    }
    lex->token = XML_TOK_NAME;
}

static void lex_string(xml_lexer_t* lex) {
    char quote = advance(lex); /* consume opening quote */
    token_clear(lex);

    while (peek(lex) != '\0' && peek(lex) != quote) {
        char c = peek(lex);
        if (c == '&') {
            advance(lex);
            if (!decode_entity(lex)) return;
        } else {
            token_append(lex, advance(lex));
        }
    }

    if (peek(lex) != quote) {
        set_error(lex, "Unterminated string");
        return;
    }
    advance(lex); /* consume closing quote */
    lex->token = XML_TOK_STRING;
}

static void skip_comment(xml_lexer_t* lex) {
    /* Already consumed '<!--' */
    while (lex->pos + 2 < lex->len) {
        if (peek(lex) == '-' && peek_n(lex, 1) == '-' && peek_n(lex, 2) == '>') {
            advance(lex); advance(lex); advance(lex);
            lex->token = XML_TOK_COMMENT;
            return;
        }
        advance(lex);
    }
    set_error(lex, "Unterminated comment");
}

static void skip_pi(xml_lexer_t* lex) {
    /* Already consumed '<?' */
    while (lex->pos + 1 < lex->len) {
        if (peek(lex) == '?' && peek_n(lex, 1) == '>') {
            advance(lex); advance(lex);
            lex->token = XML_TOK_PI;
            return;
        }
        advance(lex);
    }
    set_error(lex, "Unterminated processing instruction");
}

static void lex_text(xml_lexer_t* lex) {
    token_clear(lex);
    while (peek(lex) != '\0' && peek(lex) != '<') {
        char c = peek(lex);
        if (c == '&') {
            advance(lex);
            if (!decode_entity(lex)) return;
        } else {
            token_append(lex, advance(lex));
        }
    }

    /* Trim trailing whitespace */
    while (lex->token_sb.len > 0 &&
           isspace((unsigned char)lex->token_sb.buf[lex->token_sb.len - 1])) {
        lex->token_sb.len--;
    }
    if (lex->token_sb.buf) lex->token_sb.buf[lex->token_sb.len] = '\0';

    lex->token = XML_TOK_TEXT;
}

/* State: inside a tag (after '<' or '</') */
static void next_token_in_tag(xml_lexer_t* lex) {
    skip_whitespace(lex);

    char c = peek(lex);

    if (c == '\0') {
        lex->token = XML_TOK_EOF;
    } else if (c == '>') {
        advance(lex);
        lex->token = XML_TOK_TAG_END;
    } else if (c == '/' && peek_n(lex, 1) == '>') {
        advance(lex); advance(lex);
        lex->token = XML_TOK_SELF_CLOSE;
    } else if (c == '=') {
        advance(lex);
        lex->token = XML_TOK_EQUALS;
    } else if (c == '"' || c == '\'') {
        lex_string(lex);
    } else if (is_name_start_char(c)) {
        lex_name(lex);
    } else {
        set_error(lex, "Unexpected character in tag");
    }
}

/* State: outside tags (content) */
static void next_token_content(xml_lexer_t* lex) {
    /* Skip leading whitespace in content */
    skip_whitespace(lex);

    char c = peek(lex);

    if (c == '\0') {
        lex->token = XML_TOK_EOF;
    } else if (c == '<') {
        advance(lex);
        c = peek(lex);

        if (c == '/') {
            advance(lex);
            lex->token = XML_TOK_CLOSE_TAG;
        } else if (c == '!') {
            /* Comment or DOCTYPE */
            if (peek_n(lex, 1) == '-' && peek_n(lex, 2) == '-') {
                advance(lex); advance(lex); advance(lex);
                skip_comment(lex);
            } else {
                /* Skip DOCTYPE or other declarations */
                int depth = 1;
                while (depth > 0 && peek(lex) != '\0') {
                    c = advance(lex);
                    if (c == '<') depth++;
                    else if (c == '>') depth--;
                }
                lex->token = XML_TOK_COMMENT; /* Treat as comment */
            }
        } else if (c == '?') {
            advance(lex);
            skip_pi(lex);
        } else {
            lex->token = XML_TOK_OPEN_TAG;
        }
    } else {
        lex_text(lex);
    }
}

/* ============================================================================
 * PARSER
 * ============================================================================ */

typedef struct {
    xml_lexer_t lex;
    openmc_xml_doc_t* doc;
    openmc_xml_element_t* current;
    bool in_tag;
} xml_parser_t;

static openmc_xml_element_t* new_element(xml_parser_t* p, const char* tag_name) {
    openmc_xml_element_t* elem = (openmc_xml_element_t*)arena_alloc(
        &p->doc->arena, sizeof(openmc_xml_element_t));
    if (!elem) return NULL;

    memset(elem, 0, sizeof(*elem));
    elem->tag_name = arena_strdup(&p->doc->arena, tag_name);

    return elem;
}

static bool add_attr(xml_parser_t* p, openmc_xml_element_t* elem,
                     const char* name, const char* value) {
    if (elem->attr_count >= elem->attr_capacity) {
        size_t new_cap = elem->attr_capacity ? elem->attr_capacity * 2 : 4;
        openmc_xml_attr_t* new_attrs = (openmc_xml_attr_t*)arena_alloc(
            &p->doc->arena, new_cap * sizeof(openmc_xml_attr_t));
        if (!new_attrs) return false;

        if (elem->attrs && elem->attr_count > 0) {
            memcpy(new_attrs, elem->attrs, elem->attr_count * sizeof(openmc_xml_attr_t));
        }
        elem->attrs = new_attrs;
        elem->attr_capacity = new_cap;
    }

    openmc_xml_attr_t* attr = &elem->attrs[elem->attr_count++];
    attr->name = arena_strdup(&p->doc->arena, name);
    attr->value = arena_strdup(&p->doc->arena, value);

    return true;
}

static bool add_child(xml_parser_t* p, openmc_xml_element_t* parent,
                      openmc_xml_element_t* child) {
    if (parent->child_count >= parent->child_capacity) {
        size_t new_cap = parent->child_capacity ? parent->child_capacity * 2 : 4;
        openmc_xml_element_t** new_children = (openmc_xml_element_t**)arena_alloc(
            &p->doc->arena, new_cap * sizeof(openmc_xml_element_t*));
        if (!new_children) return false;

        if (parent->children && parent->child_count > 0) {
            memcpy(new_children, parent->children,
                   parent->child_count * sizeof(openmc_xml_element_t*));
        }
        parent->children = new_children;
        parent->child_capacity = new_cap;
    }

    parent->children[parent->child_count++] = child;
    child->parent = parent;

    return true;
}

static void parser_error(xml_parser_t* p, const char* msg) {
    snprintf(p->doc->error_msg, sizeof(p->doc->error_msg),
             "Line %d, col %d: %s", p->lex.line, p->lex.col, msg);
    p->doc->error_line = p->lex.line;
    p->doc->error_col = p->lex.col;
}

static bool parse_element(xml_parser_t* p);

static bool parse_content(xml_parser_t* p, openmc_xml_element_t* parent) {
    while (true) {
        next_token_content(&p->lex);

        if (p->lex.has_error) {
            parser_error(p, p->lex.error_msg);
            return false;
        }

        if (p->lex.token == XML_TOK_EOF) {
            parser_error(p, "Unexpected end of file");
            return false;
        }

        if (p->lex.token == XML_TOK_COMMENT || p->lex.token == XML_TOK_PI) {
            continue; /* Skip comments and PIs */
        }

        if (p->lex.token == XML_TOK_CLOSE_TAG) {
            /* End of current element */
            return true;
        }

        if (p->lex.token == XML_TOK_TEXT) {
            /* Text content */
            if (p->lex.token_sb.len > 0) {
                parent->text_content = arena_strdup(&p->doc->arena, p->lex.token_sb.buf);
            }
            continue;
        }

        if (p->lex.token == XML_TOK_OPEN_TAG) {
            /* Child element */
            p->current = parent;
            if (!parse_element(p)) return false;
            continue;
        }

        parser_error(p, "Unexpected token in content");
        return false;
    }
}

static bool parse_element(xml_parser_t* p) {
    /* Token should be OPEN_TAG, need to read tag name */
    next_token_in_tag(&p->lex);

    if (p->lex.token != XML_TOK_NAME) {
        parser_error(p, "Expected element name");
        return false;
    }

    openmc_xml_element_t* elem = new_element(p, p->lex.token_sb.buf);
    if (!elem) {
        parser_error(p, "Out of memory");
        return false;
    }

    /* Add to parent or set as root */
    if (p->current) {
        if (!add_child(p, p->current, elem)) {
            parser_error(p, "Out of memory");
            return false;
        }
    } else {
        p->doc->root = elem;
    }

    /* Parse attributes */
    while (true) {
        next_token_in_tag(&p->lex);

        if (p->lex.has_error) {
            parser_error(p, p->lex.error_msg);
            return false;
        }

        if (p->lex.token == XML_TOK_TAG_END) {
            /* Element has content */
            if (!parse_content(p, elem)) return false;

            /* Should be at close tag name now */
            next_token_in_tag(&p->lex);
            if (p->lex.token != XML_TOK_NAME ||
                strcmp(p->lex.token_sb.buf, elem->tag_name) != 0) {
                parser_error(p, "Mismatched close tag");
                return false;
            }

            next_token_in_tag(&p->lex);
            if (p->lex.token != XML_TOK_TAG_END) {
                parser_error(p, "Expected '>'");
                return false;
            }
            return true;
        }

        if (p->lex.token == XML_TOK_SELF_CLOSE) {
            /* Self-closing element */
            return true;
        }

        if (p->lex.token == XML_TOK_NAME) {
            /* Attribute */
            char* attr_name = arena_strdup(&p->doc->arena, p->lex.token_sb.buf);

            next_token_in_tag(&p->lex);
            if (p->lex.token != XML_TOK_EQUALS) {
                parser_error(p, "Expected '=' after attribute name");
                return false;
            }

            next_token_in_tag(&p->lex);
            if (p->lex.token != XML_TOK_STRING) {
                parser_error(p, "Expected attribute value");
                return false;
            }

            if (!add_attr(p, elem, attr_name, p->lex.token_sb.buf)) {
                parser_error(p, "Out of memory");
                return false;
            }
            continue;
        }

        parser_error(p, "Unexpected token in tag");
        return false;
    }
}

static bool parse_document(xml_parser_t* p) {
    /* Skip leading whitespace, comments, PIs, XML declaration */
    while (true) {
        next_token_content(&p->lex);

        if (p->lex.has_error) {
            parser_error(p, p->lex.error_msg);
            return false;
        }

        if (p->lex.token == XML_TOK_EOF) {
            parser_error(p, "Empty document");
            return false;
        }

        if (p->lex.token == XML_TOK_COMMENT || p->lex.token == XML_TOK_PI) {
            continue;
        }

        if (p->lex.token == XML_TOK_OPEN_TAG) {
            p->current = NULL;
            return parse_element(p);
        }

        if (p->lex.token == XML_TOK_TEXT && p->lex.token_sb.len == 0) {
            continue; /* Empty text (whitespace only) */
        }

        parser_error(p, "Expected root element");
        return false;
    }
}

/* ============================================================================
 * PUBLIC API
 * ============================================================================ */

openmc_xml_doc_t* openmc_xml_parse_string(const char* xml_string, size_t length) {
    if (!xml_string) return NULL;

    if (length == 0) length = strlen(xml_string);

    openmc_xml_doc_t* doc = (openmc_xml_doc_t*)calloc(1, sizeof(openmc_xml_doc_t));
    if (!doc) return NULL;

    if (!arena_init(&doc->arena)) {
        free(doc);
        return NULL;
    }

    xml_parser_t parser;
    lexer_init(&parser.lex, xml_string, length, &doc->arena);
    parser.doc = doc;
    parser.current = NULL;
    parser.in_tag = false;

    if (!parse_document(&parser)) {
        /* Error already set in doc */
        return doc; /* Return doc so caller can get error message */
    }

    return doc;
}

openmc_xml_doc_t* openmc_xml_parse_file(const char* filename) {
    if (!filename) return NULL;

    FILE* f = fopen(filename, "rb");
    if (!f) {
        openmc_xml_doc_t* doc = (openmc_xml_doc_t*)calloc(1, sizeof(openmc_xml_doc_t));
        if (doc) {
            snprintf(doc->error_msg, sizeof(doc->error_msg),
                     "Failed to open file: %s", filename);
        }
        return doc;
    }

    /* Get file size */
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    fseek(f, 0, SEEK_SET);

    if (size <= 0) {
        fclose(f);
        openmc_xml_doc_t* doc = (openmc_xml_doc_t*)calloc(1, sizeof(openmc_xml_doc_t));
        if (doc) {
            snprintf(doc->error_msg, sizeof(doc->error_msg),
                     "Empty or invalid file: %s", filename);
        }
        return doc;
    }

    char* content = (char*)malloc(size + 1);
    if (!content) {
        fclose(f);
        return NULL;
    }

    size_t read = fread(content, 1, size, f);
    fclose(f);
    content[read] = '\0';

    openmc_xml_doc_t* doc = openmc_xml_parse_string(content, read);
    free(content);

    return doc;
}

void openmc_xml_doc_free(openmc_xml_doc_t* doc) {
    if (!doc) return;
    arena_free(&doc->arena);
    free(doc);
}

const char* openmc_xml_get_error(const openmc_xml_doc_t* doc) {
    if (!doc) return "NULL document";
    if (doc->error_msg[0]) return doc->error_msg;
    return NULL;
}

/* ============================================================================
 * QUERY API
 * ============================================================================ */

openmc_xml_element_t* openmc_xml_find_child(const openmc_xml_element_t* elem,
                                             const char* tag_name) {
    if (!elem || !tag_name) return NULL;

    for (size_t i = 0; i < elem->child_count; i++) {
        if (strcmp(elem->children[i]->tag_name, tag_name) == 0) {
            return elem->children[i];
        }
    }
    return NULL;
}

size_t openmc_xml_find_children(const openmc_xml_element_t* elem,
                                 const char* tag_name,
                                 openmc_xml_element_t** out_children,
                                 size_t max_count) {
    if (!elem || !tag_name || !out_children) return 0;

    size_t found = 0;
    for (size_t i = 0; i < elem->child_count && found < max_count; i++) {
        if (strcmp(elem->children[i]->tag_name, tag_name) == 0) {
            out_children[found++] = elem->children[i];
        }
    }
    return found;
}

size_t openmc_xml_count_children(const openmc_xml_element_t* elem,
                                  const char* tag_name) {
    if (!elem || !tag_name) return 0;

    size_t count = 0;
    for (size_t i = 0; i < elem->child_count; i++) {
        if (strcmp(elem->children[i]->tag_name, tag_name) == 0) {
            count++;
        }
    }
    return count;
}

const char* openmc_xml_get_attr(const openmc_xml_element_t* elem,
                                 const char* attr_name) {
    if (!elem || !attr_name) return NULL;

    for (size_t i = 0; i < elem->attr_count; i++) {
        if (strcmp(elem->attrs[i].name, attr_name) == 0) {
            return elem->attrs[i].value;
        }
    }
    return NULL;
}

int openmc_xml_get_attr_int(const openmc_xml_element_t* elem,
                             const char* attr_name,
                             int default_value) {
    const char* val = openmc_xml_get_attr(elem, attr_name);
    if (!val) return default_value;
    return (int)strtol(val, NULL, 10);
}

double openmc_xml_get_attr_double(const openmc_xml_element_t* elem,
                                   const char* attr_name,
                                   double default_value) {
    const char* val = openmc_xml_get_attr(elem, attr_name);
    if (!val) return default_value;
    return strtod(val, NULL);
}

size_t openmc_parse_doubles(const char* str, double* out_values, size_t max_values) {
    if (!str || !out_values || max_values == 0) return 0;

    size_t count = 0;
    const char* p = str;

    while (*p && count < max_values) {
        /* Skip whitespace and commas */
        while (*p && (isspace((unsigned char)*p) || *p == ',')) p++;
        if (!*p) break;

        char* end;
        double val = strtod(p, &end);
        if (end == p) break; /* No valid number */

        out_values[count++] = val;
        p = end;
    }

    return count;
}

size_t openmc_parse_ints(const char* str, int* out_values, size_t max_values) {
    if (!str || !out_values || max_values == 0) return 0;

    size_t count = 0;
    const char* p = str;

    while (*p && count < max_values) {
        while (*p && (isspace((unsigned char)*p) || *p == ',')) p++;
        if (!*p) break;

        char* end;
        long val = strtol(p, &end, 10);
        if (end == p) break;

        out_values[count++] = (int)val;
        p = end;
    }

    return count;
}
