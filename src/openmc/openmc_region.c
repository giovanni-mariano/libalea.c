// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file openmc_region.c
 * @brief Parser for OpenMC region expressions
 *
 * Grammar (in order of precedence, lowest to highest):
 *   region     := union_expr
 *   union_expr := inter_expr ('|' inter_expr)*
 *   inter_expr := primary (primary)*
 *   primary    := '(' region ')' | '~' primary | surface_ref
 *   surface_ref := ['-']? INTEGER
 *
 * Intersection is implicit (space-separated), union uses '|', complement uses '~'.
 */

#include "openmc_region.h"
#include "core/alea_ops.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* ============================================================================
 * LEXER
 * ============================================================================ */

typedef enum {
    REG_TOK_EOF,
    REG_TOK_ERROR,
    REG_TOK_SURF_INSIDE,    /* -N (negative half-space) */
    REG_TOK_SURF_OUTSIDE,   /* N (positive half-space) */
    REG_TOK_UNION,          /* | */
    REG_TOK_COMPLEMENT,     /* ~ */
    REG_TOK_LPAREN,         /* ( */
    REG_TOK_RPAREN          /* ) */
} region_token_t;

typedef struct {
    const char* src;
    size_t len;
    size_t pos;

    region_token_t token;
    int surface_id;         /* For SURF_INSIDE/SURF_OUTSIDE */

    char error_msg[128];
} region_lexer_t;

static void reg_lex_init(region_lexer_t* lex, const char* src) {
    lex->src = src;
    lex->len = src ? strlen(src) : 0;
    lex->pos = 0;
    lex->token = REG_TOK_EOF;
    lex->surface_id = 0;
    lex->error_msg[0] = '\0';
}

static char reg_peek(region_lexer_t* lex) {
    if (lex->pos >= lex->len) return '\0';
    return lex->src[lex->pos];
}

static char reg_advance(region_lexer_t* lex) {
    if (lex->pos >= lex->len) return '\0';
    return lex->src[lex->pos++];
}

static void reg_skip_ws(region_lexer_t* lex) {
    while (isspace((unsigned char)reg_peek(lex))) {
        reg_advance(lex);
    }
}

static void reg_next_token(region_lexer_t* lex) {
    reg_skip_ws(lex);

    char c = reg_peek(lex);

    if (c == '\0') {
        lex->token = REG_TOK_EOF;
        return;
    }

    if (c == '|') {
        reg_advance(lex);
        lex->token = REG_TOK_UNION;
        return;
    }

    if (c == '~') {
        reg_advance(lex);
        lex->token = REG_TOK_COMPLEMENT;
        return;
    }

    if (c == '(') {
        reg_advance(lex);
        lex->token = REG_TOK_LPAREN;
        return;
    }

    if (c == ')') {
        reg_advance(lex);
        lex->token = REG_TOK_RPAREN;
        return;
    }

    /* Surface reference: optional '-' followed by integer */
    if (c == '-' || isdigit((unsigned char)c)) {
        bool negative = false;
        if (c == '-') {
            negative = true;
            reg_advance(lex);
            reg_skip_ws(lex);
        }

        if (!isdigit((unsigned char)reg_peek(lex))) {
            snprintf(lex->error_msg, sizeof(lex->error_msg),
                     "Expected surface number after '-' at position %zu", lex->pos);
            lex->token = REG_TOK_ERROR;
            return;
        }

        /* Parse integer */
        int id = 0;
        while (isdigit((unsigned char)reg_peek(lex))) {
            id = id * 10 + (reg_advance(lex) - '0');
        }

        lex->surface_id = id;
        lex->token = negative ? REG_TOK_SURF_INSIDE : REG_TOK_SURF_OUTSIDE;
        return;
    }

    snprintf(lex->error_msg, sizeof(lex->error_msg),
             "Unexpected character '%c' at position %zu", c, lex->pos);
    lex->token = REG_TOK_ERROR;
}

/* ============================================================================
 * PARSER
 * ============================================================================ */

typedef struct {
    region_lexer_t lex;
    openmc_region_ctx_t* ctx;
} region_parser_t;

static void reg_parser_init(region_parser_t* p, openmc_region_ctx_t* ctx, const char* expr) {
    reg_lex_init(&p->lex, expr);
    p->ctx = ctx;
    reg_next_token(&p->lex);
}

static alea_node_id_t parse_union_expr(region_parser_t* p);
static alea_node_id_t parse_inter_expr(region_parser_t* p);
static alea_node_id_t parse_primary(region_parser_t* p);

/**
 * Get the pre-created node for a surface reference.
 * Returns the appropriate pos/neg sense node directly — no new nodes created.
 */
static alea_node_id_t get_surface_node(region_parser_t* p, int surface_id, bool inside) {
    openmc_region_ctx_t* ctx = p->ctx;

    /* Look up in surface map */
    if (surface_id < 0 || (size_t)surface_id >= ctx->surface_map_size) {
        snprintf(ctx->error_msg, sizeof(ctx->error_msg),
                 "Surface %d not found", surface_id);
        return ALEA_NODE_ID_INVALID;
    }

    const openmc_surface_node_pair_t* pair = &ctx->surface_node_map[surface_id];
    if (pair->pos_node == ALEA_NODE_ID_INVALID) {
        snprintf(ctx->error_msg, sizeof(ctx->error_msg),
                 "Surface %d not registered", surface_id);
        return ALEA_NODE_ID_INVALID;
    }

    return inside ? pair->neg_node : pair->pos_node;
}

/**
 * primary := '(' region ')' | '~' primary | surface_ref
 */
static alea_node_id_t parse_primary(region_parser_t* p) {
    region_token_t tok = p->lex.token;

    if (tok == REG_TOK_LPAREN) {
        reg_next_token(&p->lex);
        alea_node_id_t node = parse_union_expr(p);
        if (node == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;

        if (p->lex.token != REG_TOK_RPAREN) {
            snprintf(p->ctx->error_msg, sizeof(p->ctx->error_msg),
                     "Expected ')' at position %zu", p->lex.pos);
            return ALEA_NODE_ID_INVALID;
        }
        reg_next_token(&p->lex);
        return node;
    }

    if (tok == REG_TOK_COMPLEMENT) {
        reg_next_token(&p->lex);
        alea_node_id_t inner = parse_primary(p);
        if (inner == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;
        return alea_create_complement(p->ctx->sys, inner);
    }

    if (tok == REG_TOK_SURF_INSIDE) {
        int id = p->lex.surface_id;
        reg_next_token(&p->lex);
        return get_surface_node(p, id, true);
    }

    if (tok == REG_TOK_SURF_OUTSIDE) {
        int id = p->lex.surface_id;
        reg_next_token(&p->lex);
        return get_surface_node(p, id, false);
    }

    if (tok == REG_TOK_ERROR) {
        snprintf(p->ctx->error_msg, sizeof(p->ctx->error_msg),
                 "%s", p->lex.error_msg);
    } else {
        snprintf(p->ctx->error_msg, sizeof(p->ctx->error_msg),
                 "Unexpected token at position %zu", p->lex.pos);
    }
    return ALEA_NODE_ID_INVALID;
}

/**
 * Check if current token can start a primary expression
 */
static bool can_start_primary(region_token_t tok) {
    return tok == REG_TOK_LPAREN ||
           tok == REG_TOK_COMPLEMENT ||
           tok == REG_TOK_SURF_INSIDE ||
           tok == REG_TOK_SURF_OUTSIDE;
}

/**
 * inter_expr := primary (primary)*
 * Space-separated terms are implicitly intersected.
 */
static alea_node_id_t parse_inter_expr(region_parser_t* p) {
    alea_node_id_t left = parse_primary(p);
    if (left == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;

    while (can_start_primary(p->lex.token)) {
        alea_node_id_t right = parse_primary(p);
        if (right == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;
        left = alea_create_intersection(p->ctx->sys, left, right);
        if (left == ALEA_NODE_ID_INVALID) {
            snprintf(p->ctx->error_msg, sizeof(p->ctx->error_msg),
                     "Failed to create intersection");
            return ALEA_NODE_ID_INVALID;
        }
    }

    return left;
}

/**
 * union_expr := inter_expr ('|' inter_expr)*
 */
static alea_node_id_t parse_union_expr(region_parser_t* p) {
    alea_node_id_t left = parse_inter_expr(p);
    if (left == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;

    while (p->lex.token == REG_TOK_UNION) {
        reg_next_token(&p->lex);
        alea_node_id_t right = parse_inter_expr(p);
        if (right == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;
        left = alea_create_union(p->ctx->sys, left, right);
        if (left == ALEA_NODE_ID_INVALID) {
            snprintf(p->ctx->error_msg, sizeof(p->ctx->error_msg),
                     "Failed to create union");
            return ALEA_NODE_ID_INVALID;
        }
    }

    return left;
}

/* ============================================================================
 * PUBLIC API
 * ============================================================================ */

void openmc_region_ctx_init(openmc_region_ctx_t* ctx, alea_system_t* sys, arena_t* arena) {
    memset(ctx, 0, sizeof(*ctx));
    ctx->sys = sys;
    ctx->arena = arena;
    ctx->surface_node_map = NULL;
    ctx->surface_map_size = 0;
}

int openmc_region_register_surface(openmc_region_ctx_t* ctx, int surface_id,
                                   alea_node_id_t pos_node, alea_node_id_t neg_node) {
    if (!ctx || surface_id < 0) return -1;

    /* Grow map if needed */
    size_t needed = (size_t)surface_id + 1;
    if (needed > ctx->surface_map_size) {
        size_t new_size = needed * 2;
        if (new_size < 64) new_size = 64;

        openmc_surface_node_pair_t* new_map;
        if (ctx->arena) {
            new_map = (openmc_surface_node_pair_t*)arena_alloc(
                ctx->arena, new_size * sizeof(openmc_surface_node_pair_t));
        } else {
            new_map = (openmc_surface_node_pair_t*)realloc(
                ctx->surface_node_map, new_size * sizeof(openmc_surface_node_pair_t));
        }
        if (!new_map) return -1;

        /* Initialize new entries as unregistered */
        for (size_t i = ctx->surface_map_size; i < new_size; i++) {
            new_map[i].pos_node = ALEA_NODE_ID_INVALID;
            new_map[i].neg_node = ALEA_NODE_ID_INVALID;
        }

        /* Copy old data if using arena (realloc already copies) */
        if (ctx->arena && ctx->surface_node_map) {
            memcpy(new_map, ctx->surface_node_map,
                   ctx->surface_map_size * sizeof(openmc_surface_node_pair_t));
        }

        ctx->surface_node_map = new_map;
        ctx->surface_map_size = new_size;
    }

    ctx->surface_node_map[surface_id].pos_node = pos_node;
    ctx->surface_node_map[surface_id].neg_node = neg_node;
    return 0;
}

alea_node_id_t openmc_parse_region(openmc_region_ctx_t* ctx, const char* region_expr) {
    if (!ctx || !ctx->sys) {
        if (ctx) strcpy(ctx->error_msg, "Invalid context");
        return ALEA_NODE_ID_INVALID;
    }

    if (!region_expr || !*region_expr) {
        strcpy(ctx->error_msg, "Empty region expression");
        return ALEA_NODE_ID_INVALID;
    }

    region_parser_t parser;
    reg_parser_init(&parser, ctx, region_expr);

    alea_node_id_t root = parse_union_expr(&parser);

    if (root != ALEA_NODE_ID_INVALID && parser.lex.token != REG_TOK_EOF) {
        snprintf(ctx->error_msg, sizeof(ctx->error_msg),
                 "Unexpected content after expression at position %zu", parser.lex.pos);
        return ALEA_NODE_ID_INVALID;
    }

    return root;
}

const char* openmc_region_get_error(const openmc_region_ctx_t* ctx) {
    if (!ctx) return "NULL context";
    if (ctx->error_msg[0]) return ctx->error_msg;
    return NULL;
}
