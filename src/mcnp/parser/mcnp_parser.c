// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "mcnp_parser.h"
#include "mcnp_lexer.h"
#include "util/alea_log.h"
#include "util/str_builder.h"
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>

// Platform-specific includes for memory-mapped files
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#endif
#include "util/compat.h"

/* ========================================================================== */
/* INTERNAL STATE AND HELPERS                                                 */
/* ========================================================================== */

// Represents the current parsing block (cells, surfaces, or data)
typedef enum {
    STATE_TITLE_CARD,
    STATE_CELL_BLOCK,
    STATE_SURFACE_BLOCK,
    STATE_DATA_BLOCK,
    STATE_DONE
} parser_state_t;

// A structure to hold the memory-mapped file data
typedef struct {
    const char* content; // Pointer to the file content in memory
    size_t size;         // Size of the file

#ifdef _WIN32
    HANDLE hFile;
    HANDLE hMapping;
#else
    int fd;
#endif
} mapped_file_t;

// Forward declarations for internal parsing functions
static int parse_cell_card(mcnp_context_t* ctx, const char* line, size_t len);
static int parse_surface_card(mcnp_context_t* ctx, const char* line, size_t len);
static int parse_material_card(mcnp_context_t* ctx, const char* line, size_t len);
static int parse_transform_card(mcnp_context_t* ctx, const char* line, size_t len);

/* ========================================================================== */
/* FILE MEMORY-MAPPING (PLATFORM-SPECIFIC)                                    */
/* ========================================================================== */

static int map_file(const char* filepath, mapped_file_t* mf) {
#ifdef _WIN32
    mf->hFile = CreateFileA(filepath, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (mf->hFile == INVALID_HANDLE_VALUE) return 0;

    mf->size = GetFileSize(mf->hFile, NULL);
    if (mf->size == INVALID_FILE_SIZE) { CloseHandle(mf->hFile); return 0; }

    mf->hMapping = CreateFileMapping(mf->hFile, NULL, PAGE_READONLY, 0, 0, NULL);
    if (mf->hMapping == NULL) { CloseHandle(mf->hFile); return 0; }

    mf->content = (const char*)MapViewOfFile(mf->hMapping, FILE_MAP_READ, 0, 0, 0);
    if (mf->content == NULL) { CloseHandle(mf->hMapping); CloseHandle(mf->hFile); return 0; }
#else
    mf->fd = open(filepath, O_RDONLY);
    if (mf->fd == -1) return 0;

    struct stat s;
    if (fstat(mf->fd, &s) == -1) { close(mf->fd); return 0; }
    mf->size = s.st_size;

    mf->content = (const char*)mmap(NULL, mf->size, PROT_READ, MAP_PRIVATE, mf->fd, 0);
    if (mf->content == MAP_FAILED) { close(mf->fd); return 0; }
#endif
    return 1;
}

static void unmap_file(mapped_file_t* mf) {
#ifdef _WIN32
    UnmapViewOfFile(mf->content);
    CloseHandle(mf->hMapping);
    CloseHandle(mf->hFile);
#else
    munmap((void*)mf->content, mf->size);
    close(mf->fd);
#endif
}

/* ========================================================================== */
/* DYNAMIC ARRAY MANAGEMENT (USING THE ARENA)                                 */
/* ========================================================================== */

// Generic function to ensure capacity in our dynamic pointer arrays
static int ensure_capacity(mcnp_context_t* ctx, void*** array_ptr, size_t* capacity, size_t element_size) {
    if (*capacity == 0) {
        *capacity = 64; // Initial capacity
        *array_ptr = arena_alloc(&ctx->arena, *capacity * element_size);
    } else {
        size_t new_capacity = *capacity * 2;
        void** new_array = arena_alloc(&ctx->arena, new_capacity * element_size);
        if (new_array) {
            memcpy(new_array, *array_ptr, *capacity * element_size);
            *array_ptr = new_array;
            *capacity = new_capacity;
        } else {
            ALEA_LOG_ERROR("Out of memory growing array (capacity %zu)", *capacity);
            return 0;
        }
    }
    return *array_ptr != NULL;
}

/**
 * @brief Creates a new MCNP context with a dedicated memory arena.
 *
 * @param filename The name of the file being parsed (can be NULL).
 * @return A pointer to the new context, or NULL on failure.
 */
mcnp_context_t* mcnp_context_create(const char* filename) {
    // The context itself is a small allocation
    mcnp_context_t* ctx = (mcnp_context_t*)malloc(sizeof(mcnp_context_t));
    if (!ctx) return NULL;
    memset(ctx, 0, sizeof(mcnp_context_t));

    // Initialize the main arena where ALL data will live
    if (!arena_init(&ctx->arena)) {
        free(ctx);
        return NULL;
    }

    // Initialize a scratch arena for temporary per-card allocations.
    // 4 MB initial block is enough for any single logical card.
    if (!arena_init_with_size(&ctx->scratch, 4 * 1024 * 1024)) {
        arena_free(&ctx->arena);
        free(ctx);
        return NULL;
    }

    // Allocate the initial filename and title card from the new arena
    if (filename) {
        ctx->source_filename = arena_strdup(&ctx->arena, filename);
    }
    ctx->title_card = arena_strdup(&ctx->arena, "Default Title Card - To be replaced by parser");
    //ctx->parse_time = time(NULL);

    // Note: The arrays for cells/surfaces/materials will be allocated
    // from the arena on-demand by the parser.
    
    ALEA_LOG_INFO("Created MCNP context for file: %s\n", filename ? filename : "unknown");
   
    return ctx;
}

/**
 * @brief Destroys the MCNP context and frees all associated memory.
 *
 * Thanks to the arena, this is incredibly simple and fast.
 * @param ctx The context to destroy.
 */
void mcnp_context_destroy(mcnp_context_t* ctx) {
    if (!ctx) return;
    
    ALEA_LOG_INFO("Destroying MCNP context for file: %s\n", ctx->source_filename ? ctx->source_filename : "unknown");
    
    // This one call frees the memory for the arena buffer, which contains
    // ALL cells, surfaces, materials, strings, and comments. No loops needed.
    arena_free(&ctx->arena);
    arena_free(&ctx->scratch);

    // Free the context struct itself.
    free(ctx);
}

/* ========================================================================== */
/* MAIN PARSER LOGIC                                                          */
/* ========================================================================== */

int mcnp_parse_file(const char* filename, mcnp_context_t** out_context) {
    mapped_file_t mf = {0};
    if (!map_file(filename, &mf)) {
        ALEA_LOG_ERROR("ERROR: Could not open or map file: %s\n", filename);
        return 0;
    }

    mcnp_context_t* ctx = mcnp_context_create(filename);
    if (!ctx) {
        unmap_file(&mf);
        return 0;
    }

    parser_state_t state = STATE_TITLE_CARD;
    const char* current = mf.content;
    const char* end = mf.content + mf.size;
    int line_num = 0;
    int had_fatal_error = 0;

    // --- String builders backed by scratch arena (reclaimed between cards) ---
    str_builder_t line_sb;
    str_builder_init(&line_sb, &ctx->scratch, 4096);

    str_builder_t comment_sb;       // Comments before the current card
    str_builder_init(&comment_sb, &ctx->scratch, 1024);
    str_builder_t peek_comment_sb;  // Comments consumed during look-ahead (for next card)
    str_builder_init(&peek_comment_sb, &ctx->scratch, 512);

    while (current < end) {
        const char* line_start = current;
        const char* line_end = strchr(line_start, '\n');
        if (line_end == NULL) line_end = end;
        line_num++;


        const char* trimmed_start = line_start;
        while (trimmed_start < line_end && isspace((unsigned char)*trimmed_start)) {
            trimmed_start++;
        }

        int is_blank = (trimmed_start == line_end);
        int is_comment = (trimmed_start < line_end && tolower((unsigned char)*trimmed_start) == 'c' && isspace((unsigned char)*(trimmed_start + 1)));


        if (state == STATE_TITLE_CARD) {
            if (!is_blank) {
                size_t title_len = line_end - trimmed_start;
                const char* title_end = line_end;
                mcnp_lexer_trim_trailing_whitespace(&title_end, trimmed_start);
                title_len = title_end - trimmed_start;

                ctx->title_card = (char*)arena_alloc(&ctx->arena, title_len + 1);
                memcpy(ctx->title_card, trimmed_start, title_len);
                ctx->title_card[title_len] = '\0';
                state = STATE_CELL_BLOCK;
            }
            current = line_end + 1;
            continue;
        }

        if (is_blank) {
            if (state == STATE_CELL_BLOCK) state = STATE_SURFACE_BLOCK;
            else if (state == STATE_SURFACE_BLOCK) state = STATE_DATA_BLOCK;
            str_builder_reset(&comment_sb);
            str_builder_reset(&peek_comment_sb);
            current = line_end + 1;
            continue;
        }
        if (is_comment) {
            // Accumulate "C" comment lines (preserve full line from trimmed_start)
            if (state == STATE_CELL_BLOCK) {
                str_builder_write(&comment_sb, trimmed_start, line_end - trimmed_start);
                str_builder_putc(&comment_sb, '\n');
            }
            current = line_end + 1;
            continue;
        }

        size_t initial_len = line_end - line_start;

        // Track the last inline $ comment across all physical lines of this card
        const char* last_dollar_start = NULL;
        const char* last_dollar_end = NULL;

        // Strip inline comment from this line (if any)
        const char* comment_pos = (const char*)memchr(line_start, '$', initial_len);
        if (comment_pos) {
            initial_len = comment_pos - line_start;
            last_dollar_start = comment_pos + 1;
            last_dollar_end = line_end;
        }

        str_builder_reset(&line_sb);
        str_builder_write(&line_sb, line_start, initial_len);
        
        const char* peek_pos = line_end + 1;

        while (peek_pos < end) {
            const char* next_line_start = peek_pos;
            const char* next_line_end = strchr(next_line_start, '\n');
            if (next_line_end == NULL) next_line_end = end;

            // Classify the peeked line
            const char* peek_trimmed_start = next_line_start;
            while(peek_trimmed_start < next_line_end && isspace((unsigned char)*peek_trimmed_start)) {
                peek_trimmed_start++;
            }
            (void)(peek_trimmed_start == next_line_end);  /* peek_is_blank - reserved for future use */
            int peek_is_comment = (peek_trimmed_start < next_line_end && tolower((unsigned char)*peek_trimmed_start) == 'c' && isspace((unsigned char)*(peek_trimmed_start+1)));

            // If the peeked line is a comment, accumulate for the *next* card
            // and continue looking for continuations of the current card.
            if (peek_is_comment) {
                if (state == STATE_CELL_BLOCK) {
                    str_builder_write(&peek_comment_sb, peek_trimmed_start,
                                      next_line_end - peek_trimmed_start);
                    str_builder_putc(&peek_comment_sb, '\n');
                }
                peek_pos = next_line_end + 1;
                line_num++;
                continue;
            }
                        
            // Now, check if the non-comment is a continuation
            const char* prev_line_logical_end = line_sb.buf + line_sb.len;
            mcnp_lexer_trim_trailing_whitespace(&prev_line_logical_end, line_sb.buf);
            int ends_with_ampersand = (prev_line_logical_end > line_sb.buf && *(prev_line_logical_end-1) == '&');
            if (ends_with_ampersand) {
                line_sb.len = (prev_line_logical_end - line_sb.buf) - 1;
            }

            int indent = 0;
            const char* p = next_line_start;
            while(indent < 5 && p < next_line_end && isspace((unsigned char)*p)) {
                indent++; p++;
            }
            int is_continuation = (indent == 5);

            if (is_continuation || ends_with_ampersand) {
                // Comments peeked before this continuation are mid-card; discard them
                str_builder_reset(&peek_comment_sb);

                size_t continuation_len = next_line_end - next_line_start;

                // Strip inline comment from continuation line (if any)
                const char* cont_comment_pos = (const char*)memchr(next_line_start, '$', continuation_len);
                if (cont_comment_pos) {
                    continuation_len = cont_comment_pos - next_line_start;
                    last_dollar_start = cont_comment_pos + 1;
                    last_dollar_end = next_line_end;
                } else {
                    // This continuation has no $, so it becomes the "last line"
                    // and clears any previous $ capture
                    last_dollar_start = NULL;
                    last_dollar_end = NULL;
                }

                // Append a space to replace the newline/&
                str_builder_putc(&line_sb, ' ');

                // Find start of actual content on the continued line
                const char* content_start = next_line_start;
                while(content_start < next_line_end && isspace((unsigned char)*content_start)) content_start++;
                continuation_len = next_line_end - content_start;

                // Re-check if comment marker appears after leading whitespace
                if (content_start < next_line_end) {
                    const char* adj_comment_pos = (const char*)memchr(content_start, '$', continuation_len);
                    if (adj_comment_pos) {
                        continuation_len = adj_comment_pos - content_start;
                    }
                }

                str_builder_write(&line_sb, content_start, continuation_len);

                peek_pos = next_line_end + 1;
                line_num++;
            } else { break; }
        }
        
        str_builder_finish(&line_sb);

        int success = 1;
        const char* logical_line_ptr = line_sb.buf;
        while(isspace((unsigned char)*logical_line_ptr)) logical_line_ptr++;

        switch (state) {
            case STATE_CELL_BLOCK: if (!parse_cell_card(ctx, logical_line_ptr, strlen(logical_line_ptr))) success = 0; break;
            case STATE_SURFACE_BLOCK: if (!parse_surface_card(ctx, logical_line_ptr, strlen(logical_line_ptr))) success = 0; break;
            case STATE_DATA_BLOCK: {
                char c = tolower((unsigned char)*logical_line_ptr);
                char c2 = logical_line_ptr[0];
                /* Material cards: Mn */
                if (c == 'm' && isdigit((unsigned char)logical_line_ptr[1])) {
                    if (!parse_material_card(ctx, logical_line_ptr, strlen(logical_line_ptr))) success = 0;
                }
                /* Transform cards: TRn or *TRn */
                else if ((c == 't' && tolower((unsigned char)logical_line_ptr[1]) == 'r') ||
                         (c2 == '*' && tolower((unsigned char)logical_line_ptr[1]) == 't' && 
                          tolower((unsigned char)logical_line_ptr[2]) == 'r')) {
                    if (!parse_transform_card(ctx, logical_line_ptr, strlen(logical_line_ptr))) success = 0;
                }
                break;
            }
            default: break;
        }

        // Attach accumulated comments to the just-parsed cell
        if (state == STATE_CELL_BLOCK && success && ctx->cell_count > 0) {
            mcnp_cell_t* last_cell = ctx->cells[ctx->cell_count - 1];
            if (comment_sb.len > 0) {
                str_builder_finish(&comment_sb);
                last_cell->comments = arena_strdup(&ctx->arena, comment_sb.buf);
            }
            if (last_dollar_start) {
                // Trim leading/trailing whitespace from the $ comment text
                while (last_dollar_start < last_dollar_end && isspace((unsigned char)*last_dollar_start))
                    last_dollar_start++;
                const char* dollar_trim_end = last_dollar_end;
                while (dollar_trim_end > last_dollar_start && isspace((unsigned char)*(dollar_trim_end - 1)))
                    dollar_trim_end--;
                if (dollar_trim_end > last_dollar_start) {
                    size_t dlen = dollar_trim_end - last_dollar_start;
                    char* buf = (char*)arena_alloc(&ctx->arena, dlen + 1);
                    memcpy(buf, last_dollar_start, dlen);
                    buf[dlen] = '\0';
                    last_cell->inline_comment = buf;
                }
            }
        }
        str_builder_reset(&comment_sb);
        // Transfer any peeked comments (consumed during look-ahead) for the next card
        if (peek_comment_sb.len > 0) {
            str_builder_write(&comment_sb, peek_comment_sb.buf, peek_comment_sb.len);
            str_builder_reset(&peek_comment_sb);
        }

        // Reset the scratch arena to reclaim per-card temporaries (geo_sb, param_sb,
        // line_copy, etc.).  Peek comments now in comment_sb must survive the reset.
        {
            char  peek_static[4096];
            char* peek_save = NULL;
            size_t peek_len = comment_sb.len;
            if (peek_len > 0) {
                if (peek_len < sizeof(peek_static)) {
                    memcpy(peek_static, comment_sb.buf, peek_len + 1);
                    peek_save = peek_static;
                } else {
                    peek_save = (char*)malloc(peek_len + 1);
                    if (peek_save) memcpy(peek_save, comment_sb.buf, peek_len + 1);
                    else peek_len = 0; /* malloc failed: drop comments rather than crash */
                }
            }

            arena_reset(&ctx->scratch);
            str_builder_init(&line_sb, &ctx->scratch, 4096);
            str_builder_init(&comment_sb, &ctx->scratch, 1024);
            str_builder_init(&peek_comment_sb, &ctx->scratch, 512);

            if (peek_len > 0 && peek_save) {
                str_builder_write(&comment_sb, peek_save, peek_len);
                if (peek_save != peek_static) free(peek_save);
            }
        }

        if (!success) {
            ALEA_LOG_ERROR("ERROR: Parsing failed on logical line starting near line %d.\n", line_num);
            had_fatal_error = 1;
            break;
        }
        current = peek_pos;
    }

    unmap_file(&mf);
    if (had_fatal_error) {
        mcnp_context_destroy(ctx);
        *out_context = NULL;
        return 0;
    }
    *out_context = ctx;
    ALEA_LOG_INFO("Parsing complete. Cells: %zu, Surfaces: %zu, Materials: %zu\n", ctx->cell_count, ctx->surface_count, ctx->material_count);
    return 1;
}

/* ========================================================================== */
/* CARD-SPECIFIC PARSING ROUTINES                                             */
/* ========================================================================== */

static int parse_cell_card(mcnp_context_t* ctx, const char* line, size_t len) {
    if (ctx->cell_count >= ctx->cell_capacity) {
        if (!ensure_capacity(ctx, (void***)&ctx->cells, &ctx->cell_capacity, sizeof(mcnp_cell_t*))) return 0;
    }

    //mcnp_cell_t* cell = arena_alloc_struct(&ctx->arena, mcnp_cell_t);
    mcnp_cell_t* cell = (mcnp_cell_t*)arena_alloc(&ctx->arena, sizeof(mcnp_cell_t));
    if (!cell) { ALEA_LOG_ERROR("Out of memory allocating cell"); return 0; }
    memset(cell, 0, sizeof(*cell));
    
    (void)len; /* len was used for line_copy; lexer reads via pointer only */

    // Comments have already been stripped during line assembly.
    // The lexer only advances the cursor pointer without modifying the string,
    // so we can read directly from the caller's buffer — no copy needed.
    const char* cursor = line;
    const char* token_start = NULL;
    size_t token_len;
    char token_buf[128];

    // --- ID, Material, Density using lexer ---
    token_len = mcnp_lexer_get_next_token(&cursor, &token_start);
    if (token_len == 0) { ALEA_LOG_ERROR("Truncated cell card: no cell ID"); return 0; }
    if (token_len >= sizeof(token_buf)) { ALEA_LOG_WARN("Cell card: token too long (%zu)", token_len); return 0; }
    memcpy(token_buf, token_start, token_len); token_buf[token_len] = '\0';
    cell->cell_id = strtol(token_buf, NULL, 10);

    token_len = mcnp_lexer_get_next_token(&cursor, &token_start);
    if (token_len == 0) { ALEA_LOG_ERROR("Cell %d: missing material ID", cell->cell_id); return 0; }
    if (token_len >= sizeof(token_buf)) { ALEA_LOG_WARN("Cell %d: material token too long (%zu)", cell->cell_id, token_len); return 0; }
    memcpy(token_buf, token_start, token_len); token_buf[token_len] = '\0';

    // Detect LIKE BUT cells: "cell_id LIKE ref_id BUT ..."
    // "LIKE" appears where material_id would normally be
    int is_like_card = (token_len == 4 && alea_strncasecmp(token_buf, "LIKE", 4) == 0);

    if (!is_like_card) {
        cell->material_id = strtol(token_buf, NULL, 10);
    } else {
        cell->material_id = 0;
    }

    if (!is_like_card && cell->material_id != 0) {
        token_len = mcnp_lexer_get_next_token(&cursor, &token_start);
        if (token_len == 0) { ALEA_LOG_ERROR("Cell %d: missing density for material %d", cell->cell_id, cell->material_id); return 0; }
        if (token_len >= sizeof(token_buf)) { ALEA_LOG_WARN("Cell %d: density token too long (%zu)", cell->cell_id, token_len); return 0; }
        memcpy(token_buf, token_start, token_len); token_buf[token_len] = '\0';
        cell->density = strtod(token_buf, NULL);
    } else {
        cell->density = 0.0;
    }

    // --- Separate Geometry from Parameters using lexer ---
    str_builder_t geo_sb, param_sb;
    str_builder_init(&geo_sb, &ctx->scratch, 1024);
    str_builder_init(&param_sb, &ctx->scratch, 1024);
    int in_parameters = 0;  // Flag: once we hit a parameter keyword, everything is parameters

    // For LIKE cells, prepend "LIKE" to geometry since it was consumed as the material token
    if (is_like_card) {
        str_builder_write(&geo_sb, "LIKE", 4);
    }

    // Loop through all remaining tokens using lexer
    while ((token_len = mcnp_lexer_get_next_token(&cursor, &token_start)) > 0) {
        // Check if this token is a parameter keyword
        if (!in_parameters && mcnp_lexer_is_parameter_keyword(token_start, token_len)) {
            in_parameters = 1;  // Switch to parameter mode permanently
        }

        str_builder_t* sb = in_parameters ? &param_sb : &geo_sb;
        if (sb->len > 0) str_builder_putc(sb, ' ');
        str_builder_write(sb, token_start, token_len);
    }

    str_builder_finish(&geo_sb);
    str_builder_finish(&param_sb);

    // Finalize strings in the arena
    if (geo_sb.len > 0) {
        cell->geometry_definition = arena_strdup(&ctx->arena, geo_sb.buf);
    }
    if (param_sb.len > 0) {
        cell->parameters = arena_strdup(&ctx->arena, param_sb.buf);
    }

    ctx->cells[ctx->cell_count++] = cell;
    return 1;
}

static int parse_surface_card(mcnp_context_t* ctx, const char* line, size_t len) {
    if (ctx->surface_count >= ctx->surface_capacity) {
        if (!ensure_capacity(ctx, (void***)&ctx->surfaces, &ctx->surface_capacity, sizeof(mcnp_surface_t*))) return 0;
    }

    //mcnp_surface_t* surf = arena_alloc_struct(&ctx->arena, mcnp_surface_t);
    mcnp_surface_t* surf = (mcnp_surface_t*)arena_alloc(&ctx->arena, sizeof(mcnp_surface_t));
    if (!surf) { ALEA_LOG_ERROR("Out of memory allocating surface"); return 0; }
    
    // Set defaults
    surf->boundary_type = SURFACE_TRANSMISSIVE;
    surf->transform_id = 0;
    surf->periodic_surface_id = 0;

    const char* comment_start_ptr = memchr(line, '$', len);
    if (comment_start_ptr) {
        surf->comments = arena_strdup(&ctx->arena, comment_start_ptr + 1);
        len = comment_start_ptr - line;
    }

    char* line_copy = (char*)arena_alloc(&ctx->scratch, len + 1);
    if (!line_copy) { ALEA_LOG_ERROR("Out of memory for surface line copy"); return 0; }
    memcpy(line_copy, line, len);
    line_copy[len] = '\0';

    const char* cursor = line_copy;
    const char* token_start = NULL;
    size_t token_len;
    char token_buf[128];

    // --- Token 1: J (Surface ID with optional prefix) using lexer ---
    token_len = mcnp_lexer_get_next_token(&cursor, &token_start);
    if (token_len == 0) { ALEA_LOG_ERROR("Truncated surface card: no surface ID"); return 0; }

    // Check for prefixes
    if (*token_start == '*') {
        surf->boundary_type = SURFACE_REFLECTIVE;
        token_start++; token_len--;
    } else if (*token_start == '+') {
        surf->boundary_type = SURFACE_WHITE;
        token_start++; token_len--;
    }
    
    if (token_len == 0 || token_len >= sizeof(token_buf)) {
        ALEA_LOG_ERROR("Surface card: invalid surface ID token (len=%zu)", token_len);
        return 0;
    }
    memcpy(token_buf, token_start, token_len);
    token_buf[token_len] = '\0';
    surf->surface_id = strtol(token_buf, NULL, 10);
    
    // --- Token 2: Could be N (optional) or A (mnemonic) ---
    const char* next_token_start_pos = cursor;
    token_len = mcnp_lexer_get_next_token(&cursor, &token_start);
    if (token_len == 0) { ALEA_LOG_ERROR("Surface %d: missing type mnemonic", surf->surface_id); return 0; }

    // Try to parse as number
    char* endptr;
    if (token_len >= sizeof(token_buf)) { ALEA_LOG_WARN("Surface %d: token too long (%zu)", surf->surface_id, token_len); return 0; }
    memcpy(token_buf, token_start, token_len);
    token_buf[token_len] = '\0';
    long n_value = strtol(token_buf, &endptr, 10);

    if (*endptr == '\0') {
        // Was a valid number (N parameter)
        if (n_value > 0) {
            surf->transform_id = (int)n_value;
        } else if (n_value < 0) {
            surf->boundary_type = SURFACE_PERIODIC;
            surf->periodic_surface_id = (int)(-n_value);
        }
        // Next token must be mnemonic A
        token_len = mcnp_lexer_get_next_token(&cursor, &token_start);
        if (token_len == 0) {
            ALEA_LOG_ERROR("Surface %d: missing type mnemonic", surf->surface_id);
            return 0;
        }
        surf->mnemonic = (char*)arena_alloc(&ctx->arena, token_len + 1);
        memcpy(surf->mnemonic, token_start, token_len);
        surf->mnemonic[token_len] = '\0';
    } else {
        // Was NOT a number, so it IS the mnemonic A
        cursor = next_token_start_pos; // Rewind
        token_len = mcnp_lexer_get_next_token(&cursor, &token_start);
        if (token_len == 0) {
            ALEA_LOG_ERROR("Surface %d: missing type mnemonic", surf->surface_id);
            return 0;
        }
        surf->mnemonic = (char*)arena_alloc(&ctx->arena, token_len + 1);
        memcpy(surf->mnemonic, token_start, token_len);
        surf->mnemonic[token_len] = '\0';
    }
    
    // --- Remaining text: coefficients ---
    while (*cursor && isspace((unsigned char)*cursor)) cursor++;
    if (*cursor) {
        surf->coefficients_str = arena_strdup(&ctx->arena, cursor);
    }
    
    ctx->surfaces[ctx->surface_count++] = surf;
    return 1;
}

static int parse_material_card(mcnp_context_t* ctx, const char* line, size_t len) {
    (void)len;  /* Currently unused, kept for API consistency */
    if (ctx->material_count >= ctx->material_capacity) {
        if (!ensure_capacity(ctx, (void***)&ctx->materials, &ctx->material_capacity, sizeof(mcnp_material_t*))) return 0;
    }
    
    //mcnp_material_t* mat = arena_alloc_struct(&ctx->arena, mcnp_material_t);
    mcnp_material_t* mat = (mcnp_material_t*)arena_alloc(&ctx->arena, sizeof(mcnp_material_t));
    if (!mat) { ALEA_LOG_ERROR("Out of memory allocating material"); return 0; }
    
    //mat->definition = arena_strdup(&ctx->arena, line);

    memset(mat, 0, sizeof(mcnp_material_t));
    
    /* Parse material ID from "Mn" or "mn" */
    const char* p = line;
    if (tolower((unsigned char)*p) == 'm') p++;
    mat->material_id = (int)strtol(p, (char**)&p, 10);
    
    /* Skip whitespace */
    while (*p && isspace((unsigned char)*p)) p++;
    
    /* Rest is the definition */
    if (*p) {
        mat->definition = arena_strdup(&ctx->arena, p);
    }

    ctx->materials[ctx->material_count++] = mat;
    return 1;
}


static int parse_transform_card(mcnp_context_t* ctx, const char* line, size_t len) {
    (void)len;  /* Currently unused, kept for API consistency */
    if (ctx->transform_count >= ctx->transform_capacity) {
        if (!ensure_capacity(ctx, (void***)&ctx->transforms, &ctx->transform_capacity, sizeof(mcnp_transform_t*))) return 0;
    }
    
    mcnp_transform_t* tr = (mcnp_transform_t*)arena_alloc(&ctx->arena, sizeof(mcnp_transform_t));
    if (!tr) { ALEA_LOG_ERROR("Out of memory allocating transform"); return 0; }
    
    memset(tr, 0, sizeof(mcnp_transform_t));
    
    const char* p = line;
    
    /* Check for * prefix (degrees) */
    if (*p == '*') {
        tr->is_star = 1;
        p++;
    }
    
    /* Skip "TR" or "tr" */
    if (tolower((unsigned char)*p) == 't') p++;
    if (tolower((unsigned char)*p) == 'r') p++;
    
    /* Parse transform ID */
    tr->transform_id = (int)strtol(p, (char**)&p, 10);
    
    /* Skip whitespace */
    while (*p && isspace((unsigned char)*p)) p++;
    
    /* Rest is the definition */
    if (*p) {
        tr->definition = arena_strdup(&ctx->arena, p);
    }

    ctx->transforms[ctx->transform_count++] = tr;
    return 1;
}
