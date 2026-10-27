// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "xml_dom.h"
#include "alea_alloc.h"
#include "expat/expat_config.h"
#include "expat.h"

#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define XML_CHUNK_SIZE 65536

typedef struct {
    XML_Parser parser;
    alea_xml_dom_doc_t* doc;
    alea_xml_dom_element_t* current;
    alea_xml_dom_limits_t limits;
    size_t depth, nodes, allocated, input_bytes;
} dom_parser_t;

alea_xml_dom_limits_t alea_xml_dom_default_limits(void) {
    return (alea_xml_dom_limits_t){
        ALEA_XML_MAX_DEPTH, ALEA_XML_MAX_NODES, ALEA_XML_MAX_ATTRIBUTES,
        ALEA_XML_MAX_INPUT_BYTES, ALEA_XML_MAX_DOM_BYTES, ALEA_XML_MAX_TEXT_BYTES
    };
}

static void dom_error(dom_parser_t* p, const char* message) {
    if (p->doc->error_msg[0]) return;
    p->doc->error_line = (int)XML_GetCurrentLineNumber64(p->parser);
    p->doc->error_col = (int)XML_GetCurrentColumnNumber64(p->parser) + 1;
    snprintf(p->doc->error_msg, sizeof(p->doc->error_msg), "Line %d, col %d: %s",
             p->doc->error_line, p->doc->error_col, message);
    (void)XML_StopParser(p->parser, XML_FALSE);
}

static void* dom_alloc(dom_parser_t* p, size_t size) {
    if (size > p->limits.max_dom_bytes - p->allocated) {
        dom_error(p, "XML document memory limit exceeded");
        return NULL;
    }
    void* result = arena_alloc(&p->doc->arena, size);
    if (!result) {
        dom_error(p, "Out of memory building XML document");
        return NULL;
    }
    p->allocated += size;
    return result;
}

static char* dom_copy(dom_parser_t* p, const char* text) {
    size_t length = strlen(text);
    if (length == SIZE_MAX) { dom_error(p, "XML string is too long"); return NULL; }
    char* copy = dom_alloc(p, length + 1);
    if (copy) memcpy(copy, text, length + 1);
    return copy;
}

static void XMLCALL start_element(void* context, const XML_Char* name,
                                  const XML_Char** attributes) {
    dom_parser_t* p = context;
    if (p->doc->error_msg[0]) return;
    if (p->depth >= p->limits.max_depth) { dom_error(p, "XML nesting limit exceeded"); return; }
    if (p->nodes >= p->limits.max_nodes) { dom_error(p, "XML element limit exceeded"); return; }
    if (strchr(name, ':')) { dom_error(p, "XML namespaces are not supported"); return; }
    size_t count = 0;
    while (attributes[2 * count]) {
        const char* attr = attributes[2 * count];
        if (strchr(attr, ':') || strcmp(attr, "xmlns") == 0) {
            dom_error(p, "XML namespaces are not supported"); return;
        }
        if (count >= p->limits.max_attributes || count >= SIZE_MAX / 2 - 1) {
            dom_error(p, "XML attribute limit exceeded"); return;
        }
        count++;
    }
    if (count > SIZE_MAX / sizeof(alea_xml_dom_attr_t)) {
        dom_error(p, "XML attributes are too large"); return;
    }
    alea_xml_dom_element_t* elem = dom_alloc(p, sizeof(*elem));
    if (!elem) return;
    elem->tag_name = dom_copy(p, name);
    if (!elem->tag_name) return;
    if (count) {
        elem->attrs = dom_alloc(p, count * sizeof(*elem->attrs));
        if (!elem->attrs) return;
        for (size_t i = 0; i < count; i++) {
            elem->attrs[i].name = dom_copy(p, attributes[2 * i]);
            if (!elem->attrs[i].name) return;
            elem->attrs[i].value = dom_copy(p, attributes[2 * i + 1]);
            if (!elem->attrs[i].value) return;
        }
        elem->attr_count = elem->attr_capacity = count;
    }
    if (p->current) {
        alea_xml_dom_element_t* parent = p->current;
        if (parent->child_count == parent->child_capacity) {
            size_t capacity = parent->child_capacity ? parent->child_capacity * 2 : 4;
            if (capacity < parent->child_capacity || capacity > SIZE_MAX / sizeof(elem)) {
                dom_error(p, "XML child list is too large"); return;
            }
            alea_xml_dom_element_t** children = dom_alloc(p, capacity * sizeof(*children));
            if (!children) return;
            if (parent->child_count)
                memcpy(children, parent->children, parent->child_count * sizeof(*children));
            parent->children = children;
            parent->child_capacity = capacity;
        }
        parent->children[parent->child_count++] = elem;
        elem->parent = parent;
    } else {
        p->doc->root = elem;
    }
    p->current = elem;
    p->nodes++;
    p->depth++;
}

static void XMLCALL character_data(void* context, const XML_Char* text, int length) {
    dom_parser_t* p = context;
    alea_xml_dom_element_t* elem = p->current;
    if (!elem || p->doc->error_msg[0] || length <= 0) return;
    size_t size = (size_t)length;
    if (size > p->limits.max_text_bytes - elem->text_length ||
        size >= SIZE_MAX - elem->text_length) {
        dom_error(p, "XML element text limit exceeded"); return;
    }
    size_t needed = elem->text_length + size + 1;
    if (needed > elem->text_capacity) {
        size_t capacity = elem->text_capacity ? elem->text_capacity : 64;
        while (capacity < needed) {
            if (capacity > SIZE_MAX / 2) { capacity = needed; break; }
            capacity *= 2;
        }
        char* grown = dom_alloc(p, capacity);
        if (!grown) return;
        if (elem->text_length) memcpy(grown, elem->text_content, elem->text_length);
        elem->text_content = grown;
        elem->text_capacity = capacity;
    }
    memcpy(elem->text_content + elem->text_length, text, size);
    elem->text_length += size;
    elem->text_content[elem->text_length] = '\0';
}

static void XMLCALL end_element(void* context, const XML_Char* name) {
    (void)name;
    dom_parser_t* p = context;
    if (p->doc->error_msg[0] || !p->current) return;
    alea_xml_dom_element_t* elem = p->current;
    /* Match importer expectations for numeric lists, while preserving all
     * interior text across Expat callbacks, entity references and comments. */
    if (elem->text_content) {
        size_t first = 0, last = elem->text_length;
        while (first < last && isspace((unsigned char)elem->text_content[first])) first++;
        while (last > first && isspace((unsigned char)elem->text_content[last - 1])) last--;
        memmove(elem->text_content, elem->text_content + first, last - first);
        elem->text_length = last - first;
        elem->text_content[elem->text_length] = '\0';
    }
    p->current = elem->parent;
    p->depth--;
}

static void XMLCALL reject_doctype(void* context, const XML_Char* name,
                                   const XML_Char* system_id, const XML_Char* public_id,
                                   int internal_subset) {
    (void)name; (void)system_id; (void)public_id; (void)internal_subset;
    dom_error(context, "DTD and entity declarations are not allowed");
}

static void XMLCALL reject_entity(void* context, const XML_Char* name, int parameter) {
    (void)name; (void)parameter;
    dom_error(context, "Entity references are not allowed");
}

static int XMLCALL reject_external(XML_Parser parser, const XML_Char* context,
                                    const XML_Char* base, const XML_Char* system_id,
                                    const XML_Char* public_id) {
    (void)context; (void)base; (void)system_id; (void)public_id;
    dom_error(XML_GetUserData(parser), "External entities are not allowed");
    return XML_STATUS_ERROR;
}

static int parser_init(dom_parser_t* p, const alea_xml_dom_limits_t* limits) {
    memset(p, 0, sizeof(*p));
    p->limits = limits ? *limits : alea_xml_dom_default_limits();
    p->doc = alea_calloc(1, sizeof(*p->doc));
    if (!p->doc) return 0;
    if (!arena_init(&p->doc->arena)) {
        snprintf(p->doc->error_msg, sizeof(p->doc->error_msg), "Out of memory initializing XML document");
        return 0;
    }
    const XML_Memory_Handling_Suite memory = {alea_malloc, alea_realloc, free};
    p->parser = XML_ParserCreate_MM(NULL, &memory, NULL);
    if (!p->parser) {
        snprintf(p->doc->error_msg, sizeof(p->doc->error_msg), "Out of memory initializing XML parser");
        return 0;
    }
    XML_SetUserData(p->parser, p);
    XML_SetElementHandler(p->parser, start_element, end_element);
    XML_SetCharacterDataHandler(p->parser, character_data);
    XML_SetStartDoctypeDeclHandler(p->parser, reject_doctype);
    XML_SetSkippedEntityHandler(p->parser, reject_entity);
    XML_SetExternalEntityRefHandler(p->parser, reject_external);
    return 1;
}

static int parse_chunk(dom_parser_t* p, const char* bytes, size_t length, int final) {
    if (length > p->limits.max_input_bytes - p->input_bytes) {
        dom_error(p, "XML input size limit exceeded"); return 0;
    }
    p->input_bytes += length;
    if (XML_Parse(p->parser, bytes, (int)length, final) == XML_STATUS_ERROR) {
        if (!p->doc->error_msg[0]) dom_error(p, XML_ErrorString(XML_GetErrorCode(p->parser)));
        return 0;
    }
    return 1;
}

alea_xml_dom_doc_t* alea_xml_dom_parse_string_with_limits(
        const char* xml, size_t length, const alea_xml_dom_limits_t* limits) {
    if (!xml) return NULL;
    if (length == 0) length = strlen(xml);
    dom_parser_t p;
    if (!parser_init(&p, limits)) return p.doc;
    if (length > p.limits.max_input_bytes) {
        dom_error(&p, "XML input size limit exceeded");
    } else {
        size_t offset = 0;
        do {
            size_t count = length - offset;
            if (count > XML_CHUNK_SIZE) count = XML_CHUNK_SIZE;
            int final = count == length - offset;
            if (!parse_chunk(&p, xml + offset, count, final)) break;
            offset += count;
            if (final) break;
        } while (offset < length);
    }
    XML_ParserFree(p.parser);
    return p.doc;
}

alea_xml_dom_doc_t* alea_xml_dom_parse_string(const char* xml, size_t length) {
    return alea_xml_dom_parse_string_with_limits(xml, length, NULL);
}

alea_xml_dom_doc_t* alea_xml_dom_parse_file(const char* filename) {
    if (!filename) return NULL;
    dom_parser_t p;
    if (!parser_init(&p, NULL)) return p.doc;
    FILE* file = fopen(filename, "rb");
    if (!file) {
        dom_error(&p, "Failed to open XML file");
    } else {
        char buffer[XML_CHUNK_SIZE];
        for (;;) {
            size_t count = fread(buffer, 1, sizeof(buffer), file);
            if (ferror(file)) { dom_error(&p, "Failed to read XML file"); break; }
            int final = feof(file) != 0;
            if (!parse_chunk(&p, buffer, count, final) || final) break;
        }
        if (fclose(file) != 0 && !p.doc->error_msg[0]) dom_error(&p, "Failed to close XML file");
    }
    XML_ParserFree(p.parser);
    return p.doc;
}

void alea_xml_dom_doc_free(alea_xml_dom_doc_t* doc) {
    if (!doc) return;
    arena_free(&doc->arena);
    free(doc);
}

const char* alea_xml_dom_get_error(const alea_xml_dom_doc_t* doc) {
    if (!doc) return "NULL document";
    return doc->error_msg[0] ? doc->error_msg : NULL;
}

alea_xml_dom_element_t* alea_xml_dom_find_child(const alea_xml_dom_element_t* elem,
                                             const char* tag_name) {
    if (!elem || !tag_name) return NULL;

    for (size_t i = 0; i < elem->child_count; i++) {
        if (strcmp(elem->children[i]->tag_name, tag_name) == 0) {
            return elem->children[i];
        }
    }
    return NULL;
}

size_t alea_xml_dom_find_children(const alea_xml_dom_element_t* elem,
                                 const char* tag_name,
                                 alea_xml_dom_element_t** out_children,
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

size_t alea_xml_dom_count_children(const alea_xml_dom_element_t* elem,
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

const char* alea_xml_dom_get_attr(const alea_xml_dom_element_t* elem,
                                 const char* attr_name) {
    if (!elem || !attr_name) return NULL;

    for (size_t i = 0; i < elem->attr_count; i++) {
        if (strcmp(elem->attrs[i].name, attr_name) == 0) {
            return elem->attrs[i].value;
        }
    }
    return NULL;
}

int alea_xml_dom_get_attr_int(const alea_xml_dom_element_t* elem,
                             const char* attr_name,
                             int default_value) {
    const char* val = alea_xml_dom_get_attr(elem, attr_name);
    if (!val) return default_value;
    return (int)strtol(val, NULL, 10);
}

double alea_xml_dom_get_attr_double(const alea_xml_dom_element_t* elem,
                                   const char* attr_name,
                                   double default_value) {
    const char* val = alea_xml_dom_get_attr(elem, attr_name);
    if (!val) return default_value;
    return strtod(val, NULL);
}

size_t alea_xml_dom_parse_doubles(const char* str, double* out_values, size_t max_values) {
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

size_t alea_xml_dom_parse_ints(const char* str, int* out_values, size_t max_values) {
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
