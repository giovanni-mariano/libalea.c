// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_xml/xml_dom.h
 * @brief Minimal XML parser for native ALEA model files
 *
 * Custom XML parser following project style (no third-party deps).
 * Handles elements, attributes, text content, and basic XML escapes.
 * Namespaces, CDATA, DTDs, and external entities are not accepted.
 */

#ifndef ALEA_XML_DOM_H
#define ALEA_XML_DOM_H

#include <stddef.h>
#include <stdbool.h>
#include "util/arena.h"


/* ============================================================================
 * XML ATTRIBUTE
 * ============================================================================ */

typedef struct {
    char* name;
    char* value;
} alea_xml_dom_attr_t;

/* ============================================================================
 * XML ELEMENT
 * ============================================================================ */

typedef struct alea_xml_dom_element {
    char* tag_name;
    alea_xml_dom_attr_t* attrs;
    size_t attr_count;
    size_t attr_capacity;
    struct alea_xml_dom_element** children;
    size_t child_count;
    size_t child_capacity;
    char* text_content;
    struct alea_xml_dom_element* parent;
} alea_xml_dom_element_t;

/* ============================================================================
 * XML DOCUMENT
 * ============================================================================ */

typedef struct {
    alea_xml_dom_element_t* root;
    arena_t arena;
    char error_msg[512];
    int error_line;
    int error_col;
} alea_xml_dom_doc_t;

/* ============================================================================
 * PARSER API
 * ============================================================================ */

/**
 * @brief Parse XML from file
 * @param filename Path to XML file
 * @return Parsed document or NULL on error
 */
alea_xml_dom_doc_t* alea_xml_dom_parse_file(const char* filename);

/**
 * @brief Parse XML from string
 * @param xml_string XML content
 * @param length String length (0 = null-terminated)
 * @return Parsed document or NULL on error
 */
alea_xml_dom_doc_t* alea_xml_dom_parse_string(const char* xml_string, size_t length);

/**
 * @brief Free document and all associated memory
 */
void alea_xml_dom_doc_free(alea_xml_dom_doc_t* doc);

/**
 * @brief Get error message from failed parse
 */
const char* alea_xml_dom_get_error(const alea_xml_dom_doc_t* doc);

/* ============================================================================
 * ELEMENT QUERY API
 * ============================================================================ */

/**
 * @brief Find first child element by tag name
 * @param elem Parent element
 * @param tag_name Tag to find
 * @return First matching child or NULL
 */
alea_xml_dom_element_t* alea_xml_dom_find_child(const alea_xml_dom_element_t* elem,
                                             const char* tag_name);

/**
 * @brief Find all children with given tag name
 * @param elem Parent element
 * @param tag_name Tag to find
 * @param out_children Output array (must be pre-allocated)
 * @param max_count Maximum children to return
 * @return Number of children found
 */
size_t alea_xml_dom_find_children(const alea_xml_dom_element_t* elem,
                                 const char* tag_name,
                                 alea_xml_dom_element_t** out_children,
                                 size_t max_count);

/**
 * @brief Count children with given tag name
 */
size_t alea_xml_dom_count_children(const alea_xml_dom_element_t* elem,
                                  const char* tag_name);

/**
 * @brief Get attribute value by name
 * @return Attribute value or NULL if not found
 */
const char* alea_xml_dom_get_attr(const alea_xml_dom_element_t* elem,
                                 const char* attr_name);

/**
 * @brief Get attribute as integer
 * @param elem Element
 * @param attr_name Attribute name
 * @param default_value Value if attribute not found
 * @return Integer value or default
 */
int alea_xml_dom_get_attr_int(const alea_xml_dom_element_t* elem,
                             const char* attr_name,
                             int default_value);

/**
 * @brief Get attribute as double
 */
double alea_xml_dom_get_attr_double(const alea_xml_dom_element_t* elem,
                                   const char* attr_name,
                                   double default_value);

/**
 * @brief Parse space-separated doubles from string
 * @param str Input string
 * @param out_values Output array
 * @param max_values Maximum values to parse
 * @return Number of values parsed
 */
size_t alea_xml_dom_parse_doubles(const char* str, double* out_values, size_t max_values);

/**
 * @brief Parse space-separated integers from string
 */
size_t alea_xml_dom_parse_ints(const char* str, int* out_values, size_t max_values);


#endif /* ALEA_XML_DOM_H */
