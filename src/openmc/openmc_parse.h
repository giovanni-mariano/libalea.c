// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file openmc_parse.h
 * @brief Minimal XML parser for OpenMC geometry files
 *
 * OpenMC compatibility interface for the shared Expat-backed XML tree.
 * DTDs, external entities, and namespaces are rejected. CDATA is text.
 */

#ifndef OPENMC_PARSE_H
#define OPENMC_PARSE_H

#include <stddef.h>
#include <stdbool.h>
#include "util/arena.h"


#include "util/xml_dom.h"
typedef alea_xml_dom_attr_t openmc_xml_attr_t;
typedef alea_xml_dom_element_t openmc_xml_element_t;
typedef alea_xml_dom_doc_t openmc_xml_doc_t;

/* ============================================================================
 * PARSER API
 * ============================================================================ */

/**
 * @brief Parse XML from file
 * @param filename Path to XML file
 * @return Document (inspect get_error), or NULL on allocation/argument failure
 */
openmc_xml_doc_t* openmc_xml_parse_file(const char* filename);

/**
 * @brief Parse XML from string
 * @param xml_string XML content
 * @param length String length (0 = null-terminated)
 * @return Document (inspect get_error), or NULL on allocation/argument failure
 */
openmc_xml_doc_t* openmc_xml_parse_string(const char* xml_string, size_t length);

/**
 * @brief Free document and all associated memory
 */
void openmc_xml_doc_free(openmc_xml_doc_t* doc);

/**
 * @brief Get error message from failed parse
 */
const char* openmc_xml_get_error(const openmc_xml_doc_t* doc);

/* ============================================================================
 * ELEMENT QUERY API
 * ============================================================================ */

/**
 * @brief Find first child element by tag name
 * @param elem Parent element
 * @param tag_name Tag to find
 * @return First matching child or NULL
 */
openmc_xml_element_t* openmc_xml_find_child(const openmc_xml_element_t* elem,
                                             const char* tag_name);

/**
 * @brief Find all children with given tag name
 * @param elem Parent element
 * @param tag_name Tag to find
 * @param out_children Output array (must be pre-allocated)
 * @param max_count Maximum children to return
 * @return Number of children found
 */
size_t openmc_xml_find_children(const openmc_xml_element_t* elem,
                                 const char* tag_name,
                                 openmc_xml_element_t** out_children,
                                 size_t max_count);

/**
 * @brief Count children with given tag name
 */
size_t openmc_xml_count_children(const openmc_xml_element_t* elem,
                                  const char* tag_name);

/**
 * @brief Get attribute value by name
 * @return Attribute value or NULL if not found
 */
const char* openmc_xml_get_attr(const openmc_xml_element_t* elem,
                                 const char* attr_name);

/**
 * @brief Get attribute as integer
 * @param elem Element
 * @param attr_name Attribute name
 * @param default_value Value if attribute not found
 * @return Integer value or default
 */
int openmc_xml_get_attr_int(const openmc_xml_element_t* elem,
                             const char* attr_name,
                             int default_value);

/**
 * @brief Get attribute as double
 */
double openmc_xml_get_attr_double(const openmc_xml_element_t* elem,
                                   const char* attr_name,
                                   double default_value);

/**
 * @brief Parse space-separated doubles from string
 * @param str Input string
 * @param out_values Output array
 * @param max_values Maximum values to parse
 * @return Number of values parsed
 */
size_t openmc_parse_doubles(const char* str, double* out_values, size_t max_values);

/**
 * @brief Parse space-separated integers from string
 */
size_t openmc_parse_ints(const char* str, int* out_values, size_t max_values);


#endif /* OPENMC_PARSE_H */
