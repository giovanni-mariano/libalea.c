// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "openmc_parse.h"

openmc_xml_doc_t* openmc_xml_parse_string(const char* xml, size_t length) {
    return alea_xml_dom_parse_string(xml, length);
}
openmc_xml_doc_t* openmc_xml_parse_file(const char* filename) {
    return alea_xml_dom_parse_file(filename);
}
void openmc_xml_doc_free(openmc_xml_doc_t* doc) { alea_xml_dom_doc_free(doc); }
const char* openmc_xml_get_error(const openmc_xml_doc_t* doc) {
    return alea_xml_dom_get_error(doc);
}
openmc_xml_element_t* openmc_xml_find_child(const openmc_xml_element_t* elem, const char* tag) {
    return alea_xml_dom_find_child(elem, tag);
}
size_t openmc_xml_find_children(const openmc_xml_element_t* elem, const char* tag,
                                openmc_xml_element_t** out, size_t count) {
    return alea_xml_dom_find_children(elem, tag, out, count);
}
size_t openmc_xml_count_children(const openmc_xml_element_t* elem, const char* tag) {
    return alea_xml_dom_count_children(elem, tag);
}
const char* openmc_xml_get_attr(const openmc_xml_element_t* elem, const char* attr) {
    return alea_xml_dom_get_attr(elem, attr);
}
int openmc_xml_get_attr_int(const openmc_xml_element_t* elem, const char* attr, int fallback) {
    return alea_xml_dom_get_attr_int(elem, attr, fallback);
}
double openmc_xml_get_attr_double(const openmc_xml_element_t* elem, const char* attr, double fallback) {
    return alea_xml_dom_get_attr_double(elem, attr, fallback);
}
size_t openmc_parse_doubles(const char* text, double* values, size_t count) {
    return alea_xml_dom_parse_doubles(text, values, count);
}
size_t openmc_parse_ints(const char* text, int* values, size_t count) {
    return alea_xml_dom_parse_ints(text, values, count);
}
