// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_test.h"
#include "util/xml_dom.h"
#include "util/alea_alloc.h"
#include "openmc/openmc_parse.h"
#include "alea_xml.h"
#include "alea_openmc.h"
#include <stdlib.h>
#include <string.h>

static char* nested_xml(size_t depth) {
    char* text = malloc(depth * 7 + 1);
    if (!text) return NULL;
    for (size_t i = 0; i < depth; i++) memcpy(text + 3 * i, "<a>", 3);
    for (size_t i = 0; i < depth; i++) memcpy(text + 3 * depth + 4 * i, "</a>", 4);
    text[depth * 7] = '\0';
    return text;
}

TEST(xml_nesting_limits_apply_to_both_importers) {
    char* text = nested_xml(ALEA_XML_MAX_DEPTH);
    ASSERT_NOT_NULL(text);
    alea_xml_dom_doc_t* doc = alea_xml_dom_parse_string(text, 0);
    ASSERT_NOT_NULL(doc); ASSERT_NULL(alea_xml_dom_get_error(doc));
    alea_xml_dom_doc_free(doc); free(text);
    text = nested_xml(100000);
    ASSERT_NOT_NULL(text);
    doc = alea_xml_dom_parse_string(text, 700000);
    ASSERT_NOT_NULL(doc);
    ASSERT_NOT_NULL(strstr(alea_xml_dom_get_error(doc), "nesting"));
    alea_xml_dom_doc_free(doc);
    openmc_xml_doc_t* openmc = openmc_xml_parse_string(text, 700000);
    ASSERT_NOT_NULL(openmc);
    ASSERT_NOT_NULL(strstr(openmc_xml_get_error(openmc), "nesting"));
    openmc_xml_doc_free(openmc);
    ASSERT_NULL(alea_xml_load_string(text, 700000));
    ASSERT_NULL(openmc_load_string(text, 700000));
    free(text);
}

TEST(xml_rejects_malformed_documents_and_external_content) {
    static const char* const invalid[] = {
        "", "<a>", "<a></b>", "<a/><b/>", "<a/>garbage",
        "<a x='1' x='2'/>", "<a>&unknown;</a>", "<a>&#0;</a>",
        "<a>&#xD800;</a>", "<a>&#x110000;</a>", "<a>\xff</a>",
        "<!DOCTYPE a><a/>",
        "<!DOCTYPE a SYSTEM 'file:///nonexistent'><a/>",
        "<!DOCTYPE a [<!ENTITY e SYSTEM 'file:///nonexistent'>]><a>&e;</a>",
        "<!DOCTYPE a [<!ENTITY e 'value'>]><a>&e;</a>",
        "<a xmlns='urn:test'/>", "<a xmlns:x='urn:test'><x:b/></a>"
    };
    for (size_t i = 0; i < sizeof(invalid) / sizeof(*invalid); i++) {
        alea_xml_dom_doc_t* doc = alea_xml_dom_parse_string(invalid[i], 0);
        ASSERT_NOT_NULL(doc); ASSERT_NOT_NULL(alea_xml_dom_get_error(doc));
        alea_xml_dom_doc_free(doc);
        openmc_xml_doc_t* other = openmc_xml_parse_string(invalid[i], 0);
        ASSERT_NOT_NULL(other); ASSERT_NOT_NULL(openmc_xml_get_error(other));
        openmc_xml_doc_free(other);
    }
    const char embedded_null[] = "<a/>\0<a/>";
    alea_xml_dom_doc_t* doc = alea_xml_dom_parse_string(embedded_null, sizeof(embedded_null) - 1);
    ASSERT_NOT_NULL(doc); ASSERT_NOT_NULL(alea_xml_dom_get_error(doc));
    alea_xml_dom_doc_free(doc);
}

TEST(xml_preserves_unicode_entities_and_split_text) {
    const char* text = "<?xml version='1.0'?><root title='A &amp; &#x3b1; &#10;'> 12 &amp;<!--split--> 34<![CDATA[ <56>]]> &#x1f600; </root>";
    alea_xml_dom_doc_t* doc = alea_xml_dom_parse_string(text, 0);
    ASSERT_NOT_NULL(doc); ASSERT_NULL(alea_xml_dom_get_error(doc));
    ASSERT_STR_EQ(alea_xml_dom_get_attr(doc->root, "title"), "A & \xce\xb1 \n");
    ASSERT_STR_EQ(doc->root->text_content, "12 & 34 <56> \xf0\x9f\x98\x80");
    alea_xml_dom_doc_free(doc);

    size_t length = 200000;
    char* large = malloc(length + 8);
    ASSERT_NOT_NULL(large);
    memcpy(large, "<a>", 3); memset(large + 3, 'x', length);
    memcpy(large + 3 + length, "</a>", 5);
    doc = alea_xml_dom_parse_string(large, length + 7);
    ASSERT_NOT_NULL(doc); ASSERT_NULL(alea_xml_dom_get_error(doc));
    ASSERT_EQ(strlen(doc->root->text_content), length);
    ASSERT_EQ(doc->root->text_content[length - 1], 'x');
    alea_xml_dom_doc_free(doc); free(large);
}

TEST(xml_document_resource_limits) {
    alea_xml_dom_limits_t limits = alea_xml_dom_default_limits();
    limits.max_nodes = 2;
    alea_xml_dom_doc_t* doc = alea_xml_dom_parse_string_with_limits("<a><b/><c/></a>", 0, &limits);
    ASSERT_NOT_NULL(doc); ASSERT_NOT_NULL(alea_xml_dom_get_error(doc));
    alea_xml_dom_doc_free(doc);
    limits = alea_xml_dom_default_limits(); limits.max_attributes = 1;
    doc = alea_xml_dom_parse_string_with_limits("<a x='1' y='2'/>", 0, &limits);
    ASSERT_NOT_NULL(doc); ASSERT_NOT_NULL(alea_xml_dom_get_error(doc));
    alea_xml_dom_doc_free(doc);
    limits = alea_xml_dom_default_limits(); limits.max_text_bytes = 3;
    doc = alea_xml_dom_parse_string_with_limits("<a>abc</a>", 0, &limits);
    ASSERT_NOT_NULL(doc); ASSERT_NULL(alea_xml_dom_get_error(doc));
    alea_xml_dom_doc_free(doc);
    doc = alea_xml_dom_parse_string_with_limits("<a>abcd</a>", 0, &limits);
    ASSERT_NOT_NULL(doc); ASSERT_NOT_NULL(alea_xml_dom_get_error(doc));
    alea_xml_dom_doc_free(doc);
    limits = alea_xml_dom_default_limits(); limits.max_input_bytes = 3;
    doc = alea_xml_dom_parse_string_with_limits("<a/>", 0, &limits);
    ASSERT_NOT_NULL(doc); ASSERT_NOT_NULL(alea_xml_dom_get_error(doc));
    alea_xml_dom_doc_free(doc);
    limits = alea_xml_dom_default_limits(); limits.max_dom_bytes = 1;
    doc = alea_xml_dom_parse_string_with_limits("<a/>", 0, &limits);
    ASSERT_NOT_NULL(doc); ASSERT_NOT_NULL(alea_xml_dom_get_error(doc));
    alea_xml_dom_doc_free(doc);
}

TEST(xml_allocation_failures_clean_up_and_allow_retry) {
    const char* text = "<root a='value'><child>text &amp; data</child><child/></root>";
    int succeeded = 0;
    for (long failure = 0; failure < 200; failure++) {
        alea_alloc_fail_after(failure);
        alea_xml_dom_doc_t* doc = alea_xml_dom_parse_string(text, 0);
        alea_alloc_fail_after(-1);
        succeeded = doc && !alea_xml_dom_get_error(doc);
        alea_xml_dom_doc_free(doc);
        if (succeeded) break;
    }
    ASSERT(succeeded);
    alea_xml_dom_doc_t* retry = alea_xml_dom_parse_string(text, 0);
    ASSERT_NOT_NULL(retry); ASSERT_NULL(alea_xml_dom_get_error(retry));
    alea_xml_dom_doc_free(retry);
}

TEST(xml_file_parser_consumes_the_complete_document) {
    const char* path = "test_xml_dom.tmp";
    FILE* file = fopen(path, "wb"); ASSERT_NOT_NULL(file);
    fputs("<a/>\n<!-- trailing comment -->", file); ASSERT_EQ(fclose(file), 0);
    alea_xml_dom_doc_t* doc = alea_xml_dom_parse_file(path);
    ASSERT_NOT_NULL(doc); ASSERT_NULL(alea_xml_dom_get_error(doc));
    alea_xml_dom_doc_free(doc);
    file = fopen(path, "ab"); ASSERT_NOT_NULL(file);
    fputs("<b/>", file); ASSERT_EQ(fclose(file), 0);
    doc = alea_xml_dom_parse_file(path);
    ASSERT_NOT_NULL(doc); ASSERT_NOT_NULL(alea_xml_dom_get_error(doc));
    alea_xml_dom_doc_free(doc); remove(path);
    doc = alea_xml_dom_parse_file(path);
    ASSERT_NOT_NULL(doc); ASSERT_NOT_NULL(alea_xml_dom_get_error(doc));
    alea_xml_dom_doc_free(doc);
}

TEST(xml_importers_reject_errors_after_a_valid_root) {
    const char* openmc = "<geometry><surface id='1' type='sphere' coeffs='0 0 0 1'/>"
                         "<cell id='1' material='void' region='-1'/></geometry>";
    const char* native = "<alea version='1' length_units='cm'><surfaces>"
                         "<surface id='1' type='sphere' coeffs='0 0 0 1'/></surfaces>"
                         "<cells><cell id='1' material='0' region='-1'/></cells></alea>";
    openmc_model_t* model = openmc_load_string(openmc, 0);
    ASSERT_NOT_NULL(model); openmc_model_destroy(model);
    alea_model_t* other = alea_xml_load_string(native, 0);
    ASSERT_NOT_NULL(other); alea_model_destroy(other);
    char malformed[1024];
    snprintf(malformed, sizeof(malformed), "%s<extra/>", openmc);
    ASSERT_NULL(openmc_load_string(malformed, 0));
    ASSERT_EQ(alea_error_code(), ALEA_ERR_PARSE_ERROR);
    snprintf(malformed, sizeof(malformed), "%s<extra/>", native);
    ASSERT_NULL(alea_xml_load_string(malformed, 0));
    const char* path = "test_xml_invalid_tail.tmp";
    FILE* file = fopen(path, "wb"); ASSERT_NOT_NULL(file);
    fprintf(file, "%s<extra/>", openmc); ASSERT_EQ(fclose(file), 0);
    ASSERT_NULL(openmc_load(path));
    remove(path);
}

TEST_MAIN()
