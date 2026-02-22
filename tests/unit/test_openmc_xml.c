// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>



#include "util/arena.h"
#include "openmc/openmc_xml.h"

int main(int argc, char** argv)
{
  const char* out_path = (argc > 1) ? argv[1] : "materials.xml";

  arena_t arena;
  arena_init(&arena);

  openmc_xml_t x;
  /* 80-column wrap, 2-space indent, pretty-print on */
  openmc_xml_init_ex(&x, &arena, 1024, 80, /*indent*/ 2, /*pretty*/ true);
  
  /* Document header */
  if (!openmc_xml_start_document(&x, "1.0", "utf-8")) goto fail;
  
  /* <materials> */
  if (!openmc_xml_start_element(&x, "materials")) goto fail;
  if (!openmc_xml_end_start_tag(&x, false)) goto fail;

    /* <material id="1" name="UO2"> */
    if (!openmc_xml_start_element(&x, "material")) goto fail;
    if (!openmc_xml_attribute_i(&x, "id", 1)) goto fail;
    if (!openmc_xml_attribute(&x, "name", "UO2")) goto fail;
    if (!openmc_xml_end_start_tag(&x, false)) goto fail;

      /* <density value="10.97" units="g/cm3" /> */
      if (!openmc_xml_start_element(&x, "density")) goto fail;
      if (!openmc_xml_attribute_f(&x, "value", 10.97, 6, false)) goto fail;
      if (!openmc_xml_attribute(&x, "units", "g/cm3")) goto fail;
      if (!openmc_xml_end_start_tag(&x, true)) goto fail;

      /* Long temperatures list as text content */
      if (!openmc_xml_start_element(&x, "temperature")) goto fail;
      if (!openmc_xml_attribute(&x, "units", "K")) goto fail;
      if (!openmc_xml_end_start_tag(&x, false)) goto fail;
      double temps[24];
      for (int i = 0; i < 24; ++i) temps[i] = 300.0 + 25.0 * i;
      if (!openmc_xml_numbers(&x, temps, 24, /*precision*/ 6, /*scientific*/ false, /*max_per_line*/ 8)) goto fail;

      if (!openmc_xml_end_element(&x, "temperature")) goto fail;

      /* Two nuclides */
      if (!openmc_xml_start_element(&x, "nuclide")) goto fail;
      if (!openmc_xml_attribute(&x, "name", "U235")) goto fail;
      if (!openmc_xml_attribute_f(&x, "ao", 0.04, 6, false)) goto fail;
      if (!openmc_xml_end_start_tag(&x, true)) goto fail;

      if (!openmc_xml_start_element(&x, "nuclide")) goto fail;
      if (!openmc_xml_attribute(&x, "name", "U238")) goto fail;
      if (!openmc_xml_attribute_f(&x, "ao", 0.96, 6, false)) goto fail;
      if (!openmc_xml_end_start_tag(&x, true)) goto fail;

    if (!openmc_xml_end_element(&x, "material")) goto fail;

  if (!openmc_xml_end_element(&x, "materials")) goto fail;
  /* Write to file */
  {
    FILE* fp = fopen(out_path, "wb");
    if (!fp) {
      fprintf(stderr, "Error: could not open output file '%s'\n", out_path);
      goto fail;
    }
    if (!openmc_xml_write(&x, fp)) {
      fprintf(stderr, "Error: openmc_xml_write failed\n");
      fclose(fp);
      goto fail;
    }
    if (fclose(fp) != 0) {
      fprintf(stderr, "Warning: failed to close '%s' cleanly\n", out_path);
    }
  }

  fprintf(stderr, "Wrote OpenMC XML to: %s\n", out_path);
  arena_free(&arena);
  return 0;

fail:
  fprintf(stderr, "Builder ended in error state or an operation failed.\n");
  arena_free(&arena);
  return 1;
}
