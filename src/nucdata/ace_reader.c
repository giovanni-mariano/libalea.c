// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file ace_reader.c
 * @brief ACE file reader — Type 1 (ASCII) and Type 2 (binary)
 *
 * ACE Type 1 (ASCII) layout:
 *   Line 1:  ZAID  AWR  temperature  date  (+ comment to col 80)
 *   Line 2:  comment (cols 1-70), MAT (cols 71-80)
 *   Lines 3-12:  IZ(i), AW(i) pairs — 4 pairs per line, 10 lines
 *   Lines 13-14: NXS(1..16), 8 per line
 *   Lines 15-16: JXS(1..32), 8 per line
 *   Remaining:   XSS array, 4 values per line
 *
 * The 'address' from xsdir is the 1-based line number where the table starts.
 */

#include "alea_nucdata.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define ACE_LINE_MAX 256

/** Skip to a specific 1-based line number in a file */
static alea_error_t skip_to_line(FILE* fp, int target_line) {
    char buf[ACE_LINE_MAX];
    for (int i = 1; i < target_line; i++) {
        if (!fgets(buf, sizeof(buf), fp))
            return ALEA_ERR_FILE_READ;
    }
    return ALEA_OK;
}

/**
 * Read ACE Type 1 (ASCII) table.
 *
 * The file pointer should be positioned at the start of the table
 * (line 1 of the ACE header).
 */
static alea_error_t ace_read_type1(FILE* fp, alea_nuc_ace_table_t* table) {
    char line[ACE_LINE_MAX];

    /* --- Line 1: ZAID, AWR, temperature, date, comment --- */
    if (!fgets(line, sizeof(line), fp)) return ALEA_ERR_FILE_READ;

    /* ZAID is first token, AWR second, temperature third, date fourth */
    char date_buf[12] = {0};
    int n = sscanf(line, "%23s %lf %lf %11s",
                   table->zaid, &table->awr, &table->temperature, date_buf);
    if (n < 3) return ALEA_ERR_PARSE_ERROR;
    memcpy(table->date, date_buf, sizeof(table->date) - 1);
    table->date[sizeof(table->date) - 1] = '\0';

    /* Rest of line 1 after the 4 fields is the comment */
    /* Find position after 4th token for comment extraction */
    {
        const char* p = line;
        int tokens = 0;
        while (tokens < 4 && *p) {
            while (*p && isspace((unsigned char)*p)) p++;
            if (!*p) break;
            while (*p && !isspace((unsigned char)*p)) p++;
            tokens++;
        }
        while (*p && isspace((unsigned char)*p)) p++;
        if (*p) {
            size_t clen = strlen(p);
            if (clen > 0 && p[clen-1] == '\n') clen--;
            if (clen >= sizeof(table->comment)) clen = sizeof(table->comment) - 1;
            memcpy(table->comment, p, clen);
            table->comment[clen] = '\0';
        }
    }

    /* --- Line 2: comment/source line --- */
    if (!fgets(line, sizeof(line), fp)) return ALEA_ERR_FILE_READ;
    /* Usually contains processing info, we skip it */

    /*
     * Detect old-style vs new-style ACE format:
     *   Old-style (NJOY99):  10 IZ/AW lines (2 pairs/line) + 2 NXS + 2 JXS = 14 lines
     *   New-style (NJOY2016): 4 IZ/AW lines (4 pairs/line) + 2 NXS + 4 JXS = 10 lines
     *
     * Detection: read line 3, count how many integer-float pairs it has.
     * Old-style has 2 pairs per line; new-style has 4 pairs per line.
     */

    /* Read first IZ/AW line to detect format */
    if (!fgets(line, sizeof(line), fp)) return ALEA_ERR_FILE_READ;

    int pairs_on_first_line = 0;
    {
        const char* p = line;
        while (1) {
            char* end;
            strtol(p, &end, 10);
            if (end == p) break;
            p = end;
            strtod(p, &end);
            if (end == p) break;
            p = end;
            pairs_on_first_line++;
        }
    }

    /* New-style: 4 pairs/line → 4 IZ/AW lines total */
    /* Old-style: 2 pairs/line → 10 IZ/AW lines total (or could be 4 pairs too) */
    int iz_aw_lines;
    if (pairs_on_first_line >= 4)
        iz_aw_lines = 4;    /* new-style */
    else
        iz_aw_lines = 10;   /* old-style */

    /* Parse first IZ/AW line (already read) */
    {
        const char* p = line;
        for (int col = 0; col < pairs_on_first_line && col < 4; col++) {
            char* end;
            long iz = strtol(p, &end, 10);
            if (end == p) break;
            p = end;
            double aw = strtod(p, &end);
            if (end == p) break;
            p = end;
            table->iz[col] = (int)iz;
            table->aw[col] = aw;
        }
    }

    /* Read remaining IZ/AW lines */
    for (int row = 1; row < iz_aw_lines; row++) {
        if (!fgets(line, sizeof(line), fp)) return ALEA_ERR_FILE_READ;
        const char* p = line;
        for (int col = 0; col < pairs_on_first_line; col++) {
            int idx = row * pairs_on_first_line + col;
            if (idx < 16) {
                char* end;
                long iz = strtol(p, &end, 10);
                if (end == p) break;
                p = end;
                double aw = strtod(p, &end);
                if (end == p) break;
                p = end;
                table->iz[idx] = (int)iz;
                table->aw[idx] = aw;
            }
        }
    }

    /* --- NXS(1..16), 8 per line --- */
    for (int row = 0; row < 2; row++) {
        if (!fgets(line, sizeof(line), fp)) return ALEA_ERR_FILE_READ;
        const char* p = line;
        for (int col = 0; col < 8; col++) {
            char* end;
            long val = strtol(p, &end, 10);
            if (end == p) return ALEA_ERR_PARSE_ERROR;
            p = end;
            table->nxs[row * 8 + col] = (int)val;
        }
    }

    /* --- JXS(1..32), 8 per line --- */
    for (int row = 0; row < 4; row++) {
        if (!fgets(line, sizeof(line), fp)) return ALEA_ERR_FILE_READ;
        const char* p = line;
        for (int col = 0; col < 8; col++) {
            char* end;
            long val = strtol(p, &end, 10);
            if (end == p) return ALEA_ERR_PARSE_ERROR;
            p = end;
            table->jxs[row * 8 + col] = (int)val;
        }
    }

    /* --- XSS data array --- */
    table->xss_length = table->nxs[0]; /* NXS[1] in Fortran 1-based = nxs[0] in C */
    if (table->xss_length <= 0 || table->xss_length > 500000000) {
        /* If table length is invalid, the format detection may have been wrong.
         * Try the other format: if we guessed old-style (10 lines), it might be
         * new-style (4 lines) with only 2 pairs on the first line, or vice versa. */
        return ALEA_ERR_PARSE_ERROR;
    }

    table->xss = malloc((size_t)table->xss_length * sizeof(double));
    if (!table->xss) return ALEA_ERR_OUT_OF_MEMORY;

    int read = 0;
    while (read < table->xss_length) {
        if (!fgets(line, sizeof(line), fp)) {
            free(table->xss);
            table->xss = NULL;
            return ALEA_ERR_FILE_READ;
        }
        const char* p = line;
        while (read < table->xss_length) {
            char* end;
            double val = strtod(p, &end);
            if (end == p) break;
            table->xss[read++] = val;
            p = end;
        }
    }

    return ALEA_OK;
}

/**
 * Read ACE Type 2 (binary) table.
 */
static alea_error_t ace_read_type2(FILE* fp, int address, alea_nuc_ace_table_t* table) {
    /* Seek to byte offset */
    if (fseek(fp, (long)address, SEEK_SET) != 0)
        return ALEA_ERR_FILE_READ;

    /*
     * Binary ACE format (platform-native endianness):
     *   Record 1: ZAID(10A1), AWR(1F), TZ(1F), HD(10A1), HK(70A1), HM(10A1)
     *             then IZ(16I), AW(16F), NXS(16I), JXS(32I)
     *   Record 2: XSS(NXS[1] doubles)
     *
     * Fortran unformatted records have 4-byte length markers before and after.
     */

    /* Skip Fortran record marker */
    int32_t rec_len;
    if (fread(&rec_len, 4, 1, fp) != 1) return ALEA_ERR_FILE_READ;

    /* ZAID: 10 characters */
    char zaid_buf[11] = {0};
    if (fread(zaid_buf, 1, 10, fp) != 10) return ALEA_ERR_FILE_READ;
    /* Trim trailing spaces */
    for (int i = 9; i >= 0 && zaid_buf[i] == ' '; i--) zaid_buf[i] = '\0';
    strncpy(table->zaid, zaid_buf, sizeof(table->zaid) - 1);

    /* AWR, temperature */
    if (fread(&table->awr, sizeof(double), 1, fp) != 1) return ALEA_ERR_FILE_READ;
    if (fread(&table->temperature, sizeof(double), 1, fp) != 1) return ALEA_ERR_FILE_READ;

    /* HD: date (10 chars) */
    char hd[11] = {0};
    if (fread(hd, 1, 10, fp) != 10) return ALEA_ERR_FILE_READ;
    strncpy(table->date, hd, sizeof(table->date) - 1);

    /* HK: comment (70 chars) */
    char hk[71] = {0};
    if (fread(hk, 1, 70, fp) != 70) return ALEA_ERR_FILE_READ;
    strncpy(table->comment, hk, sizeof(table->comment) - 1);

    /* HM: 10 chars (source identifier, skip) */
    char hm[11];
    if (fread(hm, 1, 10, fp) != 10) return ALEA_ERR_FILE_READ;

    /* IZ(16), AW(16) */
    int32_t iz32[16];
    if (fread(iz32, sizeof(int32_t), 16, fp) != 16) return ALEA_ERR_FILE_READ;
    for (int i = 0; i < 16; i++) table->iz[i] = iz32[i];
    if (fread(table->aw, sizeof(double), 16, fp) != 16) return ALEA_ERR_FILE_READ;

    /* NXS(16) */
    int32_t nxs32[16];
    if (fread(nxs32, sizeof(int32_t), 16, fp) != 16) return ALEA_ERR_FILE_READ;
    for (int i = 0; i < 16; i++) table->nxs[i] = nxs32[i];

    /* JXS(32) */
    int32_t jxs32[32];
    if (fread(jxs32, sizeof(int32_t), 32, fp) != 32) return ALEA_ERR_FILE_READ;
    for (int i = 0; i < 32; i++) table->jxs[i] = jxs32[i];

    /* End-of-record marker */
    if (fread(&rec_len, 4, 1, fp) != 1) return ALEA_ERR_FILE_READ;

    /* Record 2: XSS array */
    if (fread(&rec_len, 4, 1, fp) != 1) return ALEA_ERR_FILE_READ;

    table->xss_length = table->nxs[0];
    if (table->xss_length <= 0) return ALEA_ERR_PARSE_ERROR;

    table->xss = malloc((size_t)table->xss_length * sizeof(double));
    if (!table->xss) return ALEA_ERR_OUT_OF_MEMORY;

    size_t nread = fread(table->xss, sizeof(double), (size_t)table->xss_length, fp);
    if ((int)nread != table->xss_length) {
        free(table->xss);
        table->xss = NULL;
        return ALEA_ERR_FILE_READ;
    }

    return ALEA_OK;
}

alea_error_t alea_nuc_ace_read(const char* path, int address, int file_type,
                          alea_nuc_ace_table_t* table) {
    if (!path || !table) return ALEA_ERR_NULL_ARG;

    memset(table, 0, sizeof(*table));

    FILE* fp = fopen(path, file_type == 2 ? "rb" : "r");
    if (!fp) return ALEA_ERR_FILE_NOT_FOUND;

    alea_error_t err;

    if (file_type == 2) {
        err = ace_read_type2(fp, address, table);
    } else {
        /* Type 1 ASCII: address is 1-based line number */
        if (address > 1) {
            err = skip_to_line(fp, address);
            if (err != ALEA_OK) { fclose(fp); return err; }
        }
        err = ace_read_type1(fp, table);
    }

    fclose(fp);

    if (err == ALEA_OK) {
        /* Determine table type from ZAID suffix */
        alea_nuc_parse_zaid(table->zaid, NULL, NULL, NULL, &table->type);
    }

    return err;
}

void alea_nuc_ace_free(alea_nuc_ace_table_t* table) {
    if (!table) return;
    free(table->xss);
    table->xss = NULL;
    table->xss_length = 0;
}
