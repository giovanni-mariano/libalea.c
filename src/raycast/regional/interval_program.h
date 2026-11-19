// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Short-circuit decisions for bounded-size regional expressions. Each jump
 * either finishes the program or steps back to a preceding instruction. */
#define QUERY_PROGRAM UINT32_C(0x60000000)
#define QUERY_PLANE UINT32_C(0x08000000)
#define QUERY_QUADRIC UINT32_C(0x10000000)
#define QUERY_OPERAND UINT32_C(0x07ffffff)
typedef struct {
    uint32_t operand;
    uint16_t yes, no;
} query_decision_t;
ALEA_VEC_DEFINE(query_decisions, uint32_t);
ALEA_VEC_DEFINE(wide_decisions, query_decision_t);

enum { QUERY_DECISION_LIMIT=768, QUERY_ACCEPT=UINT16_MAX };
