// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_DIGEST_H
#define ALEA_DIGEST_H
#include <stddef.h>
#include <stdint.h>
typedef struct {uint32_t h[8];uint64_t bytes;unsigned char block[64];size_t used;} alea_sha256_t;
void alea_sha256_init(alea_sha256_t* s);
void alea_sha256_update(alea_sha256_t* s,const void* bytes,size_t n);
void alea_sha256_final(alea_sha256_t* s,unsigned char digest[32]);
void alea_digest_u64(alea_sha256_t* s,uint64_t v);
int alea_digest_double(alea_sha256_t* s,double v);
void alea_digest_string(alea_sha256_t* s,const char* v);
#endif
