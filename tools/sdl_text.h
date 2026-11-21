// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_SDL_TEXT_H
#define ALEA_SDL_TEXT_H
#include <SDL.h>
#include <string.h>

/* Small built-in bitmap alphabet: no font installation or SDL_ttf required. */
static void alea_ui_text(SDL_Renderer* renderer, int x, int y, const char* text) {
    static const unsigned char letters[][7] = {
        {14,17,17,31,17,17,17},{30,17,17,30,17,17,30},
        {14,17,16,16,16,17,14},{30,17,17,17,17,17,30},
        {31,16,16,30,16,16,31},{31,16,16,30,16,16,16},
        {14,17,16,23,17,17,15},{17,17,17,31,17,17,17},
        {14,4,4,4,4,4,14},{7,2,2,2,18,18,12},
        {17,18,20,24,20,18,17},{16,16,16,16,16,16,31},
        {17,27,21,21,17,17,17},{17,25,25,21,19,19,17},
        {14,17,17,17,17,17,14},{30,17,17,30,16,16,16},
        {14,17,17,17,21,18,13},{30,17,17,30,20,18,17},
        {15,16,16,14,1,1,30},{31,4,4,4,4,4,4},
        {17,17,17,17,17,17,14},{17,17,17,17,17,10,4},
        {17,17,17,21,21,21,10},{17,17,10,4,10,17,17},
        {17,17,10,4,4,4,4},{31,1,2,4,8,16,31},
        {0,0,14,1,15,17,15},{16,16,30,17,17,17,30},
        {0,0,14,17,16,17,14},{1,1,15,17,17,17,15},
        {0,0,14,17,31,16,14},{6,8,8,28,8,8,8},
        {0,0,15,17,15,1,14},{16,16,30,17,17,17,17},
        {4,0,12,4,4,4,14},{2,0,6,2,2,18,12},
        {16,16,18,20,24,20,18},{12,4,4,4,4,4,14},
        {0,0,26,21,21,21,21},{0,0,30,17,17,17,17},
        {0,0,14,17,17,17,14},{0,0,30,17,30,16,16},
        {0,0,15,17,15,1,1},{0,0,22,25,16,16,16},
        {0,0,15,16,14,1,30},{8,8,28,8,8,9,6},
        {0,0,17,17,17,19,13},{0,0,17,17,17,10,4},
        {0,0,17,17,21,21,10},{0,0,17,10,4,10,17},
        {0,0,17,17,15,1,14},{0,0,31,2,4,8,31}
    };
    static const unsigned char digits[][7] = {
        {14,17,19,21,25,17,14},{4,12,4,4,4,4,14},
        {14,17,1,2,4,8,31},{30,1,1,14,1,1,30},
        {2,6,10,18,31,2,2},{31,16,16,30,1,1,30},
        {14,16,16,30,17,17,14},{31,1,2,4,8,8,8},
        {14,17,17,14,17,17,14},{14,17,17,15,1,1,14}
    };
    int start = x;
    for (; *text; ++text, x += 12) {
        unsigned char glyph[7] = {0};
        if (*text == '\n') { y += 24; x = start - 12; continue; }
        if (*text >= 'A' && *text <= 'Z') memcpy(glyph, letters[*text-'A'], 7);
        else if (*text >= 'a' && *text <= 'z') memcpy(glyph, letters[26+*text-'a'], 7);
        else if (*text >= '0' && *text <= '9') memcpy(glyph, digits[*text-'0'], 7);
        else if (*text == '[' || *text == ']') {
            glyph[0]=glyph[6]=14;
            for (int i=1;i<6;++i) glyph[i]=*text=='['?8:2;
        }
        else if (*text == '-') glyph[3] = 14;
        else if (*text == '+') { glyph[1]=4; glyph[2]=4; glyph[3]=31; glyph[4]=4; glyph[5]=4; }
        else if (*text == '.') glyph[6] = 4;
        else if (*text == ':') { glyph[2]=4; glyph[5]=4; }
        else if (*text == '/') { for (int i=0;i<7;++i) glyph[i]=(unsigned char)(1 << (i*4/6)); }
        else if (*text == '|') { for (int i=0;i<7;++i) glyph[i]=4; }
        for (int row = 0; row < 7; ++row)
            for (int col = 0; col < 5; ++col)
                if (glyph[row] & (16 >> col)) {
                    SDL_Rect pixel = {x+col*2, y+row*2, 2, 2};
                    SDL_RenderFillRect(renderer, &pixel);
                }
    }
}

#endif
