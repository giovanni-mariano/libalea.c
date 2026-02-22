// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_elements_data.c
 * @brief Complete natural isotopic composition data for all elements
 *
 * Data source: NIST Atomic Weights and Isotopic Compositions
 * https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
 *
 * This file contains:
 * - Natural isotope abundances for all 118 elements
 * - Standard atomic weights
 * - Element symbols and names
 *
 * For radioactive elements without stable isotopes, the most significant
 * isotopes are listed with zero abundance.
 */

#include "alea_materials.h"

/* ============================================================================
 * ISOTOPE DATA BY ELEMENT
 * ============================================================================
 * Format: {mass_number, abundance, atomic_mass}
 * Abundance is 0.0 for radioactive isotopes with no natural occurrence
 */

/* Z=1 Hydrogen */
static const alea_isotope_t isotopes_H[] = {
    {1, 0.999885, 1.00782503223},
    {2, 0.000115, 2.01410177812}
};

/* Z=2 Helium */
static const alea_isotope_t isotopes_He[] = {
    {3, 0.00000134, 3.0160293201},
    {4, 0.99999866, 4.00260325413}
};

/* Z=3 Lithium */
static const alea_isotope_t isotopes_Li[] = {
    {6, 0.0759, 6.0151228874},
    {7, 0.9241, 7.0160034366}
};

/* Z=4 Beryllium */
static const alea_isotope_t isotopes_Be[] = {
    {9, 1.0, 9.012183065}
};

/* Z=5 Boron */
static const alea_isotope_t isotopes_B[] = {
    {10, 0.199, 10.01293695},
    {11, 0.801, 11.00930536}
};

/* Z=6 Carbon */
static const alea_isotope_t isotopes_C[] = {
    {12, 0.9893, 12.0000000},
    {13, 0.0107, 13.00335483507}
};

/* Z=7 Nitrogen */
static const alea_isotope_t isotopes_N[] = {
    {14, 0.99636, 14.00307400443},
    {15, 0.00364, 15.00010889888}
};

/* Z=8 Oxygen */
static const alea_isotope_t isotopes_O[] = {
    {16, 0.99757, 15.99491461957},
    {17, 0.00038, 16.99913175650},
    {18, 0.00205, 17.99915961286}
};

/* Z=9 Fluorine */
static const alea_isotope_t isotopes_F[] = {
    {19, 1.0, 18.99840316273}
};

/* Z=10 Neon */
static const alea_isotope_t isotopes_Ne[] = {
    {20, 0.9048, 19.9924401762},
    {21, 0.0027, 20.993846685},
    {22, 0.0925, 21.991385114}
};

/* Z=11 Sodium */
static const alea_isotope_t isotopes_Na[] = {
    {23, 1.0, 22.9897692820}
};

/* Z=12 Magnesium */
static const alea_isotope_t isotopes_Mg[] = {
    {24, 0.7899, 23.985041697},
    {25, 0.1000, 24.985836976},
    {26, 0.1101, 25.982592968}
};

/* Z=13 Aluminum */
static const alea_isotope_t isotopes_Al[] = {
    {27, 1.0, 26.98153853}
};

/* Z=14 Silicon */
static const alea_isotope_t isotopes_Si[] = {
    {28, 0.92223, 27.97692653465},
    {29, 0.04685, 28.97649466490},
    {30, 0.03092, 29.973770136}
};

/* Z=15 Phosphorus */
static const alea_isotope_t isotopes_P[] = {
    {31, 1.0, 30.97376199842}
};

/* Z=16 Sulfur */
static const alea_isotope_t isotopes_S[] = {
    {32, 0.9499, 31.9720711744},
    {33, 0.0075, 32.9714589098},
    {34, 0.0425, 33.967867004},
    {36, 0.0001, 35.96708071}
};

/* Z=17 Chlorine */
static const alea_isotope_t isotopes_Cl[] = {
    {35, 0.7576, 34.968852682},
    {37, 0.2424, 36.965902602}
};

/* Z=18 Argon */
static const alea_isotope_t isotopes_Ar[] = {
    {36, 0.003336, 35.967545105},
    {38, 0.000629, 37.96273211},
    {40, 0.996035, 39.9623831237}
};

/* Z=19 Potassium */
static const alea_isotope_t isotopes_K[] = {
    {39, 0.932581, 38.9637064864},
    {40, 0.000117, 39.963998166},
    {41, 0.067302, 40.9618252579}
};

/* Z=20 Calcium */
static const alea_isotope_t isotopes_Ca[] = {
    {40, 0.96941, 39.962590863},
    {42, 0.00647, 41.95861783},
    {43, 0.00135, 42.95876644},
    {44, 0.02086, 43.95548156},
    {46, 0.00004, 45.9536890},
    {48, 0.00187, 47.95252276}
};

/* Z=21 Scandium */
static const alea_isotope_t isotopes_Sc[] = {
    {45, 1.0, 44.95590828}
};

/* Z=22 Titanium */
static const alea_isotope_t isotopes_Ti[] = {
    {46, 0.0825, 45.95262772},
    {47, 0.0744, 46.95175879},
    {48, 0.7372, 47.94794198},
    {49, 0.0541, 48.94786568},
    {50, 0.0518, 49.94478689}
};

/* Z=23 Vanadium */
static const alea_isotope_t isotopes_V[] = {
    {50, 0.00250, 49.94715601},
    {51, 0.99750, 50.94395704}
};

/* Z=24 Chromium */
static const alea_isotope_t isotopes_Cr[] = {
    {50, 0.04345, 49.94604183},
    {52, 0.83789, 51.94050623},
    {53, 0.09501, 52.94064815},
    {54, 0.02365, 53.93887916}
};

/* Z=25 Manganese */
static const alea_isotope_t isotopes_Mn[] = {
    {55, 1.0, 54.93804391}
};

/* Z=26 Iron */
static const alea_isotope_t isotopes_Fe[] = {
    {54, 0.05845, 53.93960899},
    {56, 0.91754, 55.93493633},
    {57, 0.02119, 56.93539284},
    {58, 0.00282, 57.93327443}
};

/* Z=27 Cobalt */
static const alea_isotope_t isotopes_Co[] = {
    {59, 1.0, 58.93319429}
};

/* Z=28 Nickel */
static const alea_isotope_t isotopes_Ni[] = {
    {58, 0.68077, 57.93534241},
    {60, 0.26223, 59.93078588},
    {61, 0.011399, 60.93105557},
    {62, 0.036346, 61.92834537},
    {64, 0.009255, 63.92796682}
};

/* Z=29 Copper */
static const alea_isotope_t isotopes_Cu[] = {
    {63, 0.6915, 62.92959772},
    {65, 0.3085, 64.92778970}
};

/* Z=30 Zinc */
static const alea_isotope_t isotopes_Zn[] = {
    {64, 0.4917, 63.92914201},
    {66, 0.2773, 65.92603381},
    {67, 0.0404, 66.92712775},
    {68, 0.1845, 67.92484455},
    {70, 0.0061, 69.9253192}
};

/* Z=31 Gallium */
static const alea_isotope_t isotopes_Ga[] = {
    {69, 0.60108, 68.9255735},
    {71, 0.39892, 70.92470258}
};

/* Z=32 Germanium */
static const alea_isotope_t isotopes_Ge[] = {
    {70, 0.2057, 69.92424875},
    {72, 0.2745, 71.922075826},
    {73, 0.0775, 72.923458956},
    {74, 0.3650, 73.921177761},
    {76, 0.0773, 75.921402726}
};

/* Z=33 Arsenic */
static const alea_isotope_t isotopes_As[] = {
    {75, 1.0, 74.92159457}
};

/* Z=34 Selenium */
static const alea_isotope_t isotopes_Se[] = {
    {74, 0.0089, 73.922475934},
    {76, 0.0937, 75.919213704},
    {77, 0.0763, 76.919914154},
    {78, 0.2377, 77.91730928},
    {80, 0.4961, 79.9165218},
    {82, 0.0873, 81.9166995}
};

/* Z=35 Bromine */
static const alea_isotope_t isotopes_Br[] = {
    {79, 0.5069, 78.9183376},
    {81, 0.4931, 80.9162897}
};

/* Z=36 Krypton */
static const alea_isotope_t isotopes_Kr[] = {
    {78, 0.00355, 77.92036494},
    {80, 0.02286, 79.91637808},
    {82, 0.11593, 81.91348273},
    {83, 0.11500, 82.91412716},
    {84, 0.56987, 83.9114977282},
    {86, 0.17279, 85.9106106269}
};

/* Z=37 Rubidium */
static const alea_isotope_t isotopes_Rb[] = {
    {85, 0.7217, 84.9117897379},
    {87, 0.2783, 86.9091805310}
};

/* Z=38 Strontium */
static const alea_isotope_t isotopes_Sr[] = {
    {84, 0.0056, 83.9134191},
    {86, 0.0986, 85.9092606},
    {87, 0.0700, 86.9088775},
    {88, 0.8258, 87.9056125}
};

/* Z=39 Yttrium */
static const alea_isotope_t isotopes_Y[] = {
    {89, 1.0, 88.9058403}
};

/* Z=40 Zirconium */
static const alea_isotope_t isotopes_Zr[] = {
    {90, 0.5145, 89.9046977},
    {91, 0.1122, 90.9056396},
    {92, 0.1715, 91.9050347},
    {94, 0.1738, 93.9063108},
    {96, 0.0280, 95.9082714}
};

/* Z=41 Niobium */
static const alea_isotope_t isotopes_Nb[] = {
    {93, 1.0, 92.9063730}
};

/* Z=42 Molybdenum */
static const alea_isotope_t isotopes_Mo[] = {
    {92, 0.1453, 91.90680796},
    {94, 0.0915, 93.90508490},
    {95, 0.1584, 94.90583877},
    {96, 0.1667, 95.90467612},
    {97, 0.0960, 96.90601812},
    {98, 0.2439, 97.90540482},
    {100, 0.0982, 99.9074718}
};

/* Z=43 Technetium - no stable isotopes */
static const alea_isotope_t isotopes_Tc[] = {
    {99, 0.0, 98.9062508}
};

/* Z=44 Ruthenium */
static const alea_isotope_t isotopes_Ru[] = {
    {96, 0.0554, 95.90759025},
    {98, 0.0187, 97.9052868},
    {99, 0.1276, 98.9059341},
    {100, 0.1260, 99.9042143},
    {101, 0.1706, 100.9055769},
    {102, 0.3155, 101.9043441},
    {104, 0.1862, 103.9054275}
};

/* Z=45 Rhodium */
static const alea_isotope_t isotopes_Rh[] = {
    {103, 1.0, 102.9054980}
};

/* Z=46 Palladium */
static const alea_isotope_t isotopes_Pd[] = {
    {102, 0.0102, 101.9056022},
    {104, 0.1114, 103.9040305},
    {105, 0.2233, 104.9050796},
    {106, 0.2733, 105.9034804},
    {108, 0.2646, 107.9038916},
    {110, 0.1172, 109.90517220}
};

/* Z=47 Silver */
static const alea_isotope_t isotopes_Ag[] = {
    {107, 0.51839, 106.9050916},
    {109, 0.48161, 108.9047553}
};

/* Z=48 Cadmium */
static const alea_isotope_t isotopes_Cd[] = {
    {106, 0.0125, 105.9064599},
    {108, 0.0089, 107.9041834},
    {110, 0.1249, 109.90300661},
    {111, 0.1280, 110.90418287},
    {112, 0.2413, 111.90276287},
    {113, 0.1222, 112.90440813},
    {114, 0.2873, 113.90336509},
    {116, 0.0749, 115.90476315}
};

/* Z=49 Indium */
static const alea_isotope_t isotopes_In[] = {
    {113, 0.0429, 112.90406184},
    {115, 0.9571, 114.903878776}
};

/* Z=50 Tin */
static const alea_isotope_t isotopes_Sn[] = {
    {112, 0.0097, 111.90482387},
    {114, 0.0066, 113.9027827},
    {115, 0.0034, 114.903344699},
    {116, 0.1454, 115.90174280},
    {117, 0.0768, 116.90295398},
    {118, 0.2422, 117.90160657},
    {119, 0.0859, 118.90331117},
    {120, 0.3258, 119.90220163},
    {122, 0.0463, 121.9034438},
    {124, 0.0579, 123.9052766}
};

/* Z=51 Antimony */
static const alea_isotope_t isotopes_Sb[] = {
    {121, 0.5721, 120.9038120},
    {123, 0.4279, 122.9042132}
};

/* Z=52 Tellurium */
static const alea_isotope_t isotopes_Te[] = {
    {120, 0.0009, 119.9040593},
    {122, 0.0255, 121.9030435},
    {123, 0.0089, 122.9042698},
    {124, 0.0474, 123.9028171},
    {125, 0.0707, 124.9044299},
    {126, 0.1884, 125.9033109},
    {128, 0.3174, 127.90446128},
    {130, 0.3408, 129.906222748}
};

/* Z=53 Iodine */
static const alea_isotope_t isotopes_I[] = {
    {127, 1.0, 126.9044719}
};

/* Z=54 Xenon */
static const alea_isotope_t isotopes_Xe[] = {
    {124, 0.000952, 123.9058920},
    {126, 0.000890, 125.9042983},
    {128, 0.019102, 127.9035310},
    {129, 0.264006, 128.9047808611},
    {130, 0.040710, 129.903509349},
    {131, 0.212324, 130.90508406},
    {132, 0.269086, 131.9041550856},
    {134, 0.104357, 133.90539466},
    {136, 0.088573, 135.907214484}
};

/* Z=55 Cesium */
static const alea_isotope_t isotopes_Cs[] = {
    {133, 1.0, 132.9054519610}
};

/* Z=56 Barium */
static const alea_isotope_t isotopes_Ba[] = {
    {130, 0.00106, 129.9063207},
    {132, 0.00101, 131.9050611},
    {134, 0.02417, 133.90450818},
    {135, 0.06592, 134.90568838},
    {136, 0.07854, 135.90457573},
    {137, 0.11232, 136.90582714},
    {138, 0.71698, 137.90524700}
};

/* Z=57 Lanthanum */
static const alea_isotope_t isotopes_La[] = {
    {138, 0.0008881, 137.9071149},
    {139, 0.9991119, 138.9063563}
};

/* Z=58 Cerium */
static const alea_isotope_t isotopes_Ce[] = {
    {136, 0.00185, 135.90712921},
    {138, 0.00251, 137.905991},
    {140, 0.88450, 139.9054431},
    {142, 0.11114, 141.9092504}
};

/* Z=59 Praseodymium */
static const alea_isotope_t isotopes_Pr[] = {
    {141, 1.0, 140.9076576}
};

/* Z=60 Neodymium */
static const alea_isotope_t isotopes_Nd[] = {
    {142, 0.27152, 141.9077290},
    {143, 0.12174, 142.9098200},
    {144, 0.23798, 143.9100930},
    {145, 0.08293, 144.9125793},
    {146, 0.17189, 145.9131226},
    {148, 0.05756, 147.9168993},
    {150, 0.05638, 149.9209022}
};

/* Z=61 Promethium - no stable isotopes */
static const alea_isotope_t isotopes_Pm[] = {
    {147, 0.0, 146.9151450}
};

/* Z=62 Samarium */
static const alea_isotope_t isotopes_Sm[] = {
    {144, 0.0307, 143.9120065},
    {147, 0.1499, 146.9149044},
    {148, 0.1124, 147.9148292},
    {149, 0.1382, 148.9171921},
    {150, 0.0738, 149.9172829},
    {152, 0.2675, 151.9197397},
    {154, 0.2275, 153.9222169}
};

/* Z=63 Europium */
static const alea_isotope_t isotopes_Eu[] = {
    {151, 0.4781, 150.9198578},
    {153, 0.5219, 152.9212380}
};

/* Z=64 Gadolinium */
static const alea_isotope_t isotopes_Gd[] = {
    {152, 0.0020, 151.9197995},
    {154, 0.0218, 153.9208741},
    {155, 0.1480, 154.9226305},
    {156, 0.2047, 155.9221312},
    {157, 0.1565, 156.9239686},
    {158, 0.2484, 157.9241123},
    {160, 0.2186, 159.9270624}
};

/* Z=65 Terbium */
static const alea_isotope_t isotopes_Tb[] = {
    {159, 1.0, 158.9253547}
};

/* Z=66 Dysprosium */
static const alea_isotope_t isotopes_Dy[] = {
    {156, 0.00056, 155.9242847},
    {158, 0.00095, 157.9244159},
    {160, 0.02329, 159.9252046},
    {161, 0.18889, 160.9269405},
    {162, 0.25475, 161.9268056},
    {163, 0.24896, 162.9287383},
    {164, 0.28260, 163.9291819}
};

/* Z=67 Holmium */
static const alea_isotope_t isotopes_Ho[] = {
    {165, 1.0, 164.9303288}
};

/* Z=68 Erbium */
static const alea_isotope_t isotopes_Er[] = {
    {162, 0.00139, 161.9287884},
    {164, 0.01601, 163.9292088},
    {166, 0.33503, 165.9302995},
    {167, 0.22869, 166.9320546},
    {168, 0.26978, 167.9323767},
    {170, 0.14910, 169.9354702}
};

/* Z=69 Thulium */
static const alea_isotope_t isotopes_Tm[] = {
    {169, 1.0, 168.9342179}
};

/* Z=70 Ytterbium */
static const alea_isotope_t isotopes_Yb[] = {
    {168, 0.00123, 167.9338896},
    {170, 0.02982, 169.9347664},
    {171, 0.1409, 170.9363302},
    {172, 0.2168, 171.9363859},
    {173, 0.16103, 172.9382151},
    {174, 0.32026, 173.9388664},
    {176, 0.12996, 175.9425764}
};

/* Z=71 Lutetium */
static const alea_isotope_t isotopes_Lu[] = {
    {175, 0.97401, 174.9407752},
    {176, 0.02599, 175.9426897}
};

/* Z=72 Hafnium */
static const alea_isotope_t isotopes_Hf[] = {
    {174, 0.0016, 173.9400461},
    {176, 0.0526, 175.9414076},
    {177, 0.1860, 176.9432277},
    {178, 0.2728, 177.9437058},
    {179, 0.1362, 178.9458232},
    {180, 0.3508, 179.9465570}
};

/* Z=73 Tantalum */
static const alea_isotope_t isotopes_Ta[] = {
    {180, 0.0001201, 179.9474648},
    {181, 0.9998799, 180.9479958}
};

/* Z=74 Tungsten */
static const alea_isotope_t isotopes_W[] = {
    {180, 0.0012, 179.9467108},
    {182, 0.2650, 181.94820394},
    {183, 0.1431, 182.95022275},
    {184, 0.3064, 183.95093092},
    {186, 0.2843, 185.9543628}
};

/* Z=75 Rhenium */
static const alea_isotope_t isotopes_Re[] = {
    {185, 0.3740, 184.9529545},
    {187, 0.6260, 186.9557501}
};

/* Z=76 Osmium */
static const alea_isotope_t isotopes_Os[] = {
    {184, 0.0002, 183.9524885},
    {186, 0.0159, 185.9538350},
    {187, 0.0196, 186.9557474},
    {188, 0.1324, 187.9558352},
    {189, 0.1615, 188.9581442},
    {190, 0.2626, 189.9584437},
    {192, 0.4078, 191.9614770}
};

/* Z=77 Iridium */
static const alea_isotope_t isotopes_Ir[] = {
    {191, 0.373, 190.9605893},
    {193, 0.627, 192.9629216}
};

/* Z=78 Platinum */
static const alea_isotope_t isotopes_Pt[] = {
    {190, 0.00012, 189.9599297},
    {192, 0.00782, 191.9610387},
    {194, 0.3286, 193.9626809},
    {195, 0.3378, 194.9647917},
    {196, 0.2521, 195.96495209},
    {198, 0.07356, 197.9678949}
};

/* Z=79 Gold */
static const alea_isotope_t isotopes_Au[] = {
    {197, 1.0, 196.96656879}
};

/* Z=80 Mercury */
static const alea_isotope_t isotopes_Hg[] = {
    {196, 0.0015, 195.9658326},
    {198, 0.0997, 197.96676860},
    {199, 0.1687, 198.96828064},
    {200, 0.2310, 199.96832659},
    {201, 0.1318, 200.97030284},
    {202, 0.2986, 201.97064340},
    {204, 0.0687, 203.97349398}
};

/* Z=81 Thallium */
static const alea_isotope_t isotopes_Tl[] = {
    {203, 0.2952, 202.9723446},
    {205, 0.7048, 204.9744278}
};

/* Z=82 Lead */
static const alea_isotope_t isotopes_Pb[] = {
    {204, 0.014, 203.9730440},
    {206, 0.241, 205.9744657},
    {207, 0.221, 206.9758973},
    {208, 0.524, 207.9766525}
};

/* Z=83 Bismuth */
static const alea_isotope_t isotopes_Bi[] = {
    {209, 1.0, 208.9803991}
};

/* Z=84 Polonium - no stable isotopes */
static const alea_isotope_t isotopes_Po[] = {
    {210, 0.0, 209.9828741}
};

/* Z=85 Astatine - no stable isotopes */
static const alea_isotope_t isotopes_At[] = {
    {210, 0.0, 209.9871479}
};

/* Z=86 Radon - no stable isotopes */
static const alea_isotope_t isotopes_Rn[] = {
    {222, 0.0, 222.0175782}
};

/* Z=87 Francium - no stable isotopes */
static const alea_isotope_t isotopes_Fr[] = {
    {223, 0.0, 223.0197360}
};

/* Z=88 Radium - no stable isotopes */
static const alea_isotope_t isotopes_Ra[] = {
    {226, 0.0, 226.0254103}
};

/* Z=89 Actinium - no stable isotopes */
static const alea_isotope_t isotopes_Ac[] = {
    {227, 0.0, 227.0277523}
};

/* Z=90 Thorium */
static const alea_isotope_t isotopes_Th[] = {
    {232, 1.0, 232.0380558}
};

/* Z=91 Protactinium */
static const alea_isotope_t isotopes_Pa[] = {
    {231, 1.0, 231.0358842}
};

/* Z=92 Uranium */
static const alea_isotope_t isotopes_U[] = {
    {234, 0.000054, 234.0409523},
    {235, 0.007204, 235.0439301},
    {238, 0.992742, 238.0507884}
};

/* Z=93 Neptunium - no stable isotopes */
static const alea_isotope_t isotopes_Np[] = {
    {237, 0.0, 237.0481736}
};

/* Z=94 Plutonium - no stable isotopes */
static const alea_isotope_t isotopes_Pu[] = {
    {239, 0.0, 239.0521636},
    {240, 0.0, 240.0538138},
    {241, 0.0, 241.0568517},
    {242, 0.0, 242.0587428}
};

/* Z=95 Americium - no stable isotopes */
static const alea_isotope_t isotopes_Am[] = {
    {241, 0.0, 241.0568293},
    {243, 0.0, 243.0613813}
};

/* Z=96 Curium - no stable isotopes */
static const alea_isotope_t isotopes_Cm[] = {
    {244, 0.0, 244.0627528},
    {245, 0.0, 245.0654915},
    {246, 0.0, 246.0672238}
};

/* Z=97 Berkelium - no stable isotopes */
static const alea_isotope_t isotopes_Bk[] = {
    {249, 0.0, 249.0749877}
};

/* Z=98 Californium - no stable isotopes */
static const alea_isotope_t isotopes_Cf[] = {
    {252, 0.0, 252.0816272}
};

/* Z=99 Einsteinium - no stable isotopes */
static const alea_isotope_t isotopes_Es[] = {
    {252, 0.0, 252.082980}
};

/* Z=100 Fermium - no stable isotopes */
static const alea_isotope_t isotopes_Fm[] = {
    {257, 0.0, 257.0951061}
};

/* Z=101 Mendelevium - no stable isotopes */
static const alea_isotope_t isotopes_Md[] = {
    {258, 0.0, 258.0984315}
};

/* Z=102 Nobelium - no stable isotopes */
static const alea_isotope_t isotopes_No[] = {
    {259, 0.0, 259.10103}
};

/* Z=103 Lawrencium - no stable isotopes */
static const alea_isotope_t isotopes_Lr[] = {
    {262, 0.0, 262.10961}
};

/* Z=104 Rutherfordium - no stable isotopes */
static const alea_isotope_t isotopes_Rf[] = {
    {267, 0.0, 267.12179}
};

/* Z=105 Dubnium - no stable isotopes */
static const alea_isotope_t isotopes_Db[] = {
    {268, 0.0, 268.12567}
};

/* Z=106 Seaborgium - no stable isotopes */
static const alea_isotope_t isotopes_Sg[] = {
    {271, 0.0, 271.13393}
};

/* Z=107 Bohrium - no stable isotopes */
static const alea_isotope_t isotopes_Bh[] = {
    {272, 0.0, 272.13826}
};

/* Z=108 Hassium - no stable isotopes */
static const alea_isotope_t isotopes_Hs[] = {
    {270, 0.0, 270.13429}
};

/* Z=109 Meitnerium - no stable isotopes */
static const alea_isotope_t isotopes_Mt[] = {
    {276, 0.0, 276.15159}
};

/* Z=110 Darmstadtium - no stable isotopes */
static const alea_isotope_t isotopes_Ds[] = {
    {281, 0.0, 281.16451}
};

/* Z=111 Roentgenium - no stable isotopes */
static const alea_isotope_t isotopes_Rg[] = {
    {280, 0.0, 280.16514}
};

/* Z=112 Copernicium - no stable isotopes */
static const alea_isotope_t isotopes_Cn[] = {
    {285, 0.0, 285.17712}
};

/* Z=113 Nihonium - no stable isotopes */
static const alea_isotope_t isotopes_Nh[] = {
    {284, 0.0, 284.17873}
};

/* Z=114 Flerovium - no stable isotopes */
static const alea_isotope_t isotopes_Fl[] = {
    {289, 0.0, 289.19042}
};

/* Z=115 Moscovium - no stable isotopes */
static const alea_isotope_t isotopes_Mc[] = {
    {288, 0.0, 288.19274}
};

/* Z=116 Livermorium - no stable isotopes */
static const alea_isotope_t isotopes_Lv[] = {
    {293, 0.0, 293.20449}
};

/* Z=117 Tennessine - no stable isotopes */
static const alea_isotope_t isotopes_Ts[] = {
    {292, 0.0, 292.20746}
};

/* Z=118 Oganesson - no stable isotopes */
static const alea_isotope_t isotopes_Og[] = {
    {294, 0.0, 294.21392}
};


/* ============================================================================
 * ELEMENT TABLE
 * ============================================================================ */

#define ELEM(z, sym, name, weight, iso_array) \
    {z, sym, name, weight, (alea_isotope_t*)iso_array, sizeof(iso_array)/sizeof(iso_array[0])}

const alea_element_t g_elements[119] = {
    {0, "", "", 0.0, NULL, 0},  /* Placeholder for index 0 */
    ELEM(1,   "H",  "Hydrogen",      1.008,       isotopes_H),
    ELEM(2,   "He", "Helium",        4.002602,    isotopes_He),
    ELEM(3,   "Li", "Lithium",       6.94,        isotopes_Li),
    ELEM(4,   "Be", "Beryllium",     9.0121831,   isotopes_Be),
    ELEM(5,   "B",  "Boron",         10.81,       isotopes_B),
    ELEM(6,   "C",  "Carbon",        12.011,      isotopes_C),
    ELEM(7,   "N",  "Nitrogen",      14.007,      isotopes_N),
    ELEM(8,   "O",  "Oxygen",        15.999,      isotopes_O),
    ELEM(9,   "F",  "Fluorine",      18.998403163, isotopes_F),
    ELEM(10,  "Ne", "Neon",          20.1797,     isotopes_Ne),
    ELEM(11,  "Na", "Sodium",        22.98976928, isotopes_Na),
    ELEM(12,  "Mg", "Magnesium",     24.305,      isotopes_Mg),
    ELEM(13,  "Al", "Aluminum",      26.9815385,  isotopes_Al),
    ELEM(14,  "Si", "Silicon",       28.085,      isotopes_Si),
    ELEM(15,  "P",  "Phosphorus",    30.973761998, isotopes_P),
    ELEM(16,  "S",  "Sulfur",        32.06,       isotopes_S),
    ELEM(17,  "Cl", "Chlorine",      35.45,       isotopes_Cl),
    ELEM(18,  "Ar", "Argon",         39.948,      isotopes_Ar),
    ELEM(19,  "K",  "Potassium",     39.0983,     isotopes_K),
    ELEM(20,  "Ca", "Calcium",       40.078,      isotopes_Ca),
    ELEM(21,  "Sc", "Scandium",      44.955908,   isotopes_Sc),
    ELEM(22,  "Ti", "Titanium",      47.867,      isotopes_Ti),
    ELEM(23,  "V",  "Vanadium",      50.9415,     isotopes_V),
    ELEM(24,  "Cr", "Chromium",      51.9961,     isotopes_Cr),
    ELEM(25,  "Mn", "Manganese",     54.938044,   isotopes_Mn),
    ELEM(26,  "Fe", "Iron",          55.845,      isotopes_Fe),
    ELEM(27,  "Co", "Cobalt",        58.933194,   isotopes_Co),
    ELEM(28,  "Ni", "Nickel",        58.6934,     isotopes_Ni),
    ELEM(29,  "Cu", "Copper",        63.546,      isotopes_Cu),
    ELEM(30,  "Zn", "Zinc",          65.38,       isotopes_Zn),
    ELEM(31,  "Ga", "Gallium",       69.723,      isotopes_Ga),
    ELEM(32,  "Ge", "Germanium",     72.630,      isotopes_Ge),
    ELEM(33,  "As", "Arsenic",       74.921595,   isotopes_As),
    ELEM(34,  "Se", "Selenium",      78.971,      isotopes_Se),
    ELEM(35,  "Br", "Bromine",       79.904,      isotopes_Br),
    ELEM(36,  "Kr", "Krypton",       83.798,      isotopes_Kr),
    ELEM(37,  "Rb", "Rubidium",      85.4678,     isotopes_Rb),
    ELEM(38,  "Sr", "Strontium",     87.62,       isotopes_Sr),
    ELEM(39,  "Y",  "Yttrium",       88.90584,    isotopes_Y),
    ELEM(40,  "Zr", "Zirconium",     91.224,      isotopes_Zr),
    ELEM(41,  "Nb", "Niobium",       92.90637,    isotopes_Nb),
    ELEM(42,  "Mo", "Molybdenum",    95.95,       isotopes_Mo),
    ELEM(43,  "Tc", "Technetium",    98.0,        isotopes_Tc),
    ELEM(44,  "Ru", "Ruthenium",     101.07,      isotopes_Ru),
    ELEM(45,  "Rh", "Rhodium",       102.90550,   isotopes_Rh),
    ELEM(46,  "Pd", "Palladium",     106.42,      isotopes_Pd),
    ELEM(47,  "Ag", "Silver",        107.8682,    isotopes_Ag),
    ELEM(48,  "Cd", "Cadmium",       112.414,     isotopes_Cd),
    ELEM(49,  "In", "Indium",        114.818,     isotopes_In),
    ELEM(50,  "Sn", "Tin",           118.710,     isotopes_Sn),
    ELEM(51,  "Sb", "Antimony",      121.760,     isotopes_Sb),
    ELEM(52,  "Te", "Tellurium",     127.60,      isotopes_Te),
    ELEM(53,  "I",  "Iodine",        126.90447,   isotopes_I),
    ELEM(54,  "Xe", "Xenon",         131.293,     isotopes_Xe),
    ELEM(55,  "Cs", "Cesium",        132.90545196, isotopes_Cs),
    ELEM(56,  "Ba", "Barium",        137.327,     isotopes_Ba),
    ELEM(57,  "La", "Lanthanum",     138.90547,   isotopes_La),
    ELEM(58,  "Ce", "Cerium",        140.116,     isotopes_Ce),
    ELEM(59,  "Pr", "Praseodymium",  140.90766,   isotopes_Pr),
    ELEM(60,  "Nd", "Neodymium",     144.242,     isotopes_Nd),
    ELEM(61,  "Pm", "Promethium",    145.0,       isotopes_Pm),
    ELEM(62,  "Sm", "Samarium",      150.36,      isotopes_Sm),
    ELEM(63,  "Eu", "Europium",      151.964,     isotopes_Eu),
    ELEM(64,  "Gd", "Gadolinium",    157.25,      isotopes_Gd),
    ELEM(65,  "Tb", "Terbium",       158.92535,   isotopes_Tb),
    ELEM(66,  "Dy", "Dysprosium",    162.500,     isotopes_Dy),
    ELEM(67,  "Ho", "Holmium",       164.93033,   isotopes_Ho),
    ELEM(68,  "Er", "Erbium",        167.259,     isotopes_Er),
    ELEM(69,  "Tm", "Thulium",       168.93422,   isotopes_Tm),
    ELEM(70,  "Yb", "Ytterbium",     173.045,     isotopes_Yb),
    ELEM(71,  "Lu", "Lutetium",      174.9668,    isotopes_Lu),
    ELEM(72,  "Hf", "Hafnium",       178.49,      isotopes_Hf),
    ELEM(73,  "Ta", "Tantalum",      180.94788,   isotopes_Ta),
    ELEM(74,  "W",  "Tungsten",      183.84,      isotopes_W),
    ELEM(75,  "Re", "Rhenium",       186.207,     isotopes_Re),
    ELEM(76,  "Os", "Osmium",        190.23,      isotopes_Os),
    ELEM(77,  "Ir", "Iridium",       192.217,     isotopes_Ir),
    ELEM(78,  "Pt", "Platinum",      195.084,     isotopes_Pt),
    ELEM(79,  "Au", "Gold",          196.966569,  isotopes_Au),
    ELEM(80,  "Hg", "Mercury",       200.592,     isotopes_Hg),
    ELEM(81,  "Tl", "Thallium",      204.38,      isotopes_Tl),
    ELEM(82,  "Pb", "Lead",          207.2,       isotopes_Pb),
    ELEM(83,  "Bi", "Bismuth",       208.98040,   isotopes_Bi),
    ELEM(84,  "Po", "Polonium",      209.0,       isotopes_Po),
    ELEM(85,  "At", "Astatine",      210.0,       isotopes_At),
    ELEM(86,  "Rn", "Radon",         222.0,       isotopes_Rn),
    ELEM(87,  "Fr", "Francium",      223.0,       isotopes_Fr),
    ELEM(88,  "Ra", "Radium",        226.0,       isotopes_Ra),
    ELEM(89,  "Ac", "Actinium",      227.0,       isotopes_Ac),
    ELEM(90,  "Th", "Thorium",       232.0377,    isotopes_Th),
    ELEM(91,  "Pa", "Protactinium",  231.03588,   isotopes_Pa),
    ELEM(92,  "U",  "Uranium",       238.02891,   isotopes_U),
    ELEM(93,  "Np", "Neptunium",     237.0,       isotopes_Np),
    ELEM(94,  "Pu", "Plutonium",     244.0,       isotopes_Pu),
    ELEM(95,  "Am", "Americium",     243.0,       isotopes_Am),
    ELEM(96,  "Cm", "Curium",        247.0,       isotopes_Cm),
    ELEM(97,  "Bk", "Berkelium",     247.0,       isotopes_Bk),
    ELEM(98,  "Cf", "Californium",   251.0,       isotopes_Cf),
    ELEM(99,  "Es", "Einsteinium",   252.0,       isotopes_Es),
    ELEM(100, "Fm", "Fermium",       257.0,       isotopes_Fm),
    ELEM(101, "Md", "Mendelevium",   258.0,       isotopes_Md),
    ELEM(102, "No", "Nobelium",      259.0,       isotopes_No),
    ELEM(103, "Lr", "Lawrencium",    266.0,       isotopes_Lr),
    ELEM(104, "Rf", "Rutherfordium", 267.0,       isotopes_Rf),
    ELEM(105, "Db", "Dubnium",       268.0,       isotopes_Db),
    ELEM(106, "Sg", "Seaborgium",    269.0,       isotopes_Sg),
    ELEM(107, "Bh", "Bohrium",       270.0,       isotopes_Bh),
    ELEM(108, "Hs", "Hassium",       269.0,       isotopes_Hs),
    ELEM(109, "Mt", "Meitnerium",    278.0,       isotopes_Mt),
    ELEM(110, "Ds", "Darmstadtium",  281.0,       isotopes_Ds),
    ELEM(111, "Rg", "Roentgenium",   282.0,       isotopes_Rg),
    ELEM(112, "Cn", "Copernicium",   285.0,       isotopes_Cn),
    ELEM(113, "Nh", "Nihonium",      286.0,       isotopes_Nh),
    ELEM(114, "Fl", "Flerovium",     289.0,       isotopes_Fl),
    ELEM(115, "Mc", "Moscovium",     290.0,       isotopes_Mc),
    ELEM(116, "Lv", "Livermorium",   293.0,       isotopes_Lv),
    ELEM(117, "Ts", "Tennessine",    294.0,       isotopes_Ts),
    ELEM(118, "Og", "Oganesson",     294.0,       isotopes_Og)
};

const size_t g_elements_count = 119;

#undef ELEM
