/*
n  j  k  p

n≤7: no confidence interval possible.

n≥73   ≈⌊0.50n−1.288√n⌋   ≈⌈0.50n+1+1.288√n⌉  0.990
 */

std::vector<std::tuple<int, int, int, double>> table_50_99 = {

    {8, 1, 8, 0.992},    {9, 1, 9, 0.996},    {10, 1, 10, 0.998},
    {11, 1, 11, 0.999},  {12, 2, 11, 0.994},  {13, 2, 12, 0.997},
    {14, 2, 12, 0.993},  {15, 3, 13, 0.993},  {16, 3, 14, 0.996},
    {17, 3, 15, 0.998},  {18, 4, 15, 0.992},  {19, 4, 16, 0.996},
    {20, 4, 16, 0.993},  {21, 5, 17, 0.993},  {22, 5, 18, 0.996},
    {23, 5, 19, 0.997},  {24, 6, 19, 0.993},  {25, 6, 20, 0.996},
    {26, 7, 20, 0.991},  {27, 7, 21, 0.994},  {28, 7, 21, 0.992},
    {29, 8, 22, 0.992},  {30, 8, 23, 0.995},  {31, 8, 24, 0.997},
    {32, 9, 24, 0.993},  {33, 9, 25, 0.995},  {34, 10, 25, 0.991},
    {35, 10, 26, 0.994}, {36, 10, 26, 0.992}, {37, 11, 27, 0.992},
    {38, 11, 27, 0.991}, {39, 12, 28, 0.991}, {40, 12, 29, 0.994},
    {41, 12, 30, 0.996}, {42, 13, 30, 0.992}, {43, 13, 31, 0.995},
    {44, 14, 31, 0.990}, {45, 14, 32, 0.993}, {46, 15, 33, 0.992},
    {47, 15, 33, 0.992}, {48, 15, 33, 0.991}, {49, 16, 34, 0.991},
    {50, 16, 35, 0.993}, {51, 16, 36, 0.995}, {52, 17, 36, 0.992},
    {53, 17, 37, 0.995}, {54, 18, 37, 0.991}, {55, 18, 38, 0.994},
    {56, 18, 38, 0.992}, {57, 19, 39, 0.992}, {58, 20, 40, 0.991},
    {59, 20, 40, 0.991}, {60, 20, 40, 0.990}, {61, 21, 41, 0.990},
    {62, 21, 42, 0.993}, {63, 21, 43, 0.995}, {64, 22, 43, 0.992},
    {65, 22, 44, 0.994}, {66, 23, 44, 0.991}, {67, 23, 45, 0.993},
    {68, 23, 45, 0.992}, {69, 24, 46, 0.992}, {70, 24, 46, 0.991},
    {71, 25, 47, 0.991}, {72, 25, 47, 0.990}};