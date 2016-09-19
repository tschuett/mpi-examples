/*

quantile 95° confidence level 95°
n  j  k  p

n≤58: no confidence interval possible.

n≥ 124      ≈⌊0.95n−0.427√n⌋   ≈⌈0.95n+1+0.427√n⌉       0.950
*/

std::vector<std::tuple<int, int, int, double>> table_95_95 = {
    {59, 50, 59, 0.951},    {60, 52, 60, 0.951},    {61, 53, 61, 0.953},
    {62, 54, 62, 0.955},    {63, 55, 63, 0.957},    {64, 56, 64, 0.958},
    {65, 57, 65, 0.959},    {66, 58, 66, 0.961},    {67, 59, 67, 0.962},
    {68, 60, 68, 0.963},    {69, 61, 69, 0.964},    {70, 62, 70, 0.964},
    {71, 63, 71, 0.965},    {72, 64, 72, 0.965},    {73, 65, 73, 0.966},
    {74, 66, 74, 0.966},    {75, 67, 75, 0.966},    {76, 68, 76, 0.966},
    {77, 69, 77, 0.966},    {78, 70, 78, 0.966},    {79, 71, 79, 0.966},
    {80, 72, 80, 0.965},    {81, 73, 81, 0.964},    {82, 74, 82, 0.964},
    {83, 75, 83, 0.963},    {84, 76, 84, 0.962},    {85, 77, 85, 0.961},
    {86, 78, 86, 0.960},    {87, 79, 87, 0.959},    {88, 80, 88, 0.957},
    {89, 81, 89, 0.956},    {90, 82, 90, 0.954},    {91, 83, 91, 0.952},
    {92, 84, 92, 0.950},    {93, 84, 93, 0.974},    {94, 85, 94, 0.973},
    {95, 86, 95, 0.972},    {96, 87, 96, 0.971},    {97, 88, 97, 0.970},
    {98, 89, 98, 0.969},    {99, 90, 99, 0.967},    {100, 91, 100, 0.966},
    {101, 91, 100, 0.952},  {102, 92, 101, 0.953},  {103, 93, 102, 0.953},
    {104, 94, 103, 0.954},  {105, 95, 104, 0.954},  {106, 96, 105, 0.954},
    {107, 97, 106, 0.954},  {108, 98, 107, 0.954},  {109, 99, 108, 0.954},
    {110, 100, 109, 0.954}, {111, 101, 110, 0.954}, {112, 102, 111, 0.953},
    {113, 103, 112, 0.953}, {114, 104, 113, 0.952}, {115, 105, 114, 0.951},
    {116, 106, 115, 0.950}, {117, 107, 117, 0.965}, {118, 108, 118, 0.963},
    {119, 109, 119, 0.961}, {120, 110, 120, 0.959}, {121, 110, 120, 0.967},
    {122, 111, 121, 0.966}, {123, 112, 122, 0.966}};