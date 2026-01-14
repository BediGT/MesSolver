#pragma once
#include <vector>
#include <iomanip>

const int ELEMENT_POINTS = 4;
const int ELEMENT_DIMENSION = 2;

const std::streamsize PRINT_PRECISION = 3L;
const std::streamsize PRINT_WIDTH = 9L;

std::vector<double> getGaussWeights(int n);
std::vector<double> getGaussValues(int n);