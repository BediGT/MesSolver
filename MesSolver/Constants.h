#pragma once
#include <vector>

const int ELEMENT_POINTS = 4;
const int ELEMENT_DIMENSION = 2;

std::vector<double> getGaussWeights(int n);
std::vector<double> getGaussValues(int n);