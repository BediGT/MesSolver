#pragma once
#include <vector>

class Integrater
{
public:
	static double Integrate1D(const std::vector<double>& values, double dx, int n);
	static std::vector<std::vector<double>> Integrate2D(const std::vector<std::vector<double>>& valuesA, const std::vector<std::vector<double>>& valuesB, std::vector<double> dz, int n);
};

