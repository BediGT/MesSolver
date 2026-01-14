#pragma once
#include <vector>

namespace TestHelpers
{
	bool AreMatricesEqual(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B, double eps);
    bool AreVectorsEqual(const std::vector<double>& A, const std::vector<double>& B, double eps);
};

