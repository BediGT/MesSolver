#pragma once
#include <vector>

#include "UniversalElement.h"
#include "Constants.h"

class Jacobian
{
public:
    double matrix[ELEMENT_DIMENSION][ELEMENT_DIMENSION]{};
    double matrixTransposed[ELEMENT_DIMENSION][ELEMENT_DIMENSION]{};
    double determinant{};

    std::vector<double> m_dN_dx{};
	std::vector<double> m_dN_dy{};

    Jacobian() = default;
    Jacobian(int i, const std::vector<Node*>& vertices);

    void Print() const;
};
