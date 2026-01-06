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

    std::vector<double> m_dN_dx{ ELEMENT_POINTS };
	std::vector<double> m_dN_dy{ ELEMENT_POINTS };

    Jacobian() = default;
    Jacobian(int i, const Node& node);

    void Print() const;
};
