#include "Jacobian.h"
#include <iostream>

#include "Node.h"

Jacobian::Jacobian(int i, const Node& node)
{
    UniversalElement* universalElement = UniversalElement::Get();

    for (int k = 0; k < ELEMENT_POINTS; ++k)
    {
        matrix[0][0] += universalElement->m_dN_dKsi[i][k] * node.m_x;
        matrix[0][1] += universalElement->m_dN_dKsi[i][k] * node.m_y;
        matrix[1][0] += universalElement->m_dN_dEta[i][k] * node.m_x;
        matrix[1][1] += universalElement->m_dN_dEta[i][k] * node.m_y;
    }

    determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

    matrixTransposed[1][1] = matrix[0][0] / determinant;
    matrixTransposed[1][0] = -matrix[1][0] / determinant;
    matrixTransposed[0][0] = matrix[1][1] / determinant;
    matrixTransposed[0][1] = -matrix[0][1] / determinant;

    for (int j = 0; j < ELEMENT_POINTS; ++j)
    {
        m_dN_dx.at(j) = (matrixTransposed[0][0] * universalElement->m_dN_dKsi[i][j] + matrixTransposed[0][1] * universalElement->m_dN_dEta[i][j]);
        m_dN_dy.at(j) = (matrixTransposed[1][0] * universalElement->m_dN_dKsi[i][j] + matrixTransposed[1][1] * universalElement->m_dN_dEta[i][j]);
    }
}

void Jacobian::Print() const
{
    std::cout << "Jacobian Matrix:\n";
    std::cout << "[" << matrix[0][0] << ", " << matrix[0][1] << "]\n";
    std::cout << "[" << matrix[1][0] << ", " << matrix[1][1] << "]\n";
    std::cout << "Determinant: " << determinant << "\n\n";

    std::cout << "dN/dx\n";
    for (auto val : m_dN_dx)
        std::cout << val << " ";
    std::cout << '\n';

    std::cout << "dN/dy\n";
    for (auto val : m_dN_dy)
        std::cout << val << " ";
    std::cout << '\n\n\n';
}