#include "Jacobian.h"
#include <iostream>

#include "Node.h"

Jacobian::Jacobian(int i, const std::vector<Node*>& vertices)
{
    if (UniversalElement* universalElement = UniversalElement::Get())
    {
        for (int k = 0; k < ELEMENT_POINTS; ++k)
        {
            matrix[0][1] += universalElement->m_dN_dKsi[i][k] * vertices.at(k)->m_y;
            matrix[1][0] += universalElement->m_dN_dEta[i][k] * vertices.at(k)->m_x;
            matrix[1][1] += universalElement->m_dN_dEta[i][k] * vertices.at(k)->m_y;
            matrix[0][0] += universalElement->m_dN_dKsi[i][k] * vertices.at(k)->m_x;
        }

        determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

        matrixTransposed[1][1] = matrix[0][0] / determinant;
        matrixTransposed[1][0] = -matrix[1][0] / determinant;
        matrixTransposed[0][0] = matrix[1][1] / determinant;
        matrixTransposed[0][1] = -matrix[0][1] / determinant;

        for (int j = 0; j < ELEMENT_POINTS; ++j)
        {
            m_dN_dx.emplace_back(matrixTransposed[0][0] * universalElement->m_dN_dKsi.at(i).at(j) + matrixTransposed[0][1] * universalElement->m_dN_dEta.at(i).at(j));
            m_dN_dy.emplace_back(matrixTransposed[1][0] * universalElement->m_dN_dKsi.at(i).at(j) + matrixTransposed[1][1] * universalElement->m_dN_dEta.at(i).at(j));
        }
    }
}

void Jacobian::Print() const
{
    std::cout << "\nJacobian Matrix:\n";
    std::cout << "[" << matrix[0][0] << ", " << matrix[0][1] << "]\n";
    std::cout << "[" << matrix[1][0] << ", " << matrix[1][1] << "]\n";
    std::cout << "Determinant: " << determinant << "\n\n";

    std::cout << "dN/dx\n";
    for (const auto& val : m_dN_dx)
        std::cout << val << " ";
    std::cout << '\n';

    std::cout << "dN/dy\n";
    for (const auto& val : m_dN_dy)
        std::cout << val << " ";
    std::cout << "\n";
}