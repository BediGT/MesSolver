#include "Element.h"

#include "Jacobian.h"
#include "Node.h"
#include "Edge.h"
#include "Helpers.h"
#include "ShapeFunctions.h"

#include <iostream>
#include <iomanip>

Element::Element(Node* arg1, Node* arg2, Node* arg3, Node* arg4)
{
	m_vertices.push_back(arg1);
	m_vertices.push_back(arg2);
	m_vertices.push_back(arg3);
	m_vertices.push_back(arg4);
}

void Element::CalculateJacobians()
{
    if (auto universalElement = UniversalElement::Get())
    {
		m_jacobians.clear();
        for (int i = 0; i < universalElement->GetIntegralPointsNumber(); ++i)
            m_jacobians.emplace_back(i, m_vertices);
    }
}

void Element::CalculateMatrixH(double conductivity)
{
    if (auto universalElement = UniversalElement::Get())
    {
        std::vector<double> weights = getGaussWeights(universalElement->GetGaussNumber());
        int n = universalElement->GetGaussNumber();
        int integralPointIndex = 0;

        m_matrixH.clear();
        m_matrixH.resize(ELEMENT_POINTS, std::vector<double>(ELEMENT_POINTS, 0.0));
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                std::vector<std::vector<double>> partOfMatrixH = MatrixUtils::MultiplyVectors(m_jacobians.at(integralPointIndex).m_dN_dx, m_jacobians.at(integralPointIndex).m_dN_dx);
                MatrixUtils::AddTo(partOfMatrixH, MatrixUtils::MultiplyVectors(m_jacobians.at(integralPointIndex).m_dN_dy, m_jacobians.at(integralPointIndex).m_dN_dy));
                MatrixUtils::Scale(partOfMatrixH, m_jacobians.at(integralPointIndex).determinant * weights.at(i) * weights.at(j) * conductivity);

                MatrixUtils::AddTo(m_matrixH, partOfMatrixH);

                integralPointIndex++;
            }
        }
	}
}

void Element::CalculateBoundryConditionsH(double alfa, double ambientTemp)
{
    m_vectorP.clear();
    m_vectorP.resize(ELEMENT_POINTS, 0.0);

    m_edges.clear();
    for (int i = 0; i < ELEMENT_POINTS; ++i)
        m_edges.emplace_back(m_vertices.at(i), m_vertices.at((i + 1) % ELEMENT_POINTS), static_cast<EEdgeAligment>(i));

    for (auto& edge : m_edges)
    {
        if (edge.m_bBoundaryCondition)
        {
            edge.CalculateMatrixHbcAndVectorP(alfa, ambientTemp);
            MatrixUtils::AddTo(m_matrixH, edge.m_matrixHbc);
            MatrixUtils::AddTo(m_vectorP, edge.m_vectorP);
        }
    }
}

void Element::CalculateMatrixC(double specificHeat, double density)
{
    if (auto universalElement = UniversalElement::Get())
    {
        std::vector<double> weights = getGaussWeights(universalElement->GetGaussNumber());
        std::vector<double> values = getGaussValues(universalElement->GetGaussNumber());
        int n = universalElement->GetGaussNumber();
        int integralPointIndex = 0;

        m_matrixC.clear();
        m_matrixC.resize(ELEMENT_POINTS, std::vector<double>(ELEMENT_POINTS, 0.0));
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                std::vector<double> shapeFunctionsValue(ELEMENT_POINTS);
                for (int k = 0; k < ELEMENT_POINTS; ++k)
					shapeFunctionsValue.at(k) = (shapeFunction(k, values.at(j), values.at(i)));

                std::vector<std::vector<double>> partOfMatrixC = MatrixUtils::MultiplyVectors(shapeFunctionsValue, shapeFunctionsValue);
                MatrixUtils::Scale(partOfMatrixC, m_jacobians.at(integralPointIndex).determinant * weights.at(i) * weights.at(j) * specificHeat * density);

                MatrixUtils::AddTo(m_matrixC, partOfMatrixC);

                integralPointIndex++;
            }
        }
    }
}

void Element::PrintIds() const
{
    for (auto node : m_vertices)
        std::cout << std::setw(2) << node->m_id << " ";
}

void Element::PrintJacobians() const
{
    std::cout << "\nElement Nodes Ids: ";
    for (auto node : m_vertices)
        std::cout << std::setw(2) << node->m_id << " ";
    std::cout << "\n";

    int i = 0;
    for (const auto& jacobian : m_jacobians)
    {
		std::cout << "\nJacobian " << i << ":";
        jacobian.Print();
        i++;
    }
}