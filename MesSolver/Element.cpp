#include "Element.h"

#include "Jacobian.h"
#include "Node.h"
#include "Edge.h"
#include "Integrater.h"
#include "Helpers.h"

#include <iostream>
#include <iomanip>

static double shapeFunction(int index, double ksi, double eta)
{
    switch (index)
    {
        // 0 - N1, 1 - N2, 2 - N3, 3 - N4
        case 0:
        {
            return 0.25 * (1 - ksi) * (1 - eta);
            break;
        }
        case 1:
        {
            return 0.25 * (1 + ksi) * (1 - eta);
            break;
        }
        case 2:
        {
            return 0.25 * (1 + ksi) * (1 + eta);
            break;
        }
        case 3:
        {
            return 0.25 * (1 - ksi) * (1 + eta);
            break;
        }
        default:
        {
            return 0.0;
            break;
        }
    }
}

Element::Element(Node* arg1, Node* arg2, Node* arg3, Node* arg4)
{
	m_nodes.push_back(arg1);
	m_nodes.push_back(arg2);
	m_nodes.push_back(arg3);
	m_nodes.push_back(arg4);

    for (int i = 0; i < ELEMENT_POINTS; ++i)
        m_edges.emplace_back(m_nodes[i], m_nodes[(i + 1) % ELEMENT_POINTS]);

	m_matrixH.resize(ELEMENT_POINTS, std::vector<double>(ELEMENT_POINTS, 0.0));
    m_boundryConditionH.resize(ELEMENT_POINTS, std::vector<double>(ELEMENT_POINTS, 0.0));
    m_vectorP.resize(ELEMENT_POINTS, 0.0);
    m_matrixC.resize(ELEMENT_POINTS, std::vector<double>(ELEMENT_POINTS, 0.0));
}

void Element::CalculateJacobians()
{
    if (auto universalElement = UniversalElement::Get())
    {
        for (int i = 0; i < universalElement->GetIntegralPointsNumber(); ++i)
            m_jacobians.emplace_back(i, *m_nodes[i]);
    }
}

void Element::CalculateMatrixH(double conductivity)
{
    /*std::vector<double> weights = { 1.0, 1.0 };

    for (int i = 0; i < ELEMENT_POINTS; ++i)
    {
        for (int j = 0; j < ELEMENT_POINTS; ++j)
        {
            m_matrixH[i][j] =  (m_jacobians[0].m_dN_dx.at(i) * m_jacobians[0].m_dN_dx.at(j) + m_jacobians[0].m_dN_dy[i] * m_jacobians[0].m_dN_dy[j]) *
                m_jacobians[0].determinant * weights.at(0) * weights.at(0);

            m_matrixH[i][j] += (m_jacobians[1].m_dN_dx.at(i) * m_jacobians[1].m_dN_dx.at(j) + m_jacobians[1].m_dN_dy[i] * m_jacobians[1].m_dN_dy[j]) *
                m_jacobians[1].determinant * weights.at(1) * weights.at(0);

            m_matrixH[i][j] += (m_jacobians[2].m_dN_dx.at(i) * m_jacobians[2].m_dN_dx.at(j) + m_jacobians[2].m_dN_dy[i] * m_jacobians[2].m_dN_dy[j]) *
                m_jacobians[2].determinant * weights.at(0) * weights.at(1);

            m_matrixH[i][j] += (m_jacobians[3].m_dN_dx.at(i) * m_jacobians[3].m_dN_dx.at(j) + m_jacobians[3].m_dN_dy[i] * m_jacobians[3].m_dN_dy[j]) * 
               m_jacobians[3].determinant * weights.at(1) * weights.at(1);

            m_matrixH[i][j] *= conductivity;
        }
    }*/

    if (auto universalElement = UniversalElement::Get())
    {
        std::vector<double> weights = getGaussWeights(universalElement->GetGaussNumber());

        int n = universalElement->GetGaussNumber();
        int integralPointIndex = 0;

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

void Element::CalculateBoundryConditionsH(std::vector<Node>& nodes, double alfa, double ambientTemp)
{
    std::vector<double> gaussNodes = { -sqrt(1.0 / 3.0), sqrt(1.0 / 3.0) };
    std::vector<double> weights = { 1.0, 1.0 };

    std::vector<std::pair<double, double>> edgePoints = {
        { gaussNodes.at(0), -1.0 },
        { gaussNodes.at(1), -1.0 },
        { 1.0, gaussNodes.at(0) },
        { 1.0, gaussNodes.at(1) },
        { gaussNodes.at(1), 1.0 },
        { gaussNodes.at(0), 1.0 },
        { -1.0, gaussNodes.at(1) },
        { -1.0, gaussNodes.at(0) }
    };

    for (int k = 0; k < ELEMENT_POINTS; ++k)
    {
        int nextK = k + 1;
        if (nextK == ELEMENT_POINTS)
            nextK = 0;

        if (nodes.at(m_nodesIds[k]).bBoundaryCondition && nodes.at(m_nodesIds[nextK]).bBoundaryCondition)
        {
            std::vector<double> firstPointShapeFunctions(4, 0.0);
            firstPointShapeFunctions.at(k) += shapeFunction(k, edgePoints.at(2 * k).first, edgePoints.at(2 * k).second);
            firstPointShapeFunctions.at(nextK) += shapeFunction(nextK, edgePoints.at(2 * k).first, edgePoints.at(2 * k).second);

            std::vector<double> secondPointShapeFunctions(4, 0.0);
            secondPointShapeFunctions.at(k) += shapeFunction(k, edgePoints.at(2 * k + 1).first, edgePoints.at(2 * k + 1).second);
            secondPointShapeFunctions.at(nextK) += shapeFunction(nextK, edgePoints.at(2 * k + 1).first, edgePoints.at(2 * k + 1).second);

            double detJ = 0.5 * sqrt(
                pow(nodes.at(m_nodesIds[nextK]).first - nodes.at(m_nodesIds[k]).first, 2) +
                pow(nodes.at(m_nodesIds[nextK]).second - nodes.at(m_nodesIds[k]).second, 2));

            for (int i = 0; i < ELEMENT_POINTS; ++i)
            {
                for (int j = 0; j < ELEMENT_POINTS; ++j)
                {
                    m_boundryConditionH[i][j] += firstPointShapeFunctions.at(i) * firstPointShapeFunctions.at(j) * alfa * weights.at(0) * detJ;
                    m_boundryConditionH[i][j] += secondPointShapeFunctions.at(i) * secondPointShapeFunctions.at(j) * alfa * weights.at(1) * detJ;
                }

                m_vectorP[i] += firstPointShapeFunctions.at(i) * alfa * weights.at(0) * detJ * ambientTemp;
                m_vectorP[i] += secondPointShapeFunctions.at(i) * alfa * weights.at(1) * detJ * ambientTemp;
            }
        }
    }
}

void Element::CalculateMatrixC(double specificHeat, double density)
{
    double a = sqrt(1.0 / 3.0);
    std::vector<std::pair<double, double>> gaussPoints = {
        { -a, -a },
        { a, -a },
        { a, a },
        { -a, a }
    };

    // I don't multiply by weights here, because they are all 1.0 for 2 point Gauss quadrature,
    // but in case of changing number of points, this should be updated
    for (int i = 0; i < ELEMENT_POINTS; ++i)
    {
        for (int j = 0; j < ELEMENT_POINTS; ++j)
        {
            for (int k = 0; k < ELEMENT_POINTS; ++k)
                m_matrixC[i][j] += (shapeFunction(i, gaussPoints.at(k).first, gaussPoints.at(k).second) * shapeFunction(j, gaussPoints.at(k).first, gaussPoints.at(k).second)) * m_jacobians[k].determinant;

            m_matrixC[i][j] *= specificHeat * density;
        }
    }
}

void Element::PrintIds() const
{
    for (auto node : m_nodes)
        std::cout << std::setw(2) << node->m_id << " ";
}

void Element::PrintJacobians() const
{
    std::cout << "Element Nodes Ids: ";
    for (auto node : m_nodes)
        std::cout << std::setw(2) << node->m_id << " ";
    std::cout << "\n\n";

    for (int i = 0; i < sizeof(m_jacobians) / sizeof(Jacobian); ++i)
    {
        std::cout << i + 1 << " ";
        m_jacobians[i].Print();
    }
}