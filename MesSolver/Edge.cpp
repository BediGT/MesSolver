#include "Edge.h"
#include "Node.h"
#include "ShapeFunctions.h"
#include "Constants.h"
#include "Helpers.h"
#include "UniversalElement.h"

Edge::Edge(Node* firstVertice, Node* secondVertice, EEdgeAligment aligment)
	: m_firstVertice(firstVertice), m_secondVertice(secondVertice), m_alignment(aligment)
{
	m_bBoundaryCondition = m_firstVertice->m_bBoundaryCondition && m_secondVertice->m_bBoundaryCondition;

	m_matrixHbc.resize(ELEMENT_POINTS, std::vector<double>(ELEMENT_POINTS, 0.0));
	m_vectorP.resize(ELEMENT_POINTS, 0.0);
}

void Edge::CalculateMatrixHbcAndVectorP(double alfa, double ambientTemp)
{
	if (auto universalElement = UniversalElement::Get())
	{
		int gaussNumber = universalElement->GetGaussNumber();
		std::vector<double> values = getGaussValues(gaussNumber);
		std::vector<double> weights = getGaussWeights(gaussNumber);
		double detJ = 0.5 * Node::GetDistance(*m_firstVertice, *m_secondVertice);

		for (int i = 0; i < gaussNumber; ++i)
		{
			std::vector<double> shapeFunctionValues(ELEMENT_POINTS, 0.0);

			switch (m_alignment)
			{
				case eBottom:
				{
					shapeFunctionValues.at(0) = shapeFunction(0, values.at(i), -1.0);
					shapeFunctionValues.at(1) = shapeFunction(1, values.at(i), -1.0);
					break;
				}
				case eRight:
				{
					shapeFunctionValues.at(1) = shapeFunction(1, 1.0, values.at(i));
					shapeFunctionValues.at(2) = shapeFunction(2, 1.0, values.at(i));
					break;
				}
				case eTop:
				{
					shapeFunctionValues.at(2) = shapeFunction(2, values.at(i), 1.0);
					shapeFunctionValues.at(3) = shapeFunction(3, values.at(i), 1.0);
					break;
				}
				case eLeft:
				{
					shapeFunctionValues.at(3) = shapeFunction(3, -1.0, values.at(i));
					shapeFunctionValues.at(0) = shapeFunction(0, -1.0, values.at(i));
					break;
				}
				default:
					break;
			}

			std::vector<std::vector<double>> shapeFunctionsMatrix = MatrixUtils::MultiplyVectors(shapeFunctionValues, shapeFunctionValues);
			MatrixUtils::Scale(shapeFunctionsMatrix, weights.at(i) * detJ * alfa);

			MatrixUtils::AddTo(m_matrixHbc, shapeFunctionsMatrix);

			MatrixUtils::Scale(shapeFunctionValues, weights.at(i) * detJ * alfa * ambientTemp);
			MatrixUtils::AddTo(m_vectorP, shapeFunctionValues);
		}
	}
}
