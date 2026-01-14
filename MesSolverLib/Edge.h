#pragma once
#include <vector>

class Node;

enum EEdgeAligment
{
	eBottom = 0,
	eRight = 1,
	eTop = 2,
	eLeft = 3
};

class Edge
{
public:
	Node* m_firstVertice = nullptr;
	Node* m_secondVertice = nullptr;
	bool m_bBoundaryCondition = false;
	EEdgeAligment m_alignment = eBottom;

	std::vector<std::vector<double>> m_matrixHbc{};
	std::vector<double> m_vectorP{};

	Edge() = default;

	Edge(Node* firstNode, Node* secondNode, EEdgeAligment aligment);

	void CalculateMatrixHbcAndVectorP(double alfa, double ambientTemp);
};

