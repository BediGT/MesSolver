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
	Node* m_firstNode = nullptr;
	Node* m_secondNode = nullptr;
	bool m_bBoundaryCondition = false;
	EEdgeAligment m_alignment = eBottom;

	Edge() = default;

	Edge(Node* firstNode, Node* secondNode, EEdgeAligment aligment);

	std::vector<std::vector<double>> GetPartOfHbc() const;
};

