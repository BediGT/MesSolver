#pragma once
class Node;

class Edge
{
public:
	Node* m_firstNode = nullptr;
	Node* m_secondNode = nullptr;
	bool m_bBoundaryCondition = false;

	Edge() = default;

	Edge(Node* firstNode, Node* secondNode);
};

