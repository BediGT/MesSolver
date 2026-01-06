#include "Edge.h"
#include "Node.h"

Edge::Edge(Node* firstNode, Node* secondNode)
	: m_firstNode(firstNode), m_secondNode(secondNode)
{
	m_bBoundaryCondition = m_firstNode->m_bBoundaryCondition && m_secondNode->m_bBoundaryCondition;
}