#include "Node.h"
#include <cmath>

Node::Node(int id, double x, double y, bool bBoundaryCondition)
	: m_id(id), m_x(x), m_y(y), m_bBoundaryCondition(bBoundaryCondition)
{
}

double Node::GetDistance(const Node& first, const Node& second)
{
	return sqrt(pow(second.m_x - first.m_x, 2.0) + pow(second.m_y - first.m_y, 2.0));
}
