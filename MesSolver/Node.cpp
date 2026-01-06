#include "Node.h"

Node::Node(int id, double x, double y, bool bBoundaryCondition)
	: m_id(id), m_x(x), m_y(y), m_bBoundaryCondition(bBoundaryCondition)
{
}
