#pragma once
class Node
{
public:
	int m_id = -1;
	double m_x = 0.0;
	double m_y = 0.0;
	bool m_bBoundaryCondition = false;

	Node(int id, double x, double y, bool bBoundaryCondition = false);

	static double GetDistance(const Node& first, const Node& second);
};

