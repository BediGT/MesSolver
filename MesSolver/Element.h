#pragma once
#include "Constants.h"
#include "Jacobian.h"
#include "Edge.h"
#include <vector>

class Node;

class Element
{
public:
    std::vector<Node*> m_nodes{};
    std::vector<Edge> m_edges{ ELEMENT_POINTS };
    std::vector<Jacobian> m_jacobians{};

    std::vector<std::vector<double>> m_matrixH{};
    std::vector<std::vector<double>> m_boundryConditionH{};
    std::vector<double> m_vectorP{};
    std::vector<std::vector<double>> m_matrixC{};

    Element(Node* arg1, Node* arg2, Node* arg3, Node* arg4);

    void CalculateJacobians();
    void CalculateMatrixH(double conductivity);
    void CalculateBoundryConditionsH(std::vector<Node>& nodes, double alfa, double ambientTemp);
    void CalculateMatrixC(double specificHeat, double density);

    void PrintIds() const;
    void PrintJacobians() const;
};
