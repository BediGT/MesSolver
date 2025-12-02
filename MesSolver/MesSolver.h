#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include <functional>

const int ElementPoints = 4;

struct GlobalData
{
    int SimulationTime = 0;
    int SimulationStepTime = 0;
    int Conductivity = 0;
    int Alfa = 0;
    int AmbientTemperature = 0;
    int InitialTemp = 0;
    int Density = 0;
    int SpecificHeat = 0;
    int NodesNumber = 0;
    int ElementsNumber = 0;
};

struct Node
{
    double x = 0.0;
    double y = 0.0;
    bool bBoundaryCondition = false;

    Node(double argX, double argY);
};

struct UniversalElement
{
    double dN_dKsi[ElementPoints][4] = { 0.0 };
    double dN_dEta[ElementPoints][4] = { 0.0 };
};

struct Jacobian
{
    double matrix[2][2] = { 0.0 };
    double matrixTransposed[2][2] = { 0.0 };
    double determinant = 0.0;

    double dN_dx[4] = { 0.0 };
    double dN_dy[4] = { 0.0 };

    void CalculateValue(const UniversalElement& universalElement, int i, int* nodesIds, std::vector<Node>& nodes);

    void Print() const;
};

struct Element
{
    int nodesIds[4] = { 0 };
    Jacobian jacobians[ElementPoints];
    double matrixH[4][4] = { 0 };
    double boundryConditionH[4][4] = { 0 };
    double vectorP[4] = { 0 };
    double matrixC[4][4] = { 0 };

    Element(int arg1, int arg2, int arg3, int arg4);

    void CalculateJacobians(std::vector<Node>& nodes, UniversalElement& universalElement);
    void CalculateMatrixH(double conductivity);
    void CalculateBoundryConditionsH(std::vector<Node>& nodes, double alfa, double ambientTemp);
    void CalculateMatrixC(double specificHeat, double density);

    void PrintJacobians() const;
};

enum Section {
    eHeader,
    eNode,
    eElement,
    eBoundaryCondition
};

class MesSolver
{
    GlobalData m_globalData;
    std::vector<Node> m_nodes;
    std::vector<Element> m_elements;
    std::vector<int> m_boundaryConditions = { 0 };
    UniversalElement universalElement;
    std::vector<std::vector<double>> m_globalH = {}; // nodes.size() x nodes.size()
    std::vector<double> m_globalP = { 0.0 }; // nodes.size()
    std::vector<double> m_nodesTemperatures = { 0.0 }; // nodes.size()
    std::vector<std::vector<double>> m_globalC = {}; // nodes.size() x nodes.size()

public:
    MesSolver();

    void LoadData(const std::string& strFileName); // Is also initializer
    void PrintData();

    void CalculateDeritatives();
    void PrintDerivatives() const;

    void PrintJacobiansForElements();

    void CalculateSolution();

	void StartTimeSimulation();
};




