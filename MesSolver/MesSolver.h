#pragma once
#include <vector>
#include <string>
#include <utility>

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

    Node(double argX, double argY)
        : x(argX), y(argY)
    {}
};

struct UniversalElement
{
    std::vector<std::vector<double>> dN_dKsi{};
    std::vector<std::vector<double>> dN_dEta{};

    UniversalElement() = default;

    UniversalElement(int nodesNumber)
    {
        dN_dKsi.resize(nodesNumber, std::vector<double>(ElementPoints, 0.0));
        dN_dEta.resize(nodesNumber, std::vector<double>(ElementPoints, 0.0));
    }
};

struct Jacobian
{
    double matrix[2][2]{};
    double matrixTransposed[2][2]{};
    double determinant{};

    double dN_dx[ElementPoints]{};
    double dN_dy[ElementPoints]{};

    void CalculateValue(const UniversalElement& universalElement, int i, int* nodesIds, std::vector<Node>& nodes);

    void Print() const;
};

struct Element
{
    int nodesIds[4]{};
    Jacobian jacobians[ElementPoints];
    double matrixH[4][4]{};
    double boundryConditionH[4][4]{};
    double vectorP[4]{};
    double matrixC[4][4]{};

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
    GlobalData m_globalData{};
    std::vector<Node> m_nodes{};
    std::vector<Element> m_elements{};
    std::vector<int> m_boundaryConditions{};

    std::vector<double> m_gaussValues{};
	std::vector<double> m_gaussWeights{};
    std::vector<std::pair<double, double>> m_gaussNodes{};
    std::vector<std::pair<double, double>> m_gaussBoundryNodes{};

    UniversalElement m_universalElement{};
    std::vector<std::vector<double>> m_globalH{};
    std::vector<double> m_globalP{};
    std::vector<double> m_nodesTemperatures{};
    std::vector<std::vector<double>> m_globalC{};

public:
    MesSolver(const std::string& strFileName, int gaussPoints);
    void LoadData(const std::string& strFileName);

    void CalculateDeritatives();

    void CalculateSolution();
	void StartTimeSimulation();

    void PrintData();
    void PrintDerivatives() const;
    void PrintJacobiansForElements();
};




