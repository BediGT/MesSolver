#pragma once
#include "Element.h"
#include "Node.h"
#include <vector>
#include <string>

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

    std::vector<std::vector<double>> m_globalH{};
    std::vector<double> m_globalP{};
    std::vector<double> m_nodesTemperatures{};
    std::vector<std::vector<double>> m_globalC{};

public:
    MesSolver(const std::string& strFileName, int nGaussPoints);
    void LoadData(const std::string& strFileName);

    void CalculateSolution();
	void StartTimeSimulation();

    void PrintData();
    void PrintJacobiansForElements();
};
