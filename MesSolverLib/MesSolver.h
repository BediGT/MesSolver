#pragma once
#include "Element.h"
#include "Node.h"
#include <vector>
#include <string>

struct GlobalData
{
    double SimulationTime = 0;
    double SimulationStepTime = 0;
    double Conductivity = 0;
    double Alfa = 0;
    double AmbientTemperature = 0;
    double InitialTemp = 0;
    double Density = 0;
    double SpecificHeat = 0;
    int NodesNumber = 0;
    int ElementsNumber = 0;
};

enum ESection {
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

    std::vector<std::vector<double>> m_globalH{};
    std::vector<double> m_globalP{};
    std::vector<std::vector<double>> m_globalC{};

    bool m_bShouldPrint = false;

public:
    MesSolver(int nGaussPoints, bool bShouldPrint = false);
    bool LoadData(const std::string& strFileName);

    void SetShouldPrint(bool bShouldPrint);
	void SetGaussNumber(int nGaussPoints);

    const std::vector<std::vector<double>>& GetMatrixH();
    const std::vector<double>& GetVectorP();
    const std::vector<std::vector<double>>& GetMatrixC();

    void CalculateTimeIndependentVariables();
    std::vector<std::vector<double>> SimulateWithTime();

    void PrintData();
    void PrintJacobiansForElements();
};
