#include "MesSolver.h"

#include "Helpers.h"

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <iostream>

#pragma region static Funtions
static void printMatrix(const double matrix[4][4])
{
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j) {
            std::cout << std::setw(10) << matrix[i][j] << " ";
        }
        std::cout << "\n";
    }
}

static void printMatrix(const std::vector<std::vector<double>> matrix)
{
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            std::cout << std::setw(12) << std::fixed << std::setprecision(6)
                << matrix[i][j] << " ";
        }
        std::cout << "\n";
    }
}

static void printVector(const std::vector<double> vector)
{
    for (size_t i = 0; i < vector.size(); ++i)
        std::cout << std::setw(12) << std::fixed << std::setprecision(6) << vector.at(i) << " ";
    std::cout << "\n";
}
#pragma endregion

MesSolver::MesSolver(const std::string& strFileName, int nGaussPoints)
    : m_nGaussPoints(nGaussPoints)
{
	LoadData(strFileName);

    m_globalH.resize(m_nodes.size(), std::vector<double>(m_nodes.size(), 0.0));

    for (auto& index : m_boundaryConditions)
        m_nodes.at(index).m_bBoundaryCondition = true;

    m_globalP.resize(m_nodes.size(), 0.0);

    m_nodesTemperatures.resize(m_nodes.size(), 0.0);

    m_globalC.resize(m_nodes.size(), std::vector<double>(m_nodes.size(), 0.0));
}

void MesSolver::LoadData(const std::string& strFileName)
{
    std::ifstream file(strFileName);
    if (!file.is_open())
    {
        std::cerr << "Cannot open file!" << std::endl;
        return;
    }

    std::string line;
    Section section = eHeader;

    while (std::getline(file, line))
    {
        if (line.empty()) continue;
        if (line.find("*Node") != std::string::npos)
        {
            section = eNode;
            continue;
        }
        if (line.find("*Element") != std::string::npos)
        {
            section = eElement;
            continue;
        }
        if (line.find("*BC") != std::string::npos)
        {
            section = eBoundaryCondition;
            continue;
        }

        std::istringstream iss(line);
        switch (section)
        {
        case eHeader:
        {
            std::string key;
            iss >> key;
            if (key == "SimulationTime") iss >> m_globalData.SimulationTime;
            else if (key == "SimulationStepTime") iss >> m_globalData.SimulationStepTime;
            else if (key == "Conductivity") iss >> m_globalData.Conductivity;
            else if (key == "Alfa") iss >> m_globalData.Alfa;
            else if (key == "Tot") iss >> m_globalData.AmbientTemperature;
            else if (key == "InitialTemp") iss >> m_globalData.InitialTemp;
            else if (key == "Density") iss >> m_globalData.Density;
            else if (key == "SpecificHeat") iss >> m_globalData.SpecificHeat;
            else if (key == "Nodes")
            {
                std::string dummy;
                iss >> dummy >> m_globalData.NodesNumber;
            }
            else if (key == "Elements")
            {
                std::string dummy;
                iss >> dummy >> m_globalData.ElementsNumber;
            }

            break;
        }
        case eNode:
        {
            int id;
            double x, y;
            char comma;
            iss >> id >> comma >> x >> comma >> y;
            m_nodes.push_back(Node(id, x, y));
            break;
        }
        case eElement:
        {
            int id, n1, n2, n3, n4;
            char comma;
            iss >> id >> comma >> n1 >> comma >> n2 >> comma >> n3 >> comma >> n4;
            m_elements.push_back(Element(&m_nodes.at(n1 - 1), &m_nodes.at(n2 - 1), &m_nodes.at(n3 - 1), &m_nodes.at(n4 - 1)));
            break;
        }
        case eBoundaryCondition:
        {
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream bcStream(line);
            int val;
            m_boundaryConditions.clear();
            while (bcStream >> val)
            {
                m_boundaryConditions.push_back(val - 1);
            }
            break;
        }
        }
    }

    file.close();

    PrintData();
}

void MesSolver::CalculateSolution()
{
	UniversalElement::Init(m_nGaussPoints);

    for (auto& element : m_elements)
    {
        element.CalculateJacobians();
        element.CalculateMatrixH(m_globalData.Conductivity);
        element.CalculateBoundryConditionsH(m_nodes, m_globalData.Alfa, m_globalData.AmbientTemperature);
		element.CalculateMatrixC(m_globalData.SpecificHeat, m_globalData.Density);
    }

    //PrintJacobiansForElements();

    for (auto& element : m_elements)
    {
        for (int i = 0; i < ELEMENT_POINTS; ++i)
        {
            for (int j = 0; j < ELEMENT_POINTS; ++j)
            {
                m_globalH.at(element.m_nodes[i]->m_id).at(element.m_nodes[j]->m_id) += element.m_matrixH[i][j];
                m_globalH.at(element.m_nodes[i]->m_id).at(element.m_nodes[j]->m_id) += element.m_matrixH[i][j];
                m_globalH.at(element.m_nodes[i]->m_id).at(element.m_nodes[j]->m_id) += element.m_boundryConditionH[i][j];
                m_globalC.at(element.m_nodes[i]->m_id).at(element.m_nodes[j]->m_id) += element.m_matrixC[i][j];
            }

            m_globalP.at(element.m_nodes[i]->m_id) += element.m_vectorP[i];
        }
    }

    m_nodesTemperatures = SolveLinearSystem(m_globalH, m_globalP);

    std::cout << "\nGlobal Matrix H:\n";
    printMatrix(m_globalH);

    std::cout << "\nGlobal Vector P:\n";
    printVector(m_globalP);

    std::cout << "\nNodes Temperatures:\n";
    printVector(m_nodesTemperatures);

    std::cout << "\nGlobal Matrix C:\n";
    printMatrix(m_globalC);
}

void MesSolver::StartTimeSimulation()
{
    CalculateSolution();

	std::vector<double> initialTemps(m_nodes.size(), m_globalData.InitialTemp);
    for (int i = 0; i < m_globalData.SimulationTime; i += m_globalData.SimulationStepTime)
    {
		double stepTime = static_cast<double>(m_globalData.SimulationStepTime);
        std::vector<std::vector<double>> scaledC = MatrixUtils::Scale(m_globalC, 1.0 / stepTime);
            
        initialTemps = SolveLinearSystem(MatrixUtils::Add(m_globalH, scaledC), MatrixUtils::Add(m_globalP, MatrixUtils::Multiply(scaledC, initialTemps)));

		std::cout << "\nTime: " << i + m_globalData.SimulationStepTime << " s\n";
		//printVector(initialTemps);
        std::cout << "\nMin: " << *std::min_element(initialTemps.begin(), initialTemps.end()) << "\n";
        std::cout << "\nMax: " << *std::max_element(initialTemps.begin(), initialTemps.end()) << "\n";
        std::cout << "\n\n";
    }
}

#pragma region PrintFunctions
void MesSolver::PrintData()
{
    const int precision = 6;

    std::cout.precision(precision);

    std::cout << "Global Data:\n";
    std::cout << "\tSimulation Time:      " << m_globalData.SimulationTime << " s" << "\n";
    std::cout << "\tSimulation Step Time: " << m_globalData.SimulationStepTime << " s" << "\n";
    std::cout << "\tConductivity:         " << m_globalData.Conductivity << "\n";
    std::cout << "\tAlfa:                 " << m_globalData.Alfa << "\n";
    std::cout << "\tAmbient Temperature:  " << m_globalData.AmbientTemperature << "\n";
    std::cout << "\tInitial Temperature:  " << m_globalData.InitialTemp << "\n";
    std::cout << "\tDensity:              " << m_globalData.Density << "\n";
    std::cout << "\tSpecificHeat:         " << m_globalData.SpecificHeat << "\n";
    std::cout << "\tNodes:                " << m_globalData.NodesNumber << '\n';
    std::cout << "\tElements:             " << m_globalData.ElementsNumber << "\n";

    std::cout << "\nNodes:\n";
    for (size_t i = 0; i < m_nodes.size(); ++i)
    {
        std::cout << std::left << "\tID: " << std::setw(2) << i << std::fixed << std::right <<
            ", x: " << std::setprecision(precision) << std::setw(precision + 3) << m_nodes.at(i).m_x <<
            ", y: " << std::setprecision(precision) << std::setw(precision + 3) << m_nodes.at(i).m_y << "\n";
    }

    std::cout << "\nElements:\n";
    for (size_t i = 0; i < m_elements.size(); ++i)
    {
        const auto& el = m_elements.at(i);
        std::cout << std::left << "\tID: " << std::setw(2) << i << ", Nodes: ";
		el.PrintIds();
        std::cout << "\n";
    }

    std::cout << "\nBoundary conditions:\n\t";
    for (int b : m_boundaryConditions)
    {
        std::cout << b << " ";
    }
    std::cout << "\n";
}

void MesSolver::PrintJacobiansForElements()
{
    for (const auto& element : m_elements)
    {
        element.PrintJacobians();
        std::cout << "\n";
    }
}
#pragma endregion
