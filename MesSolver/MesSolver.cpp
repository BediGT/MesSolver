#include "MesSolver.h"
#include "Helpers.h"

static double shapeFunctionDer_Ksi(int index, double eta)
{
    switch (index)
    {
        // 0 - N1, 1 - N2, 2 - N3, 3 - N4
    case 0:
    {
        return -0.25 * (1 - eta);
        break;
    }
    case 1:
    {
        return 0.25 * (1 - eta);
        break;
    }
    case 2:
    {
        return 0.25 * (1 + eta);
        break;
    }
    case 3:
    {
        return -0.25 * (1 + eta);
        break;
    }
    default:
    {
        return 0.0;
        break;
    }
    }
}

static double shapeFunctionDer_Eta(int index, double ksi)
{
    switch (index)
    {
        // 0 - N1, 1 - N2, 2 - N3, 3 - N4
    case 0:
    {
        return -0.25 * (1 - ksi);
        break;
    }
    case 1:
    {
        return -0.25 * (1 + ksi);
        break;
    }
    case 2:
    {
        return 0.25 * (1 + ksi);
        break;
    }
    case 3:
    {
        return 0.25 * (1 - ksi);
        break;
    }
    default:
    {
        return 0.0;
        break;
    }
    }
}

static double shapeFunction(int index, double ksi, double eta)
{
    switch (index)
    {
        // 0 - N1, 1 - N2, 2 - N3, 3 - N4
    case 0:
    {
        return 0.25 * (1 - ksi) * (1 - eta);
        break;
    }
    case 1:
    {
        return 0.25 * (1 + ksi) * (1 - eta);
        break;
    }
    case 2:
    {
        return 0.25 * (1 + ksi) * (1 + eta);
        break;
    }
    case 3:
    {
        return 0.25 * (1 - ksi) * (1 + eta);
        break;
    }
    default:
    {
        return 0.0;
        break;
    }
    }
}

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

Node::Node(double argX, double argY)
    : x(argX), y(argY)
{
}

void Jacobian::CalculateValue(const UniversalElement& universalElement, int i, int* nodesIds, std::vector<Node>& nodes)
{
    for (int k = 0; k < ElementPoints; ++k)
    {
        matrix[0][0] += universalElement.dN_dKsi[i][k] * nodes.at(nodesIds[k]).x;
        matrix[0][1] += universalElement.dN_dKsi[i][k] * nodes.at(nodesIds[k]).y;
        matrix[1][0] += universalElement.dN_dEta[i][k] * nodes.at(nodesIds[k]).x;
        matrix[1][1] += universalElement.dN_dEta[i][k] * nodes.at(nodesIds[k]).y;
    }

    determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

    matrixTransposed[1][1] = matrix[0][0] / determinant;
    matrixTransposed[1][0] = -matrix[1][0] / determinant;
    matrixTransposed[0][0] = matrix[1][1] / determinant;
    matrixTransposed[0][1] = -matrix[0][1] / determinant;

    for (int j = 0; j < ElementPoints; ++j)
    {
        dN_dx[j] = (matrixTransposed[0][0] * universalElement.dN_dKsi[i][j] + matrixTransposed[0][1] * universalElement.dN_dEta[i][j]);
        dN_dy[j] = (matrixTransposed[1][0] * universalElement.dN_dKsi[i][j] + matrixTransposed[1][1] * universalElement.dN_dEta[i][j]);
    }
}

void Jacobian::Print() const
{
    std::cout << "Jacobian Matrix:\n";
    std::cout << "[" << matrix[0][0] << ", " << matrix[0][1] << "]\n";
    std::cout << "[" << matrix[1][0] << ", " << matrix[1][1] << "]\n";
    std::cout << "Determinant: " << determinant << "\n\n";

    std::cout << "dN/dx\n";
    std::cout << dN_dx[0] << ", " << dN_dx[1] << ", " << dN_dx[2] << ", " << dN_dx[3] << "\n";
    std::cout << "dN/dy\n";
    std::cout << dN_dy[0] << ", " << dN_dy[1] << ", " << dN_dy[2] << ", " << dN_dy[3] << "\n\n\n";
}

Element::Element(int arg1, int arg2, int arg3, int arg4)
{
    nodesIds[0] = arg1;
    nodesIds[1] = arg2;
    nodesIds[2] = arg3;
    nodesIds[3] = arg4;
}

void Element::CalculateJacobians(std::vector<Node>& nodes, UniversalElement& universalElement)
{
    for (int i = 0; i < ElementPoints; ++i)
        jacobians[i].CalculateValue(universalElement, i, nodesIds, nodes);
}

void Element::CalculateMatrixH(double conductivity)
{
    std::vector<double> weights = { 1.0, 1.0 };

    for (int i = 0; i < ElementPoints; ++i)
    {
        for (int j = 0; j < ElementPoints; ++j)
        {
            matrixH[i][j] = (jacobians[0].dN_dx[i] * jacobians[0].dN_dx[j] + jacobians[0].dN_dy[i] * jacobians[0].dN_dy[j]) * jacobians[0].determinant * weights.at(0) * weights.at(0);
            matrixH[i][j] += (jacobians[1].dN_dx[i] * jacobians[1].dN_dx[j] + jacobians[1].dN_dy[i] * jacobians[1].dN_dy[j]) * jacobians[1].determinant * weights.at(1) * weights.at(0);
            matrixH[i][j] += (jacobians[2].dN_dx[i] * jacobians[2].dN_dx[j] + jacobians[2].dN_dy[i] * jacobians[2].dN_dy[j]) * jacobians[2].determinant * weights.at(0) * weights.at(1);
            matrixH[i][j] += (jacobians[3].dN_dx[i] * jacobians[3].dN_dx[j] + jacobians[3].dN_dy[i] * jacobians[3].dN_dy[j]) * jacobians[3].determinant * weights.at(1) * weights.at(1);
            matrixH[i][j] *= conductivity;
        }
    }
}

void Element::CalculateBoundryConditionsH(std::vector<Node>& nodes, double alfa, double ambientTemp)
{
    std::vector<double> gaussNodes = { -sqrt(1.0 / 3.0), sqrt(1.0 / 3.0) };
    std::vector<double> weights = { 1.0, 1.0 };

    std::vector<std::pair<double, double>> edgePoints = {
        { gaussNodes.at(0), -1.0 },
        { gaussNodes.at(1), -1.0 },
        { 1.0, gaussNodes.at(0) },
        { 1.0, gaussNodes.at(1) },
        { gaussNodes.at(1), 1.0 },
        { gaussNodes.at(0), 1.0 },
        { -1.0, gaussNodes.at(1) },
        { -1.0, gaussNodes.at(0) }
    };

    for (int k = 0; k < ElementPoints; ++k)
    {
        int nextK = k + 1;
        if (nextK == ElementPoints)
            nextK = 0;

        if (nodes.at(nodesIds[k]).bBoundaryCondition && nodes.at(nodesIds[nextK]).bBoundaryCondition)
        {
            std::vector<double> firstPointShapeFunctions(4, 0.0);
            firstPointShapeFunctions.at(k) += shapeFunction(k, edgePoints.at(2 * k).first, edgePoints.at(2 * k).second);
            firstPointShapeFunctions.at(nextK) += shapeFunction(nextK, edgePoints.at(2 * k).first, edgePoints.at(2 * k).second);

            std::vector<double> secondPointShapeFunctions(4, 0.0);
            secondPointShapeFunctions.at(k) += shapeFunction(k, edgePoints.at(2 * k + 1).first, edgePoints.at(2 * k + 1).second);
            secondPointShapeFunctions.at(nextK) += shapeFunction(nextK, edgePoints.at(2 * k + 1).first, edgePoints.at(2 * k + 1).second);

            double detJ = 0.5 * sqrt(
                pow(nodes.at(nodesIds[nextK]).x - nodes.at(nodesIds[k]).x, 2) +
                pow(nodes.at(nodesIds[nextK]).y - nodes.at(nodesIds[k]).y, 2));

            for (int i = 0; i < ElementPoints; ++i)
            {
                for (int j = 0; j < ElementPoints; ++j)
                {
                    boundryConditionH[i][j] += firstPointShapeFunctions.at(i) * firstPointShapeFunctions.at(j) * alfa * weights.at(0) * detJ;
                    boundryConditionH[i][j] += secondPointShapeFunctions.at(i) * secondPointShapeFunctions.at(j) * alfa * weights.at(1) * detJ;
                }

                vectorP[i] += firstPointShapeFunctions.at(i) * alfa * weights.at(0) * detJ * ambientTemp;
                vectorP[i] += secondPointShapeFunctions.at(i) * alfa * weights.at(1) * detJ * ambientTemp;
            }
        }
    }
}

void Element::CalculateMatrixC(double specificHeat, double density)
{
    double a = sqrt(1.0 / 3.0);
    std::vector<std::pair<double, double>> gaussPoints = {
        { -a, -a },
        { a, -a },
        { a, a },
        { -a, a }
    };

	// I don't multiply by weights here, because they are all 1.0 for 2 point Gauss quadrature,
    // but in case of changing number of points, this should be updated
    for (int i = 0; i < ElementPoints; ++i)
    {
        for (int j = 0; j < ElementPoints; ++j)
        {
            for (int k = 0; k < ElementPoints; ++k)
                matrixC[i][j] += (shapeFunction(i, gaussPoints.at(k).first, gaussPoints.at(k).second) * shapeFunction(j, gaussPoints.at(k).first, gaussPoints.at(k).second)) * jacobians[k].determinant;
                
            matrixC[i][j] *= specificHeat * density;
        }
    }
}

void Element::PrintJacobians() const
{
    std::cout << "Element Nodes Ids: " << nodesIds[0] << " " << nodesIds[1] << " " << nodesIds[2] << " " << nodesIds[3] << "\n\n";
    for (int i = 0; i < sizeof(jacobians) / sizeof(Jacobian); ++i)
    {
        std::cout << i + 1 << " ";
        jacobians[i].Print();
    }
}

MesSolver::MesSolver()
{

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
            m_nodes.push_back(Node(x, y));
            break;
        }
        case eElement:
        {
            int id, n1, n2, n3, n4;
            char comma;
            iss >> id >> comma >> n1 >> comma >> n2 >> comma >> n3 >> comma >> n4;
            m_elements.push_back(Element(n1 - 1, n2 - 1, n3 - 1, n4 - 1)); // Indexing from 0
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

    // TO INIT
    m_globalH.resize(m_nodes.size());
    for (auto& row : m_globalH)
        row.resize(m_nodes.size(), 0.0);

    for (int index : m_boundaryConditions)
        m_nodes.at(index).bBoundaryCondition = true;

    m_globalP.resize(m_nodes.size(), 0.0);

    m_nodesTemperatures.resize(m_nodes.size(), 0.0);

    m_globalC.resize(m_nodes.size());
    for (auto& row : m_globalC)
        row.resize(m_nodes.size(), 0.0);
}

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
            ", x: " << std::setprecision(precision) << std::setw(precision + 3) << m_nodes.at(i).x <<
            ", y: " << std::setprecision(precision) << std::setw(precision + 3) << m_nodes.at(i).y << "\n";
    }

    std::cout << "\nElements:\n";
    for (size_t i = 0; i < m_elements.size(); ++i)
    {
        const auto& el = m_elements[i];
        std::cout << std::left << "\tID: " << std::setw(2) << i << ", Nodes: ";
        for (int nid : el.nodesIds)
        {
            std::cout << std::setw(2) << nid << " ";
        }
        std::cout << "\n";
    }


    std::cout << "\nBoundary conditions:\n\t";
    for (int b : m_boundaryConditions)
    {
        std::cout << b << " ";
    }
    std::cout << "\n";
}

void MesSolver::CalculateDeritatives()
{
    const double a = 1.0 / std::sqrt(3.0);
    std::vector<std::pair<double, double>> gaussPoints =
    {
        { -a, -a },
        {  +a, -a },
        {  +a,  +a },
        { -a,  +a }
    };

    for (int i = 0; i < ElementPoints; ++i)
    {
        for (int j = 0; j < ElementPoints; ++j)
        {
            universalElement.dN_dKsi[i][j] = shapeFunctionDer_Ksi(j, gaussPoints.at(i).second);
			universalElement.dN_dEta[i][j] = shapeFunctionDer_Eta(j, gaussPoints.at(i).first);
        }
	}
}

void MesSolver::PrintDerivatives() const
{
    std::vector<double> gaussNodes = { -sqrt(1.0 / 3.0), sqrt(1.0 / 3.0) };

    std::cout << "Shape Functions Derivatives at Integration Points:\n";
    std::cout << "dN/dKsi" << '\n';
    for (int i = 0; i < ElementPoints; ++i)
    {
        std::cout << universalElement.dN_dKsi[i][0] << " ";
        std::cout << universalElement.dN_dKsi[i][1] << " ";
        std::cout << universalElement.dN_dKsi[i][2] << " ";
        std::cout << universalElement.dN_dKsi[i][3] << '\n';
    }

    std::cout << "dN/dEta" << '\n';
    for (int i = 0; i < ElementPoints; ++i)
    {
        std::cout << universalElement.dN_dEta[i][0] << " ";
        std::cout << universalElement.dN_dEta[i][1] << " ";
        std::cout << universalElement.dN_dEta[i][2] << " ";
        std::cout << universalElement.dN_dEta[i][3] << '\n';
    }
}

void MesSolver::PrintJacobiansForElements()
{
    for (const auto& element : m_elements)
    {
        element.PrintJacobians();
        std::cout << "\n";
    }
}

void MesSolver::CalculateSolution()
{
    CalculateDeritatives();
    //PrintDerivatives();

    for (auto& element : m_elements)
    {
        element.CalculateJacobians(m_nodes, universalElement);
        element.CalculateMatrixH(m_globalData.Conductivity);
        element.CalculateBoundryConditionsH(m_nodes, m_globalData.Alfa, m_globalData.AmbientTemperature);
		element.CalculateMatrixC(m_globalData.SpecificHeat, m_globalData.Density);
    }

    //PrintJacobiansForElements();

    // Calculate Globals
    for (auto& element : m_elements)
    {
        for (int i = 0; i < ElementPoints; ++i)
        {
            for (int j = 0; j < ElementPoints; ++j)
            {
                m_globalH.at(element.nodesIds[i]).at(element.nodesIds[j]) += element.matrixH[i][j];
                m_globalH.at(element.nodesIds[i]).at(element.nodesIds[j]) += element.boundryConditionH[i][j];
                m_globalC.at(element.nodesIds[i]).at(element.nodesIds[j]) += element.matrixC[i][j];
            }

            m_globalP.at(element.nodesIds[i]) += element.vectorP[i];
        }
    }

    m_nodesTemperatures = solveLinearSystem(m_globalH, m_globalP);

    /*std::cout << "\nGlobal Matrix H:\n";
    printMatrix(m_globalH);

    std::cout << "\nGlobal Vector P:\n";
    printVector(m_globalP);

    std::cout << "\nNodes Temperatures:\n";
    printVector(m_nodesTemperatures);

    std::cout << "\nGlobal Matrix C:\n";
    printMatrix(m_globalC);*/
}

void MesSolver::StartTimeSimulation()
{
    CalculateSolution();

	std::vector<double> initialTemps(m_nodes.size(), m_globalData.InitialTemp);
    for (int i = 0; i < m_globalData.SimulationTime; i += m_globalData.SimulationStepTime)
    {
		double stepTime = static_cast<double>(m_globalData.SimulationStepTime);
        std::vector<std::vector<double>> scaledC = MatrixUtils::Scale(m_globalC, 1.0 / stepTime);
            
        initialTemps = solveLinearSystem(MatrixUtils::Add(m_globalH, scaledC), MatrixUtils::Add(m_globalP, MatrixUtils::Multiply(scaledC, initialTemps)));

		std::cout << "\nTime: " << i + m_globalData.SimulationStepTime << " s\n";
		//printVector(initialTemps);
        std::cout << "\nMin: " << *std::min_element(initialTemps.begin(), initialTemps.end()) << "\n";
        std::cout << "\nMax: " << *std::max_element(initialTemps.begin(), initialTemps.end()) << "\n";
        std::cout << "\n\n";
    }
}

//static void GetGaussValues(int n, std::vector<double>& nodes, std::vector<double>& weights)
//{
//    nodes.clear();
//    weights.clear();
//
//    switch (n)
//    {
//    case 1:
//        nodes = { 0.0 };
//        weights = { 2.0 };
//        break;
//    case 2:
//        nodes = { -sqrt(1.0 / 3.0), sqrt(1.0 / 3.0) };
//        weights = { 1.0, 1.0 };
//        break;
//    case 3:
//        nodes = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };
//        weights = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
//        break;
//    case 4:
//        nodes = { -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)) };
//        weights = { (18.0 - sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0 };
//        break;
//    case 5:
//        nodes = { -1.0 / 3.0 * sqrt(5.0 + 2 * sqrt(10.0 / 7.0)), -1.0 / 3.0 * sqrt(5.0 - 2 * sqrt(10.0 / 7.0)), 0.0, 1.0 / 3.0 * sqrt(5.0 - 2 * sqrt(10.0 / 7.0)), 1.0 / 3.0 * sqrt(5.0 + 2 * sqrt(10.0 / 7.0)) };
//        weights = { (322.0 - 13.0 * sqrt(70.0)) / 900.0, (322.0 + 13.0 * sqrt(70.0)) / 900.0, 128.0 / 225.0, (322.0 + 13.0 * sqrt(70.0)) / 900.0, (322.0 - 13.0 * sqrt(70.0)) / 900.0 };
//        break;
//    default:
//        return;
//    }
//}
