#include "UniversalElement.h"

#include "Node.h"
#include "Constants.h"

#include <iostream>

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

UniversalElement* UniversalElement::m_instance = nullptr;

UniversalElement::UniversalElement(int n)
    : m_nGauss(n)
{
    m_dN_dKsi.resize(m_nGauss * m_nGauss, std::vector<double>(ELEMENT_POINTS, 0.0));
    m_dN_dEta.resize(m_nGauss * m_nGauss, std::vector<double>(ELEMENT_POINTS, 0.0));

    std::vector<double> gaussValues = getGaussValues(m_nGauss);

    int integralPointIndex = 0;
    for (int iEta = 0; iEta < m_nGauss; ++iEta)
    {
        for (int iKsi = 0; iKsi < m_nGauss; ++iKsi)
        {
            for (int i = 0; i < ELEMENT_POINTS; ++i)
            {
                m_dN_dKsi.at(integralPointIndex).at(i) = shapeFunctionDer_Ksi(i, gaussValues.at(iEta));
                m_dN_dEta.at(integralPointIndex).at(i) = shapeFunctionDer_Eta(i, gaussValues.at(iKsi));
            }
            integralPointIndex++;
        }
    }
}

void UniversalElement::Init(int n)
{
    if (m_instance)
		delete m_instance;

    m_instance = new UniversalElement(n);
}

UniversalElement* const UniversalElement::Get()
{
	return m_instance;
}

int UniversalElement::GetGaussNumber() const
{
    return m_nGauss;
}

int UniversalElement::GetIntegralPointsNumber() const
{
    return m_nGauss * m_nGauss;
}

void UniversalElement::Print() const
{
    std::cout << "Shape Functions Derivatives at Integration Points:\n";
    std::cout << "dN/dKsi" << '\n';
    for (int i = 0; i < m_nGauss * m_nGauss; ++i)
    {
        for (int k = 0; k < ELEMENT_POINTS; ++k)
            std::cout << m_dN_dKsi[i][k] << " ";
        std::cout << "\n";
    }


    std::cout << "dN/dEta" << '\n';
    for (int i = 0; i < m_nGauss * m_nGauss; ++i)
    {
        for (int k = 0; k < ELEMENT_POINTS; ++k)
            std::cout << m_dN_dEta[i][k] << " ";
        std::cout << "\n";
    }
}
