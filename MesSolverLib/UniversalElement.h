#pragma once
#include <vector>

class Node;

class UniversalElement
{
    int m_nGauss = 0;

public:
    std::vector<std::vector<double>> m_dN_dKsi{};
    std::vector<std::vector<double>> m_dN_dEta{};

private:
    static UniversalElement* m_instance;

    UniversalElement(int n);
	~UniversalElement() = default;

    UniversalElement(const UniversalElement&) = delete;
    UniversalElement& operator=(const UniversalElement&) = delete;
    UniversalElement(UniversalElement&&) = delete;
    UniversalElement& operator=(UniversalElement&&) = delete;

public:
    static void Init(int n);
    static UniversalElement* const Get();

    int GetGaussNumber() const;
    int GetIntegralPointsNumber() const;

    void Print() const;
};