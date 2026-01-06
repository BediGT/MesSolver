#pragma once

#include <vector>

std::vector<double> SolveLinearSystem(std::vector<std::vector<double>> A, std::vector<double> b);

struct MatrixUtils
{
    static std::vector<std::vector<double>> Add(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B);

    static void AddTo(std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B);

    static std::vector<double> Add(const std::vector<double>& a, const std::vector<double>& b);

    static void Scale(std::vector<std::vector<double>>& A, double scalar);

    static std::vector<double> Multiply(const std::vector<std::vector<double>>& A, const std::vector<double>& x);

    static std::vector<double> Subtract(const std::vector<double>& a, const std::vector<double>& b);

    static std::vector<std::vector<double>> MultiplyVectors(const std::vector<double>& a, const std::vector<double>& b);
};
