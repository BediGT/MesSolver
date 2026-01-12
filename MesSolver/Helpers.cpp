#include "Helpers.h"
#include <stdexcept>
#include <cmath>

std::vector<double> SolveLinearSystem(std::vector<std::vector<double>> A, std::vector<double> b)
{
    int n = (int)A.size();
    if (n == 0 || (int)A[0].size() != n || (int)b.size() != n)
        throw std::runtime_error("Matrix/vector size mismatch");

    for (int i = 0; i < n; ++i)
    {
        int maxRow = i;
        for (int k = i + 1; k < n; ++k)
            if (std::fabs(A[k][i]) > std::fabs(A[maxRow][i])) maxRow = k;

        std::swap(A[i], A[maxRow]);
        std::swap(b[i], b[maxRow]);

        double piv = A[i][i];
        if (std::fabs(piv) < 1e-14)
            throw std::runtime_error("Singular or ill-conditioned matrix");

        for (int k = i + 1; k < n; ++k) {
            double factor = A[k][i] / piv;
            for (int j = i; j < n; ++j) A[k][j] -= factor * A[i][j];
            b[k] -= factor * b[i];
        }
    }

    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i)
    {
        double sum = b[i];
        for (int j = i + 1; j < n; ++j) sum -= A[i][j] * x[j];
        x[i] = sum / A[i][i];
    }
    return x;
}

std::vector<std::vector<double>> MatrixUtils::Add(
    const std::vector<std::vector<double>>& A,
    const std::vector<std::vector<double>>& B)
{
    if (A.size() != B.size() || A[0].size() != B[0].size())
        return {};

    std::vector<std::vector<double>> C(A.size(), std::vector<double>(A[0].size(), 0.0));

    for (size_t i = 0; i < A.size(); ++i)
        for (size_t j = 0; j < A[i].size(); ++j)
            C[i][j] = A[i][j] + B[i][j];

    return C;
}

void MatrixUtils::AddTo(std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B)
{
    if (A.size() != B.size() || A[0].size() != B[0].size())
        return;

    for (size_t i = 0; i < A.size(); ++i)
        for (size_t j = 0; j < A[i].size(); ++j)
            A[i][j] += B[i][j];
}

void MatrixUtils::AddTo(std::vector<double>& A, const std::vector<double>& B)
{
    if (A.size() != B.size())
        return;

    for (size_t i = 0; i < A.size(); ++i)
            A[i] += B[i];
}

std::vector<double> MatrixUtils::Add(
    const std::vector<double>& a,
    const std::vector<double>& b)
{
    if (a.size() != b.size())
        throw std::invalid_argument("Wektory musz? mie? ten sam rozmiar");

    std::vector<double> result(a.size(), 0.0);

    for (size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] + b[i];

    return result;
}

void MatrixUtils::Scale(std::vector<std::vector<double>>& A, double scalar)
{
    for (size_t i = 0; i < A.size(); ++i)
        for (size_t j = 0; j < A[i].size(); ++j)
            A[i][j] *= scalar;
}

void MatrixUtils::Scale(std::vector<double>& A, double scalar)
{
    for (size_t i = 0; i < A.size(); ++i)
        A[i] *= scalar;
}

std::vector<double> MatrixUtils::Multiply(
    const std::vector<std::vector<double>>& A,
    const std::vector<double>& x)
{
    if (A[0].size() != x.size())
        throw std::invalid_argument("Rozmiar wektora musi odpowiada? liczbie kolumn macierzy");

    std::vector<double> y(A.size(), 0.0);

    for (size_t i = 0; i < A.size(); ++i)
        for (size_t j = 0; j < A[i].size(); ++j)
            y[i] += A[i][j] * x[j];

    return y;
}

std::vector<double> MatrixUtils::Subtract(
    const std::vector<double>& a,
    const std::vector<double>& b)
{
    if (a.size() != b.size())
        throw std::invalid_argument("Wektory musz? mie? ten sam rozmiar");

    std::vector<double> result(a.size(), 0.0);

    for (size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] - b[i];

    return result;
}

std::vector<std::vector<double>> MatrixUtils::MultiplyVectors(const std::vector<double>& a, const std::vector<double>& b)
{
    std::vector<std::vector<double>> result(a.size(), std::vector<double>(b.size(), 0.0));

    for (size_t i = 0; i < a.size(); ++i)
    {
        for (size_t j = 0; j < b.size(); ++j)
        {
            result[i][j] = a[i] * b[j];
        }
    }

    return result;
}
