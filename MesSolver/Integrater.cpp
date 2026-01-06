#include "Integrater.h"
#include "Constants.h"
#include "Helpers.h"

double Integrater::Integrate1D(const std::vector<double>& values, double dx, int n)
{
    std::vector<double> weights = getGaussWeights(n);

    double result = 0.0;

    for (int i = 0; i < n; ++i)
    {
        result += weights.at(i) * values.at(i) * dx;
    }

    return result;
}

std::vector<std::vector<double>> Integrater::Integrate2D(const std::vector<std::vector<double>>& valuesA, const std::vector<std::vector<double>>& valuesB, std::vector<double> dz, int n)
{
    std::vector<double> weights = getGaussWeights(n);
    std::vector<std::vector<double>> result{ valuesA[0].size(), std::vector<double>(valuesA[0].size()) };

	int numberOfIntegralPoints = valuesA.size();
	int integralPointIndex = 0;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            std::vector<std::vector<double>> matrix = MatrixUtils::MultiplyVectors(valuesA.at(integralPointIndex), valuesA.at(integralPointIndex));
            matrix = MatrixUtils::Add(matrix, MatrixUtils::MultiplyVectors(valuesB.at(integralPointIndex), valuesB.at(integralPointIndex)));
            matrix = MatrixUtils::Scale(matrix, dz.at(integralPointIndex) * weights.at(i) * weights.at(j));

            result = MatrixUtils::Add(result, matrix);

            integralPointIndex++;
        }
    }

    return result;
}
