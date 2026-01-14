#include "Constants.h"

std::vector<double> getGaussWeights(int n)
{
    std::vector<double> weights{};

    switch (n)
    {
        case 1:
            weights = { 2.0 };
            break;
        case 2:
            weights = { 1.0, 1.0 };
            break;
        case 3:
            weights = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
            break;
        case 4:
            weights = {
                (18.0 - sqrt(30.0)) / 36.0,
                (18.0 + sqrt(30.0)) / 36.0,
                (18.0 + sqrt(30.0)) / 36.0,
                (18.0 - sqrt(30.0)) / 36.0
            };
            break;
        case 5:
            weights = {
                (322.0 - 13.0 * sqrt(70.0)) / 900.0,
                (322.0 + 13.0 * sqrt(70.0)) / 900.0,
                128.0 / 225.0,
                (322.0 + 13.0 * sqrt(70.0)) / 900.0,
                (322.0 - 13.0 * sqrt(70.0)) / 900.0
            };
            break;
        default:
            break;
    }

    return weights;
}

std::vector<double> getGaussValues(int n)
{
    std::vector<double> values{};

    switch (n)
    {
        case 1:
            values = { 0.0 };
            break;
        case 2:
            values = { -sqrt(1.0 / 3.0),  sqrt(1.0 / 3.0) };
            break;
        case 3:
            values = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };
            break;
        case 4:
            values = {
                -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)),
                -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)),
                 sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)),
                 sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0))
            };
            break;
        case 5:
            values = {
                -1.0 / 3.0 * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)),
                -1.0 / 3.0 * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)),
                 0.0,
                 1.0 / 3.0 * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)),
                 1.0 / 3.0 * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0))
            };
            break;
        default:
            break;
    }

    return values;
}
