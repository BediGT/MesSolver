#include "ShapeFunctions.h"

double shapeFunction(int index, double ksi, double eta)
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

double shapeFunctionDer_Ksi(int index, double eta)
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

double shapeFunctionDer_Eta(int index, double ksi)
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
