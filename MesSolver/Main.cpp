#include <vector>
#include "MesSolver.h"

int main()
{
	int nGaussPointsNumber = 4;

	MesSolver solver("Test1_4_4.txt", nGaussPointsNumber);
	//MesSolver solver("Test2_4_4_MixGrid.txt", nGaussPointsNumber);
	//MesSolver solver("Test3_31_31_kwadrat.txt", nGaussPointsNumber);

	//solver.CalculateSolution();

	solver.StartTimeSimulation();

    return 0;
}

