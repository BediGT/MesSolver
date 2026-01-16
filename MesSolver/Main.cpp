#include "../MesSolverLib/MesSolver.h"
#include <vector>

int main()
{
	int nGaussPointsNumber = 3;
	bool bShouldPrint = true;

	MesSolver solver(nGaussPointsNumber, bShouldPrint);
	solver.LoadData("..\\Grids\\Test1_4_4.txt");
	//solver.LoadData("..\\Grids\\Test2_4_4_MixGrid.txt");
	//solver.LoadData("..\\Grids\\Test3_31_31_kwadrat.txt");

	//solver.CalculateSolution();

	solver.SimulateWithTime();

    return 0;
}

