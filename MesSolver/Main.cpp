#include "../MesSolverLib/MesSolver.h"
#include <vector>

int main()
{
	int nGaussPointsNumber = 4;
	bool bShouldPrint = false;

	MesSolver solver(nGaussPointsNumber, bShouldPrint);
	//solver.LoadData("..\\Grids\\Test1_4_4.txt");
	solver.LoadData("..\\Grids\\Test2_4_4_MixGrid.txt");
	//solver.LoadData("..\\Grids\\Test3_31_31_kwadrat.txt");

	solver.Simulate();

    return 0;
}

