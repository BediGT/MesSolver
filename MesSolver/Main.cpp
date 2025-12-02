#include <vector>
#include "MesSolver.h"

int main()
{
	MesSolver solver;

	//solver.LoadData("Test1_4_4.txt");
	//solver.LoadData("Test2_4_4_MixGrid.txt");
	solver.LoadData("Test3_31_31_kwadrat.txt");
	solver.PrintData();

	//solver.CalculateSolution();

	solver.StartTimeSimulation();

    return 0;
}

