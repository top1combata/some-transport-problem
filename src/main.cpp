#include "Solver.h"
#include "ProblemFunctional.h"
#include "DualFunction.h"
#include <fstream>
#include <iostream>


int main(int argc, char** argv)
{

    auto Y = RayleighRitz(20, 1e-6);
    std::cout << "\nFunctional value " << std::defaultfloat <<  ProblemFunctional(Y) << '\n';

    using namespace ElementaryDualFunctions;
    using namespace ProblemConstants;

    constexpr int numPoints = 100;
    double x = startPoint;
    double delta = (endPoint - startPoint) / (numPoints - 1);

    std::ofstream fout("data.txt");
    for (int i = 0; i < numPoints; i++, x += delta)
        fout << x << ' ' << Y.function()(x) << '\n';

    return 0;
}