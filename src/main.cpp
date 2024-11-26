#include "Solver.h"
#include "ProblemFunctional.h"
#include "DualFunction.h"
#include <fstream>
#include <iostream>
#include <iomanip>


int main(int argc, char** argv)
{
    if (argc != 4)
    {
        std::cout << "usage : " << argv[0] << " n eps step\n";
        std::cout << "n - basis size\n";
        std::cout << "eps - tolerance treshold\n";
        std::cout << "step - step in gradient descent\n";
        return 1;
    }
    unsigned n   = std::stoul(argv[1]);
    double eps   = std::stod(argv[2]);
    double step  = std::stod(argv[3]);

    auto Y = RayleighRitz(n, eps, step);
    std::cout << "\nFunctional value " << std::defaultfloat << std::setprecision(10) << ProblemFunctional(Y) << '\n';

    using namespace ProblemConstants;

    constexpr int numPoints = 100;
    double x = startPoint;
    double delta = (endPoint - startPoint) / (numPoints - 1);

    std::ofstream fout("data.txt");
    for (int i = 0; i < numPoints; i++, x += delta)
        fout << x << ' ' << Y.function()(x) << '\n';

    return 0;
}