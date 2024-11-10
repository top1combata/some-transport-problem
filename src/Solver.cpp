#include "Solver.h"
#include "ProblemFunctional.h"
#include <numeric>
#include <iostream>
#include <cassert>
#include <Eigen/Dense>

using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;


DualFunction<double> getY(const Vector& coeffs, unsigned n);
Vector gradient(Vector coeffs, double);


DualFunction<double> RayleighRitz(unsigned n, double eps)
{
    Vector coeffs = Vector::Zero(n+1, 1);

    double error = std::numeric_limits<double>::max();
    double step = 2;
    constexpr double ratio = 1.05;

    while (error > eps)
    {
        Vector grad = gradient(coeffs, step/10);
        error = grad.norm();
        std::cout << "Error : " << std::scientific << error << '\n';

        coeffs += -grad * (step / error);

        step /= ratio;
    }

    return getY(coeffs, n);
}


DualFunction<double> BasisFunction(unsigned k)
{
    using namespace ProblemConstants;
    using namespace ElementaryDualFunctions;
    using std::numbers::pi;

    return Sin(pi / (endPoint - startPoint) * (X - DualFunction<double>(startPoint)));
}


DualFunction<double> getY(const Vector& coeffs, unsigned n)
{
    using namespace ProblemConstants;
    using namespace ElementaryDualFunctions;

    double bias = (startValue * endPoint - endValue * startPoint) / (endPoint - startPoint);
    double slope = (endValue - startValue) / (endPoint - startPoint);

    DualFunction<double> Y = DualFunction<double>(bias) + slope * X;

    for (unsigned k = 0; k <= n; k++)
        Y = Y + coeffs[k] * BasisFunction(k);

    return Y;
}


Vector gradient(Vector coeffs, double dC)
{
    unsigned n = coeffs.size() - 1;
    Vector grad = Vector(coeffs.size(), 1);

    double functionalValue = ProblemFunctional(getY(coeffs, n));

    for (int k = 0; k < grad.size(); k++)
    {
        double pivot = coeffs[k];
        coeffs[k] += dC;
        double newValue = ProblemFunctional(getY(coeffs, n));
        coeffs[k] = pivot;
        
        grad[k] = (newValue - functionalValue) / dC;
    }

    return grad;
}