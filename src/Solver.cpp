#include "Solver.h"
#include "ProblemFunctional.h"
#include <numeric>
#include <iostream>
#include <cassert>
#include <Eigen/Dense>

using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;


DualFunction<double> getY(const Vector& coeffs);
Vector gradient(const Vector& coeffs, double);


DualFunction<double> RayleighRitz(unsigned n, double eps, double step)
{
    assert(eps > 0);
    assert(step > 0);
    assert(n != 0);

    Vector coeffs = Vector::Zero(n, 1);

    double error = std::numeric_limits<double>::max();
    Vector velocity = Vector::Zero(n, 1);
    constexpr double beta = 0.9;

    for (int iterNum = 0; error > eps; iterNum++)
    {
        velocity = beta * velocity + step * gradient(coeffs - beta * velocity, eps);
        error = velocity.norm();
        coeffs -= velocity;
        std::cout << "Error : " << std::scientific << error << '\n';
    }

    return getY(coeffs);
}


DualFunction<double> basisFunction(unsigned k)
{
    using namespace ProblemConstants;
    using namespace ElementaryDualFunctions;
    using std::numbers::pi;

    return Sin(k * pi / (endPoint - startPoint) * (X - DualFunction<double>(startPoint)));

    // DualFunction<double> base = (X - DualFunction<double>(startPoint)) * (DualFunction<double>(endPoint) - X);
    // return base.pow(k+1);
}


DualFunction<double> getY(const Vector& coeffs)
{
    unsigned n = coeffs.size();
    assert(n != 0);

    using namespace ProblemConstants;
    using namespace ElementaryDualFunctions;

    double bias = (startValue * endPoint - endValue * startPoint) / (endPoint - startPoint);
    double slope = (endValue - startValue) / (endPoint - startPoint);

    // line connecting boundary points
    DualFunction<double> Y = DualFunction<double>(bias) + slope * X;

    for (unsigned k = 1; k <= n; k++)
        Y = Y + coeffs[k-1] * basisFunction(k);

    return Y;
}


Vector gradient(const Vector& coeffs, double dC)
{
    unsigned n = coeffs.size();
    assert(n != 0);

    Vector grad = Vector(n, 1);

    for (int k = 0; k < n; k++)
    {
        Vector delta = Vector::Zero(n, 1);
        delta[k] = dC;
        
        // second order scheme
        // [f(x + h) - f(x - h)] / 2h
        double nextValue = ProblemFunctional(getY(coeffs + delta));
        double prevValue = ProblemFunctional(getY(coeffs - delta));
        grad[k] = (nextValue - prevValue) / (2 * dC);
    }

    return grad;
}