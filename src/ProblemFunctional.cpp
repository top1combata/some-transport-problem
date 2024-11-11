#include "ProblemFunctional.h"
#include "workspace.hpp"
#include <cmath>


double ProblemFunctional(DualFunction<double> Y)
{
    using namespace ElementaryDualFunctions;
    using namespace ProblemConstants;

    DualFunction<double> Z = Sin(5*X) * Sin(Y);

    std::function<double(double)> curve = [Y, Z](double x)
    {
        double dy = Y.derivative()(x);
        double dz = Z.derivative()(x);
        return std::sqrt(1 + dy*dy + dz*dz);
    };

    std::function<double(double)> integrand = [curve](double x)
    {
        return beta.function()(x) * curve(x);
    };

    static Workspace<double> workspace(128, 5);

    constexpr double epsabs = 1e-7;
    double abserr, integral1, integral2;

    workspace.qag(curve, startPoint, endPoint, epsabs, 0, integral1, abserr);
    workspace.qag(integrand, startPoint, endPoint, epsabs, 0, integral2, abserr);

    double FunctionalValue = alpha / 2 * integral1 * integral1 + integral2;
    return FunctionalValue;
}