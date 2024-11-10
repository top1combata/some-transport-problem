#pragma once

#include "DualFunction.h"

namespace ProblemConstants
{
    // interval [startpoint, endpoint]
    constexpr double startPoint = 0.03;
    constexpr double endPoint = 1;

    // boundary values for problem
    constexpr double startValue = 0;
    constexpr double endValue = 1;
    
    constexpr double alpha = 0.1;

    using namespace ElementaryDualFunctions;
    static const DualFunction<double> beta(0.5);
}

double ProblemFunctional(DualFunction<double> y);