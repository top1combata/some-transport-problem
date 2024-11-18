#include "Solver.h"
#include "ProblemFunctional.h"
#include <numeric>
#include <iostream>
#include <cassert>
#include <Eigen/Dense>

using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;


DualFunction<double> getY(const Vector& coeffs);
Vector gradient(const Vector& coeffs, double);
double goldenRatioSearch(const std::function<double(double)>& f, double l, double r, double eps);


DualFunction<double> RayleighRitz(unsigned n, double eps, double step)
{
    assert(eps > 0);
    assert(step > 0);
    assert(n != 0);

    Vector coeffs = Vector::Zero(n, 1);
    Vector prevGrad = gradient(coeffs, step * eps);
    Vector grad  = prevGrad;
    Vector conjugate = -prevGrad;

    std::function<double(double)> f = [&coeffs, &conjugate](double x)
    {
        try
        {
            return ProblemFunctional(getY(coeffs + x * conjugate));
        }
        catch (const char* message)
        {
            std::cout << "ERROR : " << message << '\n';
            std::exit(0);
        }
    };

    double error = std::numeric_limits<double>::max();

    for (int iterNum = 0; error > eps; iterNum++)
    {
        prevGrad = grad;
        grad = gradient(coeffs, step * eps);

        // Polak–Ribiere
        double beta = (iterNum ? grad.dot(grad - prevGrad) / prevGrad.dot(prevGrad) : 0);
        beta = std::max(beta, 0.0);

        conjugate = -grad + beta * conjugate;
        double alpha = goldenRatioSearch(f, 0, step, eps * step);

        coeffs += alpha * conjugate;
        // error = (coeffs - prevCoeffs).norm();
        error = grad.norm();
        std::cout << "Residual " << std::scientific << error << '\n';
    }

    return getY(coeffs);
}


double goldenRatioSearch(const std::function<double(double)>& f, double l, double r, double eps)
{
    using std::numbers::phi;

    struct point
    {
        double x, y;
    };

    // left.x < middleLeft.x < middleRight.x < right.x
    point left, right, middleLeft, middleRight;
    left.x        = l;
    right.x       = r;
    middleLeft.x  = (r - l) / (phi*phi) + l;
    middleRight.x = right.x - middleLeft.x + left.x;
    
    left.y        = f(left.x);
    middleLeft.y  = f(middleLeft.x);
    middleRight.y = f(middleRight.x);
    right.y       = f(right.x);

    while (right.x - left.x > eps)
    {
        if (middleLeft.y < middleRight.y)
        {
            right = middleRight;
            middleRight = middleLeft;

            middleLeft.x = left.x + right.x - middleRight.x;
            middleLeft.y = f(middleLeft.x);
        }
        else
        {
            left = middleLeft;
            middleLeft = middleRight;

            middleRight.x = right.x - middleLeft.x + left.x;
            middleRight.y = f(middleRight.x);
        }
    }
    
    return (left.x + right.x) / 2;
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