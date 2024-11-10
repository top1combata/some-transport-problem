#pragma once

#include <functional>
#include <concepts>
#include <cmath>


template<class Real>
class DualFunction
{
public:

    using FunctionType = std::function<Real(Real)>;

    DualFunction(FunctionType _function, FunctionType _derivative) :
        m_function(_function),
        m_derivative(_derivative) 
        {}

    template<class Scalar>
    explicit DualFunction(Scalar scalar) 
    requires requires (Scalar scalar){static_cast<Real>(scalar);}
        :
        m_function([scalar](Real){return Real(scalar);}),
        m_derivative([](Real){return Real(0);})
        {}

    const FunctionType& function() const
    {
        return m_function;
    }

    const FunctionType& derivative() const
    {
        return m_derivative;
    }


    DualFunction operator()(const DualFunction& g) const
    {
        auto function = [f = this->function(), g = g.function()](Real x)
        {
            return f(g(x));
        };
        auto derivative = [df = this->derivative(), g = g.function(), dg = g.derivative()](Real x)
        {
            return df(g(x)) * dg(x);
        };
        return DualFunction(function, derivative);
    };

private:

    FunctionType m_function;
    FunctionType m_derivative;
};


template<class Real>
DualFunction<Real> operator+(const DualFunction<Real>& f, const DualFunction<Real>& g)
{
    auto function = [f = f.function(), g = g.function()](Real x)
    {
        return f(x) + g(x);
    };

    auto derivative = [df = f.derivative(), dg = g.derivative()](Real x)
    {
        return df(x) + dg(x);
    };
    return DualFunction<Real>(function, derivative);
}


template<class Real>
DualFunction<Real> operator-(const DualFunction<Real>& f, const DualFunction<Real>& g)
{
    auto function = [f = f.function(), g = g.function()](Real x)
    {
        return f(x) - g(x);
    };

    auto derivative = [df = f.derivative(), dg = g.derivative()](Real x)
    {
        return df(x) - dg(x);
    };
    return DualFunction<Real>(function, derivative);
}


template<class Real>
DualFunction<Real> operator*(const DualFunction<Real>& f, const DualFunction<Real>& g)
{
    auto function = [f = f.function(), g = g.function()](Real x)
    {
        return f(x) * g(x);
    };

    auto derivative = [f = f.function(), g = g.function(), df = f.derivative(), dg = g.derivative()](Real x)
    {
        return f(x) * dg(x) + df(x) * g(x);
    };
    return DualFunction<Real>(function, derivative);
}


template<class Real>
DualFunction<Real> operator/(const DualFunction<Real>& f, const DualFunction<Real>& g)
{
    auto function = [f = f.function(), g = g.function()](Real x)
    {
        return f(x) / g(x);
    };

    auto derivative = [f = f.function(), g = g.function(), df = f.derivative(), dg = g.derivative()](Real x)
    {
        Real gx = g(x);
        return (df(x) * gx - f(x) * dg(x)) / gx / gx;
    };
    return DualFunction<Real>(function, derivative);
}


template<class Real, class Scalar>
requires requires (Scalar scalar){static_cast<Real>(scalar);}
DualFunction<Real> operator*(const DualFunction<Real>& f, Scalar scalar)
{
    auto function = [f = f.function(), scalar](Real x)
    {
        return scalar * f(x);
    };

    auto derivative = [df = f.derivative(), scalar](Real x)
    {
        return scalar * df(x);
    };
    return DualFunction<Real>(function, derivative);
}


template<class Real, class Scalar>
DualFunction<Real> operator*(Scalar scalar, const DualFunction<Real>& f)
requires requires (Scalar scalar){static_cast<Real>(scalar);}
{
    return f * scalar;
}


template<class Real>
DualFunction<Real> operator-(const DualFunction<Real>& f)
{
    return (-1) * f;
}

namespace ElementaryDualFunctions
{
    static DualFunction<double> X([](double x){return x;}, [](double){return 1;});
    static DualFunction<double> Sin([](double x){return std::sin(x);}, [](double x){return std::cos(x);});
    static DualFunction<double> Cos([](double x){return std::cos(x);}, [](double x){return -std::sin(x);});
    static DualFunction<double> Exp([](double x){return std::exp(x);}, [](double x){return std::exp(x);});
    static DualFunction<double> Log([](double x){return std::log(x);}, [](double x){return 1/x;});
    static DualFunction<double> Sqrt([](double x){return std::sqrt(x);}, [](double x){return 1 / sqrt(x) / 2;});
}