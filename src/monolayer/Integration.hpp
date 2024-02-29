#include <cstdlib>
#include <iostream>

#ifndef __INTEGRATION_H__
#define __INTEGRATION_H__

struct Quadrature {
    // Abstract base class for elementary quadrature algorithms.
    long n; // Current level of refinement.

    virtual double next() = 0;
    // Returns the value of the integral at the nth stage of refinement.
    // The function next() must be defined in the derived class.
};

template <typename T> struct Integrand {
    T& func;
    void** args;
    Integrand(T& f) : func(f), args(nullptr) {}
};

template <class T> struct Trapzd : Quadrature {
    double a, b, s; // Limits of integration and current value of integral.
    Integrand<T>& integrand;

    // Trapzd(){};

    // func is function or functor to be integrated between limits: a and b
    Trapzd(Integrand<T>& integrand, const double aa, const double bb)
        : a(aa), b(bb), integrand(integrand) {
        n = 0;
    }

    // Returns the nth stage of refinement of the extended trapezoidal rule.
    // On the first call (n = 1), the routine returns the crudest estimate
    // of integral of f x / dx in [a,b]. Subsequent calls set n=2,3,... and
    // improve the accuracy by adding 2n - 2 additional interior points.
    double next() {
        double lower, upper, tnm, sum, del;
        long it, j;
        n++;

        if (n == 1) {
            return (s = (b - a) *
                        (0.25 * integrand.func(a, integrand.args) +
                         0.25 * integrand.func(b, integrand.args) +
                         0.5 * integrand.func((a + b) / 2.0, integrand.args)));
        } else {
            for (it = 1, j = 1; j < n - 1; j++) {
                it <<= 1;
            }
            tnm = it;
            // This is the spacing of the points to be added.
            del = (b - a) / tnm;

            for (sum = 0.0, j = 0; j < tnm; j++) {
                lower = a + j * del;
                upper = lower + del;
                sum +=
                    0.25 * integrand.func(lower, integrand.args) +
                    0.25 * integrand.func(upper, integrand.args) +
                    0.5 * integrand.func((lower + upper) / 2.0, integrand.args);
            }
            // This replaces s by its refined value.
            // s = 0.5 * (s + (b - a) * sum / tnm);
            s = (b - a) * sum / tnm;
            return s;
        }
    }
};

template <typename T>
double qtrap(Integrand<T>& func, const double a, const double b,
             const double eps = 1.0e-10) {
    // Returns the integral of the function or functor func from a to b.
    // The constants EPS can be set to the desired fractional accuracy and
    // JMAX so that 2 to the power JMAX-1 is the maximum allowed number of
    // steps. longegration is performed by the trapezoidal rule.

    const long JMAX = 20;
    double s, olds = 0.0; // Initial value of olds is arbitrary.

    Trapzd<T> t(func, a, b);

    for (long j = 0; j < JMAX; j++) {
        s = t.next();

        if (j > 5) // Avoid spurious early convergence.
        {
            if (abs(s - olds) < eps * abs(olds) || (s == 0.0 && olds == 0.0)) {
                return s;
            }
        }
        olds = s;
    }
#ifdef DEBUG
    std::cerr << "Too many steps in routine qtrap" << std::endl;
    std::cerr << "Last value: " << olds << std::endl;
#endif
    // throw("Too many steps in routine qtrap");
    return s;
}

template <typename T>
double quad(T& function, double lower, double upper, void** args) {
    Integrand<T> f(function); // Functor func here has no parameters.
    f.args = args;
    return qtrap(f, lower, upper, 1e-6);
}
#endif
