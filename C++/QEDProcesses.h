// QEDProcesses.h
#ifndef QEDPROCESSES_H
#define QEDPROCESSES_H

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>


class QED{
public:
    QED();

    double thomasFermiWavelength(int atomic_number);
    double i1Function(double delta, double screening_length, double q = 1.0);
    double i2Function(double delta, double screening_length, double q = 1.0);
    double nonRelativisticBremsstrahlungCrossSection(double gamma1, double k, double atomic_number);
    double relBremmstrahlungDifferentialCrossSection(double atomic_number, double k, double gamma1);
    double coulombCorrectionTerm(double atomic_number, int n = 10);
    double ultraRelBremmstrahlungDifferentialCrossSection(double atomic_number, double k,double gamma1);
    double bremmstrahlungCrossSection(double atomic_number, double gamma, double lower_bound, double upper_bound);
    double diffBetheHeitlerCrossSection(double gamma_p, double k, double atomic_number);
    double betheHeitlerCrossSection(double k, double atomic_number, double lower_bound, double upper_bound); // Add this line
};


#endif // QEDPROCESSES_H
