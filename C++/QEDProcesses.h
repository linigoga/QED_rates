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
    double betheHeitlerCrossSection(double k, double atomic_number, double lower_bound, double upper_bound);
    std::tuple<double,double,double,double,double> calculateCCoefficient(double gamma);
    double nonRelativisticCoulombTridentDiffCrossSection(double atomic_number, double gamma, double positron_energy);
    double relativisticCoulombTridentDiffCrossSection(double atomic_number, double gamma, double positron_energy);
    double totalCoulombTridentCrossSectionLowerLimit(double atomic_number, double gamma);
    double totalCoulombTridentCrossSectionUpperLimit(double atomic_number, double gamma);
    double totalCoulombTridentCrossSection(double atomic_number, double gamma);
    double differentialCoulombTridentCrossSectionLowerLimit(double atomic_number, double gamma, double positron_energy);    
    double differentialCoulombTridentCrossSectionUpperLimit(double atomic_number, double gamma, double positron_energy);
    double differentialCoulombTridentCrossSection(double atomic_number, double gamma, double positron_energy);

};


#endif // QEDPROCESSES_H
