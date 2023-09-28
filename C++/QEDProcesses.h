// QEDProcesses.h
#ifndef QEDPROCESSES_H
#define QEDPROCESSES_H

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>
#include <map>


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

class QEDBlackburn{
    public:
    QEDBlackburn();

    double auxilaryFunctionT(double chi);
    double auxilaryFunctionR(double x);
    double auxilaryFunctionG(double x);
    double calculateChi(double gamma_e, double a0, double omega_0, double phi, double n);
    double calculatePairProductionProbability(double a0, double gamma_e, double omega_0, double omega, double n);
    double calculateRadiatedEnergy(double gamma_e, double a0, double omega_0, double chi_c);
    double calculateCriticalChiRR(double gamma_e,double a0, double omega_0, double chi);
    double criticalChi(double gamma_e, double a0, double omega_0, double n);
    double criticalChiModulus(double gamma_e, double a0, double omega_0, double n);
    double criticalFrequency(double gamma_e, double a0, double omega_0, double chi_c);
    double calculateVariance(double gamma_e, double ao, double omega_0, double n, double chi_c);
    double calculateCriticalPhase(double gamma_e, double a0, double omega_0, double n);
    double calculateCriticalPhaseModulus(double gamma_e, double a0, double omega_0, double n);
    double correctionFactorFhe(double n, double phi_c);
    double calculatePhotonEnergySpectrum(double gamma_e, double a0, double omega_0, double omega, double n);
    double positronYield(double gamma_e, double a0, double omega_0,double n);
    double calculateChiAsFunctionOfPhi(double gamma_e, double a0, double omega_0, double phi, double n);
    double calculateGammaAsFunctionOfPhi(double gamma_e, double a0, double omega_0, double phi, double n);
  //  double calculateProbabilityDensity(double gamma_e,double a0, double omega_0, double n);
};

class QEDReconstructionMethods{
    public:
    QEDReconstructionMethods();

    double particleBinningAmaro(std::map<std::string,float> cache, std::string beam_type);
    double make3dGaussDistribution(double x, double y, double z, double a0, double waist, double wavelength);
    double calculatePositronsProducedFromBeam(double gamma_e, double a0, double omega_0,double delta, double r_, double waist, double ne,double n);
};

#endif // QEDPROCESSES_H
