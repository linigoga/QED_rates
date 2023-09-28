#include <pybind11/pybind11.h>
#include "QEDProcesses.h"

namespace py = pybind11;

PYBIND11_MODULE(QEDProcesses, m) {
    py::class_<QED>(m, "QED")
        .def(py::init<>())
        .def("thomas_fermi_wavelength", &QED::thomasFermiWavelength, py::arg("atomic_number"))
        .def("i1_function", &QED::i1Function, py::arg("delta"), py::arg("screening_length"), py::arg("q"))
        .def("i2_function", &QED::i2Function, py::arg("delta"), py::arg("screening_length"), py::arg("q"))
        .def("non_relativistic_bremsstrahlung_cross_section", &QED::nonRelativisticBremsstrahlungCrossSection, py::arg("gamma1"), py::arg("k"), py::arg("atomic_number"))
        .def("rel_bremmstrahlung_differential_cross_section", &QED::relBremmstrahlungDifferentialCrossSection, py::arg("atomic_number"), py::arg("k"), py::arg("gamma1"))
        .def("coulomb_correction_term", &QED::coulombCorrectionTerm, py::arg("atomic_number"), py::arg("n"))
        .def("ultra_rel_bremmstrahlung_differential_cross_section", &QED::ultraRelBremmstrahlungDifferentialCrossSection, py::arg("atomic_number"), py::arg("k"), py::arg("gamma1"))
        .def("bremmstrahlung_cross_section", &QED::bremmstrahlungCrossSection, py::arg("atomic_number"), py::arg("gamma"), py::arg("lower_bound"), py::arg("upper_bound"))
        .def("diff_bethe_heitler_cross_section", &QED::diffBetheHeitlerCrossSection, py::arg("gamma_p"), py::arg("k"), py::arg("atomic_number"))
        .def("bethe_heitler_cross_section", &QED::betheHeitlerCrossSection, py::arg("k"), py::arg("atomic_number"), py::arg("lower_bound"), py::arg("upper_bound"))
        .def("calculate_c_coefficient", &QED::calculateCCoefficient, py::arg("gamma"))
        .def("non_relativistic_coulomb_trident_diff_cross_section", &QED::nonRelativisticCoulombTridentDiffCrossSection, py::arg("atomic_number"), py::arg("gamma"), py::arg("positron_energy"))
        .def("relativistic_coulomb_trident_diff_cross_section", &QED::relativisticCoulombTridentDiffCrossSection, py::arg("atomic_number"), py::arg("gamma"), py::arg("positron_energy"))
        .def("total_coulomb_trident_cross_section_lower_limit", &QED::totalCoulombTridentCrossSectionLowerLimit, py::arg("atomic_number"), py::arg("gamma"))
        .def("total_coulomb_trident_cross_section_upper_limit", &QED::totalCoulombTridentCrossSectionUpperLimit, py::arg("atomic_number"), py::arg("gamma"))
        .def("total_coulomb_trident_cross_section", &QED::totalCoulombTridentCrossSection, py::arg("atomic_number"), py::arg("gamma"))
        .def("differential_coulomb_trident_cross_section_lower_limit", &QED::differentialCoulombTridentCrossSectionLowerLimit, py::arg("atomic_number"), py::arg("gamma"), py::arg("positron_energy"))
        .def("differential_coulomb_trident_cross_section_upper_limit", &QED::differentialCoulombTridentCrossSectionUpperLimit, py::arg("atomic_number"), py::arg("gamma"), py::arg("positron_energy"))
        .def("differential_coulomb_trident_cross_section", &QED::differentialCoulombTridentCrossSection, py::arg("atomic_number"), py::arg("gamma"), py::arg("positron_energy"));
    
    py::class_<QEDBlackburn>(m, "QEDBlackburn")
        .def(py::init<>())
        .def("auxilary_function_t", &QEDBlackburn::auxilaryFunctionT, py::arg("chi"))
        .def("auxilary_function_r", &QEDBlackburn::auxilaryFunctionR, py::arg("x"))
        .def("auxilary_function_g", &QEDBlackburn::auxilaryFunctionG, py::arg("x"))
        .def("calculate_chi", &QEDBlackburn::calculateChi, py::arg("gamma_e"), py::arg("a0"), py::arg("omega_0"), py::arg("phi"), py::arg("n"))
        .def("calculate_pair_production_probability", &QEDBlackburn::calculatePairProductionProbability, py::arg("a0"), py::arg("gamma_e"), py::arg("omega_0"), py::arg("omega"), py::arg("n"))
        .def("calculate_radiated_energy", &QEDBlackburn::calculateRadiatedEnergy, py::arg("gamma_e"), py::arg("a0"), py::arg("omega_0"), py::arg("chi_c"))
        .def("calculate_critical_chi_rr", &QEDBlackburn::calculateCriticalChiRR, py::arg("gamma_e"), py::arg("a0"), py::arg("omega_0"), py::arg("chi"))
        .def("critical_chi", &QEDBlackburn::criticalChi, py::arg("gamma_e"), py::arg("a0"), py::arg("omega_0"), py::arg("n"))
        .def("critical_chi_modulus", &QEDBlackburn::criticalChiModulus, py::arg("gamma_e"), py::arg("a0"), py::arg("omega_0"), py::arg("n"))
        .def("critical_frequency", &QEDBlackburn::criticalFrequency, py::arg("gamma_e"), py::arg("a0"), py::arg("omega_0"), py::arg("chi_c"))
        .def("calculate_critical_phase", &QEDBlackburn::calculateCriticalPhase, py::arg("gamma_e"), py::arg("a0"), py::arg("omega_0"), py::arg("n"))
        .def("calculate_critical_phase_modulus", &QEDBlackburn::calculateCriticalPhaseModulus, py::arg("gamma_e"), py::arg("a0"), py::arg("omega_0"), py::arg("n"))
        .def("correction_factor_fhe", &QEDBlackburn::correctionFactorFhe, py::arg("n"), py::arg("phi_c"))
        .def("calculate_photon_energy_spectrum", &QEDBlackburn::calculatePhotonEnergySpectrum, py::arg("gamma_e"), py::arg("a0"), py::arg("omega_0"), py::arg("omega"), py::arg("n"))
        .def("positron_yield", &QEDBlackburn::positronYield, py::arg("gamma_e"), py::arg("a0"), py::arg("omega_0"), py::arg("n"))
        .def("calculate_chi_as_function_of_phi", &QEDBlackburn::calculateChiAsFunctionOfPhi, py::arg("gamma_e"), py::arg("a0"), py::arg("omega_0"), py::arg("phi"), py::arg("n"))
        .def("calculate_gamma_as_function_of_phi", &QEDBlackburn::calculateGammaAsFunctionOfPhi, py::arg("gamma_e"), py::arg("a0"), py::arg("omega_0"), py::arg("phi"), py::arg("n"))
        .def("calculate_variance", &QEDBlackburn::calculateVariance, py::arg("gamma_e"), py::arg("ao"), py::arg("omega_0"), py::arg("n"), py::arg("chi_c"));

        
    py::class_<QEDReconstructionMethods>(m, "QEDReconstructionMethods")
        .def(py::init<>())
        .def("particle_binning_amaro", &QEDReconstructionMethods::particleBinningAmaro, py::arg("cache"), py::arg("beam_type"))
        .def("make_3d_gauss_distribution", &QEDReconstructionMethods::make3dGaussDistribution)
        .def("calculate_positrons_produced_from_beam", &QEDReconstructionMethods::calculatePositronsProducedFromBeam);
}