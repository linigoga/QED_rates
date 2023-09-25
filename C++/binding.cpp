#include <pybind11/pybind11.h>
#include "QEDProcesses.h"

namespace py = pybind11;

PYBIND11_MODULE(QEDProcesses, m) {
    py::class_<QED>(m, "QED")
        .def(py::init<>())
        .def("thomasFermiWavelength", &QED::thomasFermiWavelength, py::arg("atomic_number"))
        .def("i1Function", &QED::i1Function, py::arg("delta"), py::arg("screening_length"), py::arg("q"))
        .def("i2Function", &QED::i2Function, py::arg("delta"), py::arg("screening_length"), py::arg("q"))
        .def("nonRelativisticBremsstrahlungCrossSection", &QED::nonRelativisticBremsstrahlungCrossSection, py::arg("gamma1"), py::arg("k"), py::arg("atomic_number"))
        .def("relBremmstrahlungDifferentialCrossSection", &QED::relBremmstrahlungDifferentialCrossSection, py::arg("atomic_number"), py::arg("k"), py::arg("gamma1"))
        .def("coulombCorrectionTerm", &QED::coulombCorrectionTerm, py::arg("atomic_number"), py::arg("n"))
        .def("ultraRelBremmstrahlungDifferentialCrossSection", &QED::ultraRelBremmstrahlungDifferentialCrossSection, py::arg("atomic_number"), py::arg("k"), py::arg("gamma1"))
        .def("bremmstrahlungCrossSection", &QED::bremmstrahlungCrossSection, py::arg("atomic_number"), py::arg("gamma"), py::arg("lower_bound"), py::arg("upper_bound"))
        .def("diffBetheHeitlerCrossSection", &QED::diffBetheHeitlerCrossSection, py::arg("gamma_p"), py::arg("k"), py::arg("atomic_number"))
        .def("betheHeitlerCrossSection", &QED::betheHeitlerCrossSection, py::arg("k"), py::arg("atomic_number"), py::arg("lower_bound"), py::arg("upper_bound"))
        .def("calculateCCoefficient", &QED::calculateCCoefficient, py::arg("gamma"))
        .def("nonRelativisticCoulombTridentDiffCrossSection", &QED::nonRelativisticCoulombTridentDiffCrossSection, py::arg("atomic_number"), py::arg("gamma"), py::arg("positron_energy"))
        .def("relativisticCoulombTridentDiffCrossSection", &QED::relativisticCoulombTridentDiffCrossSection, py::arg("atomic_number"), py::arg("gamma"), py::arg("positron_energy"))
        .def("totalCoulombTridentCrossSectionLowerLimit", &QED::totalCoulombTridentCrossSectionLowerLimit, py::arg("atomic_number"), py::arg("gamma"))
        .def("totalCoulombTridentCrossSectionUpperLimit", &QED::totalCoulombTridentCrossSectionUpperLimit, py::arg("atomic_number"), py::arg("gamma"))
        .def("totalCoulombTridentCrossSection", &QED::totalCoulombTridentCrossSection, py::arg("atomic_number"), py::arg("gamma"))
        .def("differentialCoulombTridentCrossSectionLowerLimit", &QED::differentialCoulombTridentCrossSectionLowerLimit, py::arg("atomic_number"), py::arg("gamma"), py::arg("positron_energy"))
        .def("differentialCoulombTridentCrossSectionUpperLimit", &QED::differentialCoulombTridentCrossSectionUpperLimit, py::arg("atomic_number"), py::arg("gamma"), py::arg("positron_energy"))
        .def("differentialCoulombTridentCrossSection", &QED::differentialCoulombTridentCrossSection, py::arg("atomic_number"), py::arg("gamma"), py::arg("positron_energy"));
}
