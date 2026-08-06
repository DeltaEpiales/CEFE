#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <memory>

#include "cefe/geometry/CausalDiamond.hpp"
#include "cefe/geometry/Metric.hpp"
#include "cefe/core/EntropicFieldEngine.hpp"
#include "cefe/core/QuantumWaveEngine.hpp"
#include "cefe/core/LatticeQFTEngine.hpp"
#include "../../dynamic_tensors/DynamicTensorEngine.hpp"

namespace py = pybind11;

PYBIND11_MODULE(cefe_core, m) {
    m.doc() = "Causal Entropic Field Engine core C++ extension";

    // Bind BoundaryType Enum
    py::enum_<cefe::core::BoundaryType>(m, "BoundaryType")
        .value("DIRICHLET", cefe::core::BoundaryType::DIRICHLET)
        .value("NEUMANN", cefe::core::BoundaryType::NEUMANN)
        .export_values();

    // Bind Metric Base Class
    py::class_<cefe::geometry::Metric, std::shared_ptr<cefe::geometry::Metric>>(m, "Metric");

    // Bind MinkowskiMetric
    py::class_<cefe::geometry::MinkowskiMetric, cefe::geometry::Metric, std::shared_ptr<cefe::geometry::MinkowskiMetric>>(m, "MinkowskiMetric")
        .def(py::init<>());

    // Bind SchwarzschildMetric
    py::class_<cefe::geometry::SchwarzschildMetric, cefe::geometry::Metric, std::shared_ptr<cefe::geometry::SchwarzschildMetric>>(m, "SchwarzschildMetric")
        .def(py::init<double>(), py::arg("Rs"));

    // Bind AntiDeSitterMetric
    py::class_<cefe::geometry::AntiDeSitterMetric, cefe::geometry::Metric, std::shared_ptr<cefe::geometry::AntiDeSitterMetric>>(m, "AntiDeSitterMetric")
        .def(py::init<double>(), py::arg("L"));

    // Bind FLRWMetric
    py::class_<cefe::geometry::FLRWMetric, cefe::geometry::Metric, std::shared_ptr<cefe::geometry::FLRWMetric>>(m, "FLRWMetric")
        .def(py::init<>());

    // Bind KerrMetric
    py::class_<cefe::geometry::KerrMetric, cefe::geometry::Metric, std::shared_ptr<cefe::geometry::KerrMetric>>(m, "KerrMetric")
        .def(py::init<double, double>(), py::arg("M"), py::arg("a"));

    // Bind GridPoint struct
    py::class_<cefe::geometry::GridPoint>(m, "GridPoint")
        .def_readonly("t", &cefe::geometry::GridPoint::t)
        .def_readonly("x", &cefe::geometry::GridPoint::x)
        .def_readonly("y", &cefe::geometry::GridPoint::y)
        .def_readonly("z", &cefe::geometry::GridPoint::z)
        .def_readonly("is_boundary", &cefe::geometry::GridPoint::is_boundary)
        .def_readonly("index", &cefe::geometry::GridPoint::index);

    // Bind CausalDiamondGrid
    py::class_<cefe::geometry::CausalDiamondGrid, std::shared_ptr<cefe::geometry::CausalDiamondGrid>>(m, "CausalDiamondGrid")
        .def(py::init<double, double, std::shared_ptr<cefe::geometry::Metric>>(), 
             py::arg("radius"), py::arg("spacing"), py::arg("metric") = std::make_shared<cefe::geometry::MinkowskiMetric>())
        .def("get_radius", &cefe::geometry::CausalDiamondGrid::get_radius)
        .def("get_spacing", &cefe::geometry::CausalDiamondGrid::get_spacing)
        .def("get_dof", &cefe::geometry::CausalDiamondGrid::get_dof)
        .def("get_area_bound", &cefe::geometry::CausalDiamondGrid::get_area_bound, py::arg("subregion_radius"))
        .def("get_points", &cefe::geometry::CausalDiamondGrid::get_points);

    // Bind EntropicFieldEngine
    py::class_<cefe::core::EntropicFieldEngine>(m, "EntropicFieldEngine")
        .def(py::init<std::shared_ptr<cefe::geometry::CausalDiamondGrid>>(), py::arg("grid"))
        .def("set_mass", &cefe::core::EntropicFieldEngine::set_mass, py::arg("m"))
        .def("set_boundary_condition", &cefe::core::EntropicFieldEngine::set_boundary_condition, py::arg("bc"))
        .def("set_interaction_coupling", &cefe::core::EntropicFieldEngine::set_interaction_coupling, py::arg("lambda"))
        .def("evolve_to_slice", &cefe::core::EntropicFieldEngine::evolve_to_slice, py::arg("target_t"), py::arg("temperature") = 0.0)
        .def("evolve_to_slice_sparse", &cefe::core::EntropicFieldEngine::evolve_to_slice_sparse, py::arg("target_t"), py::arg("num_modes"), py::arg("temperature") = 0.0)
        .def("compute_entanglement_entropy", &cefe::core::EntropicFieldEngine::compute_entanglement_entropy, py::arg("subregion_radius"))
        .def("compute_entanglement_entropy_indices", &cefe::core::EntropicFieldEngine::compute_entanglement_entropy_indices, py::arg("subregion_indices"))
        .def("initialize_field_state", &cefe::core::EntropicFieldEngine::initialize_field_state, py::arg("target_t"), py::arg("temperature") = 0.0)
        .def("step_forward", &cefe::core::EntropicFieldEngine::step_forward, py::arg("dt"))
        .def("get_field_energy", &cefe::core::EntropicFieldEngine::get_field_energy)
        .def("get_state_size", &cefe::core::EntropicFieldEngine::get_state_size)
        .def("get_covariance_C", &cefe::core::EntropicFieldEngine::get_covariance_C)
        .def("get_covariance_P", &cefe::core::EntropicFieldEngine::get_covariance_P)
        .def("get_slice_indices", &cefe::core::EntropicFieldEngine::get_slice_indices);

    // Bind DynamicTensorEngine
    py::class_<cefe::dynamic_tensors::DynamicTensorEngine>(m, "DynamicTensorEngine")
        .def(py::init<std::shared_ptr<cefe::geometry::CausalDiamondGrid>>(), py::arg("grid"))
        .def("compute_stress_energy_tensor", &cefe::dynamic_tensors::DynamicTensorEngine::compute_stress_energy_tensor, py::arg("field_phi"), py::arg("mass"))
        .def("linearized_einstein_update", &cefe::dynamic_tensors::DynamicTensorEngine::linearized_einstein_update, py::arg("time_step"), py::arg("G_constant") = 1.0);

    // Bind QuantumWaveEngine
    py::class_<cefe::core::QuantumWaveEngine>(m, "QuantumWaveEngine")
        .def(py::init<std::shared_ptr<cefe::geometry::CausalDiamondGrid>>(), py::arg("grid"))
        .def("set_mass", &cefe::core::QuantumWaveEngine::set_mass, py::arg("m"))
        .def("set_hbar", &cefe::core::QuantumWaveEngine::set_hbar, py::arg("h"))
        .def("initialize_state", &cefe::core::QuantumWaveEngine::initialize_state, py::arg("target_t"))
        .def("set_gaussian_packet", &cefe::core::QuantumWaveEngine::set_gaussian_packet, py::arg("x0"), py::arg("y0"), py::arg("sigma"), py::arg("px"), py::arg("py"))
        .def("set_double_slit_potential", &cefe::core::QuantumWaveEngine::set_double_slit_potential, py::arg("slit_width"), py::arg("slit_separation"), py::arg("barrier_thickness"), py::arg("barrier_x"))
        .def("set_atomic_potential", &cefe::core::QuantumWaveEngine::set_atomic_potential, py::arg("Z"), py::arg("softening"))
        .def("build_hamiltonian", &cefe::core::QuantumWaveEngine::build_hamiltonian, py::arg("target_t"))
        .def("prepare_crank_nicolson", &cefe::core::QuantumWaveEngine::prepare_crank_nicolson, py::call_guard<py::gil_scoped_release>(), py::arg("dt"))
        .def("step_forward", &cefe::core::QuantumWaveEngine::step_forward, py::call_guard<py::gil_scoped_release>())
        .def("get_total_probability", &cefe::core::QuantumWaveEngine::get_total_probability)
        .def("get_probability_density", &cefe::core::QuantumWaveEngine::get_probability_density)
        .def("get_slice_indices", &cefe::core::QuantumWaveEngine::get_slice_indices);

    // Bind LatticeQFTEngine
    py::class_<cefe::core::LatticeQFTEngine>(m, "LatticeQFTEngine")
        .def(py::init<std::shared_ptr<cefe::geometry::CausalDiamondGrid>>())
        .def("set_mass", &cefe::core::LatticeQFTEngine::set_mass)
        .def("set_coupling", &cefe::core::LatticeQFTEngine::set_coupling)
        .def("initialize_state", &cefe::core::LatticeQFTEngine::initialize_state)
        .def("set_gaussian_packet", &cefe::core::LatticeQFTEngine::set_gaussian_packet)
        .def("set_opposing_packets", &cefe::core::LatticeQFTEngine::set_opposing_packets)
        .def("build_laplacian", &cefe::core::LatticeQFTEngine::build_laplacian)
        .def("step_forward", &cefe::core::LatticeQFTEngine::step_forward, py::call_guard<py::gil_scoped_release>())
        .def("get_field_amplitude", &cefe::core::LatticeQFTEngine::get_field_amplitude)
        .def("get_total_energy", &cefe::core::LatticeQFTEngine::get_total_energy);
}
