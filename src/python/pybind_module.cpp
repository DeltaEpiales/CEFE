#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <memory>

#include "cefe/geometry/CausalDiamond.hpp"
#include "cefe/geometry/Metric.hpp"
#include "cefe/core/EntropicFieldEngine.hpp"
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
        .def("get_state_size", &cefe::core::EntropicFieldEngine::get_state_size);

    // Bind DynamicTensorEngine
    py::class_<cefe::dynamic_tensors::DynamicTensorEngine>(m, "DynamicTensorEngine")
        .def(py::init<std::shared_ptr<cefe::geometry::CausalDiamondGrid>>(), py::arg("grid"))
        .def("compute_stress_energy_tensor", &cefe::dynamic_tensors::DynamicTensorEngine::compute_stress_energy_tensor, py::arg("field_phi"), py::arg("mass"))
        .def("linearized_einstein_update", &cefe::dynamic_tensors::DynamicTensorEngine::linearized_einstein_update, py::arg("time_step"), py::arg("G_constant") = 1.0);
}
