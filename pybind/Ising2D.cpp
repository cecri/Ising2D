#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Ising2D.hpp"

namespace py=pybind11;

PYBIND11_MODULE(Ising2D, m) {
	py::class_<Ising2D>(m, "Ising2D")
		.def(py::init<uint32_t, uint32_t>())
		.def_property_readonly("nrows", &Ising2D::get_nrows)
		.def_property_readonly("ncols", &Ising2D::get_ncols)
		.def_property_readonly("size", &Ising2D::size)
		.def("energy", &Ising2D::energy)
		.def("to_idx", &Ising2D::to_idx)
		.def("to_coord", &Ising2D::to_coord)
		.def("neighbors", &Ising2D::neighbors)
		.def("magnetization", &Ising2D::magnetization)
		.def("all_neighbors", &Ising2D::all_neighbors);

	py::class_<WolffSampler>(m, "WolffSampler")
		.def(py::init<const Ising2D&, double>())
		.def("set_seed", &WolffSampler::set_seed)
		.def("randomize_conf", &WolffSampler::randomize_conf)
		.def_property_readonly("conf", &WolffSampler::get_conf)
		.def("sweep", &WolffSampler::sweep);
}
