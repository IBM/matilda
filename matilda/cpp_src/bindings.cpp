#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/numpy.h"
#include "pybind11/buffer_info.h"
#include <pybind11/operators.h>
#include <pybind11/iostream.h>
#include "SparseMatrix.hpp"
#include "ModularInt.hpp"
#include "FilteredSimplicialComplex.hpp"
#include "HomologyComputer.hpp"
#include "RowVector.hpp"
#include "Real.hpp"

namespace py = pybind11;
using namespace matilda;

class PyFilteredSimplicialComplex : public FilteredSimplicialComplex
{
    using FilteredSimplicialComplex::construct_vietoris_from_metric;

public:
    void py_construct_vietoris_from_metric(py::array_t<float> numpy_array,
                                           int_fast64_t dimension,
                                           float diameter)
    {
        py::scoped_ostream_redirect stream(
            std::cout,                                // std::ostream&
            py::module_::import("sys").attr("stdout") // Python output
        );
        auto buffer = numpy_array.request();
        Matrix matrix((float *)buffer.ptr, buffer.shape[0], buffer.shape[1]);
        construct_vietoris_from_metric(matrix, dimension, diameter);
    }
};

template <typename Coefficient, typename Row>
class PyPersistentHomologyComputer : public PersistentHomologyComputer<Coefficient, Row>
{
    using PersistentHomologyComputer<Coefficient, Row>::compute_persistent_homology;
    using PersistentHomologyComputer<Coefficient, Row>::compute_persistent_homology_no_representatives;

public:
    std::map<int_fast64_t, std::map<int_fast64_t, std::map<int_fast64_t, int_fast64_t>>> py_persistent_cycles;

    std::map<int_fast64_t, std::map<int_fast64_t, float>> py_boundary_matrix;
    std::map<int_fast64_t, std::map<int_fast64_t, float>> py_reduced_boundary_matrix;

    void py_compute_persistent_homology(PyFilteredSimplicialComplex fsc,
                                        int_fast64_t top_degree,
                                        bool verbose,
                                        int_fast64_t modulus_param = 3)
    {
        constexpr int_fast64_t step_size = 10000;
        std::vector<int_fast64_t> endpoints;
        int_fast64_t i = 0;
        for (; i < fsc.simplices.size(); i += step_size)
        {
            endpoints.push_back(i);
        }
        if (i > fsc.simplices.size() - 1)
        {
            endpoints.push_back(fsc.simplices.size());
        }
        for (i = 1; i < endpoints.size(); ++i)
        {
            if (PyErr_CheckSignals() != 0)
                throw py::error_already_set();
            compute_persistent_homology(fsc, top_degree, verbose, true, endpoints[i - 1], endpoints[i], modulus_param);
        }

        // convert sparse boundary matrix into a pybind friendly dict
        for (auto &x : this->boundary_matrix.rows)
        {
            for (auto it = x.second.begin(); it != x.second.end(); ++it)
            {
                py_boundary_matrix[x.first][it->first] = it->second.value;
            }
        }

        // convert sparse reduced boundary matrix into a pybind friendly dict
        for (auto &x : this->store_matrix.rows)
        {
            for (auto it = x.second.begin(); it != x.second.end(); ++it)
            {
                py_reduced_boundary_matrix[x.first][it->first] = it->second.value;
            }
        }

        // convert persistent_cycles into a pybind friendly nested dict
        for (auto &x : this->persistent_cycles)
        {
            for (auto &y : x.second)
            {
                for (auto it = y.second.begin(); it != y.second.end(); ++it)
                {
                    py_persistent_cycles[x.first][y.first][it->first] = it->second.value;
                }
            }
        }
    }
    void py_compute_persistent_homology_no_representatives(PyFilteredSimplicialComplex fsc,
                                                           int_fast64_t top_degree,
                                                           bool verbose,
                                                           int_fast64_t modulus_param = 3)
    {
        constexpr int_fast64_t step_size = 10000;
        std::vector<int_fast64_t> endpoints;
        int_fast64_t i = 0;
        for (; i < fsc.simplices.size(); i += step_size)
        {
            endpoints.push_back(i);
        }
        if (i > fsc.simplices.size() - 1)
        {
            endpoints.push_back(fsc.simplices.size());
        }
        for (i = 1; i < endpoints.size(); ++i)
        {
            if (PyErr_CheckSignals() != 0)
                throw py::error_already_set();
            compute_persistent_homology_no_representatives(fsc, top_degree, verbose, true, endpoints[i - 1], endpoints[i], modulus_param);
        }
    }
    void py_print_persistent_cycles()
    {
        py::scoped_ostream_redirect stream(
            std::cout,
            py::module_::import("sys").attr("stdout"));
    }
};

PYBIND11_MODULE(matildacpp, m)
{
    m.doc() = "matilda C++ backend for homology computations";

    py::class_<PyFilteredSimplicialComplex>(m, "FilteredSimplicialComplex", py::dynamic_attr())
        .def(py::init<>())
        .def_readwrite("simplices", &PyFilteredSimplicialComplex::simplices)
        .def_readwrite("simplices_indices", &PyFilteredSimplicialComplex::simplices_indices)
        .def_readwrite("appears_at", &PyFilteredSimplicialComplex::appears_at)
        .def_readwrite("dimension", &PyFilteredSimplicialComplex::dimension)
        .def("construct_vietoris_from_metric",
             &PyFilteredSimplicialComplex::py_construct_vietoris_from_metric,
             py::arg("matrix"),
             py::arg("dimension"),
             py::arg("diameter"));

    py::class_<PyPersistentHomologyComputer<ModularInt, RowVector<ModularInt>>>(m, "PersistentHomologyComputerMod", py::dynamic_attr())
        .def(py::init<>())
        .def_readwrite("bars", &PyPersistentHomologyComputer<ModularInt, RowVector<ModularInt>>::bars)

        .def_readwrite("proto_bars", &PyPersistentHomologyComputer<ModularInt, RowVector<ModularInt>>::proto_bars)

        .def_readwrite("persistent_cycles", &PyPersistentHomologyComputer<ModularInt, RowVector<ModularInt>>::py_persistent_cycles)

        // naive attempt
        .def_readwrite("boundary_matrix", &PyPersistentHomologyComputer<ModularInt, RowVector<ModularInt>>::py_boundary_matrix)
        .def_readwrite("reduced_boundary_matrix", &PyPersistentHomologyComputer<ModularInt, RowVector<ModularInt>>::py_reduced_boundary_matrix)

        .def("compute_persistent_homology",
             &PyPersistentHomologyComputer<ModularInt, RowVector<ModularInt>>::py_compute_persistent_homology,
             py::arg("fsc"),
             py::arg("top_degree"),
             py::arg("verbose"),
             py::arg("modulus_param"))
        .def("compute_persistent_homology_no_representatives",
             &PyPersistentHomologyComputer<ModularInt, RowVector<ModularInt>>::py_compute_persistent_homology_no_representatives,
             py::arg("fsc"),
             py::arg("top_degree"),
             py::arg("verbose"),
             py::arg("modulus_param"));

    py::class_<PyPersistentHomologyComputer<Real, RowVector<Real>>>(m, "PersistentHomologyComputerReal", py::dynamic_attr())
        .def(py::init<>())
        .def_readwrite("bars", &PyPersistentHomologyComputer<Real, RowVector<Real>>::bars)

        .def_readwrite("proto_bars", &PyPersistentHomologyComputer<Real, RowVector<Real>>::proto_bars)

        .def_readwrite("persistent_cycles", &PyPersistentHomologyComputer<Real, RowVector<Real>>::py_persistent_cycles)

        // naive attempt
        .def_readwrite("boundary_matrix", &PyPersistentHomologyComputer<Real, RowVector<Real>>::py_boundary_matrix)
        .def_readwrite("reduced_boundary_matrix", &PyPersistentHomologyComputer<Real, RowVector<Real>>::py_reduced_boundary_matrix)

        .def("compute_persistent_homology",
             &PyPersistentHomologyComputer<Real, RowVector<Real>>::py_compute_persistent_homology,
             py::arg("fsc"),
             py::arg("top_degree"),
             py::arg("verbose"),
             py::arg("modulus_param"))
        .def("compute_persistent_homology_no_representatives",
             &PyPersistentHomologyComputer<Real, RowVector<Real>>::py_compute_persistent_homology_no_representatives,
             py::arg("fsc"),
             py::arg("top_degree"),
             py::arg("verbose"),
             py::arg("modulus_param"));

    py::class_<Matrix>(m, "Matrix")
        .def(py::init([](py::array_t<float> numpy_matrix,
                         uint_fast64_t n_rows,
                         uint_fast64_t n_columns)
                      { return std::unique_ptr<Matrix>(new Matrix((float *)numpy_matrix.request().ptr,
                                                                  n_rows,
                                                                  n_columns)); }));
}
