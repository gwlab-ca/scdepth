// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#include "bindings.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>
#include <string>
#include <utility>
#include <pybind11/stl.h>
#include "goverlaps.hpp"

namespace py = pybind11;
using namespace scdepth;


template <typename T>
py::array move_vec_1d(std::vector<T>&& v) {
    auto *vec = new std::vector<T>(std::move(v));

    py::capsule owner(vec, [](void *p) {
        delete reinterpret_cast<std::vector<T>*>(p);
    });

    return py::array(
        py::dtype::of<T>(),
        { static_cast<py::ssize_t>(vec->size()) },
        { static_cast<py::ssize_t>(sizeof(T)) },
        vec->data(),
        owner
    );
}

void scdepth::bind_govs(py::module_ &m) {

    m.def(
        "calculate_gene_overlaps", [](const std::string& baseline_file,
                                            const std::vector<std::string>& files) -> py::dict {
        GeneOverlap go;
        if(!go.init(baseline_file, files)){
            throw std::runtime_error("GeneOverlap::init failed");
        }
        {
            py::gil_scoped_release release;
            if(!go.calculate()){
                throw std::runtime_error("GeneOverlap::calculate failed");
            }
        }

        py::dict out;

        out["counts"]   = move_vec_1d<uint32_t>(std::move(go.counts));
        out["counts_1"] = move_vec_1d<uint32_t>(std::move(go.counts_1));
        out["counts_2"] = move_vec_1d<uint32_t>(std::move(go.counts_2));
        out["counts_3p"] = move_vec_1d<uint32_t>(std::move(go.counts_3p));

        for (uint32_t i = 0; i < go.steps; ++i) {
            auto& m = go.mean_outs.at(i);

            out[py::str("stability_{}").format(i)]     = move_vec_1d<double>(std::move(m.stability));
            out[py::str("stability_1_{}").format(i)]   = move_vec_1d<double>(std::move(m.stability_1));
            out[py::str("stability_2_{}").format(i)]   = move_vec_1d<double>(std::move(m.stability_2));
            out[py::str("stability_3p_{}").format(i)]   = move_vec_1d<double>(std::move(m.stability_3p));
        }

        out["barcodes"] = py::int_(go.barcodes);
        out["genes"]    = py::int_(go.genes);
        out["steps"]    = py::int_(go.steps);
        out["seeds"]    = py::int_(go.seeds);

        return out;
    });
}
