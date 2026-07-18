// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "bindings.hpp"
#include "counter.hpp"

using namespace scdepth;
namespace py = pybind11;

void scdepth::bind_barcode_counter(py::module_& m) {
    py::class_<BarcodeCounter>(
        m, "BarcodeCounter",
        "Caches read data from a scRNA-seq sam/bam/cram file")
        .def(py::init<>())

        .def("set_count_parameters", &BarcodeCounter::set_count_parameters,
             py::arg("min_gene") = 0.95, py::arg("min_gene_bases") = 40u,
             py::arg("min_exonic") = 0.95, py::arg("min_intronic") = 15u,
             py::arg("min_qual")              = 255,
             py::arg("discard_unknown_juncs") = false,
             py::arg("probes") = false, "Configure read filtering")

        .def("init", &BarcodeCounter::init, py::arg("lib_string"),
             py::arg("fwd"), py::arg("barcode_tag"), py::arg("barcode_re"),
             py::arg("umi_tag"), py::arg("sample_tag") = "",
             py::arg("samples")       = std::vector<std::string>{},
             py::arg("random_hex_re") = "", py::arg("random_hex_value") = "",
             py::arg("barcode_length") = 0u, py::arg("umi_length") = 0u,
             "Configure the scRNA-seq library properties")

        .def("prepare_bam", &BarcodeCounter::prepare_bam, py::arg("gtf"),
             py::arg("bam"), py::arg("out"), py::arg("threads") = 1,
             py::arg("max_tags")     = static_cast<size_t>(200000000),
             py::arg("max_tag_frac") = 0.95,
             py::call_guard<py::gil_scoped_release>{},
             "Prepare the bam for processing")

        .def("process_reads", &BarcodeCounter::process_reads, py::arg("chunk"),
             py::call_guard<py::gil_scoped_release>{},
             "Process up to `chunk` reads and return the number processed")

        .def("finish", &BarcodeCounter::finish,
             py::call_guard<py::gil_scoped_release>{},
             "Finalize and flush outputs")

        .def("done", &BarcodeCounter::done, "True when processing is complete")

        .def("total_reads", &BarcodeCounter::total_reads,
             "Total Reads Processed")

        .def("countable_reads", &BarcodeCounter::countable_reads,
             "Countable Reads Processed")

        .def("__repr__",
             [](const BarcodeCounter&) { return "<BarcodeCounter>"; });
}
