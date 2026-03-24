// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#include "bindings.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>
#include "visbins.hpp"

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

void scdepth::bind_vis_bins(py::module_ &m) {
    m.def(
        "aggregate_visium_bins",
        [](const Downsampler &ds,
           uint32_t step,
           uint32_t row_div,
           uint32_t col_div) -> py::tuple {

            VisiumBin out;
            {
                py::gil_scoped_release release;
                out = scdepth::aggregate_visium_bins(ds, step, row_div, col_div);
            }

            py::dict d, b;
            if(out.invalid) py::make_tuple(b, d);

            d["invalid"] = out.invalid;
            d["n_rows"] = py::int_(out.n_rows);
            d["n_cols"] = py::int_(out.n_cols);
            d["n_bins"] = py::int_(out.n_bins);
            d["step"] = step;
            d["row_div"] = row_div;
            d["col_div"] = col_div;
            d["row_total"] = out.total_rows;
            d["col_total"] = out.total_cols;


            // 1D arrays length n_bins
            b["reads"] = move_vec_1d<uint32_t>(std::move(out.reads));
            b["molecules"] = move_vec_1d<uint32_t>(std::move(out.mols));
            b["genes"] = move_vec_1d<uint32_t>(std::move(out.genes));
            if(!out.mt_mols.empty()) b["mt_molecules"] = move_vec_1d<uint32_t>(std::move(out.mt_mols));
            if(!out.mod_mols.empty()) b["mod_molecules"] = move_vec_1d<uint32_t>(std::move(out.mod_mols));
            b["in_tissue"] = move_vec_1d<uint32_t>(std::move(out.in_tissue));
            b["countable"] = move_vec_1d<uint32_t>(std::move(out.countable));
            b["total"] = move_vec_1d<uint32_t>(std::move(out.total));
            b["n_barcodes"] = move_vec_1d<uint32_t>(std::move(out.n_barcodes));

            // Convenience for pandas (also 1D, length n_bins)
            b["array_row"] = move_vec_1d<uint32_t>(std::move(out.bin_rows));
            b["array_col"] = move_vec_1d<uint32_t>(std::move(out.bin_cols));

            return py::make_tuple(b, d);
        },
        py::arg("downsampler"),
        py::arg("step"),
        py::arg("row_div"),
        py::arg("col_div"),
        "Aggregate per-barcode outputs into Visium HD bins. Returns a dict of 1d numpy arrays."
    );

    /*
    py::class_<RowColHists>(m, "RowColHists").def(py::init<>());
    py::class_<RowColSums>(m, "RowColSums").def(py::init<>());

    m.def("aggregate_visium_hists",
        [](const Downsampler& ds, RowColHists& out, uint32_t step,
            unsigned int min_molecules, bool in_tissue) -> py::dict
        {
            {
                py::gil_scoped_release release;
                scdepth::aggregate_visium_hists(ds, out, step, min_molecules, in_tissue);
            }

            const ssize_t H = (ssize_t)out.max_hist - 1;

            py::dict d;
            py::object base = py::cast(&out, py::return_value_policy::reference);

            d["rows"] = make_2d_u32(out.rows, (ssize_t)out.total_rows, H, base);
            d["cols"] = make_2d_u32(out.cols, (ssize_t)out.total_cols, H, base);

            d["row_reads"] = make_1d_u32(out.row_reads, base);
            d["row_molecules"] = make_1d_u32(out.row_molecules, base);
            d["row_support"] = make_1d_u32(out.row_support, base);
            d["row_countable"] = make_1d_u32(out.row_countable, base);
            d["row_total"] = make_1d_u32(out.row_total, base);

            d["col_reads"] = make_1d_u32(out.col_reads, base);
            d["col_molecules"] = make_1d_u32(out.col_molecules, base);
            d["col_support"] = make_1d_u32(out.col_support, base);
            d["col_countable"] = make_1d_u32(out.col_countable, base);
            d["col_total"] = make_1d_u32(out.col_total, base);

            d["row_total"] = out.total_rows;
            d["col_total"] = out.total_cols;
            d["max_hist"] = out.max_hist;

            if(ds.output.has_mt){
                d["col_mt_molecules"] = make_1d_u32(out.col_mt_molecules, base);
                d["row_mt_molecules"] = make_1d_u32(out.row_mt_molecules, base);
            }

            return d;
        },
        py::arg("downsampler"),
        py::arg("out"),
        py::arg("step"),
        py::arg("min_molecules"),
        py::arg("in_tissue"),
        //py::arg("permute") = false,
        //py::arg("frac") = 0.0,
        //py::arg("seed") = 0,
        "Aggregate barcode molecule count hists by row and col.  The last last histogram bin(ie the tail) is not included."
    );

    m.def("aggregate_visium",
        [](const Downsampler& ds, RowColSums& out, uint32_t step,
            unsigned int min_molecules, bool in_tissue, bool permute, double frac, uint32_t seed) -> py::dict
        {
            {
                py::gil_scoped_release release;
                scdepth::aggregate_visium(ds, out, step, min_molecules, in_tissue, permute, frac, seed);
            }

            py::dict d;
            py::object base = py::cast(&out, py::return_value_policy::reference);

            d["row_reads"] = make_1d_u32(out.row_reads, base);
            d["row_molecules"] = make_1d_u32(out.row_molecules, base);
            d["row_support"] = make_1d_u32(out.row_support, base);
            d["row_lost"] = make_1d_u32(out.row_lost, base);

            d["col_reads"] = make_1d_u32(out.col_reads, base);
            d["col_molecules"] = make_1d_u32(out.col_molecules, base);
            d["col_support"] = make_1d_u32(out.col_support, base);
            d["col_lost"] = make_1d_u32(out.col_lost, base);

            d["row_total"] = out.total_rows;
            d["col_total"] = out.total_cols;
            d["max_hist"] = out.max_hist;

            return d;
        },
        py::arg("downsampler"),
        py::arg("out"),
        py::arg("step"),
        py::arg("min_molecules"),
        py::arg("in_tissue"),
        py::arg("permute") = false,
        py::arg("frac") = 0.0,
        py::arg("seed") = 0,
        "Aggregate barcode molecule counts by row and col with optional permuting or bootstrapping."
    );
    */

}
