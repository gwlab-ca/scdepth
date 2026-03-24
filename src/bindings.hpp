// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "dsdata.hpp"

namespace py = pybind11;
namespace scdepth{

py::object sparse_to_csr(SparseMatrix &mat, size_t n_rows, size_t n_cols, py::handle owner);

py::array_t<uint32_t> make_1d_u32(std::vector<uint32_t> &v, py::handle owner);
py::array_t<uint64_t> make_1d_u64(std::vector<uint64_t> &v, py::handle owner);

py::array_t<uint32_t> make_2d_u32(std::vector<uint32_t> &v, size_t dim0, size_t dim1, py::handle owner);
py::array_t<uint64_t> make_2d_u64(std::vector<uint64_t> &v, size_t dim0, size_t dim1, py::handle owner);

py::array make_3d_u64(std::vector<uint64_t> &v, size_t dim0, size_t dim1, size_t dim2, py::object owner);
py::array make_3d_u32(std::vector<uint32_t> &v, size_t dim0, size_t dim1, size_t dim2, py::object owner);


void bind_downsampler(pybind11::module_ &m);
void bind_barcode_counter(pybind11::module_ &m);
void bind_vis_bins(pybind11::module_& m);
void bind_govs(pybind11::module_& m);

};
