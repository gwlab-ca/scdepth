// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "bindings.hpp"

using namespace scdepth;
namespace py = pybind11;


PYBIND11_MODULE(_bindings, m) {
    m.doc() = "scdepth.bindings";
    scdepth::bind_barcode_counter(m);
    scdepth::bind_downsampler(m);
    scdepth::bind_vis_bins(m);
    scdepth::bind_govs(m);
}

py::array_t<uint32_t> scdepth::make_2d_u32(std::vector<uint32_t> &v,
                                         size_t dim0,        // e.g. fracs
                                         size_t dim1,        // e.g. barcodes
                                         py::handle owner)
{
    using ssize = py::ssize_t;
    if (v.size() != dim0 * dim1) {
        throw std::runtime_error("Size mismatch in make_2d_u32");
    }

    return py::array_t<uint32_t>(
        { static_cast<ssize>(dim0),
          static_cast<ssize>(dim1) },
        { static_cast<ssize>(dim1 * sizeof(uint32_t)),  // row stride
          static_cast<ssize>(sizeof(uint32_t)) },       // col stride
        v.data(),
        owner
    );
}

py::array_t<uint64_t> scdepth::make_2d_u64(std::vector<uint64_t> &v,
                                         size_t dim0,        // e.g. fracs
                                         size_t dim1,        // e.g. barcodes
                                         py::handle owner)
{
    using ssize = py::ssize_t;
    if (v.size() != dim0 * dim1) {
        throw std::runtime_error("Size mismatch in make_2d_u64");
    }

    return py::array_t<uint64_t>(
        { static_cast<ssize>(dim0),
          static_cast<ssize>(dim1) },
        { static_cast<ssize>(dim1 * sizeof(uint64_t)),  // row stride
          static_cast<ssize>(sizeof(uint64_t)) },       // col stride
        v.data(),
        owner
    );
}

py::array scdepth::make_3d_u64(std::vector<uint64_t> &v,
                            size_t dim0, size_t dim1, size_t dim2,
                            py::object owner) {
    return py::array(
        py::buffer_info(
            v.data(),
            sizeof(uint64_t),
            py::format_descriptor<uint64_t>::format(),
            3,
            { (ssize_t)dim0, (ssize_t)dim1, (ssize_t)dim2 },
            { (ssize_t)(dim1 * dim2 * sizeof(uint64_t)),
              (ssize_t)(dim2 * sizeof(uint64_t)),
              (ssize_t)sizeof(uint64_t) }
        ),
        owner
    );
}

py::array scdepth::make_3d_u32(std::vector<uint32_t> &v,
                            size_t dim0, size_t dim1, size_t dim2,
                            py::object owner) {
    return py::array(
        py::buffer_info(
            v.data(),
            sizeof(uint32_t),
            py::format_descriptor<uint32_t>::format(),
            3,
            { (ssize_t)dim0, (ssize_t)dim1, (ssize_t)dim2 },
            { (ssize_t)(dim1 * dim2 * sizeof(uint32_t)),
              (ssize_t)(dim2 * sizeof(uint32_t)),
              (ssize_t)sizeof(uint32_t) }
        ),
        owner
    );
}

py::object scdepth::sparse_to_csr(SparseMatrix &mat,
                                size_t n_rows,
                                size_t n_cols,
                                py::handle owner)
{
    namespace py = pybind11;
    using ssize = py::ssize_t;

    const ssize nnz = static_cast<ssize>(mat.indices.size());
    const ssize n_indptr = static_cast<ssize>(mat.indptr.size());

    if (mat.data.size() != static_cast<size_t>(nnz)) {
        throw std::runtime_error("SparseMatrix data/indices size mismatch");
    }

    // 1D numpy arrays (zero-copy, lifetime tied to owner)
    py::array_t<uint32_t> data(
        { nnz },
        { static_cast<ssize>(sizeof(uint32_t)) },
        mat.data.data(),
        owner
    );

    py::array_t<uint32_t> indices(
        { nnz },
        { static_cast<ssize>(sizeof(uint32_t)) },
        mat.indices.data(),
        owner
    );

    py::array_t<uint32_t> indptr(
        { n_indptr },
        { static_cast<ssize>(sizeof(uint32_t)) },
        mat.indptr.data(),
        owner
    );

    auto sparse = py::module_::import("scipy.sparse");

    py::tuple triplet(3);
    triplet[0] = data;
    triplet[1] = indices;
    triplet[2] = indptr;

    // shape = (n_rows, n_cols)
    py::tuple shape(2);
    shape[0] = py::int_(n_rows);
    shape[1] = py::int_(n_cols);

    py::object csr = sparse.attr("csr_matrix")(triplet, py::arg("shape") = shape);
    return csr;
}

py::array_t<uint32_t> scdepth::make_1d_u32(std::vector<uint32_t> &v,
                                         py::handle owner)
{
    using ssize = py::ssize_t;
    return py::array_t<uint32_t>(
        { static_cast<ssize>(v.size()) },
        { static_cast<ssize>(sizeof(uint32_t)) },
        v.data(),
        owner
    );
}

py::array_t<uint64_t> scdepth::make_1d_u64(std::vector<uint64_t> &v,
                                         py::handle owner)
{
    using ssize = py::ssize_t;
    return py::array_t<uint64_t>(
        { static_cast<ssize>(v.size()) },
        { static_cast<ssize>(sizeof(uint64_t)) },
        v.data(),
        owner
    );
}



