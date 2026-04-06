// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "bindings.hpp"
#include "visbins.hpp"
#include "downsampler.hpp"

using namespace scdepth;
namespace py = pybind11;

py::dict barcodes_to_cols(const Downsampler& ds){
    auto& v = ds.barcodes;
    const py::ssize_t n = static_cast<py::ssize_t>(v.size());

    py::array_t<uint64_t> offset(n);
    py::array_t<uint64_t> seed(n);
    py::array_t<uint32_t> index(n);
    py::array_t<uint32_t> total(n);
    py::array_t<uint32_t> countable(n);
    py::array_t<uint32_t> raw_molecules(n);

    auto off = offset.mutable_unchecked<1>();
    auto sed = seed.mutable_unchecked<1>();
    auto idx = index.mutable_unchecked<1>();
    auto tot = total.mutable_unchecked<1>();
    auto cnt = countable.mutable_unchecked<1>();
    auto raw = raw_molecules.mutable_unchecked<1>();

    py::array_t<py::object> barcode_obj(n);

    auto bptr = static_cast<py::object*>(barcode_obj.mutable_data());

    for(py::ssize_t i = 0; i < n; ++i){
        const auto& bc = v[(size_t)i];

        // construct Python string and store into object array
        bptr[i] = py::str(bc.barcode);

        off(i) = bc.offset;
        sed(i) = bc.seed;
        idx(i) = bc.index;
        tot(i) = bc.total;
        cnt(i) = bc.countable;
        raw(i) = bc.raw_molecules;
    }

    py::dict d;
    d["barcode"] = std::move(barcode_obj);
    d["offset"] = std::move(offset);
    d["seed"] = std::move(seed);
    d["index"] = std::move(index);
    d["total"] = std::move(total);
    d["countable"] = std::move(countable);
    d["raw_molecules"] = std::move(raw_molecules);
    return d;
}

void scdepth::bind_downsampler(py::module_ &m) {

    py::class_<Downsampler, std::unique_ptr<Downsampler>>(m, "Downsampler")
        .def(py::init<>())

        .def("init",
             &Downsampler::init,
             py::arg("prefix"),
             py::arg("mt_prefix") = "",
             py::arg("mt_file") = "",
             py::arg("mod_file") = "",
             py::arg("exclude_file") = "",

             py::arg("max_hist") = 50,
             py::arg("build_matrices") = false,
             py::arg("calc_sau") = false,
             py::call_guard<py::gil_scoped_release>{},
             "Initialize the downsampler from a cached prefix")

        .def("init_visium",
             &Downsampler::init_visium,
             py::arg("rows"),
             py::arg("cols"),
             py::arg("in_tissue"),
             py::arg("countable"),
             py::arg("total"),
             py::arg("total_rows"),
             py::arg("total_cols"),
             py::call_guard<py::gil_scoped_release>{},
             "Initialize the downsampler to create binned tiles/barcodes for visium HD data")

        .def("downsample",
             &Downsampler::downsample,
             py::arg("fracs"),
             py::arg("umi_len"),
             py::arg("seed"),
             py::arg("threads") = 1,
             py::arg("aggregate_only") = false,
             py::arg("umi_mode") = "directed",
             py::arg("correct_multi_umis") = true,
             py::call_guard<py::gil_scoped_release>{},
             "Run downsampling for the given barcodes (optional) and fractions")

        .def("reset_visium", &Downsampler::reset_visium, 
             py::call_guard<py::gil_scoped_release>{},
             "Reset visium tile aggregation information")

        .def("clear_output", &Downsampler::clear_output, 
             py::call_guard<py::gil_scoped_release>{},
             "Clear output memory")

        .def("write_gene_mats", &Downsampler::write_gene_mats, 
             py::arg("output"),
             py::arg("idx"),
             py::arg("bin_div") = 0,
             py::call_guard<py::gil_scoped_release>{},
             "Write the gene detection matrices to a file for comparisons (optionally bin for visium HD data)")

        .def("write_gene_baseline", &Downsampler::write_gene_baseline, 
             py::arg("output"),
             py::arg("idx"),
             py::arg("step") = 0,
             py::arg("bin_div") = 0,
             py::call_guard<py::gil_scoped_release>{},
             "Write the gene detection full count matrix to a file for comparisons (optionally bin for visium HD data)")

        .def_property_readonly("fracs",
                [](Downsampler &self) {
                    // std::vector<double> -> Python list[float]
                    return self.output.fracs;
                },
                "Downsampling fractions actually used")

        //Scalars --------------------------------------------------------------------------- 
        .def_property_readonly("total_fracs",
            [](Downsampler &self) {
                return self.output.steps;
            },
            "Number of downsampling fractions")

        .def_property_readonly("total_barcodes",
            [](Downsampler &self) {
                return self.output.barcodes;
            },
            "Number of raw barcodes in the results")

        .def_property_readonly("total_genes",
            [](Downsampler &self) {
                return self.output.genes;
            },
            "Number of genes in the results")

        .def_property_readonly("max_hist",
            [](Downsampler &self) {
                return self.output.max_hist;
            },
            "Maximum histogram bin (histograms have size max_hist+1)")

        .def_property_readonly("seed",
            [](Downsampler &self) {
                return self.output.seed;
            },
            "Downsampling seed")

        .def_property_readonly("build_mats",
            [](Downsampler &self) {
                return self.output.build_mats;
            },
            "Whether sparse matrices were built")

        .def_property_readonly("has_MT",
            [](Downsampler &self) {
                return self.output.has_mt;
            },
            "A MT gene list was included")

        .def_property_readonly("has_module",
            [](Downsampler &self) {
                return self.output.has_mod;
            },
            "A module gene list was included")

        .def_property_readonly("has_excluded",
            [](Downsampler &self) {
                return self.output.has_exc;
            },
            "An excluded gene list was included")

        .def_property_readonly("has_sau",
            [](Downsampler &self) {
                return self.output.calc_sau;
            },
            "Has spliced/ambig/unspliced quantifications")

        .def_property_readonly("aggregate_only",
            [](Downsampler &self) {
                return self.output.aggregate_only;
            },
            "Whether sparse matrices were built")

        .def_property_readonly("barcodes",
            [](Downsampler& self)->py::dict{
                return barcodes_to_cols(self);
            },
            "Barcode columns as a dict-of-arrays (DataFrame-friendly)")

        //aggregate counts  --------------------------------------------------------------------------- 
        .def_property_readonly("spliced_reads",
            [](Downsampler &self) {
                auto &res = self.output;
                return make_1d_u64(res.spliced.reads, py::cast(&self));
            },
            "Number of spliced reads")

        .def_property_readonly("unspliced_reads",
            [](Downsampler &self) {
                auto &res = self.output;
                return make_1d_u64(res.unspliced.reads, py::cast(&self));
            },
            "Number of unspliced reads")

        .def_property_readonly("ambiguous_reads",
            [](Downsampler &self) {
                auto &res = self.output;
                return make_1d_u64(res.ambiguous.reads, py::cast(&self));
            },
            "Number of ambiguous reads")

        .def_property_readonly("total_reads",
            [](Downsampler &self) {
                auto &res = self.output;
                return make_1d_u64(res.total.reads, py::cast(&self));
            },
            "Number of total reads")

        .def_property_readonly("spliced_molecules",
            [](Downsampler &self) {
                auto &res = self.output;
                return make_1d_u64(res.spliced.molecules, py::cast(&self));
            },
            "Number of spliced molecules")

        .def_property_readonly("unspliced_molecules",
            [](Downsampler &self) {
                auto &res = self.output;
                return make_1d_u64(res.unspliced.molecules, py::cast(&self));
            },
            "Number of unspliced molecules")

        .def_property_readonly("ambiguous_molecules",
            [](Downsampler &self) {
                auto &res = self.output;
                return make_1d_u64(res.ambiguous.molecules, py::cast(&self));
            },
            "Number of ambiguous molecules")

        .def_property_readonly("total_molecules",
            [](Downsampler &self) {
                auto &res = self.output;
                return make_1d_u64(res.total.molecules, py::cast(&self));
            },
            "Number of total molecules")


        //UMI correction / ambig umi filtering / exclude filtering -----------------------------------------------
        .def_property_readonly("reads_discarded",
            [](Downsampler &self) {
                auto &res = self.output;
                return make_1d_u64(res.reads_discarded, py::cast(&self));
            },
            "Number of reads lost to ambiguous UMI/gene mappings")

        .def_property_readonly("molecules_discarded",
            [](Downsampler &self) {
                auto &res = self.output;
                return make_1d_u64(res.mols_discarded, py::cast(&self));
            },
            "Number of molecules lost to ambiguous UMI/gene mappings")

        .def_property_readonly("reads_excluded",
            [](Downsampler &self) {
                auto &res = self.output;
                return make_1d_u64(res.reads_excluded, py::cast(&self));
            },
            "Number of reads lost to excluded gene filter")

        .def_property_readonly("molecules_ambig",
            [](Downsampler &self) {
                auto &res = self.output;
                return make_1d_u64(res.mols_ambig, py::cast(&self));
            },
            "Number of molecule subgraphs mapping to at least 2 genes.")

        //barcode count matrices  --------------------------------------------------------------------------- 
        .def_property_readonly("total_bc_reads",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.total.bc_reads,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Total reads per barcode shape (steps x barcodes)")

        .def_property_readonly("total_bc_mols",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.total.bc_mols,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Total molecules per barcode shape (steps x barcodes)")

        .def_property_readonly("total_bc_genes",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.total.bc_genes,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Total genes per barcode shape (steps x barcodes)")

        .def_property_readonly("total_bc_mt",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.total.bc_mt_mols,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Total MT mols per barcode shape (steps x barcodes)")

        .def_property_readonly("total_bc_mod",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.total.bc_mod_mols,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Total module mols per barcode shape (steps x barcodes)")

        .def_property_readonly("spliced_bc_reads",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.spliced.bc_reads,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Spliced reads per barcode shape (steps x barcodes)")

        .def_property_readonly("spliced_bc_mols",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.spliced.bc_mols,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Spliced molecules per barcode shape (steps x barcodes)")

        .def_property_readonly("spliced_bc_genes",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.spliced.bc_genes,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Spliced genes per barcode shape (steps x barcodes)")

        .def_property_readonly("spliced_bc_mt",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.spliced.bc_mt_mols,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Spliced MT mols per barcode shape (steps x barcodes)")

        .def_property_readonly("spliced_bc_mod",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.spliced.bc_mod_mols,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Spliced module mols per barcode shape (steps x barcodes)")

        .def_property_readonly("unspliced_bc_reads",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.unspliced.bc_reads,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Unspliced reads per barcode shape (steps x barcodes)")

        .def_property_readonly("unspliced_bc_mols",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.unspliced.bc_mols,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Unspliced molecules per barcode shape (steps x barcodes)")

        .def_property_readonly("unspliced_bc_genes",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.unspliced.bc_genes,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Unspliced genes per barcode shape (steps x barcodes)")

        .def_property_readonly("unspliced_bc_mt",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.unspliced.bc_mt_mols,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Unspliced MT mols per barcode shape (steps x barcodes)")

        .def_property_readonly("unspliced_bc_mod",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.unspliced.bc_mod_mols,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Unspliced module mols per barcode shape (steps x barcodes)")

        .def_property_readonly("ambiguous_bc_reads",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.ambiguous.bc_reads,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Ambiguous reads per barcode shape (steps x barcodes)")

        .def_property_readonly("ambiguous_bc_mols",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.ambiguous.bc_mols,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Ambiguous molecules per barcode shape (steps x barcodes)")

        .def_property_readonly("ambiguous_bc_genes",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.ambiguous.bc_genes,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Ambiguous genes per barcode shape (steps x barcodes)")

        .def_property_readonly("ambiguous_bc_mt",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.ambiguous.bc_mt_mols,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Ambiguous MT mols per barcode shape (steps x barcodes)")

        .def_property_readonly("ambiguous_bc_mod",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.ambiguous.bc_mod_mols,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Ambiguous module mols per barcode shape (steps x barcodes)")

        //barcode umi correction metrics  --------------------------------------------------------------------------- 
        .def_property_readonly("bc_discarded_reads",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.bc_dis_reads,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Reads lost from multi-gene UMIs and excluded gene lists (steps x barcodes)")

        .def_property_readonly("bc_excluded_reads",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.bc_exc_reads,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Reads lost from excluded read lists (steps x barcodes)")

        .def_property_readonly("bc_discarded_moleculess",
            [](Downsampler &self) -> py::object {
                if (self.output.aggregate_only) {
                    return py::none();
                }
                auto &res = self.output;
                return make_2d_u32(res.bc_dis_mols,
                                res.steps,
                                res.barcodes,
                                py::cast(&self));
            },
            "Molecules lost from multi-gene UMIs (steps x barcodes)")

        /*
        //aggregate reads per molecule histograms  --------------------------------------------------------------------------- 
        .def_property_readonly("barcode_mhist",
            [](Downsampler &self) {
                auto &res = self.output;
                return make_3d_u32(res.cell_mhist,
                                res.steps,
                                res.barcodes,
                                res.max_hist,
                                py::cast(&self));
            },
            "Histogram of per barcode total reads per molecule (steps x barcodes x max_hist)")
        */
        .def_property_readonly("total_mhist",
            [](Downsampler &self) {
                auto &res = self.output;
                return make_2d_u32(res.total.mhist,
                                res.steps,
                                res.max_hist,
                                py::cast(&self));
            },
            "Histogram of total reads per molecule (steps x max_hist)")

        .def_property_readonly("spliced_mhist",
            [](Downsampler &self) {
                auto &res = self.output;
                return make_2d_u32(res.spliced.mhist,
                                res.steps,
                                res.max_hist,
                                py::cast(&self));
            },
            "Histogram of spliced reads per molecule (steps x max_hist)")

        .def_property_readonly("unspliced_mhist",
            [](Downsampler &self) {
                auto &res = self.output;
                return make_2d_u32(res.unspliced.mhist,
                                res.steps,
                                res.max_hist,
                                py::cast(&self));
            },
            "Histogram of unspliced reads per molecule (steps x max_hist)")

        .def_property_readonly("ambiguous_mhist",
            [](Downsampler &self) {
                auto &res = self.output;
                return make_2d_u32(res.ambiguous.mhist,
                                res.steps,
                                res.max_hist,
                                py::cast(&self));
            },
            "Histogram of ambiguous reads per MT molecule (steps x max_hist)")

        //sparse matrices  --------------------------------------------------------------------------- 
        .def("total_csr",
            [](Downsampler &self, size_t step) {
                auto &res = self.output;
                if (!res.build_mats)
                    throw std::runtime_error("build_mats was false");
                if (step >= res.steps)
                    throw std::out_of_range("step out of range");
                py::print(
                    "Test",
                    "indptr", res.total_mat[step].indptr.size(),
                    "indices", res.total_mat[step].indices.size(),
                    "data", res.total_mat[step].data.size(),
                    "barcodes", res.barcodes,
                    "genes", res.genes
                );
                return sparse_to_csr(res.total_mat[step],
                                    res.barcodes,
                                    res.genes,
                                    py::cast(&self));
            },
            py::arg("step"))

        .def("total_csr_bin",
            [](Downsampler &self, size_t step, uint32_t bin_div) {
                auto &res = self.output;
                if(!res.build_mats)
                    throw std::runtime_error("build_mats was false");
                if(step >= res.steps)
                    throw std::out_of_range("step out of range");
                if(!res.has_visium)
                    throw std::out_of_range("must be visium HD data");
                bin_gene_counts(self, step, bin_div, res.binned_mat);
                return sparse_to_csr(res.binned_mat,
                                    res.binned_mat.indptr.size() - 1,
                                    res.genes,
                                    py::cast(&self), true);
            },
            py::arg("step"),
            py::arg("bin_div"))

        .def("spliced_csr",
            [](Downsampler &self, size_t step) {
                auto &res = self.output;
                if (!res.build_mats)
                    throw std::runtime_error("build_mats was false");
                if (step >= res.steps)
                    throw std::out_of_range("step out of range");
                return sparse_to_csr(res.spliced_mat[step],
                                    res.barcodes,
                                    res.genes,
                                    py::cast(&self));
            },
            py::arg("step"))

        .def("ambiguous_csr",
            [](Downsampler &self, size_t step) {
                auto &res = self.output;
                if (!res.build_mats)
                    throw std::runtime_error("build_mats was false");
                if (step >= res.steps)
                    throw std::out_of_range("step out of range");
                return sparse_to_csr(res.ambiguous_mat[step],
                                    res.barcodes,
                                    res.genes,
                                    py::cast(&self));
            },
            py::arg("step"))

        .def("unspliced_csr",
            [](Downsampler &self, size_t step) {
                auto &res = self.output;
                if (!res.build_mats)
                    throw std::runtime_error("build_mats was false");
                if (step >= res.steps)
                    throw std::out_of_range("step out of range");
                return sparse_to_csr(res.unspliced_mat[step],
                                    res.barcodes,
                                    res.genes,
                                    py::cast(&self));
            },
            py::arg("step"))

        ;
}
