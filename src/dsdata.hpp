// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#pragma once

#include <cassert>
#include <cstdint>
#include <string>
#include <vector>

namespace scdepth {

struct SparseMatrix {
    SparseMatrix() { indptr.push_back(0); }
    void clear() {
        indptr.clear();
        indptr.push_back(0);
        indices.clear();
        data.clear();
    }

    bool empty() const { return indptr.size() < 2; }

    void merge(const SparseMatrix& rhs);
    // barcodes x genes
    std::vector<uint32_t> indptr;
    std::vector<uint32_t> indices;
    std::vector<uint32_t> data;
};

struct GeneCountMatrix {
    GeneCountMatrix() { indptr.push_back(0); }
    void clear() {
        indptr.clear();
        indptr.push_back(0);
        indices.clear();
    }

    bool empty() const { return indptr.size() < 2; }

    void                  merge(const GeneCountMatrix& rhs);
    std::vector<uint32_t> indptr;
    std::vector<uint32_t> indices;
};

template <typename T>
inline void sum_vectors(std::vector<T>& a, const std::vector<T>& b,
                        bool optional_rhs) {
    if(optional_rhs && b.empty())
        return;
    assert(a.size() == b.size());
    for(size_t i = 0; i < a.size(); i++) {
        a[i] += b[i];
    }
}

struct DownsampleResultsType {
    void reset(const std::vector<double>& fracs, size_t barcodes,
               size_t                 max_hist,
               /*size_t genes,*/ bool alloc_barcodes, bool has_mt, bool has_exc,
               size_t samples);

    void merge(const DownsampleResultsType& rhs);

    std::vector<uint32_t> bc_reads;     // fracs x barcodes
    std::vector<uint32_t> bc_mols;      // fracs x barcodes
    std::vector<uint32_t> bc_genes;     // fracs x barcodes
    std::vector<uint32_t> bc_mt_mols;   // convenience MT molecular counts
    std::vector<uint32_t> bc_mod_mols;  // convenience module molecular counts
    std::vector<uint32_t> mhist;        // reads per molecule
    std::vector<uint32_t> sample_mhist; // reads per molecule per sample

    std::vector<uint64_t> molecules;        // fracs
    std::vector<uint64_t> reads;            // fracs
    std::vector<uint64_t> sample_molecules; // samples x fracs
    std::vector<uint64_t> sample_reads;     // samples x fracs
    // std::vector<uint64_t> gene_molecules;   // fracs x genes
    // std::vector<uint64_t> gene_reads;       // fracs x genes
};

struct DownsampleResultsLocal {
    void reset(const std::vector<double>& fracs, size_t barcodes, size_t genes,
               size_t max_hist, bool build_mats, bool has_visium, bool calc_sau,
               bool aggregate_only, size_t samples);

    std::vector<double>   fracs;
    std::vector<uint64_t> reads_discarded;
    std::vector<uint64_t> reads_excluded;
    std::vector<uint64_t> mols_discarded;

    std::vector<GeneCountMatrix> gcounts;
    std::vector<uint64_t>        spliced_molecules; // fracs
    std::vector<uint64_t>        spliced_reads;     // fracs
    std::vector<uint32_t>        spliced_mhist;     // reads per molecule
    std::vector<uint32_t>        spliced_sample_mhist;

    std::vector<uint64_t> unspliced_molecules; // fracs
    std::vector<uint64_t> unspliced_reads;     // fracs
    std::vector<uint32_t> unspliced_mhist;     // reads per molecule
    std::vector<uint32_t> unspliced_sample_mhist;

    std::vector<uint64_t> ambiguous_molecules; // fracs
    std::vector<uint64_t> ambiguous_reads;     // fracs
    std::vector<uint32_t> ambiguous_mhist;     // reads per molecule
    std::vector<uint32_t> ambiguous_sample_mhist;

    std::vector<uint64_t> total_molecules; // fracs
    std::vector<uint64_t> total_reads;     // fracs
    std::vector<uint32_t> total_mhist;     // reads per molecule
    std::vector<uint32_t> total_sample_mhist;

    std::vector<SparseMatrix> total_mat;
    std::vector<SparseMatrix> spliced_mat;
    std::vector<SparseMatrix> ambiguous_mat;
    std::vector<SparseMatrix> unspliced_mat;

    std::vector<uint64_t> sample_reads_discarded;     // samples x fracs
    std::vector<uint64_t> sample_reads_excluded;      // samples x fracs
    std::vector<uint64_t> sample_mols_discarded;      // samples x fracs
    std::vector<uint64_t> sample_total_molecules;     // samples x fracs
    std::vector<uint64_t> sample_total_reads;         // samples x fracs
    std::vector<uint64_t> sample_spliced_molecules;   // samples x fracs
    std::vector<uint64_t> sample_spliced_reads;       // samples x fracs
    std::vector<uint64_t> sample_unspliced_molecules; // samples x fracs
    std::vector<uint64_t> sample_unspliced_reads;     // samples x fracs
    std::vector<uint64_t> sample_ambiguous_molecules; // samples x fracs
    std::vector<uint64_t> sample_ambiguous_reads;     // samples x fracs

    size_t   steps      = 0;
    size_t   barcodes   = 0;
    size_t   genes      = 0;
    size_t   samples    = 0;
    uint32_t max_hist   = 0;
    bool     build_mats = false;
    bool     has_visium = false;
    bool     has_mt     = false;
    bool     has_mod    = false;
    bool     has_exc    = false;
    bool     calc_sau   = false;
};

struct DownsampleResults {
    void reset(const std::vector<double>& fracs, size_t barcodes, size_t genes,
               size_t max_hist, bool build_mats, bool calc_sau,
               bool alloc_barcodes, bool has_visium, size_t samples);

    void merge(const DownsampleResultsLocal& rhs);

    /*
    size_t get_hidx_start(size_t f, size_t b) const{
        return (f * barcodes + b) * max_hist;
    }

    size_t get_hidx(size_t f, size_t b, size_t m) const{
        return get_hidx_start(f, b) + m;
    }

    void set_cell_mhist(size_t f, size_t b, const std::vector<uint32_t> &
    mhist){ size_t hs = get_hidx_start(f, b); assert(mhist.size() == max_hist);
        for(size_t i = 0; i < max_hist; i++){
            cell_mhist[hs + i] = mhist[i];
        }
    }
    */

    std::vector<double>   fracs;
    std::vector<uint32_t> bc_dis_reads; // (fracs x n_barcodes)
    std::vector<uint32_t> bc_dis_mols;  // (fracs x n_barcodes)
    std::vector<uint32_t> bc_exc_reads; // (fracs x n_barcodes)

    std::vector<uint64_t> reads_discarded;
    std::vector<uint64_t> reads_excluded;
    std::vector<uint64_t> mols_discarded;
    std::vector<uint64_t> sample_reads_discarded;
    std::vector<uint64_t> sample_reads_excluded;
    std::vector<uint64_t> sample_mols_discarded;

    std::vector<GeneCountMatrix> gcounts;
    std::string                  primer_mode;

    // std::vector<uint32_t>               cell_mhist;

    DownsampleResultsType spliced;   // not used if calc_sau
    DownsampleResultsType ambiguous; // not used if calc_sau
    DownsampleResultsType unspliced; // not used if total only
    DownsampleResultsType total;     // always populated

    std::vector<SparseMatrix> total_mat;     // optional count matrices
    std::vector<SparseMatrix> spliced_mat;   // optional count matrices
    std::vector<SparseMatrix> ambiguous_mat; // optional count matrices
    std::vector<SparseMatrix> unspliced_mat; // optional count matrices
    SparseMatrix              binned_mat;

    size_t   steps          = 0;
    size_t   barcodes       = 0;
    size_t   genes          = 0;
    size_t   samples        = 0;
    uint64_t seed           = 0;
    uint32_t max_hist       = 0;
    bool     build_mats     = false;
    bool     aggregate_only = false;
    bool     has_visium     = false;
    bool     has_mt         = false;
    bool     has_mod        = false;
    bool     has_exc        = false;
    bool     calc_sau       = false;
};

} // namespace scdepth
