// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#pragma once

#include "downsampler.hpp"
#include <cassert>
#include <cstdint>
#include <vector>

namespace scdepth {

struct VisiumBin {
    std::vector<uint32_t> reads;
    std::vector<uint32_t> mols;
    std::vector<uint32_t> genes;
    std::vector<uint32_t> mt_mols;
    std::vector<uint32_t> mod_mols;
    std::vector<uint32_t> in_tissue;
    std::vector<uint32_t> countable;
    std::vector<uint32_t> total;
    std::vector<uint32_t> n_barcodes;

    std::vector<uint32_t> bin_rows;
    std::vector<uint32_t> bin_cols;

    size_t n_rows;
    size_t n_cols;
    size_t total_rows;
    size_t total_cols;
    size_t n_bins;
    size_t max_hist;
    bool   invalid = false;
};

/*
struct RowColHists{
    void resize(uint32_t barcodes, uint32_t total_rows, uint32_t total_cols,
uint32_t max_hist, bool has_mt){ this->total_rows = total_rows; this->total_cols
= total_cols; this->max_hist = max_hist; this->barcodes = barcodes;

        const size_t H = max_hist - 1;
        rows.assign(size_t(total_rows) * H, 0);
        cols.assign(size_t(total_cols) * H, 0);
        row_molecules.assign(total_rows, 0);
        row_reads.assign(total_rows, 0);
        row_support.assign(total_rows, 0);
        row_total.assign(total_rows, 0);
        row_countable.assign(total_rows, 0);

        col_molecules.assign(total_cols, 0);
        col_total.assign(total_cols, 0);
        col_countable.assign(total_cols, 0);
        col_reads.assign(total_cols, 0);
        col_support.assign(total_cols, 0);
        seed = 0; pidx.clear();
        if(has_mt){
            col_mt_molecules.assign(total_cols, 0);
            row_mt_molecules.assign(total_cols, 0);
        }
    }

    uint32_t              step;
    uint32_t              barcodes;
    uint32_t              total_rows;
    uint32_t              total_cols;
    uint32_t              max_hist;
    uint32_t              seed;
    std::vector<uint32_t> pidx;
    std::vector<uint32_t> rows;
    std::vector<uint32_t> cols;
    std::vector<uint32_t> row_reads;
    std::vector<uint32_t> row_molecules;
    std::vector<uint32_t> row_mt_molecules;
    std::vector<uint32_t> row_total;
    std::vector<uint32_t> row_countable;
    std::vector<uint32_t> row_support;
    std::vector<uint32_t> col_reads;
    std::vector<uint32_t> col_mt_molecules;
    std::vector<uint32_t> col_molecules;
    std::vector<uint32_t> col_support;
    std::vector<uint32_t> col_total;
    std::vector<uint32_t> col_countable;
};

struct RowColSums{
    void resize(uint32_t barcodes, uint32_t total_rows, uint32_t total_cols){
        this->total_rows = total_rows;
        this->total_cols = total_cols;
        this->barcodes = barcodes;

        row_molecules.assign(total_rows, 0);
        row_reads.assign(total_rows, 0);
        row_support.assign(total_rows, 0);
        row_lost.assign(total_rows, 0);
        col_molecules.assign(total_cols, 0);
        col_reads.assign(total_cols, 0);
        col_support.assign(total_cols, 0);
        col_lost.assign(total_cols, 0);
        seed = 0; pidx.clear();
    }

    uint32_t              step;
    uint32_t              barcodes;
    uint32_t              total_rows;
    uint32_t              total_cols;
    uint32_t              max_hist;
    uint32_t              seed;
    std::vector<uint32_t> pidx;
    std::vector<uint32_t> row_reads;
    std::vector<uint32_t> row_molecules;
    std::vector<uint32_t> row_support;
    std::vector<uint32_t> row_lost;
    std::vector<uint32_t> col_reads;
    std::vector<uint32_t> col_molecules;
    std::vector<uint32_t> col_support;
    std::vector<uint32_t> col_lost;
};
*/

VisiumBin aggregate_visium_bins(const Downsampler& ds, uint32_t step,
                                uint32_t row_div, uint32_t col_div);

struct BinnedMap {
    void                build(const Downsampler& ds, uint32_t bin_div);
    size_t              trows;
    size_t              tcols;
    size_t              N;
    std::vector<size_t> bc2bin;
    std::vector<size_t> bin_ptr;
    std::vector<size_t> bin_barcodes;

    uint32_t bin_div;
};

bool make_filt_gbin_mat(const BinnedMap& bm, const GeneCountMatrix& gmat,
                        GeneCountMatrix& go, const std::vector<uint32_t>& idx);

bool make_filt_cbin_mat(const BinnedMap& bm, const SparseMatrix& gmat,
                        SparseMatrix& go, const std::vector<uint32_t>& idx);

void bin_gene_counts(const Downsampler& ds, uint32_t step, uint32_t bin_div,
                     SparseMatrix& out);

// void aggregate_visium_hists(const Downsampler & ds, RowColHists & out,
// uint32_t step, unsigned int min_molecules,
//                             bool in_tissue); /*, bool permute, double frac,
//                             uint32_t seed);*/

// void aggregate_visium(const Downsampler & ds, RowColSums & out, uint32_t
// step, unsigned int min_molecules,
//                             bool in_tissue, bool permute, double frac,
//                             uint32_t seed);

}; // namespace scdepth
