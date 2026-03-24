// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#include "visbins.hpp"
#include "downsampler.hpp"
#include <iostream>
#include <numeric>

using namespace scdepth;

void bin_genes(const GeneCountMatrix & gmat, const std::vector<uint64_t> & barcode2bin, 
        size_t n_bins, std::vector<uint32_t> & genes){
    std::vector<uint32_t> gidx;
    size_t n_barcodes = barcode2bin.size();
    genes.assign(n_bins, 0);

    std::vector<size_t> bin_ptr(n_bins + 1, 0);
    for(size_t i = 0; i < n_barcodes; i++) bin_ptr[barcode2bin[i] + 1]++;
    for(size_t i = 1; i < bin_ptr.size(); i++) bin_ptr[i] += bin_ptr[i - 1];

    std::vector<size_t> write_ptr = bin_ptr;
    std::vector<size_t> bin_barcodes(bin_ptr.back());

    for(size_t b = 0; b < n_barcodes; ++b) {
        bin_barcodes[write_ptr[barcode2bin[b]]++] = b;
    }

    //for each fraction
    const size_t n_rows = gmat.indptr.size() - 1;
    for(size_t bin = 0; bin < n_bins; ++bin){
        size_t begin = bin_ptr[bin], end = bin_ptr[bin + 1];
        gidx.clear();
        for(size_t i = begin; i < end; ++i){
            size_t barcode = bin_barcodes[i];
            if(barcode >= n_rows){
                //std::cerr << "error barcode " << barcode << " is greater than " << n_rows << "\n";
                continue;

            }
            //merge gene ids from each barcode in the bin
            auto g_start = std::next(gmat.indices.begin(), gmat.indptr[barcode]);
            auto g_end = std::next(gmat.indices.begin(), gmat.indptr[barcode + 1]);
            gidx.insert(gidx.end(), g_start, g_end);
        }
        std::sort(gidx.begin(), gidx.end());
        auto uend = std::unique(gidx.begin(), gidx.end()); 
        /*
        std::cerr << "n_bins = " << n_bins << " n_barcodes = " << n_barcodes 
            << " raw_count = " << gidx.size()
            << " count = " << std::distance(gidx.begin(), uend) << "\n";
        */
        //std::cout  n_bins = " << n_bins << " bin = " << bin << " size = " << bc_genes.size() << "\n";
        genes[bin] = std::distance(gidx.begin(), uend);
    }
}

void scdepth::bin_gene_counts(const Downsampler & ds, uint32_t step, uint32_t bin_div, SparseMatrix & out){
    auto & dso = ds.output;
    auto trows = (ds.total_rows + bin_div -1) / bin_div;
    auto tcols = (ds.total_cols + bin_div -1) / bin_div;
    auto n_bins = trows * tcols;
    std::vector<uint64_t> bc2bin;
    bc2bin.assign(dso.barcodes, 0);
    for(size_t i = 0; i < dso.barcodes; i++){
        size_t r = ds.bin_pos[i].r / bin_div;
        size_t c = ds.bin_pos[i].c / bin_div;
        size_t ii = r * tcols + c;
        bc2bin[i] = ii;
    }
    size_t n_barcodes = bc2bin.size();
    std::vector<size_t> bin_ptr(n_bins + 1, 0);
    for(size_t i = 0; i < n_barcodes; i++) bin_ptr[bc2bin[i] + 1]++;
    for(size_t i = 1; i < bin_ptr.size(); i++) bin_ptr[i] += bin_ptr[i - 1];

    std::vector<size_t> write_ptr = bin_ptr;
    std::vector<size_t> bin_barcodes(bin_ptr.back());

    for(size_t b = 0; b < n_barcodes; ++b) {
        bin_barcodes[write_ptr[bc2bin[b]]++] = b;
    }

    out.clear();
    out.indptr.reserve(n_bins + 1);
    std::vector<std::pair<uint32_t, uint32_t>> pairs;
    auto & mat = ds.output.total_mat[step];

    const size_t n_rows = mat.indptr.size() - 1;
    for(size_t bin = 0; bin < n_bins; ++bin){
        size_t begin = bin_ptr[bin], end = bin_ptr[bin + 1];
        pairs.clear();
        for(size_t i = begin; i < end; ++i){
            size_t barcode = bin_barcodes[i];
            if(barcode >= n_rows){
                continue;

            }

            auto g_start = std::next(mat.indices.begin(), mat.indptr[barcode]);
            auto g_end = std::next(mat.indices.begin(), mat.indptr[barcode + 1]);
            auto d_start = std::next(mat.data.begin(), mat.indptr[barcode]);
            while(g_start != g_end){
                pairs.push_back({*g_start, *d_start});
                g_start++;
                d_start++;
            }
        }

        std::sort(pairs.begin(), pairs.end());
        auto pstart = pairs.begin();
        while(pstart != pairs.end()){
            auto pit = std::next(pstart);
            while(pit != pairs.end() && pit->first == pstart->first){
                pstart->second += pit->second;
                pit++;
            }
            out.indices.push_back(pstart->first);
            out.data.push_back(pstart->second);
            pstart = pit;
        }
        out.indptr.push_back(static_cast<uint32_t>(out.indices.size()));
    }
}


void BinnedMap::build(const Downsampler & ds, uint32_t bin_div){
    auto & dso = ds.output;
    this->bin_div = bin_div;
    this->trows = (ds.total_rows + bin_div -1) / bin_div;
    this->tcols = (ds.total_cols + bin_div -1) / bin_div;
    this->N = trows * tcols;
    std::vector<uint64_t> bc2bin;
    bc2bin.assign(dso.barcodes, 0);
    for(size_t i = 0; i < dso.barcodes; i++){
        size_t r = ds.bin_pos[i].r / bin_div;
        size_t c = ds.bin_pos[i].c / bin_div;
        size_t ii = r * tcols + c;
        bc2bin[i] = ii;
    }
    size_t n_barcodes = bc2bin.size();
    bin_ptr.assign(N + 1, 0);
    for(size_t i = 0; i < n_barcodes; i++) bin_ptr[bc2bin[i] + 1]++;
    for(size_t i = 1; i < bin_ptr.size(); i++) bin_ptr[i] += bin_ptr[i - 1];

    std::vector<size_t> write_ptr = bin_ptr;
    bin_barcodes.assign(bin_ptr.back(), 0);

    for(size_t b = 0; b < n_barcodes; ++b) {
        bin_barcodes[write_ptr[bc2bin[b]]++] = b;
    }
}

bool scdepth::make_filt_gbin_mat(const BinnedMap & bm, const GeneCountMatrix & gmat, GeneCountMatrix & go, 
        const std::vector<uint32_t> & idx){

    const size_t n_rows = gmat.indptr.size() - 1;
    auto & bin_ptr = bm.bin_ptr;
    auto & bin_barcodes = bm.bin_barcodes;
    std::vector<uint32_t> gidx;
    /*
    std::cout << "gene count mat" 
        << " indptr = " << gmat.indptr.size()
        << " indices = " << gmat.indices.size()
        << "\n";
    */

    for(auto bin : idx){
        size_t begin = bin_ptr[bin], end = bin_ptr[bin + 1];
        gidx.clear();
        for(size_t i = begin; i < end; ++i){
            size_t barcode = bin_barcodes[i];
            if(barcode >= n_rows){
                continue;

            }
            auto g_start = std::next(gmat.indices.begin(), gmat.indptr[barcode]);
            auto g_end = std::next(gmat.indices.begin(), gmat.indptr[barcode + 1]);
            gidx.insert(gidx.end(), g_start, g_end);
        }
        std::sort(gidx.begin(), gidx.end());
        auto uend = std::unique(gidx.begin(), gidx.end()); 
        go.indices.insert(go.indices.end(), gidx.begin(), uend);
        go.indptr.push_back(go.indices.size());
    }
    return true;
}

bool scdepth::make_filt_cbin_mat(const BinnedMap & bm, const SparseMatrix & gmat, SparseMatrix & go, 
        const std::vector<uint32_t> & idx){

    const size_t n_rows = gmat.indptr.size() - 1;
    auto & bin_ptr = bm.bin_ptr;
    auto & bin_barcodes = bm.bin_barcodes;
    std::vector<std::pair<uint32_t, uint32_t>> gidx;

    for(auto bin : idx){
        size_t begin = bin_ptr[bin], end = bin_ptr[bin + 1];
        gidx.clear();
        for(size_t i = begin; i < end; ++i){
            size_t barcode = bin_barcodes[i];
            if(barcode >= n_rows){
                continue;
            }

            for(size_t j = gmat.indptr[barcode]; j < gmat.indptr[barcode + 1]; j++){
                gidx.push_back({gmat.indices[j], gmat.data[j]});
            }
        }

        std::sort(gidx.begin(), gidx.end());
        auto it = gidx.begin();
        while(it != gidx.end()){
            auto c = it->first;
            auto s = 0;
            while(it != gidx.end() && it->first == c){
                s += it->second;
                it++;
            }
            go.indices.push_back(c);
            go.data.push_back(s);
        }
        go.indptr.push_back(go.data.size());
    }
    return true;
}


VisiumBin scdepth::aggregate_visium_bins(const Downsampler & ds, uint32_t step, uint32_t row_div, uint32_t col_div){
    VisiumBin out;
    out.invalid = true;
    auto & dso = ds.output;
    if(ds.bin_pos.empty() || ds.barcodes.size() != ds.bin_pos.size()){
        std::cerr << "[error] the bin positions must be specified and the same size as `barcdoes`";
        return out;
    }
    if(dso.aggregate_only){
        std::cerr << "[error] cannot bin data with aggregate_only = false";
        return out;
    }
    if(col_div < 2 && row_div < 2){
        std::cerr << "Bin div must be at least 2 for one dimension";
        return out;
    }

    if(step >= dso.steps){
        std::cerr << "[error] invalid step\n";
        return out;
    }

    std::vector<uint64_t> bc2bin;

    size_t trows = (ds.total_rows + row_div -1) / row_div;
    size_t tcols = (ds.total_cols + col_div -1) / col_div;

    size_t N = trows * tcols;

    out.reads.assign(N, 0);
    out.mols.assign(N, 0);
    out.genes.assign(N, 0);

    if(dso.has_mt) out.mt_mols.assign(N, 0);
    if(dso.has_mod) out.mod_mols.assign(N, 0);
    out.in_tissue.assign(N, 0);
    out.countable.assign(N, 0);
    out.total.assign(N, 0);
    out.n_barcodes.assign(N, 0);

    out.n_rows = trows;
    out.n_cols = tcols;
    out.total_rows = ds.total_rows;
    out.total_cols = ds.total_cols;

    out.n_bins = N;
    out.bin_rows.resize(N);
    out.bin_cols.resize(N);
    for(size_t r = 0; r < trows; r++){
        for(size_t c = 0; c < tcols; c++){
            size_t idx = r * tcols + c;
            out.bin_rows[idx] = r;
            out.bin_cols[idx] = c;
        }
    }

    bc2bin.assign(dso.barcodes, 0);

    for(size_t i = 0; i < dso.barcodes; i++){
        size_t r = ds.bin_pos[i].r / row_div;
        size_t c = ds.bin_pos[i].c / col_div;

        //if(r >= trows || c >= tcols) continue;
        size_t idx = r * tcols + c;
        bc2bin[i] = idx;
        size_t didx = step * dso.barcodes + i;
        out.reads[idx] += dso.total.bc_reads[didx];
        out.mols[idx] += dso.total.bc_mols[didx];

        if(dso.has_mt) out.mt_mols[idx] += dso.total.bc_mt_mols[didx];
        if(dso.has_mod) out.mod_mols[idx] += dso.total.bc_mod_mols[didx];
        out.in_tissue[idx] += (ds.bin_pos[i].in_tissue > 0);
        out.countable[idx] += ds.bin_pos[i].countable;
        out.total[idx] += ds.bin_pos[i].total;
        out.n_barcodes[idx]++;
    }

    bin_genes(ds.output.gcounts[step], bc2bin, N, out.genes);
    out.invalid = false;

    return out;
}
/*
void scdepth::aggregate_visium_hists(const Downsampler & ds, RowColHists & out, uint32_t step,
        unsigned int min_molecules, bool in_tissue)
{
    auto & dso = ds.output;
    out.resize(dso.barcodes, ds.total_rows, ds.total_cols, dso.max_hist, dso.has_mt);
    const size_t H = out.max_hist - 1;
    uint32_t * rows = out.rows.data();
    uint32_t * cols = out.cols.data();
    const uint32_t * mh = ds.output.cell_mhist.data();
    for(size_t i = 0; i < dso.barcodes; i++){
        size_t r = ds.bin_pos[i].r;
        size_t c = ds.bin_pos[i].c;
        //size_t p = out.pidx.empty() ? i : out.pidx[i];
        size_t p = i;
        size_t mols = ds.output.total.bc_mols[size_t(step) * dso.barcodes + p];
        size_t reads = ds.output.total.bc_reads[size_t(step) * dso.barcodes + p];
        if(mols < min_molecules || (in_tissue && !ds.bin_pos[i].in_tissue)){
            continue;
        }
        out.row_molecules[r] += mols;
        out.row_reads[r] += reads;
        out.row_support[r]++;
        out.row_total[r] += ds.barcodes[p].total;
        out.row_countable[r] += ds.barcodes[p].countable;

        out.col_molecules[c] += mols;
        out.col_reads[c] += reads;
        out.col_support[c]++;
        out.col_total[c] += ds.barcodes[p].total;
        out.col_countable[c] += ds.barcodes[p].countable;
        if(dso.has_mt){
            size_t mt = ds.output.total.bc_mt_mols[size_t(step) * dso.barcodes + p];
            out.row_mt_molecules[r] += mt;
            out.col_mt_molecules[c] += mt;
        }

        const size_t in_idx  = ds.output.get_hidx_start(step, p);
        const uint32_t * in_ptr  = mh + in_idx;

        uint32_t* row_ptr = rows + r * H;
        uint32_t* col_ptr = cols + c * H;

        for (size_t m = 0; m < H; ++m) {
            const auto v = in_ptr[m];
            row_ptr[m] += v;
            col_ptr[m] += v;
        }
    }
}

void scdepth::aggregate_visium(const Downsampler & ds, RowColSums & out, uint32_t step,
        unsigned int min_molecules, bool in_tissue, bool permute, double frac, uint32_t seed)
{
    auto & dso = ds.output;
    out.resize(dso.barcodes, ds.total_rows, ds.total_cols);
    if(permute){
        out.pidx.resize(out.barcodes);
        std::iota(out.pidx.begin(), out.pidx.end(), 0);
        pcg64 rng(seed);
        std::shuffle(out.pidx.begin(), out.pidx.end(), rng);
    }else{
        out.pidx.clear();
        out.seed = 0;
    }
    pcg64 rng(seed);
    for(size_t i = 0; i < dso.barcodes; i++){
        size_t r = ds.bin_pos[i].r;
        size_t c = ds.bin_pos[i].c;
        size_t p = out.pidx.empty() ? i : out.pidx[i];
        size_t mols = ds.output.total.bc_mols[size_t(step) * dso.barcodes + p];
        size_t reads = ds.output.total.bc_reads[size_t(step) * dso.barcodes + p];
        if(mols < min_molecules || (in_tissue && !ds.bin_pos[i].in_tissue)){
            continue;
        }
        if (frac > 0.0 && frac < 1.0) {
            double u = std::generate_canonical<double, 53>(rng);
            if (u >= frac) {
                out.row_lost[r]++;
                out.col_lost[c]++;
                continue;
            }
        }

        out.row_molecules[r] += mols;
        out.row_reads[r] += reads;
        out.row_support[r]++;
        out.col_molecules[c] += mols;
        out.col_reads[c] += reads;
        out.col_support[c]++;
    }
}
*/
