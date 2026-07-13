// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#pragma once


#include "tags.hpp"
#include <cstdint>
#include <vector>
#include <string>
#include <cassert>
#include "dsdata.hpp"
#include "dsmisc.hpp"

namespace scdepth{

using gtl::flat_hash_map;
class Downsampler{
    public:
        using smap = flat_hash_map<std::string, uint32_t>;
        Downsampler() = default;
        Downsampler(const Downsampler&) = delete;
        Downsampler& operator=(const Downsampler&) = delete;

        Downsampler(Downsampler&&) noexcept = default;
        Downsampler& operator=(Downsampler&&) noexcept = default;

        bool init(const std::string & prefix, const std::string & mt_prefix, const std::string & mt_file, 
                const std::string & mod_file, const std::string & exclude_file,  
                size_t max_hist = 40, bool build_matrices = false, bool calc_sau = false);

        bool init_visium(const std::vector<uint32_t> & rows, const std::vector<uint32_t> & cols, 
                std::vector<uint32_t> & in_tissue, std::vector<uint32_t> & countable_reads, 
                std::vector<uint32_t> & total_reads, size_t total_rows, size_t total_cols);

        bool reset_visium();

        void clear_output() {
            output = DownsampleResults();
        }

        bool downsample(std::vector<double> & fracs, 
                uint32_t umi_len, uint64_t seed, unsigned int threads = 1, bool aggregate_only = false,
                const std::string & umi_mode = "directed", bool correct_multi_umis = true,
                const std::string & primer_mode = "merge",
                const std::vector<uint32_t> & barcode2sample = {});

        bool write_gene_mats(const std::string & out, const std::vector<uint32_t> & idx,
                uint32_t bin_div = 0);

        bool write_gene_baseline(const std::string & out, const std::vector<uint32_t> & idx,
                uint32_t step = 0, uint32_t bin_div = 0);

        DownsampleResults           output;

        std::vector<BarcodePos>     bin_pos;
        std::vector<BarcodeCount>   barcodes;
        std::vector<bool>           mt_filter;
        std::vector<bool>           exclude_filter;
        std::vector<bool>           mod_filter;
        std::vector<uint32_t>       barcode2sample;
        std::string                 prefix;
        std::string                 primer_mode_str = "";

        uint32_t                    bin_div = 0;
        uint32_t                    total_rows = 0;
        uint32_t                    total_cols = 0;
        uint32_t                    samples = 0;

        PrimerMode                  primer_mode = PrimerMode::Merge;

        bool                        has_visium = false;
        bool                        aggregate_only = false;
        bool                        umi_directed = true;
        bool                        umi_multi = false;
        bool                        has_mt = false;
        bool                        has_mod = false;
        bool                        primer_filter = false;

    private:
        std::vector<size_t>         chunks_;

        uint64_t                    seed_ = 42;
        uint64_t                    countable_= 0;
        uint64_t                    raw_molecules_= 0;
        size_t                      max_hist_= 40;
        size_t                      genes_= 0;
        uint32_t                    umi_len_= 0;
        int                         threads_= 1;
        bool                        build_mats_= false;
        bool                        calc_sau_ = false;
        bool                        has_exc_ = false;
};

}
