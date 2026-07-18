// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#pragma once

#include "dsdata.hpp"
#include "htslib/bgzf.h"
#include <cassert>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace scdepth {

struct GeneOverlapStep {
    void reset(size_t barcodes, size_t seeds) {
        size_t N = barcodes * seeds;

        stability.assign(N, 0.0);
        stability_1.assign(N, 0.0);
        stability_2.assign(N, 0.0);
        stability_3p.assign(N, 0.0);
    }

    std::vector<double> stability;
    std::vector<double> stability_1;
    std::vector<double> stability_2;
    std::vector<double> stability_3p;
};

struct GeneOverlapMean {
    void reset(size_t barcodes) {
        stability.assign(barcodes, 0.0);
        stability_1.assign(barcodes, 0.0);
        stability_2.assign(barcodes, 0.0);
        stability_3p.assign(barcodes, 0.0);
    }

    std::vector<double> stability;
    std::vector<double> stability_1;
    std::vector<double> stability_2;
    std::vector<double> stability_3p;
};

class GeneOverlap {
public:
    GeneOverlap()                              = default;
    GeneOverlap(const GeneOverlap&)            = delete;
    GeneOverlap& operator=(const GeneOverlap&) = delete;

    GeneOverlap(GeneOverlap&&) noexcept            = default;
    GeneOverlap& operator=(GeneOverlap&&) noexcept = default;

    bool init(const std::string&              baseline_file,
              const std::vector<std::string>& files);
    bool calculate();

    std::vector<uint32_t> counts;
    std::vector<uint32_t> counts_1;
    std::vector<uint32_t> counts_2;
    std::vector<uint32_t> counts_3p;

    std::vector<GeneOverlapMean> mean_outs;
    uint32_t                     steps;
    uint32_t                     genes;
    uint32_t                     barcodes;
    uint32_t                     seeds;

private:
    bool read_baseline_(const std::string& baseline_file);
    bool calc_barcode_(uint32_t barcode);

    struct SeedMats {
        SeedMats() {}

        ~SeedMats() {
            if(bin != nullptr)
                bgzf_close(bin);
            bin = nullptr;
        }

        bool init(const std::string& file);
        bool next();

        double frac() { return fracs[idx]; }

        std::string         file;
        GeneCountMatrix     mat;
        std::vector<double> fracs;
        size_t              idx = 0;
        uint64_t            seed;
        uint32_t            genes;
        uint32_t            barcodes;
        BGZF*               bin = nullptr;
    };

    std::vector<SeedMats> seeds_;
    SparseMatrix          baseline_;
    GeneOverlapStep       sout_;
};

} // namespace scdepth
