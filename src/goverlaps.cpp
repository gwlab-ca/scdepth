// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#include "goverlaps.hpp"
#include "dsmisc.hpp"
#include <iostream>
#include <cmath>
#include <iostream>

using namespace scdepth;

static inline double mean_over_seeds(const std::vector<double>& v, size_t off, uint32_t seeds){
    double sum = 0.0;
    for(uint32_t i = 0; i < seeds; i++) sum += v[off + i];
    return (seeds ? sum / double(seeds) : 0.0);
}

bool GeneOverlap::calculate(){

    counts.assign(barcodes, 0);
    counts_1.assign(barcodes, 0);
    counts_2.assign(barcodes, 0);
    counts_3p.assign(barcodes, 0);

    mean_outs.clear();
    mean_outs.resize(steps);
    for(auto & o : mean_outs) o.reset(barcodes);

    //std::cout << "Counting barcodes\n";
    for(size_t barcode = 0; barcode < barcodes; barcode++){
        const size_t base0 = baseline_.indptr[barcode];
        const size_t base1 = baseline_.indptr[barcode + 1];
        const uint32_t* b_dat  = baseline_.data.data();

        uint32_t B  = uint32_t(base1 - base0);
        uint32_t B1 = 0, B2 = 0, B3p = 0;
        for(size_t p = base0; p < base1; ++p){
            uint32_t c = b_dat[p];
            if(c == 1) B1++;
            else if(c == 2) B2++;
            else if(c >= 3) ++B3p;
        }
        counts[barcode] = B;
        counts_1[barcode] = B1;
        counts_2[barcode] = B2;
        counts_3p[barcode] = B3p;
    }

    for(size_t step = 0; step < steps; step++){
        //std::cout << "Step " << (step + 1) << " / " << steps << "\n";
        sout_.reset(barcodes, seeds);
        for(auto & s : seeds_){
            if(!s.next()){
                std::cerr << "[error] Error loading next seed for file " << s.file << "\n";
            }
        }
        for(size_t b = 0; b < barcodes; b++){
            calc_barcode_(b);
        }

        auto& m = mean_outs[step];


        for(uint32_t bc = 0; bc < barcodes; bc++){
            const size_t off = size_t(bc) * seeds;

            m.stability[bc] = mean_over_seeds(sout_.stability, off, seeds);
            m.stability_1[bc] = mean_over_seeds(sout_.stability_1, off, seeds);
            m.stability_2[bc] = mean_over_seeds(sout_.stability_2, off, seeds);
            m.stability_3p[bc] = mean_over_seeds(sout_.stability_3p, off, seeds);
        }
    }
    return true;
}

bool GeneOverlap::calc_barcode_(uint32_t barcode){
    const size_t base0 = baseline_.indptr[barcode];
    const size_t base1 = baseline_.indptr[barcode + 1];

    const uint32_t* b_idx  = baseline_.indices.data();
    const uint32_t* b_dat  = baseline_.data.data();

    const size_t seeds = seeds_.size();
    const size_t row_offset = size_t(barcode) * seeds;
    uint32_t B = counts[barcode];
    uint32_t B1 = counts_1[barcode];
    uint32_t B2 = counts_2[barcode];
    uint32_t B3p = counts_3p[barcode];

    for(size_t i = 0; i < seeds; i++){
        auto & s = seeds_[i];
        const size_t down0 = s.mat.indptr[barcode];
        const size_t down1 = s.mat.indptr[barcode + 1];

        const uint32_t* d_idx = s.mat.indices.data();

        uint32_t I = 0, I1 = 0, I2 = 0, I3p = 0;

        size_t p = base0;
        size_t q = down0;
        while(p < base1 && q < down1){
            uint32_t bi = b_idx[p];
            uint32_t di = d_idx[q];

            if(bi == di){
                I++;
                uint32_t c = b_dat[p];
                if(c == 1) I1++;
                else if(c == 2) I2++;
                else if(c >= 3) I3p++;
                p++; q++;
            } else if (bi < di){
                p++;
            } else {
                q++;
            }
        }

        const size_t k = row_offset + i;

        sout_.stability[k]    = (B  ? double(I)  / double(B)  : 0.0);
        sout_.stability_1[k]  = (B1 ? double(I1) / double(B1) : 0.0);
        sout_.stability_2[k]  = (B2 ? double(I2) / double(B2) : 0.0);
        sout_.stability_3p[k] = (B3p? double(I3p)/ double(B3p): 0.0);
    }

    return true;
}

bool GeneOverlap::SeedMats::init(const std::string & file){
    idx = 0;
    this->file = file;
    bin = bgzf_open(file.c_str(), "r");
    if(bin == nullptr){
        std::cerr << "[error] Error opening seed " << file << "\n";
        return false;
    }

    if(bgzf_read(bin, reinterpret_cast<void*>(&seed), sizeof(seed)) <= 0){
        std::cerr << "[error] Error reading seed " << file << "\n";
        bgzf_close(bin);
        return false;
    }

    if(bgzf_read(bin, reinterpret_cast<void*>(&genes), sizeof(genes)) <= 0){
        std::cerr << "[error] Error reading genes " << file << "\n";
        bgzf_close(bin);
        return false;
    }

    if(bgzf_read(bin, reinterpret_cast<void*>(&barcodes), sizeof(barcodes)) <= 0){
        std::cerr << "[error] Error reading barcodes " << file << "\n";
        bgzf_close(bin);
        return false;
    }

    if(!bz_read_vector(bin, fracs)){
        std::cerr << "[error] Error reading fracs " << file << "\n";
        bgzf_close(bin);
        return false;
    };

    return true;
}

bool GeneOverlap::SeedMats::next(){
    if(idx >= fracs.size()) return false;
    if(!bz_read_vector(bin, mat.indptr)){
        std::cerr << "[error] Error reading indptr " << file << "\n";
        bgzf_close(bin);
        return false;
    }

    if(!bz_read_vector(bin, mat.indices)){
        std::cerr << "[error] Error reading indices " << file << "\n";
        bgzf_close(bin);
        return false;
    }
    return true;
}

bool vectors_equal(const std::vector<double>& a, const std::vector<double>& b, double eps = 1e-10){
    if (a.size() != b.size()) return false;

    for(size_t i = 0; i < a.size(); ++i)
        if(std::fabs(a[i] - b[i]) > eps) return false;

    return true;
}

bool GeneOverlap::init(const std::string & baseline_file, const std::vector<std::string> & files){
    if(files.empty())    return false;
    if(files.size() < 2) return false;

    if(!read_baseline_(baseline_file)) return false;
    seeds_.resize(files.size());
    for(size_t i = 0; i < files.size(); i++){
        if(!seeds_[i].init(files[i])){
            std::cerr << "[error] Error opening " << files[i] << "\n";
        }
    }
    steps = seeds_[0].fracs.size();
    seeds = seeds_.size();

    for(size_t i = 1; i < seeds_.size(); i++){
        if(seeds_[i].barcodes != barcodes){
            std::cerr << "[error] barcodes for seed and baseline (" << barcodes << " vs " << seeds_[i].barcodes << ")\n";
            seeds_.clear();
            return false;
        }
        if(seeds_[i].genes != genes){
            std::cerr << "[error] genes for seed and baseline (" << genes << " vs " << seeds_[i].genes << ")\n";
            seeds_.clear();
            return false;
        }

        if(!vectors_equal(seeds_[i].fracs, seeds_[0].fracs)){
            std::cerr << "[error] fractions for each seed are not equal (" << seeds_[0].file << " vs " << seeds_[i].file << ")\n";
            seeds_.clear();
            return false;
        }
    }
    return true;
}

bool GeneOverlap::read_baseline_(const std::string & baseline_file){
    BGZF * bin = bgzf_open(baseline_file.c_str(), "r");
    baseline_.clear();
    if(bin == nullptr){
        return false;
    }

    if(bgzf_read(bin, reinterpret_cast<void*>(&genes), sizeof(genes)) <= 0){
        std::cerr << "[error] Error reading genes " << baseline_file << "\n";
        bgzf_close(bin);
        return false;
    }

    if(bgzf_read(bin, reinterpret_cast<void*>(&barcodes), sizeof(barcodes)) <= 0){
        std::cerr << "[error] Error reading barcodes " << baseline_file << "\n";
        bgzf_close(bin);
        return false;
    }

    if(!bz_read_vector(bin, baseline_.indptr)){
        std::cerr << "[error] Error reading indptr " << baseline_file << "\n";
        bgzf_close(bin);
        return false;
    }

    if(!bz_read_vector(bin, baseline_.indices)){
        std::cerr << "[error] Error reading indices " << baseline_file << "\n";
        bgzf_close(bin);
        return false;
    }

    if(!bz_read_vector(bin, baseline_.data)){
        std::cerr << "[error] Error reading data " << baseline_file << "\n";
        bgzf_close(bin);
        return false;
    }

    return true;

}
