// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#pragma once

#include "gtl/phmap.hpp"
#include "tags.hpp"
#include <cstdint>
#include <vector>

namespace scdepth {
using gtl::flat_hash_map;

struct CorrUMITag {
    void make_tag(uint32_t gene, uint32_t umi) {
        this->gene = gene;
        this->umi  = umi;
    }

    bool operator<(const CorrUMITag& rhs) const {
        return std::tie(gene, umi) < std::tie(rhs.gene, rhs.umi);
    }

    uint32_t total() const {
        return uint32_t(spliced) + uint32_t(unspliced) + uint32_t(ambiguous);
    }

    static constexpr uint16_t U16MAX = std::numeric_limits<uint16_t>::max();
    static constexpr uint8_t  F_SPLICED_OVF   = 1u << 0;
    static constexpr uint8_t  F_UNSPLICED_OVF = 1u << 1;
    static constexpr uint8_t  F_AMBIG_OVF     = 1u << 2;
    static constexpr uint8_t  F_INVALID       = 1u << 3;
    static constexpr uint8_t  F_NRAW_OVF      = 1u << 4;
    static constexpr uint8_t  F_RANDOM_HEX    = 1u << 5;

    void inc_spliced(uint16_t v) { inc_sat(spliced, v, F_SPLICED_OVF); }
    void inc_unspliced(uint16_t v) { inc_sat(unspliced, v, F_UNSPLICED_OVF); }
    void inc_ambiguous(uint16_t v) { inc_sat(ambiguous, v, F_AMBIG_OVF); }

    void inc_n_umis(uint8_t c = 1) {
        uint16_t v = uint16_t(n_umis) + c;
        if(v > 255) {
            n_umis = 255;
            flags |= F_NRAW_OVF;
        } else {
            n_umis = uint8_t(v);
        }
    }

    void set_invalid() { flags |= F_INVALID; }

    bool invalid() const { return (flags & F_INVALID) != 0; }
    bool spliced_overflow() const { return (flags & F_SPLICED_OVF) != 0; }
    bool unspliced_overflow() const { return (flags & F_UNSPLICED_OVF) != 0; }
    bool ambiguous_overflow() const { return (flags & F_AMBIG_OVF) != 0; }
    bool n_umi_overflow() const { return (flags & F_NRAW_OVF) != 0; }
    bool is_random_hex() const { return (flags & F_RANDOM_HEX) != 0; }

    uint32_t gene      = 0;
    uint32_t umi       = 0;
    uint16_t spliced   = 0;
    uint16_t unspliced = 0;
    uint16_t ambiguous = 0;
    uint8_t  flags     = 0;
    uint8_t  n_umis    = 0;

    void inc_sat(uint16_t& x, uint16_t v, uint16_t flag_bit) {
        uint32_t sum = uint32_t(x) + uint32_t(v);
        if(sum <= U16MAX) {
            x = uint16_t(sum);
        } else {
            x = U16MAX;
            flags |= flag_bit;
        }
    }
};

struct UMIResults {
    void reset() { umis.clear(); }

    std::vector<CorrUMITag> umis;
    uint32_t                dis_reads = 0;
    uint32_t                dis_mols  = 0;
    uint32_t                dis_ambig = 0;

    // uint32_t                   ed0_cross = 0;
    // uint32_t                   ed1_cross = 0;
    // uint32_t                   ed1_same = 0;
};

struct UMINode {
    uint32_t umi       = 0;
    uint32_t spliced   = 0;
    uint32_t unspliced = 0;
    uint32_t ambiguous = 0;
    uint32_t flags     = 0;

    bool is_random_hex() const { return (flags & F_RANDOM_HEX) != 0; }

    uint32_t total() const { return spliced + unspliced + ambiguous; }
    static constexpr uint8_t F_RANDOM_HEX = 1u << 5;
    std::vector<uint16_t>    edges;
};

class UMIDirectional {

public:
    void downsample(double frac, uint32_t umi_len, const BarcodeCount& bc,
                    const BarcodeData& bd, UMIResults& out, bool directional,
                    bool       correct_multi_umis,
                    PrimerMode primer_mode); //, bool profile_singles);

private:
    template <typename T> void build_graph_brute_(T build_edge);
    template <typename T> void build_graph_hash_(T build_edge);
    void                       correct_umis_(uint32_t gene, UMIResults& out);
    void correct_multi_umis_(UMIResults& out); //, bool profile_singles);
    // void profile_singles_(UMIResults & out);

    template <typename T> struct CacheVect {
        void resize_clear(size_t n) {
            c = n;
            if(v.size() < c)
                v.resize(c);
            for(size_t i = 0; i < c; i++) {
                v[i].clear();
            }
        }

        void push_back() {
            if(c == v.size())
                v.emplace_back();
            else
                v[c].clear();
            c++;
        }

        void clear() { c = 0; }

        size_t size() const { return c; }

        std::vector<T>& back() { return v[c - 1]; }

        const std::vector<T>& back() const { return v[c - 1]; }

        std::vector<T>& operator[](size_t i) { return v[i]; }

        const std::vector<T>& operator[](size_t i) const { return v[i]; }

        size_t                      c = 0;
        std::vector<std::vector<T>> v;
    };

    std::vector<UMINode> nodes_;
    std::vector<UMITag>  tags_;
    // CacheVect<uint32_t>                         rev_;
    CacheVect<uint32_t>                            comps_full_;
    std::vector<int>                               visited_;
    std::vector<uint32_t>                          stack_;
    std::vector<uint32_t>                          idx_;
    flat_hash_map<uint32_t, uint32_t>              umi_map_;
    flat_hash_map<uint32_t, std::vector<uint32_t>> ed_umi_map_;
    // std::vector<CorrUMITag>                       final_umis_;
    size_t   n_nodes_ = 0;
    uint32_t umi_len_ = 0;
};

} // namespace scdepth
