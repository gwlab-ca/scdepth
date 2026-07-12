// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#pragma once

#include <cstdint>
#include "gtl/phmap.hpp"
#include "htslib/bgzf.h"


namespace scdepth{

using gtl::flat_hash_map;

struct GenePos{
    std::string ref;
    std::string gid;
    std::string gname;
    uint32_t    lft;
    uint32_t    rgt;
};

struct RawTag{
    static constexpr uint32_t F_UNSPLICED  = 1u << 0;
    static constexpr uint32_t F_SPLICED    = 1u << 1;
    static constexpr uint32_t F_RANDOM_HEX = 1u << 2;

    void make_tag(uint32_t barcode, uint32_t gene, uint32_t umi, bool spliced, bool unspliced, bool random_hex){
        this->barcode = barcode;
        this->gene = gene;
        this->umi = umi;
        pack_flags(spliced, unspliced, random_hex);
    }

    void pack_flags(bool spliced, bool unspliced, bool random_hex = false){
        flag = 0;
        if (spliced) {
            flag |= F_SPLICED;
        }

        if (unspliced) {
            flag |= F_UNSPLICED;
        }

        if (random_hex) {
            flag |= F_RANDOM_HEX;
        }
    }

    bool is_spliced() const {
        return (flag & F_SPLICED) != 0;
    }

    bool is_unspliced() const {
        return (flag & F_UNSPLICED) != 0;
    }

    bool is_random_hex() const {
        return (flag & F_RANDOM_HEX) != 0;
    }

    uint32_t unpack_random_hex() const {
        return is_random_hex() ? 1u : 0u;
    }

    uint32_t unpack_spliced() const {
        return is_spliced() ? 1u : 0u;
    }

    uint32_t unpack_unspliced() const {
        return is_unspliced() ? 1u : 0u;
    }

    bool operator<(const RawTag & rhs) const{
        return std::tie(barcode, gene, umi, flag) 
            < std::tie(rhs.barcode, rhs.gene, rhs.umi, rhs.flag);
    }

    bool same_umi(const RawTag & rhs) const {
        bool rh1 = is_random_hex();
        bool rh2 = rhs.is_random_hex();
        return std::tie(barcode, gene, umi, rh1) 
            == std::tie(rhs.barcode, rhs.gene, rhs.umi, rh2);
    }

    bool operator==(const RawTag & rhs) const{
        return std::tie(barcode, gene, umi, flag) 
            == std::tie(rhs.barcode, rhs.gene, rhs.umi, rhs.flag);
    }

    uint32_t barcode;
    uint32_t gene;
    uint32_t umi;
    uint32_t flag;
};

struct BarcodeCount{
    bool operator<(const BarcodeCount & rhs) const{
        return index < rhs.index;
    }

    std::string barcode;
    uint64_t    offset = 0;
    uint64_t    seed = 0;
    uint32_t    index = 0;
    uint32_t    total = 0;
    uint32_t    countable = 0;
    uint32_t    raw_molecules = 0;
    uint32_t    total_data_bytes = 0;
    std::string sample;
    uint32_t    random_hex = 0;
    uint32_t    poly_a = 0;
    bool        has_qc = false;
};

struct BarcodeHeader{
    void make_header(uint32_t barcode, uint32_t n_count_tags = 0, uint32_t block_size = 0){
        this->barcode = barcode;
        this->n_count_tags = n_count_tags;
        this->block_size = block_size;
    }
    uint32_t barcode;
    uint32_t n_count_tags = 0; // how many counttag entries
    uint32_t block_size = 0;
};

static_assert(std::is_trivially_copyable_v<BarcodeHeader>);
static_assert(std::is_standard_layout_v<BarcodeHeader>);
static_assert(sizeof(BarcodeHeader) == 12, "BarcodeHeader size mismatch");

struct UMITag {
    void make_tag(uint32_t gene, uint32_t umi){
        this->gene = gene;
        this->umi  = umi;
    }

    bool operator<(const UMITag& rhs) const{
        bool rh1 = is_random_hex(), rh2 = rhs.is_random_hex();
        return std::tie(gene, umi, rh1) < std::tie(rhs.gene, rhs.umi, rh2);
    }

    bool operator==(const UMITag& rhs) const{
        bool rh1 = is_random_hex(), rh2 = rhs.is_random_hex();
        return std::tie(gene, umi, rh1) == std::tie(rhs.gene, rhs.umi, rh2);
    }

    uint32_t total() const{
        return uint32_t(spliced) + uint32_t(unspliced) + uint32_t(ambiguous);
    }

    static constexpr uint16_t U16MAX = std::numeric_limits<uint16_t>::max();
    static constexpr uint16_t F_SPLICED_OVF   = 1u << 0;
    static constexpr uint16_t F_UNSPLICED_OVF = 1u << 1;
    static constexpr uint16_t F_AMBIG_OVF     = 1u << 2;
    static constexpr uint16_t F_INVALID       = 1u << 3;
    static constexpr uint16_t F_RANDOM_HEX    = 1u << 4;

    void inc_spliced(uint16_t v)   { inc_sat(spliced,   v, F_SPLICED_OVF); }
    void inc_unspliced(uint16_t v) { inc_sat(unspliced, v, F_UNSPLICED_OVF); }
    void inc_ambiguous(uint16_t v) { inc_sat(ambiguous, v, F_AMBIG_OVF); }

    void set_invalid(){ flags |= F_INVALID; }

    bool invalid() const            { return (flags & F_INVALID) != 0; }
    bool spliced_overflow() const   { return (flags & F_SPLICED_OVF) != 0; }
    bool unspliced_overflow() const { return (flags & F_UNSPLICED_OVF) != 0; }
    bool ambiguous_overflow() const { return (flags & F_AMBIG_OVF) != 0; }
    bool is_random_hex() const      { return (flags & F_RANDOM_HEX) != 0; }

    void set_random_hex(bool value) {
        if (value) {
            flags |= F_RANDOM_HEX;
        } else {
            flags &= ~F_RANDOM_HEX;
        }
    }

    bool same_molecule(const UMITag& rhs) const {
        return gene == rhs.gene &&
            umi == rhs.umi &&
            is_random_hex() == rhs.is_random_hex();
    }


    uint32_t gene = 0;
    uint32_t umi = 0;
    uint16_t spliced = 0;
    uint16_t unspliced = 0;
    uint16_t ambiguous = 0;
    uint16_t flags = 0;

    void inc_sat(uint16_t& x, uint16_t v, uint16_t flag_bit){
        uint32_t sum = uint32_t(x) + uint32_t(v);
        if (sum <= U16MAX) {
            x = uint16_t(sum);
        } else {
            x = U16MAX;
            flags |= flag_bit;
        }
    }
};

static_assert(std::is_trivially_copyable_v<UMITag>);
static_assert(std::is_standard_layout_v<UMITag>);
static_assert(sizeof(UMITag) == 16, "unexpected padding on CountTag");

struct BarcodeData{
    BarcodeHeader               bc;
    std::vector<UMITag>         counts;
};

enum StrandMode {
    TAG_UNKNOWN = 0,
    TAG_FWD = 1, // Tag maps to the forward strand // For 10X 3' R2 maps to the forward strand of the transcript
    TAG_REV = 2  // Tag maps to the reverse strand // For 10X 5' R2 maps to the reverse strand of the transcript
};

enum class PrimerMode {
    Merge,
    PolyAOnly,
    RandomHexOnly
};

struct ADNA4 {
    static const std::size_t size_ = 2;
    static const char * alphabet_str_;
    static const char ltable_[256];
    static const char rtable_[4];
};

inline std::pair<uint32_t, bool> seq2int(const char * seq, uint32_t length){
    uint32_t code = 0;
    for(size_t i = 0; i < length; i++){
        char c = ADNA4::ltable_[static_cast<size_t>(seq[i])];
        if(c < 0) return {0, false};
        code = (code << ADNA4::size_) | c;
    }
    return {code, true};
}

inline std::string int2seq(uint32_t val, uint32_t N){
    std::string seq("", N);
    //cout << "val = " << val << "\n";
    uint32_t mask = (1 << ADNA4::size_) - 1;
    for(size_t i = 0; i < N; i++){
        //cout << "  i = " << i << " val = " << val << " " << (val & 0x3) << "\n";
        seq[N - i - 1] = ADNA4::alphabet_str_[val & mask];
        val >>= ADNA4::size_;
    }
    return seq;
}

long int read_barcode(BGZF * fin, BarcodeData & data);
long int read_barcode_debug(BGZF * fin, BarcodeData & data, const std::vector<BarcodeCount> & bcs);
std::vector<BarcodeCount> read_barcode_index(const std::string & fin);
}
