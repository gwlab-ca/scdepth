// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#pragma once

#include "gtf.hpp"
#include "gtl/phmap.hpp"
#include "htslib/sam.h"
#include "shards.hpp"
#include "tags.hpp"
#include <regex>
#include <string>
#include <vector>

namespace scdepth {
using gtl::flat_hash_map;
using gtl::flat_hash_set;

enum TagSummary {
    TOTAL_READS     = 0,
    SPLICED_READS   = 1,
    AMBIGUOUS_READS = 2,
    UNSPLICED_READS = 3,
    COUNTABLE_READS = 4,
    LOW_QUALITY     = 5,
    BAD_TAGS        = 6,
    NO_GENE         = 7,
    TOTAL_FIELDS    = 8,
};

static constexpr std::array<std::string_view, TagSummary::TOTAL_FIELDS>
    TagNames = {"total_reads",     "spliced_reads",   "ambiguous_reads",
                "unspliced_reads", "countable_reads", "low_quality",
                "bad_tags",        "no_gene"};

struct TagCounter {
    using CountArray = std::array<uint64_t, TagSummary::TOTAL_FIELDS>;
    CountArray merged_counts{};
    CountArray polyA_counts{};
    CountArray random_hex_counts{};

    void inc(TagSummary f, bool has_random_hex, bool is_random_hex) {
        size_t idx = static_cast<size_t>(f);
        merged_counts[idx]++;
        inc_pa_hex(f, has_random_hex, is_random_hex);
    }

    void inc_pa_hex(TagSummary f, bool has_random_hex, bool is_random_hex) {
        size_t idx = static_cast<size_t>(f);
        if(has_random_hex) {
            if(is_random_hex)
                random_hex_counts[idx]++;
            else
                polyA_counts[idx]++;
        }
    }
};

class BarcodeCounter {
public:
    using hmap = flat_hash_map<std::string, uint32_t>;
    using cmap = flat_hash_map<std::string, TagCounter>;
    using smap = flat_hash_set<std::string>;
    ~BarcodeCounter() {
        if(bh_ != nullptr) {
            sam_hdr_destroy(bh_);
            bh_ = nullptr;
        }
        if(bf_ != nullptr) {
            sam_close(bf_);
            bf_ = nullptr;
        }
    }

    void set_count_parameters(double       min_gene              = 0.95,
                              unsigned int min_gene_bases        = 40,
                              double       min_exonic            = 0.95,
                              unsigned int min_intronic          = 15,
                              uint8_t      min_qual              = 255,
                              bool         discard_unknown_juncs = false,
                              bool         probes                = false);

    void init(const std::string& lib_string, bool fwd,
              const std::string& barcode_tag, const std::string& barcode_re,
              const std::string& umi_tag, const std::string& sample_tag = "",
              const std::vector<std::string>& samples          = {},
              const std::string&              random_hex_re    = "",
              const std::string&              random_hex_value = "",
              uint32_t barcode_length = 0, uint32_t umi_length = 0);

    bool prepare_bam(const std::string& gtf, const std::string& bam,
                     const std::string& out, int threads = 1,
                     size_t max_tags = 200000000, double max_tag_frac = 0.95);

    size_t process_reads(size_t chunk);

    bool done() const { return done_; }

    size_t total_reads() const {
        return full_.merged_counts[TagSummary::TOTAL_READS];
    }

    size_t countable_reads() const {
        return full_.merged_counts[TagSummary::COUNTABLE_READS];
    }

    bool finish();

private:
    char read2strand_(const bam1_t* b) const {
        bool rev = bam_is_rev(b);
        if((strand_ == StrandMode::TAG_FWD && !rev) ||
           (strand_ == StrandMode::TAG_REV && rev))
            return '+';
        return '-';
    }

    std::vector<BarcodeCount> barcode_counts_;
    std::vector<Exon>         read_blocks_;
    Genes                     genes_;
    RefTrees                  gtrees_;
    TagShards                 shards_;
    std::regex                barcode_re_;
    std::regex                random_hex_re_;
    hmap                      barcode_map_;
    hmap                      gene_map_;
    std::string               bam_file_;
    std::string               out_file_;
    std::string               lib_string_;
    std::string               barcode_regex_str_;
    std::string               random_hex_regex_str_;
    std::string               random_hex_value_;
    TagCounter                full_;
    cmap                      ssum_;

    std::vector<std::string> samples_;
    smap                     sample_set_;
    char                     barcode_tag_[2];
    char                     umi_tag_[2];
    char                     sample_tag_[2];
    samFile*                 bf_                    = nullptr;
    sam_hdr_t*               bh_                    = nullptr;
    double                   min_gene_              = 0.95;
    double                   min_exonic_            = 0.95;
    unsigned int             min_gene_bases_        = 40;
    unsigned int             min_intronic_          = 15;
    uint32_t                 barcode_length_        = 0;
    uint32_t                 umi_length_            = 0;
    StrandMode               strand_                = StrandMode::TAG_UNKNOWN;
    uint8_t                  min_qual_              = 255;
    bool                     done_                  = true;
    bool                     visium_hd_             = false;
    bool                     discard_unknown_juncs_ = false;
    bool                     has_probes_            = false;
};

unsigned int block_overlaps(const std::vector<Exon>& blocks,
                            const std::vector<Exon>& exons);
unsigned int check_junctions(const std::vector<Exon>& blocks,
                             const Gene::Juncs&       junctions);
unsigned int get_read_blocks(bam1_t* b, std::vector<Exon>& blocks);
} // namespace scdepth
