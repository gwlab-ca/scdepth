// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#pragma once

#include <string>
#include <vector>
#include <regex>
#include "gtl/phmap.hpp"
#include "htslib/sam.h"
#include "tags.hpp"
#include "shards.hpp"
#include "gtf.hpp"

namespace scdepth{
using gtl::flat_hash_map;
using gtl::flat_hash_set;

struct TagSummary{
    uint64_t                     low_quality = 0;
    uint64_t                     bad_tags = 0;
    uint64_t                     no_gene = 0;
    uint64_t                     countable_reads = 0;
    uint64_t                     total_reads = 0;
    uint64_t                     spliced_reads = 0;
    uint64_t                     ambiguous_reads = 0;
    uint64_t                     unspliced_reads = 0;
};

class BarcodeCounter{
    public:
        using hmap = flat_hash_map<std::string, uint32_t>;
        using smap = flat_hash_set<std::string>;
        using tsmap = flat_hash_map<std::string, TagSummary>;
        ~BarcodeCounter(){
            if(bh_ != nullptr){
                sam_hdr_destroy(bh_);
                bh_ = nullptr;
            }
            if(bf_ != nullptr){
                sam_close(bf_);
                bf_ = nullptr;
            }
        }

        void set_count_parameters(double min_gene = 0.95, unsigned int min_gene_bases = 40, 
                double min_exonic = 0.95, unsigned int min_intronic = 15, uint8_t min_qual = 255,
                bool discard_unknown_juncs = false, bool probes = false);

        void init(const std::string & lib_string, bool fwd, const std::string & barcode_tag, 
                const std::string & barcode_re, const std::string & umi_tag, 
                const std::string & sample_tag = "", const std::vector<std::string> & samples = {},
                const std::string & random_hex_re = "", const std::string & random_hex_value = "",
                uint32_t barcode_length = 0, uint32_t umi_length = 0);

        bool prepare_bam(const std::string & gtf, const std::string & bam, const std::string & out,
                int threads = 1, size_t max_tags = 200000000, double max_tag_frac = 0.95);

        size_t process_reads(size_t chunk);

        bool done() const {
            return done_;
        }

        size_t total_reads() const {
            return full_.total_reads;
        }

        size_t countable_reads() const {
            return full_.countable_reads;
        }

        bool finish();
    private:
        char read2strand_(const bam1_t * b) const{
            bool rev = bam_is_rev(b);
            if((strand_ == StrandMode::TAG_FWD && !rev) || 
                    (strand_ == StrandMode::TAG_REV && rev)) 
                return '+';
            return '-';
        }

        std::vector<BarcodeCount>    barcode_counts_;
        std::vector<Exon>            read_blocks_;
        Genes                        genes_;
        RefTrees                     gtrees_;
        TagShards                    shards_;
        std::regex                   barcode_re_;
        std::regex                   random_hex_re_;
        hmap                         barcode_map_;
        hmap                         gene_map_;
        std::string                  bam_file_;
        std::string                  out_file_;
        std::string                  lib_string_;
        std::string                  barcode_regex_str_;
        std::string                  random_hex_regex_str_;
        std::string                  random_hex_value_;
        TagSummary                   full_;
        tsmap                        ssum_;

        std::vector<std::string>     samples_;
        smap                         sample_set_;
        char                         barcode_tag_[2];
        char                         umi_tag_[2];
        char                         sample_tag_[2];
        samFile                    * bf_ = nullptr;
        sam_hdr_t                  * bh_ = nullptr;
        double                       min_gene_ = 0.95;
        double                       min_exonic_ = 0.95;
        unsigned int                 min_gene_bases_ = 40;
        unsigned int                 min_intronic_ = 15;
        uint32_t                     barcode_length_ = 0;
        uint32_t                     umi_length_ = 0;
        StrandMode                   strand_ = StrandMode::TAG_UNKNOWN;
        uint8_t                      min_qual_ = 255;
        bool                         done_ = true;
        bool                         visium_hd_ = false;
        bool                         discard_unknown_juncs_ = false;
        bool                         has_probes_ = false;
};

unsigned int block_overlaps(const std::vector<Exon> & blocks, const std::vector<Exon> & exons);
unsigned int check_junctions(const std::vector<Exon> & blocks, const Gene::Juncs & junctions);
unsigned int get_read_blocks(bam1_t * b, std::vector<Exon> & blocks);
}
