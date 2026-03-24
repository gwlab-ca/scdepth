// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#pragma once

#include <vector>
#include <cstdint>
#include "tags.hpp"
#include "gtf.hpp"

namespace scdepth{

class TagShards{
    public:
        TagShards()
        {
            reset();
        }

        void reset(){
            shard_count_ = 0;
            tid_ = 0;
            bidx_ = 0;
            cutoff_ = max_elems_ * write_frac_;
            merged_.clear();
            tags_.clear();
            tags_.reserve(max_elems_);
            tags_written_ = 0;
        }

        void init(const std::string & out_prefix, size_t max_elems, double write_frac){
            out_prefix_ = out_prefix;
            max_elems_ = max_elems;
            write_frac_ = write_frac;
        }

        bool write();
        bool merge_shards(std::vector<BarcodeCount> & bc_counts);

        void add_tag(int tid, int lft, const RawTag & tag){
            bool past_break = update_break_(tid, lft);
            if(tags_.size() >= cutoff_ && past_break){
                write();
            }
            tags_.push_back(tag);
        }

        size_t tags_written() const {
            return tags_written_;
        }

        size_t tags_cached() const {
            return tags_.size();
        }

        GeneBreaks & gene_breaks() {
            return breaks_;
        }

        const GeneBreaks & gene_breaks() const {
            return breaks_;
        }

        size_t size() const {
            return shard_files_.size();
        }

        size_t raw_molecules() const {
            return raw_molecules_;
        }

    private:
        struct ShardFile{
            ShardFile(const std::string & file) : file(file){
            }

            ~ShardFile(){
                close();
            }

            void close() {
                if(fin != nullptr) bgzf_close(fin);
                fin = nullptr;
            }

            uint32_t barcode() const {
                return next.bc.barcode;
            }

            BGZF                * fin = nullptr;
            std::string           file;
            BarcodeData           next;
            bool                  done = true;
        };

        void collapse_(std::vector<UMITag> & counts);
        bool merge_();

        void close_all_(){
            for(auto & f : shard_files_){
                if(f.fin != nullptr) bgzf_close(f.fin);
                f.fin = nullptr;
            }
        }

        long int write_one_(BGZF * bout, const std::string & fname, const BarcodeHeader & bc,
            std::vector<UMITag>::const_iterator & mit);

        long int write_one_(BGZF * bout, const std::string & fname, const BarcodeData & d);

        bool update_break_(int tid, int lft){
            bool past_break = false;
            if(tid != tid_){
                tid_ = tid;
                bidx_ = 0;
                past_break = true;
            }
            while(bidx_ < breaks_[tid].size() && static_cast<int>(breaks_[tid][bidx_]) < lft){
                past_break = true;
                bidx_++;
            }
            return past_break;
        }

        std::vector<RawTag>          tags_;
        std::vector<BarcodeHeader>   barcodes_;
        std::vector<UMITag>          merged_;
        std::vector<ShardFile>       shard_files_;
        std::vector<size_t>          sort_idx_;
        std::vector<UMITag>          sort_counts_;
        std::string                  out_prefix_;
        GeneBreaks                   breaks_;
        double                       write_frac_ = 0.95;
        size_t                       max_elems_ = 10000000;
        size_t                       cutoff_ = 0;
        size_t                       raw_molecules_ = 0;
        size_t                       shard_count_ = 0;
        size_t                       tags_written_ = 0;
        int                          tid_ = -1;
        uint32_t                     bidx_ = 0;
};

}
