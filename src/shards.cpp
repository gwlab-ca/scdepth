// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#include "shards.hpp"
#include "htslib/bgzf.h"
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>

using namespace scdepth;

void TagShards::collapse_(std::vector<UMITag>& counts) {
    size_t out = 0;
    for(size_t i = 0; i < counts.size();) {
        UMITag acc = counts[i];
        size_t j   = i + 1;
        while(j < counts.size() && acc.same_molecule(counts[j])) {
            acc.inc_spliced(counts[j].spliced);
            acc.inc_unspliced(counts[j].unspliced);
            acc.inc_ambiguous(counts[j].ambiguous);
            acc.flags |= counts[j].flags;
            j++;
        }
        counts[out++] = acc;
        i             = j;
    }
    counts.resize(out);
}

bool TagShards::merge_() {
    if(tags_.empty())
        return true;
    std::sort(tags_.begin(), tags_.end());
    auto start = tags_.begin();
    barcodes_.push_back(BarcodeHeader());
    barcodes_.back().make_header(start->barcode, 0, 0);

    while(start != tags_.end()) {
        merged_.push_back(UMITag());
        merged_.back().make_tag(start->gene, start->umi);
        barcodes_.back().block_size += sizeof(UMITag);
        barcodes_.back().n_count_tags++;
        if(start->unpack_spliced()) {
            merged_.back().inc_spliced(1);
        } else if(start->unpack_unspliced()) {
            merged_.back().inc_unspliced(1);
        } else {
            merged_.back().inc_ambiguous(1);
        }
        if(start->is_random_hex()) {
            merged_.back().set_random_hex(true);
        }
        auto it = std::next(start);
        while(it != tags_.end() && start->same_umi(*it)) {
            if(it->unpack_spliced()) {
                merged_.back().inc_spliced(1);
            } else if(it->unpack_unspliced()) {
                merged_.back().inc_unspliced(1);
            } else {
                merged_.back().inc_ambiguous(1);
            }
            it++;
        }
        if(it != tags_.end() && it->barcode != barcodes_.back().barcode) {
            barcodes_.push_back(BarcodeHeader());
            barcodes_.back().make_header(it->barcode, 0, 0);
        }
        start = it;
    }
    return true;
}

bool TagShards::write() {
    std::stringstream outss;
    if(!merge_()) {
        std::cerr << "[error] Error merge call failed unexpectedly\n";
        return false;
    }
    if(merged_.empty()) {
        return true;
    }
    outss << out_prefix_ << ".shard_tags." << std::setfill('0') << std::setw(4)
          << shard_count_ << ".gz";
    BGZF* bout = bgzf_open(outss.str().c_str(), "w");
    if(bout == nullptr) {
        std::cerr << "[error] Error opening shard tag file " << outss.str()
                  << "\n";
        return false;
    }
    shard_files_.push_back(ShardFile(outss.str()));
    std::vector<UMITag>::const_iterator mit = merged_.cbegin();
    for(auto& bc : barcodes_) {
        if((write_one_(bout, shard_files_.back().file, bc, mit) < 0)) {
            bgzf_close(bout);
            return false;
        }
    }

    if(mit != merged_.end()) {
        std::cerr << "[error] merged had leftover CountTags after final "
                     "barcode write\n";
        bgzf_close(bout);
        return false;
    }

    if(bgzf_flush(bout) < 0) {
        std::cerr << "[error] Error flushing the shard tag file " << outss.str()
                  << "\n";
        bgzf_close(bout);
        return false;
    }

    if(bgzf_close(bout) < 0) {
        std::cerr << "[error] Error closing the shard tag file " << outss.str()
                  << "\n";
        return false;
    }

    shard_count_++;
    tags_written_ += tags_.size();
    tags_.clear();
    barcodes_.clear();
    merged_.clear();
    return true;
}

long int TagShards::write_one_(BGZF* bout, const std::string& fname,
                               const BarcodeHeader&                 bc,
                               std::vector<UMITag>::const_iterator& mit) {
    if(bgzf_write(bout, reinterpret_cast<const void*>(&bc),
                  sizeof(BarcodeHeader)) <= 0) {
        std::cerr << "[error] Error writing to shard file " << fname << "\n";
        bgzf_close(bout);
        return -1;
    };

    long int written = 0;

    for(size_t i = 0; i < bc.n_count_tags; i++) {
        if(mit == merged_.end()) {
            std::cerr << "[error] Unexpected end of UMIs\n";
            bgzf_close(bout);
            return -1;
        }
        const auto& m = *mit;
        if(bgzf_write(bout, reinterpret_cast<const void*>(&(m)),
                      sizeof(UMITag)) <= 0) {
            std::cerr << "[error] Error writing to shard file " << fname
                      << "\n";
            bgzf_close(bout);
            return -1;
        };
        written += sizeof(UMITag);
        mit++;
    }
    if(written != bc.block_size) {
        std::cerr << "[error] block_size mismatch for barcode = " << bc.barcode
                  << " tag count = " << bc.n_count_tags
                  << " expected written = "
                  << (sizeof(UMITag) * bc.n_count_tags)
                  << " header = " << bc.block_size << " written = " << written
                  << "\n";
        bgzf_close(bout);
        return -1;
    }
    // std::cerr << "written = " << (written + sizeof(BarcodeHeader)) << "\n";
    return written + sizeof(BarcodeHeader);
}

long int TagShards::write_one_(BGZF* bout, const std::string& fname,
                               const BarcodeData& d) {
    std::vector<UMITag>::const_iterator mit = d.counts.cbegin();
    long int written = write_one_(bout, fname, d.bc, mit);
    if(written < 0) {
        return -1;
    }

    if(mit != d.counts.end()) {
        std::cerr << "Error writing record for barcode " << d.bc.barcode
                  << " count iterator is not at the end\n";
        return -1;
    }
    return written;
}

bool TagShards::merge_shards(std::vector<BarcodeCount>& bc_counts) {
    // check for in_memory
    bool in_memory = !tags_.empty();
    if(in_memory && !merge_())
        return false;
    std::cerr << "Merging " << size() << " shard files";
    if(in_memory)
        std::cerr << " and in memory tags\n";
    std::vector<UMITag>::const_iterator        mit = merged_.cbegin();
    std::vector<BarcodeHeader>::const_iterator bit = barcodes_.cbegin();
    uint32_t bc_min = std::numeric_limits<uint32_t>::max();
    raw_molecules_  = 0;

    std::string out_file = out_prefix_ + "_tags.gz";
    BGZF*       bout     = bgzf_open(out_file.c_str(), "w");
    if(bout == nullptr) {
        std::cerr << "[error] Error opening output file " << out_file << "\n";
        return false;
    }
    if(in_memory) {
        bc_min = barcodes_.front().barcode;
    }

    for(auto& f : shard_files_) {
        f.fin = bgzf_open(f.file.c_str(), "r");
        if(f.fin == nullptr) {
            std::cerr << "[error] Error opening shard file " << f.file << "\n";
            close_all_();
            return false;
        }
        auto ret = read_barcode(f.fin, f.next);
        if(ret <= 0) {
            std::cerr << "[error] Error reading barcode from shard file "
                      << f.file << "\n";
            close_all_();
            bgzf_close(bout);
            return false;
        }
        f.done = false;
        bc_min = std::min(bc_min, f.barcode());
    }

    if(bc_min == std::numeric_limits<uint32_t>::max()) {
        std::cerr << "[error] No barcode's found in shard files";
        close_all_();
        bgzf_close(bout);
        return false;
    }

    BarcodeData merged;
    while(bc_min < std::numeric_limits<uint32_t>::max()) {
        uint32_t bc_next = std::numeric_limits<uint32_t>::max();
        merged.counts.clear();
        merged.bc.n_count_tags = 0;
        merged.bc.block_size   = 0;
        merged.bc.barcode      = bc_min;
        // shard files first then in memory
        // std::cout << "merging barcode = " << bc_min << " detected = " <<
        // barcodes_detected << " / " << bc_counts.size() << "\n";
        for(auto& f : shard_files_) {
            if(f.done)
                continue;
            // std::cout << "  shard = " << bit->barcode << "\n";
            if(f.barcode() == bc_min) {
                // merge
                // merged.bc.n_count_tags += f.next.bc.n_count_tags;
                // merged.bc.block_size += f.next.bc.block_size;
                merged.counts.insert(merged.counts.end(), f.next.counts.begin(),
                                     f.next.counts.end());

                auto ret = read_barcode(f.fin, f.next);
                if(ret < 0) {
                    std::cerr
                        << "[error] Error reading barcode from shard tag file "
                        << f.file << "\n";
                    close_all_();
                    bgzf_close(bout);
                    return false;
                } else if(ret == 0) {
                    f.done = true;
                    f.close();
                } else {
                    bc_next = std::min(bc_next, f.barcode());
                }
            } else {
                bc_next = std::min(bc_next, f.barcode());
            }
        }

        // read in memory barcode
        if(in_memory && bit != barcodes_.end()) {
            // std::cout << "  in memory = " << bit->barcode << "\n";
            if(bit->barcode == bc_min) {
                // merged.bc.n_count_tags += bit->n_count_tags;
                // merged.bc.block_size += bit->block_size;
                for(size_t i = 0; i < bit->n_count_tags; i++) {
                    if(mit == merged_.end()) {
                        std::cerr << "[error] Error merging in memory barcode "
                                     "count (unexpected end)"
                                  << bit->barcode << "\n";
                        close_all_();
                        bgzf_close(bout);
                        return false;
                    }
                    merged.counts.push_back(*mit);
                    mit++;
                }
                bit++;
            }
            if(bit != barcodes_.end())
                bc_next = std::min(bc_next, bit->barcode);
        }
        std::sort(merged.counts.begin(), merged.counts.end());
        collapse_(merged.counts);
        merged.bc.n_count_tags = merged.counts.size();
        merged.bc.block_size   = merged.counts.size() * sizeof(UMITag);

        auto     offset  = bgzf_tell(bout);
        long int written = write_one_(bout, out_file, merged);
        // std::cerr << "written = " << written << "\n";
        if(written < 0) {
            std::cerr << "[error] Error writing merged record\n";
            close_all_();
            bgzf_close(bout);
            return false;
        }

        bc_counts[merged.bc.barcode].offset           = offset;
        bc_counts[merged.bc.barcode].total_data_bytes = written;
        bc_counts[merged.bc.barcode].raw_molecules    = merged.counts.size();
        raw_molecules_ += merged.counts.size();
        long int expected = merged.bc.block_size + sizeof(BarcodeHeader);
        if(written != expected) {
            std::cerr << "[error] Bytes mismatch written = " << written
                      << " expected = " << expected
                      << " block size = " << merged.bc.block_size << "\n";
            return false;
        }
        bc_min = bc_next;
        // std::cerr << "    bc_next = " << bc_next << "\n";
        // write the barcode setup index etc
        // check gene order
    }

    for(auto& b : bc_counts) {
        if(b.countable == 0) {
            b.offset = bgzf_tell(bout);
        }
    }

    bool error = false;
    if(in_memory) {
        if(bit != barcodes_.end()) {
            std::cerr
                << "[error] Error in memory barcodes were not fully merged\n";
            std::cerr << "  unmerged barcode = " << bit->barcode << "\n";
            error = true;
        }
        if(mit != merged_.end()) {
            std::cerr << "[error] Error in memory UMIs were not fully merged\n";
            error = true;
        }
    }

    for(auto& f : shard_files_) {
        if(!f.done) {
            error = true;
            std::cerr << "[error] Error shard " << f.file
                      << " barcodes were not fully merged\n";
        }
    }
    close_all_();

    if(error) {
        return false;
    }

    if(bgzf_flush(bout) < 0) {
        std::cerr << "[error] Error flushing " << out_file << "\n";
        bgzf_close(bout);
        return false;
    }

    if(bgzf_close(bout) < 0) {
        std::cerr << "[error] Error closing " << out_file << "\n";
        return false;
    }
    // clean up
    for(auto& f : shard_files_) {
        try {
            if(!std::filesystem::remove(f.file)) {
                std::cerr << "[error] Error removing shard file " << f.file
                          << "\n";
            }
        } catch(const std::filesystem::filesystem_error& err) {
            std::cout << "[error] File system error removing file " << f.file
                      << "  error = " << err.what() << '\n';
        }
    }
    /*
    for(size_t i = 0; i < bc_counts.size(); i++){
        const auto & b = bc_counts[i];
        std::cerr << "barcode = " << i << " total = " << b.total << " countable
    = " << b.countable << " offset = " << b.offset << "\n";

    }
    */

    return true;
}
