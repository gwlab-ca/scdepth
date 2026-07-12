// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#include "dsthread.hpp"
#include "htslib/bgzf.h"
#include "downsampler.hpp"
#include <iostream>

using namespace scdepth;

DownsamplerThread::DownsamplerThread(Downsampler & data, uint32_t umi_len, 
        size_t genes, size_t chunk_start, size_t chunk_end) 
    : data_(data), cstart_(chunk_start), cend_(chunk_end), umi_len_(umi_len)
{
    out.reset(data.output.fracs, data.output.barcodes, genes, data.output.max_hist, 
            data.output.build_mats, data.has_visium, data.output.calc_sau,
            data.aggregate_only, data.samples);
}

DownsamplerThread::~DownsamplerThread(){
    if (thread_.joinable()) {
        thread_.join();
    }
}

void DownsamplerThread::start(){
    thread_ = std::thread(std::ref(*this));
}

void DownsamplerThread::join() {
    if(thread_.joinable())
        thread_.join();
}

void DownsamplerThread::operator()(){
    std::string tagf = data_.prefix + "_tags.gz";
    BGZF * fin = bgzf_open(tagf.c_str(), "r");
    if(fin == nullptr){
        std::cerr << "[error] Error opening tag file " << tagf << "\n";
        error = true;
        return;
    }

    UMIDirectional corrector;
    UMIResults     bco;
    for(size_t i = cstart_; i < cend_; i++){
        auto & bc = data_.barcodes[i];
        if(bc.countable == 0) continue;
        if(bgzf_seek(fin, data_.barcodes[i].offset, SEEK_SET) < 0){
            std::cerr << "[error] Error seeking to " << data_.barcodes[cstart_].offset << " for " << tagf << "\n";
            error = true;
            bgzf_close(fin);
            return;
        }
        break;
    }

    //std::cout << " starting thread " << cstart_ << " - " << cend_ << " use qc = " << data_.use_qc << "\n";
    //Find the first valid barcode to seek to
    for(size_t i = cstart_; i < cend_; i++){
        //get the first barcode
        auto & bc = data_.barcodes[i];
        processed++;
        if(bc.countable == 0){
            if(out.build_mats){
                for(size_t i = 0; i < data_.output.fracs.size(); i++){
                    if(out.calc_sau){
                        out.spliced_mat[i].indptr.push_back(out.spliced_mat[i].data.size());
                        out.unspliced_mat[i].indptr.push_back(out.unspliced_mat[i].data.size());
                        out.ambiguous_mat[i].indptr.push_back(out.ambiguous_mat[i].data.size());
                        if(out.spliced_mat[i].indptr.size() != (bc.index - cstart_ + 2)){
                            std::cerr << "[error] indptr size = " << out.spliced_mat[i].indptr.size() << " c = " << i << " cstart = " << cstart_ << " index = " << bc.index << "\n";
                        }
                    }
                    out.total_mat[i].indptr.push_back(out.total_mat[i].data.size());
                    if(out.total_mat[i].indptr.size() != (bc.index - cstart_ + 2)){
                        std::cerr << "[error] indptr size = " << out.total_mat[i].indptr.size() << " c = " << i << " cstart = " << cstart_ << " index = " << bc.index << "\n";
                    }
                }
            }

            if(!data_.aggregate_only){
                for(size_t i = 0; i < data_.output.fracs.size(); i++){
                    out.gcounts[i].indptr.push_back(out.gcounts[i].indices.size());
                }
            }
            continue;
        }
        if(read_barcode(fin, b_) <= 0){
            std::cerr << "[error] Error reading barcode c = " << i << " (" << cstart_ << " - " << cend_ << ") "
                << bc.barcode << " index = " << bc.index << " bgzf_tell = " << bgzf_tell(fin) 
                << " offset = " << bc.offset << "\n";
            error = true;
            bgzf_close(fin);
            return;
        }

        //std::cout << "Barcode read" << std::endl;
        /*
        std::cout << "barcode idx = " << b_.bc.barcode << "\n";
        std::cout << "    raw counts = " << b_.counts.size() << "\n";
        for(auto & c : b_.counts){
            std::cerr << "    tag  gene = " << c.gene << " reads = " << c.total() 
                << " umi = " << c.umi << " " << int2seq(c.umi, umi_len_) << "\n";
        }
        */

        auto & o = data_.output;
        for(size_t i = 0; i < o.fracs.size(); i++){
            corrector.downsample(o.fracs[i], umi_len_, bc, b_, bco,
                    data_.umi_directed, data_.umi_multi, data_.primer_mode); //, data_.umi_singletons);

            //o.total_singles[i] += bco.total_singles;
            //o.ed0_cross[i] += bco.ed0_cross;
            //o.ed1_cross[i] += bco.ed1_cross;
            //o.ed1_same[i] += bco.ed1_same;
            merge_barcode_(i, bc, bco);

            /*
            std::cout << "  barcode = " << bc.index 
                << " spliced indptr size = " << out.spliced_mat[i].indptr.size()
                << " unspliced indptr size = " << out.unspliced_mat[i].indptr.size()
                << " ambiguous indptr size = " << out.ambiguous_mat[i].indptr.size()
                << "\n";
            */
            /*
            if(out.spliced_mat[i].indptr.size() != (bc.index - cstart_ + 1)){
                std::cerr << "[error] indptr size = " << out.spliced_mat[i].indptr.size() << " c = " << c << " cstart = " << cstart_ << " index = " << bc.index << "\n";
            }
            if(out.unspliced_mat[i].indptr.size() != (bc.index - cstart_ + 1)){
                std::cerr << "[error] indptr size = " << out.unspliced_mat[i].indptr.size() << " c = " << c << " cstart = " << cstart_ << " index = " << bc.index << "\n";
            }
            if(out.ambiguous_mat[i].indptr.size() != (bc.index - cstart_ + 1)){
                std::cerr << "[error] indptr size = " << out.ambiguous_mat[i].indptr.size() << " c = " << c << " cstart = " << cstart_ << " index = " << bc.index << "\n";
            }
            */
        }

    }
}

void DownsamplerThread::merge_barcode_(size_t i, const BarcodeCount & bc, const UMIResults & bco){
    //std::cout << "don't forget to keep track of the number of reads lost due to ambiguous umi/genes\n";
    //merge the stuff
    size_t hidx = out.max_hist * i;
    size_t bidx = out.barcodes * i + bc.index;
    bool bmat = out.build_mats;

    uint32_t bc_excluded_reads = 0;
    //bc_mhist_.assign(data.max_hist, 0);
    auto it = bco.umis.begin();
    uint32_t bc_spliced_mols = 0, bc_unspliced_mols = 0, bc_ambiguous_mols = 0, bc_total_mols = 0;
    uint32_t bc_spliced_reads = 0, bc_unspliced_reads = 0, bc_ambiguous_reads = 0, bc_total_reads = 0;
    while(it != bco.umis.end()){
        uint32_t spliced_mols = 0, unspliced_mols = 0, ambiguous_mols = 0, total_mols = 0;
        uint32_t spliced_reads = 0, unspliced_reads = 0, ambiguous_reads = 0, total_reads = 0;
        auto start = it;
        bool is_mt = data_.has_mt && data_.mt_filter[start->gene];
        bool is_mod = data_.has_mod && data_.mod_filter[start->gene];

        while(it != bco.umis.end() && start->gene == it->gene){
            auto & u = *it;
            // build matrices
            if(u.invalid()){
                it++;
                continue;
            }
            if(data_.exclude_filter[it->gene]){
                bc_excluded_reads += u.total();
                it++;
                continue;
            }

            uint32_t ut = u.total();
            uint32_t mc = std::min(ut, out.max_hist);
            if(mc == 0) continue;
            //1 will be in index 0
            mc--;
            uint32_t m = mc + hidx;
            /*
            if(u.n_umis > 1) {
                bc_merged_umis++;
                bc_extra_umis += uint32_t(u.n_umis - 1);
            }
            if(ut == 1) bc_singleton_mols++;
            */

            total_reads += ut;
            total_mols++;
            out.total_mhist[m]++;
            if(data_.samples > 0){
                uint32_t sample = data_.barcode2sample[bc.index];
                size_t idx = (static_cast<size_t>(sample) * out.steps + i) * out.max_hist + mc;
                out.total_sample_mhist[idx]++;
                if(u.spliced > 0){
                    out.spliced_sample_mhist[idx]++;
                }else if(u.unspliced > 0){
                    out.unspliced_sample_mhist[idx]++;
                }else{
                    out.ambiguous_sample_mhist[idx]++;
                }
            }
            //bc_mhist_[mc]++;
            if(u.spliced > 0){
                out.spliced_mhist[m]++;
                spliced_mols++;
                spliced_reads += ut;
            }else if(u.unspliced > 0){
                out.unspliced_mhist[m]++;
                unspliced_mols++;
                unspliced_reads += ut;
            }else{
                out.ambiguous_mhist[m]++;
                ambiguous_mols++;
                ambiguous_reads += ut;
            }
            it++;
        }

        //update barcode totals
        bc_spliced_mols += spliced_mols;
        bc_unspliced_mols += unspliced_mols;
        bc_ambiguous_mols += ambiguous_mols;
        bc_total_mols += total_mols;
        bc_spliced_reads += spliced_reads;
        bc_unspliced_reads += unspliced_reads;
        bc_ambiguous_reads += ambiguous_reads;
        bc_total_reads += total_reads;

        //stuff that depends on a given gene
        if(!data_.aggregate_only){
            data_.output.total.bc_genes[bidx] += (total_mols > 0);

            if(is_mt){
                data_.output.total.bc_mt_mols[bidx] += total_mols;
                if(data_.output.calc_sau){
                    data_.output.spliced.bc_mt_mols[bidx] += spliced_mols;
                    data_.output.unspliced.bc_mt_mols[bidx] += unspliced_mols;
                    data_.output.ambiguous.bc_mt_mols[bidx] += ambiguous_mols;
                }
            }
            if(is_mod){
                data_.output.total.bc_mod_mols[bidx] += total_mols;
                if(data_.output.calc_sau){
                    data_.output.spliced.bc_mod_mols[bidx] += spliced_mols;
                    data_.output.unspliced.bc_mod_mols[bidx] += unspliced_mols;
                    data_.output.ambiguous.bc_mod_mols[bidx] += ambiguous_mols;
                }
            }
            if(data_.output.calc_sau){
                data_.output.spliced.bc_genes[bidx] += (spliced_mols > 0);
                data_.output.unspliced.bc_genes[bidx] += (unspliced_mols > 0);
                data_.output.ambiguous.bc_genes[bidx] += (ambiguous_mols > 0);
            }
            //gene matrix tracking
            if(total_mols > 0){
                out.gcounts[i].indices.push_back(start->gene);
            }
        }

        if(bmat){
            if(out.calc_sau){
                if(spliced_mols > 0){
                    out.spliced_mat[i].indices.push_back(start->gene);
                    out.spliced_mat[i].data.push_back(spliced_mols);
                }
                if(unspliced_mols > 0){
                    out.unspliced_mat[i].indices.push_back(start->gene);
                    out.unspliced_mat[i].data.push_back(unspliced_mols);
                }
                if(ambiguous_mols > 0){
                    out.ambiguous_mat[i].indices.push_back(start->gene);
                    out.ambiguous_mat[i].data.push_back(ambiguous_mols);
                }
            }
            out.total_mat[i].indices.push_back(start->gene);
            out.total_mat[i].data.push_back(total_mols);
        }
    }

    out.reads_discarded[i] += bco.dis_reads;
    out.mols_discarded[i] += bco.dis_mols;
    out.reads_excluded[i] += bc_excluded_reads;
    if(data_.samples > 0){
        uint32_t sample = data_.barcode2sample[bc.index];
        size_t idx = (static_cast<size_t>(sample) * out.steps + i);
        out.sample_reads_discarded[idx] += bco.dis_reads;
        out.sample_mols_discarded[idx] += bco.dis_mols;
        out.sample_reads_excluded[idx] += bc_excluded_reads;
    }
    if(bc_total_mols > 0){
        out.total_reads[i] += bc_total_reads;
        out.total_molecules[i] += bc_total_mols;
        if(!data_.aggregate_only){
            //data_.output.set_cell_mhist(i, bc.index, bc_mhist_);
            data_.output.total.bc_reads[bidx] += bc_total_reads;
            data_.output.total.bc_mols[bidx] += bc_total_mols;
        }
        out.spliced_reads[i] += bc_spliced_reads;
        out.spliced_molecules[i] += bc_spliced_mols;
        out.unspliced_reads[i] += bc_unspliced_reads;
        out.unspliced_molecules[i] += bc_unspliced_mols;
        out.ambiguous_reads[i] += bc_ambiguous_reads;
        out.ambiguous_molecules[i] += bc_ambiguous_mols;
        if(data_.samples > 0){
            uint32_t sample = data_.barcode2sample[bc.index];
            size_t idx = (static_cast<size_t>(sample) * out.steps + i);
            out.sample_total_reads[idx] += bc_total_reads;
            out.sample_total_molecules[idx] += bc_total_mols;
            out.sample_spliced_reads[idx] += bc_spliced_reads;
            out.sample_spliced_molecules[idx] += bc_spliced_mols;
            out.sample_unspliced_reads[idx] += bc_unspliced_reads;
            out.sample_unspliced_molecules[idx] += bc_unspliced_mols;
            out.sample_ambiguous_reads[idx] += bc_ambiguous_reads;
            out.sample_ambiguous_molecules[idx] += bc_ambiguous_mols;
        }
    }

    if(data_.output.calc_sau){
        if(!data_.aggregate_only){
            data_.output.spliced.bc_reads[bidx] += bc_spliced_reads;
            data_.output.spliced.bc_mols[bidx] += bc_spliced_mols;

            data_.output.unspliced.bc_reads[bidx] += bc_unspliced_reads;
            data_.output.unspliced.bc_mols[bidx] += bc_unspliced_mols;

            data_.output.ambiguous.bc_reads[bidx] += bc_ambiguous_reads;
            data_.output.ambiguous.bc_mols[bidx] += bc_ambiguous_mols;
        }
    }

    if(!data_.aggregate_only){
        data_.output.bc_dis_reads[bidx] += bco.dis_reads;
        data_.output.bc_dis_mols[bidx] += bco.dis_mols;
        data_.output.bc_exc_reads[bidx] += bc_excluded_reads;
        out.gcounts[i].indptr.push_back(out.gcounts[i].indices.size());
    }

    //finish the matrices for this barcode
    if(bmat){
        if(out.calc_sau){
            out.spliced_mat[i].indptr.push_back(out.spliced_mat[i].data.size());
            out.unspliced_mat[i].indptr.push_back(out.unspliced_mat[i].data.size());
            out.ambiguous_mat[i].indptr.push_back(out.ambiguous_mat[i].data.size());
        }
        out.total_mat[i].indptr.push_back(out.total_mat[i].data.size());
    }
}
