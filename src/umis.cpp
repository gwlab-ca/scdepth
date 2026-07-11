// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#include "umis.hpp"
#include "pcg_random.hpp"
#include <random>

using namespace scdepth;

struct DirectedEdges{
    inline void operator()(size_t i, size_t j, UMINode & n1, UMINode & n2) const noexcept{
        if(n1.total() >= (n2.total() * 2 - 1)){
            n1.edges.push_back(j);
        }
        if(n2.total() >= (n1.total() * 2 - 1)){
            n2.edges.push_back(i);
        }
    }
};

void UMIDirectional::downsample(double frac, uint32_t umi_len,
        const BarcodeCount & bc, const BarcodeData & bd, UMIResults & out,
        bool directional, bool correct_multi_umis, PrimerMode primer_mode){ //, bool profile_singles){

    umi_len_ = umi_len;
    pcg64 rng(bc.seed);
    using Binom = std::binomial_distribution<uint32_t>;
    Binom binom;
    auto it = bd.counts.begin();
    /*
    size_t tot = 0, otot = 0;
    */
    //std::cout << "  Downsampling barcode " << bc.barcode << " frac = " << frac << " raw mols = " << bc.raw_molecules  << "\n";
    //std::cout << "  check counts = " << bd.counts.size() << "\n";
    out.umis.clear();
    out.dis_reads = 0;
    out.dis_ambig = 0;
    out.dis_mols = 0;

    //size_t lost = 0, processed = 0;
    while(it != bd.counts.end()){
        auto start = it;
        //std::cout << "    downsampling gene " << start->gene << " index = " << std::distance(bd.counts.begin(), start) << " / " << bd.counts.size() << "\n";
        n_nodes_ = 0;
        //s  nodes = 0;
        while(it != bd.counts.end() && start->gene == it->gene){
            if((primer_mode == PrimerMode::PolyAOnly && it->is_random_hex()) || 
                    (primer_mode == PrimerMode::RandomHexOnly && !it->is_random_hex()))
            {
               it++;
               continue;
            }

            uint32_t spliced = 0, unspliced = 0, ambiguous = 0;
            if(frac < 1.0){
                if(it->spliced > 0){
                    Binom::param_type params(it->spliced, frac);
                    binom.param(params);
                    spliced = binom(rng);
                }
                if(it->unspliced > 0){
                    Binom::param_type params(it->unspliced, frac);
                    binom.param(params);
                    unspliced = binom(rng);
                }
                if(it->ambiguous > 0){
                    Binom::param_type params(it->ambiguous, frac);
                    binom.param(params);
                    ambiguous = binom(rng);
                }
            }else{
                spliced = it->spliced;
                unspliced = it->unspliced;
                ambiguous = it->ambiguous;
            }
            if((spliced + ambiguous + unspliced) > 0){
                if(n_nodes_ >= nodes_.size()) nodes_.resize(n_nodes_ + 5);
                nodes_[n_nodes_].spliced = spliced;
                nodes_[n_nodes_].unspliced = unspliced;
                nodes_[n_nodes_].ambiguous = ambiguous;
                nodes_[n_nodes_].flags |= it->flags;
                nodes_[n_nodes_].umi = it->umi;
                nodes_[n_nodes_].edges.clear();
                n_nodes_++;
            }
            /*std::cout << "    n_nodes_ = " << n_nodes_ << " gene = " << start->gene 
                << " umi = " << int2seq(it->umi, umi_len_)
                << " spliced = " << spliced << " unspliced = " << unspliced << " ambiguous = " << ambiguous << "\n";*/
            it++;
        }

        //std::cout << "    downsampled gene " << start->gene << " number of nodes = " << n_nodes_  << "\n";
        if(n_nodes_ == 0) continue;
        if(n_nodes_ < 500){
            if(directional) build_graph_brute_(DirectedEdges{});
        }else{
            if(directional) build_graph_hash_(DirectedEdges{});
        }
        correct_umis_(start->gene, out);
        //lost += (n_nodes_ - out.umis.size());
        //std::cout << "    after gene = " << start->gene << " raw = " << n_nodes_ 
        //    << " corrected = " << out.umis.size() << " lost so far = " << lost << "\n";
        /*
        if(final_umis_.size() != n_nodes_){
            std::cout << "    after gene = " << start->gene << " raw = " << n_nodes_ 
                << " corrected = " << final_umis_.size() << " lost so far = " << lost << "\n";
        }
        lost += (n_nodes_ - final_umis_.size());
        processed += n_nodes_;
        */
        //out.umis.insert(out.umis.end(), final_umis_.begin(), final_umis_.end());
    }


    //std::cout << "correct multi umis\n";
    //std::cout << "    mols before umi correction = " << out.umis.size() << "\n";
    if(correct_multi_umis){
        correct_multi_umis_(out); //, profile_singles);
    }
    /*
    size_t dcmols = 0, rlost = 0; 
    for(auto & c : out.umis){
        if(!c.invalid()){
            dcmols++;
        }else{
            rlost += c.total();
        }
    }
    std::cout << "    mols after umi correction = " << dcmols << "\n";
    */
}

template <typename T>
void UMIDirectional::build_graph_brute_(T build_edge){
    static constexpr uint32_t mask = 0x55555555u;
    for(size_t i = 0; i < n_nodes_; i++){
        auto & n1 = nodes_[i];
        for(size_t j = i + 1; j < n_nodes_; j++){
            auto & n2 = nodes_[j];
            //no links between different primers
            if(n1.is_random_hex() != n2.is_random_hex()) continue;
            auto mm = n1.umi ^ n2.umi;
            auto cnt = __builtin_popcount((mm | (mm >> 1)) & mask);
            if(cnt != 1) continue;
            build_edge(i, j, n1, n2);
        }
    }
}

template <typename T>
void UMIDirectional::build_graph_hash_(T build_edge){
    umi_map_.clear();
    umi_map_.reserve(2 * n_nodes_);
    for(size_t i = 0; i < n_nodes_; i++){
        umi_map_[nodes_[i].umi] = i;
    }
    for(size_t i = 0; i < n_nodes_; i++){
        auto & n1 = nodes_[i];
        for(uint32_t p = 0; p < umi_len_; p++){
            uint32_t s = 2 * p;
            uint32_t m = ~(0b11u << s); 
            uint32_t b = (n1.umi >> s) & 0b11u;
            for(uint32_t nb = 0; nb < 4; nb++){
                if(b == nb) continue;
                uint32_t nu = (n1.umi & m) | (nb << s);
                auto it = umi_map_.find(nu);
                if(it == umi_map_.end()) continue;
                auto j = it->second;
                if(j <= i) continue;
                auto & n2 = nodes_[j];
                if(n1.is_random_hex() != n2.is_random_hex()) continue;
                build_edge(i, j, n1, n2);
            }
        }
    }
}

void UMIDirectional::correct_umis_(uint32_t gene, UMIResults & out){
    idx_.resize(n_nodes_);
    std::iota(idx_.begin(), idx_.end(), 0);
    std::sort(idx_.begin(), idx_.end(),
        [this](uint32_t a, uint32_t b) {
            if (nodes_[a].total() != nodes_[b].total())
                return nodes_[a].total() > nodes_[b].total();
            return nodes_[a].umi < nodes_[b].umi;
        }
    );

    visited_.assign(n_nodes_, -1);
    comps_full_.clear();
    stack_.clear();
    for(auto i : idx_){
        if(visited_[i] != -1) continue;
        size_t cid = comps_full_.size();
        comps_full_.push_back();
        auto & c = comps_full_.back();
        stack_.push_back(i);
        visited_[i] = cid;
        while(!stack_.empty()){
            uint32_t n = stack_.back();
            stack_.pop_back();
            c.push_back(n);
            for(auto u : nodes_[n].edges){
                if(visited_[u] == -1){
                    visited_[u] = cid;
                    stack_.push_back(u);
                }
            }
            /*
            if(weakly_){
                for(auto u : rev_[n]){
                    if(visited_[u] == -1){
                        visited_[u] = cid;
                        stack_.push_back(u);
                    }
                }
            }
            */
        }
    }

    size_t start = out.umis.size();
    out.umis.resize(start + comps_full_.size());
    //if(n_nodes_ >= 100) std::cout << "  umi components for gene = " << gene << " = " << comps_full_.size() << "\n";
    for(size_t i = 0; i < comps_full_.size(); i++){
        auto & c = comps_full_[i];
        if(c.empty()) continue;
        //if(n_nodes_ >= 100) std::cout << "    component = " << i << " size = " << c.size() << "\n";
        std::sort(c.begin(), c.end(),
            [this](uint32_t a, uint32_t b) {
                if (nodes_[a].total() != nodes_[b].total())
                    return nodes_[a].total() > nodes_[b].total();
                return nodes_[a].umi < nodes_[b].umi;
            }
        );

        size_t idx = i + start;
        out.umis[idx].umi = nodes_[c.front()].umi;
        out.umis[idx].flags = 0;
        out.umis[idx].gene = gene;
        out.umis[idx].spliced = 0;
        out.umis[idx].unspliced = 0;
        out.umis[idx].ambiguous = 0;
        out.umis[idx].n_umis = 0;
        for(auto j : c){
            out.umis[idx].inc_spliced(nodes_[j].spliced);
            out.umis[idx].inc_unspliced(nodes_[j].unspliced);
            out.umis[idx].inc_ambiguous(nodes_[j].ambiguous);
            out.umis[idx].inc_n_umis();
            out.umis[idx].flags |= nodes_[j].flags;
        }

    }
}

void UMIDirectional::correct_multi_umis_(UMIResults & out){ //, bool profile_singles){
    std::sort(out.umis.begin(), out.umis.end(), 
        [](const CorrUMITag & c1, const CorrUMITag & c2){
            if(c1.umi == c2.umi) return c1.total() > c2.total();
            return c1.umi < c2.umi;
        }
    );

    /*
    std::cout << "  correcting for ambiguous umis = " << out.umis.size() << "\n";
    for(auto & u : out.umis){
        std::cout << "    corrected umi = " << u.umi << " count = " << u.count << "\n";
    }
    */
    auto start = out.umis.begin();
    while(start != out.umis.end()){
        auto it = std::next(start);
        bool all_bad = false, has_bad = false;
        while(it != out.umis.end() && start->umi == it->umi){
            if(it->total() == start->total()) all_bad = true;
            out.dis_reads += it->total();
            out.dis_mols++;
            has_bad = true;
            it->set_invalid();
            it++;
        }
        if(all_bad){
            out.dis_reads += start->total();  // only if you define ties as “discard all”
            out.dis_mols++;
            start->set_invalid();
        }
        if(has_bad) out.dis_ambig++;
        start = it;
    }

    std::sort(out.umis.begin(), out.umis.end());
}
