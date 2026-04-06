// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#include "dsdata.hpp"
#include "dsmisc.hpp"
#include "downsampler.hpp"
#include "htslib/bgzf.h"
#include <iostream>

using namespace scdepth;


void SparseMatrix::merge(const SparseMatrix & rhs){
    if(rhs.empty()) return;
    if(empty()){
        indptr = rhs.indptr;
        indices = rhs.indices;
        data = rhs.data;
        //std::cout << "Assign sparse indptr = " << indptr.size() << " data = " << data.size() << " indices = " << indices.size() << " back = " << indptr.back() << "\n";
        return;
    }

    /*
    std::cout << "Merge sparse indptr = " << indptr.size() << " data = " << data.size() << " indices = " << indices.size() << " back = " << indptr.back() << "\n";
    std::cout << "         rhs indptr = " << rhs.indptr.size() << " data = " << rhs.data.size() << " indices = " << rhs.indices.size() << " back = " << rhs.indptr.back() << "\n";
    */
    assert(indptr.back() == data.size());
    assert(rhs.indptr.back() == rhs.data.size());
    indices.reserve(indices.size() + rhs.indices.size());
    data.reserve(data.size() + rhs.data.size());
    //both matrices have + 1 elements
    indptr.reserve(indptr.size() + rhs.indptr.size() - 1);

    size_t osize = indptr.size();
    uint32_t delta = indptr.back();
    indices.insert(indices.end(), rhs.indices.begin(), rhs.indices.end());
    data.insert(data.end(), rhs.data.begin(), rhs.data.end());
    //skip the first one becuase indptr is always N+1 elements
    indptr.insert(indptr.end(), std::next(rhs.indptr.begin()), rhs.indptr.end());
    //move one after the first element from rhs
    //      |last from this| +1 from this/rhs.indptr[0] | rhs indptr[1] |
    //index | osize  - 2   |            osize - 1       | osize         |
    //must add delta to indptr[1:] 
    for(auto it = std::next(indptr.begin(), osize); it != indptr.end(); it++){
        (*it) += delta;
    }
    //std::cout << "         after indptr = " << indptr.size() << " data = " << data.size() << " indices = " << indices.size() << " back = " << indptr.back() << "\n";
}

void GeneCountMatrix::merge(const GeneCountMatrix & rhs){
    if(rhs.empty()) return;
    if(empty()){
        indptr = rhs.indptr;
        indices = rhs.indices;
        //std::cout << "Assign sparse indptr = " << indptr.size() << " indices = " << indices.size() << " back = " << indptr.back() << "\n";
        return;
    }

    //std::cout << "Merge sparse indptr = " << indptr.size() << " indices = " << indices.size() << " back = " << indptr.back() << "\n";
    //std::cout << "         rhs indptr = " << rhs.indptr.size() << " indices = " << rhs.indices.size() << " back = " << rhs.indptr.back() << "\n";
    indices.reserve(indices.size() + rhs.indices.size());
    //both matrices have + 1 elements
    indptr.reserve(indptr.size() + rhs.indptr.size() - 1);

    size_t osize = indptr.size();
    uint32_t delta = indptr.back();
    indices.insert(indices.end(), rhs.indices.begin(), rhs.indices.end());
    //skip the first one becuase indptr is always N+1 elements
    indptr.insert(indptr.end(), std::next(rhs.indptr.begin()), rhs.indptr.end());
    //move one after the first element from rhs
    //      |last from this| +1 from this/rhs.indptr[0] | rhs indptr[1] |
    //index | osize  - 2   |            osize - 1       | osize         |
    //must add delta to indptr[1:] 
    for(auto it = std::next(indptr.begin(), osize); it != indptr.end(); it++){
        (*it) += delta;
    }
    //std::cout << "         after indptr = " << indptr.size() << " indices = " << indices.size() << " back = " << indptr.back() << "\n";
}

void DownsampleResultsType::reset(const std::vector<double> & fracs, size_t barcodes, size_t max_hist, 
        /*size_t genes,*/ bool alloc_barcodes, bool has_mt, bool has_mod)
{
    mhist.assign(fracs.size() * max_hist, 0);
    reads.assign(fracs.size(), 0);
    molecules.assign(fracs.size(), 0);

    if(alloc_barcodes){
        size_t N = barcodes * fracs.size();
        bc_reads.assign(N, 0);
        bc_mols.assign(N, 0);
        bc_genes.assign(N, 0);

        if(has_mt)  bc_mt_mols.assign(N, 0);
        else        bc_mt_mols.resize(0);
        if(has_mod) bc_mod_mols.assign(N, 0);
        else        bc_mod_mols.resize(0);
    }else{
        bc_reads.resize(0);
        bc_mols.resize(0);
        bc_genes.resize(0);
        bc_mt_mols.resize(0);
        bc_mod_mols.resize(0);
    }

}

void DownsampleResultsType::merge(const DownsampleResultsType & rhs){
    sum_vectors(mhist, rhs.mhist, false);
    sum_vectors(reads, rhs.reads, false);
    sum_vectors(molecules, rhs.molecules, false);
}


void DownsampleResultsLocal::reset(const std::vector<double> & fracs, size_t barcodes, size_t genes, size_t max_hist, 
        bool build_mats, bool has_visium, bool calc_sau, bool aggregate_only){

    this->fracs = fracs;
    this->steps = fracs.size();
    this->barcodes = barcodes;
    this->max_hist = max_hist;
    this->genes = genes;

    reads_discarded.resize(fracs.size());
    reads_excluded.resize(fracs.size());
    mols_discarded.resize(fracs.size());

    spliced_reads.resize(fracs.size());
    unspliced_reads.resize(fracs.size());
    ambiguous_reads.resize(fracs.size());

    spliced_molecules.resize(fracs.size());
    unspliced_molecules.resize(fracs.size());
    ambiguous_molecules.resize(fracs.size());

    spliced_mhist.resize(fracs.size() * max_hist);
    unspliced_mhist.resize(fracs.size() * max_hist);
    ambiguous_mhist.resize(fracs.size() * max_hist);

    total_reads.resize(fracs.size());
    total_molecules.resize(fracs.size());
    total_mhist.resize(fracs.size() * max_hist);

    this->calc_sau = calc_sau;
    this->build_mats = build_mats;
    this->has_visium = has_visium;

    if(!aggregate_only){
        gcounts.assign(fracs.size(), GeneCountMatrix());
    }else{
        gcounts.resize(0);
    }

    if(build_mats){
        total_mat.assign(fracs.size(), SparseMatrix());
        if(calc_sau){
            spliced_mat.assign(fracs.size(), SparseMatrix());
            ambiguous_mat.assign(fracs.size(), SparseMatrix());
            unspliced_mat.assign(fracs.size(), SparseMatrix());
        }
    }else{
        total_mat.resize(0);
        spliced_mat.resize(0);
        ambiguous_mat.resize(0);
        unspliced_mat.resize(0);
    }
}

void DownsampleResults::reset(const std::vector<double> & fracs, size_t barcodes, size_t genes, size_t max_hist, 
        bool build_mats, bool calc_sau, bool alloc_barcodes, bool has_visium){

    this->fracs = fracs;
    this->steps = fracs.size();
    this->barcodes = barcodes;
    this->max_hist = max_hist;
    this->genes = genes;

    reads_discarded.resize(fracs.size());
    reads_excluded.resize(fracs.size());
    mols_discarded.resize(fracs.size());
    mols_ambig.resize(fracs.size());

    //total_singles.resize(fracs.size());
    //ed0_cross.resize(fracs.size());
    //ed1_cross.resize(fracs.size());
    //ed1_same.resize(fracs.size());

    if(alloc_barcodes){
        size_t N = barcodes * fracs.size();
        bc_dis_reads.assign(N, 0);
        bc_dis_mols.assign(N, 0);
        bc_exc_reads.assign(N, 0);
        //cell_mhist.assign(N * max_hist, 0);
    }else{
        bc_exc_reads.resize(0);
        bc_dis_reads.resize(0);
        bc_dis_mols.resize(0);
        //cell_mhist.resize(0);
    }

    spliced.reset(fracs, barcodes, max_hist, (alloc_barcodes && calc_sau), has_mt, has_exc);
    ambiguous.reset(fracs, barcodes, max_hist, (alloc_barcodes && calc_sau), has_mt, has_exc);
    unspliced.reset(fracs, barcodes, max_hist, (alloc_barcodes && calc_sau), has_mt, has_exc);

    total.reset(fracs, barcodes, max_hist, alloc_barcodes, has_mt, has_exc);
    this->calc_sau = calc_sau;
    this->build_mats = build_mats;
    this->has_visium = has_visium;

    if(!aggregate_only){
        gcounts.assign(fracs.size(), GeneCountMatrix());
    }else{
        gcounts.resize(0);
    }
    if(build_mats){
        total_mat.assign(fracs.size(), SparseMatrix());
        if(calc_sau){
            spliced_mat.assign(fracs.size(), SparseMatrix());
            ambiguous_mat.assign(fracs.size(), SparseMatrix());
            unspliced_mat.assign(fracs.size(), SparseMatrix());
        }
    }else{
        total_mat.resize(0);
        spliced_mat.resize(0);
        ambiguous_mat.resize(0);
        unspliced_mat.resize(0);
    }
}

void DownsampleResults::merge(const DownsampleResultsLocal & rhs){
    assert(fracs.size() == rhs.fracs.size());

    sum_vectors(spliced.mhist, rhs.spliced_mhist, false);
    sum_vectors(spliced.reads, rhs.spliced_reads, false);
    sum_vectors(spliced.molecules, rhs.spliced_molecules, false);

    sum_vectors(unspliced.mhist, rhs.unspliced_mhist, false);
    sum_vectors(unspliced.reads, rhs.unspliced_reads, false);
    sum_vectors(unspliced.molecules, rhs.unspliced_molecules, false);

    sum_vectors(ambiguous.mhist, rhs.ambiguous_mhist, false);
    sum_vectors(ambiguous.reads, rhs.ambiguous_reads, false);
    sum_vectors(ambiguous.molecules, rhs.ambiguous_molecules, false);

    sum_vectors(total.mhist, rhs.total_mhist, false);
    sum_vectors(total.reads, rhs.total_reads, false);
    sum_vectors(total.molecules, rhs.total_molecules, false);

    sum_vectors(reads_discarded, rhs.reads_discarded, false);
    sum_vectors(reads_excluded, rhs.reads_excluded, false);
    sum_vectors(mols_discarded, rhs.mols_discarded, false);

    if(build_mats){
        if(calc_sau){
            //matrices have been built
            assert(spliced_mat.size() == rhs.spliced_mat.size());
            for(size_t i = 0; i < spliced_mat.size(); i++){
                spliced_mat[i].merge(rhs.spliced_mat[i]);
            }
            assert(ambiguous_mat.size() == rhs.ambiguous_mat.size());
            for(size_t i = 0; i < ambiguous_mat.size(); i++){
                ambiguous_mat[i].merge(rhs.ambiguous_mat[i]);
            }
            assert(unspliced_mat.size() == rhs.unspliced_mat.size());
            for(size_t i = 0; i < unspliced_mat.size(); i++){
                unspliced_mat[i].merge(rhs.unspliced_mat[i]);
            }
        }
        assert(total_mat.size() == rhs.total_mat.size());
        for(size_t i = 0; i < total_mat.size(); i++){
            total_mat[i].merge(rhs.total_mat[i]);
        }
    }
    if(!aggregate_only){
        assert(gcounts.size() == rhs.gcounts.size());
        for(size_t i = 0; i < gcounts.size(); i++){
            gcounts[i].merge(rhs.gcounts[i]);
        }
    }
}

