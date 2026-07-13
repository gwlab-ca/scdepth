// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#include "downsampler.hpp"
#include "pcg_random.hpp"
#include <numeric>
#include <random>
#include <iostream>
#include "gtf.hpp"
#include "dsthread.hpp"
#include "visbins.hpp"

using namespace scdepth;

bool Downsampler::init(const std::string & prefix, 
        const std::string & mt_prefix, const std::string & mt_file, 
        const std::string & mod_file, const std::string & exclude_file, 
        size_t max_hist, bool build_matrices, bool calc_sau){

    this->prefix = prefix;
    this->build_mats_ = build_matrices;
    this->calc_sau_ = calc_sau;
    max_hist_ = max_hist;
    barcodes = read_barcode_index(prefix + "_barcode_index.txt.gz");
    if(barcodes.empty()){
        std::cerr << "[error] barcode index is empty";
        return false;
    }
    if(!barcode2sample.empty() && barcode2sample.size() != barcodes.size()){
        std::cerr << "barcode2sample size mismatch " << barcode2sample.size() << " vs expected = " << barcodes.size() << "\n";
        return false;
    }

    auto genes = read_gtf2mapping(prefix + "_genes.txt.gz");
    if(genes.empty()){
        std::cerr << "[error] gene mapping file is empty";
        return false;
    }

    mt_filter.resize(genes.size());
    if(!mt_prefix.empty() && !mt_file.empty()){
        std::cerr << "[error] cannot provide MT prefix and MT file";
        return false;
    }else if(!mt_prefix.empty()){
        for(size_t i = 0; i < genes.size(); i++){
            if(genes[i].gene_name.compare(0, mt_prefix.size(), mt_prefix) == 0){
                mt_filter[genes[i].gidx] = true;
            }
        }
        this->has_mt = true;
    }else if(!mt_file.empty()){
        this->has_mt = true;
        set_gene_filter(genes, mt_file, mt_filter);
    }

    exclude_filter.resize(genes.size());
    if(!exclude_file.empty()){
        has_exc_ = true;
        set_gene_filter(genes, exclude_file, mod_filter, false);
    }

    mod_filter.resize(genes.size());
    if(!mod_file.empty()){
        this->has_mod = true;
        set_gene_filter(genes, mod_file, mod_filter);
    }

    genes_ = genes.size();
    countable_ = std::accumulate(barcodes.begin(), barcodes.end(), 0, 
        [](uint32_t sum, const BarcodeCount & bc){ return sum + bc.countable;
    });

    raw_molecules_ = std::accumulate(barcodes.begin(), barcodes.end(), 0, 
        [](uint32_t sum, const BarcodeCount & bc){ return sum + bc.raw_molecules; 
    });

    return true;
}

bool Downsampler::init_visium(const std::vector<uint32_t> & rows, const std::vector<uint32_t> & cols, 
        std::vector<uint32_t> & in_tissue, std::vector<uint32_t> & countable_reads, 
        std::vector<uint32_t> & total_reads, size_t total_rows, size_t total_cols){

    if(barcodes.empty()){
        std::cerr << "Barcodes are not loaded `init` must be run first";
        return false;
    }

    if(barcodes.size()!= rows.size() 
        || barcodes.size() != cols.size() 
        || in_tissue.size() != barcodes.size()
        || countable_reads.size() != barcodes.size()
        || total_reads.size() != barcodes.size())
    {

        std::cerr << "barcodes/rows/cols all must be the same size for rows/cols and empty or the same size for countable_reads/total_reads/in_tissue"
            << " barcodes = " << barcodes.size() 
            << " rows = " << rows.size() 
            << " cols = " << cols.size() 
            << " in_tissue = " << in_tissue.size() 
            << " countable_reads = " << countable_reads.size() 
            << " total_reads = " << total_reads.size() 
            << "\n";

        return false;
    }

    has_visium = true;
    bin_pos.resize(rows.size());
    //uint32_t tot_rows = 0, tot_cols = 0;
    for(size_t i = 0; i < rows.size(); i++){
        bin_pos[i].b = i;
        bin_pos[i].r = rows[i];
        bin_pos[i].c = cols[i];
        bin_pos[i].in_tissue = (in_tissue[i] > 0);
        bin_pos[i].countable = countable_reads[i];
        bin_pos[i].total = total_reads[i];
    }

    this->total_rows = total_rows;
    this->total_cols = total_cols;

    return true;
}

bool Downsampler::reset_visium(){
    total_rows = 0;
    total_cols = 0;
    bin_pos.clear();
    bin_pos.resize(0);
    has_visium = false;
    return true;
}

uint32_t sum_sparse(size_t row, const SparseMatrix & s){
    uint32_t ret = 0;
    for(size_t i = s.indptr[row]; i < s.indptr[row + 1]; i++){
        ret += s.data[i];
    }
    return ret;
}

bool Downsampler::downsample(std::vector<double> & fracs, 
        uint32_t umi_len, uint64_t seed, unsigned int threads, bool aggregate_only,
        const std::string & umi_mode, bool correct_multi_umis,
        const std::string & primer_mode, const std::vector<uint32_t> & barcode2sample){
        //bool profile_umi_singletons){
    chunks_.clear();
    if(fracs.empty()){
        std::cerr << "[error] fractions list is empty\n";
        return false;
    }
    for(auto & f : fracs){
        if(f > 1.0){
            std::cerr << "[error] fractions must be >0 and <= 1.0 frac = " << f << "\n";
            return false;
        }
    }
    this->barcode2sample = barcode2sample;
    samples = 0;
    for(auto s : barcode2sample){
        samples = std::max(s + 1, samples);
    }
    //std::cout << "total samples = " << samples << std::endl;
    this->aggregate_only = aggregate_only;
    if(threads < 1) threads = 1;
    if(threads > barcodes.size()) threads = barcodes.size();

    //umi_singletons = profile_umi_singletons;
    umi_multi = correct_multi_umis;
    if(umi_mode == "directed"){
        umi_directed = true;
    }else if(umi_mode == "none"){
        umi_directed = false;
    }else{
        std::cerr << "[error] invalid umi resolution mode\n";
        return false;
    }

    if(primer_mode == "merge"){
        this->primer_mode = PrimerMode::Merge;
    }else if(primer_mode == "polyA"){
        this->primer_mode = PrimerMode::PolyAOnly;
    }else if(primer_mode == "random_hex"){
        this->primer_mode = PrimerMode::RandomHexOnly;
    }else{
        std::cerr << "[error] invalid primer mode\n";
        return false;
    }

    clear_output();
    output.has_mt = has_mt;
    output.has_mod = has_mod;
    output.has_exc= has_exc_;
    output.primer_mode = primer_mode;
    has_visium = false;
    if(total_rows > 0 && total_cols > 0){
        has_visium = true;
    }
    if(!has_visium) reset_visium();
    output.reset(fracs, barcodes.size(), genes_, max_hist_, build_mats_, this->calc_sau_, !aggregate_only, has_visium, samples);
    output.aggregate_only = aggregate_only;
    //if visium output is specified get it ready
    if(total_rows > 0 && total_cols > 0){
        has_visium = true;
    }
    //std::cout << "total rows = " << total_rows << " cols = " << total_cols << " use qc = " << use_qc << "\n";
    output.has_visium = has_visium;
    umi_len_ = umi_len;
    uint64_t mols_per_thread = raw_molecules_ / threads;
    uint64_t ctotal = 0;
    chunks_.push_back(0);
    for(size_t i = 0; i < barcodes.size(); i++){
        size_t bcount = i - chunks_.back();
        if(chunks_.size() < threads && bcount > 0 && 
                (ctotal + barcodes[i].raw_molecules) > mols_per_thread){
            chunks_.push_back(i);
            ctotal = 0;
        }
        ctotal += barcodes[i].raw_molecules;
    }
    if((chunks_.size() == 1) || chunks_.back() != barcodes.size()){
        chunks_.push_back(barcodes.size());
    }
    /*
    for(size_t i = 0; i < chunks_.size() -1; i++){
        std::cout << "Chunk = " << chunks_[i] << " - " << chunks_[i + 1] << "\n";
    }
    */
    output.seed = seed;
    pcg64 rng(seed);
    std::uniform_int_distribution<uint64_t> dist;
    for(auto & bc : barcodes){
        bc.seed = dist(rng);
    }

    if((chunks_.size() - 1) < threads){
        threads = chunks_.size() - 1;
    }
    //std::cout << " threads = " << threads << " chunks = " << chunks_.size() << "\n";

    std::vector<DownsamplerThread*> workers;
    for (size_t i = 0; i < threads; i++) {
        workers.push_back(new DownsamplerThread(*this, umi_len_, genes_, chunks_[i], chunks_[i + 1]));
    }

    for(size_t i = 1; i < threads; i++){
        workers[i]->start();
    }
    (*workers[0])();
    for(size_t i = 1; i < threads; i++){
        (*workers[i]).join();
    }

    bool error = false;
    for(size_t i = 0; i < workers.size(); i++){
        auto & w = workers[i];
        //handle error
        if(w->error){
            std::cerr << "[error] with worker " << workers[i] << " chunk = " << chunks_[i] << " - " << chunks_[i + 1] << "\n";
            error = true;
        }

    }
    for(size_t i = 0; i < workers.size(); i++){
        if(!error){
            output.merge(workers[i]->out);
        }
        delete workers[i];
        workers[i] = nullptr;
    }
    workers.resize(0);

    /*
    if(build_mats_ && calc_sau_){
        for(size_t i = 0; i < fracs.size(); i++){
            for(size_t j = 0; j < barcodes.size(); j++){
                size_t bidx = i * barcodes.size() + barcodes[j].index;
                size_t expected = output.spliced.bc_mols[bidx];
                size_t mat_sum = sum_sparse(barcodes[j].index, output.spliced_mat[i]);
                size_t mgenes = output.spliced_mat[i].indptr[j + 1] - output.spliced_mat[i].indptr[j];
                if(expected != mat_sum){
                    std::cout << "sparse spliced matrix error for frac = " << i << " j = " << j
                        << " spliced mat barcode= " << barcodes[j].index 
                        << " expected = " << expected 
                        << " matrix sum = " << mat_sum
                        << " delta = " << ((long int)expected - (long int)mat_sum)
                        << " reads = " << output.spliced.bc_reads[bidx]
                        << " genes = " << output.spliced.bc_genes[bidx]
                        << " mgenes = " << mgenes
                        << " delta = " << ((long int)output.spliced.bc_genes[bidx] - (long int)mgenes)
                        << "\n";
                    return true;
                }

                expected = output.unspliced.bc_mols[bidx];
                mat_sum = sum_sparse(barcodes[j].index, output.unspliced_mat[i]);
                mgenes = output.unspliced_mat[i].indptr[j + 1] - output.unspliced_mat[i].indptr[j];
                if(expected != mat_sum){
                    std::cout << "sparse unspliced matrix error for frac = " << i << " j = " << j
                        << " unspliced mat barcode= " << barcodes[j].index 
                        << " expected = " << expected 
                        << " matrix sum = " << mat_sum
                        << " delta = " << ((long int)expected - (long int)mat_sum)
                        << " reads = " << output.unspliced.bc_reads[bidx]
                        << " genes = " << output.unspliced.bc_genes[bidx]
                        << " mgenes = " << mgenes
                        << " delta = " << ((long int)output.unspliced.bc_genes[bidx] - (long int)mgenes)
                        << "\n";
                    return true;
                }

                expected = output.ambiguous.bc_mols[bidx];
                mat_sum = sum_sparse(barcodes[j].index, output.ambiguous_mat[i]);
                mgenes = output.ambiguous_mat[i].indptr[j + 1] - output.ambiguous_mat[i].indptr[j];
                if(expected != mat_sum){
                    std::cout << "sparse ambiguous matrix error for frac = " << i << " j = " << j
                        << " ambiguous mat barcode= " << barcodes[j].index 
                        << " expected = " << expected 
                        << " matrix sum = " << mat_sum
                        << " delta = " << ((long int)expected - (long int)mat_sum)
                        << " reads = " << output.ambiguous.bc_reads[bidx]
                        << " genes = " << output.ambiguous.bc_genes[bidx]
                        << " mgenes = " << mgenes
                        << " delta = " << ((long int)output.ambiguous.bc_genes[bidx] - (long int)mgenes)
                        << "\n";
                    return true;
                }
            }
        }
    }
    */
    return !error;
}

bool Downsampler::write_gene_mats(const std::string & out, const std::vector<uint32_t> & idx,
        uint32_t bin_div){

    if(bin_div > 0 && !has_visium){
        std::cerr << "[error] Cannot use a bin divisor without init_visium\n";
        return false;
    }

    if(idx.size() >= output.barcodes){
        std::cerr << "[error] idx size exceeds the number of barcodes\n";
        return false;
    }

    if(!std::is_sorted(idx.begin(), idx.end())){
        std::cerr << "[error] idx must be sorted in asecending order\n";
        return false;
    }

    if(idx.back() >= output.barcodes){
        std::cerr << "[error] idx contains elements that exceed the number of barcodes\n";
        return false;
    }

    BGZF * bout= bgzf_open(out.c_str(), "w");
    if(bout == nullptr){
        std::cerr << "[error] Error opening output file\n";
        return false;
    }

    if(bgzf_write(bout, reinterpret_cast<const void*>(&output.seed), sizeof(uint64_t)) <= 0){
        std::cerr << "[error] Error writing seed\n";
        bgzf_close(bout);
        return false;
    };

    uint32_t N = output.genes;
    if(bgzf_write(bout, reinterpret_cast<const void*>(&N), sizeof(uint32_t)) <= 0){
        std::cerr << "[error] Error writing gene count\n";
        bgzf_close(bout);
        return false;
    };

    N = idx.size();
    if(bgzf_write(bout, reinterpret_cast<const void*>(&N), sizeof(uint32_t)) <= 0){
        std::cerr << "[error] Error writing barcode count\n";
        bgzf_close(bout);
        return false;
    };

    if(!bz_write_vector(bout, output.fracs)){
        std::cerr << "[error] Error writing fracs\n";
        bgzf_close(bout);
        return false;
    };

    GeneCountMatrix go;
    BinnedMap bm;
    if(bin_div > 0){
        bm.build(*this, bin_div);
    }
    for(auto & g : output.gcounts){
        go.clear();
        if(bin_div > 0){
            make_filt_gbin_mat(bm, g, go, idx);
        }else{
            for(auto & i : idx){
                go.indices.insert(go.indices.end(), std::next(g.indices.begin(), g.indptr[i]), std::next(g.indices.begin(), g.indptr[i + 1]));
                go.indptr.push_back(go.indices.size());
            }
        }

        if(!bz_write_vector(bout, go.indptr)){
            std::cerr << "[error] Error writing indptr\n";
            bgzf_close(bout);
            return false;

        }
        if(!bz_write_vector(bout, go.indices)){
            std::cerr << "[error] Error writing indices\n";
            bgzf_close(bout);
            return false;
        }
    }

    bgzf_close(bout);

    return true;
}

bool Downsampler::write_gene_baseline(const std::string & out, const std::vector<uint32_t> & idx,
        uint32_t step, uint32_t bin_div){

    if(!build_mats_){
        std::cerr << "[error] Build mats must be called for the baseline\n";
        return false;
    }

    if(bin_div > 0 && !has_visium){
        std::cerr << "[error] Cannot use a bin divisor without init_visium\n";
        return false;
    }

    if(step >= output.fracs.size()){
        std::cerr << "[error] Step cannot be larger than the number of fractions\n";
        return false;
    }

    if(idx.size() > output.barcodes){
        std::cerr << "[error] idx size exceeds the number of barcodes\n";
        return false;
    }

    if(!std::is_sorted(idx.begin(), idx.end())){
        std::cerr << "[error] idx must be sorted in asecending order\n";
        return false;
    }

    if(idx.back() >= output.barcodes){
        std::cerr << "[error] idx contains elements that exceed the number of barcodes\n";
        return false;
    }

    BGZF * bout= bgzf_open(out.c_str(), "w");
    if(bout == nullptr){
        std::cerr << "[error] Error opening output file\n";
        return false;
    }

    uint32_t N = output.genes;
    if(bgzf_write(bout, reinterpret_cast<const void*>(&N), sizeof(uint32_t)) <= 0){
        std::cerr << "[error] Error writing genes\n";
        bgzf_close(bout);
        return false;
    };

    N = idx.size();
    if(bgzf_write(bout, reinterpret_cast<const void*>(&N), sizeof(uint32_t)) <= 0){
        std::cerr << "[error] Error writing barcode count\n";
        bgzf_close(bout);
        return false;
    };

    auto & g = output.total_mat[step];

    SparseMatrix go;
    if(bin_div > 0){
        BinnedMap bm;
        bm.build(*this, bin_div);
        make_filt_cbin_mat(bm, g, go, idx);
    }else{
        for(auto & i : idx){
            go.indices.insert(go.indices.end(), std::next(g.indices.begin(), g.indptr[i]), 
                    std::next(g.indices.begin(), g.indptr[i + 1]));
            go.data.insert(go.data.end(), std::next(g.data.begin(), g.indptr[i]), 
                    std::next(g.data.begin(), g.indptr[i + 1]));
            go.indptr.push_back(go.indices.size());
        }
    }


    if(!bz_write_vector(bout, go.indptr)){
        std::cerr << "[error] Error writing indptr\n";
        bgzf_close(bout);
        return false;

    }
    if(!bz_write_vector(bout, go.indices)){
        std::cerr << "[error] Error writing indices\n";
        bgzf_close(bout);
        return false;
    }

    if(!bz_write_vector(bout, go.data)){
        std::cerr << "[error] Error writing data\n";
        bgzf_close(bout);
        return false;
    }

    bgzf_close(bout);
    return true;
}
