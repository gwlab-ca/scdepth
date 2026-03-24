// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#include "downsampler.hpp"
#include <string>
#include <iostream>
#include <string>
#include <vector>
#include "goverlaps.hpp"
#include <numeric>

using namespace scdepth;


int main(int argc, char * argv[]){
    Downsampler ds;

    std::string prefix = argv[1];
    unsigned int threads = std::stoi(argv[2]);
    unsigned int umi_len = std::stoi(argv[3]);

    ds.init(prefix, "MT-", "", "", "", 50, true, false);

    //std::vector<double> fracs = {1.0, 0.75, 0.5, 0.25, 0.1};
    //std::vector<double> fracs = {0.1, 0.25, 0.5, 0.75, 1.0};
    std::vector<double> fracs = {1.0};
    ds.downsample(fracs, umi_len, 42, threads, true);

    /*
    std::cout << "Gene counts check\n";
    for(auto & g : ds.output.gcounts){
        std::cout << "  " << g.indptr.size() << " / " << g.indices.size() << "\n";
    }

    std::vector<uint32_t> idx(ds.barcodes.size(), 0);
    std::iota(idx.begin(), idx.end(), 0);
    std::cout << "write baseline\n";
    std::cout << "Check " << idx.back() << " sz = " << ds.barcodes.size() << "\n";
    ds.write_gene_baseline("test.gz", idx, 0);
    //auto bc1 = parse_bc_db(argv[3]);
    //auto bc2 = parse_bc_db(argv[4]);
    //uint32_t umi_len = std::stoi(argv[3]);
    //

    std::string prefix = argv[1];

    std::string full = prefix + "_full.gz";
    std::vector<std::string> seeds;
    for(size_t i = 0; i < 10; i++){
        seeds.push_back(prefix + "_seed_" + std::to_string(i) + ".gz");
    }

    GeneOverlap go;
    if(!go.init(full, seeds)){
        std::cout << "init error\n";
    }
    if(!go.calculate()){
        std::cout << "calculate error\n";
    }
    */

    return 0;
}
