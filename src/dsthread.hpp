// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#pragma once
#include "dsdata.hpp"
#include <thread>
#include "umis.hpp"

namespace scdepth{

class Downsampler;

class DownsamplerThread{
    public:
        DownsamplerThread(Downsampler & data, uint32_t umi_len, 
                size_t genes, size_t chunk_start, size_t chunk_end);

        ~DownsamplerThread();
        DownsamplerThread & operator=(const DownsamplerThread &) = delete;
        DownsamplerThread(const DownsamplerThread &) = delete;
        DownsamplerThread & operator=(const DownsamplerThread &&) = delete;
        DownsamplerThread(const DownsamplerThread &&) = delete;

        void start();
        void join();
        void operator()();

        bool                   error = false;
        DownsampleResultsLocal out;
        size_t                 processed = 0;

    private:
        void                    merge_barcode_(size_t i, const BarcodeCount & bc, const UMIResults & bco);
        std::thread             thread_;
        BarcodeData             b_;
        Downsampler           & data_;
        std::vector<uint32_t>   bc_mhist_;
        size_t                  cstart_ = 0;
        size_t                  cend_ = 0;
        uint32_t                umi_len_ = 0;
};


}
