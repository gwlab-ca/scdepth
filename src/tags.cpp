// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#include "tags.hpp"
#include <charconv>
#include <iostream>
#include <sstream>
#include <string_view>
#include <zlib.h>

using namespace scdepth;

const char* scdepth::ADNA4::alphabet_str_ = "ACGT";
const char  scdepth::ADNA4::ltable_[256]  = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, 0,  -1, 1,  -1, -1, -1, 2,  -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, 3,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, 0,  -1, 1,  -1, -1, -1, 2,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, 3,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1};
const char scdepth::ADNA4::rtable_[4] = {3, 2, 1, 0};
long int   scdepth::read_barcode_debug(BGZF* fin, BarcodeData& data,
                                       const std::vector<BarcodeCount>& bcs) {
    std::cout << "Read Next offset = " << bgzf_utell(fin)
              << " number of barcodes = " << bcs.size() << "\n";
    std::cout << "  Counts Check = " << data.counts.size()
              << " Expected = " << data.bc.n_count_tags
              << " Capacity = " << data.counts.capacity();
    if(!data.counts.empty())
        std::cout << " ptr = " << &data.counts[0];
    std::cout << "\n";
    data.counts.clear();
    std::cout << "  bc check = " << &data.bc << " size = " << sizeof(data.bc)
              << "\n";
    size_t bytes_read = 0;
    auto   ret        = bgzf_read(fin, reinterpret_cast<void*>(&data.bc),
                                    sizeof(BarcodeHeader));
    std::cout << "  Read barcode ret = " << ret
              << " barcode = " << data.bc.barcode
              << " size = " << data.bc.block_size
              << " counts = " << data.bc.n_count_tags << "\n";
    std::cout << "  Expected offset = " << bcs[data.bc.barcode].offset << "\n";
    if(ret < 0)
        return -1;
    else if(ret == 0)
        return 0;
    else if(ret != sizeof(BarcodeHeader))
        return -1;
    std::cout << "Calling resize\n";
    std::cout << "  Counts Check = " << data.counts.size()
              << " Expected = " << data.bc.n_count_tags
              << " Capacity = " << data.counts.capacity();
    data.counts.resize(data.bc.n_count_tags);

    if(data.bc.barcode >= bcs.size()) {
        std::cout << "[error] bgzf_read error barcode is larger than the "
                       "barcode index\n";
        return -1;
    }
    auto& d = bcs[data.bc.barcode];
    std::cout << "  index barcode = " << d.barcode << " offset = " << d.offset
              << " bytes = " << d.total_data_bytes << " reads = " << d.countable
              << "\n";
    std::cout << "    bgzf_utell = " << bgzf_utell(fin) << "\n";
    size_t midx      = 0;
    size_t umi_count = 0;
    for(size_t i = 0; i < data.bc.n_count_tags; i++) {
        std::cout << "    Reading i = " << i << " midx = " << midx
                  << " counts size = " << data.counts.size() << "\n";
        auto ret = bgzf_read(fin, reinterpret_cast<void*>(&data.counts[midx]),
                               sizeof(UMITag));
        std::cout << "    Read UMI gene = " << data.counts[midx].gene
                  << " umi = " << data.counts[midx].umi << " bytes = " << ret
                  << "\n";
        // std::cout << "    Read UMI gene = " << data.counts[midx].gene
        //     << " umi = " << data.counts[midx].umi << " nfrags = " <<
        //     data.counts[midx].n_frags << "\n";
        if(ret <= 0) {
            std::cout
                << "[error] bgzf_read error or EOF while reading CountTag\n";
            return -1;
        } else if(ret != sizeof(UMITag)) {
            std::cout << "[error] truncated CountTag\n";
            return -1;
        }
        bytes_read += ret;
        umi_count++;
        midx++;
    }
    std::cout << "  bytes read = " << bytes_read
              << " expected = " << d.total_data_bytes
              << " umi count = " << umi_count << " raw = " << d.raw_molecules
              << "\n";
    if(umi_count != d.raw_molecules) {
        std::cout << "[error] barcode = " << data.bc.barcode
                  << " expected molecules = " << d.raw_molecules
                  << " read molecules = " << data.counts.size() << "\n";
    }
    if(bytes_read != data.bc.block_size) {
        std::cout << "[error] barcode = " << data.bc.barcode
                  << " Bytes read = " << bytes_read
                  << " block size = " << data.bc.block_size << "\n";
        return -1;
    }
    std::cout << "  Counts Confirm = " << data.counts.size()
              << " Expected = " << data.bc.n_count_tags
              << " Capacity = " << data.counts.capacity();
    if(!data.counts.empty())
        std::cout << " ptr = " << &data.counts[0];
    std::cout << "\n";
    std::cout << "Done\n";
    return bytes_read + sizeof(BarcodeHeader);
}

long int scdepth::read_barcode(BGZF* fin, BarcodeData& data) {
    data.counts.clear();
    size_t bytes_read = 0;
    auto   ret        = bgzf_read(fin, reinterpret_cast<void*>(&data.bc),
                                  sizeof(BarcodeHeader));
    if(ret < 0)
        return -1;
    else if(ret == 0)
        return 0;
    else if(ret != sizeof(BarcodeHeader))
        return -1;
    data.counts.resize(data.bc.n_count_tags);
    // std::cout << "Read barcode " << data.bc.barcode << " size = " <<
    // data.bc.block_size << " counts = " << data.bc.n_count_tags << "\n";
    size_t midx = 0;
    for(size_t i = 0; i < data.bc.n_count_tags; i++) {
        auto ret = bgzf_read(fin, reinterpret_cast<void*>(&data.counts[midx]),
                             sizeof(UMITag));
        // std::cout << "  Read count gene = " << data.counts[midx].gene
        //     << " umi = " << data.counts[midx].umi << " nfrags = " <<
        //     data.counts[midx].n_frags << "\n";
        if(ret <= 0) {
            std::cerr
                << "[error] bgzf_read error or EOF while reading CountTag\n";
            return -1;
        } else if(ret != sizeof(UMITag)) {
            std::cerr << "[error] truncated CountTag\n";
            return -1;
        }
        bytes_read += ret;
        midx++;
    }
    if(bytes_read != data.bc.block_size) {
        std::cerr << "[error] barcode = " << data.bc.barcode
                  << " Bytes read = " << bytes_read
                  << " block size = " << data.bc.block_size << "\n";
        return -1;
    }
    return bytes_read;
}

// whitespace/tab-delimited tokenizer
static inline bool next_tok(std::string_view& sv, std::string_view& tok) {
    // skip leading spaces/tabs
    size_t i = 0;
    while(i < sv.size() && (sv[i] == '\t' || sv[i] == ' '))
        ++i;
    sv.remove_prefix(i);
    if(sv.empty())
        return false;

    size_t j = 0;
    while(j < sv.size() && sv[j] != '\t' && sv[j] != ' ')
        ++j;
    tok = sv.substr(0, j);
    sv.remove_prefix(j);
    return true;
}

template <class UInt>
static inline bool parse_uint(std::string_view sv, UInt& out) {
    const char* b = sv.data();
    const char* e = b + sv.size();
    auto        r = std::from_chars(b, e, out, 10);
    return r.ec == std::errc{} && r.ptr == e;
}

inline bool parse_barcode_line(const std::string& line, BarcodeCount& bc) {
    std::string_view sv(line), tok;

    // barcode
    if(!next_tok(sv, tok))
        return false;
    bc.barcode.assign(tok.data(), tok.size());

    if(!next_tok(sv, tok) || !parse_uint<uint32_t>(tok, bc.index))
        return false;
    if(!next_tok(sv, tok) || !parse_uint<uint32_t>(tok, bc.total))
        return false;
    if(!next_tok(sv, tok) || !parse_uint<uint32_t>(tok, bc.countable))
        return false;
    if(!next_tok(sv, tok) || !parse_uint<uint32_t>(tok, bc.raw_molecules))
        return false;
    if(!next_tok(sv, tok) || !parse_uint<uint64_t>(tok, bc.offset))
        return false;
    if(!next_tok(sv, tok) || !parse_uint<uint32_t>(tok, bc.total_data_bytes))
        return false;
    bc.seed       = 0;
    bc.random_hex = 0;
    bc.poly_a     = 0;
    bc.sample.clear();

    // Optional only in newer files
    if(next_tok(sv, tok)) {
        bc.sample.assign(tok.data(), tok.size());
        if(!next_tok(sv, tok) || !parse_uint<uint32_t>(tok, bc.random_hex))
            return false;
        if(!next_tok(sv, tok) || !parse_uint<uint32_t>(tok, bc.poly_a))
            return false;
    }
    return true;
}

std::vector<BarcodeCount> scdepth::read_barcode_index(const std::string& fin) {
    std::vector<BarcodeCount> ret;
    gzFile                    gz = gzopen(fin.c_str(), "rb");
    if(!gz)
        return ret;

    gzbuffer(gz, 4 << 20); // 1MB buffer (try 4MB too)

    std::string line;
    line.resize(4096); // starting capacity
    char* buf = line.data();

    // read header line
    if(!gzgets(gz, buf, (int)line.size())) {
        gzclose(gz);
        return ret;
    }

    while(true) {
        char* p = gzgets(gz, buf, (int)line.size());
        if(!p)
            break;

        while(strchr(buf, '\n') == nullptr && !gzeof(gz)) {
            size_t old = line.size();
            line.resize(old * 2);
            buf = line.data();
            break;
        }

        // strip newline
        size_t n = strlen(buf);
        if(n && (buf[n - 1] == '\n' || buf[n - 1] == '\r'))
            buf[--n] = '\0';
        if(n == 0)
            continue;

        BarcodeCount bc;
        if(!parse_barcode_line(buf, bc)) {
            ret.clear();
            gzclose(gz);
            return ret;
        }
        ret.push_back(std::move(bc));
    }

    gzclose(gz);
    std::sort(ret.begin(), ret.end());
    return ret;
}
