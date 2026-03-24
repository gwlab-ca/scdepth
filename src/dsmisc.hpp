
// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#pragma once


#include <cstdint>
#include <cassert>
#include <atomic>
#include <vector>
#include "htslib/bgzf.h"

namespace scdepth{

struct SpinLock {
    std::atomic_flag f = ATOMIC_FLAG_INIT;
    void lock() noexcept {
        while (f.test_and_set(std::memory_order_acquire)) { /* spin */ }
    }
    void unlock() noexcept {
        f.clear(std::memory_order_release);
    }

    void reset() noexcept {
        f.clear(std::memory_order_relaxed);
    }
};

struct BarcodePos{
    uint32_t b = 0;
    uint32_t r = 0;
    uint32_t c = 0;
    uint32_t in_tissue = 0;
    uint32_t countable = 0;
    uint32_t total = 0;

    bool operator<(const BarcodePos & rhs){
        return (r == rhs.r) ? (c < rhs.c) : (r < rhs.r);
    }
};

template <typename T>
inline bool bz_write_vector(BGZF * bout, const std::vector<T> & v){
    uint64_t sz = v.size();
    if(bgzf_write(bout, reinterpret_cast<const void*>(&sz), sizeof(sz)) <= 0){
        return false;
    };

    if(bgzf_write(bout, reinterpret_cast<const void*>(&v[0]), sizeof(T) * sz) <= 0){
        return false;
    };
    return true;
}

template <typename T>
inline bool bz_read_vector(BGZF * bin, std::vector<T> & v){
    uint64_t sz = 0;

    if(bgzf_read(bin, reinterpret_cast<void*>(&sz), sizeof(sz)) <= 0){
        return false;
    }

    v.resize(sz);

    if(sz > 0){
        if(bgzf_read(bin, reinterpret_cast<void*>(&v[0]), sizeof(T) * sz) <= 0){
            return false;
        }
    }

    return true;
}

};
