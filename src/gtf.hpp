// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#pragma once

#include "gtl/phmap.hpp"
#include "htslib/sam.h"
#include "iitree.hpp"
#include <cstdint>
#include <set>
#include <string>
#include <vector>

using gtl::flat_hash_map;
using gtl::flat_hash_set;

namespace scdepth {

struct Exon {
    Exon() {}

    Exon(uint32_t lft, uint32_t rgt) : lft(lft), rgt(rgt) {}

    bool operator<(const Exon& rhs) const {
        return std::tie(lft, rgt) < std::tie(rhs.lft, rhs.rgt);
    }

    bool operator==(const Exon& rhs) const {
        return std::tie(lft, rgt) == std::tie(rhs.lft, rhs.rgt);
    }

    bool overlaps(const Exon& rhs) const {
        return lft <= rhs.rgt && rhs.lft <= rgt;
    }

    bool overlaps(uint32_t rlft, uint32_t rrgt) const {
        return lft <= rrgt && rlft <= rgt;
    }

    struct HashFunction {
        size_t operator()(const Exon& ex) const {
            size_t lfth = std::hash<int>()(ex.lft);
            size_t rgth = std::hash<int>()(ex.rgt) << 1;
            return lfth ^ rgth;
        }
    };

    uint32_t lft;
    uint32_t rgt;
};

struct Gene {

    using Juncs = flat_hash_set<Exon, Exon::HashFunction>;
    Gene() {}

    Gene(std::string gene_id, int tid, uint32_t lft, uint32_t rgt, char strand)
        : gene_id(gene_id), tid(tid), lft(lft), rgt(rgt), strand(strand) {}
    bool operator<(const Gene& rhs) const {
        return std::tie(tid, lft, rgt, gene_id) <
               std::tie(rhs.tid, rhs.lft, rhs.rgt, rhs.gene_id);
    }

    std::string           gene_id;
    std::string           gene_name;
    std::vector<Exon>     exons;
    Juncs                 junctions;
    std::set<std::string> transcripts;
    int                   tid;
    uint32_t              gidx;
    uint32_t              lft;
    uint32_t              rgt;
    char                  strand;
};

struct GeneIDEntry {
    std::string gene_id;
    std::string gene_name;
    uint32_t    gidx = 0;
};

using Genes      = std::vector<Gene>;
using GeneIDMap  = std::vector<GeneIDEntry>;
using GeneTree   = IITree<long int, uint32_t>;
using RefTrees   = std::vector<GeneTree>;
using GeneBreaks = std::vector<std::vector<uint32_t>>;

Genes                    gtf2genes(std::string gtf, sam_hdr_t* hdr);
RefTrees                 genes2trees(const Genes& genes, unsigned int n_chroms);
std::vector<std::string> genes2ids(const Genes& genes);
GeneBreaks               gene_breaks(const Genes& genes, unsigned int n_chroms,
                                     uint32_t delta = 5000);
bool      gtf2mapping(const Genes& genes, const std::string& out);
bool      gtf2mapping(const std::string& gtf, const std::string& out);
GeneIDMap read_gtf2mapping(const std::string& gtf_map);
void      set_gene_filter(const GeneIDMap& genes, const std::string& mod_file,
                          std::vector<bool>& out, bool warn = true);

} // namespace scdepth
