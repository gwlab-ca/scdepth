// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#include "gtf.hpp"
#include "gzstream.hpp"
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>

using namespace scdepth;

template <typename O>
Genes gtf2genes_parse(const std::string& gtf, sam_hdr_t* hdr) {
    std::regex                 words_regex{R"res((\S+)\s+"?([^;"]+)"?;)res"};
    std::array<std::string, 9> fields;
    std::string                line, gene_id, tx_id, gene_name;
    flat_hash_map<std::string, int>               tid_map;
    flat_hash_map<std::string, Gene>              gene_map;
    flat_hash_map<std::string, std::vector<Exon>> tx_map;
    Genes                                         genes;

    O in(gtf.c_str());
    if(!in) {
        std::cerr << "[error] Couldn't open " << gtf << "\n";
        return genes;
    }
    size_t line_no = 0;
    while(getline(in, line)) {
        line_no++;
        if(line.empty() || line[0] == '#') {
            continue;
        }
        bool               parse_error = false;
        std::istringstream ss(line);
        for(size_t i = 0; i < fields.size(); i++) {
            if(!std::getline(ss, fields[i], '\t')) {
                parse_error = true;
                break;
            }
        }

        if(parse_error) {
            std::cerr << "[error] Error malformed line on line number "
                      << line_no << " [" << line << "]\n";
            genes.clear();
            return genes;
        }

        gene_id.clear();
        gene_name.clear();
        tx_id.clear();

        if(fields[2] != "gene" && fields[2] != "exon")
            continue;
        int tid = -1;
        if(hdr != nullptr) {
            tid = bam_name2id(hdr, fields[0].c_str());
        } else {
            // make arbitrary tid map based on fields
            auto res = tid_map.insert({fields[0], tid_map.size()});
            tid      = res.first->second;
        }

        if(tid < 0)
            continue;
        for(auto it = std::sregex_iterator(fields[8].begin(), fields[8].end(),
                                           words_regex);
            it != std::sregex_iterator(); it++) {
            const std::smatch& m = *it;
            if(m[1].str() == "gene_id") {
                gene_id = m[2].str();
            } else if(m[1].str() == "gene_name") {
                gene_name = m[2].str();
            } else if(m[1].str() == "transcript_id") {
                tx_id = m[2].str();
            }
        }

        if(gene_id.empty()) {
            std::cerr << "[error] Error missing gene_id on line number "
                      << line_no << " [" << line << "]\n";
            genes.clear();
            return genes;
        }

        if(gene_name.empty()) {
            gene_name = gene_id;
        }

        auto  ret = gene_map.insert(std::make_pair(gene_id, Gene()));
        auto& g   = ret.first->second;
        if(fields[2] == "gene") {
            g.gene_id   = gene_id;
            g.gene_name = gene_name;
            g.tid       = tid;
            g.lft       = std::stoi(fields[3]) - 1;
            g.rgt       = std::stoi(fields[4]) - 1;
            g.strand    = fields[6][0];
        } else if(fields[2] == "exon" && !tx_id.empty()) {
            auto tret =
                tx_map.insert(std::make_pair(tx_id, std::vector<Exon>()));
            tret.first->second.push_back(
                Exon(std::stoi(fields[3]) - 1, std::stoi(fields[4]) - 1));
            g.transcripts.insert(tx_id);
        }
    }

    for(auto& p : gene_map) {
        if(p.second.gene_id.empty()) {
            std::cerr << "[error] Missing gene record for " << p.first << "\n";
            genes.clear();
            return genes;
        }
        auto& g = p.second;

        for(auto& tid : g.transcripts) {
            auto res = tx_map[tid];
            std::sort(res.begin(), res.end());
            g.exons.push_back(res.front());
            for(size_t i = 1; i < res.size(); i++) {
                g.junctions.insert(Exon(res[i - 1].rgt, res[i].lft));
                g.exons.push_back(res[i]);
            }
        }
        // size_t N = g.exons.size();
        if(g.exons.size() > 1) {
            std::sort(g.exons.begin(), g.exons.end());
            auto res = g.exons.begin();
            auto it  = res;
            while(++it != g.exons.end()) {
                // IF they don't overlap add after res
                if(!res->overlaps(*it) && ++res != it) {
                    *res = std::move(*it);
                } else {
                    // res->lft = std::min(res->lft, it->lft);
                    res->rgt = std::max(res->rgt, it->rgt);
                }
            }
            g.exons.erase(++res, g.exons.end());
        }

        genes.push_back(g);
        g.transcripts.clear();
    }
    gene_map.clear();
    std::sort(genes.begin(), genes.end());
    /*
    std::unordered_map<std::string, size_t> gid_sort;
    {

        std::vector<std::string> gids;
        for(auto & g : genes){
            gids.push_back(g.gene_id);
        }
        std::sort(gids.begin(), gids.end());
        for(size_t i = 0; i < gids.size(); i++){
            gid_sort[gids[i]] = i;
        }
    }
    */
    for(size_t i = 0; i < genes.size(); i++) {
        genes[i].gidx = i; // gid_sort[genes[i].gene_id];
        /*
        std::cout << "  tid = " << genes[i].tid << " pos = " << genes[i].lft <<
        " - " << genes[i].rgt
            << " id = " << genes[i].gene_id << " gidx = " << genes[i].gidx
            << " exons = " << genes[i].exons.size() << " junctions = " <<
        genes[i].junctions.size()
            << "\n";
        */
    }
    return genes;
}

scdepth::Genes scdepth::gtf2genes(std::string gtf, sam_hdr_t* hdr) {
    std::string gzp = ".gz";
    if(gtf.size() > 3 && std::equal(gzp.rbegin(), gzp.rend(), gtf.rbegin())) {
        return gtf2genes_parse<scdepth::gzifstream>(gtf, hdr);
    } else {
        return gtf2genes_parse<std::ifstream>(gtf, hdr);
    }
}

scdepth::RefTrees scdepth::genes2trees(const Genes& genes,
                                       unsigned int n_chroms) {
    RefTrees trees;
    trees.resize(n_chroms);
    for(size_t i = 0; i < genes.size(); i++) {
        auto& g = genes[i];
        trees[g.tid].add(g.lft, g.rgt, i);
    }
    for(size_t i = 0; i < trees.size(); i++) {
        trees[i].index();
        // std::cout << "  tid = " << i << " index size = " << trees[i].size()
        // << "\n";
    }
    return trees;
}

std::vector<std::string> scdepth::genes2ids(const Genes& genes) {
    std::vector<std::string> gene_ids;
    gene_ids.resize(genes.size());
    for(auto& g : genes) {
        gene_ids[g.gidx] = g.gene_id;
    }
    return gene_ids;
}

scdepth::GeneBreaks scdepth::gene_breaks(const Genes& genes,
                                         unsigned int n_chroms,
                                         uint32_t     delta) {
    GeneBreaks ret;
    ret.resize(n_chroms);
    auto it = genes.begin();
    // std::cout << "Number of genes = " << genes.size() << "\n";
    while(it != genes.end()) {
        std::vector<uint32_t> breaks;
        auto                  start = it;
        auto                  crgt  = start->rgt;
        // std::cout << "Processing " << start->tid << "\n";
        while(it != genes.end() && it->tid == start->tid) {
            if((crgt + delta) < it->lft) {
                // std::cout << "    breakpoint " << crgt << "\n";
                breaks.push_back(crgt + delta / 2);
                crgt = it->rgt;
            } else {
                crgt = std::max(it->rgt, crgt);
            }
            // std::cout << "  " << it->gene_id << " " << it->lft << " - " <<
            // it->rgt << "\n";;
            it++;
        }
        ret[start->tid] = breaks;
    }
    // std::cout << "Done gene breaks\n";
    return ret;
}

template <typename O>
bool gene_mapping_write(const Genes& genes, const std::string& out) {
    O of(out.c_str());
    if(!of) {
        std::cerr << "[error] Error opening " << out << " for writing\n";
        return false;
    }
    of << "gene_index\tgene_id\tgene_name\n";
    for(size_t i = 0; i < genes.size(); i++) {
        of << i << "\t" << genes[i].gene_id << "\t" << genes[i].gene_name
           << "\n";
    }
    return true;
}

bool scdepth::gtf2mapping(const Genes& genes, const std::string& out) {
    std::string gzp = ".gz";
    if(out.size() > 3 && std::equal(gzp.rbegin(), gzp.rend(), out.rbegin())) {
        return gene_mapping_write<scdepth::gzofstream>(genes, out);
    } else {
        return gene_mapping_write<std::ofstream>(genes, out);
    }
}

bool scdepth::gtf2mapping(const std::string& gtf, const std::string& out) {
    Genes genes = gtf2genes(gtf, nullptr);
    return gtf2mapping(genes, out);
}

GeneIDMap scdepth::read_gtf2mapping(const std::string& gtf_map) {
    std::string line;
    GeneIDMap   ret;
    gzifstream  in(gtf_map.c_str());
    if(!in) {
        std::cerr << "[error] could not open the gene id map " << gtf_map
                  << "\n";
        return ret;
    }

    if(!std::getline(in, line)) {
        std::cerr << "[error] file is empty: " << gtf_map << "\n";
        return ret;
    }

    while(std::getline(in, line)) {
        if(line.empty()) {
            continue;
        }
        std::istringstream iss(line);
        GeneIDEntry        gene;
        if(!(iss >> gene.gidx >> gene.gene_id >> gene.gene_name)) {
            std::cerr << "[error] failed to parse gtf mapping line from "
                      << gtf_map << " line :\n"
                      << line << "\n";
            ret.clear();
            return ret;
        }
        ret.push_back(gene);
    }
    return ret;
}

void scdepth::set_gene_filter(const GeneIDMap&   genes,
                              const std::string& mod_file,
                              std::vector<bool>& out, bool warn) {
    out.resize(genes.size());
    std::unordered_map<std::string, uint32_t> gidx_map;
    for(auto& g : genes) {
        gidx_map[g.gene_id] = g.gidx;
    }
    std::string   line;
    std::ifstream in(mod_file.c_str());
    if(!in) {
        std::cerr << "Error opening " << mod_file << "\n";
        exit(1);
    }

    while(std::getline(in, line)) {
        if(line.empty() || line[0] == '#') {
            continue;
        }
        auto it = gidx_map.find(line);
        if(it != gidx_map.end()) {
            out[it->second] = true;
        } else if(warn) {
            std::cerr << "[warning] " << mod_file << " gene_id " << line
                      << " is not found in the gtf index\n";
        }
    }
}
