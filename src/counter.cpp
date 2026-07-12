// SPDX-License-Identifier: MIT
// Copyright (c) Gavin W. Wilson

#include "counter.hpp"
#include <iomanip>
#include <iostream>
#include <fstream>
#include "gzstream.hpp"
using namespace scdepth;

std::pair<int, int> extract_coords(const std::string& s) {
    // s_002um_00320_02575-1
    auto first = s.find('_', s.find("um_"));
    auto second = s.find('_', first + 1);
    auto dash = s.find('-', second + 1);

    if (first == std::string::npos ||
        second == std::string::npos ||
        dash == std::string::npos) {
        throw std::runtime_error("Unexpected barcode format: " + s);
    }

    int x = std::stoi(s.substr(first + 1, second - first - 1));
    int y = std::stoi(s.substr(second + 1, dash - second - 1));

    return {x, y};
}

void BarcodeCounter::init(const std::string & lib_string, bool fwd, const std::string & barcode_tag, 
        const std::string & barcode_re, const std::string & umi_tag,
        const std::string & sample_tag, const std::vector<std::string> & samples,
        const std::string & random_hex_re, const std::string & random_hex_value,
        uint32_t barcode_length, uint32_t umi_length)
{
    if(fwd) strand_ = StrandMode::TAG_FWD; 
    else    strand_ = StrandMode::TAG_REV;
    barcode_tag_[0] = barcode_tag[0]; barcode_tag_[1] = barcode_tag[1]; 
    umi_tag_[0] = umi_tag[0]; umi_tag_[1] = umi_tag[1]; 
    lib_string_ = lib_string;
    barcode_length_ = barcode_length;
    umi_length_ = umi_length;
    barcode_regex_str_ = barcode_re;
    random_hex_regex_str_ = random_hex_re;
    random_hex_value_ = random_hex_value;
    if(samples.empty() != sample_tag.empty()){
        std::cerr << "sample ids and sample_tag must both be empty or specified\n";
        exit(1);
    }
    sample_tag_[0] = sample_tag[0]; sample_tag_[1] = sample_tag[1]; 
    samples_ = samples;
    std::cout << "samples = " << samples.size() << " tag = " << sample_tag << "\n";
    for(auto & s : samples){
        sample_set_.insert(s);
    }
    //std::cout << "barcode tag = " << barcode_tag << " len = " << barcode_length << " regex = " << barcode_regex_str_ << "\n";
    try{
        barcode_re_ = std::regex(
                barcode_re,
                std::regex::ECMAScript | std::regex::optimize
                );
    }catch(const std::regex_error & e){
        std::cerr << "Error compiling barcode regex " << e.what() << " code = " << e.code() 
            << " regex = " << barcode_re 
            << "\n";
        exit(1);
    }

    if(!random_hex_regex_str_.empty()){
        try{
            random_hex_re_ = std::regex(
                    random_hex_re,
                    std::regex::ECMAScript | std::regex::optimize
                    );
        }catch(const std::regex_error & e){
            std::cerr << "Error compiling random_hex_re_ regex " << e.what() << " code = " << e.code() 
                << " regex = " << random_hex_re 
                << "\n";
            exit(1);
        }
    }
}

void BarcodeCounter::set_count_parameters(double min_gene, unsigned int min_gene_bases, 
        double min_exonic, unsigned int min_intronic, uint8_t min_qual,
        bool discard_unknown_juncs, bool probes){

    min_gene_ = min_gene;
    min_exonic_ = min_exonic;
    min_gene_bases_ = min_gene_bases;
    min_intronic_ = min_intronic;
    min_qual_ = min_qual;
    discard_unknown_juncs_ = discard_unknown_juncs;
    has_probes_ = probes;
}

bool BarcodeCounter::prepare_bam(const std::string & gtf, const std::string & in, const std::string & out,
        int threads, size_t max_tags, double max_tag_frac){

    bf_ = sam_open(in.c_str(), "r");
    if(bf_ == nullptr){
        std::cerr << "Error opening sam/bam/cram header\n";
        return false;
    }
    bh_ = sam_hdr_read(bf_);
    if(bh_ == nullptr){
        std::cerr << "Error reading sam/bam/cram header\n";
        sam_close(bf_);
        return false;
    }

    hts_set_threads(bf_, threads);

    genes_ = gtf2genes(gtf, bh_);
    if(genes_.empty()){
        std::cerr << "Error loading genes/exons from gtf\n";
        return false;
    }
    if(has_probes_){
        for(size_t i = 0; i < genes_.size(); i++){
            auto & g = genes_[i];
            if(g.gidx != i){
                std::cerr << "mismatch between i " << i << " gidx = " << g.gidx << "\n";
            }
            gene_map_[g.gene_id] = g.gidx;
        }
    }
    gtrees_ = genes2trees(genes_, bh_->n_targets);
    shards_.init(out, max_tags, max_tag_frac);
    shards_.gene_breaks() = gene_breaks(genes_, bh_->n_targets);
    done_ = false;
    bam_file_ = in;
    out_file_ = out;
    barcode_counts_.clear();
    full_ = TagCounter();
    ssum_.clear();
    shards_.reset();
    std::cout.imbue(std::locale(""));
    return true;
}

unsigned int scdepth::get_read_blocks(bam1_t * b, std::vector<Exon> & blocks){
    blocks.clear();
    int lft = b->core.pos;
    int rgt = lft;
    auto cig = bam_get_cigar(b);
    unsigned int rbases = 0;
    for(uint32_t i = 0; i < b->core.n_cigar; i++){
        auto op = bam_cigar_op(cig[i]);
        auto len = bam_cigar_oplen(cig[i]);
        switch(op){
            case BAM_CMATCH: case BAM_CEQUAL: case BAM_CDIFF:
                rbases += len;
                rgt += len;
                break;
            case BAM_CDEL:
                rbases += len;
                rgt += len;
                break;
            case BAM_CREF_SKIP:
                blocks.push_back(Exon(lft, rgt - 1));
                lft = rgt + len;
                rgt = lft;
                break;
            default:
                break;
        }
    }
    blocks.push_back(Exon(lft, rgt - 1));
    return rbases;
}

unsigned int scdepth::block_overlaps(const std::vector<Exon> & blocks, const std::vector<Exon> & exons){
    unsigned int bases = 0;
    auto it = blocks.begin();
    auto eit = exons.begin();

    while(eit != exons.end() && it != blocks.end()){
        if(eit->rgt < it->lft){
            eit++;   //exon before block
        }else if(it->rgt < eit->lft){
            it++;    //block before exon
        }else if(it->overlaps(*eit)){
            bases += std::min(it->rgt, eit->rgt) - std::max(it->lft, eit->lft) + 1;
            if(it->rgt < eit->rgt){
                it++; // finished the block
            }else if(eit->rgt < it->rgt){
                eit++; // finished the exon
            }else{
                eit++; it++; //finished both
            }
        }else{
            std::cout << "    This should never happen block = " << it->lft << " - " << it->rgt << " exon = " << eit->lft << " - " << eit->rgt << "\n";
            exit(1);
        }
    }
    return bases;
}

unsigned int scdepth::check_junctions(const std::vector<Exon> & blocks, const Gene::Juncs & junctions){
    unsigned int found = 0;
    for(size_t i = 1; i < blocks.size(); i++){
        Exon junc(blocks[i-1].rgt, blocks[i].lft);
        if(junctions.find(junc) != junctions.end()){
            found++;
        }
    }
    return found;
}

bool extract_flag_from_qname(std::string_view qname, const std::regex& pattern, std::string_view true_value, bool & error){
    std::match_results<std::string_view::const_iterator> match;

    if (!std::regex_search(qname.begin(), qname.end(), match, pattern) ||
        match.size() < 2)
    {
        error = true;
        return false;
    }

    error = false;
    return std::string_view(match[1].first, match[1].length()) == true_value;
}

size_t BarcodeCounter::process_reads(size_t chunk){
    if(done_) return 0;
    size_t processed = 0;
    bam1_t * rec = bam_init1();
    std::vector<size_t> overlaps;
    std::string barcode;
    uint32_t umi = 0;
    bool umi_okay = false;
    bool track_primer_type = !random_hex_regex_str_.empty();
    //std::vector<Exon> tmp;
    while(processed < chunk){
        auto ret = sam_read1(bf_, bh_, rec);
        TagCounter * samp_sum = nullptr;
        if(ret < -1){
            std::cerr << "[error] Error processing the sam/bam/cram file\n";
            done_ = true;
            return 0;
        }else if(ret == -1){
            done_ = true;
            break;
        }
        processed++;
        full_.inc(TagSummary::TOTAL_READS, false, false);

        auto ptr = bam_aux_get(rec, barcode_tag_);
        char * BC = ptr == NULL ? NULL : bam_aux2Z(ptr);
        char * sample = nullptr;
        ptr = bam_aux_get(rec, umi_tag_);
        char * UMI = ptr == NULL ? NULL : bam_aux2Z(ptr);
        if(UMI == NULL || BC == NULL){
            full_.inc(TagSummary::BAD_TAGS, false, false);
            continue;
        }

        auto cb_len = strlen(BC);
        auto umi_len = strlen(UMI);
        int random_hex = 0;

        if(umi_length_ == 0){
            umi_length_ = umi_len;
        }
        if((barcode_length_ != 0 && cb_len != barcode_length_) || umi_len != umi_length_){
            full_.inc(TagSummary::BAD_TAGS, false, false);
            //std::cout << "bad tag length cb = " << cb_len << " vs " << barcode_length_ 
            //    << " umi = " << umi_len << " vs " << umi_length_ << "\n";
            continue;
        }

        if(!std::regex_match(BC, barcode_re_)){
            full_.inc(TagSummary::BAD_TAGS, false, false);
            //std::cout << "bad CB RE match CB = " << BC << "\n";
            continue;
        }

        if(!sample_set_.empty()){
            ptr = bam_aux_get(rec, sample_tag_);
            sample = ptr == NULL ? NULL : bam_aux2Z(ptr);
            if(sample == NULL || sample_set_.find(sample) == sample_set_.end()){
                full_.inc(TagSummary::BAD_TAGS, false, false);
                continue;
            }
            barcode.assign(sample, strlen(sample));
            barcode += '_';
            barcode.append(BC, BC + strlen(BC));
            samp_sum = &ssum_[sample];
            samp_sum->inc(TagSummary::TOTAL_READS, false, false);

        }else{
            barcode.assign(BC, BC + strlen(BC));
        }


        if(track_primer_type){
            bool hex_error = false;
            random_hex = extract_flag_from_qname(
                bam_get_qname(rec),
                random_hex_re_,
                random_hex_value_,
                hex_error
            );
            if(hex_error){
                full_.inc(BAD_TAGS, false, false);
                if(samp_sum != nullptr) samp_sum->inc(TagSummary::BAD_TAGS, false, false);
                continue;
            }
            full_.inc_pa_hex(TOTAL_READS, true, random_hex);
            //increment the total reads since they weren't counted above
            if(samp_sum != nullptr) samp_sum->inc_pa_hex(TagSummary::TOTAL_READS, true, random_hex);
        }


        std::tie(umi, umi_okay) = seq2int(UMI, umi_length_);
        if(!umi_okay){
            //std::cout << "bad UMI UMI = " << UMI << "\n";
            full_.inc(TagSummary::BAD_TAGS, track_primer_type, random_hex);
            if(samp_sum) samp_sum->inc(TagSummary::BAD_TAGS, track_primer_type, random_hex);
            continue;
        }

        auto res= barcode_map_.insert({barcode, barcode_map_.size()});
        size_t bidx = res.first->second;
        if(barcode_counts_.size() <= bidx) barcode_counts_.resize(bidx + 1);
        barcode_counts_[bidx].total++;
        if(sample && barcode_counts_[bidx].sample.empty()){
            barcode_counts_[bidx].sample = sample;
        }

        if(rec->core.tid < 0 || rec->core.flag & BAM_FUNMAP || rec->core.qual < min_qual_){
            full_.inc(TagSummary::LOW_QUALITY, track_primer_type, random_hex);
            if(samp_sum) samp_sum->inc(TagSummary::LOW_QUALITY, track_primer_type, random_hex);
            continue;
        }
        auto ref = gtrees_[rec->core.tid]; 
        //ref not in gtf
        uint32_t lft = rec->core.pos;
        uint32_t rgt = bam_endpos(rec) - 1;
        char xs = read2strand_(rec);
        //std::cout << "Read overlaps = " << overlaps.size() << " xs = " << xs << " pos = " << lft << " - " << rgt << "\n";

        if(has_probes_){
            auto ptr = bam_aux_get(rec, "GX");
            char * GX = ptr == NULL ? NULL : bam_aux2Z(ptr);
            uint32_t gidx = 0;
            //std::cout << "check GX = " << GX << "\n";
            if(GX == NULL || strchr(GX, ';') != NULL){
                full_.inc(TagSummary::NO_GENE, track_primer_type, random_hex);
                if(samp_sum) samp_sum->inc(TagSummary::NO_GENE, track_primer_type, random_hex);
                continue;
            }else{
                auto git = gene_map_.find(GX);
                if(git != gene_map_.end()){
                    gidx = git->second;
                }else{
                    full_.inc(TagSummary::NO_GENE, track_primer_type, random_hex);
                    if(samp_sum) samp_sum->inc(TagSummary::NO_GENE, track_primer_type, random_hex);
                    continue;
                }
            }
            RawTag tag;
            tag.make_tag(bidx, gidx, umi, 0, 0, random_hex);
            full_.inc(TagSummary::AMBIGUOUS_READS, track_primer_type, random_hex);
            full_.inc(TagSummary::COUNTABLE_READS, track_primer_type, random_hex);
            if(samp_sum) {
                samp_sum->inc(TagSummary::COUNTABLE_READS, track_primer_type, random_hex);
                samp_sum->inc(TagSummary::AMBIGUOUS_READS, track_primer_type, random_hex);
            }
            barcode_counts_[bidx].countable++;
            if(random_hex) barcode_counts_[bidx].random_hex++;
            else           barcode_counts_[bidx].poly_a++;
            shards_.add_tag(rec->core.tid, rec->core.pos, tag);
            continue;
        }

        if(!ref.overlap(lft, rgt, overlaps)) {
            full_.inc(TagSummary::NO_GENE, track_primer_type, random_hex);
            if(samp_sum) samp_sum->inc(TagSummary::NO_GENE, track_primer_type, random_hex);
            continue;
        }
        //std::cout << "check overlaps\n";
        size_t rbases = get_read_blocks(rec, read_blocks_);
        //if(read_blocks_.size() < 2) continue;

        size_t sgidx = 0;
        size_t agidx = 0;
        size_t ugidx = 0;
        size_t scount = 0, ucount = 0, acount = 0;
        //std::cout << "Read blocks = " << read_blocks_.size() << "\n";
        for(auto o : overlaps){
            size_t idx = o;
            idx = ref.data(o);
            const Gene & g = genes_[idx];
            //std::cout << "  Potential gene " << g.gene_id << " strand = " << g.strand << " pos = " << g.lft << " - " << g.rgt << "\n";
            /* 
            if(read_blocks_.size() > 1){
                std::cout << " Junctions passed = " << check_junctions(read_blocks_, g->junctions);
            }
            std::cout << " junctions: ";
            for(auto & e : g->junctions){
                tmp.push_back(e);
            }
            std::sort(tmp.begin(), tmp.end());
            for(auto & e : tmp){
                std::cout << " (" << e.lft << " - " << e.rgt << ") ";
            }
            std::cout << " exons: ";
            for(auto & e : g->exons){
                std::cout << " (" << e.lft << " - " << e.rgt << ") ";
            }
            std::cout << "\n";
            */
            if(g.strand != xs) continue;

            unsigned int gbases = 0;
            unsigned int in_gene = 0;
            //std::cout << "Checking block overlap\n";
            for(auto & b : read_blocks_){
                if(b.overlaps(g.lft, g.rgt)){
                    gbases += (std::min(b.rgt, g.rgt) - std::max(b.lft, g.lft)) + 1;
                    in_gene++;
                }
            }
            if(gbases < min_gene_bases_ || (1.0 * gbases / rbases) < min_gene_ 
                    || in_gene != read_blocks_.size()) continue;

            //std::cout << "Checking junctions\n";
            unsigned int jcnt = 0;
            if(read_blocks_.size() > 1){
                jcnt = check_junctions(read_blocks_, g.junctions);
                if(discard_unknown_juncs_ && jcnt != (read_blocks_.size() -1)){
                    continue;
                }
            }
            //std::cout << "Checking exon overlaps\n";
            auto eoverlaps = block_overlaps(read_blocks_, g.exons);
            if(eoverlaps > gbases){
                std::cerr << "[error] PROBLEM this should not happen!!!\n";
                exit(1);
            }
            auto ioverlaps = gbases - eoverlaps;
            double efrac = 1.0 * eoverlaps / gbases;
            /*
            std::cout << std::fixed << std::setprecision(3) 
                << "    gid = " << g.gidx << " gene_id = " << g.gene_id << " gbases = " << gbases 
                << " in_gene = " << in_gene << " gene frac = " << (1.0 * gbases / rbases)
                << " ioverlaps = " << ioverlaps << " blocks = " << read_blocks_.size() << " efrac = " << min_exonic_
                << "\n";
            */
            if(ioverlaps >= min_intronic_ && read_blocks_.size() == 1){
                ucount++;
                ugidx = g.gidx;
            }else if(efrac >= min_exonic_){
                if(jcnt > 0){
                    scount++;
                    sgidx = g.gidx;
                }else{
                    acount++;
                    agidx = g.gidx;
                }
            }
        }
        RawTag tag;
        if((scount + acount) == 1){
            if(random_hex) barcode_counts_[bidx].random_hex++;
            else           barcode_counts_[bidx].poly_a++;
            barcode_counts_[bidx].countable++;
            full_.inc(TagSummary::COUNTABLE_READS, track_primer_type, random_hex);
            if(samp_sum) samp_sum->inc(TagSummary::COUNTABLE_READS, track_primer_type, random_hex);
            if(scount == 1){
                tag.make_tag(bidx, sgidx, umi, 1, 0, random_hex);
                //std::cout << " spliced = " << sgidx;
                full_.inc(TagSummary::SPLICED_READS, track_primer_type, random_hex);
                if(samp_sum) samp_sum->inc(TagSummary::SPLICED_READS, track_primer_type, random_hex);
            }else{
                tag.make_tag(bidx, agidx, umi, 0, 0, random_hex);
                //std::cout << " ambiguous = " << agidx;
                full_.inc(TagSummary::AMBIGUOUS_READS, track_primer_type, random_hex);
                if(samp_sum) samp_sum->inc(TagSummary::AMBIGUOUS_READS, track_primer_type, random_hex);
            }
        }else if((scount + acount) == 0 && ucount == 1){
            tag.make_tag(bidx, ugidx, umi, 0, 1, random_hex);
            full_.inc(TagSummary::UNSPLICED_READS, track_primer_type, random_hex);
            //std::cout << " unspliced = " << ugidx;
            full_.inc(TagSummary::COUNTABLE_READS, track_primer_type, random_hex);
            if(samp_sum) {
                samp_sum->inc(TagSummary::COUNTABLE_READS, track_primer_type, random_hex);
                samp_sum->inc(TagSummary::UNSPLICED_READS, track_primer_type, random_hex);
            }
            barcode_counts_[bidx].countable++;
            if(random_hex) barcode_counts_[bidx].random_hex++;
            else           barcode_counts_[bidx].poly_a++;
        }else{
            full_.inc(TagSummary::NO_GENE, track_primer_type, random_hex);
            if(samp_sum) samp_sum->inc(TagSummary::NO_GENE, track_primer_type, random_hex);
            continue;
        }

        shards_.add_tag(rec->core.tid, lft, tag);
    }

    if(done_){
        sam_close(bf_);
        bam_hdr_destroy(bh_);
        bf_ = nullptr;
        bh_ = nullptr;
    }
    bam_destroy1(rec);
    return processed;
}

bool BarcodeCounter::finish(){
    if(!done_){
        std::cerr << "[warn] The bam file has unprocessed reads!\n";
    }

    if(!shards_.merge_shards(barcode_counts_)){
        std::cerr << "[errror] Error merging " << shards_.size() << " shard files\n";
        return false;
    }
    std::cout << "finish samples = " << samples_.size() << "\n";

    {
        gzofstream gzo(out_file_ + "_barcode_index.txt.gz");
        gzo << "barcode\tbarcode_idx\ttotal_reads\tcountable_reads\traw_molecules\toffset\ttotal_bytes\tsample\tcountable_random_hex\tcountable_poly_A\n";

        std::vector<std::string> bseqs;
        bseqs.resize(barcode_map_.size());
        for(auto & p : barcode_map_){
            bseqs[p.second] = p.first;
        }
        for(auto  & b : bseqs){
            size_t bidx = barcode_map_[b];
            gzo << b 
                << "\t" << bidx 
                << "\t" << barcode_counts_[bidx].total 
                << "\t" << barcode_counts_[bidx].countable 
                << "\t" << barcode_counts_[bidx].raw_molecules 
                << "\t" << barcode_counts_[bidx].offset 
                << "\t" << barcode_counts_[bidx].total_data_bytes 
                << "\t" << barcode_counts_[bidx].sample 
                << "\t" << barcode_counts_[bidx].random_hex 
                << "\t" << barcode_counts_[bidx].poly_a 
                << "\n";
        }
    }

    {
        std::ofstream os(out_file_ + "_summary.txt");
        os << "key\ttype\tvalue\n";
        for(size_t i = 0; i < TagSummary::TOTAL_FIELDS; i++){
            os << TagNames[i] << "\treads\t" << full_.merged_counts[i] << "\n";
        }
        if(!random_hex_regex_str_.empty()){
            for(size_t i = 0; i < TagSummary::TOTAL_FIELDS; i++){
                os << "polyA_" << TagNames[i] << "\treads\t" << full_.polyA_counts[i] << "\n";
            }
            for(size_t i = 0; i < TagSummary::TOTAL_FIELDS; i++){
                os << "random_hex_" << TagNames[i] << "\treads\t" << full_.random_hex_counts[i] << "\n";
            }
        }
        os << "raw_molecules\treads\t" <<  shards_.raw_molecules() << "\n";
        os << "bam\tparam\t" << bam_file_ << "\n";
        os << "min_gene\tparam\t" << std::fixed << std::setprecision(4) << min_gene_ << "\n";
        os << "min_exonic\tparam\t" << std::fixed << std::setprecision(4) << min_exonic_ << "\n";
        os << "min_gene_bases\tparam\t" << min_gene_bases_ << "\n";
        os << "min_intronic\tparam\t" << min_intronic_ << "\n";
        os << "min_qual\tparam\t" << static_cast<int>(min_qual_) << "\n";
        os << "discard_unknown_juncs\tparam\t" << static_cast<int>(discard_unknown_juncs_) << "\n";
        os << "library_string\tlibrary\t" << lib_string_ << "\n";
        os << "barcode_tag\tlibrary\t" << barcode_tag_[0] << barcode_tag_[1] << "\n";
        os << "barcode_length\tlibrary\t" << barcode_length_ << "\n";
        os << "barcode_regex\tlibrary\t" << barcode_regex_str_ << "\n";
        os << "umi_tag\tlibrary\t" << umi_tag_[0] << umi_tag_[1] << "\n";
        os << "umi_length\tlibrary\t" << umi_length_ << "\n";
        os << "strand_mode\tlibrary\t" << (strand_ == StrandMode::TAG_FWD ? "fwd" : "rev") << "\n";
        os << "barcodes_detected\tlibrary\t" << barcode_counts_.size() << "\n";
    }

    std::cout << "samples = " << samples_.size() << " ssum size = " << ssum_.size() << "\n";
    if(!samples_.empty()){
        std::ofstream os(out_file_ + "_sample_summary.txt");
        os << "key\ttype\tvalue\n";
        for(auto & s : samples_){
            const TagCounter & ss = ssum_[s];
            for(size_t i = 0; i < TagSummary::TOTAL_FIELDS; i++){
                os << s << "_" << TagNames[i] << "\treads\t" << ss.merged_counts[i] << "\n";
            }
            if(!random_hex_regex_str_.empty()){
                for(size_t i = 0; i < TagSummary::TOTAL_FIELDS; i++){
                    os << s << "_" << "polyA_" << TagNames[i] << "\treads\t" << ss.polyA_counts[i] << "\n";
                }
                for(size_t i = 0; i < TagSummary::TOTAL_FIELDS; i++){
                    os << s << "_" << "random_hex_" << TagNames[i] << "\treads\t" << ss.random_hex_counts[i] << "\n";
                }
            }
        }
    }

    gtf2mapping(genes_, out_file_ + "_genes.txt.gz");

    return true;

}
