#include "parblockaligner.hpp"
#include "../external/threadpool/include/scheduler.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#define EXTENSION_LENGTH 10
#define MAX_ARGS 100
#define MAX_CIGAR_LEN 10000
#define _DEBUG_BLOCK 0

extern std::string tmp_path_prefix;
extern char** seq, **comment, *thread_str;
extern int block_sequence_type, threads_num;
extern std::mt19937 rand_gen; // not thread-safe random generator
extern item_int first_end;
extern int_ threshold;

#pragma GCC push_options
#pragma GCC optimize ("O0")
int Frag_dispatch_with_lock(const char* File, item_int& sequence_id, item_int& ll, item_int& rr, gap_tree &this_tree, std::mutex *lock)
{
    if(sequence_id.size() <= 1) return 0;
    thread_local FILE* f;
    f = fopen(File, "w");
    if (f == NULL) { fprintf(stderr, "Error: Cannot open temp file %s. Program will exit.\n", File); exit(1); }
    for (size_t i = 0; i < sequence_id.size(); ++i)
        WriteSeq_SSW_Frag(f, ll[i], rr[i], sequence_id[i]);
    fclose(f);
    std::string outputfile, backupfile;
    outputfile += File;
    outputfile += ".out";
    backupfile = outputfile + ".bak";
#if _DEBUG
    f = fopen(backupfile.c_str(), "w");
    if (f == NULL) { fprintf(stderr, "Error: Cannot open backup temp file %s. Program will exit.\n", File); exit(1); }
    for (size_t i = 0; i < sequence_id.size(); ++i)
        WriteSeq_SSW_Frag(f, ll[i], rr[i], sequence_id[i]);
    fclose(f);
#endif
    lock->lock();
    if (Fragalignment(File, outputfile.c_str())) { fprintf(stderr, "Error: can not execute alignment. Program will exit.\n"); exit(1); }
    lock->unlock();
#if _DEBUG
    if (Filedelete(File) || Filemove(outputfile.c_str(), File)) { fprintf(stderr, "Error: can not move file to source name. Program will exit.\n"); exit(1); }
    if (read_gaps_from_file(File, sequence_id, ll, rr, this_tree)) { fprintf(stderr, "Error: can not read sequence info. Program will exit.\n"); exit(1); }
#else
    if (read_gaps_from_file(outputfile.c_str(), sequence_id, ll, rr, this_tree)) { fprintf(stderr, "Error: can not read sequence info. Program will exit.\n"); exit(1); }
    if (Filedelete(File) || Filedelete(outputfile.c_str())) { fprintf(stderr, "Warning: can not remove alignment file.\n"); }
#endif
    return 0;
}

int Center_block_dispatch_with_lock(const char* File, center_int& center_info, gap_tree& this_tree, std::mutex *lock)
{
    thread_local FILE* f;
    f = fopen(File, "w");
    if (f == NULL) { fprintf(stderr, "Error: Cannot open temp file %s. Program will exit.\n", File); exit(1); }
    for (auto& item : center_info)
        WriteSeq_SSW_Frag(f, item.second.first, item.second.second, item.first);
    fclose(f);
    std::string outputfile, backupfile;
    outputfile += File;
    outputfile += ".out";
    backupfile = outputfile + ".bak";
#if _DEBUG
    f = fopen(backupfile.c_str(), "w");
    if (f == NULL) { fprintf(stderr, "Error: Cannot open backup temp file %s. Program will exit.\n", File); exit(1); }
    for (auto& item : center_info)
        WriteSeq_SSW_Frag(f, item.second.first, item.second.second, item.first);
    fclose(f);
#endif
    lock->lock();
    if (Fragalignment(File, outputfile.c_str())) { fprintf(stderr, "Error: can not execute alignment. Program will exit.\n"); exit(1); }
    lock->unlock();
#if _DEBUG
    if (Filedelete(File) || Filemove(outputfile.c_str(), File)) { fprintf(stderr, "Error: can not move file to source name. Program will exit.\n"); exit(1); }
    if (read_gaps_from_file(File, center_info, this_tree)) { fprintf(stderr, "Error: can not read sequence info. Program will exit.\n"); exit(1); }
#else
    if (read_gaps_from_file(outputfile.c_str(), center_info, this_tree)) { fprintf(stderr, "Error: can not read sequence info. Program will exit.\n"); exit(1); }
    if (Filedelete(File) || Filedelete(outputfile.c_str())) { fprintf(stderr, "Warning: can not remove alignment file.\n"); }
#endif
    return 0;
}

int TwoProfile_dispatch_with_lock(item_int &dispatch_ids, item_int &other_ids, item_int &ll, item_int &rr, map_int &mp,
                                  gap_tree &this_tree, const int depth, const unsigned int part, std::mutex *lock)
{
    // write File name
    std::string File1, File2, File3;
    File1 = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_"; File3 = File2 = File1;
    File1 += "A.fasta"; File2 += "B.fasta"; File3 += "res.fasta";

    // Write common sequences into file1
    TwoProfile_writeseq(File1.c_str(), dispatch_ids, ll, rr, mp, this_tree);

    // Write common sequences into file2
    TwoProfile_writeseq(File2.c_str(), other_ids,    ll, rr, mp, this_tree);

    // call two_profile_align
    lock->lock();
    if(TwoProfileAlignment(File1.c_str(), File2.c_str(), File3.c_str()))
    { 
        fprintf(stderr, "Error: can not two profile alignment. Program will exit.\n"); 
        exit(1); 
    };
    lock->unlock();
    // read result from two profiles
    if (read_gaps_from_merged_profile(File3.c_str(), dispatch_ids, other_ids, ll, rr, mp, this_tree))
    {
        fprintf(stderr, "Error: can not read sequence info. Program will exit.\n");
        exit(1);
    }
    return 0;
}

#pragma GCC pop_options

/*
 The find similarity function.
 * Use ssw program for finding the items in seq_list called need_align_seq_id, and made need_align_seq_id-th sequence segment
   in [need_align_center - (template_seq_len + ext_size), need_align_center + (template_seq_len + ext_size)]
 * modifing args:
 - real_other_ids: a set, contains failed to find ids;
 - center_block:   a map, first is the id of successfully to find, second is tuple, contains start place and end place
 */
void Findsim_and_extension_with_lock(item_int& seq_list, set_int& not_same_len_list, set_int& real_ids, int_ need_align_center,
    item_int& all_ll, item_int& all_rr, map_int& ll_rr_map,
    int_ template_seq_id, int_ template_seq_begin, int_ template_seq_end, center_int &center_block, set_int &real_other_ids,
    gap_tree &gap_tree_main, int depth, unsigned int part, std::mutex *lock)
{
    FILE* f_seq;
    // calling ssw function
    int8_t* table, * mat, * ref_int;
    int32_t nmatrix, ref_len;

    if (block_sequence_type == 1) // Protein
    {
        table = aa_table;
        mat = blosum62;
        nmatrix = 24;
    }
    else if (block_sequence_type == 2) //DNA
    {
        table = nt_table;
        mat = mat_nt;
        nmatrix = 5;
    }
    else { fprintf(stderr, "Error: sequence type is error. Program will exit.\n"); exit(1); }

    // first time for ref, without gaps in this system
    ref_len = template_seq_end - template_seq_begin;
    ref_int = (int8_t*)malloc(sizeof(int8_t) * (size_t)(ref_len + 1));
    if (ref_int == NULL) { fprintf(stderr, "\nError: Can not alloc sequence. Program will exit\n"); exit(1); }
    for (size_t i = 0; i < ref_len; ++i) ref_int[i] = table[(size_t)seq[template_seq_id][template_seq_begin + i]];

    // make alignment using ssw
    s_profile* profile = ssw_init(ref_int, ref_len, mat, nmatrix, 2);
    int gapopen = 3, gapext = 1;
    bool foundsw = false;
    int_ ssw_len_half = ref_len >> 1;
#if 0
#pragma loop(hint_parallel(12))
#elif defined(_OPENMP)
#pragma omp parallel for num_threads(threads_num)
#endif
    for (int i = 0; i < seq_list.size(); ++i)
    {
        int_* real_id = &seq_list[i];
        int_ now_id = ll_rr_map[*real_id];
        if (all_ll[now_id] == all_rr[now_id]) // no length
        {
            // processed before, so no need to process
            continue;
        }
        int_ ll = all_ll[now_id], rr = all_rr[now_id];
#if _DEBUG_BLOCK
        assert(ll <= rr && all_ll[now_id] <= ll && rr <= all_rr[now_id]); // may be has no length?
        fprintf(stderr, "Info: SSW Program Border = [%d, %d), find similar sequence in [%d, %d), sequence_id = %d\n",
                  all_ll[now_id], all_rr[now_id], ll, rr, *real_id);
#endif
        int8_t* seq_ssw = (int8_t*)malloc(sizeof(int8_t) * (rr - ll));
        if (seq_ssw == NULL) { fprintf(stderr, "Error: SSW can not alloc enough memory. Program will exit.\n"); exit(1); }
        for (int_ i = ll; i < rr; ++i) seq_ssw[i - ll] = table[(size_t)seq[*real_id][i]];
        s_align* alignment_result = ssw_align(profile, seq_ssw, rr - ll, gapopen, gapext, 1, 0, 0, 15);
        if (alignment_result->ref_begin1 == -1 || alignment_result->read_begin1 == -1 || alignment_result->cigarLen == 0)
        {
#if _DEBUG
            fprintf(stderr, "Warning: not found in sequence %d\n", *real_id);
#endif
#ifdef _OPENMP
#pragma omp critical
#endif
            real_other_ids.insert(*real_id);
            align_destroy(alignment_result);
            free(seq_ssw);
            continue;
        }
#ifdef _OPENMP
#pragma omp critical
#endif
        center_block[*real_id] = std::make_pair(ll + alignment_result->ref_begin1, ll + alignment_result->ref_end1 + 1);
#if _DEBUG
        assert(ll + alignment_result->ref_end1 + 1 <= all_rr[now_id]);
#endif
        foundsw = true;
        align_destroy(alignment_result);
        free(seq_ssw);
    }
    init_destroy(profile);
    free(ref_int);
    if(foundsw || not_same_len_list.size()) // maybe not found
    {
        std::string depth_part, seq_file_name, out_file_name;
        depth_part = tmp_path_prefix; depth_part += std::to_string(depth); depth_part += "_"; depth_part += std::to_string(part);
        // write center sequence line and common sequences and ids for output and align
        item_int vll, vrr, seq_frags_ids;
        vll.reserve(seq_list.size() + 1); vrr.reserve(seq_list.size() + 1); seq_frags_ids.reserve(seq_list.size() + 1);
        seq_frags_ids.emplace_back(template_seq_id);
        vll.emplace_back(template_seq_begin); vrr.emplace_back(template_seq_end);
        for (auto& real_id : seq_list)
        {
            if (real_other_ids.find(real_id) != real_other_ids.end()) continue;
            auto& [ll, rr] = center_block[real_id];
            vll.emplace_back(ll); vrr.emplace_back(rr); seq_frags_ids.emplace_back(real_id);
        }
        // add sequences which has not same block length
        for (auto &real_id: not_same_len_list)
        {
            auto& [ll, rr] = center_block[real_id];
            vll.emplace_back(ll); vrr.emplace_back(rr); seq_frags_ids.emplace_back(real_id);
        }
#if _DEBUG
        seq_file_name = depth_part; seq_file_name += "SSW_RAW.fasta";
        f_seq = fopen(seq_file_name.c_str(), "w");
        if(f_seq == NULL) { fprintf(stderr, "Can not open file %s (common file). Program will exit.\n", seq_file_name.c_str()); exit(1); }
        for (int i = 0; i < seq_frags_ids.size(); ++ i)
            WriteSeq_SSW_Frag(f_seq, vll[i], vrr[i], seq_frags_ids[i]);
        fclose(f_seq);
#endif

        int_ template_seq_len = template_seq_end - template_seq_begin;
        out_file_name = depth_part + ".SSW_Align_out.fasta";
        Frag_dispatch_with_lock(out_file_name.c_str(), seq_frags_ids, vll, vrr, gap_tree_main, lock);

        // the result has been written in template_seq_id, let us write it!
        auto &gap_template = gap_tree_main.get_gaps(template_seq_id);
        for(auto &real_id: real_ids)
        {
            if(real_id == template_seq_id) continue; // this sequence is written by frag
            auto& gap_now = gap_tree_main.get_gaps(real_id);
            auto &[this_ll, this_rr] = center_block[real_id];
            int_ this_gaps, delta;
            for(delta = -1; delta < template_seq_len; ++delta)
            {
                this_gaps = gap_template.get_gaps(template_seq_begin + delta);
                gap_now.add_gaps(this_ll + delta, this_gaps);
            }
        }
    }
    else
    {
        if(! foundsw) fprintf(stderr, "Warning: not found vaild part for alignment.\n");
    }
}

class decompose_Task
{
public:
    decompose_Task(item_int& sequence_ids, block_irreg& center_block, std::vector<block_irreg>& other_blocks, item_int& ll, item_int& rr, int depth,
                   unsigned int part, map_int& cluster_mp_ids, gap_tree& this_tree, std::mutex *lock, std::string& final_name) :
        sequence_ids(sequence_ids), center_block(center_block), other_blocks(other_blocks), ll(ll), rr(rr), depth(depth), part(part),
        cluster_mp_ids(cluster_mp_ids), this_tree(this_tree), lock(lock), final_name(final_name) {}
    item_int& sequence_ids;
    block_irreg& center_block;
    std::vector<block_irreg>& other_blocks;
    item_int& ll, &rr;
    int depth;
    unsigned int part;
    map_int& cluster_mp_ids;
    gap_tree& this_tree;
    std::mutex *lock;
    std::string& final_name;
};

class TaskRunnerSSW: public staccato::task<TaskRunnerSSW>
{
public:
    TaskRunnerSSW(decompose_Task task_) : task_(task_) {}
    void execute();
private:
    decompose_Task task_;
};


unsigned int rand_gen_lock()
{
    static std::mutex _;
    std::lock_guard<std::mutex> __(_);
    return rand_gen();
}

/*
 * The decopmose function. This function will decompose seqeuences into *three* parts like this.
 * |||||||||||||||||||||||.||---------------------------||.|||||||||||||||||||||
 * |||||||-----------|||||.||---------------------------||.|||||-----------|||||
 * |||||||-*********-|||||.||---------------------------||.|||||-*********-|||||
 * |||||||-*********-|||||.||---------------------------||.|||||-*********-|||||
 * |||||||-*********-|||||.||---------------------------||.|||||-*********-|||||
 * |||||||-----------|||||.||---------------------------||.|||||-----------|||||
 * |||||||||||||||||||||||.||---------------------------||.|||||||||||||||||||||
 * .............................................................................
 * |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 * |||||||-----------|||||||||||||||||||||||||||||||||||||||||||-----------|||||
 * |||||||-*********-|||||||||||||||||||||||||||||||||||||||||||-|||||||||-|||||
 * |||||||-*********-|||||||||||||||||||||||||||||||||||||||||||-|||||||||-|||||
 * |||||||-----------|||||||||||||||||||||||||||||||||||||||||||-----------|||||
 * with SSW
 * Return the merged file name. (may not need to return?)
 */

void TaskRunnerSSW::execute()
{
    auto &[sequence_ids, center_block, other_blocks, ll, rr, depth, part, cluster_mp_ids, this_tree, lock, final_name] = task_;
    if(center_block.container.size() <= 0)
    {
        final_name = "";
        return;
    }
    // Part 0: sort the blocks
    std::sort(other_blocks.begin(), other_blocks.end(), std::greater<block_irreg>());
    // Part 1: left_blocks = [], right_blocks = [], down_blocks = []
    std::vector <block_irreg> left_blocks, right_blocks, down_blocks;
    item_int other_ids, * real_other_ids = nullptr, *real_dispatch_lr = nullptr;
    item_int ll_dispatch_ll, ll_dispatch_rr, rr_dispatch_ll, rr_dispatch_rr, other_dispatch_ll, other_dispatch_rr;
    set_int other_set, center_items, all_set, zero_set;
    center_int center_block_info;
    map_int mp_ids; // not in other blocks: real_id->id
    block_irreg* left_block, * right_block, * down_block;
    std::string file_left, file_right, file_down;
    gap_load_int gap_preload_in_main;
    int_ center_block_len = center_block.length, max_len = -1, min_len = INT_MAX;
    gap_tree center_tree;
    bool use_center_tree = false;
    // Part 2: check the number in center_block is same as the cluster
    for (size_t i = 0; i < sequence_ids.size(); ++i)
    {
        mp_ids[sequence_ids[i]] = i;
        if(ll[i] == rr[i]) zero_set.insert(sequence_ids[i]);
        else all_set.insert(sequence_ids[i]);
    }
    for (auto& item : center_block.container)
    {
        auto& [real_id, block_ll, delta] = item;
        center_block_info[real_id] = std::make_pair(block_ll, block_ll + center_block_len - delta);
        max_len = (std::max)(center_block_len - delta, max_len);
        min_len = (std::min)(center_block_len - delta, min_len);
    }
    if(max_len <= 0)
    {
        fprintf(stderr, "Warning: this block is invalid. Program will return.\n");
        std::string final_file_name = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_final.fasta";
        Frag_dispatch_with_lock(final_file_name.c_str(), sequence_ids, ll, rr, this_tree, lock);
        final_name = "";
        return ; // can not find invalid block
    }
    if (sequence_ids.size() != center_block_info.size())
    {
        assert(sequence_ids.size() > center_block_info.size());
        // Part 2.1: determine how many sequences are not in clusters, and try to find the similar parts by extending sequences into 10 bps
        int_ valid_beg, valid_end, valid_seq_id = -1; // [valid_beg, valid_end)
        fprintf(stderr, "Info: this block is not contained all sequences. Will change this block information.\n");
        set_int not_same_len_seq, tmp_; // center block sequence may have different length
        for(auto& item: center_block_info)
            if(max_len == item.second.second - item.second.first)
            {
                center_items.insert(item.first);
                if(valid_seq_id != -1) continue;
                valid_seq_id = item.first;
                std::tie(valid_beg, valid_end) = item.second;
            }
            else
            {
                not_same_len_seq.insert(item.first);
            }
        std::set_difference(all_set.begin(), all_set.end(), center_items.begin(), center_items.end(), std::inserter(tmp_, tmp_.begin()));
        std::set_difference(tmp_.begin(), tmp_.end(), not_same_len_seq.begin(), not_same_len_seq.end(), std::back_inserter(other_ids));
        all_set.clear(); tmp_.clear();
        use_center_tree = true;
        // Part 2.2: try to use s-w find sequence. if not found, insert these sequences into down_blocks list
        Findsim_and_extension_with_lock(other_ids, not_same_len_seq, center_items, center_block.center, ll, rr, mp_ids, valid_seq_id, valid_beg, valid_end,
                                        center_block_info, other_set, center_tree, depth, part, lock);
        // Part 2.3: determine how many sequences will be in dispatch mode
        real_dispatch_lr = new item_int;
        // Part 2.4: determine real_dispatch_lr = center_block_info.first
        real_dispatch_lr -> reserve(center_block_info.size());
        for(auto &item: center_block_info) real_dispatch_lr->emplace_back(item.first);
#if _DEBUG
        write_gaps_to_file_center((tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + ".ssw.fasta").c_str(), center_block_info, this_tree);
#endif
    }
    else
    {
        real_dispatch_lr = new item_int;
        for (auto& item : all_set) real_dispatch_lr->emplace_back(item);
        all_set.clear();
        if (min_len != max_len)
        {
            // block is irregular, align it!
            fprintf(stderr, "Info: this block is not same (although contains all sequences). Program will align this block.\n");
            std::string center_file_name = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_block.fasta";
            use_center_tree = true;
            Center_block_dispatch_with_lock(center_file_name.c_str(), center_block_info, center_tree, lock);
        }
    }
    // Part 3: determine other_blocks is in left or right or down, or just simply, drop it
    for (auto& block_item : other_blocks)
    {
        left_block = new block_irreg;
        right_block = new block_irreg;
        if (other_set.size()) down_block = new block_irreg;
        else down_block = nullptr;
        // check the place for every block
        for (auto& sequence_item : block_item.container)
        {
            auto &[this_id, this_ll, delta] = sequence_item;
            auto this_rr = this_ll + block_item.length - delta;
            auto mp_place = mp_ids[this_id];
            auto &border_ll = ll[mp_place], &border_rr = rr[mp_place];
            if(zero_set.find(this_id) != zero_set.end()) continue; // ignore
            if(down_block && other_set.find(this_id) != other_set.end()) // will be in other block
            {
                if(border_ll <= this_ll && this_rr <= border_rr && this_ll + threshold < this_rr)
                {
#if _DEBUG_BLOCK
                    fprintf(stderr, "Border = [%d, %d); ", border_ll, border_rr);
                    fprintf(stderr, "dispatch to other = [%d, %d)\n", this_ll, this_rr);
#endif
                    down_block->container.emplace_back(sequence_item);
                }
            }
            else
            {
                auto& [centerll, centerrr] = center_block_info[this_id];
                auto ll_border_ll = border_ll, ll_border_rr = centerll, rr_border_ll = centerrr, rr_border_rr = border_rr;
#if _DEBUG_BLOCK
                fprintf(stderr, "Center block = [%d, %d), will dispatch block = [%d, %d)\n", centerll, centerrr, this_ll, this_rr);
                fprintf(stderr, "Border = [%d, %d); left border = [%d, %d); right border = [%d, %d)\n",
                        border_ll, border_rr, ll_border_ll, ll_border_rr, rr_border_ll, rr_border_rr);
#endif
                if (centerll <= this_ll && this_rr <= centerrr) continue; // dropped
                else if (this_ll <= centerll && centerrr <= this_rr) // across the block
                {
#if _DEBUG_BLOCK
                    fprintf(stderr,
                            "Warning: sequence %d has condition:\ncenter block = [%d, %d), this block = [%d, %d)\n",
                            this_id, centerll, centerrr, this_ll, this_rr);
#endif
                    sequence_irreg sequence_item_backup;
                    // insert it to left block
                    if(ll_border_ll <= this_ll && this_ll + block_item.length - (this_rr - centerll) <= ll_border_rr &&
                       this_ll + threshold < this_ll + block_item.length - (this_rr - centerll))
                    {
                        sequence_item_backup = sequence_item;
                        std::get<2>(sequence_item_backup) += this_rr - centerll;
                        left_block->container.emplace_back(sequence_item_backup);
#if _DEBUG_BLOCK
                        fprintf(stderr, "Left = [%d, %d);\n", this_ll, this_ll + block_item.length - (this_rr - centerll));
#endif
                    }
                    // insert it to right block
                    if(rr_border_ll <= centerrr && center_block.length + this_ll <= rr_border_rr &&
                       centerrr + threshold < center_block.length + this_ll)
                    {
                        sequence_item_backup = sequence_item;
                        std::get<1>(sequence_item_backup) = centerrr;
                        std::get<2>(sequence_item_backup) += (centerrr - this_ll);
                        right_block->container.emplace_back(sequence_item_backup);
#if _DEBUG_BLOCK
                        fprintf(stderr, "Right = [%d, %d)\n", centerrr, centerrr + block_item.length - (centerrr - this_ll));
#endif
                    }
                }
                else if (this_ll <= centerll && centerll <= this_rr && this_rr <= centerrr)
                    // will be in left block, and end place is centerll
                {
                    if(ll_border_ll <= this_ll && centerll <= ll_border_rr && this_ll + threshold < centerll)
                    {
                        delta += (this_rr - centerll);
                        left_block->container.emplace_back(sequence_item);
#if _DEBUG_BLOCK
                        fprintf(stderr, "Border = [%d, %d); make changes ", ll_border_ll, ll_border_rr);
                        fprintf(stderr, "dispatch to left = [%d, %d); ", this_ll, this_ll + block_item.length - delta);
                        fprintf(stderr, "check = %d\n", centerll == this_ll + block_item.length - delta);
#endif
                    }
                }
                else if (centerll <= this_ll && this_ll <= centerrr && centerrr <= this_rr)
                    // will be in right block, and start place is centerrr
                {
                    if(rr_border_ll <= centerrr && this_rr <= rr_border_rr && centerrr + threshold < this_rr)
                    {
                        delta += (centerrr - this_ll);
                        this_ll = centerrr;
                        right_block->container.emplace_back(sequence_item);
#if _DEBUG_BLOCK
                        fprintf(stderr, "Border = [%d, %d); make changes ", rr_border_ll, rr_border_rr);
                        fprintf(stderr, "dispatch to right = [%d, %d); ", this_ll, this_ll + block_item.length - delta);
                        fprintf(stderr, "check = %d\n", this_rr == this_ll + block_item.length - delta);
#endif
                    }
                }
                else if (this_rr <= centerll) // will be in left block
                {
                    if(ll_border_ll <= this_ll && this_rr <= ll_border_rr)
                    {
                        left_block->container.emplace_back(sequence_item);
#if _DEBUG_BLOCK
                        fprintf(stderr, "Border = [%d, %d); ", ll_border_ll, ll_border_rr);
                        fprintf(stderr, "dispatch to left = [%d, %d)\n", this_ll, this_rr);
#endif
                    }
                }
                else if (centerrr <= this_ll) // will be in right block
                {
                    if(rr_border_ll <= this_ll && this_rr <= rr_border_rr)
                    {
                        right_block->container.emplace_back(sequence_item);
#if _DEBUG_BLOCK
                        fprintf(stderr, "Border = [%d, %d); ", rr_border_ll, rr_border_rr);
                        fprintf(stderr, "dispatch to right = [%d, %d)\n", this_ll, this_rr);
#endif
                    }
                }
                else
                {
                    fprintf(stderr,
                           "Designer Error: not considered condition: center = [%d, %d), this = [%d, %d)\n",
                            centerll, centerrr, this_ll, this_rr);
                    assert(false);
                }
            }
        }
        // insert left block to left blocks
        if (left_block->container.size() >= 2)
        {
            left_block->length = block_item.length;
            left_block->center = find_center(left_block->container);
            left_block->calc_area();
            left_blocks.emplace_back(std::move(*left_block));
        }
        else
        {
            left_block->freecontainer();
            delete left_block;
        }
        // insert right block to right blocks
        if (right_block->container.size() >= 2)
        {
            right_block->length = block_item.length;
            right_block->center = find_center(right_block->container);
            right_block->calc_area();
            right_blocks.emplace_back(std::move(*right_block));
        }
        else
        {
            right_block->freecontainer();
            delete right_block;
        }
        // insert down block to down blocks
        if (other_set.size())
        {
            if (down_block->container.size() >= 2)
            {
                down_block->length = block_item.length;
                down_block->center = find_center(down_block->container);
                down_block->calc_area();
                down_blocks.emplace_back(std::move(*down_block));
            }
            else
            {
                down_block->freecontainer();
                delete down_block;
            }
        }
    }
    // Part 4: decompose three parts: left_blocks, right_blocks, other_blocks
    this_tree.alloc_lr();
    // Please sort id and center block
    // Part 4.1: if can not be decomposed, call fragalign for making alignment
    // determine the dispatch_ll and dispatch_rr, and call for sub process
    block_irreg left_center, right_center, other_center;
    unsigned int rand_left, rand_right, rand_down;
    bool calc_left = true, calc_right = true, calc_down = true;
    if (left_blocks.size())
    {
        calc_left = false;
        std::sort(left_blocks.begin(), left_blocks.end());
        left_center = left_blocks.back();
        left_blocks.pop_back(); // remove last one
        ll_dispatch_ll.resize(real_dispatch_lr->size());
        ll_dispatch_rr.resize(real_dispatch_lr->size());
        for (size_t i = 0; i < real_dispatch_lr->size(); ++i)
        {
            ll_dispatch_rr[i] = center_block_info[real_dispatch_lr->at(i)].first;
            ll_dispatch_ll[i] = ll[mp_ids[real_dispatch_lr->at(i)]];
        }
        // dispatch to left next-level
#if _DEBUG
        rand_left = rand_gen_lock();
        fprintf(stderr, "Recurse: left decompose %d %u, dispatch to %u\n", depth, part, rand_left);
        spawn(new(child()) TaskRunnerSSW(decompose_Task(*real_dispatch_lr, std::ref(left_center), std::ref(left_blocks), std::ref(ll_dispatch_ll), std::ref(ll_dispatch_rr), depth + 1, rand_left,  std::ref(cluster_mp_ids), *this_tree.left_tree, std::ref(lock), std::ref(file_left))));
#else
        spawn(new(child()) TaskRunnerSSW(decompose_Task(*real_dispatch_lr, std::ref(left_center), std::ref(left_blocks), std::ref(ll_dispatch_ll), std::ref(ll_dispatch_rr), depth + 1, rand_gen_lock(), std::ref(cluster_mp_ids), *this_tree.left_tree, std::ref(lock), std::ref(file_left))));
#endif
    }
    else
    {
        file_left = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_left.fasta";
        ll_dispatch_ll.resize(real_dispatch_lr->size());
        ll_dispatch_rr.resize(real_dispatch_lr->size());
        for (size_t i = 0; i < real_dispatch_lr->size(); ++i)
        {
            ll_dispatch_ll[i] = ll[mp_ids[real_dispatch_lr->at(i)]];
            ll_dispatch_rr[i] = center_block_info[real_dispatch_lr->at(i)].first;
        }
    }

    if (right_blocks.size())
    {
        calc_right = false;
        std::sort(right_blocks.begin(), right_blocks.end());
        right_center = right_blocks.back();
        right_blocks.pop_back(); // remove last one
        rr_dispatch_ll.resize(real_dispatch_lr->size());
        rr_dispatch_rr.resize(real_dispatch_lr->size());
        for (size_t i = 0; i < real_dispatch_lr->size(); ++i)
        {
            rr_dispatch_ll[i] = center_block_info[real_dispatch_lr->at(i)].second;
            rr_dispatch_rr[i] = rr[mp_ids[real_dispatch_lr->at(i)]];
        }
#if _DEBUG
        rand_right = rand_gen_lock();
        fprintf(stderr, "Recurse: right decompose %d %u, dispatch to %u\n", depth, part, rand_right);
        spawn(new(child()) TaskRunnerSSW(decompose_Task(*real_dispatch_lr, std::ref(right_center), std::ref(right_blocks), std::ref(rr_dispatch_ll), std::ref(rr_dispatch_rr), depth + 1, rand_right, std::ref(cluster_mp_ids), *this_tree.right_tree, std::ref(lock), std::ref(file_right))));
#else
        spawn(new(child()) TaskRunnerSSW(decompose_Task(*real_dispatch_lr, std::ref(right_center), std::ref(right_blocks), std::ref(rr_dispatch_ll), std::ref(rr_dispatch_rr), depth + 1, rand_gen_lock(), std::ref(cluster_mp_ids), *this_tree.right_tree, std::ref(lock), std::ref(file_right))));
#endif
    }
    else
    {
        file_right = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_right.fasta";
        rr_dispatch_ll.resize(real_dispatch_lr->size());
        rr_dispatch_rr.resize(real_dispatch_lr->size());
        for (size_t i = 0; i < real_dispatch_lr->size(); ++i)
        {
            rr_dispatch_ll[i] = center_block_info[real_dispatch_lr->at(i)].second;
            rr_dispatch_rr[i] = rr[mp_ids[real_dispatch_lr->at(i)]];
        }
    }

    if (other_set.size())
    {
        real_other_ids = new item_int;
        real_other_ids->insert(real_other_ids->end(), other_set.begin(), other_set.end());
        if(down_blocks.size())
        {
            calc_down = false;
            this_tree.alloc_down();
            std::sort(down_blocks.begin(), down_blocks.end());
            other_center = down_blocks.back();
            down_blocks.pop_back(); // remove last one
            other_dispatch_ll.resize(real_other_ids->size());
            other_dispatch_rr.resize(real_other_ids->size());
            for (size_t i = 0; i < real_other_ids->size(); ++i)
            {
                other_dispatch_ll[i] = ll[mp_ids[real_other_ids -> at(i)]];
                other_dispatch_rr[i] = rr[mp_ids[real_other_ids -> at(i)]];
            }
#if _DEBUG
            rand_down = rand_gen_lock();
            fprintf(stderr, "Recurse: down decompose %d %u, dispatch to %u\n", depth, part, rand_down);
            spawn(new(child()) TaskRunnerSSW(decompose_Task(*real_other_ids, std::ref(other_center), std::ref(down_blocks), std::ref(other_dispatch_ll), std::ref(other_dispatch_rr), depth + 1, rand_down,  std::ref(cluster_mp_ids), *this_tree.down_tree, std::ref(lock), std::ref(file_down))));
#else
            spawn(new(child()) TaskRunnerSSW(decompose_Task(*real_other_ids, std::ref(other_center), std::ref(down_blocks), std::ref(other_dispatch_ll), std::ref(other_dispatch_rr), depth + 1, rand_gen_lock(), std::ref(cluster_mp_ids), *this_tree.down_tree, std::ref(lock), std::ref(file_down))));
#endif
        }
        else
        {
            this_tree.alloc_down();
            file_down = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_down.fasta";
            other_dispatch_ll.resize(real_other_ids->size());
            other_dispatch_rr.resize(real_other_ids->size());
            for (size_t i = 0; i < real_other_ids->size(); ++i)
            {
                other_dispatch_ll[i] = ll[mp_ids[real_other_ids -> at(i)]];
                other_dispatch_rr[i] = rr[mp_ids[real_other_ids -> at(i)]];
            }
        }
    }
    else calc_down = false;

    if(calc_left)
    {
#if _DEBUG
        fprintf(stderr, "Recurse: dispatch_left %d %u\n", depth, part);
#endif
        Frag_dispatch_with_lock((file_left).c_str(), *real_dispatch_lr, ll_dispatch_ll, ll_dispatch_rr, *this_tree.left_tree, lock);
    }
    if(calc_right)
    {
#if _DEBUG
        fprintf(stderr, "Recurse: dispatch_right %d %u\n", depth, part);
#endif
        Frag_dispatch_with_lock((file_right).c_str(), *real_dispatch_lr, rr_dispatch_ll, rr_dispatch_rr, *this_tree.right_tree, lock);
    }
    if(calc_down)
    {
#if _DEBUG
        fprintf(stderr, "Recurse: dispatch_down %d %u\n", depth, part);
#endif
        Frag_dispatch_with_lock((file_down).c_str(), *real_other_ids, other_dispatch_ll, other_dispatch_rr, *this_tree.down_tree, lock);
    }

    // merge left and right gaps
    wait();
    this_tree.merge_lr();
    if (use_center_tree)
    {
        this_tree.merge(center_tree);
    }
#if _DEBUG
    write_gaps_to_file_with_lr((tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_lr.fasta").c_str(),
        *real_dispatch_lr, ll_dispatch_ll, rr_dispatch_rr, ll_dispatch_rr, rr_dispatch_ll, this_tree);
    assert(aligned(*real_dispatch_lr, ll_dispatch_ll, rr_dispatch_rr, this_tree));
#endif

    refinement(*real_dispatch_lr, ll_dispatch_rr, this_tree);
    refinement(*real_dispatch_lr, rr_dispatch_ll, this_tree);

#if _DEBUG
    write_gaps_to_file_with_lr((tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_lr_refined.fasta").c_str(),
        *real_dispatch_lr, ll_dispatch_ll, rr_dispatch_rr, ll_dispatch_rr, rr_dispatch_ll, this_tree);
    assert(aligned(*real_dispatch_lr, ll_dispatch_ll, rr_dispatch_rr, this_tree));
#endif

    // merge other gaps
    if(other_set.size())
    {    
        this_tree.merge_down();
        // merge real_dispatch and real_other_ids
#if _DEBUG
        fprintf(stderr, "Recurse: dispatch_down merging %d %u\n", depth, part);
#endif
        TwoProfile_dispatch_with_lock(*real_dispatch_lr, *real_other_ids, ll, rr, mp_ids, this_tree, depth, part, lock);
        delete real_other_ids;
    }
    // calculate the length for adding 0-length sequences
    int_ sequence_length = 0, has_length_id = 0;
    while(ll[has_length_id] == rr[has_length_id]) ++ has_length_id;
    sequence_length = rr[has_length_id] - ll[has_length_id] + this_tree.this_tree[sequence_ids[has_length_id]].gap_sum();
    for(auto &zero_id: zero_set) this_tree.add_gaps(zero_id, mp_ids[sequence_ids[zero_id]], sequence_length);
    // calculate done
    final_name = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_final.fasta";
#if _DEBUG
    // Part 5: merge and print to file
    write_gaps_to_file_in_dispatch(final_name.c_str(), sequence_ids, ll, rr, this_tree); // may be useless in program, can be deprecated
    assert(aligned(sequence_ids, ll, rr, this_tree));
#endif
    // clean up pointer items
    if(real_dispatch_lr != &sequence_ids) delete real_dispatch_lr;

}

class TaskRunnerNoSSW: public staccato::task<TaskRunnerNoSSW>
{
public:
    TaskRunnerNoSSW(decompose_Task task_) : task_(task_) {}
    void execute();
private:
    decompose_Task task_;
};

/*
 * The decopmose function. This function will decompose seqeuences into *three* parts like this.
 * |||||||||||||||||||||||.||---------------------------||.|||||||||||||||||||||
 * |||||||-----------|||||.||---------------------------||.|||||-----------|||||
 * |||||||-*********-|||||.||---------------------------||.|||||-*********-|||||
 * |||||||-*********-|||||.||---------------------------||.|||||-*********-|||||
 * |||||||-*********-|||||.||---------------------------||.|||||-*********-|||||
 * |||||||-----------|||||.||---------------------------||.|||||-----------|||||
 * |||||||||||||||||||||||.||---------------------------||.|||||||||||||||||||||
 * .............................................................................
 * |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 * |||||||-----------|||||||||||||||||||||||||||||||||||||||||||-----------|||||
 * |||||||-*********-|||||||||||||||||||||||||||||||||||||||||||-|||||||||-|||||
 * |||||||-*********-|||||||||||||||||||||||||||||||||||||||||||-|||||||||-|||||
 * |||||||-----------|||||||||||||||||||||||||||||||||||||||||||-----------|||||
 */

void TaskRunnerNoSSW::execute()
{
    auto &[sequence_ids, center_block, other_blocks, ll, rr, depth, part, cluster_mp_ids, this_tree, lock, final_name] = task_;
    if(center_block.container.size() <= 0)
    {
        final_name = "";
        return;
    }
    // Part 0: sort the blocks
    std::sort(other_blocks.begin(), other_blocks.end(), std::greater<block_irreg>());
    // Part 1: left_blocks = [], right_blocks = [], down_blocks = []
    std::vector <block_irreg> left_blocks, right_blocks, down_blocks;
    item_int other_ids, * real_other_ids = nullptr, *real_dispatch_lr = nullptr;
    item_int ll_dispatch_ll, ll_dispatch_rr, rr_dispatch_ll, rr_dispatch_rr, other_dispatch_ll, other_dispatch_rr;
    set_int other_set, center_items, all_set, zero_set;
    center_int center_block_info;
    map_int mp_ids; // not in other blocks: real_id->id
    block_irreg* left_block, * right_block, * down_block;
    std::string file_left, file_right, file_down;
    gap_load_int gap_preload_in_main;
    int_ center_block_len = center_block.length, max_len = -1, min_len = INT_MAX;
    gap_tree center_tree;
    bool use_center_tree = false;
    // Part 2: check the number in center_block is same as the cluster
    for (size_t i = 0; i < sequence_ids.size(); ++i)
    {
        mp_ids[sequence_ids[i]] = i;
        if(ll[i] == rr[i]) zero_set.insert(sequence_ids[i]);
        else all_set.insert(sequence_ids[i]);
    }
    for (auto &item: center_block.container)
    {
        auto &[real_id, block_ll, delta] = item;
        center_block_info[real_id] = std::make_pair(block_ll, block_ll + center_block_len - delta);
        max_len = (std::max)(center_block_len - delta, max_len);
        min_len = (std::min)(center_block_len - delta, min_len);
    }
    if(max_len <= 0)
    {
        fprintf(stderr, "Warning: this block is invalid. Program will return.\n");
        std::string final_file_name = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_final.fasta";
        Frag_dispatch_with_lock(final_file_name.c_str(), sequence_ids, ll, rr, this_tree, lock);
        final_name = "";
        return ; // can not find invalid block
    }
    if (min_len != max_len)
    {
        // block is irregular, align it!
        fprintf(stderr, "Info: this block is not same (although contains all sequences). Program will align this block.\n");
        std::string center_file_name = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_block.fasta";
        use_center_tree = true;
        Center_block_dispatch_with_lock(center_file_name.c_str(), center_block_info, center_tree, lock);
    }
    if (sequence_ids.size() != center_block_info.size())
    {
        assert(sequence_ids.size() > center_block_info.size());
        for (auto& item : center_block_info) center_items.insert(item.first);
        std::set_difference(all_set.begin(), all_set.end(), center_items.begin(), center_items.end(), std::inserter(other_set, other_set.begin()));
        all_set.clear();
        // NO SSW!
        // Part 2.3: determine how many sequences will be in dispatch mode
        real_dispatch_lr = new item_int;
        // Part 2.4: determine real_dispatch_lr = center_block_info.first
        real_dispatch_lr->reserve(center_block_info.size());
        for(auto &item: center_block_info) real_dispatch_lr->emplace_back(item.first);
    }
    else
    {
        real_dispatch_lr = new item_int;
        for (auto& item : all_set) real_dispatch_lr->emplace_back(item);
        all_set.clear();
    }
    // Part 3: determine other_blocks is in left or right or down, or just simply, drop it
    for (auto& block_item : other_blocks)
    {
        left_block = new block_irreg;
        right_block = new block_irreg;
        if (other_set.size()) down_block = new block_irreg;
        else down_block = nullptr;
        // check the place for every block
        for (auto& sequence_item : block_item.container)
        {
            auto &[this_id, this_ll, delta] = sequence_item;
            auto this_rr = this_ll + block_item.length - delta;
            auto mp_place = mp_ids[this_id];
            auto &border_ll = ll[mp_place], &border_rr = rr[mp_place];
            if(zero_set.find(this_id) != zero_set.end()) continue;
            if(down_block && other_set.find(this_id) != other_set.end()) // will be in other block
            {
                if(border_ll <= this_ll && this_rr <= border_rr && this_ll + threshold < this_rr)
                {
#if _DEBUG_BLOCK
                    fprintf(stderr, "Border = [%d, %d); ", border_ll, border_rr);
                    fprintf(stderr, "dispatch to other = [%d, %d)\n", this_ll, this_rr);
#endif
                    down_block->container.emplace_back(sequence_item);
                }
            }
            else
            {
                auto& [centerll, centerrr] = center_block_info[this_id];
                auto ll_border_ll = border_ll, ll_border_rr = centerll, rr_border_ll = centerrr, rr_border_rr = border_rr;
#if _DEBUG_BLOCK
                fprintf(stderr, "Center block = [%d, %d), will dispatch block = [%d, %d)\n", centerll, centerrr, this_ll, this_rr);
                fprintf(stderr, "Border = [%d, %d); left border = [%d, %d); right border = [%d, %d)\n",
                        border_ll, border_rr, ll_border_ll, ll_border_rr, rr_border_ll, rr_border_rr);
#endif
                if (centerll <= this_ll && this_rr <= centerrr) continue; // dropped
                else if (this_ll <= centerll && centerrr <= this_rr) // across the block
                {
#if _DEBUG_BLOCK
                    fprintf(stderr,
                            "Warning: sequence %d has condition:\ncenter block = [%d, %d), this block = [%d, %d)\n",
                            this_id, centerll, centerrr, this_ll, this_rr);
#endif
                    sequence_irreg sequence_item_backup;
                    // insert it to left block
                    if(ll_border_ll <= this_ll && this_ll + block_item.length - (this_rr - centerll) <= ll_border_rr &&
                       this_ll + threshold < this_ll + block_item.length - (this_rr - centerll))
                    {
                        sequence_item_backup = sequence_item;
                        std::get<2>(sequence_item_backup) += this_rr - centerll;
                        left_block->container.emplace_back(sequence_item_backup);
#if _DEBUG_BLOCK
                        fprintf(stderr, "Left = [%d, %d);\n", this_ll, this_ll + block_item.length - (this_rr - centerll));
#endif
                    }
                    // insert it to right block
                    if(rr_border_ll <= centerrr && center_block.length + this_ll <= rr_border_rr &&
                       centerrr + threshold < center_block.length + this_ll)
                    {
                        sequence_item_backup = sequence_item;
                        std::get<1>(sequence_item_backup) = centerrr;
                        std::get<2>(sequence_item_backup) += (centerrr - this_ll);
                        right_block->container.emplace_back(sequence_item_backup);
#if _DEBUG_BLOCK
                        fprintf(stderr, "Right = [%d, %d)\n", centerrr, centerrr + block_item.length - (centerrr - this_ll));
#endif
                    }
                }
                else if (this_ll <= centerll && centerll <= this_rr && this_rr <= centerrr)
                    // will be in left block, and end place is centerll
                {
                    if(ll_border_ll <= this_ll && centerll <= ll_border_rr && this_ll + threshold < centerll)
                    {
                        delta += (this_rr - centerll);
                        left_block->container.emplace_back(sequence_item);
#if _DEBUG_BLOCK
                        fprintf(stderr, "Border = [%d, %d); make changes ", ll_border_ll, ll_border_rr);
                        fprintf(stderr, "dispatch to left = [%d, %d); ", this_ll, this_ll + block_item.length - delta);
                        fprintf(stderr, "check = %d\n", centerll == this_ll + block_item.length - delta);
#endif
                    }
                }
                else if (centerll <= this_ll && this_ll <= centerrr && centerrr <= this_rr)
                    // will be in right block, and start place is centerrr
                {
                    if(rr_border_ll <= centerrr && this_rr <= rr_border_rr && centerrr + threshold < this_rr)
                    {
                        delta += (centerrr - this_ll);
                        this_ll = centerrr;
                        right_block->container.emplace_back(sequence_item);
#if _DEBUG_BLOCK
                        fprintf(stderr, "Border = [%d, %d); make changes ", rr_border_ll, rr_border_rr);
                        fprintf(stderr, "dispatch to right = [%d, %d); ", this_ll, this_ll + block_item.length - delta);
                        fprintf(stderr, "check = %d\n", this_rr == this_ll + block_item.length - delta);
#endif
                    }
                }
                else if (this_rr <= centerll) // will be in left block
                {
                    if(ll_border_ll <= this_ll && this_rr <= ll_border_rr)
                    {
                        left_block->container.emplace_back(sequence_item);
#if _DEBUG_BLOCK
                        fprintf(stderr, "Border = [%d, %d); ", ll_border_ll, ll_border_rr);
                        fprintf(stderr, "dispatch to left = [%d, %d)\n", this_ll, this_rr);
#endif
                    }
                }
                else if (centerrr <= this_ll) // will be in right block
                {
                    if(rr_border_ll <= this_ll && this_rr <= rr_border_rr)
                    {
                        right_block->container.emplace_back(sequence_item);
#if _DEBUG_BLOCK
                        fprintf(stderr, "Border = [%d, %d); ", rr_border_ll, rr_border_rr);
                        fprintf(stderr, "dispatch to right = [%d, %d)\n", this_ll, this_rr);
#endif
                    }
                }
                else
                {
                    fprintf(stderr,
                           "Designer Error: not considered condition: center = [%d, %d), this = [%d, %d)\n",
                            centerll, centerrr, this_ll, this_rr);
                    assert(false);
                }
            }
        }
        // insert left block to left blocks
        if (left_block->container.size() >= 2)
        {
            left_block->length = block_item.length;
            left_block->center = find_center(left_block->container);
            left_block->calc_area();
            left_blocks.emplace_back(std::move(*left_block));
        }
        else
        {
            left_block->freecontainer();
            delete left_block;
        }
        // insert right block to right blocks
        if (right_block->container.size() >= 2)
        {
            right_block->length = block_item.length;
            right_block->center = find_center(right_block->container);
            right_block->calc_area();
            right_blocks.emplace_back(std::move(*right_block));
        }
        else
        {
            right_block->freecontainer();
            delete right_block;
        }
        // insert down block to down blocks
        if (other_set.size())
        {
            if (down_block->container.size() >= 2)
            {
                down_block->length = block_item.length;
                down_block->center = find_center(down_block->container);
                down_block->calc_area();
                down_blocks.emplace_back(std::move(*down_block));
            }
            else
            {
                down_block->freecontainer();
                delete down_block;
            }
        }
    }
    // Part 4: decompose three parts: left_blocks, right_blocks, other_blocks
    this_tree.alloc_lr();
    // Please sort id and center block
    // Part 4.1: if can not be decomposed, call fragalign for making alignment
    // determine the dispatch_ll and dispatch_rr, and call for sub process
    block_irreg left_center, right_center, other_center;
    unsigned int rand_left, rand_right, rand_down;
    bool calc_left = true, calc_right = true, calc_down = true;
    if (left_blocks.size())
    {
        calc_left = false;
        std::sort(left_blocks.begin(), left_blocks.end());
        left_center = left_blocks.back();
        left_blocks.pop_back(); // remove last one
        ll_dispatch_ll.resize(real_dispatch_lr->size());
        ll_dispatch_rr.resize(real_dispatch_lr->size());
        for (size_t i = 0; i < real_dispatch_lr->size(); ++i)
        {
            ll_dispatch_rr[i] = center_block_info[real_dispatch_lr->at(i)].first;
            ll_dispatch_ll[i] = ll[mp_ids[real_dispatch_lr->at(i)]];
        }
        // dispatch to left next-level
#if _DEBUG
        rand_left = rand_gen_lock();
        fprintf(stderr, "Recurse: left decompose %d %u, dispatch to %u\n", depth, part, rand_left);
        spawn(new(child()) TaskRunnerNoSSW(decompose_Task(*real_dispatch_lr, std::ref(left_center), std::ref(left_blocks), std::ref(ll_dispatch_ll), std::ref(ll_dispatch_rr), depth + 1, rand_left,  std::ref(cluster_mp_ids), *this_tree.left_tree, std::ref(lock), std::ref(file_left))));
#else
        spawn(new(child()) TaskRunnerNoSSW(decompose_Task(*real_dispatch_lr, std::ref(left_center), std::ref(left_blocks), std::ref(ll_dispatch_ll), std::ref(ll_dispatch_rr), depth + 1, rand_gen_lock(), std::ref(cluster_mp_ids), *this_tree.left_tree, std::ref(lock), std::ref(file_left))));
#endif
    }
    else
    {
        file_left = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_left.fasta";
        ll_dispatch_ll.resize(real_dispatch_lr->size());
        ll_dispatch_rr.resize(real_dispatch_lr->size());
        for (size_t i = 0; i < real_dispatch_lr->size(); ++i)
        {
            ll_dispatch_ll[i] = ll[mp_ids[real_dispatch_lr->at(i)]];
            ll_dispatch_rr[i] = center_block_info[real_dispatch_lr->at(i)].first;
        }
    }

    if (right_blocks.size())
    {
        calc_right = false;
        std::sort(right_blocks.begin(), right_blocks.end());
        right_center = right_blocks.back();
        right_blocks.pop_back(); // remove last one
        rr_dispatch_ll.resize(real_dispatch_lr->size());
        rr_dispatch_rr.resize(real_dispatch_lr->size());
        for (size_t i = 0; i < real_dispatch_lr->size(); ++i)
        {
            rr_dispatch_ll[i] = center_block_info[real_dispatch_lr->at(i)].second;
            rr_dispatch_rr[i] = rr[mp_ids[real_dispatch_lr->at(i)]];
        }
#if _DEBUG
        rand_right = rand_gen_lock();
        fprintf(stderr, "Recurse: right decompose %d %u, dispatch to %u\n", depth, part, rand_right);
        spawn(new(child()) TaskRunnerNoSSW(decompose_Task(*real_dispatch_lr, std::ref(right_center), std::ref(right_blocks), std::ref(rr_dispatch_ll), std::ref(rr_dispatch_rr), depth + 1, rand_right, std::ref(cluster_mp_ids), *this_tree.right_tree, std::ref(lock), std::ref(file_right))));
#else
        spawn(new(child()) TaskRunnerNoSSW(decompose_Task(*real_dispatch_lr, std::ref(right_center), std::ref(right_blocks), std::ref(rr_dispatch_ll), std::ref(rr_dispatch_rr), depth + 1, rand_gen_lock(), std::ref(cluster_mp_ids), *this_tree.right_tree, std::ref(lock), std::ref(file_right))));
#endif
    }
    else
    {
        file_right = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_right.fasta";
        rr_dispatch_ll.resize(real_dispatch_lr->size());
        rr_dispatch_rr.resize(real_dispatch_lr->size());
        for (size_t i = 0; i < real_dispatch_lr->size(); ++i)
        {
            rr_dispatch_ll[i] = center_block_info[real_dispatch_lr->at(i)].second;
            rr_dispatch_rr[i] = rr[mp_ids[real_dispatch_lr->at(i)]];
        }
    }

    if (other_set.size())
    {
        real_other_ids = new item_int;
        real_other_ids->insert(real_other_ids->end(), other_set.begin(), other_set.end());
        if(down_blocks.size())
        {
            calc_down = false;
            this_tree.alloc_down();
            std::sort(down_blocks.begin(), down_blocks.end());
            other_center = down_blocks.back();
            down_blocks.pop_back(); // remove last one
            other_dispatch_ll.resize(real_other_ids->size());
            other_dispatch_rr.resize(real_other_ids->size());
            for (size_t i = 0; i < real_other_ids->size(); ++i)
            {
                other_dispatch_ll[i] = ll[mp_ids[real_other_ids -> at(i)]];
                other_dispatch_rr[i] = rr[mp_ids[real_other_ids -> at(i)]];
            }
#if _DEBUG
            rand_down = rand_gen_lock();
            fprintf(stderr, "Recurse: down decompose %d %u, dispatch to %u\n", depth, part, rand_down);
            spawn(new(child()) TaskRunnerNoSSW(decompose_Task(*real_other_ids, std::ref(other_center), std::ref(down_blocks), std::ref(other_dispatch_ll), std::ref(other_dispatch_rr), depth + 1, rand_down,  std::ref(cluster_mp_ids), *this_tree.down_tree, std::ref(lock), std::ref(file_down))));
#else
            spawn(new(child()) TaskRunnerNoSSW(decompose_Task(*real_other_ids, std::ref(other_center), std::ref(down_blocks), std::ref(other_dispatch_ll), std::ref(other_dispatch_rr), depth + 1, rand_gen_lock(), std::ref(cluster_mp_ids), *this_tree.down_tree, std::ref(lock), std::ref(file_down))));
#endif
        }
        else
        {
            this_tree.alloc_down();
            file_down = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_down.fasta";
            other_dispatch_ll.resize(real_other_ids->size());
            other_dispatch_rr.resize(real_other_ids->size());
            for (size_t i = 0; i < real_other_ids->size(); ++i)
            {
                other_dispatch_ll[i] = ll[mp_ids[real_other_ids -> at(i)]];
                other_dispatch_rr[i] = rr[mp_ids[real_other_ids -> at(i)]];
            }
        }
    }
    else calc_down = false;

    if(calc_left)
    {
#if _DEBUG
        fprintf(stderr, "Recurse: dispatch_left %d %u\n", depth, part);
#endif
        Frag_dispatch_with_lock((file_left).c_str(), *real_dispatch_lr, ll_dispatch_ll, ll_dispatch_rr, *this_tree.left_tree, lock);
    }
    if(calc_right)
    {
#if _DEBUG
        fprintf(stderr, "Recurse: dispatch_right %d %u\n", depth, part);
#endif
        Frag_dispatch_with_lock((file_right).c_str(), *real_dispatch_lr, rr_dispatch_ll, rr_dispatch_rr, *this_tree.right_tree, lock);
    }
    if(calc_down)
    {
#if _DEBUG
        fprintf(stderr, "Recurse: dispatch_down %d %u\n", depth, part);
#endif
        Frag_dispatch_with_lock((file_down).c_str(), *real_other_ids, other_dispatch_ll, other_dispatch_rr, *this_tree.down_tree, lock);
    }

    // merge left and right gaps
    wait();
    this_tree.merge_lr();
    if (use_center_tree)
    {
        this_tree.merge(center_tree);
    }
#if _DEBUG
    write_gaps_to_file_with_lr((tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_lr.fasta").c_str(),
        *real_dispatch_lr, ll_dispatch_ll, rr_dispatch_rr, ll_dispatch_rr, rr_dispatch_ll, this_tree);
    assert(aligned(*real_dispatch_lr, ll_dispatch_ll, rr_dispatch_rr, this_tree));
#endif

    refinement(*real_dispatch_lr, ll_dispatch_rr, this_tree);
    refinement(*real_dispatch_lr, rr_dispatch_ll, this_tree);

#if _DEBUG
    write_gaps_to_file_with_lr((tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_lr_refined.fasta").c_str(),
        *real_dispatch_lr, ll_dispatch_ll, rr_dispatch_rr, ll_dispatch_rr, rr_dispatch_ll, this_tree);
    assert(aligned(*real_dispatch_lr, ll_dispatch_ll, rr_dispatch_rr, this_tree));
#endif

    // merge other gaps
    if(other_set.size())
    {    
        this_tree.merge_down();
        // merge real_dispatch and real_other_ids
#if _DEBUG
        fprintf(stderr, "Recurse: dispatch_down merging %d %u\n", depth, part);
#endif
        TwoProfile_dispatch_with_lock(*real_dispatch_lr, *real_other_ids, ll, rr, mp_ids, this_tree, depth, part, lock);
        delete real_other_ids;
    }
    // calculate the length for adding 0-length sequences
    int_ sequence_length = 0, has_length_id = 0;
    while(ll[has_length_id] == rr[has_length_id]) ++ has_length_id;
    sequence_length = rr[has_length_id] - ll[has_length_id] + this_tree.this_tree[sequence_ids[has_length_id]].gap_sum();
    for(auto &zero_id: zero_set) this_tree.add_gaps(zero_id, mp_ids[sequence_ids[zero_id]], sequence_length);
    // calculate done

    final_name = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_final.fasta";
#if _DEBUG
    // Part 5: merge and print to file
    write_gaps_to_file_in_dispatch(final_name.c_str(), sequence_ids, ll, rr, this_tree); // may be useless in program, can be deprecated
    assert(aligned(sequence_ids, ll, rr, this_tree));
#endif
    // clean up pointer items
    if(real_dispatch_lr != &sequence_ids) delete real_dispatch_lr;
}

std::string decompose_with_SSW(item_int& sequence_ids, block_irreg& center_block, std::vector <block_irreg>& other_blocks, item_int& ll, item_int& rr, map_int &cluster_mp_ids, gap_tree &this_tree, int threads)
{
    staccato::scheduler<TaskRunnerSSW> ____(3, threads);
    std::string final_file_name;
    std::mutex lock;
    ____.spawn(new(____.root()) TaskRunnerSSW(decompose_Task(sequence_ids, center_block, other_blocks, ll, rr, 0, 0u, cluster_mp_ids, std::ref(this_tree), &lock, final_file_name)));
    ____.wait();
    return final_file_name;
}

std::string decompose_without_SSW(item_int& sequence_ids, block_irreg& center_block, std::vector <block_irreg>& other_blocks, item_int& ll, item_int& rr, map_int &cluster_mp_ids, gap_tree &this_tree, int threads)
{
    staccato::scheduler<TaskRunnerNoSSW> ____(3, threads);
    std::string final_file_name;
    std::mutex lock;
    ____.spawn(new(____.root()) TaskRunnerNoSSW(decompose_Task(sequence_ids, center_block, other_blocks, ll, rr, 0, 0u, cluster_mp_ids, std::ref(this_tree), &lock, final_file_name)));
    ____.wait();
    return final_file_name;

}
