#include "libblockaligner.hpp"
#include "irreg_part.hpp"
//#define _OPENMP 1
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
extern item_int first_end;
extern std::mt19937 rand_gen;
extern int_ threshold;

#pragma GCC push_options
#pragma GCC optimize ("O0")
void WriteSeq_SSW_Frag(FILE* f, const int_& ll, const int_& rr, const int_& real_id)
{
#if _DEBUG
    fprintf(f, "> %d %s %d %d\n", real_id, comment[real_id], ll, rr);
#else
    fprintf(f, "> %d %s\n", real_id, comment[real_id]);
#endif
    for (int_ i = ll; i < rr; ++i)
    {
#if 0
        auto this_gap = gap_stores[cluster_id].get_gaps(i - 1);
        for (int_ gap = 0; gap < this_gap; ++gap) fputc('-', f);
#endif
        fputc(seq[real_id][i], f);
    }
    fputc('\n', f);
}

int Frag_dispatch(const char* File, item_int& sequence_id, item_int& ll, item_int& rr, gap_tree &this_tree)
{
    FILE* f = fopen(File, "w");
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
    if (Fragalignment(File, outputfile.c_str())) { fprintf(stderr, "Error: can not execute alignment. Program will exit.\n"); exit(1); }
#if _DEBUG
    if (Filedelete(File) || Filemove(outputfile.c_str(), File)) { fprintf(stderr, "Error: can not move file to source name. Program will exit.\n"); exit(1); }
    if (read_gaps_from_file(File, sequence_id, ll, rr, this_tree)) { fprintf(stderr, "Error: can not read sequence info. Program will exit.\n"); exit(1); }
#else
    if (read_gaps_from_file(outputfile.c_str(), sequence_id, ll, rr, this_tree)) { fprintf(stderr, "Error: can not read sequence info. Program will exit.\n"); exit(1); }
    if (Filedelete(File) || Filedelete(outputfile.c_str())) { fprintf(stderr, "Warning: can not remove alignment file.\n"); }
#endif
    return 0;
}

int Center_block_dispatch(const char* File, center_int& center_info, gap_tree& this_tree)
{
    FILE* f = fopen(File, "w");
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
    if (Fragalignment(File, outputfile.c_str())) { fprintf(stderr, "Error: can not execute alignment. Program will exit.\n"); exit(1); }
#if _DEBUG
    if (Filedelete(File) || Filemove(outputfile.c_str(), File)) { fprintf(stderr, "Error: can not move file to source name. Program will exit.\n"); exit(1); }
    if (read_gaps_from_file(File, center_info, this_tree)) { fprintf(stderr, "Error: can not read sequence info. Program will exit.\n"); exit(1); }
#else
    if (read_gaps_from_file(outputfile.c_str(), center_info, this_tree)) { fprintf(stderr, "Error: can not read sequence info. Program will exit.\n"); exit(1); }
    if (Filedelete(File) || Filedelete(outputfile.c_str())) { fprintf(stderr, "Warning: can not remove alignment file.\n"); }
#endif
    return 0;
}


int Frag_dispatch_full_seq(const char* File, const item_int& seq_id)
{
    FILE* f = fopen(File, "w");
    if (f == NULL) { fprintf(stderr, "Error: can not open file %s. Program will exit.\n", File); exit(1); }
    for(auto &id: seq_id)
    {
        fprintf(f, "> %s\n", comment[id]);
        fputs(seq[id], f);
        fputc('\n', f);
    }
    fclose(f);
    std::string outputfile;
    outputfile += File;
    outputfile += ".out";
    if (Fragalignment(File, outputfile.c_str())) { fprintf(stderr, "Error: can not execute alignment. Program will exit.\n"); exit(1); }
    if (Filedelete(File) || Filemove(outputfile.c_str(), File)) { fprintf(stderr, "Error: can not move file to source name. Program will exit.\n"); exit(1); }
    return 0;
}

void TwoProfile_writeseq(const char* name, item_int &ids, item_int &ll, item_int &rr, map_int &mp, gap_tree &this_tree)
{
    int_ this_gap;
    FILE *f = fopen(name, "w");
    if(f == NULL) { fprintf(stderr, "Error: can not make temp file %s. Program will exit.\n", name); exit(1); }
    for(auto &sequence_id: ids)
    {
        auto &real_place = mp[sequence_id];
        auto &lll = ll[real_place], &rrr = rr[real_place];
        write_seq(sequence_id, lll, rrr, this_tree.get_gaps(sequence_id), f);
    }
    fclose(f);
}

int TwoProfile_dispatch(item_int &dispatch_ids, item_int &other_ids, item_int &ll, item_int &rr, map_int &mp, 
                        gap_tree &this_tree, const int depth, const unsigned int part)
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
    if(TwoProfileAlignment(File1.c_str(), File2.c_str(), File3.c_str())) { fprintf(stderr, "Error: can not two profile alignment. Program will exit.\n"); exit(1); };

    // read result from two profiles
    if (read_gaps_from_merged_profile(File3.c_str(), dispatch_ids, other_ids, ll, rr, mp, this_tree))
    {
        fprintf(stderr, "Error: can not read sequence info. Program will exit.\n");
        exit(1);
    }

    Filedelete(File1.c_str());
    Filedelete(File2.c_str());
    Filedelete(File3.c_str());
    return 0;
}

int TwoProfile_dispatch(const std::string& File1, const std::string& File2)
{
    // call two_profile_align
    std::string File3 = File1 + ".tmp";
    TwoProfileAlignment(File1.c_str(), File2.c_str(), File3.c_str());
    Filemove(File3.c_str(), File1.c_str());
    return 0;
}

int Profiles_dispatch(std::vector<item_int> &files, std::vector<std::string> &file_names, const char* out_name, std::string &center_file_name, std::string &file_name)
{
    int_ cluster_id, max_id, max_len;
    FILE *f_center = fopen(center_file_name.c_str(), "wb");
    if(f_center == NULL) { fprintf(stderr, "Error: can not write center file. Program will exit.\n"); exit(1); }
    for(auto &item: files)
    {
        // choose the longest sequence in profile-profile alignment
        max_id = item[0];
        max_len = first_end[max_id];
        for(auto &seq: item)
        {
            if(first_end[seq] > max_len)
            {
                max_id = seq;
                max_len = first_end[seq];
            }
        }
        fprintf(f_center, "> %d %s\n", max_id, comment[max_id]);
        for(size_t i = 0; seq[max_id][i]; ++ i) fputc(seq[max_id][i], f_center);
        fputc('\n', f_center);
    }
    fclose(f_center);

    FILE *f_list = fopen(file_name.c_str(), "wb");
    if(f_list == NULL) { fprintf(stderr, "Error: can not write list file. Program will exit.\n"); exit(1); }
    for(auto &seq_name: file_names)
    {
        fputs(seq_name.c_str(), f_list);
        fputc('\n', f_list);
    }
    fclose(f_list);
    extern aligner_config align_profiles_config_;
    if(align_profiles_config_.profile2_prefix) return WMSAProfileProfilealignment(center_file_name.c_str(), file_name.c_str(), out_name);
    else return MAFFTProfileProfilealignment(file_name.c_str(), out_name);
}
#pragma GCC pop_options

#if _DEBUG
bool aligned(item_int &seq_ids, item_int &ll, item_int &rr, gap_tree &gap_tree_main)
{
    if(seq_ids.size() < 1) return true;
    int_ seq_len = rr[0] - ll[0] + gap_tree_main.this_tree[seq_ids[0]].gap_sum();
    for(int_ i = 1, this_seq_len; i < seq_ids.size(); ++ i)
    {
        this_seq_len = rr[i] - ll[i] + gap_tree_main.this_tree[seq_ids[i]].gap_sum();
        if(this_seq_len != seq_len)
        {
            fprintf(stderr, "Error: not aligned in sequence %d (in [%d, %d), inserted %d gaps)\n", seq_ids[i], ll[i], rr[i], gap_tree_main.this_tree[seq_ids[i]].gap_sum());
            fprintf(stderr, "- Sequence %d: %d\n- Sequence %d: %d\n", seq_ids[0], seq_len, seq_ids[i], this_seq_len);
            return false;
        }
    }
    return true;
}
#endif

/*
 The find similarity function.
 * Use ssw program for finding the items in seq_list called need_align_seq_id, and made need_align_seq_id-th sequence segment
   in [need_align_center - (template_seq_len + ext_size), need_align_center + (template_seq_len + ext_size)]
 * modifing args:
 - real_other_ids: a set, contains failed to find ids;
 - center_block:   a map, first is the id of successfully to find, second is tuple, contains start place and end place
 */
int8_t blosum50[] = {
    //  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
        5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -2, -1, -1, -5,	// A
       -2,  7, -1, -2, -4,  1,  0, -3,  0, -4, -3,  3, -2, -3, -3, -1, -1, -3, -1, -3, -1,  0, -1, -5,	// R
       -1, -1,  7,  2, -2,  0,  0,  0,  1, -3, -4,  0, -2, -4, -2,  1,  0, -4, -2, -3,  5,  0, -1, -5,	// N
       -2, -2,  2,  8, -4,  0,  2, -1, -1, -4, -4, -1, -4, -5, -1,  0, -1, -5, -3, -4,  6,  1, -1, -5,	// D
       -1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -3, -3, -1, -5,	// C
       -1,  1,  0,  0, -3,  7,  2, -2,  1, -3, -2,  2,  0, -4, -1,  0, -1, -1, -1, -3,  0,  4, -1, -5,	// Q
       -1,  0,  0,  2, -3,  2,  6, -3,  0, -4, -3,  1, -2, -3, -1, -1, -1, -3, -2, -3,  1,  5, -1, -5,	// E
        0, -3,  0, -1, -3, -2, -3,  8, -2, -4, -4, -2, -3, -4, -2,  0, -2, -3, -3, -4, -1, -2, -1, -5,	// G
       -2,  0,  1, -1, -3,  1,  0, -2, 10, -4, -3,  0, -1, -1, -2, -1, -2, -3,  2, -4,  0,  0, -1, -5,	// H
       -1, -4, -3, -4, -2, -3, -4, -4, -4,  5,  2, -3,  2,  0, -3, -3, -1, -3, -1,  4, -4, -3, -1, -5,	// I
       -2, -3, -4, -4, -2, -2, -3, -4, -3,  2,  5, -3,  3,  1, -4, -3, -1, -2, -1,  1, -4, -3, -1, -5,	// L
       -1,  3,  0, -1, -3,  2,  1, -2,  0, -3, -3,  6, -2, -4, -1,  0, -1, -3, -2, -3,  0,  1, -1, -5,	// K
       -1, -2, -2, -4, -2,  0, -2, -3, -1,  2,  3, -2,  7,  0, -3, -2, -1, -1,  0,  1, -3, -1, -1, -5,	// M
       -3, -3, -4, -5, -2, -4, -3, -4, -1,  0,  1, -4,  0,  8, -4, -3, -2,  1,  4, -1, -4, -4, -1, -5,	// F
       -1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3, -2, -1, -1, -5,	// P
        1, -1,  1,  0, -1,  0, -1,  0, -1, -3, -3,  0, -2, -3, -1,  5,  2, -4, -2, -2,  0,  0, -1, -5,	// S
        0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  2,  5, -3, -2,  0,  0, -1, -1, -5, 	// T
       -3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1,  1, -4, -4, -3, 15,  2, -3, -5, -2, -1, -5, 	// W
       -2, -1, -2, -3, -3, -1, -2, -3,  2, -1, -1, -2,  0,  4, -3, -2, -2,  2,  8, -1, -3, -2, -1, -5, 	// Y
        0, -3, -3, -4, -1, -3, -3, -4, -4,  4,  1, -3,  1, -1, -3, -2,  0, -3, -1,  5, -3, -3, -1, -5, 	// V
       -2, -1,  5,  6, -3,  0,  1, -1,  0, -4, -4,  0, -3, -4, -2,  0,  0, -5, -3, -3,  6,  1, -1, -5, 	// B
       -1,  0,  0,  1, -3,  4,  5, -2,  0, -3, -3,  1, -1, -4, -1,  0, -1, -2, -2, -3,  1,  5, -1, -5, 	// Z
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5, 	// X
       -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1 	// *
};

int8_t blosum62[] = {
    //  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
        4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4,  // A
       -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4,  // R
       -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4,  // N
       -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4,  // D
        0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4,  // C
       -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4,  // Q
       -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4,  // E
        0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4,  // G
       -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4,  // H
       -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4,  // I
       -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4,  // L
       -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4,  // K
       -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4,  // M
       -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4,  // F
       -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4,  // P
        1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4,  // S
        0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4,  // T
       -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4,  // W
       -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4,  // Y
        0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4,  // V
       -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4,  // B
       -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4,  // Z
        0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4,  // X
       -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1   // *
};

int8_t mat_nt[] = {
        2, -2, -2, -2, 0,
       -2,  2, -2, -2, 0,
       -2, -2,  2, -2, 0,
       -2, -2, -2,  2, 0,
        0,  0,  0,  0, 0
};


/* This table is used to transform amino acid letters into numbers. */
int8_t aa_table[128] = {
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
    14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
    23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
    14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23
};

/* This table is used to transform nucleotide letters into numbers. */
int8_t nt_table[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

void Findsim_and_extension(item_int& seq_list, set_int& not_same_len_list, set_int& real_ids, int_ need_align_center,
    item_int& all_ll, item_int& all_rr, map_int& ll_rr_map,
    int_ template_seq_id, int_ template_seq_begin, int_ template_seq_end, center_int &center_block, set_int &real_other_ids,
    gap_tree &gap_tree_main, int depth, unsigned int part)
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
#if _DEBUG
            fprintf(stderr, "Warning: sequence %d has no length\n", *real_id);
#endif
#ifdef _OPENMP
#pragma omp critical
#endif
            real_other_ids.insert(*real_id);
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
        seq_file_name = depth_part; seq_file_name += "SSW_Align.fasta";
        f_seq = fopen(seq_file_name.c_str(), "w");
        if(f_seq == NULL) { fprintf(stderr, "Can not open file %s (common file). Program will exit.\n", seq_file_name.c_str()); exit(1); }
        // write center sequence line
        WriteSeq_SSW_Frag(f_seq, template_seq_begin, template_seq_end, template_seq_id);

        // write common sequences and ids for output
        item_int vll, vrr, seq_frags_ids;
        vll.reserve(seq_list.size() + 1); vrr.reserve(seq_list.size() + 1); seq_frags_ids.reserve(seq_list.size() + 1);
        seq_frags_ids.emplace_back(template_seq_id);
        vll.emplace_back(template_seq_begin); vrr.emplace_back(template_seq_end);
        for (auto& real_id : seq_list)
        {
            if (real_other_ids.find(real_id) != real_other_ids.end()) continue;
            auto& [ll, rr] = center_block[real_id];
            vll.emplace_back(ll); vrr.emplace_back(rr); seq_frags_ids.emplace_back(real_id);
            WriteSeq_SSW_Frag(f_seq, ll, rr, real_id);
        }
        // add sequences which has not same block length
        for (auto &real_id: not_same_len_list)
        {
            auto& [ll, rr] = center_block[real_id];
            vll.emplace_back(ll); vrr.emplace_back(rr); seq_frags_ids.emplace_back(real_id);
            WriteSeq_SSW_Frag(f_seq, ll, rr, real_id);
        }
        fclose(f_seq);

        int_ template_seq_len = template_seq_end - template_seq_begin;
        out_file_name = depth_part + ".out.fasta";
        Frag_dispatch(out_file_name.c_str(), seq_frags_ids, vll, vrr, gap_tree_main);

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
#if 0
        out_file_name = depth_part + ".out_SSW.fasta";
        write_gaps_to_file(out_file_name.c_str(), center_block, true, gap_tree_main);
#endif
    }
    else 
    {
        if(! foundsw) fprintf(stderr, "Warning: not found vaild part for alignment.\n");
    }
}

void refinement(item_int &seq_ids, item_int& border, gap_tree& this_tree)
{
    int_ gaps = this_tree.get_gaps(seq_ids[0], border[0] - 1);
    if (gaps == 0)
    {
        fprintf(stderr, "Warning: sequence %d in place %d has 0 gaps. No need to refine.\n", seq_ids[0], border[0] - 1);
        return;
    }
    for (int i = 1; i < seq_ids.size(); ++i)
    {
        gaps = (std::min)(gaps, this_tree.get_gaps(seq_ids[i], border[i] - 1));
        if (gaps == 0)
        {
            fprintf(stderr, "Warning: sequence %d in place %d has 0 gaps. No need to refine.\n", seq_ids[i], border[i] - 1);
            return;
        }
    }
    for (int i = 0; i < seq_ids.size(); ++i)
    {
#if _DEBUG
        int_ now_ = this_tree.get_gaps(seq_ids[i], border[i] - 1);
#endif
        this_tree.add_gaps(seq_ids[i], border[i] - 1, -gaps);
#if _DEBUG
        assert(now_ == this_tree.get_gaps(seq_ids[i], border[i] - 1) + gaps);
#endif
    }
    fprintf(stderr, "Info: refined %d gaps.\n", gaps);
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
 * Return the merged file name. (may not need to return?)
 */
std::string decompose_without_SSW(item_int& sequence_ids,
                                  block_irreg& center_block,
                                  std::vector <block_irreg>& other_blocks,
                                  item_int& ll, item_int& rr,
                                  const int depth, const unsigned int part, map_int &cluster_mp_ids, gap_tree &this_tree)
{
#if _DEBUG
    write_raw_to_file((tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_raw.fasta").c_str(), sequence_ids, ll, rr);
#endif
    if(center_block.container.size() <= 0) return "";
    // Part 0: sort the blocks
    std::sort(other_blocks.begin(), other_blocks.end(), std::greater<block_irreg>());
    // Part 1: left_blocks = [], right_blocks = [], down_blocks = []
    std::vector <block_irreg> left_blocks, right_blocks, down_blocks;
    item_int other_ids, * real_other_ids = nullptr, *real_dispatch_lr = nullptr;
    item_int ll_dispatch_ll, ll_dispatch_rr, rr_dispatch_ll, rr_dispatch_rr, other_dispatch_ll, other_dispatch_rr;
    set_int other_set, center_items, all_set;
    center_int center_block_info;
    map_int mp_ids; // not in other blocks: real_id->id
    block_irreg* left_block, * right_block, * down_block;
    std::string file_left, file_right, file_down, final_file;
    gap_load_int gap_preload_in_main;
    int_ center_block_len = center_block.length, max_len = -1, min_len = INT_MAX;
    gap_tree center_tree;
    bool use_center_tree = false;
    // Part 2: check the number in center_block is same as the cluster
    for (size_t i = 0; i < sequence_ids.size(); ++i) mp_ids[sequence_ids[i]] = i, all_set.insert(sequence_ids[i]);
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
        Frag_dispatch(final_file_name.c_str(), sequence_ids, ll, rr, this_tree);
        return ""; // can not find invalid block
    }
    if (min_len != max_len)
    {
        // block is irregular, align it!
        fprintf(stderr, "Info: this block is not same (although contains all sequences). Program will align this block.\n");
        std::string center_file_name = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_block.fasta";
        use_center_tree = true;
        Center_block_dispatch(center_file_name.c_str(), center_block_info, center_tree);
    }
    if (sequence_ids.size() != center_block_info.size())
    {
        assert(sequence_ids.size() > center_block_info.size());
        for (auto& item : center_block_info) center_items.insert(item.first);
        std::set_difference(all_set.begin(), all_set.end(), center_items.begin(), center_items.end(), std::inserter(other_set, other_set.begin()));
        all_set.clear();
        // Part 2.3: determine how many sequences will be in dispatch mode
        real_dispatch_lr = new item_int;
        // Part 2.4: determine real_dispatch_lr = center_block_info.first
        real_dispatch_lr -> reserve(center_block_info.size());
        for(auto &item: center_block_info) real_dispatch_lr -> emplace_back(item.first);
    }
    else
    {
        all_set.clear();
        real_dispatch_lr = &sequence_ids;
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
    if (left_blocks.size())
    {
        std::sort(left_blocks.begin(), left_blocks.end());
        block_irreg left_center = left_blocks.back();
        left_blocks.pop_back(); // remove last one
        ll_dispatch_ll.resize(real_dispatch_lr -> size());
        ll_dispatch_rr.resize(real_dispatch_lr -> size());
        for (size_t i = 0; i < real_dispatch_lr -> size(); ++i)
        {
            ll_dispatch_rr[i] = center_block_info[real_dispatch_lr -> at(i)].first;
            ll_dispatch_ll[i] = ll[mp_ids[real_dispatch_lr -> at(i)]];
        }
        // dispatch to left next-level
#if _DEBUG
        unsigned int rand_left = rand_gen();
        fprintf(stderr, "Recurse: left decompose %d %u, dispatch to %u\n", depth, part, rand_left);
        file_left = decompose_without_SSW(*real_dispatch_lr, left_center, left_blocks, ll_dispatch_ll, ll_dispatch_rr, depth + 1, rand_left, cluster_mp_ids, *this_tree.left_tree);
#else
        file_left = decompose_without_SSW(*real_dispatch_lr, left_center, left_blocks, ll_dispatch_ll, ll_dispatch_rr, depth + 1, rand_gen(), cluster_mp_ids, *this_tree.left_tree);
#endif
    }
    else
    {
        file_left = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_left.fasta";
        ll_dispatch_ll.resize(real_dispatch_lr -> size());
        ll_dispatch_rr.resize(real_dispatch_lr -> size());
        for (size_t i = 0; i < real_dispatch_lr -> size(); ++i)
        {
            ll_dispatch_ll[i] = ll[mp_ids[real_dispatch_lr -> at(i)]];
            ll_dispatch_rr[i] = center_block_info[real_dispatch_lr -> at(i)].first;
        }
#if _DEBUG
        fprintf(stderr, "Recurse: dispatch_left %d %u\n", depth, part);
#endif
        Frag_dispatch((file_left).c_str(), *real_dispatch_lr, ll_dispatch_ll, ll_dispatch_rr, *this_tree.left_tree);
    }

    if (right_blocks.size())
    {
        std::sort(right_blocks.begin(), right_blocks.end());
        block_irreg right_center = right_blocks.back();
        right_blocks.pop_back(); // remove last one
        rr_dispatch_ll.resize(real_dispatch_lr -> size());
        rr_dispatch_rr.resize(real_dispatch_lr -> size());
        for (size_t i = 0; i < real_dispatch_lr -> size(); ++i)
        {
            rr_dispatch_ll[i] = center_block_info[real_dispatch_lr -> at(i)].second;
            rr_dispatch_rr[i] = rr[mp_ids[real_dispatch_lr -> at(i)]];
        }
#if _DEBUG
        unsigned int rand_right = rand_gen();
        fprintf(stderr, "Recurse: right decompose %d %u, dispatch to %u\n", depth, part, rand_right);
        file_right = decompose_without_SSW(*real_dispatch_lr, right_center, right_blocks, rr_dispatch_ll, rr_dispatch_rr, depth + 1, rand_right, cluster_mp_ids, *this_tree.right_tree);
#else
        file_right = decompose_without_SSW(*real_dispatch_lr, right_center, right_blocks, rr_dispatch_ll, rr_dispatch_rr, depth + 1, rand_gen(), cluster_mp_ids, *this_tree.right_tree);
#endif
    }
    else
    {
        file_right = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_right.fasta";
        rr_dispatch_ll.resize(real_dispatch_lr -> size());
        rr_dispatch_rr.resize(real_dispatch_lr -> size());
        for (size_t i = 0; i < real_dispatch_lr -> size(); ++i)
        {
            rr_dispatch_ll[i] = center_block_info[real_dispatch_lr -> at(i)].second;
            rr_dispatch_rr[i] = rr[mp_ids[real_dispatch_lr -> at(i)]];
        }
#if _DEBUG
        fprintf(stderr, "Recurse: dispatch_right %d %u\n", depth, part);
#endif
        Frag_dispatch((file_right).c_str(), *real_dispatch_lr, rr_dispatch_ll, rr_dispatch_rr, *this_tree.right_tree);
    }

    // merge left and right gaps
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

    if (other_set.size())
    {
        real_other_ids = new item_int;
        real_other_ids->insert(real_other_ids->end(), other_set.begin(), other_set.end());
        if(down_blocks.size())
        {
            this_tree.alloc_down();
            std::sort(down_blocks.begin(), down_blocks.end());
            block_irreg other_center = down_blocks.back();
            down_blocks.pop_back(); // remove last one
            other_dispatch_ll.resize(real_other_ids->size());
            other_dispatch_rr.resize(real_other_ids->size());
            for (size_t i = 0; i < real_other_ids->size(); ++i)
            {
                other_dispatch_ll[i] = ll[mp_ids[real_other_ids -> at(i)]];
                other_dispatch_rr[i] = rr[mp_ids[real_other_ids -> at(i)]];
            }
#if _DEBUG
            unsigned int rand_down = rand_gen();
            fprintf(stderr, "Recurse: down decompose %d %u, dispatch to %u\n", depth, part, rand_down);
            file_down = decompose_without_SSW(*real_other_ids, other_center, down_blocks, other_dispatch_ll, other_dispatch_rr, depth + 1, rand_down, cluster_mp_ids, *this_tree.down_tree);
#else
            file_down = decompose_without_SSW(*real_other_ids, other_center, down_blocks, other_dispatch_ll, other_dispatch_rr, depth + 1, rand_gen(), cluster_mp_ids, *this_tree.down_tree);
#endif
            // merge other gaps
            this_tree.merge_down();
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
#if _DEBUG
            fprintf(stderr, "Recurse: dispatch_down %d %u\n", depth, part);
#endif
            Frag_dispatch((file_down).c_str(), *real_other_ids, other_dispatch_ll, other_dispatch_rr, *this_tree.down_tree);
            // merge other gaps
            this_tree.merge_down();
        }
#if _DEBUG
        fprintf(stderr, "Recurse: dispatch_down merging %d %u\n", depth, part);
#endif
        // merge real_dispatch and real_other_ids
        TwoProfile_dispatch(*real_dispatch_lr, *real_other_ids, ll, rr, mp_ids, this_tree, depth, part);
        delete real_other_ids;
    }
    final_file = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_final.fasta";
#if _DEBUG
    // Part 5: merge and print to file
    write_gaps_to_file_in_dispatch(final_file.c_str(), sequence_ids, ll, rr, this_tree); // may be useless in program, can be deprecated
    assert(aligned(sequence_ids, ll, rr, this_tree));
#endif
    // clean up pointer items
    if(real_dispatch_lr != &sequence_ids) delete real_dispatch_lr;
    return final_file;
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
std::string decompose_with_SSW(item_int& sequence_ids,
                               block_irreg& center_block,
                               std::vector <block_irreg>& other_blocks,
                               item_int& ll, item_int& rr,
                               const int depth, const unsigned int part, map_int &cluster_mp_ids, gap_tree &this_tree)
{
#if _DEBUG
    write_raw_to_file((tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_raw.fasta").c_str(), sequence_ids, ll, rr);
#endif
    if(center_block.container.size() <= 0) return "";
    // Part 0: sort the blocks
    std::sort(other_blocks.begin(), other_blocks.end(), std::greater<block_irreg>());
    // Part 1: left_blocks = [], right_blocks = [], down_blocks = []
    std::vector <block_irreg> left_blocks, right_blocks, down_blocks;
    item_int other_ids, * real_other_ids = nullptr, *real_dispatch_lr = nullptr;
    item_int ll_dispatch_ll, ll_dispatch_rr, rr_dispatch_ll, rr_dispatch_rr, other_dispatch_ll, other_dispatch_rr;
    set_int other_set, center_items, all_set;
    center_int center_block_info;
    map_int mp_ids; // not in other blocks: real_id->id
    block_irreg* left_block, * right_block, * down_block;
    std::string file_left, file_right, file_down, final_file;
    gap_load_int gap_preload_in_main;
    int_ center_block_len = center_block.length, max_len = -1, min_len = INT_MAX;
    gap_tree center_tree;
    bool use_center_tree = false;
    // Part 2: check the number in center_block is same as the cluster
    for (size_t i = 0; i < sequence_ids.size(); ++i) mp_ids[sequence_ids[i]] = i, all_set.insert(sequence_ids[i]);
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
        Frag_dispatch(final_file_name.c_str(), sequence_ids, ll, rr, this_tree);
        return ""; // can not find invalid block
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
        Findsim_and_extension(other_ids, not_same_len_seq, center_items, center_block.center, ll, rr, mp_ids, valid_seq_id, valid_beg, valid_end, 
                              center_block_info, other_set, center_tree, depth, part);
        // Part 2.3: determine how many sequences will be in dispatch mode
        real_dispatch_lr = new item_int;
        // Part 2.4: determine real_dispatch_lr = center_block_info.first
        real_dispatch_lr -> reserve(center_block_info.size());
        for(auto &item: center_block_info) real_dispatch_lr -> emplace_back(item.first);
#if _DEBUG
        write_gaps_to_file_center((tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + ".ssw.fasta").c_str(), center_block_info, this_tree);
#endif
    }
    else
    {
        all_set.clear();
        real_dispatch_lr = &sequence_ids;
        if (min_len != max_len)
        {
            // block is irregular, align it!
            fprintf(stderr, "Info: this block is not same (although contains all sequences). Program will align this block.\n");
            std::string center_file_name = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_block.fasta";
            use_center_tree = true;
            Center_block_dispatch(center_file_name.c_str(), center_block_info, center_tree);
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
    if (left_blocks.size())
    {
        std::sort(left_blocks.begin(), left_blocks.end());
        block_irreg left_center = left_blocks.back();
        left_blocks.pop_back(); // remove last one
        ll_dispatch_ll.resize(real_dispatch_lr -> size());
        ll_dispatch_rr.resize(real_dispatch_lr -> size());
        for (size_t i = 0; i < real_dispatch_lr -> size(); ++i)
        {
            ll_dispatch_rr[i] = center_block_info[real_dispatch_lr -> at(i)].first;
            ll_dispatch_ll[i] = ll[mp_ids[real_dispatch_lr -> at(i)]];
        }
        // dispatch to left next-level
#if _DEBUG
        unsigned int rand_left = rand_gen();
        fprintf(stderr, "Recurse: left decompose %d %u, dispatch to %u\n", depth, part, rand_left);
        file_left = decompose_with_SSW(*real_dispatch_lr, left_center, left_blocks, ll_dispatch_ll, ll_dispatch_rr, depth + 1, rand_left, cluster_mp_ids, *this_tree.left_tree);
#else
        file_left = decompose_with_SSW(*real_dispatch_lr, left_center, left_blocks, ll_dispatch_ll, ll_dispatch_rr, depth + 1, rand_gen(), cluster_mp_ids, *this_tree.left_tree);
#endif
    }
    else
    {
        file_left = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_left.fasta";
        ll_dispatch_ll.resize(real_dispatch_lr -> size());
        ll_dispatch_rr.resize(real_dispatch_lr -> size());
        for (size_t i = 0; i < real_dispatch_lr -> size(); ++i)
        {
            ll_dispatch_ll[i] = ll[mp_ids[real_dispatch_lr -> at(i)]];
            ll_dispatch_rr[i] = center_block_info[real_dispatch_lr -> at(i)].first;
        }
#if _DEBUG
        fprintf(stderr, "Recurse: dispatch_left %d %u\n", depth, part);
#endif
        Frag_dispatch((file_left).c_str(), *real_dispatch_lr, ll_dispatch_ll, ll_dispatch_rr, *this_tree.left_tree);
    }

    if (right_blocks.size())
    {
        std::sort(right_blocks.begin(), right_blocks.end());
        block_irreg right_center = right_blocks.back();
        right_blocks.pop_back(); // remove last one
        rr_dispatch_ll.resize(real_dispatch_lr -> size());
        rr_dispatch_rr.resize(real_dispatch_lr -> size());
        for (size_t i = 0; i < real_dispatch_lr -> size(); ++i)
        {
            rr_dispatch_ll[i] = center_block_info[real_dispatch_lr -> at(i)].second;
            rr_dispatch_rr[i] = rr[mp_ids[real_dispatch_lr -> at(i)]];
        }
#if _DEBUG
        unsigned int rand_right = rand_gen();
        fprintf(stderr, "Recurse: right decompose %d %u, dispatch to %u\n", depth, part, rand_right);
        file_right = decompose_with_SSW(*real_dispatch_lr, right_center, right_blocks, rr_dispatch_ll, rr_dispatch_rr, depth + 1, rand_right, cluster_mp_ids, *this_tree.right_tree);
#else
        file_right = decompose_with_SSW(*real_dispatch_lr, right_center, right_blocks, rr_dispatch_ll, rr_dispatch_rr, depth + 1, rand_gen(), cluster_mp_ids, *this_tree.right_tree);
#endif
    }
    else
    {
        file_right = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_right.fasta";
        rr_dispatch_ll.resize(real_dispatch_lr -> size());
        rr_dispatch_rr.resize(real_dispatch_lr -> size());
        for (size_t i = 0; i < real_dispatch_lr -> size(); ++i)
        {
            rr_dispatch_ll[i] = center_block_info[real_dispatch_lr -> at(i)].second;
            rr_dispatch_rr[i] = rr[mp_ids[real_dispatch_lr -> at(i)]];
        }
#if _DEBUG
        fprintf(stderr, "Recurse: dispatch_right %d %u\n", depth, part);
#endif
        Frag_dispatch((file_right).c_str(), *real_dispatch_lr, rr_dispatch_ll, rr_dispatch_rr, *this_tree.right_tree);
    }

    // merge left and right gaps
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

    if (other_set.size())
    {
        real_other_ids = new item_int;
        real_other_ids->insert(real_other_ids->end(), other_set.begin(), other_set.end());
        if(down_blocks.size())
        {
            this_tree.alloc_down();
            std::sort(down_blocks.begin(), down_blocks.end());
            block_irreg other_center = down_blocks.back();
            down_blocks.pop_back(); // remove last one
            other_dispatch_ll.resize(real_other_ids->size());
            other_dispatch_rr.resize(real_other_ids->size());
            for (size_t i = 0; i < real_other_ids->size(); ++i)
            {
                other_dispatch_ll[i] = ll[mp_ids[real_other_ids -> at(i)]];
                other_dispatch_rr[i] = rr[mp_ids[real_other_ids -> at(i)]];
            }
#if _DEBUG
            unsigned int rand_down = rand_gen();
            fprintf(stderr, "Recurse: down decompose %d %u, dispatch to %u\n", depth, part, rand_down);
            file_down = decompose_with_SSW(*real_other_ids, other_center, down_blocks, other_dispatch_ll, other_dispatch_rr, depth + 1, rand_down, cluster_mp_ids, *this_tree.down_tree);
#else
            file_down = decompose_with_SSW(*real_other_ids, other_center, down_blocks, other_dispatch_ll, other_dispatch_rr, depth + 1, rand_gen(), cluster_mp_ids, *this_tree.down_tree);
#endif
            // merge other gaps
            this_tree.merge_down();
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
#if _DEBUG
            fprintf(stderr, "Recurse: dispatch_down %d %u\n", depth, part);
#endif
            Frag_dispatch((file_down).c_str(), *real_other_ids, other_dispatch_ll, other_dispatch_rr, *this_tree.down_tree);
            // merge other gaps
            this_tree.merge_down();
        }
#if _DEBUG
        fprintf(stderr, "Recurse: dispatch_down merging %d %u\n", depth, part);
#endif
        // merge real_dispatch and real_other_ids
        TwoProfile_dispatch(*real_dispatch_lr, *real_other_ids, ll, rr, mp_ids, this_tree, depth, part);
        delete real_other_ids;
    }
    final_file = tmp_path_prefix + std::to_string(depth) + "_" + std::to_string(part) + "_final.fasta";
#if _DEBUG
    // Part 5: merge and print to file
    write_gaps_to_file_in_dispatch(final_file.c_str(), sequence_ids, ll, rr, this_tree); // may be useless in program, can be deprecated
    assert(aligned(sequence_ids, ll, rr, this_tree));
#endif
    // clean up pointer items
    if(real_dispatch_lr != &sequence_ids) delete real_dispatch_lr;
    return final_file;
}
