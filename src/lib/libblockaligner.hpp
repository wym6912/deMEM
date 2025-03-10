#ifndef __BLOCKALIGN__
#define __BLOCKALIGN__

#include "ext_prog_call.hpp"
#include "seqio.hpp"
#include "irreg_part.hpp"
#include "libfile_v2.hpp"
#include "../external/ssw/src/ssw_cpp.hpp"
#include "../external/ssw/src/ssw.h"
#include <stdexcept>
#include <functional>
#include <random>
#include <bitset>

typedef std::pair<int_, int_> ssw_end_beg;
typedef std::pair<item_int, std::map<int_, ssw_end_beg> > findsim_args;

std::string decompose_with_SSW(item_int& sequence_ids, block_irreg& center_block, std::vector <block_irreg>& other_blocks, item_int& ll, item_int& rr, const int depth, const unsigned int part, map_int &cluster_mp_ids, gap_tree &this_tree);
std::string decompose_without_SSW(item_int& sequence_ids, block_irreg& center_block, std::vector <block_irreg>& other_blocks, item_int& ll, item_int& rr, const int depth, const unsigned int part, map_int &cluster_mp_ids, gap_tree &this_tree);
int Profiles_dispatch(std::vector<item_int> &files, std::vector<std::string> &file_names, const char* out_name, std::string &center_file_name, std::string &file_name);
int TwoProfile_dispatch(const std::string& File1, const std::string& File2);
int Frag_dispatch_full_seq(const char* File, const item_int& seq_id);

// used in parblockaligner.cpp
void WriteSeq_SSW_Frag(FILE* f, const int_& ll, const int_& rr, const int_& real_id);
void TwoProfile_writeseq(const char* name, item_int &ids, item_int &ll, item_int &rr, map_int &mp, gap_tree &this_tree);
void refinement(item_int &seq_ids, item_int& border, gap_tree& this_tree);
extern int8_t blosum50[], blosum62[], mat_nt[], aa_table[], nt_table[];

#if _DEBUG
bool aligned(item_int &seq_ids, item_int &ll, item_int &rr, gap_tree &gap_tree_main);
#endif

#endif