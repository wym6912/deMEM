#ifndef __GAPIO__
#define __GAPIO__
#include "gaps.hpp"
#include "../external/kseq/kseq.h"
#if _DEBUG
#include <cassert>
#include <cstring>
#endif

#if (defined(__linux__) || defined(__APPLE__))
#include <unistd.h>
#else
#include "../include/unistd.h"
#endif

KSEQ_INIT(int, read)

void read_seq();

// 0 is ok, other is bad
int read_gaps_from_file(const char *file_name, item_int& sequence_id, item_int& ll, item_int& rr, gap_tree &this_tree);
int read_gaps_from_file(const char* file_name, center_int& center_info, gap_tree& this_tree);
int write_gaps_to_file(const char *file_name, item_int& sequence_id, item_int& ll, item_int& rr, bool first, gap_tree &this_tree);
int write_gaps_to_file(const char* file_name, center_int &center_info, bool first, gap_tree &this_tree);
int write_gaps_to_file_with_lr(const char* file_name, item_int& sequence_id, item_int& ll, item_int& rr, item_int& ll_border, item_int& rr_border, gap_tree& this_tree);
int write_gaps_to_file_in_dispatch(const char *file_name, item_int& sequence_id, item_int& ll, item_int& rr, gap_tree &this_tree);
int write_gaps_to_file_center(const char *file_name, center_int &center_block, gap_tree &this_tree);
int write_raw_to_file(const char *file_name, item_int& sequence_id, item_int& ll, item_int& rr);
int read_gaps_from_merged_profile(const char* file_name, item_int& sequence_id1, item_int& sequence_id2, item_int& ll, item_int& rr, map_int& mp, gap_tree &this_tree);
void write_seq(int_ &seq_id, int_ &ll, int_ &rr, gap_system &gap_load, FILE *f);
inline void write_preload_seq(int_ &seq_id, int_ &ll, int_ &rr, int_ &border_l, int_ &border_r, pair_int &gap_border, gap_system &gap_load, FILE *f);

int gap_system_init(const int& sequence_total); // 0 is ok, other is bad
int gap_system_reinit(item_int& v, int size);
int gap_system_free();

#endif
