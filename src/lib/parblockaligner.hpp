#ifndef __PARBLOCKALIGN__
#define __PARBLOCKALIGN__

#include "libblockaligner.hpp"

std::string decompose_with_SSW(item_int& sequence_ids, block_irreg& center_block, std::vector <block_irreg>& other_blocks, item_int& ll, item_int& rr, map_int &cluster_mp_ids, gap_tree &this_tree, int threads);
std::string decompose_without_SSW(item_int& sequence_ids, block_irreg& center_block, std::vector <block_irreg>& other_blocks, item_int& ll, item_int& rr, map_int &cluster_mp_ids, gap_tree &this_tree, int threads);

#endif