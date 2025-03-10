#ifndef __GAPS__
#define __GAPS__

#include <vector>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <map>
#include <set>
#include <exception>
#include <utility>
#include <climits>

#if (defined(_WIN32) || defined(_WIN64))
#define LONG_LONG_MAX LLONG_MAX
#endif

typedef int int_;
typedef std::vector<int_> item_int;
typedef std::set<int_> set_int;
typedef std::unordered_map<int_, int_> map_int;
typedef std::pair<int_, int_> pair_int;

typedef map_int gap_item;
typedef std::unordered_map<int_, pair_int> gap_load_int;
typedef gap_load_int center_int;

typedef struct gap_system
{
	gap_item gaps; // gap place is behind of sequence
	int_ max_cap, sum_gaps;
	gap_system();
	gap_system(int_ seq_size);
	~gap_system();
	      void add_gaps(const int_ real_place, const int_ gaps);
	      void set_gaps(const int_ real_place, const int_ gaps);
	const int_ get_gaps(const int_ real_place);
	      int  deinit  (const int_ seq_size);
		  int  deinit  ();
	const int_ gap_sum ();
}gaps;

typedef struct gap_tree
{
	std::unordered_map<int, gaps> this_tree;
	gap_tree *left_tree, *right_tree, *down_tree;
	gap_tree();
	~gap_tree();
	      void  add_gaps(const int_ seq_id, const int_ real_place, const int_ gaps);
	const int_  get_gaps(const int_ seq_id, const int_ real_place);
	      gaps& get_gaps(const int_ seq_id);
		  void  set_gaps(const int_ seq_id, const int_ real_place, const int_ gaps);
	      int   deinit  (const int_ seq_id, const int_ seq_size);
	      int   alloc_lr();
		  int   alloc_down();
		  int   merge_lr();
		  int   merge_down();
		  int   merge(gap_tree* tree);
		  int   merge(gap_tree& tree);

} gap_tree;

#endif