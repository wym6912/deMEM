#include "gaps.hpp"

gap_system::gap_system() 
{
	sum_gaps = 0;
	max_cap = sizeof(int_) == 4 ? INT_MAX : LONG_LONG_MAX;
}

gap_system::gap_system(int_ seq_size)
{
	if(seq_size < 0)
	{
		fprintf(stderr, "Error: initial gap system size (= %d) is less than 0. Program will exit.\n", seq_size);
		exit(1);
	}
	max_cap = seq_size + 1;
	sum_gaps = 0;
}

gap_system::~gap_system()
{
	gaps.clear();
}

void gap_system::add_gaps(const int_ real_place, const int_ gaps)
{
	int_ tree_place = real_place + 1;
	if (tree_place < 0)
	{
		fprintf(stderr, "Error: real place %d is less than 0. Program will exit.\n", tree_place);
		exit(1);
	}
	else if (tree_place >= this->max_cap)
	{
		max_cap = tree_place + 1;
		fprintf(stderr, "Warning: will extend the gap size to %d.\n", this->max_cap);
	}
	if (gaps == 0) return;
	this->gaps[tree_place] += gaps;
	sum_gaps += gaps;
}

const int_ gap_system::get_gaps(const int_ real_place)
{
	int_ gap_place = real_place + 1;
	if (gap_place < 0 || gap_place >= this->max_cap)
	{
		fprintf(stderr, "Error: real place %d is invaild. Program will exit.\n", gap_place);
#if _DEBUG
		throw "out_of_range";
#endif
		exit(1);
	}
	const auto& it = this->gaps.find(gap_place);
	if (it == this->gaps.end() || it->first != gap_place) return 0;
	return it->second;
}

int gap_system::deinit(const int_ seq_size)
{
	if (seq_size < 0)
	{
		fprintf(stderr, "Error: initial gap system size (= %d) is less than 0. Program will exit.\n", seq_size);
		exit(1);
	}
	gaps.clear();
	max_cap = seq_size + 1;
	sum_gaps = 0;
	return 0;
}

int gap_system::deinit()
{
	gaps.clear();
	sum_gaps = 0;
	max_cap = sizeof(int_) == 4 ? INT_MAX : LONG_LONG_MAX;
	return 0;
}

const int_ gap_system::gap_sum()
{
	return sum_gaps;
}

void gap_system::set_gaps(const int_ real_place, const int_ gaps)
{
	int_ tree_place = real_place + 1, gap_before = 0;
	if (tree_place < 0)
	{
		fprintf(stderr, "Error: real place %d is less than 0. Program will exit.\n", tree_place);
		exit(1);
	}
	else if (tree_place >= this->max_cap)
	{
		max_cap = tree_place + 1;
		fprintf(stderr, "Warning: will extend the gap size to %d.\n", this->max_cap);
	}
	const auto &it = this->gaps.find(tree_place);
	if(it == this->gaps.end()) gap_before = 0;
	else gap_before = it->second;
	this->gaps[tree_place] = gaps;
	sum_gaps += gaps - gap_before;
}

void gap_tree::add_gaps(const int_ seq_id, const int_ real_place, const int_ gaps)
{
	int_ tree_place = real_place + 1;
	if (tree_place < 0)
	{
		fprintf(stderr, "Error: real place %d is less than 0. Program will exit.\n", tree_place);
		exit(1);
	}
	else if (tree_place >= this->this_tree[seq_id].max_cap)
	{
		this_tree[seq_id].max_cap = tree_place + 1;
		fprintf(stderr, "Warning: will extend the gap size to %d.\n", this->this_tree[seq_id].max_cap);
	}
	if (gaps == 0) return;
	this->this_tree[seq_id].gaps[tree_place] += gaps;
	this->this_tree[seq_id].sum_gaps += gaps;
}

const int_ gap_tree::get_gaps(const int_ seq_id, const int_ real_place)
{
	return this_tree[seq_id].get_gaps(real_place);
}

gap_system& gap_tree::get_gaps(const int_ seq_id)
{
	return this_tree[seq_id];
}

void gap_tree::set_gaps(const int_ seq_id, const int_ real_place, const int_ gaps)
{
	int_ tree_place = real_place + 1, gap_before = 0;
	if (tree_place < 0)
	{
		fprintf(stderr, "Error: real place %d is less than 0. Program will exit.\n", tree_place);
		exit(1);
	}
	else if (tree_place >= this->this_tree[seq_id].max_cap)
	{
		this_tree[seq_id].max_cap = tree_place + 1;
		fprintf(stderr, "Warning: will extend the gap size to %d.\n", this->this_tree[seq_id].max_cap);
	}
	const auto &it = this->this_tree[seq_id].gaps.find(tree_place);
	if(it == this->this_tree[seq_id].gaps.end()) gap_before = 0;
	else gap_before = it->second;
	this->this_tree[seq_id].gaps[tree_place] = gaps;
	this_tree[seq_id].sum_gaps += gaps - gap_before;
}

int gap_tree::deinit(const int_ seq_id, const int_ seq_size)
{
	return this_tree.at(seq_id).deinit(seq_size);
}

int gap_tree::alloc_lr()
{
	left_tree  = new gap_tree;
	right_tree = new gap_tree;
	return 0;
}

int gap_tree::alloc_down()
{
	down_tree = new gap_tree;
	return 0;
}

gap_tree::gap_tree()
{
	left_tree = right_tree = down_tree = nullptr;
}

gap_tree::~gap_tree()
{
	if(left_tree)  delete left_tree;
	if(right_tree) delete right_tree;
	if(down_tree)  delete down_tree;
	this_tree.clear();
}

int gap_tree::merge(gap_tree* tree)
{
	if (tree == nullptr) return 1;
	for(auto &tree_node: tree->this_tree)
	{
		const auto &seq_id = tree_node.first;
		for(auto &place_with_gaps: tree_node.second.gaps)
		{
			this_tree[seq_id].gaps[place_with_gaps.first] += place_with_gaps.second;
			this_tree[seq_id].sum_gaps += place_with_gaps.second;
		}
	}
	delete tree;
	return 0;
}

int gap_tree::merge(gap_tree& tree)
{
	for (auto& tree_node : tree.this_tree)
	{
		const auto& seq_id = tree_node.first;
		for (auto& place_with_gaps : tree_node.second.gaps)
		{
			this_tree[seq_id].gaps[place_with_gaps.first] += place_with_gaps.second;
			this_tree[seq_id].sum_gaps += place_with_gaps.second;
		}
	}
	return 0;
}


int gap_tree::merge_lr()
{
	int val = merge(left_tree) || merge(right_tree);
	left_tree = right_tree = nullptr;
	return val;
}

int gap_tree::merge_down()
{
	int val = merge(down_tree);
	down_tree = nullptr;
	return val;
}
