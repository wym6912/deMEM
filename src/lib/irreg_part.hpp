#ifndef __IRREG_PART__
#define __IRREG_PART__

#include "gaps.hpp"

#include <tuple>
#include <vector>
#include <map>
#include <string>
#include <cassert>
#include <algorithm>

typedef std::tuple<int_, int_, int_> sequence_irreg;

typedef struct block_irreg
{
	int_ length, center;
	long long area;
	/* element tuple: first is sequence id, second is sequence start place, third is extension length
	   extension length = (block length - element length)
	 */
	std::vector <sequence_irreg> container;
	block_irreg()
	{
		this->length = 0;
		this->center = 0;
		this->area = 0;
	}

	bool operator < (const block_irreg& b) const noexcept;
	bool operator > (const block_irreg& b) const noexcept;

	void calc_area();
	void freecontainer();
	void sortcontainer();
} block_irreg;

typedef struct cluster_irreg
{
	std::vector <int_> sequences; // sequence ids, only changed once
	block_irreg* common_block; // only one block in blocks
	std::vector <block_irreg> sub_blocks; // many blocks in sequences
	cluster_irreg()
	{
		this->common_block = nullptr;
	}
	cluster_irreg(cluster_irreg&& c) noexcept : sequences(std::move(c.sequences)), common_block(c.common_block), sub_blocks(std::move(c.sub_blocks)) {}
	cluster_irreg(const cluster_irreg& c) = default;
} cluster_irreg;

template <typename T> void insertion_sort(T* data, size_t T_size);
template <typename _VectorType, typename _ItemType> void insertion_sort(_VectorType data);

int_ find_center(const std::vector <sequence_irreg>& sequences);

#endif