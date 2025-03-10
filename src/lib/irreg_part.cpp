#include "irreg_part.hpp"

template <typename T>
void insertion_sort(T* data, size_t sz)
{
	for (size_t i = 1, j; i < sz; i++) {
		T key = data[i];
		j = i - 1;
		while ((j >= 0) && (key < data[j])) {
			data[j + 1] = data[j];
			--j;
		}
		data[j + 1] = key;
	}
}

template <typename _VectorType, typename _ItemType>
void insertion_sort(_VectorType data)
{
	for (size_t i = 1, j; i < data.size(); i++) {
		_ItemType key = data[i];
		j = i - 1;
		while ((j >= 0) && (key < data[j])) {
			data[j + 1] = data[j];
			--j;
		}
		data[j + 1] = key;
	}
}

void block_irreg::calc_area()
{
	this->area = 0;
	for (size_t i = 0; i < this->container.size(); ++i)
	{
		this->area += this->length
			- (long long)std::abs(std::get<1>(this->container[i]) - this->center)
			- std::get<2>(this->container[i]); // this will be modified by calling SSW
	}
}

void block_irreg::freecontainer()
{
	this->container.clear();
}

void block_irreg::sortcontainer()
{
	insertion_sort<std::vector <sequence_irreg>, sequence_irreg>(this->container);
}

bool block_irreg::operator < (const block_irreg& b) const noexcept
{
	// first compare with area (length * numbers), next compare the length
	if(this->area != b.area) return this->area < b.area;
	return this->length < b.length;
}

bool block_irreg::operator > (const block_irreg& b) const noexcept
{
	// first compare with area (length * numbers), next compare the length
	if (this->area != b.area) return this->area > b.area;
	return this->length > b.length;
}

int_ find_center(const std::vector <sequence_irreg>& sequences)
{
	std::map<int_, int_> appeal_time;
	std::vector <int_> max_data;
	int_ max_time = 0;
	long long sum = 0;
	appeal_time.clear();
	max_data.clear();
	for (size_t i = 0; i < sequences.size(); ++i)
	{
		++appeal_time[std::get<1>(sequences[i])];
		sum += std::get<1>(sequences[i]) - std::get<2>(sequences[i]);
	}
	for (std::map<int_, int_>::iterator it = appeal_time.begin(); it != appeal_time.end(); ++it) max_time = (std::max)(it->second, max_time);
	for (std::map<int_, int_>::iterator it = appeal_time.begin(); it != appeal_time.end(); ++it)
		if (max_time == it->second)
			max_data.push_back(it->first);
	std::sort(max_data.begin(), max_data.end());
	if (max_data.size() & 1) return max_data[max_data.size() >> 1];
	else
	{
		assert((max_data.size() != 0));
		sum /= max_data.size();
		return (std::abs(sum - max_data[max_data.size() >> 1]) > std::abs(sum - max_data[(max_data.size() >> 1) - 1])) ?
			max_data[(max_data.size() >> 1) - 1] : max_data[max_data.size() >> 1];
	}
}