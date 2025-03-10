#include "seqio.hpp"

extern char* c_input, * c_output, * c_seq, * thread_str;
extern char** seq, ** comment;
extern item_int first_end, this_beg, this_end;
extern int threads_num, block_sequence_type;
extern int_ seq_sum;

extern gap_system* final_gap_stores;

#if _DEBUG
bool simple_char(char res, char des)
{
	switch (res)
	{
		case 'A':
			if (des == 'V' || des == 'D' || des == 'H' || des == 'W' || des == 'M' || des == 'R' || des == 'N') return true;
			break;
		case 'C':
			if (des == 'V' || des == 'B' || des == 'H' || des == 'S' || des == 'M' || des == 'Y' || des == 'N') return true;
			break;
		case 'G':
			if (des == 'V' || des == 'D' || des == 'B' || des == 'S' || des == 'K' || des == 'R' || des == 'N') return true;
			break;
		case 'T':
			if (des == 'D' || des == 'B' || des == 'H' || des == 'W' || des == 'K' || des == 'Y' || des == 'N') return true;
			break;
	}
	return false;
}

bool is_same_char(char ch1, char ch2)
{
	if (ch1 == ch2) return true;
	if (block_sequence_type == 2)
	{
		if (ch1 == 'A' || ch1 == 'T' || ch1 == 'C' || ch1 == 'G')
		{
			return simple_char(ch1, ch2);
		}
		else if (ch2 == 'A' || ch2 == 'T' || ch2 == 'C' || ch2 == 'G')
		{
			return simple_char(ch2, ch1);
		}
	}
	return false;
}
#endif

void read_seq()
{
	fprintf(stderr, "Reading sequence file...");
	FILE* seq_point = fopen(c_seq, "r");
	if (seq_point == NULL) { fprintf(stderr, "Error: file %s cannot open. Program will exit\n", c_seq); exit(1); }
	kseq_t* seq_file = kseq_init(fileno(seq_point));
	seq_sum = 0;
	int max_alloc = 20000;
	seq = (char**)malloc(sizeof(char*) * max_alloc);
	comment = (char**)malloc(sizeof(char*) * max_alloc);
	first_end.reserve(max_alloc);
	if (seq == NULL || comment == NULL) { fprintf(stderr, "Error: Cannot allocate enough space on reading sequences. Program will exit.\n"); exit(1); }
	char** tmp_alloc;
	while (kseq_read(seq_file) >= 0)
	{
		if (max_alloc < seq_sum)
		{
			max_alloc <<= 1;
			if (max_alloc < 0) { fprintf(stderr, "Error: Sequence is too many. Program will exit.\n"); exit(1); }
			first_end.reserve(max_alloc);
			tmp_alloc = seq;
			seq = (char**)calloc(max_alloc, sizeof(char*));
			if (seq == NULL) { fprintf(stderr, "Error: Cannot allocate enough space on reading sequences. Program will exit.\n"); exit(1); }
			memcpy(seq, tmp_alloc, sizeof(char*) * (max_alloc >> 1));
			free(tmp_alloc);
			tmp_alloc = comment;
			comment = (char**)calloc(max_alloc, sizeof(char*));
			if (comment == NULL) { fprintf(stderr, "Error: Cannot allocate enough space on reading sequences. Program will exit.\n"); exit(1); }
			memcpy(comment, tmp_alloc, sizeof(char*) * (max_alloc >> 1));
			free(tmp_alloc);
		}
		seq[seq_sum] = (char*)malloc(sizeof(char) * (seq_file->seq.l + 1));
		if (seq[seq_sum] == NULL) { fprintf(stderr, "Error: Cannot allocate enough space on reading sequences. Program will exit.\n"); exit(1); }
		strncpy(seq[seq_sum], seq_file->seq.s, seq_file->seq.l);
		seq[seq_sum][seq_file->seq.l] = 0;
		first_end.emplace_back(seq_file->seq.l);
		comment[seq_sum] = (char*)malloc(sizeof(char) * (seq_file->comment.l + seq_file->name.l + 1));
		if (comment[seq_sum] == NULL) { fprintf(stderr, "Error: Cannot allocate enough space on reading sequences. Program will exit.\n"); exit(1); }
		strncpy(comment[seq_sum], seq_file->name.s, seq_file->name.l);
		strncpy(comment[seq_sum] + seq_file->name.l, seq_file->comment.s, seq_file->comment.l);
		comment[seq_sum][seq_file->comment.l + seq_file->name.l] = 0;
		seq_sum++;
	}
	kseq_destroy(seq_file);
	fclose(seq_point);
	fprintf(stderr, "done.\n");
}

int read_gaps_from_file(const char* file_name, item_int& sequence_id, item_int& ll, item_int& rr, gap_tree &this_tree)
{
	// find gaps use the sequence_id order
	FILE* file = fopen(file_name, "r");
	if (file == NULL) { fprintf(stderr, "Error: can not open file %s in reading gaps. Program will exit.\n", file_name); exit(1); }
	kseq_t* seq_file = kseq_init(fileno(file));
	int_ spaces, now_id;
	for (size_t i = 0; i < sequence_id.size(); ++i)
	{
		kseq_read(seq_file);
		spaces = 0;
		now_id = -1;
		const auto& real_id = sequence_id[i];
		for (size_t j = 0; j < seq_file->seq.l; ++j)
		{
			if (seq_file->seq.s[j] != '-')
			{
				if(spaces) this_tree.add_gaps(real_id, ll[i] + now_id, spaces);
				spaces = 0;
				++ now_id;
#if _DEBUG
				assert(is_same_char(toupper(seq_file->seq.s[j]), toupper(seq[real_id][ll[i] + now_id])));
#endif
			}
			else ++spaces;
		}
		if (spaces) this_tree.add_gaps(real_id, ll[i] + now_id, spaces);
	}
	kseq_destroy(seq_file);
	fclose(file);
	return 0;
}

int read_gaps_from_file(const char* file_name, center_int& center_info, gap_tree& this_tree)
{
	// find gaps use the sequence_id order
	FILE* file = fopen(file_name, "r");
	if (file == NULL) { fprintf(stderr, "Error: can not open file %s in reading gaps. Program will exit.\n", file_name); exit(1); }
	kseq_t* seq_file = kseq_init(fileno(file));
	int_ spaces, now_id;
	for (auto &item: center_info)
	{
		kseq_read(seq_file);
		spaces = 0;
		now_id = -1;
		const auto& real_id = item.first;
		for (size_t j = 0; j < seq_file->seq.l; ++j)
		{
			if (seq_file->seq.s[j] != '-')
			{
				if (spaces) this_tree.add_gaps(real_id, item.second.first + now_id, spaces);
				spaces = 0;
				++now_id;
#if _DEBUG
				assert(is_same_char(toupper(seq_file->seq.s[j]), toupper(seq[real_id][item.second.first + now_id])));
#endif
			}
			else ++spaces;
		}
		if (spaces) this_tree.add_gaps(real_id, item.second.first + now_id, spaces);
	}
	kseq_destroy(seq_file);
	fclose(file);
	return 0;
}

int read_gaps_from_merged_profile(const char* file_name, item_int& sequence_id1, item_int& sequence_id2, item_int& ll, item_int& rr, map_int& mp, gap_tree &this_tree)
{
	// find gaps use the sequence_id order
	FILE* file = fopen(file_name, "r");
	if (file == NULL) { fprintf(stderr, "Error: can not open file %s in reading gaps. Program will exit.\n", file_name); exit(1); }
	kseq_t* seq_file = kseq_init(fileno(file));
	int_ spaces, now_id;
	for (size_t i = 0; i < sequence_id1.size(); ++i)
	{
		kseq_read(seq_file);
		spaces = 0;
		now_id = -1;
		const auto& real_id = sequence_id1[i], &map_place = mp[real_id];
		const auto& start_place = ll[map_place], &end_place = rr[map_place];
		for (size_t j = 0; j < seq_file->seq.l; ++j)
		{
			if (seq_file->seq.s[j] != '-')
			{
				this_tree.set_gaps(real_id, start_place + now_id, spaces);
				spaces = 0;
				++ now_id;
#if _DEBUG
				assert(is_same_char(toupper(seq_file->seq.s[j]), toupper(seq[real_id][start_place + now_id])));
#endif
			}
			else ++spaces;
		}
		if (spaces) 
		{
			// must be rr[mp[i]] - 1
			this_tree.set_gaps(real_id, start_place + now_id, spaces);
		}
	}
	for (size_t i = 0; i < sequence_id2.size(); ++i)
	{
		kseq_read(seq_file);
		spaces = 0;
		now_id = -1;
		const auto& real_id = sequence_id2[i], &map_place = mp[real_id];
		const auto& start_place = ll[map_place], &end_place = rr[map_place];
		for (size_t j = 0; j < seq_file->seq.l; ++j)
		{
			if (seq_file->seq.s[j] != '-')
			{
				this_tree.set_gaps(real_id, start_place + now_id, spaces);
				spaces = 0;
				++ now_id;
#if _DEBUG
				assert(is_same_char(toupper(seq_file->seq.s[j]), toupper(seq[real_id][start_place + now_id])));
#endif
			}
			else ++spaces;
		}
		if (spaces) 
		{
			// must be rr[mp[i]] - 1
			this_tree.set_gaps(real_id, start_place + now_id, spaces);
		}
	}
	kseq_destroy(seq_file);
	fclose(file);
	return 0;
}

int write_gaps_to_file(const char* file_name, center_int &center_info, bool first, gap_tree &this_tree)
{
	// find gaps use the sequence_id order
	FILE* file = fopen(file_name, "w");
	if (file == NULL) { fprintf(stderr, "Error: can not open file %s in writing gaps. Program will exit.\n", file_name); exit(1); }
	int_ spaces;
	for (auto &block: center_info)
	{
		const auto& real_id = block.first;
		const auto &[lli, rri] = block.second;
		fputc('>', file);
		fputs(comment[real_id], file);
		fputc('\n', file);
		auto& this_seq_gap = this_tree.get_gaps(real_id);
		if (first)
		{
			spaces = this_seq_gap.get_gaps(lli - 1);
			for (int_ space = 0; space < spaces; ++space) fputc('-', file);
		}
		for (int_ j = lli; j < rri; ++j)
		{
			spaces = this_seq_gap.get_gaps(j);
			fputc(seq[real_id][j], file);
#if _DEBUG
			assert(spaces >= 0);
#endif
			for (int_ space = 0; space < spaces; ++space) fputc('-', file);
		}
		fputc('\n', file);
	}
	fclose(file);
	return 0;
}

int write_gaps_to_file(const char* file_name, item_int& sequence_id, item_int& ll, item_int& rr, bool first, gap_tree &this_tree)
{
	// find gaps use the sequence_id order
	FILE* file = fopen(file_name, "w");
	if (file == NULL) { fprintf(stderr, "Error: can not open file %s in writing gaps. Program will exit.\n", file_name); exit(1); }
	int_ spaces;
	for (size_t i = 0; i < sequence_id.size(); ++i)
	{
		auto &real_id = sequence_id[i], &lli = ll[i], &rri = rr[i];
		fputc('>', file);
		fputs(comment[real_id], file);
		fputc('\n', file);
		auto &this_seq_gap = this_tree.get_gaps(real_id);
		if (first)
		{
			spaces = this_seq_gap.get_gaps(lli - 1);
			for (int_ space = 0; space < spaces; ++space) fputc('-', file);
		}
		for (int_ j = lli; j < rri; ++j)
		{
			spaces = this_seq_gap.get_gaps(j);
			fputc(seq[real_id][j], file);
#if _DEBUG
			assert(spaces >= 0);
#endif
			for (int_ space = 0; space < spaces; ++space) fputc('-', file);
		}
		fputc('\n', file);
	}
	fclose(file);
	return 0;
}

int write_gaps_to_file_with_lr(const char* file_name, item_int& sequence_id, item_int& ll, item_int& rr, item_int& ll_border, item_int& rr_border, gap_tree& this_tree)
{
	// find gaps use the sequence_id order
	FILE* file = fopen(file_name, "w");
	if (file == NULL) { fprintf(stderr, "Error: can not open file %s in writing gaps. Program will exit.\n", file_name); exit(1); }
	int_ spaces;
	for (size_t i = 0; i < sequence_id.size(); ++i)
	{
		auto& real_id = sequence_id[i], & lli = ll[i], & rri = rr[i], & llbi = ll_border[i], & rrbi = rr_border[i];
		auto& this_seq_gap = this_tree.get_gaps(real_id);
		fprintf(file, "> %d %s %d %d left_border = %d (with %d gaps), right_border = %d (with %d gaps), len_sum = %d\n", 
			    real_id, comment[real_id], lli, rri, llbi, this_seq_gap.get_gaps(llbi - 1), rrbi, this_seq_gap.get_gaps(rrbi - 1),
			    this_seq_gap.gap_sum() + rri - lli);
		//continue;
		spaces = this_seq_gap.get_gaps(lli - 1);
		for (int_ space = 0; space < spaces; ++space) fputc('-', file);
		for (int_ j = lli; j < rri; ++j)
		{
			spaces = this_seq_gap.get_gaps(j);
			fputc(seq[real_id][j], file);
#if _DEBUG
			assert(spaces >= 0);
#endif
			for (int_ space = 0; space < spaces; ++space) fputc('-', file);
		}
		fputc('\n', file);
	}
	fclose(file);
	return 0;
}


int write_gaps_to_file_in_dispatch(const char *file_name, item_int& sequence_id, item_int& ll, item_int& rr, gap_tree &this_tree)
{
	// find gaps use the sequence_id order
	FILE* file = fopen(file_name, "w");
	if (file == NULL) { fprintf(stderr, "Error: can not open file %s in writing gaps. Program will exit.\n", file_name); exit(1); }
	int_ spaces;
	for (size_t i = 0; i < sequence_id.size(); ++i)
	{
		auto &real_id = sequence_id[i], &lli = ll[i], &rri = rr[i];
		write_seq(real_id, lli, rri, this_tree.get_gaps(real_id), file);
	}
	fclose(file);
	return 0;
}

int write_raw_to_file(const char *file_name, item_int& sequence_id, item_int& ll, item_int& rr)
{
	FILE* file = fopen(file_name, "w");
	if (file == NULL) { fprintf(stderr, "Error: can not open file %s in writing gaps. Program will exit.\n", file_name); exit(1); }
	for (size_t i = 0; i < sequence_id.size(); ++i)
	{
		auto &seq_id = sequence_id[i];
		fputc('>', file);
		fputs(comment[seq_id], file);
#if _DEBUG
		assert(ll[i] <= rr[i]);
#endif
		fprintf(file, " %d %d %d", sequence_id[i], ll[i], rr[i]);
		fputc('\n', file);
		for (int_ j = ll[i]; j < rr[i]; ++j) fputc(seq[seq_id][j], file);
		fputc('\n', file);
	}
	fclose(file);
	return 0;
}

int write_gaps_to_file_center(const char *file_name, center_int &center_block, gap_tree &this_tree)
{
	FILE* file = fopen(file_name, "w");
	if (file == NULL) { fprintf(stderr, "Error: can not open file %s in writing gaps. Program will exit.\n", file_name); exit(1); }
	int_ spaces;
	for (auto& item : center_block)
	{
		auto& [real_id, ____] = item;
		int_ real_id_ = real_id;
		auto &[ll, rr] = ____;
		write_seq(real_id_, ll, rr, this_tree.get_gaps(real_id), file);
	}
	fclose(file);
	return 0;
}

int gap_system_init(const int &sequence_total)
{
	final_gap_stores = new gaps[sequence_total];
	return 0;
}

int gap_system_free()
{
	delete[] final_gap_stores;
	final_gap_stores = nullptr;
	return 0;
}

int gap_system_reinit(item_int& v, int size)
{
	if (final_gap_stores == nullptr) gap_system_init(size);
	for (size_t i = 0; i < size; ++i) final_gap_stores[i].deinit(v[i] + 1);
	return 0;
}

void write_seq(int_ &seq_id, int_ &ll, int_ &rr, gap_system &gap_load, FILE *f)
{
#if _DEBUG
    fprintf(f, "> %d %s %d %d\n", seq_id, comment[seq_id], ll, rr);
#else
    fprintf(f, "> %s\n", comment[seq_id]);
#endif
	int_ spaces = gap_load.get_gaps(ll - 1);
#if _DEBUG
	assert(spaces >= 0);
#endif
	for (int_ space = 0; space < spaces; ++space) fputc('-', f);
	for (int_ j = ll; j < rr; ++j)
	{
		spaces = gap_load.get_gaps(j);
		fputc(seq[seq_id][j], f);
#if _DEBUG
		assert(spaces >= 0);
#endif
		for (int_ space = 0; space < spaces; ++space) fputc('-', f);
	}
	fputc('\n', f);
#if _DEBUG
	fflush(f);
#endif
}

inline void write_preload_seq(int_ &seq_id, int_ &ll, int_ &rr, int_ &border_l, int_ &border_r, pair_int &gap_border, gap_system &gap_load, FILE *f)
{
	auto [pre_ll, pre_rr] = gap_border;
	if (ll != border_l) pre_ll = 0;
	if (rr != border_r) pre_rr = 0;
#if _DEBUG
    fprintf(f, "> %d %s %d %d\n", seq_id, comment[seq_id], ll, rr);
#else
    fprintf(f, "> %s\n", comment[seq_id]);
#endif
	int_ spaces = gap_load.get_gaps(ll - 1) - pre_ll;
#if _DEBUG
	assert(spaces >= 0);
#endif
	for (int_ space = 0; space < spaces; ++space) fputc('-', f);
	for (int_ j = ll; j < rr - 1; ++j)
	{
		spaces = gap_load.get_gaps(j);
		fputc(seq[seq_id][j], f);
#if _DEBUG
		assert(spaces >= 0);
#endif
		for (int_ space = 0; space < spaces; ++space) fputc('-', f);
	}
	if (ll != rr) 
	{
		fputc(seq[seq_id][rr - 1], f);
		spaces = gap_load.get_gaps(rr - 1) - pre_rr;
#if _DEBUG
		assert(spaces >= 0);
#endif
		for (int_ space = 0; space < spaces; ++space) fputc('-', f);
	}
	fputc('\n', f);
#if _DEBUG
	fflush(f);
#endif
}