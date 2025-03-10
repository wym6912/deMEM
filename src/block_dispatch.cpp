#include <cstdio>
#include <cstdlib>
#include <string>
#include <chrono>

#if (defined(__linux__) || defined(__APPLE__))
#include <getopt.h>
#define _atoi64(val) strtoll(val, NULL, 10)
#else
#include "include/getopt9/include/getopt.h"
#include "include/unistd.h"
#endif

#include "lib/parblockaligner.hpp"
#include "lib/libblockaligner.hpp"
#include "lib/version.h"
#include "lib/libfile_v2.hpp"
#include "lib/ext_prog_config.hpp"

char* c_input = NULL, * c_output = NULL, * c_seq = NULL, * thread_str = NULL;
char** seq = NULL, **comment = NULL;
item_int first_end, this_beg, this_end, real_size, not_in_any_cluster;
int_ seq_sum, cluster_sum, threshold;
int not_in_any_cluster_sz, block_sequence_type, threads_num;

json config_json;

std::vector <cluster_irreg> blockclusters;
cluster_irreg only_one_cluster;

std::string path_prefix, this_path, this_moved_path, tmp_path_prefix, cluster_path_prefix;
std::vector <std::string> final_files;
std::vector <item_int> seq_ids;
std::mt19937 rand_gen;

map_int cluster_map;

gap_system* final_gap_stores;
gap_tree gap_main;

void version()
{
    fprintf(stderr, "blockaligner version %d.%d.%d.%d\n", VER_MAJOR, VER_MINOR, VER_RELEASE_BLOCK, VER_BUILD);
    exit(0);
}

void usage(char name[])
{
    fprintf(stderr, "\n\tUsage: %s [options]\n\n", name);
    fprintf(stderr, "Make alignment by a user specified aligner\n");
    fprintf(stderr, "Available options:\n");
    fprintf(stderr, "\t--in        FILE      input file name (Required)\n");
    fprintf(stderr, "\t--seq       FILE      sequence file name (Required)\n");
    fprintf(stderr, "\t--out       FILE      output file name (Required)\n");
    fprintf(stderr, "\t--json      FILE      JSON config file name (Required)\n");
    fprintf(stderr, "\t--help                print help message\n");
    fprintf(stderr, "\t--version             show program version\n");
    fprintf(stderr, "Example:\n\t%s --in seq.mem.unq --seq seq.fasta --out seq_out.fasta --json settings.json\n", name);
}


void get_args(int argc, char* argv[])
{
    extern char* optarg;
    extern int optind;
    int c;

    config_json["thread"] = 1;
    config_json["path"] = "swap/";
    config_json["high_sim"] = false;
    config_json["threshold_block"] = threshold = 10;

#if ( defined(_WIN32) || defined(_WIN64) )
    std::string consts[] = { "in", "out", "help", "version", "seq", "json" };
#endif

    while (1)
    {
        //int this_option_optind = optind ? optind : 1;
        int option_index = 0;
#if ( defined(_WIN32) || defined(_WIN64) )
        static struct option long_options[] =
        {
            {(char*)consts[0].c_str(),        required_argument, 0, 'f'},
            {(char*)consts[1].c_str(),        required_argument, 0, 'o'},
            {(char*)consts[2].c_str(),        no_argument,       0, 'h'},
            {(char*)consts[3].c_str(),        no_argument,       0, 'v'},
            {(char*)consts[4].c_str(),        required_argument, 0, 's'},
            {(char*)consts[5].c_str(),        required_argument, 0, 'j'},
            {0,                               0,                 0,  0 }
        };
#else
        static struct option long_options[] =
        {
            {"in",      required_argument, 0, 'f'},
            {"out",     required_argument, 0, 'o'},
            {"help",    no_argument,       0, 'h'},
            {"version", no_argument,       0, 'v'},
            {"seq",     required_argument, 0, 's'},
            {"json",    required_argument, 0, 'j'},
            {0,         0,                 0,  0 }
        };
#endif
        c = getopt_long(argc, argv, "", long_options, &option_index);
        if (c == -1) break;
        switch (c)
        {
        case 0:
            fprintf(stderr, "Not supported: option %s", long_options[option_index].name);
            if (optarg) fprintf(stderr, " with arg %s", optarg);
            fprintf(stderr, "\n");
            break;
        case 'h':
            usage(argv[0]);
            version();
            break;
        case 'f':
            c_input = argv[optind - 1];
            break;
        case 'o':
            c_output = argv[optind - 1];
            break;
        case 's':
            c_seq = argv[optind - 1];
            break;
        case 'j':
            json_read(argv[optind - 1], config_json);
            fprintf(stderr, "JSON = %s\n", config_json.dump().c_str());
            break;
        case 'v':
            version();
        }
    }

    if (c_input == NULL || c_output == NULL || c_seq == NULL || config_json.empty())
    {
        fprintf(stderr, "ERROR: Please determine input/output/sequence/json file. use %s --help for more information.\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    if(config_json["all_mafft"] == true)
    {
        config_json["align_method"] = 1;   if(config_json.find("mafft_align_method") == config_json.end())   config_json["mafft_align_method"] = 3;
        config_json["merge_method"] = 1;   if(config_json.find("mafft_merge_method") == config_json.end())   config_json["mafft_merge_method"] = 0;
        config_json["cluster_method"] = 1; if(config_json.find("mafft_cluster_method") == config_json.end()) config_json["mafft_cluster_method"] = 3;
    }
    if(config_json["all_wmsa"] == true)
    {
        config_json["align_method"] = 0;   if(config_json.find("wmsa_align_method") == config_json.end())   config_json["wmsa_align_method"] = 0;
        config_json["merge_method"] = 0;   if(config_json.find("wmsa_merge_method") == config_json.end())   config_json["wmsa_merge_method"] = 2;
        config_json["cluster_method"] = 0; if(config_json.find("wmsa_cluster_method") == config_json.end()) config_json["wmsa_cluster_method"] = 2;
    }
    if(config_json["all_abpoa"] == true)
    {
        config_json["align_method"] = 3;
        config_json["merge_method"] = 3;
        config_json["cluster_method"] = 3; if(config_json.find("abpoa_cluster_method") == config_json.end()) config_json["abpoa_cluster_method"] = 0;
    }

    threshold = config_json["threshold_block"];
    block_sequence_type = config_json["seq_type"];
    path_prefix = config_json["path"];
    thread_str = strdup(std::to_string((int)config_json["thread"]).c_str());
    threads_num = (int)config_json["thread"];
    fprintf(stderr, "Use %s threads\n", thread_str);
    print_program_method(config_json);
}

void read_clusters()
{
    fprintf(stderr, "Reading clusters...");
    FILE* file = fopen(c_input, "r");
    if (file == NULL) { fprintf(stderr, "Error: can not open input file %s. Please check this file!\n", c_input); exit(1); }
    cluster_irreg* this_cluster;
	std::vector <sequence_irreg> sequences;
	sequence_irreg sequence;
    block_irreg* this_mem;
    bool *seq_in_cluster = new bool[seq_sum + 1]{}; not_in_any_cluster_sz = seq_sum;
	int block_length, block_size; bool end = false;
	do
	{
		this_cluster = new cluster_irreg;
		do
		{		
			if(fscanf(file, "\n%d:%d", &block_length, &block_size) != 2)
			{
				end = true;
				break;
			}
			for(int i = 0; i < block_size; ++ i)
			{
				fscanf_noreturn(file, 2, "%d %d", &std::get<0>(sequence), &std::get<1>(sequence));
                std::get<2>(sequence) = 0;
				sequences.emplace_back(std::move(sequence));
			}
            if(block_size)
            {
                this_mem = new block_irreg;
                this_mem->container = std::move(sequences);
                this_mem->length = block_length;
                this_mem->center = find_center(this_mem->container);
                this_mem->calc_area();
                this_cluster->sub_blocks.emplace_back(std::move(*this_mem));
                // fprintf(stderr, "Reading block with this_mem = %p and block_length = %d, block_size = %d\n", this_mem, block_length, block_size);
                delete this_mem;
                sequences.clear();
            }
		} while (block_length && block_size);
        if(! end)
        {
		    // sequence tag
		    fscanf_noreturn(file, 1, "%d", &block_size);
		    this_cluster->sequences.resize(block_size);
		    for(int i = 0; i < block_size; ++ i) 
            {
                fscanf_noreturn(file, 1, "%d", &(this_cluster->sequences[i]));
                seq_in_cluster[this_cluster->sequences[i]] = true; 
                -- not_in_any_cluster_sz;
            }
		    blockclusters.emplace_back(std::move(*this_cluster));
            fprintf(stderr, "Reading seq...done.\n");
        }
        delete this_cluster;
	} while (! end);
    fprintf(stderr, "done.\n");
    // TODO: parallelize
    for(auto &cluster: blockclusters)
    {
        sort(cluster.sub_blocks.begin(), cluster.sub_blocks.end());
        cluster.common_block = new block_irreg;
        *cluster.common_block = std::move(cluster.sub_blocks.back());
        cluster.sub_blocks.pop_back();
#if 0
        fprintf(stderr, "cluster has %llu sequences\n", cluster.sequences.size());
        if(cluster.sub_blocks.size())
        {
            fprintf(stderr, "%lld %lld\n", cluster.common_block->area, cluster.sub_blocks.back().area);
            fprintf(stderr, "center block has %llu sequences and %llu sub-blocks\n", cluster.common_block->container.size(), cluster.sub_blocks.size());
        }
        else
        {
            fprintf(stderr, "%lld\n", cluster.common_block->area);
        }
#endif
    }
    fprintf(stderr, "done. Now calculating sequence ids not in any cluster...\n");
    // calcluate the block not used
    not_in_any_cluster.clear();
    not_in_any_cluster.reserve(not_in_any_cluster_sz);
    for(int i = 0; i < seq_sum; ++ i)
    {
        if(! seq_in_cluster[i]) not_in_any_cluster.emplace_back(i);
    }
    if(not_in_any_cluster.size() != not_in_any_cluster_sz)
    {
        fprintf(stderr, "Error: cluster has redunant sequences. Program will exit. size = %llu, expected %llu, found %d sequences\nID =", 
                not_in_any_cluster.size(), not_in_any_cluster_sz, seq_sum);
        for(auto &seqid: not_in_any_cluster) fprintf(stderr, " %d", seqid);
        fputc('\n', stderr);
        exit(EXIT_FAILURE);
    }
    fclose(file);
}

void freeall()
{
    for (int_ i = 0; i < seq_sum; ++i) free(seq[i]);
    for (int i = 0; i < blockclusters.size(); ++i)
    {
        if (blockclusters[i].common_block) delete blockclusters[i].common_block;
    }
}

int main(int argc, char *argv[])
{
    get_args(argc, argv);
    aligner_method_init(config_json);
    read_seq();
    read_clusters();
#if ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__)) )
    std::string cwd_prefix;
    char* cwd = getcurrdir();
    cwd_prefix = cwd;
    free(cwd);
    path_prefix = cwd_prefix + "/" + path_prefix;
    cwd_prefix.clear();
#endif
    tmp_path_prefix = path_prefix + "tmp/";
    cluster_path_prefix = path_prefix + "clstr/";
    if (makedir(path_prefix.c_str()))
    {
        fprintf(stderr, "Warning: can not make dir %s. Program will rewrite any files in this folder.\n", path_prefix.c_str());
        // fprintf(stderr, "Error: can not make dir %s. Program will exit.\n", path_prefix.c_str()); exit(1);
    }
    if (makedir(tmp_path_prefix.c_str()) || makedir(cluster_path_prefix.c_str()))
    {
        fprintf(stderr, "Warning: can not make dir %s or %s. Program will rewrite any files in these folders.\n",
                         tmp_path_prefix.c_str(), cluster_path_prefix.c_str());
    }
    final_files.reserve(blockclusters.size());
    fprintf(stderr, "Making alignment for clusters...");
    rand_gen.seed(std::random_device{}());
    std::chrono::time_point<std::chrono::system_clock> time_start = std::chrono::system_clock::now(), time_end;
    if(config_json.find("highsim") != config_json.end() && ! config_json["highsim"])
    {
        fprintf(stderr, "Treat as clusters.\n");
        for (int i = 0; i < blockclusters.size(); ++ i)
        {
            this_beg.clear(); this_beg.assign(blockclusters[i].sequences.size(), 0);
            this_end.clear(); this_end.resize(this_beg.size());
            cluster_sum += this_beg.size();
            cluster_map.clear();
            fprintf(stderr, "Cluster %d has %llu sequences\n", i, this_beg.size());
            for (size_t j = 0; j < this_end.size(); ++j) this_end[j] = first_end[blockclusters[i].sequences[j]];
            // make clustermap
            for (size_t j = 0; j < blockclusters[i].sequences.size(); ++ j) cluster_map.emplace(blockclusters[i].sequences[j], j);
            // reinit gap system
            gap_system_reinit(this_end, (const int_)this_beg.size());
#if 1
            this_path = decompose_without_SSW(blockclusters[i].sequences, *(blockclusters[i].common_block), blockclusters[i].sub_blocks, this_beg, this_end, cluster_map, gap_main, threads_num);
#else
            this_path = decompose_without_SSW(blockclusters[i].sequences, *(blockclusters[i].common_block), blockclusters[i].sub_blocks, this_beg, this_end, 0, 0, cluster_map, gap_main);
#endif
            this_moved_path = cluster_path_prefix + "cluster_" + std::to_string(i) + ".fasta";
            write_gaps_to_file(this_moved_path.c_str(), blockclusters[i].sequences, this_beg, this_end, true, gap_main);
            // Filemove(this_path.c_str(), this_moved_path.c_str());
            final_files.emplace_back(std::move(this_moved_path));
            seq_ids.emplace_back(std::move(blockclusters[i].sequences));
            gap_system_free();
            gap_main.this_tree.clear();
#if (defined(_WIN32) || defined(_WIN64))
            removedir((tmp_path_prefix + "\0\0").c_str());
#else
            removedir(tmp_path_prefix.c_str());
#endif
            makedir(tmp_path_prefix.c_str());
        }
        fprintf(stderr, "Done.\n");
        fprintf(stderr, "Making profile-profile alignment...\n");
        if (final_files.size() == 0)
        {
            fprintf(stderr, "Warning: not found any clusters. There is no need to make profile-profile alignment.\n");
        }
        else if(final_files.size() == 1)
        {
            fprintf(stderr, "Info: only one cluster found. There is no need to make profile-profile alignment.\n");
            Filemove(final_files[0].c_str(), c_output);
        }
        else if(final_files.size() == 2)
        {
            fprintf(stderr, "Info: found two clusters. Use two profiles alignment method.\n");
            TwoProfile_dispatch(final_files[0], final_files[1]);
            Filemove(final_files[0].c_str(), c_output);
        }
        else
        {
            std::string center_file_name = tmp_path_prefix + "center.fasta";
            std::string file_name = tmp_path_prefix + "file_list.txt";
            if(Profiles_dispatch(seq_ids, final_files, c_output, center_file_name, file_name))
            {
                fprintf(stderr, "Error: Can not execute profile-profile alignment. Program will exit.\n");
                exit(1);
            }
        }
        fprintf(stderr, "Done.\n");
        if(not_in_any_cluster_sz)
        {
            fprintf(stderr, "Warning: clusters are not contained all sequences. Making alignment for these sequences\n");
            std::string frag_name = path_prefix + "frag.fasta", out_name_ = c_output;
            Frag_dispatch_full_seq(frag_name.c_str(), not_in_any_cluster);
            if(not_in_any_cluster_sz != seq_sum) TwoProfile_dispatch(out_name_, frag_name);
            else if(Filemove(frag_name.c_str(), out_name_.c_str())) { fprintf(stderr, "Error: can not move final file. Program will exit.\n"); exit(1); }
        }
    }
    else
    {
        fprintf(stderr, "Treat as one cluster.\n");
        only_one_cluster.sequences.reserve(seq_sum);
        only_one_cluster.common_block = nullptr;
        for(auto &cluster_: blockclusters)
        {
            only_one_cluster.sequences.insert(only_one_cluster.sequences.end(), cluster_.sequences.begin(), cluster_.sequences.end());
            if (cluster_.common_block) only_one_cluster.sub_blocks.push_back(*cluster_.common_block);
            cluster_.common_block->freecontainer();
            cluster_.common_block = nullptr;
            only_one_cluster.sub_blocks.insert(only_one_cluster.sub_blocks.end(), cluster_.sub_blocks.begin(), cluster_.sub_blocks.end());
        }
        only_one_cluster.sequences.insert(only_one_cluster.sequences.end(), not_in_any_cluster.begin(), not_in_any_cluster.end());
        blockclusters.clear();
        if(only_one_cluster.sub_blocks.size())
        {
            fprintf(stderr, "Found blocks from MEMs. Aligning...\n");
            std::sort(only_one_cluster.sequences.begin(), only_one_cluster.sequences.end());
            std::sort(only_one_cluster.sub_blocks.begin(), only_one_cluster.sub_blocks.end());
            gap_system_reinit(first_end, seq_sum);
            item_int first_beg(seq_sum, 0);
            block_irreg center_block = only_one_cluster.sub_blocks.back();
            only_one_cluster.sub_blocks.pop_back();
            for(int_ i = 0; i < seq_sum; ++ i) cluster_map[i] = i;
#if 1
            this_path = decompose_with_SSW(only_one_cluster.sequences, center_block, only_one_cluster.sub_blocks, first_beg, first_end, cluster_map, gap_main, threads_num);
#else
            this_path = decompose_with_SSW(only_one_cluster.sequences, center_block, only_one_cluster.sub_blocks, first_beg, first_end, 0, 0, cluster_map, gap_main);
#endif
            write_gaps_to_file(c_output, only_one_cluster.sequences, first_beg, first_end, true, gap_main);
            fprintf(stderr, "Done.\n");
        }
        else
        {
            fprintf(stderr, "Info: no MEMs. Aligning directly\n");
            not_in_any_cluster.clear();
            for (int_ i = 0; i < seq_sum; ++i) not_in_any_cluster.emplace_back(i);
            Frag_dispatch_full_seq(c_output, not_in_any_cluster);
        }
    }
    time_end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = time_end - time_start;
#if (defined(_WIN32) || defined(_WIN64))
    if (removedir((tmp_path_prefix + "\0\0").c_str()) || removedir((cluster_path_prefix + "\0\0").c_str()))
#else
    if (removedir(tmp_path_prefix.c_str()) || removedir(cluster_path_prefix.c_str()))
#endif
        fprintf(stderr, "Error: can not remove dir %s or %s. Please remove it manually.\n", tmp_path_prefix.c_str(), cluster_path_prefix.c_str());
#if (defined(_WIN32) || defined(_WIN64))
    if (removedir((path_prefix + "\0\0").c_str()))
#else
    if (removedir(path_prefix.c_str()))
#endif
        fprintf(stderr, "Warning: folder %s can not be removed. Please remove it manually.\n", path_prefix.c_str());
    freeall();
    aligner_method_destroy();
    fprintf(stderr, "Block dispatch time elapsed: %lfs\n", elapsed.count());
    return 0;
}