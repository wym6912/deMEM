#ifndef __EXT_PROG_CONFIG__
#define __EXT_PROG_CONFIG__

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>

#include "../include/json.hpp"

#if defined(__linux__)
#include <linux/limits.h>
#include <unistd.h>
#endif

typedef struct aligner_config
{
    bool in_arg, out_to_stdout, single_thread, has_two_profile_prefix, change_penalty;
    char* prog_path;
#if( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__) ) )
    char* prog_name;
#endif
    char* in_arg_str;
    char* out_arg_str;
    char* DNA, *Protein;
    char* Thread_first;
    char* other_args;
    char* profile1_prefix, *profile2_prefix;
    char* gapopen, *gapext;
} aligner_config;

void aligner_method_destroy();
void print_method(int method);
void print_mafft_method(int mafft_method);
void print_wmsa_align_method(int wmsa_method);
void print_mafft_merge_method(int mafft_method);
void print_wmsa_merge_method(int wmsa_method);
void print_abpoa_cluster_method(int abpoa_method);
char* expand_prog_path(const char *ext_program_path);

using json = nlohmann::json;
void json_read(const char* file_name, json& config_json);
void print_json(std::string file_name, const json& config_json);
void print_json_sample(const char* file_name, const json& config_json);
void print_program_method(const json& config_json);
void aligner_method_init(const json& config_json);

#endif