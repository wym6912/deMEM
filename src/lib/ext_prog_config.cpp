#include "ext_prog_config.hpp"

aligner_config align_seqs_config_, align_two_profile_config_, align_profiles_config_;

char* expand_prog_path(const char *ext_program_path)
{
    std::string program_path;
#if (defined(_WIN32) || defined(_WIN64))
    program_path = _pgmptr;
    program_path.erase(program_path.find_last_of('\\') + 1, program_path.size());
#else
    char progpath[PATH_MAX];
    readlink("/proc/self/exe", progpath, PATH_MAX);
    program_path = progpath;
    program_path.erase(program_path.find_last_of('/') + 1, program_path.size());
#endif
    char *return_path = nullptr;
    program_path += ext_program_path;
    return_path = strdup(program_path.c_str());
    return return_path;
}

char* mafft_method_determine(int mafft_method)
{
    switch (mafft_method)
    {
    case 0:
        return "--nofft --retree 1";
    case 1:
        return "--nofft --retree 2";
    case 2:
        return "--nofft --maxiterate 1000";
    case 3:
        return "--fft --retree 1";
    case 4:
        return "--fft --retree 2";
    case 5:
        return "--fft --maxiterate 1000";
    case 6:
        return "--globalpair --maxiterate 1000";
    case 7:
        return "--localpair --maxiterate 1000";
    }
}

char* wmsa_method_merge_method_determine(int wmsa_method)
{
    switch (wmsa_method)
    {
    case 0:
        return "-F -N";
    case 1:
        return "-G -N";
    case 2:
        return "-F -A";
    case 3:
        return "-G -A";
    }
}

void print_wmsa_align_method(int wmsa_method)
{
    switch (wmsa_method)
    {
    case 0:
        fprintf(stderr, " FFT+K-band\n");
        break;
    case 1:
        fprintf(stderr, " Only K-band\n");
        break;
    default:
        fprintf(stderr, " Not determined. \nError: not determined WMSA method. Program will exit.\n");
        exit(1);
    }
}

void print_mafft_method(int mafft_method)
{
    switch (mafft_method)
    {
    case 0:
        fprintf(stderr, "NW-NS-1\n");
        break;
    case 1:
        fprintf(stderr, "NW-NS-2\n");
        break;
    case 2:
        fprintf(stderr, "NW-NS-i\n");
        break;
    case 3:
        fprintf(stderr, "FFT-NS-1\n");
        break;
    case 4:
        fprintf(stderr, "FFT-NS-2\n");
        break;
    case 5:
        fprintf(stderr, "FFT-NS-i\n");
        break;
    case 6:
        fprintf(stderr, "G-INS-i\n");
        break;
    case 7:
        fprintf(stderr, "L-INS-i\n");
        break;
    default:
        fprintf(stderr, "Not determined. \nError: not determined MAFFT method. Program will exit.\n");
        exit(1);
    }
}

void print_method(int method)
{
    switch(method)
    {
    case 0:
        fprintf(stderr, " WMSA\n");
        break;
    case 1:
        fprintf(stderr, " MAFFT\n");
        break;
    case 2:
        fprintf(stderr, " HAlign 3\n");
        break;
    case 3:
        fprintf(stderr, " abPOA\n");
        break;
    default:
        fprintf(stderr, " Not determined. \nError: not determined method. Please check your argument and run again. Program will exit.\n");
        exit(1);
    }
}

void print_wmsa_merge_method(int wmsa_method)
{
    switch (wmsa_method)
    {
    case 0:
        fprintf(stderr, " FFT+K-band with simple gap mode\n");
        break;
    case 1:
        fprintf(stderr, " K-band with simple gap mode\n");
        break;
    case 2:
        fprintf(stderr, " FFT+K-band with gap open and extension mode\n");
        break;
    case 3:
        fprintf(stderr, " K-band with gap open and extension mode\n");
        break;
    default:
        fprintf(stderr, " Not determined. \nError: not determined WMSA method. Program will exit.\n");
        exit(1);
    }
}

void print_abpoa_cluster_method(int abpoa_method)
{
    switch (abpoa_method)
    {
    case 0:
        fprintf(stderr, " use center sequences make pre-align\n");
        break;
    case 1:
        fprintf(stderr, " directly align clusters\n");
        break;
    default:
        fprintf(stderr, " Not determined. \nError: not determined abPOA cluster method. Program will exit.\n");
        exit(1);
    }
}

void print_mafft_merge_method(int mafft_merge_method)
{
    switch (mafft_merge_method)
    {
    case 0:
        fprintf(stderr, " NW\n");
        break;
    case 1:
        fprintf(stderr, " FFT\n");
        break;
    default:
        fprintf(stderr, "Not determined. \nError: not determined MAFFT method. Program will exit.\n");
        exit(1);
    }
}

void print_program_method(const json& config_json)
{
    fprintf(stderr, "Alignment method:"); print_method(config_json["align_method"]);
    if(config_json["align_method"] == 0) { fprintf(stderr, "Alignment WMSA method:"); print_wmsa_align_method(config_json["wmsa_align_method"]);}
    else if(config_json["align_method"] == 1) { fprintf(stderr, "Alignment MAFFT method:"); print_mafft_method(config_json["mafft_align_method"]); }
    fprintf(stderr, "Merge method:"); print_method(config_json["merge_method"]);
    if(config_json["merge_method"] == 0) { fprintf(stderr, "Merge WMSA method:"); print_wmsa_merge_method(config_json["wmsa_merge_method"]); }
    else if(config_json["merge_method"] == 1) { fprintf(stderr, "Merge MAFFT method:"); print_mafft_merge_method(config_json["mafft_merge_method"]); }
    fprintf(stderr, "Cluster merge method:"); print_method(config_json["cluster_method"]);
    if(config_json["cluster_method"] == 0) { fprintf(stderr, "Cluster merge WMSA method:"); print_wmsa_merge_method(config_json["wmsa_cluster_method"]); }
    else if(config_json["cluster_method"] == 1) { fprintf(stderr, "Cluster merge MAFFT method:"); print_mafft_method(config_json["mafft_cluster_method"]); }
}

void json_read(const char* file, json& config_json)
{
    try
    {
        std::ifstream s(file);
        s >> config_json;
        s.close();
    }
    catch (...)
    {
        fprintf(stderr, "Error: can not process json file %s. Program will exit.\n", file);
        exit(1);
    }
}

void print_json(std::string file_name, const json& config_json)
{
    std::ofstream s(file_name);
    s << config_json;
    s.flush();
    s.close();
}

void print_json_sample(const char* file_name, const json& config_json)
{
    std::ofstream s(file_name);
    s << std::setw(4) << config_json;
    s.flush();
    s.close();
}

void HAlign_seqs_config()
{
#ifdef _DEBUG
#if (defined(_WIN32) || defined(_WIN64))
    align_seqs_config_.prog_path = expand_prog_path("..\\..\\src\\external\\subhalign\\bin\\halign.bat");
#else
    align_seqs_config_.prog_path = expand_prog_path("/src/external/subhalign/bin/halign");
#endif
#else
#if (defined(_WIN32) || defined(_WIN64))
    align_seqs_config_.prog_path = expand_prog_path("\\halign.bat");
#else
    align_seqs_config_.prog_path = "halign";
#endif
#endif
#if ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__) ) )
    align_seqs_config_.prog_name = expand_prog_path("/halign");
#endif
    align_seqs_config_.in_arg = false; align_seqs_config_.in_arg_str = nullptr;
    align_seqs_config_.DNA = align_seqs_config_.Protein = nullptr;
    align_seqs_config_.gapopen = align_seqs_config_.gapext = nullptr;
    align_seqs_config_.change_penalty = false;
    align_seqs_config_.out_to_stdout = false; align_seqs_config_.out_arg_str = "-o";
    align_seqs_config_.single_thread = false; align_seqs_config_.Thread_first = "-t";
    align_seqs_config_.other_args = nullptr;
}

void MAFFT_seqs_config(const json &config_json)
{
#ifdef _DEBUG
#if (defined(_WIN32) || defined(_WIN64))
    align_seqs_config_.prog_path = expand_prog_path("..\\..\\src\\external\\submafft\\windows\\mafft.bat");
#else
    align_seqs_config_.prog_path = expand_prog_path("/src/external/submafft/linux/bin/mafft");
#endif
#else
#if (defined(_WIN32) || defined(_WIN64))
    align_seqs_config_.prog_path = expand_prog_path("\\..\\mafft.bat");
#else
    align_seqs_config_.prog_path = expand_prog_path("/mafft");
#endif
#endif
#if ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__) ) )
    align_seqs_config_.prog_name = "mafft";
#endif
    align_seqs_config_.in_arg = false; align_seqs_config_.in_arg_str = nullptr;
    align_seqs_config_.DNA = "--nuc"; align_seqs_config_.Protein = "--amino";
    align_seqs_config_.out_to_stdout = true; align_seqs_config_.out_arg_str = nullptr;
    align_seqs_config_.single_thread = false; align_seqs_config_.Thread_first = "--thread";
#if 1
    align_seqs_config_.other_args = mafft_method_determine(config_json["mafft_align_method"]);
#else
    std::string mafft_method_len = mafft_method_determine(config_json["mafft_align_method"]);
    align_seqs_config_.other_args = (char*)malloc((mafft_method_len.length() + 2) * sizeof(char));
    if (align_seqs_config_.other_args == nullptr) { fprintf(stderr, "Error: can not alloc space for alignment. Program will exit.\n"); exit(1); }
    mafft_method_len.copy(align_seqs_config_.other_args, 0, mafft_method_len.length());
    align_seqs_config_.other_args[mafft_method_len.length()] = 0;
#endif
    if(config_json.find("align_gapopen") != config_json.end() || config_json.find("align_gapext") != config_json.end())
    {
        align_seqs_config_.change_penalty = true;
        if(config_json.find("align_gapopen") != config_json.end())
            align_seqs_config_.gapopen = strdup(("--op " + std::to_string((double)config_json["align_gapopen"])).c_str());
        else align_seqs_config_.gapopen = nullptr;
        if(config_json.find("align_gapext") != config_json.end())
            align_seqs_config_.gapext = strdup(("--ep " + std::to_string((double)config_json["align_gapext"])).c_str());
        else align_seqs_config_.gapext = nullptr;
    }
    else align_seqs_config_.change_penalty = false;
}

void WMSA_seqs_config(const json &config_json)
{
#ifdef _DEBUG
#if (defined(_WIN32) || defined(_WIN64))
    align_seqs_config_.prog_path = expand_prog_path("..\\Release\\fragalign.exe");
#else
    align_seqs_config_.prog_path = expand_prog_path("/fragalign");
#endif
#else
#if (defined(_WIN32) || defined(_WIN64))
    align_seqs_config_.prog_path = expand_prog_path("\\fragalign.exe");
#else
    align_seqs_config_.prog_path = expand_prog_path("/fragalign");
#endif
#endif
#if ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__) ) )
    align_seqs_config_.prog_name = "fragalign";
#endif
    align_seqs_config_.in_arg = true; align_seqs_config_.in_arg_str = "-i";
    align_seqs_config_.DNA = "-D"; align_seqs_config_.Protein = "-P";
    align_seqs_config_.out_to_stdout = true; align_seqs_config_.out_arg_str = nullptr;
    align_seqs_config_.single_thread = false; align_seqs_config_.Thread_first = "-T";
    if(config_json["wmsa_align_method"] == 0) align_seqs_config_.other_args = "-F";
    else align_seqs_config_.other_args = "-G";
    if(config_json.find("align_gapopen") != config_json.end() || config_json.find("align_gapext") != config_json.end())
    {
        align_seqs_config_.change_penalty = true;
        if(config_json.find("align_gapopen") != config_json.end())
            align_seqs_config_.gapopen = strdup(("-f " + std::to_string((double)config_json["align_gapopen"])).c_str());
        else align_seqs_config_.gapopen = nullptr;
        if(config_json.find("align_gapext") != config_json.end())
            align_seqs_config_.gapext = strdup(("-g " + std::to_string((double)config_json["align_gapext"])).c_str());
        else align_seqs_config_.gapext = nullptr;
    }
    else align_seqs_config_.change_penalty = false;
}

void abPOA_seqs_config(const json &config_json)
{
#ifdef _DEBUG
#if (defined(_WIN32) || defined(_WIN64))
    align_seqs_config_.prog_path = expand_prog_path("..\\Release\\abPOA.exe");
#else
    align_seqs_config_.prog_path = expand_prog_path("/abpoa");
#endif
#else
#if (defined(_WIN32) || defined(_WIN64))
    align_seqs_config_.prog_path = expand_prog_path("\\abPOA.exe");
#else
    align_seqs_config_.prog_path = expand_prog_path("/abpoa");
#endif
#endif
#if ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__) ) )
    align_seqs_config_.prog_name = "abPOA";
#endif
    align_seqs_config_.in_arg = false; align_seqs_config_.in_arg_str = nullptr;
    align_seqs_config_.DNA = ""; align_seqs_config_.Protein = "-c";
    align_seqs_config_.out_to_stdout = false; align_seqs_config_.out_arg_str = "-o";
    align_seqs_config_.single_thread = true;
    align_seqs_config_.other_args = "-r 1";
    if(config_json.find("align_gapopen") != config_json.end() || config_json.find("align_gapext") != config_json.end())
    {
        align_seqs_config_.change_penalty = true;
        if(config_json.find("align_gapopen") != config_json.end())
            align_seqs_config_.gapopen = strdup(("-O "+ std::to_string((double)config_json["align_gapopen"])).c_str());
        else align_seqs_config_.gapopen = nullptr;
        if(config_json.find("align_gapext") != config_json.end())
            align_seqs_config_.gapext = strdup(("-E " + std::to_string((double)config_json["align_gapext"])).c_str());
        else align_seqs_config_.gapext = nullptr;
    }
    else align_seqs_config_.change_penalty = false;
}

void MAFFT_two_profile_config(const json &config_json)
{
#ifdef _DEBUG
#if (defined(_WIN32) || defined(_WIN64))
    align_two_profile_config_.prog_path = expand_prog_path("..\\..\\src\\external\\submafft\\windows\\usr\\lib\\mafft\\mafft-profile.exe");
#else
    align_two_profile_config_.prog_path = expand_prog_path("/src/external/submafft/linux/bin/mafft-profile");
#endif
#else
#if (defined(_WIN32) || defined(_WIN64))
    align_two_profile_config_.prog_path = expand_prog_path("..\\usr\\lib\\mafft\\mafft-profile.exe");
#else
    align_two_profile_config_.prog_path = expand_prog_path("/mafft-profile");
#endif
#endif
#if ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__) ) )
    align_two_profile_config_.prog_name = "mafft-profile";
#endif
    align_two_profile_config_.in_arg = false; align_two_profile_config_.in_arg_str = nullptr;
    align_two_profile_config_.DNA = "-D"; align_two_profile_config_.Protein = "-P";
    align_two_profile_config_.out_to_stdout = true; align_two_profile_config_.out_arg_str = nullptr;
    align_two_profile_config_.single_thread = true; align_two_profile_config_.Thread_first = nullptr;
    align_two_profile_config_.has_two_profile_prefix = false;
    align_two_profile_config_.profile1_prefix = nullptr;
    align_two_profile_config_.profile2_prefix = nullptr;
    if(config_json["mafft_merge_method"] == 1) align_two_profile_config_.other_args = (char*)"-F";
    else align_two_profile_config_.other_args = nullptr;
    if(config_json.find("merge_gapopen") != config_json.end() || config_json.find("merge_gapext") != config_json.end())
    {
        align_two_profile_config_.change_penalty = true;
        if(config_json.find("merge_gapopen") != config_json.end())
            align_two_profile_config_.gapopen = strdup(("--op " + std::to_string((double)config_json["merge_gapopen"])).c_str());
        else align_two_profile_config_.gapopen = nullptr;
        if(config_json.find("merge_gapext") != config_json.end())
            align_two_profile_config_.gapext = strdup(("--ep " + std::to_string((double)config_json["merge_gapext"])).c_str());
        else align_two_profile_config_.gapext = nullptr;
    }
    else align_two_profile_config_.change_penalty = false;

}

void WMSA_two_profile_config(const json &config_json)
{
#ifdef _DEBUG
#if (defined(_WIN32) || defined(_WIN64))
    align_two_profile_config_.prog_path = expand_prog_path("\\..\\Release\\profile_two_align.exe");
#else
    align_two_profile_config_.prog_path = expand_prog_path("/profile_two_align");
#endif
#else
#if (defined(_WIN32) || defined(_WIN64))
    align_two_profile_config_.prog_path = expand_prog_path("\\profile_two_align.exe");
#else
    align_two_profile_config_.prog_path = expand_prog_path("/profile_two_align");
#endif
#endif
#if ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__) ) )
    align_two_profile_config_.prog_name = "profile_two_align";
#endif
    align_two_profile_config_.in_arg = true; align_two_profile_config_.in_arg_str = "-i";
    align_two_profile_config_.DNA = "-D"; align_two_profile_config_.Protein = "-P";
    align_two_profile_config_.out_to_stdout = false; align_two_profile_config_.out_arg_str = nullptr;
    align_two_profile_config_.single_thread = true; align_two_profile_config_.Thread_first = nullptr;
    align_two_profile_config_.has_two_profile_prefix = true;
    align_two_profile_config_.profile1_prefix = "-p";
    align_two_profile_config_.profile2_prefix = "-q";
    align_two_profile_config_.other_args = wmsa_method_merge_method_determine(config_json["wmsa_merge_method"]);
    if(config_json.find("merge_gapopen") != config_json.end() || config_json.find("merge_gapext") != config_json.end())
    {
        align_two_profile_config_.change_penalty = true;
        if(config_json.find("merge_gapopen") != config_json.end())
            align_two_profile_config_.gapopen = strdup(("-f " + std::to_string((double)config_json["merge_gapopen"])).c_str());
        else align_two_profile_config_.gapopen = nullptr;
        if(config_json.find("merge_gapext") != config_json.end())
            align_two_profile_config_.gapext = strdup(("-g " + std::to_string((double)config_json["merge_gapext"])).c_str());
        else align_two_profile_config_.gapext = nullptr;
    }
    else align_two_profile_config_.change_penalty = false;
}

void abPOA_two_profile_config(const json &config_json)
{
#ifdef _DEBUG
#if (defined(_WIN32) || defined(_WIN64))
    align_two_profile_config_.prog_path = expand_prog_path("\\..\\Release\\abPOA.exe");
#else
    align_two_profile_config_.prog_path = expand_prog_path("abpoa");
#endif
#else
#if (defined(_WIN32) || defined(_WIN64))
    align_two_profile_config_.prog_path = expand_prog_path("abPOA.exe");
#else
    align_two_profile_config_.prog_path = expand_prog_path("/abpoa");
#endif
#endif
#if ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__) ) )
    align_two_profile_config_.prog_name = "abPOA_two_profile";
#endif
    align_two_profile_config_.in_arg = false; align_two_profile_config_.in_arg_str = nullptr;
    align_two_profile_config_.DNA = ""; align_two_profile_config_.Protein = "-c";
    align_two_profile_config_.out_to_stdout = true; align_two_profile_config_.out_arg_str = nullptr;
    align_two_profile_config_.single_thread = true; align_two_profile_config_.Thread_first = nullptr;
    align_two_profile_config_.has_two_profile_prefix = true;
    align_two_profile_config_.profile1_prefix = "-i";
    align_two_profile_config_.profile2_prefix = "";
    align_two_profile_config_.other_args = "-r 1";
    if(config_json.find("merge_gapopen") != config_json.end() || config_json.find("merge_gapext") != config_json.end())
    {
        align_two_profile_config_.change_penalty = true;
        if(config_json.find("merge_gapopen") != config_json.end())
            align_two_profile_config_.gapopen = strdup(("-O " + std::to_string((double)config_json["merge_gapopen"])).c_str());
        else align_two_profile_config_.gapopen = nullptr;
        if(config_json.find("merge_gapext") != config_json.end())
            align_two_profile_config_.gapext = strdup(("-E " + std::to_string((double)config_json["merge_gapext"])).c_str());
        else align_two_profile_config_.gapext = nullptr;
    }
    else align_two_profile_config_.change_penalty = false;
}

void WMSA_profiles_config(const json &config_json)
{
#ifdef _DEBUG
#if (defined(_WIN32) || defined(_WIN64))
    align_profiles_config_.prog_path = expand_prog_path("\\..\\Release\\profilealign.exe");
#else
    align_profiles_config_.prog_path = expand_prog_path("/profilealign");
#endif
#else
#if (defined(_WIN32) || defined(_WIN64))
    align_profiles_config_.prog_path = expand_prog_path("\\profilealign.exe");
#else
    align_profiles_config_.prog_path = expand_prog_path("/profilealign");
#endif
#endif
#if ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__) ) )
    align_profiles_config_.prog_name = "profilealign";
#endif
    align_profiles_config_.in_arg = false; align_profiles_config_.in_arg_str = nullptr;
    align_profiles_config_.DNA = "-D"; align_profiles_config_.Protein = "-P";
    align_profiles_config_.out_to_stdout = true; align_profiles_config_.out_arg_str = nullptr;
    align_profiles_config_.single_thread = false; align_profiles_config_.Thread_first = "-T";
    align_profiles_config_.has_two_profile_prefix = true;
    align_profiles_config_.profile1_prefix = "-p";
    align_profiles_config_.profile2_prefix = "-i";
    align_profiles_config_.other_args = wmsa_method_merge_method_determine(config_json["wmsa_cluster_method"]);
    if(config_json.find("cluster_gapopen") != config_json.end() || config_json.find("cluster_gapext") != config_json.end())
    {
        align_profiles_config_.change_penalty = true;
        if(config_json.find("cluster_gapopen") != config_json.end())
            align_profiles_config_.gapopen = strdup(("-f " + std::to_string((double)config_json["cluster_gapopen"])).c_str());
        else align_profiles_config_.gapopen = nullptr;
        if(config_json.find("merge_gapext") != config_json.end())
            align_profiles_config_.gapext = strdup(("-g " + std::to_string((double)config_json["merge_gapext"])).c_str());
        else align_profiles_config_.gapext = nullptr;
    }
    else align_profiles_config_.change_penalty = false;
}

void MAFFT_profiles_config(const json &config_json)
{
#ifdef _DEBUG
#if (defined(_WIN32) || defined(_WIN64))
    align_profiles_config_.prog_path = expand_prog_path("..\\..\\src\\external\\submafft\\windows\\profiles_wrapper.exe");
#else
    align_profiles_config_.prog_path = expand_prog_path("/profiles_wrapper");
#endif
#else
#if (defined(_WIN32) || defined(_WIN64))
    align_profiles_config_.prog_path = expand_prog_path("..\\profiles_wrapper.exe");
#else
    align_profiles_config_.prog_path = expand_prog_path("/profiles_wrapper");
#endif
#endif
#if ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__) ) )
    align_profiles_config_.prog_name = "mafft_profiles_align";
#endif
    align_profiles_config_.in_arg = false; align_profiles_config_.in_arg_str = "--type";
    align_profiles_config_.DNA = "2"; align_profiles_config_.Protein = "1";
    align_profiles_config_.out_to_stdout = false; align_profiles_config_.out_arg_str = "--out";
    align_profiles_config_.single_thread = false; align_profiles_config_.Thread_first = "--thread";
    align_profiles_config_.has_two_profile_prefix = true;
    align_profiles_config_.profile1_prefix = "--list";
    align_profiles_config_.profile2_prefix = "";
    align_profiles_config_.other_args = strdup(("--method " + std::to_string((int)config_json["mafft_cluster_method"])).c_str());
    if(config_json.find("cluster_gapopen") != config_json.end() || config_json.find("cluster_gapext") != config_json.end())
    {
        align_profiles_config_.change_penalty = true;
        if(config_json.find("cluster_gapopen") != config_json.end())
            align_profiles_config_.gapopen = strdup(("-f " + std::to_string((double)config_json["cluster_gapopen"])).c_str());
        else align_profiles_config_.gapopen = nullptr;
        if(config_json.find("cluster_gapext") != config_json.end())
            align_profiles_config_.gapext = strdup(("-g " + std::to_string((double)config_json["cluster_gapext"])).c_str());
        else align_profiles_config_.gapext = nullptr;
    }
    else align_profiles_config_.change_penalty = false;
}

void abPOA_profiles_config(const json &config_json)
{
#ifdef _DEBUG
#if (defined(_WIN32) || defined(_WIN64))
    align_profiles_config_.prog_path = expand_prog_path("\\..\\Release\\abPOA_profile.exe");
#else
    align_profiles_config_.prog_path = expand_prog_path("abpoa_profile");
#endif
#else
#if (defined(_WIN32) || defined(_WIN64))
    align_profiles_config_.prog_path = expand_prog_path("abPOA_profile.exe");
#else
    align_profiles_config_.prog_path = expand_prog_path("/abpoa_profile");
#endif
#endif
#if ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__) ) )
    align_profiles_config_.prog_name = "abPOA_two_profile";
#endif
    align_profiles_config_.in_arg = false; align_profiles_config_.in_arg_str = nullptr;
    align_profiles_config_.DNA = ""; align_profiles_config_.Protein = "-c";
    align_profiles_config_.out_to_stdout = true; align_profiles_config_.out_arg_str = nullptr;
    align_profiles_config_.single_thread = true; align_profiles_config_.Thread_first = nullptr;
    align_profiles_config_.has_two_profile_prefix = true;
    align_profiles_config_.profile1_prefix = "-C";
    align_profiles_config_.profile2_prefix = "-l";
    if(config_json["abpoa_cluster_method"] == 1) align_profiles_config_.other_args = nullptr;
    else align_profiles_config_.other_args = "-p";
    if(config_json.find("merge_gapopen") != config_json.end() || config_json.find("merge_gapext") != config_json.end())
    {
        align_profiles_config_.change_penalty = true;
        if(config_json.find("merge_gapopen") != config_json.end())
            align_profiles_config_.gapopen = strdup(("-O " + std::to_string((double)config_json["merge_gapopen"])).c_str());
        else align_profiles_config_.gapopen = nullptr;
        if(config_json.find("merge_gapext") != config_json.end())
            align_profiles_config_.gapext = strdup(("-E " + std::to_string((double)config_json["merge_gapext"])).c_str());
        else align_profiles_config_.gapext = nullptr;
    }
    else align_profiles_config_.change_penalty = false;
}

#if _DEBUG
void config_show()
{
    /* Show in Windows mode */
    fprintf(stderr, "Windows style arguments:\n");
    /* Part 1: Fragalign */
    std::string arg_string, InFile = "<in>", OutFile = "<out>", thread_str = "<Threads>", block_sequence_dim = "/";
	arg_string += align_seqs_config_.prog_path;
	arg_string += " ";
    if (align_seqs_config_.Protein && align_seqs_config_.DNA)
    {
        arg_string += align_seqs_config_.Protein;
        arg_string += block_sequence_dim;
        arg_string += align_seqs_config_.DNA;
        arg_string += " ";
    }
	if(! align_seqs_config_.single_thread)
	{
		arg_string += align_seqs_config_.Thread_first;
		arg_string += " ";
		arg_string += thread_str;

	}
	arg_string += " ";
    if(align_seqs_config_.other_args)
	{
		arg_string += align_seqs_config_.other_args;
		arg_string += " ";
	}
    if (align_seqs_config_.out_to_stdout)
    {
        if (align_seqs_config_.in_arg)
        {
            arg_string += align_seqs_config_.in_arg_str;
            arg_string += " ";
        }
        arg_string += InFile;
        arg_string += " > ";
        arg_string += OutFile;
    }
    else // fit HAlign-3
    {
        arg_string += align_seqs_config_.out_arg_str;
        arg_string += " ";
        arg_string += OutFile;
        arg_string += " ";
        if (align_seqs_config_.in_arg)
        {
            arg_string += align_seqs_config_.in_arg_str;
            arg_string += " ";
        }
        arg_string += InFile;
        arg_string += " ";
    }
    fprintf(stderr, "Fragalign args = \n%s\n", arg_string.c_str());
    arg_string.clear();
    /* Part 2: two-profile align */
    std::string Profile1 = "<Profile 1 path>", Profile2 = "<Profile 2 path>";
    arg_string += align_two_profile_config_.prog_path;
	arg_string += " ";
	arg_string += align_two_profile_config_.Protein;
    arg_string += block_sequence_dim;
	arg_string += align_two_profile_config_.DNA;
	arg_string += " ";
	if(align_two_profile_config_.has_two_profile_prefix) arg_string += align_two_profile_config_.profile1_prefix;
	arg_string += " ";
	arg_string += Profile1;
	arg_string += " ";
	if(align_two_profile_config_.has_two_profile_prefix) arg_string += align_two_profile_config_.profile2_prefix;
	arg_string += " ";
	arg_string += Profile2;
	arg_string += " ";
	if(align_two_profile_config_.other_args) arg_string += align_two_profile_config_.other_args;
	arg_string += " ";
    if(align_two_profile_config_.out_to_stdout)
    {
        arg_string += "> " + OutFile;
    }
    else arg_string += "\nWill print to " + Profile1;
    fprintf(stderr, "two-profile args = \n%s\n", arg_string.c_str());
    arg_string.clear();
    /* Part 3: profiles align */
    std::string CenterFile = "<Center File>", ListFile = "<List File>";
    arg_string += align_profiles_config_.prog_path;
	arg_string += " ";
	if(align_profiles_config_.has_two_profile_prefix) arg_string += align_profiles_config_.profile1_prefix;
	arg_string += " ";
	arg_string += CenterFile;
	arg_string += " ";
	if(align_profiles_config_.has_two_profile_prefix) arg_string += align_profiles_config_.profile2_prefix;
	arg_string += " ";
	arg_string += ListFile;
	arg_string += " ";
    arg_string += align_profiles_config_.Protein;
    arg_string += block_sequence_dim;
	arg_string += align_profiles_config_.DNA;
	arg_string += " ";
	if(! align_profiles_config_.single_thread)
	{
		arg_string += align_profiles_config_.Thread_first;
		arg_string += " ";
		arg_string += thread_str;
		arg_string += " ";
	}
    fprintf(stderr, "profiles args = \n%s\n", arg_string.c_str());
    /* Show in Linux mode */
    //exit(0);
}
#endif

void aligner_method_init(const json& config_json)
{
    /* seqs method init */
    if(config_json["align_method"] == 0) WMSA_seqs_config(config_json);
    else if(config_json["align_method"] == 1) MAFFT_seqs_config(config_json);
    else if(config_json["align_method"] == 2) HAlign_seqs_config();
    else if(config_json["align_method"] == 3) abPOA_seqs_config(config_json);

    /* two profiles method init */
    if(config_json["merge_method"] == 0) WMSA_two_profile_config(config_json);
    else if(config_json["merge_method"] == 1) MAFFT_two_profile_config(config_json);
    else if(config_json["merge_method"] == 3) abPOA_two_profile_config(config_json);

    /* profiles method init */
    if(config_json["cluster_method"] == 0) WMSA_profiles_config(config_json);
    else if(config_json["cluster_method"] == 1) MAFFT_profiles_config(config_json);
    else if(config_json["cluster_method"] == 3) abPOA_profiles_config(config_json);

#if _DEBUG
    config_show();
#endif
}

void aligner_method_destroy()
{
    if(align_seqs_config_.gapopen) free(align_seqs_config_.gapopen);
    if(align_seqs_config_.gapext)  free(align_seqs_config_.gapext);
    free(align_seqs_config_.prog_path);
    if(align_two_profile_config_.gapopen) free(align_two_profile_config_.gapopen);
    if(align_two_profile_config_.gapext)  free(align_two_profile_config_.gapext);
    free(align_two_profile_config_.prog_path);
    if(align_profiles_config_.gapopen) free(align_profiles_config_.gapopen);
    if(align_profiles_config_.gapext)  free(align_profiles_config_.gapext);
    free(align_profiles_config_.prog_path);
}
