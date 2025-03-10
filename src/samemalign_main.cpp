#include "lib/libfile_v2.hpp"
#include "include/getopt9/include/getopt.h"
#include "lib/ext_prog_config.hpp"
#include "lib/version.h"
#include "include/json.hpp"
#include <string>
#include <cstdlib>
#include <chrono>
#include <random>

#if defined(_WIN32) || defined(_WIN64)
#include <io.h>
#include "include/dirent.h"
#else
#include <unistd.h>
#include <dirent.h>
#endif

#define MAX_ARGS 100

char* c_input = NULL, * c_output = NULL, * tmp_dir = NULL;
char* json_name = NULL;
std::string program_path, json_name_;
// threshold_type: false: ratio, true: threshold number
bool removechr, no_verbose, json_save, mem_regen;
char* no_verbose_s = NULL;
json config_json, tmp_json;
std::mt19937 rand_gen;

void usage(char name[])
{
    fprintf(stderr, "\n\tUsage: %s (%llu-bit) [options]\n\n", name, sizeof(int*) << 3);
    fprintf(stderr, "Use deBruijn graph find MEM, and make alignment by a specified aligner like WMSA/MAFFT/HAlign 3\n");
    fprintf(stderr, "Available options: \n");
    fprintf(stderr, "File process: \n");
    fprintf(stderr, "\t--in         FILE      input sequence file name (Required)\n");
    fprintf(stderr, "\t--out        FILE      output file name (Required)\n");
    fprintf(stderr, "\t--removechr            remove characters in sequence\n");
    fprintf(stderr, "Config JSON file: (Recommended)\n");
    fprintf(stderr, "\t--json       FILE      JSON config, contain \"Conditions\", \"MEM arugments\" and \"Alignment arugments\" configurations\n");
    fprintf(stderr, "\t--jsonsample FILE      output JSON config sample to FILE and exit\n");
    fprintf(stderr, "\t--keepjson             keep JSON file when finish\n");
    fprintf(stderr, "Conditions: \n");
    fprintf(stderr, "\t--thread     N         use N threads to run program (default: 1)\n");
    fprintf(stderr, "\t--type       X         sequence type: X = 1 is Protein, X = 2 is DNA (Required)\n"); // X = 0 is common string
    fprintf(stderr, "MEM arugments: \n");
    fprintf(stderr, "\t--threshold  t         the minimal length of MEM (default = 10)\n");
    fprintf(stderr, "\t--memregen             re-generate MEM (default: false)\n");
    // fprintf(stderr, "\t--ratio      r         the minimal ratio length of MEMs, the ratio must be in (0, 1] (not included 0.0)\n");
    fprintf(stderr, "Alignment arugments: \n");
    fprintf(stderr, "\t--highsim              this file is high similarity, will treat all sequences into one cluster; will ignore --cluster,\n");
    fprintf(stderr, "\t                       --wmsac and --mafftc\n");
    fprintf(stderr, "\t--path       PATH      temp file will be in path, default = \"swap/\" (only support first-level directory)\n");
    fprintf(stderr, "\t--threshold2 t         the minimal length of blockaligner (default = 10)\n");
    fprintf(stderr, "\t                       not recommended for t < 10\n");
    fprintf(stderr, "\t--align      A         alignment method: A = 0 is WMSA, A = 1 is MAFFT, A = 2 is HAlign, A = 3 is abPOA (default: 0)\n");
    fprintf(stderr, "\t--wmsaa      w         WMSA align method: w = 0 use FFT+K-band align, w = 1 only use K-band align (default: 0)\n");
    fprintf(stderr, "\t--maffta     m         MAFFT align method: m = 0 is NW-NS-1, m = 1 is NW-NS-2, m = 2 is NW-NS-i, m = 3 is FFT-NS-1,\n");
    fprintf(stderr, "\t                                           m = 4 is FFT-NS-2, m = 5 is FFT-NS-i, m = 6 is G-INS-i, m = 7 is L-INS-i\n");
    fprintf(stderr, "\t                                           (default: 3)\n");
    fprintf(stderr, "\t--gapopena   a         use user-defined gap open/gap penalty in align for WMSA, MAFFT and abPOA\n");
    fprintf(stderr, "\t--gapexta    a         use user-defined gap extension in align for MAFFT and abPOA\n");
    fprintf(stderr, "\t--merge      M         merge two profiles method: M = 0 is WMSA, M = 1 is MAFFT, M = 2 is abPOA (default: 0)\n");
    fprintf(stderr, "\t--wmsam      w         WMSA merge method: w = 0 use FFT+K-band with no gap open calculation, \n");
    fprintf(stderr, "\t                                          w = 1 use K-band with no gap open calculation,\n");
    fprintf(stderr, "\t                                          w = 2 use FFT+K-band with gap open calculation,\n");
    fprintf(stderr, "\t                                          w = 3 use K-band with gap open calculation (default: 0)\n");
    fprintf(stderr, "\t--mafftm     m         MAFFT two profiles align method: m = 0 is NW-NS, m = 1 is FFT-NS (default: 1)\n");
    fprintf(stderr, "\t--gapopenm   m         use user-defined gap open/gap penalty in merge two profiles for WMSA, MAFFT and abPOA\n");
    fprintf(stderr, "\t--gapextm    m         use user-defined gap extension in merge two profiles in WMSA, MAFFT and abPOA\n");
    fprintf(stderr, "\t--cluster    C         merge clusters method: C = 0 is WMSA, C = 1 is MAFFT, C = 2 is abPOA (default: 0)\n");
    fprintf(stderr, "\t--wmsac      w         WMSA cluster merge method: w = 0 use FFT+K-band with no gap open calculation, \n");
    fprintf(stderr, "\t                                                  w = 1 use K-band with no gap open calculation,\n");
    fprintf(stderr, "\t                                                  w = 2 use FFT+K-band with gap open calculation,\n");
    fprintf(stderr, "\t                                                  w = 3 use K-band with gap open calculation (default: 1)\n");
    fprintf(stderr, "\t--mafftc     m         MAFFT cluster merge method: m = 0 is NW-NS-1, m = 1 is NW-NS-2, m = 2 is NW-NS-i, m = 3 is FFT-NS-1,\n");
    fprintf(stderr, "\t                                                   m = 4 is FFT-NS-2, m = 5 is FFT-NS-i, m = 6 is G-INS-i, m = 7 is L-INS-i\n");
    fprintf(stderr, "\t                                                   (default: 3)\n");
    fprintf(stderr, "\t--abpoac     p         abPOA cluster merge method: p = 0 use center sequences make pre-align and insert clusters,\n");
    fprintf(stderr, "\t                                                   p = 1 directly align clusters (default: 0)\n");
    fprintf(stderr, "\t--gapopenc   c         use user-defined gap open/gap penalty in merge clusters for WMSA, MAFFT and abPOA\n");
    fprintf(stderr, "\t--gapextc    c         use user-defined gap extension in merge clusters for WMSA, MAFFT and abPOA\n");
    fprintf(stderr, "\t--allmafft             all method use MAFFT\n");
    fprintf(stderr, "\t--allwmsa              all method use WMSA\n");
    fprintf(stderr, "\t--allabpoa             all method use abPOA\n");
    fprintf(stderr, "Others: \n");
    fprintf(stderr, "\t--noverbose            not print the info for subprograms\n");
    fprintf(stderr, "\t--help                 print help message\n");
    fprintf(stderr, "\t--version              show program version\n");
    fprintf(stderr, "Hint: If read json file and use same arugments, program will omit the arugment settings.\n");
    fprintf(stderr, "Example:\n\t%s --in seq.fasta --out seq_out.fasta\n", name);
}

void version()
{
#if (defined(_WIN32) || defined(_WIN64))
    program_path = _pgmptr;
    program_path.erase(program_path.find_last_of('\\') + 1, program_path.size());
#else
    char progpath[PATH_MAX];
    readlink("/proc/self/exe", progpath, PATH_MAX);
    program_path = progpath;
    program_path.erase(program_path.find_last_of('/') + 1, program_path.size());
#endif
    fprintf(stderr, "SAMEMAlign version %d.%d.%d.%d\n", VER_MAJOR, VER_MINOR, VER_RELEASE_MAIN, VER_BUILD);
    fprintf(stderr, "Submodules version:\n");
#if defined(_WIN32) || defined(_WIN64)
    std::string prog_args = program_path + ".\\splitmem.exe --version";
    // splitMEM
    system_spawn_Win((char*)prog_args.c_str(), NULL, NULL, NULL);
    // block dispatch
    prog_args = program_path + ".\\block_dispatch.exe --version";
    system_spawn_Win((char*)prog_args.c_str(), NULL, NULL, NULL);
    // mafft
    prog_args = program_path + "..\\usr\\lib\\mafft\\version.exe";
    fprintf(stderr, "MAFFT version: ");
    system_spawn_Win((char*)prog_args.c_str(), NULL, NULL, NULL);
    // abPOA
    prog_args = program_path + ".\\abpoa.exe --version";
    fprintf(stderr, "\nabPOA version: ");
    system_spawn_Win((char*)prog_args.c_str(), NULL, NULL, NULL);
    // subWMSA
    fprintf(stderr, "WMSA components version:\n");
    prog_args = program_path + ".\\fragalign.exe -v";
    system_spawn_Win((char*)prog_args.c_str(), NULL, NULL, NULL);
    prog_args = program_path + ".\\profile_two_align.exe -v";
    system_spawn_Win((char*)prog_args.c_str(), NULL, NULL, NULL);
    prog_args = program_path + ".\\profilealign.exe -v";
    system_spawn_Win((char*)prog_args.c_str(), NULL, NULL, NULL);

#else
    // splitMEM
    char prog_name[program_path.length() + 50];
    const char** argv = (const char**)malloc(sizeof(char*) * 3);
    if (argv == NULL) { fprintf(stderr, "ERROR: can not allocate enough memory for arugment. Program will exit.\n"); exit(1); }
    sprintf(prog_name, "%ssplitmem", program_path.c_str());
    argv[0] = "splitmem";
    argv[1] = "--version";
    argv[2] = NULL;
    system_spawn((const char*)prog_name, argv, NULL, NULL, NULL);
    // block dispatch
    sprintf(prog_name, "%sblock_dispatch", program_path.c_str());
    argv[0] = "block_dispatch";
    argv[1] = "--version";
    system_spawn((const char*)prog_name, argv, NULL, NULL, NULL);
    // mafft
    sprintf(prog_name, "%s../libexec/mafft/version", program_path.c_str());
    argv[0] = "mafft";
    argv[1] = NULL;
    fprintf(stderr, "MAFFT version: ");
    system_spawn((const char*)prog_name, argv, NULL, NULL, NULL);
    // abPOA
    sprintf(prog_name, "%sabpoa", program_path.c_str());
    argv[0] = "abpoa";
    argv[1] = "--version";
    fprintf(stderr, "\nabPOA version: ");
    system_spawn((const char*)prog_name, argv, NULL, NULL, NULL);
    // subWMSA
    fprintf(stderr, "WMSA components version:\n");
    sprintf(prog_name, "%sfragalign", program_path.c_str());
    argv[0] = "fragalign";
    argv[1] = "-v";
    system_spawn((const char*)prog_name, argv, NULL, NULL, NULL);
    sprintf(prog_name, "%sprofile_two_align", program_path.c_str());
    argv[0] = "profile_two_align";
    system_spawn((const char*)prog_name, argv, NULL, NULL, NULL);
    sprintf(prog_name, "%sprofilealign", program_path.c_str());
    argv[0] = "profilealign";
    system_spawn((const char*)prog_name, argv, NULL, NULL, NULL);

#endif
    exit(0);
}


bool checkallprogram()
{
#if (defined(_WIN32) || defined(_WIN64))
    program_path = _pgmptr;
    program_path.erase(program_path.find_last_of('\\') + 1, program_path.size());
#else
    char progpath[PATH_MAX];
    readlink("/proc/self/exe", progpath, PATH_MAX);
    program_path = progpath;
    program_path.erase(program_path.find_last_of('/') + 1, program_path.size());
#endif
    //fprintf(stderr, "Program Directory = %s\n", program_path.c_str());
#if (defined(_WIN32) || defined(_WIN64))
#ifdef _DEBUG
    if(access((program_path + "..\\Release\\profilealign.exe").c_str(), 0)) { fprintf(stderr, "Error: cannot found profilealign.exe."); return false; }
    if(access((program_path + "..\\Release\\profile_two_align.exe").c_str(), 0)) { fprintf(stderr, "Error: cannot found profile_two_align.exe."); return false; }
    if(access((program_path + "..\\Release\\fragalign.exe").c_str(), 0)) { fprintf(stderr, "Error: cannot found fragalign.exe."); return false; }
    if(access((program_path + "..\\Release\\abPOA.exe").c_str(), 0)) { fprintf(stderr, "Error: cannot found abPOA.exe."); return false; }
    if(access((program_path + "..\\Release\\abPOA_profile.exe").c_str(), 0)) { fprintf(stderr, "Error: cannot found abPOA_profile.exe."); return false; }
#else
    if(access((program_path + "profilealign.exe").c_str(), 0)) { fprintf(stderr, "Error: cannot found profilealign.exe."); return false; }
    if(access((program_path + "profile_two_align.exe").c_str(), 0)) { fprintf(stderr, "Error: cannot found profile_two_align.exe."); return false; }
    if(access((program_path + "fragalign.exe").c_str(), 0)) {{ fprintf(stderr, "Error: cannot found fragalign.exe."); return false; }}
    if(access((program_path + "abPOA.exe").c_str(), 0)) {{ fprintf(stderr, "Error: cannot found abPOA.exe."); return false; }}
    if(access((program_path + "abPOA_profile.exe").c_str(), 0)) {{ fprintf(stderr, "Error: cannot found abPOA_profile.exe."); return false; }}
#endif // defined(_DEBUG)
    if(access((program_path + "splitmem.exe").c_str(), 0)) { fprintf(stderr, "Error: cannot found splitmem.exe."); return false; }
    if(access((program_path + "block_dispatch.exe").c_str(), 0)) { fprintf(stderr, "Error: cannot found blockalign.exe."); return false; }
    if(access((program_path + "seq_remove_char.exe").c_str(), 0)) { fprintf(stderr, "Error: cannot found seq_remove_char.exe."); return false; }
#else
    if(access((program_path + "profilealign").c_str(), 0)) { fprintf(stderr, "Error: cannot found profilealign."); return false; }
    if(access((program_path + "profile_two_align").c_str(), 0)) { fprintf(stderr, "Error: cannot found profile_two_align."); return false; }
    if(access((program_path + "fragalign").c_str(), 0)) { fprintf(stderr, "Error: cannot found fragalign."); return false; }
    if(access((program_path + "abpoa").c_str(), 0)) { fprintf(stderr, "Error: cannot found abpoa."); return false; }
    if(access((program_path + "abpoa_profile").c_str(), 0)) { fprintf(stderr, "Error: cannot found abpoa_profile."); return false; }
    if(access((program_path + "splitmem").c_str(), 0)) { fprintf(stderr, "Error: cannot found splitmem."); return false; }
    if(access((program_path + "block_dispatch").c_str(), 0)) { fprintf(stderr, "Error: cannot found blockalign."); return false; }
    if(access((program_path + "seq_remove_char").c_str(), 0)) { fprintf(stderr, "Error: cannot found seq_remove_char."); return false; }
#endif
    return true;
}

void get_args(int argc, char* argv[])
{
    extern char* optarg;
    extern int optind;
    int c, threads_num;
    no_verbose = false;
    removechr = false;
    json_save = false;
    mem_regen = false;

    config_json["seq_type"] = 0;
    config_json["highsim"] = false;
    config_json["all_mafft"] = false;
    config_json["all_wmsa"] = false;
    config_json["all_abpoa"] = false;
    config_json["align_method"] = 0;
    config_json["wmsa_align_method"] = 0;
    config_json["mafft_align_method"] = 3;
    config_json["merge_method"] = 0;
    config_json["wmsa_merge_method"] = 3;
    config_json["mafft_merge_method"] = 1;
    config_json["cluster_method"] = 0;
    config_json["wmsa_cluster_method"] = 0;
    config_json["mafft_cluster_method"] = 3;
    config_json["abpoa_cluster_method"] = 0;
    config_json["threshold"] = 10;
    config_json["threshold_block"] = 10;
    config_json["path"] = "swap/"; tmp_dir = "swap/";

#if ( defined(_WIN32) || defined(_WIN64) )
    std::string consts[] = { "in", "out", "help", "version", "thread", "type", "threshold", "ratio", "noverbose", "threshold2", "highsim", "align",
                             "maffta", "merge", "mafftm", "cluster", "mafftc", "allmafft", "allwmsa", "removechr", "wmsaa", "wmsam", "wmsac", "json",
                             "jsonsample", "gapopena", "gapexta", "gapopenm", "gapextm", "gapopenc", "gapextc", "path", "keepjson", "allabpoa",
                             "abpoac", "memregen" };
#endif

    while (1)
    {
        //int this_option_optind = optind ? optind : 1;
        int option_index = 0;
#if ( defined(_WIN32) || defined(_WIN64) )
        static struct option long_options[] =
        {
            {(char*)consts[0].c_str(),  required_argument, 0, 'f'},
            {(char*)consts[1].c_str(),  required_argument, 0, 'o'},
            {(char*)consts[2].c_str(),  no_argument,       0, 'h'},
            {(char*)consts[3].c_str(),  no_argument,       0, 'v'},
            {(char*)consts[4].c_str(),  required_argument, 0, 't'},
            {(char*)consts[5].c_str(),  required_argument, 0, 'T'},
            {(char*)consts[6].c_str(),  required_argument, 0, 'x'},
            {(char*)consts[7].c_str(),  required_argument, 0, 'r'},
            {(char*)consts[8].c_str(),  no_argument,       0, 'n'},
            {(char*)consts[9].c_str(),  required_argument, 0, 's'},
            {(char*)consts[10].c_str(), no_argument,       0, 'H'},
            {(char*)consts[11].c_str(), required_argument, 0, 'A'},
            {(char*)consts[12].c_str(), required_argument, 0, 'a'},
            {(char*)consts[13].c_str(), required_argument, 0, 'M'},
            {(char*)consts[14].c_str(), required_argument, 0, 'm'},
            {(char*)consts[15].c_str(), required_argument, 0, 'C'},
            {(char*)consts[16].c_str(), required_argument, 0, 'c'},
            {(char*)consts[17].c_str(), no_argument,       0, '['},
            {(char*)consts[18].c_str(), no_argument,       0, ']'},
            {(char*)consts[19].c_str(), no_argument,       0, ':'},
            {(char*)consts[20].c_str(), required_argument, 0, 'l'},
            {(char*)consts[21].c_str(), required_argument, 0, 'e'},
            {(char*)consts[22].c_str(), required_argument, 0, 'u'},
            {(char*)consts[23].c_str(), required_argument, 0, 'J'},
            {(char*)consts[24].c_str(), required_argument, 0, 'j'},
            {(char*)consts[25].c_str(), required_argument, 0, '1'},
            {(char*)consts[26].c_str(), required_argument, 0, '!'},
            {(char*)consts[27].c_str(), required_argument, 0, '2'},
            {(char*)consts[28].c_str(), required_argument, 0, '@'},
            {(char*)consts[29].c_str(), required_argument, 0, '3'},
            {(char*)consts[30].c_str(), required_argument, 0, '#'},
            {(char*)consts[31].c_str(), required_argument, 0, 'p'},
            {(char*)consts[32].c_str(), no_argument,       0, 'K'},
            {(char*)consts[33].c_str(), no_argument,       0, ';'},
            {(char*)consts[34].c_str(), required_argument, 0, 'P'},
            {(char*)consts[35].c_str(), no_argument,       0, '0'},
            {0,                         0,                 0,  0 }
        };
#else
        static struct option long_options[] =
        {
            {"in",         required_argument, 0, 'f'},
            {"out",        required_argument, 0, 'o'},
            {"help",       no_argument,       0, 'h'},
            {"version",    no_argument,       0, 'v'},
            {"thread",     required_argument, 0, 't'},
            {"type",       required_argument, 0, 'T'},
            {"threshold",  required_argument, 0, 'x'},
            {"ratio",      required_argument, 0, 'r'},
            {"noverbose",  no_argument,       0, 'n'},
            {"threshold2", required_argument, 0, 's'},
            {"highsim",    no_argument,       0, 'H'},
            {"align",      required_argument, 0, 'A'},
            {"maffta",     required_argument, 0, 'a'},
            {"merge",      required_argument, 0, 'M'},
            {"mafftm",     required_argument, 0, 'm'},
            {"cluster",    required_argument, 0, 'C'},
            {"mafftc",     required_argument, 0, 'c'},
            {"allmafft",   no_argument,       0, '['},
            {"allwmsa",    no_argument,       0, ']'},
            {"removechr",  no_argument,       0, ':'},
            {"wmsaa",      required_argument, 0, 'l'},
            {"wmsam",      required_argument, 0, 'e'},
            {"wmsac",      required_argument, 0, 'u'},
            {"json",       required_argument, 0, 'J'},
            {"jsonsample", required_argument, 0, 'j'},
            {"gapopena",   required_argument, 0, '1'},
            {"gapexta",    required_argument, 0, '!'},
            {"gapopenm",   required_argument, 0, '2'},
            {"gapextm",    required_argument, 0, '@'},
            {"gapopenc",   required_argument, 0, '3'},
            {"gapextc",    required_argument, 0, '#'},
            {"path",       required_argument, 0, 'p'},
            {"keepjson",   no_argument,       0, 'K'},
            {"allabpoa",   no_argument,       0, ';'},
            {"abpoac",     required_argument, 0, 'P'},
            {"memregen",   no_argument,       0, '0'},
            {0,            0,                 0,  0 }
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
        case 't':
            threads_num = atoi(argv[optind - 1]);
            if(! threads_num || threads_num < 0) { fprintf(stderr, "Invaild threads: %s\n", argv[optind - 1]); break; }
            fprintf(stderr, "Use %d thread(s) for all programs.\n", threads_num);
            config_json["thread"] = threads_num;
            break;
        case 'T':
            config_json["seq_type"] = atoi(argv[optind - 1]);
            break;
        case 'x':
            if(config_json.find("threhsold") != config_json.end())
                fprintf(stderr, "Warning: will omit the setting: %s %s\n", "--threshold", config_json["threshold"]);
            if(config_json.find("ratio") != config_json.end())
                fprintf(stderr, "Warning: will omit the setting: %s %s\n", "--ratio", config_json["ratio"]);
            config_json.erase("threshold"); config_json.erase("ratio");
            config_json["threshold"] = atoi(argv[optind - 1]);
            break;
        case 'r':
#if 1
            fprintf(stderr, "Warning: currently not supported this arugment. Please use --threshold instead.\n");
            break;
#else
            if(config_json.find("threhsold") != config_json.end())
                fprintf(stderr, "Warning: will omit the setting: %s %s\n", "--threshold", config_json["threshold"]);
            if(config_json.find("ratio") != config_json.end())
                fprintf(stderr, "Warning: will omit the setting: %s %s\n", "--ratio", config_json["ratio"]);
            config_json.erase("threshold"); config_json.erase("ratio");
            config_json["ratio"] = atoi(argv[optind - 1]);
            break;
#endif
        case 'n':
            no_verbose = true;
            break;
        case 'v':
            version();
        case 's':
            config_json["threshold_block"] = atoi(argv[optind - 1]);
            break;
        case 'H':
            config_json["highsim"] = true;
            break;
        case 'A':
            config_json["align_method"] = atoi(argv[optind - 1]);
            fprintf(stderr, "Alignment method:"); print_method(config_json["align_method"]);
            break;
        case 'a':
            config_json["mafft_align_method"] = atoi(argv[optind - 1]);
            fprintf(stderr, "Alignment MAFFT method:"); print_mafft_method(config_json["mafft_align_method"]);
            break;
        case 'l':
            config_json["wmsa_align_method"] = atoi(argv[optind - 1]);
            fprintf(stderr, "Alignment WMSA method:"); print_wmsa_align_method(config_json["wmsa_align_method"]);
            break;
        case '1':
            config_json["align_gapopen"] = atof(argv[optind - 1]);
            fprintf(stderr, "Gap open in align: %f\n", config_json["align_gapopen"]);
            break;
        case '!':
            config_json["align_gapext"] = atof(argv[optind - 1]);
            fprintf(stderr, "Gap extension in align: %f\n", config_json["align_gapext"]);
            break;
        case 'M':
            config_json["merge_method"] = atoi(argv[optind - 1]);
            if(config_json["merge_method"] == 2) // it means method = abPOA
            {
                config_json["merge_method"] = 3;
                // abPOA
            }
            fprintf(stderr, "Merge method:"); print_method(config_json["merge_method"]);
            break;
        case 'm':
            config_json["mafft_merge_method"] = atoi(argv[optind - 1]);
            fprintf(stderr, "Merge MAFFT method:"); print_mafft_merge_method(config_json["mafft_merge_method"]);
            break;
        case 'e':
            config_json["wmsa_merge_method"] = atoi(argv[optind - 1]);
            fprintf(stderr, "Merge WMSA method:"); print_wmsa_merge_method(config_json["wmsa_merge_method"]);
            break;
        case '2':
            config_json["merge_gapopen"] = atof(argv[optind - 1]);
            fprintf(stderr, "Gap open in merge: %f\n", config_json["merge_gapopen"]);
            break;
        case '@':
            config_json["merge_gapext"] = atof(argv[optind - 1]);
            fprintf(stderr, "Gap extension in merge: %f\n", config_json["merge_gapext"]);
            break;
        case 'C':
            config_json["cluster_method"] = atoi(argv[optind - 1]);
            if(atoi(argv[optind - 1]) == 2)  // it means method = abPOA
            {
                config_json["merge_method"] = 3;
                // abPOA
            }
            fprintf(stderr, "Cluster merge method:"); print_method(config_json["cluster_method"]);
            break;
        case 'c':
            config_json["mafft_cluster_method"] = atoi(argv[optind - 1]);
            fprintf(stderr, "Cluster merge MAFFT method:"); print_mafft_method(config_json["mafft_cluster_method"]);
            break;
        case 'u':
            config_json["wmsa_cluster_method"] = atoi(argv[optind - 1]);
            fprintf(stderr, "Cluster merge WMSA method:"); print_wmsa_merge_method(config_json["wmsa_cluster_method"]);
            break;
        case 'P':
            config_json["abpoa_cluster_method"] = atoi(argv[optind - 1]);
            fprintf(stderr, "Cluster merge abPOA method:"); print_abpoa_cluster_method(config_json["abpoa_cluster_method"]);
            break;
        case '3':
            config_json["cluster_gapopen"] = atof(argv[optind - 1]);
            fprintf(stderr, "Gap cluster in merge: %f\n", config_json["cluster_gapopen"]);
            break;
        case '#':
            config_json["cluster_gapext"] = atof(argv[optind - 1]);
            fprintf(stderr, "Gap extension in cluster: %f\n", config_json["cluster_gapext"]);
            break;
        case '[':
            config_json["all_mafft"] = true;
            if (config_json.find("all_wmsa") != config_json.end() && config_json["all_wmsa"])
            {
                fprintf(stderr, "Warning: conflict arguments: --allmafft vs --allwmsa. Program will use --allmafft.\n");
                config_json["all_wmsa"] = false;
            }
            if (config_json.find("all_abpoa") != config_json.end() && config_json["all_abpoa"])
            {
                fprintf(stderr, "Warning: conflict arguments: --allmafft vs --allabpoa. Program will use --allmafft.\n");
                config_json["all_abpoa"] = false;
            }
            break;
        case ']':
            config_json["all_wmsa"] = true;
            if (config_json.find("all_mafft") != config_json.end() && config_json["all_mafft"])
            {
                fprintf(stderr, "Warning: conflict arguments: --allwmsa vs --allmafft. Program will use --allwmsa.\n");
                config_json["all_mafft"] = false;
            }
            if (config_json.find("all_abpoa") != config_json.end() && config_json["all_abpoa"])
            {
                fprintf(stderr, "Warning: conflict arguments: --allwmsa vs --allabpoa. Program will use --allwmsa.\n");
                config_json["all_abpoa"] = false;
            }
            break;
        case ';':
            config_json["all_abpoa"] = true;
            if (config_json.find("all_wmsa") != config_json.end() && config_json["all_wmsa"])
            {
                fprintf(stderr, "Warning: conflict arguments: --allabpoa vs --allwmsa. Program will use --allabpoa.\n");
                config_json["all_wmsa"] = false;
            }
            if (config_json.find("all_mafft") != config_json.end() && config_json["all_mafft"])
            {
                fprintf(stderr, "Warning: conflict arguments: --allabpoa vs --allmafft. Program will use --allabpoa.\n");
                config_json["all_mafft"] = false;
            }
            break;
        case ':':
            removechr = true;
            break;
        case 'J':
            json_name = argv[optind - 1];
            fprintf(stderr, "Info: Program will read config from file %s, will omit settings.\n", json_name);
            break;
        case 'j':
            config_json["seq_type"] = 1; // sample
            config_json["align_gapopen"] = 100.0;
            config_json["align_gapext"] = 0.0;
            config_json["merge_gapopen"] = 100.0;
            config_json["merge_gapext"] = 0.0;
            config_json["cluster_gapopen"] = 100.0;
            config_json["cluster_gapext"] = 0.0;
            print_json_sample(argv[optind - 1], config_json);
            exit(0);
            break;
        case 'K':
            json_save = true;
            break;
        case 'p':
            config_json["path"] = argv[optind - 1];
            tmp_dir = argv[optind - 1];
            break;
        case '0':
            mem_regen = true;
            fprintf(stderr, "Info: regen MEM.\n");
            break;
        }
    }

    if (json_name)
    {
        tmp_json = json();
        json_read(json_name, tmp_json);
        fprintf(stderr, "JSON = %s\n", config_json.dump().c_str());
        config_json.merge_patch(tmp_json);
    }

    if(config_json["all_wmsa"] == true)
    {
        config_json["align_method"] = 0;   if(config_json.find("wmsa_align_method") == config_json.end())   config_json["wmsa_align_method"] = 0;
        config_json["merge_method"] = 0;   if(config_json.find("wmsa_merge_method") == config_json.end())   config_json["wmsa_merge_method"] = 2;
        config_json["cluster_method"] = 0; if(config_json.find("wmsa_cluster_method") == config_json.end()) config_json["wmsa_cluster_method"] = 2;
    }

    if(config_json["all_mafft"] == true)
    {
        config_json["align_method"] = 1;   if(config_json.find("mafft_align_method") == config_json.end())   config_json["mafft_align_method"] = 3;
        config_json["merge_method"] = 1;   if(config_json.find("mafft_merge_method") == config_json.end())   config_json["mafft_merge_method"] = 0;
        config_json["cluster_method"] = 1; if(config_json.find("mafft_cluster_method") == config_json.end()) config_json["mafft_cluster_method"] = 3;
    }

    if (config_json.find("thread") == config_json.end())
    {
        config_json["thread"] = 1;
        fprintf(stderr, "Info: Use 1 thread for all programs.\n");
    }
    if (c_input == NULL || c_output == NULL)
    {
        fprintf(stderr, "ERROR: Please determine input/output file. Use %s --help for more information.\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if (!(config_json["seq_type"] == 1 || config_json["seq_type"] == 2))
    {
        fprintf(stderr, "ERROR: Please determine sequence type by adding arugment: Protein: --type 1, DNA: --type 2. Use %s --help for more information.\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if(no_verbose)
    {
#if defined(_WIN32) || defined(_WIN64)
        no_verbose_s = "nul\0";
#else
        no_verbose_s = "/dev/null\0";
#endif
    }
    // fprintf(stderr, "Blockaligner Alignment arugments:\n");
    // print_program_method(config_json);
    // fprintf(stderr, "JSON = %s", config_json.dump().c_str());
}

int removecharcall()
{
    std::string file_name = c_input;
    file_name += ".tmp";
    Filemove(c_input, file_name.c_str());
    int exitval;
#if defined(_WIN32) || defined(_WIN64)
    std::string prog = program_path + ".\\seq_remove_char.exe ";
    prog += file_name + " ";
    prog += c_input;
    exitval = system_spawn_Win((char*)prog.c_str(), NULL, NULL, NULL);
    exitval |= Filedelete(file_name.c_str());
#else
    char prog_name[program_path.length() + 30];
    const char** argv = (const char**)malloc(4 * sizeof(char *));
    if(argv == NULL) { fprintf(stderr, "ERROR: can not allocate enough memory for removing char. Program will exit.\n"); exit(1); }
    argv[0] = "seq_remove_char"; argv[1] = file_name.c_str(); argv[2] = (const char*)c_input; argv[3] = NULL;
    sprintf(prog_name, "%sseq_remove_char\0", program_path.c_str());
    exitval = system_spawn((const char*)prog_name, argv, NULL, NULL, NULL);
    exitval |= Filedelete(file_name.c_str());
#endif
    return exitval;
}

int splitmem_call()
{
    makedir(tmp_dir);
#if defined(_WIN32) || defined(_WIN64)
    std::string prog_args = program_path + ".\\splitmem.exe -i ";
    prog_args += c_input;
    prog_args += " -o ";
    prog_args += c_input; prog_args += ".mem";
    prog_args += " -d "; prog_args += tmp_dir;
    if(config_json.find("threshold") != config_json.end())
        prog_args += " -k " + std::to_string((int)config_json["threshold"]);
    if(mem_regen) prog_args += " -n ";
    if(! no_verbose) return system_spawn_Win((char *)prog_args.c_str(), NULL, NULL, NULL);
    else return system_spawn_Win((char *)prog_args.c_str(), NULL, NULL, no_verbose_s);
#else
    char prog_name[program_path.length() + 15];
    int exitval;
    sprintf(prog_name, "%ssplitmem\0", program_path.c_str());
    const char **argv = (const char**)malloc(sizeof(char *) * MAX_ARGS);
    if(argv == NULL) { fprintf(stderr, "ERROR: can not allocate enough memory for arugment. Program will exit.\n"); exit(1); }
    std::string prog_out = c_input; prog_out += ".mem";
    argv[0] = "splitmem";
    argv[1] = "-i"; argv[2] = c_input;
    argv[3] = "-o"; argv[4] = prog_out.c_str();
    argv[5] = "-d"; argv[6] = tmp_dir;
    if(config_json.find("threshold") != config_json.end())
    {
        argv[7] = "-k";
        argv[8] = strdup(std::to_string((int)config_json["threshold"]).c_str());
        if(mem_regen) argv[9] = "-n", argv[10] = NULL;
        else argv[9] = NULL;
    }
    else
    {
        if(mem_regen) argv[7] = "-n";
        else argv[7] = NULL;
        argv[8] = NULL;
    }
    if(! no_verbose) exitval = system_spawn((const char*)prog_name, argv, NULL, NULL, NULL);
    else exitval = system_spawn((const char*)prog_name, argv, NULL, NULL, no_verbose_s);
    if(argv[8]) free((void*)argv[8]);
    return exitval;
#endif
}

int blockalign_call()
{
    rand_gen.seed(std::random_device{}());
    json_name_ = std::to_string(rand_gen()) + ".json";
    fprintf(stderr, "Info: config will print to json file %s\n", json_name_.c_str());
    print_json(json_name_, config_json);
#if defined(_WIN32) || defined(_WIN64)
    std::string prog_args = program_path + ".\\block_dispatch.exe --in ";
    prog_args += c_input; prog_args += ".mem";
    prog_args += " --seq "; prog_args += c_input;
    prog_args += " --out "; prog_args += c_output;
    prog_args += " --json "; prog_args += json_name_;
    if (! no_verbose) return system_spawn_Win((char *)prog_args.c_str(), NULL, NULL, NULL);
    else return system_spawn_Win((char *)prog_args.c_str(), NULL, NULL, no_verbose_s);
#else
    char prog_name[program_path.length() + 20];
    sprintf(prog_name, "%sblock_dispatch", program_path.c_str());
    const char **argv = (const char**)malloc(sizeof(char *) * MAX_ARGS);
    int now_argv;
    if(argv == NULL) { fprintf(stderr, "ERROR: can not allocate enough memory for arugment. Program will exit.\n"); exit(1); }
    std::string file_mem = c_input;
    file_mem += ".mem";
    argv[0] = "block_dispatch";
    argv[1] = "--in";   argv[2] = (char*)file_mem.c_str();
    argv[3] = "--out";  argv[4] = c_output;
    argv[5] = "--seq";  argv[6] = c_input;
    argv[7] = "--json"; argv[8] = (char*)json_name_.c_str();
    argv[9] = NULL;
    if(! no_verbose) return system_spawn((const char*)prog_name, argv, NULL, NULL, NULL);
    else return system_spawn((const char*)prog_name, argv, NULL, NULL, no_verbose_s);
#endif
}

int main(int argc, char* argv[])
{
    if(! checkallprogram()) { fprintf(stderr, "\nFATAL ERROR: cannot found compoment. Program will exit.\n"); exit(-1); }
    get_args(argc, argv);
    if (removechr)
    {
        fprintf(stderr, "Info: detecting and removing characters...\n");
        if (removecharcall()) { fprintf(stderr, "Error: can not execute seq_remove_char. Program will exit.\n"); exit(1); }
        fprintf(stderr, "Done. ");
    }
    std::chrono::steady_clock::time_point t_start = std::chrono::steady_clock::now(), t_MEM, t_align;
    fprintf(stderr, "Finding MEMs...\n");
    // Step 1: call findmem
    if(splitmem_call())    { fprintf(stderr, "Error: can not execute splitmem. Program will exit.\n"); exit(1); }
    else fprintf(stderr, "Done. Now making alignment by clusters...\n");
    t_MEM = std::chrono::steady_clock::now();
    // Step 2: call blockaligner
    if(blockalign_call()) { fprintf(stderr, "Error: can not execute blockaligner. Program will exit.\n"); exit(1); }
    else fprintf(stderr, "Done. \n");
    // json is generated by program and clean it
    fprintf(stderr, "Cleaning json temp file %s... ", json_name_.c_str());
    if(Filedelete(json_name_.c_str())) fprintf(stderr, "failed. Please remove it manually.\n");
    else fprintf(stderr, "done.\n");
    t_align = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_all, time_MEM, time_align;
    time_all   = std::chrono::duration_cast<std::chrono::duration<double>>(t_align - t_start);
    time_MEM   = std::chrono::duration_cast<std::chrono::duration<double>>(t_MEM - t_start);
    time_align = std::chrono::duration_cast<std::chrono::duration<double>>(t_align - t_MEM);
    fprintf(stderr,
            "\n\nTime usage: \n- All: %lf s\n- MEM: %lf s\n- alignment: %lf s\n\n\n",
            time_all.count(), time_MEM.count(), time_align.count());
    version();
}