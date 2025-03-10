#include "ext_prog_call.hpp"

extern aligner_config align_seqs_config_, align_two_profile_config_, align_profiles_config_;
extern char* thread_str;
extern int block_sequence_type;

#if ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__)) )
const char delim[2] = " ";
inline void string_to_argv(std::string &str, wordexp_t &p, const char *_func_)
{
#if 0
	fprintf(stderr, "Arg string = %s, called by %s\n", str.c_str(), _func_);
#endif
    if(wordexp(str.c_str(), &p, 0)) { fprintf(stderr, "Error: cannot express string %s. Program will exit\n", str.c_str()); exit(1); }
}
#endif

// Profile-Profile Alignment. 0 is ok, other is bad.
// It will call wmsa profile-profile alignment program.
int WMSAProfileProfilealignment(const char *CenterFile, const char* ListFile, const char *FinalFile)
{
	std::string arg_string;
	arg_string += align_profiles_config_.prog_path;
	arg_string += " ";
	if(align_profiles_config_.has_two_profile_prefix) arg_string += align_profiles_config_.profile1_prefix;
	arg_string += " \"";
	arg_string += CenterFile;
	arg_string += "\" ";
	if(align_profiles_config_.has_two_profile_prefix) arg_string += align_profiles_config_.profile2_prefix;
	arg_string += " \"";
	arg_string += ListFile;
	arg_string += "\" ";
	if(block_sequence_type == 1) { if(align_profiles_config_.Protein) arg_string += align_profiles_config_.Protein; }
	else if(block_sequence_type == 2) { if(align_profiles_config_.DNA) arg_string += align_profiles_config_.DNA; }
	else { fprintf(stderr, "Error: sequence type is error. Program will exit.\n"); exit(1); }
	arg_string += " ";
	if(align_profiles_config_.change_penalty)
	{
		if(align_profiles_config_.gapopen)
		{
			arg_string += align_profiles_config_.gapopen;
			arg_string += " ";
		}
		if(align_profiles_config_.gapext)
		{
			arg_string += align_profiles_config_.gapext;
			arg_string += " ";
		}
	}
	if(! align_profiles_config_.single_thread)
	{
		arg_string += align_profiles_config_.Thread_first;
		arg_string += " ";
		arg_string += thread_str;
		arg_string += " ";
	}
#if ( defined(_WIN32) || defined(_WIN64) )
	return system_spawn_Win((char *)arg_string.c_str(), NULL, FinalFile, NULL);
#elif ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__)) )
	wordexp_t p;
	string_to_argv(arg_string, p, __func__);
	int exit_val = system_spawn_wordexp(align_profiles_config_.prog_path, &p, NULL, FinalFile, NULL);
	wordfree(&p);
	return exit_val;
#endif
}

// Profile-Profile Alignment. 0 is ok, other is bad.
// It will call wmsa profile-profile alignment program.
int MAFFTProfileProfilealignment(const char* ListFile, const char *FinalFile)
{
	std::string arg_string;
	arg_string += align_profiles_config_.prog_path;
	arg_string += " ";
	if (align_profiles_config_.other_args) arg_string += align_profiles_config_.other_args;
	arg_string += " ";
	if(align_profiles_config_.has_two_profile_prefix) arg_string += align_profiles_config_.profile1_prefix;
	arg_string += " \"";
	arg_string += ListFile;
	arg_string += "\" ";
	arg_string += align_profiles_config_.in_arg_str; // --type
	arg_string += " ";
	if(block_sequence_type == 1) arg_string += align_profiles_config_.Protein;
	else if(block_sequence_type == 2) arg_string += align_profiles_config_.DNA;
	else { fprintf(stderr, "Error: sequence type is error. Program will exit.\n"); exit(1); }
	arg_string += " ";
	if(align_profiles_config_.change_penalty)
	{
		if(align_profiles_config_.gapopen)
		{
			arg_string += align_profiles_config_.gapopen;
			arg_string += " ";
		}
		if(align_profiles_config_.gapext)
		{
			arg_string += align_profiles_config_.gapext;
			arg_string += " ";
		}
	}
	if(! align_profiles_config_.single_thread)
	{
		arg_string += align_profiles_config_.Thread_first;
		arg_string += " ";
		arg_string += thread_str;
		arg_string += " ";
	}
	arg_string += align_profiles_config_.out_arg_str;
	arg_string += " \"";
	arg_string += FinalFile;
	arg_string += "\"";
#if ( defined(_WIN32) || defined(_WIN64) )
	return system_spawn_Win((char *)arg_string.c_str(), NULL, NULL, NULL);
#elif ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__)) )
	wordexp_t p;
	string_to_argv(arg_string, p, __func__);
	int exit_val = system_spawn_wordexp(align_profiles_config_.prog_path, &p, NULL, FinalFile, NULL);
	wordfree(&p);
	return exit_val;
#endif
}

// Fragment Alignment. 0 is ok, other is bad.
// It will call fragment alignment subprogram.
int Fragalignment(const char* InFile, const char *OutFile)
{
	// fragalign.exe / fragalign must be in the same folder with program
	const char* real_stdout = NULL;
	std::string arg_string;
	arg_string += align_seqs_config_.prog_path;
	arg_string += " ";
	if (align_seqs_config_.other_args)
	{
		arg_string += align_seqs_config_.other_args;
		arg_string += " ";
	}
	if (block_sequence_type == 1) { if(align_seqs_config_.Protein) arg_string += align_seqs_config_.Protein; }
	else if (block_sequence_type == 2) { if(align_seqs_config_.DNA) arg_string += align_seqs_config_.DNA; }
	else { fprintf(stderr, "Error: sequence type is error. Program will exit.\n"); exit(1); }
	arg_string += " ";
	if(! align_seqs_config_.single_thread)
	{
		arg_string += align_seqs_config_.Thread_first;
		arg_string += " ";
		arg_string += thread_str;
	}
	arg_string += " ";
	if(align_seqs_config_.change_penalty)
	{
		if(align_seqs_config_.gapopen)
		{
			arg_string += align_seqs_config_.gapopen;
			arg_string += " ";
		}
		if(align_seqs_config_.gapext)
		{
			arg_string += align_seqs_config_.gapext;
			arg_string += " ";
		}
	}
	arg_string += " ";
	if (align_seqs_config_.out_to_stdout)
	{
		if (align_seqs_config_.in_arg)
		{
			arg_string += align_seqs_config_.in_arg_str;
			arg_string += " ";
		}
		arg_string += "\"";
		arg_string += InFile;
		real_stdout = OutFile;
		arg_string += "\"";
	}
	else // fit HAlign-3
	{
		arg_string += align_seqs_config_.out_arg_str;
		arg_string += " \"";
		arg_string += OutFile;
		arg_string += "\" ";
		if (align_seqs_config_.in_arg)
		{
			arg_string += align_seqs_config_.in_arg_str;
			arg_string += " ";
		}
		arg_string += "\"";
		arg_string += InFile;
		arg_string += "\"";
	}
#if ( defined(_WIN32) || defined(_WIN64) )
	return system_spawn_Win((char *)arg_string.c_str(), NULL, real_stdout, NULL);
#else
	wordexp_t p;
	string_to_argv(arg_string, p, __func__);
	int exit_val = system_spawn_wordexp(align_seqs_config_.prog_path, &p, NULL, real_stdout, NULL);
	wordfree(&p);
	return exit_val;
#endif
	}

// Two profiles alignment. 0 is ok, other is bad.
// It will call 2-profile alignment subprogram.
int TwoProfileAlignment(const char *Profile1, const char *Profile2, const char* OutFile)
{
	// profile_two_align.exe / profile_two_align must be in the same folder with program
	int exitval;
	std::string arg_string;
	arg_string += align_two_profile_config_.prog_path;
	arg_string += " ";
    if (block_sequence_type == 1) { if(align_two_profile_config_.Protein) arg_string += align_two_profile_config_.Protein; }
	else if (block_sequence_type == 2) { if(align_two_profile_config_.DNA) arg_string += align_two_profile_config_.DNA; }
	else { fprintf(stderr, "Error: sequence type is error. Program will exit.\n"); exit(1); }
	arg_string += " ";
	if(align_two_profile_config_.other_args) arg_string += align_two_profile_config_.other_args;
	arg_string += " ";
	if(align_two_profile_config_.has_two_profile_prefix) arg_string += align_two_profile_config_.profile1_prefix;
	arg_string += " \"";
	arg_string += Profile1;
	arg_string += "\" ";
	if(align_two_profile_config_.has_two_profile_prefix) arg_string += align_two_profile_config_.profile2_prefix;
	arg_string += " \"";
	arg_string += Profile2;
	arg_string += "\" ";
	if(align_two_profile_config_.change_penalty)
	{
		if(align_two_profile_config_.gapopen)
		{
			arg_string += align_two_profile_config_.gapopen;
			arg_string += " ";
		}
		if(align_two_profile_config_.gapext)
		{
			arg_string += align_two_profile_config_.gapext;
			arg_string += " ";
		}
	}
	arg_string += " ";
#if ( defined(_WIN32) || defined(_WIN64) )
	if(align_two_profile_config_.out_to_stdout) return system_spawn_Win((char*)arg_string.c_str(), NULL, OutFile, NULL);
	else
	{
		exitval = system_spawn_Win((char*)arg_string.c_str(), NULL, NULL, NULL);
		exitval |= Filemove(Profile1, OutFile);
		return exitval;
	}
#elif ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__)) )
	wordexp_t p;
	string_to_argv(arg_string, p, __func__);
	if(align_two_profile_config_.out_to_stdout)
		exitval = system_spawn_wordexp(align_two_profile_config_.prog_path, &p, NULL, OutFile, NULL);
	else
	{
		exitval = system_spawn_wordexp(align_two_profile_config_.prog_path, &p, NULL, NULL, NULL);
		exitval |= Filemove(Profile1, OutFile);
	}
	wordfree(&p);
	return exitval;
#endif
}
