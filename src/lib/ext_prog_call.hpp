#ifndef __EXT_PROG_CALL__
#define __EXT_PROG_CALL__

#define MAX_ARGS 100

extern "C"
{
#include "program_spawn.h"
}

#include "ext_prog_config.hpp"
#include "libfile_v2.hpp"
#include <string>
#include <vector>

int TwoProfileAlignment(const char *Profile1, const char *Profile2, const char* OutFile);
int WMSAProfileProfilealignment(const char *CenterFile, const char* ListFile, const char *FinalFile);
int MAFFTProfileProfilealignment(const char* ListFile, const char *FinalFile);
int Fragalignment(const char* InFile, const char *OutFile);

#endif