#ifndef __FILE_CONTROL__
#define __FILE_CONTROL__

#include "program_spawn.h"

#include <cstdarg>

void fscanf_noreturn(FILE* _Stream, const int should_be_value, const char* _Format, ...);

char* getcurrdir();
int makedir(const char* str);
int removedir(const char* str);

int Filecopy(const char* src, const char* dest, bool append);
int Filecopy(const char* src, const char* dest);
int Filedelete(const char* file);
int Filemove(const char* src, const char* dest);
void Filemerge(const char* src, const char* File1, const char* File2, const char* File3);

#endif