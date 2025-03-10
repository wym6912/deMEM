#ifndef PROGRAM_SPAWN__
#define PROGRAM_SPAWN__
#ifdef __cplusplus
extern "C"
{
#endif //__cplusplus

#if ( defined(_WIN32) || defined(_WIN64) )
#include <process.h>
#include <windows.h>
#include <io.h>
#include <direct.h>
#elif ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__)) )
#include <spawn.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <wait.h>
#include <errno.h>
#include <fcntl.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <wordexp.h>
#endif

#include <stdio.h>

#if ( defined(_WIN32) || defined(_WIN64) ) // Windows System Call Functions
// Call subprograms, return the subprogram exit value (make the program and args in a string)
int system_spawn_Win(char* prog_with_arg, const char* stdin_stream, const char* stdout_stream, const char* stderr_stream);
#elif ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__)) ) // Linux / MacOS Call Functions (posix_spawn mode)
// Call subprograms, return the subprogram exit value
int system_spawn(const char* prog, const char** argv, const char* stdin_stream, const char* stdout_stream, const char* stderr_stream);
int system_spawn_wordexp(const char* prog, wordexp_t *wt, const char* stdin_stream, const char* stdout_stream, const char* stderr_stream);
#endif

#ifdef __cplusplus
}
#endif // __cplusplus

#endif