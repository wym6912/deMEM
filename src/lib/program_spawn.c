#ifdef __cplusplus
extern "C"
{
#endif

#include "program_spawn.h"

#define DEBUG_ARGS 0

#if (defined(_WIN32) || defined(_WIN64))

#define CMD_MAX_LEN 1024

int system_spawn_Win(char* prog_with_arg, const char* stdin_stream, const char* stdout_stream, const char* stderr_stream)
{
    STARTUPINFO si;
    PROCESS_INFORMATION pi;
    DWORD exit_value;
    SECURITY_ATTRIBUTES sa;
    sa.nLength = sizeof(sa);
    sa.lpSecurityDescriptor = NULL;
    sa.bInheritHandle = TRUE;

    ZeroMemory(&si, sizeof(si));
    si.cb = sizeof(si);
    ZeroMemory(&pi, sizeof(pi));
    
    LPCSTR szCmd = "%SystemRoot%\\System32\\cmd.exe /C ", szWritableCmd;
    szWritableCmd = (LPCSTR)calloc(CMD_MAX_LEN, sizeof(CHAR));
    if (szWritableCmd == NULL) { fprintf(stderr, "ERROR: cannot run command. Program will exit\n"); exit(1); }
    ExpandEnvironmentStrings(szCmd, szWritableCmd, CMD_MAX_LEN);
    strcat_s(szWritableCmd, strlen(szWritableCmd) + strlen(prog_with_arg) + 1, prog_with_arg);

    if (stdin_stream != NULL)
    {
        strcat_s(szWritableCmd, strlen(szWritableCmd) + 4, " < ");
        strcat_s(szWritableCmd, strlen(szWritableCmd) + strlen(stdin_stream) + 1, stdin_stream);
    }
    else si.hStdInput = NULL;

    if (stdout_stream != NULL)
    {
        strcat_s(szWritableCmd, strlen(szWritableCmd) + 4, " > ");
        strcat_s(szWritableCmd, strlen(szWritableCmd) + strlen(stdout_stream) + 1, stdout_stream);
    }
    else si.hStdOutput = NULL;

    if (stderr_stream != NULL)
    {
        strcat_s(szWritableCmd, strlen(szWritableCmd) + 4, " 2>");
        strcat_s(szWritableCmd, strlen(szWritableCmd) + strlen(stderr_stream) + 1, stderr_stream);
    }
    else si.hStdError = NULL;

    if(stderr_stream && strcmp(stderr_stream, "nul")) fprintf(stderr, "Run command: %s\n", szWritableCmd);

    // Start the child process. 
    if (!CreateProcess(NULL, szWritableCmd, NULL, NULL, FALSE, 0, NULL, NULL, &si, &pi))
    {
        fprintf(stderr, "CreateProcess failed (%d): ", GetLastError());
        FormatMessage(
            FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
            NULL,
            GetLastError(),
            MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), 
            (LPTSTR)&szWritableCmd,
            0,
            NULL
        );
        fprintf(stderr, "%s", szWritableCmd);
        return -1;
    }

    // Wait until child process exits.
    do
    {
        WaitForSingleObject(pi.hProcess, INFINITE);
        GetExitCodeProcess(pi.hProcess, &exit_value);
    } while (exit_value == STILL_ACTIVE);

    // Close process and thread handles. 
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);

    return (int)exit_value;
}
#elif ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__) ) )
#pragma GCC push_options
#pragma GCC optimize("O0")
int system_spawn(const char* prog, const char** argv, const char *stdin_stream, const char *stdout_stream, const char *stderr_stream)
{
    posix_spawn_file_actions_t action;
    posix_spawn_file_actions_init(&action);
    if(stdin_stream != NULL) posix_spawn_file_actions_addopen(&action, STDIN_FILENO, stdin_stream, O_RDONLY | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if(stdout_stream != NULL) posix_spawn_file_actions_addopen(&action, STDOUT_FILENO, stdout_stream, O_CREAT | O_WRONLY | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if(stderr_stream != NULL) posix_spawn_file_actions_addopen(&action, STDERR_FILENO, stderr_stream, O_WRONLY | O_APPEND, 0);
    pid_t pid;
    int status = posix_spawn(&pid, prog, &action, NULL, argv, NULL), s = 0;
#if DEBUG_ARGS
    while (argv[s]) { fputs(argv[s++], stderr); fputc(' ', stderr); }
    fputc('\n', stderr);
#endif
    // stop until the subprogram ended
    if (status)
    {
        fprintf(stderr, "Error on posix_spawn: %s, calling %s, args = \n", strerror(status), prog);
        while(argv[s]) { fputs(argv[s++], stderr); fputc(' ', stderr); }
        fputc('\n', stderr);
    }
    else
    {
#if 0
        printf("PID of child: %ld, called by %s\n", (long)pid, prog);
#endif
        do 
        {
            s = waitpid(pid, &status, WUNTRACED | WCONTINUED);
#if 0
            if (s == -1)
                fprintf(stderr, "waitpid");

            fprintf(stderr, "Child status: ");
            if (WIFEXITED(status))
            {
                fprintf(stderr, "Normally exited, status = %d\n", WEXITSTATUS(status));
            }
            else if (WIFSIGNALED(status))
            {
                fprintf(stderr, "Killed by signal %d\n", WTERMSIG(status));
            }
            else if (WIFSTOPPED(status))
            {
                fprintf(stderr, "Stopped by signal %d\n", WSTOPSIG(status));
            }
#endif
        } while (!WIFEXITED(status) && !WIFSIGNALED(status));
    }
    posix_spawn_file_actions_destroy(&action);
    return status;
}
#pragma GCC pop_options

int system_spawn_wordexp(const char* prog, wordexp_t *wt, const char* stdin_stream, const char* stdout_stream, const char* stderr_stream)
{
    posix_spawn_file_actions_t action;
    posix_spawn_file_actions_init(&action);
    if(stdin_stream != NULL) posix_spawn_file_actions_addopen(&action, STDIN_FILENO, stdin_stream, O_RDONLY | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if(stdout_stream != NULL) posix_spawn_file_actions_addopen(&action, STDOUT_FILENO, stdout_stream, O_CREAT | O_WRONLY | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if(stderr_stream != NULL) posix_spawn_file_actions_addopen(&action, STDERR_FILENO, stderr_stream, O_WRONLY | O_APPEND, 0);
    pid_t pid;
    int status = posix_spawn(&pid, prog, &action, NULL, (const char**)wt->we_wordv, NULL), s;
#if DEBUG_ARGS
    for(s = 0; s < wt->we_wordc; ++ s) { fputs(wt->we_wordv[s], stderr); fputc(' ', stderr); } fputc('\n', stderr);
#endif
    // stop until the subprogram ended
    if (status)
    {
        fprintf(stderr, "Error on posix_spawn: %s, calling %s, args = \n", strerror(status), prog);
        for(s = 0; s < wt->we_wordc; ++ s) { fputs(wt->we_wordv[s], stderr); fputc(' ', stderr); } fputc('\n', stderr);
    }
    else
    {
        printf("PID of child: %ld, called by %s\n", (long)pid, prog);
        do 
        {
            s = waitpid(pid, &status, WUNTRACED | WCONTINUED);
            if (s == -1)
                fprintf(stderr, "waitpid");

            fprintf(stderr, "Child status: ");
            if (WIFEXITED(status))
            {
                fprintf(stderr, "Normally exited, status = %d\n", WEXITSTATUS(status));
            }
            else if (WIFSIGNALED(status))
            {
                fprintf(stderr, "Killed by signal %d\n", WTERMSIG(status));
            }
            else if (WIFSTOPPED(status))
            {
                fprintf(stderr, "Stopped by signal %d\n", WSTOPSIG(status));
            }
        } while (!WIFEXITED(status) && !WIFSIGNALED(status));
    }
    posix_spawn_file_actions_destroy(&action);
    return status;
}
#endif

#ifdef __cplusplus
}
#endif
