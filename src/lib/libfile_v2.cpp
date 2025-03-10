#include "libfile_v2.hpp"

void Filecopy(FILE* src, FILE* des)
{
	if (src == NULL || des == NULL) { fprintf(stderr, "Error: File can not open when copying files"); exit(1); }
	char c;
	while ((c = fgetc(src)) != EOF) fputc(c, des);
}

int Filecopy(const char* src, const char* dest, bool append)
{
	FILE* f1, * f2;
	f1 = fopen(src,  "r");
	f2 = fopen(dest, append ? "w+" : "w");
	Filecopy(f1, f2);
	fclose(f1);
	fclose(f2);
	return 0;
}

int Filecopy(const char* src, const char* dest)
{
#if ( defined(_WIN32) || defined(_WIN64) )
	return !CopyFile(src, dest, true);
#elif ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__)) )
	int fdr = open(src, O_RDONLY);
	int fdw = open(dest, O_WRONLY | O_CREAT, 0600);

	if (fdr == -1 || fdw == -1)
	{
		perror("file open error");
	}

	char buff[256] = { 0 };
	int len = 0;

	while ((len = read(fdr, buff, 256)) > 0)
	{
		write(fdw, buff, len);
	}

	close(fdr);
	close(fdw);
	return 0;
#endif

}

int Filedelete(const char* file)
{
#if ( defined(_WIN32) || defined(_WIN64) )
	return !DeleteFile(file);
#elif ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__)) )
	return remove(file);
#endif
	
}

int Filemove(const char* src, const char* dest)
{
#if ( defined(_WIN32) || defined(_WIN64) )
	return !MoveFileEx(src, dest, MOVEFILE_REPLACE_EXISTING);
#elif ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__)) )
	return rename(src, dest);
#endif

}

void Filemerge(const char* src, const char* File1, const char* File2, const char* File3)
{
	// need modify: the order is File1-File2-File3, I want to use the order is same as sequence id
	FILE* f_read, * f1, * f2, * f3;
	f_read = fopen(src, "w");
	f1 = fopen(File1, "r");
	f2 = fopen(File2, "r");
	f3 = fopen(File3, "r");
	if (f1 == NULL || f2 == NULL || f3 == NULL || f_read == NULL)
	{
		fprintf(stderr, "Error: can not open file. Program will exit.\n");
		exit(1);
	}
	Filecopy(f1, f_read);
	Filecopy(f2, f_read);
	Filecopy(f3, f_read);
	fclose(f1);
	fclose(f2);
	fclose(f3);
	fclose(f_read);
}

int makedir(const char *str)
{
#if ( defined(_WIN32) || defined(_WIN64) )
	return !CreateDirectory(str, NULL);
#elif ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__)) )
	return mkdir(str, 0777);
#endif
}

int removedir(const char* str)
{
#if ( defined(_WIN32) || defined(_WIN64) )
	SHFILEOPSTRUCT FileOp = {0};
	FileOp.wFunc = FO_DELETE;
	FileOp.pFrom = (PCZZSTR)str;
	FileOp.pTo = NULL;
	FileOp.fFlags = FOF_NOCONFIRMATION;
	return SHFileOperation(&FileOp);
#elif ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__)) )
    char* argv[5];
    char *cmddd = (char*)malloc((strlen(str) + 10) * sizeof(char));
    int exitval = 0;
    if(cmddd == NULL) { fprintf(stderr, "Error: can not execute removedir.\n"); return 1; }
    sprintf(cmddd, "rm -rf %s", str);
    argv[0] = "sh"; argv[1] = "-c"; argv[2] = cmddd; argv[3] = NULL;
    exitval = system_spawn("/bin/sh", (const char**)argv, NULL, NULL, NULL);
    free(cmddd);
    return exitval;
#endif
}

char* getcurrdir()
{
#if ( defined(_WIN32) || defined(_WIN64) )
    return _getcwd(NULL, 0);
#elif ( defined(__linux__) || (defined(__APPLE__) && defined(__MACH__)) )
    return getcwd(NULL, 0);
#endif
}

void fscanf_noreturn(FILE* _Stream, const int should_be_value, const char* _Format, ...)
{
	int _Result;
	va_list _ArgList;
	va_start(_ArgList, _Format);
	_Result = vfscanf(_Stream, _Format, _ArgList);
	va_end(_ArgList);
	if (_Result != should_be_value) { fprintf(stderr, "Error: input format error. Program will exit.\n"); exit(1); }
}