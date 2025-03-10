/*
This file is part of the Mingw32 package.
unistd.h maps (roughly) to io.h
*/
#ifndef _UNISTD_H
#define _UNISTD_H

#include <io.h>
#include <process.h>
#include <windows.h>
#include <limits.h>

#ifndef ssize_t
#define ssize_t long int
#endif

#ifndef SSIZE_MAX
#define SSIZE_MAX SIZE_MAX
#endif

/*
* This code is by https://github.com/ivanrad/getline - it works!
*/

#ifndef getdelim
inline ssize_t getdelim(char** lineptr, size_t* n, int delim, FILE* stream) {
    char* cur_pos, * new_lineptr;
    size_t new_lineptr_len;
    int c;

    if (lineptr == NULL || n == NULL || stream == NULL) {
        errno = EINVAL;
        return -1;
    }

    if (*lineptr == NULL) {
        *n = 128; /* init len */
        if ((*lineptr = (char*)malloc(*n)) == NULL) {
            errno = ENOMEM;
            return -1;
        }
    }

    cur_pos = *lineptr;
    for (;;) {
        c = getc(stream);

        if (ferror(stream) || (c == EOF && cur_pos == *lineptr))
            return -1;

        if (c == EOF)
            break;

        if ((*lineptr + *n - cur_pos) < 2) {
            if (SSIZE_MAX / 2 < *n) {
#ifdef EOVERFLOW
                errno = EOVERFLOW;
#else
                errno = ERANGE; /* no EOVERFLOW defined */
#endif
                return -1;
            }
            new_lineptr_len = *n * 2;

            if ((new_lineptr = (char*)realloc(*lineptr, new_lineptr_len)) == NULL) {
                errno = ENOMEM;
                return -1;
            }
            cur_pos = new_lineptr + (cur_pos - *lineptr);
            *lineptr = new_lineptr;
            *n = new_lineptr_len;
        }

        *cur_pos++ = (char)c;

        if (c == delim)
            break;
    }

    *cur_pos = '\0';
    return (ssize_t)(cur_pos - *lineptr);
}

#endif

#ifndef getline

inline ssize_t getline(char** lineptr, size_t* n, FILE* stream) {
    return getdelim(lineptr, n, '\n', stream);
}

#endif
#endif /* _UNISTD_H */
