//! ***********************************************************************
//!
//!   Copyright (C) 2018 Robert Farmer & The MESA Team
//!
//!   This program is free software: you can redistribute it and/or modify
//!   it under the terms of the GNU Lesser General Public License
//!   as published by the Free Software Foundation,
//!   either version 3 of the License, or (at your option) any later version.
//!
//!   This program is distributed in the hope that it will be useful,
//!   but WITHOUT ANY WARRANTY; without even the implied warranty of
//!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//!   See the GNU Lesser General Public License for more details.
//!
//!   You should have received a copy of the GNU Lesser General Public License
//!   along with this program. If not, see <https://www.gnu.org/licenses/>.
//!
//! ***********************************************************************

#include "utils_c_system.h"
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h> /* mkdir(2) */
#include <sys/types.h>
#include <unistd.h>

#define PATHLEN 4096
static const int SUCCESS = 0;
static const char TEMPLATE[] = ".temp-XXXXXX";

/* Makes a single directory at path (mkdir path) */
int c_mkdir(const char *restrict path) {

    // If folder already exists return, don't relay on EEXIST as sometimes
    // EACCES happens first
    if (is_dir(path) == 1) {
        return SUCCESS;
    }

    if (mkdir(path, S_IRWXU) != SUCCESS) {
        if (errno != EEXIST) {
            printf("MKDIR error on '%s' Errno %d :: %s\n", path, errno,
                   strerror(errno));
            return -1;
        }
    }
    return SUCCESS;
}

/* Makes a directory at path, potentially making needed parent directories
 * (mkdir -p path) */
int c_mkdir_p(const char *restrict path) {
    char _path[PATHLEN];
    char *p;

    // If folder already exists return, don't relay on EEXIST as sometimes
    // EACCES happens first
    if (is_dir(path) == 1) {
        return SUCCESS;
    }

    if (strlen(path) >= PATHLEN) {
        return -1;
    }

    strcpy(_path, path);

    /* Iterate the string */
    for (p = _path + 1; *p; p++) {
        if (*p == '/') {
            /* Temporarily truncate */
            *p = '\0';

            if (c_mkdir(_path) != SUCCESS) {
                return -1;
            }

            *p = '/';
        }
    }
    if (c_mkdir(_path) != SUCCESS) {
        return -1;
    }

    return SUCCESS;
}

int is_dir(const char *restrict path) {
    DIR *dir = opendir(path);
    if (dir) {
        closedir(dir);
        return 1;
    } else {
        return 0;
    }
}

/* Moves file src to dest (mv src dest)
 * if src and dest are on the same filesystem, then mv the files
 * if not, copy src to a temp file on same filesystem as dest then mv the temp
 * file to dest, this prevents the filesystem ever being in a inconsistent
 * state. dest either does not exist or does exist and is "full", thus we
 * prevent problems with other programs reading the file before its finished
 * being written.
 */

int c_mv(const char *restrict src, const char *restrict dest) {
    char dest_temp[PATHLEN];

    if (rename(src, dest) == SUCCESS) {
        return SUCCESS;
    }

    if (errno != EXDEV) {
        return -1;
    }

    if (strlen(dest) >= PATHLEN) {
        return -1;
    }

    strcpy(dest_temp, dest);

    size_t i = strlen(dest);
    while (i) {
        if (dest_temp[i - 1] == '/') {
            break;
        }

        i--;
    }

    if (strlen(TEMPLATE) >= PATHLEN - i) {
        return -1;
    }

    strcpy(dest_temp + i, TEMPLATE);

    int fd = mkstemp(dest_temp);

    if (fd == -1) {
        return -1;
    }

    close(fd);

    if (c_cp(src, dest_temp) != SUCCESS) {
        remove(dest_temp);
        return -1;
    }

    if (rename(dest_temp, dest) != SUCCESS) {
        remove(dest_temp);
        return -1;
    }

    remove(src);

    return SUCCESS;
}

/* copies file src to dest (cp src dest) */
int c_cp(const char *restrict src, const char *restrict dest) {
    const int BUF_SIZE = 4096;
    int fd_src, fd_dest;
    char buf[BUF_SIZE];
    ssize_t numRead;
    const int dbg = false;

    if (dbg) {
        printf("cp %s to %s\n", src, dest);
    }

    fd_src = open(src, O_RDONLY);
    if (fd_src < 0) {
        return -1;
    }

    fd_dest = open(dest, O_CREAT | O_WRONLY | O_TRUNC,
                   S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);
    if (fd_dest < 0) {
        goto error;
    }

    while ((numRead = read(fd_src, buf, BUF_SIZE)) > 0) {
        if (write(fd_dest, buf, numRead) != numRead) {
            goto error;
        }
    }

    if (numRead == -1) {
        goto error;
    }

    close(fd_src);
    close(fd_dest);

    return SUCCESS;

error:
    close(fd_src);
    if (fd_dest >= 0) {
        close(fd_dest);
    }
    return -2;
}
