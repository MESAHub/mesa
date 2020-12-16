//! ***********************************************************************
//!
//!   Copyright (C) 2018 Robert Farmer, Bill Paxton
//!
//!   MESA is free software; you can use it and/or modify
//!   it under the combined terms and restrictions of the MESA MANIFESTO
//!   and the GNU General Library Public License as published
//!   by the Free Software Foundation; either version 2 of the License,
//!   or (at your option) any later version.
//!
//!   You should have received a copy of the MESA MANIFESTO along with
//!   this software; if not, it is available at the mesa website:
//!   http://mesa.sourceforge.net/
//!
//!   MESA is distributed in the hope that it will be useful,
//!   but WITHOUT ANY WARRANTY; without even the implied warranty of
//!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//!   See the GNU Library General Public License for more details.
//!
//!   You should have received a copy of the GNU Library General Public License
//!   along with this software; if not, write to the Free Software
//!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//!
//! ***********************************************************************

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>   /* mkdir(2) */
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <libgen.h>
#include <sys/types.h>
#include <stdbool.h>
#include <dirent.h>
#include "utils_c_system.h"

static const int SUCCESS=0;
static const int PATHLEN=4096;
static const char TEMPLATE[]="/.temp-XXXXXX";
static const int LEN_TEMP=strlen(TEMPLATE);

/* Makes a single directory at path (mkdir path) */
int c_mkdir(const char * restrict path) {
    
    // If folder already exists return, don't relay on EEXIST as sometimes EACCES happens first
     if(is_dir(path)==1) 
        return SUCCESS;
    
    if (mkdir(path, S_IRWXU) != SUCCESS) {
        if (errno != EEXIST){
            printf("MKDIR error on '%s' Errno %d :: %s\n",path, errno, strerror(errno));
            return -1; 
         }
    }
    return SUCCESS;
}

/* Makes a directory at path, potentially making needed parent directories (mkdir -p path) */
int c_mkdir_p(const char * restrict path) {
    char _path[PATHLEN];
    char *p; 
    
    // If folder already exists return, don't relay on EEXIST as sometimes EACCES happens first
    if(is_dir(path)==1) 
        return SUCCESS;

    strncpy(_path, path,PATHLEN);
    
    /* Iterate the string */
    for (p = _path + 1; *p; p++) {
        if (*p == '/') {
            /* Temporarily truncate */
            *p = '\0';

            if(c_mkdir(_path) != SUCCESS)
                return -1;

            *p = '/';
        }
    }   
    if(c_mkdir(_path) != SUCCESS)
        return -1;

    return SUCCESS;
}


int is_dir(const char * restrict path) {
    DIR* dir = opendir(path);
    if(dir){
        closedir(dir);
        return 1;
    } else {
        return 0;
    }
    
}


/* Moves file src to dest (mv src dest)
 * if src and dest are on the same filesystem, then mv the files
 * if not, copy src to a temp file on same filesystem as dest then mv the temp
 * file to dest, this prevents the filesystem ever being in a inconsistent state.
 * dest either does not exist or does exist and is "full", thus we prevent problems
 * with other programs reading the file before its finished being written.
 */

int c_mv(const char * restrict src, const char * restrict dest) {
    char realSrc[PATHLEN], realDest[PATHLEN];
    char *tmp;
    char destDir[PATHLEN], srcFile[PATHLEN], temp[PATHLEN], copy[PATHLEN];
    struct stat srcStat, destStat, destDirStat;
    int remove_dest, fTemp, index;
    bool dbg = false;
    
    remove_dest = 0;
    /* Get real paths for files */
    realpath(src, realSrc);
        
    /* Cant use realpath as file might not exist yet */
    strncpy(realDest,dest,PATHLEN);
        
    if (stat(realDest,&destStat) == 0)
        if (!S_ISREG(destStat.st_mode))
            strncat(realDest,"/",PATHLEN - strlen(realDest) - 1);
        
    if(dbg){
        printf("mv %s to %s\n",src, dest);
        printf("mv %s to %s\n",realSrc, realDest);
    }
    
    /* Get folder where output will go */
    strncpy(copy, realDest, PATHLEN);
    tmp = dirname(copy);   
    strncpy(destDir, tmp, PATHLEN-1);
    
    /* If destination folder does not end in a / add one */
    index = strlen(destDir);
    if(index > PATHLEN-2)
        return -2;
    
    if(strcmp(&destDir[index],"/") != 0){
        destDir[index+1]  = '/';
        destDir[index+2] = '\0';
    }
    
    /* Get basename of input file */
    strncpy(copy, realSrc, PATHLEN);
    tmp = basename(copy);
    strncpy(srcFile, tmp, PATHLEN-1);
    
    if(dbg)
        printf("+ %s %s %s %s\n",srcFile, realSrc, destDir, realDest);
    
    /* Does input file exist? */ 
    if (stat(realSrc,&srcStat) < 0)
        return -1;
    
    /* Does output exist? */
    if (stat(realDest,&destStat) == 0){
        if(dbg)
            printf("** %s %s\n", realDest, srcFile);
            
        /* Check if folder*/
        if (S_ISREG(destStat.st_mode)){
            /* If a file allready, mark for removal*/
            remove_dest = 1;
        } else {
            /* is folder, append src filename to foldername */
            strcat(realDest, srcFile);
        }
    } 
    if(dbg)
        printf("*** %s %s\n",destDir, realDest);
        
    /* Get stat for output folder */
    if(stat(destDir, &destDirStat))
        return -2;
    
    /* is src and dest on the same file system? */
    if (srcStat.st_dev == destDirStat.st_dev){
        /* Delete output file if it exists */
        if(remove_dest)
            remove(realDest);
            
        if(dbg)
            printf("%s %s\n", realSrc, realDest);
        /* Move src to dest */
        if(rename(realSrc, realDest) != SUCCESS)
            return -3; 
            
    } else{
        /* Copy src to a temp file on dest then move to final location */
        
        /* Build temporay filename */
        strncpy(temp, destDir, PATHLEN - LEN_TEMP);
        strncat(temp, TEMPLATE, PATHLEN - strlen(temp) - 1);
        
        /* Make temp file */
        fTemp = mkstemp(temp);
        if (fTemp < 0)
            goto error;
        /* Just need the file made dont need the file descriptor */
        close(fTemp);
            
        /* Copy data to temp file */
        if (c_cp(realSrc, temp) != SUCCESS)
            goto error; 
            
        /* Delete output file if it exists */
        if(remove_dest)
            remove(realDest);
                        
        /* Move from temp to final dest */
        if(rename(temp, realDest) != SUCCESS)
            goto error;
        
        /* Remove orignal */
        remove(realSrc);
    }
    
    return SUCCESS;
    
    error:
        remove(temp);
        if(fTemp >= 0)
            close(fTemp);
        return -8;
    
}


/* copies file src to dest (cp src dest) */
int c_cp(const char * restrict src, const char * restrict dest)
{
    const int BUF_SIZE=4096;
    int fd_src, fd_dest;
    char buf[BUF_SIZE];
    ssize_t numRead;
    const int dbg=false;

    if(dbg)
        printf("cp %s to %s\n",src,dest);

    fd_src = open(src, O_RDONLY);
    if (fd_src < 0)
        return -1;

    fd_dest = open(dest, O_CREAT | O_WRONLY | O_TRUNC,S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP |S_IROTH);
    if (fd_dest < 0)
        goto error;

    while ((numRead = read(fd_src, buf, BUF_SIZE)) > 0){
        if (write(fd_dest, buf, numRead) != numRead)
            goto error;
    }
    
    if(numRead == -1)
        goto error;
        
    close(fd_src);
    close(fd_dest);
        
    return SUCCESS;

    error:
        close(fd_src);
        if (fd_dest >= 0)
            close(fd_dest);
        return -2;
}
