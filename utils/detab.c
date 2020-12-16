/* detab.c */

#include <stdio.h>
#include <stdlib.h> 
#include <stdbool.h> 
#include <math.h>
#include <string.h>
#include <unistd.h>

#define MAXSIZE 250000
#define TABSPACE 3

   char buffer[MAXSIZE], outbuff[MAXSIZE*TABSPACE];

   static void Detab_File(char *fname) {
      FILE *file;
      int num_chars, i, j, col, num_tabs;
      char c;
      if ((file=fopen(fname,"r"))==NULL) { printf("cannot open %s to read\n",fname); return; }
      num_chars = fread(buffer, sizeof(char), MAXSIZE, file);
      fclose(file);
      if (num_chars == 0) { printf("%s is empty or unreadable", fname); return; }
      if (num_chars >= MAXSIZE) { printf("%s is too large.  max size is %i\n", fname, MAXSIZE); return; }
      num_tabs = 0;
      for (i=0, j=0, col=0; i < num_chars; i++) {
         c = buffer[i];
         if (c == '\t') {
            num_tabs++;
            do { outbuff[j++] = ' '; col++; } while (col % TABSPACE != 0);
            }
         else {
            if (c == '\n') col = 0;
            else col++;
            outbuff[j++] = c;
            }
         }
      if ((file = fopen(fname,"w"))==NULL) { printf("cannot open %s to write\n",fname); return; }
      i = fwrite(outbuff, sizeof(char), j, file);
      fclose(file);
      if (i != j) printf("error in writing %s\n", fname);
      else printf("new %s has %i characters versus %i in old with %i tabs\n", fname, i, num_chars, num_tabs);
      }

   int main(int argc, char *argv[]) {
      int i;
      printf("argc=%i\n", argc);
      for (i = 1; i < argc; i++) {
         printf("%s\n",argv[i]);
         Detab_File(argv[i]);
         }
      return 0;
      }
