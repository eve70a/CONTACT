/* c_test_caddon.c - example program using the CONTACT library (module 1) from C

   Copyright 2008-2023 by Vtech CMCC.
 
   Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main()
{
   int reid, cpid;
   int ver, ierr;
   char *c_outpath; 
   int len_outpath;
#ifdef _WIN32
__declspec(dllimport) :: cntc_initialize
__declspec(dllimport) :: cntc_finalize
#endif

   for (reid = 1; reid<=2; reid++) {
      printf("test_caddon: starting for reid=%d\n", reid);
      ierr = 15031968;
      c_outpath = "obj/";
      len_outpath = strlen(c_outpath);
      cntc_initialize(&reid, &ver, &ierr, c_outpath, &len_outpath);
      printf("test_caddon: obtained ver=%d, ierr=%d\n", ver, ierr);
   }

   for (reid = 1; reid<=2; reid++) {
      printf("test_caddon: finalizing for reid=%d\n", reid);
      cntc_finalize(reid);
   }
   printf("test_caddon: done, exiting.\n");
}

