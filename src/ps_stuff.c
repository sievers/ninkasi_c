#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>


#include "ninkasi.h"
#include "ps_stuff.h"



/*--------------------------------------------------------------------------------*/
void *psAlloc(size_t nelem)
{
  return malloc(nelem);
}
/*--------------------------------------------------------------------------------*/
void psFree(void *vec)
{
  if (vec)
    free(vec);
  //else
  // fprintf(stderr,"Warning - skipping free of null vector.\n");
}

/*--------------------------------------------------------------------------------*/
void *psRealloc(void *ptr, size_t size)
{
  return realloc(ptr,size);
}

/*--------------------------------------------------------------------------------*/
int psFindTag(char **tags, int ntag, char *tag)
{
  if (ntag==0)
    return NO_TAG;
  for (int i=0;i<ntag;i++) {
    int len1=strlen(tag);
    int len2=strlen(tags[i]);
    if (len1==len2)
      if (strncmp(tag,tags[i],len1)==0)
	return i;
  }
  return NO_TAG;

}
/*--------------------------------------------------------------------------------*/
void psAddTag(char ***tags, int *ntag, int **level, char *tag) 
{
  *tags=(char **)realloc(*tags,sizeof(char *)*(*ntag+1));
  *level=(int *)realloc(*level,sizeof(int)*(*ntag+1));

  int len=strlen(tag);
  char *tag_copy=(char *)malloc(sizeof(char)*(len+2));


  assert(tag_copy);
  strncpy(tag_copy,tag,len+1);
  (*tags)[*ntag]=tag_copy;
  (*level)[*ntag]=0; //default trace level is 0 unles set otherwise.
  (*ntag)++;
}
/*--------------------------------------------------------------------------------*/
int psTraceLevel(char *tag, int level, bool set_level)
//hold staticness etc. for pstrace.should be thread-safe
{
  static char **tags;
  static int *levels; 
  static int ntag=0;
  int curtag;
#pragma omp critical (PSTRACE)
  {
    curtag=psFindTag(tags,ntag,tag);
    if (curtag==NO_TAG) {
      psAddTag(&tags,&ntag,&levels, tag);
      curtag=ntag-1;
    }
    if (set_level)
      levels[curtag]=level;	
    
  } 
  return levels[curtag];
}
/*--------------------------------------------------------------------------------*/
void psTraceSetLevel(char *tag, int level)
{
  psTraceLevel(tag,level,true);
}

/*--------------------------------------------------------------------------------*/
int psTraceGetLevel(char *tag)
{
  return psTraceLevel(tag, 0, false);
}
/*--------------------------------------------------------------------------------*/
void psTrace(char *tag, int tracelevel, char *format, ...)
{
  va_list args;
  va_start(args,format);
  if (tracelevel<=psTraceGetLevel(tag))
    vfprintf(stdout,format,args);
  va_end(args);
  
}

/*--------------------------------------------------------------------------------*/
actData **psAllocMatrix(int n, int m)
{
  return matrix(n,m);
}

/*--------------------------------------------------------------------------------*/
void psFreeAndZero(void **p)
{
  if ((*p)!=NULL) {
    psFree(*p);
    *p=NULL;
  }
    
}
/*--------------------------------------------------------------------------------*/
void psError(int code, bool whatever, char *format,...)
{
  va_list args;
  va_start(args,format);
  vfprintf(stdout,format,args);
  

}
