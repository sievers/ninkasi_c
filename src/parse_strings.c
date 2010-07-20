#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "parse_strings.h"

#define NOT_FOUND -1

//this is the maximum length of a single line when reading from stdin/a file
#define PARSE_STRINGS_MAX_STRING_LEN 2048

int isnum_float(char c)
     /*return 1 if c is 0..9, -, or . */
{

  return( (isdigit(c)||(c=='.')||(c=='-')));

}

/*-------------------------------------------------------------------*/
int isnum_int(char c)
     /*return 1 if c is 0..9 or -*/
{

  return( (isdigit(c)||(c=='-')));

}

/*-------------------------------------------------------------------*/
#if 0
char **scmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	char **m;

	/* allocate pointers to rows */
	m=(char **) malloc((size_t)((nrow+NR_END)*sizeof(char*)));
	assert(m);
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(char *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(char)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}
#endif
/*-------------------------------------------------------------------*/
char *scvector(long nrl, long nrh)
{
  long i;
  char *cvec;
  cvec = (char *)malloc((size_t)(nrh-nrl+1)*sizeof(char));
  assert(cvec);

  
  return cvec-nrl;
}

/*-------------------------------------------------------------------*/
char **scmatrix_vector(long nrl, long nrh)
{
  long i;
  char **cvec;
  cvec = (char **)malloc((size_t)(nrh-nrl+1)*sizeof(char *));
  assert(cvec);
  
  return cvec-nrl;
}


/*-------------------------------------------------------------------*/

int is_char_in_string(char c, char *delim)
{
  int i,len;
  len=strlen(delim);

  for (i=0;i<len;i++)
    if (c==delim[i])
      return 1;

  return 0;
}

/*-------------------------------------------------------------------*/
void cut_trailing_newline(char *string_in)
{
  if (string_in[strlen(string_in)-1]=='\n')
    string_in[strlen(string_in)-1]='\0';
}

/*-------------------------------------------------------------------*/
#if 1
void free_argv(int argc, char **argv)
//free an argv where each line was allocated separately.
{
  //don't free if argv is empty
  if (argv) {
    for (int i=0;i<argc;i++)
      free(argv[i]);
  }
  free(argv);
}
#else
void free_argv(char **argv)
{
  free(argv[0]);
  free(argv);
}
#endif

/*--------------------------------------------------------------------------------*/
char *copy_string(char *str)
{
  int len=strlen(str);
  char *mystr=(char *)malloc(sizeof(char)*(len+1));
  strcpy(mystr,str);
  return mystr;
}

/*-------------------------------------------------------------------*/
char **create_argv_new(char *line_in, int *narg, char *delim_in)
{
  char *line=copy_string(line_in);
  char *delim=copy_string(delim_in);
  //printf("strings are copies.\n");
  int ntok=0;
  char *saveptr;
  char **argv=NULL;
  //printf("calling strtok_r.\n");
  char *tok=(char *)strtok_r(line,delim,&saveptr);
  assert(tok);
  //printf("called.\n");
  while (tok) {    
    //printf("reallocing.\n");
    argv=(char **)realloc(argv,(ntok+1)*sizeof(char *));
    //printf("realloced.\n");
    assert(argv);
    //printf("strlen is %ld\n",strlen(tok));
    argv[ntok]=copy_string(tok);
    ntok++;

    tok=(char *)strtok_r(NULL,delim,&saveptr);
  }
  
  free(line);
  free(delim);
  *narg=ntok;
  return argv;
}

  
  
  



/*-------------------------------------------------------------------*/
char **create_argv(char *line, int *narg, char *delim)
{
  char *cur_ptr, *head_ptr,**local_argv;
  int i,j,maxlen,i_max,totlen,curlen,ntok,last_good;

  totlen=strlen(line);
  /*fprintf(stderr,"Line length is %d\n",totlen);*/
  head_ptr=scvector(0,strlen(line)+2);
  strcpy(head_ptr,line);
  cur_ptr=head_ptr;
  /*replace all delimiter strings with the null string*/
  for (i=0;i<totlen;i++)
    {
      if (is_char_in_string(head_ptr[i],delim))
	head_ptr[i]='\0';
    }
  maxlen=0;
  ntok=0;
  last_good=0;
  curlen=0;
  /*do we start with a string or a delimiter*/
  if (head_ptr[0])
    {
      ntok++;
      last_good=1;
      curlen=1;
    }
  for (i=1;i<totlen;i++)
    {
      /*fprintf(stderr,"%3d %3d %3d %3d %c  %6d\n",i,curlen,ntok,head_ptr[i],head_ptr[i],maxlen);*/
      if (head_ptr[i])
	{
	  curlen++;
	  last_good=1;
	}
      else
	{
	  if (last_good)
	    {
 	      if (curlen>maxlen)
		maxlen=curlen;
	    }
	  curlen=0;
	  last_good=0;
	}
      if ((head_ptr[i])&&(!head_ptr[i-1]))
	ntok++;

    }
  local_argv=scmatrix_vector(0,ntok);
  i=0;
  if (head_ptr)
    {
      local_argv[i]=head_ptr;
      i++;
    }
  if (head_ptr[0])
    {
      local_argv[i]=head_ptr;
      i++;
    }
  for (j=1;j<totlen;j++)
    {
      if ((head_ptr[j])&&(!head_ptr[j-1]))
	{
	  local_argv[i]=head_ptr+j;
	  i++;
	}
    }
  
  *narg=ntok;
  return local_argv;



}
      
  
/*--------------------------------------------------------------------------------*/
int find_tok(int argc, char **argv, char *to_find) 
/* Look for a token inside of argv matching to_find, and return a pointer
   to it.  It is suffient for the first strlen(to_find) characters to match. */
{
  int len=strlen(to_find);
  for (int i=0;i<argc; i++) {
    int ll=strlen(argv[i]);
    if (ll>=len) {
      if (strncmp(argv[i],to_find,len)==0)
	return i; 
    }   
  }
  return NOT_FOUND;
}

/*--------------------------------------------------------------------------------*/
char *find_argument(int argc, char **argv, char *to_find, int *found_list) 
{
  int i=find_tok(argc,argv,to_find);
  if (i==-1) 
    return NULL;
  
  found_list[i]=1;
  //check if we match the token exactly
  if (strlen(argv[i])==strlen(to_find)) { 
    if (i<argc-1) {
      found_list[i+1]=1;
      return argv[i+1];
    }
    else {
      fprintf(stderr,"Warning - token %s requested with an argument, but none was supplied.\n",to_find);
      return NULL;	
    }
  }
  else {  //the option is later in the same string
    return &(argv[i][strlen(to_find)]);
  }
  
}
/*--------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
int get_command_line_int(int argc, char *argv[], char *match, int *val, int *found_list)
     /*find integer value in command line marked by *match. */
{
  int i;
  char *value_start;


  i=1;
  while ((i<=argc)&&(strncmp(argv[i],match,strlen(match))!=0))
    {
      /*fprintf(stderr,"Didn't find match in %d %s\n",i,argv[i]);*/
      i++;
    }

  /*if we didn't find the match string, return 1. */
  if (i==argc+1)
    return 1;
	
  found_list[i]=1;
  /*check to see if token is longer than the match string.
    if so, assume value is pasted on to end of string.
  */

  /*fprintf(stderr,"Found a match %s\n",argv[i]);*/
  if (strlen(argv[i])>strlen(match))
    value_start=argv[i]+strlen(match);
  else
    {
      /*check to see if the match string is the very end of the input.
	If so, fail to find a value.*/
      if (i==argc)
	return 1;

      /*otherwise, value is contained in the next string line.*/
      value_start=argv[i+1]; 
      if (isnum_int(value_start[0]))
	found_list[i+1]=1;
    }

  /*fprintf(stderr,"value_start is %s\n",value_start);*/
  if (isnum_int(value_start[0]))
    {
      /*if yes, we have found an integer at the start of the string.*/
      *val=atoi(value_start);
    }
  else
    {
      /*we found a non-numeric character where we expected the value.*/
      return 2;
    }

  return 0; /*we got one!*/

}

/*-------------------------------------------------------------------*/

int get_command_line_long(int argc, char *argv[], char *match, long *val, int *found_list)
     /*find integer value in command line marked by *match. */
{
  int i;
  char *value_start;


  i=1;
  while ((i<=argc)&&(strncmp(argv[i],match,strlen(match))!=0))
    {
      /*fprintf(stderr,"Didn't find match in %d %s\n",i,argv[i]);*/
      i++;
    }

  /*if we didn't find the match string, return 1. */
  if (i==argc+1)
    return 1;
	
  found_list[i]=1;
  /*check to see if token is longer than the match string.
    if so, assume value is pasted on to end of string.
  */

  /*fprintf(stderr,"Found a match %s\n",argv[i]);*/
  if (strlen(argv[i])>strlen(match))
    value_start=argv[i]+strlen(match);
  else
    {
      /*check to see if the match string is the very end of the input.
	If so, fail to find a value.*/
      if (i==argc)
	return 1;

      /*otherwise, value is contained in the next string line.*/
      value_start=argv[i+1];
      if (isnum_int(value_start[0]))
	found_list[i+1]=1;
    }

  /*fprintf(stderr,"value_start is %s\n",value_start);*/
  if (isnum_int(value_start[0]))
    {
      /*if yes, we have found an integer at the start of the string.*/
      *val=atoi(value_start);
    }
  else
    {
      /*we found a non-numeric character where we expected the value.*/
      return 2;
    }
  return 0; /*we got one!*/
}


/*-------------------------------------------------------------------*/

int get_command_line_float(int argc, char *argv[], char *match, float *val, int *found_list)
     /*find integer value in command line marked by *match. */
{
  int i;
  char *value_start;


  i=1;
  while ((i<=argc)&&(strncmp(argv[i],match,strlen(match))!=0))
    i++;
  
  /*if we didn't find the match string, return 1. */
  if (i==argc+1)
    return 1;
	
  found_list[i]=1;

  /*check to see if token is longer than the match string.
    if so, assume value is pasted on to end of string.
  */
  if (strlen(argv[i])>strlen(match))
    value_start=argv[i]+strlen(match);
  else
    {
      /*check to see if the match string is the very end of the input.
	If so, fail to find a value.*/
      if (i==argc)
	return 1;

      /*otherwise, value is contained in the next string line.*/
      value_start=argv[i+1];
      if (isnum_float(value_start[0]))
	found_list[i+1]=1;
    }
  
  /*fprintf(stderr,"value_start is %s\n",value_start);*/
  if (isnum_float(value_start[0]))
    {
      /*if yes, we have found an integer at the start of the string.*/
      *val=atof(value_start);
    }
  else
    {
      /*we found a non-numeric character where we expected the value.*/
      return 2;
    }
  return 0; /*we got one!*/
}



/*-------------------------------------------------------------------*/

int get_command_line_double(int argc, char *argv[], char *match, double *val, int *found_list)
     /*find integer value in command line marked by *match. */
{
  int i;
  char *value_start;


  i=1;
  while ((i<=argc)&&(strncmp(argv[i],match,strlen(match))!=0))
    i++;
  
  /*if we didn't find the match string, return 1. */
  if (i==argc+1)
    return 1;
	
  found_list[i]=1;

  /*check to see if token is longer than the match string.
    if so, assume value is pasted on to end of string.
  */
  if (strlen(argv[i])>strlen(match))
    value_start=argv[i]+strlen(match);
  else
    {
      /*check to see if the match string is the very end of the input.
	If so, fail to find a value.*/
      if (i==argc)
	return 1;

      /*otherwise, value is contained in the next string line.*/
      value_start=argv[i+1];
      if (isnum_float(value_start[0]))
	found_list[i+1]=1;
    }
  
  /*fprintf(stderr,"value_start is %s\n",value_start);*/
  if (isnum_float(value_start[0]))
    {
      /*if yes, we have found an integer at the start of the string.*/
      *val=atof(value_start);
    }
  else
    {
      /*we found a non-numeric character where we expected the value.*/
      return 2;
    }
  return 0; /*we got one!*/
}




/*-------------------------------------------------------------------*/

int get_command_line_string(int argc, char *argv[], char *match, char *val, int *found_list)
     /*find integer value in command line marked by *match. */
{
  int i;
  char *value_start;


  i=1;
  while ((i<=argc)&&(strncmp(argv[i],match,strlen(match))!=0))
    i++;
  
  /*if we didn't find the match string, return 1. */
  if (i==argc+1)
    return 1;
	
  found_list[i]=1;

  /*check to see if token is longer than the match string.
    if so, assume value is pasted on to end of string.
  */
  if (strlen(argv[i])>strlen(match))
    value_start=argv[i]+strlen(match);
  else
    {
      /*check to see if the match string is the very end of the input.
	If so, fail to find a value.*/
      if (i==argc)
	return 1;

      /*otherwise, value is contained in the next string line.*/
      value_start=argv[i+1];
      found_list[i+1]=1;
    }
  
  strcpy(val,value_start);

  return 0; /*we got one!*/
}


/*-------------------------------------------------------------------*/
int exists_in_command_line(int argc, char *argv[], char *match, int *found_list)
     /*find integer value in command line marked by *match. */
{
#if 1
  int i=find_tok(argc,argv,match);
  if (i==NOT_FOUND) {
    return 0;
  }
  else {
    found_list[i]=1;
    return 1;
  }
#else
  int i;
  char *value_start;

  i=1;
  while ((i<=argc)&&(strncmp(argv[i],match,strlen(match))!=0))
    i++;
  
  /*if we didn't find the match string, return 1. */
  if (i==argc+1)
    return 1;
	
  found_list[i]=1;

  return 0; /*we have entry in command line.*/
#endif  
}

/*-------------------------------------------------------------------*/
int where_in_command_line(int argc, char *argv[], char *match)
     /*find integer value in command line marked by *match. */
{
  int i;
  char *value_start;

  i=1;
  while ((i<=argc)&&(strncmp(argv[i],match,strlen(match))!=0))
    i++;
  
  /*if we didn't find the match string, return 1. */
  if (i==argc+1)
    return -1;
  return i; /*we have entry in command line.*/
}

/*-------------------------------------------------------------------*/
int exists_in_string(char *string_in, char *tag)
     /*Check to see if tag exists in string_in*/
{
  int len1, len2, i;
  len1=strlen(string_in);
  len2=strlen(tag);

  for (i=0;i<=len1-len2;i++)
    {
      /*if (strncmp(&string_in[i],tag,len2)==0)
	return 1;*/
      /*return the starting position of the string*/
      if (strncmp(&string_in[i],tag,len2)==0)
	return i+1;

    }
  return 0;
}

/*--------------------------------------------------------------------------------*/
char *whitespace_pad_char(char *str, char *delim)
//add whitespace around character delim, return a malloced copy of the string.
{
  int ndelim=0;
  int len=strlen(str);
  for (int i=0;i<len;i++)
    if (strchr(delim,str[i]))
      ndelim++;
  if (ndelim==0)
    return copy_string(str);
  
  
  int totlen=len+2+2*ndelim;
  char *dest=(char *)malloc(sizeof(char)*(totlen));
  int j=0;
  for (int i=0;i<len;i++) {
    if (strchr(delim,str[i])) {
      dest[j]=' ';
      j++;
      dest[j]=str[i];
      j++;
      dest[j]=' ';
      j++;
    }
    else {
      dest[j]=str[i];
      j++;
    }
  }
  dest[j]='\0';
  return dest;
}

/*-------------------------------------------------------------------*/
int expand_tokens(char *string_out, char *string_in,char *to_replace, char **replacement_strings, int n_to_replace, int max_len)
{
  /*
    routine to replace all characters that exist in to_replace with their 
    expansions in replacement_strings.  Replacement happens in string_in
    and is put in string_out.  maximum allowable length is in max_len.
  */

  int i,j,k,n,needs_replacing,hit_end;

  j=0;
  hit_end=0;
  for (i=0;i<strlen(string_in);i++)
    {
      needs_replacing=0;
      for (k=0;k<n_to_replace;k++)
	{
	  if(string_in[i]==to_replace[k])
	    {
	      needs_replacing=1;
	      for (n=0;n<strlen(replacement_strings[k]);n++)
		{
		  j++;
		  if (j==max_len-1)
		    {
		      string_out[j]='\0';
		      return 1;
		    }
		  string_out[j]=replacement_strings[k][n];
		}
	    }
	}
      if (needs_replacing==0)
	{
	  if (j==max_len-1)
	    {
	      string_out[j]='\0';
	      return 1;
	    }
	  j++;
	  string_out[j]=string_in[i];
	}
    }
  j++;
  string_out[j]='\0';
  return 0;
}


/*--------------------------------------------------------------------------------*/
void lengthen_string(char **cur_string,int *len_in, int n_to_add)
{
  char *new_string;
  int len;

  len=*len_in;
  new_string=scvector(0,len+n_to_add);
  strncpy(new_string,*cur_string,len);

  if (len>0)
    free(*cur_string);
  *len_in = len+n_to_add;
  *cur_string=new_string;
}

/*--------------------------------------------------------------------------------*/

char *read_all_stdin(const char *fname)
{
  /*pass in null to read stdin*/
  FILE *instream;
  if (fname) {
    instream=fopen(fname,"r");
    assert(instream);
  }
  else
    instream=stdin;
  char line_in[PARSE_STRINGS_MAX_STRING_LEN],*big_line,*cur_spot;
  int line_len,total_len,i,finished;

  total_len=1;
  big_line=scvector(0,1);
  //sprintf(big_line,"");
  big_line[0]='\0';
  
  finished=0;
  while ((fgets(line_in,PARSE_STRINGS_MAX_STRING_LEN-1,instream))&&(finished==0))
    {
      /*fprintf(stderr,"line_in is %s\n",line_in);*/
      if (strncmp(line_in,"#finished",strlen("#finished")-1)==0)
	finished=1;
      for (i=0;i<(int)strlen(line_in);i++)
	if ((line_in[i]=='#')||(line_in[i]=='!')||(line_in[i]=='\n'))
	  line_in[i]='\0';      
      /*fprintf(stderr,"line_in is .%s.\n",line_in);*/
      lengthen_string(&big_line,&total_len,strlen(line_in)+3);
      sprintf(big_line,"%s %s",big_line,line_in);
    }
  /*if we read from a file, close it here.*/
  if (fname)
    fclose(instream);
  return big_line;
}
/*--------------------------------------------------------------------------------*/
int char_set_from_argv(int argc, char **argv,char *tag, char *stoptag, int *found_list, char *destination, int len)
{
  int i,which_argc,nguess;
  which_argc=where_in_command_line(argc,argv,tag);
  if (which_argc<1)
    {
      fprintf(stderr,"Tag %s not found in char_set_from_argv.\n",tag);
      return -1;
    }
  nguess=where_in_command_line(argc-which_argc,&argv[which_argc],stoptag)-1;
  if (nguess<0)
    nguess=argc-which_argc;
  
  found_list[which_argc]=1;
  for (i=0;i<nguess;i++)
    {
      sprintf(destination+i*len*sizeof(char),"%s",argv[which_argc+i+1]);
      found_list[which_argc+i+1]=1;
    }
  return nguess;
}

/*--------------------------------------------------------------------------------*/
int double_set_from_argv(int argc, char **argv,char *tag, char *stoptag, int *found_list, double *destination, int len)
{
  int i,which_argc,nguess,nuse;
  which_argc=where_in_command_line(argc,argv,tag);
  if (which_argc<1)
    {
      fprintf(stderr,"Tag %s not found in char_set_from_argv.\n",tag);
      return -1;
    }
  nguess=where_in_command_line(argc-which_argc,&argv[which_argc],stoptag)-1;
  if (nguess<0)
    nguess=argc-which_argc;
  if (nguess>len)
    {
      fprintf(stderr,"Warning - have %d guess tags in %s, but only space for %d.  Truncating....\n",nguess,tag,len);
      nguess=len;
    }
  found_list[which_argc]=1;
  for (i=0;i<nguess;i++)
    {
      destination[i]=atof(argv[which_argc+i+1]);
      found_list[which_argc+i+1]=1;
    }
  return nguess;
}
/*--------------------------------------------------------------------------------*/
int int_set_from_argv(int argc, char **argv,char *tag, char *stoptag, int *found_list, int *destination, int len)
{
  int i,which_argc,nguess,nuse;
  which_argc=where_in_command_line(argc,argv,tag);
  if (which_argc<1)
    {
      fprintf(stderr,"Tag %s not found in char_set_from_argv.\n",tag);
      return -1;
    }
  nguess=where_in_command_line(argc-which_argc,&argv[which_argc],stoptag)-1;
  if (nguess<0)
    nguess=argc-which_argc;
  if (nguess>len)
    {
      fprintf(stderr,"Warning - have %d guess tags in %s, but only space for %d.  Truncating....\n",nguess,tag,len);
      nguess=len;
    }
  found_list[which_argc]=1;
  for (i=0;i<nguess;i++)
    {
      destination[i]=atoi(argv[which_argc+i+1]);
      found_list[which_argc+i+1]=1;
    }
  return nguess;
}

/*--------------------------------------------------------------------------------*/
char *find_first_delim_in_string(char *str, char *delim)
//find the first occurence of any character in delim inside of str, and return a pointer
//to that location in str.  If none exist, return null.
{
  
  int len=strlen(str);
  for (int i=0;i<len;i++) {
    if (strchr(delim,str[i]))
      return str+i;
  }
  return NULL;

}
/*--------------------------------------------------------------------------------*/
char *argv_subset_to_string(int argc, char **argv,int first, int last)
{
  
  assert(first>=0);
  assert(first<argc);
  assert(last>=0);
  assert(last<argc);
  assert(first<=last);
  int totlen=0;
  for (int i=first;i<=last;i++)
    totlen+=strlen(argv[i])+2;
  char *c=(char *)malloc(sizeof(char)*totlen);
  c[0]='\0';
  for (int i=first;i<=last;i++)
    sprintf(c,"%s %s",c,argv[i]);
  return c;
}
/*--------------------------------------------------------------------------------*/
char *get_list_as_string(int argc, char **argv, char *tag, char *start_tag, char *stop_tag, int *found_list)
//pull a list out of an argv/argc delimited by start_tag and stop_tag.  If
//null is passed in for those guys, [] brackets are default.  If the next character after the tag is not a 
//delimiter, then return the next argument in a string copy.
{
  if (!start_tag)
    start_tag="[";
  if (!stop_tag)
    stop_tag="]";
  
  char *start=find_argument(argc, argv, tag, found_list);
  if (start==NULL)
    return NULL;
  printf("start is %s\n",start);

  if (!strchr(start_tag,start[0])) {
    //printf("didn't find list beginning in %s and #%c#\n",tag,start[0]);
    return copy_string(start);
  }
  int first_tok=find_tok(argc,argv,tag);
  
  int last_tok;
  for (int i=first_tok;i<argc;i++) {
    if (find_first_delim_in_string(argv[i],stop_tag)) {
      last_tok=i;
      break;
    }
    if (i==argc-1) {
      fprintf(stderr,"warning - unclosed loop detected in get_list_as_string with requested tag %s and delimiters %s %s\n",tag,start_tag,stop_tag);
      return NULL;
    }
  }
  char *crude_line=argv_subset_to_string(argc,argv,first_tok,last_tok);
  for (int i=first_tok;i<=last_tok;i++)
    found_list[i]=1;


  char *start_ptr=find_first_delim_in_string(crude_line,start_tag);
  assert(start_ptr);
  if (strlen(start_ptr)==0) {
    fprintf(stderr,"Warning - blank loop detected in get_list_as_string  with requested tag %s and delimiters %s %s\n",tag,start_tag,stop_tag);
    free(crude_line);
    return NULL;
  }
  
  start_ptr++;  //move to the first non-delimiter character
  char *stop_ptr=find_first_delim_in_string(start_ptr,stop_tag);
  assert(stop_ptr);  //we should have already failed if this doesn't exist
  stop_ptr[0]='\0';  //null terminate at this point;
  char *good_line=copy_string(start_ptr);
  free(crude_line);
  return good_line;
  
}
/*--------------------------------------------------------------------------------*/
char **get_list_from_argv(int argc, char **argv, char *tag, int *nfound, int *found_list)
{
  char *myptr=get_list_as_string(argc,argv,tag,NULL,NULL,found_list);
  if (myptr==NULL) {
    *nfound=0;
    return NULL;
  }
  char **myargv=create_argv_new(myptr,nfound," \n");
  free(myptr);
  return myargv;
}
/*--------------------------------------------------------------------------------*/
int *get_int_list_from_argv(int argc, char **argv, char *tag, int *nfound, int *found_list)
//pull an integer list out of a parameter file.  Going to use default delimeters here.
//int argc, char **argv, char *tag, char *start_tag, char *stop_tag, int *found_list)
{

  char delim=':';
  char *myptr=get_list_as_string(argc,argv,tag,NULL,NULL,found_list);
  if (!myptr) {
    *nfound=0;
    return NULL;
  }

  printf("my line is %s\n",myptr);
  char *new_ptr=whitespace_pad_char(myptr,&delim);
  free(myptr);

  int myargc;
  char **myargv=create_argv_new(new_ptr,&myargc," ");
  free(new_ptr);
  if (myargc==0) {
    *nfound=0;
    return NULL;
  }



  assert(strchr(myargv[0],delim)==NULL);  //we shouldn't start with a delimeter.  This would be undefined.
  assert(strchr(myargv[myargc-1],delim)==NULL);  //we shouldn't end with one, either.
  int nelem=0;
  for (int i=0;i<myargc;i++) {
    nelem++;
    if (strchr(myargv[i],delim)) {
      int upper=atoi(myargv[i+1]);
      int lower=atoi(myargv[i-1]);
      assert(upper>=lower);
      nelem+=upper-lower-2;  //yeah, I think this math is correct
    }
  }
  int *vals=(int *)malloc(sizeof(int)*nelem);
  int cur=0;
  for (int i=0;i<myargc;i++) {
    if (strchr(myargv[i],delim)) {

      assert(strchr(myargv[i+1],delim)==NULL);  //don't support two colons in a row.
      int jmax=atoi(myargv[i+1]);
      int j=vals[cur-1];
      if (jmax>j) {
	for (j=j+1;j<jmax;j++) {
	  //printf("j and cur are %d %d\n",j,cur);
	  vals[cur]=j;
	  cur++;
	}
      }
      else
	i++;
    }
    else {
      vals[cur]=atoi(myargv[i]);
      cur++;    
    }    
  }
  free_argv(myargc,myargv);
  *nfound=nelem;
  return vals;
}
  

