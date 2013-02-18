
/// \file mbCuts.c

/// Contains methods for manipulating a set of time periods that are
/// to be cut from further analysis.  Cuts can be array-wide (global)
/// or detector-specific.


//#include "act.h"
#include <ninkasi.h>
#include "mbCuts.h"
#include <string.h>  // For memset()
#include <assert.h>

#include <omp.h>

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Local (static) method declarations
////////////////////////////////////////////////////////////////////////////////////////////////////////////

static void CutListFree(mbCutList *list);
void CutsFree(mbCuts *cuts);
static inline int MAX(int a, int b);
static void cutsCombine(mbCutList *list);
static void cutListExtend(mbCutList *list, int first, int last);
static bool isThisListAlwaysCut(const mbCutList *list);
static bool isThisListNeverCut(const mbCutList *list);

static mbCutList *cutListOr(const mbCutList *a, const mbCutList *b);



static inline int MAX(int a, int b)
{
  if (a>b)
    return a;
  else
    return b;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Public functions on the structures mbCuts,  mbCutList and mbSingleCut.
////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Allocate a mbSingleCut structure.
/// Also registers the deallocator.
/// \param first  The first data index to be cut.
/// \param last   The last data index to be cut.




void CutsFree(mbCuts *cuts)
{
  for (int r=0; r<cuts->nrow; r++)
    for (int c=0; c<cuts->ncol; c++)
      psFree(cuts->detCuts[r][c]);
  psFree(cuts->detCuts[0]);
  psFree(cuts->detCuts);
  psFree(cuts->globalCuts);
  //psFree(cuts->tod);
}

/*--------------------------------------------------------------------------------*/


int mbCutsRead(mbCuts *c, const char *filename)
{

  FILE *f;
  int nrow, ncol;
  int MAX_LINE = 200000; // global cut could be v. long with turnaround cuts
  char line[MAX_LINE];
  int first, last;
  int nread = 0;

  f = fopen( filename, "r" );
  if (f == NULL) {
    psTrace("moby", 2, "mbCutsRead: Failed to open cut file: %s", filename);
    return nread;
  }

  fgets( line, MAX_LINE, f );
  sscanf( line, "%d %d\n", &nrow, &ncol );
  //printf("nread and ncol are %d %d\n",nrow,ncol);
  nread++;
  if ( c->nrow != nrow || c->ncol != ncol ){
    //psTrace("moby", 2, "nrow and ncol in %s does not fit cuts object", filename);
    fprintf(stderr,"nrow and ncol in %s does not fit cuts object\n", filename);
    return nread;
  }

  { // read in the global cuts
    fgets( line, MAX_LINE, f );
    char *tok = NULL;
    tok = strtok( line, " " );
    int nnread = 0;
    if (tok != NULL) {
      while ( *tok != '\n' ){
        nnread = sscanf( tok, "(%d,%d)", &first, &last );
        if (nnread != 2) {
          psTrace( "moby", 3, "Bad format in global cuts: %s\n", 
                   filename );
          return nnread;
        }
	
        mbCutsExtendGlobal( c, first, last);
        tok = strtok( NULL, " " );
      }
    }
  }
  //printf("reading detector cuts now...\n");

  { // now read in the per detector cuts
    char *tok = NULL;
    int nnread = 0;
    int row, col;
    printf("starting per-detector cuts.\n");
    while( fgets( line, MAX_LINE, f ) != NULL ) {
      nnread = sscanf( line, "r%dc%d: ", &row, &col );
      //printf("parsing %d %d %s\n",row,col,line);
      if (nnread != 2) {
        psTrace( "moby", 3, "Bad cutfile format : %s :\n%s\n",
                 filename, line );
        return nnread;
      }
      tok = strtok( line, " " );
      tok = strtok( NULL, " " ); //skip the row/col designator
      if ( tok == NULL ) continue;
      while ((tok!=NULL)&&( *tok != '\n' )){
        nnread = sscanf( tok, "(%d,%d)", &first, &last );
        if (nnread != 2) {
          psTrace( "moby", 3, "Bad format in global cuts: %s :\n%s\n",
                   filename, line );
          return nnread;
        }
	//printf("extending cuts on %d %d with %d %d\n",row,col,first,last);
        mbCutsExtend( c, first, last, row, col);
        tok = strtok( NULL, " " );
      }
    }
  }
  fclose(f);
  printf("finished cuts read.\n");
  return nread;
}
/*--------------------------------------------------------------------------------*/
void mbCutsDecimateOneCut(mbCutList *cut)
{
  if (cut==NULL)
    return;
  if (isThisListAlwaysCut(cut))
    return;
  int ncut=cut->ncuts;
  if (ncut==0) {
    //fprintf(stderr,"weirdness in mbCutsDecimateOneCut.\n");
    return;
  }
  if (cut==NULL)
    return;

  mbSingleCut *head=cut->head;
  int ncur=0;
  while (head!=NULL) {
    ncur++;
    head->indexFirst= (head->indexFirst)/2;
    head->indexLast= (head->indexLast+1)/2;
    head=head->next;
    if (ncur==ncut)
      if (head) {
	fprintf(stderr,"Cuts have a problem.\n");
	return;
      }
  }
  if (ncur!=ncut)
    fprintf(stderr,"mismatch in mbCutsDecimateOneCut, %d vs %d\n",ncur,ncut);
}
/*--------------------------------------------------------------------------------*/
void mbCutsDecimate(mbCuts *cuts)
{
  //fprintf(stderr,"Doing globals.\n");
  mbCutsDecimateOneCut(cuts->globalCuts);

  for (int row=0;row<cuts->nrow;row++)
    for (int col=0;col<cuts->ncol;col++) {
      if (!mbCutsIsAlwaysCut(cuts,row,col)) {
	//fprintf(stderr,"doing %d  %d with %d cuts.\n",row,col,cuts->detCuts[row][col]);
	mbCutsDecimateOneCut(cuts->detCuts[row][col]);
      }

    }
}
/*--------------------------------------------------------------------------------*/

/// Extend the list of cuts for all detectors
/// \param cuts   The cuts object to be modified.
/// \param first  The first data index to be cut.
/// \param last   The last data index to be cut.
 
void mbCutsExtendGlobal(mbCuts *cuts, int first, int last)
{
  if (cuts == NULL)
    return;

  if (first < 0)
    first = 0;
  if (last < first)
    return;

  cutListExtend(cuts->globalCuts, first, last);
}



/// Extend the list of cuts for a single detector (or create list, if none).
/// \param cuts   The cuts object to be modified.
/// \param first  The first data index to be cut.
/// \param last   The last data index to be cut.
/// \param row    The detector row number.
/// \param col    The detector column number.
 
void mbCutsExtend(mbCuts *cuts, int first, int last, int row, int col)
{
  if (cuts == NULL ||
      row < 0 || row >= cuts->nrow ||
      col < 0 || col >= cuts->ncol)
    return;

  if (first < 0)
    first = 0;
  if (last < first)
    return;

  mbCutList *list = cuts->detCuts[row][col];
  if (list == NULL) {
    list = mbCutListAlloc();
    cuts->detCuts[row][col] = list;
  }
  cutListExtend(list, first, last);
}



/// Extend any mbCutList by adding a new cut on index [first,last] both inclusive.
/// Takes care of sorting the new cut into the appropriate place in the list.
/// \param list   The mbCutList to be extended.
/// \param first  The first data index to be cut.
/// \param last   The last data index to be cut.

static void cutListExtend(mbCutList *list, int first, int last)
{
  assert (first <= last);
  mbSingleCut *newCut = mbSingleCutAlloc(first, last);

  // There is no list yet.  This is going to be easy.
  if (list->head == NULL) {
    list->head = newCut;
    list->ncuts = 1;
    return;
  }

  // Add the current cut to an existing list, ensuring that it goes AFTER the
  // last existing list element with the same or lower firstIndex.
  // Do not worry if newCut overlaps any existing cuts: cutsCombine will fix that.
  mbSingleCut *tail = list->head;
  while (tail->next != NULL)
    tail=tail->next;
  int firstFirst = list->head->indexFirst;
  int lastFirst = tail->indexFirst;

  // Add the newCut at the tail of the list
  if (first >= lastFirst) {       
    tail->next = newCut;

    // Add the newCut at the head of the list
  } else if (first < firstFirst) { 
    newCut->next = list->head;
    list->head = newCut;

    // Add the newCut somewhere in the list not at the head or tail.
  } else { 
    mbSingleCut *c = list->head;
    while (first >= c->next->indexFirst)
      c = c->next;
    newCut->next = c->next;
    c->next = newCut;
  }

  // Combine overlapping cuts into a single cut, where possible.
  list->ncuts++;
  cutsCombine(list);
}



/// Allocate a mbCutList structure.
/// Also registers the deallocator.

mbCutList *mbCutListAlloc()
{
  mbCutList *list = psAlloc(sizeof(mbCutList));
  //psMemSetDeallocator(list, (psFreeFunc)CutListFree);
    
  list->head = NULL;
  list->ncuts = 0;

  return list;
}




/// If any consecutive cuts in a mbCutList are continguous,combine them into one.
/// This assumes that the indexFirst values in the list are pre-sorted, and it 
/// assures that they remain that way.
/// \param list  The mbCutList to be fixed.

static void cutsCombine(mbCutList *list)
{
  // No need to check empty lists.
  if (list == NULL || list->head == 0)
    return;

  mbSingleCut *ca, *cb;
  ca=list->head; 
  cb=ca->next;
  while (cb != NULL) {
    // If 1st sample after 1st cut (ca) is included in the next cut (cb), then extend ca to end at
    // the later of ca or cb.
    if (ca->indexLast+1 >= cb->indexFirst ||
        ca->indexLast == INT_MAX) {
      ca->indexLast = MAX(ca->indexLast, cb->indexLast);
      ca->next = cb->next;
      psFree(cb);

      cb = ca->next;
      list->ncuts --;
      continue;
    }
    ca = ca->next;
    cb = ca->next;
  }
}



/// Allocate a mbSingleCut structure.
/// Also registers the deallocator.
/// \param first  The first data index to be cut.
/// \param last   The last data index to be cut.

mbSingleCut *mbSingleCutAlloc(int first, int last)
{
  mbSingleCut *cut = psAlloc(sizeof(mbSingleCut));

  cut->next = NULL;
  cut->indexFirst = first;
  cut->indexLast = last;
    
  return cut;
}



/// Allocate a mbCuts structure.
/// Allocates the per-detector mbCutList objects but gives them empty cuts lists.
/// Also registers the deallocator.

mbCuts *mbCutsAlloc(int nrow,    ///< Number of rows
		    int ncol)    ///< Number of columns
{
  assert (nrow > 0);
  assert (ncol > 0);

  mbCuts *cuts = psAlloc(sizeof(mbCuts));
  //psMemSetDeallocator(cuts, (psFreeFunc)CutsFree);

  cuts->nrow = nrow;
  cuts->ncol = ncol;

  cuts->detCuts = psAlloc(nrow*sizeof(mbCutList **));
  cuts->detCuts[0] = psAlloc(nrow*ncol*sizeof(mbCutList *));
  for (int i=1; i<nrow; i++)
    cuts->detCuts[i] = cuts->detCuts[0] + i*ncol;
  for (int r=0; r<nrow; r++)
    for (int c=0; c<ncol; c++)
      cuts->detCuts[r][c] = NULL;

  cuts->globalCuts = mbCutListAlloc();
  cuts->globalCuts->head=NULL;
  //cuts->tod = NULL;
  
  omp_init_lock(&(cuts->cutlock));
  return cuts;  
}


/// Return whether no data are cut.
/// Recognize this state by having a cut list with no entries.
/// \param list  The cut list to be checked.

static bool isThisListNeverCut(const mbCutList *list)
{
  if (list == NULL || list->head == NULL)
    return true;
  return false;
}
	


/// Return whether a single detector is cut for the entire TOD.
/// \param list  The cut list to be checked.
/// \param row   The detector row number.
/// \param col   The detector column number.

bool mbCutsIsAlwaysCut(const mbCuts *cuts, int row, int col)
{
  // Say it is cut if arguments are invalid.
  if (cuts == NULL)
    return true;
  if (row < 0 || row >= cuts->nrow) return true;
  if (col < 0 || col >= cuts->ncol) return true;

  // For valid arguments, find status.
  if (isThisListAlwaysCut(cuts->globalCuts))
    return true;
  if (isThisListNeverCut(cuts->detCuts[row][col]))
    return false;
  return (isThisListAlwaysCut(cuts->detCuts[row][col]));
}



/// Return whether all data are cut.
/// Recognize this state by having a cut encompassing all nonnegative integers.
/// \param list  The cut list to be checked.

static bool isThisListAlwaysCut(const mbCutList *list)
{
  if (list == NULL || list->head == NULL)
    return false;
  if (list->head->indexFirst <= 0 &&
      list->head->indexLast == INT_MAX)
    return true;
  return false;
}
	


static void
mbUncutFree( mbUncut *uncut )
{
    psFree( uncut->indexFirst );
    psFree( uncut->indexLast );
}


/*--------------------------------------------------------------------------------*/

///Invert (i.e. get cuts from uncuts) an mbUncut object for a detector
mbUncut *mbCutsInvertUncut(mbUncut *uncut, int ndata)
{
  if (uncut==NULL)
    return NULL;  //don't know what to do with empty uncuts
  

  mbUncut *cut = psAlloc( sizeof(mbUncut) );
  cut->indexFirst=(int *)malloc(sizeof(int)*(uncut->nregions+2));
  cut->indexLast=(int *)malloc(sizeof(int)*(uncut->nregions+2));
  cut->nregions=0;
  //printf("uncut nregions is %d\n",uncut->nregions);
  //for (int i=0;i<uncut->nregions;i++)
    //printf("region %d is %5d %5d\n",i,uncut->indexFirst[i],uncut->indexLast[i]);
  if (uncut->nregions==0) {
    cut->nregions=1;
    cut->indexFirst[0]=0;
    cut->indexLast[0]=ndata;
    return cut;
  }

  if (uncut->indexFirst >0) {
    cut->indexFirst[cut->nregions]=0;
    cut->indexLast[cut->nregions]=uncut->indexFirst[0]-1;
    cut->nregions++;
  }
  for (int i=0;i<uncut->nregions-1;i++) {
    cut->indexFirst[cut->nregions]=uncut->indexLast[i];
    cut->indexLast[cut->nregions]=uncut->indexFirst[i+1];
    cut->nregions++;    
  }
  if (uncut->indexLast[uncut->nregions-1]<ndata) {
    cut->indexFirst[cut->nregions]=uncut->indexLast[uncut->nregions-1];
    cut->indexLast[cut->nregions]=ndata;
    cut->nregions++;
  }
  return cut;
  
}

///Get an mbUncut object for a detector
mbUncut *
mbCutsGetUncut( 
        const mbCuts *cuts, ///< mbCuts from which to draw the mbUncut
        int row,            ///< detector row
        int col,            ///< detector column
        int min,            ///< minimum data sample number to include in the mbUncut
        int max             ///< maximum data sample number to include in the mbUncut 
        )
{
    if ( row < 0 || row >= cuts->nrow ||
         col < 0 || col >= cuts->ncol )
        return NULL;

    mbUncut *uncut = psAlloc( sizeof(mbUncut) );
    //psMemSetDeallocator( uncut, (psFreeFunc) mbUncutFree );

    if ( cuts == NULL || mbCutsIsNeverCut(cuts, row, col) )
    {
        uncut->indexFirst = psAlloc( sizeof(int) );
        uncut->indexLast = psAlloc( sizeof(int) );
        uncut->nregions = 1;
        uncut->indexFirst[0] = min;
        uncut->indexLast[0] = max;
        return uncut;
    }
    else if ( mbCutsIsAlwaysCut(cuts, row, col) )
    {
        uncut->nregions = 0;
        uncut->indexFirst = NULL;
        uncut->indexLast = NULL;
        return uncut;
    }

    mbCutList *list = cutListOr(cuts->globalCuts, cuts->detCuts[row][col]);
    assert (list->head != NULL);

    uncut->indexFirst = psAlloc( (list->ncuts+1)*sizeof(int) );
    uncut->indexLast = psAlloc( (list->ncuts+1)*sizeof(int) );
    uncut->nregions = 1;
    uncut->indexFirst[0] = min;
    uncut->indexLast[0] = max;

    for ( mbSingleCut *c = list->head; c != NULL; c=c->next )
    {
        // ignore cuts outside of [min,max]
        if ( c->indexLast < min )
            continue;
        if ( c->indexFirst > max )
            break;

        int i = uncut->nregions - 1;

        if ( c->indexFirst <= min ) // cut intersects min
            uncut->indexFirst[i] = c->indexLast + 1;
        else if ( c->indexLast >= max ) // cut intersects max
            uncut->indexLast[i] = c->indexFirst - 1;
        else // cut contained in (min,max)
        {
            uncut->nregions++;
            uncut->indexLast[i] = c->indexFirst - 1;
            uncut->indexFirst[i+1] = c->indexLast + 1;
            uncut->indexLast[i+1] = max;
        }
    }

    psFree(list);
    return uncut;
}


/// Return whether a single detector is cut for none of the TOD.
/// \param list  The cut list to be checked.
/// \param row   The detector row number.
/// \param col   The detector column number.

bool mbCutsIsNeverCut(const mbCuts *cuts, int row, int col)
{
  // Say it is cut if arguments are invalid.
  if (cuts == NULL)
    return false;
  if (row < 0 || row >= cuts->nrow) return false;
  if (col < 0 || col >= cuts->ncol) return false;

  // For valid arguments, find status.
  if (isThisListNeverCut(cuts->globalCuts) &&
      isThisListNeverCut(cuts->detCuts[row][col]))
    return true;
  return false;
}


static mbCutList *cutListOr(const mbCutList *a, const mbCutList *b)
{
  mbCutList *or = mbCutListAlloc();
  mbSingleCut *cut;
  if (!isThisListNeverCut( a )){
    for (cut = a->head; cut != NULL; cut=cut->next)
      cutListExtend(or, cut->indexFirst, cut->indexLast);
  }
  if (!isThisListNeverCut( b )){
    for (cut = b->head; cut != NULL; cut=cut->next)
        cutListExtend(or, cut->indexFirst, cut->indexLast);
  }
  return or;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Methods to change cut status
////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Set a single detector to be always cut.
/// \param cuts    The cut list to be changed.
/// \param row     Row number for the cut detector.
/// \param col     Column number for the cut detector.
/// \return True on success, false if error.

bool mbCutsSetAlwaysCut(mbCuts *cuts, int row, int col)
{
  if (cuts == NULL ||
      row < 0 || row >= cuts->nrow ||
      col < 0 || col >= cuts->ncol)
    return false;
  omp_set_lock(&(cuts->cutlock));
  mbCutsExtend(cuts, 0, INT_MAX, row, col);
  omp_unset_lock(&(cuts->cutlock));
  return true;
  
}

/*--------------------------------------------------------------------------------*/

int mbCutsWrite(const mbCuts *cuts, const char *filename)
{
  FILE *f;
  int nline = 0;

  f = fopen( filename, "w" );
  if (f == NULL) {
    psTrace("moby", 2, "Failed to open cut file: %s", filename);
    return nline;
  }

  fprintf( f, "%d %d\n", cuts->nrow, cuts->ncol );

  for ( mbSingleCut *ct = cuts->globalCuts->head; ct != NULL; ct = ct->next ){
    fprintf( f, "(%d,%d) ", ct->indexFirst, ct->indexLast );
  }

  nline += 2;

  fputc('\n', f);

  for ( int c = 0; c < cuts->ncol; c++ ) {
    for ( int r = 0; r < cuts->nrow; r++) {
      if ( isThisListNeverCut( cuts->detCuts[r][c] ) ) continue;
      fprintf( f, "r%2.2dc%2.2d: ", r, c );
      for ( mbSingleCut *ct = cuts->detCuts[r][c]->head; ct != NULL; ct = ct->next ){
        fprintf( f, "(%d,%d) ", ct->indexFirst, ct->indexLast );
      }
      fputc('\n', f);
      nline++;
    }
  }

  fclose(f);
  return nline;
}






/// This routine takes a list of mbCuts and returns an mbCuts
/// that is a combination of all those cuts

mbCuts *mbCutsAllocFromCuts( const mbCuts **cutsList, int ncuts )
{
  if ( cutsList == NULL ) return NULL;
  if ( cutsList[0] == NULL ) return NULL;
  mbCuts *uberCuts;
  uberCuts = mbCutsAlloc( cutsList[0]->nrow, cutsList[0]->ncol );
  for ( int i = 0; i < ncuts; i++ )
    {
      mbSingleCut *sc;
      const mbCuts *cuts = cutsList[i];

      // First global cuts
      if (cuts->globalCuts != NULL)
	{
	  for ( sc = cuts->globalCuts->head; sc != NULL; sc = sc->next )
	    {
	      mbCutsExtendGlobal( uberCuts, sc->indexFirst, sc->indexLast );
	    }
	}

      // now cuts for individual detectors
      for ( int r = 0; r < cuts->nrow; r++ ) {
	for ( int c = 0; c < cuts->ncol; c++ ) {
	  if ( cuts->detCuts[r][c] != NULL ) {
	    for ( sc = cuts->detCuts[r][c]->head; sc != NULL;
		  sc = sc->next ) {
	      mbCutsExtend( uberCuts, sc->indexFirst, sc->indexLast, r, c );
	    }
	  }
	}
      }
    }

  return uberCuts;
}


