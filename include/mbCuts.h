#if !defined(MB_CUTS_H)
#define MB_CUTS_H 1

#include <ps_stuff.h>
#include <omp.h>

//#define HAVE_MBTOD

#ifdef HAVE_MBTOD
#include "mbTOD.h"
#endif

/// \file mbCuts.h
/// The public interface for the mbCuts object and its methods.
/// The object allows us to remove from analysis periods of data for
/// single detectors or the entire array.


/// Represents a single range of data indices to be cut.
/// Stored as a linked list of cuts.

typedef struct mbSingleCut {
  struct mbSingleCut *next;///< Pointer to next mbCut in the list.
  int indexFirst;        ///< First index to be cut.
  int indexLast;         ///< Last index to be cut.
} mbSingleCut;



/// Represents a full set of cuts for a detector (or array), with 0 or more ranges of cut data.
/// The cuts are stored as a linked list.  
/// This object contains the head and length of the list.
/// Having no cuts at all is represented by head==NULL.

typedef struct {
  mbSingleCut *head;     ///< Pointer to linked list of cuts (or NULL if no cuts).
  int ncuts;             ///< Number of cuts in the list.
} mbCutList;


typedef struct
{
    int nregions;
    int *indexFirst;
    int *indexLast;
}
mbUncut;

typedef struct
{
  int nregions;
  int *nparams;
  int *starting_param;
  actData ***precon;
}
mbCutFitParams;

/// Represents all cuts for a time-ordered-data set.
/// All cuts includes both global and detector-specific cuts.
/// To be useful, it should be associated with a TOD.

typedef struct {
  int nrow;              ///< Number of detector rows.
  int ncol;              ///< Number of detector columns.
  mbCutList ***detCuts;  ///< Ptr to 2D array of per-detector cutLists.
  mbCutList *globalCuts; ///< The global cutList.

  omp_lock_t cutlock;
} mbCuts;




// Allocators and other setup methods (deallocators are private)
// Note: perhaps we want these first two to be private?
mbSingleCut *mbSingleCutAlloc(int first, int last);
mbCutList *mbCutListAlloc();
mbCuts *mbCutsAlloc(int nrow, int ncol);
mbCuts *mbCutsAllocFromFile( const char *filename );
mbCuts *mbCutsAllocFromCuts( const mbCuts **cutsList, int ncuts );


void CutsFree(mbCuts *cuts);



// Add a range to one of the cuts lists.
void mbCutsExtendGlobal(mbCuts *cuts, int first, int last);
void mbCutsExtend(mbCuts *cuts, int first, int last, int row, int col);
int mbCutsExtendByArray( mbCuts *cuts, int row, int col, const int *array, int ndata );


// Set cut status
bool mbCutsSetAlwaysCut(mbCuts *cuts, int row, int col);
bool mbCutsSetAlwaysCutList(mbCuts *cuts, int ndets,  const int *rowlist, const int *collist);
bool mbCutsSetNeverCut(mbCuts *cuts, int row, int col);
bool mbCutsSetNeverCutList(mbCuts *cuts, int ndets,  const int *rowlist, const int *collist);
void mbCutsBuffer( mbCuts *cuts, int size, int ndata );


// Inquire about cut status
bool mbCutsIsAlwaysCut(const mbCuts *cuts, int row, int col);
bool mbCutsIsNeverCut(const mbCuts *cuts, int row, int col);
bool mbCutsIsCut(const mbCuts *cuts, int row, int col, int index);
int mbCutsListAlwaysCut(const mbCuts *cuts, int **rowlist, int **collist);
#ifdef HAVE_MBTOD
int mbCutsGoodSampleList(const mbCuts *cuts, char **goodlist, int row, int col, mbTOD *tod);
#endif
int mbCutsGetNCut(const mbCuts *cuts, int row, int col );
mbCutList *mbGetCutList( mbCuts *cuts, int row, int col );

mbUncut *
mbCutsGetUncut( const mbCuts *cuts, int row, int col, int min, int max );
mbUncut *mbCutsInvertUncut(mbUncut *uncut, int ndata);
int    mbUncutGetNumberOfRegions( mbUncut *uncut );
void   mbUncutGetRegionLimits( mbUncut *uncut, int i, int *min, int *max );

// I/O
int mbCutsWrite(const mbCuts *cuts, const char *filename );
int mbCutsRead( mbCuts *cuts, const char *filename );

// Specific Cuts
#ifdef HAVE_MBTOD
double mbCutsBuzzCut( mbTOD* tod, int ndata, float max_rms, float kill_frac );
double mbCutsRFCut( mbTOD* tod );
#endif

void mbCutsDecimate(mbCuts *cuts);


#endif // defined(MB_CUTS_H)


