/*! \file memory.c
 *
 *  \brief Memory allocation and deallocation subroutines
 */

#include "fasp.h"

#if DLMALLOC

#include "dlmalloc.h"

#elif NEDMALLOC

#include "nedmalloc.h"

#ifdef __cplusplus
extern "C" {
#endif
    void * nedcalloc(size_t no, size_t size);
    void * nedrealloc(void *mem, size_t size);
    void   nedfree(void *mem);
#ifdef __cplusplus
}
#endif

#endif

/*---------------------------------*/
/*--      Global Variables       --*/
/*---------------------------------*/

unsigned INT total_alloc_mem   = 0;        //! Total allocated memory amount
unsigned INT total_alloc_count = 0;        //! Total number of allocations
const    INT Million           = 1048576;  //! 1M = 1024*1024

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void * fasp_mem_calloc (LONGLONG size, INT type)
 *
 * \brief Allocate, initiate, and check memory
 *
 * \param size    Number of memory blocks
 * \param type    Size of memory blocks
 *
 * \return        Void pointer to the allocated memory
 *
 * \author Chensong Zhang
 * \date   2010/08/12
 *
 * Modified by Chunsheng Feng on 12/20/2013
 * Modified by Chunsheng Feng on 07/23/2013
 * Modified by Chunsheng Feng on 07/30/2013
 * Modified by Chensong Zhang on 07/30/2013: print error if failed
 */
void * fasp_mem_calloc (LONGLONG  size,
                        INT       type)
{
    const LONGLONG tsize = size*type;
    
    void * mem = NULL;
    
#if DEBUG_MODE > 2
    printf("### DEBUG: Trying to allocate %.3lfMB RAM!\n", (REAL)tsize/Million);
#endif
    
    if ( tsize > 0 ) {
        
#if DLMALLOC
        mem = dlcalloc(size,type);
#elif NEDMALLOC
        mem = nedcalloc(size,type);
#else
        mem = calloc(size,type);
#endif
        
#if CHMEM_MODE
        total_alloc_mem += tsize;
        total_alloc_count++;
#endif
    }
    
    if ( mem == NULL ) {
        printf("### WARNING: Cannot allocate %.3lfMB RAM!\n", (REAL)tsize/Million);
    }
    
    return mem;
}

/**
 * \fn void * fasp_mem_realloc (void * oldmem, LONGLONG type)
 *
 * \brief Reallocate, initiate, and check memory
 *
 * \param oldmem  Pointer to the existing mem block
 * \param type    Size of memory blocks
 *
 * \return        Void pointer to the reallocated memory
 *
 * \author Chensong Zhang
 * \date   2010/08/12
 *
 * Modified by Chunsheng Feng on 07/23/2013
 * Modified by Chensong Zhang on 07/30/2013: print error if failed
 */
void * fasp_mem_realloc (void     *oldmem,
                         LONGLONG  tsize)
{
    void * mem = NULL;

#if DEBUG_MODE > 2
    printf("### DEBUG: Trying to allocate %.3lfMB RAM!\n", (REAL)tsize/Million);
#endif
    
    if ( tsize > 0 ) {
        
#if DLMALLOC
        mem = dlrealloc(oldmem,tsize);
#elif NEDMALLOC
        mem = nedrealloc(oldmem,tsize);
#else
        mem = realloc(oldmem,tsize);
#endif
        
    }
    
    if ( mem == NULL ) {
        printf("### WARNING: Cannot allocate %.3lfMB RAM!\n", (REAL)tsize/Million);
    }
    
    return mem;
}

/**
 * \fn void fasp_mem_free (void *mem)
 *
 * \brief Free up previous allocated memory body
 *
 * \param mem   Pointer to the memory body need to be freed
 *
 * \return      NULL pointer
 *
 * \author Chensong Zhang
 * \date   2010/12/24
 */
void fasp_mem_free (void *mem)
{
    if ( mem ) {
#if DLMALLOC
        dlfree(mem);
#elif NEDMALLOC
        nedfree(mem);
#else
        free(mem);
#endif
        
#if CHMEM_MODE
        total_alloc_count--;
#endif
    }
}

/**
 * \fn void fasp_mem_usage ()
 *
 * \brief Show total allocated memory currently
 *
 * \author Chensong Zhang
 * \date   2010/08/12
 */
void fasp_mem_usage ()
{
#if CHMEM_MODE
    printf("### DEBUG: Number of alloc = %d, allocating memory = %.3fMB.\n",
           total_alloc_count, (REAL)total_alloc_mem/Million);
#endif
}

/**
 * \fn SHORT fasp_mem_check (void * ptr, const char *message, const INT ERR)
 *
 * \brief Check wether a point is null or not.
 *
 * \param ptr       Void pointer to be checked
 * \param message   Error message to print
 * \param ERR       Integer error code
 *
 * \return          FASP_SUCCESS or error code
 *
 * \author Chensong Zhang
 * \date   11/16/2009
 */
SHORT fasp_mem_check (void        *ptr,
                      const char  *message,
                      INT          ERR)
{
    if ( ptr == NULL ) {
        printf("### ERROR: %s", message);
        return ERR;
    }
    else {
        return FASP_SUCCESS;
    }
}

/**
 * \fn SHORT fasp_mem_iludata_check (ILU_data *iludata)
 *
 * \brief Check wether a ILU_data has enough work space
 *
 * \param iludata    Pointer to be cheked
 *
 * \return           FASP_SUCCESS if success, else ERROR (negative value)
 *
 * \author Xiaozhe Hu, Chensong Zhang
 * \date   11/27/09
 */
SHORT fasp_mem_iludata_check (ILU_data *iludata)
{
    const INT memneed = 2*iludata->row; // estimated memory usage
    
    if ( iludata->nwork >= memneed ) {
        return FASP_SUCCESS;
    }
    else {
        printf("### ERROR: ILU needs %d RAM, only %d allocated!\n",
               memneed, iludata->nwork);
        return ERROR_ALLOC_MEM;
    }
}

/**
 * \fn SHORT fasp_mem_dcsr_check (dCSRmat *A)
 *
 * \brief Check wether a dCSRmat A has sucessfully allocated memory
 *
 * \param A   Pointer to be cheked
 *
 * \return    FASP_SUCCESS if success, else ERROR message (negative value)
 *
 * \author Xiaozhe Hu
 * \date   11/27/09
 */
SHORT fasp_mem_dcsr_check (dCSRmat *A)
{
    if ( (A->IA == NULL) || (A->JA == NULL) || (A->val == NULL) ) {
        printf("### ERROR: Something is wrong with the matrix!\n");
        return ERROR_ALLOC_MEM;
    }
    else {
        return FASP_SUCCESS;
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
