/*! \file memory.c
 *  \brief Memory allocation and deallocation
 *
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

unsigned INT total_alloc_mem   = 0; // Total allocated memory amount
unsigned INT total_alloc_count = 0; // Total number of allocations

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void * fasp_mem_calloc (LONG size, INT type)
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
 * Modified by Chunsheng Feng on 07/23/2013
 * Modified by Chunsheng Feng on 07/30/2013
 */
void * fasp_mem_calloc (LONG size, 
                        INT type)
{
    const LONG tsize = size*type;    
    void * mem = NULL;
    
#if DEBUG_MODE
    printf("### DEBUG: Trying to allocate %.3fKB RAM!\n", tsize/1024.0);
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
#endif
    }
    
    if ( mem == NULL ) {
        printf("### ERROR: Fail to allocate %.3fKB RAM!\n", tsize/1024.0);
//        exit(ERROR_ALLOC_MEM);
    }
    
#if CHMEM_MODE
    total_alloc_count++;
#endif
    
    return mem;
}

/**
 * \fn void * fasp_mem_realloc (void * oldmem, LONG type)
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
 */
void * fasp_mem_realloc (void * oldmem, 
                         LONG tsize)
{
#if DLMALLOC
    void * mem = dlrealloc(oldmem,tsize);
#elif NEDMALLOC
    void * mem = nedrealloc(oldmem, tsize);
#else
    void * mem = realloc(oldmem,tsize);
#endif
    
    if (mem==NULL) {
        printf("### ERROR: Fail to allocate %.3fKB RAM!\n", tsize/1024.0);    
        exit(ERROR_ALLOC_MEM);
    }

#if CHMEM_MODE    
    total_alloc_mem += tsize;
#endif

    return mem;
}

/**
 * \fn void fasp_mem_free (void* mem)
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
void fasp_mem_free (void * mem)
{
    if (mem) {
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
    printf("### DEBUG: Total number of alloc = %d, allocated memory %.3fMB.\n",
           total_alloc_count,total_alloc_mem/1e+6);
#endif
}    

/**
 * \fn SHORT fasp_mem_check (void * ptr, char *message, const INT ERR)
 *
 * \brief Check wether a point is null or not. 
 *
 * \param ptr       Void pointer to be checked
 * \param message   Error message to print
 * \param ERR       Integer error code
 *
 * \return          SUCCESS or error code
 *
 * \author Chensong Zhang
 * \date   11/16/2009
 */
SHORT fasp_mem_check (void *ptr, 
                      char *message, 
                      INT ERR)
{    
    if (ptr==NULL) {
        printf("### ERROR: %s", message);
        return ERR;
    }
    
    return SUCCESS;
} 

/**
 * \fn SHORT fasp_mem_iludata_check (ILU_data *iludata)
 *
 * \brief Check wether a ILU_data has enough work space
 *
 * \param iludata    Pointer to be cheked
 *
 * \return           SUCCESS if success, else ERROR (negative value) 
 *
 * \author Xiaozhe Hu, Chensong Zhang
 * \date   11/27/09
 */
SHORT fasp_mem_iludata_check (ILU_data *iludata)
{    
    const INT memneed = 2*iludata->row; // estimated memory usage
    
    if (iludata->nwork >= memneed) {
        return SUCCESS;
    }
    else {
        printf("### ERROR: ILU needs %d memory, only %d available!!!\n", 
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
 * \return    SUCCESS if success, else ERROR message (negative value) 
 *
 * \author Xiaozhe Hu
 * \date   11/27/09
 */
SHORT fasp_mem_dcsr_check (dCSRmat *A)
{
    if ((A->IA == NULL) || (A->JA == NULL) || (A->val == NULL)) {
        printf("### ERROR: dCSRmat has not sucessfully allocate memory!!!\n");
        return ERROR_ALLOC_MEM;
    }
    else {
        return SUCCESS;
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
