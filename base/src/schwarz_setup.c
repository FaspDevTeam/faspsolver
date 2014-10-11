/*! \file schwarz_setup.c
 *
 *  \brief Setup phase for the Schwarz methods
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"
#include "forts_ns.h"

#include "mg_util.inl"

static void schwarz_levels (INT, dCSRmat *, INT *, INT *, INT *, INT *, INT);

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn fasp_schwarz_get_block_matrix (Schwarz_data *schwarz, INT nblk, INT *iblock,
 *                                    INT *jblock, INT *mask)
 *
 * \brief Form schwarz partition data
 *
 * \param schwarz Pointer to the shcwarz data
 * \param nblk    Number of partitions
 * \param iblock  Pointer to number of vertices on each level
 * \param jblock  Pointer to vertices of each level
 * \param mask    Pointer to flag array
 *
 * \author Zheng Li, Chensong Zhang
 * \date   2014/09/29
 */
void fasp_schwarz_get_block_matrix (Schwarz_data *schwarz,
                                    INT nblk,
                                    INT *iblock,
                                    INT *jblock,
                                    INT *mask)
{
    INT i, j, iblk, ki, kj, kij, kil, is, ibl0, ibl1, nloc, iaa, iab;
    INT maxbs = 0, count, nnz;
    
    dCSRmat A = schwarz->A;
    dCSRmat *blk = schwarz->blk_data;
    
    INT  *ia  = A.IA;
    INT  *ja  = A.JA;
    REAL *val = A.val;
    
    // get maximal block size
    for (is=0; is<nblk; ++is) {
        ibl0 = iblock[is];
        ibl1 = iblock[is+1];
        nloc = ibl1-ibl0;
        maxbs = MAX(maxbs, nloc);
    }
    
    // allocate memory for each sub_block's right hand
    schwarz->xloc1   = fasp_dvec_create(maxbs);
    schwarz->rhsloc1 = fasp_dvec_create(maxbs);
    
    for (is=0; is<nblk; ++is) {
        ibl0 = iblock[is];
        ibl1 = iblock[is+1];
        nloc = ibl1-ibl0;
        count = 0;
        for (i=0; i<nloc; ++i ) {
            iblk = ibl0 + i;
            ki   = jblock[iblk];
            iaa  = ia[ki]-1;
            iab  = ia[ki+1]-1;
            count += iab - iaa;
            mask[ki] = i+1;
        }
        
        blk[is] = fasp_dcsr_create(nloc, nloc, count);
        blk[is].IA[0] = 0;
        nnz = 0;
        
        for (i=0; i<nloc; ++i) {
            iblk = ibl0 + i;
            ki = jblock[iblk];
            iaa = ia[ki]-1;
            iab = ia[ki+1]-1;
            for (kij = iaa; kij<iab; ++kij) {
                kj = ja[kij]-1;
                j  = mask[kj];
                if(j != 0) {
                    blk[is].JA[nnz] = j-1;
                    blk[is].val[nnz] = val[kij];
                    nnz ++;
                }
            }
            blk[is].IA[i+1] = nnz;
        }
        
        blk[is].nnz = nnz;
        
        // zero the mask so that everyting is as it was
        for (i=0; i<nloc; ++i) {
            iblk = ibl0 + i;
            ki   = jblock[iblk];
            mask[ki] = 0;
        }
    }
}

/**
 * \fn INT fasp_schwarz_setup (Schwarz_data *schwarz, INT mmsize,
 *                             INT maxlev, INT schwarz_type)
 *
 * \brief Setup phase for the Schwarz methods
 *
 * \param schwarz        Pointer to the shcwarz data
 * \param mmsize         Max block size
 * \param maxlev         Max number of levels
 * \param schwarz_type   Type of the Schwarz method
 *
 * \return               FASP_SUCCESS if succeed
 *
 * \author Ludmil, Xiaozhe Hu
 * \date   03/22/2011
 *
 * Modified by Zheng Li on 10/09/2014
 */
INT fasp_schwarz_setup (Schwarz_data *schwarz,
                        Schwarz_param *param)
{
    // information about A
    dCSRmat A = schwarz->A;
    INT n   = A.row;
    INT *ia = A.IA;
    INT *ja = A.JA;
    REAL *a = A.val;
    
    INT  block_solver = param->schwarz_blksolver;
    INT  maxlev = param->schwarz_maxlvl;
    
    // local variables
    INT n1 = n,i;
    INT inroot = -10, nsizei = -10, nsizeall = -10, nlvl = 0;
    INT maxbs=0;
    INT *jb=NULL;
    ivector MIS;
    
    // data for schwarz method
    INT nblk;
    INT *iblock = NULL, *jblock = NULL, *mask = NULL, *maxa = NULL;
    REAL *au = NULL, *al = NULL, *rhsloc = NULL;
    INT memt = 0;
    
    // return
    INT flag = 0;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif
    
    // allocate memory
    maxa    = (INT *)fasp_mem_calloc(n1,sizeof(INT));
    mask    = (INT *)fasp_mem_calloc(n1,sizeof(INT));
    iblock  = (INT *)fasp_mem_calloc(n1,sizeof(INT));
    jblock  = (INT *)fasp_mem_calloc(n1,sizeof(INT));
    
    nsizeall=0;
    memset(mask,   0, sizeof(INT)*n1);
    memset(iblock, 0, sizeof(INT)*n1);
    memset(maxa,   0, sizeof(INT)*n1);
    
    maxa[0]=0;
    
    MIS = fasp_sparse_MIS(&A);
    
    /*-------------------------------------------*/
    // find the blocks
    /*-------------------------------------------*/
    
    // first pass: do a maxlev level sets out for each node
    for ( i = 0; i < MIS.row; i++ ) {
        inroot = MIS.val[i];
        schwarz_levels(inroot,&A,mask,&nlvl,maxa,jblock,maxlev);
        nsizei=maxa[nlvl];
        nsizeall+=nsizei;
    }
    
#if DEBUG_MODE
    fprintf(stdout,"### DEBUG: nsizeall is: %d\n",nsizeall);
#endif
    
    /* We only calculated the size of this up to here. So we can reallocate jblock */
    jblock = (INT *)fasp_mem_realloc(jblock,(nsizeall+n)*sizeof(INT));
    
    // second pass: redo the same again, but this time we store in jblock
    maxa[0]=0;
    iblock[0]=0;
    nsizeall=0;
    jb=jblock;
    for (i=0;i<MIS.row;i++) {
        inroot = MIS.val[i];
        schwarz_levels(inroot,&A,mask,&nlvl,maxa,jb,maxlev);
        nsizei=maxa[nlvl];
        iblock[i+1]=iblock[i]+nsizei;
        nsizeall+=nsizei;
        jb+=nsizei;
    }
    nblk = MIS.row;
    
#if DEBUG_MODE
    fprintf(stdout,"### DEBUG: nsizall is: %d %d\n",nsizeall,iblock[nblk]);
    fprintf(stdout,"### DEBUG: nnz is: %d\n",ia[n]);
#endif
    
    /*-------------------------------------------*/
    //  LU decomposition of blocks
    /*-------------------------------------------*/
    n1 = (nblk + iblock[nblk]);
    
    memset(mask, 0, sizeof(INT)*n);
    
    schwarz->blk_data = (dCSRmat*)fasp_mem_calloc(nblk, sizeof(dCSRmat));
    
    fasp_schwarz_get_block_matrix(schwarz, nblk, iblock, jblock, mask);
    
#if DEBUG_MODE
    fprintf(stdout,"### DEBUG: Number of nonzeroes for LU=%d maxbs=%d\n",
            memt, maxbs);
#endif
    
    dCSRmat *blk = schwarz->blk_data;
    
    // Setup for each block solver
    switch (block_solver) {
            
#if WITH_MUMPS
        case SOLVER_MUMPS: {
            /* use MUMPS direct solver on each block */
            Mumps_data *mumps = (Mumps_data*)fasp_mem_calloc(nblk, sizeof(Mumps_data));
            for (i=0; i<nblk; ++i)
                mumps[i] = fasp_mumps_factorize(&blk[i], NULL, NULL, PRINT_NONE);
            schwarz->mumps = mumps;
            
            break;
        }
#endif
            
#if WITH_UMFPACK
        case SOLVER_UMFPACK: {
            /* use UMFPACK direct solver on each block */
            void **numeric	= (void**)fasp_mem_calloc(nblk, sizeof(void*));
            dCSRmat Ac_tran;
            for (i=0; i<nblk; ++i) {
                fasp_dcsr_trans(&blk[i], &Ac_tran);
                fasp_dcsr_sort(&Ac_tran);
                fasp_dcsr_cp(&Ac_tran, &blk[i]);
                numeric[i] = fasp_umfpack_factorize(&blk[i], 0);
            }
            schwarz->numeric = numeric;
            fasp_dcsr_free(&Ac_tran);
            
            break;
        }
#endif
            
        default: {
            /* need to do nothing for iterative methods */
        }
    }
    
#if DEBUG_MODE
    fprintf(stdout,"### DEBUG: n = %d, #blocks = %d, max block size = %d\n",
            n, nblk, maxbs);
#endif
    
    /*-------------------------------------------*/
    //  return
    /*-------------------------------------------*/
    schwarz->nblk   = nblk;
    schwarz->iblock = iblock;
    schwarz->jblock = jblock;
    schwarz->mask   = mask;
    schwarz->schwarz_type = param->schwarz_type;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return flag;
}

/**
 * \fn void fasp_dcsr_schwarz_smoother (Schwarz_data *schwarz, dvector *x, dvector *b)
 *
 * \brief Schwarz smoother
 *
 * \param schwarz Pointer to the shcwarz data
 * \param x       Pointer to solution vector
 * \param b       Pointer to right hand
 *
 * \author Zheng Li, Chensong Zhang
 * \date   2014/10/5
 */
void fasp_dcsr_schwarz_forward_smoother (Schwarz_data  *schwarz,
                                         Schwarz_param *param,
                                         dvector       *x,
                                         dvector       *b)
{
    INT i, j, iblk, ki, kj, kij, is, ibl0, ibl1, nloc, iaa, iab;
    
    // Schwarz partition
    INT  nblk = schwarz->nblk;
    dCSRmat *blk = schwarz->blk_data;
    INT  *iblock = schwarz->iblock;
    INT  *jblock = schwarz->jblock;
    INT  *mask   = schwarz->mask;
    INT  block_solver = param->schwarz_blksolver;
    
    // Schwarz data
    dCSRmat A = schwarz->A;
    INT *ia = A.IA;
    INT *ja = A.JA;
    REAL *val = A.val;
    
    // Local solution and right hand vectors
    dvector rhs = schwarz->rhsloc1;
    dvector u   = schwarz->xloc1;
    
#if WITH_UMFPACK
    void **numeric = schwarz->numeric;
#endif
    
#if WITH_MUMPS
    Mumps_data *mumps = schwarz->mumps;
#endif
    
    for (is=0; is<nblk; ++is) {
        // Form the right hand of eack block
        ibl0 = iblock[is];
        ibl1 = iblock[is+1];
        nloc = ibl1-ibl0;
        for (i=0; i<nloc; ++i ) {
            iblk = ibl0 + i;
            ki   = jblock[iblk];
            mask[ki] = i+1;
        }
        
        for (i=0; i<nloc; ++i) {
            iblk = ibl0 + i;
            ki = jblock[iblk];
            rhs.val[i] = b->val[ki];
            iaa = ia[ki]-1;
            iab = ia[ki+1]-1;
            for (kij = iaa; kij<iab; ++kij) {
                kj = ja[kij]-1;
                j  = mask[kj];
                if(j == 0) {
                    rhs.val[i] -= val[kij]*x->val[kj];
                }
            }
        }
        
        // Solve each block
        switch (block_solver) {
                
#if WITH_MUMPS
            case SOLVER_MUMPS: {
                /* use MUMPS direct solver on each block */
                fasp_mumps_solve(&blk[is], &rhs, &u, mumps[is], 0);
                break;
            }
#endif
                
#if WITH_UMFPACK
            case SOLVER_UMFPACK: {
                /* use UMFPACK direct solver on each block */
                fasp_umfpack_solve(&blk[is], &rhs, &u, numeric[is], 0);
                break;
            }
#endif
                /* use iterative solver on each block */
            default:
                u.row = blk[is].row;
                rhs.row = blk[is].row;
                fasp_dvec_set(u.row, &u, 0);
                fasp_solver_dcsr_pvgmres(&blk[is], &rhs, &u, NULL, 1e-8, 100, 20, 1, 0);
        }
        
        //zero the mask so that everyting is as it was
        for (i=0; i<nloc; ++i) {
            iblk = ibl0 + i;
            ki   = jblock[iblk];
            mask[ki] = 0;
            x->val[ki] = u.val[i];
        }
    }
}

/**
 * \fn void fasp_dcsr_schwarz_smoother (Schwarz_data *schwarz, dvector *x, dvector *b)
 *
 * \brief Schwarz smoother
 *
 * \param schwarz Pointer to the shcwarz data
 * \param x       Pointer to solution vector
 * \param b       Pointer to right hand
 *
 * \author Zheng Li, Chensong Zhang
 * \date   2014/10/5
 */
void fasp_dcsr_schwarz_backward_smoother (Schwarz_data *schwarz,
                                          Schwarz_param *param,
                                          dvector *x,
                                          dvector *b)
{
    INT i, j, iblk, ki, kj, kij, is, ibl0, ibl1, nloc, iaa, iab;
    
    // Schwarz partition
    INT  nblk = schwarz->nblk;
    dCSRmat *blk = schwarz->blk_data;
    INT  *iblock = schwarz->iblock;
    INT  *jblock = schwarz->jblock;
    INT  *mask   = schwarz->mask;
    INT  block_solver = param->schwarz_blksolver;
    
    // Schwarz data
    dCSRmat A = schwarz->A;
    INT *ia = A.IA;
    INT *ja = A.JA;
    REAL *val = A.val;
    
    // Local solution and right hand vectors
    dvector rhs = schwarz->rhsloc1;
    dvector u   = schwarz->xloc1;
    
#if WITH_UMFPACK
    void **numeric = schwarz->numeric;
#endif
    
#if WITH_MUMPS
    Mumps_data *mumps = schwarz->mumps;
#endif
    
    for (is=nblk-1; is>=0; --is) {
        // Form the right hand of eack block
        ibl0 = iblock[is];
        ibl1 = iblock[is+1];
        nloc = ibl1-ibl0;
        for (i=0; i<nloc; ++i ) {
            iblk = ibl0 + i;
            ki   = jblock[iblk];
            mask[ki] = i+1;
        }
        
        for (i=0; i<nloc; ++i) {
            iblk = ibl0 + i;
            ki = jblock[iblk];
            rhs.val[i] = b->val[ki];
            //iaa = ia[ki];
            //iab = ia[ki+1];
            iaa = ia[ki]-1;
            iab = ia[ki+1]-1;
            for (kij = iaa; kij<iab; ++kij) {
                kj = ja[kij]-1;
                //kj = ja[kij];
                j  = mask[kj];
                if(j == 0) {
                    rhs.val[i] -= val[kij]*x->val[kj];
                }
            }
        }
        
        // Solve each block
        switch (block_solver) {
                
#if WITH_MUMPS
            case SOLVER_MUMPS: {
                /* use MUMPS direct solver on each block */
                fasp_mumps_solve(&blk[is], &rhs, &u, mumps[is], 0);
                break;
            }
#endif
                
#if WITH_UMFPACK
            case SOLVER_UMFPACK: {
                /* use UMFPACK direct solver on each block */
                fasp_umfpack_solve(&blk[is], &rhs, &u, numeric[is], 0);
                break;
            }
#endif
            default:
                /* use iterative solver on each block */
                rhs.row = blk[is].row;
                u.row   = blk[is].row;
                fasp_dvec_set(u.row, &u, 0);
                fasp_solver_dcsr_pvgmres (&blk[is], &rhs, &u, NULL, 1e-8, 100, 20, 1, 0);
        }
        
        //zero the mask so that everyting is as it was
        for (i=0; i<nloc; ++i) {
            iblk = ibl0 + i;
            ki   = jblock[iblk];
            mask[ki] = 0;
            x->val[ki] = u.val[i];
        }
    }
}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static void schwarz_levels (INT inroot, dCSRmat *A, INT *mask, INT *nlvl,
 *                                 INT *iblock, INT *jblock, INT maxlev)
 *
 * \brief Form the level hierarchy of input root node
 *
 * \param inroot  Root node
 * \param A       Pointer to CSR matrix
 * \param mask    Pointer to flag array
 * \param nlvl    The number of levels to expand from root node
 * \param iblock  Pointer to vertices number of each level
 * \param jblock  Pointer to vertices of each level
 * \param maxlev  The maximal number of levels to expand from root node
 *
 * \author Zheng Li
 * \date   2014/09/29
 */
static void schwarz_levels (INT inroot,
                            dCSRmat *A,
                            INT *mask,
                            INT *nlvl,
                            INT *iblock,
                            INT *jblock,
                            INT maxlev)
{
    INT *ia = A->IA;
    INT *ja = A->JA;
    INT nnz = A->nnz;
    INT i, j, lvl, lbegin, lvlend, nsize, node;
    INT jstrt, jstop, nbr, lvsize;
    
    // This is diagonal
    if (ia[inroot+1]-ia[inroot] <= 1) {
        lvl = 0;
        iblock[lvl] = 0;
        jblock[iblock[lvl]] = inroot;
        lvl ++;
        iblock[lvl] = 1;
    }
    else {
        // input node as root node (level 0)
        lvl = 0;
        jblock[0] = inroot;
        lvlend = 0;
        nsize  = 1;
        // mark root node
        mask[inroot] = 1;
        
        lvsize = nnz;
        
        // start to form the level hierarchy for root node(root, level0, level 1,...)
        while (lvsize > 0 && lvl < maxlev) {
            lbegin = lvlend;
            lvlend = nsize;
            iblock[lvl] = lbegin;
            lvl ++;
            for(i=lbegin; i<lvlend; ++i) {
                node = jblock[i];
                jstrt = ia[node]-1;
                jstop = ia[node+1]-1;
                for (j = jstrt; j<jstop; ++j) {
                    nbr = ja[j]-1;
                    if (mask[nbr] == 0) {
                        jblock[nsize] = nbr;
                        mask[nbr] = lvl;
                        nsize ++;
                    }
                }
            }
            lvsize = nsize - lvlend;
        }
        
        iblock[lvl] = nsize;
        
        // reset mask array
        for (i = 0; i< nsize; ++i) {
            node = jblock[i];
            mask[node] = 0;
        }
    }
    
    *nlvl = lvl;
}

#if 0 // TODO: Need to remove this! --Chensong
/**
 * \fn INT fasp_schwarz_setup (Schwarz_data *schwarz, INT mmsize,
 *                             INT maxlev, INT schwarz_type)
 *
 * \brief Setup phase for the Schwarz methods
 *
 * \param schwarz        Pointer to the shcwarz data
 * \param mmsize         Max block size
 * \param maxlev         Max number of levels
 * \param schwarz_type   Type of the Schwarz method
 *
 * \return               FASP_SUCCESS if succeed
 *
 * \author Ludmil, Xiaozhe Hu
 * \date   03/22/2011
 *
 * Modified by Chunsheng Feng, Zheng Li on 08/28/2012
 */
static INT fasp_schwarz_setup (Schwarz_data *schwarz,
                               INT mmsize,
                               INT maxlev,
                               INT schwarz_type)
{
    // information about A
    dCSRmat A = schwarz->A;
    INT n   = A.row;
    INT *ia = A.IA;
    INT *ja = A.JA;
    REAL *a = A.val;
    
    // local variables
    INT n1=n+1,i;
    INT inroot=-10,nsizei=-10,nsizeall=-10,nlvl=0;
    INT maxbs=0;
    INT *jb=NULL;
    ivector MIS;
    
    // data for schwarz method
    INT nblk;
    INT *iblock=NULL, *jblock=NULL, *mask=NULL, *maxa=NULL;
    REAL *au=NULL, *al=NULL, *rhsloc=NULL;
    INT memt=0;
    
    // return
    INT flag = 0;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif
    
    // allocate memory
    maxa    = (INT *)fasp_mem_calloc(n1,sizeof(INT));
    mask    = (INT *)fasp_mem_calloc(n1,sizeof(INT));
    iblock  = (INT *)fasp_mem_calloc(n1,sizeof(INT));
    jblock  = (INT *)fasp_mem_calloc(n1,sizeof(INT));
    
    nsizeall=0;
    memset(mask,   0, sizeof(INT)*n1);
    memset(iblock, 0, sizeof(INT)*n1);
    memset(maxa,   0, sizeof(INT)*n1);
    
    maxa[0]=1;
    
    MIS = fasp_sparse_MIS(&A);
    
    /*-------------------------------------------*/
    // find the blocks
    /*-------------------------------------------*/
    // first pass
    for ( i = 0; i < MIS.row; i++ ) {
        // for each node do a maxlev level sets out
        inroot = MIS.val[i]+1;
        levels_(&inroot,ia,ja,mask,&nlvl,maxa,jblock,&maxlev);
        nsizei=maxa[nlvl]-1;
        nsizeall+=nsizei;
    }
    
#if DEBUG_MODE
    fprintf(stdout,"### DEBUG: nsizeall is: %d\n",nsizeall);
#endif
    
    /* We only calculated the size of this up to here. So we can reallocate jblock */
    jblock  = (INT *)fasp_mem_realloc(jblock,(nsizeall+n)*sizeof(INT));
    
    // second pass
    /* we redo the same again, but this time we store in jblock */
    maxa[0]=1;
    iblock[0]=1;
    nsizeall=0;
    jb=jblock;
    for (i=0;i<MIS.row;i++) {
        inroot = MIS.val[i]+1;
        levels_(&inroot,ia,ja,mask,&nlvl,maxa,jb,&maxlev);
        nsizei=maxa[nlvl]-1;
        iblock[i+1]=iblock[i]+nsizei;
        nsizeall+=nsizei;
        jb+=nsizei;
    }
    
    nblk = MIS.row;
    
#if DEBUG_MODE
    fprintf(stdout,"### DEBUG: nsizall is: %d %d\n",nsizeall,iblock[nblk]-1);
    fprintf(stdout,"### DEBUG: nnz is: %d\n",ia[n]-1);
#endif
    
    /*-------------------------------------------*/
    //  LU decomposition of blocks
    /*-------------------------------------------*/
    n1 = (nblk + iblock[nblk]-1);
    
    maxa = (INT *)fasp_mem_realloc(maxa,n1*sizeof(INT));
    
    memset(maxa, 0, sizeof(INT)*n1);
    memset(mask, 0, sizeof(INT)*n);
    
    // first estimate the memroy we need.
    mxfrm2_(&n,ia,ja,&nblk,iblock,jblock,mask,maxa,&memt,&maxbs);
    
#if DEBUG_MODE
    fprintf(stdout,"### DEBUG: Number of nonzeroes for LU=%d maxbs=%d\n",
            memt, maxbs);
#endif
    
    // allocate the memory
    al     = (REAL *)fasp_mem_calloc(memt,  sizeof(REAL));
    au     = (REAL *)fasp_mem_calloc(memt,  sizeof(REAL));
    rhsloc = (REAL *)fasp_mem_calloc(maxbs, sizeof(REAL));
    
    //  LU decomposition
    sky2ns_(&n,ia,ja,a,&nblk,iblock,jblock,mask,maxa,au,al);
    
#if DEBUG_MODE
    fprintf(stdout,"### DEBUG: n = %d, #blocks = %d, max block size = %d\n",
            n, nblk, maxbs);
#endif
    
    /*-------------------------------------------*/
    //  return
    /*-------------------------------------------*/
    schwarz->nblk = nblk;
    schwarz->iblock = iblock;
    schwarz->jblock = jblock;
    schwarz->rhsloc = rhsloc;
    schwarz->au = au;
    schwarz->al = al;
    
    schwarz->schwarz_type = schwarz_type;
    schwarz->memt = memt;
    schwarz->mask = mask;
    schwarz->maxa = maxa;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return flag;
}

#endif

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
