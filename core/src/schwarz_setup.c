/*!
 *  \file schwarz_setup.c
 *  \brief Setup phase for the Schwarz methods
 *
 *  Created by Hu Xiaozhe on 3/22/11.
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"
#include "forts_ns.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_schwarz_setup(Schwarz_data *schwarz, INT mmsize, INT maxlev, INT schwarz_type)
 *
 * \brief Setup phase for the Schwarz methods
 *
 * \param schwarz        Pointer to the shcwarz data
 * \param mmsize         Max block size
 * \param maxlev         Max number of levels
 * \param schwarz_type   Type of the Schwarz method
 *
 * \return               SUCCESS if succeed
 *
 * \author Ludmil, Xiaozhe Hu
 * \date   03/22/2011 
 */
INT fasp_schwarz_setup (Schwarz_data *schwarz, 
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
    printf("### DEBUG: fasp_schwarz_setup ...... [Start]\n");
#endif

	// allocate memory
	maxa    = (INT *)fasp_mem_calloc(n1,sizeof(INT));
	mask    = (INT *)fasp_mem_calloc(n1,sizeof(INT));
	iblock  = (INT *)fasp_mem_calloc(n1,sizeof(INT));
	jblock  = (INT *)fasp_mem_calloc(n1,sizeof(INT));
	
	nsizeall=0;
	
	for (i=0;i<n1;i++) {
        mask[i]=0;
        iblock[i]=0;
        maxa[i]=0;
	}
	
	maxa[0]=1;
	
	MIS = fasp_sparse_MIS(&A);
	
	/*-------------------------------------------*/
	//! find the blocks
	/*-------------------------------------------*/
	// first pass
	for (i=0;i<MIS.row;i++) {
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
	//! LU decomposition of blocks
	/*-------------------------------------------*/
	n1 = (nblk + iblock[nblk]-1);

	maxa = (INT *)fasp_mem_realloc(maxa,n1*sizeof(INT));
	
	for (i=0;i<n1;i++) maxa[i]=0;

	for (i=0;i<n; i++) mask[i]=0;
	
	// first estimate the memroy we need.
	mxfrm2_(&n,ia,ja,&nblk,iblock,jblock,mask,maxa,&memt,&maxbs);
	
#if DEBUG_MODE
	fprintf(stdout,"### DEBUG: Number of nonzeroes for LU=%d maxbs=%d\n",memt, maxbs);
#endif
	
	// allocate the memory
	al     = (REAL *)fasp_mem_calloc(memt,  sizeof(REAL));
	au     = (REAL *)fasp_mem_calloc(memt,  sizeof(REAL));
	rhsloc = (REAL *)fasp_mem_calloc(maxbs, sizeof(REAL));
    
	//  LU decomposition
    sky2ns_(&n,ia,ja,a,&nblk,iblock,jblock,mask,maxa,au,al);
	  
	printf("Schwarz setup succeeded: matrix size = %d, #blocks = %d, max block size =%d\n", 
           n, nblk, maxbs);
      
	/*-------------------------------------------*/  
	//! return
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
    printf("### DEBUG: fasp_schwarz_setup ...... [Finish]\n");
#endif
	
	return flag;	
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
