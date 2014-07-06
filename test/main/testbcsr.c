/*! \file testbcsr.c
 *  \brief The main test function for FASP solvers -- BCSR format
 */

#include "fasp.h"
#include "fasp_functs.h"

/**
 * \fn int main (int argc, const char * argv[])
 *
 * \brief This is the main function for a few simple tests.
 *
 * \author Xiaozhe Hu
 * \date   04/07/2014
 * 
 */
int main (int argc, const char * argv[]) 
{
    block_dCSRmat Abcsr;
	dvector b, uh;
    
    dCSRmat A;
    dCSRmat *A_diag;
    
    dBSRmat Absr;
    
    INT i,j;
    block_dCSRmat Aibcsr;
    INT NumLayers;
    dCSRmat *local_A;
    ivector *local_index;
    
	int status=FASP_SUCCESS;
	
	// Step 0. Set parameters
	input_param     inpar;  // parameters from input files
	itsolver_param  itpar;  // parameters for itsolver
	AMG_param       amgpar; // parameters for AMG
	ILU_param       ilupar; // parameters for ILU
    
    // Set solver parameters: use ./ini/bsr.dat
    fasp_param_set(argc, argv, &inpar);
    fasp_param_init(&inpar, &itpar, &amgpar, &ilupar, NULL);
    
    // Set local parameters
	const int print_level   = inpar.print_level;
	const int problem_num   = inpar.problem_num;
	const int itsolver_type = inpar.solver_type;
	const int precond_type  = inpar.precond_type;
	const int output_type   = inpar.output_type;
    
    // Set output devices
    if (output_type) {
		char *outputfile = "out/test.out";
		printf("Redirecting outputs to file: %s ...\n", outputfile);
		freopen(outputfile,"w",stdout); // open a file for stdout
	}
    
    printf("Test Problem %d\n", problem_num);
    
	// Step 1. Input stiffness matrix and right-hand side
	char filename1[512], *datafile1;
	char filename2[512], *datafile2;
    char filename3[512], *datafile3;
	char filename4[512], *datafile4;
    char filename5[512], *datafile5;
	char filename6[512], *datafile6;
    char filename7[512], *datafile7;
	char filename8[512], *datafile8;
    char filename9[512], *datafile9;
	char filename10[512], *datafile10;
	
	strncpy(filename1,inpar.workdir,128);
	strncpy(filename2,inpar.workdir,128);
    strncpy(filename3,inpar.workdir,128);
	strncpy(filename4,inpar.workdir,128);
    strncpy(filename5,inpar.workdir,128);
	strncpy(filename6,inpar.workdir,128);
    strncpy(filename7,inpar.workdir,128);
	strncpy(filename8,inpar.workdir,128);
    strncpy(filename9,inpar.workdir,128);
	strncpy(filename10,inpar.workdir,128);
    
    
    // Default test problem from black-oil benchmark: SPE01
	if (problem_num == 10) {
        
        strncpy(filename1,inpar.workdir,128);
        datafile1="/pnp-data/A.dat";
        strcat(filename1,datafile1);

        fasp_dcoo_shift_read(filename1, &A);
        
        // form block CSR matrix
        INT row = A.row/3;
		ivector phi_idx;
        ivector n_idx;
        ivector p_idx;
        
        fasp_ivec_alloc(row, &phi_idx);
        fasp_ivec_alloc(row, &n_idx);
        fasp_ivec_alloc(row, &p_idx);
        
        printf("row = %d\n",row);
        for (i=0; i<row; i++){
            phi_idx.val[i] = 3*i;
            n_idx.val[i] = 3*i+1;
            p_idx.val[i] = 3*i+2;
        }

        // Assemble the matrix in block dCSR format
		Abcsr.brow = 3; Abcsr.bcol = 3;
		Abcsr.blocks = (dCSRmat **)calloc(9, sizeof(dCSRmat *));
        for (i=0; i<9 ;i++) {
            Abcsr.blocks[i] = (dCSRmat *)fasp_mem_calloc(1, sizeof(dCSRmat));
        }
		      
		// A11
        fasp_dcsr_getblk(&A, phi_idx.val, phi_idx.val, row, row, Abcsr.blocks[0]);
        // A12
        fasp_dcsr_getblk(&A, phi_idx.val, n_idx.val, row, row, Abcsr.blocks[1]);
        // A13
        fasp_dcsr_getblk(&A, phi_idx.val, p_idx.val, row, row, Abcsr.blocks[2]);
        // A21
        fasp_dcsr_getblk(&A, n_idx.val, phi_idx.val, row, row, Abcsr.blocks[3]);
        // A22
        fasp_dcsr_getblk(&A, n_idx.val, n_idx.val, row, row, Abcsr.blocks[4]);
        // A23
        fasp_dcsr_getblk(&A, n_idx.val, p_idx.val, row, row, Abcsr.blocks[5]);
        // A31
        fasp_dcsr_getblk(&A, p_idx.val, phi_idx.val, row, row, Abcsr.blocks[6]);
        // A32
        fasp_dcsr_getblk(&A, p_idx.val, n_idx.val, row, row, Abcsr.blocks[7]);
        // A33
        fasp_dcsr_getblk(&A, p_idx.val, p_idx.val, row, row, Abcsr.blocks[8]);
        
        // form right hand side
        dvector b_temp;
        strncpy(filename2,inpar.workdir,128);
        datafile2="/pnp-data/rhs.dat"; strcat(filename2,datafile2);
        fasp_dvec_read(filename2, &b_temp);
        
        
        fasp_dvec_alloc(b_temp.row, &b);
        for (i=0; i<row; i++){
            
            b.val[i]        = b_temp.val[3*i];
            b.val[row+i]    = b_temp.val[3*i+1];
            b.val[2*row+i]  = b_temp.val[3*i+2];
            
        }
        
        //Absr = fasp_format_dcsr_dbsr(&A, 3);
        
    }
    
    // test problem for sweeping method
    else if (problem_num == 20) {
        
        printf("-------------------------\n");
        printf("Read matrix A\n");
        printf("-------------------------\n");
        
        strncpy(filename1,inpar.workdir,128);
        datafile1="/5layers/A.dat";
        strcat(filename1,datafile1);
        
        fasp_dcoo_shift_read(filename1, &A);
        
        
        printf("-------------------------\n");
        printf("Read right hand size b\n");
        printf("-------------------------\n");
        dvector b_temp;
        strncpy(filename2,inpar.workdir,128);
        datafile2="/5layers/b.dat";
        strcat(filename2,datafile2);
        
        fasp_dvec_read(filename2, &b_temp);
        
        NumLayers = 5;
        INT l;
        
        dCSRmat Ai;
        
        // -----------------------------------------
        // read in the matrix for the preconditioner
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Read preconditioner Ai\n");
        printf("-------------------------\n");
    
        strncpy(filename1,inpar.workdir,128);
        datafile1="/5layers/Ai.dat";
        strcat(filename1,datafile1);
        
        fasp_dcoo_shift_read(filename1, &Ai);
        
        // -----------------------------------------
        // read in the index for forming the blocks
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Read gloabl index\n");
        printf("-------------------------\n");
        ivector *global_idx = (ivector *)fasp_mem_calloc(NumLayers,sizeof(ivector));
        {
        
        // layer 0
        strncpy(filename3,inpar.workdir,128);
        datafile3="/5layers/global_idx_0.dat";
        strcat(filename3,datafile3);
        
        fasp_ivecind_read(filename3, &global_idx[0]);
        
        // layer 1
        strncpy(filename4,inpar.workdir,128);
        datafile4="/5layers/global_idx_1.dat";
        strcat(filename4,datafile4);
        
        fasp_ivecind_read(filename4, &global_idx[1]);
        
        // layer 2
        strncpy(filename5,inpar.workdir,128);
        datafile5="/5layers/global_idx_2.dat";
        strcat(filename5,datafile5);
        
        fasp_ivecind_read(filename5, &global_idx[2]);
        
        // layer 3
        strncpy(filename6,inpar.workdir,128);
        datafile6="/5layers/global_idx_3.dat";
        strcat(filename6,datafile6);
        
        fasp_ivecind_read(filename6, &global_idx[3]);
        
        // layer 4
        strncpy(filename7,inpar.workdir,128);
        datafile7="/5layers/global_idx_4.dat";
        strcat(filename7,datafile7);
        
        fasp_ivecind_read(filename7, &global_idx[4]);
        }
        
        // ------------------------------------------------
        // form the correct index (real and imaginary part)
        // ------------------------------------------------
        printf("-------------------------\n");
        printf("Modify global index\n");
        printf("-------------------------\n");
        
        ivector *global_index = (ivector *)fasp_mem_calloc(NumLayers,sizeof(ivector));
        
        for (l=0; l<NumLayers; l++){
            
            fasp_ivec_alloc(2*global_idx[l].row, &global_index[l]);
        
            for (i=0; i<global_idx[l].row; i++){
            
                global_index[l].val[i] = global_idx[l].val[i];
                global_index[l].val[i+global_idx[l].row] = global_idx[l].val[i] + (A.row/2);
            
            }
        
        }
        
        // -----------------------------------------
        // form the tridiagonal matrix of A
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Tridiagonalize A\n");
        printf("-------------------------\n");
        
        Abcsr.brow = NumLayers; Abcsr.bcol = NumLayers;
		Abcsr.blocks = (dCSRmat **)calloc(NumLayers*NumLayers, sizeof(dCSRmat *));
        for (i=0; i<NumLayers*NumLayers; i++) {
            Abcsr.blocks[i] = (dCSRmat *)fasp_mem_calloc(1, sizeof(dCSRmat));
        }
        
        for (i=0; i<NumLayers; i++){
            
            for (j=0; j<NumLayers; j++){
                
                if (j==i) {
                    fasp_dcsr_getblk(&A, global_index[i].val, global_index[j].val, global_index[i].row, global_index[j].row, Abcsr.blocks[i*NumLayers+j]);
                }
                else if (j == (i-1)) {
                    fasp_dcsr_getblk(&A, global_index[i].val, global_index[j].val, global_index[i].row, global_index[j].row, Abcsr.blocks[i*NumLayers+j]);
                }
                else if (j == (i+1)) {
                    fasp_dcsr_getblk(&A, global_index[i].val, global_index[j].val, global_index[i].row, global_index[j].row, Abcsr.blocks[i*NumLayers+j]);
                }
                else {
                    Abcsr.blocks[i*NumLayers+j] = NULL;
                }
                    
            }
            
        }

        // -----------------------------------------
        // form b corresponding to the tridiagonal matrix
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Modify b for tridiagonlaized A\n");
        printf("-------------------------\n");
        fasp_dvec_alloc(b_temp.row, &b);
        
        INT start=0;
        for (l=0; l<NumLayers; l++){
            
            for (i=0; i<global_index[l].row; i++)
            {
                b.val[start+i] = b_temp.val[global_index[l].val[i]];
            }
            
            start = start + global_index[l].row;
        }
        
        // -----------------------------------------
        // for tridiagonal matrix of Ai
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Tridiagonalize Ai\n");
        printf("-------------------------\n");
        
        Aibcsr.brow = NumLayers; Aibcsr.bcol = NumLayers;
		Aibcsr.blocks = (dCSRmat **)calloc(NumLayers*NumLayers, sizeof(dCSRmat *));
        for (i=0; i<NumLayers*NumLayers; i++) {
            Aibcsr.blocks[i] = (dCSRmat *)fasp_mem_calloc(1, sizeof(dCSRmat));
        }
        
        for (i=0; i<NumLayers; i++){
            
            for (j=0; j<NumLayers; j++){
                
                if (j==i) {
                    fasp_dcsr_getblk(&Ai, global_index[i].val, global_index[j].val, global_index[i].row, global_index[j].row, Aibcsr.blocks[i*NumLayers+j]);
                }
                else if (j == (i-1)) {
                    fasp_dcsr_getblk(&Ai, global_index[i].val, global_index[j].val, global_index[i].row, global_index[j].row, Aibcsr.blocks[i*NumLayers+j]);
                }
                else if (j == (i+1)) {
                    fasp_dcsr_getblk(&Ai, global_index[i].val, global_index[j].val, global_index[i].row, global_index[j].row, Aibcsr.blocks[i*NumLayers+j]);
                }
                else {
                    Aibcsr.blocks[i*NumLayers+j] = NULL;
                }
                
            }
            
        }
        
        // -----------------------------------------
        // read in local A (Schur Complement approximation)
        // -----------------------------------------
        printf("----------------------------------------\n");
        printf("Read Schur complements approximations\n");
        printf("----------------------------------------\n");
        local_A = (dCSRmat *)fasp_mem_calloc(NumLayers, sizeof(dCSRmat));

        // first level is just A11 block
        local_A[0] = fasp_dcsr_create (Aibcsr.blocks[0]->row, Aibcsr.blocks[0]->col, Aibcsr.blocks[0]->nnz);
        fasp_dcsr_cp(Aibcsr.blocks[0], &(local_A[0]));
        
        // other levels
        {
        // layer 1
        strncpy(filename1,inpar.workdir,128);
        datafile1="/5layers/local_A_0.dat";
        strcat(filename1,datafile1);
        
        fasp_dcoo_shift_read(filename1, &(local_A[1]));
        
        // layer 2
        strncpy(filename2,inpar.workdir,128);
        datafile2="/5layers/local_A_1.dat";
        strcat(filename2,datafile2);
        
        fasp_dcoo_shift_read(filename2, &(local_A[2]));
        
        // layer 3
        strncpy(filename3,inpar.workdir,128);
        datafile3="/5layers/local_A_2.dat";
        strcat(filename3,datafile3);
        
        fasp_dcoo_shift_read(filename3, &(local_A[3]));
        
        // layer 4
        strncpy(filename4,inpar.workdir,128);
        datafile4="/5layers/local_A_3.dat";
        strcat(filename4,datafile4);
        
        fasp_dcoo_shift_read(filename4, &(local_A[4]));
        }
        
        // -----------------------------------------
        // read in local pointer for Schur complement
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Read local pointer\n");
        printf("-------------------------\n");
        ivector *local_ptr = (ivector *)fasp_mem_calloc(NumLayers,sizeof(ivector));
        
        {
            
            // layer 0 (no need to read)
            
            // layer 1
            strncpy(filename1,inpar.workdir,128);
            datafile1="/5layers/local_ptr_0.dat";
            strcat(filename1,datafile1);
            
            fasp_ivecind_read(filename1, &local_ptr[1]);
            
            // layer 2
            strncpy(filename2,inpar.workdir,128);
            datafile2="/5layers/local_ptr_1.dat";
            strcat(filename2,datafile2);
            
            fasp_ivecind_read(filename2, &local_ptr[2]);
            
            // layer 3
            strncpy(filename3,inpar.workdir,128);
            datafile3="/5layers/local_ptr_2.dat";
            strcat(filename3,datafile3);
            
            fasp_ivecind_read(filename3, &local_ptr[3]);
            
            // layer 4
            strncpy(filename4,inpar.workdir,128);
            datafile4="/5layers/local_ptr_3.dat";
            strcat(filename4,datafile4);
            
            fasp_ivecind_read(filename4, &local_ptr[4]);
        }
        
        // -----------------------------------------
        // generate local index for Schur complement
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Generate local index\n");
        printf("-------------------------\n");
        local_index = (ivector *)fasp_mem_calloc(NumLayers,sizeof(ivector));
        
        {
            INT local_size;
            INT local_A_size;
            
            // layer 0
            local_index[0] = fasp_ivec_create(Aibcsr.blocks[0]->row);
            
            for (i=0; i<Aibcsr.blocks[0]->row; i++) local_index[0].val[i] = i;
            
            // other layers
            for (l=1; l<NumLayers; l++){
                
                // compute local size
                local_size = local_ptr[l].val[1] - local_ptr[l].val[0];
                local_A_size = (local_A[l].row)/2;
                
                // allocate
                local_index[l] = fasp_ivec_create(local_size*2);
                
                // generate
                for (i=0; i<local_size; i++){
                    
                    local_index[l].val[i] = local_ptr[l].val[0] + i;
                    local_index[l].val[i+local_size] = local_index[l].val[i]+local_A_size;
                    
                }

            }
            
        }
        
        // -----------------------------------------
        // cleaning
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Cleaning \n");
        printf("-------------------------\n");
        
        fasp_dcsr_free(&A);
        fasp_dvec_free(&b_temp);
        fasp_dcsr_free(&Ai);
        
        for(l=0;l<NumLayers;l++) fasp_ivec_free(&(global_idx[l]));
        for(l=0;l<NumLayers;l++) fasp_ivec_free(&(local_ptr[l]));
        
        //------------------------------------
        // check
        //------------------------------------
        //fasp_dcsr_write_coo("A54.dat", Aibcsr.blocks[23]);
        //fasp_ivec_write("local_index_4.dat", &local_index[4]);
        
        
}
    
    else if (problem_num == 21) {
        
        printf("-------------------------\n");
        printf("Read matrix A\n");
        printf("-------------------------\n");
        
        strncpy(filename1,inpar.workdir,128);
        datafile1="/4layers/A.dat";
        strcat(filename1,datafile1);
        
        fasp_dcoo_shift_read(filename1, &A);
        
        
        printf("-------------------------\n");
        printf("Read right hand size b\n");
        printf("-------------------------\n");
        dvector b_temp;
        strncpy(filename2,inpar.workdir,128);
        datafile2="/4layers/b.dat";
        strcat(filename2,datafile2);
        
        fasp_dvec_read(filename2, &b_temp);
        
        NumLayers = 4;
        INT l;
        
        dCSRmat Ai;
        
        // -----------------------------------------
        // read in the matrix for the preconditioner
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Read preconditioner Ai\n");
        printf("-------------------------\n");
        
        strncpy(filename1,inpar.workdir,128);
        datafile1="/4layers/Ai.dat";
        strcat(filename1,datafile1);
        
        fasp_dcoo_shift_read(filename1, &Ai);
        
        // -----------------------------------------
        // read in the index for forming the blocks
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Read gloabl index\n");
        printf("-------------------------\n");
        ivector *global_idx = (ivector *)fasp_mem_calloc(NumLayers,sizeof(ivector));
        {
            
            // layer 0
            strncpy(filename3,inpar.workdir,128);
            datafile3="/4layers/global_idx_0.dat";
            strcat(filename3,datafile3);
            
            fasp_ivecind_read(filename3, &global_idx[0]);
            
            // layer 1
            strncpy(filename4,inpar.workdir,128);
            datafile4="/4layers/global_idx_1.dat";
            strcat(filename4,datafile4);
            
            fasp_ivecind_read(filename4, &global_idx[1]);
            
            // layer 2
            strncpy(filename5,inpar.workdir,128);
            datafile5="/4layers/global_idx_2.dat";
            strcat(filename5,datafile5);
            
            fasp_ivecind_read(filename5, &global_idx[2]);
            
            // layer 3
            strncpy(filename6,inpar.workdir,128);
            datafile6="/4layers/global_idx_3.dat";
            strcat(filename6,datafile6);
            
            fasp_ivecind_read(filename6, &global_idx[3]);
            
        }
        
        // ------------------------------------------------
        // form the correct index (real and imaginary part)
        // ------------------------------------------------
        printf("-------------------------\n");
        printf("Modify global index\n");
        printf("-------------------------\n");
        
        ivector *global_index = (ivector *)fasp_mem_calloc(NumLayers,sizeof(ivector));
        
        for (l=0; l<NumLayers; l++){
            
            fasp_ivec_alloc(2*global_idx[l].row, &global_index[l]);
            
            for (i=0; i<global_idx[l].row; i++){
                
                global_index[l].val[i] = global_idx[l].val[i];
                global_index[l].val[i+global_idx[l].row] = global_idx[l].val[i] + (A.row/2);
                
            }
            
        }
        
        // -----------------------------------------
        // form the tridiagonal matrix of A
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Tridiagonalize A\n");
        printf("-------------------------\n");
        
        Abcsr.brow = NumLayers; Abcsr.bcol = NumLayers;
		Abcsr.blocks = (dCSRmat **)calloc(NumLayers*NumLayers, sizeof(dCSRmat *));
        for (i=0; i<NumLayers*NumLayers; i++) {
            Abcsr.blocks[i] = (dCSRmat *)fasp_mem_calloc(1, sizeof(dCSRmat));
        }
        
        for (i=0; i<NumLayers; i++){
            
            for (j=0; j<NumLayers; j++){
                
                if (j==i) {
                    fasp_dcsr_getblk(&A, global_index[i].val, global_index[j].val, global_index[i].row, global_index[j].row, Abcsr.blocks[i*NumLayers+j]);
                }
                else if (j == (i-1)) {
                    fasp_dcsr_getblk(&A, global_index[i].val, global_index[j].val, global_index[i].row, global_index[j].row, Abcsr.blocks[i*NumLayers+j]);
                }
                else if (j == (i+1)) {
                    fasp_dcsr_getblk(&A, global_index[i].val, global_index[j].val, global_index[i].row, global_index[j].row, Abcsr.blocks[i*NumLayers+j]);
                }
                else {
                    Abcsr.blocks[i*NumLayers+j] = NULL;
                }
                
            }
            
        }
        
        // -----------------------------------------
        // form b corresponding to the tridiagonal matrix
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Modify b for tridiagonlaized A\n");
        printf("-------------------------\n");
        fasp_dvec_alloc(b_temp.row, &b);
        
        INT start=0;
        for (l=0; l<NumLayers; l++){
            
            for (i=0; i<global_index[l].row; i++)
            {
                b.val[start+i] = b_temp.val[global_index[l].val[i]];
            }
            
            start = start + global_index[l].row;
        }
        
        // -----------------------------------------
        // for tridiagonal matrix of Ai
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Tridiagonalize Ai\n");
        printf("-------------------------\n");
        
        Aibcsr.brow = NumLayers; Aibcsr.bcol = NumLayers;
		Aibcsr.blocks = (dCSRmat **)calloc(NumLayers*NumLayers, sizeof(dCSRmat *));
        for (i=0; i<NumLayers*NumLayers; i++) {
            Aibcsr.blocks[i] = (dCSRmat *)fasp_mem_calloc(1, sizeof(dCSRmat));
        }
        
        for (i=0; i<NumLayers; i++){
            
            for (j=0; j<NumLayers; j++){
                
                if (j==i) {
                    fasp_dcsr_getblk(&Ai, global_index[i].val, global_index[j].val, global_index[i].row, global_index[j].row, Aibcsr.blocks[i*NumLayers+j]);
                }
                else if (j == (i-1)) {
                    fasp_dcsr_getblk(&Ai, global_index[i].val, global_index[j].val, global_index[i].row, global_index[j].row, Aibcsr.blocks[i*NumLayers+j]);
                }
                else if (j == (i+1)) {
                    fasp_dcsr_getblk(&Ai, global_index[i].val, global_index[j].val, global_index[i].row, global_index[j].row, Aibcsr.blocks[i*NumLayers+j]);
                }
                else {
                    Aibcsr.blocks[i*NumLayers+j] = NULL;
                }
                
            }
            
        }
        
        // -----------------------------------------
        // read in local A (Schur Complement approximation)
        // -----------------------------------------
        printf("----------------------------------------\n");
        printf("Read Schur complements approximations\n");
        printf("----------------------------------------\n");
        local_A = (dCSRmat *)fasp_mem_calloc(NumLayers, sizeof(dCSRmat));
        
        // first level is just A11 block
        local_A[0] = fasp_dcsr_create (Aibcsr.blocks[0]->row, Aibcsr.blocks[0]->col, Aibcsr.blocks[0]->nnz);
        fasp_dcsr_cp(Aibcsr.blocks[0], &(local_A[0]));
        
        // other levels
        {
            // layer 1
            strncpy(filename1,inpar.workdir,128);
            datafile1="/4layers/local_A_0.dat";
            strcat(filename1,datafile1);
            
            fasp_dcoo_shift_read(filename1, &(local_A[1]));
            
            // layer 2
            strncpy(filename2,inpar.workdir,128);
            datafile2="/4layers/local_A_1.dat";
            strcat(filename2,datafile2);
            
            fasp_dcoo_shift_read(filename2, &(local_A[2]));
            
            // layer 3
            strncpy(filename3,inpar.workdir,128);
            datafile3="/4layers/local_A_2.dat";
            strcat(filename3,datafile3);
            
            fasp_dcoo_shift_read(filename3, &(local_A[3]));
            
        }
        
        // -----------------------------------------
        // read in local pointer for Schur complement
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Read local pointer\n");
        printf("-------------------------\n");
        ivector *local_ptr = (ivector *)fasp_mem_calloc(NumLayers,sizeof(ivector));
        
        {
            
            // layer 0 (no need to read)
            
            // layer 1
            strncpy(filename1,inpar.workdir,128);
            datafile1="/4layers/local_ptr_0.dat";
            strcat(filename1,datafile1);
            
            fasp_ivecind_read(filename1, &local_ptr[1]);
            
            // layer 2
            strncpy(filename2,inpar.workdir,128);
            datafile2="/4layers/local_ptr_1.dat";
            strcat(filename2,datafile2);
            
            fasp_ivecind_read(filename2, &local_ptr[2]);
            
            // layer 3
            strncpy(filename3,inpar.workdir,128);
            datafile3="/4layers/local_ptr_2.dat";
            strcat(filename3,datafile3);
            
            fasp_ivecind_read(filename3, &local_ptr[3]);
            
        }
        
        // -----------------------------------------
        // generate local index for Schur complement
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Generate local index\n");
        printf("-------------------------\n");
        local_index = (ivector *)fasp_mem_calloc(NumLayers,sizeof(ivector));
        
        {
            INT local_size;
            INT local_A_size;
            
            // layer 0
            local_index[0] = fasp_ivec_create(Aibcsr.blocks[0]->row);
            
            for (i=0; i<Aibcsr.blocks[0]->row; i++) local_index[0].val[i] = i;
            
            // other layers
            for (l=1; l<NumLayers; l++){
                
                // compute local size
                local_size = local_ptr[l].val[1] - local_ptr[l].val[0];
                local_A_size = (local_A[l].row)/2;
                
                // allocate
                local_index[l] = fasp_ivec_create(local_size*2);
                
                // generate
                for (i=0; i<local_size; i++){
                    
                    local_index[l].val[i] = local_ptr[l].val[0] + i;
                    local_index[l].val[i+local_size] = local_index[l].val[i]+local_A_size;
                    
                }
                
            }
            
        }
        
        // -----------------------------------------
        // cleaning
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Cleaning \n");
        printf("-------------------------\n");
        
        fasp_dcsr_free(&A);
        fasp_dvec_free(&b_temp);
        fasp_dcsr_free(&Ai);
        
        for(l=0;l<NumLayers;l++) fasp_ivec_free(&(global_idx[l]));
        for(l=0;l<NumLayers;l++) fasp_ivec_free(&(local_ptr[l]));
        
        //------------------------------------
        // check
        //------------------------------------
        //fasp_dcsr_write_coo("A54.dat", Aibcsr.blocks[23]);
        //fasp_ivec_write("local_index_4.dat", &local_index[4]);
        
        
    }
    
    else if (problem_num == 22) {
        
        printf("-------------------------\n");
        printf("Read matrix A\n");
        printf("-------------------------\n");
        
        strncpy(filename1,inpar.workdir,128);
        datafile1="/3layers/A.dat";
        strcat(filename1,datafile1);
        
        fasp_dcoo_shift_read(filename1, &A);
        
        
        printf("-------------------------\n");
        printf("Read right hand size b\n");
        printf("-------------------------\n");
        dvector b_temp;
        strncpy(filename2,inpar.workdir,128);
        datafile2="/3layers/b.dat";
        strcat(filename2,datafile2);
        
        fasp_dvec_read(filename2, &b_temp);
        
        NumLayers = 3;
        INT l;
        
        dCSRmat Ai;
        
        // -----------------------------------------
        // read in the matrix for the preconditioner
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Read preconditioner Ai\n");
        printf("-------------------------\n");
        
        strncpy(filename1,inpar.workdir,128);
        datafile1="/3layers/Ai.dat";
        strcat(filename1,datafile1);
        
        fasp_dcoo_shift_read(filename1, &Ai);
        
        // -----------------------------------------
        // read in the index for forming the blocks
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Read gloabl index\n");
        printf("-------------------------\n");
        ivector *global_idx = (ivector *)fasp_mem_calloc(NumLayers,sizeof(ivector));
        {
            
            // layer 0
            strncpy(filename3,inpar.workdir,128);
            datafile3="/3layers/global_idx_0.dat";
            strcat(filename3,datafile3);
            
            fasp_ivecind_read(filename3, &global_idx[0]);
            
            // layer 1
            strncpy(filename4,inpar.workdir,128);
            datafile4="/3layers/global_idx_1.dat";
            strcat(filename4,datafile4);
            
            fasp_ivecind_read(filename4, &global_idx[1]);
            
            // layer 2
            strncpy(filename5,inpar.workdir,128);
            datafile5="/3layers/global_idx_2.dat";
            strcat(filename5,datafile5);
            
            fasp_ivecind_read(filename5, &global_idx[2]);
            
            
        }
        
        // ------------------------------------------------
        // form the correct index (real and imaginary part)
        // ------------------------------------------------
        printf("-------------------------\n");
        printf("Modify global index\n");
        printf("-------------------------\n");
        
        ivector *global_index = (ivector *)fasp_mem_calloc(NumLayers,sizeof(ivector));
        
        for (l=0; l<NumLayers; l++){
            
            fasp_ivec_alloc(2*global_idx[l].row, &global_index[l]);
            
            for (i=0; i<global_idx[l].row; i++){
                
                global_index[l].val[i] = global_idx[l].val[i];
                global_index[l].val[i+global_idx[l].row] = global_idx[l].val[i] + (A.row/2);
                
            }
            
        }
        
        // -----------------------------------------
        // form the tridiagonal matrix of A
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Tridiagonalize A\n");
        printf("-------------------------\n");
        
        Abcsr.brow = NumLayers; Abcsr.bcol = NumLayers;
		Abcsr.blocks = (dCSRmat **)calloc(NumLayers*NumLayers, sizeof(dCSRmat *));
        for (i=0; i<NumLayers*NumLayers; i++) {
            Abcsr.blocks[i] = (dCSRmat *)fasp_mem_calloc(1, sizeof(dCSRmat));
        }
        
        for (i=0; i<NumLayers; i++){
            
            for (j=0; j<NumLayers; j++){
                
                if (j==i) {
                    fasp_dcsr_getblk(&A, global_index[i].val, global_index[j].val, global_index[i].row, global_index[j].row, Abcsr.blocks[i*NumLayers+j]);
                }
                else if (j == (i-1)) {
                    fasp_dcsr_getblk(&A, global_index[i].val, global_index[j].val, global_index[i].row, global_index[j].row, Abcsr.blocks[i*NumLayers+j]);
                }
                else if (j == (i+1)) {
                    fasp_dcsr_getblk(&A, global_index[i].val, global_index[j].val, global_index[i].row, global_index[j].row, Abcsr.blocks[i*NumLayers+j]);
                }
                else {
                    Abcsr.blocks[i*NumLayers+j] = NULL;
                }
                
            }
            
        }
        
        // -----------------------------------------
        // form b corresponding to the tridiagonal matrix
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Modify b for tridiagonlaized A\n");
        printf("-------------------------\n");
        fasp_dvec_alloc(b_temp.row, &b);
        
        INT start=0;
        for (l=0; l<NumLayers; l++){
            
            for (i=0; i<global_index[l].row; i++)
            {
                b.val[start+i] = b_temp.val[global_index[l].val[i]];
            }
            
            start = start + global_index[l].row;
        }
        
        // -----------------------------------------
        // for tridiagonal matrix of Ai
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Tridiagonalize Ai\n");
        printf("-------------------------\n");
        
        Aibcsr.brow = NumLayers; Aibcsr.bcol = NumLayers;
		Aibcsr.blocks = (dCSRmat **)calloc(NumLayers*NumLayers, sizeof(dCSRmat *));
        for (i=0; i<NumLayers*NumLayers; i++) {
            Aibcsr.blocks[i] = (dCSRmat *)fasp_mem_calloc(1, sizeof(dCSRmat));
        }
        
        for (i=0; i<NumLayers; i++){
            
            for (j=0; j<NumLayers; j++){
                
                if (j==i) {
                    fasp_dcsr_getblk(&Ai, global_index[i].val, global_index[j].val, global_index[i].row, global_index[j].row, Aibcsr.blocks[i*NumLayers+j]);
                }
                else if (j == (i-1)) {
                    fasp_dcsr_getblk(&Ai, global_index[i].val, global_index[j].val, global_index[i].row, global_index[j].row, Aibcsr.blocks[i*NumLayers+j]);
                }
                else if (j == (i+1)) {
                    fasp_dcsr_getblk(&Ai, global_index[i].val, global_index[j].val, global_index[i].row, global_index[j].row, Aibcsr.blocks[i*NumLayers+j]);
                }
                else {
                    Aibcsr.blocks[i*NumLayers+j] = NULL;
                }
                
            }
            
        }
        
        // -----------------------------------------
        // read in local A (Schur Complement approximation)
        // -----------------------------------------
        printf("----------------------------------------\n");
        printf("Read Schur complements approximations\n");
        printf("----------------------------------------\n");
        local_A = (dCSRmat *)fasp_mem_calloc(NumLayers, sizeof(dCSRmat));
        
        // first level is just A11 block
        local_A[0] = fasp_dcsr_create (Aibcsr.blocks[0]->row, Aibcsr.blocks[0]->col, Aibcsr.blocks[0]->nnz);
        fasp_dcsr_cp(Aibcsr.blocks[0], &(local_A[0]));
        
        // other levels
        {
            // layer 1
            strncpy(filename1,inpar.workdir,128);
            datafile1="/3layers/local_A_0.dat";
            strcat(filename1,datafile1);
            
            fasp_dcoo_shift_read(filename1, &(local_A[1]));
            
            // layer 2
            strncpy(filename2,inpar.workdir,128);
            datafile2="/3layers/local_A_1.dat";
            strcat(filename2,datafile2);
            
            fasp_dcoo_shift_read(filename2, &(local_A[2]));
            
            
        }
        
        // -----------------------------------------
        // read in local pointer for Schur complement
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Read local pointer\n");
        printf("-------------------------\n");
        ivector *local_ptr = (ivector *)fasp_mem_calloc(NumLayers,sizeof(ivector));
        
        {
            
            // layer 0 (no need to read)
            
            // layer 1
            strncpy(filename1,inpar.workdir,128);
            datafile1="/3layers/local_ptr_0.dat";
            strcat(filename1,datafile1);
            
            fasp_ivecind_read(filename1, &local_ptr[1]);
            
            // layer 2
            strncpy(filename2,inpar.workdir,128);
            datafile2="/3layers/local_ptr_1.dat";
            strcat(filename2,datafile2);
            
            fasp_ivecind_read(filename2, &local_ptr[2]);
            
           
            
        }
        
        // -----------------------------------------
        // generate local index for Schur complement
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Generate local index\n");
        printf("-------------------------\n");
        local_index = (ivector *)fasp_mem_calloc(NumLayers,sizeof(ivector));
        
        {
            INT local_size;
            INT local_A_size;
            
            // layer 0
            local_index[0] = fasp_ivec_create(Aibcsr.blocks[0]->row);
            
            for (i=0; i<Aibcsr.blocks[0]->row; i++) local_index[0].val[i] = i;
            
            // other layers
            for (l=1; l<NumLayers; l++){
                
                // compute local size
                local_size = local_ptr[l].val[1] - local_ptr[l].val[0];
                local_A_size = (local_A[l].row)/2;
                
                // allocate
                local_index[l] = fasp_ivec_create(local_size*2);
                
                // generate
                for (i=0; i<local_size; i++){
                    
                    local_index[l].val[i] = local_ptr[l].val[0] + i;
                    local_index[l].val[i+local_size] = local_index[l].val[i]+local_A_size;
                    
                }
                
            }
            
        }
        
        // -----------------------------------------
        // cleaning
        // -----------------------------------------
        printf("-------------------------\n");
        printf("Cleaning \n");
        printf("-------------------------\n");
        
        fasp_dcsr_free(&A);
        fasp_dvec_free(&b_temp);
        fasp_dcsr_free(&Ai);
        
        for(l=0;l<NumLayers;l++) fasp_ivec_free(&(global_idx[l]));
        for(l=0;l<NumLayers;l++) fasp_ivec_free(&(local_ptr[l]));
        
        //------------------------------------
        // check
        //------------------------------------
        //fasp_dcsr_write_coo("A54.dat", Aibcsr.blocks[23]);
        //fasp_ivec_write("local_index_4.dat", &local_index[4]);
        
        
    }
    
    // test problem for MHD benchmark
    else if (problem_num == 30) {
        
        //----------------------------------
        // read in data
        //----------------------------------
        
        // read in A
        strncpy(filename1,inpar.workdir,128);
        datafile1="/MHD/t0.05_h0.5/A.dat";
        strcat(filename1,datafile1);
        
        fasp_dcoo_read(filename1, &A);
        
        // read in b
        strncpy(filename2,inpar.workdir,128);
        datafile2="/MHD/t0.05_h0.5/b.dat";
        strcat(filename2,datafile2);
        
        fasp_dvec_read(filename2, &b);
        
        //---------------------------------
        // read in diagonal block data
        //---------------------------------
        A_diag = (dCSRmat *)fasp_mem_calloc(4, sizeof(dCSRmat));
        
        // read in A_U
        strncpy(filename3,inpar.workdir,128);
        datafile3="/MHD/t0.05_h0.5/H1H1.dat";
        strcat(filename3, datafile3);
        
        fasp_dcoo_read(filename3, &(A_diag[0]));
        
        // read in M_P
        dCSRmat M_P;
        strncpy(filename4,inpar.workdir,128);
        datafile4="/MHD/t0.05_h0.5/L2L2.dat";
        strcat(filename4,datafile4);
        
        fasp_dcoo_read(filename4, &(A_diag[1]));
        
        // read in A_E
        dCSRmat A_E;
        strncpy(filename5,inpar.workdir,128);
        datafile5="/MHD/t0.05_h0.5/Hcurl.dat";
        strcat(filename5,datafile5);
        
        fasp_dcoo_read(filename5, &(A_diag[2]));
        
        // read in M_B
        dCSRmat M_B;
        strncpy(filename6,inpar.workdir,128);
        datafile6="/MHD/t0.05_h0.5/Hdiv.dat";
        strcat(filename6,datafile6);
        
        fasp_dcoo_read(filename6, &(A_diag[3]));
        
        // form block CSR matrix
		ivector U_idx;
        ivector P_idx;
        ivector E_idx;
        ivector B_idx;
        
        fasp_ivec_alloc(A_diag[0].row, &U_idx);
        fasp_ivec_alloc(A_diag[1].row, &P_idx);
        fasp_ivec_alloc(A_diag[2].row, &E_idx);
        fasp_ivec_alloc(A_diag[3].row, &B_idx);
        
        for (i=0; i<A_diag[0].row; i++){
            U_idx.val[i] = i;
        }
        
        for (i=0; i<A_diag[1].row; i++){
            P_idx.val[i] = A_diag[0].row+i;
        }
    
        for (i=0; i<A_diag[2].row; i++){
            E_idx.val[i] = A_diag[0].row+A_diag[1].row+i;
        }
        
        for (i=0; i<A_diag[3].row; i++){
            B_idx.val[i] = A_diag[0].row+A_diag[1].row+A_diag[2].row+i;
        }
        
        // Assemble the matrix in block dCSR format
		Abcsr.brow = 4; Abcsr.bcol = 4;
		Abcsr.blocks = (dCSRmat **)calloc(16, sizeof(dCSRmat *));
        for (i=0; i<16 ;i++) {
            Abcsr.blocks[i] = (dCSRmat *)fasp_mem_calloc(1, sizeof(dCSRmat));
        }
        
		// A11
        fasp_dcsr_getblk(&A, U_idx.val, U_idx.val, U_idx.row, U_idx.row, Abcsr.blocks[0]);
        // A12
        fasp_dcsr_getblk(&A, U_idx.val, P_idx.val, U_idx.row, P_idx.row, Abcsr.blocks[1]);
        // A13
        fasp_dcsr_getblk(&A, U_idx.val, E_idx.val, U_idx.row, E_idx.row, Abcsr.blocks[2]);
        // A14
        fasp_dcsr_getblk(&A, U_idx.val, B_idx.val, U_idx.row, B_idx.row, Abcsr.blocks[3]);
        // A21
        fasp_dcsr_getblk(&A, P_idx.val, U_idx.val, P_idx.row, U_idx.row, Abcsr.blocks[4]);
        // A22
        fasp_dcsr_getblk(&A, P_idx.val, P_idx.val, P_idx.row, P_idx.row, Abcsr.blocks[5]);
        // A23
        fasp_dcsr_getblk(&A, P_idx.val, E_idx.val, P_idx.row, E_idx.row, Abcsr.blocks[6]);
        // A24
        fasp_dcsr_getblk(&A, P_idx.val, B_idx.val, P_idx.row, B_idx.row, Abcsr.blocks[7]);
        // A31
        fasp_dcsr_getblk(&A, E_idx.val, U_idx.val, E_idx.row, U_idx.row, Abcsr.blocks[8]);
        // A32
        fasp_dcsr_getblk(&A, E_idx.val, P_idx.val, E_idx.row, P_idx.row, Abcsr.blocks[9]);
        // A33
        fasp_dcsr_getblk(&A, E_idx.val, E_idx.val, E_idx.row, E_idx.row, Abcsr.blocks[10]);
        // A34
        fasp_dcsr_getblk(&A, E_idx.val, B_idx.val, E_idx.row, B_idx.row, Abcsr.blocks[11]);
        // A41
        fasp_dcsr_getblk(&A, B_idx.val, U_idx.val, B_idx.row, U_idx.row, Abcsr.blocks[12]);
        // A42
        fasp_dcsr_getblk(&A, B_idx.val, P_idx.val, B_idx.row, P_idx.row, Abcsr.blocks[13]);
        // A43
        fasp_dcsr_getblk(&A, B_idx.val, E_idx.val, B_idx.row, E_idx.row, Abcsr.blocks[14]);
        // A44
        fasp_dcsr_getblk(&A, B_idx.val, B_idx.val, B_idx.row, B_idx.row, Abcsr.blocks[15]);
        
        // use A0 block in the preconditioner
        fasp_dcsr_cp(Abcsr.blocks[0], &(A_diag[0]));
        
    }
    
    else if (problem_num == 31) {
        
        //----------------------------------
        // read in data
        //----------------------------------
        
        // read in A
        strncpy(filename1,inpar.workdir,128);
        datafile1="/MHD/t0.05_h0.25/A.dat";
        strcat(filename1,datafile1);
        
        fasp_dcoo_read(filename1, &A);
        
        // read in b
        strncpy(filename2,inpar.workdir,128);
        datafile2="/MHD/t0.05_h0.25/b.dat";
        strcat(filename2,datafile2);
        
        fasp_dvec_read(filename2, &b);
        
        //---------------------------------
        // read in diagonal block data
        //---------------------------------
        A_diag = (dCSRmat *)fasp_mem_calloc(4, sizeof(dCSRmat));
        
        // read in A_U
        strncpy(filename3,inpar.workdir,128);
        datafile3="/MHD/t0.05_h0.25/H1H1.dat";
        strcat(filename3, datafile3);
        
        fasp_dcoo_read(filename3, &(A_diag[0]));
        
        // read in M_P
        dCSRmat M_P;
        strncpy(filename4,inpar.workdir,128);
        datafile4="/MHD/t0.05_h0.25/L2L2.dat";
        strcat(filename4,datafile4);
        
        fasp_dcoo_read(filename4, &(A_diag[1]));
        
        // read in A_E
        dCSRmat A_E;
        strncpy(filename5,inpar.workdir,128);
        datafile5="/MHD/t0.05_h0.25/Hcurl.dat";
        strcat(filename5,datafile5);
        
        fasp_dcoo_read(filename5, &(A_diag[2]));
        
        // read in M_B
        dCSRmat M_B;
        strncpy(filename6,inpar.workdir,128);
        datafile6="/MHD/t0.05_h0.25/Hdiv.dat";
        strcat(filename6,datafile6);
        
        fasp_dcoo_read(filename6, &(A_diag[3]));
        
        // form block CSR matrix
		ivector U_idx;
        ivector P_idx;
        ivector E_idx;
        ivector B_idx;
        
        fasp_ivec_alloc(A_diag[0].row, &U_idx);
        fasp_ivec_alloc(A_diag[1].row, &P_idx);
        fasp_ivec_alloc(A_diag[2].row, &E_idx);
        fasp_ivec_alloc(A_diag[3].row, &B_idx);
        
        for (i=0; i<A_diag[0].row; i++){
            U_idx.val[i] = i;
        }
        
        for (i=0; i<A_diag[1].row; i++){
            P_idx.val[i] = A_diag[0].row+i;
        }
        
        for (i=0; i<A_diag[2].row; i++){
            E_idx.val[i] = A_diag[0].row+A_diag[1].row+i;
        }
        
        for (i=0; i<A_diag[3].row; i++){
            B_idx.val[i] = A_diag[0].row+A_diag[1].row+A_diag[2].row+i;
        }
        
        // Assemble the matrix in block dCSR format
		Abcsr.brow = 4; Abcsr.bcol = 4;
		Abcsr.blocks = (dCSRmat **)calloc(16, sizeof(dCSRmat *));
        for (i=0; i<16 ;i++) {
            Abcsr.blocks[i] = (dCSRmat *)fasp_mem_calloc(1, sizeof(dCSRmat));
        }
        
		// A11
        fasp_dcsr_getblk(&A, U_idx.val, U_idx.val, U_idx.row, U_idx.row, Abcsr.blocks[0]);
        // A12
        fasp_dcsr_getblk(&A, U_idx.val, P_idx.val, U_idx.row, P_idx.row, Abcsr.blocks[1]);
        // A13
        fasp_dcsr_getblk(&A, U_idx.val, E_idx.val, U_idx.row, E_idx.row, Abcsr.blocks[2]);
        // A14
        fasp_dcsr_getblk(&A, U_idx.val, B_idx.val, U_idx.row, B_idx.row, Abcsr.blocks[3]);
        // A21
        fasp_dcsr_getblk(&A, P_idx.val, U_idx.val, P_idx.row, U_idx.row, Abcsr.blocks[4]);
        // A22
        fasp_dcsr_getblk(&A, P_idx.val, P_idx.val, P_idx.row, P_idx.row, Abcsr.blocks[5]);
        // A23
        fasp_dcsr_getblk(&A, P_idx.val, E_idx.val, P_idx.row, E_idx.row, Abcsr.blocks[6]);
        // A24
        fasp_dcsr_getblk(&A, P_idx.val, B_idx.val, P_idx.row, B_idx.row, Abcsr.blocks[7]);
        // A31
        fasp_dcsr_getblk(&A, E_idx.val, U_idx.val, E_idx.row, U_idx.row, Abcsr.blocks[8]);
        // A32
        fasp_dcsr_getblk(&A, E_idx.val, P_idx.val, E_idx.row, P_idx.row, Abcsr.blocks[9]);
        // A33
        fasp_dcsr_getblk(&A, E_idx.val, E_idx.val, E_idx.row, E_idx.row, Abcsr.blocks[10]);
        // A34
        fasp_dcsr_getblk(&A, E_idx.val, B_idx.val, E_idx.row, B_idx.row, Abcsr.blocks[11]);
        // A41
        fasp_dcsr_getblk(&A, B_idx.val, U_idx.val, B_idx.row, U_idx.row, Abcsr.blocks[12]);
        // A42
        fasp_dcsr_getblk(&A, B_idx.val, P_idx.val, B_idx.row, P_idx.row, Abcsr.blocks[13]);
        // A43
        fasp_dcsr_getblk(&A, B_idx.val, E_idx.val, B_idx.row, E_idx.row, Abcsr.blocks[14]);
        // A44
        fasp_dcsr_getblk(&A, B_idx.val, B_idx.val, B_idx.row, B_idx.row, Abcsr.blocks[15]);
        
        // use A0 block in the preconditioner
        fasp_dcsr_cp(Abcsr.blocks[0], &(A_diag[0]));
        
    }
    
    else if (problem_num == 32) {
        
        //----------------------------------
        // read in data
        //----------------------------------
        
        // read in A
        strncpy(filename1,inpar.workdir,128);
        datafile1="/MHD/t0.05_h0.18/A.dat";
        strcat(filename1,datafile1);
        
        fasp_dcoo_read(filename1, &A);
        
        // read in b
        strncpy(filename2,inpar.workdir,128);
        datafile2="/MHD/t0.05_h0.18/b.dat";
        strcat(filename2,datafile2);
        
        fasp_dvec_read(filename2, &b);
        
        //---------------------------------
        // read in diagonal block data
        //---------------------------------
        A_diag = (dCSRmat *)fasp_mem_calloc(4, sizeof(dCSRmat));
        
        // read in A_U
        strncpy(filename3,inpar.workdir,128);
        datafile3="/MHD/t0.05_h0.18/H1H1.dat";
        strcat(filename3, datafile3);
        
        fasp_dcoo_read(filename3, &(A_diag[0]));
        
        // read in M_P
        dCSRmat M_P;
        strncpy(filename4,inpar.workdir,128);
        datafile4="/MHD/t0.05_h0.18/L2L2.dat";
        strcat(filename4,datafile4);
        
        fasp_dcoo_read(filename4, &(A_diag[1]));
        
        // read in A_E
        dCSRmat A_E;
        strncpy(filename5,inpar.workdir,128);
        datafile5="/MHD/t0.05_h0.18/Hcurl.dat";
        strcat(filename5,datafile5);
        
        fasp_dcoo_read(filename5, &(A_diag[2]));
        
        // read in M_B
        dCSRmat M_B;
        strncpy(filename6,inpar.workdir,128);
        datafile6="/MHD/t0.05_h0.18/Hdiv.dat";
        strcat(filename6,datafile6);
        
        fasp_dcoo_read(filename6, &(A_diag[3]));
        
        // form block CSR matrix
		ivector U_idx;
        ivector P_idx;
        ivector E_idx;
        ivector B_idx;
        
        fasp_ivec_alloc(A_diag[0].row, &U_idx);
        fasp_ivec_alloc(A_diag[1].row, &P_idx);
        fasp_ivec_alloc(A_diag[2].row, &E_idx);
        fasp_ivec_alloc(A_diag[3].row, &B_idx);
        
        for (i=0; i<A_diag[0].row; i++){
            U_idx.val[i] = i;
        }
        
        for (i=0; i<A_diag[1].row; i++){
            P_idx.val[i] = A_diag[0].row+i;
        }
        
        for (i=0; i<A_diag[2].row; i++){
            E_idx.val[i] = A_diag[0].row+A_diag[1].row+i;
        }
        
        for (i=0; i<A_diag[3].row; i++){
            B_idx.val[i] = A_diag[0].row+A_diag[1].row+A_diag[2].row+i;
        }
        
        // Assemble the matrix in block dCSR format
		Abcsr.brow = 4; Abcsr.bcol = 4;
		Abcsr.blocks = (dCSRmat **)calloc(16, sizeof(dCSRmat *));
        for (i=0; i<16 ;i++) {
            Abcsr.blocks[i] = (dCSRmat *)fasp_mem_calloc(1, sizeof(dCSRmat));
        }
        
		// A11
        fasp_dcsr_getblk(&A, U_idx.val, U_idx.val, U_idx.row, U_idx.row, Abcsr.blocks[0]);
        // A12
        fasp_dcsr_getblk(&A, U_idx.val, P_idx.val, U_idx.row, P_idx.row, Abcsr.blocks[1]);
        // A13
        fasp_dcsr_getblk(&A, U_idx.val, E_idx.val, U_idx.row, E_idx.row, Abcsr.blocks[2]);
        // A14
        fasp_dcsr_getblk(&A, U_idx.val, B_idx.val, U_idx.row, B_idx.row, Abcsr.blocks[3]);
        // A21
        fasp_dcsr_getblk(&A, P_idx.val, U_idx.val, P_idx.row, U_idx.row, Abcsr.blocks[4]);
        // A22
        fasp_dcsr_getblk(&A, P_idx.val, P_idx.val, P_idx.row, P_idx.row, Abcsr.blocks[5]);
        // A23
        fasp_dcsr_getblk(&A, P_idx.val, E_idx.val, P_idx.row, E_idx.row, Abcsr.blocks[6]);
        // A24
        fasp_dcsr_getblk(&A, P_idx.val, B_idx.val, P_idx.row, B_idx.row, Abcsr.blocks[7]);
        // A31
        fasp_dcsr_getblk(&A, E_idx.val, U_idx.val, E_idx.row, U_idx.row, Abcsr.blocks[8]);
        // A32
        fasp_dcsr_getblk(&A, E_idx.val, P_idx.val, E_idx.row, P_idx.row, Abcsr.blocks[9]);
        // A33
        fasp_dcsr_getblk(&A, E_idx.val, E_idx.val, E_idx.row, E_idx.row, Abcsr.blocks[10]);
        // A34
        fasp_dcsr_getblk(&A, E_idx.val, B_idx.val, E_idx.row, B_idx.row, Abcsr.blocks[11]);
        // A41
        fasp_dcsr_getblk(&A, B_idx.val, U_idx.val, B_idx.row, U_idx.row, Abcsr.blocks[12]);
        // A42
        fasp_dcsr_getblk(&A, B_idx.val, P_idx.val, B_idx.row, P_idx.row, Abcsr.blocks[13]);
        // A43
        fasp_dcsr_getblk(&A, B_idx.val, E_idx.val, B_idx.row, E_idx.row, Abcsr.blocks[14]);
        // A44
        fasp_dcsr_getblk(&A, B_idx.val, B_idx.val, B_idx.row, B_idx.row, Abcsr.blocks[15]);
        
        // use A0 block in the preconditioner
        fasp_dcsr_cp(Abcsr.blocks[0], &(A_diag[0]));
        
    }
    
    else {
		printf("### ERROR: Unrecognized problem number %d\n", problem_num);
		return ERROR_INPUT_PAR;
	}
    
	// Step 2. Solve the system
    
    // Print out solver parameters
    if (print_level>PRINT_NONE) fasp_param_solver_print(&itpar);
    
    // Set initial guess
    fasp_dvec_alloc(b.row, &uh); 
    fasp_dvec_set(b.row, &uh, 0.0);
    
    // Preconditioned Krylov methods
    if ( itsolver_type > 0 && itsolver_type < 20 ) {
        
		// Using no preconditioner for Krylov iterative methods
		if (precond_type == PREC_NULL) {
            status = fasp_solver_bdcsr_krylov(&Abcsr, &b, &uh, &itpar);
		}	
        
		// Using diag(A) as preconditioner for Krylov iterative methods
		else if (precond_type >= 20 &&  precond_type < 30) {
            
            if (Abcsr.brow == 3) {
                status = fasp_solver_bdcsr_krylov_block_3(&Abcsr, &b, &uh, &itpar, &amgpar);
            }
            else if (Abcsr.brow == 4) {
                status = fasp_solver_bdcsr_krylov_block_4(&Abcsr, &b, &uh, &itpar, &amgpar, A_diag);
            }
            else {
                //status = fasp_solver_bdcsr_krylov_block(&Abcsr, &b, &uh, &itpar, &amgpar);
            }
		}
        
        // sweeping preconditioners
        else if (precond_type >= 30 && precond_type < 40) {
            status = fasp_solver_bdcsr_krylov_sweeping(&Abcsr, &b, &uh, &itpar, NumLayers, &Aibcsr, local_A, local_index);
        }
        
		else {
			printf("### ERROR: Wrong preconditioner type %d!!!\n", precond_type);		
			exit(ERROR_SOLVER_PRECTYPE);
		}
        
	}
    
	else {
		printf("### ERROR: Wrong solver type %d!!!\n", itsolver_type);		
		status = ERROR_SOLVER_TYPE;
        goto FINISHED;
	}
    	
	if (status<0) {
		printf("\n### ERROR: Solver failed! Exit status = %d.\n\n", status);
	}
    
    if (output_type) fclose (stdout);
    
 FINISHED:
    // Clean up memory
    fasp_bdcsr_free(&Abcsr);
    //fasp_dbsr_free(&Absr);
	fasp_dvec_free(&b);
	fasp_dvec_free(&uh);
    
	return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
