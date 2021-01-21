#define MAXCHARLEN 128

/*---------------------------------*/
/*--      Define Baseline        --*/
/*---------------------------------*/
typedef struct baseline* Baseline;
struct baseline{
    int num; // number of baseline problems.
    char **prob; // name of baseline problems, num x MAXCHARLEN.
    int *callnums; // number of call per baseline problems.
};

Baseline CreateBaseline(int n)
{
    Baseline bl = (Baseline) malloc(sizeof(struct baseline));
    bl->num = n;
    bl->callnums = (int *)malloc(sizeof(int)*n);
    bl->prob = (char **)malloc(sizeof(char *)*n);
    int i;
    for(i = 0; i < n; i++)
    {
        bl->prob[i] = (char *)malloc(MAXCHARLEN * sizeof(char));
    }
    return bl;
}

void FreeBaseline(Baseline bl)
{
    int i;
    free(bl->callnums); bl->callnums = NULL;
    for(i = 0; i < bl->num; i++)
    {
        free(bl->prob[i]);
    }
    free(bl->prob); bl->prob = NULL;
}


void PrintBaseline(Baseline bl)
{
    int i;
    printf("Baseline->num = %d\n", bl->num);
    for(i = 0; i < bl->num; i++)
    {
        printf("Baseline->prob[%d] = %s\n", i, bl->prob[i]);
    }
    for(i = 0; i < bl->num; i++)
    {
        printf("Baseline->callnums[%d] = %d\n", i, bl->callnums[i]);
    }
}



/*---------------------------------*/
/*--      Define Problem         --*/
/*---------------------------------*/
typedef struct problem* Problem;
struct problem{
    int num; // number of problem problems.
    char **prob; // name of problem problems, num x MAXCHARLEN.
};

Problem CreateProblem(int n)
{
    Problem pb = (Problem) malloc(sizeof(struct problem));
    pb->num = n;
    pb->prob = (char **)malloc(sizeof(char *)*n);
    int i;
    for(i = 0; i < n; i++)
    {
        pb->prob[i] = (char *)malloc(MAXCHARLEN * sizeof(char));
    }
    return pb;
}

void FreeProblem(Problem pb)
{
    int i;
    for(i = 0; i < pb->num; i++)
    {
        free(pb->prob[i]);
    }
    free(pb->prob); pb->prob = NULL;
}


void PrintProblem(Problem pb)
{
    int i;
    printf("Problem->num = %d\n", pb->num);
    for(i = 0; i < pb->num; i++)
    {
        printf("Problem->prob[%d] = %s\n", i, pb->prob[i]);
    }
}



/*---------------------------------*/
/*--      Define Algorithm         --*/
/*---------------------------------*/
typedef struct algorithm* Algorithm;
struct algorithm{
    int num; // number of algorithms.
    char **para; // name of algorithms, num x MAXCHARLEN.
};

Algorithm CreateAlgorithm(int n)
{
    Algorithm ag = (Algorithm) malloc(sizeof(struct algorithm));
    ag->num = n;
    ag->para = (char **)malloc(sizeof(char *)*n);
    int i;
    for(i = 0; i < n; i++)
    {
        ag->para[i] = (char *)malloc(MAXCHARLEN * sizeof(char));
    }
    return ag;
}

void FreeAlgorithm(Algorithm ag)
{
    int i;
    for(i = 0; i < ag->num; i++)
    {
        free(ag->para[i]);
    }
    free(ag->para); ag->para = NULL;
}


void PrintAlgorithm(Algorithm ag)
{
    int i;
    printf("Algorithm->num = %d\n", ag->num);
    for(i = 0; i < ag->num; i++)
    {
        printf("Algorithm->para[%d] = %s\n", i, ag->para[i]);
    }
}


/*---------------------------------*/
/*--       ReadInputFile         --*/
/*---------------------------------*/
int ReadInputFile(const char *filename, Baseline *blOut, Problem *pbOut, Algorithm *agOut)
{
    Baseline bl;
    Problem  pb;
    Algorithm ag;
    FILE *fpReadInput = fopen(filename, "r");
    int numBl, numPb, numAg;
    char buffer[512], bufTemp[128];
    int isBaseline = 0, isProblem = 0, isAlgorithm = 0;
    int baselineID = 0, baselineNum = 0;
    int probID;
    if(!fpReadInput)
    {
        printf("%s file does not exist!!!\n", filename);
        return -1;
    }
    
    fscanf(fpReadInput, "%d %d %d\n", &numBl, &numPb, &numAg);
    // printf("numBl = %d, numPb = %d, numAg = %d\n", numBl, numPb, numAg);
    bl = CreateBaseline(numBl);
    pb = CreateProblem(numPb);
    ag = CreateAlgorithm(numAg);

    while (!feof(fpReadInput))
    {
        fscanf(fpReadInput, "%s", buffer);
        // printf("buffer = %s\n", buffer);
        if(strcmp(buffer, "Baseline")==0){
            // printf("buffer = %s\n", buffer);
            isBaseline = 1;
            isProblem = 0;
            isAlgorithm = 0;
            fscanf(fpReadInput, "%*[^\n]");
            continue;
        }
        if(strcmp(buffer, "Problem")==0){
            // printf("buffer = %s\n", buffer);
            isBaseline = 0;
            isProblem = 1;
            isAlgorithm = 0;
            fscanf(fpReadInput, "%*[^\n]");
            continue;
        }
        if(strcmp(buffer, "Algorithm")==0){
            // printf("buffer = %s\n", buffer);
            isBaseline = 0;
            isProblem = 0;
            isAlgorithm = 1;
            fscanf(fpReadInput, "%*[^\n]");
            continue;
        }

        if(strcmp(buffer, "/")==0){
            fscanf(fpReadInput, "%*[^\n]");
        }else{ 
            // Read Baseline
            if(isBaseline){
                fscanf(fpReadInput, "%s %d\n", bufTemp, &baselineNum);
                baselineID = atoi(buffer);
                strcpy(bl->prob[baselineID - 1], bufTemp);
                bl->callnums[baselineID - 1] = baselineNum;
                // printf("baselineID = %d, baseline_prob = %s, baselineCount = %d\n", baselineID, bufTemp, baselineNum);
            }
            // Read Problem
            if(isProblem){
                fscanf(fpReadInput, "%s\n", bufTemp);
                probID = atoi(buffer);
                strcpy(pb->prob[probID - 1], bufTemp);
                // printf("probID = %d, buffer = %s\n", probID, bufTemp);
            }
            // Read Algorithm
            if(isAlgorithm){
                fscanf(fpReadInput, "%s\n", bufTemp);
                probID = atoi(buffer);
                strcpy(ag->para[probID - 1], bufTemp);
                // printf("algID = %d, buffer = %s\n", probID, bufTemp);
            }
        }      
    }
    
    // return
    *blOut = bl;
    *pbOut = pb;
    *agOut = ag;
    return 1;
} 


#if 0
int ReadInputFile(const char *filename)
{
    FILE *fpReadInput = fopen(filename, "r");
    char buffer[512], bufTemp[128];
    int isBaseline = 0, isProblem = 0, isAlgorithm = 0;
    int baselineID = 0, baselineNum = 0;
    int probID;
    if(!fpReadInput)
    {
        printf("%s file does not exist!!!\n", filename);
        return -1;
    }
    while (!feof(fpReadInput))
    {
        fscanf(fpReadInput, "%s", buffer);
        // printf("buffer = %s\n", buffer);
        if(strcmp(buffer, "Baseline")==0){
            // printf("buffer = %s\n", buffer);
            isBaseline = 1;
            isProblem = 0;
            isAlgorithm = 0;
            fscanf(fpReadInput, "%*[^\n]");
            continue;
        }
        if(strcmp(buffer, "Problem")==0){
            // printf("buffer = %s\n", buffer);
            isBaseline = 0;
            isProblem = 1;
            isAlgorithm = 0;
            fscanf(fpReadInput, "%*[^\n]");
            continue;
        }
        if(strcmp(buffer, "Algorithm")==0){
            // printf("buffer = %s\n", buffer);
            isBaseline = 0;
            isProblem = 0;
            isAlgorithm = 1;
            fscanf(fpReadInput, "%*[^\n]");
            continue;
        }

        if(strcmp(buffer, "/")==0){
            fscanf(fpReadInput, "%*[^\n]");
        }else{ 
            // Read Baseline
            if(isBaseline){
                fscanf(fpReadInput, "%s %d\n", bufTemp, &baselineNum);
                baselineID = atoi(buffer);
                printf("baselineID = %d, baseline_prob = %s, baselineCount = %d\n", baselineID, bufTemp, baselineNum);
            }
            // Read Problem
            if(isProblem){
                fscanf(fpReadInput, "%s\n", bufTemp);
                probID = atoi(buffer);
                printf("probID = %d, buffer = %s\n", probID, bufTemp);
            }
            // Read Algorithm
            if(isAlgorithm){
                fscanf(fpReadInput, "%s\n", bufTemp);
                probID = atoi(buffer);
                printf("algID = %d, buffer = %s\n", probID, bufTemp);
            }
        }      
    }
    return 1;
} 
#endif




void fasp_param_set_from_file (const int argc, const char *argv[], input_param *iniparam)
{
    int      arg_index   = 1;
    int      print_usage = FALSE;
    SHORT    status      = FASP_SUCCESS;

    // Option 1. set default input parameters
    fasp_param_input_init(iniparam);
    
    // fasp_param_input(file, iniparam);

    while ( arg_index < argc ) {

        if ( strcmp(argv[arg_index], "-help") == 0 ) {
            print_usage = TRUE; break;
        }

        // Option 2. Get parameters from an ini file
        else if ( strcmp(argv[arg_index], "-ini") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Missing ini filename! [%s]\n", __FUNCTION__);
                print_usage = TRUE; break;
            }
            strcpy(iniparam->inifile, argv[arg_index]);
            fasp_param_input(iniparam->inifile,iniparam);
            if ( ++arg_index >= argc ) break;
        }

        // Option 3. Get parameters from command line input
        else if ( strcmp(argv[arg_index], "-print") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting print level (from 0 to 10).\n");
                print_usage = TRUE; break;
            }
            iniparam->print_level = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-output") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting output type (0 or 1).\n");
                print_usage = TRUE; break;
            }
            iniparam->output_type = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-solver") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting solver type.\n");
                print_usage = TRUE; break;
            }
            iniparam->solver_type = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-precond") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting preconditioner type.\n");
                print_usage = TRUE; break;
            }
            iniparam->precond_type = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-maxit") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting max number of iterations.\n");
                print_usage = TRUE; break;
            }
            iniparam->itsolver_maxit = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-tol") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting tolerance for itsolver.\n");
                print_usage = TRUE; break;
            }
            iniparam->itsolver_tol = atof(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-amgmaxit") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting max num of iterations for AMG.\n");
                print_usage = TRUE; break;
            }
            iniparam->AMG_maxit = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-amgtol") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting tolerance for AMG.\n");
                print_usage = TRUE; break;
            }
            iniparam->AMG_tol = atof(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-amgtype") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting AMG type (1, 2, 3).\n");
                print_usage = TRUE; break;
            }
            iniparam->AMG_type = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-amgcycle") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting AMG cycle type (1, 2, 3, 12, 21).\n");
                print_usage = TRUE; break;
            }
            iniparam->AMG_cycle_type = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-amgcoarsening") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting AMG coarsening type.\n");
                print_usage = TRUE; break;
            }
            iniparam->AMG_coarsening_type = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-amginterplation") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting AMG interpolation type.\n");
                print_usage = TRUE; break;
            }
            iniparam->AMG_interpolation_type = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-amgsmoother") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting AMG smoother type.\n");
                print_usage = TRUE; break;
            }
            iniparam->AMG_smoother = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-amgsthreshold") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting AMG strong threshold.\n");
                print_usage = TRUE; break;
            }
            iniparam->AMG_strong_threshold = atof(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-amgscouple") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting AMG strong coupled threshold.\n");
                print_usage = TRUE; break;
            }
            iniparam->AMG_strong_coupled = atof(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else {
            print_usage = TRUE;
            break;
        }

    }

    if ( print_usage ) {

        printf("FASP command line options:\n");

        printf("================================================================\n");
        printf("  -ini              [CharValue] : Ini file name\n");
        printf("  -print            [IntValue]  : Print level\n");
        printf("  -output           [IntValue]  : Output to screen or a log file\n");
        printf("  -solver           [IntValue]  : Solver type\n");
        printf("  -precond          [IntValue]  : Preconditioner type\n");
        printf("  -maxit            [IntValue]  : Max number of iterations\n");
        printf("  -tol              [RealValue] : Tolerance for iterative solvers\n");
        printf("  -amgmaxit         [IntValue]  : Max number of AMG iterations\n");
        printf("  -amgtol           [RealValue] : Tolerance for AMG methods\n");
        printf("  -amgtype          [IntValue]  : AMG type\n");
        printf("  -amgcycle         [IntValue]  : AMG cycle type\n");
        printf("  -amgcoarsening    [IntValue]  : AMG coarsening type\n");
        printf("  -amginterpolation [IntValue]  : AMG interpolation type\n");
        printf("  -amgsmoother      [IntValue]  : AMG smoother type\n");
        printf("  -amgsthreshold    [RealValue] : AMG strong threshold\n");
        printf("  -amgscoupled      [RealValue] : AMG strong coupled threshold\n");
        printf("  -help                         : Brief help messages\n");

        exit(ERROR_INPUT_PAR);

    }

    // sanity checks
    status = fasp_param_check(iniparam);

    // if meet unexpected input, stop the program
    fasp_chkerr(status, __FUNCTION__);

}