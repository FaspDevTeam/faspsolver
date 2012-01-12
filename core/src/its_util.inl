/*! \file its_util.inl
 *  \brief Utilies for iterative solvers
 */

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

#define ITS_ZEROSOL printf("### WARNING: Iteration stopped due to the solution is close to zero!\n")

#define ITS_RESTART printf("### WARNING: Iteration restarted due to stagnation!\n")

#define ITS_STAGGED printf("### WARNING: Iteration stopped due to staggnation!\n")

#define ITS_ZEROTOL printf("### WARNING: The tolerence might be too small!\n")

#define ITS_REALRES(relres) printf("### WARNING: The actual relative residual = %e!\n",(relres))

#define ITS_COMPRES(relres) printf("### WARNING: The computed relative residual = %e!\n",(relres))

inline static void ITS_FINAL (const INT iter, const INT MaxIt, const REAL relres) 
{
    if (iter>MaxIt){
        printf("Maximal iteration %d reached with relative residual %e.\n", MaxIt, relres);
    }
    else if (iter >= 0)
        printf("Number of iterations = %d with relative residual %e.\n", iter, relres);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
