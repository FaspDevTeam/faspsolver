/*! \file doxygen.h
 *  \brief Main page for Doygen documentation
 */
 
/** \mainpage FASP Solver Package
 *
 *
 * The purpose of this project is to construct a pool of problems and solvers for 
 * these problems. 
 * It is a community of matrices and solvers; hereafter we refer it as the Community. 
 * For a problem (the Community contains problem descriptions), there are a bunch of 
 * solvers (the Community contains algorithm descriptions and source codes) for this 
 * problem. A solver could be applied to different problems (the Community compares 
 * performance and convergence of the solver for various of problems)
 *
 * 
 * \section sec_multigrid Our goal and strategy 
 *
 * We can think of the Community in a Multigrid + Capitalism point of view.
 *
 * \subsection ss_fine Fine grid (free market stage)
 *
 * (1) Collect problems and solvers. Allow duplications, for example same solution 
 * algorithm, but different implementation or different programming languages. 
 *
 * (2) Do not enforce much regulation. Allow the market to be FREE. 
 * 
 * (3) Keep all the record: problem description, solver code, test results, etc. 
 * Do not throw anything away.
 *
 * (4) We are currently at the initial stage. We try to find a minimal set of 
 * standard or rules. And then we let the market to evolve freely.
 *
 *
 * \subsection ss_coarse Coarse grid (state capitalism stage)
 *
 * (1)	As the market evolves, we might find at certain time when the market is 
 * out-of-control. This basically means the "free market" is very successful. And 
 * now we need to give more restrict standard or regulation.
 *
 * (2) We still keep everything before standardized; just a new branch containing 
 * standardized solvers. Different versions marked with different colors for 
 * convenience.
 * 
 * (3) At certain stage, we might write professional level codes for the solvers 
 * and form a package. We are now working toward this goal. 
 *
 * (4)	Host an annual user-developer meeting and set up a solver competition. Start 
 * a journal for computational softwares. 
 *
 *
 * 
 * \section sec_guide General guidelines 
 *
 * 1. For the moment, this project is still restricted to our own group only. 
 * It does not have to be user-friendly to people outside of the group. But at 
 * least, for people in the group, the problems should be easy to understand 
 * and the solvers should be easy to compile and run. 
 *
 * 2. Write as much comments as possible for easy maintenance and usage. 
 * See \ref page_comment for more detailed instructions on how to add comment 
 * using Doxygen. 
 *
 * 3. We mainly interested in the problems with limited information (mostly 
 * algebraic) only (maybe with some mesh information and/or local stiffness 
 * matrices); See \ref page_steps for solving steps.  
 * 
 * 4. We are mathematicians, not software engineers. We are going to do it as
 * mathematicians without worrying about the implementation quality too much. 
 * We don't need to build an industry quality software package. It is for 
 * academic (research and education) purpose. Version control is still very 
 * important for our collaborative work; see \ref page_comment for the HG 
 * settings.
 *
 *
 */
 
 
/** 
 * \page page_cvs Obtaining FASP 
 *
 * In order to coordinate our collaboration, we use HG as our version control tool. 
 * If you are willing to contribute to the project, you can send 
 * Chensong Zhang multigrid@me.com an email.
 *
 * HG (Mecurial) is easy to use and it is available at all platforms. For Linux and 
 * OS X users, set up the following environment variables.
 *
 * Windows version is supposed to be fool-proof and we don't include an instruction 
 * for Windows version here. 
 *
 * We use SSH for security. But we don't use passwd access, instead everybody should 
 * send me you public SSH authorization code. If you already have one, it should locate 
 * in your home directory at
 *
 * ~/.ssh/id_rsa.pub
 *
 * If you don't have one in .ssh, you can create one by 
 * 
 * $ ssh-keygen -t rsa
 *
 * The program will ask you a few questions, just hit on RETURN for every question. 
 *
 * If you don't even have .ssh direction, you need to create on 
 * 
 * $ mkdir ~/.ssh
 * 
 * When you find this ***.pub file, please send it to me. So I can give you access to 
 * the HG server without passwd. Now we are all set to go. 
 *
 * Now just do 
 * 
 * $ hg clone ssh://psucuda\@ganymede.screms.math.psu.edu/fasp 
 * 
 * will give you the current version of the package. 
 * 
 */
 

/**
 * \page page_readme Installation and Building
 *
 * \section build Installation -- Building The Library:
 * 
 * This is a simple instruction on building and testing. For more details, please refer to 
 * the README files and the user's manual. 
 *
 * To compile, you need a Fortran and a C compiler. By default, we use gfortan and gcc, 
 * respectively; see core/Makefile and core/make.in for details.
 *
 * First, if you are building FASP for the first time, do setup and generate library and 
 * output directories:
 * >> ./bin/setup.sh
 *
 * which setup lib/ and output/ directories if they are missing. 
 *
 * And then you can make the FASP library in the "core/" directory: 
 *
 * >> make
 *
 * which makes the libfasp.a static library.
 *
 * If you wish to see the detailed usage of "make" or need any help, please 
 *
 * >> make help
 * 
 *
 * \section test Test:
 *
 * By running the executable test.ex, we can get numerical results for different test examples. 
 *
 * >> ./test.ex
 *
 * is the terminal command for reading ini/input.dat and run test problems. There is also a regression test tool:
 *
 * >> ./regression.ex 
 *
 * \section input Input Parameters: 
 *
 * And test.ex reads parameters from ini/input.dat, where you can choose:
 *
 * Solver type ( 0--5 ):  
 *		- 0 AMG as itertive solver;
 *		- 1 CG method;
 *		- 2 BiCGstab method;
 *		- 3 MinRes method;
 *		- 4 GMRES method;
 *		- 5 Variable-restart GMRES method.
 *
 * Precond type ( 0--3 ): 
 *		- 0 No preconditioners;
 *		- 1 Diagonal preconditioner;
 *		- 2 AMG preconditioner: classical AMG or smooth aggregation AMG;
 *		- 3 ILU preconditioner: ILUk and ILUt.
 *
 */ 

 
/**
 * \page page_steps Solving Test Problems
 *
 * To solve a problem (some times we only have the coefficent matrix and the 
 * right-hand side), we can do the following five steps: 
 *
 * \section step1 Step 1: obtain a problem
 *
 * Get a discrete PDE from somewhere: Usually reading from a disk file or passing 
 * by other programs.
 *  
 * \section step2 Step 2: check matrix property
 *
 * Run a sequence of tests to see: whether this matrix is symmetric, positive 
 * definite, sparse, etc. 
 *   
 * \section step3 Step 3: select a solver
 *
 * Using some artificial intelligence, we pick a solver which is suitable for 
 * the problem we are dealing with. This is the brain of this whole project. 
 *   
 * \section step4 Step 4: solve the system
 *
 * Once the solver has been chosen, we solve the system. 
 *   
 * \section step5 Step 5: post processing
 *   
 * Graphical or text output.
 *
 */
 
 
/**
 * \page page_algorithms Algorithms
 *
 * Here we list all the algorithms we implemented and tested. 
 *
 * 1. Sparse BLAS and other generic operations
 *
 * 2. Finite element discretizations
 *
 * (i) P1, P2, P3, P4 finite element basis functions
 *
 * (ii) General assembling framework for unstructured grids
 *
 * 3. Iterative solvers
 *
 * (i) CG
 *
 * (ii) BiCGStab
 *
 * (iii) GMRES
 *
 * (iv) MINRES
 * 
 * 4. Preconditioners
 *
 * (i) Diagonal preconditioner
 *
 * (ii) Classical AMG preconditioners: Ruge-Stuben | Energy-Minimization
 *
 * (iii) Smoothed aggregation AMG preconditioners
 *
 * (iv) Incomplete decomposition preconditioners: ILUt | ILUk
 *
 */
 
 
/**
 * \page developers Developers
 *
 * Names are listed in alphabetical order:
 *
 *      -Brannick, James (Penn State University, USA)
 *
 *		-Feng, Chunsheng (Xiangtan University, China)
 *
 *      -Hu, Xiaozhe (Peking University, China)
 * 
 *      -Shu, Shi (Xiangtan University, China)
 *
 *      -Xu, Jinchao (Penn State University, USA)
 *
 *      -Zhang, Chensong (Chinese Academy of Sciences, China)
 *
 *      -Zhou, Zhiyang (Xiangtan University, China)
 *
 *      -Zikatanov, Ludmil (Penn State Univeristy, USA)
 *
 * Former developers:
 * 
 *      -Huang, Xuehai (Shanghai Jiaotong University, China)
 *
 *      -Zhang, Shiquan (Sichuan University, China)
 *
 *      -Zhang, Shuo (Chinese Academy of Sciences, China)
 *
 */
   
/**
 * \page page_comment Documentation with Doxygen
 *
 * We use Doxygen as our automatically documentation generator which will make our 
 * future maintainance minimized. You can obtain the software (Windows, Linux and 
 * OS X) as well as its manual on the official website 
 * 
 * http://www.doxygen.org
 *
 * For an oridinary user, Doxygen is completely trivial to use. We only need to use 
 * some special marker in the usual comment as we put in c-files. 
 * See test.c and matvecio.c for examples.  
 *
 */
 
/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
 
