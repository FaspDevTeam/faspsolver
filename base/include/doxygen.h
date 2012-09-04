/*! \file doxygen.h
 *  \brief Main page for Doygen documentation
 */
 
/** \mainpage FASP Solver Package
 *
 *
 * ## Introduction
 *
 * Over the last few decades, researchers have expended significant effort on developing 
 * efficient iterative methods for solving discretized partial differential equations 
 * (PDEs). Though these efforts have yielded many mathematically optimal solvers such as 
 * the multigrid method, the unfortunate reality is that multigrid methods have not been 
 * much used in practical applications. This marked gap between theory and practice is 
 * mainly due to the fragility of traditional multigrid (MG) methodology and the complexity
 * of its implementation. We aim to develop techniques and the corresponding software that 
 * will narrow this gap, specifically by developing mathematically optimal solvers that are 
 * robust and easy to use in practice. 
 * 
 * We believe that there is no one-size-for-all solution method for discrete linear systems 
 * from different applications. And, efficient iterative solvers can be constructed by taking 
 * the properties of PDEs and discretizations into account. In this project, we plan to 
 * construct a pool of discrete problems arising from partial differential equations (PDEs) 
 * or PDE systems and efficient linear solvers for these problems. We mainly utilize the 
 * methodology of Auxiliary Space Preconditioning (ASP) to construct efficient linear solvers. 
 * Due to this reason, this software package is called Fast Auxiliary Space Preconditioning
 * or FASP for short.  
 *
 * FASP contains the kernel part and several applications (ranging from fluid dynamics to 
 * reservoir simulation). The kernel part is open-source and licensed under GNU Lesser General
 * Public License or LGPL version 3.0 or later. Some of the applications contain contributions
 * from and owned partially by other parties.
 *
 *
 * ## How to obtain FASP
 *
 * In order to coordinate our collaboration, we use HG as our version control tool. 
 * HG (Mecurial) is easy to use and it is available at all platforms. For Linux and 
 * OS X users, set up the following environment variables.
 *
 * To obtain the FASP package, just do 
 * 
 * $ hg clone https://faspusers@bitbucket.org/fasp/faspsolver 
 * 
 * will give you the current version of the package. 
 *
 * > For the moment, FASP is under alpha testing. If you wish to obtain a current version of
 * > FASP or you have any questions, feel free to contact us at faspdev@gmail.com.
 * 
 *
 * ## Installation and Building
 *
 * This is a simple instruction on building and testing. For more details, please refer to 
 * the README files and the user's guide in "doc/". 
 *
 * To compile, you need a Fortran and a C compiler. By default, we use gfortan and gcc, 
 * respectively; see base/Makefile and base/make.in for details.
 *
 * First, you can make the FASP library in the "base/" directory: 
 *
 * $ make
 *
 * which makes the libfasp.a static library.
 *
 * If you wish to see the detailed usage of "make" or need any help, please 
 *
 * $ make help
 *
 */ 


/**
 * \page developers Developers
 *
 * Project Leader: 
 *
 * - Xu, Jinchao (Penn State University, USA)
 *
 * Developers (in alphabetic order):
 *
 * - Feng, Chunsheng (Xiangtan University, China)
 *
 * - Hu, Xiaozhe (Penn State University, USA)
 * 
 * - Zhang, Chensong (Chinese Academy of Sciences, China)
 *
 * - Zikatanov, Ludmil (Penn State Univeristy, USA)
 *
 * With contributions from (in alphabetic order):
 *
 * - Brannick, James (Penn State University, USA)
 * 
 * - Cao, Fei (Penn State University, USA)
 *
 * - Huang, Feiteng (Sichuang University, China)
 *
 * - Huang, Xuehai (Shanghai Jiaotong University, China)
 *
 * - Li, Zheng (Xiangtan University, China)
 *
 * - Qiao, Changhe (Penn State University, USA)
 *
 * - Shu, Shi (Xiangtan University, China)
 *
 * - Wang, Ziteng (Peking University, China)
 *
 * - Yang, Kai (Penn State University, USA)
 *      
 * - Yue, Xiaoqiang (Xiangtan University, China)
 *
 * - Zhang, Shiquan (Sichuan University, China)
 *
 * - Zhang, Shuo (Chinese Academy of Sciences, China)
 *
 * - Zhang, Weifeng (Kunming University of Science and Technology, China)
 *
 * - Zhou, Zhiyang (Xiangtan University, China)
 *
 * Project Coordinator:
 * 
 * - Zhang, Chensong (Chinese Academy of Sciences, China)
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
 *
 */
 
/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
 
