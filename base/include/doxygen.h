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
 * > For the moment, FASP is under alpha testing. If you wish to obtain a current version of
 * > FASP or you have any questions, feel free to contact us at faspdev@gmail.com.
 * 
 * This software distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
 * See the GNU Lesser General Public License for more details.
 *
 */
 
/**
 * \page download Download
 * \brief Download page for the FASP package
 *
 * ## How to obtain FASP
 *
 * The most updated release version of FASP can be downloaded from
 *
 * > http://fasp.sourceforge.net/download/faspsolver.zip
 *
 * We use HG (Mecurial) as our version control tool. HG is easy to use and it is available at
 * all platforms. For people who is interested in the developer version, you can obtain the
 * FASP package with hg: 
 * 
 * > $ hg clone https://faspusers@bitbucket.org/fasp/faspsolver 
 * 
 * will give you the current version of the package. 
 *
 * ## Installation and Building
 *
 * This is a simple instruction on building and testing. For more details, please refer to 
 * the README files and the short
 * <a href="http://fasp.sourceforge.net/download/userguide.pdf">\bf User's Guide</a> in "doc/". 
 *
 * To compile, you need a Fortran and a C compiler.  First, you can type:
 *
 * > $ make config
 * 
 * which will config the environment automatically. And, then, you can need to type:
 *
 * > $ make install
 *
 * which will make the FASP shared static library and install to PREFIX/.
 *
 * If you wish to see the detailed usage of "make" or need any help, please type:
 *
 * > $ make help
 *
 * Tutorial examples can be found in "tutorial/".
 *
 */ 


/**
 * \page developers Developers
 * \brief FASP developers and contributors
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
 * - Li, Zheng (Xiangtan University, China)
 *
 * - Shu, Shi (Xiangtan University, China)
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
 * - Qiao, Changhe (Penn State University, USA)
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
 * \page doxygen_comment Doxygen
 * \brief Documentation Auto-Generation
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
 
