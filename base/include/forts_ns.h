#ifdef __cplusplus 
extern "C" {void cut0_(INT *n,INT *ia,INT *ja,REAL *a,INT *iaw,INT *jaw,INT *jblk,INT *iblk,INT *nblk,INT *lwork1,INT *lwork2,INT *lwork3,INT *msize);}
extern "C" {void chsize_(REAL *a,REAL *b,REAL *tol,INT *imin);}
extern "C" {void shift_(INT *nxadj,INT *nadj,INT *n);}
extern "C" {void dfs_(INT *n,INT *ia,INT *ja,INT *nblk,INT *iblk,INT *jblk,INT *lowlink,INT *iedge,INT *numb);}
extern "C" {void permat_(INT *iord,INT *ia,INT *ja,REAL *an,INT *n,INT *m,INT *iat,INT *jat,REAL *ant);}
extern "C" {void pervec_(INT *iord,REAL *u1,REAL *u2,INT *n);}
extern "C" {void perback_(INT *iord,REAL *u1,REAL *u2,INT *n);}
extern "C" {void perm0_(INT *iord,INT *ia,INT *ja,REAL *an,INT *n,INT *m,INT *iat,INT *jat,REAL *ant);}
extern "C" {void icopyv_(INT *iu,INT *iv,INT *n);}
extern "C" {void mxfrm2_(INT *n,INT *ia,INT *ja,INT *nblk,INT *iblock,INT *jblock,INT *mask,INT *maxa,INT *memt,INT *maxbs);}
extern "C" {void sky2ns_(INT *n,INT *ia,INT *ja,REAL *a,INT *nblk,INT *iblock,INT *jblock,INT *mask,INT *maxa,REAL *au,REAL *al);}
extern "C" {void fbgs2ns_(INT *n,INT *ia,INT *ja,REAL *a,REAL *x,REAL *b,INT *nblk,INT *iblock,INT *jblock,INT *mask,INT *maxa,REAL *au,REAL *al,REAL *rhsloc,INT *memt);}
extern "C" {void bbgs2ns_(INT *n,INT *ia,INT *ja,REAL *a,REAL *x,REAL *b,INT *nblk,INT *iblock,INT *jblock,INT *mask,INT *maxa,REAL *au,REAL *al,REAL *rhsloc,INT *memt);}
extern "C" {void doluns_(REAL *au,REAL *al,INT *maxa,INT *nn);}
extern "C" {void sluns_(REAL *au,REAL *al,REAL *v,INT *maxa,INT *nn);}
extern "C" {void dolu_(REAL *a,INT *maxa,INT *nn);}
extern "C" {void slvlu_(REAL *a,REAL *v,INT *maxa,INT *nn);}
extern "C" {void ijacrs_(INT *ln,INT *ia,INT *ja,REAL *a,INT *n,INT *nnz,INT *ir,INT *ic,REAL *aij);}
extern "C" {void sympat_(INT *ln,INT *ia,INT *ja,INT *n,INT *ir,INT *ic,REAL *aij);}
extern "C" {void levels_(INT *inroot,INT *ia,INT *ja,INT *mask,INT *nlvl,INT *iblock,INT *jblock,INT *maxlev);}
#else
extern void cut0_(INT *n,INT *ia,INT *ja,REAL *a,INT *iaw,INT *jaw,INT *jblk,INT *iblk,INT *nblk,INT *lwork1,INT *lwork2,INT *lwork3,INT *msize);
extern void chsize_(REAL *a,REAL *b,REAL *tol,INT *imin);
extern void shift_(INT *nxadj,INT *nadj,INT *n);
extern void dfs_(INT *n,INT *ia,INT *ja,INT *nblk,INT *iblk,INT *jblk,INT *lowlink,INT *iedge,INT *numb);
extern void permat_(INT *iord,INT *ia,INT *ja,REAL *an,INT *n,INT *m,INT *iat,INT *jat,REAL *ant);
extern void pervec_(INT *iord,REAL *u1,REAL *u2,INT *n);
extern void perback_(INT *iord,REAL *u1,REAL *u2,INT *n);
extern void perm0_(INT *iord,INT *ia,INT *ja,REAL *an,INT *n,INT *m,INT *iat,INT *jat,REAL *ant);
extern void icopyv_(INT *iu,INT *iv,INT *n);
extern void mxfrm2_(INT *n,INT *ia,INT *ja,INT *nblk,INT *iblock,INT *jblock,INT *mask,INT *maxa,INT *memt,INT *maxbs);
extern void sky2ns_(INT *n,INT *ia,INT *ja,REAL *a,INT *nblk,INT *iblock,INT *jblock,INT *mask,INT *maxa,REAL *au,REAL *al);
extern void fbgs2ns_(INT *n,INT *ia,INT *ja,REAL *a,REAL *x,REAL *b,INT *nblk,INT *iblock,INT *jblock,INT *mask,INT *maxa,REAL *au,REAL *al,REAL *rhsloc,INT *memt);
extern void bbgs2ns_(INT *n,INT *ia,INT *ja,REAL *a,REAL *x,REAL *b,INT *nblk,INT *iblock,INT *jblock,INT *mask,INT *maxa,REAL *au,REAL *al,REAL *rhsloc,INT *memt);
extern void doluns_(REAL *au,REAL *al,INT *maxa,INT *nn);
extern void sluns_(REAL *au,REAL *al,REAL *v,INT *maxa,INT *nn);
extern void dolu_(REAL *a,INT *maxa,INT *nn);
extern void slvlu_(REAL *a,REAL *v,INT *maxa,INT *nn);
extern void ijacrs_(INT *ln,INT *ia,INT *ja,REAL *a,INT *n,INT *nnz,INT *ir,INT *ic,REAL *aij);
extern void sympat_(INT *ln,INT *ia,INT *ja,INT *n,INT *ir,INT *ic,REAL *aij);
extern void levels_(INT *inroot,INT *ia,INT *ja,INT *mask,INT *nlvl,INT *iblock,INT *jblock,INT *maxlev);
#endif