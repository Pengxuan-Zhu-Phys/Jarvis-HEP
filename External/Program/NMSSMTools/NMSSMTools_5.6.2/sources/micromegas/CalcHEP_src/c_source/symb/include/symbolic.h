#ifndef __SYMBOLIC_
#define __SYMBOLIC_

extern void  calcallproc(void);

extern int sqDiagList( int **sd, int nCore);

extern void  calcWithFork(int np, int * diag,int fi);

extern void updateMenuQ(void);

extern int writeVertexCode(char*pathToModels,int Model,int N,char**field,char*label);

#endif
