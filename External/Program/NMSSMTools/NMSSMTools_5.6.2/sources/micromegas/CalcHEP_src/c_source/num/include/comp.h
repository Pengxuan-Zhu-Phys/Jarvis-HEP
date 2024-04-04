#ifndef __COMPOSITES__
#define __COMPOSITES__

#include "file_scr.h"
#include "../../../include/VandP_size.h"
extern int fillCompositeArray(void);
extern int rdrcomp_(FILE *);
extern int wrtcomp_(FILE *);


extern table compTab;
extern int nComps,nCompParts[60];
extern char compName[60][P_NAME_SIZE], compParts[60][60][P_NAME_SIZE];
 

#endif
