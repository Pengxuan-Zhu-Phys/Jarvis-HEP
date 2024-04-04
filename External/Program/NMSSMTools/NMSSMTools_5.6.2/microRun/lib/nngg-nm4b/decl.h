* decl.h
* these declarations are included "everywhere"
* this file is part of FormCalc
* last modified 2 Sep 12 th


#ifndef DECL_H
#define DECL_H

* declarations for the whole file (e.g. preprocessor defs)

#include "types.h"
#include "user.h"

#define MW2 MW**2
#define Finite 1D0

#else

* declarations for every subroutine

#include "user.h"
#include "util.h"
#include "looptools.h"
#include "renconst.h"

#endif

