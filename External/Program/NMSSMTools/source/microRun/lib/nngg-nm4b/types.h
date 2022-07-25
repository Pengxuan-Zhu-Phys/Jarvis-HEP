* types.h
* real-based type declarations
* this file is part of FormCalc
* last modified 17 Jul 12 th


#ifndef TYPES_H
#define TYPES_H

#define RealType double precision
#define ComplexType double complex
#define Re DBLE
#define Im DIMAG
#define Conjugate DCONJG
#define ToComplex DCMPLX

#define Sq(c) Re((c)*Conjugate(c))

#define marker double precision

#endif

