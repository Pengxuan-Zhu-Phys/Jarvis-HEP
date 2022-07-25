#define _LONG_          to use long double type in numerical calculations
//#define _QUADICC_     quadrouple precision for icc compiler. It has to be accompany with compiler option '-Qoption,cpp,--extended_float_type'  
//#define _QUADGCC_     quadrouple precision for gcc   https://gcc.gnu.org/onlinedocs/gcc-4.7.2/libquadmath.pdf

#ifndef __NTYPE__
#define __NTYPE__

#include<math.h>
#include<complex.h>

#ifdef _LONG_  

#define REAL    long double
#define COMPLEX long double complex

#define nTypePre 
#define nTypePost l


#else 

#ifdef _QUADICC_

#define REAL _Quad
#define COMPLEX  complex _Quad
//http://software.intel.com/en-us/forums/topic/289725


#define nTypePre __
#define nTypePost q

#else
#ifdef _QUADGCC_

#include  <quadmath.h>

#define REAL __float128
#define COMPLEX  __complex128

// https://gcc.gnu.org/onlinedocs/gcc-4.7.2/libquadmath.pdf

#define nTypePre 
#define nTypePost q

#else 

#define REAL double
#define COMPLEX double complex

#define nTypePre
#define nTypePost

#endif
#endif
#endif

#define nTypeSSS_(S1,S2,S3) S1##S2##S3
#define nTypeSSS(S1,S2,S3) nTypeSSS_(S1,S2,S3)


#define Fabs    nTypeSSS(nTypePre, fabs,  nTypePost)    
#define Sqrt    nTypeSSS(nTypePre, sqrt,  nTypePost)    
#define Pow     nTypeSSS(nTypePre, pow,   nTypePost)     
#define Exp     nTypeSSS(nTypePre, exp,   nTypePost)     
#define Log     nTypeSSS(nTypePre, log,   nTypePost)     
#define log10   nTypeSSS(nTypePre, log10, nTypePost)   
#define Sin     nTypeSSS(nTypePre, sin,   nTypePost)     
#define Cos     nTypeSSS(nTypePre, cos,   nTypePost)     
#define Tan     nTypeSSS(nTypePre, tan,   nTypePost)     
#define Asin    nTypeSSS(nTypePre, asin,  nTypePost)    
#define Acos    nTypeSSS(nTypePre, acos,  nTypePost)    
#define Atan    nTypeSSS(nTypePre, atan,  nTypePost)    
#define Atan2   nTypeSSS(nTypePre, atan2, nTypePost)   
#define Sinh    nTypeSSS(nTypePre, sinh,  nTypePost)    
#define Cosh    nTypeSSS(nTypePre, cosh,  nTypePost)    
#define Tanh    nTypeSSS(nTypePre, tanh,  nTypePost)    
#define Asinh   nTypeSSS(nTypePre, asinh, nTypePost)   
#define Acosh   nTypeSSS(nTypePre, acosh, nTypePost)   
#define Atanh   nTypeSSS(nTypePre, atanh, nTypePost)   
#define Fabs    nTypeSSS(nTypePre, fabs,  nTypePost)

#define Creal   nTypeSSS(nTypePre, creal, nTypePost)   
#define Cimag   nTypeSSS(nTypePre, cimag, nTypePost)   
#define Carg    nTypeSSS(nTypePre, carg,  nTypePost)    
#define Cabs    nTypeSSS(nTypePre, cabs,  nTypePost)    
#define Conj    nTypeSSS(nTypePre, conj,  nTypePost)    
#define Cproj   nTypeSSS(nTypePre, cproj, nTypePost)   

#define Csqrt   nTypeSSS(nTypePre, csqrt, nTypePost)   
#define Cpow    nTypeSSS(nTypePre, cpow,  nTypePost)    
#define Clog    nTypeSSS(nTypePre, clog,  nTypePost)    
#define Clog10  nTypeSSS(nTypePre, clog10,nTypePost)  
#define Cexp    nTypeSSS(nTypePre, cexp,  nTypePost)    
#define Csin    nTypeSSS(nTypePre, csin,  nTypePost)    
#define Ccos    nTypeSSS(nTypePre, ccos,  nTypePost)    
#define Ctan    nTypeSSS(nTypePre, ctan,  nTypePost)    
#define Casin   nTypeSSS(nTypePre, casin, nTypePost)   
#define Cacos   nTypeSSS(nTypePre, cacos, nTypePost)   
#define Catan   nTypeSSS(nTypePre, catan, nTypePost)   
#define Csinh   nTypeSSS(nTypePre, csinh, nTypePost)   
#define Ccosh   nTypeSSS(nTypePre, ccosh, nTypePost)   
#define Ctanh   nTypeSSS(nTypePre, ctanh, nTypePost)   
#define Casinh  nTypeSSS(nTypePre, casinh,nTypePost)  
#define Cacosh  nTypeSSS(nTypePre, cacosh,nTypePost)  
#define Catanh  nTypeSSS(nTypePre, catanh,nTypePost)


#endif
