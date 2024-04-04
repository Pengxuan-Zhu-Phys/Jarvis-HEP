
/*   Terms library 
 *   by A. Semenov
 *   1995 - 1998
 */
 
#ifndef TERMS_H
#define TERMS_H

#define TERMTYPE(t) (t>>56)

#ifdef MTHREAD
#include <pthread.h>
extern pthread_key_t TermsKey;
#endif

typedef signed long int Atom, Integer, Float,
		Atomic, Functor, Term, Compound, List, Label;
		
#include "atoms.h"

union LISTpoint 
	{
	struct {
		Term value;
		List  nextl;
		union LISTpoint *nextp;
		} d;
	struct 	{
		List  l;
		union LISTpoint *next;
		} p;
	};

	
typedef struct
	{
		int h_off, h_sz;
		List *ori_list;
		
		List prev_l;
		List cur_l;
		union LISTpoint *prev_p;
		union LISTpoint *cur_p;
	} ListIterator;
		
		
	
/*extern union LISTpoint *LISTS_next_free;*/
extern union LISTpoint *LISTS_buffers[400096];
extern int LISTS_blocks, LISTS_allocs, LISTS_free;


#define is_empty_list(l) (l==0)
#define LIST_pointer(l) (LISTS_buffers[(l&0xfffff0000)/0x10000]+(l&0xffff))
#define ListFirst(l) (LIST_pointer(l)->d.value)
#define ListTail(l) (LIST_pointer(l)->d.nextl)
#define ListConcat(l1,l2) (LIST_pointer(l1)->d.nextl=l2, \
	 LIST_pointer(l1)->d.nextp=LIST_pointer(l2))

#ifdef USE_INLINE

#define __inline inline static
#include "lists.ci"
#undef __inline

#else

List NewList(void);
/*int  is_empty_list(List);
Term ListFirst(List);
List ListTail(List);*/
List AppendFirst(List,Term);
List AppendLast(List,Term);
int  ListLength(List l);
void ChangeList(List l, Term nv);
int 	is_list(Term a);

#endif

union COMPOpoint 
	{
	struct {
		Atomic a[3];
		Atomic next_a;
		union COMPOpoint *next_p;
		} data;
	struct 	{
		union COMPOpoint *next;
		Atomic this_a;
		} p;
	};
	


/*extern union COMPOpoint *COMPO_next_free;*/
extern union COMPOpoint *COMPO_buffers[400096];
extern int COMPO_blocks, COMPO_tall, COMPO_tfree;

extern long int LABEL_count;

extern char *ATOM_buffers[256];
extern int  ATOM_buffill[256], ATOM_count, ATOM_stcount, ATOM_tscount;


#ifdef USE_INLINE

#define __inline inline static
#include "data.ci"
#undef __inline

#else

Integer	NewInteger(long val);
long 	IntegerValue(Integer i);

Functor NewFunctor(Atom name, int arity);
/*Atom 	FunctorName(Functor f);
int 	FunctorArity(Functor f);*/


#endif

#define CompoundFunctor(c) ((COMPO_buffers[(c&0xfffff0000)/0x10000]+(c&0xffff))->data.a[0])
#define CompoundArg1(c) ((COMPO_buffers[(c&0xfffff0000)/0x10000]+(c&0xffff))->data.a[1])
#define CompoundArg2(c) ((COMPO_buffers[(c&0xfffff0000)/0x10000]+(c&0xffff))->data.a[2])

#define FunctorName(f) ((f & 0xffffff) + (1L<<56))
#define FunctorArity(f) ((int)(((f>>48)&0xff)-1))

#define CompoundName(c) FunctorName(CompoundFunctor(c))
#define CompoundArity(c) FunctorArity(CompoundFunctor(c))

#ifdef MTHOPT

List CopyList2(List, List *, int);
void FreeAtomic2(Atomic, int);
Term CopyTerm2(Term, int);

#define MTHOPTDECL int *mindp=pthread_getspecific(TermsKey), mind = *mindp;

#else

#define CopyList2(a,b,c) CopyList(a,b)
#define FreeAtomic2(a,b) FreeAtomic(a)
#define CopyTerm2(a,b) CopyTerm(a)

#define MTHOPTDECL

#endif


extern char *OutputDirectory;


	/*   File data.c  */

int 	is_atomic(Term), is_compound(Term), is_functor(Atomic),
	is_atom(Atomic), is_float(Atomic), is_integer(Atomic);

void 	InitAtoms(void);

Atom 	NewAtom(char *s, int len);
char 	*AtomValue(Atom);
void    SetAtomProperty(Atom a, Atom type, Term value);
Term    GetAtomProperty(Atom a, Atom type);
void    RemoveAtomProperty(Atom a, Atom type);
Term    GetProperties(Term,Term);
Term 	SetProperty(Term,Term);
List    AtomPropertiesList(Atom);


typedef struct { double r, i;} cmplx;
Float   NewComplex(cmplx);
Float 	NewFloat(double d);
double 	FloatValue(Float);
cmplx  ComplexValue(Float);
void   SetFloat(Float, double, double);
double ImaginaryValue(Float);

Compound NewCompound(Functor f);
/*Functor CompoundFunctor(Compound c);
Term    CompoundArg1(Compound c);
Term    CompoundArg2(Compound c);*/
Term    CompoundArgN(Compound c, int arg);
Term    ConsumeCompoundArg(Compound c, int arg);
void    SetCompoundArg(Compound c, int arg, Term t);
void    SetCompoundName(Compound c, Atom n);
Compound MakeCompound(Atom name, int Arity);
Compound MakeCompound1(Atom name, Term arg1);
Compound MakeCompound2(Atom name, Term arg1, Term arg2);
Compound MakeCompound3(Atom name, Term arg1, Term arg2, Term arg3);
/*Atom	CompoundName(Compound c);
int     CompoundArity(Compound c);
*/
		
Atomic  NewLabel(void);
long int 	LabelValue(Atomic l);

void 	FreeAtomic(Atomic);
Term	CopyTerm(Term);
int 	EqualTerms(Term,Term);
char 	AtomicType(Atomic);
int 	is_atomic(Atomic a);
int 	is_float(Atomic a);
int 	is_integer(Atomic a);
int 	is_atom(Atomic a);
int 	is_functor(Atomic a);
int 	is_compound(Term a);
int     is_label(Atomic a);

void 	AtomStatistics(void);
void 	AtomStat1(void);
long    TermMemory(void);

	/* File lists.c */

List SortedList(List li, int (*cmp)(Term,Term));
List ListSplit(List, List, List *);
List CutFromList(List l, List c);
void FreeList(List l);
void RemoveList(List l);
Term ListNth(List l, int arg);
List ListNthList(List l, int arg);
void InsertList(List l, Term ins);
void DumpList(List l);
void ListStatistics(void);
void ListStat1(void);
long ListMemory(void);
List ConcatList(List, List);
List CopyList(List, List *);
int  ListMember(List, Term);
List IntersectList(List, List);
List MakeList1(Term);
List MakeList2(Term,Term);
List MakeList3(Term,Term,Term);
List MakeList4(Term,Term,Term,Term);


	/*  File lexic.c  */

int  	SetInputFile(char *name);
int     SetInputFileM(char *line);
void 	CloseInputFile(void);
Atomic 	ReadAtomic(void);
void  	UnreadAtomic(Atomic);
char 	*CurrentInputFile(void);
int   	CurrentInputLine(void);
void 	WritePrompt(int );
void 	SkipInputLine(void);
int 	IsTermInput(void);

	/* File write.c  */

int  DisplayTerm(Term);
int  WriteTerm(Term);
int  fWriteTerm(FILE *, Term);
int  sWriteTerm(char *, Term);
int  fDisplayTerm(FILE *, Term);
void DumpTerm(Term);

extern int AlwaysBracets, WideWriting, NoQuotes, fortr_digi;


	/* File itrio.c */

int  itrSetOut(char *);
void itrOut(Term);
void itrCloseOut(void);
int  itrSetIn(char *);
Term itrIn(void);
void itrCloseIn(void);

	/* File  ops.c  */

void SetOperator(Atom name, Atom assoc, int prior);
Compound  GetOperator(Atom name, Atom type, Atom *assoc, int *prior);
Compound  CurrentOperator(int no);

#define OP_NONE		0
/*
#define OP_PREFIX	1
#define OP_INFIX	2
#define OP_POSTFIX	3
*/
	/*	File  parse.c */

Term ReadTerm(void);
Term ReadTermM(char *line);
Term ReadDisplTerm(void);

extern int ExitOnError;


#endif

