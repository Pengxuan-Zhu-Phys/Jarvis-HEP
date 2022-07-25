#include <string.h>
#include <stdlib.h>
#include "lanhep.h"

#define TBSIZE 65500

#ifdef MTHREAD
#include <pthread.h>

extern pthread_key_t TermsKey;
static pthread_mutex_t mtx  = PTHREAD_MUTEX_INITIALIZER;
#endif


List CopyList2(List, List *, int);
void FreeAtomic2(Atomic a, int m);
/*#define DATA_DBG*/

/*********

  1 - atoms
  2-30 functors
  31- floats
  32-47 - compounds
  48-63 - lists
  64 - * labels
  
  ---- new notation
  
  1 - atoms
  2 - functors
  3 - floats
  4 - compounds
  5 - lists
  6 - labels
***********/

union fpoint
	{
	struct { double d,i;} v;
	struct 	{
		union fpoint *next;
		Float f;
		} p;
	};



static union fpoint *f_next_free=NULL;
static union fpoint *fbuffers[256];
static int f_blocks=0, f_dtf=0, f_ftd=0;

Float NewFloat(double d)
	{
	union fpoint *p;
	Float f;
	if(f_next_free==NULL)
		{
		int i;
		if(f_blocks==256)
			{  puts("Internal error (too much floats or memory lack)."); exit(0); }
		p=fbuffers[f_blocks]=(union fpoint*)malloc(sizeof(union fpoint)*1000);
		if(p==NULL)
			{  puts("Internal error (lack space for float)."); exit(0); }
		for(i=0;i<1000;i++)
			{
			p[i].p.next=f_next_free;
			f_next_free=p+i;
			p[i].p.f=(3L<<56)+f_blocks*0x10000+i;
			}
		f_blocks++;
		}
	p=f_next_free;
	f_next_free=f_next_free->p.next;
	f=p->p.f;
	p->v.d=d;
	p->v.i=0.0;
	f_dtf++;
	return f;
	}


Float NewComplex(cmplx cx)
{
  Float f=NewFloat(0.0);
  SetFloat(f,cx.r,cx.i);
  return f;
}

double FloatValue(Float f)
	{
	union fpoint *p;
	int bno,bpos;
	bno=f&0xff0000;
	bno/=0x10000;
	bpos=f&0xffff;
	p=fbuffers[bno]+bpos;
	return p->v.d;
	}

cmplx ComplexValue(Float f)
{
  cmplx ret;
  ret.r=FloatValue(f);
  ret.i=ImaginaryValue(f);
  return ret;
}

void SetFloat(Float f, double r, double i)
{
  union fpoint *p;
	int bno,bpos;
	bno=f&0xff0000;
	bno/=0x10000;
	bpos=f&0xffff;
	p=fbuffers[bno]+bpos;
	p->v.d=r;
	p->v.i=i;
}

double ImaginaryValue(Float f)
{
  union fpoint *p;
	int bno,bpos;
	bno=f&0xff0000;
	bno/=0x10000;
	bpos=f&0xffff;
	p=fbuffers[bno]+bpos;
	return p->v.i;
}

#ifdef MTHREAD
union COMPOpoint *COMPO_next_free[17]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
#else
union COMPOpoint *COMPO_next_free=NULL;
#endif
union COMPOpoint *COMPO_buffers[400096];
int COMPO_blocks=0, COMPO_tall=0, COMPO_tfree=0;
long int LABEL_count = 0;


/*
#ifndef USE_INLINE

#define __inline
#include "data.ci"
#undef __inline

#endif
*/

#define __inline

__inline Integer NewInteger(long l)
	{
	Integer i;
	i=l-2223372036854775807L;
	if(i>=0)
		{  puts("Internal error (integer out of range)."); 
		/*abort();*/
		   return NewInteger(0); }
	return i;
	}

__inline long IntegerValue(Integer i)
	{
	return i+2223372036854775807L;
	}

/*static int newfu[30];*/


__inline Functor NewFunctor(Atom name, int arity)
	{
	if(arity<1 || arity >250)
		{  puts("Internal error (arity out of range)."); exit(0); }
	/*newfu[arity]++;*/
	return (name&0xffffff) + (((long)arity+1)<<48)+ (2L<<56);
	}

/*
__inline Atom 	FunctorName(Functor f)
	{
	return (f & 0xffffff) + 0x1000000;
	}

__inline int 	FunctorArity(Functor f)
	{
	return (f/0x1000000)-1;
	}
*/

__inline void COMPO_free(Compound c)
	{
	union COMPOpoint *p;
	long int bno,bpos;
#ifdef MTHREAD
	int *ind;
	ind=pthread_getspecific(TermsKey);
#endif
	bno=c&0xfffff0000;
	bno/=0x10000;
	bpos=c&0xffff;
	p=COMPO_buffers[bno]+bpos;
	p->p.this_a=c;
	p->data.next_a=1;
#ifdef MTHREAD
	p->p.next=COMPO_next_free[*ind];
	COMPO_next_free[*ind]=p;
#else
	p->p.next=COMPO_next_free;
	COMPO_next_free=p;
#endif
	COMPO_tfree++;
	return;
	}
#ifdef MTHOPT
__inline void COMPO_free2(Compound c, int ind)
	{
	union COMPOpoint *p;
	int bno,bpos;
	bno=c&0xfffff0000;
	bno/=0x10000;
	bpos=c&0xffff;
	p=COMPO_buffers[bno]+bpos;
	p->p.this_a=c;
	p->data.next_a=1;
#ifdef MTHREAD
	p->p.next=COMPO_next_free[ind];
	COMPO_next_free[ind]=p;
#else
	p->p.next=COMPO_next_free;
	COMPO_next_free=p;
#endif
	COMPO_tfree++;
	return;
	}
#endif
	
__inline union COMPOpoint *COMPO_alloc(Compound *cco)
	{
	union COMPOpoint *p;
	Compound co;
	int cb;
#ifdef MTHREAD
	int *ind;
	ind=pthread_getspecific(TermsKey);
#endif	

#ifdef MTHREAD
	if(COMPO_next_free[*ind]==0)
#else
	if(COMPO_next_free==NULL)
#endif
		{
		unsigned int i;
		if(COMPO_blocks==400096)
			{  puts("Internal error (too much terms)."); exit(0); }
		i=sizeof(union COMPOpoint)*TBSIZE;

#ifdef MTHREAD	
		pthread_mutex_lock(&mtx);	
#endif
		p=COMPO_buffers[COMPO_blocks]=(union COMPOpoint*)malloc(i);
		cb=COMPO_blocks;
		COMPO_blocks++;
#ifdef MTHREAD			
		pthread_mutex_unlock(&mtx);
#endif
		if(p==NULL)
			{  puts("Internal error (lack space for terms)."); exit(0); }
		for(i=0;i<TBSIZE;i++)
			{
#ifdef MTHREAD
			p[i].p.next=COMPO_next_free[*ind];
			COMPO_next_free[*ind]=p+i;
#else
			p[i].p.next=COMPO_next_free;
			COMPO_next_free=p+i;
#endif
			p[i].p.this_a=(4L<<56)+cb*0x10000+i;
			}
		
		}
#ifdef MTHREAD
	p=COMPO_next_free[*ind];
	COMPO_next_free[*ind]=COMPO_next_free[*ind]->p.next;
#else
	p=COMPO_next_free;
	COMPO_next_free=COMPO_next_free->p.next;
#endif
	co=p->p.this_a;
	*cco=co;
	COMPO_tall++;
	return p;
	}

#ifdef MTHOPT	
__inline union COMPOpoint *COMPO_alloc2(Compound *cco, int ind)
	{
	union COMPOpoint *p;
	Compound co;
	int cb;


#ifdef MTHREAD
	if(COMPO_next_free[ind]==0)
#else
	if(COMPO_next_free==NULL)
#endif
		{
		int i;
		if(COMPO_blocks==400096)
			{  puts("Internal error (too much terms)."); exit(0); }
		i=sizeof(union COMPOpoint)*TBSIZE;

#ifdef MTHREAD	
		pthread_mutex_lock(&mtx);	
#endif
		p=COMPO_buffers[COMPO_blocks]=(union COMPOpoint*)malloc(i);
		cb=COMPO_blocks;
		COMPO_blocks++;
#ifdef MTHREAD			
		pthread_mutex_unlock(&mtx);
#endif
		if(p==NULL)
			{  puts("Internal error (lack space for terms)."); exit(0); }
		for(i=0;i<TBSIZE;i++)
			{
#ifdef MTHREAD
			p[i].p.next=COMPO_next_free[ind];
			COMPO_next_free[ind]=p+i;
#else
			p[i].p.next=COMPO_next_free;
			COMPO_next_free=p+i;
#endif
			p[i].p.this_a=(4L<<56)+cb*0x10000+i;
			}
		
		}
#ifdef MTHREAD
	p=COMPO_next_free[ind];
	COMPO_next_free[ind]=COMPO_next_free[ind]->p.next;
#else
	p=COMPO_next_free;
	COMPO_next_free=COMPO_next_free->p.next;
#endif
	co=p->p.this_a;
	*cco=co;
	COMPO_tall++;
	return p;
	}
#endif
	
__inline Compound  NewCompound(Functor fu)
	{
	union COMPOpoint *p;
	Compound c;
	int arity;

	arity=FunctorArity(fu);
	p=COMPO_alloc(&c);
	p->data.a[0]=fu;
	p->data.a[1]=p->data.a[2]=0;
	p->data.next_a=0; p->data.next_p=NULL;
	arity-=2;
	while(arity>0)
		{
		Compound c1;
		union COMPOpoint *p1;
		p1=COMPO_alloc(&c1);
		p->data.next_a=c1;
		p->data.next_p=p1;
		p1->data.a[0]=p1->data.a[1]=p1->data.a[2]=p1->data.next_a=0;
		p1->data.next_p=NULL;
		p=p1;
		arity-=3;
		}
	return c;
	}

#ifdef MTHOPT
__inline Compound  NewCompound2(Functor fu, int m)
	{
	union COMPOpoint *p;
	Compound c;
	int arity;

	arity=FunctorArity(fu);
	p=COMPO_alloc2(&c,m);
	p->data.a[0]=fu;
	p->data.a[1]=p->data.a[2]=0;
	p->data.next_a=0; p->data.next_p=NULL;
	arity-=2;
	while(arity>0)
		{
		Compound c1;
		union COMPOpoint *p1;
		p1=COMPO_alloc2(&c1,m);
		p->data.next_a=c1;
		p->data.next_p=p1;
		p1->data.a[0]=p1->data.a[1]=p1->data.a[2]=p1->data.next_a=0;
		p1->data.next_p=NULL;
		p=p1;
		arity-=3;
		}
	return c;
	}
#endif

/*
__inline Functor CompoundFunctor(Compound c)
	{
	union COMPOpoint *p;
	int bno,bpos;
	bno=c&0xfffff0000;
	bno/=0x10000;
	bpos=c&0xffff;
	p=COMPO_buffers[bno]+bpos;
	return p->data.a[0];
	}
*/

__inline void SetCompoundName(Compound c, Atom name)
	{
	union COMPOpoint *p;
	long int bno,bpos;
	int ar;
	Functor f;
	bno=c&0xfffff0000;
	bno/=0x10000;
	bpos=c&0xffff;
	p=COMPO_buffers[bno]+bpos;
	f=p->data.a[0];
	ar=FunctorArity(f);
	p->data.a[0]=NewFunctor(name,ar);
	}
/*
__inline Term CompoundArg1(Compound c)
	{
	union COMPOpoint *p;
	int bno,bpos;
	bno=c&0xfffff0000;
	bno/=0x10000;
	bpos=c&0xffff;
	p=COMPO_buffers[bno]+bpos;
	return p->data.a[1];
	}

__inline Term CompoundArg2(Compound c)
	{
	union COMPOpoint *p;
	int bno,bpos;
	bno=c&0xfffff0000;
	bno/=0x10000;
	bpos=c&0xffff;
	p=COMPO_buffers[bno]+bpos;
	return p->data.a[2];
	}
*/
__inline Term CompoundArgN(Compound c, int arg)
	{
	union COMPOpoint *p;
	long int bno,bpos;
	bno=c&0xfffff0000;
	bno/=0x10000;
	bpos=c&0xffff;
	p=COMPO_buffers[bno]+bpos;
	while(arg>2)
		{ p=p->data.next_p; arg-=3; }
	return p->data.a[arg];
	}

__inline Term ConsumeCompoundArg(Compound c, int arg)
	{
	union COMPOpoint *p;
	long int bno,bpos;
	Term ret;
	bno=c&0xfffff0000;
	bno/=0x10000;
	bpos=c&0xffff;
	p=COMPO_buffers[bno]+bpos;

	while(arg>2)
		{ p=p->data.next_p; arg-=3; }
	ret=p->data.a[arg];
	p->data.a[arg]=0;
	return ret;
	}


__inline void SetCompoundArg(Compound c, int arg, Term t)
	{
	union COMPOpoint *p;
	long int bno,bpos;
	bno=c&0xfffff0000;
	bno/=0x10000;
	bpos=c&0xffff;
	p=COMPO_buffers[bno]+bpos;

	while(arg>2)
		{ p=p->data.next_p; arg-=3; }
	/*if(p->data.a[arg]!=0) FreeAtomic(p->data.a[arg]);*/
	p->data.a[arg]=t;
	}


__inline Compound MakeCompound(Atom name, int arity)
	{
	return NewCompound(NewFunctor(name,arity));
	}

__inline Compound MakeCompound1(Atom name, Term arg1)
    {
    Term tt;
    tt=NewCompound(NewFunctor(name,1));
    SetCompoundArg(tt,1,arg1);
    return tt;
    }

__inline Compound MakeCompound2(Atom name, Term arg1, Term arg2)
    {
    Term tt;
    tt=NewCompound(NewFunctor(name,2));
    SetCompoundArg(tt,1,arg1);
    SetCompoundArg(tt,2,arg2);
    return tt;
    }

__inline Compound MakeCompound3(Atom name, Term arg1, Term arg2, Term arg3)
    {
    Term tt;
    tt=NewCompound(NewFunctor(name,3));
    SetCompoundArg(tt,1,arg1);
    SetCompoundArg(tt,2,arg2);
    SetCompoundArg(tt,3,arg3);
    return tt;
    }

/*
__inline Atom CompoundName(Compound c)
	{
	return FunctorName(CompoundFunctor(c));
	}

__inline int CompoundArity(Compound c)
	{
	return FunctorArity(CompoundFunctor(c));
	}
*/
			 
__inline int is_atomic(Atomic a)
	{
	return is_atom(a) || is_integer(a) || is_float(a) || 
		is_label(a) || a==0;
	}

__inline int is_float(Atomic a)
	{
	return TERMTYPE(a)==3;
	}

__inline int is_integer(Atomic a)
	{
	return a<0;
	}

__inline int is_atom(Atomic a)
	{
	return TERMTYPE(a)==1;
	}

__inline int is_functor(Atomic a)
	{
	return TERMTYPE(a)==2;
	}

__inline int is_compound(Term a)
	{
	return TERMTYPE(a)==4;
	}

__inline int is_label(Term a)
	{
	return TERMTYPE(a)==6;
	}


__inline Atomic NewLabel(void)
	{
	return ((++LABEL_count)&0xffffffffffffff)+(6L<<56);
	}

__inline long   LabelValue(Atomic l)
	{
	return l-(6L<<56);
	}


void FreeCompound(Compound c)
	{
	union COMPOpoint *p,*ps,*p1;
	long int bno,bpos,arity;
	Compound c1;
	bno=c&0xfffff0000;
	bno/=0x10000;
	bpos=c&0xffff;
	p=COMPO_buffers[bno]+bpos;
	ps=p;
	arity=FunctorArity(p->data.a[0]);


	do
		{
		p1=p->data.next_p;
		c1=p->data.next_a;
		if(p!=ps) FreeAtomic(p->data.a[0]);
		FreeAtomic(p->data.a[1]);
		FreeAtomic(p->data.a[2]);
		COMPO_free(c);
		p=p1;
		c=c1;
		}   while(p!=NULL && c!=0);

	}
	
#ifdef MTHOPT
void FreeCompound2(Compound c, int m)
	{
	union COMPOpoint *p,*ps,*p1;
	int bno,bpos,arity;
	Compound c1;
	bno=c&0xfffff0000;
	bno/=0x10000;
	bpos=c&0xffff;
	p=COMPO_buffers[bno]+bpos;
	ps=p;
	arity=FunctorArity(p->data.a[0]);


	do
		{
		p1=p->data.next_p;
		c1=p->data.next_a;
		if(p!=ps) FreeAtomic2(p->data.a[0],m);
		FreeAtomic2(p->data.a[1],m);
		FreeAtomic2(p->data.a[2],m);
		COMPO_free2(c,m);
		p=p1;
		c=c1;
		}   while(p!=NULL && c!=0);

	}
#endif

void FreeAtomic(Atomic a)
	{
	long type;
	if(a<0) return;
	type = TERMTYPE(a);
	if(type==3)
		{
		union fpoint *p;
		int bno,bpos;
		bno=a&0xff0000;
		bno/=0x10000;
		bpos=a&0xffff;
		p=fbuffers[bno]+bpos;
		p->p.f=a;
		p->p.next=f_next_free;
		f_next_free=p;
		f_ftd++;
		return;
		}
	if(type==4)
		{
		FreeCompound(a);
		return;
		}
	if(type==5)
		{
		FreeList(a);
		return;
		}

	}
	
#ifdef MTHOPT	
void FreeAtomic2(Atomic a, int m)
	{
	long type;
	if(a<0) return;
	type = TERMTYPE(a);
	if(type==3)
		{
		union fpoint *p;
		int bno,bpos;
		bno=a&0xff0000;
		bno/=0x10000;
		bpos=a&0xffff;
		p=fbuffers[bno]+bpos;
		p->p.f=a;
		p->p.next=f_next_free;
		f_next_free=p;
		f_ftd++;
		return;
		}
	if(type==4)
		{
		FreeCompound2(a,m);
		return;
		}
	if(type==5)
		{
		FreeList2(a,m);
		return;
		}

	}
#endif
	
Term CopyTerm(Term t)
     {
     char ty;
     ty=AtomicType(t);
     switch(ty)
         {
         case 'i':
         case 'a':
         case 'e':
         case '?':
         case 'u':
         case 'L':
              return t;
         case 'f':
              {
              double f;
              f=FloatValue(t);
              return NewFloat(f);
              }
         case 'c':
              {
              Term tt;
              int i,ar;
              ar=CompoundArity(t);
              tt=NewCompound(NewFunctor(CompoundName(t),ar));
              for(i=1;i<=ar;i++)
                  SetCompoundArg(tt,i,CopyTerm(CompoundArgN(t,i)));
              return tt;
              }
         case 'l':
              {
              List li,ll;
              li=CopyList(t,0);
              /*ll=t;
              while(!is_empty_list(ll))
                  {
                  li=AppendLast(li,CopyTerm(ListFirst(ll)));
                  ll=ListTail(ll);
                  }*/
              return li;
              }
         }
     return 0;
     }
	 
#ifdef MTHOPT
Term CopyTerm2(Term t, int m)
     {
     char ty;
     ty=AtomicType(t);
     switch(ty)
         {
         case 'i':
         case 'a':
         case 'e':
         case '?':
         case 'u':
         case 'L':
              return t;
         case 'f':
              {
              double f;
              f=FloatValue(t);
              return NewFloat(f);
              }
         case 'c':
              {
              Term tt;
              int i,ar;
              ar=CompoundArity(t);
              tt=NewCompound2(NewFunctor(CompoundName(t),ar),m);
              for(i=1;i<=ar;i++)
                  SetCompoundArg(tt,i,CopyTerm2(CompoundArgN(t,i),m));
              return tt;
              }
         case 'l':
              {
              List li,ll;
              li=CopyList2(t,0,m);
              /*ll=t;
              while(!is_empty_list(ll))
                  {
                  li=AppendLast(li,CopyTerm(ListFirst(ll)));
                  ll=ListTail(ll);
                  }*/
              return li;
              }
         }
     return 0;
     }
	 
#endif
	 
char AtomicType(Atomic a)
	{
	long type;
	if(a<0) return 'i';
	type = TERMTYPE(a);
	if(type==1)
		return 'a';
	if(type==3)	/* 0xff */
		return 'f';
	if(type==4)   /* 0xfe */
		return 'c';
	if(type==5)  /*  0xfd */ 
		return 'l';
	if(type==6)  /*  0xfc */
		return 'L';
	if(type==2)
		return 'u';
	if(a==0)
		return 'e';
	return '?';
	}


void AtomStatistics(void)
	{
	/*int i;*/
	printf("Atoms: %d       , %d a->s  , %d s->a\n",ATOM_count,ATOM_tscount,
				ATOM_stcount);
	printf("Terms: %d blocks use %ld Kbytes, %d allocs ( %d remain) \n",
		COMPO_blocks,COMPO_blocks*(int)sizeof(union COMPOpoint)*64,
		COMPO_tall,COMPO_tall-COMPO_tfree);
	printf("Float: %d blocks, %d allocs, %d frees\n",f_blocks,f_dtf,f_ftd);
	printf("Labels:          %d allocs\n",LABEL_count);
	/*for(i=1;i<30;i++)
		printf("arity %2d: %d allocs\n",i,newfu[i]);*/
	}

void AtomStat1(void)
	{
	printf("%d atoms, ",ATOM_count);
	printf("%d terms (%d blocks), ",COMPO_tall-COMPO_tfree,COMPO_blocks);
	}

long TermMemory(void)
	{
	return (long)COMPO_blocks*(long)sizeof(union COMPOpoint)*64;
	}

int EqualTerms(Term t1, Term t2)
	{
	if(t1==t2)
		return 1;
	if(is_compound(t1) && is_compound(t2) && 
		CompoundFunctor(t1) == CompoundFunctor(t2))
		{
		int i,a;
		a=(int)CompoundArity(t1);
		for(i=1;i<=a;i++)
			if(!EqualTerms(CompoundArgN(t1,i),CompoundArgN(t2,i)))
				return 0;
		return 1;
		}
	if(is_list(t1) && is_list(t2) && ListLength(t1)==ListLength(t2))
		{
		List l1,l2;
		l1=t1;
		l2=t2;
		while(!is_empty_list(l1))
			{
			if(!EqualTerms(ListFirst(l1),ListFirst(l2)))
				return 0;
			l1=ListTail(l1);
			l2=ListTail(l2);
			}
		return 1;
		}
	return 0;
	}


