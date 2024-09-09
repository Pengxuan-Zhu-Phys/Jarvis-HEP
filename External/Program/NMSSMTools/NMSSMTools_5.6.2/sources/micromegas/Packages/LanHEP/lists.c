#include <stdlib.h>
#include "lanhep.h"

#define TBSIZE 65500

#ifdef MTHREAD
#include <pthread.h>

extern pthread_key_t TermsKey;
static pthread_mutex_t mtx = PTHREAD_MUTEX_INITIALIZER;
#endif

#ifdef MTHREAD
union LISTpoint *LISTS_next_free[17]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};	
#else
union LISTpoint *LISTS_next_free=NULL;	
#endif
union LISTpoint *LISTS_buffers[400096];
int LISTS_blocks=0, LISTS_allocs=0, LISTS_free=0;

Term CopyTerm2(Term, int);
void FreeAtomic2(Atomic a, int m);
void RemoveList2(List li, int m);
/*
#ifndef USE_INLINE
#define __inline
#include "lists.ci"
#undef __inline
#endif
*/

/************************************/
/*********** File lists.ci **********/
/************************************/

#define __inline

/*
__inline int is_empty_list(List l)
	{
	return l==0;
	}

	
__inline union LISTpoint *LIST_pointer(List l)
	{
	register int bno,bpos;
	if(is_empty_list(l))
		return NULL;
	bno=l&0xfff0000;
	bno/=0x10000;
	bpos=l&0xffff;
	return LISTS_buffers[bno]+bpos;
	}	
*/
		
__inline void LIST_free(List l, union LISTpoint *p)
	{
#ifdef MTHREAD
	int *ind;
	ind=pthread_getspecific(TermsKey);
	p->p.next=LISTS_next_free[*ind];
	p->p.l=l;
	LISTS_next_free[*ind]=p;
	LISTS_free++;
#else
	p->p.next=LISTS_next_free;
	p->p.l=l;
	LISTS_next_free=p;
	LISTS_free++;
#endif
	}

#ifdef MTHOPT
__inline void LIST_free2(List l, union LISTpoint *p, int ind)
	{
#ifdef MTHREAD
	p->p.next=LISTS_next_free[ind];
	p->p.l=l;
	LISTS_next_free[ind]=p;
	LISTS_free++;
#else
	p->p.next=LISTS_next_free;
	p->p.l=l;
	LISTS_next_free=p;
	LISTS_free++;
#endif
	}
#endif
		
__inline union LISTpoint *LIST_alloc(List *li)
	{
	union LISTpoint *p;
	int lb;
#ifdef MTHREAD
	int *ind;
	ind=pthread_getspecific(TermsKey);
#endif
	
#ifdef MTHREAD
	if(LISTS_next_free[*ind]==NULL)
#else
	if(LISTS_next_free==NULL)
#endif
		{
		int i;
		if(LISTS_blocks==400096)
			{  puts("Internal error (too much lists)."); exit(0); }
		i=sizeof(union LISTpoint)*TBSIZE;

#ifdef MTHREAD	
		pthread_mutex_lock(&mtx);	
#endif		
		p=LISTS_buffers[LISTS_blocks]=(union LISTpoint*)malloc((unsigned)i);
		lb=LISTS_blocks;
		LISTS_blocks++;
#ifdef MTHREAD			
		pthread_mutex_unlock(&mtx);
#endif		
		if(p==NULL)
			{  printf("Internal error (lack space for lists, %d blocks of %d bytes).\n",
					LISTS_blocks, i); exit(0); }
		for(i=0;i<TBSIZE;i++)
			{
#ifdef MTHREAD
			p[i].p.next=LISTS_next_free[*ind];
			LISTS_next_free[*ind]=p+i;
#else			
			p[i].p.next=LISTS_next_free;
			LISTS_next_free=p+i;
#endif	
			p[i].p.l=(5L<<56)+lb*0x10000+i;
			}
		
		}
#ifdef MTHREAD
	p=LISTS_next_free[*ind];
	LISTS_next_free[*ind]=p->p.next;
#else	
	p=LISTS_next_free;
	LISTS_next_free=p->p.next;
#endif	
	*li=p->p.l;
	p->d.value=0;
	p->d.nextp=NULL;
	p->d.nextl=0;
	LISTS_allocs++;

	return p;
	}

#ifdef MTHOPT	
__inline union LISTpoint *LIST_alloc2(List *li, int ind)
	{
	union LISTpoint *p;
	int lb;

	
#ifdef MTHREAD
	if(LISTS_next_free[ind]==NULL)
#else
	if(LISTS_next_free==NULL)
#endif
		{
		int i;
		if(LISTS_blocks==400096)
			{  puts("Internal error (too much lists)."); exit(0); }
		i=sizeof(union LISTpoint)*TBSIZE;

#ifdef MTHREAD	
		pthread_mutex_lock(&mtx);	
#endif		
		p=LISTS_buffers[LISTS_blocks]=(union LISTpoint*)malloc(i);
		lb=LISTS_blocks;
		LISTS_blocks++;
#ifdef MTHREAD			
		pthread_mutex_unlock(&mtx);
#endif		
		if(p==NULL)
			{  printf("Internal error (lack space for lists, %d blocks of %d bytes).\n",
					LISTS_blocks, i); exit(0); }
		for(i=0;i<TBSIZE;i++)
			{
#ifdef MTHREAD
			p[i].p.next=LISTS_next_free[ind];
			LISTS_next_free[ind]=p+i;
#else			
			p[i].p.next=LISTS_next_free;
			LISTS_next_free=p+i;
#endif	
			p[i].p.l=(5L<<56)+lb*0x10000+i;
			}
		
		}
#ifdef MTHREAD
	p=LISTS_next_free[ind];
	LISTS_next_free[ind]=p->p.next;
#else	
	p=LISTS_next_free;
	LISTS_next_free=p->p.next;
#endif	
	*li=p->p.l;
	p->d.value=0;
	p->d.nextp=NULL;
	p->d.nextl=0;
	LISTS_allocs++;

	return p;
	}
#endif	
	
__inline List NewList(void)
	{
	return 0;
	}
	
/*	
__inline Term ListFirst(List l)
	{
	union LISTpoint *p;
	p=LIST_pointer(l);
	return p->d.value;
	}
	
__inline List ListTail(List l)
	{
	union LISTpoint *p;
	p=LIST_pointer(l);
	return p->d.nextl;
	}	
*/
__inline List AppendFirst(List li,Term t)
	{
	union LISTpoint *p1, *p2;
	List l1;
	p1=LIST_alloc(&l1);
	p1->d.value=t;
	if(li==0)
		return l1;
	p2=LIST_pointer(li);
	p1->d.nextl=li;
	p1->d.nextp=p2;
	return l1;
	}
	
__inline List AppendLast(List li,Term t)
	{
	union LISTpoint *p1, *p2;
	List l1;
	p1=LIST_alloc(&l1);
	p1->d.value=t;
	if(li==0)
		return l1;
	p2=LIST_pointer(li);
	while(p2->d.nextp!=NULL)
		p2=p2->d.nextp;
	p2->d.nextl=l1;
	p2->d.nextp=p1;
	return li;
	}
#ifdef MTHOPT
__inline List AppendFirst2(List li,Term t, int m)
	{
	union LISTpoint *p1, *p2;
	List l1;
	p1=LIST_alloc2(&l1,m);
	p1->d.value=t;
	if(li==0)
		return l1;
	p2=LIST_pointer(li);
	p1->d.nextl=li;
	p1->d.nextp=p2;
	return l1;
	}
	
__inline List AppendLast2(List li,Term t, int m)
	{
	union LISTpoint *p1, *p2;
	List l1;
	p1=LIST_alloc2(&l1,m);
	p1->d.value=t;
	if(li==0)
		return l1;
	p2=LIST_pointer(li);
	while(p2->d.nextp!=NULL)
		p2=p2->d.nextp;
	p2->d.nextl=l1;
	p2->d.nextp=p1;
	return li;
	}
#endif
	
__inline int ListLength(List l)
	{
	union LISTpoint *p;
	int len;
	len=0;
	if(l==0) return 0;
	p=LIST_pointer(l);
	while(p!=NULL)
		{
		p=p->d.nextp;
		len++;
		}
	return len;
	}


__inline void ChangeList(List li, Term nv)
     {
     union LISTpoint *p;
     p=LIST_pointer(li);
     /*FreeAtomic(p->d.value);*/
     p->d.value=nv;
     }

__inline int is_list(Term a)
	{
	return TERMTYPE(a) == 5;
	}



	
Term ListNth(List l, int arg)
	{
	union LISTpoint *p;
	p=LIST_pointer(l);
	while(arg>1)
		{
		p=p->d.nextp;
		arg--;
		}
	return p->d.value;
	}
	
List ListNthList(List l, int arg)
	{
	union LISTpoint *p;
	if(arg<2)
		return l;
	p=LIST_pointer(l);
	while(arg>2)
		{
		p=p->d.nextp;
		arg--;
		}
	return p->d.nextl;
	}

	
List ConcatList(List li, List l1)
	{
	union LISTpoint *p1, *p2;
	if(li==0)
		return l1;
	if(l1==0)
		return li;
	p1=LIST_pointer(l1);
	p2=LIST_pointer(li);
	while(p2->d.nextp!=NULL)
		p2=p2->d.nextp;
	p2->d.nextl=l1;
	p2->d.nextp=p1;
	return li;
	}

List CopyList(List l, List *e)
{
	List ret=0, ee=0;
	
	for(;l;l=ListTail(l))
	{
		Term t=CopyTerm(ListFirst(l));
		if(ret==0)
		{
			ret=AppendFirst(ret,t); ee=ret;
		}
		else
		{
			ee=AppendLast(ee,t); ee=ListTail(ee);
		}
	}
	if(e)
		*e=ee;
	return ret;
}
#ifdef MTHOPT
List CopyList2(List l, List *e, int m)
{
	List ret=0, ee=0;
	
	for(;l;l=ListTail(l))
	{
		Term t=CopyTerm2(ListFirst(l),m);
		if(ret==0)
		{
			ret=AppendFirst2(ret,t,m); ee=ret;
		}
		else
		{
			ee=AppendLast2(ee,t,m); ee=ListTail(ee);
		}
	}
	if(e)
		*e=ee;
	return ret;
}
#endif

List ListSplit(List list, List cut, List *rest) 	
	{
	union LISTpoint *pi,*pcut, *pl;
	
	pi=LIST_pointer(list);
	pcut=LIST_pointer(cut);
	
	if(list==cut)
		{
		*rest=pi->d.nextl;
		pi->d.nextp=NULL;
		pi->d.nextl=0;
		return 0;
		}
		
	pl=pi;
	while(pl->d.nextl!=cut)
		pl=pl->d.nextp;
	pcut=pl->d.nextp;
	pl->d.nextp=NULL;
	pl->d.nextl=0;
	*rest=pcut->d.nextl;
	pcut->d.nextp=NULL;
	pcut->d.nextl=0;
	return list;
	}
/*	
void FreeList(List li)
	{
	union LISTpoint *l,*l1;
	List li1;
	l=LIST_pointer(li);
	while(l!=NULL)
		{
		l1=l->d.nextp;
		li1=l->d.nextl;
		FreeAtomic(l->d.value);
		LIST_free(li,l);
		l=l1;
		li=li1;
		}
	}
*/

void FreeList(List li)
	{
	union LISTpoint *l;
	l=LIST_pointer(li);
	while(l!=NULL)
		{
		FreeAtomic(l->d.value);
		l=l->d.nextp;
		}
	RemoveList(li);   
	}
#ifdef MTHOPT
void FreeList2(List li, int m)
	{
	union LISTpoint *l;
	l=LIST_pointer(li);
	while(l!=NULL)
		{
		FreeAtomic2(l->d.value,m);
		l=l->d.nextp;
		}
	RemoveList2(li,m);   
	}
#endif
		
void RemoveList(List li)
	{
	union LISTpoint *l,*l1;
	List li1;
	if(li==0) return;
	l=LIST_pointer(li);
	while(l!=NULL)
		{
		l1=l->d.nextp;
		li1=l->d.nextl;
		LIST_free(li,l);
		l=l1;
		li=li1;
		}

	}	
	
#ifdef MTHOPT
void RemoveList2(List li, int m)
	{
	union LISTpoint *l,*l1;
	List li1;
	if(li==0) return;
	l=LIST_pointer(li);
	while(l!=NULL)
		{
		l1=l->d.nextp;
		li1=l->d.nextl;
		LIST_free2(li,l,m);
		l=l1;
		li=li1;
		}

	}	
#endif

List CutFromList(List list, List cut)
     {
     union LISTpoint *pi,*pcut, *pl;
     List rest;
	pi=LIST_pointer(list);
	pcut=LIST_pointer(cut);

	if(list==cut)
		{
		rest=pi->d.nextl;
		pi->d.nextp=NULL;
		pi->d.nextl=0;
        FreeList(list);
		return rest;
		}
		
	pl=pi;
	while(pl->d.nextl!=cut)
		pl=pl->d.nextp;
	pcut=pl->d.nextp;
	pl->d.nextp=pcut->d.nextp;
	pl->d.nextl=pcut->d.nextl;

	pcut->d.nextp=NULL;
	pcut->d.nextl=0;
    FreeList(cut);
	return list;
	}


void InsertList(List l,Term ins)
	{
	union LISTpoint *p,*p1;
	List l1;
	p=LIST_pointer(l);
	p1=LIST_alloc(&l1);
	p1->d.value=ins;
	p1->d.nextl=p->d.nextl;
	p1->d.nextp=p->d.nextp;
	p->d.nextl=l1;
	p->d.nextp=p1;
	}
	



void ListStatistics(void)
	{
	printf("Lists: %d blocks use %ld Kbytes, %d allocs ( %d remain)\n",
		LISTS_blocks,LISTS_blocks*(int)sizeof(union LISTpoint)*64,
		LISTS_allocs,LISTS_allocs-LISTS_free);
	}

long ListMemory(void)
	{
	return (long)LISTS_blocks*(long)sizeof(union LISTpoint)*64;
	}

void ListStat1(void)
	{
	printf("%d lists (%d blocks)\n",LISTS_allocs-LISTS_free,LISTS_blocks);
	}





void DumpList(List l)
	{
	printf("\t[\n\t");
	while(!is_empty_list(l))
		{
		WriteTerm(ListFirst(l));
		printf("\n\t");
		l=ListTail(l);
		}
	printf("]\n");
	}


int ListMember(List l, Term t)
	{
	int i=1;
	while(!is_empty_list(l))
		{
		if(ListFirst(l)==t)
			return i;
		i++;
		l=ListTail(l);
		}
	return 0;
	}

List IntersectList(List l1, List l2)
	{
	List li;
	li=NewList();
	while(!is_empty_list(l1))
		{
		if(ListMember(l2,ListFirst(l1)))
			li=AppendLast(li,CopyTerm(ListFirst(l1)));
		l1=ListTail(l1);
		}
	return li;
	}


List SortedList(List li, int (*cmp)(Term,Term))
	{
	List l,l1,l2,li1;
	li1=li;
	if(is_empty_list(li))
		return NewList();
	l=AppendFirst(NewList(),ListFirst(li));
	li=ListTail(li);
	while(!is_empty_list(li))
		{
		Term t;
		t=ListFirst(li);
		if(cmp(t,ListFirst(l))<=0)
			{
			l=AppendFirst(l,t);
			goto cnt;
			}
		l1=l;
		l2=ListTail(l);
		while(!is_empty_list(l2))
			{
			if(cmp(t,ListFirst(l2))<=0)
				{
				InsertList(l1,t);
				goto cnt;
				}
			l1=l2;
			l2=ListTail(l2);
			}
		InsertList(l1,t);
	cnt:
		li=ListTail(li);
		}
	RemoveList(li1);
	return l;
	}

List MakeList1(Term t1)
	{
	register List l;
	l=AppendFirst(NewList(),t1);
	return l;
	}

	
List MakeList2(Term t1, Term t2)
	{
	register List l;
	l=AppendFirst(NewList(),t2);
	l=AppendFirst(l,t1);
	return l;
	}
	
List MakeList3(Term t1, Term t2, Term t3)
	{
	register List l;
	l=AppendFirst(NewList(),t3);
	l=AppendFirst(l,t2);
	l=AppendFirst(l,t1);
	return l;
	}
	
List MakeList4(Term t1, Term t2, Term t3, Term t4)
	{
	register List l;
	l=AppendFirst(NewList(),t4);
	l=AppendFirst(l,t3);
	l=AppendFirst(l,t2);
	l=AppendFirst(l,t1);
	return l;
	}
		

	
/*
static char cbuf[100];


static int t_cmp(Term t1, Term t2)
	{
	char c1,c2;
	c1=AtomicType(t1);
	c2=AtomicType(t2);
	if(c1<c2) return -1;
	if(c1>c2) return 1;
	switch(c1)
		{
		case 'i':
			if(t1==t2) return 0;
			if(IntegerValue(t1)<IntegerValue(t2))
				return -1;
			return 1;
		case 'a':
			return strcmp(AtomValue(t1),AtomValue(t2));
		case 'f':
			if(FloatValue(t1)==FloatValue(t2)) return 0;
			if(FloatValue(t1)<FloatValue(t2))
				return -1;
			return 1;
		default: return 0;
		}
	return 0;
	}	
	
static void displlist(List l)
	{
	Atomic last_a=0;
	int reps=0;
	while(!is_empty_list(l))
		{
		Atomic a;
		a=ListFirst(l);
		if(last_a==0)
			{
			last_a=a;
			reps=1;
			l=ListTail(l);
			continue;
			}
		if(a==last_a)
			{
			reps++;
			l=ListTail(l);
			continue;
			}
		DisplayTerm(last_a);
		printf("\t\t:%d\n",reps);
		reps=1;
		last_a=a;
		l=ListTail(l);
		}
	
	DisplayTerm(last_a);
	printf("\t\t:%d\n",reps);
	}

*/
