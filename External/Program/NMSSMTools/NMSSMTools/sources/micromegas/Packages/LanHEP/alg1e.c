#include <setjmp.h>
#include "lanhep.h"

extern jmp_buf alg1_jmp_buf;

static Term SetInd_e(Term t, Term *ir, Term *is);

static void set_null(List al)
{
	List il=0,l1,l2;
	for(l1=al;l1;l1=ListTail(l1))
		if(CompoundArg2(ListFirst(l1)))
		{
			il=CompoundArg2(ListFirst(l1));
			break;
		}
		
	if(is_empty_list(il))
		return;
	
	for(l1=al;l1;l1=ListTail(l1))
		if(is_empty_list(CompoundArg2(ListFirst(l1))))
		{
			Term m,w,i;
			if(CompoundArg1(ListFirst(l1)))
				return;
			
			i=CopyTerm(il);
			for(l2=i;l2;l2=ListTail(l2))
				SetCompoundArg(ListFirst(l2),2,NewLabel());
			
			w=MakeCompound2(OPR_WILD,CopyTerm(i),0);
			m=MakeCompound(A_MTERM,4);
			SetCompoundArg(m,1,NewInteger(1));
			SetCompoundArg(m,2,NewInteger(1));
			SetCompoundArg(m,3,MakeList1(w));
			
			SetCompoundArg(ListFirst(l1),1,MakeList1(m));
			SetCompoundArg(ListFirst(l1),2,i);
		}

}

Term alg1_mk_wild(Term f, List *rll, List *sll)
	{
	Term t,gr;
	List l;
/*	WriteTerm(f); puts("");*/
	t=ConsumeCompoundArg(f,1);
	FreeAtomic(f);
	if(!is_compound(t))
		{
		ErrorInfo(320);
		printf("bad expression in array");
		longjmp(alg1_jmp_buf,1);
		}
	if(CompoundName(t)==OPR_SECO)
		{
		f=ConsumeCompoundArg(t,1);
		gr=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		f=CommaToList(f);
		}
	else
		{
		f=CommaToList(t);
		gr=MakeCompound1(OPR_WILD,NewInteger(ListLength(f)));
		}
	if(ListLength(f)==1)
		{
		ErrorInfo(321);
		printf("1 component in array\n");
		longjmp(alg1_jmp_buf,1);
		}
	gr=SpecToRepr(gr);
	if(gr==0)
		longjmp(alg1_jmp_buf,1);
	
	l=f;
	while(!is_empty_list(l))
		{
		ChangeList(l,ExprTo1(ListFirst(l)));
		l=ListTail(l);
		}
	
	set_null(f);
		
	SetInd_e(f,rll,sll);
	
	t=MakeCompound(OPR_WILD,3);
	SetCompoundArg(t,1,gr);
	SetCompoundArg(t,3,NewInteger(ListLength(f)));
	*rll=AppendLast(*rll,t);
/*	puts("    making wild :");
	WriteTerm(*rll); puts("	: indices");
	DumpList(f);	*/
	return f;
	
	}


static List resort_ind_e(List ii, List oind)
	{
	List ret;
	ret=NewList();
	while(!is_empty_list(oind))
		{
		Term t;
		t=ListFirst(oind);
		if(is_atom(CompoundArg2(t)))
			{
			List l;
			l=ii;
			while(!is_empty_list(l))
				{
				Term t1;
				t1=ListFirst(l);
				if(CompoundArg2(t1)==CompoundArg2(t))
					{
					if(!EqualTerms(CompoundArg1(t),CompoundArg1(t1)))
						{
						ErrorInfo(322);
						printf(" misused index '%s'\n",AtomValue(CompoundArg2(t1)));
						longjmp(alg1_jmp_buf,1);
						}
					ChangeList(l,0);
					ret=AppendLast(ret,t1);
					ii=CutFromList(ii,l);
					goto fi;
					}
				l=ListTail(l);
				}
			return 0;
			}
		else
			{
			List l;
			l=ii;
			while(!is_empty_list(l))
				{
				Term t1;
				t1=ListFirst(l);
				if(is_label(CompoundArg2(t1)) &&
					EqualTerms(CompoundArg1(t),CompoundArg1(t1)))
						{
						ChangeList(l,0);
						ret=AppendLast(ret,t1);
						ii=CutFromList(ii,l);
						goto fi;
						}
				l=ListTail(l);
				}
			return 0;
			}
			
		fi:	
		oind=ListTail(oind);
		}
		
	return ret;
	}	
				

		
		
static Term SetInd_e(Term t, Term *ir, Term *is)
	{
	List li,oind,si;
	int first=1, fil=0;
	li=t;
	oind=si=NewList();
	while(!is_empty_list(li))
		{
		Term ct,ii;
		ct=ListFirst(li);
		ii=ConsumeCompoundArg(ct,2);
		if(first)
			{
			oind=CopyTerm(ii);
			SetCompoundArg(ct,2,ii);
			fil=ListLength(oind);
			first=0;
			}
		else
			{
			if(ListLength(ii)!=fil)
				{
				ErrorInfo(323);
				printf("expressions in array have different indices.\n");
				longjmp(alg1_jmp_buf,1);
				}
			if(fil!=0)
				{
				ii=resort_ind_e(ii,oind);
				if(ii==0)
					{
					ErrorInfo(323);
					printf("expressions in array  have different indices.\n");
					longjmp(alg1_jmp_buf,1);
					}
				SetCompoundArg(ct,2,ii);
				}
			}
			
		li=ListTail(li);
		
		}
		
	li=oind;
	while(!is_empty_list(li))
		{
		Term lll;
		lll=CompoundArg2(ListFirst(li));
		if(is_atom(lll))
			si=AppendLast(si,lll);
		SetCompoundArg(ListFirst(li),2,0);
		li=ListTail(li);
		}
		
	*is=ConcatList(si,*is);
	*ir=oind;
	return t;
	}
	
static int infno[100];
	
static void sel_inf(Term w, int o)
{
	List l1,l2,l3;
	for(l1=CompoundArg2(w);l1;l1=ListTail(l1))
	{
		Term a1=ListFirst(l1);
		int ch=0;
		for(l2=CompoundArg1(a1);l2;l2=ListTail(l2))
		{
			int po=0;
			for(l3=CompoundArgN(ListFirst(l2),3);l3;l3=ListTail(l3))
			{
				Term prp=0;
				if(CompoundName(ListFirst(l3))==OPR_PARAMETER)
					prp=GetAtomProperty(CompoundArg2(ListFirst(l3)),
							A_INFINITESIMAL);
				if(prp && IntegerValue(CompoundArg1(prp))>0)
					po+=(int)IntegerValue(CompoundArg1(prp));
			}
			if(po>9)
				po=9;
			if(o==-1)
				infno[po]++;
			else if(po!=o)
			{
				FreeAtomic(ListFirst(l2));
				ChangeList(l2,0);
				ch++;
			}
		}
		if(ch)
		{
			l2=ConsumeCompoundArg(a1,1);
			do
			{
				for(l3=l2;l3;l3=ListTail(l3))
				if(ListFirst(l3)==0)
				{
					l2=CutFromList(l2,l3);
					break;
				}
			} while(l3);
			SetCompoundArg(a1,1,l2);
					
		}
	}
}

void alg1_inf_wild(Term sub)
{
	List l,l1,nl;
	Term w;
	int i,no=0;
	return;
	l=CompoundArg1(sub);
	if(ListLength(l)!=1)
		return;
	l1=CompoundArgN(ListFirst(l),3);
	if(l1==0 || ListLength(l1)!=1 || !is_compound(ListFirst(l1))
			|| CompoundName(ListFirst(l1))!=OPR_WILD)
		return;
	w=ListFirst(l1);
	for(l=CompoundArg1(w);l;l=ListTail(l))
		if(CompoundName(ListFirst(l))==OPR_WILD)
			no++;
	if(no!=1)
		return;
	
	for(i=0;i<10;i++)
		infno[i]=0;
	sel_inf(w,-1);
	no=0;
	for(i=1;i<10;i++)
		if(infno[i])
			no++;
	if(no==0)
		return;
	
	nl=0;
	for(i=1;i<10;i++)
	{
		Term nw;
		if(infno[i]==0)
			continue;
		nw=CopyTerm(ListFirst(CompoundArg1(sub)));
		sel_inf(ListFirst(CompoundArgN(nw,3)),i);
		AppendLast(CompoundArgN(nw,3),
				MakeCompound2(A_INFINITESIMAL,0,NewInteger(i)));
		nl=AppendLast(nl,nw);
	}
	sel_inf(w,0);
	ConcatList(CompoundArg1(sub),nl);
	

	/*DumpList(CompoundArg1(sub));
	puts("\n");*/	
	
}
	
