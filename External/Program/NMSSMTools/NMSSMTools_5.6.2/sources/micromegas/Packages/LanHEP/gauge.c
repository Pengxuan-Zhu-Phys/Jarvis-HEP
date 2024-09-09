#include <string.h>
#include "lanhep.h"

Term ProcGauge(Term t, Term ind)
	{
	Term fld,grp,cnst;
	if(!is_compound(t) || CompoundArity(t)!=3)
		{
		ErrorInfo(224);
		printf("wrong arguments number in call ");
		WriteTerm(t);
		puts(".");
		return 0;
		}
	fld=ConsumeCompoundArg(t,1);
	grp=ConsumeCompoundArg(t,2);
	cnst=ConsumeCompoundArg(t,3);
	FreeAtomic(t);
	
	
	
	return 0;
	}
	
	
Term ProcGenGauge(Term t, Term ind)
	{
	
	
	return 0;
	}
	
static int ord_ind(List pind, Term a1, Atom p)
{
	List l1, l2, l3, l4;
	
	l1=ConsumeCompoundArg(a1,2);
	if(ListLength(l1)!=ListLength(pind))
	{
		ErrorInfo(702);
		printf("brst_transform: expression for particle %s with different indices\n",
				AtomValue(p));
		return 0;
	}
	
	if(pind==0)
		return 1;
	
	l2=NewList();
	
	for(l3=pind;l3;l3=ListTail(l3))
	{
		Term t;
		t=CompoundArg1(ListFirst(l3));
		for(l4=l1;l4;l4=ListTail(l4))
		{
			if(EqualTerms(CompoundArg1(ListFirst(l4)),t))
				break;
		}
		if(l4==0)
		{
			ErrorInfo(702);
			printf("brst_transform: expression for particle %s with different indices\n",
				AtomValue(p));
			return 0;
		}
		l2=AppendLast(l2,ListFirst(l4));
		ChangeList(l4,0);
		l1=CutFromList(l1,l4);
	}
	
	SetCompoundArg(a1,2,l2);
	return 1;
	
}
	
Term ProcBRSTTransf(Term t, Term ind)
{
	Term t1;
	Atom tp;
	
	if(!is_compound(t) || CompoundArity(t)!=1)
	{
		ErrorInfo(225);
		puts("brst_transform: wrong syntax");
		FreeAtomic(t);
		return 0;
	}
	
	tp=CompoundName(t);
	if(tp!=A_BRST_TRANSF && tp!=A_BRSTI_TRANSF)
	{
		ErrorInfo(225);
		puts("brst_transform: should be 'brst_transform' or 'brsti_transform'.");
		FreeAtomic(t);
		return 0;
	}
	
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	t=CommaToList(t1);
	
	for(t1=t;t1;t1=ListTail(t1))
	{
		Term t2, a1;
		List pind;
		t2=ListFirst(t1);
		
		if(!is_compound(t2)||CompoundArity(t2)!=2||CompoundName(t2)!=OPR_RARROW)
		{
			ErrorInfo(701);
			printf("brst_transform: '");
			WriteTerm(t2);
			puts("' is illegal");
			continue;
		}
		
		if(!is_atom(CompoundArg1(t2)))
		{
			ErrorInfo(701);
			printf("brst_transform: '");
			WriteTerm(t2);
			puts("' is illegal");
			continue;
		}
		
		if(!is_particle(CompoundArg1(t2),&pind))
		{
			ErrorInfo(701);
			printf("brst_transform: '");
			WriteTerm(CompoundArg1(t2));
			puts("' is not a particle.");
			continue;
		}
		
		a1=ExprTo1(ConsumeCompoundArg(t2,2));
		if(a1==0)
			continue;
		
		if(ord_ind(pind,a1,CompoundArg1(t2))==0)
			continue;
		
/*		printf("prt %s ",AtomValue(CompoundArg1(t2)));
		WriteTerm(pind);puts("");
		WriteTerm(a1);
		puts("");
*/		
		SetAtomProperty(CompoundArg1(t2),tp,a1);
		
	}
	
	
	return 0;
}

static Term sub_brst(Term m1, int pos, Term a1, int chs)
{
	List l1;
	List l2;
	
	if(chs<0)
	{
		SetCompoundArg(m1,1,NewInteger(-IntegerValue(CompoundArg1(m1))));
	}
	
	l1=ListNthList(CompoundArgN(m1,3),pos);
	l2=ConsumeCompoundArg(ListFirst(l1),1);
	FreeAtomic(ListFirst(l1));
	ChangeList(l1,MakeCompound2(A_ALG1,l2,a1));
	return m1;
}

Term alg1_proc_brst(Term a1, Atom tp)
{
	List l1, l2, l3;
	
	
	l1=ConsumeCompoundArg(a1,1);
	
	l2=NewList();
	
	for(l3=l1;l3;l3=ListTail(l3))
	{
		Term m1;
		List l4,l5,l6;
		int pos=0;
		int sn=1, chs=0;
		Term mm;
		
		m1=ListFirst(l3);
		
		/* deriv -> derivp */
		
		l4=ConsumeCompoundArg(m1,3);
		l5=NewList();
		mm=0;
		for(l6=l4;l6;l6=ListTail(l6))
		{
			Term r;
			r=ListFirst(l6);
			l5=AppendLast(l5,r);
			if(CompoundName(r)==OPR_SPECIAL &&
					CompoundArg2(r)==A_MOMENT)
				mm=r;
			if(mm && CompoundName(r)==OPR_FIELD)
			{
				SetCompoundArg(mm,2,A_MOMENT_S);
				mm=0;
				l5=AppendLast(l5,MakeCompound2(OPR_SPECIAL,0,A_MOMENT_E));
			}
		}
		RemoveList(l4);
		SetCompoundArg(m1,3,l5);
		
		for(l4=CompoundArgN(m1,3);l4;l4=ListTail(l4))
		{
			Term sa;
			
			pos++;
			if(chs)
			{
				sn=-sn;
				chs=0;
			}
			
			if(CompoundName(ListFirst(l4))!=OPR_FIELD)
				continue;
			if(CompoundArg2(ListFirst(l4))==A_VEV)
				continue;
			sa=GetAtomProperty(CompoundArg2(ListFirst(l4)),PROP_TYPE);
			if(is_compound(sa) && CompoundName(sa)==OPR_FIELD &&
				(CompoundArg2(sa)==NewInteger(2)||CompoundArg2(sa)==NewInteger(2)))
				sn=-sn;
			sa=GetAtomProperty(CompoundArg2(ListFirst(l4)),tp);
			if(sa==NewInteger(0))
				continue;
			if(sa==0)
			{
				printf("brst: transform for '%s' is not defined.\n",
						AtomValue(CompoundArg2(ListFirst(l4))));
				SetAtomProperty(CompoundArg2(ListFirst(l4)),tp,NewInteger(0));
				continue;
			}
			
			l2=AppendLast(l2,sub_brst(CopyTerm(m1),pos,CopyTerm(sa),sn));
		}
		
	}
	
	FreeAtomic(l1);
	l2=SetIntAlgs(l2);
	SetCompoundArg(a1,1,l2);
	
	return a1;
}

		
Term ProcBRST(Term t, Term ind)
{
	Atom tp;
	char *s;
	Term t1;
	
	if(!is_compound(t) || CompoundArity(t)!=1)
	{
		ErrorInfo(717);
		puts("brst: wrong syntax");
		return NewInteger(0);
	}
	
	t1=ConsumeCompoundArg(t,1);
	
	s=AtomValue(CompoundName(t));
	if(strcmp(s,"brst")==0)
		tp=A_BRST_TRANSF;
	else if(strcmp(s,"brsti")==0)
		tp=A_BRSTI_TRANSF;
	else
	{
		ErrorInfo(715);
		puts("brst: should be 'brst' or 'brsti'.");
		return NewInteger(0);
	}
	
	FreeAtomic(t);
	t=t1;
	
	t=ExprTo1(t);
	if(t==0)
		return NewInteger(0);
	
	t=alg1_proc_brst(t,tp);
	
	if(CompoundArg1(t)==0)
	{
		FreeAtomic(t);
		return NewInteger(0);
	}
	else
	{
		if(ind)
			t=il_to_caret(t,ind);
		
		return t;
	}
	return 0;
}


