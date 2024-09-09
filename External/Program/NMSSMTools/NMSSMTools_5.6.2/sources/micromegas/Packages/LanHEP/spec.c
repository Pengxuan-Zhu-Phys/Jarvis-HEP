#include "lanhep.h"


static List specs=0;


static List mk_spec_ind(Term t)
	{
    List ret,l;
    Term t1,t2;
    ret=NewList();
    t=l=CommaToList(t);
    while(!is_empty_list(l))
    	{
    	t1=SpecToRepr(ListFirst(l));
    	if(t1==0)
    		{
    		FreeAtomic(ret);
    		RemoveList(t);
    		return 0;
    		}
    	t2=MakeCompound2(A_I,t1,0);
    	ret=AppendLast(ret,t2);
    	l=ListTail(l);
		}
	RemoveList(t);
    return ret;
    }

	


static void proc_specs(Term t)
	{
    Term ilist;
	if(is_compound(t) && FunctorName(CompoundFunctor(t))==OPR_COMMA)
		{
		Term t1,t2;
		t1=ConsumeCompoundArg(t,1);
		t2=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		proc_specs(t1);
		proc_specs(t2);
		return;
		}
		
   	if(is_atom(t))
		{
		ReportRedefined(t,"special");
        SetAtomProperty(t,PROP_TYPE,OPR_SPECIAL);
		specs=AppendLast(specs,t);
		return;
		}

    if(!is_compound(t) || FunctorName(CompoundFunctor(t))!=OPR_COLON)
		{
		if(!IsTermInput())
			printf("File \"%s\", line %d: ",CurrentInputFile(),
				CurrentInputLine());
		printf("Semantic error: unexpected \'");
		WriteTerm(t);
		printf("\' in \'special\' statement.\n");
		FreeAtomic(t);
		return;
		}
		
    if(!is_atom(CompoundArg1(t)))
		{
		if(!IsTermInput())
			printf("File \"%s\", line %d: ",CurrentInputFile(),
				CurrentInputLine());
		printf("Semantic error: unexpected \'");
		WriteTerm(CompoundArg1(t));
		printf("\' in \'special\' statement.\n");
		FreeAtomic(t);
		return;
		}
	if(is_compound(CompoundArg2(t))&&CompoundName(CompoundArg2(t))==OPR_SECO)
		{
		Integer tp=CompoundArg2(CompoundArg2(t));
		if(!is_integer(tp)||IntegerValue(tp)*IntegerValue(tp)!=1)
			{
			ErrorInfo(0);
			printf("special: symmetry type must be 1 or -1.\n");
			FreeAtomic(t);
			return;
			}
		SetCompoundArg(t,2,ConsumeCompoundArg(CompoundArg2(t),1));
		SetAtomProperty(CompoundArg1(t),A_SYMM,tp);
		}		
    ilist=mk_spec_ind(ConsumeCompoundArg(t,2));
    if(ilist==0)
   		{
   		FreeAtomic(t);
		return;
		}
	specs=AppendLast(specs,CompoundArg1(t));
	ReportRedefined(CompoundArg1(t),"special");
    SetAtomProperty(CompoundArg1(t),PROP_INDEX,ilist);
    SetAtomProperty(CompoundArg1(t),PROP_TYPE,OPR_SPECIAL);
    color_check_spec(CompoundArg1(t),ilist);
   /* WriteTerm(ilist); puts(""); */
    FreeAtomic(t);
	return;
	}



Term ProcessSpecial(Term t, Term ind)
	{
	Term arg;
	char *s;
	arg=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	if(is_atom(arg))
		{
		s=AtomValue(arg);
		if(s[0]=='?' && s[1]==0)
			{
			List li;
			li=specs;
			while(!is_empty_list(li))
				{
				Term t1;
				t1=ListFirst(li);
				WriteTerm(t1);
				printf(":\t");
                WriteTerm(GetAtomProperty(t1,PROP_INDEX));
                puts("");
				li=ListTail(li);
				}
			return 0;
			}
		}
	proc_specs(arg);
	return 0;
	}

int is_special(Atom p, List *ind)
	{
	Term prt;
    prt=GetAtomProperty(p,PROP_TYPE);
    if(prt!=OPR_SPECIAL)
        return 0;
    if(ind!=NULL)
        *ind=GetAtomProperty(p,PROP_INDEX);
    return 1;
	}

void ClearSpecial(Atom sp)
	{
	List l;
	l=specs;
	while(!is_empty_list(l))
		{
		if(ListFirst(l)==sp)
			{
			specs=CutFromList(specs,l);
			break;
			}
		l=ListTail(l);
		}
	RemoveAtomProperty(sp,PROP_TYPE);
	RemoveAtomProperty(sp,PROP_INDEX);
	}

static Term ListToComma(List li)
	{
	if(ListLength(li)==1)
		return ListFirst(li);
	return MakeCompound2(OPR_COMMA,ListFirst(li),ListToComma(ListTail(li)));
	}

Term ProcDelta(Term t, Term ind)
	{
	Term ip, ip1;
	Term mm, dd;
	Term a1;
	
	if(ind!=0 && ListLength(ind)!=2)
		{
		ErrorInfo(329);
		printf("wrong indices number at 'delta' function.\n");
		return 0;
		}
	ip=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	
	if(is_integer(ip) && IntegerValue(ip)>1 && IntegerValue(ip)<50)
		{
		int i,j;
		List l1;
		Term ff;
		l1=NewList();

		for(i=1;i<=IntegerValue(ip);i++)
			{
			Term sua,sual;
			sual=NewList();
			for(j=1;j<=IntegerValue(ip);j++)
				sual=AppendLast(sual,NewInteger(i==j));
			sua=ListToComma(sual);
			FreeAtomic(sual);
			sua=MakeCompound1(A_FBRACET,sua);
			if(ind)
				sua=MakeCompound2(OPR_CARET,sua,ListFirst(ind));
			l1=AppendLast(l1,sua);
			}

		ff=ListToComma(l1);
		RemoveList(l1);
		ff=MakeCompound1(A_FBRACET,ff);
		if(ind)
			ff=MakeCompound2(OPR_CARET,ff,ListNth(ind,2));
		a1=ExprTo1(ff);
		}
	else
		{	
		ip=SpecToRepr(ip);
		ip1=CopyTerm(ip);
		SetCompoundArg(ip1,1,CompoundArg2(ip));
		SetCompoundArg(ip1,2,CompoundArg1(ip));
		ip=MakeCompound2(A_I,ip,NewLabel());
		ip1=MakeCompound2(A_I,ip1,NewLabel());
/*		WriteTerm(ip);
		WriteTerm(ip1);
		puts("");*/
		dd=MakeCompound2(OPR_SPECIAL,MakeList2(CopyTerm(ip),CopyTerm(ip1)),A_DELTA);
		mm=MakeCompound(A_MTERM,4);
		SetCompoundArg(mm,1,NewInteger(1));
		SetCompoundArg(mm,2,NewInteger(1));
		SetCompoundArg(mm,3,MakeList1(dd));
		a1=MakeCompound2(A_ALG1,MakeList1(mm),MakeList2(ip,ip1));
		}
	if(ind==0)
		return a1;
	else
		return MakeCompound2(OPR_CARET,a1,
				MakeCompound2(OPR_CARET,ListFirst(ind),ListNth(ind,2)));
	}
	
