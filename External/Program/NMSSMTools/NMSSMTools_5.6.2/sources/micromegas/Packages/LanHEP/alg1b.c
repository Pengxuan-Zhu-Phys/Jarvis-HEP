#include <setjmp.h>
#include "lanhep.h"

extern jmp_buf alg1_jmp_buf;

int opSetGpm=1;

Term SplitIndices(Term t, List *ilist)
	{
	Term rt,tt,t1=0;
	rt=ConsumeCompoundArg(t,1);
	tt=ConsumeCompoundArg(t,2);	
	FreeAtomic(t);
	while(is_compound(tt))
		{
		t1=CompoundName(tt);
		if(t1!=OPR_CARET && t1!=OPR_USCORE)
			{
			if(!IsTermInput())
				printf("File \"%s\", line %d: ",CurrentInputFile(),
					CurrentInputLine());
			printf("Semantic error: \'");
			WriteTerm(tt);
			printf("\' is not appropriate index.\n");
			FreeAtomic(tt); return 0;
			}		
		t1=ConsumeCompoundArg(tt,1);
		if(!is_atom(t1) && !is_integer(t1))
			{
			if(!IsTermInput())
				printf("File \"%s\", line %d: ",CurrentInputFile(),
					CurrentInputLine());
			printf("Semantic error: \'");
			WriteTerm(t1);
			printf("\' is not appropriate index.\n");
			FreeAtomic(tt);FreeAtomic(t1); return 0;
			}
		*ilist=AppendLast(*ilist,t1);
		t1=ConsumeCompoundArg(tt,2);
		FreeAtomic(tt);
		tt=t1;
		}
	if(!is_atom(tt) && !is_integer(tt))
		{
		if(!IsTermInput())
			printf("File \"%s\", line %d: ",CurrentInputFile(),
				CurrentInputLine());
		printf("Semantic error: \'");
		WriteTerm(tt);
		printf("\' is not appropriate index.\n");
		FreeAtomic(tt); return 0;
		}
	*ilist=AppendLast(*ilist,tt);
	return rt;
	}


static int equal_repres(Term f1, Term f2)
	{
	if(CompoundName(f1)==CompoundName(f2) &&
		(EqualTerms(CompoundArg1(f1),CompoundArg1(f2))
		|| EqualTerms(CompoundArg1(f1),CompoundArg2(f2))) )
				return 1;
	return 0;
	}


static int inst_in(Term a, List ind)
	{
	int i=0;
	while(!is_empty_list(ind))
		{
		Term il;
		il=CompoundArg1(ListFirst(ind));
		if(equal_repres(il,a))
			i++;
		ind=ListTail(ind);
		}
	return i;
	}
	
static int must_skip(int realno, List ind)
	{
	int res,mustbe;
	List l1;
	mustbe=ListLength(ind);
	if(realno==mustbe)
		return 0;
	if(realno>mustbe)
		return -1;
	l1=DefaultIndex;
	res=1;
	while(!is_empty_list(l1))
		{
		mustbe-=inst_in(ListFirst(l1),ind);
		if(realno==mustbe)
			return res;
		if(realno>mustbe)
			return -1;
		l1=ListTail(l1);
		res++;
		}
	return -1;
	}


static int should_skip(int no, Term a)
	{
	List l;
	int i=1;
	l=DefaultIndex;
	while(!is_empty_list(l))
		{
		Term gg;
		gg=ListFirst(l);
		if(equal_repres(gg,a))
			return 1;
		if(i==no)
			return 0;
		i++;
		l=ListTail(l);
		}
	return 0;
	}

static Term cc_particle(Term t1, List *ind)
	{
	Term t, prt;
	t=ConsumeCompoundArg(t1,1);
	FreeAtomic(t1);
	prt=GetAtomProperty(t,PROP_TYPE);
    	if( !(is_compound(prt) && CompoundName(prt)==OPR_PARTICLE))
		{
		ErrorInfo(216);
		printf(" cc(\'");WriteTerm(t);printf("\') is undefined.\n");
		longjmp(alg1_jmp_buf,1);
        }
	t1=t;
/*	if(CompoundArg1(prt)==t)
		t=CompoundArg2(prt);
	else
		t=CompoundArg1(prt);*/
	if(ind!=NULL)
		*ind=CopyTerm(GetAtomProperty(t,PROP_INDEX));
	if(!is_empty_list(*ind) && CompoundName(CompoundArg1(ListFirst(*ind)))==A_LORENTZ)
		{
		Term tt, in1,in2;
		tt=CompoundArg1(ListFirst(*ind));
		in1=ConsumeCompoundArg(tt,1);
		in2=ConsumeCompoundArg(tt,2);
		SetCompoundArg(tt,1,in2);
		SetCompoundArg(tt,2,in1);
		}
	/*
	WriteTerm(*ind); puts("");
	*/
	return t1;
	}

	
Term AtomicTo1(Term t, Term ind)
	{
	List rl;
	Term ret;
	Term t1;
	Atom ttype;
	
	if(is_compound(t) && CompoundName(t)==A_ALG1)
	{
		List l;
		ttype=A_ALG1;
		rl=CopyTerm(CompoundArg2(t));
		for(l=rl;l;l=ListTail(l))
			SetCompoundArg(ListFirst(l),2,0);
		goto cnt;
	}
	
	if(is_compound(t) && CompoundName(t)==A_FBRACET)
		{
		t=alg1_mk_wild(t,&rl,&ind);
		ttype=OPR_WILD;
		goto cnt;
		}
	/*
	if(is_compound(t) && CompoundName(t)==A_CC)
		{
		t=cc_particle(t,&rl);
		ttype=OPR_FIELD;
		goto cnt;
		}
	*/
	
	if(is_let(t,&rl))
		{
		ttype=OPR_LET;
		goto cnt;
		}
	
	if(is_parameter(t) || (is_compound(t) && (CompoundName(t)==A_COS ||
						CompoundName(t)==A_SIN)) )
		{
		rl=NewList();
		ttype=OPR_PARAMETER;
		goto cnt;
		}
		
	if(is_particle(t,&rl))
		{
		ttype=OPR_FIELD;
		goto cnt;
		}
		
	if(is_special(t,&rl))
		{
		ttype=OPR_SPECIAL;
		goto cnt;
		}
	
	ErrorInfo(301);
	printf(" \'%s\' undefined object.\n",AtomValue(t));
	FreeAtomic(ind);
	longjmp(alg1_jmp_buf,1);
cnt:
	ret=MakeCompound(A_MTERM,4);
	SetCompoundArg(ret,1,NewInteger(1));
	SetCompoundArg(ret,2,NewInteger(1));
	t1=0;
	if(!is_empty_list(rl))
		{
		Term ttt;
		int skipped;
		rl=CopyTerm(rl);
		if(is_empty_list(ind))
			{
			SetCompoundArg(ret,3,AppendFirst(NewList(),MakeCompound2(ttype,rl,t)));
			return AppendFirst(NewList(),ret);
			}
		skipped=must_skip(ListLength(ind),rl);
		if(skipped==-1)
			{
			ErrorInfo(302);
			printf(" can not set indices ");
			WriteTerm(ind);
			printf(" to \'%s\'.\n",AtomValue(t));
			FreeAtomic(ind); FreeAtomic(rl);
			longjmp(alg1_jmp_buf,1);
			}
		
		t1=rl;	
		ttt=ind;
		while(!is_empty_list(rl))
			{
			if(skipped!=0 && should_skip(skipped,CompoundArg1(ListFirst(rl))))
				{ rl=ListTail(rl); continue; }
			SetCompoundArg(ListFirst(rl),2,ListFirst(ttt));
			rl=ListTail(rl);
			ttt=ListTail(ttt);
			}
		FreeAtomic(ind);
		}
		
	if(ttype==OPR_PARAMETER && is_compound(t))
	{
	Term mmm=MakeCompound(ttype,CompoundArity(t)+2);
	SetCompoundArg(mmm,1,t1);
	SetCompoundArg(mmm,2,CompoundName(t));
	SetCompoundArg(mmm,3,CompoundArg1(t));
	if(CompoundArity(t)==2)
		{
		Term m=ExprTo1(ConsumeCompoundArg(t,2));
		m=CompoundArg1(m);
		if(ListLength(m)!=1)
			{
			ErrorInfo(0);
			puts("bad expression in sin/cos.");
			longjmp(alg1_jmp_buf,1);
			}
		m=ListFirst(m);
		SetCompoundArg(mmm,4,m);
		}
	SetCompoundArg(ret,3,AppendFirst(NewList(),mmm));
	}
	else
	SetCompoundArg(ret,3,AppendFirst(NewList(),MakeCompound2(ttype,t1,t)));
	return AppendFirst(NewList(),ret);
	
	}
	
static void mult_no(List e1, int no)
	{
	while(!is_empty_list(e1))
		{
		long int v;
		v=IntegerValue(CompoundArg1(ListFirst(e1)));
		v*=no;
		SetCompoundArg(ListFirst(e1),1,NewInteger(v));
		e1=ListTail(e1);
		}
	}
	
		
	
static Term multiply(Term t1, Term t2)
	{
	Term ret;
	long int n1,n2,d1,d2,num,den,cf;
	List l1,l2;
	ret=MakeCompound(A_MTERM,4);
	n1=IntegerValue(CompoundArg1(t1));
	n2=IntegerValue(CompoundArg1(t2));
	d1=IntegerValue(CompoundArg2(t1));
	d2=IntegerValue(CompoundArg2(t2));
	num=n1*n2;
	den=d1*d2;
	if(den<0) den=-den;
	cf=gcf(num,den);
	num/=cf;
	den/=cf;
	if(d1<0 && d2<0)
		{
		num=-num;
		}
	else
		if((d1<0 && d2>0) || (d1>0 && d2<0))
			{
			den=-den;
			}
	SetCompoundArg(ret,1,NewInteger(num));
	SetCompoundArg(ret,2,NewInteger(den));
	l1=ConsumeCompoundArg(t1,3);
	l2=ConsumeCompoundArg(t2,3);
	SetCompoundArg(ret,3,ConcatList(l1,l2));
	l1=ConsumeCompoundArg(t1,4);
	l2=ConsumeCompoundArg(t2,4);
	SetCompoundArg(ret,4,ConcatList(l1,l2));
	FreeAtomic(t1); FreeAtomic(t2);
	return ret;
	}
	
static Term multiply_l(Term t1,Term t2)
	{
	List l1,l2,ll;
	ll=NewList();
	l1=t1;
	while(!is_empty_list(l1))
		{
		l2=t2;
		while(!is_empty_list(l2))
			{
			Term q1,q2;
			q1=CopyTerm(ListFirst(l1));
			q2=CopyTerm(ListFirst(l2));
			ll=AppendLast(ll,multiply(q1,q2));
			l2=ListTail(l2);
			}
		l1=ListTail(l1);
		}
	FreeAtomic(t1);
	FreeAtomic(t2);
	return ll;
	}
	
static void invert_term(Term t)
	{
	
		{
		long int i1,i2;
		i1=IntegerValue(CompoundArg1(t));
		i2=IntegerValue(CompoundArg2(t));
		if(i1<0 && i2>0)
			{
			i1=-i1;
			i2=-i2;
			goto ee;
			}
		if(i1>0 && i2<0)
			{
			i1=-i1;
			goto ee;
			}
		if(i1<0 && i2<0)
			{
			i2=-i2;
			goto ee;
			}
	ee:
		SetCompoundArg(t,1,NewInteger(i2));
		SetCompoundArg(t,2,NewInteger(i1));
		}
		
		{
		Term t1,t2;
		t1=ConsumeCompoundArg(t,3);
		t2=ConsumeCompoundArg(t,4);
		SetCompoundArg(t,3,t2);
		SetCompoundArg(t,4,t1);
		}
		
	}


Term ExpandTerm(Term t)
	{
	List il=0;
	Term res;
	/*puts("");
	WriteTerm(t);
	puts("\n");*/
	if(is_compound(t) && (CompoundName(t)==OPR_USCORE || CompoundName(t)==OPR_CARET))
		{
		t=SplitIndices(t,&il);
		if(t==0)
			longjmp(alg1_jmp_buf,1);
		}
	if(is_function(t,NULL))
		return ExpandTerm(CallFunction(t,il));
	if(is_float(t))
		{
		ErrorInfo(303);
		printf(" illegal use of floating point number %f.\n",FloatValue(t));
		FreeAtomic(t);
		longjmp(alg1_jmp_buf,1);
		}
	if(is_integer(t))
		{
		if(IntegerValue(t)==0)
			return 0;
		
		if(il)
			{
			ErrorInfo(304);
			puts("Integer can't have indices.\n");
			longjmp(alg1_jmp_buf,1);
			}
			
		res=MakeCompound(A_MTERM,4);
		SetCompoundArg(res,1,t);
		SetCompoundArg(res,2,NewInteger(1));
		return AppendFirst(0,res);
		}
		
	if(t==A_I)
		{
		res=MakeCompound(A_MTERM,4);
		SetCompoundArg(res,1,NewInteger(1));
		SetCompoundArg(res,2,NewInteger(-1));
		return AppendFirst(0,res);
		}
		
	if(is_atom(t) || (is_compound(t) && CompoundName(t)==A_FBRACET) ||
		(is_compound(t) && CompoundName(t)==A_ALG1 && CompoundArity(t)==2) ||
		(is_compound(t) && (CompoundName(t)==A_SIN || CompoundName(t)==A_COS)
		 	&& CompoundArity(t)<3 && is_integer(CompoundArg1(t))))
		 return AtomicTo1(t,il);
	
	if(!is_empty_list(il))
		{
		ErrorInfo(307);
		printf(" bad indices ");
		WriteTerm(il);
		printf(" in expression.\n");
		longjmp(alg1_jmp_buf,1);
		}


	if(opSetGpm && is_compound(t) && CompoundArity(t)==2 && CompoundArg1(t)==NewInteger(1)
			&& is_atom(CompoundArg2(t))
			&& (CompoundArg2(t)==A_GAMMA5 || GetAtomProperty(CompoundArg2(t),A_GAMMA5))
			&& (CompoundName(t)==OPR_PLUS || CompoundName(t)==OPR_MINUS))
	{
		Term ret,sp;
		ret=MakeCompound(A_MTERM,4);
		SetCompoundArg(ret,1,NewInteger(2));
		SetCompoundArg(ret,2,NewInteger(1));
		sp=MakeCompound(OPR_SPECIAL,2);
		if(CompoundName(t)==OPR_PLUS)
			SetCompoundArg(sp,2,A_GAMMAP);
		else
			SetCompoundArg(sp,2,A_GAMMAM);
		SetCompoundArg(sp,1,CopyTerm(GetAtomProperty(A_GAMMA5,PROP_INDEX)));
		SetCompoundArg(ret,3,AppendLast(NewList(),sp));
		return AppendLast(NewList(),ret);
	}
		

	
	
	if(is_compound(t) && CompoundName(t)==OPR_PLUS && CompoundArity(t)==2)
		{
		Term t1,t2;
		t1=ConsumeCompoundArg(t,1);
		t2=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		t1=ExpandTerm(t1);
		t2=ExpandTerm(t2);
		return ConcatList(t1,t2);
		}
		
	if(is_compound(t) && CompoundName(t)==OPR_MINUS && CompoundArity(t)==2)
		{
		Term t1,t2;
		t1=ConsumeCompoundArg(t,1);
		t2=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		t1=ExpandTerm(t1);
		t2=ExpandTerm(t2);
		mult_no(t2,-1);
		return ConcatList(t1,t2);
		}
		
	if(is_compound(t) && CompoundName(t)==OPR_MINUS && CompoundArity(t)==1)
		{
		Term t1;
		t1=ConsumeCompoundArg(t,1);
		FreeAtomic(t);
		t1=ExpandTerm(t1);
		mult_no(t1,-1);
		return t1;
		}
	
	if(is_compound(t) && CompoundName(t)==OPR_MLT && CompoundArity(t)==2)
		{
		Term t1,t2;
		t1=ConsumeCompoundArg(t,1);
		t2=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		t1=ExpandTerm(t1);
		t2=ExpandTerm(t2);
		return multiply_l(t1,t2);
		}
		
	if(is_compound(t) && CompoundName(t)==OPR_DIV && CompoundArity(t)==2)
		{
		Term t1,t2,t2s;
		t1=ConsumeCompoundArg(t,1);
		t2=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		t1=ExpandTerm(t1);
		t2s=CopyTerm(t2);
		t2=ExpandTerm(t2);
		if(ListLength(t2)!=1)
			{
			ErrorInfo(306);
			printf(" cannot divide by \'");
			WriteTerm(t2s); printf("\'.\n");
			FreeAtomic(t1); FreeAtomic(t2); FreeAtomic(t2s);
			longjmp(alg1_jmp_buf,1);
			}
		FreeAtomic(t2s);
		invert_term(ListFirst(t2));
		return multiply_l(t1,t2);
		}
		
	if(is_compound(t) && CompoundName(t)==OPR_POW && CompoundArity(t)==2)
		{
		Term t1,t2,ret;
		int i,pp;
		t1=ConsumeCompoundArg(t,1);
		t2=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		if(!is_integer(t2) || IntegerValue(t2)<1)
			{
			ErrorInfo(307);
			printf(" illegal power \'");
			WriteTerm(t2); printf("\'.\n");
			FreeAtomic(t1); FreeAtomic(t2);
			longjmp(alg1_jmp_buf,1);
			}
		pp=(int)IntegerValue(t2);
		t1=ExpandTerm(t1);
		ret=CopyTerm(t1);
		for(i=2;i<=pp;i++)
			{
			Term tmp;
			tmp=CopyTerm(t1);
			ret=multiply_l(ret,tmp);
			}
		FreeAtomic(t1);
		return ret;
		}
	
	ErrorInfo(309);
	printf("bad expression ");
	WriteTerm(t);
	puts("");
	longjmp(alg1_jmp_buf,1);	
	
	}
