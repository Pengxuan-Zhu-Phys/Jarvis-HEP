#include "lanhep.h"
#include <string.h>

static int allow_transf_lets=1;

int block_transf=0;

Term ExprTo1kl(Term);

static void c_mlt(long int n1, long int d1, long int n2, long int d2, long int *n, long int *d)
	{
	*n=n1*n2;
	*d=d1*d2;
	if(d1<0 && d2<0)
		(*n)=-(*n);
	}	 

static int only_prm(Term a1)
{
	List l1,l2;
	for(l1=CompoundArg1(a1);l1;l1=ListTail(l1))
	{
		if(IntegerValue(CompoundArg2(ListFirst(l1)))<0)
			return 0;
		for(l2=CompoundArgN(ListFirst(l1),3);l2;l2=ListTail(l2))
		{

			if(CompoundName(ListFirst(l2))==OPR_PARAMETER)
				continue;
			if(CompoundName(ListFirst(l2))==OPR_WILD)
			{
				List l3;
				for(l3=CompoundArg2(ListFirst(l2));l3;l3=ListTail(l3))
					if(!only_prm(ListFirst(l3)))
						return 0;
				continue;
			}
			return 0;
		}
	}
	
	return 1;
}

	
Term alg1_inv_alg(Term a1)
	{
	List l, l1, lm;
	Term m2_b, m2_b2;
	List m2_i;
	l=CompoundArg1(a1);
	if(CompoundArg2(a1))
		return 0;
	
	m2_i=m2_b=0;
	
	for(l1=l;l1;l1=ListTail(l1))
		{
		Term m1;
		List l2;
		int  ino;
		ino=0;
		m1=ListFirst(l1);
		for(l2=CompoundArgN(m1,3);l2;l2=ListTail(l2))
			{
			Term prp=0;
			if(CompoundName(ListFirst(l2))==OPR_PARAMETER)
				prp=GetAtomProperty(CompoundArg2(ListFirst(l2)),A_INFINITESIMAL);
			else
				{
				RemoveList(m2_i);
				return 0;
				}
			if(prp && IntegerValue(CompoundArg1(prp))>0)
				ino+=(int)IntegerValue(CompoundArg1(prp));
			}
		if(!ino)
			{
			if(m2_b)
				{
				RemoveList(m2_i);
				return 0;
				}
			m2_b=m1;
			}
		else
			m2_i=AppendLast(m2_i,m1);
		}
	if(m2_b==0)
		return 0;
	m2_b=CopyTerm(m2_b);
	l=ConsumeCompoundArg(m2_b,3);
	l1=ConsumeCompoundArg(m2_b,4);
	SetCompoundArg(m2_b,3,l1);
	SetCompoundArg(m2_b,4,l);
		
		{
		long int n, ns, d, ds;
		n=IntegerValue(CompoundArg1(m2_b));
		d=IntegerValue(CompoundArg2(m2_b));
		ns=ds=1;
		if(n<0)
			{
			ns=-1;
			n=-n;
			}
		if(d<0)
			{
			ds=-1;
			d=-d;
			}
		SetCompoundArg(m2_b,1,NewInteger(d*ns));
		SetCompoundArg(m2_b,2,NewInteger(n*ds));
		}
	
	lm=AppendLast(NewList(),m2_b);
	if(m2_i)
		{
		m2_b2=CopyTerm(m2_b);
		l=ConsumeCompoundArg(m2_b2,3);
		l=ConcatList(l,CopyTerm(l));
		SetCompoundArg(m2_b2,3,l);
		l=ConsumeCompoundArg(m2_b2,4);
		l=ConcatList(l,CopyTerm(l));
		SetCompoundArg(m2_b2,4,l);
			{
			long int n,d,n1,d1;
			n=IntegerValue(CompoundArg1(m2_b2));
			d=IntegerValue(CompoundArg2(m2_b2));
			c_mlt(n,d,-n,d,&n1,&d1);
			SetCompoundArg(m2_b2,1,NewInteger(n1));
			SetCompoundArg(m2_b2,2,NewInteger(d1));
			}
		for(l=m2_i;l;l=ListTail(l))
			{
			Term a1;
			a1=MakeCompound(A_MTERM,4);
			l1=ConcatList(CopyTerm(CompoundArgN(m2_b2,3)),
						  CopyTerm(CompoundArgN(ListFirst(l),3)));
			SetCompoundArg(a1,3,l1);
			l1=ConcatList(CopyTerm(CompoundArgN(m2_b2,4)),
						  CopyTerm(CompoundArgN(ListFirst(l),4)));
			SetCompoundArg(a1,4,l1);
				{
				long int n,d;
				c_mlt(	IntegerValue(CompoundArg1(m2_b2)),
						IntegerValue(CompoundArg2(m2_b2)),
						IntegerValue(CompoundArg1(ListFirst(l))),
						IntegerValue(CompoundArg2(ListFirst(l))),
						&n,&d);
				SetCompoundArg(a1,1,NewInteger(n));
				SetCompoundArg(a1,2,NewInteger(d));
				}
			lm=AppendLast(lm,a1);
			
			}
		RemoveList(m2_i);
		FreeAtomic(m2_b2);
		}
/*    DumpList(lm);*/	
	return MakeCompound2(A_ALG1,lm,0);		
	
	}	
	

static void alg1_set_cos0(Term a1)
	{
	List l;
	Term c0;
	int scno=0;
	for(l=CompoundArg1(a1);l;l=ListTail(l))
		{
		List l2;
		for(l2=CompoundArgN(ListFirst(l),3);l2;l2=ListTail(l2))
			{
			Term t=ListFirst(l2);
			if(CompoundName(t)==OPR_PARAMETER && 
			 (CompoundArg2(t)==A_SIN || CompoundArg2(t)==A_COS))
			 	scno++;
			if(CompoundName(t)==OPR_WILD)
				{
				List l3;
				for(l3=CompoundArg2(t);l3;l3=ListTail(l3))
					alg1_set_cos0(ListFirst(l3));
				}
			}
		}
	if(scno==0)
		return;
	
	c0=MakeCompound3(OPR_PARAMETER,0,A_COS,NewInteger(0));
	
	for(l=CompoundArg1(a1);l;l=ListTail(l))
		{
		List l2;
		scno=0;
		for(l2=CompoundArgN(ListFirst(l),3);l2;l2=ListTail(l2))
			{
			Term t=ListFirst(l2);
			if(CompoundName(t)==OPR_PARAMETER && 
			 (CompoundArg2(t)==A_SIN || CompoundArg2(t)==A_COS))
			 	scno++;
			}
		if(scno)
			continue;
		l2=ConsumeCompoundArg(ListFirst(l),3);
		l2=AppendLast(l2,CopyTerm(c0));
		SetCompoundArg(ListFirst(l),3,l2);
		}
	FreeAtomic(c0);
	/*WriteTerm(a1);puts("\n");*/
	
	}
	
extern void alg1_inf_wild(Term);

Term ProcLet(Term t, Term ind)
	{
	Term t1,nm,sub,kl=0;
	List il,ol,l1;	
	int transf_fl=0;
	Atom anti1=0, anti2=0;
	
	ol=il=NewList();
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	
	
	if(is_compound(t1) && CompoundArity(t1)==2 && CompoundName(t1)==OPR_COMMA)
		{
		Term a1,a2;
		a1=ConsumeCompoundArg(t1,1);
		a2=ConsumeCompoundArg(t1,2);
		FreeAtomic(t1);
		ProcLet(MakeCompound1(OPR_LET,a1),0);
		ProcLet(MakeCompound1(OPR_LET,a2),0);
		return 0;
		}
		
		
	if(!is_compound(t1) || CompoundArity(t1)!=2 || 
	 (CompoundName(t1)!=OPR_EQSIGN && CompoundName(t1)!=OPR_RARROW))
		{
		ErrorInfo(325);
		printf("bad argument in let statement\n");
		return 0;
		}
		
	if(CompoundName(t1)==OPR_EQSIGN && is_atom(CompoundArg1(t1)) &&
		is_atom(CompoundArg2(t1)) &&
		(CompoundArg2(t1)==A_GAMMA5 || GetAtomProperty(CompoundArg2(t1),
		A_GAMMA5)))
		{
		SetAtomProperty(CompoundArg1(t1),A_GAMMA5,NewInteger(1));
		}
		
	if(CompoundName(t1)==OPR_RARROW)
		transf_fl=1;
	
	if(CompoundName(t1)==OPR_EQSIGN && is_atom(CompoundArg1(t1)) &&
			is_compound(CompoundArg2(t1)) && 
			CompoundName(CompoundArg2(t1))==A_ANTI &&
			is_atom(CompoundArg1(CompoundArg2(t1))))
		{
		anti1=CompoundArg1(t1);
		anti2=CompoundArg1(CompoundArg2(t1));
		}
		
	nm=ConsumeCompoundArg(t1,1);
	sub=ConsumeCompoundArg(t1,2);
	FreeAtomic(t1);
	
	if(transf_fl)
		allow_transf_lets=0;
	
	if(is_compound(nm) && (CompoundName(nm)==OPR_USCORE ||
								 CompoundName(nm)==OPR_CARET))
		nm=SplitIndices(nm,&il);
	
	if(!is_atom(nm))
		{
		ErrorInfo(325);
		printf("bad left argument in let call\n");
		FreeAtomic(sub);
		FreeAtomic(t1);
		return 0;
		}
		
	if(transf_fl)
		{
		if(!is_parameter(nm) && !is_particle(nm,NULL))
			{
			ErrorInfo(728);
			printf("Unknown object '%s'.\n",AtomValue(nm));
			return 0;
			}
		}
	
	
	if(GetAtomProperty(nm,A_KEEP_LETS))
	{
		Term prp;
		kl=ExprTo1kl(sub);
		if(kl==0)
			return 0;
		sub=CopyTerm(kl);
		prp=GetAtomProperty(nm,A_KEEP_LETS);
		SetCompoundArg(prp,1,kl);
	}
	else
	{
	
	sub=ExprTo1(sub);
	if(sub==0)
		return 0;
	alg1_set_cos0(sub);
	alg1_inf_wild(sub);
	}
	
	if(sub==0)
		return 0;
	if(transf_fl)
		allow_transf_lets=1;
			
	t1=CopyTerm(CompoundArg2(sub));

	
	if(is_empty_list(il))
		{
		l1=t1;
		while(!is_empty_list(l1))
			{
			Term t;
			t=CompoundArg2(ListFirst(l1));
			if(!is_label(t))
				{
				ErrorInfo(326);
				printf("unbalanced index '");
				WriteTerm(t);
				printf("' in let statement\n");
				FreeAtomic(sub);
				FreeAtomic(t1);
				FreeAtomic(ol);
				return 0;
				}
			ol=AppendLast(ol,t);
			l1=ListTail(l1);
			}
		/* mk_let(nm,sub,ol,t1); */
		}
	else
		{
		List t2;
		t2=0;
		if(ListLength(il)!=ListLength(t1))
				{
				ErrorInfo(327);
				printf("distinct indices number ");
				printf("in let statement.\n");
				FreeAtomic(sub);
				FreeAtomic(t1);
				return 0;
				}			
		l1=t1;
		while(!is_empty_list(l1))
			{
			Term t;
			t=CompoundArg2(ListFirst(l1));
			if(!is_atom(t) || !ListMember(il,t))
				{
				ErrorInfo(326);
				printf("unbalanced index '");
				WriteTerm(t);
				printf("' in let statement\n");
				FreeAtomic(sub);
				FreeAtomic(t1);
				FreeAtomic(ol);
				return 0;
				}
			ol=AppendLast(ol,t);
			l1=ListTail(l1);
			}
		l1=il;
		while(!is_empty_list(l1))
			{
			List l2;
			if(!ListMember(ol,ListFirst(l1)))
				{
				ErrorInfo(326);
				printf("unbalanced index '");
				WriteTerm(t);
				printf("' in let statement\n");
				FreeAtomic(sub);
				FreeAtomic(t1);
				FreeAtomic(ol);
				return 0;
				}
				
			for(l2=t1;l2;l2=ListTail(l2))
				if(ListFirst(l1)==CompoundArg2(ListFirst(l2)))
					{
					t2=AppendLast(t2,ListFirst(l2));
					break;
					}
				
			l1=ListTail(l1);
			}
		RemoveList(t1);
		t1=t2;
		FreeAtomic(ol);
		ol=il;		
		}	
	
	/* mk_let(nm,sub,ol,t1); */
	
	
	l1=t1;
	while(!is_empty_list(l1))
		{
		SetCompoundArg(ListFirst(l1),2,0);
		l1=ListTail(l1);
		}
		
	/*WriteTerm(sub); puts(""); getchar();*/

		
	t=MakeCompound(OPR_LET,5);
	SetCompoundArg(t,1,sub);
	SetCompoundArg(t,2,ol);
	SetCompoundArg(t,3,alg1_inv_alg(sub));
	
	
	if(transf_fl==0)
		{
		Term tt;
		tt=alg1_guess_mpl(sub);
		if(tt)
		{
			SetCompoundArg(t,4,NewInteger(1));
			SetCompoundArg(t,5,tt);
		}
		else
		{
			int tp;
			tt=alg1_guess_mtr(sub, &tp);
			if(tt)
			{
				SetCompoundArg(t,4,NewInteger(tp));
				SetCompoundArg(t,5,tt);
			}
		}
		
		ReportRedefined(nm,"let-substitution");
		SetAtomProperty(nm,PROP_INDEX,t1);
		SetAtomProperty(nm,PROP_TYPE,t);
		alg1_let_cw(nm);
		if(anti1 && anti2)
			{
			SetAtomProperty(anti1,A_ANTI,anti2);
			SetAtomProperty(anti2,A_ANTI,anti1);
			}
		else
			if(only_prm(CompoundArg1(t)))
				SetAtomProperty(nm,A_ANTI,nm);
		return 0;
		}
		
	l1=GetAtomProperty(nm,PROP_INDEX);
	if(!EqualTerms(l1,t1))
		{
		ErrorInfo(729);
		puts("transformed object has other indices types");
		return 0;
		}
	if(GetAtomProperty(nm,OPR_LET))
		{
		WarningInfo(0);
		printf("Warning: transformation rule for '%s' is redefined.\n",
			AtomValue(nm));
		}
	SetAtomProperty(nm,OPR_LET,t);
	return 0;
	}



int is_let(Term t, List *l)
	{
	Term a;
	if(!block_transf && allow_transf_lets && GetAtomProperty(t,OPR_LET))
		{
		if(l!=NULL)
			*l=GetAtomProperty(t,PROP_INDEX);
		return 1;
		}
	a=GetAtomProperty(t,PROP_TYPE);
	if(!is_compound(a) || CompoundName(a)!=OPR_LET)
		return 0;
	if(l!=NULL)
		*l=GetAtomProperty(t,PROP_INDEX);
	return 1;
	}
	
void ClearLet(Atom ll)
	{
	RemoveAtomProperty(ll,PROP_TYPE);
	RemoveAtomProperty(ll,PROP_INDEX);
	}
		


Term ProcKeepLets(Term t, Term ind)
{
	Term t1;
	List l;
	
	if(!is_compound(t) || CompoundArity(t)!=1)
	{
		ErrorInfo(331);
		printf("bad syntax in 'keep_lets' statement.\n");
		return 0;
	}
	
	t1=CommaToList(ConsumeCompoundArg(t,1));
	FreeAtomic(t);
	
	for(l=t1;l;l=ListTail(l))
	{
		Atom p;
		p=ListFirst(l);
		if(!is_atom(p))
		{
			ErrorInfo(332);
			printf("unexpected '");
			WriteTerm(p);
			printf("' in 'keep_lets' statement.\n");
			continue;
		}
		if(GetAtomProperty(p,PROP_TYPE))
		{
			ErrorInfo(333);
			printf("'keep_lets': object '%s' is already defined.\n",
					AtomValue(p));
			continue;
		}
		SetAtomProperty(p,A_KEEP_LETS,MakeCompound1(A_KEEP_LETS,0));
	}
	
	return 0;
}

int prmcmp(Term p1, Term p2)
{
	return strcmp(AtomValue(CompoundArg2(p1)),AtomValue(CompoundArg2(p2)));
}

void alg1_opt_let(Term a1)
{
	List l,l1,l2,lm=ConsumeCompoundArg(a1,1);
	List nl=0;
	Label mylbl=NewLabel();
	
	/*printf("%d -> ",ListLength(lm));*/
	
	for(l=lm;l;l=ListTail(l))
	{
		Term m1=ListFirst(l);
		long int n,d,g;
		List ln=ConsumeCompoundArg(m1,3);
		List ld=ConsumeCompoundArg(m1,4);
		Term f1=0, f2=0, w=0;
		n=IntegerValue(CompoundArg1(m1));
		d=IntegerValue(CompoundArg2(m1));
				
	rpt1:
		for(l1=ld;l1;l1=ListTail(l1))
		{
			int n=ListMember(ln,ListFirst(l1));
			if(n)
			{
				ld=CutFromList(ld,l1);
				l1=ListNthList(ln,n);
				ln=CutFromList(ln,l1);
				goto rpt1;
			}
			if(CompoundArg2(ListFirst(l1))==A_SQRT2)
			{
				Term t=ListFirst(l1);
				ChangeList(l1,0);
				ld=CutFromList(ld,l1);
				ln=AppendFirst(ln,t);
				d*=2;
				goto rpt1;
			}
		}
	rpt2:
		for(l1=ln;l1;l1=ListTail(l1))
			if(CompoundArg2(ListFirst(l1))==A_SQRT2)
			{
				List l2;
				for(l2=ListTail(l1);l2;l2=ListTail(l2))
					if(CompoundArg2(ListFirst(l2))==A_SQRT2)
					{
						ln=CutFromList(ln,l1);
						ln=CutFromList(ln,l2);
						n*=2;
						goto rpt2;
					}
				}
		g=gcf(n,d);
		n/=g;
		d/=g;
		
		for(l1=ln;l1;l1=l1?ListTail(l1):0)
		{
			Term sp=ListFirst(l1);
			if(CompoundName(sp)==OPR_PARAMETER)
				continue;
			if(CompoundName(sp)==OPR_FIELD)
			{
				if(f1==0)
					f1=l1;
				else if(f2==0)
					f2=l1;
				else
				{
					ErrorInfo(1280);
					puts("bad let-sub: >2 prtc");
					return;
				}
				continue;
			}
			if(CompoundName(sp)==OPR_WILD)
			{
				if(w)
				{
					ErrorInfo(1281);
					puts("bad let-sub: unk spec");
					return;
				}
				w=l1;
				continue;
			}
			ErrorInfo(110);
			WriteTerm(sp);puts(" : bad let: unk stuff");
			return;
		}
		
		if(f2==0)
		{
			ErrorInfo(112);
			puts("bad let-subs: not 2 prtc");
			return;
		}
		
		l1=f1;f1=ListFirst(f1);ChangeList(l1,0);ln=CutFromList(ln,l1);
		l1=f2;f2=ListFirst(f2);ChangeList(l1,0);ln=CutFromList(ln,l1);
		if(w)
		{l1=w;w=ListFirst(w);ChangeList(l1,0);ln=CutFromList(ln,l1);}
		
		ln=SortedList(ln,prmcmp);
		ld=SortedList(ld,prmcmp);
		
		if((GetAtomProperty(CompoundArg2(f1),A_ANTI)==CompoundArg2(f1) ||
				GetAtomProperty(CompoundArg2(f2),A_ANTI)==CompoundArg2(f2))
			&& strcmp(AtomValue(CompoundArg2(f1)),AtomValue(CompoundArg2(f2)))>0)
		{
			Term tmp=f1;
			f1=f2;
			f2=tmp;
		}
		
		if(ListLength(CompoundArg1(f1))==1 && ListLength(CompoundArg1(f1))==1)
		{
			Term i1=ListFirst(CompoundArg1(f1));
			Term i2=ListFirst(CompoundArg1(f2));
			if(CompoundName(CompoundArg1(i1))!=A_COLOR ||
				CompoundName(CompoundArg1(i2))!=A_COLOR ||
					CompoundArg2(i1)!=CompoundArg2(i2))
			{
				ErrorInfo(233);
				puts("bad color stru in let-sub.");
				return;
			}
			SetCompoundArg(i1,2,mylbl);
			SetCompoundArg(i2,2,mylbl);
		}
		
		ln=AppendLast(ln,f1);
		ln=AppendLast(ln,f2);
		if(w) ln=AppendLast(ln,w);
		
		SetCompoundArg(m1,1,NewInteger(n));
		SetCompoundArg(m1,2,NewInteger(d));
		SetCompoundArg(m1,3,ln);
		SetCompoundArg(m1,4,ld);
		
		for(l1=nl;l1;l1=ListTail(l1))
		{
			if(EqualTerms(CompoundArgN(ListFirst(l1),3),ln) &&
					EqualTerms(CompoundArgN(ListFirst(l1),4),ld) &&
					IntegerValue(CompoundArg2(ListFirst(l1)))*d>0)
			{
				Term om1=ListFirst(l1);
				long int n1=IntegerValue(CompoundArg1(om1));
				long int d1=IntegerValue(CompoundArg2(om1));
				long int rn,rd,cc=0;
				if(d1<0)
				{
					d=-d;d1=-d1;cc=1;
				}
				rn=n*d1+n1*d;
				rd=d*d1;
				if(rn==0)
				{
					nl=CutFromList(nl,l1);
					break;
				}
				g=gcf(rn,rd);
				rn/=g; rd/=g;
				if(cc) rd=-rd;
				SetCompoundArg(om1,1,NewInteger(rn));
				SetCompoundArg(om1,2,NewInteger(rd));
				FreeAtomic(m1);
				break;
			}
		}
		
		if(l1==0)
			nl=AppendLast(nl,m1);
		
	}
	
	RemoveList(lm);
	/*DumpList(nl);*/
	
	/*printf("%d\n",ListLength(nl));*/
	SetCompoundArg(a1,1,nl);
}

Term ProcOptLets(Term ls, Term ind)
{
	List l=ConsumeCompoundArg(ls,1);
	/*WriteTerm(l);puts("");*/
	FreeAtomic(ls);
	ls=CommaToList(l);
	/*WriteTerm(ls);puts("");*/
	for(l=ls;l;l=ListTail(l))
	{
		Term prp=GetAtomProperty(ListFirst(l),PROP_TYPE);
		if(!prp || CompoundName(prp)!=OPR_LET)
		{
			ErrorInfo(900);
			WriteTerm(ListFirst(l));
			puts(" is not a let subst.");
		}
		
		alg1_opt_let(CompoundArg1(prp));
		
	}
	return 0;
}
	
