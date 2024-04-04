#include <setjmp.h>
#include "lanhep.h"

extern jmp_buf alg1_jmp_buf;


static void repl_ind(Term a2, List oi, List ni)
	{

	a2=CompoundArg1(a2);
	while(!is_empty_list(a2))
		{
		List t;
		t=CompoundArgN(ListFirst(a2),3);
		while(!is_empty_list(t))
			{
			List t1;
			t1=CompoundArg1(ListFirst(t));

			while(!is_empty_list(t1))
				{
				Term a;
				List loi,lni;
				a=CompoundArg2(ListFirst(t1));
				loi=oi;
				lni=ni;
				while(!is_empty_list(loi))
					{
					if(a==CompoundArg2(ListFirst(loi)))
						{
						SetCompoundArg(ListFirst(t1),2,
								CompoundArg2(ListFirst(lni)));
						break;
						}
					loi=ListTail(loi);
					lni=ListTail(lni);
					}
				t1=ListTail(t1);
				}

			t=ListTail(t);
			}
		a2=ListTail(a2);
		}				
	return ;
	}			
		

static List mk_let(Term m1, List cut, Term a1)
	{
	List l;
	List lb,le;
	List sd;
	long int num1,den1;
	
	l=ConsumeCompoundArg(a1,1);
	FreeAtomic(a1);
	a1=l;
	
	
	num1=IntegerValue(ConsumeCompoundArg(m1,1));
	den1=IntegerValue(ConsumeCompoundArg(m1,2));
	l=ConsumeCompoundArg(m1,3);
	sd=ConsumeCompoundArg(m1,4);
	lb=ListSplit(l,cut,&le);
	FreeAtomic(cut);
	FreeAtomic(m1);
	
	l=a1;
	while(!is_empty_list(l))
		{
		long int n1,n2,d1,d2,num,den,cf;
		List lb1,le1,lm;
		m1=ListFirst(l);
		lm=ConsumeCompoundArg(m1,3);
		lb1=CopyTerm(lb);
		le1=CopyTerm(le);
		lb1=ConcatList(lb1,lm);
		lb1=ConcatList(lb1,le1);
		SetCompoundArg(m1,3,lb1);
		
		lm=ConsumeCompoundArg(m1,4);
		lm=ConcatList(lm,CopyTerm(sd));
		SetCompoundArg(m1,4,lm);
		
		n1=num1;
		d1=den1;
		n2=IntegerValue(CompoundArg1(m1));
		d2=IntegerValue(CompoundArg2(m1));
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
		SetCompoundArg(m1,1,NewInteger(num));
		SetCompoundArg(m1,2,NewInteger(den));
		l=ListTail(l);
		}
	FreeAtomic(lb);
	FreeAtomic(le);
	FreeAtomic(sd);
	return a1;
	}



static List mk_let_d(Term m1,  Term a1)
	{
	List l;
	List lb;
	List sd;
	long int num1,den1;


	l=ConsumeCompoundArg(a1,1);
	FreeAtomic(a1);
	a1=l;
	
	num1=IntegerValue(ConsumeCompoundArg(m1,1));
	den1=IntegerValue(ConsumeCompoundArg(m1,2));
	l=ConsumeCompoundArg(m1,3);
	sd=ConsumeCompoundArg(m1,4);
	lb=l;
	FreeAtomic(m1);
	
	l=a1;
	while(!is_empty_list(l))
		{
		long int n1,n2,d1,d2,num,den,cf;
		List lb1,lm;
		m1=ListFirst(l);
		lm=ConsumeCompoundArg(m1,3);
		lb1=CopyTerm(lb);
		lb1=ConcatList(lb1,lm);
		SetCompoundArg(m1,3,lb1);
		
		lm=ConsumeCompoundArg(m1,4);
		lm=ConcatList(lm,CopyTerm(sd));
		SetCompoundArg(m1,4,lm);
		
		n1=num1;
		d1=den1;
		n2=IntegerValue(CompoundArg1(m1));
		d2=IntegerValue(CompoundArg2(m1));
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
		SetCompoundArg(m1,1,NewInteger(num));
		SetCompoundArg(m1,2,NewInteger(den));
		l=ListTail(l);
		}
	FreeAtomic(lb);
	FreeAtomic(sd);
	return a1;
	}


static List s_l_1(Term m1)
	{
	List l,l1;
	
	l=CompoundArgN(m1,3);

	while(!is_empty_list(l))
		{
		Term t1;
		t1=ListFirst(l);
		if(CompoundName(t1)==A_ALG1)
			{
			Term sub,a1,ila,ill;
			
			a1=sub=ConsumeCompoundArg(t1,2);
			ill=ConsumeCompoundArg(t1,1);
			FreeAtomic(t1);
			ChangeList(l,0);
					
			ila=CopyTerm(CompoundArg2(sub));
			
			repl_ind(a1,ila,ill);

			return SetIntAlgs(mk_let(m1,l,a1));
			
			}
		l=ListTail(l);
		}
	
	
	l=CompoundArgN(m1,4);
	while(!is_empty_list(l))
		{
		Term t1;
		t1=ListFirst(l);
		if(CompoundName(t1)==A_ALG1)
			{
			Term sub, a1;
			sub=alg1_inv_alg(CompoundArg2(t1));
			if(sub==0)
			{
				ErrorInfo(375);
				puts("polynomial in denominator");
				longjmp(alg1_jmp_buf,1);
			}
			
			l1=ConsumeCompoundArg(m1,4);
			l1=CutFromList(l1,l);
			SetCompoundArg(m1,4,l1);
			
			a1=sub;	
			
			return SetIntAlgs(mk_let_d(m1,a1));
			
			}
		l=ListTail(l);
		}
		
	

	
	return AppendFirst(NewList(),m1);
	}

static Term alg1_df_df(Term m, int pos, List il)
{
	List nl,l1,l2;
	Term obj;
	
	nl=ConsumeCompoundArg(m,3);
	
	l1=ListNthList(nl,pos);
	
	obj=ListFirst(l1);
	ChangeList(l1,0);
	nl=CutFromList(nl,l1);
	
	l1=ConsumeCompoundArg(obj,1);
	FreeAtomic(obj);
	obj=l1;
	
	for(l1=obj,l2=il;l1;l1=ListTail(l1),l2=ListTail(l2))
		nl=AppendLast(nl,MakeCompound2(OPR_SPECIAL,
				MakeList2(ListFirst(l1),CopyTerm(ListFirst(l2))),A_DELTA));
	
	SetCompoundArg(m,3,nl);
	RemoveList(obj);
	
	return m;
}
	
			
static Term alg1_df(Term a1, Atom p)
{
	Term ind1, ind2,l1,l2,mo,mn=0;
	ind1=ConsumeCompoundArg(a1,2);
	ind2=CopyTerm(GetAtomProperty(p,PROP_INDEX));
	
	for(l1=ind2;l1;l1=ListTail(l1))
	{
		Label l;
		Term g, r1,r2;
		l=NewLabel();
		SetCompoundArg(ListFirst(l1),2,l);
		g=CompoundArg1(ListFirst(l1));
		r1=ConsumeCompoundArg(g,1);
		r2=ConsumeCompoundArg(g,2);
		SetCompoundArg(g,1,r2);
		SetCompoundArg(g,2,r1);
	}
	
	mo=ConsumeCompoundArg(a1,1);
	
	for(l1=mo;l1;l1=ListTail(l1))
	{
		Term mt;
		int i;
		mt=ListFirst(l1);
		for(l2=CompoundArgN(mt,4);l2;l2=ListTail(l2))
			if(CompoundArg2(ListFirst(l2))==p)
			{
				ErrorInfo(371);
				printf("df: object '%s' in denominator.\n",AtomValue(p));
				FreeAtomic(mo);
				return a1;
			}
		
		for(l2=CompoundArgN(mt,3),i=1;l2;l2=ListTail(l2),i++)
			if(CompoundArg2(ListFirst(l2))==p)
				mn=AppendLast(mn,alg1_df_df(CopyTerm(mt),i,ind2));
				
		
	}
	
	FreeAtomic(mo);
	
	SetCompoundArg(a1,1,mn);
	SetCompoundArg(a1,2,ConcatList(ind1,ind2));

	return a1;
}			

void renewlab(Term);

void alg1_kl_to_ia(List al)
{
	for(;al;al=ListTail(al))
	{
		List l1;
		for(l1=CompoundArgN(ListFirst(al),3);l1;l1=ListTail(l1))
		{
			Term prp;
			if(CompoundName(ListFirst(l1))==OPR_LET &&
				(prp=GetAtomProperty(CompoundArg2(ListFirst(l1)),A_KEEP_LETS))
				&& CompoundArg1(prp))
			{
				List il,kll;
				kll=CopyTerm(CompoundArg1(prp));
				renewlab(kll);
				il=ConsumeCompoundArg(ListFirst(l1),1);
				FreeAtomic(ListFirst(l1));
				ChangeList(l1,MakeCompound2(A_ALG1,il,kll));
			}
		}
	}
}


List SetIntAlgs(List l)
	{
	List l1,lr;
	
	
	lr=NewList();
	l1=l;
	while(!is_empty_list(l1))
		{
		lr=ConcatList(lr,s_l_1(ListFirst(l1)));
		l1=ListTail(l1);
		}
	RemoveList(l);	
	return lr;
	}
	
	

Term ProcDF(Term t, Term ind)
{
	
	Term t1, ti;
	int ano, i;
	
	if(!is_compound(t) || CompoundArity(t)<2 || !is_atom(CompoundArg1(t)) )
	{
		ErrorInfo(367);
		printf("wrong arguments in 'df' call.\n");
		return 0;
	}
	
	ano=CompoundArity(t)-1;
	
	t1=CompoundArg1(t);
	
	for(i=0;i<ano;i++)
	if(!is_atom(CompoundArgN(t,i+2)))
	{
		ErrorInfo(367);
		printf("wrong arguments in 'df' call.\n");
		return 0;
	}
	
	for(i=-1;i<ano;i++)
	if(GetAtomProperty(CompoundArgN(t,i+2),PROP_TYPE)==0)
	{
		ErrorInfo(368);
		printf("undefined object '%s'\n",AtomValue(CompoundArgN(t,i+2)));
		return 0;
	}
		
	
	ti=GetAtomProperty(t1,A_KEEP_LETS);
	if(ti==0 || CompoundArg1(ti)==0)
	{
		ErrorInfo(369);
		printf("'df': object '%s' must be declared in keep_lets statement.\n",
				AtomValue(t1));
		return 0;
	}
	
	
	ti=CopyTerm(ti);
	renewlab(ti);
	ti=CompoundArg1(ti);
	
	for(i=0;i<ano;i++)
		ti=alg1_df(ti,CompoundArgN(t,i+2));
	
	FreeAtomic(t);
	
	
	if(CompoundArg1(ti)==0)
	{
		FreeAtomic(ti);
		return NewInteger(0);
	}
	else
	{
		if(ind)
			ti=il_to_caret(ti,ind);
		
		return ti;
	}
	
}

static int has_ued=0, has_uedl=0;

static Term deriv5th=0;

Term ProcUED5th(Term t, Term ind)
{
	Term t1, ali;
	int d5th0=0;
	
	if(!is_compound(t) || CompoundArity(t)!=1)
	{
		ErrorInfo(0);
		printf("wrong syntax in ued_5th statement.\n");
		return 0;
	}
	
	t1=CommaToList(ConsumeCompoundArg(t,1));
	FreeAtomic(t);
	t=t1;
	
	for(t1=t;t1;t1=ListTail(t1))
	{
		Term v=ListFirst(t1);
		Term f, l;
		List il1=0,il2;
		if(!is_compound(v) || CompoundArity(v)!=2 || 
			!is_atom(CompoundArg1(v)))
		{
			ErrorInfo(0);
			printf("ued_5th: wrong expression ");WriteTerm(v);puts("");
			continue;
		}
		
		f=CompoundArg1(v);
		l=CompoundArg2(v);

		if(l==NewInteger(0))
		{
			SetAtomProperty(f,A_DERIV5,0);
			if(f==A_DERIV)
				{
				SetAtomProperty(A_MOMENT,A_DERIV5,0);
				deriv5th=0;
				d5th0=1;
				}
			continue;
		}
		
		has_ued++;
		if(f==A_DERIV5)
		{
			f=A_MOMENT;
			l=MakeCompound2(OPR_MLT,l,A_I);
			l=MakeCompound2(OPR_MLT,l,A_DERIV5);
		}
		
		if(f==A_DERIV)
		{
			deriv5th=CopyTerm(l);
		}
		
		
		if(f!=A_MOMENT && !is_particle(f,&il1) && !is_let(f,&il1))		
		{
			ErrorInfo(0);
			printf("ued_5th: is not a particle: ");WriteTerm(f);puts("");
			continue;
		}
		
		if(is_let(f,0) && f!=A_DERIV)
			has_uedl=1;
		
		if(il1)
			{
			il1=CopyTerm(il1);
			il1=CutFromList(il1,il1);
			}
		
		if(!is_atom(l) || !is_let(l,0))
		{
			Atom nl;
			char cbuf[32];
			sprintf(cbuf,"__5th__%05d",has_ued);
			nl=NewAtom(cbuf,0);
			l=MakeCompound1(OPR_LET,MakeCompound2(OPR_EQSIGN,nl,l));
			/*WriteTerm(l);puts("");*/
			CallFunction(l,0);
			l=nl;
			il2=GetAtomProperty(nl,PROP_INDEX);
			if(!EqualTerms(il1,il2))
			{
				ErrorInfo(0);
				WriteTerm(il1);puts("");
				WriteTerm(il2);puts("");
				printf("ued_5th: wrong indices for '%s'.\n",AtomValue(f));
				continue;	
			}
			il2=GetAtomProperty(nl,PROP_TYPE);
			if(il2==0 || CompoundName(il2)!=OPR_LET)
			{
				ErrorInfo(0);
				printf("ued_5th: rule for '%s' is not set.\n",AtomValue(f));
				continue;	
			}
		}
		SetAtomProperty(f,A_DERIV5,l);
		ali=MakeCompound1(OPR_ALIAS,MakeCompound2(OPR_EQSIGN,
			MakeCompound2(OPR_CARET,f,NewInteger(5)),l));
		CallFunction(ali,0);
		if(il1 && ListLength(il1)==1)
		{
		ali=MakeCompound1(OPR_ALIAS,MakeCompound2(OPR_EQSIGN,
			MakeCompound2(OPR_CARET,f,
					MakeCompound2(OPR_CARET,A_I,NewInteger(5))),
					MakeCompound2(OPR_CARET,l,A_I)));
		CallFunction(ali,0);
		
		}
		
	}

	
	if(GetAtomProperty(A_MOMENT,A_DERIV5)==0)
		{
		Term nl,l;
		char cbuf[32];
		if(deriv5th==0)
			{
			if(d5th0) return 0;
			WarningInfo(0);
			puts("ued_5th: the rule for 'deriv' must be set.");
			return 0;
			}
		
		l=MakeCompound2(OPR_MLT,A_I,deriv5th);
		sprintf(cbuf,"__5th__%05d",++has_ued);
		nl=NewAtom(cbuf,0);
		l=MakeCompound1(OPR_LET,MakeCompound2(OPR_EQSIGN,nl,l));
		CallFunction(l,0);
		SetAtomProperty(A_MOMENT,A_DERIV5,nl);
		}
	
	return 0;
	
}

void alg1_let5th(Term a1)
{
	List l1;
	
	if(a1==0 || has_uedl==0)
		return;
	
	List newt=0;
	
	for(l1=a1;l1;l1=ListTail(l1))
	{
		List pl=CompoundArgN(ListFirst(l1),3);
		List l2, ul=0;
		Label la;
		int i1,i2;
		for(l2=pl,i1=1;l2;l2=ListTail(l2),i1++)
		{
			List l3,l4,l5;
			if(CompoundName(ListFirst(l2))!=OPR_LET)
				continue;
			if(GetAtomProperty(CompoundArg2(ListFirst(l2)),A_DERIV5)==0)
				continue;
			if(is_particle(CompoundArg2(ListFirst(l2)),0))
				continue;
			if(GetAtomProperty(CompoundArg2(ListFirst(l2)),OPR_SCALAR)==NewInteger(0))
				continue;
			
			/*DumpList(pl);puts("-->");*/

			l3=CompoundArg1(ListFirst(l2));
			if(CompoundName(CompoundArg1(ListFirst(l3)))!=A_LORENTZ ||
				   CompoundArg1(CompoundArg1(ListFirst(l3)))!=NewInteger(2))
				{
				ErrorInfo(0);printf("ued_5th: '%s' is not vector field.\n",
						AtomValue(CompoundArg2(ListFirst(l2))));
				l3=0;
				}
			if(l3)
			{
				la=CompoundArg2(ListFirst(l3));
				List mlc;
				if(ListMember(ul,la))
					continue;
				for(l4=pl,i2=1;l4;l4=ListTail(l4),i2++)
				{
					if(l4==l2)
						continue;
					for(l5=CompoundArg1(ListFirst(l4));l5;l5=ListTail(l5))
						if(CompoundArg2(ListFirst(l5))==la)
							break;
					if(l5)
						break;
				}
				if(l4==0)
					{ErrorInfo(0);puts("Internal error c5thgen.");
					continue;}
				ul=AppendLast(ul,la);
				
				if(CompoundName(ListFirst(l4))==OPR_SPECIAL &&
					(CompoundArg2(ListFirst(l4))!=A_GAMMA &&
					 CompoundArg2(ListFirst(l4))!=A_MOMENT))
					{
					WarningInfo(0);
					WriteTerm(CompoundArg2(ListFirst(l4)));
					puts(": 5th component is not defined.");
					continue;
					}
				if(CompoundName(ListFirst(l2))==OPR_FIELD &&
					GetAtomProperty(CompoundArg2(ListFirst(l2)),A_DERIV5)==0)
					continue;
				if(CompoundName(ListFirst(l4))==OPR_FIELD &&
					GetAtomProperty(CompoundArg2(ListFirst(l4)),A_DERIV5)==0)
					continue;
				mlc=CopyTerm(ListFirst(l1));

				SetCompoundArg(mlc,1,NewInteger(
						-IntegerValue(CompoundArg1(mlc))));

				l3=ListNthList(CompoundArgN(mlc,3),i1);
				l4=ConsumeCompoundArg(ListFirst(l3),1);
				for(l5=l4;l5;l5=ListTail(l5))
					if(CompoundArg2(ListFirst(l5))==la)
						{
						l4=CutFromList(l4,l5);
						break;
						}
				SetCompoundArg(ListFirst(l3),1,l4);
				
				{
					Term prp;
					prp=GetAtomProperty(CompoundArg2(ListFirst(l3)),A_DERIV5);
					if(prp==0)
					{
						ErrorInfo(0);
						printf("ued_5th for '");
						WriteTerm(CompoundArg2(ListFirst(l3)));
						puts("' is not defined.");
						continue;
					}
					SetCompoundName(ListFirst(l3),OPR_LET);
					SetCompoundArg(ListFirst(l3),2,prp);
				}
				l3=ListNthList(CompoundArgN(mlc,3),i2);
				l4=ConsumeCompoundArg(ListFirst(l3),1);
				for(l5=l4;l5;l5=ListTail(l5))
					if(CompoundArg2(ListFirst(l5))==la)
						{
						l4=CutFromList(l4,l5);
						break;
						}
				SetCompoundArg(ListFirst(l3),1,l4);
				if(CompoundArg2(ListFirst(l3))==A_GAMMA)
				{
					SetCompoundArg(ListFirst(l3),2,A_GAMMA5);
					SetCompoundArg(mlc,2,NewInteger(
						-IntegerValue(CompoundArg2(mlc))));
					SetCompoundArg(mlc,1,NewInteger(
						-IntegerValue(CompoundArg1(mlc))));
				}
				else
				{
					Term prp;
					prp=GetAtomProperty(CompoundArg2(ListFirst(l3)),A_DERIV5);
					if(prp==0)
					{
						ErrorInfo(0);
						printf("ued_5th for '");
						WriteTerm(CompoundArg2(ListFirst(l3)));
						puts("' is not defined.");
						continue;
					}
					SetCompoundName(ListFirst(l3),OPR_LET);
					SetCompoundArg(ListFirst(l3),2,prp);
				}
				/*newt=AppendLast(newt,mlc);*/
				/*puts("new is generated");*/
				/*printf("--> ");DumpList(CompoundArgN(mlc,3));*/
				/*mlc=MakeList1(mlc);
				mlc=SetLets(mlc);*/
				/*DumpList(mlc);	*/			
				if(ListLength(ul)==1)
					AppendLast(a1,mlc);
				else
					newt=AppendLast(newt,mlc);
			}
		}
		if(ul) FreeAtomic(ul);

	}

	ConcatList(a1,newt);

}
	
void alg1_add5th(Term a1)
{
	List l1,ml;
		List newt=0;
		for(l1=CompoundArg1(a1);l1;l1=ListTail(l1))
		{
			List pl=CompoundArgN(ListFirst(l1),3);
			List l2, ul=0;
			Label la;
			int i1,i2;
			/*DumpList(pl);*/
			for(l2=pl,i1=1;l2;l2=ListTail(l2),i1++)
			{
				List l3,l4,l5;
				for(l3=CompoundArg1(ListFirst(l2));l3;l3=ListTail(l3))
					if(CompoundName(CompoundArg1(ListFirst(l3)))==A_LORENTZ &&
					   CompoundArg1(CompoundArg1(ListFirst(l3)))==NewInteger(2))
					   break;
				if(l3)
				{
					la=CompoundArg2(ListFirst(l3));
					List mlc;
					if(ListMember(ul,la))
						continue;
					for(l4=ListTail(l2),i2=i1+1;l4;l4=ListTail(l4),i2++)
					{
						for(l5=CompoundArg1(ListFirst(l4));l5;l5=ListTail(l5))
							if(CompoundArg2(ListFirst(l5))==la)
								break;
						if(l5)
							break;
					}
					if(l4==0)
						{ErrorInfo(0);puts("Internal error c5thgen.");
						continue;}
					ul=AppendLast(ul,la);
					if(CompoundName(ListFirst(l2))==OPR_SPECIAL &&
						(CompoundArg2(ListFirst(l2))!=A_GAMMA &&
						 CompoundArg2(ListFirst(l2))!=A_MOMENT))
						{
						WarningInfo(0);
						WriteTerm(CompoundArg2(ListFirst(l2)));
						puts(": 5th component is not defined.");
						continue;
						}
					if(CompoundName(ListFirst(l4))==OPR_SPECIAL &&
						(CompoundArg2(ListFirst(l4))!=A_GAMMA &&
						 CompoundArg2(ListFirst(l4))!=A_MOMENT))
						{
						WarningInfo(0);
						WriteTerm(CompoundArg2(ListFirst(l4)));
						puts(": 5th component is not defined.");
						continue;
						}
					if(CompoundName(ListFirst(l2))==OPR_FIELD &&
						GetAtomProperty(CompoundArg2(ListFirst(l2)),A_DERIV5)==0)
						continue;
					if(CompoundName(ListFirst(l4))==OPR_FIELD &&
						GetAtomProperty(CompoundArg2(ListFirst(l4)),A_DERIV5)==0)
						continue;
					mlc=CopyTerm(ListFirst(l1));
					
					SetCompoundArg(mlc,1,NewInteger(
							-IntegerValue(CompoundArg1(mlc))));
							
					l3=ListNthList(CompoundArgN(mlc,3),i1);
					l4=ConsumeCompoundArg(ListFirst(l3),1);
					for(l5=l4;l5;l5=ListTail(l5))
						if(CompoundArg2(ListFirst(l5))==la)
							{
							l4=CutFromList(l4,l5);
							break;
							}
					SetCompoundArg(ListFirst(l3),1,l4);
					if(CompoundArg2(ListFirst(l3))==A_GAMMA)
					{
						SetCompoundArg(ListFirst(l3),2,A_GAMMA5);
						SetCompoundArg(mlc,2,NewInteger(
							-IntegerValue(CompoundArg2(mlc))));
					SetCompoundArg(mlc,1,NewInteger(
							-IntegerValue(CompoundArg1(mlc))));
					
							
					}
					else
					{
						Term prp;
						prp=GetAtomProperty(CompoundArg2(ListFirst(l3)),A_DERIV5);
						if(prp==0)
						{
							ErrorInfo(0);
							printf("ued_5th for '");
							WriteTerm(CompoundArg2(ListFirst(l3)));
							puts("' is not defined.");
							continue;
						}
						SetCompoundName(ListFirst(l3),OPR_LET);
						SetCompoundArg(ListFirst(l3),2,prp);
					}
					l3=ListNthList(CompoundArgN(mlc,3),i2);
					l4=ConsumeCompoundArg(ListFirst(l3),1);
					for(l5=l4;l5;l5=ListTail(l5))
						if(CompoundArg2(ListFirst(l5))==la)
							{
							l4=CutFromList(l4,l5);
							break;
							}
					SetCompoundArg(ListFirst(l3),1,l4);
					if(CompoundArg2(ListFirst(l3))==A_GAMMA)
					{
						SetCompoundArg(ListFirst(l3),2,A_GAMMA5);
						SetCompoundArg(mlc,2,NewInteger(
							-IntegerValue(CompoundArg2(mlc))));
					SetCompoundArg(mlc,1,NewInteger(
							-IntegerValue(CompoundArg1(mlc))));
					puts("2");
					}
					else
					{
						Term prp;
						prp=GetAtomProperty(CompoundArg2(ListFirst(l3)),A_DERIV5);
						if(prp==0)
						{
							ErrorInfo(0);
							printf("ued_5th for '");
							WriteTerm(CompoundArg2(ListFirst(l3)));
							puts("' is not defined.");
							continue;
						}
						SetCompoundName(ListFirst(l3),OPR_LET);
						SetCompoundArg(ListFirst(l3),2,prp);
					}
					/*newt=AppendLast(newt,mlc);*/
					/*puts("new is generated");*/
					/*DumpList(CompoundArgN(mlc,3));*/
					mlc=MakeList1(mlc);
					mlc=SetLets(mlc);
					/*DumpList(mlc);	*/			
					if(ListLength(ul)==1)
						ConcatList(CompoundArg1(a1),mlc);
					else
						newt=ConcatList(newt,mlc);
				}
			}
			if(ul) FreeAtomic(ul);
		
		}
	
		ConcatList(CompoundArg1(a1),newt);
	
}

void alg1_rem_sincos(Term a1)
{
	List l1,ml;
	
	/*alg1_dump(a1);*/
	
	if(has_ued) /* generate 5th components */
	{
	alg1_add5th(a1);
	}
	
	
	for(l1=CompoundArg1(a1);l1;l1=ListTail(l1))
	{
		List scl=0;
		List pl=ConsumeCompoundArg(ListFirst(l1),3);
		List l2;
		int df=1;
		int nf=0;
		
/*		DumpList(pl);*/
		
		rpt0:
		
		for(l2=pl;l2;l2=ListTail(l2))
		{
			if(CompoundArg2(ListFirst(l2))==A_DERIV5)
			{
				List l3;
				for(l3=ListTail(l2);l3;l3=ListTail(l3))
					if(CompoundArg2(ListFirst(l3))==A_COS ||
						CompoundArg2(ListFirst(l3))==A_SIN)
						break;
				if(l3==0)
				{
				SetCompoundArg(ListFirst(l1),1,NewInteger(0));
				WarningInfo(0);puts("deriv5 does not followed by field.");
				SetCompoundArg(ListFirst(l1),3,pl);
				}
				else
				{
					Atom prm=0;
					if(CompoundArg2(ListFirst(l3))==A_COS)
					{
					SetCompoundArg(ListFirst(l3),2,A_SIN);
					SetCompoundArg(ListFirst(l1),1,
						NewInteger(-IntegerValue(CompoundArg1(ListFirst(l1)))*
						IntegerValue(CompoundArgN(ListFirst(l3),3))));
					}
					else
					{
					SetCompoundArg(ListFirst(l3),2,A_COS);
					SetCompoundArg(ListFirst(l1),1,
						NewInteger(IntegerValue(CompoundArg1(ListFirst(l1)))*
						IntegerValue(CompoundArgN(ListFirst(l3),3))));
					}
					if(CompoundArity(ListFirst(l3))==4)
						prm=CompoundArgN(ListFirst(l3),4);
					if(prm)
						{
						/*WriteTerm(prm);puts("");*/
						if(CompoundArg1(prm)!=NewInteger(1))
							SetCompoundArg(ListFirst(l1),1,
							NewInteger(IntegerValue(CompoundArg1(ListFirst(l1)))
							*IntegerValue(CompoundArg1(prm))));
						if(CompoundArg2(prm)!=NewInteger(1))
							SetCompoundArg(ListFirst(l1),2,
							NewInteger(IntegerValue(CompoundArg2(ListFirst(l1)))
							*IntegerValue(CompoundArg2(prm))));
						if(CompoundArgN(prm,3))
							ConcatList(pl,CopyTerm(CompoundArgN(prm,3)));
						if(CompoundArgN(prm,4))
							{
							List dl=ConsumeCompoundArg(ListFirst(l1),4);
							dl=ConcatList(dl,CopyTerm(CompoundArgN(prm,4)));
							SetCompoundArg(ListFirst(l1),4,dl);
							}
						
						}
					
				}
				pl=CutFromList(pl,l2);
				goto rpt0;
			}
		}
		
		
		/*DumpList(pl);*/
		
		rpt:
		
		for(l2=pl;l2;l2=ListTail(l2))
		{
			Term t=ListFirst(l2);
			if(CompoundName(t)==OPR_PARAMETER && 
				(CompoundArg2(t)==A_COS || CompoundArg2(t)==A_SIN))
				{
					Term ii=CompoundArgN(t,3);
					if(!is_integer(ii))
					{
						puts("internal error (remsincos)");
						break;
					}
					scl=AppendLast(scl,MakeCompound1(CompoundArg2(t),ii));
					pl=CutFromList(pl,l2);
					goto rpt;
				}
		}
		
		if(scl==0 || CompoundArg1(ListFirst(l1))==NewInteger(0))
		{
			SetCompoundArg(ListFirst(l1),3,pl);
			continue;
		}
		/*WriteTerm(scl);puts("");*/
		
		scl=AppendFirst(0,MakeCompound1(OPR_PLUS,scl));
		while(ListLength(CompoundArg1(ListFirst(scl)))>1)
		{
			List nl=0;
			df*=2;
			/*DumpList(scl);*/
			for(l2=scl;l2;l2=ListTail(l2))
			{
				Term t1,t2;
				List r;
				Atom plus=CompoundName(ListFirst(l2));
				Atom minus=(plus==OPR_PLUS?OPR_MINUS:OPR_PLUS);
				r=ConsumeCompoundArg(ListFirst(l2),1);
				t1=ListFirst(r);
				ChangeList(r,0);r=CutFromList(r,r);
				t2=ListFirst(r);
				ChangeList(r,0);r=CutFromList(r,r);
				if(CompoundName(t1)==A_COS && CompoundName(t2)==A_COS)
				{
					nl=AppendLast(nl,MakeCompound1(plus,AppendFirst(
						CopyTerm(r),MakeCompound1(A_COS,NewInteger(
						IntegerValue(CompoundArg1(t1))-
						IntegerValue(CompoundArg1(t2)))))));
					nl=AppendLast(nl,MakeCompound1(plus,AppendFirst(
						          r,MakeCompound1(A_COS,NewInteger(
						IntegerValue(CompoundArg1(t1))+
						IntegerValue(CompoundArg1(t2)))))));
				}
				if(CompoundName(t1)==A_SIN && CompoundName(t2)==A_SIN)
				{
					nl=AppendLast(nl,MakeCompound1(plus,AppendFirst(
						CopyTerm(r),MakeCompound1(A_COS,NewInteger(
						IntegerValue(CompoundArg1(t1))-
						IntegerValue(CompoundArg1(t2)))))));
					nl=AppendLast(nl,MakeCompound1(minus,AppendFirst(
						          r,MakeCompound1(A_COS,NewInteger(
						IntegerValue(CompoundArg1(t1))+
						IntegerValue(CompoundArg1(t2)))))));
				}
				if(CompoundName(t1)==A_SIN && CompoundName(t2)==A_COS)
				{
					nl=AppendLast(nl,MakeCompound1(plus,AppendFirst(
						CopyTerm(r),MakeCompound1(A_SIN,NewInteger(
						IntegerValue(CompoundArg1(t1))-
						IntegerValue(CompoundArg1(t2)))))));
					nl=AppendLast(nl,MakeCompound1(plus,AppendFirst(
						          r,MakeCompound1(A_SIN,NewInteger(
						IntegerValue(CompoundArg1(t1))+
						IntegerValue(CompoundArg1(t2)))))));
				}
				if(CompoundName(t1)==A_COS && CompoundName(t2)==A_SIN)
				{
					nl=AppendLast(nl,MakeCompound1(minus,AppendFirst(
						CopyTerm(r),MakeCompound1(A_SIN,NewInteger(
						IntegerValue(CompoundArg1(t1))-
						IntegerValue(CompoundArg1(t2)))))));
					nl=AppendLast(nl,MakeCompound1(plus,AppendFirst(
						          r,MakeCompound1(A_SIN,NewInteger(
						IntegerValue(CompoundArg1(t1))+
						IntegerValue(CompoundArg1(t2)))))));
				}
				
				FreeAtomic(t1);
				FreeAtomic(t2);
			}
			
			FreeAtomic(scl);
			
			scl=nl;
		}
			
		
		SetCompoundArg(ListFirst(l1),3,pl);
		
		for(l2=scl;l2;l2=ListTail(l2))
		{
			Term t1=ListFirst(l2);
			Term t2=ListFirst(CompoundArg1(t1));
			if(CompoundName(t1)==OPR_PLUS && CompoundName(t2)==A_COS 
					&& CompoundArg1(t2)==NewInteger(0))
				nf++;
			if(CompoundName(t1)==OPR_MINUS && CompoundName(t2)==A_COS 
					&& CompoundArg1(t2)==NewInteger(0))
				nf--;
		}
		
		/*printf("%d/%d\n",nf,df);*/
		
		if(nf==0)
			SetCompoundArg(ListFirst(l1),1,0);
		else
		{
			long int n=IntegerValue(CompoundArg1(ListFirst(l1))),
				d=IntegerValue(CompoundArg2(ListFirst(l1))),
				gc;
			n*=nf;
			d*=df;
			gc=gcf(n,d);
			n/=gc;
			d/=gc;
			SetCompoundArg(ListFirst(l1),1,NewInteger(n));
			SetCompoundArg(ListFirst(l1),2,NewInteger(d));
		}
			
	/*WriteTerm(ListFirst(l1));puts("");*/
	
	}
	
	ml=ConsumeCompoundArg(a1,1);
/*	
rpt2:
	for(l1=ml;l1;l1=ListTail(l1))
		if(CompoundArg1(ListFirst(l1))==0)
			{
			ml=CutFromList(ml,l1);
			goto rpt2;
			}
*/
	
	{
	List ml2=0, mle=0;
	for(l1=ml;l1;l1=ListTail(l1))
		if(CompoundArg1(ListFirst(l1))==0)
			FreeAtomic(ListFirst(l1));
		else if(ml2==0)
			{
			ml2=AppendLast(ml2,ListFirst(l1)); mle=ml2;
			}
		else
			{
			AppendLast(mle,ListFirst(l1));
			mle=ListTail(mle);
			}
	RemoveList(ml);
	ml=ml2;
	}
			
	

	
	SetCompoundArg(a1,1,ml);

}


