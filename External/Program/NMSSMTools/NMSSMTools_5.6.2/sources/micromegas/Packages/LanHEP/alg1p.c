#include <string.h>
#include "lanhep.h"

static int prtcmp(Term p1, Term p2)
	{
	Term pp1,pp2;
	pp1=CompoundArg1(p1);
	pp2=CompoundArg1(p2);
	if(is_atom(pp1) && is_atom(pp2))
		return strcmp(AtomValue(pp1),AtomValue(pp2));
	if(is_integer(pp1) && is_integer(pp2))
		return (int)(IntegerValue(pp1)-IntegerValue(pp2));
	if(is_integer(pp1))
		return -1;
	return 1;
	}

		
static List mn_mlt(List l, Atom par, int pw)
	{
	List l1;
	
	for(l1=ListTail(l);l1;l1=ListTail(l1))
		{
		Term pp;
		pp=ListFirst(l1);
		if(CompoundArg1(pp)==par)
			{
			int n;
			n=pw+(int)IntegerValue(CompoundArg2(pp));
			if(n==0)
				return CutFromList(l,l1);
			if(par==A_SQRT2)
			{
				long int n,d,c;
				n=IntegerValue(CompoundArg1(ListFirst(l)));
				d=IntegerValue(CompoundArg2(ListFirst(l)));
				if(pw>0)
					n*=2;
				else
					d*=2;
				c=gcf(n,d);
				n/=c;
				d/=c;
				SetCompoundArg(ListFirst(l),1,NewInteger(n));
				SetCompoundArg(ListFirst(l),2,NewInteger(d));
				return CutFromList(l,l1);
			}
			SetCompoundArg(pp,2,NewInteger(n));
			return l;
			}
		}
	l=AppendLast(l,MakeCompound2(OPR_POW,par,NewInteger(pw)));
	l=SortedList(l,prtcmp);
	return l;
	} 		


		
static int is_comp(Term a1)
{
	List l;
	for(l=CompoundArg1(a1);l;l=ListTail(l))
	{
		List l2;
		int fno=0;
		if(CompoundArgN(ListFirst(l),4))
			return 0;
		for(l2=CompoundArgN(ListFirst(l),3);l2;l2=ListTail(l2))
		{
			if(CompoundName(ListFirst(l2))==OPR_PARAMETER)
				continue;
			if(CompoundName(ListFirst(l2))==OPR_FIELD)
			{
				fno++;
				continue;
			}
			return 0;
		}
		if(fno!=1)
			return 0;
	}
	return 1;
}


static Term mk_mpl_el(Term a1)
{
	List ret, l;
	ret=NewList();
	for(l=CompoundArg1(a1);l;l=ListTail(l))
	{
		Atom prt=0;
		List prms=0;
		List l1;
		long int num, den;
		num=IntegerValue(CompoundArg1(ListFirst(l)));
		den=IntegerValue(CompoundArg2(ListFirst(l)));
		
		prms=AppendFirst(prms,MakeCompound2(OPR_DIV,
				NewInteger(num),NewInteger(den)));
		
		for(l1=CompoundArgN(ListFirst(l),3);l1;l1=ListTail(l1))
		{
			if(CompoundName(ListFirst(l1))==OPR_FIELD)
				prt=CompoundArg2(ListFirst(l1));
			else
				prms=mn_mlt(prms,CompoundArg2(ListFirst(l1)),1);
		}
		if(prt==0)
		{
			puts("Internal error (a1p1)");
			exit(0);
		}
		ret=AppendLast(ret,MakeCompound1(prt,prms));
	}
	
	return ret;
}


Term alg1_guess_mpl(Term a1)
{
	List ml, l;
	Term ret;
	int i,ii;
	ml=CompoundArg1(a1);
		
	if(ListLength(ml)!=1)
		return 0;
	ml=ListFirst(ml);
	
	if(CompoundArg1(ml)!=NewInteger(1) || CompoundArg2(ml)!=NewInteger(1))
		return 0;
	if(CompoundArgN(ml,4))
		return 0;
	ml=CompoundArgN(ml,3);
	if(ListLength(ml)!=1)
		return 0;
	ml=ListFirst(ml);
	if(CompoundName(ml)!=OPR_WILD)
		return 0;


	for(l=CompoundArg2(ml);l;l=ListTail(l))
		if(!is_comp(ListFirst(l)))
			return 0;
	
	ii=ListLength(CompoundArg2(ml));
	ret=MakeCompound(A_I,ii);
	
	for(i=1,l=CompoundArg2(ml);i<=ii;i++,l=ListTail(l))
		SetCompoundArg(ret,i,mk_mpl_el(ListFirst(l)));
	
	
	return ret;

}


static void remake_2(Term a1, List sub, Label ind1, Label ind2, int d1, int d2)
{
	List tind1, tind2;
	List l1,l2;
	List r1,r2;
	Term ww,mm;

	
	tind1=MakeCompound(OPR_WILD,3);
	SetCompoundArg(tind1,1,
			MakeCompound2(OPR_WILD,NewInteger(d1),NewInteger(d1)));
	SetCompoundArg(tind1,2,ind1);
	SetCompoundArg(tind1,3,NewInteger(d1));
	
	tind2=MakeCompound(OPR_WILD,3);
	SetCompoundArg(tind2,1,
			MakeCompound2(OPR_WILD,NewInteger(d2),NewInteger(d2)));
	SetCompoundArg(tind2,2,ind2);
	SetCompoundArg(tind2,3,NewInteger(d2));
	
	r1=0;
	
	for(l1=sub;l1;l1=ListTail(l1))
	{
		Term mm, ww, aa;
		r2=0;
		
		for(l2=ListFirst(l1);l2;l2=ListTail(l2))
		{
			List l3;
			Term aa1;
			Term aa1ml=0;
			
			for(l3=ListFirst(l2);l3;l3=ListTail(l3))
			{
				List ln=0,ld=0,l4;	
				Term m1;
				m1=MakeCompound(A_MTERM,4);
				SetCompoundArg(m1,1,CompoundArg1(ListFirst(ListFirst(l3))));
				SetCompoundArg(m1,2,CompoundArg2(ListFirst(ListFirst(l3))));
				for(l4=ListTail(ListFirst(l3));l4;l4=ListTail(l4))
				{
					Atom prm;
					int pw,i;
					Term pt;
					prm=CompoundArg1(ListFirst(l4));
					pw=(int)IntegerValue(CompoundArg2(ListFirst(l4)));

					if(pw>0)
					{
						for(i=0;i<pw;i++)
						{
							pt=MakeCompound2(OPR_PARAMETER,0,prm);
							ln=AppendLast(ln,pt);
						}
					}
					else
					{
						for(i=0;i<-pw;i++)
						{
							pt=MakeCompound2(OPR_PARAMETER,0,prm);
							ld=AppendLast(ld,pt);
						}
					}
				}
				SetCompoundArg(m1,3,ln);
				SetCompoundArg(m1,4,ld);
				aa1ml=AppendLast(aa1ml,m1);
			}
			aa1=MakeCompound2(A_ALG1,aa1ml,0);
			r2=AppendLast(r2,aa1);
		}
		
		ww=MakeCompound2(OPR_WILD,MakeList1(CopyTerm(tind1)),r2);
		mm=MakeCompound(A_MTERM,4);
		SetCompoundArg(mm,1,NewInteger(1));
		SetCompoundArg(mm,2,NewInteger(1));
		SetCompoundArg(mm,3,MakeList1(ww));
		aa=MakeCompound2(A_ALG1,MakeList1(mm),MakeList1(CopyTerm(tind1)));
		r1=AppendLast(r1,aa);
		
	}
	
	ww=MakeCompound2(OPR_WILD,MakeList2(CopyTerm(tind1),CopyTerm(tind2)),r1);
	mm=MakeCompound(A_MTERM,4);
	SetCompoundArg(mm,1,NewInteger(1));
	SetCompoundArg(mm,2,NewInteger(1));
	SetCompoundArg(mm,3,MakeList1(ww));
/*	
	printf("remake2: %d\n",EqualTerms(mm,ListFirst(CompoundArg1(a1))));
*/	
	FreeAtomic(ConsumeCompoundArg(a1,1));
	SetCompoundArg(a1,1,MakeList1(mm));
		
}

static void remake_3(Term a1, List sub, Label ind1, Label ind2, Label ind3,
		int d1, int d2, int d3)
{
	List tind1, tind2, tind3;
	List l1,l2,l22;
	List r1,r2,r3;
	Term ww,mm;

	tind1=MakeCompound(OPR_WILD,3);
	SetCompoundArg(tind1,1,
			MakeCompound2(OPR_WILD,NewInteger(d1),NewInteger(d1)));
	SetCompoundArg(tind1,2,ind1);
	SetCompoundArg(tind1,3,NewInteger(d1));
	
	tind2=MakeCompound(OPR_WILD,3);
	SetCompoundArg(tind2,1,
			MakeCompound2(OPR_WILD,NewInteger(d2),NewInteger(d2)));
	SetCompoundArg(tind2,2,ind2);
	SetCompoundArg(tind2,3,NewInteger(d2));
	
	tind3=MakeCompound(OPR_WILD,3);
	SetCompoundArg(tind3,1,
			MakeCompound2(OPR_WILD,NewInteger(d3),NewInteger(d3)));
	SetCompoundArg(tind3,2,ind3);
	SetCompoundArg(tind3,3,NewInteger(d3));
	
	r1=0;
	
	for(l1=sub;l1;l1=ListTail(l1))
	{
		Term mm, ww, aa;
		r2=0;
		
		for(l2=ListFirst(l1);l2;l2=ListTail(l2))
		{
			Term mm, ww, aa;
			r3=0;
			
			for(l22=ListFirst(l2);l22;l22=ListTail(l22))
			{
				List l3;
				Term aa1;
				Term aa1ml=0;

				for(l3=ListFirst(l22);l3;l3=ListTail(l3))
				{
					List ln=0,ld=0,l4;	
					Term m1;
					m1=MakeCompound(A_MTERM,4);
					SetCompoundArg(m1,1,CompoundArg1(ListFirst(ListFirst(l3))));
					SetCompoundArg(m1,2,CompoundArg2(ListFirst(ListFirst(l3))));
					for(l4=ListTail(ListFirst(l3));l4;l4=ListTail(l4))
					{
						Atom prm;
						int pw,i;
						Term pt;
						prm=CompoundArg1(ListFirst(l4));
						pw=(int)IntegerValue(CompoundArg2(ListFirst(l4)));

						if(pw>0)
						{
							for(i=0;i<pw;i++)
							{
								pt=MakeCompound2(OPR_PARAMETER,0,prm);
								ln=AppendLast(ln,pt);
							}
						}
						else
						{
							for(i=0;i<-pw;i++)
							{
								pt=MakeCompound2(OPR_PARAMETER,0,prm);
								ld=AppendLast(ld,pt);
							}
						}
					}
					SetCompoundArg(m1,3,ln);
					SetCompoundArg(m1,4,ld);
					aa1ml=AppendLast(aa1ml,m1);
				}
				
				aa1=MakeCompound2(A_ALG1,aa1ml,0);
				r3=AppendLast(r3,aa1);

			}
			ww=MakeCompound2(OPR_WILD,MakeList1(CopyTerm(tind1)),r3);
			mm=MakeCompound(A_MTERM,4);
			SetCompoundArg(mm,1,NewInteger(1));
			SetCompoundArg(mm,2,NewInteger(1));
			SetCompoundArg(mm,3,MakeList1(ww));
			aa=MakeCompound2(A_ALG1,MakeList1(mm),MakeList1(CopyTerm(tind1)));
			r2=AppendLast(r2,aa);
			
		}
		
		ww=MakeCompound2(OPR_WILD,
				MakeList2(CopyTerm(tind1),CopyTerm(tind2)),r2);
		mm=MakeCompound(A_MTERM,4);
		SetCompoundArg(mm,1,NewInteger(1));
		SetCompoundArg(mm,2,NewInteger(1));
		SetCompoundArg(mm,3,MakeList1(ww));
		aa=MakeCompound2(A_ALG1,MakeList1(mm),
				MakeList2(CopyTerm(tind1),CopyTerm(tind2)));
		r1=AppendLast(r1,aa);
		
	}
	
	ww=MakeCompound2(OPR_WILD,
			MakeList3(CopyTerm(tind1),CopyTerm(tind2),CopyTerm(tind3)),r1);
	mm=MakeCompound(A_MTERM,4);
	SetCompoundArg(mm,1,NewInteger(1));
	SetCompoundArg(mm,2,NewInteger(1));
	SetCompoundArg(mm,3,MakeList1(ww));
/*	
	printf("remake3: %d\n",EqualTerms(mm,ListFirst(CompoundArg1(a1))));
*/
	FreeAtomic(ConsumeCompoundArg(a1,1));
	SetCompoundArg(a1,1,MakeList1(mm));
		
}

static int mtr_has_prm, mtr_has_sum;

static List a12pl_add(List sum, List m)
{
	List l;
	for(l=sum;l;l=ListTail(l))
	{
		Term l1;
		long int n1,n2,d1,d2,n,d,c;
		
		l1=ListFirst(l);
		if(!EqualTerms(ListTail(l1),ListTail(m)))
			continue;
		if(IntegerValue(CompoundArg2(ListFirst(l1))) *
				IntegerValue(CompoundArg2(ListFirst(m)))<0)
			continue;
		
		n1=IntegerValue(CompoundArg1(ListFirst(l1)));
		d1=IntegerValue(CompoundArg2(ListFirst(l1)));
		n2=IntegerValue(CompoundArg1(ListFirst(m)));
		d2=IntegerValue(CompoundArg2(ListFirst(m)));
		
		n=n1*d2+n2*d1;
		if(n==0)
		{
			FreeAtomic(m);
			return CutFromList(sum,l);
		}
		
		d=d1*d2;
		
		if(d1<0)
			d=-d,n=-n;
		
		c=gcf(n,d);
		n/=c;
		d/=c;
		SetCompoundArg(ListFirst(l1),1,NewInteger(n));
		SetCompoundArg(ListFirst(l1),2,NewInteger(d));
		FreeAtomic(m);
		return sum;
	}
		
	return AppendLast(sum,m);
}


static List a12pl_mk1(Term m1)
{
	List ret=0,l;
	
	ret=AppendFirst(ret,
			MakeCompound2(OPR_DIV,CompoundArg1(m1),CompoundArg2(m1)));
	
	for(l=CompoundArgN(m1,3);l;l=ListTail(l))
		ret=mn_mlt(ret,CompoundArg2(ListFirst(l)),1);
	for(l=CompoundArgN(m1,4);l;l=ListTail(l))
		ret=mn_mlt(ret,CompoundArg2(ListFirst(l)),-1);
	
	return ret;
}


static List a1_2_pl(Term a1)
{
	List l1;
	List ret=0;
	
	for(l1=CompoundArg1(a1);l1;l1=ListTail(l1))
		ret=a12pl_add(ret,a12pl_mk1(ListFirst(l1)));
	
	return ret;
}


static int is_matr(Term t)
{
	if(CompoundName(t)==A_ALG1)
	{
		List l;
		for(l=CompoundArg1(t);l;l=ListTail(l))
		{
			List l1;
			for(l1=CompoundArgN(ListFirst(l),4);l1;l1=ListTail(l1))
				if(CompoundName(ListFirst(l1))!=OPR_PARAMETER)
					return 0;
			
			for(l1=CompoundArgN(ListFirst(l),3);l1;l1=ListTail(l1))
			{
				if(CompoundName(ListFirst(l1))==OPR_PARAMETER)
					continue;
				if(CompoundName(ListFirst(l1))==OPR_WILD &&
						is_matr(ListFirst(l1)))
					continue;
				return 0;
			}
		}
		return 1;
	}
	
	if(CompoundName(t)==OPR_WILD)
	{
		List l;
		for(l=CompoundArg2(t);l;l=ListTail(l))
			if(!is_matr(ListFirst(l)))
				return 0;
		return 1;
		
	}
	

	return 0;
}

extern void sub_ind_alg(Term a1, Atomic ind, int i);
		
Term alg1_guess_mtr(Term a1, int *tp)
{
	Term a11;
	long int dd,d1,d2,d3=0;
	Atomic ind1, ind2, ind3=0;
	int i,j,k;
	
	dd=ListLength(CompoundArg2(a1));
	if(dd!=2 && dd!=3)
		return 0;
	
	if(!is_matr(a1))
		return 0;
	
	
	d1=IntegerValue(CompoundArgN(ListFirst(CompoundArg2(a1)),3));
	d2=IntegerValue(CompoundArgN(ListFirst(ListTail(CompoundArg2(a1))),3));
	if(dd==3)
		d3=IntegerValue(CompoundArgN(
				ListFirst(ListTail(ListTail(CompoundArg2(a1)))),3));
	
	ind1=CompoundArg2(ListFirst(CompoundArg2(a1)));
	ind2=CompoundArg2(ListFirst(ListTail(CompoundArg2(a1))));
	if(dd==3)
		ind3=CompoundArg2(ListFirst(ListTail(ListTail(CompoundArg2(a1)))));
	
	mtr_has_prm=0;
	mtr_has_sum=0;
	
	if(dd==2)
	{
	List l1=0,l2;
	

		for(j=1;j<=d2;j++)
		{
			l2=0;
			for(i=1;i<=d1;i++)
			{
				Term res;

				a11=CopyTerm(a1);

				sub_ind_alg(a11,ind1,i);
				sub_ind_alg(a11,ind2,j);

				res=a1_2_pl(a11);

				if(ListLength(res)>1)
					mtr_has_sum=1;
				if(res && ListLength(ListFirst(res))>1)
					mtr_has_prm=1;

				FreeAtomic(a11);
				
				l2=AppendLast(l2,res);
			}
			l1=AppendLast(l1,l2);
		}

		if(mtr_has_sum)
			mtr_has_prm=1;
		
		remake_2(a1,l1,ind1,ind2,(int)d1,(int)d2);
/*
		WriteTerm(a1);
		puts("\n");
*/		
		*tp=2+mtr_has_prm+4*mtr_has_sum;
		return l1;
	}
	else
	{
		List l1=0,l2,l3;
		

		for(k=1;k<=d3;k++)
		{
			l2=0;
			for(j=1;j<=d2;j++)
			{
				l3=0;
				for(i=1;i<=d1;i++)
				{
					Term res;

					a11=CopyTerm(a1);

					sub_ind_alg(a11,ind1,i);
					sub_ind_alg(a11,ind2,j);
					sub_ind_alg(a11,ind3,k);

					res=a1_2_pl(a11);

					if(ListLength(res)>1)
						mtr_has_sum=1;
					if(res && ListLength(ListFirst(res))>1)
						mtr_has_prm=1;

					FreeAtomic(a11);
					l3=AppendLast(l3,res);
				}
				l2=AppendLast(l2,l3);
			}
			l1=AppendLast(l1,l2);
		}

		if(mtr_has_sum)
			mtr_has_prm=1;
		
		remake_3(a1,l1,ind1,ind2,ind3,(int)d1,(int)d2,(int)d3);
		
		*tp=2+mtr_has_prm+4*mtr_has_sum;
		return l1;
	}
	

	return 0;
}
