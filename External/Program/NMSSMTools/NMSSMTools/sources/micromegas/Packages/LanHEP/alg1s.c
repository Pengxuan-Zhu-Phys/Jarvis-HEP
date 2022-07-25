#include "lanhep.h"
#include <string.h>

int opSplitCol1=2;

extern int VerbMode;

static int alg1_cw_poly(Term a1);
static int alg1_cw_mono(Term m1);

static int alg1_sw_poly(Term a1);
static int alg1_sw_mono(Term m1);

static int alg1_fw_poly(Term a1);
static int alg1_fw_mono(Term m1);

extern int VerbMode;

static int acmp(Atom a, Atom b)
{
	return strcmp(AtomValue(a), AtomValue(b));
}

static List mk_ind(Term inc, Term ins, Term inv1, Term inv2)
{
	List ret=0;
	if(ins)
		ret=MakeList1(CopyTerm(ins));
	if(inv1)
		ret=AppendLast(ret,CopyTerm(inv1));
	if(inv2)
		ret=AppendLast(ret,CopyTerm(inv2));
	if(inc)
		ret=AppendLast(ret,CopyTerm(inc));
	return ret;
}

static Term mk_im_f0(int apcol, Term inc, Term inv1, Term inv2)
{
	Atom c;
	Term t;
	t=MakeCompound(A_I,3);
	if(inv1==0 && inv2==0)
		SetCompoundArg(t,1,NewInteger(0));
	else if(inv2==0)
		SetCompoundArg(t,1,NewInteger(1));
	else
		SetCompoundArg(t,1,NewInteger(2));
	SetCompoundArg(t,2,apcol?NewInteger(8):NewInteger(1));
	SetCompoundArg(t,3,NewInteger(0));
	c=ProcImPrt(t,0);
	return MakeCompound2(OPR_FIELD,mk_ind(inc,0,inv1,inv2),c);
}

static Term mk_im_f(int apdim, int apcol, Term inw, Term inc,
		Term inv1, Term inv2)
{
	Term a,m;
	List al=0;
	int i;
	if(apdim==0)
		return mk_im_f0(apcol,inc, inv1, inv2);
	
	for(i=0;i<apdim;i++)
	{
		m=MakeCompound(A_MTERM,4);
		SetCompoundArg(m,1,NewInteger(1));
		SetCompoundArg(m,2,NewInteger(1));
		SetCompoundArg(m,3,MakeList1(mk_im_f0(apcol,inc,inv1,inv2)));
		a=MakeCompound2(A_ALG1,MakeList1(m),mk_ind(inc,0,inv1,inv2));
		al=AppendLast(al,a);
	}
	
	m=MakeCompound(OPR_WILD,2);
	SetCompoundArg(m,2,al);
	al=AppendLast(mk_ind(inc,0,inv1,inv2),CopyTerm(inw));
	SetCompoundArg(m,1,al);
	return m;
}

static void mk_im_f_20(int apcol, Term inc1, Term inc2, Term *ap1, Term *ap2)
{
	Atom c,cc;
	Term t;
	t=MakeCompound(A_I,3);
	SetCompoundArg(t,1,NewInteger(0));
	SetCompoundArg(t,2,apcol?NewInteger(apcol):NewInteger(1));
	SetCompoundArg(t,3,NewInteger(1));
	c=ProcImPrt(t,0);
	cc=/*ProcImPrt(A_I,0);*/ GetAtomProperty(c,A_ANTI);
	
	*ap1=MakeCompound2(OPR_FIELD,apcol?MakeList1(CopyTerm(inc1)):0,c);
	*ap2=MakeCompound2(OPR_FIELD,apcol?MakeList1(CopyTerm(inc2)):0,cc);
}

static void mk_im_f_2(int apdim, int apcol, Term inw1, Term inc1,
		Term inw2, Term inc2, Term *ap1, Term *ap2)
{
	Term a,m;
	List al1=0, al2=0;
	int i;
	
	if(apcol<0)
	{
		Term *tmp;
		Term tmp1;
		tmp=ap1;
		ap1=ap2;
		ap2=tmp;
		tmp1=inc1;
		inc1=inc2;
		inc2=tmp1;
		apcol=-apcol;
	}
	
	if(apdim==0)
	{
		mk_im_f_20(apcol, inc1, inc2, ap1, ap2);
		return;
	}
	
	for(i=0;i<apdim;i++)
	{
		Term f1,f2;
		mk_im_f_20(apcol, inc1, inc2, &f1, &f2);
		m=MakeCompound(A_MTERM,4);
		SetCompoundArg(m,1,NewInteger(1));
		SetCompoundArg(m,2,NewInteger(1));
		SetCompoundArg(m,3,MakeList1(f1));
		a=MakeCompound2(A_ALG1,MakeList1(m),0);
		if(apcol)
			SetCompoundArg(a,2,MakeList1(CopyTerm(inc1)));
		al1=AppendLast(al1,a);
		m=MakeCompound(A_MTERM,4);
		SetCompoundArg(m,1,NewInteger(1));
		SetCompoundArg(m,2,NewInteger(1));
		SetCompoundArg(m,3,MakeList1(f2));
		a=MakeCompound2(A_ALG1,MakeList1(m),0);
		if(apcol)
			SetCompoundArg(a,2,MakeList1(CopyTerm(inc2)));
		al2=AppendLast(al2,a);
	}
	
	m=MakeCompound(OPR_WILD,2);
	SetCompoundArg(m,2,al1);
	if(apcol)
		al1=MakeList1(CopyTerm(inc1));
	else
		al1=NewList();
	al1=AppendLast(al1,CopyTerm(inw1));
	SetCompoundArg(m,1,al1);
	*ap1=m;
	
	m=MakeCompound(OPR_WILD,2);
	SetCompoundArg(m,2,al2);
	if(apcol)
		al1=MakeList1(CopyTerm(inc2));
	else
		al1=NewList();
	al1=AppendLast(al1,CopyTerm(inw2));
	SetCompoundArg(m,1,al1);
	*ap2=m;
	
	return ;
}

List alg1_spl_col2(List a1)
{
	
	List ar=0, al=0, ltc=0;
	List l1,l2;
	
	
	for(l1=a1;l1;l1=ListTail(l1))
	{

		long int nu=0, de=0, apdim=0, apcol=0, need_sq2;	
		Term inw=0, inc=0, inv=0, inv2=0, ins=0;
		List in1=0, in2=0;
		Atom Maux;
		Atom le1, le2;
		Term le;
		
		if(CompoundArgN(ListFirst(l1),4))
			continue;
		
		l2=CompoundArgN(ListFirst(l1),3);
		if(ListLength(l2)!=2)
			continue;
		if(CompoundName(ListFirst(l2))!=OPR_LET || 
				CompoundName(ListFirst(ListTail(l2)))!=OPR_LET)
			continue;

		le1=CompoundArg2(ListFirst(l2));
		le2=CompoundArg2(ListFirst(ListTail(l2)));

		if(GetAtomProperty(le1,A_COLOR)!=NewInteger(2))
			continue;
		if(GetAtomProperty(le1,OPR_FIELD)!=NewInteger(2))
			continue;
		if(GetAtomProperty(le2,A_COLOR)!=NewInteger(2))
			continue;
		if(GetAtomProperty(le2,OPR_FIELD)!=NewInteger(2))
			continue;
		if(le1!=le2 && le1!=GetAtomProperty(l2,A_ANTI))
			continue;

		nu=IntegerValue(CompoundArg1(ListFirst(l1)));
		de=IntegerValue(CompoundArg2(ListFirst(l1)));
		in1=CopyTerm(CompoundArg1(ListFirst(l2)));
		in2=CopyTerm(CompoundArg1(ListFirst(ListTail(l2))));

		
		if((nu!=1 && nu!=-1) || (de!=1 && de!=2 && de!=4))
			continue;

		for(l2=in1;l2;l2=ListTail(l2))
		{
			if(CompoundName(ListFirst(l2))==OPR_WILD)
			{
				if(apdim)
					break;
				apdim=IntegerValue(CompoundArgN(ListFirst(l2),3));
				inw=ListFirst(l2);
			}
			else if(CompoundName(CompoundArg1(ListFirst(l2)))==A_COLOR)
			{
				if(apcol)
					break;
				if( CompoundArg1(CompoundArg1(ListFirst(l2)))!=
					CompoundArg2(CompoundArg1(ListFirst(l2))))
					apcol=3;
				else
					apcol=8;
				inc=ListFirst(l2);
			}
			else if(CompoundName(CompoundArg1(ListFirst(l2)))==A_LORENTZ && 
					CompoundArg1(CompoundArg1(ListFirst(l2)))!=
					CompoundArg2(CompoundArg1(ListFirst(l2))))
			{
				if(ins)
					break;
				ins=ListFirst(l2);
			}
			else if(CompoundName(CompoundArg1(ListFirst(l2)))==A_LORENTZ && 
					CompoundArg1(CompoundArg1(ListFirst(l2)))==
					CompoundArg2(CompoundArg1(ListFirst(l2))))
			{
				if(inv && inv2)
					break;
				if(!inv) inv=ListFirst(l2);
				else inv2=ListFirst(l2);
			}
			else
			{
				break;
			}
		}
		
		if(l2)
			continue;
		if(ins && inv2)
			continue;
		if(le1!=le2 || ins)
			continue;
	
		ltc=AppendLast(ltc,l1);
	
	


		
		if(nu==-1)
		{
			nu=1;
			if(de<0)
			{
				de=-de;
				nu=-nu;
			}
			else
				de=-de;
		}

		if(de==4 || de==-4) 
			de/=2;
		else if(de==2 || de==-2)
		{
			de/=2;
			need_sq2=0;
		}
		else
			need_sq2=1;

		if(nu<0)
			nu=-nu;
	
		ar=mk_im_f((int)apdim,(int)apcol,inw,inc,inv,inv2);
		Maux=NewAtom("Maux",0);
	
		l2=MakeList1(ar);
		l2=AppendLast(l2,MakeCompound2(OPR_LET,CopyTerm(in1),le1));
		if(inv2==0)
			l2=AppendLast(l2,MakeCompound2(OPR_PARAMETER,0,Maux));
		if(need_sq2)
			l2=AppendLast(l2,MakeCompound2(OPR_PARAMETER,0,A_SQRT2));
		
		le=MakeCompound(A_MTERM,4);
		SetCompoundArg(le,1,NewInteger(nu));
		SetCompoundArg(le,2,NewInteger(de));
		SetCompoundArg(le,3,l2);
		
		a1=AppendLast(a1,le);
	
	
		if(VerbMode)
		{
			WarningInfo(214);
			printf("%s*%s ",AtomValue(le1),AtomValue(le2));
			puts("4-colored multiplet product splitted");
		}
	}
	

	for(l1=ltc;l1;l1=ListTail(l1))
		a1=CutFromList(a1,ListFirst(l1));
	RemoveList(ltc);

	
	return a1;
	
	
}


List alg1_spl_col(List a1)
{
	List ar=0, al=0, ltc=0;
	List l1,l2;
	long int nu=0, de=0;
	int apdim=0, apcol=0, need_sq2;
	Term inw=0, inc=0;
	List in1=0, in2=0;
	Atom Maux;
	
	if(opSplitCol1==0)
		return a1;
	
	if(opSplitCol1==-1)
	{
		for(l1=a1;l1;l1=ListTail(l1))
		{
			int cw=0;
			for(l2=CompoundArgN(ListFirst(l1),3);l2;l2=ListTail(l2))
			{
				if(CompoundName(ListFirst(l2))==OPR_FIELD &&
						GetAtomProperty(CompoundArg2(ListFirst(l2)),A_COLOR))
				{
					cw++;
					continue;
				}
				if(CompoundName(ListFirst(l2))==OPR_LET)
				{
					Term pp;
					pp=GetAtomProperty(CompoundArg2(ListFirst(l2)),A_COLOR);
					if(pp==0)
						continue;
					if(!is_integer(pp))
					{
						cw++;
						continue;
					}
					if(IntegerValue(pp)<0)
					{
						cw-=1000;
						continue;
					}
					cw+=(int)IntegerValue(pp);
					continue;
				}
			}
			if(cw>3)
				ltc=AppendLast(ltc,l1);
		}
		
		for(l1=ltc;l1;l1=ListTail(l1))
			a1=CutFromList(a1,ListFirst(l1));
		RemoveList(ltc);
		return a1;
	}
	
				
	
	for(l1=a1;l1;l1=ListTail(l1))
	{
		Atom le1, le2;
		
		if(CompoundArgN(ListFirst(l1),4))
			continue;
		
		l2=CompoundArgN(ListFirst(l1),3);
		if(ListLength(l2)!=2)
			continue;
		if(CompoundName(ListFirst(l2))!=OPR_LET || 
				CompoundName(ListFirst(ListTail(l2)))!=OPR_LET)
			continue;

		le1=CompoundArg2(ListFirst(l2));
		le2=CompoundArg2(ListFirst(ListTail(l2)));
				
		if(GetAtomProperty(le1,opSplitCol1==1?A_COLOR:OPR_SCALAR)!=
															NewInteger(2))
			continue;
		if(GetAtomProperty(le2,opSplitCol1==1?A_COLOR:OPR_SCALAR)!=
															NewInteger(2))
			continue;
	
		if(ltc==0)
		{
			if(!ListMember(al,le1))
				al=AppendLast(al,le1);
			if(!ListMember(ar,le2))
				ar=AppendLast(ar,le2);
			ltc=AppendLast(ltc,l1);
			nu=IntegerValue(CompoundArg1(ListFirst(l1)));
			de=IntegerValue(CompoundArg2(ListFirst(l1)));
			in1=CopyTerm(CompoundArg1(ListFirst(l2)));
			in2=CopyTerm(CompoundArg1(ListFirst(ListTail(l2))));
		}
		else
		{
			int bad=0;
			List l3,l4;
			if(IntegerValue(CompoundArg1(ListFirst(l1)))!=nu)
				bad=1;
			if(IntegerValue(CompoundArg2(ListFirst(l1)))!=de)
				bad=1;
			for(l3=in1,l4=CompoundArg1(ListFirst(l2));l3&&l4;
							l3=ListTail(l3),l4=ListTail(l4))
					if(!EqualTerms(CompoundArg1(ListFirst(l3)),
									CompoundArg1(ListFirst(l3))))
				bad=1;
			if(l3 || l4)
				bad=1;
			for(l3=in2,l4=CompoundArg1(ListFirst(ListTail(l2)));l3&&l4;
							l3=ListTail(l3),l4=ListTail(l4))
					if(!EqualTerms(CompoundArg1(ListFirst(l3)),
									CompoundArg1(ListFirst(l3))))
				bad=1;
			if(l3 || l4)
				bad=1;
			if(bad)
			{
				RemoveList(ltc);
				RemoveList(al);
				RemoveList(ar);
				FreeAtomic(in1);
				FreeAtomic(in2);
				a1=alg1_spl_col2(a1);	
				return a1;
			}

			if(!ListMember(al,le1))
				al=AppendLast(al,le1);
			if(!ListMember(ar,le2))
				ar=AppendLast(ar,le2);
			ltc=AppendLast(ltc,l1);
		}
	}
	
	if(is_empty_list(ltc))
	{
		a1=alg1_spl_col2(a1);	
		return a1;
	}
	
	al=SortedList(al,acmp);
	ar=SortedList(ar,acmp);
		
	if(ListLength(al)!=ListLength(ar) || !EqualTerms(al,ar) ||
			ListLength(al)*ListLength(al)!=ListLength(ltc) ||
			(nu!=1 && nu!=-1) || (de!=1 && de!=2))
	{
		RemoveList(ltc);
		RemoveList(al);
		RemoveList(ar);
		FreeAtomic(in1);
		FreeAtomic(in2);
		a1=alg1_spl_col2(a1);	
		return a1;
	}
	
	for(l1=in1;l1;l1=ListTail(l1))
	{
		if(CompoundName(ListFirst(l1))==OPR_WILD)
		{
			if(apdim)
			{
				RemoveList(ltc);
				RemoveList(al);
				RemoveList(ar);
				FreeAtomic(in1);
				FreeAtomic(in2);
				return a1;
			}
			apdim=(int)IntegerValue(CompoundArgN(ListFirst(l1),3));
			inw=ListFirst(l1);
		}
		else
		{
			if(apcol || CompoundName(CompoundArg1(ListFirst(l1)))!=A_COLOR ||
					CompoundArg1(CompoundArg1(ListFirst(l1)))!=
					CompoundArg2(CompoundArg1(ListFirst(l1))))
			{
				RemoveList(ltc);
				RemoveList(al);
				RemoveList(ar);
				FreeAtomic(in1);
				FreeAtomic(in2);
				a1=alg1_spl_col2(a1);	
				return a1;
			}
			apcol=8;
			inc=ListFirst(l1);
		}
	}
	
	for(l1=ltc;l1;l1=ListTail(l1))
		a1=CutFromList(a1,ListFirst(l1));
	RemoveList(ltc);
	RemoveList(ar);
	FreeAtomic(in2);

	for(l1=a1;l1;l1=ListTail(l1))
	{
		Atom le1, le2;
		Term mm;
		
		mm=ListFirst(l1);
		
		l2=CompoundArgN(mm,3);
		if(ListLength(l2)!=2)
			continue;
		if(CompoundName(ListFirst(l2))!=OPR_LET || 
				CompoundName(ListFirst(ListTail(l2)))!=OPR_LET)
			continue;
		le1=CompoundArg2(ListFirst(l2));
		le2=CompoundArg2(ListFirst(ListTail(l2)));
		
		if(le1==le2)
			continue;
				
		for(l2=a1;l2!=l1;l2=ListTail(l2))
		{
			Atom lee1, lee2;
			Term mm1;
			List l22;
			mm1=ListFirst(l2);

			l22=CompoundArgN(mm1,3);
			if(ListLength(l22)!=2)
				continue;
			if(CompoundName(ListFirst(l22))!=OPR_LET || 
					CompoundName(ListFirst(ListTail(l22)))!=OPR_LET)
				continue;
			lee1=CompoundArg2(ListFirst(l22));
			lee2=CompoundArg2(ListFirst(ListTail(l22)));

			if(lee1==lee2 || lee1!=le2 || lee2!=le1)
				continue;
			if(!EqualTerms(CompoundArg1(mm),CompoundArg1(mm1)))
				continue;
			if(!EqualTerms(CompoundArg2(mm),CompoundArg2(mm1)))
				continue;
			if(!EqualTerms(CompoundArgN(mm,4),CompoundArgN(mm1,4)))
				continue;
			
			break;
		}

		if(l2!=l1)
		{
			long int nu1, de1;

			a1=CutFromList(a1,l2);

			nu1=IntegerValue(CompoundArg1(ListFirst(l1)));
			de1=IntegerValue(CompoundArg2(ListFirst(l1)));
			if(de1/2==de1-de1/2)
				de1/=2;
			else
				nu1*=2;
			SetCompoundArg(ListFirst(l1),1,NewInteger(nu1));
			SetCompoundArg(ListFirst(l1),2,NewInteger(de1));
		}
		

	}
		
	if(nu==-1)
	{
		nu=1;
		if(de<0)
		{
			de=-de;
			nu=-nu;
		}
		else
			de=-de;
	}
	
/*	if(de<0)
		nu=-nu;
	de=-de;
*/	
	if(de==2 || de==-2)
	{
		de/=2;
		need_sq2=0;
	}
	else
		need_sq2=1;
	
	if(nu<0)
		nu=-nu;
	
	ar=mk_im_f(apdim,apcol,inw,inc,0,0);
	Maux=NewAtom("Maux",0);
	
	for(l1=al;l1;l1=ListTail(l1))
	{
		Term le;
		
		le=ListFirst(l1);
		l2=MakeList1(CopyTerm(ar));
		l2=AppendLast(l2,MakeCompound2(OPR_LET,CopyTerm(in1),le));
		l2=AppendLast(l2,MakeCompound2(OPR_PARAMETER,0,Maux));
		if(need_sq2)
			l2=AppendLast(l2,MakeCompound2(OPR_PARAMETER,0,A_SQRT2));
		
		le=MakeCompound(A_MTERM,4);
		SetCompoundArg(le,1,NewInteger(nu));
		SetCompoundArg(le,2,NewInteger(de));
		SetCompoundArg(le,3,l2);
		a1=AppendLast(a1,le);
	}
	
	if(VerbMode==1)
	{
		WarningInfo(214);
		puts("4 colored multiplets splitted; use '-vv' option to see result");
	} 
	
	FreeAtomic(in1);
	FreeAtomic(ar);
	RemoveList(al);
	
	return a1;
	
}

extern int allow_dfdfc;
extern List fromdfdfc;

Term ProcDFDFC(Term t, Term ind)
{
	Term t1;
	List la, lb, l1, l2, l3;
	int need_crdc=0;
	int apcol=0, apdim=0;
	Term inc1=0, inw1=0, inc2=0, inw2=0;
	Term ap1=0, ap2=0, Maux=0;
	List ret=0;
	
	
	t=ProcDF(t,0);
	if(t==0)
		return 0;
	
	t1=CopyTerm(t);
	alg1_anti(t1);
	
	for(la=0,l1=CompoundArg2(t1);l1;l1=ListTail(l1))
		la=AppendLast(la,CompoundArg2(ListFirst(l1)));
	
	for(lb=0,l1=CompoundArg1(t1);l1;l1=ListTail(l1))
	for(l2=CompoundArgN(ListFirst(l1),3);l2;l2=ListTail(l2))
	for(l3=CompoundArg1(ListFirst(l2));l3;l3=ListTail(l3))
		if(!ListMember(la,CompoundArg2(ListFirst(l3))) &&
			!ListMember(lb,CompoundArg2(ListFirst(l3))))
				lb=AppendLast(lb,CompoundArg2(ListFirst(l3)));
	
	RemoveList(la);
	for(la=0,l1=lb;l1;l1=ListTail(l1))
		la=AppendLast(la,NewLabel());
	
	for(l1=CompoundArg1(t1);l1;l1=ListTail(l1))
	for(l2=CompoundArgN(ListFirst(l1),3);l2;l2=ListTail(l2))
	for(l3=CompoundArg1(ListFirst(l2));l3;l3=ListTail(l3))
	{
		Label l;
		int pos;
		l=CompoundArg2(ListFirst(l3));
		pos=ListMember(lb,l);
		if(pos)
		{
			l=ListNth(la,pos);
			SetCompoundArg(ListFirst(l3),2,l);
		}
	}
		
		
	RemoveList(la);
	RemoveList(lb);
	
	if(opSplitCol1<1 || !allow_dfdfc)
		need_crdc=0;
	else if(opSplitCol1==1)
	{
		for(l1=CompoundArg1(t);l1;l1=ListTail(l1))
			if(alg1_cw_mono(ListFirst(l1))==2)
				need_crdc=1;
	}
	else if(opSplitCol1==2)
	{
		for(l1=CompoundArg1(t);l1;l1=ListTail(l1))
			if(alg1_sw_mono(ListFirst(l1))==2)
				need_crdc=1;
	}
	
	if(need_crdc)
	for(l1=CompoundArg2(t);l1;l1=ListTail(l1))
	{
		if(CompoundName(ListFirst(l1))==OPR_WILD)
		{
			if(apdim)
			{
				need_crdc=0;
				break;
			}
			apdim=(int)IntegerValue(CompoundArgN(ListFirst(l1),3));
			inw1=CopyTerm(ListFirst(l1));
		}
		else
		{
			char *rbuf;
			if(apcol || CompoundName(CompoundArg1(ListFirst(l1)))!=A_COLOR )
			{
				need_crdc=0;
				break;
			}
			inc1=CopyTerm(ListFirst(l1));
			rbuf=AtomValue(CompoundArg1(CompoundArg1(inc1)));
			if(rbuf[1]=='8')
				apcol=8;
			else if(rbuf[1]=='3' && rbuf[2]==0)
				apcol=3;
			else if(rbuf[1]=='3' && rbuf[2]=='b')
				apcol=-3;
			else
			{
				printf("Internal error (a1s01:'%s')\n",rbuf);
				need_crdc=0;
				break;
			}
		}
	}
	
	if(need_crdc)
	for(l1=CompoundArg2(t1);l1;l1=ListTail(l1))
	{
		if(CompoundName(ListFirst(l1))==OPR_WILD)
			inw2=CopyTerm(ListFirst(l1));
		else
			inc2=CopyTerm(ListFirst(l1));
	}
	
	for(l1=CompoundArg1(t);l1;l1=ListTail(l1))
	for(l2=CompoundArg1(t1);l2;l2=ListTail(l2))
	{
		Term m, m1, m2;
		long int n,n1,n2,d,d1,d2,c;
		
		m1=ListFirst(l1);
		m2=ListFirst(l2);
		
		if(opSplitCol1==-1 && alg1_cw_mono(m1)==2 
				&& alg1_cw_mono(m2)==2)
			continue;
		
		if(opSplitCol1==1 && need_crdc && alg1_cw_mono(m1)==2 
				&& alg1_cw_mono(m2)==2)
			continue;
		
		if(opSplitCol1==2 && need_crdc && alg1_sw_mono(m1)==2 
				&& alg1_sw_mono(m2)==2)
			continue;
		
		n1=IntegerValue(CompoundArg1(m1));
		d1=IntegerValue(CompoundArg2(m1));
		n2=IntegerValue(CompoundArg1(m2));
		d2=IntegerValue(CompoundArg2(m2));
		n=n1*n2;
		d=d1*d2;
		if(d1<0 && d2<0)
			n=-n;
		c=gcf(n,d);
		n/=c;
		d/=c;
		
		m=MakeCompound(A_MTERM,4);
		SetCompoundArg(m,1,NewInteger(n));
		SetCompoundArg(m,2,NewInteger(d));
		SetCompoundArg(m,3,ConcatList(CopyTerm(CompoundArgN(m1,3)),
				CopyTerm(CompoundArgN(m2,3))));
		SetCompoundArg(m,4,ConcatList(CopyTerm(CompoundArgN(m1,4)),
				CopyTerm(CompoundArgN(m2,4))));
		ret=AppendLast(ret,m);
	}
	

		
	if(need_crdc)
	{
		Maux=NewAtom("Maux",0);
		mk_im_f_2(apdim,apcol,inw1,inc1,inw2,inc2,&ap2,&ap1);
	
		for(l1=CompoundArg1(t);l1;l1=ListTail(l1))
		{
			Term m;
			m=ListFirst(l1);
			
			if((opSplitCol1==1 && alg1_cw_mono(m)!=2) ||
					(opSplitCol1==2 && alg1_sw_mono(m)!=2))
				continue;
			
			ChangeList(l1,0);
			
			l2=ConsumeCompoundArg(m,3);
			l2=AppendLast(l2,CopyTerm(ap1));
			l2=AppendLast(l2,MakeCompound2(OPR_PARAMETER,0,Maux));
			SetCompoundArg(m,3,l2);
/*			
			n=IntegerValue(CompoundArg1(m));
			d=IntegerValue(CompoundArg2(m));
			
			if(d<0)
				n=-n;
			d=-d;
			
			SetCompoundArg(m,1,NewInteger(n));
			SetCompoundArg(m,2,NewInteger(d));
*/			
			ret=AppendLast(ret,m);
		}
		
		for(l1=CompoundArg1(t1);l1;l1=ListTail(l1))
		{
			Term m;
			m=ListFirst(l1);
			
			if((opSplitCol1==1 && alg1_cw_mono(m)!=2) ||
					(opSplitCol1==2 && alg1_sw_mono(m)!=2))
				continue;
			
			ChangeList(l1,0);
			
			l2=ConsumeCompoundArg(m,3);
			l2=AppendLast(l2,CopyTerm(ap2));
			l2=AppendLast(l2,MakeCompound2(OPR_PARAMETER,0,Maux));
			SetCompoundArg(m,3,l2);
/*			
			n=IntegerValue(CompoundArg1(m));
			d=IntegerValue(CompoundArg2(m));
			
			if(d<0)
				n=-n;
			d=-d;
			
			SetCompoundArg(m,1,NewInteger(n));
			SetCompoundArg(m,2,NewInteger(d));
*/			
			fromdfdfc=AppendLast(fromdfdfc,m);
		}
		
		FreeAtomic(ap1);
		FreeAtomic(ap2);
		FreeAtomic(inc1);
		FreeAtomic(inc2);
		FreeAtomic(inw1);
		FreeAtomic(inw2);
	}
	
	FreeAtomic(t);
	FreeAtomic(t1);
		
	return MakeCompound2(A_ALG1,ret,0);
	
}


static int alg1_cw_mono(Term m1)
{
	int ret=0;
	
	List l1,l2;
	
	for(l1=CompoundArgN(m1,3);l1;l1=ListTail(l1))
	{
		if(CompoundName(ListFirst(l1))==OPR_PARAMETER &&
				(CompoundArg2(ListFirst(l1))==A_SIN||
				 CompoundArg2(ListFirst(l1))==A_COS))
		{
			return -1;
		}
		
		if(CompoundName(ListFirst(l1))==OPR_FIELD &&
				GetAtomProperty(CompoundArg2(ListFirst(l1)),A_COLOR))
		{
			ret++;
			continue;
		}
		
		if(CompoundName(ListFirst(l1))==OPR_WILD)
		{
			int rt1;
			List al;
			al=CompoundArg2(ListFirst(l1));
			if(is_empty_list(al))
				continue;
			rt1=alg1_cw_poly(ListFirst(al));
			if(rt1==-1)
				return -1;
			for(l2=ListTail(al);l2;l2=ListTail(l2))
				if(rt1!=alg1_cw_poly(ListFirst(l2)))
					return -1;
			ret+=rt1;
			continue;
		}
		
		if(CompoundName(ListFirst(l1))==OPR_LET)
		{
			Term prop;
			prop=GetAtomProperty(CompoundArg2(ListFirst(l1)),A_COLOR);
			if(prop && !is_integer(prop))
				ret++;
			else
				if(is_integer(prop))
				{
					if(IntegerValue(prop)==-1)
						return -1;
					ret+=(int)IntegerValue(prop);
				}
		}
	}
	
	return ret;
}

static int alg1_sw_mono(Term m1)
{
	int ret=0;
	
	List l1,l2;
	
	for(l1=CompoundArgN(m1,3);l1;l1=ListTail(l1))
	{
		if(CompoundName(ListFirst(l1))==OPR_PARAMETER &&
				(CompoundArg2(ListFirst(l1))==A_SIN||
				 CompoundArg2(ListFirst(l1))==A_COS))
		{
			return -1;
		}
		if(CompoundName(ListFirst(l1))==OPR_FIELD)
		{
			Term prp;
			if(CompoundArg2(ListFirst(l1))==A_VEV)
				return -1;
			prp=GetAtomProperty(CompoundArg2(ListFirst(l1)),PROP_INDEX);
			for(l2=prp;l2;l2=ListTail(l2))
				if(CompoundName(CompoundArg1(ListFirst(l2)))==A_LORENTZ)
					return -1;
			ret++;
			continue;
		}
		
		if(CompoundName(ListFirst(l1))==OPR_WILD)
		{
			int rt1;
			List al;
			al=CompoundArg2(ListFirst(l1));
			if(is_empty_list(al))
				continue;
			rt1=alg1_sw_poly(ListFirst(al));
			if(rt1==-1)
				return -1;
			for(l2=ListTail(al);l2;l2=ListTail(l2))
				if(rt1!=alg1_sw_poly(ListFirst(l2)))
					return -1;
			ret+=rt1;
			continue;
		}
		
		if(CompoundName(ListFirst(l1))==OPR_LET)
		{
			Term prop;
			prop=GetAtomProperty(CompoundArg2(ListFirst(l1)),OPR_SCALAR);
			if(IntegerValue(prop)==-1)
				return -1;
			ret+=(int)IntegerValue(prop);
		}
	}
	
	return ret;
}

static int alg1_fw_mono(Term m1)
{
	int ret=0;
	
	List l1,l2;
	
	for(l1=CompoundArgN(m1,3);l1;l1=ListTail(l1))
	{
		if(CompoundName(ListFirst(l1))==OPR_PARAMETER &&
				(CompoundArg2(ListFirst(l1))==A_SIN||
				 CompoundArg2(ListFirst(l1))==A_COS))
		{
			return -1;
		}
		if(CompoundName(ListFirst(l1))==OPR_FIELD)
		{
			Term prp;
			if(CompoundArg2(ListFirst(l1))==A_VEV)
				return -1;
			ret++;
			continue;
		}
		
		if(CompoundName(ListFirst(l1))==OPR_WILD)
		{
			int rt1;
			List al;
			al=CompoundArg2(ListFirst(l1));
			if(is_empty_list(al))
				continue;
			rt1=alg1_fw_poly(ListFirst(al));
			if(rt1==-1)
				return -1;
			for(l2=ListTail(al);l2;l2=ListTail(l2))
				if(rt1!=alg1_fw_poly(ListFirst(l2)))
					return -1;
			ret+=rt1;
			continue;
		}
		
		if(CompoundName(ListFirst(l1))==OPR_LET)
		{
			Term prop;
			prop=GetAtomProperty(CompoundArg2(ListFirst(l1)),OPR_FIELD);
			if(IntegerValue(prop)==-1)
				return -1;
			ret+=(int)IntegerValue(prop);
		}
	}
	
	return ret;
}


static int alg1_cw_poly(Term a1)
{
	int ret=0;
	List l1;
	
	l1=CompoundArg1(a1);
	if(is_empty_list(l1))
		return 0;
	ret=alg1_cw_mono(ListFirst(l1));
	
	for(l1=ListTail(l1);l1;l1=ListTail(l1))
		if(alg1_cw_mono(ListFirst(l1))!=ret)
			return -1;
	
	return ret;
}

static int alg1_sw_poly(Term a1)
{
	int ret=0;
	List l1;
	
	l1=CompoundArg1(a1);
	if(is_empty_list(l1))
		return 0;
	ret=alg1_sw_mono(ListFirst(l1));
	
	for(l1=ListTail(l1);l1;l1=ListTail(l1))
		if(alg1_sw_mono(ListFirst(l1))!=ret)
			return -1;
	
	return ret;
}

static int alg1_fw_poly(Term a1)
{
	int ret=0;
	List l1;
	
	l1=CompoundArg1(a1);
	if(is_empty_list(l1))
		return 0;
	ret=alg1_fw_mono(ListFirst(l1));
	
	for(l1=ListTail(l1);l1;l1=ListTail(l1))
		if(alg1_fw_mono(ListFirst(l1))!=ret)
			return -1;
	
	return ret;
}

void alg1_let_sw(Atom let)
{
	Term prp;
	int cw;
	prp=GetAtomProperty(let,PROP_TYPE);
	cw=alg1_sw_poly(CompoundArg1(prp));

/*if(cw==2) printf("scalar weight of '%s' is %d\n",AtomValue(let),cw);*/

	SetAtomProperty(let,OPR_SCALAR,NewInteger(cw));
}

void alg1_let_fw(Atom let)
{
	Term prp;
	int cw;
	prp=GetAtomProperty(let,PROP_TYPE);
	cw=alg1_fw_poly(CompoundArg1(prp));

/*if(cw==2) printf("field weight of '%s' is %d\n",AtomValue(let),cw);*/

	SetAtomProperty(let,OPR_FIELD,NewInteger(cw));
}


void alg1_let_cw(Atom let)
{
	Term prp;
	int cw;
	prp=GetAtomProperty(let,PROP_TYPE);
	cw=alg1_cw_poly(CompoundArg1(prp));

/*if(cw==2) printf("color weight of '%s' is %d\n",AtomValue(let),cw);*/

	SetAtomProperty(let,A_COLOR,NewInteger(cw));
	alg1_let_sw(let);
	alg1_let_fw(let);
}

void alg1_rem_c4(Term a1)
{
	List l;
	List ml;
	ml=ConsumeCompoundArg(a1,1);
	
rpt:
	for(l=ml;l;l=ListTail(l))
		if(alg1_cw_mono(ListFirst(l))>3)
		{
			ml=CutFromList(ml,l);
			goto rpt;
		}
		
	SetCompoundArg(a1,1,ml);
}


