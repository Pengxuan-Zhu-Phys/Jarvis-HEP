#include <string.h>
#include <stdio.h>
#include "lanhep.h"

List tri_si_co_list = 0;		
List tri_sc_2_list = 0;

extern int opUndefAngleComb;

int opNoDummies=0;
int tri_dbg_mode=0;
int opTriHeu=1;
int tri_third_angle=0;

List tri_find_sub(List);
List tri_cvt(List, Atom, Atom);
List tri_csf(List *);
void tri_set_sub(List);
void tri_red_numf(List, List, List, List);
static int mk_s2(List *);
void tri_set_rule(List, List, List);
static List list_vars(List ml, int do2);
static int red_1(int, List, List *, List *), red_2(int, List, List *, List *);
int  tri_divs(Term, int *, List *, List *, List);
int  tri_divf(Term, List *, List *, List);
int tri_opti_1(int, List *);
Atom tri_dummy_prm(List);

static int pcmp(Term p1, Term p2)
	{
	return strcmp(AtomValue(CompoundArg1(p1)),AtomValue(CompoundArg1(p2)));
	}


int tri_heu(List m2l)
{
	List l, lu, l1, l2, csf;
	
	if(opTriHeu==0)
		return 0;
	
	m2l=CopyTerm(m2l);


	if(CompoundArgN(ListFirst(m2l),3))
		for(l=m2l;l;l=ListTail(l))
			FreeAtomic(ConsumeCompoundArg(ListFirst(l),3));


	if(ListLength(m2l)==1)
	{
		List l1;
		int co;
		
		if(tri_dbg_mode)
			puts("TriHeuG: monomial");
		
		l1=CopyTerm(m2l);
		l=ConsumeCompoundArg(ListFirst(l1),2);
		co=mk_s2(&l);
		SetCompoundArg(ListFirst(l1),2,l);
		SetCompoundArg(ListFirst(m2l),1,NewInteger(co));
		SetCompoundArg(ListFirst(l1),1,NewInteger(1));
		tri_set_rule(0,m2l,l1);
		FreeAtomic(m2l);
		FreeAtomic(l1);
		return 1;
	}
	
	
	/*   First, remove common symbolic factor */
    
	csf=tri_csf(&m2l);
	
	if(ListLength(m2l)==1 && CompoundArg2(ListFirst(m2l))==0)
	{
		FreeAtomic(m2l);
		m2l=MakeCompound(A_MTERM,3);
		SetCompoundArg(m2l,1,NewInteger(1));
		SetCompoundArg(m2l,2,csf);
		m2l=MakeList1(m2l);
		if(tri_dbg_mode)
			puts("TriHeuG: expr is const after csf extraction");
		tri_set_rule(0,m2l,m2l);
		FreeAtomic(m2l);

		return 1;
	}
	
	
	l=tri_find_sub(m2l);
	if(l)
	{
		if(tri_dbg_mode)
			puts("TriHeuG: sub found after csf extraction");
		tri_set_rule(csf,m2l,l);
		FreeAtomic(m2l);
		FreeAtomic(l);
		return 1;
	}
	
	lu=list_vars(m2l,1);
	
	if(ListLength(lu)==1 && IntegerValue(ListFirst(lu))>0)
	{
		if(red_1((int)IntegerValue(ListFirst(lu)),m2l,&l1,&l2)==0)
		{
			return 0;
		}
	}

	if(ListLength(lu)==1 && IntegerValue(ListFirst(lu))<0)
	{

		if(red_2(-(int)IntegerValue(ListFirst(lu)),m2l,&l1,&l2)==0)
		{
			return 0;
		}

	}
	
	
	if(ListLength(lu)==1)
	{
		/*FreeAtomic(csf);FreeAtomic(l1);FreeAtomic(l2);*/
		
		
		/*{WriteTerm(csf);WriteTerm(l1);WriteTerm(l2);puts("");	}	*/
/*		puts("!!!");
		FreeAtomic(m2l);
		puts("!");
		FreeAtomic(l1);
		puts("!!");
		FreeAtomic(l2);*/
		tri_set_rule(csf,l1,l2);
		return 1;
		
	}
		
	/***********************************************************************
	 ***                                                                 ***
	 ***              2 unbiased angles                                  ***
	 ***                                                                 ***
	 ***********************************************************************/
	

	if(ListLength(lu)==2 && IntegerValue(ListFirst(lu))>0)
	{
		int ft, a1, a2;
		a1=(int) IntegerValue(ListFirst(ListTail(lu)));
		a2=(int) IntegerValue(ListFirst(lu));
		
		
		/**********************  additive case  **************************/
		
		if(tri_divs(
				ListNth(tri_si_co_list,a1),
				&ft,&l1,&l2,m2l))
		{
			List l1_r,l1_l,l2_r,l2_l;
			Term m2;
			if(red_1(a1,l1,&l1_l,&l1_r)==0)
				return 0;
			if(red_1(a2,l2,&l2_l,&l2_r)==0)
				return 0;

			FreeAtomic(l1);
			FreeAtomic(l2);
			
			if(ft)
			{
				m2=MakeCompound(A_MTERM,3);
				SetCompoundArg(m2,1,NewInteger(ft));
				l1=AppendLast(NewList(),m2);
			}
			else
				l1=NewList();
			l1=ConcatList(l1,l1_l);
			l1=ConcatList(l1,l2_l);
			
			if(ft)
			{
				m2=MakeCompound(A_MTERM,3);
				SetCompoundArg(m2,1,NewInteger(ft));
				l2=AppendLast(NewList(),m2);
			}
			else
				l2=NewList();
			
			l2=ConcatList(l2,l1_r);
			l2=ConcatList(l2,l2_r);
			
			if(tri_dbg_mode)
				puts("TriHeuG: additive rule (++)");
			tri_set_rule(csf,l1,l2);
			FreeAtomic(l1);
			FreeAtomic(l2);
			return 1;
		}
		
		/**************  multiplicative case *************************/
		
		if(tri_divf(
				ListNth(tri_si_co_list,a1),
				&l1,&l2,m2l))
		{
			List l1_r,l1_l,l2_r,l2_l;
			Term m2;
			if(red_1(a1,l1,&l1_l,&l1_r)==0)
				return 0;
			if(red_1(a2,l2,&l2_l,&l2_r)==0)
				return 0;

			FreeAtomic(l1);
			FreeAtomic(l2);
			
			tri_red_numf(l1_l,l1_r,0,0);
			tri_red_numf(l2_l,l2_r,0,0);
			
			m2=MakeCompound(A_MTERM,3);
			SetCompoundArg(m2,1,NewInteger(1));
			SetCompoundArg(m2,2,MakeList2(l1_l,l2_l));
			l1=AppendLast(NewList(),m2);

			m2=MakeCompound(A_MTERM,3);
			SetCompoundArg(m2,1,NewInteger(1));
			SetCompoundArg(m2,2,MakeList2(l1_r,l2_r));
			l2=AppendLast(NewList(),m2);
						
			if(tri_dbg_mode)
				puts("TriHeuG: multiplicative rule (++)");
			tri_set_rule(csf,l1,l2);
			FreeAtomic(l1);
			FreeAtomic(l2);
			return 1;
		}
		if(tri_dbg_mode)
			puts("TriHeu: pair (++) fail");
	}
	
	/***********************************************************************
	 ***                                                                 ***
	 ***              2 biased + 1 unbiased angles                       ***
	 ***                                                                 ***
	 ***********************************************************************/
	
	if(ListLength(lu)==2 && IntegerValue(ListFirst(lu))<0)
	{
		int ft, a1, a2;
		a1= (int)IntegerValue(ListFirst(ListTail(lu)));
		a2=-(int)IntegerValue(ListFirst(lu));
		
		/**************   additive case   ********************/
		
		if(tri_divs(
				ListNth(tri_si_co_list,a1),
				&ft,&l1,&l2,m2l))
		{
			List l1_r,l1_l,l2_r,l2_l;
			Term m2;
			if(red_1(a1,l1,&l1_l,&l1_r)==0)
				return 0;
			if(red_2(a2,l2,&l2_l,&l2_r)==0)
				return 0;

			FreeAtomic(l1);
			FreeAtomic(l2);
			
			if(ft)
			{
				m2=MakeCompound(A_MTERM,3);
				SetCompoundArg(m2,1,NewInteger(ft));
				l1=AppendLast(NewList(),m2);
			}
			else
				l1=NewList();
			l1=ConcatList(l1,l1_l);
			l1=ConcatList(l1,l2_l);
			
			if(ft)
			{
				m2=MakeCompound(A_MTERM,3);
				SetCompoundArg(m2,1,NewInteger(ft));
				l2=AppendLast(NewList(),m2);
			}
			else
				l2=NewList();
			
			l2=ConcatList(l2,l1_r);
			l2=ConcatList(l2,l2_r);
			
			if(tri_dbg_mode)
				puts("TriHeuG: additive rule (-+)");
			tri_set_rule(csf,l1,l2);
			FreeAtomic(l1);
			FreeAtomic(l2);
			return 1;
		}
		
		/**************   multiplicative case   ********************/
	    
		if(tri_divf(
				ListNth(tri_si_co_list,a1),
				&l2,&l1,m2l))
		{
			List l1_r,l1_l,l2_r,l2_l;
			Term m2;

			if(red_1(a1,l1,&l1_l,&l1_r)==0)
				return 0;

			if(red_2(a2,l2,&l2_l,&l2_r)==0)
				return 0;

			FreeAtomic(l1);
			FreeAtomic(l2);
			
			tri_red_numf(l1_l,l1_r,0,0);
			tri_red_numf(l2_l,l2_r,0,0);
			
			m2=MakeCompound(A_MTERM,3);
			SetCompoundArg(m2,1,NewInteger(1));
			SetCompoundArg(m2,2,MakeList2(l1_l,l2_l));
			l1=AppendLast(NewList(),m2);

			m2=MakeCompound(A_MTERM,3);
			SetCompoundArg(m2,1,NewInteger(1));
			SetCompoundArg(m2,2,MakeList2(l1_r,l2_r));
			l2=AppendLast(NewList(),m2);
						
			if(tri_dbg_mode)
				puts("TriHeuG: multiplicative rule (-+)");
			tri_set_rule(csf,l1,l2);
			FreeAtomic(l1);
			FreeAtomic(l2);
			return 1;
		}
		
		tri_third_angle=a1;
		
		if(red_2(a2,m2l,&l1,&l2))
		{
			if(tri_dbg_mode)
				puts("TriHeuG: fool (-+)");
			tri_set_rule(csf,l1,l2);
			FreeAtomic(m2l);
			FreeAtomic(l1);
			FreeAtomic(l2);
			return 1;
		}
		if(tri_dbg_mode)
			puts("TriHeu: pair (-+) fail");
	}
	
	if(tri_dbg_mode)
	{
		DumpList(m2l);
		printf("vars: ");WriteTerm(lu);puts("");
	}
	
	return 0;
}

List tri_csf(List *ml)
{
	List csf, csf1, l, lu, m2l;
	
	m2l = *ml;
	csf=CopyTerm(CompoundArg2(ListFirst(m2l)));
	
	for(l=ListTail(m2l);l;l=ListTail(l))
	{
		List l1,l2;
		for(l1=csf;l1;l1=ListTail(l1))
		{
			for(l2=CompoundArg2(ListFirst(l));l2;l2=ListTail(l2))
				if(CompoundArg1(ListFirst(l1))==CompoundArg1(ListFirst(l2)))
					break;
			if(l2)
			{
				int pt,pc;
				pt=(int)IntegerValue(CompoundArg2(ListFirst(l2)));
				pc=(int)IntegerValue(CompoundArg2(ListFirst(l1)));
				if(pt<pc)
					SetCompoundArg(ListFirst(l1),2,NewInteger(pt));
			}
			else
				SetCompoundArg(ListFirst(l1),2,NewInteger(0));
		}
	}
	
xyz1:
	for(l=csf;l;l=ListTail(l))
		if(CompoundArg2(ListFirst(l))==NewInteger(0))
		{
			csf=CutFromList(csf,l);
			goto xyz1;
		}
		
	
	for(l=m2l;l;l=ListTail(l))
	{
		List l1,l2,l3;
		l1=ConsumeCompoundArg(ListFirst(l),2);
		
		for(l2=csf;l2;l2=ListTail(l2))
			for(l3=l1;l3;l3=ListTail(l3))
				if(CompoundArg1(ListFirst(l2))==CompoundArg1(ListFirst(l3)))
				{
					int pt,pc;
					pt=(int)IntegerValue(CompoundArg2(ListFirst(l3)));
					pc=(int)IntegerValue(CompoundArg2(ListFirst(l2)));
					SetCompoundArg(ListFirst(l3),2,NewInteger(pt-pc));
					break;
				}
	xyz2:
		for(l2=l1;l2;l2=ListTail(l2))
			if(CompoundArg2(ListFirst(l2))==NewInteger(0))
			{
				l1=CutFromList(l1,l2);
				goto xyz2;
			}
			
		SetCompoundArg(ListFirst(l),2,l1);
	}
	
	lu=list_vars(m2l,0);
	
	for(l=lu;l;l=ListTail(l))
		m2l=tri_cvt(m2l,
				CompoundArg1(ListNth(tri_si_co_list,(int)IntegerValue(ListFirst(l)))),
				CompoundArg2(ListNth(tri_si_co_list,(int)IntegerValue(ListFirst(l)))));
	
	csf1=CopyTerm(CompoundArg2(ListFirst(m2l)));
	
	for(l=ListTail(m2l);l;l=ListTail(l))
	{
		List l1,l2;
		for(l1=csf1;l1;l1=ListTail(l1))
		{
			for(l2=CompoundArg2(ListFirst(l));l2;l2=ListTail(l2))
				if(CompoundArg1(ListFirst(l1))==CompoundArg1(ListFirst(l2)))
					break;
			if(l2)
			{
				int pt,pc;
				pt=(int)IntegerValue(CompoundArg2(ListFirst(l2)));
				pc=(int)IntegerValue(CompoundArg2(ListFirst(l1)));
				if(pt<pc)
					SetCompoundArg(ListFirst(l1),2,NewInteger(pt));
			}
			else
				SetCompoundArg(ListFirst(l1),2,NewInteger(0));
		}
	}
	
xyz11:
	for(l=csf1;l;l=ListTail(l))
		if(CompoundArg2(ListFirst(l))==NewInteger(0))
		{
			csf1=CutFromList(csf1,l);
			goto xyz11;
		}
		
	
	for(l=m2l;l;l=ListTail(l))
	{
		List l1,l2,l3;
		l1=ConsumeCompoundArg(ListFirst(l),2);
		
		for(l2=csf1;l2;l2=ListTail(l2))
			for(l3=l1;l3;l3=ListTail(l3))
				if(CompoundArg1(ListFirst(l2))==CompoundArg1(ListFirst(l3)))
				{
					int pt,pc;
					pt=(int)IntegerValue(CompoundArg2(ListFirst(l3)));
					pc=(int)IntegerValue(CompoundArg2(ListFirst(l2)));
					SetCompoundArg(ListFirst(l3),2,NewInteger(pt-pc));
					break;
				}
	xyz12:
		for(l2=l1;l2;l2=ListTail(l2))
			if(CompoundArg2(ListFirst(l2))==NewInteger(0))
			{
				l1=CutFromList(l1,l2);
				goto xyz12;
			}
			
		SetCompoundArg(ListFirst(l),2,l1);
	}
	
	for(l=lu;l;l=ListTail(l))
		m2l=tri_cvt(m2l,
				CompoundArg2(ListNth(tri_si_co_list,(int)IntegerValue(ListFirst(l)))),
				CompoundArg1(ListNth(tri_si_co_list,(int)IntegerValue(ListFirst(l)))));
	
	csf=ConcatList(csf,csf1);
	csf=SortedList(csf,pcmp);
	
	*ml=m2l;
	return csf;
	
}

static int isqrt(int x)
{
	int i;
	for(i=1;i*i<x;i++);
	if(i*i==x)
		return i;
	return 0;
}

static int red_1(int var_no, List ml, List *l_l, List *l_r)
{
	List l1,l2,lsq,csf,mlsv;

	l1=tri_find_sub(ml);
	
	if(l1)
	{
		if(tri_dbg_mode)
			puts("TriHeu1: exact sub found");
		*l_l=CopyTerm(ml);
		*l_r=l1;
		return 1;
	}
	
	if(ListLength(ml)==1)
	{
		List l,l1;
		int co;
		int cnum;
		
		if(tri_dbg_mode)
			puts("TriHeu1: monomial");
		
		cnum=(int)IntegerValue(CompoundArg1(ListFirst(ml)));
		l1=CopyTerm(ml);
		ml=CopyTerm(ml);
		l=ConsumeCompoundArg(ListFirst(l1),2);
		co=mk_s2(&l);
		SetCompoundArg(ListFirst(l1),2,l);
		SetCompoundArg(ListFirst(ml),1,NewInteger(cnum));
		SetCompoundArg(ListFirst(l1),1,NewInteger(cnum/co));
		*l_l=ml;
		*l_r=l1;
		
		return 1;
	}
	
	mlsv=ml;
	ml=CopyTerm(ml);
	csf=tri_csf(&ml);
	
	if(ListLength(ml)==1 && CompoundArg2(ListFirst(ml))==0)
	{
		int fa;
		fa=(int)IntegerValue(CompoundArg1(ListFirst(ml)));
		FreeAtomic(ml);
		ml=MakeCompound(A_MTERM,3);
		SetCompoundArg(ml,1,NewInteger(fa));
		SetCompoundArg(ml,2,csf);
		ml=MakeList1(ml);
		if(tri_dbg_mode)
			puts("TriHeu1: expr is const after csf extraction");
		*l_l=CopyTerm(ml);
		*l_r=CopyTerm(ml);
		FreeAtomic(ml);
		ml=ConsumeCompoundArg(ListFirst(*l_r),2);
		fa/=mk_s2(&ml);
		SetCompoundArg(ListFirst(*l_r),2,ml);
		SetCompoundArg(ListFirst(*l_r),1,NewInteger(fa));
		return 1;
	}

	FreeAtomic(ml);
	FreeAtomic(csf);
	ml=mlsv;

	lsq=0;
	if(ListLength(ml)==3)
	{
		Term m1,m2,m3;
		int i1,i2,i3,q1,q2,cf;
		
		m1=ListFirst(ml);
		m2=ListNth(ml,2);
		m3=ListNth(ml,3);
		
		if(CompoundArg2(m1)!=0 || ListLength(CompoundArg2(m2))!=1 ||
				ListLength(CompoundArg2(m3))!=1)
			goto xyz1;
		if(CompoundArg1(ListFirst(CompoundArg2(m2)))!=
				CompoundArg1(ListFirst(CompoundArg2(m3))))
			goto xyz1;
		
		if(CompoundArg2(ListFirst(CompoundArg2(m2)))!=NewInteger(2))
			goto xyz1;
		if(CompoundArg2(ListFirst(CompoundArg2(m3)))!=NewInteger(4))
			goto xyz1;
		

		i1=(int)IntegerValue(CompoundArg1(m1));
		i2=(int)IntegerValue(CompoundArg1(m2));
		i3=(int)IntegerValue(CompoundArg1(m3));
		cf=(int)gcf(i1,i2);
		cf=(int)gcf(cf,i3);
		if(i1<0)
			cf=-cf;
		i1/=cf;
		i2/=cf;
		i3/=cf;

		if(i3<0)
			goto xyz1;
		q1=isqrt(i1);
		q2=isqrt(i3);

		if(q1==0 || q2==0 || (i2!=2*q1*q2 && i2!=-2*q1*q2))
			goto xyz1;
		
		m3=CompoundArg1(ListFirst(CompoundArg2(m3)));
		m1=MakeCompound(A_MTERM,3);
		SetCompoundArg(m1,1,NewInteger(q1));
		m2=MakeCompound(A_MTERM,3);
		SetCompoundArg(m2,1,NewInteger(i2>0?q2:-q2));
		SetCompoundArg(m2,2,MakeList1(MakeCompound2(OPR_POW,m3,NewInteger(2))));
		l1=MakeList2(m1,m2);

		l2=tri_find_sub(l1);
		if(l2==0)
		{
			if(opUndefAngleComb==1)
			{
				tri_set_sub(l1);
				l2=tri_find_sub(l1);
				if(l2==0)
				{
					puts("Internal error (tri02)");
					return 0;
				}
			}
			else
			{
				Atom np;
				Term m2;
				if(tri_dbg_mode)
					puts("TriHeu1: aux sub for sqrt meth");
				np=tri_dummy_prm(CopyTerm(l1));
				m2=MakeCompound(A_MTERM,3);
				SetCompoundArg(m2,1,NewInteger(1));
				SetCompoundArg(m2,2,
						MakeList1(MakeCompound2(OPR_POW,np,NewInteger(1))));
				l2=MakeList1(m2);
				tri_set_rule(0,l1,l2);
				/*tri_set_rule(0,l1,l1);
				l2=CopyTerm(l1);*/
			}
		}
		
		m1=MakeCompound(A_MTERM,3);
		SetCompoundArg(m1,1,NewInteger(cf));
		SetCompoundArg(m1,2,MakeList2(l1,CopyTerm(l1)));
		*l_l=MakeList1(m1);
		
		m1=MakeCompound(A_MTERM,3);
		SetCompoundArg(m1,1,NewInteger(cf));
		SetCompoundArg(m1,2,MakeList2(l2,CopyTerm(l2)));
		*l_r=MakeList1(m1);
		
		if(tri_dbg_mode)
			puts("TriHeu1: sub**2 found");
		
		return 1;
	}
xyz1:
	
	if(!opUndefAngleComb)
	{
		List mlsv=CopyTerm(ml);
		if(tri_dbg_mode)
			puts("TriHeu1: leaving this as it is... (1)");
		ml=CopyTerm(ml);
		tri_opti_1(var_no,&ml);
		if(!EqualTerms(ml,mlsv))
		{
			
			if(tri_dbg_mode)
			puts("TriHeu1: opti_1 has improved the expr");
			*l_l=mlsv;
			*l_r=ml;
			return 1;
		}
		FreeAtomic(mlsv);
			
		l1=ml;
		if(ListLength(l1)>2 || opNoDummies)
			l2=CopyTerm(ml);
		else
		{
			Atom np;
			Term m2;
			int cf;
			List l3;
			l3=CopyTerm(l1);
			cf=(int)gcf(IntegerValue(CompoundArg1(ListFirst(l3))),
						IntegerValue(CompoundArg1(ListFirst(ListTail(l3)))));
			if(IntegerValue(CompoundArg1(ListFirst(l3)))<0)
				cf=-cf;
			SetCompoundArg(ListFirst(l3),1,
					NewInteger(IntegerValue(CompoundArg1(ListFirst(l3)))/cf));
			SetCompoundArg(ListFirst(ListTail(l3)),1,
					NewInteger(IntegerValue(
					CompoundArg1(ListFirst(ListTail(l3))))/cf));
			
			np=tri_dummy_prm(l3);
			m2=MakeCompound(A_MTERM,3);
			SetCompoundArg(m2,1,NewInteger(cf));
			SetCompoundArg(m2,2,
					MakeList1(MakeCompound2(OPR_POW,np,NewInteger(1))));
			l2=MakeList1(m2);
		}

		/*{
			List k1,k2;
			k1=CopyTerm(l1);
			k2=CopyTerm(l2);
			tri_set_rule(0,k1,k2);
			FreeAtomic(k1);
			FreeAtomic(k2);
		}*/

		if(tri_dbg_mode)
			puts("TriHeu1: well, ...");
		
		*l_l=l1;
		*l_r=l2;
		return 1;
	}
	
	
	if(tri_dbg_mode)
		puts("TriHeu1: leaving this as it is...(2)");
	tri_opti_1(var_no,&ml);
	l1=CopyTerm(ml);
	tri_set_sub(l1);
	l2=tri_find_sub(l1);
	if(tri_dbg_mode)
			puts("TriHeu1: well, ...");
	*l_l=l1;
	*l_r=l2;
	return 1;
}

int tri_heu_2(Term v2d, List ml, List *l_l, List *l_r);

static int red_2(int var_no, List ml, List *l_l, List *l_r)
{
	
	return tri_heu_2(ListNth(tri_sc_2_list,var_no),ml,l_l,l_r);
		
}
 

static List list_vars(List m2l, int do2)
{
	List lu, lu1, l;
	Term lu2;
	
	lu=NewList();
	for(l=m2l;l;l=ListTail(l))
	{
		List l1;
		for(l1=CompoundArg2(ListFirst(l));l1;l1=ListTail(l1))
			if(!ListMember(lu,CompoundArg1(ListFirst(l1))))
				lu=AppendLast(lu,CompoundArg1(ListFirst(l1)));
	}
	
	lu1=NewList();
	
	for(l=lu;l;l=ListTail(l))
	{
		List l1;
		int i=1;
		
		for(l1=tri_si_co_list;l1;l1=ListTail(l1))
		{
			if(ListFirst(l)==CompoundArg1(ListFirst(l1)) ||
					ListFirst(l)==CompoundArg2(ListFirst(l1)))
			{
				if(!ListMember(lu1,NewInteger(i)))
					lu1=AppendLast(lu1,NewInteger(i));
			}
			i++;
		}
				
	}

	if(!do2)
	{
		FreeAtomic(lu);
		return lu1;
	}
	
	lu2=0;
	for(l=lu1;l;l=ListTail(l))
	{
		Term sc;
		sc=ListNth(tri_si_co_list,(int)IntegerValue(ListFirst(l)));
		if(CompoundArgN(sc,4)==0)
			continue;
		if(lu2==0)
		{
			lu2=CompoundArgN(sc,4);
			continue;
		}
		if(CompoundArgN(sc,4)!=lu2)
		{
			FreeAtomic(m2l);
			return 0;
		}
	}
	
	if(lu2)
	{
		Term lu22;
		lu22=ListNth(tri_sc_2_list,(int)IntegerValue(lu2));
		if(!ListMember(lu1,CompoundArg1(lu22)) || 
				!ListMember(lu1,CompoundArg2(lu22)))
			lu2=0;
		if(CompoundArgN(lu22,3)!=NewInteger(4))
			lu2=0;
	}
	
/*	if(lu2)
	{
		Term lu22;
		Atom c1,s1,c2,s2;
		
		lu22=ListNth(tri_sc_2_list,IntegerValue(lu2));
		s1=CompoundArg1(ListNth(tri_si_co_list,IntegerValue(CompoundArg1(lu22))));
		c1=CompoundArg2(ListNth(tri_si_co_list,IntegerValue(CompoundArg1(lu22))));
		s2=CompoundArg1(ListNth(tri_si_co_list,IntegerValue(CompoundArg2(lu22))));
		c2=CompoundArg2(ListNth(tri_si_co_list,IntegerValue(CompoundArg2(lu22))));
		
		for(l=m2l;l;l=ListTail(l))
		{
			int k1=0,k2=0;
			for(l1=CompoundArg2(ListFirst(l));l1;l1=ListTail(l1))
			{
				Atom pp;
				pp=CompoundArg1(ListFirst(l1));
				if(pp==s1 || pp==c1)
					k1++;
				if(pp==s2 || pp==c2)
					k2++;
			}
			if(k1 && k2)
				break;
		}
		
		if(is_empty_list(l))
			lu2=0;
	}
*/		
	if(lu2)
	{
		Term lu22;
		
		lu22=ListNth(tri_sc_2_list,(int)IntegerValue(lu2));
		
			for(l=lu1;l;l=ListTail(l))
				if(ListFirst(l)==CompoundArg1(lu22))
				{
					lu1=CutFromList(lu1,l);
					break;
				}
			for(l=lu1;l;l=ListTail(l))
				if(ListFirst(l)==CompoundArg2(lu22))
				{
					lu1=CutFromList(lu1,l);
					break;
				}
		
	}
	
	FreeAtomic(lu);
	
	if(lu2)
		lu1=AppendFirst(lu1,NewInteger(-IntegerValue(lu2)));
	
	return lu1;
	
}


static int mk_s2(List *p)
{
	List l1,l2,l3;
	int ret=1,i;
xyz:
	
	for(l1=*p;l1;l1=ListTail(l1))
		for(l2=ListTail(l1);l2;l2=ListTail(l2))
		{
			Atom a1,a2;
			if(CompoundArg2(ListFirst(l1))!=CompoundArg2(ListFirst(l2)))
				continue;
			a1=CompoundArg1(ListFirst(l1));
			a2=CompoundArg1(ListFirst(l2));
			for(l3=tri_si_co_list;l3;l3=ListTail(l3))
				if(((CompoundArg1(ListFirst(l3))==a1 
						&& CompoundArg2(ListFirst(l3))==a2) ||
				    (CompoundArg1(ListFirst(l3))==a2 
						&& CompoundArg2(ListFirst(l3))==a1)) &&
						CompoundArgN(ListFirst(l3),3))
				{
					int po;
					po=(int)IntegerValue(CompoundArg2(ListFirst(l1)));
					if(ListLength(*p)>2)
					{
						*p=CutFromList(*p,l1);
						*p=CutFromList(*p,l2);
					}
					else
					{
						FreeAtomic(*p);
						*p=NewList();
					}
					*p=AppendFirst(*p,MakeCompound2(OPR_POW,
							CompoundArgN(ListFirst(l3),3),NewInteger(po)));
					for(i=0;i<po;i++)
						ret*=2;
					goto xyz;
				}
		}
	if(ret>1)
		*p=SortedList(*p,pcmp);
	return ret;
}

int tri_mcmp(Term m1, Term m2)
	{
	int rs;
	m1=CompoundArg2(m1);
	m2=CompoundArg2(m2);
	rs=ListLength(m1)-ListLength(m2);

	if(rs)
		return rs;

	while(!is_empty_list(m1))
		{
		int rs;
		rs=strcmp(AtomValue(CompoundArg1(ListFirst(m1))),
			AtomValue(CompoundArg1(ListFirst(m2))));
		if(rs)
			return rs;
		rs=(int)IntegerValue(CompoundArg2(ListFirst(m1))) -
			(int)IntegerValue(CompoundArg2(ListFirst(m2)));
		if(rs)
			return rs;
		m1=ListTail(m1);
		m2=ListTail(m2);
		}
	puts("eq mt");
	return 0;
	}



Term tri_m2_l_to_e(List);
static int toe_mno=0;

static Term m2_p_to_e(Term pt)
{
	if(is_list(pt))
		return tri_m2_l_to_e(pt);
	
	if(is_list(CompoundArg1(pt)) && CompoundArg2(pt)==NewInteger(1))
		return tri_m2_l_to_e(CompoundArg1(pt));
	
	if(is_list(CompoundArg1(pt)))
		return MakeCompound2(OPR_POW,tri_m2_l_to_e(CompoundArg1(pt)),CompoundArg2(pt));
	
	if(CompoundArg2(pt)==NewInteger(1))
		return CompoundArg1(pt);
	
	return CopyTerm(pt);
	
}

static Term m2_m_to_e(Term m2)
{
	int n;
	List l1;
	Term ret;
	
	n=(int)IntegerValue(CompoundArg1(m2));
	if(n<0)
		n=-n;
	
	if(n==1 && CompoundArg2(m2)==0)
		return NewInteger(1);
	
	if(n==1)
		ret=m2_p_to_e(ListFirst(CompoundArg2(m2)));
	else
		ret=NewInteger(n);
	
	if(n==1)
		l1=ListTail(CompoundArg2(m2));
	else
		l1=CompoundArg2(m2);
	
	for(;l1;l1=ListTail(l1))
		ret=MakeCompound2(OPR_MLT,ret,m2_p_to_e(ListFirst(l1)));
	
	return ret;
	
}

Term tri_m2_l_to_e(List ml)
{
	Term ret;
	List l;
	
	ret=m2_m_to_e(ListFirst(ml));
	if(IntegerValue(CompoundArg1(ListFirst(ml)))<0)
		ret=MakeCompound1(OPR_MINUS,ret);
	
	toe_mno+=ListLength(ml)-1;
	
	if(ListLength(ml)==1)
		return ret;
	
	for(l=ListTail(ml);l;l=ListTail(l))
		if(IntegerValue(CompoundArg1(ListFirst(l)))<0)
			ret=MakeCompound2(OPR_MINUS,ret,m2_m_to_e(ListFirst(l)));
		else
			ret=MakeCompound2(OPR_PLUS,ret,m2_m_to_e(ListFirst(l)));
	
	return ret;
}

void tri_red_numf(List l1, List l2, List l3, List l4)
{
	int cf;
	List l;
	
	cf=(int)IntegerValue(CompoundArg1(ListFirst(l1)));
	if(cf<0)
		cf=-cf;
	for(l=l1;l;l=ListTail(l))
		cf=(int)gcf(cf,IntegerValue(CompoundArg1(ListFirst(l))));
	for(l=l2;l;l=ListTail(l))
		cf=(int)gcf(cf,IntegerValue(CompoundArg1(ListFirst(l))));
	if(l3)
		for(l=l3;l;l=ListTail(l))
			cf=(int)gcf(cf,IntegerValue(CompoundArg1(ListFirst(l))));
	if(l4)
		for(l=l4;l;l=ListTail(l))
			cf=(int)gcf(cf,IntegerValue(CompoundArg1(ListFirst(l))));
	for(l=l1;l;l=ListTail(l))
		SetCompoundArg(ListFirst(l),1,
				NewInteger(IntegerValue(CompoundArg1(ListFirst(l)))/cf));
	for(l=l2;l;l=ListTail(l))
		SetCompoundArg(ListFirst(l),1,
				NewInteger(IntegerValue(CompoundArg1(ListFirst(l)))/cf));
	if(l3)
		for(l=l3;l;l=ListTail(l))
			SetCompoundArg(ListFirst(l),1,
				NewInteger(IntegerValue(CompoundArg1(ListFirst(l)))/cf));
	if(l4)
		for(l=l4;l;l=ListTail(l))
			SetCompoundArg(ListFirst(l),1,
				NewInteger(IntegerValue(CompoundArg1(ListFirst(l)))/cf));
}
	
void tri_set_rule(List csf, List left, List right)
{
	Term t,e1,e2;
	
	List l1;
	
	tri_red_numf(left, right, 0, 0);
	
	if(csf)
	{
		Term t;
		int nu;
		l1=CopyTerm(csf);
		nu=mk_s2(&csf);
		l1=AppendLast(l1,left);
		t=MakeCompound(A_MTERM,3);
		SetCompoundArg(t,1,NewInteger(nu));
		SetCompoundArg(t,2,l1);
		left=AppendLast(NewList(),t);

		l1=csf;
		l1=AppendLast(l1,right);
		t=MakeCompound(A_MTERM,3);
		SetCompoundArg(t,1,NewInteger(1));
		SetCompoundArg(t,2,l1);
		right=AppendLast(NewList(),t);
	}
	
	
	e1=tri_m2_l_to_e(left);
	toe_mno=0;
	e2=tri_m2_l_to_e(right);
	t=MakeCompound2(OPR_EQSIGN,e1,e2);
	if(opTriHeu==3 || ((opTriHeu==2 || opTriHeu==4) && toe_mno>2))
	{
		printf("SetAngle(");
		WriteTerm(t);
		puts(").");
	}
	t=MakeCompound1(A_I,t);
	ProcSetAngle(t,0);
}

