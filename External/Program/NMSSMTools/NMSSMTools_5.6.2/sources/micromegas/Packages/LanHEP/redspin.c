#include "lanhep.h"
#include "math.h"

/*
#define DEB
*/

static void split_m2(Term m, Term *m1, Term *m2, int *anti, int spin)
	{
	int diff=0, fff1=1, fff2=1;
	List l;
	List co1=0, co2=0;
	
	int i1,i2,n1,n2,d1,d2;
	
	*m1=MakeCompound(A_MTERM,3);
	*m2=MakeCompound(A_MTERM,3);
	
	l=CompoundArg2(m);
	
	
	while(!is_empty_list(l))
		{
		Term t;
		int i1,i2;
		t=ListFirst(l);
		i2=(int)IntegerValue(CompoundArg2(t));
		i1=i2/2;
		i2-=i1;
		if(i1!=i2)
			diff=1;
		if(i1)
			co1=AppendLast(co1,MakeCompound2(OPR_POW,
				CompoundArg1(t),NewInteger(i1)));
		if(i2)
			co2=AppendLast(co2,MakeCompound2(OPR_POW,
				CompoundArg1(t),NewInteger(i2)));
		l=ListTail(l);
		}
		
	
	
	l=CompoundArg1(m);
	i1=(int)IntegerValue(CompoundArg1(l));
	i2=(int)IntegerValue(CompoundArg2(l));
	
	if(i1<0)
		{
			List l1;
			
			for(l1=co1;l1;l1=ListTail(l1))
				if(CompoundArg1(ListFirst(l1))==A_I)
					break;
			if(l1)
			{
				co1=CutFromList(co1,l1);
				fff1=-1;
			}
			else
				co1=AppendFirst(co1,MakeCompound2(OPR_POW,A_I,NewInteger(1)));
			
			for(l1=co2;l1;l1=ListTail(l1))
				if(CompoundArg1(ListFirst(l1))==A_I)
					break;
			if(l1)
			{
				co2=CutFromList(co2,l1);
				fff2=-1;
			}
			else
				co2=AppendFirst(co2,MakeCompound2(OPR_POW,A_I,NewInteger(1)));
			i1=-i1;
		}
	
	
	n1=(int)floor(sqrt(i1)+0.5);
	d1=(int)floor(sqrt(i2)+0.5);
	
	if(n1*n1!=i1 || d1*d1!=i2)
		{
		diff=1;
		n1=i1;
		d1=i2;
		n2=1;
		d2=1;
		}
	else
		{
		n2=n1;
		d2=d1;
		}
		
	n1*=fff1;
	n2*=fff2;
		
	if(diff)
		*anti=1;
		
	if(spin<2)
		{
			List l1;
			co1=AppendLast(co1,
					MakeCompound2(OPR_POW,NewAtom("Maux",0),NewInteger(1)));
			if(spin==1)
			{
			for(l1=co1;l1;l1=ListTail(l1))
				if(CompoundArg1(ListFirst(l1))==A_I)
					break;
			if(l1)
			{
				co1=CutFromList(co1,l1);
				n1=-n1;
			}
			else
				co1=AppendFirst(co1,MakeCompound2(OPR_POW,A_I,NewInteger(1)));
			}
			
			co2=AppendLast(co2,
					MakeCompound2(OPR_POW,NewAtom("Maux",0),NewInteger(1)));
			if(spin==1)
			{
			for(l1=co2;l1;l1=ListTail(l1))
				if(CompoundArg1(ListFirst(l1))==A_I)
					break;
			if(l1)
			{
				co2=CutFromList(co2,l1);
				n2=-n2;
			}
			else
				co2=AppendFirst(co2,MakeCompound2(OPR_POW,A_I,NewInteger(1)));
			}
		}
		
	if(*anti==0)
		{
		co1=AppendLast(co1,MakeCompound2(OPR_POW,A_SQRT2,NewInteger(1)));	
		co2=AppendLast(co2,MakeCompound2(OPR_POW,A_SQRT2,NewInteger(1)));	
		}
		
			
	SetCompoundArg(*m1,2,co1);
	SetCompoundArg(*m2,2,co2);
	
	SetCompoundArg(*m1,1,MakeCompound2(OPR_DIV,NewInteger(n1),
			NewInteger(d1)));
	SetCompoundArg(*m2,1,MakeCompound2(OPR_DIV,NewInteger(n2),
			NewInteger(d2)));
	
	}


static int check_charge(List pl1, List prt1, List osp1, Term cs1,
						List pl2, List prt2, List osp2, Term cs2, List fi)
	{
	List l,l1;
	l=prt1;
	l1=prt2;
	if(osp1 || osp2) return 1;
	if(cs1 && CompoundArg1(cs1)!=CompoundArg1(cs2))
		return 1;
	while(!is_empty_list(l))
		{
		if(is_empty_list(l1))
			return 1;
		if(CompoundArg1(ListFirst(l))!=CompoundArg1(ListFirst(l1)))
			return 1;
		l=ListTail(l);
		l1=ListTail(l1);
		}
	return 0;
	}
	
static List spin_reduce_2(Term a2, List pl, List pil, List csp, List osp)
	{
	List pl1, pl2;
	List prt, prt1, prt2;
	List osp1, osp2, csp1, csp2;
	List l1,l2;
	List i1,i2;
	List co1,co2;
	Label collab;
	Term m2,m22,a22,m;	
	int anti;

/*	
	return AppendLast(NewList(),a2);
*/	
	
	prt1=prt2=pl1=pl2=osp1=osp2=csp1=csp2=NewList();
	
		
	
			
	/********* sep prtcls *************/
		
		
	prt=CompoundArg1(a2);
	
	l1=pl;
	l2=pil;
	
	while(!is_empty_list(l1))
		{
		if(CompoundArg1(ListFirst(l2))==NewInteger(-1))
			{
			collab=CompoundArg2(ListFirst(l2));
			pl1=AppendLast(pl1,ListFirst(l1));
			ChangeList(l1,0);
			pl=CutFromList(pl,l1);
			pil=CutFromList(pil,l2);
			goto cnt1;
			}
		l1=ListTail(l1);
		l2=ListTail(l2);
		}
	printf("Internal error in vertex ");
	WriteVertex(CompoundArg1(a2));
	printf(" (rs01).\n");
	return 0;

cnt1:
	
	l1=pl;
	l2=pil;
	while(!is_empty_list(l1))
		{
		if(CompoundArg2(ListFirst(l2))==collab)
			{
			pl1=AppendLast(pl1,ListFirst(l1));
			ChangeList(l1,0);
			pl=CutFromList(pl,l1);
			pil=CutFromList(pil,l2);
			goto cnt2;
			}
		l1=ListTail(l1);
		l2=ListTail(l2);
		}
	l1=csp;
	while(!is_empty_list(l1))
		{
		if(ListFirst(CompoundArg2(ListFirst(l1)))==collab)
			{
			csp1=AppendLast(csp1,ListFirst(l1));
			collab=ListNth(CompoundArg2(ListFirst(l1)),2);
			ChangeList(l1,0);
			csp=CutFromList(csp,l1);
			goto cnt1;
			}
		l1=ListTail(l1);
		}
	printf("Internal error in vertex ");
	WriteVertex(CompoundArg1(a2));
	printf(" (rs02).\n");
	return 0;
	
cnt2:	
	
	l1=pl;
	l2=pil;
	
	while(!is_empty_list(l1))
		{
		if(CompoundArg1(ListFirst(l2))==NewInteger(-1))
			{
			collab=CompoundArg2(ListFirst(l2));
			pl2=AppendLast(pl2,ListFirst(l1));
			ChangeList(l1,0);
			pl=CutFromList(pl,l1);
			pil=CutFromList(pil,l2);
			goto cnt22;
			}
		l1=ListTail(l1);
		l2=ListTail(l2);
		}
	printf("Internal error in vertex ");
	WriteVertex(CompoundArg1(a2));
	printf(" (rs05).\n");
	return 0;

cnt22:
	
	l1=pl;
	l2=pil;
	while(!is_empty_list(l1))
		{
		if(CompoundArg2(ListFirst(l2))==collab)
			{
			pl2=AppendLast(pl2,ListFirst(l1));
			ChangeList(l1,0);
			pl=CutFromList(pl,l1);
			pil=CutFromList(pil,l2);
			goto cnt3;
			}
		l1=ListTail(l1);
		l2=ListTail(l2);
		}
	l1=csp;
	while(!is_empty_list(l1))
		{
		if(ListFirst(CompoundArg2(ListFirst(l1)))==collab)
			{
			csp2=AppendLast(csp2,ListFirst(l1));
			collab=ListNth(CompoundArg2(ListFirst(l1)),2);
			ChangeList(l1,0);
			csp=CutFromList(csp,l1);
			goto cnt22;
			}
		l1=ListTail(l1);
		}
	printf("Internal error in vertex ");
	WriteVertex(CompoundArg1(a2));
	printf(" (rs03).\n");
	return 0;
	
cnt3:	
	
	
	l1=prt;
	while(!is_empty_list(l1))
		{
		Label lab;
		lab=CompoundArg2(ListFirst(l1));
		if(lab==CompoundArg1(ListFirst(pl1)) ||
		   lab==CompoundArg1(ListNth(pl1,2)) )
		   	prt1=AppendLast(prt1,ListFirst(l1));
		else
			prt2=AppendLast(prt2,ListFirst(l1));
		l1=ListTail(l1);
		}

#ifdef DEB	
	printf("prt1: "); WriteTerm(pl1); WriteTerm(prt1); puts("");
	printf("prt2: "); WriteTerm(pl2); WriteTerm(prt2); puts("");
#endif
	
	if(!is_empty_list(csp) || !is_empty_list(pl))
		{
		printf("Internal error in vertex ");
		WriteVertex(CompoundArg1(a2));
		printf(" (rs04).\n");
		return 0;
		}
	
#ifdef DEB
printf("g1: "); WriteTerm(csp1); puts("");	
printf("g2: "); WriteTerm(csp2); puts("");	
#endif	
	
	
	
	/***********  sep spec ***************/
	
	l1=osp;
	while(!is_empty_list(l1))
		{
		Term t;
		t=ListFirst(l1);
#ifdef DEB
		WriteTerm(t); 
#endif
		if(CompoundName(t)==A_MOMENT)
			{
			l2=CompoundArg1(t);
			if(l2==CompoundArg1(ListFirst(pl1)) ||
				l2==CompoundArg1(ListNth(pl1,2)))
				{
#ifdef DEB
				puts(" -> 1");
#endif
				osp1=AppendLast(osp1,t);
				goto cnt4;
				}
			if(l2==CompoundArg1(ListFirst(pl2)) ||
				l2==CompoundArg1(ListNth(pl2,2)))
				{
#ifdef DEB
				puts(" -> 2");
#endif
				osp2=AppendLast(osp2,t);
				goto cnt4;
				}
			printf("Internal error in vertex ");
			WriteVertex(CompoundArg1(a2));
			printf(" (rs06).\n");
			return 0;
			}
			
		printf("Internal error in vertex ");
		WriteVertex(CompoundArg1(a2));
		printf(" (rs07 -- "); WriteTerm(CompoundName(t));
		printf(":"); WriteTerm(CompoundArg1(t)); puts(").");
		return 0;
cnt4:	l1=ListTail(l1);
		}
		
	RemoveList(osp);
	
	
	
	/**********  free ind ***************/
	
	i1=i2=NewList();
	l1=pl1;
	while(!is_empty_list(l1))
		{
		l2=CompoundArg2(ListFirst(l1));
		while(!is_empty_list(l2))
			{
			i1=XorList(i1,ListFirst(l2));
			l2=ListTail(l2);
			}
		l1=ListTail(l1);
		}
	
	l1=osp1;
	while(!is_empty_list(l1))
		{
		l2=CompoundArg2(ListFirst(l1));
		while(!is_empty_list(l2))
			{
			i1=XorList(i1,ListFirst(l2));
			l2=ListTail(l2);
			}
		l1=ListTail(l1);
		}
		
	l1=csp1;
	while(!is_empty_list(l1))
		{
		l2=CompoundArg2(ListFirst(l1));
		while(!is_empty_list(l2))
			{
			i1=XorList(i1,ListFirst(l2));
			l2=ListTail(l2);
			}
		l1=ListTail(l1);
		}
		
	l1=pl2;
	while(!is_empty_list(l1))
		{
		l2=CompoundArg2(ListFirst(l1));
		while(!is_empty_list(l2))
			{
			i2=XorList(i2,ListFirst(l2));
			l2=ListTail(l2);
			}
		l1=ListTail(l1);
		}
	
	l1=osp2;
	while(!is_empty_list(l1))
		{
		l2=CompoundArg2(ListFirst(l1));
		while(!is_empty_list(l2))
			{
			i2=XorList(i2,ListFirst(l2));
			l2=ListTail(l2);
			}
		l1=ListTail(l1);
		}
		
	l1=csp2;
	while(!is_empty_list(l1))
		{
		l2=CompoundArg2(ListFirst(l1));
		while(!is_empty_list(l2))
			{
			i2=XorList(i2,ListFirst(l2));
			l2=ListTail(l2);
			}
		l1=ListTail(l1);
		}
	
	
	if(ListLength(i1)!=ListLength(i2) ||
		ListLength(i1)>2)
		{
		printf("Internal error in vertex ");
		WriteVertex(CompoundArg1(a2));
		printf(" (rs08).\n");
		return 0;
		}
//	FreeAtomic(i2);

#ifdef DEB
printf("im indices: ");	WriteTerm(i1); printf(" x "); WriteTerm(i2); puts("");
#endif


	
	
	/**********    splitting  ***********/

	anti=check_charge(pl1,prt1,osp1, csp1, pl2,prt2,osp2,csp2,i1);
	
	
	l1=ConsumeCompoundArg(a2,5);
	m=ListFirst(l1);
	RemoveList(l1);
	prt=ConsumeCompoundArg(a2,1);
	a22=MakeCompound(A_ALG2,5);
	co2=0;
	l1=ConsumeCompoundArg(m,3);
	RemoveList(l1);
	
	split_m2(m,&m2,&m22,&anti, ListLength(i1));
	FreeAtomic(m);
	
	l2=mk_im_field(0, ListLength(i1), anti, 
		       ConcatList(CopyTerm(prt1),CopyTerm(prt2)));
	
	pl1=ConcatList(pl1,osp1);
	pl1=ConcatList(pl1,csp1);
	co1=NewLabel();
	prt1=AppendLast(prt1,MakeCompound2(OPR_DIV,ListFirst(l2),co1));
	switch(ListLength(i1))
		{
		case 0: co2=MakeCompound2(OPR_SCALAR,co1,CopyTerm(i1)); break;
		case 1: co2=MakeCompound2(OPR_VECTOR,co1,CopyTerm(i1)); break;
		case 2: co2=MakeCompound2(OPR_TENSOR,co1,CopyTerm(i1)); break;
		}
	pl1=AppendLast(pl1,co2);
	SetCompoundArg(m2,3,pl1);
	SetCompoundArg(a2,1,prt1);
	SetCompoundArg(a2,5,AppendLast(NewList(),m2));
	
	pl2=ConcatList(pl2,osp2);
	pl2=ConcatList(pl2,csp2);
	co1=NewLabel();
	prt2=AppendLast(prt2,MakeCompound2(OPR_DIV,ListNth(l2,2),co1));
	RemoveList(l2);
	switch(ListLength(i1))
		{
		case 0: co2=MakeCompound2(OPR_SCALAR,co1,i1); break;
		case 1: co2=MakeCompound2(OPR_VECTOR,co1,i1); break;
		case 2: co2=MakeCompound2(OPR_TENSOR,co1,i1); break;
		}
	pl2=AppendLast(pl2,co2);
	SetCompoundArg(m22,3,pl2);
	SetCompoundArg(a22,1,prt2);
	SetCompoundArg(a22,5,AppendLast(NewList(),m22));
	
	
#ifdef DEB
	puts("\n"); WriteTerm(a2); puts("\n"); WriteTerm(a22); puts("");
#endif	

	sort_spin_prt(a2);
	sort_spin_prt(a22);
	
	if(anti)	
		return MakeList2(a2,a22);
	else
		{
		FreeAtomic(a22);
		return AppendLast(NewList(),a2);
		}
	
	
	
	}
	

	
	
static int lab2sp(Label lab, List pl)
	{
	while(!is_empty_list(pl))
		{
		if(CompoundArg2(ListFirst(pl))==lab)
			{
			Term prp;
			prp=GetAtomProperty(CompoundArg1(ListFirst(pl)),PROP_INDEX);
			if(prp==0) 
				return 0;
			prp=CompoundArg1(ListFirst(prp));
			if(is_compound(prp) && CompoundName(prp)==A_LORENTZ)
				{
				if(CompoundArg1(prp)==NewInteger(1))
					return 1;
				if(CompoundArg1(prp)==NewInteger(0))
					return -1;
				}
			return 0;
			}
		pl=ListTail(pl);
		}
	puts("Internal error (l2s).");
	return 0;
	}

List spinor_reduce(Term a2)
	{
	List m2, cp_l, op_l, cs_l, os_l, pi_l, l1;
	cp_l=pi_l=op_l=cs_l=os_l=NewList();

#ifdef DEB
	WriteTerm(a2); puts("");
#endif

	m2=ListFirst(CompoundArgN(a2,5));
	
	
	l1=CompoundArgN(m2,3);
	while(!is_empty_list(l1))
		{
		Term t1;
		t1=ListFirst(l1);
		if(CompoundName(t1)==OPR_TENSOR || CompoundName(t1)==OPR_VECTOR ||
		   CompoundName(t1)==OPR_SPINOR || CompoundName(t1)==OPR_SCALAR ||
		   CompoundName(t1)==OPR_SPINOR3)
		   	{
		   	int sp;
		   	sp=lab2sp(CompoundArg1(t1), CompoundArg1(a2));
		   	if(sp)
		   		{
		   		cp_l=AppendLast(cp_l,t1);
		   		pi_l=AppendLast(pi_l,MakeCompound2(OPR_DIV,NewInteger(sp),
		   			ListFirst(CompoundArg2(t1))));
		   		goto cnt1;
		   		}
		   	op_l=AppendLast(op_l,t1);
		   	goto cnt1;
		   	}
		if(CompoundName(t1)==OPR_SPECIAL && (CompoundArg1(t1)==A_GAMMA ||
				CompoundArg1(t1)==A_GAMMA5))
			{
			cs_l=AppendLast(cs_l,t1);
			goto cnt1;
			}
		os_l=AppendLast(os_l,t1);
	cnt1:
		l1=ListTail(l1);
		}
		
	if(ListLength(cp_l)==0)
		{
		if(ListLength(cs_l))
			{
			printf("Error in vertex "); WriteVertex(CompoundArg1(a2));
			printf(": illegal spinor structure.\n");
			}
		RemoveList(op_l);
		RemoveList(os_l);
		return AppendFirst(NewList(),a2);
		}

#ifdef DEB	
	printf("sp: ");WriteTerm(cp_l);printf("\n\twith ");WriteTerm(pi_l);puts("");
	printf("ss: ");WriteTerm(cs_l);puts("");
	printf("os: ");WriteTerm(os_l);puts("");	
#endif
		
	if(ListLength(cp_l)==2)
		{
		
		RemoveList(op_l);
		RemoveList(os_l);
		RemoveList(cp_l);
		RemoveList(cs_l);
		FreeAtomic(pi_l);
		return AppendFirst(NewList(),a2);
		}
	
		
	
	if(ListLength(cp_l)!=4 || op_l)
		{
		printf("Error in vertex "); WriteVertex(CompoundArg1(a2));
		printf(": illegal spinor structure.\n");
		RemoveList(op_l);
		RemoveList(os_l);
		RemoveList(cp_l);
		RemoveList(cs_l);
		FreeAtomic(pi_l);
		return AppendFirst(NewList(),a2);
		}
		
	
	
	return spin_reduce_2(a2, cp_l, pi_l, cs_l, os_l);
	
	}
	

int need_spin_rdc(Term a2)
{
	List l1;
	int sno=0;
	
	for(l1=CompoundArg1(a2);l1;l1=ListTail(l1))
		if(is_compound(CompoundArg2(ListFirst(l1))) &&
				(CompoundName(CompoundArg2(ListFirst(l1)))==OPR_SPINOR
				||CompoundName(CompoundArg2(ListFirst(l1)))==OPR_SPINOR3) )
			sno++;
	
	if(sno==0 || sno==2)
		return 0;
	
	if(sno==4)
		return 1;
	
	printf("Error in vertex ");
	WriteVertex(CompoundArg1(a2));
	printf(": illegal spinor structure.\n");
	return 0;
}

	

