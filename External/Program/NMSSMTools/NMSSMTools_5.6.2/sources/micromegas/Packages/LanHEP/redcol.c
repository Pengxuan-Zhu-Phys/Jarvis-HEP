#include "lanhep.h"
#include "math.h"

/*
#define DEB
*/

extern int is_color_eps(Atom a);
extern int NoColors;
extern int WriteColors;
int opNeutrC8=0;

void  sort_spin_prt(Term a2)
	{
	List prt, l1,l2;
	prt=CompoundArg1(a2);


	l1=prt;
	while(!is_empty_list(l1))
		{
		l2=ListTail(l1);
		while(!is_empty_list(l2))
			{
			Atom p1,p2;
			Term i1,i2;
			p1=ListFirst(l1);
			p2=ListFirst(l2);
			/*WriteTerm(p1); WriteTerm(p2); puts("prtcls");*/
			i1=GetAtomProperty(CompoundArg1(p1),PROP_INDEX);
			i2=GetAtomProperty(CompoundArg1(p2),PROP_INDEX);
			/*WriteTerm(i1); WriteTerm(i2); puts("props");*/
			if(i1 && i2)
				{
				i1=CompoundArg1(ListFirst(i1));
				i2=CompoundArg1(ListFirst(i2));
				if(is_compound(i1) && is_compound(i2) &&
					CompoundName(i1)==A_LORENTZ && CompoundName(i2)==A_LORENTZ &&
					CompoundArg1(i1)==NewInteger(1) &&
					CompoundArg1(i2)==NewInteger(0))
					{
					Term m2;
					ChangeList(l1,p2);
					ChangeList(l2,p1);
					m2=CompoundArg1(ListFirst(CompoundArgN(a2,5)));
					SetCompoundArg(m2,1,NewInteger(-IntegerValue(
						CompoundArg1(m2))));
					return;
					}
				}
			l2=ListTail(l2);
			}
		l1=ListTail(l1);
		}
	}


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
	if(osp1 || osp2) 
	{
	if(ListLength(osp1)!=ListLength(osp2))
		return 1;
	l=osp1;l1=osp2;
	for(;l;l=ListTail(l),l1=ListTail(l1))
		if(CompoundArg1(ListFirst(l))!=CompoundArg1(ListFirst(l1)))
			return 1;
	}
	l=prt1;
	l1=prt2;
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

	
int check_splitting_c2=0;
static int warnc8=0;
	
static List color_reduce_2(Term a2, List pl, List pil, List csp, List osp)
	{
	List pl1, pl2;
	List prt, prt1, prt2;
	List osp1, osp2;
	List l1,l2;
	List i1,i2;
	List co1,co2;
	Label collab;
	Term m2,m22,a22,m, csp1, csp2;
	int anti,ctype;

#ifdef DEB

	printf("color prt: "); WriteTerm(pl);
	printf(" with ");      WriteTerm(pil);
	printf("\ncol sp:    ");  WriteTerm(csp);
	printf("\noth sp:    ");  WriteTerm(osp); puts("\n");

#endif

/*
	return AppendLast(NewList(),a2);
*/

	/*******  prec csp   ***********/


	csp1=ListFirst(csp);
	csp2=ListNth(csp,2);

	co1=co2=collab=0;

	l1=CompoundArg2(csp1);
	l2=CompoundArg2(csp2);
	for(i1=1;i1<=3;i1++)
	for(i2=1;i2<=3;i2++)
		{
		if(ListNth(l1,(int)i1)==ListNth(l2,(int)i2))
			{
			if(co1)
				{
				printf("Internal error in vertex ");
				WriteVertex(CompoundArg1(a2));
				printf(" (rc21).\n");
				return AppendFirst(NewList(),a2);
				}
			co1=i1;
			co2=i2;
			collab=ListNth(l1,(int)i1);
			}
		}

	if(collab==0)
		{
		printf("Internal error in vertex ");
		WriteVertex(CompoundArg1(a2));
		printf(" (rc22).\n");
		return AppendFirst(NewList(),a2);
		}

	i1=GetAtomProperty(CompoundArg1(csp1),A_COLOR);
	i2=GetAtomProperty(CompoundArg1(csp2),A_COLOR);
	
	if(i1==A_COLOR_EPS && i2==A_COLOR_EPS)
		{
		ctype=1;
		if(!is_color_eps(CompoundArg1(csp1)))
			{
			co1=csp1;
			csp1=csp2;
			csp2=co1;
			}
		goto cnt4;
		}
	
	if(i1==A_COLOR_EPS && i2==A_COLOR_LAMBDA)
		{
		if(co2==1)
			{
			co1=csp1;
			csp1=csp2;
			csp2=co1;
			}
		ctype=1;
		goto cnt4;
		}
		
	if(i1==A_COLOR_LAMBDA && i2==A_COLOR_EPS)
		{
		if(co1==2)
			{
			co1=csp1;
			csp1=csp2;
			csp2=co1;
			}
		ctype=1;
		goto cnt4;
		}
		
		
	if(i1==A_COLOR_F && i2==A_COLOR_F)
		{
		ctype=3;
		goto cnt4;
		}
		
	if(i1==A_COLOR_LAMBDA && i2==A_COLOR_LAMBDA)
		{
		if(co1==3 && co2==3)
			{
			ctype=3;
			goto cnt4;
			}
		if(co1==1 && co2==2)
			{
			ctype=1;
			goto cnt4;
			}
		if(co1==2 && co2==1)
			{
			co1=csp1;
			csp1=csp2;
			csp2=co1;
			ctype=1;
			goto cnt4;
			}
		printf("Internal error in vertex ");
		WriteVertex(CompoundArg1(a2));
		printf(" (rc23).\n");
		return AppendFirst(NewList(),a2);
		}
	if(i1==A_COLOR_L6 && i2==A_COLOR_L6)
		{
		  
		if(co1==3 && co2==3)
			{
			ctype=3;
			goto cnt4;
			}
		if(co1==1 && co2==2)
			{
			ctype=6;
			goto cnt4;
			}
		if(co1==2 && co2==1)
			{
			co1=csp1;
			csp1=csp2;
			csp2=co1;
			ctype=6;
			goto cnt4;
			}
		printf("Internal error in vertex ");
		WriteVertex(CompoundArg1(a2));
		printf(" (rc23).\n");
		return AppendFirst(NewList(),a2);
		}
	if((i1==A_COLOR_LAMBDA && co1!=3) || (i2==A_COLOR_LAMBDA && co2!=3))
		{
		printf("Internal error in vertex ");
		WriteVertex(CompoundArg1(a2));
		printf(" (rc24).\n");
		return AppendFirst(NewList(),a2);
		}
	ctype=3;
cnt4:

/*	printf("csp1: "); WriteTerm(csp1);
	printf("\ncsp2: "); WriteTerm(csp2); WriteTerm(collab); puts("");
*/

	/********* sep prtcls *************/


	prt=CompoundArg1(a2);
	prt1=prt2=pl1=pl2=osp1=osp2=NewList();
	l1=pl;
	l2=pil;

	co1=CompoundArg2(csp1);
	co2=CompoundArg2(csp2);

	while(!is_empty_list(l1))
		{
		Label llab;
		llab=CompoundArg2(ListFirst(l2));

		if(ListMember(co1,llab))
			pl1=AppendLast(pl1,ListFirst(l1));
		else
			{
			if(!ListMember(co2,llab))
				{
				printf("Internal error in vertex ");
				WriteVertex(CompoundArg1(a2));
				printf(" (rc25).\n");
				return AppendFirst(NewList(),a2);
				}
			pl2=AppendLast(pl2,ListFirst(l1));
			}
		l1=ListTail(l1);
		l2=ListTail(l2);
		}

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
				goto cnt1;
				}
			if(l2==CompoundArg1(ListFirst(pl2)) ||
				l2==CompoundArg1(ListNth(pl2,2)))
				{
#ifdef DEB
				puts(" -> 2");
#endif
				osp2=AppendLast(osp2,t);
				goto cnt1;
				}
			printf("Internal error in vertex ");
			WriteVertex(CompoundArg1(a2));
			printf(" (rc26).\n");
			return AppendFirst(NewList(),a2);
			}
		if(CompoundName(t)==OPR_SPECIAL &&
				(CompoundArg1(t)==A_GAMMA || CompoundArg1(t)==A_GAMMA5||
				CompoundArg1(t)==A_GAMMAM ||CompoundArg1(t)==A_GAMMAP))
			{
			if(CompoundName(ListFirst(pl1))==OPR_SPINOR &&
				CompoundName(ListNth(pl1,2))==OPR_SPINOR)
				{
				if(CompoundName(ListFirst(pl2))==OPR_SPINOR ||
					CompoundName(ListNth(pl2,2))==OPR_SPINOR)
					{
					printf("Internal error in vertex ");
					WriteVertex(CompoundArg1(a2));
					printf(" (rc26).\n");
					return AppendFirst(NewList(),a2);
					}
#ifdef DEB
				puts(" -> 1");
#endif
				osp1=AppendLast(osp1,t);
				goto cnt1;
				}
			if(CompoundName(ListFirst(pl2))==OPR_SPINOR &&
				CompoundName(ListNth(pl2,2))==OPR_SPINOR)
				{
				if(CompoundName(ListFirst(pl1))==OPR_SPINOR ||
					CompoundName(ListNth(pl1,2))==OPR_SPINOR)
					{
					printf("Internal error in vertex ");
					WriteVertex(CompoundArg1(a2));
					printf(" (rc27).\n");
					return AppendFirst(NewList(),a2);
					}
#ifdef DEB
				puts(" -> 2");
#endif
				osp2=AppendLast(osp2,t);
				goto cnt1;
				}
			printf("Internal error in vertex ");
			WriteVertex(CompoundArg1(a2));
			printf(" (rc28).\n");
			return AppendFirst(NewList(),a2);
			}
		printf("Internal error in vertex ");
		WriteVertex(CompoundArg1(a2));
		printf(" (rc29 -- "); WriteTerm(CompoundName(t));
		printf(":"); WriteTerm(CompoundArg1(t)); puts(").");
		return AppendFirst(NewList(),a2);
cnt1:	l1=ListTail(l1);
		}

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

	l2=CompoundArg2(csp1);
	while(!is_empty_list(l2))
		{
		i1=XorList(i1,ListFirst(l2));
		l2=ListTail(l2);
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

	l2=CompoundArg2(csp2);
	while(!is_empty_list(l2))
		{
		i2=XorList(i2,ListFirst(l2));
		l2=ListTail(l2);
		}

	if(ListLength(i1)!=ListLength(i2))
		{
		printf("Internal error in vertex ");
		WriteVertex(CompoundArg1(a2));
		printf(" (rc2A).\n");
		return AppendFirst(NewList(),a2);
		}

	if(ListLength(i1)>3)
		{
		if(check_splitting_c2)
			return 0;
		printf("Can not split 4-color vertex ");
		WriteVertex(CompoundArg1(a2));
		printf(": too many vector indices.\n");
		/*WriteTerm(csp1);printf(",");WriteTerm(csp2);puts("");
		WriteTerm(a2);puts("");*/
		return 0 /*AppendFirst(NewList(),a2)*/ ;
		}
	FreeAtomic(i2);
	
	RemoveList(csp);
	RemoveList(osp);
	RemoveList(pl);
	FreeAtomic(pil);

	if(check_splitting_c2)
		return a2;	

/*	WriteTerm(i1);puts(": im indices");
*/
	l1=i1;
	while(!is_empty_list(l1))
		{
		if(ListFirst(l1)==collab)
			{
			i1=CutFromList(i1,l1);
			break;
			}
		l1=ListTail(l1);
		}

	i1=AppendLast(i1,collab);

/*	WriteTerm(i1);puts(": im indices");*/

	/**********    splitting  ***********/
/*
	printf("color prt1: ");
	WriteTerm(pl1); WriteTerm(prt1);
	printf("\ncolor prt2: ");
	WriteTerm(pl2); WriteTerm(prt2); puts("");
*/
	anti=check_charge(pl1,prt1,osp1, csp1, pl2,prt2,osp2,csp2,i1);

	if(opNeutrC8==0 && ctype==3 && anti)
		{
		warnc8++;
		if(warnc8==5) puts("More problematic vertices follow...");
		if(warnc8<5)
		{
		printf("Warning: ");
		WriteVertex(CompoundArg1(a2));
		printf(" splitting this vertex can make non-hermitian lagrangian.\n");
		printf("Use 'option NeutralC8=1.' ");
		printf("which however may introduce additional interaction.\n");
		}
		}

	l1=ConsumeCompoundArg(a2,5);
	m=ListFirst(l1);
	RemoveList(l1);
	prt=ConsumeCompoundArg(a2,1);
	a22=MakeCompound(A_ALG2,5);
	co2=0;
	l1=ConsumeCompoundArg(m,3);
	RemoveList(l1);

	split_m2(m,&m2,&m22,&anti, ListLength(i1)-1);
	FreeAtomic(m);
	l2=mk_im_field(ctype, ListLength(i1)-1, (ctype==3&&opNeutrC8)?0:anti, 
	  ConcatList(CopyTerm(prt1),CopyTerm(prt2)));
	
	pl1=ConcatList(pl1,osp1);
	pl1=AppendLast(pl1,csp1);
	co1=NewLabel();
	prt1=AppendLast(prt1,MakeCompound2(OPR_DIV,ListFirst(l2),co1));
	switch(ListLength(i1)-1)
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
	pl2=AppendLast(pl2,csp2);
	co1=NewLabel();
	prt2=AppendLast(prt2,MakeCompound2(OPR_DIV,ListNth(l2,2),co1));
	RemoveList(l2);
	switch(ListLength(i1)-1)
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
	printf("color prt1: ");
	WriteTerm(pl1); WriteTerm(prt1);
	printf("\ncolor prt2: ");
	WriteTerm(pl2); WriteTerm(prt2); puts("");

	puts("\n"); WriteTerm(a2); puts("\n"); WriteTerm(a22); puts("\n--2--");

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

static List color_reduce_0(Term a2, List pl, List pil, List osp)
	{
	List pl1, pl2;
	List prt, prt1, prt2;
	List osp1, osp2;
	List l1,l2;
	List i1,i2;
	List co1,co2;
	Term m2,m22,a22,m;
	int anti;


#ifdef DEB

	printf("color prt: "); WriteTerm(pl);
	printf(" with ");      WriteTerm(pil);
	printf("\noth sp:    ");  WriteTerm(osp); puts("\n");

#endif


	/*WriteTerm(a2); puts("");*/

	/*******    sep prt ***********/

	prt=CompoundArg1(a2);
	prt1=prt2=pl1=pl2=osp1=osp2=NewList();
	pl1=AppendLast(NewList(),ListFirst(pl));
	l1=ListTail(pl);
	l2=ListTail(pil);
	while(!is_empty_list(l1))
		{
		if(CompoundArg2(ListFirst(l2))==CompoundArg2(ListFirst(pil)))
			pl1=AppendLast(pl1,ListFirst(l1));
		else
			pl2=AppendLast(pl2,ListFirst(l1));
		l1=ListTail(l1);
		l2=ListTail(l2);
		}
	RemoveList(pl);
	FreeAtomic(pil);

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
printf("prt1: "); WriteTerm(prt1); WriteTerm(pl1); puts("");
printf("prt2: "); WriteTerm(prt2); WriteTerm(pl2); puts("");
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
				goto cnt1;
				}
			if(l2==CompoundArg1(ListFirst(pl2)) ||
				l2==CompoundArg1(ListNth(pl2,2)))
				{
#ifdef DEB
				puts(" -> 2");
#endif
				osp2=AppendLast(osp2,t);
				goto cnt1;
				}
			printf("Internal error in vertex ");
			WriteVertex(CompoundArg1(a2));
			printf(" (rc02).\n");
			return AppendFirst(NewList(),a2);
			}
		if(CompoundName(t)==OPR_SPECIAL &&
				(CompoundArg1(t)==A_GAMMA || CompoundArg1(t)==A_GAMMA5 ||
				CompoundArg1(t)==A_GAMMAM ||CompoundArg1(t)==A_GAMMAP ))
			{
			Label i1,i2;
			i1=ListFirst(CompoundArg2(t));
			i2=ListNth(CompoundArg2(t),2);
			if( (i1==ListFirst(CompoundArg2(ListFirst(pl1))) ||
			     i2==ListFirst(CompoundArg2(ListNth(pl1,2))) ) ||
			    (i2==ListFirst(CompoundArg2(ListFirst(pl1))) ||
			     i1==ListFirst(CompoundArg2(ListNth(pl1,2))) ))
			     {
				 if( (i1==ListFirst(CompoundArg2(ListFirst(pl2))) ||
			     i2==ListFirst(CompoundArg2(ListNth(pl2,2))) ) ||
			    (i2==ListFirst(CompoundArg2(ListFirst(pl2))) ||
			     i1==ListFirst(CompoundArg2(ListNth(pl2,2))) ))
			     	goto fail;
#ifdef DEB
				puts(" -> 1");
#endif
				osp1=AppendLast(osp1,t);
				goto cnt1;
				}

			if( (i1==ListFirst(CompoundArg2(ListFirst(pl2))) ||
			     i2==ListFirst(CompoundArg2(ListNth(pl2,2))) ) ||
			    (i2==ListFirst(CompoundArg2(ListFirst(pl2))) ||
			     i1==ListFirst(CompoundArg2(ListNth(pl2,2))) ))
			     {
				 if( (i1==ListFirst(CompoundArg2(ListFirst(pl1))) ||
			     i2==ListFirst(CompoundArg2(ListNth(pl1,2))) ) ||
			    (i2==ListFirst(CompoundArg2(ListFirst(pl1))) ||
			     i1==ListFirst(CompoundArg2(ListNth(pl1,2))) ))
				 	goto fail;
#ifdef DEB
				puts(" -> 2");
#endif
				osp2=AppendLast(osp2,t);
				goto cnt1;
				}


			printf("Internal error in vertex ");
			WriteVertex(CompoundArg1(a2));puts("");
			WriteTerm(t);puts("");
			printf(" (rc05).\n");
			return AppendFirst(NewList(),a2);
			}
fail:
		printf("Internal error in vertex ");
		WriteVertex(CompoundArg1(a2));
		printf(" (rc06 -- "); WriteTerm(CompoundName(t));
		printf(":"); WriteTerm(CompoundArg1(t)); puts(").");
		return AppendFirst(NewList(),a2);
cnt1:	l1=ListTail(l1);
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

#ifdef DEB
	WriteTerm(i1);WriteTerm(i2);puts(": im indices");
#endif

	if(ListLength(i1)!=ListLength(i2) ||
		ListLength(i1)>2)
		{
		printf("Internal error in vertex ");
		WriteVertex(CompoundArg1(a2));
		printf(" (rc04).\n");
		return AppendFirst(NewList(),a2);
		}
	FreeAtomic(i2);

	/**********    splitting  ***********/


	anti=check_charge(pl1,prt1,osp1, 0, pl2,prt2,osp2, 0, i1);


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

	l2=mk_im_field(0, ListLength(i1), anti, ConcatList(CopyTerm(prt1),CopyTerm(prt2)));

	pl1=ConcatList(pl1,osp1);
	co1=NewLabel();
	prt1=AppendLast(prt1,MakeCompound2(OPR_DIV,ListFirst(l2),co1));
	switch(ListLength(i1))
		{
		case 0: co2=MakeCompound2(OPR_SCALAR,co1,NewList()); break;
		case 1: co2=MakeCompound2(OPR_VECTOR,co1,CopyTerm(i1)); break;
		case 2: co2=MakeCompound2(OPR_TENSOR,co1,CopyTerm(i1)); break;
		}
	pl1=AppendLast(pl1,co2);
	SetCompoundArg(m2,3,pl1);
	SetCompoundArg(a2,1,prt1);
	SetCompoundArg(a2,5,AppendLast(NewList(),m2));

	pl2=ConcatList(pl2,osp2);
	co1=NewLabel();
	prt2=AppendLast(prt2,MakeCompound2(OPR_DIV,ListNth(l2,2),co1));
	RemoveList(l2);
	switch(ListLength(i1))
		{
		case 0: co2=MakeCompound2(OPR_SCALAR,co1,NewList()); break;
		case 1: co2=MakeCompound2(OPR_VECTOR,co1,i1); break;
		case 2: co2=MakeCompound2(OPR_TENSOR,co1,i1); break;
		}
	pl2=AppendLast(pl2,co2);
	SetCompoundArg(m22,3,pl2);
	SetCompoundArg(a22,1,prt2);
	SetCompoundArg(a22,5,AppendLast(NewList(),m22));

#ifdef DEB
	printf("color prt1: ");
	WriteTerm(pl1); WriteTerm(prt1);
	printf("\ncolor prt2: ");
	WriteTerm(pl2); WriteTerm(prt2); puts("");
	puts("\n"); WriteTerm(a2); puts("\n"); WriteTerm(a22); puts("--0--");
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


static int lab2col(Label lab, List pl, int *pos)
	{
	while(!is_empty_list(pl))
		{
		if(CompoundArg2(ListFirst(pl))==lab)
			{
			Term prp;
			prp=GetAtomProperty(CompoundArg1(ListFirst(pl)),A_COLOR);
			if(prp==0)
				return 0;
			*pos=(int)IntegerValue(CompoundArg2(prp));
			return (int)IntegerValue(CompoundArg1(prp));
			}
		pl=ListTail(pl);
		}
	puts("Internal error (l2c).");
	return 0;
	}

extern Atom ColorK6, ColorK6b;
	
List color_reduce(Term a2)
	{
	List m2, cp_l, op_l, cs_l, os_l, pi_l, l1,l2;
	int colc6=0, colc6b=0;
	
	cp_l=pi_l=op_l=cs_l=os_l=NewList();

/*	WriteTerm(CompoundArg1(a2)); puts("");
*/	
	m2=ListFirst(CompoundArgN(a2,5));
	l1=CompoundArgN(m2,3);
	while(!is_empty_list(l1))
		{
		Term t1;
		t1=ListFirst(l1);
		if(CompoundName(t1)==OPR_TENSOR || CompoundName(t1)==OPR_VECTOR ||
		   CompoundName(t1)==OPR_SPINOR || CompoundName(t1)==OPR_SCALAR )
		   	{
		   	int pos, col;
		   	col=lab2col(CompoundArg1(t1), CompoundArg1(a2), &pos);
		   	if(col)
		   		{
		   		cp_l=AppendLast(cp_l,t1);
		   		pi_l=AppendLast(pi_l,MakeCompound2(OPR_DIV,NewInteger(col),
		   			ListNth(CompoundArg2(t1),pos)));
		   		goto cnt1;
		   		}
		   	op_l=AppendLast(op_l,t1);
		   	goto cnt1;
		   	}
		if(CompoundName(t1)==OPR_SPECIAL && GetAtomProperty(CompoundArg1(t1),
				A_COLOR)!=0)
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
			printf(": illegal color structure (no cp, some cs).\n");
			}
		RemoveList(op_l);
		RemoveList(os_l);
		return AppendFirst(NewList(),a2);
		}

	if(NoColors && FAOutput && cs_l)
	{
		RemoveList(op_l);
		RemoveList(os_l);
		RemoveList(cp_l);
		RemoveList(cs_l);
		FreeAtomic(pi_l);
		FreeAtomic(a2);
		return NewList();
	}
		
	if(ListLength(cp_l)==1)
		{
		printf("Error in vertex "); WriteVertex(CompoundArg1(a2));
		printf(": illegal color structure (single color particle).\n");
		RemoveList(op_l);
		RemoveList(os_l);
		RemoveList(cp_l);
		RemoveList(cs_l);
		FreeAtomic(pi_l);
		return AppendFirst(NewList(),a2);
		}

	if(ListLength(cp_l)==2)
		{
		if(ListLength(cs_l))
			{
			printf("Error in vertex "); WriteVertex(CompoundArg1(a2));
			printf(": illegal color structure. (2cp, some cs)\n");
			}
		RemoveList(op_l);
		RemoveList(os_l);
		RemoveList(cp_l);
		RemoveList(cs_l);
		FreeAtomic(pi_l);
		return AppendFirst(NewList(),a2);
		}

	if(ListLength(cp_l)==3)
		{
		if(ListLength(cs_l)!=1)
			{
			printf("Error in vertex "); WriteVertex(CompoundArg1(a2));
			printf(": illegal color structure. (3cp, !1 cs)\n");
			}
		l2=CompoundArg2(ListFirst(cs_l));
		if(ListFirst(l2)==ListNth(l2,2) || ListFirst(l2)==ListNth(l2,3) ||
				ListNth(l2,2)==ListNth(l2,3))
			{
			printf("Error in vertex "); WriteVertex(CompoundArg1(a2));
			printf(": illegal color structure. (self-conv in cs)\n");
			}
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
		printf(": illegal color structure. (too many prt)\n");
		RemoveList(op_l);
		RemoveList(os_l);
		RemoveList(cp_l);
		RemoveList(cs_l);
		FreeAtomic(pi_l);
		return AppendFirst(NewList(),a2);
		}

    if(GetAtomProperty(A_COLOR,A_S_COLOR)==NewInteger(0))
        return 0;

	if(is_empty_list(cs_l))
		return color_reduce_0(a2, cp_l, pi_l, os_l);

	if(ListLength(cs_l)==2)
		return color_reduce_2(a2, cp_l, pi_l, cs_l, os_l);

  
	 for(l1=cs_l;l1;l1=ListTail(l1))
	 {
	   if(CompoundArg1(ListFirst(l1))==ColorK6)
	     colc6++;
	   if(CompoundArg1(ListFirst(l1))==ColorK6b)
	     colc6b++;
	 }
	
	if(ListLength(cs_l)!=4 || !colc6 || !colc6b)
	{
	printf("Error in vertex "); WriteVertex(CompoundArg1(a2));
	printf(": illegal color structure. (incorrect #cs)\n");
	}
	RemoveList(os_l);
	RemoveList(cp_l);
	RemoveList(cs_l);
	FreeAtomic(pi_l);


	return AppendFirst(NewList(),a2);
	}

int opSplitCol2=1;

int need_col_rdc(Term a2)
{
	
	List l1,l2;
	int cp_no, cs_no;
	
	cp_no=0;
	
	if(WriteColors || (NoColors && !FAOutput) || opSplitCol2==0)
		return 0;
	
	for(l1=CompoundArg1(a2);l1;l1=ListTail(l1))
		if(GetAtomProperty(CompoundArg1(ListFirst(l1)),A_COLOR))
			cp_no++;
	
	if(cp_no>=4)
		return 1;

	for(l1=CompoundArgN(a2,5);l1;l1=ListTail(l1))
	{
		cs_no=0;
		for(l2=CompoundArgN(ListFirst(l1),3);l2;l2=ListTail(l2))
			if(is_atom(CompoundArg1(ListFirst(l2))) &&
					GetAtomProperty(CompoundArg1(ListFirst(l2)),A_COLOR))
				cs_no++;
		if((cp_no<3 && cs_no>0) || (cp_no==3 && cs_no>1))
		{
			printf("Error in vetrex ");
			WriteVertex(CompoundArg1(a2));
			printf(": illegal color structure (%d cp, %d cs)\n",cp_no,cs_no);
			return 0;
		}
		if(cs_no && NoColors && FAOutput)
			return 1;
	}
	
	return 0;
}
	
