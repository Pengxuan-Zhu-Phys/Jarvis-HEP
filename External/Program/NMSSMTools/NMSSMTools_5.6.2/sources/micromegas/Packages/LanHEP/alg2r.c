#include "lanhep.h"

static Integer ferm_ind[8];

extern int WriteColors, opRemDotWithFerm;

int is_color_d(Atom);

static char mb_buf[32];

static Atom mk_gamma(int t, int n)
	{
	if(t==-1)
		return NewAtom("G5",0);
	if(t==-2)
		return NewAtom("(1+G5)",0);
	if(t==-3)
		return NewAtom("(1-G5)",0);
	
	sprintf(mb_buf,"G(%c%d)",t?(t==2?'M':'p'):'m',n);
	return NewAtom(mb_buf,0);
	}

static Atom mk_cnvl(int t1, int n1, int t2, int n2)
	{
	if(t1>t2)
		{
		int tmp;
		tmp=t1;
		t1=t2;
		t2=tmp;
		tmp=n1;
		n1=n2;
		n2=tmp;
		}
	if(t1==t2 && n1>n2)
		{
		int tmp;
		tmp=n1;
		n1=n2;
		n2=tmp;
		}
	sprintf(mb_buf,"%c%d.%c%d",t1?(t1==2?'M':'p'):'m',n1,
			t2?(t2==2?'M':'p'):'m',n2);

	return NewAtom(mb_buf,0);
	}
	
static Atom mk_cnvlc(int n1, int n2)
	{
	if(n1>n2)
		{
		int tmp;
		tmp=n1;
		n1=n2;
		n2=tmp;
		}
	sprintf(mb_buf,"c%d.c%d",n1,n2);

	return NewAtom(mb_buf,0);
	}

static Atom mk_color(char t, int i1, int i2, int i3)
{
	sprintf(mb_buf,"%c(c%d,c%d,c%d)",t,i1,i2,i3);
	return NewAtom(mb_buf,0);
}
		
static Atom mk_eps(Atom p1, Atom p2, Atom p3, Atom p4)
	{
	sprintf(mb_buf,"eps(%s,%s,%s,%s)",AtomValue(p1),AtomValue(p2),AtomValue(p3),AtomValue(p4));

	return NewAtom(mb_buf,0);
	}

static int check_mpg(Term m2)
{
	List l1;
	int mp=0, g=0;
	for(l1=CompoundArgN(m2,3);l1;l1=ListTail(l1))
	{
		char *a;
		a=AtomValue(ListFirst(l1));
		if(a[0]=='G' || a[3]=='G' || ferm_ind[0])
			g++;
		if(a[0]=='m' || a[0]=='M' || a[0]=='p')
			mp++;
	}
	return g?mp:0;
}


static int corr_mpg(Term m2)
{
	List l1;
	int ret=0;
	
	for(l1=CompoundArgN(m2,3);l1;l1=ListTail(l1))
	{
		char *a;
		a=AtomValue(ListFirst(l1));
		if(a[0]=='m' || a[0]=='M' || a[0]=='p')
		{
			char buf[28];
			sprintf(buf,"(G(%2.2s)*G(%2.2s)+G(%2.2s)*G(%2.2s))",
					a,a+3,a+3,a);
			ChangeList(l1,NewAtom(buf,0));
			ret++;
		}
	}
	return ret;
}


static int col_ind(List pl, Integer ii)
{
	List l2;
	int n=0;
	int ppp;
	for(l2=pl,ppp=1;l2;l2=ListTail(l2),ppp++)
		{
		Term t1=CompoundArg2(ListFirst(l2));
		if(CompoundName(t1)==OPR_SCALAR && CompoundArg1(t1) 
				&& ListFirst(CompoundArg1(t1))==ii)
			n=ppp;
		if((CompoundName(t1)==OPR_VECTOR || CompoundName(t1)==OPR_SPINOR) 
				&& ListLength(CompoundArg1(t1))==2 
				&& ListFirst(ListTail(CompoundArg1(t1)))==ii)
			n=ppp;
		if((CompoundName(t1)==OPR_TENSOR || CompoundName(t1)==OPR_SPINOR3) 
				&& ListLength(CompoundArg1(t1))==3 
				&& ListFirst(ListTail(ListTail(CompoundArg1(t1))))==ii)
			n=ppp;
		}
	return n;
}
	
	
	
static void m_reduce(Term m2, List pl, int pis, int fr)
	{
	List l,l1,lg,ls;

/*
	 WriteTerm(pl); WriteTerm(m2); printf("  fi: %d %d\n", IntegerValue(ferm_ind[0]), 
		IntegerValue(ferm_ind[1]));
*/

	l1=NewList();
	lg=NewList();
	
	ls=l=ConsumeCompoundArg(m2,3);
	
	while(!is_empty_list(l))
		{
		Term t;
		t=ListFirst(l);
		/*ChangeList(l,0);*/
		if(t==0) goto cnt;

		if(WriteColors && CompoundName(t)==OPR_SPECIAL &&
				GetAtomProperty(CompoundArg1(t),A_COLOR)==A_COLOR_LAMBDA)
		{
			int i1, i2, i3;
			i1=col_ind(pl,ListFirst(CompoundArg2(t)));
			i2=col_ind(pl,ListFirst(ListTail(CompoundArg2(t))));
			i3=col_ind(pl,ListFirst(ListTail(ListTail(CompoundArg2(t)))));
			l1=AppendLast(l1,mk_color('L',i1,i2,i3));
			goto cnt;
		}
		if(WriteColors && CompoundName(t)==OPR_SPECIAL && 
				GetAtomProperty(CompoundArg1(t),A_COLOR)==A_COLOR_F)
		{
			int i1, i2, i3;
			i1=col_ind(pl,ListFirst(CompoundArg2(t)));
			i2=col_ind(pl,ListFirst(ListTail(CompoundArg2(t))));
			i3=col_ind(pl,ListFirst(ListTail(ListTail(CompoundArg2(t)))));
			l1=AppendLast(l1,mk_color('F',i1,i2,i3));
			goto cnt;
		}

		if(WriteColors && CompoundName(t)==OPR_SPECIAL && is_color_d(CompoundArg1(t)))
		{
			int i1, i2, i3;
			i1=col_ind(pl,ListFirst(CompoundArg2(t)));
			i2=col_ind(pl,ListFirst(ListTail(CompoundArg2(t))));
			i3=col_ind(pl,ListFirst(ListTail(ListTail(CompoundArg2(t)))));
			l1=AppendLast(l1,mk_color('D',i1,i2,i3));
			goto cnt;
		}
				
		if(CompoundName(t)==OPR_SPECIAL && is_color_d(CompoundArg1(t)))
		{
			l1=AppendLast(l1,NewAtom("d_SU3",0));
			goto cnt;
		}
		
		
		if(CompoundName(t)==A_MOMENT)
			{
			List l2;
			Integer mi;
			int ppp=1;
			if(CompoundArg1(t)==0) goto cnt;
			mi=ListFirst(CompoundArg2(t));
			l2=pl;
			while(!is_empty_list(l2))
				{
				Term t1;
				t1=CompoundArg2(ListFirst(l2));
				if(CompoundName(t1)==OPR_VECTOR && ListFirst(CompoundArg1(t1))==mi)
					{
					l1=AppendLast(l1,mk_cnvl(1,(int)IntegerValue(CompoundArg1(t)),0,ppp));
					goto cnt;
					}
				if(CompoundName(t1)==OPR_TENSOR && ListFirst(CompoundArg1(t1))==mi)
					{
					l1=AppendLast(l1,mk_cnvl(1,(int)IntegerValue(CompoundArg1(t)),0,ppp));
					goto cnt;
					}
				if(CompoundName(t1)==OPR_TENSOR && ListNth(CompoundArg1(t1),2)==mi)
					{
					l1=AppendLast(l1,mk_cnvl(1,(int)IntegerValue(CompoundArg1(t)),2,ppp));
					goto cnt;
					}
				if(CompoundName(t1)==OPR_SPINOR3 && ListNth(CompoundArg1(t1),2)==mi)
					{
					l1=AppendLast(l1,mk_cnvl(1,(int)IntegerValue(CompoundArg1(t)),0,ppp));
					goto cnt;
					}
				ppp++;
				l2=ListTail(l2);
				}
			l2=ls;
			while(!is_empty_list(l2))
				{
				Term t1;
				t1=ListFirst(l2);
				if(CompoundName(t1)==OPR_SPECIAL && 
					CompoundArg1(t1)==A_GAMMA &&
					ListNth(CompoundArg2(t1),3)==mi)
					{
					Term tg;
					tg=MakeCompound(A_GAMMA,3);
					SetCompoundArg(tg,1,ListFirst(CompoundArg2(t1)));
					SetCompoundArg(tg,2,ListNth(CompoundArg2(t1),2));
					SetCompoundArg(tg,3,mk_gamma(1,(int)IntegerValue(CompoundArg1(t))));
					lg=AppendLast(lg,tg);
					ChangeList(l2,0);
					goto cnt;
					}
				l2=ListTail(l2);
				}
			l2=ListTail(l);
			while(!is_empty_list(l2))
				{
				Term t1;
				t1=ListFirst(l2);
				if(CompoundName(t1)==A_MOMENT && 
					ListFirst(CompoundArg2(t1))==mi)
					{
					l1=AppendLast(l1,mk_cnvl(1,(int)IntegerValue(CompoundArg1(t)),
											1,(int)IntegerValue(CompoundArg1(t1))));
					SetCompoundArg(t1,1,0);
					goto cnt;
					}
				if(CompoundName(t1)==OPR_SPECIAL && CompoundArg1(t1)==A_EPS_V && 
						ListMember(CompoundArg2(t1),mi))
					{
					List l3;
					for(l3=CompoundArg2(t1);l3;l3=ListTail(l3))
						{
						if(ListFirst(l3)==mi)
							ChangeList(l3,NewInteger(-IntegerValue(CompoundArg1(t))));
						}
					goto cnt;
					}
					
				l2=ListTail(l2);
				}
			printf("Internal error in vertex "); WriteVertex(pl); printf(":");
			printf(" (milf)\n");
			goto cnt;				
			}
			
		if(CompoundName(t)==OPR_SPECIAL && CompoundArg1(t)==A_DELTA)
			{
			Integer i1,i2;
			List l2;
			int p1=-1,p2=-1, pp1=-1, pp2=-1;
			int ppp=1;
			i1=ListFirst(CompoundArg2(t));
			i2=ListNth(CompoundArg2(t),2);
			l2=pl;
			while(!is_empty_list(l2))
				{
				Term t1;
				t1=CompoundArg2(ListFirst(l2));
				if(CompoundName(t1)==OPR_VECTOR && ListFirst(CompoundArg1(t1))==i1)
					{
					p1=ppp;
					pp1=0;
					}
				if(CompoundName(t1)==OPR_TENSOR && ListFirst(CompoundArg1(t1))==i1)
					{
					p1=ppp;
					pp1=0;
					}
				if(CompoundName(t1)==OPR_TENSOR && ListNth(CompoundArg1(t1),2)==i1)
					{
					p1=ppp;
					pp1=2;
					} 
				if(CompoundName(t1)==OPR_SPINOR3 && ListNth(CompoundArg1(t1),2)==i1)
					{
					p1=ppp;
					pp1=0;
					} 
				if(CompoundName(t1)==OPR_VECTOR && ListFirst(CompoundArg1(t1))==i2)
					{
					p2=ppp;
					pp2=0;
					}
				if(CompoundName(t1)==OPR_TENSOR && ListFirst(CompoundArg1(t1))==i2)
					{
					p2=ppp;
					pp2=0;
					}
				if(CompoundName(t1)==OPR_TENSOR && ListNth(CompoundArg1(t1),2)==i2)
					{
					p2=ppp;
					pp2=2;
					}
				if(CompoundName(t1)==OPR_SPINOR3 && ListNth(CompoundArg1(t1),2)==i2)
					{
					p2=ppp;
					pp2=0;
					}

				ppp++;
				l2=ListTail(l2);
				}
			if(p1!=-1 && p2!=-1)
				{
				l1=AppendLast(l1,mk_cnvl(pp1,p1,pp2,p2));
				goto cnt;
				}
			if(p1!=-1 || p2!=-1)
				{
				printf("Internal error in vertex "); WriteVertex(pl); printf(":");
				puts(" (dilf)");
				}
			if(!WriteColors)
				goto cnt;
			p1=col_ind(pl, i1);
			p2=col_ind(pl, i2);
			if(p1==0 && p2==0)
				goto cnt;
			if(p1==0 || p2==0)
				{
					printf("Internal error in vertex "); WriteVertex(pl); printf(":");
					puts(" (dilfcc)");
					goto cnt;
				}
				l1=AppendLast(l1,mk_cnvlc(p1,p2));
				goto cnt;
			}
			
		if(CompoundName(t)==OPR_SPECIAL && CompoundArg1(t)==A_EPS_V)
			{
			List il;
			for(il=CompoundArg2(t);il;il=ListTail(il))
			{
				Term i1;
				List l2;
				int ppp=1,pp1=-1,p1=-1;
				i1=ListFirst(il);
				if(is_integer(i1) &&IntegerValue(i1)<0)
				{
					char cbuf[8];
					sprintf(cbuf,"p%d",(int)-IntegerValue(i1));
					ChangeList(il,NewAtom(cbuf,0));
					continue;
				}
				for(l2=pl;l2;l2=ListTail(l2))
				{
					Term t1;
					t1=CompoundArg2(ListFirst(l2));
					if(CompoundName(t1)==OPR_VECTOR && ListFirst(CompoundArg1(t1))==i1)
						{
						p1=ppp;
						pp1=0;
						break;
						}
					if(CompoundName(t1)==OPR_TENSOR && ListFirst(CompoundArg1(t1))==i1)
						{
						p1=ppp;
						pp1=0;
						break;
						}
					if(CompoundName(t1)==OPR_TENSOR && ListNth(CompoundArg1(t1),2)==i1)
						{
						p1=ppp;
						pp1=2;
						break;
						}
					if(CompoundName(t1)==OPR_SPINOR3 && ListNth(CompoundArg1(t1),2)==i1)
						{
						p1=ppp;
						pp1=0;
						break;
						}
					ppp++;
				}
				
				if(pp1!=-1)
				{
					char cbuf[8];
					sprintf(cbuf,"%c%d",pp1?'M':'m',p1);
					ChangeList(il,NewAtom(cbuf,0));
					continue;
				}
				
				for(l2=ListTail(l);l2;l2=ListTail(l2))
				{
					Term t1;
					t1=CompoundArg1(l2);
					if(CompoundName(t1)==A_MOMENT && CompoundArg1(t1) &&
							ListFirst(CompoundArg2(t1))==i1)
					{
						char cbuf[8];
						sprintf(cbuf,"p%d",(int)IntegerValue(CompoundArg1(t1)));
						ChangeList(il,NewAtom(cbuf,0));
						SetCompoundArg(t1,1,0);
						break;
					}
				}
				if(l2)
					continue;
				printf("Internal error in vertex "); WriteVertex(pl); printf(":");
				puts(" (eilf)");
			}

			il=CompoundArg2(t);
			l1=AppendLast(l1,mk_eps(ListNth(il,1),ListNth(il,2),ListNth(il,3),ListNth(il,4)));
			goto cnt;	
			}
				
		if(CompoundName(t)==OPR_SPECIAL && CompoundArg1(t)==A_GAMMA)
			{
			Integer mi;
			List l2;
			Term tg;
			int ppp=1;
			tg=MakeCompound(A_GAMMA,3);
			SetCompoundArg(tg,1,ListFirst(CompoundArg2(t)));
			SetCompoundArg(tg,2,ListNth(CompoundArg2(t),2));
			mi=ListNth(CompoundArg2(t),3);
			l2=pl;
			while(!is_empty_list(l2))
				{
				Term t1;
				t1=CompoundArg2(ListFirst(l2));
				if(CompoundName(t1)==OPR_VECTOR && ListFirst(CompoundArg1(t1))==mi)
					{
					SetCompoundArg(tg,3,mk_gamma(0,ppp));
					lg=AppendLast(lg,tg);
					goto cnt;
					}
				if(CompoundName(t1)==OPR_TENSOR && ListFirst(CompoundArg1(t1))==mi)
					{
					SetCompoundArg(tg,3,mk_gamma(0,ppp));
					lg=AppendLast(lg,tg);
					goto cnt;
					}
				if(CompoundName(t1)==OPR_TENSOR && ListNth(CompoundArg1(t1),2)==mi)
					{
					SetCompoundArg(tg,3,mk_gamma(2,ppp));
					lg=AppendLast(lg,tg);
					goto cnt;
					}
				if(CompoundName(t1)==OPR_SPINOR3 && ListNth(CompoundArg1(t1),2)==mi)
					{
					SetCompoundArg(tg,3,mk_gamma(0,ppp));
					lg=AppendLast(lg,tg);
					goto cnt;
					}
				ppp++;
				l2=ListTail(l2);
				}
			l2=ListTail(l);
			while(!is_empty_list(l2))
				{
				Term t1;
				t1=ListFirst(l2);
				if(CompoundName(t1)==A_MOMENT && 
					ListFirst(CompoundArg2(t1))==mi)
					{
					SetCompoundArg(tg,3,mk_gamma(1,(int)IntegerValue(CompoundArg1(t1))));
					lg=AppendLast(lg,tg);
					ChangeList(l2,0);
					goto cnt;
					}
				l2=ListTail(l2);
				}
			printf("Internal error in vertex "); WriteVertex(pl); printf(":");
			printf(" (gilf)\n");
			goto cnt;				
			}
		
		if(CompoundName(t)==OPR_SPECIAL && CompoundArg1(t)==A_GAMMA5)
			{
			Term tg;
			tg=MakeCompound(A_GAMMA,3);
			SetCompoundArg(tg,1,ListFirst(CompoundArg2(t)));
			SetCompoundArg(tg,2,ListNth(CompoundArg2(t),2));
			SetCompoundArg(tg,3,mk_gamma(-1,0));
			lg=AppendLast(lg,tg);
			goto cnt;
			}
		if(CompoundName(t)==OPR_SPECIAL && CompoundArg1(t)==A_GAMMAP)
			{
			Term tg;
			tg=MakeCompound(A_GAMMA,3);
			SetCompoundArg(tg,1,ListFirst(CompoundArg2(t)));
			SetCompoundArg(tg,2,ListNth(CompoundArg2(t),2));
			SetCompoundArg(tg,3,mk_gamma(-2,0));
			lg=AppendLast(lg,tg);
			goto cnt;
			}
		if(CompoundName(t)==OPR_SPECIAL && CompoundArg1(t)==A_GAMMAM)
			{
			Term tg;
			tg=MakeCompound(A_GAMMA,3);
			SetCompoundArg(tg,1,ListFirst(CompoundArg2(t)));
			SetCompoundArg(tg,2,ListNth(CompoundArg2(t),2));
			SetCompoundArg(tg,3,mk_gamma(-3,0));
			lg=AppendLast(lg,tg);
			goto cnt;
			}
			
		if(CompoundName(t)==OPR_PARAMETER)
			{
			List pl=ConsumeCompoundArg(m2,2);
			pl=AppendLast(pl,MakeCompound2(OPR_POW,CompoundArg1(t),NewInteger(1)));
			SetCompoundArg(m2,2,pl);
			goto cnt;
			}
	cnt:
		l=ListTail(l);
		}
	
/*	
	if(g5s)
		l1=AppendLast(l1,g5s);
	if(gs)
		l1=AppendLast(l1,gs);
	WriteTerm(l1); WriteTerm(pl); puts("");
	
	WriteTerm(ls);puts("");
	WriteTerm(l1);puts("\n");
*/
				
	FreeAtomic(ls);
	
	if(fr==0 && !is_empty_list(lg))
		{
		printf("Error in vertex "); WriteVertex(pl); printf(":");
		puts(" gamma-matrices in the vertex without fermions");
		return ;
		}
	
	
	if(!is_empty_list(lg))
		{
		Label lab;
		ls=NewList();
		lab=ferm_ind[0];
	res:	
		l=lg;
		if(lab==ferm_ind[1])
			{
			if(!is_empty_list(lg))
				{
				printf("Error in vertex "); WriteVertex(pl); printf(":");
				printf(" a trace gamma^a^a^mu.\n");
				return ;
				}
			goto mkg_exi;
			}
		if(is_empty_list(lg))
			{
			printf("Error in vertex "); WriteVertex(pl); printf(":");
			printf(" bad spinor structure. ");
			puts("");
			return ;
			}
		while(!is_empty_list(l))
			{
			if(CompoundArg1(ListFirst(l))==lab)
				{
				ls=AppendLast(ls,CompoundArgN(ListFirst(l),3));
				lab=CompoundArg2(ListFirst(l));
				lg=CutFromList(lg,l);
				goto res;
				}
			l=ListTail(l);
			}
		printf("Error in vertex "); WriteVertex(pl); printf(":");
		printf(" bad spinor structure.");
		puts("");
		return ;
	mkg_exi:
		l1=ConcatList(l1,ls);
		}
		
		

			
	SetCompoundArg(m2,3,l1);

/*	WriteTerm(m2);puts("reduced\n");		*/
	}

static int prt_inds(List p)
	{
	int ret;
	ret=0;
	while(!is_empty_list(p))
		{
		ret+=ListLength(CompoundArg1(CompoundArg2(ListFirst(p))));
		p=ListTail(p);
		}
	return ret;
	}
		

static int ex2(int p)
{
	int ret=1;
	while(p) {ret*=2; p--;}
	return ret;
}

void alg2_reduce(Term a2)
	{
	
		List l,prt;
		int pis,fr=0;
		int req_mpg=0;
		
		ferm_ind[0]=0;
		
/*		printf("a2: "); WriteTerm(a2); puts("");*/
		prt=CompoundArg1(a2);
		{	/* Gluon vertices test */
		
		/*
		if(ListLength(prt)==4 && CompoundArg1(ListFirst(prt))==A_GLUON &&
			CompoundArg1(ListNth(prt,2))==A_GLUON && 
			CompoundArg1(ListNth(prt,3))==A_GLUON &&
			CompoundArg1(ListNth(prt,4))==A_GLUON)
				{
				Term m1;
				m1=ListFirst(CompoundArgN(a2,5));
				m1=CompoundArg1(ListFirst(CompoundArg2(m1)));
				SetCompoundArg(a2,3,m1);
				SetCompoundArg(a2,1,A_GLUON);
				SetCompoundArg(a2,2,NewInteger(4));
				SetCompoundArg(a2,5,
					NewAtom("m1.m3*m2.M3-m1.M3*m2.m3",0));
				l2=ListTail(l2);
				continue;
				}*/
		}
		pis=prt_inds(prt);
		l=prt;
		while(!is_empty_list(l))
			{
			if(CompoundName(CompoundArg2(ListFirst(l)))==OPR_SPINOR ||
			CompoundName(CompoundArg2(ListFirst(l)))==OPR_SPINOR3)
				{
				ferm_ind[fr]=ListFirst(CompoundArg1(CompoundArg2(ListFirst(l))));
				fr++;
				}
			l=ListTail(l);
			}
		if(fr!=0 && fr!=2)
			{
			printf("Error in vertex "); WriteVertex(prt); printf(":");
			printf(" %d fermions omitted in CompHEP output\n",fr);
			return;
			}
		l=CompoundArgN(a2,5);
		while(!is_empty_list(l))
			{
			int cno;
			m_reduce(ListFirst(l),prt,pis,fr);
			cno=check_mpg(ListFirst(l));
			if(cno>req_mpg) req_mpg=cno;
			l=ListTail(l);
			}
			
		if(req_mpg && opRemDotWithFerm)
		{
			Term nf;
			int dn;
			/*WriteVertex(CompoundArg1(a2));puts("");*/
			/*WriteTerm(a2);puts("");*/
			nf=CompoundArgN(a2,2);
			dn=(int)IntegerValue(CompoundArg2(nf));
			SetCompoundArg(nf,2,NewInteger(ex2(req_mpg)*dn));
			for(l=CompoundArgN(a2,5);l;l=ListTail(l))
			{
				int crm=corr_mpg(ListFirst(l));
				if(crm!=req_mpg)
					SetCompoundArg(ListFirst(l),1,
						NewInteger(ex2(req_mpg-crm)*
							IntegerValue(CompoundArg1(ListFirst(l)))));
			}
			/*WriteTerm(a2);puts("");*/
		}
		
			
	}		

	
	
	
