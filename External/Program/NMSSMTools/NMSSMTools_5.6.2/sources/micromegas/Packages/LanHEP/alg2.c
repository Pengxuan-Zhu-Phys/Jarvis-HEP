#include <setjmp.h>
#include <string.h>
#include "lanhep.h"


jmp_buf alg2_jmp_buf;

extern Term alg2_mterm_to_2(Term t); /* alg2a.c */
extern void alg2_norm(Term a2);      /* alg2b.c */
extern List alg2_denorm(Term a2);


static void RedSqrt2(List a2l);

extern Term alg2_rem_lambdaf(Term a2);


Term Alg1to2(Term t)
	{
	Term t1,ret=0;
	int alen, acur;
	char buf[40];
	if(t==0)
		return 0;
	
	
		
	if(setjmp(alg2_jmp_buf)!=0)
		return 0;
				
	if(!is_empty_list(CompoundArg2(t)))
		{
		ErrorInfo(324);
		printf("non-scalar lagrangian term\n");
		FreeAtomic(t);
		return 0;
		}
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
/*	DumpList(t1); puts(""); getchar();*/
	t=t1;
	alen=ListLength(t);
	acur=0;
	
	while(!is_empty_list(t))
		{
		Term t2,y,auxt=0;
		acur++;
		sprintf(buf,"mterm_to_2: %d of %d",acur,alen);
		RegisterLine(buf);
		
		t2=ListFirst(t);

		
		y=alg2_mterm_to_2(t2);	
        
		
			
		if(!y)
		{
			t=ListTail(t);
			UnregisterLine();
			continue;
		}
				
		RedSqrt2(y);
		if(TexOutput || FAOutput) auxt=alg2_rem_lambdaf(y);
		else auxt=0;
		
		alg2_norm(y);
		
		if(auxt) alg2_norm(auxt);
	
		ret=AppendFirst(ret,y);
		if(auxt) ret=AppendFirst(ret,auxt);
		
		t=ListTail(t);
		UnregisterLine();
		}
	RemoveList(t1);
/*	t1=ret;
	alen=ListLength(t1);
	acur=0;
	while(!is_empty_list(t1))
		{
		acur++;
		sprintf(buf,"varderiv: %d of %d",acur,alen);
		RegisterLine(buf);
		alg2_varderiv(ListFirst(t1));
		t1=ListTail(t1);
		UnregisterLine();
		}
*/	

	
	return ret;
	}
	
	
static void rdsq(Term dv, int *po)
	{
	long int num, den, gc;
	num=IntegerValue(CompoundArg1(dv));
	den=IntegerValue(CompoundArg2(dv));
	if(*po>1)
		while(*po>1)
			{
			num*=2;
			*po-=2;
			}
	if(*po<0)
		while(*po<0)
			{
			den*=2;
			*po+=2;
			}
	gc=gcf(num,den);
	num/=gc;
	den/=gc;
	SetCompoundArg(dv,1,NewInteger(num));
	SetCompoundArg(dv,2,NewInteger(den));
	}

static void RedSqrt2(List a2)
	{
		Term ml;

		ml=CompoundArgN(a2,5);
		while(!is_empty_list(ml))
			{
			Term m1;
			List l1,l1s,l2;
			m1=ListFirst(ml);
/*			WriteTerm(m1); puts("\n");*/
			l1s=l1=ConsumeCompoundArg(m1,2);
			l2=NewList();
			while(!is_empty_list(l1))
				{
				Term pw;
				int po;
				pw=ListFirst(l1);
				if(CompoundArg1(pw)==A_SQRT2)
					{
					po=(int)IntegerValue(CompoundArg2(pw));
					rdsq(CompoundArg1(m1),&po);
					if(po!=0)
						{
						SetCompoundArg(pw,2,NewInteger(po));
						l2=AppendLast(l2,pw);
						}
					}
				else
					l2=AppendLast(l2,pw);
				l1=ListTail(l1);
				}
			RemoveList(l1s);
			SetCompoundArg(m1,2,l2);			
						
/*			WriteTerm(m1); puts("\n"); getchar();	*/
		
			ml=ListTail(ml);
			}

	}
	
static int ncmp(Term i1, Term i2)
{
	if(IntegerValue(i1)>IntegerValue(i2))
		return 1;
	return -1;
}

void alg2_multbyi(Term a2)
{
	List l,l1,ci;
	int mn=0,in=1;
	if(is_empty_list(CompoundArgN(a2,5))||CompoundArg2(a2)==NewInteger(0))
		return;

	l=ConsumeCompoundArg(a2,3);
	for(l1=l;l1;l1=ListTail(l1))
	{
		if(CompoundArg1(ListFirst(l1))==A_I)
		{
			l=CutFromList(l,l1);
			in++;
			break;
		}
	}
	SetCompoundArg(a2,3,l);
	ci=NewList();
	for(l1=CompoundArg1(a2);l1;l1=ListTail(l1))
	{
		Atom p;
		Term prop;
		p=CompoundArg1(ListFirst(l1));
		prop=GetAtomProperty(p,A_COLOR);
		if(prop && is_compound(prop))
		{
			if(CompoundArg1(prop)==NewInteger(3))
				ci=AppendLast(ci,NewInteger(8));
			else
				ci=AppendLast(ci,NewInteger(3));
		}
	}
	ci=SortedList(ci,ncmp);
	
	if(ListLength(ci)==4 && ListFirst(ci)==NewInteger(8))
	{
		in+=2;
	}

	if(ListLength(ci)==3)
	{
		if(ListFirst(ci)==NewInteger(3))
			mn++;
		else
			in++;
	}
	FreeAtomic(ci);
	while(in>1)
		in-=2,mn++;
	if(mn==2 || mn==4)
		mn=0;
	if(in)
	{
		l=ConsumeCompoundArg(a2,3);
		l=AppendFirst(l,MakeCompound2(OPR_POW,A_I,NewInteger(1)));
		SetCompoundArg(a2,3,l);
	}
	if(mn)
		SetCompoundArg(CompoundArg2(a2),1,
				NewInteger(-IntegerValue(CompoundArg1(CompoundArg2(a2)))));
	return;
	 
}

int kill_gamma_pm=0;

void alg2_red_1pm5(Term a2)
{
	List l1,l2;
	int pmf=0;
	int chf=0;

	if(kill_gamma_pm)
	{
		alg2_kill_gpm(a2);
		return;
	}
	
	l1=ConsumeCompoundArg(a2,5);
rpt:
	for(l2=l1;l2;l2=ListTail(l2))
	{
		List l3;
		for(l3=CompoundArgN(ListFirst(l2),3);l3;l3=ListTail(l3))
		{
			Term t;
			t=CompoundArg1(ListFirst(l3));
			if(t==A_GAMMAM || t==A_GAMMAP)
			{
				List l22;
				pmf=1;
				if(t==A_GAMMAM)
					SetCompoundArg(ListFirst(l3),1,A_GAMMAP);
				else
					SetCompoundArg(ListFirst(l3),1,A_GAMMAM);
				for(l22=ListTail(l2);l22;l22=ListTail(l22))
				{
					if(EqualTerms(ListFirst(l2),ListFirst(l22)))
					{
						Term lab1, lab2;
						List l33;
						chf=1;
						l1=CutFromList(l1,l22);
						SetCompoundArg(ListFirst(l3),1,A_DELTA);
						lab1=ListFirst(CompoundArg2(ListFirst(l3)));
						lab2=ListFirst(ListTail(CompoundArg2(ListFirst(l3))));
						for(l33=CompoundArgN(ListFirst(l2),3);l33;l33=ListTail(l33))
						{
							if(l33==l3)
								continue;
							if(CompoundArg1(ListFirst(l33))==A_GAMMA && 
									ListNth(CompoundArg2(ListFirst(l33)),2)==lab1)
							{
								ChangeList(ListTail(CompoundArg2(ListFirst(l33))),lab2);
								l33=ConsumeCompoundArg(ListFirst(l2),3);
								l33=CutFromList(l33,l3);
								SetCompoundArg(ListFirst(l2),3,l33);
								goto rpt;
							}
						}
						goto rpt;
					}
				}
				SetCompoundArg(ListFirst(l3),1,t);
			}
		}
	}

		
rpt1:
	for(l2=l1;l2;l2=ListTail(l2))
	{
		List l3;
		for(l3=CompoundArgN(ListFirst(l2),3);l3;l3=ListTail(l3))
		{
			Term t;
			t=CompoundArg1(ListFirst(l3));
			if(t==A_GAMMAP)
			{
				List l22;
				SetCompoundArg(ListFirst(l3),1,A_GAMMAM);
				SetCompoundArg(ListFirst(l2),1,NewInteger(-IntegerValue(CompoundArg1(ListFirst(l2)))));
				for(l22=l1;l22;l22=ListTail(l22))
				{
					if(l2==l22)
						continue;
					if(EqualTerms(ListFirst(l2),ListFirst(l22)))
					{
						chf=1;
						SetCompoundArg(ListFirst(l2),1,NewInteger(-IntegerValue(CompoundArg1(ListFirst(l2)))));
						l1=CutFromList(l1,l22);
						SetCompoundArg(ListFirst(l3),1,A_GAMMA5);
						goto rpt1;
					}
				}
				SetCompoundArg(ListFirst(l3),1,A_GAMMAP);
				SetCompoundArg(ListFirst(l2),1,NewInteger(-IntegerValue(CompoundArg1(ListFirst(l2)))));
			}
		}
	}
	
		
	if(pmf && FAOutput==0 && TexOutput==0)
	{
		SetCompoundArg(CompoundArg2(a2),2,NewInteger(2*IntegerValue(CompoundArg2(CompoundArg2(a2)))));
		for(l2=l1;l2;l2=ListTail(l2))
		{
			List l3;
			pmf=0;
			for(l3=CompoundArgN(ListFirst(l2),3);l3;l3=ListTail(l3))
			{
				Term t;
				t=CompoundArg1(ListFirst(l3));
				if(t==A_GAMMAM || t==A_GAMMAP)
					pmf=1;
			}
			if(pmf==0)
				SetCompoundArg(ListFirst(l2),1,NewInteger(2*IntegerValue(CompoundArg1(ListFirst(l2)))));
		}
	}
	
	if(chf)
	{
		List l3;
		for(l2=l1;l2;l2=ListTail(l2))
		for(l3=ListTail(l2);l3;l3=ListTail(l3))
		{
			if(ListFirst(l2)==0 || ListFirst(l3)==0)
				continue;
			if(EqualTerms(CompoundArg2(ListFirst(l2)),
					CompoundArg2(ListFirst(l3))) &&
				EqualTerms(CompoundArgN(ListFirst(l2),3),
					CompoundArgN(ListFirst(l3),3)))
			{
				long int i1, i2;
				i1=IntegerValue(CompoundArg1(ListFirst(l2)));
				i2=IntegerValue(CompoundArg1(ListFirst(l3)));
				SetCompoundArg(ListFirst(l2),1,NewInteger(i1+i2));
				FreeAtomic(ListFirst(l3));
				ChangeList(l3,0);
				chf=0;
			}
		}
		
		if(chf==0)
		for(l2=l1;l2;l2=ListTail(l2))
			if(ListFirst(l2) && CompoundArg1(ListFirst(l2))==NewInteger(0))
			{
				FreeAtomic(ListFirst(l2));
				ChangeList(l2,0);
			}
			
		if(chf==0)
		{
	rpt3:		
			for(l2=l1;l2;l2=ListTail(l2))
				if(ListFirst(l2)==0)
				{
					l1=CutFromList(l1,l2);
					goto rpt3;
				}
		}
	}



	if(FAOutput)
	{
		for(l2=l1;l2;l2=ListTail(l2))
		{
			List l3;
			for(l3=CompoundArgN(ListFirst(l2),3);l3;l3=ListTail(l3))
				if(CompoundArg1(ListFirst(l3))==A_GAMMA5)
				{
					Term m2=CopyTerm(ListFirst(l2));
					SetCompoundArg(ListFirst(l3),1,A_GAMMAP);
					for(l3=CompoundArgN(m2,3);l3;l3=ListTail(l3))
					if( CompoundArg1(ListFirst(l3))==A_GAMMA5)
					{
						SetCompoundArg(ListFirst(l3),1,A_GAMMAM);
						break;
					}
					SetCompoundArg(m2,1,
							NewInteger(-IntegerValue(CompoundArg1(m2))));
					l1=AppendLast(l1,m2);
					break;
				}
		}
	}
	
	if(FAOutput && !UFOutput)
	{
		for(l2=l1;l2;l2=ListTail(l2))
		{
			List l3;
			List g=0, gpm=0;
			for(l3=CompoundArgN(ListFirst(l2),3);l3;l3=ListTail(l3))
			{
				if(CompoundArg1(ListFirst(l3))==A_GAMMA)
					g=l3;
				if(CompoundArg1(ListFirst(l3))==A_GAMMAP ||
						CompoundArg1(ListFirst(l3))==A_GAMMAM)
					gpm=l3;
			}
			if(g && !gpm)
				{
					long int maxi=0;
					Integer lfi;
					Term m2;
					for(l3=CompoundArgN(ListFirst(l2),3);l3;l3=ListTail(l3))
					{
						List l4;
						for(l4=CompoundArg2(ListFirst(l3));l4;l4=ListTail(l4))
							if(IntegerValue(ListFirst(l4))>maxi)
								maxi=IntegerValue(ListFirst(l4));
					}					
					lfi=ListNth(CompoundArg2(ListFirst(g)),2);
					ChangeList(ListTail(CompoundArg2(ListFirst(g))),
							NewInteger(maxi+1));
					g=AppendLast(g,MakeCompound2(OPR_SPECIAL,A_GAMMAP,
							MakeList2(NewInteger(maxi+1),lfi)));
					m2=CopyTerm(ListFirst(l2));
					for(l3=CompoundArgN(m2,3);l3;l3=ListTail(l3))
					if( CompoundArg1(ListFirst(l3))==A_GAMMAP)
					{
						SetCompoundArg(ListFirst(l3),1,A_GAMMAM);
						break;
					}					
					l1=AppendLast(l1,m2);
					continue;
				}
		}
	}

	if(FAOutput)
	{
		Integer fi[2];
		int fno=0;
		for(l2=CompoundArg1(a2);l2;l2=ListTail(l2))
			if(CompoundName(CompoundArg2(ListFirst(l2)))==OPR_SPINOR)
				fi[fno++]=ListFirst(CompoundArg1(CompoundArg2(ListFirst(l2))));
		if(fno==2)
		for(l2=l1;l2;l2=ListTail(l2))
		{
			List l3;
			for(l3=CompoundArgN(ListFirst(l2),3);l3;l3=ListTail(l3))
				if(CompoundArg1(ListFirst(l3))==A_DELTA &&
						ListFirst(CompoundArg2(ListFirst(l3)))==fi[0] &&
						ListNth(CompoundArg2(ListFirst(l3)),2)==fi[1])
				{
					Term m2;
					SetCompoundArg(ListFirst(l3),1,A_GAMMAP);
					m2=CopyTerm(ListFirst(l2));
					for(l3=CompoundArgN(m2,3);l3;l3=ListTail(l3))
					if( CompoundArg1(ListFirst(l3))==A_GAMMAP)
					{
						SetCompoundArg(ListFirst(l3),1,A_GAMMAM);
						break;
					}
					l1=AppendLast(l1,m2);
					break;
				}
		}
			
	}

	SetCompoundArg(a2,5,l1);
}


void alg2_kill_gpm(Term a2)
{
	List l1,l2;
	int pmf=0;
	l1=CompoundArgN(a2,5);

	for(l2=l1;l2;l2=ListTail(l2))
	{
		List l3;
		for(l3=CompoundArgN(ListFirst(l2),3);l3;l3=ListTail(l3))
		{
			Term t;
			t=CompoundArg1(ListFirst(l3));
			if(t==A_GAMMAM || t==A_GAMMAP)
			{
				pmf=1;
				break;
			}
		}
		if(pmf)
			break;
	}
	
	if(!pmf)
		return;
	
	l1=ConsumeCompoundArg(a2,5);
	
	if(pmf)
	{
		SetCompoundArg(CompoundArg2(a2),2,NewInteger(2*IntegerValue(CompoundArg2(CompoundArg2(a2)))));
		for(l2=l1;l2;l2=ListTail(l2))
		{
			List l3;
			pmf=0;
			for(l3=CompoundArgN(ListFirst(l2),3);l3;l3=ListTail(l3))
			{
				Term t;
				t=CompoundArg1(ListFirst(l3));
				if(t==A_GAMMAM || t==A_GAMMAP)
					pmf=1;
			}
			if(pmf==0)
				SetCompoundArg(ListFirst(l2),1,NewInteger(2*IntegerValue(CompoundArg1(ListFirst(l2)))));
		}
	}
	
	for(l2=l1;l2;l2=ListTail(l2))
	{
		List l3;
		for(l3=CompoundArgN(ListFirst(l2),3);l3;l3=ListTail(l3))
		{
			Term t;
			t=CompoundArg1(ListFirst(l3));
			if(t==A_GAMMAM || t==A_GAMMAP)
			{
				Term lab1, lab2, cp;
				List l33;

				if(t==A_GAMMAM)
					pmf=-1;
				else
					pmf=1;
				
				cp=CopyTerm(ListFirst(l2));
				
				SetCompoundArg(ListFirst(l3),1,A_DELTA);
				lab1=ListFirst(CompoundArg2(ListFirst(l3)));
				lab2=ListFirst(ListTail(CompoundArg2(ListFirst(l3))));
				
				for(l33=CompoundArgN(ListFirst(l2),3);l33;l33=ListTail(l33))
				{
					if(l33==l3)
						continue;
					if(CompoundArg1(ListFirst(l33))==A_GAMMA && 
							ListNth(CompoundArg2(ListFirst(l33)),2)==lab1)
					{
						ChangeList(ListTail(CompoundArg2(ListFirst(l33))),lab2);
						l33=ConsumeCompoundArg(ListFirst(l2),3);
						l33=CutFromList(l33,l3);
						SetCompoundArg(ListFirst(l2),3,l33);
						break;
					}
				}
				
				for(l3=CompoundArgN(cp,3);l3;l3=ListTail(l3))
				{
					Term t;
					t=CompoundArg1(ListFirst(l3));
					if(t==A_GAMMAM || t==A_GAMMAP)
						break;
				}
				if(l3==0)
				{
					puts("Internal error a2kgpmsf");
					return;
				}
				if(pmf==-1)
					SetCompoundArg(cp,1,NewInteger(-IntegerValue(CompoundArg1(cp))));
				SetCompoundArg(ListFirst(l3),1,A_GAMMA5);
				l1=AppendFirst(l1,cp);
				break;
			}
		}
	}
	
	{
		int chf=1;
		List l3;
		for(l2=l1;l2;l2=ListTail(l2))
		for(l3=ListTail(l2);l3;l3=ListTail(l3))
		{
			if(ListFirst(l2)==0 || ListFirst(l3)==0)
				continue;
			if(EqualTerms(CompoundArg2(ListFirst(l2)),
					CompoundArg2(ListFirst(l3))) &&
				EqualTerms(CompoundArgN(ListFirst(l2),3),
					CompoundArgN(ListFirst(l3),3)))
			{
				long int i1, i2;
				i1=IntegerValue(CompoundArg1(ListFirst(l2)));
				i2=IntegerValue(CompoundArg1(ListFirst(l3)));
				SetCompoundArg(ListFirst(l2),1,NewInteger(i1+i2));
				FreeAtomic(ListFirst(l3));
				ChangeList(l3,0);
				chf=0;
			}
		}
		
		if(chf==0)
		for(l2=l1;l2;l2=ListTail(l2))
			if(ListFirst(l2) && CompoundArg1(ListFirst(l2))==NewInteger(0))
			{
				FreeAtomic(ListFirst(l2));
				ChangeList(l2,0);
			}
			
		if(chf==0)
		{
	rpt3:		
			for(l2=l1;l2;l2=ListTail(l2))
				if(ListFirst(l2)==0)
				{
					l1=CutFromList(l1,l2);
					goto rpt3;
				}
		}
	}
	
	
/*	for(l1;l1;l1=ListTail(l1))
	{
		WriteTerm(ListFirst(l1));
		puts("");
	}
*/	
	SetCompoundArg(a2,5,l1);
	return;
	
}
