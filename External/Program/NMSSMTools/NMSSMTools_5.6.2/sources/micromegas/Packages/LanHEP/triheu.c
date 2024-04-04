#include <string.h>
#include <stdio.h>
#include "lanhep.h"

extern int tri_dbg_mode, tri_third_angle, opTriHeu;
extern List tri_si_co_list;
		
static Atom sa1, ca1, sa2, ca2, sa3, ca3, s2a1, s2a2, s2a3, c2a1, c2a2, c2a3;

static Atom sapb1, capb1, samb1, camb1, sapb2, capb2, samb2, camb2=0; 
static int  sapb_s, capb_s, samb_s, camb_s, c2a1_s, c2a2_s, c2a3_s;

List tri_find_sub(List);

void tri_set_rule(List, List, List);

static int prtcmp(Term p1, Term p2)
	{
	Term pp1,pp2;
	pp1=CompoundArg1(p1);
	pp2=CompoundArg1(p2);
	if(pp1==A_I)
		return -1;
	if(pp2==A_I)
		return 1;
	if(is_atom(pp1) && is_atom(pp2))
		return strcmp(AtomValue(pp1),AtomValue(pp2));
	if(is_integer(pp1) && is_integer(pp2))
		return (int)IntegerValue(pp1)-(int)IntegerValue(pp2);
	if(is_integer(pp1))
		return 1;
	return -1;
	}
	
List tri_rcos_mlt(List l, Atom par, int pw)
	{
	List l1;
	
	for(l1=l;l1;l1=ListTail(l1))
		{
		Term pp;
		pp=ListFirst(l1);
		if(CompoundArg1(pp)==par)
			{
			int n;
			n=pw+(int)IntegerValue(CompoundArg2(pp));
			SetCompoundArg(pp,2,NewInteger(n));
			return l;
			}
		}
	l=AppendLast(l,MakeCompound2(OPR_POW,par,NewInteger(pw)));
	l=SortedList(l,prtcmp);
	return l;
	} 		


List tri_rcos_add(List l, Term m2)
	{
	List ml=l;
	
/*	WriteTerm(l); printf("+"); WriteTerm(m2); puts("");*/
	while(!is_empty_list(l))
		{
		Term t;
		t=ListFirst(l);
		if(EqualTerms(CompoundArg2(t),CompoundArg2(m2)))
				{
				int n1,n2,n;
				n1=(int)IntegerValue(CompoundArg1(t));
				n2=(int)IntegerValue(CompoundArg1(m2));
				n=n1+n2;
				if(n==0)
					{
					ml=CutFromList(ml,l);
					}
				else
					{
					SetCompoundArg(t,1,NewInteger(n));
					}
				FreeAtomic(m2);
				/*printf("= "); WriteTerm(ml); puts(""); getchar();*/
				return ml;
				}
		l=ListTail(l);
		}
	return AppendLast(ml,m2);
	}
	
List tri_cvt(List ml, Atom from, Atom to)
{
	List l1,l2;

beg:

	for(l1=ml;l1;l1=ListTail(l1))
	{
		Term m1,m2;
		List lp, lp2;
		int num;
		
		m1=ListFirst(l1);
		num=(int)IntegerValue(CompoundArg1(m1));
		lp=CompoundArg2(m1);
		
		for(l2=lp;l2;l2=ListTail(l2))
		{
			Term prm;
			int  prmpw;
			
			prm=CompoundArg1(ListFirst(l2));
			prmpw=(int)IntegerValue(CompoundArg2(ListFirst(l2)));
			
			if(prm==from && prmpw==2)
			{
				lp=ConsumeCompoundArg(m1,2);
				ChangeList(l1,0);
				ml=CutFromList(ml,l1);
				
				lp=CutFromList(lp,l2);
				m2=CopyTerm(m1);
				lp2=CopyTerm(lp);
				
				SetCompoundArg(m1,2,lp);
				ml=tri_rcos_add(ml,m1);
				
				lp2=tri_rcos_mlt(lp2,to,2);
				SetCompoundArg(m2,1,NewInteger(-num));
				SetCompoundArg(m2,2,lp2);
				ml=tri_rcos_add(ml,m2);
				goto beg;
			}
			
			if(prm==from && prmpw>2)
			{
				lp=ConsumeCompoundArg(m1,2);
				ChangeList(l1,0);
				ml=CutFromList(ml,l1);
				
				lp=CutFromList(lp,l2);
				m2=CopyTerm(m1);
				lp2=CopyTerm(lp);
				
				lp=tri_rcos_mlt(lp,from,prmpw-2);
				SetCompoundArg(m1,2,lp);
				ml=tri_rcos_add(ml,m1);
				
				lp2=tri_rcos_mlt(lp2,to,2);
				lp2=tri_rcos_mlt(lp2,from,prmpw-2);
				SetCompoundArg(m2,1,NewInteger(-num));
				SetCompoundArg(m2,2,lp2);
				ml=tri_rcos_add(ml,m2);
				goto beg;
			}
		}
	}
	return ml;
}

static List pm2_2_pm(List ml)
{
	List l1,l2,l3;
	
beg:

	/*WriteTerm(ml); puts("");*/

	for(l1=ml;l1;l1=ListTail(l1))
	{
		Term m1,m2;
		List lp, lp2;
		int num;
		
		m1=ListFirst(l1);
		num=(int)IntegerValue(CompoundArg1(m1));
		lp=CompoundArg2(m1);
		
		for(l2=lp;l2;l2=ListTail(l2))
		{
			Term prm;
			int  prmpw;
			
			prm=CompoundArg1(ListFirst(l2));
			prmpw=(int)IntegerValue(CompoundArg2(ListFirst(l2)));
			
			if(prm==sapb2 && prmpw>=2)
			{
				lp=ConsumeCompoundArg(m1,2);
				ChangeList(l1,0);
				ml=CutFromList(ml,l1);
				
				lp=CutFromList(lp,l2);
				m2=CopyTerm(m1);
				lp2=CopyTerm(lp);
				
				if(prmpw>2)
					lp2=tri_rcos_mlt(lp2,sapb2,prmpw-2);
				SetCompoundArg(m2,1,NewInteger(num/2));				
				SetCompoundArg(m2,2,lp2);
				ml=tri_rcos_add(ml,m2);
				
				if(prmpw>2)
					lp=tri_rcos_mlt(lp,sapb2,prmpw-2);
				lp=tri_rcos_mlt(lp,capb1,1);
				SetCompoundArg(m1,1,NewInteger(-num/2*capb_s));
				SetCompoundArg(m1,2,lp);
				ml=tri_rcos_add(ml,m1);

				goto beg;
			}
			
			if(prm==samb2 && prmpw>=2)
			{
				lp=ConsumeCompoundArg(m1,2);
				ChangeList(l1,0);
				ml=CutFromList(ml,l1);
				
				lp=CutFromList(lp,l2);
				m2=CopyTerm(m1);
				lp2=CopyTerm(lp);
				
				if(prmpw>2)
					lp2=tri_rcos_mlt(lp2,samb2,prmpw-2);
				SetCompoundArg(m2,1,NewInteger(num/2));				
				SetCompoundArg(m2,2,lp2);
				ml=tri_rcos_add(ml,m2);
				
				if(prmpw>2)
					lp=tri_rcos_mlt(lp,samb2,prmpw-2);
				lp=tri_rcos_mlt(lp,camb1,1);
				SetCompoundArg(m1,1,NewInteger(-num/2*camb_s));
				SetCompoundArg(m1,2,lp);
				ml=tri_rcos_add(ml,m1);

				goto beg;
			}
		}
		
	}

beg1:

		
	for(l1=ml;l1;l1=ListTail(l1))
	{
		Term m1;
		List lp;
		int num;
		
		m1=ListFirst(l1);
		num=(int)IntegerValue(CompoundArg1(m1));
		lp=CompoundArg2(m1);
		
		for(l2=lp;l2 && ListTail(l2);l2=ListTail(l2))
		for(l3=ListTail(l2);l3;l3=ListTail(l3))
		{
			Term prm1,prm2;
			
			prm1=CompoundArg1(ListFirst(l2));
			prm2=CompoundArg1(ListFirst(l3));
			
			if((prm1==sapb2 && prm2==capb2) ||(prm1==capb2 && prm2==sapb2)) 
			{
				lp=ConsumeCompoundArg(m1,2);
				ChangeList(l1,0);
				ml=CutFromList(ml,l1);
				
				lp=CutFromList(lp,l2);
				lp=CutFromList(lp,l3);

				lp=tri_rcos_mlt(lp,sapb1,1);
				SetCompoundArg(m1,1,NewInteger(num/2*sapb_s));
				SetCompoundArg(m1,2,lp);
				ml=tri_rcos_add(ml,m1);

				goto beg1;
			}
			
			if((prm1==samb2 && prm2==camb2) ||(prm1==camb2 && prm2==samb2)) 
			{
				lp=ConsumeCompoundArg(m1,2);
				ChangeList(l1,0);
				ml=CutFromList(ml,l1);
				
				lp=CutFromList(lp,l2);
				lp=CutFromList(lp,l3);

				lp=tri_rcos_mlt(lp,samb1,1);
				SetCompoundArg(m1,1,NewInteger(num/2*samb_s));
				SetCompoundArg(m1,2,lp);
				ml=tri_rcos_add(ml,m1);

				goto beg1;
			}
		}
	}
	
beg2:

	for(l1=ml;l1;l1=ListTail(l1))
	{
		Term m1,m2;
		List lp, lp2;
		int num;
		
		m1=ListFirst(l1);
		num=(int)IntegerValue(CompoundArg1(m1));
		lp=CompoundArg2(m1);
		
		for(l2=lp;l2 && ListTail(l2);l2=ListTail(l2))
		for(l3=ListTail(l2);l3;l3=ListTail(l3))
		{
			Term prm1,prm2;
			
			prm1=CompoundArg1(ListFirst(l2));
			prm2=CompoundArg1(ListFirst(l3));
			
			if((prm1==capb2 && prm2==camb2) || (prm1==camb2 && prm2==capb2) ||
			   (prm1==sapb2 && prm2==samb2) || (prm1==samb2 && prm2==sapb2)) 
			{
				lp=ConsumeCompoundArg(m1,2);
				ChangeList(l1,0);
				ml=CutFromList(ml,l1);
				
				lp=CutFromList(lp,l2);
				lp=CutFromList(lp,l3);
				lp2=CopyTerm(lp);
				m2=CopyTerm(m1);
				
				lp2=tri_rcos_mlt(lp2,ca1,1);
				if(prm1==capb2 || prm1==camb2)
					SetCompoundArg(m2,1,NewInteger(num/2));
				else
					SetCompoundArg(m2,1,NewInteger(-num/2));
				SetCompoundArg(m2,2,lp2);
				ml=tri_rcos_add(ml,m2);

				lp=tri_rcos_mlt(lp,ca2,1);
				SetCompoundArg(m1,1,NewInteger(num/2));
				SetCompoundArg(m1,2,lp);
				ml=tri_rcos_add(ml,m1);

				goto beg2;
			}
			
			if((prm1==sapb2 && prm2==camb2) || (prm1==camb2 && prm2==sapb2) ||
			   (prm1==samb2 && prm2==capb2) || (prm1==capb2 && prm2==samb2))
			{
				lp=ConsumeCompoundArg(m1,2);
				ChangeList(l1,0);
				ml=CutFromList(ml,l1);
				
				lp=CutFromList(lp,l2);
				lp=CutFromList(lp,l3);
				lp2=CopyTerm(lp);
				m2=CopyTerm(m1);
				
				lp2=tri_rcos_mlt(lp2,sa1,1);
				SetCompoundArg(m2,1,NewInteger(num/2));
				SetCompoundArg(m2,2,lp2);
				ml=tri_rcos_add(ml,m2);

				lp=tri_rcos_mlt(lp,sa2,1);
				if(prm1==sapb2 || prm1==camb2)
					SetCompoundArg(m1,1,NewInteger(num/2));
				else
					SetCompoundArg(m1,1,NewInteger(-num/2));
				SetCompoundArg(m1,2,lp);
				ml=tri_rcos_add(ml,m1);

				goto beg2;
			}
		}
	}
	
	ml=tri_cvt(ml,capb1,sapb1);
	ml=tri_cvt(ml,camb1,samb1);
	return ml;
}
			

static List ab_2_pm2(List ml)
{
	List l1,l2;
	beg:


	/*WriteTerm(ml); puts("");*/

	for(l1=ml;l1;l1=ListTail(l1))
	{
		Term m1,m2;
		List lp, lp2;
		int num;
		
		m1=ListFirst(l1);
		num=(int)IntegerValue(CompoundArg1(m1));
		lp=CompoundArg2(m1);
		
		for(l2=lp;l2;l2=ListTail(l2))
		{
			Term prm;
			int  prmpw;
			
			prm=CompoundArg1(ListFirst(l2));
			prmpw=(int)IntegerValue(CompoundArg2(ListFirst(l2)));
			
			if(prm==ca1 && prmpw==1)
			{
				lp=ConsumeCompoundArg(m1,2);
				ChangeList(l1,0);
				ml=CutFromList(ml,l1);
				
				lp=CutFromList(lp,l2);
				m2=CopyTerm(m1);
				lp2=CopyTerm(lp);
								
				lp2=tri_rcos_mlt(lp2,capb2,1);
				lp2=tri_rcos_mlt(lp2,camb2,1);
				SetCompoundArg(m2,2,lp2);
				ml=tri_rcos_add(ml,m2);
				
				lp=tri_rcos_mlt(lp,sapb2,1);
				lp=tri_rcos_mlt(lp,samb2,1);
				SetCompoundArg(m1,1,NewInteger(-num));
				SetCompoundArg(m1,2,lp);
				ml=tri_rcos_add(ml,m1);

				goto beg;
			}
			
			if(prm==sa1 && prmpw==1)
			{
				lp=ConsumeCompoundArg(m1,2);
				ChangeList(l1,0);
				ml=CutFromList(ml,l1);
				
				lp=CutFromList(lp,l2);
				m2=CopyTerm(m1);
				lp2=CopyTerm(lp);
								
				lp2=tri_rcos_mlt(lp2,sapb2,1);
				lp2=tri_rcos_mlt(lp2,camb2,1);
				SetCompoundArg(m2,2,lp2);
				ml=tri_rcos_add(ml,m2);
				
				lp=tri_rcos_mlt(lp,capb2,1);
				lp=tri_rcos_mlt(lp,samb2,1);
				SetCompoundArg(m1,2,lp);
				ml=tri_rcos_add(ml,m1);

				goto beg;
			}
			
			if(prm==sa1 && prmpw>1)
			{
				lp=ConsumeCompoundArg(m1,2);
				ChangeList(l1,0);
				ml=CutFromList(ml,l1);
				
				lp=CutFromList(lp,l2);
				m2=CopyTerm(m1);
				lp2=CopyTerm(lp);
				
				lp2=tri_rcos_mlt(lp2,sa1,prmpw-1);				
				lp2=tri_rcos_mlt(lp2,sapb2,1);
				lp2=tri_rcos_mlt(lp2,camb2,1);
				SetCompoundArg(m2,2,lp2);
				ml=tri_rcos_add(ml,m2);
				
				lp=tri_rcos_mlt(lp,sa1,prmpw-1);
				lp=tri_rcos_mlt(lp,capb2,1);
				lp=tri_rcos_mlt(lp,samb2,1);
				SetCompoundArg(m1,2,lp);
				ml=tri_rcos_add(ml,m1);

				goto beg;
			}
			
			if(prm==ca2 && prmpw==1)
			{
				lp=ConsumeCompoundArg(m1,2);
				ChangeList(l1,0);
				ml=CutFromList(ml,l1);
				
				lp=CutFromList(lp,l2);
				m2=CopyTerm(m1);
				lp2=CopyTerm(lp);
								
				lp2=tri_rcos_mlt(lp2,capb2,1);
				lp2=tri_rcos_mlt(lp2,camb2,1);
				SetCompoundArg(m2,2,lp2);
				ml=tri_rcos_add(ml,m2);
				
				lp=tri_rcos_mlt(lp,sapb2,1);
				lp=tri_rcos_mlt(lp,samb2,1);
				SetCompoundArg(m1,2,lp);
				ml=tri_rcos_add(ml,m1);

				goto beg;
			}
			
			if(prm==sa2 && prmpw==1)
			{
				lp=ConsumeCompoundArg(m1,2);
				ChangeList(l1,0);
				ml=CutFromList(ml,l1);
				
				lp=CutFromList(lp,l2);
				m2=CopyTerm(m1);
				lp2=CopyTerm(lp);
								
				lp2=tri_rcos_mlt(lp2,sapb2,1);
				lp2=tri_rcos_mlt(lp2,camb2,1);
				SetCompoundArg(m2,2,lp2);
				ml=tri_rcos_add(ml,m2);
				
				lp=tri_rcos_mlt(lp,capb2,1);
				lp=tri_rcos_mlt(lp,samb2,1);
				SetCompoundArg(m1,1,NewInteger(-num));
				SetCompoundArg(m1,2,lp);
				ml=tri_rcos_add(ml,m1);

				goto beg;
			}
			
			if(prm==sa2 && prmpw>1)
			{
				lp=ConsumeCompoundArg(m1,2);
				ChangeList(l1,0);
				ml=CutFromList(ml,l1);
				
				lp=CutFromList(lp,l2);
				m2=CopyTerm(m1);
				lp2=CopyTerm(lp);
				
				lp2=tri_rcos_mlt(lp2,sa2,prmpw-1);				
				lp2=tri_rcos_mlt(lp2,sapb2,1);
				lp2=tri_rcos_mlt(lp2,camb2,1);
				SetCompoundArg(m2,2,lp2);
				ml=tri_rcos_add(ml,m2);
				
				lp=tri_rcos_mlt(lp,sa2,prmpw-1);
				lp=tri_rcos_mlt(lp,capb2,1);
				lp=tri_rcos_mlt(lp,samb2,1);
				SetCompoundArg(m1,1,NewInteger(-num));
				SetCompoundArg(m1,2,lp);
				ml=tri_rcos_add(ml,m1);

				goto beg;
			}
			
		}
	}
	

	ml=tri_cvt(ml,capb2,sapb2);
	ml=tri_cvt(ml,camb2,samb2);

	
	return ml;
} 

static List p_2_b(List ml)
	{
	List l1,l2;

beg:

	for(l1=ml;l1;l1=ListTail(l1))
	{
		Term m1,m2;
		List lp, lp2;
		int num;
		
		m1=ListFirst(l1);
		num=(int)IntegerValue(CompoundArg1(m1));
		lp=CompoundArg2(m1);
		
		for(l2=lp;l2;l2=ListTail(l2))
		{
			Term prm;
			int  prmpw;
			
			prm=CompoundArg1(ListFirst(l2));
			prmpw=(int)IntegerValue(CompoundArg2(ListFirst(l2)));
					
			if(prm==NewInteger(2) && prmpw==1)
			{
				lp=ConsumeCompoundArg(m1,2);
				ChangeList(l1,0);
				ml=CutFromList(ml,l1);
				
				lp=CutFromList(lp,l2);
				m2=CopyTerm(m1);
				lp2=CopyTerm(lp);
								
				lp2=tri_rcos_mlt(lp2,ca2,1);
				lp2=tri_rcos_mlt(lp2,ca1,1);
				SetCompoundArg(m2,2,lp2);
				ml=tri_rcos_add(ml,m2);
				
				lp=tri_rcos_mlt(lp,sa2,1);
				lp=tri_rcos_mlt(lp,sa1,1);
				SetCompoundArg(m1,1,NewInteger(-num));
				SetCompoundArg(m1,2,lp);
				ml=tri_rcos_add(ml,m1);

				goto beg;
			}
			
			if(prm==NewInteger(1) && prmpw==1)
			{
				lp=ConsumeCompoundArg(m1,2);
				ChangeList(l1,0);
				ml=CutFromList(ml,l1);
				
				lp=CutFromList(lp,l2);
				m2=CopyTerm(m1);
				lp2=CopyTerm(lp);
								
				lp2=tri_rcos_mlt(lp2,sa2,1);
				lp2=tri_rcos_mlt(lp2,ca1,1);
				SetCompoundArg(m2,2,lp2);
				ml=tri_rcos_add(ml,m2);
				
				lp=tri_rcos_mlt(lp,ca2,1);
				lp=tri_rcos_mlt(lp,sa1,1);
				SetCompoundArg(m1,2,lp);
				ml=tri_rcos_add(ml,m1);

				goto beg;
			}
			
			if(prm==NewInteger(1) && prmpw>1)
			{
				lp=ConsumeCompoundArg(m1,2);
				ChangeList(l1,0);
				ml=CutFromList(ml,l1);
				
				lp=CutFromList(lp,l2);
				m2=CopyTerm(m1);
				lp2=CopyTerm(lp);
				
				lp2=tri_rcos_mlt(lp2,NewInteger(1),prmpw-1);				
				lp2=tri_rcos_mlt(lp2,sa2,1);
				lp2=tri_rcos_mlt(lp2,ca1,1);
				SetCompoundArg(m2,2,lp2);
				ml=tri_rcos_add(ml,m2);
				
				lp=tri_rcos_mlt(lp,NewInteger(1),prmpw-1);
				lp=tri_rcos_mlt(lp,ca2,1);
				lp=tri_rcos_mlt(lp,sa1,1);
				SetCompoundArg(m1,2,lp);
				ml=tri_rcos_add(ml,m1);

				goto beg;
			}
			
			
		}
	}
	
	ml=tri_cvt(ml,ca1,sa1);
	ml=tri_cvt(ml,ca2,sa2);
	
	return ml;
	
}


	
static List b_2_p(List ml)
	{
	List l1,l2;

beg:

	for(l1=ml;l1;l1=ListTail(l1))
	{
		Term m1,m2;
		List lp, lp2;
		int num;
		
		m1=ListFirst(l1);
		num=(int)IntegerValue(CompoundArg1(m1));
		lp=CompoundArg2(m1);
		
		for(l2=lp;l2;l2=ListTail(l2))
		{
			Term prm;
			int  prmpw;
			
			prm=CompoundArg1(ListFirst(l2));
			prmpw=(int)IntegerValue(CompoundArg2(ListFirst(l2)));
					
			if(prm==ca2 && prmpw==1)
			{
				lp=ConsumeCompoundArg(m1,2);
				ChangeList(l1,0);
				ml=CutFromList(ml,l1);
				
				lp=CutFromList(lp,l2);
				m2=CopyTerm(m1);
				lp2=CopyTerm(lp);
								
				lp2=tri_rcos_mlt(lp2,NewInteger(2),1);
				lp2=tri_rcos_mlt(lp2,ca1,1);
				SetCompoundArg(m2,2,lp2);
				ml=tri_rcos_add(ml,m2);
				
				lp=tri_rcos_mlt(lp,NewInteger(1),1);
				lp=tri_rcos_mlt(lp,sa1,1);
				SetCompoundArg(m1,2,lp);
				ml=tri_rcos_add(ml,m1);

				goto beg;
			}
			
			if(prm==sa2 && prmpw==1)
			{
				lp=ConsumeCompoundArg(m1,2);
				ChangeList(l1,0);
				ml=CutFromList(ml,l1);
				
				lp=CutFromList(lp,l2);
				m2=CopyTerm(m1);
				lp2=CopyTerm(lp);
								
				lp2=tri_rcos_mlt(lp2,NewInteger(1),1);
				lp2=tri_rcos_mlt(lp2,ca1,1);
				SetCompoundArg(m2,2,lp2);
				ml=tri_rcos_add(ml,m2);
				
				lp=tri_rcos_mlt(lp,NewInteger(2),1);
				lp=tri_rcos_mlt(lp,sa1,1);
				SetCompoundArg(m1,1,NewInteger(-num));
				SetCompoundArg(m1,2,lp);
				ml=tri_rcos_add(ml,m1);

				goto beg;
			}
			
			if(prm==sa2 && prmpw>1)
			{
				lp=ConsumeCompoundArg(m1,2);
				ChangeList(l1,0);
				ml=CutFromList(ml,l1);
				
				lp=CutFromList(lp,l2);
				m2=CopyTerm(m1);
				lp2=CopyTerm(lp);
				
				lp2=tri_rcos_mlt(lp2,sa2,prmpw-1);				
				lp2=tri_rcos_mlt(lp2,NewInteger(1),1);
				lp2=tri_rcos_mlt(lp2,ca1,1);
				SetCompoundArg(m2,2,lp2);
				ml=tri_rcos_add(ml,m2);
				
				lp=tri_rcos_mlt(lp,sa2,prmpw-1);
				lp=tri_rcos_mlt(lp,NewInteger(2),1);
				lp=tri_rcos_mlt(lp,sa1,1);
				SetCompoundArg(m1,1,NewInteger(-num));
				SetCompoundArg(m1,2,lp);
				ml=tri_rcos_add(ml,m1);

				goto beg;
			}
			
			
		}
	}
	
	ml=tri_cvt(ml,ca1,sa1);
	ml=tri_cvt(ml,NewInteger(2),NewInteger(1));
	
	return ml;
	
}


static List b_2_m(List ml)
{
	List l1,l2;
	for(l1=ml;l1;l1=ListTail(l1))
		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
			if(CompoundArg1(ListFirst(l2))==sa2 &&
					(IntegerValue(CompoundArg2(ListFirst(l2)))%2)==1)
				SetCompoundArg(ListFirst(l1),1,
						NewInteger(-IntegerValue(CompoundArg1(ListFirst(l1)))));
	return b_2_p(ml);
}


static List m_2_b(List ml)
{
	List l1,l2;
	
	ml=p_2_b(ml);
	
	for(l1=ml;l1;l1=ListTail(l1))
		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
			if(CompoundArg1(ListFirst(l2))==sa2 &&
					(IntegerValue(CompoundArg2(ListFirst(l2)))%2)==1)
				SetCompoundArg(ListFirst(l1),1,
						NewInteger(-IntegerValue(CompoundArg1(ListFirst(l1)))));
	return ml;
}


static int tri_x_prm(List ml, Atomic prm)
{
	List l1,l2;
	int res=-1;
	
	for(l1=ml;l1;l1=ListTail(l1))
	{
		int mres=0;
		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
			if(CompoundArg1(ListFirst(l2))==prm)
			{
				mres=(int)IntegerValue(CompoundArg2(ListFirst(l2)));
				break;
			}
		if(res==-1 || res>mres)
			res=mres;
	}
	
	if(res==0)
		return 0;
	
	for(l1=ml;l1;l1=ListTail(l1))
	{
		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
			if(CompoundArg1(ListFirst(l2))==prm)
				break;
		if(l2)
		{
			int p;
			p=(int)IntegerValue(CompoundArg2(ListFirst(l2)))-res;
			if(p)
				SetCompoundArg(ListFirst(l2),2,NewInteger(p));
			else
			{
				List l3;
				l3=ConsumeCompoundArg(ListFirst(l1),2);
				l3=CutFromList(l3,l2);
				SetCompoundArg(ListFirst(l1),2,l3);
			}
		}
	}
	
	return res;
}


static void tri_ini_a(Term);
static void tri_c_2_e(List);

void tri_wrt_sln(List, FILE *);
void tri_red_numf(List l1, List l2, List l3, List l4);
Term tri_m2_l_to_e(List);

void tri_wrt_sln1(List l, FILE *f)
{
	List l1,l2;
	Term e;
	l1=CopyTerm(l);
	l2=CopyTerm(l);
	tri_red_numf(l1,l2,0,0);
	e=tri_m2_l_to_e(l1);
	fWriteTerm(f,e);
	FreeAtomic(e);
	FreeAtomic(l1);
	FreeAtomic(l2);
}

int tri_opti_2(Term, List *);
int tri_x_c2(List *, Atom);

int tri_heu_2(Term v2d, List ml, List *l_l, List *l_r)
{
	List l1, l2, csf, csf1;
	int p;
	int csf_s=1;
	
	tri_ini_a(v2d);

	l1=CopyTerm(ml);
	
	csf=NewList();
	
	l1=b_2_p(l1);
	p=tri_x_prm(l1,NewInteger(1));
	if(p)
	{
		csf=AppendLast(csf,MakeCompound2(OPR_POW,
					CompoundArg1(CompoundArgN(v2d,4)),NewInteger(p)));
		if(sapb_s==-1 && (p%2)==1)
			csf_s=-csf_s;
	}

	l1=tri_cvt(l1,NewInteger(1),NewInteger(2));
	p=tri_x_prm(l1,NewInteger(2));
	if(p)
	{
		csf=AppendLast(csf,MakeCompound2(OPR_POW,
					CompoundArg1(CompoundArgN(v2d,6)),NewInteger(p)));
		if(capb_s==-1 && (p%2)==1)
			csf_s=-csf_s;
	}
	l1=tri_cvt(l1,NewInteger(2),NewInteger(1));
	l1=p_2_b(l1);
	
	l1=b_2_m(l1);
	p=tri_x_prm(l1,NewInteger(1));
	if(p)
	{
		csf=AppendLast(csf,MakeCompound2(OPR_POW,
					CompoundArg1(CompoundArgN(v2d,5)),NewInteger(p)));
		if(samb_s==-1 && (p%2)==1)
			csf_s=-csf_s;
	}
	
	l1=tri_cvt(l1,NewInteger(1),NewInteger(2));
	p=tri_x_prm(l1,NewInteger(2));
	if(p)
	{
		csf=AppendLast(csf,MakeCompound2(OPR_POW,
					CompoundArg1(CompoundArgN(v2d,7)),NewInteger(p)));
		if(camb_s==-1 && (p%2)==1)
			csf_s=-csf_s;
	}
	l1=tri_cvt(l1,NewInteger(2),NewInteger(1));
	l1=m_2_b(l1);

	if(c2a1)
	{
		p=tri_x_c2(&l1,sa1);
		if(p)
		{
			csf=AppendLast(csf,MakeCompound2(OPR_POW,c2a1,NewInteger(p)));
			if(c2a1_s==-1 && (p%2)==1)
				csf_s=-csf_s;
		}
	}

	if(c2a2)
	{
		p=tri_x_c2(&l1,sa2);
		if(p)
		{
			csf=AppendLast(csf,MakeCompound2(OPR_POW,c2a2,NewInteger(p)));
			if(c2a2_s==-1 && (p%2)==1)
				csf_s=-csf_s;
		}
	}
		
	
	if(ListLength(l1)==1 && CompoundArg2(ListFirst(l1))==0)
	{
		Term m2;
		int no;
		no=(int)IntegerValue(CompoundArg1(ListFirst(l1)));
		FreeAtomic(l1);
		m2=MakeCompound(A_MTERM,3);
		SetCompoundArg(m2,1,NewInteger(no*csf_s));
		SetCompoundArg(m2,2,csf);
		l1=AppendLast(NewList(),m2);
		l2=CopyTerm(l1);
		if(tri_dbg_mode)
			puts("TriHeu2: expression is const after csf extraction");
		tri_c_2_e(l2);
		*l_l=l2;
		*l_r=l1;
		return 1;
	}
	
	l2=tri_find_sub(l1);
	if(l2)
	{
		List l3;
		if(tri_dbg_mode)
			puts("TriHeu2: sub found after csf extraction");
		l3=l1;
		l1=l2;
		l2=l3;
		goto mkres;
	}


	l2=CopyTerm(l1);
	
	if(opTriHeu==4)
	{
		printf("\na/b before opt: ");
		tri_wrt_sln1(l2,stdout);puts("");
	}
	tri_opti_2(v2d, &l2);
	if(opTriHeu==4)
	{
		printf("a/b  after opt: ");
		tri_wrt_sln1(l2,stdout);puts("");
	}
	
	l1=ab_2_pm2(l1);

	l1=pm2_2_pm(l1);
	
	if(opTriHeu==4)
	{
		printf("+/- before opt: ");
		tri_wrt_sln1(l1,stdout);puts("");
	}
	tri_opti_2(v2d, &l1);
	if(opTriHeu==4)
	{
		printf("+/-  after opt: ");
		tri_wrt_sln1(l1,stdout);puts("");
	}

	
	tri_third_angle=0;

	if(ListLength(l1)<ListLength(l2))
		FreeAtomic(l2);
	else
	{
		FreeAtomic(l1);
		l1=l2;
	}
	
	l2=CopyTerm(l1);
	tri_c_2_e(l2);
	
	if(tri_dbg_mode)
		puts("TriHeu2: last resort");
	
/*	{
		List k1,k2;
		k1=CopyTerm(l1);
		k2=CopyTerm(l2);
		tri_set_rule(0,k2,k1);
		FreeAtomic(k1);
		FreeAtomic(k2);
	}
*/
	
mkres:
	if(csf)
	{
		int nu;
		List l3;
		nu=(int)IntegerValue(CompoundArg1(ListFirst(l1)));
		for(l3=ListTail(l1);l3;l3=ListTail(l3))
			nu=(int)gcf(nu,IntegerValue(CompoundArg1(ListFirst(l3))));
		for(l3=l1;l3;l3=ListTail(l3))
			SetCompoundArg(ListFirst(l3),1,
					NewInteger(IntegerValue(CompoundArg1(ListFirst(l3)))/nu));
		for(l3=l2;l3;l3=ListTail(l3))
			SetCompoundArg(ListFirst(l3),1,
					NewInteger(IntegerValue(CompoundArg1(ListFirst(l3)))/nu));
		l3=CopyTerm(csf);
		csf1=MakeCompound(A_MTERM,3);
		SetCompoundArg(csf1,1,NewInteger(nu*csf_s));
		SetCompoundArg(csf1,2,AppendLast(csf,l1));
		*l_r=AppendLast(NewList(),csf1);
		csf=l3;
		csf1=MakeCompound(A_MTERM,3);
		SetCompoundArg(csf1,1,NewInteger(nu*csf_s));
		SetCompoundArg(csf1,2,l3);
		csf1=AppendLast(NewList(),csf1);
		tri_c_2_e(csf1);
		l3=ConsumeCompoundArg(ListFirst(csf1),2);
		l3=AppendLast(l3,l2);
		SetCompoundArg(ListFirst(csf1),2,l3);
		*l_l=csf1;
		return 1;
	}
	
	*l_l=l2;
	*l_r=l1;
	
	return 1;
	
}

static void tri_c_2_e(List ml)
{
	List l1,l2;
	
	for(l1=ml;l1;l1=ListTail(l1))
		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
		{
			Term m1,m2;
			
			if(is_atom(CompoundArg1(ListFirst(l2))) &&
					GetAtomProperty(CompoundArg1(ListFirst(l2)),A_DUMMY_PRM))
			{
				SetCompoundArg(ListFirst(l2),1,CopyTerm(
					GetAtomProperty(CompoundArg1(ListFirst(l2)),A_DUMMY_PRM)));
			}

			if(CompoundArg1(ListFirst(l2))==sapb1)
			{
				m1=MakeCompound(A_MTERM,3);
				SetCompoundArg(m1,1,NewInteger(sapb_s));
				SetCompoundArg(m1,2,MakeList2(
						MakeCompound2(OPR_POW,sa1,NewInteger(1)),
						MakeCompound2(OPR_POW,ca2,NewInteger(1))));
				
				m2=MakeCompound(A_MTERM,3);
				SetCompoundArg(m2,1,NewInteger(sapb_s));
				SetCompoundArg(m2,2,MakeList2(
						MakeCompound2(OPR_POW,ca1,NewInteger(1)),
						MakeCompound2(OPR_POW,sa2,NewInteger(1))));
				SetCompoundArg(ListFirst(l2),1,MakeList2(m1,m2));
			}
			
			if(CompoundArg1(ListFirst(l2))==capb1)
			{
				m1=MakeCompound(A_MTERM,3);
				SetCompoundArg(m1,1,NewInteger(capb_s));
				SetCompoundArg(m1,2,MakeList2(
						MakeCompound2(OPR_POW,ca1,NewInteger(1)),
						MakeCompound2(OPR_POW,ca2,NewInteger(1))));
				
				m2=MakeCompound(A_MTERM,3);
				SetCompoundArg(m2,1,NewInteger(-capb_s));
				SetCompoundArg(m2,2,MakeList2(
						MakeCompound2(OPR_POW,sa1,NewInteger(1)),
						MakeCompound2(OPR_POW,sa2,NewInteger(1))));
				if(capb_s==1)
					SetCompoundArg(ListFirst(l2),1,MakeList2(m1,m2));
				else
					SetCompoundArg(ListFirst(l2),1,MakeList2(m2,m1));
			}
			
			if(CompoundArg1(ListFirst(l2))==samb1)
			{
				m1=MakeCompound(A_MTERM,3);
				SetCompoundArg(m1,1,NewInteger(samb_s));
				SetCompoundArg(m1,2,MakeList2(
						MakeCompound2(OPR_POW,sa1,NewInteger(1)),
						MakeCompound2(OPR_POW,ca2,NewInteger(1))));
				
				m2=MakeCompound(A_MTERM,3);
				SetCompoundArg(m2,1,NewInteger(-samb_s));
				SetCompoundArg(m2,2,MakeList2(
						MakeCompound2(OPR_POW,ca1,NewInteger(1)),
						MakeCompound2(OPR_POW,sa2,NewInteger(1))));
				if(samb_s==1)
					SetCompoundArg(ListFirst(l2),1,MakeList2(m1,m2));
				else
					SetCompoundArg(ListFirst(l2),1,MakeList2(m2,m1));
			}
			
			if(CompoundArg1(ListFirst(l2))==camb1)
			{
				m1=MakeCompound(A_MTERM,3);
				SetCompoundArg(m1,1,NewInteger(camb_s));
				SetCompoundArg(m1,2,MakeList2(
						MakeCompound2(OPR_POW,ca1,NewInteger(1)),
						MakeCompound2(OPR_POW,ca2,NewInteger(1))));
				
				m2=MakeCompound(A_MTERM,3);
				SetCompoundArg(m2,1,NewInteger(camb_s));
				SetCompoundArg(m2,2,MakeList2(
						MakeCompound2(OPR_POW,sa1,NewInteger(1)),
						MakeCompound2(OPR_POW,sa2,NewInteger(1))));
				SetCompoundArg(ListFirst(l2),1,MakeList2(m1,m2));
			}
			
			if(CompoundArg1(ListFirst(l2))==s2a1)
			{
				m1=MakeCompound(A_MTERM,3);
				SetCompoundArg(m1,1,NewInteger(2));
				SetCompoundArg(m1,2,MakeList2(
						MakeCompound2(OPR_POW,sa1,NewInteger(1)),
						MakeCompound2(OPR_POW,ca1,NewInteger(1))));

				SetCompoundArg(ListFirst(l2),1,MakeList1(m1));
			}
			
			if(CompoundArg1(ListFirst(l2))==c2a1)
			{
				m1=MakeCompound(A_MTERM,3);
				SetCompoundArg(m1,1,NewInteger(c2a1_s));
				SetCompoundArg(m1,2,MakeList1(
						MakeCompound2(OPR_POW,ca1,NewInteger(2))));
				
				m2=MakeCompound(A_MTERM,3);
				SetCompoundArg(m2,1,NewInteger(-c2a1_s));
				SetCompoundArg(m2,2,MakeList1(
						MakeCompound2(OPR_POW,sa1,NewInteger(2))));
				if(c2a1_s==1)
					SetCompoundArg(ListFirst(l2),1,MakeList2(m1,m2));
				else
					SetCompoundArg(ListFirst(l2),1,MakeList2(m2,m1));
			}
			
			if(CompoundArg1(ListFirst(l2))==s2a2)
			{
				m1=MakeCompound(A_MTERM,3);
				SetCompoundArg(m1,1,NewInteger(2));
				SetCompoundArg(m1,2,MakeList2(
						MakeCompound2(OPR_POW,sa2,NewInteger(1)),
						MakeCompound2(OPR_POW,ca2,NewInteger(1))));

				SetCompoundArg(ListFirst(l2),1,MakeList1(m1));
			}
			
			if(CompoundArg1(ListFirst(l2))==c2a2)
			{
				m1=MakeCompound(A_MTERM,3);
				SetCompoundArg(m1,1,NewInteger(c2a2_s));
				SetCompoundArg(m1,2,MakeList1(
						MakeCompound2(OPR_POW,ca2,NewInteger(2))));
				
				m2=MakeCompound(A_MTERM,3);
				SetCompoundArg(m2,1,NewInteger(-c2a2_s));
				SetCompoundArg(m2,2,MakeList1(
						MakeCompound2(OPR_POW,sa2,NewInteger(2))));
				if(c2a2_s==1)
					SetCompoundArg(ListFirst(l2),1,MakeList2(m1,m2));
				else
					SetCompoundArg(ListFirst(l2),1,MakeList2(m2,m1));
			}
			
			if(CompoundArg1(ListFirst(l2))==s2a3)
			{
				m1=MakeCompound(A_MTERM,3);
				SetCompoundArg(m1,1,NewInteger(2));
				SetCompoundArg(m1,2,MakeList2(
						MakeCompound2(OPR_POW,sa3,NewInteger(1)),
						MakeCompound2(OPR_POW,ca3,NewInteger(1))));

				SetCompoundArg(ListFirst(l2),1,MakeList1(m1));
			}
			
			if(CompoundArg1(ListFirst(l2))==c2a3)
			{
				m1=MakeCompound(A_MTERM,3);
				SetCompoundArg(m1,1,NewInteger(c2a3_s));
				SetCompoundArg(m1,2,MakeList1(
						MakeCompound2(OPR_POW,ca3,NewInteger(2))));
				
				m2=MakeCompound(A_MTERM,3);
				SetCompoundArg(m2,1,NewInteger(-c2a3_s));
				SetCompoundArg(m2,2,MakeList1(
						MakeCompound2(OPR_POW,sa3,NewInteger(2))));
				if(c2a3_s==1)
					SetCompoundArg(ListFirst(l2),1,MakeList2(m1,m2));
				else
					SetCompoundArg(ListFirst(l2),1,MakeList2(m2,m1));
			}
			
		}

}

static void tri_ini_a(Term v2d)
{
	Term t;
	
	t=ListNth(tri_si_co_list,(int)IntegerValue(CompoundArg1(v2d)));
	
	sa1=CompoundArg1(t);
	ca1=CompoundArg2(t);
	s2a1=CompoundArgN(t,3);
	c2a1=CompoundArgN(t,5);
	if(c2a1)
	{
		c2a1_s=CompoundName(c2a1)==OPR_PLUS ? 1 : -1;
		c2a1=CompoundArg1(c2a1);
	}
	
	t=ListNth(tri_si_co_list,(int)IntegerValue(CompoundArg2(v2d)));
	
	sa2=CompoundArg1(t);
	ca2=CompoundArg2(t);
	s2a2=CompoundArgN(t,3);
	c2a2=CompoundArgN(t,5);
	if(c2a2)
	{
		c2a2_s=CompoundName(c2a2)==OPR_PLUS ? 1 : -1;
		c2a2=CompoundArg1(c2a2);
	}
	
	sapb1=CompoundArg1(CompoundArgN(v2d,4));
	capb1=CompoundArg1(CompoundArgN(v2d,6));
	samb1=CompoundArg1(CompoundArgN(v2d,5));
	camb1=CompoundArg1(CompoundArgN(v2d,7));
	
	sapb_s=CompoundName(CompoundArgN(v2d,4))==OPR_PLUS ? 1 : -1;
	capb_s=CompoundName(CompoundArgN(v2d,6))==OPR_PLUS ? 1 : -1;
	samb_s=CompoundName(CompoundArgN(v2d,5))==OPR_PLUS ? 1 : -1;
	camb_s=CompoundName(CompoundArgN(v2d,7))==OPR_PLUS ? 1 : -1;
	
	
	
	if(tri_third_angle)
	{
		Term t;
		t=ListNth(tri_si_co_list,tri_third_angle);
		sa3=CompoundArg1(t);
		ca3=CompoundArg2(t);
		s2a3=CompoundArgN(t,3);
		c2a3=CompoundArgN(t,5);
		if(c2a3)
		{
			c2a3_s=CompoundName(c2a3)==OPR_PLUS ? 1 : -1;
			c2a3=CompoundArg1(c2a3);
		}
	}
	else
		sa3=ca3=c2a3=c2a3_s=0;
	
	
	if(camb2==0)
	{
		sapb2=NewAtom("s+2",0);
		capb2=NewAtom("c+2",0);
		samb2=NewAtom("s-2",0);
		camb2=NewAtom("c-2",0);
	}
	
}
