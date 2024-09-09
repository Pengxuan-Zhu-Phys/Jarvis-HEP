#include <stdio.h>
#include <string.h>
#include "lanhep.h"

extern List tri_si_co_list, tri_sc_2_list;
extern int tri_third_angle, opAbbrVrt;
List tri_rcos_add(List, Term);

int tri_mcmp(Term, Term);
		
int tri_divs(Term a1, int *ft, List *l_1, List *l_2, List ml)
{
	
	List l,ll;
	
	*ft=0;
	
	for(l=ml;l;l=ListTail(l))
	{
		int c1,c2;
		
		c1=c2=0;
		
/*		if(CompoundArg2(ListFirst(l))==0)
			*ft=IntegerValue(CompoundArg1(ListFirst(l)));
*/		
		for(ll=CompoundArg2(ListFirst(l));ll;ll=ListTail(ll))
		{
			Atom p;
			p=CompoundArg1(ListFirst(ll));
			if(p==CompoundArg1(a1) || p==CompoundArg2(a1))
				c1++;
			else
				c2++;
		}
		
		if(c1 && c2)
		{
			return 0;
		}
		
	}
	
	*l_1=NewList();
	*l_2=NewList();
	
	for(l=ml;l;l=ListTail(l))
	{
		int c1,c2;
		
		c1=c2=0;
/*		
		if(CompoundArg2(ListFirst(l))==0)
			continue;
*/		
		for(ll=CompoundArg2(ListFirst(l));ll;ll=ListTail(ll))
		{
			Atom p;
			p=CompoundArg1(ListFirst(ll));
			if(p==CompoundArg1(a1) || p==CompoundArg2(a1))
				c1++;
			else
				c2++;
		}
		
		if(!c1)
			*l_2=AppendLast(*l_2,ListFirst(l));
		else
			*l_1=AppendLast(*l_1,ListFirst(l));
	}
	
	RemoveList(ml);
	
	return 1;
}

static int tri_x_c2_1(List *l, Atom s)
{
	
	List l1,l2,l3,l4,lf=0;
	int scoeff=0;
	
	l1=CopyTerm(*l);

	
	for(l2=l1;l2;l2=ListTail(l2))
		for(l3=CompoundArg2(ListFirst(l2));l3;l3=ListTail(l3))
			if(CompoundArg1(ListFirst(l3))==s && 
					IntegerValue(CompoundArg2(ListFirst(l3)))>scoeff)
				scoeff=(int)IntegerValue(CompoundArg2(ListFirst(l3)));
	if(scoeff<2)
	{
		FreeAtomic(l1);
		return 0;
	}
	
	while(scoeff>1)
	{

	xyz:
		for(l3=l1;l3;l3=ListTail(l3))
		{
			for(l4=CompoundArg2(ListFirst(l3));l4;l4=ListTail(l4))
				if(CompoundArg1(ListFirst(l4))==s && 
						CompoundArg2(ListFirst(l4))==NewInteger(scoeff))
					break;
			if(l4)
			{
				Term m2;
				int n;
				m2=ListFirst(l3);
				ChangeList(l3,0);
				l1=CutFromList(l1,l3);
				if(scoeff==2)
				{
					l3=ConsumeCompoundArg(m2,2);
					l3=CutFromList(l3,l4);
					SetCompoundArg(m2,2,l3);
				}
				else
					SetCompoundArg(ListFirst(l4),2,
						NewInteger(IntegerValue(CompoundArg2(ListFirst(l4)))-2));
				n=(int)IntegerValue(CompoundArg1(m2));
				if(n%2)
				{
					FreeAtomic(l1);
					FreeAtomic(m2);
					FreeAtomic(lf);
					return 0;
				}
				SetCompoundArg(m2,1,NewInteger(-n/2));
				lf=tri_rcos_add(lf,CopyTerm(m2));
				SetCompoundArg(m2,1,NewInteger(n/2));
				l1=tri_rcos_add(l1,m2);
				goto xyz;
			}
		}
		
		scoeff--;
	}
	
	if(!is_empty_list(l1))
	{
		FreeAtomic(l1);
		FreeAtomic(lf);
		return 0;
	}
	
	
	FreeAtomic(*l);
	*l=lf;
	
	
	return 1;
}
	

int tri_x_c2(List *l, Atom s)
{
	int ret=0;
	if(s==0)
		return 0;
	while(tri_x_c2_1(l,s))
		ret++;
	return ret;
}

int tri_divf(Term a1, List *l_1, List *l_2, List ml)
{
	
	List l1,l2,f_l;
	Atom c,s, ret1, ret2;

	s=CompoundArg1(a1);
	c=CompoundArg2(a1);
	f_l=NewList();
	
	for(l1=ml;l1;l1=ListTail(l1))
	{
		Term m2;
		List ls=0,lc=0,lp;
		int pc, ps;
		
		m2=CopyTerm(ListFirst(l1));
		lp=ConsumeCompoundArg(m2,2);
		
		for(l2=lp;l2;l2=ListTail(l2))
		{
			if(CompoundArg1(ListFirst(l2))==c)
				lc=l2;
			if(CompoundArg1(ListFirst(l2))==s)
				ls=l2;
		}
		
		if(lc==0)
			pc=0;
		else
		{
			pc=1;
			lp=CutFromList(lp,lc);
		}
		
		if(ls==0)
			ps=0;
		else
		{
			ps=(int)IntegerValue(CompoundArg2(ListFirst(ls)));
			lp=CutFromList(lp,ls);
		}
		
		SetCompoundArg(m2,2,lp);
		
		for(l2=f_l;l2;l2=ListTail(l2))
			if( CompoundArg1(ListFirst(l2))==NewInteger(ps) &&
				CompoundArg1(ListFirst(l2))==NewInteger(ps))
				{
					AppendLast(CompoundArgN(ListFirst(l2),3),m2);
					break;
				}
		
		if(is_empty_list(l2))
		{
			Term t;
			t=MakeCompound(OPR_POW,4);
			SetCompoundArg(t,1,NewInteger(ps));
			SetCompoundArg(t,2,NewInteger(pc));
			SetCompoundArg(t,3,AppendLast(NewList(),m2));
			f_l=AppendLast(f_l,t);
		}
		
	}
	
	for(l1=f_l;l1;l1=ListTail(l1))
	{
		int cf;
		List ml;
		
		ml=ConsumeCompoundArg(ListFirst(l1),3);
		
		ml=SortedList(ml,tri_mcmp);
		
		cf=(int)IntegerValue(CompoundArg1(ListFirst(ml)));
		if(cf<0)
			cf=-cf;
		for(l2=ListTail(ml);l2;l2=ListTail(l2))
			cf=(int)gcf(cf,IntegerValue(CompoundArg1(ListFirst(l2))));
		if(IntegerValue(CompoundArg1(ListFirst(ml)))<0)
			cf=-cf;
		for(l2=ml;l2;l2=ListTail(l2))
			SetCompoundArg(ListFirst(l2),1,
					NewInteger(IntegerValue(CompoundArg1(ListFirst(l2)))/cf));
		SetCompoundArg(ListFirst(l1),3,ml);
		SetCompoundArg(ListFirst(l1),4,NewInteger(cf));
	}
	
	
	for(l1=ListTail(f_l);l1;l1=ListTail(l1))
		if(!EqualTerms(CompoundArgN(ListFirst(f_l),3),
				CompoundArgN(ListFirst(l1),3)))
	{
		FreeAtomic(f_l);
		return 0;
	}
	
	
	ret1=ConsumeCompoundArg(ListFirst(f_l),3);
	ret2=NewList();
	
	for(l1=f_l;l1;l1=ListTail(l1))
	{
		Term m2;
		List pl;
		int pc, ps;
		ps=(int)IntegerValue(CompoundArg1(ListFirst(l1)));
		pc=(int)IntegerValue(CompoundArg2(ListFirst(l1)));
		pl=NewList();
		if(pc)
			pl=AppendLast(pl,MakeCompound2(OPR_POW,c,NewInteger(pc)));
		if(ps)
			pl=AppendLast(pl,MakeCompound2(OPR_POW,s,NewInteger(ps)));
		m2=MakeCompound(A_MTERM,3);
		SetCompoundArg(m2,1,NewInteger(
				IntegerValue(CompoundArgN(ListFirst(l1),4))/128));
		SetCompoundArg(m2,2,pl);
		ret2=AppendLast(ret2,m2);
	}
	
	for(l1=ret1;l1;l1=ListTail(l1))
		SetCompoundArg(ListFirst(l1),1,
					NewInteger(IntegerValue(CompoundArg1(ListFirst(l1)))*128));
/*	
	DumpList(ret1);
	DumpList(ret2);
	puts("YES!!!!!!!!!!!!!!!!!!");
	
	getchar();
*/
	
	*l_1=ret1;
	*l_2=ret2;
	FreeAtomic(f_l);
	FreeAtomic(ml);
	
	return 1;
}

List tri_rcos_mlt(List, Atom, int);
Term tri_m2_l_to_e(List);

static void rpl_pow(Term e)
{
	int a,i;
	if(!is_compound(e))
		return;
	a=CompoundArity(e);
	for(i=1;i<=a;i++)
		rpl_pow(CompoundArgN(e,i));
	if(CompoundName(e)==OPR_POW)
		SetCompoundName(e,OPR_CARET);
}

Atom tri_dummy_prm(List ml)
{
	Term e, m2_1,m2_2;
	char buf[380];
	List l;
	Atom np;
	int i;
	
		
	if(ListLength(ml)!=2)
		puts("Internal error (mkdupr)");
	
	e=tri_m2_l_to_e(ml);
	if(!(FAOutput&&opAbbrVrt) &&(ChepVersion>3 || FAOutput))
		rpl_pow(e);
	sWriteTerm(buf+1,e);
	FreeAtomic(e);
	buf[0]='(';
	for(i=0;buf[i];i++);
	buf[i]=')';
	buf[i+1]=0;
	np=NewAtom(buf,0);

	SetAtomProperty(np,PROP_TYPE,OPR_PARAMETER);
	SetAtomProperty(np,A_DUMMY_PRM,ml);
	//WriteTerm(np); printf(" -> "); WriteTerm(ml);puts("");
	buf[0]='(';
	buf[1]=0;
	m2_1=ListFirst(ml);
	m2_2=ListFirst(ListTail(ml));
	if((IntegerValue(CompoundArg1(m2_1))!=1 && 
			IntegerValue(CompoundArg1(m2_1))!=-1) ||
			CompoundArg2(m2_1)==0)
		sprintf(buf+strlen(buf),"%ld",IntegerValue(CompoundArg1(m2_1)));
	if(CompoundArg2(m2_1) && IntegerValue(CompoundArg1(m2_1))==-1)
		sprintf(buf+strlen(buf),"-");
	for(l=CompoundArg2(m2_1);l;l=ListTail(l))
	{
		Atom prm, prmn;
		prm=CompoundArg1(ListFirst(l));
		prmn=GetAtomProperty(prm,A_TEXNAME);
		if(!prmn)
			prmn=prm;
		sprintf(buf+strlen(buf)," %s ",AtomValue(prmn));
		if(IntegerValue(CompoundArg2(ListFirst(l)))>1)
			sprintf(buf+strlen(buf),"{}^%ld",IntegerValue(CompoundArg2(ListFirst(l))));
		if(TEX_set_dot && ListTail(l))
			sprintf(buf+strlen(buf)," \\cdot ");
	}
	
	if((IntegerValue(CompoundArg1(m2_2))!=1 && 
			IntegerValue(CompoundArg1(m2_2))!=-1) ||
			CompoundArg2(m2_2)==0)
		sprintf(buf+strlen(buf),"%+ld",IntegerValue(CompoundArg1(m2_2)));
	if(CompoundArg2(m2_2) && IntegerValue(CompoundArg1(m2_2))==-1)
		sprintf(buf+strlen(buf),"-");
	if(CompoundArg2(m2_2) && IntegerValue(CompoundArg1(m2_2))==1)
		sprintf(buf+strlen(buf),"+");
	
	for(l=CompoundArg2(m2_2);l;l=ListTail(l))
	{
		Atom prm, prmn;
		prm=CompoundArg1(ListFirst(l));
		prmn=GetAtomProperty(prm,A_TEXNAME);
		if(!prmn)
			prmn=prm;
		sprintf(buf+strlen(buf)," %s ",AtomValue(prmn));
		if(IntegerValue(CompoundArg2(ListFirst(l)))>1)
			sprintf(buf+strlen(buf),"{}^%ld",IntegerValue(CompoundArg2(ListFirst(l))));
		if(TEX_set_dot && ListTail(l))
			sprintf(buf+strlen(buf)," \\cdot ");
	}
	
	sprintf(buf+strlen(buf),")");

	SetAtomProperty(np,A_TEXNAME,NewAtom(buf,0));
	SetAtomProperty(np,A_TEXLENGTH,NewInteger(6));
/*puts(buf);*/
	
	return np;
}

static int tri_rep_c(List *ml, Atom s, Atom c, int c2a)
{
	List l1,l2;
	
	if(s==0 || c==0)
		return 0;
		
	for(l1=*ml;l1;l1=ListTail(l1))
	{
		List pp, pps;
		int  nu;
		Term m2;
		
		pps=CopyTerm(CompoundArg2(ListFirst(l1)));
		pp=tri_rcos_mlt(CopyTerm(pps),s,2);
		
		nu=(int)IntegerValue(CompoundArg1(ListFirst(l1)));
		
		for(l2=*ml;l2;l2=ListTail(l2))
		{
			if(l1==l2)
				continue;
			if(EqualTerms(pp,CompoundArg2(ListFirst(l2))) &&
					CompoundArg1(ListFirst(l2))==NewInteger(c2a?-2*nu:-nu))
				break;
		}
		
		FreeAtomic(pp);
		
		if(is_empty_list(l2))
		{
			FreeAtomic(pps);
			continue;
		}
		
		*ml=CutFromList(*ml,l1);
		*ml=CutFromList(*ml,l2);
		
		m2=MakeCompound(A_MTERM,3);
		SetCompoundArg(m2,1,NewInteger(c2a?c2a*nu:nu));
		pps=tri_rcos_mlt(pps,c,c2a?1:2);
		SetCompoundArg(m2,2,pps);
		*ml=AppendLast(*ml,m2);
/*puts(c2a?"+2c":"+c");*/
		return 1;
		
	}
	
	return 0;
}

static int tri_rep_c2(List *ml, Atom s, Atom c, Atom c2)
{
	List l1,l2;
	
	if(s==0 || c==0 || c2==0)
		return 0;
		
	for(l1=*ml;l1;l1=ListTail(l1))
	{
		List pp, pps;
		int  nu;
		Term m2;
		
		pps=CopyTerm(CompoundArg2(ListFirst(l1)));
		pp=tri_rcos_mlt(CopyTerm(pps),s,2);
		pp=tri_rcos_mlt(pp,c,2);
		
		nu=(int)IntegerValue(CompoundArg1(ListFirst(l1)));
		
		for(l2=*ml;l2;l2=ListTail(l2))
		{
			if(l1==l2)
				continue;
			if(EqualTerms(pp,CompoundArg2(ListFirst(l2))) &&
					CompoundArg1(ListFirst(l2))==NewInteger(-4*nu))
				break;
		}
		
		FreeAtomic(pp);
		
		if(is_empty_list(l2))
		{
			FreeAtomic(pps);
			continue;
		}
		
		*ml=CutFromList(*ml,l1);
		*ml=CutFromList(*ml,l2);
		
		m2=MakeCompound(A_MTERM,3);
		SetCompoundArg(m2,1,NewInteger(nu));
		pps=tri_rcos_mlt(pps,c2,2);
		SetCompoundArg(m2,2,pps);
		*ml=AppendLast(*ml,m2);

		return 1;
		
	}
	
	return 0;
}

static int tri_rep_cc(List *ml, Atom s1, Atom s2, Atom c1, Atom c2)
{
	List l1,l2,l3;
	
	if(s1==0 || c1==0 || s2==0 || c2==0)
		return 0;
		
	for(l1=*ml;l1;l1=ListTail(l1))
	{
		List pp, pps;
		int  nu;
		Term m2;
		
		pps=CopyTerm(CompoundArg2(ListFirst(l1)));
		nu=(int)IntegerValue(CompoundArg1(ListFirst(l1)))/2;
		if(2*nu!=IntegerValue(CompoundArg1(ListFirst(l1))))
		{
			FreeAtomic(pps);
			continue;
		}
		
		for(l2=pps;l2;l2=ListTail(l2))
			if(CompoundArg1(ListFirst(l2))==s1)
				break;
		if(is_empty_list(l2) || IntegerValue(CompoundArg2(ListFirst(l2)))<2)
		{
			FreeAtomic(pps);
			continue;
		}
		if(IntegerValue(CompoundArg2(ListFirst(l2)))>2)
			SetCompoundArg(ListFirst(l2),2,
					NewInteger(IntegerValue(CompoundArg2(ListFirst(l2)))-2));
		else
			pps=CutFromList(pps,l2);
		
		for(l2=pps;l2;l2=ListTail(l2))
			if(CompoundArg1(ListFirst(l2))==s2)
				break;
		if(is_empty_list(l2) || IntegerValue(CompoundArg2(ListFirst(l2)))<2)
		{
			FreeAtomic(pps);
			continue;
		}
		if(IntegerValue(CompoundArg2(ListFirst(l2)))>2)
			SetCompoundArg(ListFirst(l2),2,
					NewInteger(IntegerValue(CompoundArg2(ListFirst(l2)))-2));
		else
			pps=CutFromList(pps,l2);
		
		
		pp=tri_rcos_mlt(CopyTerm(pps),s1,2);
		
		for(l2=*ml;l2;l2=ListTail(l2))
		{
			if(l1==l2)
				continue;
			if(EqualTerms(pp,CompoundArg2(ListFirst(l2))) &&
					CompoundArg1(ListFirst(l2))==NewInteger(-nu))
				break;
		}
		
		FreeAtomic(pp);
		
		if(is_empty_list(l2))
		{
			FreeAtomic(pps);
			continue;
		}
		
		pp=tri_rcos_mlt(CopyTerm(pps),s2,2);
		
		for(l3=*ml;l3;l3=ListTail(l3))
		{
			if(l1==l3 || l2==l3)
				continue;
			if(EqualTerms(pp,CompoundArg2(ListFirst(l3))) &&
					CompoundArg1(ListFirst(l3))==NewInteger(-nu))
				break;
		}
		
		FreeAtomic(pp);
		
		if(is_empty_list(l3))
		{
			FreeAtomic(pps);
			continue;
		}
		
		
		
		*ml=CutFromList(*ml,l1);
		*ml=CutFromList(*ml,l2);
		*ml=CutFromList(*ml,l3);
		
		m2=MakeCompound(A_MTERM,3);
		SetCompoundArg(m2,1,NewInteger(-nu));
		pps=tri_rcos_mlt(pps,s1,2);
		pps=tri_rcos_mlt(pps,c2,2);
		SetCompoundArg(m2,2,CopyTerm(pps));
		*ml=tri_rcos_add(*ml,m2);
		
		m2=MakeCompound(A_MTERM,3);
		SetCompoundArg(m2,1,NewInteger(-nu));
		pps=tri_rcos_mlt(pps,s2,2);
		pps=tri_rcos_mlt(pps,c1,2);
		SetCompoundArg(m2,2,pps);
		*ml=tri_rcos_add(*ml,m2);
		return 1;
		
	}
	
	return 0;
}


static int tri_rep_d(List *ml, int s1, Atom i1, int s2, Atom i2,
				Atom o1, Atom o2)
{
	List l1,l2;
	
	if(i1==0 || i2==0 || o1==0 || o2==0)
		return 0;
	
	for(l1=*ml;l1;l1=ListTail(l1))
	{
		List pp, pps;
		int  nu;
		Term m2;
		
		pps=CopyTerm(CompoundArg2(ListFirst(l1)));
		nu=(int)IntegerValue(CompoundArg1(ListFirst(l1)))*s1;		
		
		for(l2=pps;l2;l2=ListTail(l2))
			if(CompoundArg1(ListFirst(l2))==i1)
				break;
		
		if(is_empty_list(l2))
		{
			FreeAtomic(pps);
			continue;
		}
		
		if(CompoundArg2(ListFirst(l2))==NewInteger(1))
			pps=CutFromList(pps,l2);
		else
			SetCompoundArg(ListFirst(l2),2,
					NewInteger(IntegerValue(CompoundArg2(ListFirst(l2)))-1));
		
		pp=tri_rcos_mlt(CopyTerm(pps),i2,1);

		
		for(l2=*ml;l2;l2=ListTail(l2))
		{
			if(l1==l2)
				continue;
			if(EqualTerms(pp,CompoundArg2(ListFirst(l2))) &&
					CompoundArg1(ListFirst(l2))==NewInteger(nu*s2))
				break;
		}
		
		FreeAtomic(pp);
		
		if(is_empty_list(l2))
		{
			FreeAtomic(pps);
			continue;
		}
		
		*ml=CutFromList(*ml,l1);
		*ml=CutFromList(*ml,l2);
		
		m2=MakeCompound(A_MTERM,3);
		SetCompoundArg(m2,1,NewInteger(2*nu));
		pps=tri_rcos_mlt(pps,o1,1);
		pps=tri_rcos_mlt(pps,o2,1);
		SetCompoundArg(m2,2,pps);
		*ml=AppendLast(*ml,m2);
/*puts("+d");*/
		return 1;
		
	}
	
	return 0;
}

static int tri_rep_e(List *ml, int s1, Atom i11, Atom i12, 
				int s2, Atom i21, Atom i22, Atom o)
{
	List l1,l2,l3;
	
	if(o==0 || i11==0 || i12==0 || i21==0 || i22==0)
		return 0;
	
	for(l1=*ml;l1;l1=ListTail(l1))
	{
		List pp, pps;
		int  nu;
		Term m2;
		
		pps=CopyTerm(CompoundArg2(ListFirst(l1)));
		nu=(int)IntegerValue(CompoundArg1(ListFirst(l1)))*s1;		
		
		for(l2=pps;l2;l2=ListTail(l2))
			if(CompoundArg1(ListFirst(l2))==i11)
				break;
		
		for(l3=pps;l3;l3=ListTail(l3))
			if(CompoundArg1(ListFirst(l3))==i12)
				break;
		
		if(is_empty_list(l2) || is_empty_list(l3))
		{
			FreeAtomic(pps);
			continue;
		}
		
		if(CompoundArg2(ListFirst(l2))==NewInteger(1))
			pps=CutFromList(pps,l2);
		else
			SetCompoundArg(ListFirst(l2),2,
					NewInteger(IntegerValue(CompoundArg2(ListFirst(l2)))-1));
		
		if(CompoundArg2(ListFirst(l3))==NewInteger(1))
			pps=CutFromList(pps,l3);
		else
			SetCompoundArg(ListFirst(l3),2,
					NewInteger(IntegerValue(CompoundArg2(ListFirst(l3)))-1));
		
		pp=tri_rcos_mlt(CopyTerm(pps),i21,1);
		pp=tri_rcos_mlt(pp,i22,1);

		
		for(l2=*ml;l2;l2=ListTail(l2))
		{
			if(l1==l2)
				continue;
			if(EqualTerms(pp,CompoundArg2(ListFirst(l2))) &&
					CompoundArg1(ListFirst(l2))==NewInteger(nu*s2))
				break;
		}
		
		FreeAtomic(pp);
		
		if(is_empty_list(l2))
		{
			FreeAtomic(pps);
			continue;
		}
		
		*ml=CutFromList(*ml,l1);
		*ml=CutFromList(*ml,l2);
		
		m2=MakeCompound(A_MTERM,3);
		SetCompoundArg(m2,1,NewInteger(nu));
		pps=tri_rcos_mlt(pps,o,1);
		SetCompoundArg(m2,2,pps);
		*ml=AppendLast(*ml,m2);
/*puts("+e");*/
		return 1;
		
	}
	
	return 0;
}

static int tri_rep_f(List *ml, Atom sp, Atom sm, Atom cp, Atom cm, 
			int sn1, int sn2, int has1, Atom o)
{
	List l1,k1,k2,k3,k4,k5=0;
	
	if(o==0 || sm==0 || sp==0 || cm==0 || cp==0)
		return 0;
	
	for(l1=*ml;l1;l1=ListTail(l1))
	{
		List pp, pps;
		int  nu;
		Term m2;
		
		pps=CopyTerm(CompoundArg2(ListFirst(l1)));
		nu=(int)IntegerValue(CompoundArg1(ListFirst(l1)))*sn2/2;		
		
		for(k1=pps;k1;k1=ListTail(k1))
			if(CompoundArg1(ListFirst(k1))==cp)
				break;
		
		for(k2=pps;k2;k2=ListTail(k2))
			if(CompoundArg1(ListFirst(k2))==cm)
				break;
		
		for(k3=pps;k3;k3=ListTail(k3))
			if(CompoundArg1(ListFirst(k3))==sp)
				break;
		
		for(k4=pps;k4;k4=ListTail(k4))
			if(CompoundArg1(ListFirst(k4))==sm)
				break;
		
		if(is_empty_list(k1) || is_empty_list(k2) || 
				is_empty_list(k3) || is_empty_list(k4))
		{
			FreeAtomic(pps);
			continue;
		}
		
		if(CompoundArg2(ListFirst(k1))==NewInteger(1))
			pps=CutFromList(pps,k1);
		else
			SetCompoundArg(ListFirst(k1),2,
					NewInteger(IntegerValue(CompoundArg2(ListFirst(k1)))-1));
		
		if(CompoundArg2(ListFirst(k2))==NewInteger(1))
			pps=CutFromList(pps,k2);
		else
			SetCompoundArg(ListFirst(k2),2,
					NewInteger(IntegerValue(CompoundArg2(ListFirst(k2)))-1));
		
		if(CompoundArg2(ListFirst(k3))==NewInteger(1))
			pps=CutFromList(pps,k3);
		else
			SetCompoundArg(ListFirst(k3),2,
					NewInteger(IntegerValue(CompoundArg2(ListFirst(k3)))-1));
		
		if(CompoundArg2(ListFirst(k4))==NewInteger(1))
			pps=CutFromList(pps,k4);
		else
			SetCompoundArg(ListFirst(k4),2,
					NewInteger(IntegerValue(CompoundArg2(ListFirst(k4)))-1));
		
		k1=l1;
		
		
		pp=tri_rcos_mlt(CopyTerm(pps),sp,2);
		pp=tri_rcos_mlt(pp,sm,2);
		
		for(k2=*ml;k2;k2=ListTail(k2))
		{
			if(k1==k2)
				continue;
			if(EqualTerms(pp,CompoundArg2(ListFirst(k2))) &&
					CompoundArg1(ListFirst(k2))==NewInteger(2*sn1*nu))
				break;
		}
		
		FreeAtomic(pp);
		
		if(is_empty_list(k2))
		{
			FreeAtomic(pps);
			continue;
		}
		
		pp=tri_rcos_mlt(CopyTerm(pps),sm,2);
		
		for(k3=*ml;k3;k3=ListTail(k3))
		{
			if(k1==k3 || k2==k3)
				continue;
			if(EqualTerms(pp,CompoundArg2(ListFirst(k3))) &&
					CompoundArg1(ListFirst(k3))==NewInteger(-sn1*nu))
				break;
		}
		
		FreeAtomic(pp);
		
		if(is_empty_list(k3))
		{
			FreeAtomic(pps);
			continue;
		}
		
		pp=tri_rcos_mlt(CopyTerm(pps),sp,2);
		
		for(k4=*ml;k4;k4=ListTail(k4))
		{
			if(k1==k4 || k2==k4 || k3==k4)
				continue;
			if(EqualTerms(pp,CompoundArg2(ListFirst(k4))) &&
					CompoundArg1(ListFirst(k4))==NewInteger(-sn1*nu))
				break;
		}
		
		FreeAtomic(pp);
		
		if(is_empty_list(k4))
		{
			FreeAtomic(pps);
			continue;
		}
		
		
		if(has1)
		{
			for(k5=*ml;k5;k5=ListTail(k5))
			{
				if(k1==k5 || k2==k5 || k3==k5 || k4==k5)
					continue;
				if(EqualTerms(pps,CompoundArg2(ListFirst(k5))) &&
						CompoundArg1(ListFirst(k5))==NewInteger(nu))
					break;
			}

			if(is_empty_list(k5))
			{
				FreeAtomic(pps);
				continue;
			}
		}
		
		*ml=CutFromList(*ml,k1);
		*ml=CutFromList(*ml,k2);
		*ml=CutFromList(*ml,k3);
		*ml=CutFromList(*ml,k4);
		if(has1)
			*ml=CutFromList(*ml,k5);
		
		m2=MakeCompound(A_MTERM,3);
		SetCompoundArg(m2,1,NewInteger(nu));
		pps=tri_rcos_mlt(pps,o,2);
		SetCompoundArg(m2,2,pps);
		*ml=AppendLast(*ml,m2);
/*puts("+f");*/
		return 1;
		
	}
	
	return 0;
}

static int tri_rep_g(List *ml, Atom sp, Atom sm, Atom cp, Atom cm, 
			int sp_s, int sm_s, int cp_s, int cm_s,
			int sn1, int sn2, Atom o1, Atom o2)
{
	List l1,k1,k2,k3,k4;
	
	if(o1==0 || o2==0 || sm==0 || sp==0 || cm==0 || cp==0)
		return 0;
	
	for(l1=*ml;l1;l1=ListTail(l1))
	{
		List pp, pps;
		int  nu;
		Term m2;
		
		pps=CopyTerm(CompoundArg2(ListFirst(l1)));
		nu=(int)IntegerValue(CompoundArg1(ListFirst(l1)))*sp_s*cp_s;		
		
		for(k1=pps;k1;k1=ListTail(k1))
			if(CompoundArg1(ListFirst(k1))==cp)
				break;
				
		for(k3=pps;k3;k3=ListTail(k3))
			if(CompoundArg1(ListFirst(k3))==sp)
				break;
		
		
		if(is_empty_list(k1) || is_empty_list(k3))
		{
			FreeAtomic(pps);
			continue;
		}
		
		if(CompoundArg2(ListFirst(k1))==NewInteger(1))
			pps=CutFromList(pps,k1);
		else
			SetCompoundArg(ListFirst(k1),2,
					NewInteger(IntegerValue(CompoundArg2(ListFirst(k1)))-1));
				
		if(CompoundArg2(ListFirst(k3))==NewInteger(1))
			pps=CutFromList(pps,k3);
		else
			SetCompoundArg(ListFirst(k3),2,
					NewInteger(IntegerValue(CompoundArg2(ListFirst(k3)))-1));
		
		k1=l1;
		
		
		pp=tri_rcos_mlt(CopyTerm(pps),sm,1);
		pp=tri_rcos_mlt(pp,cm,1);
		
		for(k2=*ml;k2;k2=ListTail(k2))
		{
			if(k1==k2)
				continue;
			if(EqualTerms(pp,CompoundArg2(ListFirst(k2))) &&
					CompoundArg1(ListFirst(k2))==NewInteger(sn1*nu*sm_s*cm_s))
				break;
		}
		
		FreeAtomic(pp);
		
		if(is_empty_list(k2))
		{
			FreeAtomic(pps);
			continue;
		}
		
		pp=tri_rcos_mlt(CopyTerm(pps),sm,2);
		pp=tri_rcos_mlt(pp,sp,1);
		pp=tri_rcos_mlt(pp,cp,1);
		
		for(k3=*ml;k3;k3=ListTail(k3))
		{
			if(k1==k3 || k2==k3)
				continue;
			if(EqualTerms(pp,CompoundArg2(ListFirst(k3))) &&
					CompoundArg1(ListFirst(k3))==NewInteger(-2*nu*sp_s*cp_s))
				break;
		}
		
		FreeAtomic(pp);
		
		if(is_empty_list(k3))
		{
			FreeAtomic(pps);
			continue;
		}
		
		pp=tri_rcos_mlt(CopyTerm(pps),sp,2);
		pp=tri_rcos_mlt(pp,sm,1);
		pp=tri_rcos_mlt(pp,cm,1);
		
		for(k4=*ml;k4;k4=ListTail(k4))
		{
			if(k1==k4 || k2==k4 || k3==k4)
				continue;
			if(EqualTerms(pp,CompoundArg2(ListFirst(k4))) &&
					CompoundArg1(ListFirst(k4))==NewInteger(-2*sn1*nu*sm_s*cm_s))
				break;
		}
		
		FreeAtomic(pp);
		
		if(is_empty_list(k4))
		{
			FreeAtomic(pps);
			continue;
		}
				
		*ml=CutFromList(*ml,k1);
		*ml=CutFromList(*ml,k2);
		*ml=CutFromList(*ml,k3);
		*ml=CutFromList(*ml,k4);
		
		m2=MakeCompound(A_MTERM,3);
		SetCompoundArg(m2,1,NewInteger(sn2*nu));
		pps=tri_rcos_mlt(pps,o1,1);
		pps=tri_rcos_mlt(pps,o2,1);
		SetCompoundArg(m2,2,pps);
		*ml=AppendLast(*ml,m2);
/*puts("+g");*/
		return 1;
		
	}
	
	return 0;
}

static int tri_rep_h(List *ml, Atom i1, Atom i2, Atom i3, Atom i4, 
			int s_1, int s_2,  int has1, Atom o1, Atom o2, int s_o)
{
	List l1,k1,k2,k3=0;
	
	if(i1==0 || i2==0 || i3==0 || i4==0 || o1==0 || o2==0)
		return 0;
	
	for(l1=*ml;l1;l1=ListTail(l1))
	{
		List pp, pps;
		int  nu;
		Term m2;
		
		pps=CopyTerm(CompoundArg2(ListFirst(l1)));
		nu=(int)IntegerValue(CompoundArg1(ListFirst(l1)))*s_1;
		
		for(k1=pps;k1;k1=ListTail(k1))
			if(CompoundArg1(ListFirst(k1))==i1)
				break;
				
		for(k2=pps;k2;k2=ListTail(k2))
			if(CompoundArg1(ListFirst(k2))==i2)
				break;
		
		if(is_empty_list(k1) || is_empty_list(k2))
		{
			FreeAtomic(pps);
			continue;
		}
		
		if(i1==i2 && IntegerValue(CompoundArg2(ListFirst(k1)))<2)
		{
			FreeAtomic(pps);
			continue;
		}



		if(i1==i2)
		{
			if(CompoundArg2(ListFirst(k1))==NewInteger(2))
				pps=CutFromList(pps,k1);
			else
				SetCompoundArg(ListFirst(k1),2,
					NewInteger(IntegerValue(CompoundArg2(ListFirst(k1)))-2));
		}
		else
		{
			if(CompoundArg2(ListFirst(k1))==NewInteger(1))
				pps=CutFromList(pps,k1);
			else
				SetCompoundArg(ListFirst(k1),2,
					NewInteger(IntegerValue(CompoundArg2(ListFirst(k1)))-1));
			
			if(CompoundArg2(ListFirst(k2))==NewInteger(1))
				pps=CutFromList(pps,k2);
			else
				SetCompoundArg(ListFirst(k2),2,
					NewInteger(IntegerValue(CompoundArg2(ListFirst(k2)))-1));
		}
		
		k1=l1;
		
		if(i3==i4)
			pp=tri_rcos_mlt(CopyTerm(pps),i3,2);
		else
		{
			pp=tri_rcos_mlt(CopyTerm(pps),i3,1);
			pp=tri_rcos_mlt(pp,i4,1);
		}
		
				
		for(k2=*ml;k2;k2=ListTail(k2))
		{
			if(k1==k2)
				continue;
			if(EqualTerms(pp,CompoundArg2(ListFirst(k2))) &&
					CompoundArg1(ListFirst(k2))==NewInteger(nu*s_2))
				break;
		}
		
		FreeAtomic(pp);
		
		if(is_empty_list(k2))
		{
			FreeAtomic(pps);
			continue;
		}
		
		if(has1)
		{
			for(k3=*ml;k3;k3=ListTail(k3))
			{
				if(k1==k3 || k2==k3)
					continue;
				if(EqualTerms(pps,CompoundArg2(ListFirst(k3))) &&
						CompoundArg1(ListFirst(k3))==NewInteger(nu))
					break;
			}

			if(is_empty_list(k3))
			{
				FreeAtomic(pps);
				continue;
			}
		}
		
		*ml=CutFromList(*ml,k1);
		*ml=CutFromList(*ml,k2);
		if(has1)
			*ml=CutFromList(*ml,k3);
		
		m2=MakeCompound(A_MTERM,3);
		SetCompoundArg(m2,1,NewInteger(nu*s_o));
		pps=tri_rcos_mlt(pps,o1,1);
		pps=tri_rcos_mlt(pps,o2,1);
		SetCompoundArg(m2,2,pps);
		*ml=AppendLast(*ml,m2);
/*puts("+h");*/
		return 1;
		
	}

	return 0;
}


int tri_rep_cs3(List *ml, Atom s1, Atom s2, Atom c1, Atom c2,
		Atom sapb, Atom samb, Atom capb, Atom camb,
		int sapbs, int sambs, int capbs, int cambs);

int tri_opti_2(Term v2d, List *m2l)
{
	
	Atom sa1, sa2, sa3, ca1, ca2, ca3;
	Atom s2a1, s2a2, s2a3, c2a1, c2a2, c2a3;
	Atom sapb, capb, samb, camb;
	int sapb_s, capb_s, samb_s, camb_s, c2a1_s=1, c2a2_s=1, c2a3_s=1;
	
	int chs;
	
	if(ListLength(*m2l)==1)
		return 1;
	
	{
		int i;
		i=(int)IntegerValue(CompoundArg1(v2d));

		sa1=CompoundArg1(ListNth(tri_si_co_list,i));
		ca1=CompoundArg2(ListNth(tri_si_co_list,i));
		s2a1=CompoundArgN(ListNth(tri_si_co_list,i),3);
		c2a1=CompoundArgN(ListNth(tri_si_co_list,i),5);
		if(c2a1)
		{
			c2a1_s=CompoundName(c2a1)==OPR_PLUS?1:-1;
			c2a1=CompoundArg1(c2a1);
		}

		i=(int)IntegerValue(CompoundArg2(v2d));

		sa2=CompoundArg1(ListNth(tri_si_co_list,i));
		ca2=CompoundArg2(ListNth(tri_si_co_list,i));
		s2a2=CompoundArgN(ListNth(tri_si_co_list,i),3);
		c2a2=CompoundArgN(ListNth(tri_si_co_list,i),5);
		if(c2a2)
		{
			c2a2_s=CompoundName(c2a2)==OPR_PLUS?1:-1;
			c2a2=CompoundArg1(c2a2);
		}
		
		if(tri_third_angle)
		{
			i=tri_third_angle;
			sa3=CompoundArg1(ListNth(tri_si_co_list,i));
			ca3=CompoundArg2(ListNth(tri_si_co_list,i));
			s2a3=CompoundArgN(ListNth(tri_si_co_list,i),3);
			c2a3=CompoundArgN(ListNth(tri_si_co_list,i),5);
			if(c2a3)
			{
				c2a3_s=CompoundName(c2a3)==OPR_PLUS?1:-1;
				c2a3=CompoundArg1(c2a3);
			}
		}
		else
			sa3=ca3=s2a3=c2a3=c2a3_s=0;
		
		sapb=CompoundArg1(CompoundArgN(v2d,4));
		capb=CompoundArg1(CompoundArgN(v2d,6));
		samb=CompoundArg1(CompoundArgN(v2d,5));
		camb=CompoundArg1(CompoundArgN(v2d,7));

		sapb_s=CompoundName(CompoundArgN(v2d,4))==OPR_PLUS ? 1 : -1;
		capb_s=CompoundName(CompoundArgN(v2d,6))==OPR_PLUS ? 1 : -1;
		samb_s=CompoundName(CompoundArgN(v2d,5))==OPR_PLUS ? 1 : -1;
		camb_s=CompoundName(CompoundArgN(v2d,7))==OPR_PLUS ? 1 : -1;
	}
	
	
	tri_rep_cs3(m2l, sa1, sa2, ca1, ca2, sapb, samb, capb, camb,
				sapb_s, samb_s, capb_s, camb_s);

	
	do
	{
		chs=0;
				
		chs|=tri_rep_f(m2l, sapb, samb, capb, camb, 1, 
				-samb_s*sapb_s*camb_s*capb_s, 1, c2a1);
		chs|=tri_rep_f(m2l, sapb, samb, capb, camb, 1, 
				 samb_s*sapb_s*camb_s*capb_s, 1, c2a2);
		chs|=tri_rep_f(m2l, sapb, samb, capb, camb,-1, 
				 samb_s*sapb_s*camb_s*capb_s, 0, s2a1);
		chs|=tri_rep_f(m2l, sapb, samb, capb, camb,-1, 
				-samb_s*sapb_s*camb_s*capb_s, 0, s2a2);
		
		chs|=tri_rep_g(m2l, sapb, samb, capb, camb,
				sapb_s, samb_s, capb_s, camb_s, 1,c2a1_s,s2a1,c2a1);
		chs|=tri_rep_g(m2l, sapb, samb, capb, camb,
				sapb_s, samb_s, capb_s, camb_s,-1,c2a2_s,s2a2,c2a2);
		
		chs|=tri_rep_h(m2l, sapb, sapb, samb, samb, -1, -1, 1,
								c2a1,c2a2,c2a1_s*c2a2_s);
		chs|=tri_rep_h(m2l, sapb, sapb, samb, samb,  1, -1, 0,
								s2a1,s2a2,1);
		chs|=tri_rep_h(m2l, capb, sapb, camb, samb,  
				capb_s*sapb_s, -camb_s*samb_s, 0, c2a1,s2a2,c2a1_s);
		chs|=tri_rep_h(m2l, capb, sapb, camb, samb,  
				capb_s*sapb_s,  camb_s*samb_s, 0, c2a2,s2a1,c2a2_s);
		
		chs|=tri_rep_e(m2l, capb_s*camb_s,capb,camb,
					-sapb_s*samb_s,sapb,samb,c2a1);
		chs|=tri_rep_e(m2l, capb_s*camb_s,camb,capb,
					 sapb_s*samb_s,sapb,samb,c2a2);
		chs|=tri_rep_e(m2l, capb_s*samb_s,capb,samb,
					 camb_s*sapb_s,camb,sapb,s2a1);
		chs|=tri_rep_e(m2l,-capb_s*samb_s,capb,samb,
					 camb_s*sapb_s,camb,sapb,s2a2);
		
		chs|=tri_rep_d(m2l, camb_s,camb,-capb_s,capb,sa1,sa2);
		chs|=tri_rep_d(m2l, samb_s,samb, sapb_s,sapb,sa1,ca2);
		chs|=tri_rep_d(m2l,-samb_s,samb, sapb_s,sapb,ca1,sa2);
		chs|=tri_rep_d(m2l, camb_s,camb, capb_s,capb,ca1,ca2);
		
		if(chs)
			continue;
		
		chs|=tri_rep_c(m2l,sa1,ca1,0);
		chs|=tri_rep_c(m2l,sa2,ca2,0);
		chs|=tri_rep_c(m2l,sa3,ca3,0);
		chs|=tri_rep_c(m2l,ca1,sa1,0);
		chs|=tri_rep_c(m2l,ca2,sa2,0);
		chs|=tri_rep_c(m2l,ca3,sa3,0);
		chs|=tri_rep_c(m2l,sapb,capb,0);
		chs|=tri_rep_c(m2l,samb,camb,0);
		
		chs|=tri_rep_cc(m2l,sa1,sa1,ca1,ca2);
		chs|=tri_rep_cc(m2l,sapb,samb,capb,camb);
		
		if(chs)
			continue;
		
		chs|=tri_rep_c(m2l,sa1,c2a1,c2a1_s);
		chs|=tri_rep_c(m2l,sa2,c2a2,c2a2_s);
		chs|=tri_rep_c(m2l,sa3,c2a3,c2a3_s);
		
		chs|=tri_rep_c2(m2l,sa1,ca1,c2a1);
		chs|=tri_rep_c2(m2l,sa2,ca2,c2a2);
		chs|=tri_rep_c2(m2l,sa3,ca3,c2a3);
		
	} while(chs);
	
	
	return ListLength(*m2l);
}


int tri_opti_1(int v, List *m2l)
{
	Atom ca, sa, c2a, s2a, c2a_s=1;
	int chs;
	
	sa=CompoundArg1(ListNth(tri_si_co_list,v));
	ca=CompoundArg2(ListNth(tri_si_co_list,v));
	s2a=CompoundArgN(ListNth(tri_si_co_list,v),3);
	c2a=CompoundArgN(ListNth(tri_si_co_list,v),5);
	if(c2a)
	{
		c2a_s=CompoundName(c2a)==OPR_PLUS?1:-1;
		c2a=CompoundArg1(c2a);
	}

	do
	{
		chs=0;
				
		chs|=tri_rep_c(m2l,sa,ca,0);
		chs|=tri_rep_c2(m2l,sa,ca,c2a);
/*		chs|=tri_rep_c(m2l,sa,c2a,c2a_s);*/
		
	} while(chs);
	
	
	return ListLength(*m2l);
}

static int pcmp(Term p1, Term p2)
        {
                if(CompoundArg1(p1)==A_I)
                        return -1;
                if(CompoundArg1(p2)==A_I)
                        return 1;
                return strcmp(AtomValue(CompoundArg1(p1)),
                                AtomValue(CompoundArg1(p2)));
        }
		
/*
static int pcmp(Term p1, Term p2)
	{
	return strcmp(AtomValue(CompoundArg1(p1)),AtomValue(CompoundArg1(p2)));
	}
*/
		
void alg2_red_comsico(Term a2)
{
	List csf;
	List l1,l2,l3;
	int rden=1, rnum=1, i, c=0;
	
	csf=ConsumeCompoundArg(a2,3);
	if(csf==0)
		return;
	
/*
	WriteTerm(csf);WriteTerm(CompoundArg2(a2));
	puts("");
*/
			
xyz:
	
	for(l1=csf;l1;l1=ListTail(l1))
		for(l2=ListTail(l1);l2;l2=ListTail(l2))
		{
			Atom a1,a2;
			int po1,po2;
			po1=(int)IntegerValue(CompoundArg2(ListFirst(l1)));
			po2=(int)IntegerValue(CompoundArg2(ListFirst(l2)));
			if(po1!=po2 && po1!=-po2)
				continue;
			a1=CompoundArg1(ListFirst(l1));
			a2=CompoundArg1(ListFirst(l2));
			for(l3=tri_si_co_list;l3;l3=ListTail(l3))
				if(((CompoundArg1(ListFirst(l3))==a1 
						&& CompoundArg2(ListFirst(l3))==a2) ||
				    (CompoundArg1(ListFirst(l3))==a2 
						&& CompoundArg2(ListFirst(l3))==a1)))
				{ 
					if(po1==po2 && CompoundArgN(ListFirst(l3),3))
					{
						if(ListLength(csf)>2)
						{
							csf=CutFromList(csf,l1);
							csf=CutFromList(csf,l2);
						}
						else
						{
							FreeAtomic(csf);
							csf=NewList();
						}
						csf=AppendFirst(csf,MakeCompound2(OPR_POW,
								CompoundArgN(ListFirst(l3),3),NewInteger(po1)));
						if(po1>0)
							for(i=0;i<po1;i++)
								rden*=2;
						else
							for(i=0;i>po1;i--)
								rnum*=2;
						c++;
						goto xyz;
					}
					if(po1==-po2 && CompoundArgN(ListFirst(l3),6))
					{
						int po;
						if(a1==CompoundArg1(ListFirst(l3)))
							po=po1;
						else
							po=po2;
						if(ListLength(csf)>2)
						{
							csf=CutFromList(csf,l1);
							csf=CutFromList(csf,l2);
						}
						else
						{
							FreeAtomic(csf);
							csf=NewList();
						}
						csf=AppendFirst(csf,MakeCompound2(OPR_POW,
								CompoundArgN(ListFirst(l3),6),NewInteger(po)));
						c++;
						goto xyz;
					}
				}
		}
	
	csf=SortedList(csf,pcmp);
	if(rnum>1 || rden>1)
	{
		int num, den;
		num=(int)IntegerValue(CompoundArg1(CompoundArg2(a2)));
		den=(int)IntegerValue(CompoundArg2(CompoundArg2(a2)));
		den*=rden;
		num*=rnum;
		rnum=(int)gcf(num,den);
		num/=rnum;
		den/=rnum;
		SetCompoundArg(CompoundArg2(a2),1,NewInteger(num));
		SetCompoundArg(CompoundArg2(a2),2,NewInteger(den));
	}
/*
	if(c)
	{
	WriteTerm(csf);WriteTerm(CompoundArg2(a2));
	puts(" -- redcomsico --- done");
	}
*/
	SetCompoundArg(a2,3,csf);
	return;
}
