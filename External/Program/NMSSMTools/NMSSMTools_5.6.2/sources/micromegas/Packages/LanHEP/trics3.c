#include <stdio.h>
#include "lanhep.h"

Atom tri_dummy_prm(List);
List tri_rcos_mlt(List, Atom, int);

    /* 1 --- sa*cb**3 + ca*sb**3    = sa*cb - sa*cb*sb**2 + ca*sb**3
	 * 2 --- sa*cb**3 - ca*sb**3    = sa*cb - sa*cb*sb**2 - ca*sb**3
	 * 3 --- ca*cb**3 + sa*sb**3    = ca*cb - ca*cb*sb**2 + sa*sb**3
	 * 4 --- ca*cb**3 - sa*sb**3    = ca*cb - ca*cb*sb**2 + sa*sb**3 */

static int tri_find_pure_ab(List *ml, Atom s1, Atom s2, Atom c1, Atom c2)
{
	
	Term m2, m2_1=0, m2_2=0, m2_3=0;
	List l1,l2;
	Atom p11=0, p12=0, p21=0, p22=0, p23=0, p31=0, p32=0;
	Atom pc1=0, pc2=0, pc3=0;
	int tp=0, rs=0;
	int comf1=0, comf=0;
		
	
	if(ListLength(*ml)!=3)
		return 0;
	
	for(l1=(*ml);l1;l1=ListTail(l1))
	{
		m2=ListFirst(l1);
		l2=CompoundArg2(m2);
		if(ListLength(l2)==3)
		{
			if(m2_2)
				return 0;
			m2_2=m2;
			continue;
		}
		
		if(ListLength(l2)!=2)
			return 0;
		if(CompoundArg2(ListFirst(l2))==NewInteger(1) &&
				CompoundArg2(ListFirst(ListTail(l2)))==NewInteger(1))
		{
			if(m2_1)
				return 0;
			m2_1=m2;
			continue;
		}
		
		if((CompoundArg2(ListFirst(l2))==NewInteger(1) &&
				CompoundArg2(ListFirst(ListTail(l2)))==NewInteger(3)) ||
			(CompoundArg2(ListFirst(l2))==NewInteger(3) &&
				CompoundArg2(ListFirst(ListTail(l2)))==NewInteger(1)))
		{
			if(m2_3)
				return 0;
			m2_3=m2;
			continue;
		}
		return 0;
	}
	
	if(m2_1 ==0 || m2_2==0 || m2_3==0)
		return 0;
	
	comf=(int)IntegerValue(CompoundArg1(m2_1));
	comf1=(int)IntegerValue(CompoundArg1(m2_3));
	if(IntegerValue(CompoundArg1(m2_2))!=-comf || (comf!=comf1 && comf!=-comf1))
		return 0;
	
	
	l1=CompoundArg2(m2_3);
	if(CompoundArg2(ListFirst(l1))==NewInteger(1) &&
			CompoundArg2(ListFirst(ListTail(l1)))==NewInteger(3))
	{
		p31=CompoundArg1(ListFirst(l1));
		p32=CompoundArg1(ListFirst(ListTail(l1)));
	}
	if(CompoundArg2(ListFirst(l1))==NewInteger(3) &&
			CompoundArg2(ListFirst(ListTail(l1)))==NewInteger(1))
	{
		p32=CompoundArg1(ListFirst(l1));
		p31=CompoundArg1(ListFirst(ListTail(l1)));
	}
	if(p31==0 || p32==0)
		return 0;
	
	if(p31==c1 && p32==s2)
	{
		if(comf==comf1)
			tp=1;
		else
			tp=2;
		pc1=s1;
		pc2=c2;
		pc3=s2;
	}
	if(p31==c2 && p32==s1)
	{
		if(comf==comf1)
			tp=1;
		else
			tp=2;
		rs=1;
		pc1=s2;
		pc2=c1;
		pc3=s1;
	}
	if(p31==s1 && p32==s2)
	{
		if(comf==comf1)
			tp=3;
		else
			tp=4;
		pc1=c1;
		pc2=c2;
		pc3=s2;
	}
	if(p31==s2 && p32==s1)
	{
		if(comf==comf1)
			tp=3;
		else
			tp=4;
		rs=1;
		pc1=c2;
		pc2=c1;
		pc3=s1;
	}
	
	if(tp==0)
		return 0;
		
	if(rs)
	{
		Atom tmp;
		tmp=s1;
		s1=s2;
		s2=tmp;
		tmp=c1;
		c1=c2;
		c2=tmp;
	}

/*	
	puts("");WriteTerm(*ml);printf(" t=%d r=%d\n\n",tp,rs);
*/
	
	l1=CompoundArg2(m2_1);
	p11=CompoundArg1(ListFirst(l1));
	p12=CompoundArg1(ListFirst(ListTail(l1)));
	
	if((p11!=pc1 && p11!=pc2) || (p12!=pc1 && p12!=pc2))
		return 0;
	
	
	for(l1=CompoundArg2(m2_2);l1;l1=ListTail(l1))
	{
		if(CompoundArg2(ListFirst(l1))==NewInteger(2))
		{
			if(p21)
				return 0;
			p21=CompoundArg1(ListFirst(l1));
			continue;
		}
		if(CompoundArg2(ListFirst(l1))!=NewInteger(1))
			return 0;
		if(p22==0)
		{
			p22=CompoundArg1(ListFirst(l1));
			continue;
		}
		if(p23)
			return 0;
		p23=CompoundArg1(ListFirst(l1));
	}
	
	if(p21!=pc3)
		return 0;
	if((p22!=pc1 && p22!=pc2) || (p23!=pc1 && p23!=pc2))
		return 0;
	
		
	m2_1=MakeCompound(A_MTERM,3);
	SetCompoundArg(m2_1,1,NewInteger(1));
	SetCompoundArg(m2_1,2,MakeList2(
			MakeCompound2(OPR_POW,tp<3?s1:c1,NewInteger(1)),
			MakeCompound2(OPR_POW,        c2,NewInteger(3))));

	m2_2=MakeCompound(A_MTERM,3);
	SetCompoundArg(m2_2,1,NewInteger((tp%2)?1:-1));
	SetCompoundArg(m2_2,2,MakeList2(
			MakeCompound2(OPR_POW,tp<3?c1:s1,NewInteger(1)),
			MakeCompound2(OPR_POW,        s2,NewInteger(3))));

	l2=MakeList2(m2_1,m2_2);

	m2=MakeCompound(A_MTERM,3);
	SetCompoundArg(m2,1,NewInteger(comf));
	SetCompoundArg(m2,2,MakeList1(
			MakeCompound2(OPR_POW,tri_dummy_prm(l2),NewInteger(1))));

	FreeAtomic(*ml);
	*ml=NewList();
	*ml=AppendLast(*ml,m2);
	return 1;
}

		
static List find_3_mono_pm(List ml, Atom sp, Atom sm, Atom cp, Atom cm, int type)
{
	/* 1 --- sa*cb**3 + ca*sb**3  *
	 * 2 --- sa*cb**3 - ca*sb**3  *
	 * 3 --- ca*cb**3 + sa*sb**3  *
	 * 4 --- ca*cb**3 - sa*sb**3  */
	
	List ret, l1, pp;
	
	ret=NewList();

	pp=0;
	pp=tri_rcos_mlt(pp,sp,1);
	pp=tri_rcos_mlt(pp,sm,1);
	pp=tri_rcos_mlt(pp,cp,1);
	pp=tri_rcos_mlt(pp,cm,1);
	for(l1=ml;l1;l1=ListTail(l1))
		if(EqualTerms(pp,CompoundArg2(ListFirst(l1))))
			break;
	FreeAtomic(pp);
	if(is_empty_list(l1))
		goto fail;
	ret=AppendLast(ret,l1);

	pp=0;
	pp=tri_rcos_mlt(pp,sp,(type==1||type==4)?1:3);
	pp=tri_rcos_mlt(pp,sm,(type==1||type==4)?3:1);
	pp=tri_rcos_mlt(pp,cp,1);
	pp=tri_rcos_mlt(pp,cm,1);
	for(l1=ml;l1;l1=ListTail(l1))
		if(EqualTerms(pp,CompoundArg2(ListFirst(l1))))
			break;
	FreeAtomic(pp);
	if(is_empty_list(l1))
		goto fail;
	ret=AppendLast(ret,l1);

	pp=0;
	pp=tri_rcos_mlt(pp,sp,2);
	for(l1=ml;l1;l1=ListTail(l1))
		if(EqualTerms(pp,CompoundArg2(ListFirst(l1))))
			break;
	FreeAtomic(pp);
	if(is_empty_list(l1))
		goto fail;
	ret=AppendLast(ret,l1);

	pp=0;
	pp=tri_rcos_mlt(pp,sm,2);
	for(l1=ml;l1;l1=ListTail(l1))
		if(EqualTerms(pp,CompoundArg2(ListFirst(l1))))
			break;
	FreeAtomic(pp);
	if(is_empty_list(l1))
		goto fail;
	ret=AppendLast(ret,l1);

	pp=0;
	pp=tri_rcos_mlt(pp,(type==1||type==4)?sm:sp,4);
	for(l1=ml;l1;l1=ListTail(l1))
		if(EqualTerms(pp,CompoundArg2(ListFirst(l1))))
			break;
	FreeAtomic(pp);
	if(is_empty_list(l1))
		goto fail;
	ret=AppendLast(ret,l1);

	pp=0;
	pp=tri_rcos_mlt(pp,sp,(type==1||type==4)?2:4);
	pp=tri_rcos_mlt(pp,sm,(type==1||type==4)?4:2);
	for(l1=ml;l1;l1=ListTail(l1))
		if(EqualTerms(pp,CompoundArg2(ListFirst(l1))))
			break;
	FreeAtomic(pp);
	if(is_empty_list(l1))
		goto fail;
	ret=AppendLast(ret,l1);

	pp=0;
	pp=tri_rcos_mlt(pp,sp,2);
	pp=tri_rcos_mlt(pp,sm,2);
	for(l1=ml;l1;l1=ListTail(l1))
		if(EqualTerms(pp,CompoundArg2(ListFirst(l1))))
			break;
	FreeAtomic(pp);
	if(is_empty_list(l1))
		goto fail;
	ret=AppendLast(ret,l1);

	if(type<3)
		return ret;

	pp=0;
	for(l1=ml;l1;l1=ListTail(l1))
		if(EqualTerms(pp,CompoundArg2(ListFirst(l1))))
			break;
	if(is_empty_list(l1))
		goto fail;
	ret=AppendLast(ret,l1);

	return ret;
	
fail:
	
	RemoveList(ret);
	return 0;
}
		
		
static List find_3_mono_ab(List ml, Atom s1, Atom s2, Atom c1, Atom c2, int type)
{
	/* 1 --- sa*cb**3 + ca*sb**3  *
	 * 2 --- sa*cb**3 - ca*sb**3  *
	 * 3 --- ca*cb**3 + sa*sb**3  *
	 * 4 --- ca*cb**3 - sa*sb**3  */
	
	List ret, l1, pp;
	
	ret=NewList();

	pp=0;
	pp=tri_rcos_mlt(pp,s1,1);
	pp=tri_rcos_mlt(pp,s2,3);
	pp=tri_rcos_mlt(pp,c1,1);
	pp=tri_rcos_mlt(pp,c2,1);
	for(l1=ml;l1;l1=ListTail(l1))
		if(EqualTerms(pp,CompoundArg2(ListFirst(l1))))
			break;
	FreeAtomic(pp);
	if(is_empty_list(l1))
		goto fail;
	ret=AppendLast(ret,l1);

	pp=0;
	pp=tri_rcos_mlt(pp,s1,1);
	pp=tri_rcos_mlt(pp,s2,5);
	pp=tri_rcos_mlt(pp,c1,1);
	pp=tri_rcos_mlt(pp,c2,1);
	for(l1=ml;l1;l1=ListTail(l1))
		if(EqualTerms(pp,CompoundArg2(ListFirst(l1))))
			break;
	FreeAtomic(pp);
	if(is_empty_list(l1))
		goto fail;
	ret=AppendLast(ret,l1);

	pp=0;
	pp=tri_rcos_mlt(pp,s1,2);
	for(l1=ml;l1;l1=ListTail(l1))
		if(EqualTerms(pp,CompoundArg2(ListFirst(l1))))
			break;
	FreeAtomic(pp);
	if(is_empty_list(l1))
		goto fail;
	ret=AppendLast(ret,l1);

	pp=0;
	pp=tri_rcos_mlt(pp,s1,2);
	pp=tri_rcos_mlt(pp,s2,2);
	for(l1=ml;l1;l1=ListTail(l1))
		if(EqualTerms(pp,CompoundArg2(ListFirst(l1))))
			break;
	FreeAtomic(pp);
	if(is_empty_list(l1))
		goto fail;
	ret=AppendLast(ret,l1);

	pp=0;
	pp=tri_rcos_mlt(pp,s1,2);
	pp=tri_rcos_mlt(pp,s2,4);
	for(l1=ml;l1;l1=ListTail(l1))
		if(EqualTerms(pp,CompoundArg2(ListFirst(l1))))
			break;
	FreeAtomic(pp);
	if(is_empty_list(l1))
		goto fail;
	ret=AppendLast(ret,l1);

	pp=0;
	pp=tri_rcos_mlt(pp,s1,2);
	pp=tri_rcos_mlt(pp,s2,6);
	for(l1=ml;l1;l1=ListTail(l1))
		if(EqualTerms(pp,CompoundArg2(ListFirst(l1))))
			break;
	FreeAtomic(pp);
	if(is_empty_list(l1))
		goto fail;
	ret=AppendLast(ret,l1);

	pp=0;
	pp=tri_rcos_mlt(pp,s2,6);
	for(l1=ml;l1;l1=ListTail(l1))
		if(EqualTerms(pp,CompoundArg2(ListFirst(l1))))
			break;
	FreeAtomic(pp);
	if(is_empty_list(l1))
		goto fail;
	ret=AppendLast(ret,l1);

	if(type<3)
		return ret;

	pp=0;
	pp=tri_rcos_mlt(pp,s2,4);
	for(l1=ml;l1;l1=ListTail(l1))
		if(EqualTerms(pp,CompoundArg2(ListFirst(l1))))
			break;
	FreeAtomic(pp);
	if(is_empty_list(l1))
		goto fail;
	ret=AppendLast(ret,l1);

	pp=0;
	pp=tri_rcos_mlt(pp,s2,2);
	for(l1=ml;l1;l1=ListTail(l1))
		if(EqualTerms(pp,CompoundArg2(ListFirst(l1))))
			break;
	FreeAtomic(pp);
	if(is_empty_list(l1))
		goto fail;
	ret=AppendLast(ret,l1);

	pp=0;
	for(l1=ml;l1;l1=ListTail(l1))
		if(EqualTerms(pp,CompoundArg2(ListFirst(l1))))
			break;
	if(is_empty_list(l1))
		goto fail;
	ret=AppendLast(ret,l1);

	return ret;
	
	
fail:
	
	RemoveList(ret);
	return 0;
}

static void find_3_cnf_pm(List ll, int type, int sign, int *cnf, int *fno, 
		List *ml)
{
	List l;
	int llen, nr[8],nd[8];
	int i;
	
	llen=ListLength(ll);
	
	for(i=1,l=ll;i<=llen;i++,l=ListTail(l))
		nr[i-1]=(int)IntegerValue(CompoundArg1(ListFirst(ListFirst(l))));
	
	nd[0]=(type<3?2:4)*sign;
	nd[1]=(type<3?2:-2)*sign;
	nd[2]=type<3?1:-4;
	nd[3]=type<3?1:-4;
	nd[4]=type<3?-1:1;
	nd[5]=type<3?2:-2;
	nd[6]=type<3?1:5;
	nd[7]=4;
	
	if(ml!=NULL)
	{
		for(i=1,l=ll;i<=llen;i++,l=ListTail(l))
		{
			if((*cnf)*nd[i-1]==nr[i-1])
				*ml=CutFromList(*ml,ListFirst(l));
			else
				SetCompoundArg(ListFirst(ListFirst(l)),1,
						NewInteger(nr[i-1]-(*cnf)*nd[i-1]));
		}
		return ;
	}
	
	*fno=0;
	for(i=0;i<llen;i++)
	{
		int cc, fno1=0, j;
		cc=nr[i]/nd[i];
		for(j=0;j<llen;j++)
			if(nr[j]==nd[j]*cc)
				fno1++;
		if(fno1>(*fno))
		{
			*fno=fno1;
			*cnf=cc;
		}
	}
	
	*fno=llen-*fno;
	
}

static void find_3_cnf_ab(List ll, int type, int *cnf, int *fno, 
		List *ml)
{
	List l;
	int llen, nr[10],nd[10];
	int i;
	
	llen=ListLength(ll);
	
	for(i=1,l=ll;i<=llen;i++,l=ListTail(l))
		nr[i-1]=(int)IntegerValue(CompoundArg1(ListFirst(ListFirst(l))));
	
	nd[0]=(type%2)?2:-2;
	nd[1]=(type%2)?-2:2;
	nd[2]=type<3?1:-1;
	nd[3]=type<3?-3:3;
	nd[4]=type<3?3:-3;
	nd[5]=type<3?-2:2;
	nd[6]=type<3?1:-1;
	nd[7]=3;
	nd[8]=-3;
	nd[9]=1;
	
	
	if(ml!=NULL)
	{
		
		for(i=1,l=ll;i<=llen;i++,l=ListTail(l))
		{
			if((*cnf)*nd[i-1]==nr[i-1])
				*ml=CutFromList(*ml,ListFirst(l));
			else
				SetCompoundArg(ListFirst(ListFirst(l)),1,
						NewInteger(nr[i-1]-(*cnf)*nd[i-1]));
		}
		
		return ;
	}
	
	*fno=0;
	for(i=0;i<llen;i++)
	{
		int cc, fno1=0, j;
		cc=nr[i]/nd[i];
		for(j=0;j<llen;j++)
			if(nr[j]==nd[j]*cc)
				fno1++;
		if(fno1>(*fno))
		{
			*fno=fno1;
			*cnf=cc;
		}
	}
	
	*fno=llen-*fno;
	
}

int tri_rep_cs3(List *ml, Atom s1, Atom s2, Atom c1, Atom c2,
		Atom sapb, Atom samb, Atom capb, Atom camb,
		int sapbs, int sambs, int capbs, int cambs)
{
	List l1,l2;
	int cnf, fno;
	int has_pm=0;
	
	if(tri_find_pure_ab(ml,s1,s2,c1,c2))
		return 1;
	
	if(s1==0 || s2==0 || sapb==0 || capb==0)
		return 0;
	
	for(l1=*ml;l1;l1=ListTail(l1))
		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
	{
		Atom p;
		p=CompoundArg1(ListFirst(l2));
		if(p==sapb || p==samb || p==capb || p==camb)
		{
			has_pm++;
			break;
		}
	}
	
	if(has_pm)
	{
		int tp=0, rs=0, bfit=10, bcnf=0;
		List zl[4];
		zl[0]=zl[1]=zl[2]=zl[3]=0;
		
		if((l1=find_3_mono_pm(*ml,sapb,samb,capb,camb,3)))
		{
			
			find_3_cnf_pm(l1,3, sapbs*sambs*capbs*cambs,&cnf,&fno,0);
			if(fno<bfit)
				tp=3,rs=0,bfit=fno,bcnf=cnf;
			find_3_cnf_pm(l1,3,-sapbs*sambs*capbs*cambs,&cnf,&fno,0);
			if(fno<bfit)
				tp=3,rs=1,bfit=fno,bcnf=cnf;
			zl[2]=l1;
		}
		if((l1=find_3_mono_pm(*ml,sapb,samb,capb,camb,4)))
		{
			
			find_3_cnf_pm(l1,4, sapbs*sambs*capbs*cambs,&cnf,&fno,0);
			if(fno<bfit)
				tp=4,rs=0,bfit=fno,bcnf=cnf;
			find_3_cnf_pm(l1,4,-sapbs*sambs*capbs*cambs,&cnf,&fno,0);
			if(fno<bfit)
				tp=4,rs=1,bfit=fno,bcnf=cnf;
			zl[3]=l1;
		}

		if((l1=find_3_mono_pm(*ml,sapb,samb,capb,camb,1)))
		{
			
			find_3_cnf_pm(l1,1, sapbs*sambs*capbs*cambs,&cnf,&fno,0);
			if(fno<bfit)
				tp=1,rs=0,bfit=fno,bcnf=cnf;
			find_3_cnf_pm(l1,1,-sapbs*sambs*capbs*cambs,&cnf,&fno,0);
			if(fno<bfit)
				tp=1,rs=1,bfit=fno,bcnf=cnf;
			zl[0]=l1;
		}
		if((l1=find_3_mono_pm(*ml,sapb,samb,capb,camb,2)))
		{
			
			find_3_cnf_pm(l1,2, sapbs*sambs*capbs*cambs,&cnf,&fno,0);
			if(fno<bfit)
				tp=2,rs=0,bfit=fno,bcnf=cnf;
			find_3_cnf_pm(l1,2,-sapbs*sambs*capbs*cambs,&cnf,&fno,0);
			if(fno<bfit)
				tp=2,rs=1,bfit=fno,bcnf=cnf;
			zl[1]=l1;
		}
		
		if(bfit<5)
		{
			Term m2,m2_1,m2_2;
/*			
printf("opti_3: type=%d, revers=%d, best fit=%d\n",tp,rs,bfit);			
*/			
			if(rs)
			{
				Atom tmp;
				tmp=s1;
				s1=s2;
				s2=tmp;
				tmp=c1;
				c1=c2;
				c2=tmp;
			}
			
			l1=zl[tp-1];
			find_3_cnf_pm(l1,tp,sapbs*sambs*capbs*cambs*(rs?-1:1),
					&bcnf,&bfit,ml);
			
			m2_1=MakeCompound(A_MTERM,3);
			SetCompoundArg(m2_1,1,NewInteger(1));
			SetCompoundArg(m2_1,2,MakeList2(
					MakeCompound2(OPR_POW,tp<3?s1:c1,NewInteger(1)),
					MakeCompound2(OPR_POW,        c2,NewInteger(3))));
			
			m2_2=MakeCompound(A_MTERM,3);
			SetCompoundArg(m2_2,1,NewInteger((tp%2)?1:-1));
			SetCompoundArg(m2_2,2,MakeList2(
					MakeCompound2(OPR_POW,tp<3?c1:s1,NewInteger(1)),
					MakeCompound2(OPR_POW,        s2,NewInteger(3))));
			
			l2=MakeList2(m2_1,m2_2);
			
			m2=MakeCompound(A_MTERM,3);
			SetCompoundArg(m2,1,NewInteger(bcnf*4));
			SetCompoundArg(m2,2,MakeList1(
					MakeCompound2(OPR_POW,tri_dummy_prm(l2),NewInteger(2))));
			
			*ml=AppendLast(*ml,m2);
		}
		
		if(zl[0]) RemoveList(zl[0]);
		if(zl[1]) RemoveList(zl[1]);
		if(zl[2]) RemoveList(zl[2]);
		if(zl[3]) RemoveList(zl[3]);
		
		if(bfit<5)
			return 1;
		
		return 0;

	}
	else
	{
		int tp=0, rs=0, bfit=10, bcnf=0;
		List zl[4];
		zl[0]=zl[1]=zl[2]=zl[3]=0;
		
		if((l1=find_3_mono_ab(*ml,s1,s2,c1,c2,3)))
		{
			
			find_3_cnf_ab(l1,3,&cnf,&fno,0);
			if(fno<bfit)
				tp=3,rs=0,bfit=fno,bcnf=cnf;
			find_3_cnf_ab(l1,4,&cnf,&fno,0);
			if(fno<bfit)
				tp=4,rs=0,bfit=fno,bcnf=cnf;
			zl[2]=l1;
		}
		if((l1=find_3_mono_ab(*ml,s2,s1,c2,c1,3)))
		{
			
			find_3_cnf_ab(l1,3,&cnf,&fno,0);
			if(fno<bfit)
				tp=3,rs=1,bfit=fno,bcnf=cnf;
			find_3_cnf_ab(l1,4,&cnf,&fno,0);
			if(fno<bfit)
				tp=4,rs=1,bfit=fno,bcnf=cnf;
			zl[3]=l1;
		}

		if((l1=find_3_mono_ab(*ml,s1,s2,c1,c2,1)))
		{
			
			find_3_cnf_ab(l1,1,&cnf,&fno,0);
			if(fno<bfit)
				tp=1,rs=0,bfit=fno,bcnf=cnf;
			find_3_cnf_ab(l1,2,&cnf,&fno,0);
			if(fno<bfit)
				tp=2,rs=0,bfit=fno,bcnf=cnf;
			zl[0]=l1;
		}
		if((l1=find_3_mono_ab(*ml,s2,s1,c2,c1,1)))
		{
			
			find_3_cnf_ab(l1,1,&cnf,&fno,0);
			if(fno<bfit)
				tp=1,rs=1,bfit=fno,bcnf=cnf;
			find_3_cnf_ab(l1,2,&cnf,&fno,0);
			if(fno<bfit)
				tp=2,rs=1,bfit=fno,bcnf=cnf;
			zl[1]=l1;
		}
		
				
		if(bfit<2)
		{
			Term m2,m2_1,m2_2;
/*			
printf("opti_3: type=%d, revers=%d, best fit=%d\n",tp,rs,bfit);			
*/			
			if(rs)
			{
				Atom tmp;
				tmp=s1;
				s1=s2;
				s2=tmp;
				tmp=c1;
				c1=c2;
				c2=tmp;
			}
			
			if(tp<3)
				l1=zl[rs];
			else
				l1=zl[rs+2];
			
			find_3_cnf_ab(l1,tp,&bcnf,&bfit,ml);
			
			m2_1=MakeCompound(A_MTERM,3);
			SetCompoundArg(m2_1,1,NewInteger(1));
			SetCompoundArg(m2_1,2,MakeList2(
					MakeCompound2(OPR_POW,tp<3?s1:c1,NewInteger(1)),
					MakeCompound2(OPR_POW,        c2,NewInteger(3))));
			
			m2_2=MakeCompound(A_MTERM,3);
			SetCompoundArg(m2_2,1,NewInteger((tp%2)?1:-1));
			SetCompoundArg(m2_2,2,MakeList2(
					MakeCompound2(OPR_POW,tp<3?c1:s1,NewInteger(1)),
					MakeCompound2(OPR_POW,        s2,NewInteger(3))));
			
			l2=MakeList2(m2_1,m2_2);
			
			m2=MakeCompound(A_MTERM,3);
			SetCompoundArg(m2,1,NewInteger(bcnf));
			SetCompoundArg(m2,2,MakeList1(
					MakeCompound2(OPR_POW,tri_dummy_prm(l2),NewInteger(2))));
			
			*ml=AppendLast(*ml,m2);
		}
		
		if(zl[0]) RemoveList(zl[0]);
		if(zl[1]) RemoveList(zl[1]);
		if(zl[2]) RemoveList(zl[2]);
		if(zl[3]) RemoveList(zl[3]);
		
		if(bfit<2)
			return 1;
		
		return 0;

	}
	
	
	return 0;
}
