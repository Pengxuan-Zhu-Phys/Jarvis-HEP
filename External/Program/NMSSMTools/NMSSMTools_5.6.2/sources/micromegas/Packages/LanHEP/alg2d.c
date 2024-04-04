#include <string.h>
#include "lanhep.h"

extern int TexOutput;		

static int prtcmp(Term p1, Term p2)
	{
	if(CompoundArg1(p1)==A_I)
		return -1;
	if(CompoundArg1(p2)==A_I)
		return 1;
	return strcmp(AtomValue(CompoundArg1(p1)),AtomValue(CompoundArg1(p2)));
	}


void alg2_common_t(Term a2)
	{
	List ct, l1,l2,l3,l4, ml, gct, vil=0;
	ct=NewList();

	ml=CompoundArgN(a2,5);
	if(is_empty_list(ml) || (ListLength(ml)==1 && FAOutput==0) )
		return;
	gct=CompoundArgN(ListFirst(ml),3);
	if(is_empty_list(gct))
		return;

	if(FAOutput)
	{
		for(l1=CompoundArg1(a2);l1;l1=ListTail(l1))
			if(CompoundName(CompoundArg2(ListFirst(l1)))==OPR_VECTOR ||
					CompoundName(CompoundArg2(ListFirst(l1)))==OPR_SPINOR)
				vil=AppendLast(vil,ListFirst(CompoundArg1(
						CompoundArg2(ListFirst(l1)))));
	}
	
	l1=gct;
	while(!is_empty_list(l1))
		{
		Term spt;
		spt=ListFirst(l1);
		if(FAOutput)
		{
			Term pr;
			if(CompoundName(spt)==A_MOMENT)
				goto cnt1;
				
			if((is_atom(CompoundArg1(spt))&& 
				  (pr=GetAtomProperty(CompoundArg1(spt),A_INFINITESIMAL))))
				  goto cnt1;
				  
			if(CompoundName(spt)==OPR_SPECIAL)
			{
				l2=CompoundArg1(spt);
				if(l2==A_GAMMA||l2==A_GAMMA5||l2==A_GAMMAP||l2==A_GAMMAM)
					goto cnt1;
			}
			if(vil && CompoundName(spt)==OPR_SPECIAL && CompoundArg1(spt)==A_DELTA)
			{
				l2=CompoundArg2(spt);
				if(ListMember(vil,ListFirst(l2)))
					goto cnt1;
				if(ListMember(vil,ListFirst(ListTail(l2))))
					goto cnt1;
			}
		}
		l2=ListTail(ml);
		while(!is_empty_list(l2))
			{
			l3=CompoundArgN(ListFirst(l2),3);
			while(!is_empty_list(l3))
				{
				if(EqualTerms(spt,ListFirst(l3)))
					break;
				l3=ListTail(l3);
				}
			if(is_empty_list(l3))
				goto cnt1;
			l2=ListTail(l2);
			}
		ct=AppendLast(ct,CopyTerm(spt));
	cnt1:
		l1=ListTail(l1);
		}

	if(is_empty_list(ct))
		return;

	l1=ml;
	while(!is_empty_list(l1))
		{
		Term m2;
		m2=ListFirst(l1);
		l2=ConsumeCompoundArg(m2,3);
		l3=ct;
		while(!is_empty_list(l3))
			{
			Term rt;
			rt=ListFirst(l3);
			l4=l2;
			while(!is_empty_list(l4))
				{
				if(EqualTerms(ListFirst(l4),rt))
					{
					l2=CutFromList(l2,l4);
					goto cnt2;
					}
				l4=ListTail(l4);
				}
			puts("Internal error (a2ctf).");
		cnt2:
			l3=ListTail(l3);
			}
		SetCompoundArg(m2,3,l2);
		l1=ListTail(l1);
		}
	SetCompoundArg(a2,4,ct);

	}

int allow_sym_div=0;

static List p_m_p_s(List cm, List mo)
	{
	List l,l1;
	if(allow_sym_div)
		{
		for(l=cm;l;l=ListTail(l))
			{
			Term t=ListFirst(l),t1;
			int ic,im;
			ic=(int)IntegerValue(CompoundArg2(t));
			if(ic==0) continue;
			for(l1=mo;l1;l1=ListTail(l1))
				if(CompoundArg1(ListFirst(l1))==CompoundArg1(t))
					break;
			if(l1==0)
				{
				SetCompoundArg(t,2,NewInteger(0));
				continue;
				}
			t1=ListFirst(l1);
			im=(int)IntegerValue(CompoundArg2(t1));
			if(ic*im<0)
				{
				SetCompoundArg(t,2,NewInteger(0));
				continue;
				}
			if(ic>0 && im>0 && im<ic)
				{
				SetCompoundArg(t,2,NewInteger(im));
				continue;
				}
			if(ic<0 && im<0 && im>ic)
				{
				SetCompoundArg(t,2,NewInteger(im));
				continue;
				}
			}
		return cm;
		}
			
	l=mo;
	while(!is_empty_list(l))
		{
		Term t;
		t=ListFirst(l);
		l1=cm;
		while(!is_empty_list(l1))
			{
			Term t1;
			t1=ListFirst(l1);
			if(CompoundArg1(t1)==CompoundArg1(t))
				{
				if(IntegerValue(CompoundArg2(t))<IntegerValue(CompoundArg2(t1)))
					SetCompoundArg(t1,2,CompoundArg2(t));
				goto lcnt;
				}
			l1=ListTail(l1);
			}
		t=CopyTerm(t);
		if(IntegerValue(CompoundArg2(t))>0)
			SetCompoundArg(t,2,NewInteger(0));
		cm=AppendLast(cm,t);
	lcnt:
		l=ListTail(l);
		}
	l=cm;
	while(!is_empty_list(l))
		{
		Term t;
		t=ListFirst(l);
		l1=mo;
		while(!is_empty_list(l1))
			{
			Term t1;
			t1=ListFirst(l1);
			if(CompoundArg1(t1)==CompoundArg1(t))
				goto lcnt1;
			l1=ListTail(l1);
			}
		if(IntegerValue(CompoundArg2(t))>0)
			SetCompoundArg(t,2,NewInteger(0));
	lcnt1:
		l=ListTail(l);
		}
	return cm;
	}


void alg2_common_s(Term a2)
	{
	List m2l,l,la,l1;
	m2l=CompoundArgN(a2,5);
	l=m2l;
	if(is_empty_list(l))
		return;
	la=CopyTerm(CompoundArg2(ListFirst(l)));
	l=ListTail(l);
	while(!is_empty_list(l))
		{
		la=p_m_p_s(la,CompoundArg2(ListFirst(l)));
		l=ListTail(l);
		}
	l=la;
	while(!is_empty_list(l))
		{
		if(IntegerValue(CompoundArg2(ListFirst(l)))==0 || 
				(! (TexOutput||FAOutput) && 
				GetAtomProperty(CompoundArg1(ListFirst(l)),A_DUMMY_PRM)) ||
				(FAOutput && GetAtomProperty(CompoundArg1(ListFirst(l)),
					A_INFINITESIMAL)) )
			{
			la=CutFromList(la,l);
			l=la;
			}
		else
			l=ListTail(l);
		}
	SetCompoundArg(a2,3,la);
	l=m2l;
	while(!is_empty_list(l))
		{
		Term t;
		List l2,l3;
		t=ListFirst(l);
		l1=ConsumeCompoundArg(t,2);
		l2=la;
		while(!is_empty_list(l2))
			{
			Term t2;
			t2=ListFirst(l2);
			l3=l1;
			while(!is_empty_list(l3))
				{
				Term t3;
				t3=ListFirst(l3);
				if(CompoundArg1(t2)==CompoundArg1(t3))
					{
					int np;
					np=(int)IntegerValue(CompoundArg2(t3))-(int)IntegerValue(CompoundArg2(t2));
					if(np)
						{
						SetCompoundArg(t3,2,NewInteger(np));
						goto l2cnt;
						}
					else
						{
						l1=CutFromList(l1,l3);
						goto l2cnt;
						}
					}
				l3=ListTail(l3);
				}
			l1=AppendLast(l1,MakeCompound2(OPR_POW,CompoundArg1(t2),
						NewInteger(-IntegerValue(CompoundArg2(t2)))));
			/*puts("Internal error: csm failed");*/
		l2cnt:
			l2=ListTail(l2);
			}
		l1=SortedList(l1,prtcmp);
		SetCompoundArg(t,2,l1);
		l=ListTail(l);
		}

	}


void alg2_recommon_s(Term a2)
	{
	List m2l,l,la,l1;
	m2l=CompoundArgN(a2,5);
	l=m2l;
	if(is_empty_list(l))
		return;
	la=CopyTerm(CompoundArg2(ListFirst(l)));
	l=ListTail(l);
	while(!is_empty_list(l))
		{
		la=p_m_p_s(la,CompoundArg2(ListFirst(l)));
		l=ListTail(l);
		}
	l=la;
	while(!is_empty_list(l))
		{
		if(IntegerValue(CompoundArg2(ListFirst(l)))==0 || 
				(!(TexOutput||FAOutput) &&
				GetAtomProperty(CompoundArg1(ListFirst(l)),A_DUMMY_PRM)) ||
				(FAOutput && GetAtomProperty(CompoundArg1(ListFirst(l)),
					A_INFINITESIMAL)) )
			{
			la=CutFromList(la,l);
			l=la;
			}
		else
			l=ListTail(l);
		}
	if(is_empty_list(la))
		return;
    l1=ConsumeCompoundArg(a2,3);
    for(l=la;l;l=ListTail(l))
        {
        List l2;
        for(l2=l1;l2;l2=ListTail(l2))
            {
			if(CompoundArg1(ListFirst(l2))==CompoundArg1(ListFirst(l)))
				{
				long int b1, b2;
				b1=IntegerValue(CompoundArg2(ListFirst(l)));
				b2=IntegerValue(CompoundArg2(ListFirst(l2)));
				b1+=b2;
				if(b1==0)
					l1=CutFromList(l1,l2);
				else
					SetCompoundArg(ListFirst(l2),2,NewInteger(b1));
				break;
				}
			}
		if(is_empty_list(l2))
			l1=AppendFirst(l1,CopyTerm(ListFirst(l)));
		}

	l1=SortedList(l1,prtcmp);
	SetCompoundArg(a2,3,l1);

	l=m2l;
	while(!is_empty_list(l))
		{
		Term t;
		List l2,l3;
		t=ListFirst(l);
		l1=ConsumeCompoundArg(t,2);
		l2=la;
		while(!is_empty_list(l2))
			{
			Term t2;
			t2=ListFirst(l2);
			l3=l1;
			while(!is_empty_list(l3))
				{
				Term t3;
				t3=ListFirst(l3);
				if(CompoundArg1(t2)==CompoundArg1(t3))
					{
					long int np;
					np=IntegerValue(CompoundArg2(t3))-IntegerValue(CompoundArg2(t2));
					if(np)
						{
						SetCompoundArg(t3,2,NewInteger(np));
						goto l2cnt;
						}
					else
						{
						l1=CutFromList(l1,l3);
						goto l2cnt;
						}
					}
				l3=ListTail(l3);
				}
			l1=AppendLast(l1,MakeCompound2(OPR_POW,CompoundArg1(t2),
						NewInteger(-IntegerValue(CompoundArg2(t2)))));
			/*puts("Internal error: csm failed");*/
		l2cnt:
			l2=ListTail(l2);
			}
		l1=SortedList(l1,prtcmp);
		SetCompoundArg(t,2,l1);
		l=ListTail(l);
		}

	FreeAtomic(la);
	}


void alg2_decommon_s(Term a2)
	{
	Term cfl;
	List l;
	
	cfl=ConsumeCompoundArg(a2,3);
	if(is_empty_list(cfl))
		return;
	
	for(l=CompoundArgN(a2,5);l;l=ListTail(l))
		{
		List pfl;
		List l1,l2;
		
		pfl=ConsumeCompoundArg(ListFirst(l),2);
		for(l1=cfl;l1;l1=ListTail(l1))
			{
			Atom p;
			p=CompoundArg1(ListFirst(l1));
			for(l2=pfl;l2;l2=ListTail(l2))
				{
				if(CompoundArg1(ListFirst(l2))==p)
					{
					int pw;
					pw=(int)IntegerValue(CompoundArg2(ListFirst(l2)))
						+(int)IntegerValue(CompoundArg2(ListFirst(l1)));
					SetCompoundArg(ListFirst(l2),2,NewInteger(pw));
					break;
					}
				}
			if(is_empty_list(l2))
				pfl=AppendFirst(pfl,CopyTerm(ListFirst(l1)));
			}
	rr:
		for(l1=pfl;l1;l1=ListTail(l1))
			if(CompoundArg2(ListFirst(l1))==NewInteger(0))
				{
				pfl=CutFromList(pfl,l1);
				goto rr;
				}
		pfl=SortedList(pfl,prtcmp);
		SetCompoundArg(ListFirst(l),2,pfl);
		}
		
	FreeAtomic(cfl);
		
	}

static int mlt_list1(List l)
{
  int ret=1;
  for(;l;l=ListTail(l))
    {
      int i=(int)IntegerValue(ListFirst(l));
      i/=(int)gcf(i,ret);
      ret*=i;
    }
  return ret;
}

/*
static int mlt_list1(List l)
	{
	int ret=1;
	while(!is_empty_list(l))
		{
		if(	is_empty_list(ListTail(l)) ||
			!ListMember(ListTail(l),ListFirst(l)))
			ret*=IntegerValue(ListFirst(l));
		l=ListTail(l);
		}
	return ret;
	}
*/
static int gcf_list(List l)
	{
	int ret;
	if(is_empty_list(l))
		return 1;
	ret=(int)IntegerValue(ListFirst(l));
    if(ret<0)
        ret=-ret;
	l=ListTail(l);
	while(!is_empty_list(l))
		{
		if(ret==1)
			return 1;
		ret=(int)gcf(ret,IntegerValue(ListFirst(l)));
		l=ListTail(l);
		}
	return ret;
	}

void alg2_common_n(Term a2)
	{
	List m2l,l,l1,nl,dl;
	int cnum,cden,ti;
	m2l=CompoundArgN(a2,5);
	nl=dl=NewList();
	l=m2l;
	if(is_empty_list(l))
		return;


	while(!is_empty_list(l))
		{
		Term t;
		t=ConsumeCompoundArg(ListFirst(l),1);
		nl=AppendLast(nl,CompoundArg1(t));
		dl=AppendLast(dl,CompoundArg2(t));
		FreeAtomic(t);
		l=ListTail(l);
		}

	cnum=gcf_list(nl);
	cden=gcf_list(dl);
	if(IntegerValue(ListFirst(nl))<0)
		cnum*=-1;
	
	l=nl;
	l1=dl;
	while(!is_empty_list(l))
		{
		ChangeList(l,NewInteger(IntegerValue(ListFirst(l))/cnum));
		ChangeList(l1,NewInteger(IntegerValue(ListFirst(l1))/cden));
		l=ListTail(l);
		l1=ListTail(l1);
		}

	ti=mlt_list1(dl);

	l=nl;
	l1=dl;
	while(!is_empty_list(l))
		{
		ChangeList(l,
			NewInteger(IntegerValue(ListFirst(l))*ti/IntegerValue(ListFirst(l1))));
		l=ListTail(l);
		l1=ListTail(l1);
		}

	cden*=ti;
	ti=gcf_list(nl);
	cnum*=ti;

	l=nl;
	while(!is_empty_list(l))
		{
		ChangeList(l,NewInteger(IntegerValue(ListFirst(l))/ti));
		l=ListTail(l);
		}

	ti=(int)gcf(cnum,cden);
	cnum/=ti;
	cden/=ti;
	RemoveList(dl);

	SetCompoundArg(a2,2,MakeCompound2(OPR_DIV,NewInteger(cnum),NewInteger(cden)));
	l=m2l;
	l1=nl;

	while(!is_empty_list(l))
		{
		SetCompoundArg(ListFirst(l),1,ListFirst(l1));
		l1=ListTail(l1);
		l=ListTail(l);
		}
	RemoveList(nl);

	return ;
	}



void alg2_decommon_n(Term a2)
	{
    List l1;
	long int cnum,cden;
	Term t;

	l1=CompoundArgN(a2,5);
	if(is_empty_list(l1))
		{
		SetCompoundArg(a2,2,0);
		return;
		}

	t=ConsumeCompoundArg(a2,2);
	cnum=IntegerValue(CompoundArg1(t));
	cden=IntegerValue(CompoundArg2(t));
	FreeAtomic(t);

	while(!is_empty_list(l1))
		{
		Term t;
		long int c1,c2,c3;
		c1=cnum*IntegerValue(CompoundArg1(ListFirst(l1)));
		c2=cden;
		c3=gcf(c1,c2);
		c1/=c3;
		c2/=c3;
		t=MakeCompound2(OPR_DIV,NewInteger(c1),NewInteger(c2));
		SetCompoundArg(ListFirst(l1),1,t);
		l1=ListTail(l1);
		}


	}


void alg2_recommon_n(Term a2)
	{
    List m2l,l,l1,nl;
    long int cnum,n,d;
	m2l=CompoundArgN(a2,5);
    nl=NewList();
	l=m2l;
	if(is_empty_list(l))
		return;


	while(!is_empty_list(l))
		{
        nl=AppendLast(nl,CompoundArg1(ListFirst(l)));
		l=ListTail(l);
		}

	cnum=gcf_list(nl);
	if(IntegerValue(ListFirst(nl))<0)
		cnum*=-1;
    if(cnum==1)
        {
        RemoveList(nl);
        return;
        }

    l1=m2l;
	l=nl;
	while(!is_empty_list(l))
		{
        SetCompoundArg(ListFirst(l1),1,
            NewInteger(IntegerValue(ListFirst(l))/cnum));
		l=ListTail(l);
		l1=ListTail(l1);
		}

	RemoveList(nl);
    n=IntegerValue(CompoundArg1(CompoundArg2(a2)));
    d=IntegerValue(CompoundArg2(CompoundArg2(a2)));
    n*=cnum;
    cnum=gcf(n,d);
    n/=cnum;
    d/=cnum;
    SetCompoundArg(CompoundArg2(a2),1,NewInteger(n));
    SetCompoundArg(CompoundArg2(a2),2,NewInteger(d));

	return ;
	}

