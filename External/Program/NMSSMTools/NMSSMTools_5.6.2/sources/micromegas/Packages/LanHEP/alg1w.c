#include "lanhep.h"


static Term ListToCaret(List li)
	{
	if(ListLength(li)==1)
		return ListFirst(li);
	return MakeCompound2(OPR_CARET,ListFirst(li),ListToCaret(ListTail(li)));
	}

static Term ListToComma(List li)
	{
	if(ListLength(li)==1)
		return ListFirst(li);
	return MakeCompound2(OPR_COMMA,ListFirst(li),ListToComma(ListTail(li)));
	}

static Term s_to_expr(Term t)
	{
	List l1,l2,l3;	
	if(CompoundName(t)==OPR_WILD)
		{
		l3=l1=ConsumeCompoundArg(t,2);
		l2=NewList();
		if(l1==0)
			return NewInteger(0);
		while(!is_empty_list(l1))
			{
			l2=AppendLast(l2,Alg1ToExpr(ListFirst(l1)));
			l1=ListTail(l1);
			}
		RemoveList(l3);
		l2=ListToComma(l2);
		l2=MakeCompound1(A_FBRACET,l2);
		l3=l1=ConsumeCompoundArg(t,1);
		l1=CompoundArg2(ListNth(l1,ListLength(l1)));
		FreeAtomic(l3);
		FreeAtomic(t);
		l2=MakeCompound2(OPR_CARET,l2,l1);
		return l2;
		}
	l1=ConsumeCompoundArg(t,1);
	l2=AppendLast(NewList(),ConsumeCompoundArg(t,2));
	FreeAtomic(t);
	while(!is_empty_list(l1))
		{
		l2=AppendLast(l2,CompoundArg2(ListFirst(l1)));
		l1=ListTail(l1);
		}
	l1=ListToCaret(l2);
	RemoveList(l2);
	return l1;
	}


static Term l_to_expr(Term t)
	{
	List li;
	Term t1;
	t1=t;
	if(is_empty_list(t))
		return 0;
	li=NewList();
	while(!is_empty_list(t))
		{
		li=AppendFirst(li,s_to_expr(ListFirst(t)));
		t=ListTail(t);
		}
	RemoveList(t1);	
	t1= l2mult(li);
	RemoveList(li);
	return t1;
	}
	


static Term m_to_expr(Term t)
	{
	Integer num,den;
	int img=0;
	
	Term snum,sden;
	num=IntegerValue(ConsumeCompoundArg(t,1));
	den=IntegerValue(ConsumeCompoundArg(t,2));
	if(den<0)
		{
		img=1;
		den=-den;
		}
	snum=l_to_expr(ConsumeCompoundArg(t,3));
	sden=l_to_expr(ConsumeCompoundArg(t,4));
	FreeAtomic(t);
	if(img)
		{
		if(snum)
			snum=MakeCompound2(OPR_MLT,snum,A_I);
		else
			snum=A_I;
		}
	if(den!=1)
		{
		if(sden)
			sden=MakeCompound2(OPR_MLT,sden,NewInteger(den));
		else
			sden=NewInteger(den);
		}
	if(num!=1)
		{
		if(snum)
			snum=MakeCompound2(OPR_MLT,snum,NewInteger(num));
		else
			snum=NewInteger(num);
		}
	if(snum && sden)
		return MakeCompound2(OPR_DIV,snum,sden);
	if(snum)
		return snum;
	if(sden)
		return  MakeCompound2(OPR_DIV,NewInteger(1),sden);
	return NewInteger(1);
	}
	


Term Alg1ToExpr(Term t)
	{
	List l,l1,l2;
	if(t==0)
		return 0;
	l1=l=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	l2=NewList();
	if(is_empty_list(l))
		return NewInteger(0);
	while(!is_empty_list(l))
		{
		l2=AppendFirst(l2,m_to_expr(ListFirst(l)));
		l=ListTail(l);
		}
	RemoveList(l1);
	l1=l2plus(l2);
	RemoveList(l2);
	return l1;
	}
	
	
	
	
