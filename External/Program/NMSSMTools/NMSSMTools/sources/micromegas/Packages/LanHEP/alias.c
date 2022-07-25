#include <string.h>
#include "lanhep.h"

static List inst_ali=0;


static Term tolab(Term t, List li)
	{
	if(is_atom(t))
		{
		List l;
		l=li;
		while(!is_empty_list(l))
			{
			if(CompoundArg1(ListFirst(l))==t)
				return CompoundArg2(ListFirst(l));
			l=ListTail(l);
			}
		return t;
		}
	if(is_compound(t))
		{
		Term t1;
		int i,ar;
		t1=NewCompound(CompoundFunctor(t));
		ar=CompoundArity(t);
		for(i=1;i<=ar;i++)
			SetCompoundArg(t1,i,tolab(ConsumeCompoundArg(t,i),li));
		FreeAtomic(t);
		return t1;
		}
	if(is_list(t))
		{
		List l,l1;
		l=NewList();
		l1=t;
		while(!is_empty_list(l1))
			{
			l=AppendLast(l,tolab(ListFirst(l1),li));
			l1=ListTail(l1);
			}
		RemoveList(t);
		return l;
		}
	return t;
	}

static List plused_atoms=0;

static Term del_plus(Term t)
{
  if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==OPR_PLUS)
  {
	Term t1=ConsumeCompoundArg(t,1);
	if(!is_atom(t1))
	{
	  ErrorInfo(1);
	  printf("alias: argument of '+' must be a symbol, got ");
	  WriteTerm(t1);puts("");
	}
	else
	  plused_atoms=AppendLast(plused_atoms,t1);
	
	FreeAtomic(t);
	return t1;
  }
  
  if(is_compound(t))
  {
	int i, a=CompoundArity(t);
	for(i=1;i<=a;i++)
	{
	  Term t1=ConsumeCompoundArg(t,i);
	  SetCompoundArg(t,i,del_plus(t1));
	}
	return t;
  }
	
  if(is_list(t))
  {
	List l;
	for(l=t;l;l=ListTail(l))
	{
	  ChangeList(l,del_plus(ListFirst(l)));
	}
	return t;
  }
  
  return t;
  
}

Term SetAlias(Term al, Term sub, int ver)
	{
	List l1,l2,l3,l4;
	Term tt, alsv;
	
	alsv=CopyTerm(al);
	
	plused_atoms=0;
	al=del_plus(al);
	
	l1=AtomsInTerm(al,0);
	l2=AtomsInTerm(sub,0);
	
	
	for(l3=l1;l3;l3=ListTail(l3))
		if(AtomValue(ListFirst(l3))[0]=='_' && !ListMember(l2,ListFirst(l3)))
			l2=AppendLast(l2,ListFirst(l3));
	l3=IntersectList(l1,l2);

	
	FreeAtomic(l1);
	FreeAtomic(l2);
	
	l1=0;
	
	for(l2=l3;l2;l2=ListTail(l2))
	  if(!ListMember(plused_atoms,ListFirst(l2)))
		l1=AppendLast(l1,ListFirst(l2));
	  
	RemoveList(l3);
	RemoveList(plused_atoms);
	plused_atoms=0;
	l3=l1;
	
	while(!is_empty_list(l1))
		{
		ChangeList(l1,MakeCompound2(OPR_DIV,ListFirst(l1),NewLabel()));
		l1=ListTail(l1);
		}
	al=tolab(al,l3);
	sub=tolab(sub,l3);
	if(is_label(al))
		{
		if(ver)
			{
			ErrorInfo(211);
			printf("illegal alias\n");
			}
		FreeAtomic(l3); FreeAtomic(al); FreeAtomic(sub); return 0;
		}
	if(is_atom(al))
		l1=al;
	else
		l1=CompoundName(al);
	l2=CopyTerm(GetAtomProperty(l1,OPR_ALIAS));
	tt=MakeCompound(OPR_ALIAS,4);
	SetCompoundArg(tt,1,al);
	SetCompoundArg(tt,2,sub);
	SetCompoundArg(tt,3,l3);
	SetCompoundArg(tt,4,l4=NewLabel());
	l2=AppendLast(l2,tt);
	SetAtomProperty(l1,OPR_ALIAS,l2);
	l3=MakeCompound3(OPR_ALIAS,l1,l4,alsv);
	inst_ali=AppendLast(inst_ali,l3);
	return NewInteger(LabelValue(l4));
	}

Term InterfSetAlias(Term t,Term ind)
	{
	Term al, sub, tt;
	tt=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	if(is_compound(tt) && CompoundName(tt)==OPR_COMMA)
		{
		Term t1, t2;
		t1=ConsumeCompoundArg(tt,1);
		t2=ConsumeCompoundArg(tt,2);
		FreeAtomic(tt);
		InterfSetAlias(MakeCompound1(OPR_ALIAS,t1),ind);
		InterfSetAlias(MakeCompound1(OPR_ALIAS,t2),ind);
		return 0;
		}
	if(is_atom(tt)&&strcmp(AtomValue(tt),"?")==0)
		{
		DumpList(inst_ali);
		return 0;
		}
	if(!is_compound(tt) || CompoundArity(tt)!=2 || 
		CompoundName(tt)!=OPR_EQSIGN)
		{
		ErrorInfo(212);
		WriteTerm(tt);printf(" : wrong argument in call to SetAlias\n");
		FreeAtomic(tt);
		return 0;
		} 
	al=ConsumeCompoundArg(tt,1);
	if(!is_compound(al) && !is_atom(al))
		{
		ErrorInfo(212);
		printf("bad target for alias\n");
		FreeAtomic(al);
		return 0;
		}
	sub=ConsumeCompoundArg(tt,2);
	FreeAtomic(tt);
	if((tt=SetAlias(al,sub,1))==0)
		{
		ErrorInfo(212);
		printf("illegal alias\n");
		return 0;
		}
	return tt;
	}
	
Term do_sub(Term t, List lu)
	{
	if(is_label(t))
		{
		List l;
		l=lu;
		while(!is_empty_list(l))
			{
			Term t1;
			t1=ListFirst(l);
			if(CompoundArg1(t1)==t)
				return CopyTerm(CompoundArg2(t1));
			l=ListTail(l);
			}
		return t;
		}
	if(is_compound(t))
		{
		Term t1;
		int ac,i;
		t1=NewCompound(CompoundFunctor(t));
		ac=CompoundArity(t1);
		for(i=1;i<=ac;i++)
			SetCompoundArg(t1,i,do_sub(CompoundArgN(t,i),lu));
		return t1;
		}
	if(is_list(t))
		{
		List t1;
		t1=NewList();
		while(!is_empty_list(t))
			{
			t1=AppendLast(t1,do_sub(ListFirst(t),lu));
			t=ListTail(t);
			}
		return t1;
		}
	return CopyTerm(t);
	}
					

	
Term ali_uni(Term t, Term s, List *lu)
	{
	if(s==t)
		return 1;
	if(is_label(s))
		{
		List l;
		l=*lu;
		while(!is_empty_list(l))
			{
			Term tt;
			tt=ListFirst(l);
			if(CompoundArg1(tt)==s)
				return ali_uni(t,CompoundArg2(tt),lu);
			l=ListTail(l);
			}
		*lu=AppendFirst(*lu,MakeCompound2(OPR_DIV,s,t));
		return 1;
		}
	if(is_compound(t) && is_compound(s) && CompoundArity(s)==CompoundArity(t))
		{
		int ac,i;
		ac=CompoundArity(s);
		for(i=1;i<=ac;i++)
			if(ali_uni(CompoundArgN(t,i),CompoundArgN(s,i),lu)==0)
				return 0;
		return 1;
		}
	RemoveList(*lu);
	*lu=NewList();
	return 0;
	}
	

static int subs_in_ali;



Term proc_alias(Term t)
	{
	Atom grs;
	
	if(is_atom(t))
		{
		List prp,l;
		prp=GetAtomProperty(t,OPR_ALIAS);
		l=prp;
		while(!is_empty_list(l))
			{
			if(is_atom(CompoundArg1(ListFirst(l))))
				{
				subs_in_ali++;
				return CopyTerm(CompoundArg2(ListFirst(l)));
				}
			l=ListTail(l);
			}
		return t;
		}
	
	if(is_compound(t) && CompoundArity(t)==1 && is_atom(grs=CompoundArg1(t)) &&
			(is_particle(grs,0) || is_let(grs,0)) &&
			(GetAtomProperty(CompoundName(t),OPR_GROUP)||
			CompoundName(t)==A_LORENTZ || CompoundName(t)==OPR_WILD))
	{
		List l, lp=0, li=GetAtomProperty(grs,PROP_INDEX);
		for(l=li;l;l=ListTail(l))
		{
			List ll;
			Term sp=CompoundArg1(ListFirst(l));
			if(CompoundName(t)==A_LORENTZ && CompoundName(sp)==A_LORENTZ)
			{
				if(CompoundArg1(sp)==NewInteger(1) && CompoundArg2(sp)==NewInteger(0))
					lp=AppendLast(lp,OPR_SPINOR); else
				if(CompoundArg1(sp)==NewInteger(0) && CompoundArg2(sp)==NewInteger(1))
					lp=AppendLast(lp,A_CSPINOR); else
				if(CompoundArg1(sp)==NewInteger(2) && CompoundArg2(sp)==NewInteger(2))
					lp=AppendLast(lp,OPR_VECTOR); else
				if(CompoundArg1(sp)==NewInteger(0) && CompoundArg2(sp)==NewInteger(0))
					lp=AppendLast(lp,NewAtom("spinor2,",0)); else
				printf("Internal error: alilsp\n");
			}
			else if(CompoundName(sp)==CompoundName(t))
				lp=AppendLast(lp,CompoundArg1(sp));
		}
		if(lp==0)
			return NewInteger(1);
		if(ListLength(lp)>1)
			return lp;
		grs=ListFirst(lp);
		RemoveList(lp);
		return grs;
	}
		
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==OPR_PLUS)
	{
	  return CopyTerm(CompoundArg1(t));
	}
			
	if(is_compound(t) && GetAtomProperty(CompoundName(t),OPR_UNALIAS)==0)
		{
		List prp,l,tr;
		int ac,i;
		prp=GetAtomProperty(CompoundName(t),OPR_ALIAS);
		l=prp;
		while(!is_empty_list(l))
			{
			Term t1,lu;
			t1=ListFirst(l);
			lu=NewList();
			if(ali_uni(t,CompoundArg1(t1),&lu))
				{
				tr=do_sub(CompoundArg2(t1),lu);
				RemoveList(lu);
				subs_in_ali++;
				return tr;
				}
			l=ListTail(l);
			}
		tr=NewCompound(CompoundFunctor(t));
		ac=CompoundArity(t);
		for(i=1;i<=ac;i++)
			SetCompoundArg(tr,i,proc_alias(CompoundArgN(t,i)));
		return tr;
		} 
	
	if(is_list(t))
		{
		List tr;
		tr=NewList();
		while(!is_empty_list(t))
			{
			tr=AppendLast(tr,proc_alias(ListFirst(t)));
			t=ListTail(t);
			}
		return tr;
		}
	return CopyTerm(t);
			
	}
	
Term ProcessAlias(Term t)
	{
	 int cnt=0;
	do
		{
		Term t1;
		cnt++;
		if(cnt>15)
		{
		  ErrorInfo(100);
		  printf("possible infine loop in processing aliases.\n");
		  return 0;
		}
		subs_in_ali=0;
		t1=proc_alias(t);
		FreeAtomic(t);
		t=t1;
		}
		while(subs_in_ali);
	return t;
	}
	
Term InterfProcessAlias(Term a, Term i)
	{
	Term t1;
	t1=ConsumeCompoundArg(a,1);
	FreeAtomic(a);
	return ProcessAlias(t1);
	}


void RemoveAlias(Term lab)
	{
	List li;
	long int lv;
	if(is_compound(lab)&&CompoundName(lab)==OPR_COMMA)
		lab=CommaToList(lab);
	if(is_list(lab))
		{
		while(!is_empty_list(lab))
			{
			RemoveAlias(ListFirst(lab));
			lab=ListTail(lab);
			}
		RemoveList(lab);
		return;
		}
	if(!is_integer(lab))
		{
		for(li=inst_ali;li;li=ListTail(li))
			if(EqualTerms(lab,CompoundArgN(ListFirst(li),3)))
				break;
		if(li==0)
		{
			ErrorInfo(0);
			printf("no alias defined for ");
			WriteTerm(lab);
			puts("");
			return;
		}
		lab=NewInteger(LabelValue(CompoundArg2(ListFirst(li))));
		}
	lv=IntegerValue(lab);
	li=inst_ali;
	while(!is_empty_list(li))
		{
		Term t1;
		t1=ListFirst(li);
		if(LabelValue(CompoundArg2(t1))==lv)
			{
			Term name,pl,l1;
			name=CompoundArg1(t1);
			inst_ali=CutFromList(inst_ali,li);
			l1=pl=CopyTerm(GetAtomProperty(name,OPR_ALIAS));
			while(!is_empty_list(l1))
				{
				Term t2;
				t2=ListFirst(l1);	
				if(LabelValue(CompoundArgN(t2,4))==lv)
					{
					pl=CutFromList(pl,l1);
					if(is_empty_list(pl))
						RemoveAtomProperty(name,OPR_ALIAS);
					else
						SetAtomProperty(name,OPR_ALIAS,pl);
					return;
					}
				l1=ListTail(l1);
				}
			puts("prop lost");
			}
		li=ListTail(li);
		}
	ErrorInfo(213);
	printf(" alias #%d not found.\n",lv);
	return ;
	}
	
Term InterfRemAlias(Term t, Term ind)
	{
	Term t1;
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	RemoveAlias(t1);
	return 0;
	}

