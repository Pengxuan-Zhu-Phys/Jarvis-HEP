#include "lanhep.h"

long int gcf(long int i1, long  int i2)
	{
	if(i1==0 || i2==0)
		return 1;
	if(i1<0)
		i1=-i1;
	if(i2<0)
		i2=-i2;
	while(i1!=i2)
		if(i1>i2)
			i1-=i2;
		else
			i2-=i1;
	return i1;
	}


List CommaToList(Term t)
	{
	List li;
	li=NewList();
	while(is_compound(t) && CompoundName(t)==OPR_COMMA)
		{
		Term t1;
		li=AppendLast(li,ConsumeCompoundArg(t,1));
		t1=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		t=t1;
		}
	li=AppendLast(li,t);
	return li;
	}
	
	
List OperToList(Term t, Atom opr)
	{
	List li;
	li=NewList();
	while(is_compound(t) && CompoundName(t)==opr)
		{
		Term t1;
		li=AppendLast(li,ConsumeCompoundArg(t,1));
		t1=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		t=t1;
		}
	li=AppendLast(li,t);
	return li;
	}
	
List Oper1ToList(Term t, Atom opr)
	{
	List li;
	li=NewList();
	while(is_compound(t) && CompoundName(t)==opr)
		{
		Term t1;
		li=AppendLast(li,ConsumeCompoundArg(t,2));
		t1=ConsumeCompoundArg(t,1);
		FreeAtomic(t);
		t=t1;
		}
	li=AppendLast(li,t);
	return li;
	}

int err_cnt=0;

void ErrorInfo(int errno)
	{
	if(err_cnt++ > 10)
	{
	printf("more errors follow, exiting...\n");
	exit(1);
	}
	if(!IsTermInput())
			printf("File \"%s\", line %d: ",CurrentInputFile(),
				CurrentInputLine());
	printf("Error: ");
	}

void WarningInfo(int errno)
	{
	if(!IsTermInput())
			printf("File \"%s\", line %d: ",CurrentInputFile(),
				CurrentInputLine());
	printf("Warning: ");
	}
	
static List ait_l;
	
static void set_atom(Atom a)
	{
	List li;
	li=ait_l;
	while(!is_empty_list(li))
		{
		if(ListFirst(li)==a)
			return;
		li=ListTail(li);
		}
	ait_l=AppendFirst(ait_l,a);
	return;
	}
	
static void atoms_in_term(Term t)
	{
	if(is_atom(t))
		{
		set_atom(t);
		return;
		}
	if(is_compound(t))
		{
		int i,ac;
		ac=CompoundArity(t);
		for(i=1;i<=ac;i++)
			atoms_in_term(CompoundArgN(t,i));
		return;
		}
	if(is_list(t))
		{
		while(!is_empty_list(t))
			{
			atoms_in_term(ListFirst(t));
			t=ListTail(t);
			}
		return;
		}
	return;
	}
	
List AtomsInTerm(Term t, Term t1)
	{
	ait_l=NewList();
	atoms_in_term(t);
	return ait_l;
	}
		
Term l2plus(List li)
	{
	int le;
	le=ListLength(li);
	if(le==1)
		return ListFirst(li);
	return MakeCompound2(OPR_PLUS,l2plus(ListTail(li)),ListFirst(li));
	}
	
Term l2mult(List li)
	{
	int le;
	if(is_empty_list(li))
		return NewInteger(1);
	le=ListLength(li);
	if(le==1)
		return ListFirst(li);
	return MakeCompound2(OPR_MLT,l2mult(ListTail(li)),ListFirst(li));
	}


void WriteBlank(FILE *f, int space)
	{
	int i;
	for(i=0;i<space;i++)
		fprintf(f," ");
	}
	
	
static char *rllist[600];
static int rllen=0;
	
void RegisterLine(char *str)
	{
	rllist[rllen]=str;
	rllen++;
	if(rllen>599) rllen--;
/*	printf("Reg  line: %s\n",str);*/
	return;
	}

void UnregisterLine(void)
	{
/*	printf("Unreg line: %s\n",rllist[rllen-1]);*/
	if(rllen) rllen--;
	return;
	}

void DumpRegistered(void)
	{
	int i;
	for(i=0;i<rllen;i++)
		puts(rllist[i]);
	return;
	}

void ReportRedefined(Atom p, char *nname)
	{
	if(is_particle(p,NULL))
		{
		if(!IsTermInput())
			printf("File \"%s\", line %d: ",CurrentInputFile(),
				CurrentInputLine());
		printf("warning:\n");
		printf("redefinition of %s; it was a particle, becomes a %s.\n",
			AtomValue(p),nname);
		ClearParticle(p);
		return;
		}
		
	if(is_parameter(p))
		{
		if(!IsTermInput())
			printf("File \"%s\", line %d: ",CurrentInputFile(),
				CurrentInputLine());
		printf("warning:\n");
		printf("redefinition of %s; it was a parameter, becomes a %s.\n",
			AtomValue(p),nname);
		ClearParameter(p);
		return;
		}
		
	if(is_let(p,NULL))
		{
		if(!IsTermInput())
			printf("File \"%s\", line %d: ",CurrentInputFile(),
				CurrentInputLine());
		printf("warning:\n");
		printf("Redefinition of %s; it was a let-substitution, becomes a %s.\n",
			AtomValue(p),nname);
		ClearLet(p);
		return;
		}
		
	if(is_special(p,NULL))
		{
		if(!IsTermInput())
			printf("File \"%s\", line %d: ",CurrentInputFile(),
				CurrentInputLine());
		printf("warning:\n");
		printf("Redefinition of %s; it was a special, becomes a %s.\n",
			AtomValue(p),nname);
		ClearSpecial(p);
		return;
		}
		
		
	}
	
void WriteVertex(List pl)
	{
	printf("(");
	while(!is_empty_list(pl))
		{
		printf("%s",AtomValue(CompoundArg1(ListFirst(pl))));
		pl=ListTail(pl);
		if(pl)
			printf(", ");
		}
	printf(")");
	}
	


void ReplAtom(Term t, Atom from, Atom to)
	{
	if(is_compound(t))
		{
		int i, ar;
		ar=CompoundArity(t);
		for(i=1;i<=ar;i++)
			{
			if(CompoundArgN(t,i)==from)
				SetCompoundArg(t,i,to);
			else
				ReplAtom(CompoundArgN(t,i),from,to);
			}
		return;
		}
	if(is_list(t))
		{
		while(!is_empty_list(t))
			{
			if(ListFirst(t)==from)
				ChangeList(t,to);
			else
				ReplAtom(ListFirst(t),from,to);
			t=ListTail(t);
			}
		}
	return;
	}
	
List XorList(List l, Atomic a)
	{
	List l1;
	l1=l;
	while(!is_empty_list(l1))
		{
		if(ListFirst(l1)==a)
			return CutFromList(l,l1);
		l1=ListTail(l1);
		}
	return AppendLast(l,a);
	}

	
Term il_to_caret(Term t, List il)
{
	Term t1, t2;
	List l;
	if(is_empty_list(il))
		return t;
	
	t1=MakeCompound2(OPR_CARET,t,0);
	
	t2=t1;
	
	for(l=il;ListTail(l);l=ListTail(l))
	{
		SetCompoundArg(t2,2,MakeCompound2(OPR_CARET,ListFirst(l),0));
		t2=CompoundArg2(t2);
	}
	
	SetCompoundArg(t2,2,ListFirst(l));
	
	return t1;
}

	
