#include <string.h>
#include "lanhep.h"


/* group types: 
 *    			0 --- U(1)
 *				1 --- SU(N+1)
 *				2 --- SO(2N)
 *              3 --- Sp(2N)
 *				4 --- SO(2N+1)
 *				5 --- E(N), N=6,7,8
 *				6 --- F(4)
 *				7 --- G(2)
 */ 


/*
 *	The tag of 'group' structure:
 *		1    User name
 *		2	 type
 *		3	 rank
 *		4    name
 *      5    list of represenations
 */

static List groups =0;

Term SpecToGroup(Term t)
	{
	Term fu, ar, tt;
	int  N;
	char *name;
	int type=-1;
	if(is_atom(t))
		{
		return CopyTerm(GetAtomProperty(t,OPR_GROUP));
		}
	if(!is_compound(t) || CompoundArity(t)!=1 || !is_integer(CompoundArg1(t)))
		{
		ErrorInfo(130);
		printf("  group '"); WriteTerm(t); printf("' not supported\n");
		return 0;
		}
	fu=CompoundName(t);
	ar=CompoundArg1(t);
	N=(int)IntegerValue(ar);
	name=AtomValue(fu);
	
	if(strcmp(name,"U")==0 || strcmp(name,"u")==0)
		{
		if(N!=1)
			{
			ErrorInfo(131);
			printf("  group U(%d) not supported\n",N);
			return 0;
			}
		type=0;
		goto fi;
		}
	if(strcmp(name,"SU")==0 || strcmp(name,"su")==0)
		{
		if(N<2)
			{
			ErrorInfo(132);
			printf("  group SU(%d) not supported\n",N);
			return 0;
			}
		N-=1;
		type=1;
		goto fi;
		}
	if(strcmp(name,"A")==0 || strcmp(name,"a")==0)
		{
		if(N<1)
			{
			ErrorInfo(133);
			printf("  group A(%d) not supported\n",N);
			return 0;
			}
		type=1;
		goto fi;
		}
	if((strcmp(name,"SO")==0 || strcmp(name,"so")==0) && (N%2)==0)
		{
		if(N<1)
			{
			ErrorInfo(134);
			printf("  group SO(%d) not supported\n",N);
			return 0;
			}
		N/=2;
		type=2;
		goto fi;
		}
	if((strcmp(name,"B")==0 || strcmp(name,"b")==0) && (N%2)==0)
		{
		if(N<1)
			{
			ErrorInfo(135);
			printf("  group B(%d) not supported\n",N);
			return 0;
			}
		N/=2;
		type=2;
		goto fi;
		}
	if((strcmp(name,"SO")==0 || strcmp(name,"so")==0) && (N%2)==1)
		{
		if(N<1)
			{
			ErrorInfo(136);
			printf("  group SO(%d) not supported\n",N);
			return 0;
			}
		type=4;
		N=(N-1)/2;
		goto fi;
		}
	if((strcmp(name,"B")==0 || strcmp(name,"b")==0) && (N%2)==1)
		{
		if(N<1)
			{
			ErrorInfo(137);
			printf("  group B(%d) not supported\n",N);
			return 0;
			}
		type=4;
		N=(N-1)/2;
		goto fi;
		}
	if(strcmp(name,"Sp")==0 || strcmp(name,"sp")==0)
		{
		if(N<1 || (N%2)!=0)
			{
			ErrorInfo(138);
			printf("  group SO(%d) not supported\n",N);
			return 0;
			}
		type=3;
		N/=2;
		goto fi;
		}
	if(strcmp(name,"C")==0 || strcmp(name,"c")==0) 
		{
		if(N<1 || (N%2)!=0)
			{
			ErrorInfo(139);
			printf("  group B(%d) not supported\n",N);
			return 0;
			}
		type=3;
		N/=2;
		goto fi;
		}
	if(strcmp(name,"e")==0 || strcmp(name,"E")==0) 
		{
		if(N<6 || N>8)
			{
			ErrorInfo(140);
			printf("  group E(%d) not supported\n",N);
			return 0;
			}
		type=5;
		goto fi;
		}
	if(strcmp(name,"f")==0 || strcmp(name,"F")==0) 
		{
		if(N!=4)
			{
			ErrorInfo(141);
			printf("  group F(%d) not supported\n",N);
			return 0;
			}
		type=6;
		goto fi;
		}	
	if(strcmp(name,"g")==0 || strcmp(name,"G")==0) 
		{
		if(N!=2)
			{
			ErrorInfo(142);
			printf("  group G(%d) not supported\n",N);
			return 0;
			}
		type=7;
		goto fi;
		}
fi:
	if(type==-1)
		{
		ErrorInfo(143);
		printf("  group '"); WriteTerm(t);
		printf("'not supported\n");
		return 0;
		}
	tt=MakeCompound(OPR_GROUP,5);
	SetCompoundArg(tt,1,t);
	SetCompoundArg(tt,2,NewInteger(type));
	SetCompoundArg(tt,3,NewInteger(N));
	return tt;
	}
	
Term InterfSpecToGroup(Term t, Term ind)
	{
	Term t1;
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	return SpecToGroup(t1);
	}
	
Term SpecToRepr(Term t)	
	{
	Term gr, arg, arg1, gspec, li;
	if(t==OPR_VECTOR)
		{
		return MakeCompound2(A_LORENTZ,NewInteger(2),NewInteger(2));
		}
	if(t==OPR_SPINOR)
		{
		return MakeCompound2(A_LORENTZ,NewInteger(1),NewInteger(0));
		}
		
	if(t==A_CSPINOR)
		{
		return MakeCompound2(A_LORENTZ,NewInteger(0),NewInteger(1));
		}
		
	if(is_atom(t) && strcmp(AtomValue(t),"spinor2")==0)
		return MakeCompound2(A_LORENTZ,NewInteger(0),NewInteger(0));
		
	if(is_compound(t) && CompoundName(t)==OPR_MINUS &&
		CompoundArity(t)==1 && CompoundArg1(t)==OPR_SPINOR)
		{
		FreeAtomic(t);
		return MakeCompound2(A_LORENTZ,NewInteger(0),NewInteger(1));
		}
	if(!is_compound(t) || CompoundArity(t)!=1)
		{
		ErrorInfo(144);
		printf("Error: bad representation spec '"); WriteTerm(t);
		printf("'\n");
		return 0;
		}
	gr=CompoundName(t);
	arg=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	if(gr==OPR_WILD)
		{
		if(!is_integer(arg) || IntegerValue(arg)<1)
			{
			ErrorInfo(145);
			printf("expected dimension of 'wild' repres, got '");
			WriteTerm(arg);
			printf("'\n");
			return 0;
			}
		return MakeCompound2(OPR_WILD,arg,arg);
		}
		
	gspec=GetAtomProperty(gr,OPR_GROUP);
	if(gspec==0)
		{
		ErrorInfo(146);
		printf("unknown group '%s' in representation.\n",AtomValue(gr));
		return 0;
		}
		
	if(IntegerValue(CompoundArg2(gspec))==0)
		{		/*   U(1) processing */
		if(is_integer(arg))
			return MakeCompound1(gr,MakeCompound2(OPR_DIV,arg,NewInteger(1)));
		if(is_compound(arg) && CompoundName(arg)==OPR_DIV
			&& CompoundArity(arg)==2 && is_integer(CompoundArg1(arg))
			&& is_integer(CompoundArg2(arg)) )
			return MakeCompound1(gr,arg);
		ErrorInfo(147);
		printf("repres of U(1) group must be rational, got '");
		WriteTerm(arg);
		printf("'\n");
		FreeAtomic(arg);
		return 0;
		}
		
	if(is_atom(arg))
		{
		li=CompoundArgN(gspec,5);
		while(!is_empty_list(li))
			{
			Term t1;
			t1=ListFirst(li);
			if(arg==CompoundArg1(t1))
				return MakeCompound2(gr,CompoundArg1(t1),CompoundArg2(t1));
			if(arg==CompoundArg2(t1))
				return MakeCompound2(gr,CompoundArg2(t1),CompoundArg1(t1));
			li=ListTail(li);
			}
		
		ErrorInfo(148);
		printf("unknown representation '%s' for group '%s'\n",
			AtomValue(arg),AtomValue(gr));
		return 0;
		}
	
	
		
	
	arg=CommaToList(arg);
	li=arg;
	while(!is_empty_list(li))
		{
		if(!is_integer(ListFirst(li)))
			{
			ErrorInfo(149);
			printf("bad representation spec ");
			WriteTerm(arg);
			printf("\n");
			return 0;
			}
		li=ListTail(li);
		}
	if(ListLength(arg)!=IntegerValue(CompoundArgN(gspec,3)))
		{
		ErrorInfo(150);
		printf("wrong length of repres spec ");WriteTerm(arg);puts("");
		return 0;
		}
	arg1=NewList();
	li=arg;
	while(!is_empty_list(li))
		{
		arg1=AppendFirst(arg1,ListFirst(li));
		li=ListTail(li);
		}
	return MakeCompound2(gr,arg,arg1);;
	}
	
Term InterfSpecToRepr(Term t, Term ind)
	{
	Term t1;
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	return SpecToRepr(t1);
	}

static void proc_group(Term t)
	{
	Atom grname;
	Term grp;

	if(is_compound(t) && FunctorName(CompoundFunctor(t))==OPR_COMMA)
		{
		Term t1,t2;
		t1=ConsumeCompoundArg(t,1);
		t2=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		proc_group(t1);
		proc_group(t2);
		return;
		}
		
	if(is_atom(t))
		{
		grname=t;
		grp=MakeCompound(OPR_GROUP,5);
		SetCompoundArg(grp,1,NewAtom("generic",0));
		SetCompoundArg(grp,2,NewInteger(-1));
		SetCompoundArg(grp,3,NewInteger(-1));
		SetCompoundArg(grp,4,grname);
		SetCompoundArg(grp,5,0);
		groups=AppendLast(groups,grname);
		SetAtomProperty(grname,OPR_GROUP,grp);
		SetOperator(grname,OP_FX,800);
		return;
		}
		
	if(!is_compound(t) || FunctorName(CompoundFunctor(t))!=OPR_COLON)
		{
		ErrorInfo(151);
		printf("Semantic error: illegal expression \'");
		WriteTerm(t);
		printf("\' in 'group' statement.\n");
		FreeAtomic(t);
		return;
		}
	
	grname=ConsumeCompoundArg(t,1);
	grp=ConsumeCompoundArg(t,2);
	FreeAtomic(t);
	if(!is_atom(grname))
		{
		printf("illegal group name \'");
		WriteTerm(grname);
		printf("\'.\n");
		FreeAtomic(grname);
		FreeAtomic(grp);
		return;
		}
		
	grp=SpecToGroup(grp);
	if(grp==0)
		return;
	SetCompoundArg(grp,4,grname);
	SetCompoundArg(grp,5,0);
	groups=AppendLast(groups,grname);
	SetAtomProperty(grname,OPR_GROUP,grp);
	SetOperator(grname,OP_FX,800);
		
	return;

	}

Term ProcessGroup(Term t, Term ind)
	{
	Term arg;
	char *s;
	arg=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	if(is_atom(arg))
		{
		s=AtomValue(arg);
		if(s[0]=='?' && s[1]==0)
			{
			List li;
			WriteTerm(groups);
			puts(""); 
			li=groups;
			while(!is_empty_list(li))
				{
				Term t1;
				t1=ListFirst(li);
				WriteTerm(GetAtomProperty(t1,OPR_GROUP));
				puts("");
				li=ListTail(li);
				}
			return 0;
			}
		}
	proc_group(arg);
	return 0;
	}
	
static void proc_repres(Term t)
	{
	Term gname, rlist, gspec;
	List l1;
	if(is_compound(t) && FunctorName(CompoundFunctor(t))==OPR_COMMA)
		{
		Term t1,t2;
		t1=ConsumeCompoundArg(t,1);
		t2=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		proc_repres(t1);
		proc_repres(t2);
		return;
		}
	
	if(!is_compound(t) || CompoundArity(t)!=2 || CompoundName(t)!=OPR_COLON)
		{
		ErrorInfo(152);
		printf("bad form of representation '"); WriteTerm(t); printf("'\n");
		FreeAtomic(t);
		return ;
		}
		
	gname=ConsumeCompoundArg(t,1);
	rlist=ConsumeCompoundArg(t,2);
	FreeAtomic(t);
	if(!is_atom(gname) || (gspec=GetAtomProperty(gname,OPR_GROUP))==0)
		{
		ErrorInfo(153);
		printf("bad group name '"); WriteTerm(gname); printf("' in 'repres' call\n");
		return ;
		} 
		
	if(CompoundArg2(gspec)==NewInteger(0))
		{
		ErrorInfo(154);
		printf("group U(1) must have a charge as repres spec\n");
		return ;
		}
		
	rlist=CommaToList(rlist);
	l1=rlist;
	while(!is_empty_list(l1))
		{
		Term t1;
		t1=ListFirst(l1);
		if(is_compound(t1) && CompoundArity(t1)==2 &&
			CompoundName(t1)==OPR_DIV && is_atom(CompoundArg1(t1)) &&
			is_atom(CompoundArg2(t1)))
				{
				l1=ListTail(l1);
				continue;
				}
		if(is_atom(t1))
			{
			ChangeList(l1,MakeCompound2(OPR_DIV,t1,t1));
			l1=ListTail(l1);
			continue;
			}
		ErrorInfo(155);
		printf("bad representation name '");
		WriteTerm(t1);
		printf("'\n");
		return;
		}
	SetCompoundArg(gspec,5,ConcatList(ConsumeCompoundArg(gspec,5),rlist));
	return;
	}


Term ProcessRepres(Term t, Term ind)
	{
	Term t1;
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	proc_repres(t1);
	return 0;
	}


int is_group(Atom g, Term *grp)
	{
	Term t;
	if(!is_atom(g))
		{
		if(grp!=NULL)
			*grp=0;
		return 0;
		}
	t=GetAtomProperty(g,OPR_GROUP);
	if(grp!=NULL)
		*grp=t;
	return t!=0;
	}

void ClearGroups(void)
	{
	if(!is_empty_list(groups))
		FreeAtomic(groups);
	groups=0;
	}
	

int equal_groups(Term g1, Term g2)
	{
	Term f1,f2;
	f1=CompoundArg1(g1);
	f2=CompoundArg1(g2);
	if(CompoundArity(g1)!=CompoundArity(g2))
		return 0;
	if(CompoundArity(g1)==3 && CompoundArgN(g1,3)!=CompoundArgN(g2,3))
		return 0;
	if(CompoundName(f1)==CompoundName(f2) &&
		EqualTerms(CompoundArg1(f1),CompoundArg2(f2)) &&
		EqualTerms(CompoundArg1(f2),CompoundArg2(f1)) )
				return 1;
	return 0;
	}


void ClearGroup(Atom g)
	{
	List l;
	l=groups;
	while(!is_empty_list(l))
		{
		if(ListFirst(l)==g)
			{
			groups=CutFromList(groups,l);
			break;
			}
		l=ListTail(l);
		}
	RemoveAtomProperty(g,OPR_GROUP);
	RemoveAtomProperty(g,OP_PREFIX);
	}
