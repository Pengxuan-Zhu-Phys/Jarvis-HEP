#include "lanhep.h"

Atom OP_FX, OP_FY, OP_XF, OP_YF, OP_XFX, OP_XFY, OP_YFX;	


static List  ops_list = 0;


void SetOperator(Atom name, Atom assoc, int prior)
	{
	Compound c;
	Functor fu;
	Atom arg=0;
	if(assoc==OP_FX || assoc==OP_FY)
		arg=OP_PREFIX;
	if(assoc==OP_XFX || assoc==OP_XFY || assoc==OP_YFX)
		arg=OP_INFIX;
	if(assoc==OP_XF || assoc==OP_YF)
		arg=OP_POSTFIX;
	fu=NewFunctor(name,3);
	c=NewCompound(fu);
	SetCompoundArg(c,1,arg);
	SetCompoundArg(c,2,assoc);
	SetCompoundArg(c,3,NewInteger(prior));
    SetAtomProperty(name,arg,CopyTerm(c));
	ops_list=AppendFirst(ops_list,c);
	}
	
Compound GetOperator(Atom name, Atom type, Atom *assoc, int *prior)
	{
/*	List l;
	Functor fu;
	fu=NewFunctor(name,3);
	l=ops_list;
	while(!is_empty_list(l))
		{
		Compound c;
		c=ListFirst(l);
		if(CompoundFunctor(c)==fu && 
			CompoundArg1(c)==type)
			{
            c=GetAtomProperty(name,type);
			if(assoc!=NULL)
				*assoc=CompoundArg2(c);
			if(prior!=NULL)
				*prior=IntegerValue(CompoundArgN(c,3));
			return c;
			}
		l=ListTail(l);
		}
    return 0;*/
    Compound c;
    if(!is_atom(name))
        return 0;
    c=GetAtomProperty(name,type);
    if(c==0)
	    return 0;
    if(assoc!=NULL)
			*assoc=CompoundArg2(c);
    if(prior!=NULL)
			*prior=(int)IntegerValue(CompoundArgN(c,3));
    /*printf("/// %s: %s %ld\n",AtomValue(name),AtomValue(CompoundArg2(c)),
                IntegerValue(CompoundArgN(c,3)));*/
	return c;
	}
	
Compound CurrentOperator(int no)
	{
	List l;
	l=ops_list;
	while(!is_empty_list(l))
		{
		no--;
		if(no==0) return ListFirst(l);
		l=ListTail(l);
		}
	return 0;
	}
	
	

