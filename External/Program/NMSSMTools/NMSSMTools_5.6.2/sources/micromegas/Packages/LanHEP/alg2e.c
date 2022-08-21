#include <string.h>
#include "lanhep.h"

static List herm_matrices = 0, hermi_matrices=0;
static Term el_name[16];
Label orth_del_offdiag = 0;
Label orth_red_diag = 0;

int herm_matr_dim(Label la)
{
	List l;
	for(l=herm_matrices;l;l=ListTail(l))
		if(CompoundArg1(ListFirst(l))==la)
			return CompoundArity(CompoundArg2(ListFirst(l)));
	return 0;
}

Atom herm_matr_el(Label la, int i, int j)
{
	List l;
	for(l=herm_matrices;l;l=ListTail(l))
		if(CompoundArg1(ListFirst(l))==la)
			return CompoundArgN(CompoundArgN(CompoundArg2(ListFirst(l)),i),j);
	return 0;
}


static List rcos_add(List l, Term m2);

static int prtcmp(Term p1, Term p2)
	{
		if(CompoundArg1(p1)==A_I)
			return -1;
		if(CompoundArg1(p2)==A_I)
			return 1;
	return strcmp(AtomValue(CompoundArg1(p1)),AtomValue(CompoundArg1(p2)));
	}

static int rhcf2, rhcf11, rhcfri;
	
static int alg2_red_herm_2(List *lm, List l1, Term prop)
	{
	int line1,line2,i;
    Term m1,par,matr=0, matri=0;
    List l,l2,ltc;

   line1=(int)IntegerValue(CompoundArgN(prop,4));
    line2=(int)IntegerValue(CompoundArgN(prop,3));
    	
    l=hermi_matrices;

    while(!is_empty_list(l))
    	{
    	if(CompoundArg1(ListFirst(l))==CompoundArg1(prop))
    		{
    		matr=CompoundArg2(ListFirst(l));
			matri=CompoundArgN(ListFirst(l),3);
    		break;
    		}
    	l=ListTail(l);
    	}
    if(is_empty_list(l))
    	{
    	puts("Internal error a2ro2");
    	return 0;
    	}
    m1=CopyTerm(ListFirst(l1));
    par=l=ConsumeCompoundArg(m1,2);
    l2=NewList();
    while(!is_empty_list(l))
    	{
    	Term q,pr;
    	q=ListFirst(l);
    	pr=GetAtomProperty(CompoundArg1(q),A_HERM_MATR);
    	if(pr!=prop)
    		{
    		l2=AppendLast(l2,q);
    		ChangeList(l,0);
    		}
    	l=ListTail(l);
    	}
    FreeAtomic(par);
    par=l2;
    ltc=NewList();

    for(i=1;i<=abs((int)IntegerValue(CompoundArg2(prop)));i++)
    	{
    	Term par1;
    	par1=CompoundArgN(CompoundArgN(matr,i),line1);

    	l2=CopyTerm(par);
    	l2=AppendFirst(l2,MakeCompound2(OPR_POW,par1,NewInteger(2)));
    	l2=SortedList(l2,prtcmp);
    	SetCompoundArg(m1,2,l2);
    	l=*lm;
    	while(!is_empty_list(l))
    		{
    		if(EqualTerms(m1,ListFirst(l)))
    			{
    			ltc=AppendFirst(ltc,l);
    			break;
    			}
    		l=ListTail(l);
    		}
    	if(is_empty_list(l))
    		{
    		/*WriteTerm(m1); WriteTerm(par);
    		puts(" --- fail 1 ");*/
    		goto try2;
    		}
    	l2=ConsumeCompoundArg(m1,2);
    	FreeAtomic(l2);
    	}
		
    for(i=1;i<=abs((int)IntegerValue(CompoundArg2(prop)));i++)
    	{
    	Term par1;
    	par1=CompoundArgN(CompoundArgN(matri,i),line1);

    	l2=CopyTerm(par);
    	l2=AppendFirst(l2,MakeCompound2(OPR_POW,par1,NewInteger(2)));
    	l2=SortedList(l2,prtcmp);
    	SetCompoundArg(m1,2,l2);
    	l=*lm;
    	while(!is_empty_list(l))
    		{
    		if(EqualTerms(m1,ListFirst(l)))
    			{
    			ltc=AppendFirst(ltc,l);
    			break;
    			}
    		l=ListTail(l);
    		}
    	if(is_empty_list(l))
    		{
    		/*WriteTerm(m1); WriteTerm(par);
    		puts(" --- fail 1 ");*/
    		goto try2;
    		}
    	l2=ConsumeCompoundArg(m1,2);
    	FreeAtomic(l2);
    	}
    	
goto end;
try2:
	RemoveList(ltc);
	ltc=NewList();
	for(i=1;i<=abs((int)IntegerValue(CompoundArg2(prop)));i++)
    	{
    	Term par1;
    	par1=CompoundArgN(CompoundArgN(matr,line1),i);


    	l2=CopyTerm(par);
    	l2=AppendFirst(l2,MakeCompound2(OPR_POW,par1,NewInteger(2)));
    	l2=SortedList(l2,prtcmp);
    	SetCompoundArg(m1,2,l2);
    	l=*lm;
    	while(!is_empty_list(l))
    		{
    		if(EqualTerms(m1,ListFirst(l)))
    			{
    			ltc=AppendFirst(ltc,l);
    			break;
    			}
    		l=ListTail(l);
    		}
    	if(is_empty_list(l))
    		{
    		/*WriteTerm(m1); WriteTerm(par);
    		puts(" --- fail 2 ");*/
    		FreeAtomic(m1);
    		return 0;
    		}
    	l2=ConsumeCompoundArg(m1,2);
    	FreeAtomic(l2);
    	}
		
	for(i=1;i<=abs((int)IntegerValue(CompoundArg2(prop)));i++)
    	{
    	Term par1;
    	par1=CompoundArgN(CompoundArgN(matri,line1),i);


    	l2=CopyTerm(par);
    	l2=AppendFirst(l2,MakeCompound2(OPR_POW,par1,NewInteger(2)));
    	l2=SortedList(l2,prtcmp);
    	SetCompoundArg(m1,2,l2);
    	l=*lm;
    	while(!is_empty_list(l))
    		{
    		if(EqualTerms(m1,ListFirst(l)))
    			{
    			ltc=AppendFirst(ltc,l);
    			break;
    			}
    		l=ListTail(l);
    		}
    	if(is_empty_list(l))
    		{
   			/*WriteTerm(m1); WriteTerm(par);
    		puts(" --- fail 2 ");*/
    		FreeAtomic(m1);
    		return 0;
    		}
    	l2=ConsumeCompoundArg(m1,2);
    	FreeAtomic(l2);
    	}

end:

	SetCompoundArg(m1,2,par);
    
    l=ltc;
    while(!is_empty_list(l))
    	{
    	l2=ListFirst(l);
    	*lm=CutFromList(*lm,l2);
    	l=ListTail(l);
    	}
    RemoveList(ltc);
    *lm=rcos_add(*lm,m1);
    /*printf("11: %d & %d lines, type %d \n",line1,line2,type);*/	
	rhcf2++;
    return 1;
	
	}
	
static int alg2_red_orth_2(List *lm, List l1, Term prop)
	{
	int line1,line2,i, inv;
    Term m1,par,matr=0;
    List l,l2,ltc;

    line1=(int)IntegerValue(CompoundArgN(prop,4));
    line2=(int)IntegerValue(CompoundArgN(prop,3));
	
	if(el_name[0]==el_name[1]) inv=0;
	else if(el_name[0]==GetAtomProperty(el_name[1],A_ANTI)) inv=1;
	else {puts("Internal error (ro2ii)"); return 0;}
	
    l=herm_matrices;

    while(!is_empty_list(l))
    	{
    	if(CompoundArg1(ListFirst(l))==CompoundArg1(prop))
    		{
    		matr=CompoundArg2(ListFirst(l));
    		break;
    		}
    	l=ListTail(l);
    	}
    if(is_empty_list(l))
    	{
    	puts("Internal error a2ro2");
    	return 0;
    	}
    m1=CopyTerm(ListFirst(l1));
    par=l=ConsumeCompoundArg(m1,2);
    l2=NewList();
    while(!is_empty_list(l))
    	{
    	Term q,pr;
    	q=ListFirst(l);
    	pr=GetAtomProperty(CompoundArg1(q),A_ORTH_MATR);
    	if(pr!=prop)
    		{
    		l2=AppendLast(l2,q);
    		ChangeList(l,0);
    		}
    	l=ListTail(l);
    	}
    FreeAtomic(par);
    par=l2;
    
    ltc=NewList();

    for(i=1;i<=IntegerValue(CompoundArg2(prop));i++)
    	{
    	Term par1;
    	par1=CompoundArgN(CompoundArgN(matr,i),line1);

    	l2=CopyTerm(par);
		if(!inv)
		  l2=AppendFirst(l2,MakeCompound2(OPR_POW,par1,NewInteger(2)));
		else
		{
		  l2=AppendFirst(l2,MakeCompound2(OPR_POW,par1,NewInteger(1)));
		  l2=AppendFirst(l2,MakeCompound2(OPR_POW,GetAtomProperty(par1,A_ANTI),NewInteger(1)));
		}
    	l2=SortedList(l2,prtcmp);
    	SetCompoundArg(m1,2,l2);
    	l=*lm;
    	while(!is_empty_list(l))
    		{
    		if(EqualTerms(m1,ListFirst(l)))
    			{
    			ltc=AppendFirst(ltc,l);
    			break;
    			}
    		l=ListTail(l);
    		}
    	if(is_empty_list(l))
    		{
    		/*WriteTerm(m1); WriteTerm(par);
    		puts(" --- fail 1 ");*/
    		goto try2;
    		}
    	l2=ConsumeCompoundArg(m1,2);
    	FreeAtomic(l2);
    	}
    	
goto end;
try2:
	RemoveList(ltc);
	ltc=NewList();
	for(i=1;i<=IntegerValue(CompoundArg2(prop));i++)
    	{
    	Term par1;
    	par1=CompoundArgN(CompoundArgN(matr,line1),i);


    	l2=CopyTerm(par);
		if(!inv)
		  l2=AppendFirst(l2,MakeCompound2(OPR_POW,par1,NewInteger(2)));
		else
		{
		  l2=AppendFirst(l2,MakeCompound2(OPR_POW,par1,NewInteger(1)));
		  l2=AppendFirst(l2,MakeCompound2(OPR_POW,GetAtomProperty(par1,A_ANTI),NewInteger(1)));
		}
	
    	l2=SortedList(l2,prtcmp);
    	SetCompoundArg(m1,2,l2);
    	l=*lm;
    	while(!is_empty_list(l))
    		{
    		if(EqualTerms(m1,ListFirst(l)))
    			{
    			ltc=AppendFirst(ltc,l);
    			break;
    			}
    		l=ListTail(l);
    		}
    	if(is_empty_list(l))
    		{
   /* 		WriteTerm(m1); WriteTerm(par);
    		puts(" --- fail 2 ");*/
    		FreeAtomic(m1);
    		return 0;
    		}
    	l2=ConsumeCompoundArg(m1,2);
    	FreeAtomic(l2);
    	}

end:

	SetCompoundArg(m1,2,par);
    
    l=ltc;
    while(!is_empty_list(l))
    	{
    	l2=ListFirst(l);
    	*lm=CutFromList(*lm,l2);
    	l=ListTail(l);
    	}
    RemoveList(ltc);
    *lm=rcos_add(*lm,m1);
    /*printf("11: %d & %d lines, type %d \n",line1,line2,type);*/	
    return 1;
	
	}
	
static int alg2_red_orth_11(List *lm, List l1, Term prop1, Term prop2)
    {
    int type,line1,line2,i;
	int inv1=0, inv2=0;
    Term m1,par,matr=0;
    List l,l2,ltc;
    if(CompoundArgN(prop1,3)==CompoundArgN(prop2,3))
    	{ 
    	type=1;
    	line1=(int)IntegerValue(CompoundArgN(prop1,4));
    	line2=(int)IntegerValue(CompoundArgN(prop2,4));
    	}
    else
    	if(CompoundArgN(prop1,4)==CompoundArgN(prop2,4))
    		{
    		type=2;
    		line1=(int)IntegerValue(CompoundArgN(prop1,3));
    		line2=(int)IntegerValue(CompoundArgN(prop2,3));
    		}
    	else
    		return 0;
    if(line1==line2)
    	return 0;
    
	
	
	
    l=herm_matrices;
    while(!is_empty_list(l))
    	{
    	if(CompoundArg1(ListFirst(l))==CompoundArg1(prop1))
    		{
    		matr=CompoundArg2(ListFirst(l));
    		break;
    		}
    	l=ListTail(l);
    	}
    if(is_empty_list(l))
    	{
    	puts("Internal error a2ro11");
    	return 0;
    	}
    
   
	{
	  Atom a;
	  a=CompoundArgN(CompoundArgN(matr,(int)IntegerValue(CompoundArgN(prop1,3))),
				(int)IntegerValue(CompoundArgN(prop1,4)));
	  if(a!=el_name[0])
	  {
		if(GetAtomProperty(a,A_ANTI)==el_name[0]) inv1=1;
		else puts("Internal error (ro11i1)");
	  }
	  a=CompoundArgN(CompoundArgN(matr,(int)IntegerValue(CompoundArgN(prop2,3))),
				(int)IntegerValue(CompoundArgN(prop2,4)));
	   if(a!=el_name[1])
	  {
		if(GetAtomProperty(a,A_ANTI)==el_name[1]) inv2=1;
		else puts("Internal error (ro11i2)");
	  }
	}
    
    
    
    m1=CopyTerm(ListFirst(l1));
    par=l=ConsumeCompoundArg(m1,2);
    l2=NewList();
    while(!is_empty_list(l))
    	{
    	Term q,pr;
    	q=ListFirst(l);
    	pr=GetAtomProperty(CompoundArg1(q),A_ORTH_MATR);
    	if(pr!=prop1 && pr!=prop2)
    		{
    		l2=AppendLast(l2,q);
    		ChangeList(l,0);
    		}
    	l=ListTail(l);
    	}
    FreeAtomic(par);
    par=l2;
    
    ltc=NewList();
    
    for(i=1;i<=IntegerValue(CompoundArg2(prop1));i++)
    	{
    	Term par1,par2;
    	if(type==1)
    		{
    		par1=CompoundArgN(CompoundArgN(matr,i),line1);
    		par2=CompoundArgN(CompoundArgN(matr,i),line2);
    		}
    	else
    		{
    		par1=CompoundArgN(CompoundArgN(matr,line1),i);
    		par2=CompoundArgN(CompoundArgN(matr,line2),i);
    		}
    	if(inv1) par1=GetAtomProperty(par1,A_ANTI);
		if(inv2) par2=GetAtomProperty(par2,A_ANTI);
    	l2=CopyTerm(par);
    	l2=AppendFirst(l2,MakeCompound2(OPR_POW,par1,NewInteger(1)));
    	l2=AppendFirst(l2,MakeCompound2(OPR_POW,par2,NewInteger(1)));
    	l2=SortedList(l2,prtcmp);
    	SetCompoundArg(m1,2,l2);
    	l=*lm;
    	while(!is_empty_list(l))
    		{
    		if(EqualTerms(m1,ListFirst(l)))
    			{
    			ltc=AppendFirst(ltc,l);
    			break;
    			}
    		l=ListTail(l);
    		}
    	if(is_empty_list(l))
    		{
/*    		WriteTerm(m1); WriteTerm(par);
    		puts(" --- fail"); 
*/
    		FreeAtomic(m1);
    		return 0;
    		}
    	l2=ConsumeCompoundArg(m1,2);
    	FreeAtomic(l2);
    	}

    FreeAtomic(par);
    FreeAtomic(m1);
    
    l=ltc;
    while(!is_empty_list(l))
    	{
    	l2=ListFirst(l);
    	*lm=CutFromList(*lm,l2);
    	l=ListTail(l);
    	}
    	
    RemoveList(ltc);	

    /*printf("11: %d & %d lines, type %d \n",line1,line2,type);	*/	
    return 1;
    }
	
static int alg2_red_herm_11(List *lm, List l1, Term prop1, Term prop2)
    {
    int type,line1,line2,i;
    Term m1,par,matr=0, matri=0;
    List l,l2,ltc;
	
    if(CompoundArgN(prop1,3)==CompoundArgN(prop2,3))
    	{ 
    	type=1;
    	line1=(int)IntegerValue(CompoundArgN(prop1,4));
    	line2=(int)IntegerValue(CompoundArgN(prop2,4));
    	}
    else
    	if(CompoundArgN(prop1,4)==CompoundArgN(prop2,4))
    		{
    		type=2;
    		line1=(int)IntegerValue(CompoundArgN(prop1,3));
    		line2=(int)IntegerValue(CompoundArgN(prop2,3));
    		}
    	else
    		return 0;
    if(line1==line2)
    	return 0;
    	
    l=hermi_matrices;
    while(!is_empty_list(l))
    	{
    	if(CompoundArg1(ListFirst(l))==CompoundArg1(prop1))
    		{
    		matr=CompoundArg2(ListFirst(l));
			matri=CompoundArgN(ListFirst(l),3);
    		break;
    		}
    	l=ListTail(l);
    	}
    if(is_empty_list(l))
    	{
    	puts("Internal error a2ro11");
    	return 0;
    	}
    	
    m1=CopyTerm(ListFirst(l1));
    par=l=ConsumeCompoundArg(m1,2);
    l2=NewList();
    while(!is_empty_list(l))
    	{
    	Term q,pr;
    	q=ListFirst(l);
    	pr=GetAtomProperty(CompoundArg1(q),A_HERM_MATR);
    	if(pr!=prop1 && pr!=prop2)
    		{
    		l2=AppendLast(l2,q);
    		ChangeList(l,0);
    		}
    	l=ListTail(l);
    	}
    FreeAtomic(par);
    par=l2;
    ltc=NewList();
    
    for(i=1;i<=abs((int)IntegerValue(CompoundArg2(prop1)));i++)
    	{
    	Term par1,par2;
    	if(type==1)
    		{
    		par1=CompoundArgN(CompoundArgN(matr,i),line1);
    		par2=CompoundArgN(CompoundArgN(matr,i),line2);
    		}
    	else
    		{
    		par1=CompoundArgN(CompoundArgN(matr,line1),i);
    		par2=CompoundArgN(CompoundArgN(matr,line2),i);
    		}
    	l2=CopyTerm(par);
    	l2=AppendFirst(l2,MakeCompound2(OPR_POW,par1,NewInteger(1)));
    	l2=AppendFirst(l2,MakeCompound2(OPR_POW,par2,NewInteger(1)));
    	l2=SortedList(l2,prtcmp);
    	SetCompoundArg(m1,2,l2);
    	l=*lm;
    	while(!is_empty_list(l))
    		{
    		if(EqualTerms(m1,ListFirst(l)))
    			{
    			ltc=AppendFirst(ltc,l);
    			break;
    			}
    		l=ListTail(l);
    		}
    	if(is_empty_list(l))
    		{
/*    		WriteTerm(m1); WriteTerm(par);
    		puts(" --- fail"); 
*/
    		FreeAtomic(m1);
    		return 0;
    		}
    	l2=ConsumeCompoundArg(m1,2);
    	FreeAtomic(l2);
    	}
    for(i=1;i<=abs((int)IntegerValue(CompoundArg2(prop1)));i++)
    	{
    	Term par1,par2;
    	if(type==1)
    		{
    		par1=CompoundArgN(CompoundArgN(matri,i),line1);
    		par2=CompoundArgN(CompoundArgN(matri,i),line2);
    		}
    	else
    		{
    		par1=CompoundArgN(CompoundArgN(matri,line1),i);
    		par2=CompoundArgN(CompoundArgN(matri,line2),i);
    		}
    	l2=CopyTerm(par);
    	l2=AppendFirst(l2,MakeCompound2(OPR_POW,par1,NewInteger(1)));
    	l2=AppendFirst(l2,MakeCompound2(OPR_POW,par2,NewInteger(1)));
    	l2=SortedList(l2,prtcmp);
    	SetCompoundArg(m1,2,l2);
    	l=*lm;
    	while(!is_empty_list(l))
    		{
    		if(EqualTerms(m1,ListFirst(l)))
    			{
    			ltc=AppendFirst(ltc,l);
    			break;
    			}
    		l=ListTail(l);
    		}
    	if(is_empty_list(l))
    		{
/*    		WriteTerm(m1); WriteTerm(par);
    		puts(" --- fail"); 
*/
    		FreeAtomic(m1);
    		return 0;
    		}
    	l2=ConsumeCompoundArg(m1,2);
    	FreeAtomic(l2);
    	}

    FreeAtomic(par);
    FreeAtomic(m1);
    
    l=ltc;
    while(!is_empty_list(l))
    	{
    	l2=ListFirst(l);
    	*lm=CutFromList(*lm,l2);
    	l=ListTail(l);
    	}
    	
    RemoveList(ltc);	

    /*printf("11: %d & %d lines, type %d \n",line1,line2,type);	*/	
	rhcf11++;
    return 1;
    }
	
static int alg2_red_herm_ri(List *lm, List l1, Term prop1, Term prop2)
    {
    int type,line1,line2,i;
    Term m1,par,matr=0, matri=0;
    List l,l2,ltc;
	
	
	if(IntegerValue(CompoundArg2(prop1))<0)
	{
		Term tmp=prop1; prop1=prop2; prop2=tmp;
	}
	
    if(CompoundArgN(prop1,3)==CompoundArgN(prop2,3))
    	{ 
    	type=1;
    	line1=(int)IntegerValue(CompoundArgN(prop1,4));
    	line2=(int)IntegerValue(CompoundArgN(prop2,4));
    	}
    else
    	if(CompoundArgN(prop1,4)==CompoundArgN(prop2,4))
    		{
    		type=2;
    		line1=(int)IntegerValue(CompoundArgN(prop1,3));
    		line2=(int)IntegerValue(CompoundArgN(prop2,3));
    		}
    	else
    		return 0;
bgn:    	
    l=hermi_matrices;
    while(!is_empty_list(l))
    	{
    	if(CompoundArg1(ListFirst(l))==CompoundArg1(prop1))
    		{
    		matr=CompoundArg2(ListFirst(l));
			matri=CompoundArgN(ListFirst(l),3);
    		break;
    		}
    	l=ListTail(l);
    	}
    if(is_empty_list(l))
    	{
    	puts("Internal error a2ro11");
    	return 0;
    	}
    	
    m1=CopyTerm(ListFirst(l1));
    par=l=ConsumeCompoundArg(m1,2);
    l2=NewList();
	
    while(!is_empty_list(l))
    	{
    	Term q,pr;
    	q=ListFirst(l);
    	pr=GetAtomProperty(CompoundArg1(q),A_HERM_MATR);
    	if(pr!=prop1 && pr!=prop2)
    		{
    		l2=AppendLast(l2,q);
    		ChangeList(l,0);
    		}
    	l=ListTail(l);
    	}
    FreeAtomic(par);
    par=l2;
    ltc=NewList();
	
	/*SetCompoundArg(m1,1,NewInteger(-IntegerValue(CompoundArg1(m1))));*/
    
    for(i=1;i<=abs((int)IntegerValue(CompoundArg2(prop1)));i++)
    	{
    	Term par1,par2;
    	if(type==1)
    		{
    		par1=CompoundArgN(CompoundArgN(matr,i),line1);
    		par2=CompoundArgN(CompoundArgN(matri,i),line2);
    		}
    	else
    		{
    		par1=CompoundArgN(CompoundArgN(matr,line1),i);
    		par2=CompoundArgN(CompoundArgN(matri,line2),i);
    		}
    	l2=CopyTerm(par);
    	l2=AppendFirst(l2,MakeCompound2(OPR_POW,par1,NewInteger(1)));
    	l2=AppendFirst(l2,MakeCompound2(OPR_POW,par2,NewInteger(1)));
    	l2=SortedList(l2,prtcmp);
    	SetCompoundArg(m1,2,l2);
    	l=*lm;
    	while(!is_empty_list(l))
    		{
    		if(EqualTerms(m1,ListFirst(l)))
    			{
    			ltc=AppendFirst(ltc,l);
    			break;
    			}
    		l=ListTail(l);
    		}
    	if(is_empty_list(l))
    		{
    		/*WriteTerm(m1); WriteTerm(par);
    		puts(" --- fail"); */

    		FreeAtomic(m1);
			
			if(type==1 && line1==line2)
			{
				type=2;
				line1=(int)IntegerValue(CompoundArgN(prop1,3));
    			line2=(int)IntegerValue(CompoundArgN(prop2,3));
				goto bgn;
			}

    		return 0;
    		}
    	l2=ConsumeCompoundArg(m1,2);
    	FreeAtomic(l2);
    	}
		
	SetCompoundArg(m1,1,NewInteger(-IntegerValue(CompoundArg1(m1))));
	
    for(i=1;i<=abs((int)IntegerValue(CompoundArg2(prop1)));i++)
    	{
    	Term par1,par2;
    	if(type==1)
    		{
    		par1=CompoundArgN(CompoundArgN(matri,i),line1);
    		par2=CompoundArgN(CompoundArgN(matr,i),line2);
    		}
    	else
    		{
    		par1=CompoundArgN(CompoundArgN(matri,line1),i);
    		par2=CompoundArgN(CompoundArgN(matr,line2),i);
    		}
    	l2=CopyTerm(par);
    	l2=AppendFirst(l2,MakeCompound2(OPR_POW,par1,NewInteger(1)));
    	l2=AppendFirst(l2,MakeCompound2(OPR_POW,par2,NewInteger(1)));
    	l2=SortedList(l2,prtcmp);
    	SetCompoundArg(m1,2,l2);
    	l=*lm;
    	while(!is_empty_list(l))
    		{
    		if(EqualTerms(m1,ListFirst(l)))
    			{
    			ltc=AppendFirst(ltc,l);
    			break;
    			}
    		l=ListTail(l);
    		}
    	if(is_empty_list(l))
    		{
    		/*WriteTerm(m1); WriteTerm(par);
    		puts(" --- fail2"); */
					
			if(type==1 && line1==line2)
			{
				type=2;
				line1=(int)IntegerValue(CompoundArgN(prop1,3));
    			line2=(int)IntegerValue(CompoundArgN(prop2,3));
				goto bgn;
			}

    		FreeAtomic(m1);
    		return 0;
    		}
    	l2=ConsumeCompoundArg(m1,2);
    	FreeAtomic(l2);
    	}

    FreeAtomic(par);
    FreeAtomic(m1);
    
    l=ltc;
    while(!is_empty_list(l))
    	{
    	l2=ListFirst(l);
    	*lm=CutFromList(*lm,l2);
    	l=ListTail(l);
    	}
    	
    RemoveList(ltc);	

/*    printf("11: %d & %d lines, type %d \n",line1,line2,type);	*/
	rhcfri++;
    return 1;
    }
	
static List orth_mlt(List l, Atom par)
	{
	List l1;
	l1=l;
	while(!is_empty_list(l1))
		{
		Term pp;
		pp=ListFirst(l1);
		if(CompoundArg1(pp)==par)
			{
			int n;
			n=1+(int)IntegerValue(CompoundArg2(pp));
			SetCompoundArg(pp,2,NewInteger(n));
			return l;
			}
		l1=ListTail(l1);
		}
	l=AppendLast(l,MakeCompound2(OPR_POW,par,NewInteger(1)));
	return l;
	} 		


/* tp==1 means cols fixed rows changing and equal, tp==2 -- rows fixed */
	
static List orth4_mkprod(Term matr, int dim, int tp1, int tp2, int k1, int k2,
		int k3, int k4)
{
	List ret=0,l;
	int i,j;
	for(i=1;i<=dim;i++)
	for(j=1;j<=dim;j++)
	{
		List m=0;
		Term par1,par2, par3, par4;
    	
		if(tp1==1)
    		{
    		par1=CompoundArgN(CompoundArgN(matr,i),k1);
    		par2=CompoundArgN(CompoundArgN(matr,i),k2);
    		}
    	else
    		{
    		par1=CompoundArgN(CompoundArgN(matr,k1),i);
    		par2=CompoundArgN(CompoundArgN(matr,k2),i);
    		}
		if(tp2==1)
    		{
    		par3=CompoundArgN(CompoundArgN(matr,j),k3);
    		par4=CompoundArgN(CompoundArgN(matr,j),k4);
    		}
    	else
    		{
    		par3=CompoundArgN(CompoundArgN(matr,k3),j);
    		par4=CompoundArgN(CompoundArgN(matr,k4),j);
    		}
		m=orth_mlt(m,par1);
		m=orth_mlt(m,par2);
		m=orth_mlt(m,par3);
		m=orth_mlt(m,par4);
		m=SortedList(m,prtcmp);
		
		for(l=ret;l;l=ListTail(l))
			if(EqualTerms(CompoundArg2(ListFirst(l)),m))
			{
				FreeAtomic(m);
				SetCompoundArg(ListFirst(l),1,
						NewInteger(1+IntegerValue(CompoundArg1(ListFirst(l)))));
				break;
			}
		
		if(is_empty_list(l))
			ret=AppendLast(ret,MakeCompound2(OPR_MLT,NewInteger(1),m));
	}
	
	return ret;
}
		
static List orth3_mkprod(Term matr, int dim, int tp1, int k1, int k2,
		Atom a3)
{
	List ret=0;
	int i;
	for(i=1;i<=dim;i++)
	{
		List m=0;
		Term par1,par2;
    	
		if(tp1==1)
    		{
    		par1=CompoundArgN(CompoundArgN(matr,i),k1);
    		par2=CompoundArgN(CompoundArgN(matr,i),k2);
    		}
    	else
    		{
    		par1=CompoundArgN(CompoundArgN(matr,k1),i);
    		par2=CompoundArgN(CompoundArgN(matr,k2),i);
    		}
			
		m=orth_mlt(m,par1);
		m=orth_mlt(m,par2);
		m=orth_mlt(m,a3);
		m=SortedList(m,prtcmp);

		ret=AppendLast(ret,MakeCompound2(OPR_MLT,NewInteger(1),m));
	}
	
	return ret;
}

static int orth4_fit(List *ml, Term m1, List par, int val, List els)
{
	List ltc,par1,l1,l2;
	int coeff=0;
	int zpo=0;
	List zcf=0;
		
	ltc=NewList();
	
	zcf=CopyTerm(CompoundArg2(ListFirst(els)));
	for(l1=zcf;l1;l1=ListTail(l1))
		zpo+=(int)IntegerValue(CompoundArg2(ListFirst(l1)));
		
	for(l1=ListTail(els);l1;l1=ListTail(l1))
	{
		for(l2=zcf;l2;l2=ListTail(l2))
		{
			List l3;
			int n1,n2;
			Atom z=CompoundArg1(ListFirst(l2));
			for(l3=CompoundArg2(ListFirst(l1));l3;l3=ListTail(l3))
				if(CompoundArg1(ListFirst(l3))==z)
					break;
			n1=(int)IntegerValue(CompoundArg2(ListFirst(l2)));
			n2=l3?(int)IntegerValue(CompoundArg2(ListFirst(l3))):0;
			if(n2<n1)
				SetCompoundArg(ListFirst(l2),2,NewInteger(n2));
		}
	}
	
rpt:	
	for(l1=zcf;l1;l1=ListTail(l1))
		if(CompoundArg2(ListFirst(l1))==NewInteger(0))
			{
			zcf=CutFromList(zcf,l1);
			goto rpt;
			}
	
	if((zpo==4 && zcf) || (zpo==3 && !zcf))
	{
		puts("Internal error (orth34zcf)");
	}
	
	for(l1=els;l1;l1=ListTail(l1))
	{
		par1=ConcatList(CopyTerm(par),ConsumeCompoundArg(ListFirst(l1),2));
		par1=SortedList(par1,prtcmp);
		
		for(l2=*ml;l2;l2=ListTail(l2))
			if(EqualTerms(CompoundArg2(ListFirst(l2)),par1) &&
				EqualTerms(CompoundArgN(ListFirst(l2),3),CompoundArgN(m1,3)))
				break;
		FreeAtomic(par1);
		
		if(is_empty_list(l2))
		{
			FreeAtomic(els);
			return 0;
		}
		
		if(coeff==0)
			coeff=(int)IntegerValue(CompoundArg1(ListFirst(l2)))/
					(int)IntegerValue(CompoundArg1(ListFirst(l1)));
		if(coeff*IntegerValue(CompoundArg1(ListFirst(l1)))!=
				IntegerValue(CompoundArg1(ListFirst(l2))))
		{
			FreeAtomic(els);
			return 0;
		}
		
		ltc=AppendLast(ltc,l2);
		
	}
	
	
	FreeAtomic(els);
	for(l1=ltc;l1;l1=ListTail(l1))
		*ml=CutFromList(*ml,ListFirst(l1));
	RemoveList(ltc);
	
	if(val)
	{
		if(zcf)
			par=ConcatList(par,zcf);
		SetCompoundArg(m1,2,par);
		SetCompoundArg(m1,1,NewInteger(coeff));
		*ml=rcos_add(*ml,m1);
	}
	else
	{
		FreeAtomic(m1);
		FreeAtomic(par);
	}
	
	return 1;
}

	
static int alg2_red_orth_4(List *lm, List l1, Term *props, Term *names)
    {
    int dim;
    Term m1,par,matr=0;
    List l,l2;
	int rows[4], cols[4];
    
	rows[0]=(int)IntegerValue(CompoundArgN(props[0],3));
	rows[1]=(int)IntegerValue(CompoundArgN(props[1],3));
	rows[2]=(int)IntegerValue(CompoundArgN(props[2],3));
	rows[3]=(int)IntegerValue(CompoundArgN(props[3],3));
	
	cols[0]=(int)IntegerValue(CompoundArgN(props[0],4));
	cols[1]=(int)IntegerValue(CompoundArgN(props[1],4));
	cols[2]=(int)IntegerValue(CompoundArgN(props[2],4));
	cols[3]=(int)IntegerValue(CompoundArgN(props[3],4));
	
	dim=(int)IntegerValue(CompoundArg2(props[0]));
	
	l=herm_matrices;
    while(!is_empty_list(l))
    	{
    	if(CompoundArg1(ListFirst(l))==CompoundArg1(props[0]))
    		{
    		matr=CompoundArg2(ListFirst(l));
    		break;
    		}
    	l=ListTail(l);
    	}
		
    if(is_empty_list(l))
    	{
    	puts("Internal error a2ro114");
    	return 0;
    	}
    	
    m1=CopyTerm(ListFirst(l1));
    par=l=ConsumeCompoundArg(m1,2);
    l2=NewList();
    while(!is_empty_list(l))
    	{
    	Term q,pr;
    	q=ListFirst(l);
    	pr=CompoundArg1(q);
    	if(pr!=names[0] && pr!=names[1] && pr!=names[2] && pr!=names[3])
    		{
    		l2=AppendLast(l2,q);
    		ChangeList(l,0);
    		}
    	l=ListTail(l);
    	}
    FreeAtomic(par);
    par=l2;
    
	if(names[0]==names[3])
	{
		if(rows[0]==cols[0])
		{
			if(orth4_fit(lm,m1,par,1,orth4_mkprod(matr,dim,1,1,rows[0],
					rows[0],rows[0],rows[0])))
				return 1;
			
			if(orth4_fit(lm,m1,par,1,orth4_mkprod(matr,dim,2,2,rows[0],
					rows[0],rows[0],rows[0])))
				return 1;
			FreeAtomic(m1);
			FreeAtomic(par);
			return 0;
		}
		
		FreeAtomic(m1);
		FreeAtomic(par);
		return 0;
	}
	
	if(names[0]==names[1] && names[1]!=names[2] && names[2]==names[3])
	{
		if(orth4_fit(lm,m1,par,1,orth4_mkprod(matr,dim,1,1,cols[0],
					cols[0],cols[2],cols[2])))
			return 1;
		if(orth4_fit(lm,m1,par,1,orth4_mkprod(matr,dim,2,2,rows[0],
					rows[0],rows[2],rows[2])))
			return 1;
		
		if(orth4_fit(lm,m1,par,0,orth4_mkprod(matr,dim,1,1,cols[0],
					cols[2],cols[0],cols[2])))
			return 1;
		if(orth4_fit(lm,m1,par,0,orth4_mkprod(matr,dim,2,2,rows[0],
					rows[2],rows[0],rows[2])))
			return 1;
	}
	
	if(orth4_fit(lm,m1,par,cols[0]==cols[1] && cols[2]==cols[3],
			orth4_mkprod(matr,dim,1,1,cols[0],cols[1],cols[2],cols[3])))
		return 1;
	if(orth4_fit(lm,m1,par,cols[0]==cols[2] && cols[1]==cols[3],
			orth4_mkprod(matr,dim,1,1,cols[0],cols[2],cols[1],cols[3])))
		return 1;
	if(orth4_fit(lm,m1,par,cols[0]==cols[3] && cols[1]==cols[2],
			orth4_mkprod(matr,dim,1,1,cols[0],cols[3],cols[1],cols[2])))
		return 1;
	
	if(orth4_fit(lm,m1,par,rows[0]==rows[1] && rows[2]==rows[3],
			orth4_mkprod(matr,dim,2,2,rows[0],rows[1],rows[2],rows[3])))
		return 1;
	if(orth4_fit(lm,m1,par,rows[0]==rows[2] && rows[1]==rows[3],
			orth4_mkprod(matr,dim,2,2,rows[0],rows[2],rows[1],rows[3])))
		return 1;
	if(orth4_fit(lm,m1,par,rows[0]==rows[3] && rows[1]==rows[2],
			orth4_mkprod(matr,dim,2,2,rows[0],rows[3],rows[1],rows[2])))
		return 1;
	
	FreeAtomic(m1);
	FreeAtomic(par);
	return 0;	

    }        
	
static int alg2_red_orth_3(List *lm, List l1, Term *props, Term *names)
    {
    int dim;
    Term m1,par,matr=0;
    List l,l2;
	int rows[3], cols[3];
    
	
	rows[0]=(int)IntegerValue(CompoundArgN(props[0],3));
	rows[1]=(int)IntegerValue(CompoundArgN(props[1],3));
	rows[2]=(int)IntegerValue(CompoundArgN(props[2],3));
	
	cols[0]=(int)IntegerValue(CompoundArgN(props[0],4));
	cols[1]=(int)IntegerValue(CompoundArgN(props[1],4));
	cols[2]=(int)IntegerValue(CompoundArgN(props[2],4));
	
	dim=(int)IntegerValue(CompoundArg2(props[0]));
	
	l=herm_matrices;
    while(!is_empty_list(l))
    	{
    	if(CompoundArg1(ListFirst(l))==CompoundArg1(props[0]))
    		{
    		matr=CompoundArg2(ListFirst(l));
    		break;
    		}
    	l=ListTail(l);
    	}
		
    if(is_empty_list(l))
    	{
    	puts("Internal error a2ro113");
    	return 0;
    	}
    	
    m1=CopyTerm(ListFirst(l1));
    par=l=ConsumeCompoundArg(m1,2);
    l2=NewList();
    while(!is_empty_list(l))
    	{
    	Term q,pr;
    	q=ListFirst(l);
    	pr=CompoundArg1(q);
    	if(pr!=names[0] && pr!=names[1] && pr!=names[2])
    		{
    		l2=AppendLast(l2,q);
    		ChangeList(l,0);
    		}
    	l=ListTail(l);
    	}
    FreeAtomic(par);
    par=l2;
    	
	
	if(orth4_fit(lm,m1,par,cols[0]==cols[1],
			orth3_mkprod(matr,dim,1,cols[0],cols[1],names[2])))
		return 1;
	if(orth4_fit(lm,m1,par,cols[0]==cols[2],
			orth3_mkprod(matr,dim,1,cols[0],cols[2],names[1])))
		return 1;
	if(orth4_fit(lm,m1,par,cols[1]==cols[2],
			orth3_mkprod(matr,dim,1,cols[1],cols[2],names[0])))
		return 1;
	
	if(orth4_fit(lm,m1,par,rows[0]==rows[1],
			orth3_mkprod(matr,dim,2,rows[0],rows[1],names[2])))
		return 1;
	if(orth4_fit(lm,m1,par,rows[0]==rows[2],
			orth3_mkprod(matr,dim,2,rows[0],rows[2],names[1])))
		return 1;
	if(orth4_fit(lm,m1,par,rows[1]==rows[2],
			orth3_mkprod(matr,dim,2,rows[1],rows[2],names[0])))
		return 1;
	
	FreeAtomic(m1);
	FreeAtomic(par);
	return 0;	

    }	
	
static void alg2_red_orth_mass(Term a2)
{
	List lm,l0,l,llabs=0;
/*	WriteTerm(a2); puts(" ->"); */
	
	
	lm=ConsumeCompoundArg(a2,5);
	
	if(orth_del_offdiag)
	{
		
	b1:
		for(l=lm;l;l=ListTail(l))
		{
			List l1,l2;
			Term ss=0;
			
			l1=ConsumeCompoundArg(ListFirst(l),2);
		b2:
			for(l2=l1;l2;l2=ListTail(l2))
			{
				ss=GetAtomProperty(CompoundArg1(ListFirst(l2)),A_ORTH_MATR);
				if(ss && CompoundArg1(ss)==orth_del_offdiag)
					break;
			}
			
			if(is_empty_list(l2))
			{
				SetCompoundArg(ListFirst(l),2,l1);
				continue;
			}
			
			l1=CutFromList(l1,l2);
			if(CompoundArgN(ss,3)==CompoundArgN(ss,4))
				goto b2;
			else
			{
				FreeAtomic(l1);
				lm=CutFromList(lm,l);
				goto b1;
			}
		}
		
	}
	
	
	for(l=lm;l;l=ListTail(l))
	{
		List l1;
		for(l1=CompoundArg2(ListFirst(l));l1;l1=ListTail(l1))
		{
			Term ss;
			ss=GetAtomProperty(CompoundArg1(ListFirst(l1)),A_ORTH_MATR);
			if(ss && !ListMember(llabs,CompoundArg1(ss)))
				llabs=AppendLast(llabs,CompoundArg1(ss));
		}
	}
	
	
	for(l0=llabs;l0;l0=ListTail(l0))
	{
		Term lmsv;
		int compl=1;
		
		lmsv=CopyTerm(lm);
	
bgn:


		l=lm;
/*
		while(!is_empty_list(l))
			{
			Term m1;
			List l1,l2;
			m1=ListFirst(l);
			l1=CompoundArg2(m1);
			while(!is_empty_list(l1))
				{
				Term ss;
				if((ss=GetAtomProperty(CompoundArg1(ListFirst(l1)),A_ORTH_MATR)) && 
					CompoundArg1(ss)==ListFirst(l0) &&
					IntegerValue(CompoundArg2(ListFirst(l1)))==2)
					{
					if(alg2_red_orth_2(&lm,l,ss))
						goto bgn;
					}
				if((ss=GetAtomProperty(CompoundArg1(ListFirst(l1)),A_ORTH_MATR)) && 
					CompoundArg1(ss)==ListFirst(l0) &&
					IntegerValue(CompoundArg2(ListFirst(l1)))==1)
					{
					l2=ListTail(l1);
					while(!is_empty_list(l2))
						{
						Term ss1;
						if((ss1=GetAtomProperty(
							CompoundArg1(ListFirst(l2)),A_ORTH_MATR)) &&
							IntegerValue(CompoundArg2(ListFirst(l2)))==1 &&
							CompoundArg1(ss)==CompoundArg1(ss1))
							{
							if(alg2_red_orth_11(&lm,l,ss,ss1))
								goto bgn;
							}
						l2=ListTail(l2);
						}
					}
				l1=ListTail(l1);
				}
			l=ListTail(l);
			}
*/
		while(!is_empty_list(l))
			{
			Term m1;
			List l1;
			Term ss[16],sn[16];
			int no=0, nn;
			
			m1=ListFirst(l);
			l1=CompoundArg2(m1);
			
			while(!is_empty_list(l1))
			{
				Term ss1;
				ss1=GetAtomProperty(CompoundArg1(ListFirst(l1)),A_ORTH_MATR);
				if(ss1 && CompoundArg1(ss1)==ListFirst(l0))
					for(nn=1;nn<=IntegerValue(CompoundArg2(ListFirst(l1)));nn++)
					{
						ss[no]=ss1;
						el_name[no]=sn[no]=CompoundArg1(ListFirst(l1));
						no++;
					}
				l1=ListTail(l1);
			}
			
			
			if(no==2 && (sn[0]==sn[1] || sn[0]==GetAtomProperty(sn[1],A_ANTI))
				  && alg2_red_orth_2(&lm,l,ss[0]))
				goto bgn;
			if(no==2 && sn[0]!=sn[1] && alg2_red_orth_11(&lm,l,ss[0],ss[1]))
				goto bgn;
			
			if(no==3 && alg2_red_orth_3(&lm,l,(Term *)ss,(Term *)sn))
				goto bgn;
			
			/*if(no==4 && alg2_red_orth_4(&lm,l,(Term *)ss,(Term *)sn))
				goto bgn;*/	

			l=ListTail(l);
			}

		if(orth_red_diag!=ListFirst(l0))
		{			
			for(l=lm;l;l=ListTail(l))
			{
				List l1;
				Term ss;
				for(l1=CompoundArg2(ListFirst(l));l1;l1=ListTail(l1))
					if((ss=GetAtomProperty(CompoundArg1(ListFirst(l1)),A_ORTH_MATR))
							 && CompoundArg1(ss)==ListFirst(l0))
						compl=0;
			}

			if(compl)
				FreeAtomic(lmsv);
			else
			{
				FreeAtomic(lm);
				lm=lmsv;
			}
	   }
	   else
		   FreeAtomic(lmsv);
	}
		
	SetCompoundArg(a2,5,lm);	
	RemoveList(llabs);
	/*
	WriteTerm(a2); puts(""); getchar(); 			
	*/
}

void alg2_red_herm(Term a2)
	{
	List lm,l,l0,llabs=0;

	rhcf2=rhcf11=rhcfri=0;
	
	lm=ConsumeCompoundArg(a2,5);
	
	for(l=lm;l;l=ListTail(l))
	{
		List l1;
		for(l1=CompoundArg2(ListFirst(l));l1;l1=ListTail(l1))
		{
			Term ss;
			ss=GetAtomProperty(CompoundArg1(ListFirst(l1)),A_HERM_MATR);
			if(ss && !ListMember(llabs,CompoundArg1(ss)))
				llabs=AppendLast(llabs,CompoundArg1(ss));
		}
	}
	
	
	for(l0=llabs;l0;l0=ListTail(l0))
	{
		
	bgn:
	
		l=lm;
		while(!is_empty_list(l))
			{
			Term m1;
			List l1;
			Term ss[16],sn[16];
			int no=0, nn;
			
			m1=ListFirst(l);
			l1=CompoundArg2(m1);
			
			while(!is_empty_list(l1))
			{
				Term ss1;
				ss1=GetAtomProperty(CompoundArg1(ListFirst(l1)),A_HERM_MATR);
				if(ss1 && CompoundArg1(ss1)==ListFirst(l0))
					for(nn=1;nn<=IntegerValue(CompoundArg2(ListFirst(l1)));nn++)
					{
						ss[no]=ss1;
						sn[no]=CompoundArg1(ListFirst(l1));
						no++;
					}
				l1=ListTail(l1);
			}
			
			if(no==2 && sn[0]==sn[1] && alg2_red_herm_2(&lm,l,ss[0]))
				goto bgn;
			if(no==2 && sn[0]!=sn[1] && CompoundArg2(ss[0])==CompoundArg2(ss[1])
					&& alg2_red_herm_11(&lm,l,ss[0],ss[1]))
				goto bgn;
			if(no==2 && sn[0]!=sn[1] && CompoundArg2(ss[0])!=CompoundArg2(ss[1])
					&& alg2_red_herm_ri(&lm,l,ss[0],ss[1]))
				goto bgn;
			
			l=ListTail(l);
			}
	}
	SetCompoundArg(a2,5,lm);	
	
/*	if(rhcf2||rhcf11||rhcfri)
	{
		WriteVertex(CompoundArg1(a2));
		printf(" r2=%d r11=%d rri=%d\n",rhcf2,rhcf11,rhcfri);
	}*/

/*	WriteTerm(a2); puts(""); */

	}

void alg2_red_orth(Term a2)
	{
	List lm,l,l0,llabs=0;

/*	WriteTerm(a2); puts(" ->"); */
	
	if(hermi_matrices)
		alg2_red_herm(a2);
	
	if(ListLength(CompoundArg1(a2))==2)
	{
		alg2_red_orth_mass(a2);
		return;
	}
	
	lm=ConsumeCompoundArg(a2,5);
	
	for(l=lm;l;l=ListTail(l))
	{
		List l1;
		for(l1=CompoundArg2(ListFirst(l));l1;l1=ListTail(l1))
		{
			Term ss;
			ss=GetAtomProperty(CompoundArg1(ListFirst(l1)),A_ORTH_MATR);
			if(ss && !ListMember(llabs,CompoundArg1(ss)))
				llabs=AppendLast(llabs,CompoundArg1(ss));
		}
	}
	
	
	for(l0=llabs;l0;l0=ListTail(l0))
	{
		
	bgn:

		l=lm;
		while(!is_empty_list(l))
			{
			Term m1;
			List l1;
			Term ss[16],sn[16];
			int no=0, nn;
			
			m1=ListFirst(l);
			l1=CompoundArg2(m1);
			
			while(!is_empty_list(l1))
			{
				Term ss1;
				ss1=GetAtomProperty(CompoundArg1(ListFirst(l1)),A_ORTH_MATR);
				if(ss1 && CompoundArg1(ss1)==ListFirst(l0))
					for(nn=1;nn<=IntegerValue(CompoundArg2(ListFirst(l1)));nn++)
					{
						ss[no]=ss1;
						el_name[no]=sn[no]=CompoundArg1(ListFirst(l1));
						no++;
					}
				l1=ListTail(l1);
			}
			
			
			if(no==2 && (sn[0]==sn[1] || sn[0]==GetAtomProperty(sn[1],A_ANTI))
					&& alg2_red_orth_2(&lm,l,ss[0]))
				goto bgn;
			if(no==2 && sn[0]!=sn[1] && alg2_red_orth_11(&lm,l,ss[0],ss[1]))
				goto bgn;
			
			if(no==3)
			{
				//if(alg2_red_orth_3(&lm,l,(Term *)ss,(Term *)sn))
					//goto bgn;
			}
			
			/*if(no==4)
			{
				if(alg2_red_orth_4(&lm,l,(Term *)ss,(Term *)sn))
					goto bgn;
			}*/	

			l=ListTail(l);
			}
	}
	
	SetCompoundArg(a2,5,lm);	
	
	
/*	WriteTerm(a2); puts(""); */

	}

Term ProcHermMatr(Term t)
{
	Term t1,t2,t3,tr=0,ti=0;
	List ll,l,l1;
	int dim, i, j;
	Label lab;

	if(CompoundArity(t)==3 && is_atom(CompoundArg1(t)) &&
			is_atom(CompoundArg2(t))&&is_integer(CompoundArgN(t,3)))
		{
		char cbuf1[128],cbuf2[128];
		int p1r=-1,p2r=-1,p1i=-1,p2i=-1;
		dim=(int)IntegerValue(CompoundArgN(t,3));
		if(dim<2 || dim>9)
			{
			ErrorInfo(200);
			puts("wrong dimension in OrthMatrix statement.");
			return 0;
			}
		strcpy(cbuf1,AtomValue(CompoundArg1(t)));
		strcpy(cbuf2,AtomValue(CompoundArg2(t)));
		for(i=0;cbuf1[i];i++)
			{
			if(cbuf1[i]=='_')
				{
				if(p1r==-1)
					p1r=i;
				else if(p2r==-1)
					p2r=i;
				else
					{
					ErrorInfo(201);
					puts("HermMatrix: parameter template must have two '_'.");
					return 0;
					}
				}
			}
		if(p1r==-1 || p2r==-1)
			{
			ErrorInfo(201);
			puts("HermMatrix: parameter template must have two '_'.");
			return 0;
			}
		for(i=0;cbuf2[i];i++)
			{
			if(cbuf2[i]=='_')
				{
				if(p1i==-1)
					p1i=i;
				else if(p2i==-1)
					p2i=i;
				else
					{
					ErrorInfo(201);
					puts("HermMatrix: parameter template must have two '_'.");
					return 0;
					}
				}
			}
		if(p1i==-1 || p2i==-1)
			{
			ErrorInfo(201);
			puts("HermMatrix: parameter template must have two '_'.");
			return 0;
			}
			
		for(i=0;i<dim;i++)
		for(j=0;j<dim;j++)
			{
			Atom a;
			cbuf1[p1r]='1'+(char)i;
			cbuf1[p2r]='1'+(char)j;
			a=NewAtom(cbuf1,0);
			if(!is_parameter(a))
				{
				ErrorInfo(202);
				printf("OrthMatrix: %s is undefined.\n",cbuf1);
				return 0;
				}
			cbuf2[p1i]='1'+(char)i;
			cbuf2[p2i]='1'+(char)j;
			a=NewAtom(cbuf2,0);
			if(!is_parameter(a))
				{
				ErrorInfo(202);
				printf("OrthMatrix: %s is undefined.\n",cbuf2);
				return 0;
				}
			}
			
		lab=NewLabel();
		t1=MakeCompound(A_I,dim);
		for(i=0;i<dim;i++)
			{
			t2=MakeCompound(A_I,dim);
			for(j=0;j<dim;j++)
				{
				Atom a;
				cbuf1[p1r]='1'+(char)i;
				cbuf1[p2r]='1'+(char)j;
				a=NewAtom(cbuf1,0);
				SetCompoundArg(t2,j+1,a);
				}
			SetCompoundArg(t1,i+1,t2);
			}
		ti=MakeCompound(A_I,dim);
		for(i=0;i<dim;i++)
			{
			t2=MakeCompound(A_I,dim);
			for(j=0;j<dim;j++)
				{
				Atom a;
				cbuf2[p1i]='1'+(char)i;
				cbuf2[p2i]='1'+(char)j;
				a=NewAtom(cbuf2,0);
				SetCompoundArg(t2,j+1,a);
				}
			SetCompoundArg(ti,i+1,t2);
			}
		goto setprp1;
		}

		
	t1=ConsumeCompoundArg(t,1);
	ti=ConsumeCompoundArg(t,2);
	FreeAtomic(t);
	t=t1;
	/*WriteTerm(t);puts("");*/
	if(!is_compound(t) || CompoundName(t)!=A_FBRACET || CompoundArity(t)!=1)
		{
		ErrorInfo(220);
		puts("Illegal argument in HermMatrix call");
		FreeAtomic(t);
		return 0;
		}
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	t=t1;
	lab=NewLabel();

	ll=CommaToList(t);
	/*WriteTerm(ll); puts(" : ll");*/
	dim=ListLength(ll);
	if(dim<2)
		{
		ErrorInfo(220);
		puts("Illegal argument in OrthMatrix call");
		FreeAtomic(ll);
		return 0;
		}
	t1=MakeCompound(A_I,dim);
	l=ll;
	i=0;
	while(!is_empty_list(l))
		{
		i++;
		t=ListFirst(l);
		if(!is_compound(t) || CompoundName(t)!=A_FBRACET || CompoundArity(t)!=1)
			{
			ErrorInfo(220);
			puts("Illegal argument in OrthMatrix call");
			FreeAtomic(ll);
			return 0;
			}
		t3=ConsumeCompoundArg(t,1);
		FreeAtomic(t);
		t=t3;
		t=CommaToList(t);
		if(ListLength(t)!=dim)
			{
			ErrorInfo(220);
			puts("Illegal argument in OrthMatrix call");
			FreeAtomic(ll);
			FreeAtomic(t1);
			return 0;
			}
		ChangeList(l,t);
		j=0;
		l1=t;
		t2=MakeCompound(A_I,dim);
		while(!is_empty_list(l1))
			{
			Term tt;
			j++;
			tt=ListFirst(l1);
			if(!is_parameter(tt))
				{
				ErrorInfo(220);
				WriteTerm(tt);
				puts(" undefined in HermMatrix call");
				FreeAtomic(ll);
				FreeAtomic(t1);
				FreeAtomic(t2);
				return 0;
				}
			SetCompoundArg(t2,j,tt);
			l1=ListTail(l1);
			}
		SetCompoundArg(t1,i,t2);
		l=ListTail(l);
		}
	FreeAtomic(ll);	

setprp1:
	for(i=1;i<=dim;i++)
	for(j=1;j<=dim;j++)
		{						
		Term tt;
		tt=CompoundArgN(CompoundArgN(t1,i),j);					
		t3=MakeCompound(A_I,4);
		SetCompoundArg(t3,1,lab);
		SetCompoundArg(t3,2,NewInteger(dim));
		SetCompoundArg(t3,3,NewInteger(i));
		SetCompoundArg(t3,4,NewInteger(j));
		SetAtomProperty(tt,A_HERM_MATR,t3);
		}
	tr=t1;
	
	if(is_compound(ti)&&CompoundName(ti)==A_I)
		{
		t1=ti;
		goto setprp2;
		}
	
	t=t1=ti;
	/*WriteTerm(t);puts("");*/
	if(!is_compound(t) || CompoundName(t)!=A_FBRACET || CompoundArity(t)!=1)
		{
		ErrorInfo(220);
		puts("Illegal argument in HermMatrix call");
		FreeAtomic(t);
		return 0;
		}
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	t=t1;

	ll=CommaToList(t);
	/*WriteTerm(ll); puts(" : ll");*/
	if(dim!=ListLength(ll))
		{
		ErrorInfo(220);
		puts("Illegal argument in OrthMatrix call");
		FreeAtomic(ll);
		return 0;
		}
	t1=MakeCompound(A_I,dim);
	l=ll;
	i=0;
	
	while(!is_empty_list(l))
		{
		i++;
		t=ListFirst(l);
		if(!is_compound(t) || CompoundName(t)!=A_FBRACET || CompoundArity(t)!=1)
			{
			ErrorInfo(220);
			puts("Illegal argument in OrthMatrix call");
			FreeAtomic(ll);
			return 0;
			}
		t3=ConsumeCompoundArg(t,1);
		FreeAtomic(t);
		t=t3;
		t=CommaToList(t);
		if(ListLength(t)!=dim)
			{
			ErrorInfo(220);
			puts("Illegal argument in OrthMatrix call");
			FreeAtomic(ll);
			FreeAtomic(t1);
			return 0;
			}
		ChangeList(l,t);
		j=0;
		l1=t;
		t2=MakeCompound(A_I,dim);
		while(!is_empty_list(l1))
			{
			Term tt;
			j++;
			tt=ListFirst(l1);
			if(!is_parameter(tt))
				{
				ErrorInfo(220);
				WriteTerm(tt);
				puts(" undefined in HermMatrix call");
				FreeAtomic(ll);
				FreeAtomic(t1);
				FreeAtomic(t2);
				return 0;
				}
			SetCompoundArg(t2,j,tt);
			l1=ListTail(l1);
			}
		SetCompoundArg(t1,i,t2);
		l=ListTail(l);
		}
	FreeAtomic(ll);	
	
setprp2:
	for(i=1;i<=dim;i++)
	for(j=1;j<=dim;j++)
		{						
		Term tt;
		tt=CompoundArgN(CompoundArgN(t1,i),j);					
		t3=MakeCompound(A_I,4);
		SetCompoundArg(t3,1,lab);
		SetCompoundArg(t3,2,NewInteger(-dim));
		SetCompoundArg(t3,3,NewInteger(i));
		SetCompoundArg(t3,4,NewInteger(j));
		SetAtomProperty(tt,A_HERM_MATR,t3);
		}
	/*WriteTerm(tr);puts("");
	WriteTerm(t1);puts("");*/
	
	hermi_matrices=AppendLast(hermi_matrices,MakeCompound3(OPR_DIV,lab,tr,t1));
/*	DumpList(hermi_matrices);*/
/*	DisplayTerm(t1); puts("");
	WriteTerm(herm_matrices); puts("");*/
	return 0;
}	

Term ProcRegMatrix(Term t, Term ind)
	{
	Term t1,t2,t3;
	List ll,l,l1;
	int dim, i, j;
	Label lab;
	
	if(CompoundArity(t)==3 && is_atom(CompoundArg1(t)) &&
			is_atom(CompoundArg2(t)) && is_integer(CompoundArgN(t,3)))
		return ProcHermMatr(t);
	
	if(CompoundArity(t)==2 && is_atom(CompoundArg1(t)) &&
			is_integer(CompoundArg2(t)))
		{
		char cbuf[128];
		int p1=-1,p2=-1;
		dim=(int)IntegerValue(CompoundArg2(t));
		if(dim<2 || dim>9)
			{
			ErrorInfo(200);
			puts("wrong dimension in OrthMatrix statement.");
			return 0;
			}
		strcpy(cbuf,AtomValue(CompoundArg1(t)));
		for(i=0;cbuf[i];i++)
			{
			if(cbuf[i]=='_')
				{
				if(p1==-1)
					p1=i;
				else if(p2==-1)
					p2=i;
				else
					{
					ErrorInfo(201);
					puts("OrthMatrix: parameter template must have two '_'.");
					return 0;
					}
				}
			}
		if(p1==-1 || p2==-1)
			{
			ErrorInfo(201);
			puts("OrthMatrix: parameter template must have two '_'.");
			return 0;
			}
		for(i=0;i<dim;i++)
		for(j=0;j<dim;j++)
			{
			Atom a;
			cbuf[p1]='1'+(char)i;
			cbuf[p2]='1'+(char)j;
			a=NewAtom(cbuf,0);
			if(!is_parameter(a))
				{
				ErrorInfo(202);
				printf("OrthMatrix: %s is undefined.\n",cbuf);
				return 0;
				}
			}
		lab=NewLabel();
		t1=MakeCompound(A_I,dim);
		for(i=0;i<dim;i++)
			{
			t2=MakeCompound(A_I,dim);
			for(j=0;j<dim;j++)
				{
				Atom a;
				cbuf[p1]='1'+(char)i;
				cbuf[p2]='1'+(char)j;
				a=NewAtom(cbuf,0);
				SetCompoundArg(t2,j+1,a);
				}
			SetCompoundArg(t1,i+1,t2);
			}
		goto setprp;
		}
				
		
				
	
	if(CompoundArity(t)==2)
		return ProcHermMatr(t);
	
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	t=t1;
	if(!is_compound(t) || CompoundName(t)!=A_FBRACET || CompoundArity(t)!=1)
		{
		ErrorInfo(220);
		puts("Illegal argument in OrthMatrix call");
		FreeAtomic(t);
		return 0;
		}
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	t=t1;
	lab=NewLabel();

	ll=CommaToList(t);
	/*WriteTerm(ll); puts(" : ll");*/
	dim=ListLength(ll);
	if(dim<2)
		{
		ErrorInfo(220);
		puts("Illegal argument in OrthMatrix call");
		FreeAtomic(ll);
		return 0;
		}
	t1=MakeCompound(A_I,dim);
	l=ll;
	i=0;
	while(!is_empty_list(l))
		{
		i++;
		t=ListFirst(l);
		if(!is_compound(t) || CompoundName(t)!=A_FBRACET || CompoundArity(t)!=1)
			{
			ErrorInfo(220);
			puts("Illegal argument in OrthMatrix call");
			FreeAtomic(ll);
			return 0;
			}
		t3=ConsumeCompoundArg(t,1);
		FreeAtomic(t);
		t=t3;
		t=CommaToList(t);
		if(ListLength(t)!=dim)
			{
			ErrorInfo(220);
			puts("Illegal argument in OrthMatrix call");
			FreeAtomic(ll);
			FreeAtomic(t1);
			return 0;
			}
		ChangeList(l,t);
		j=0;
		l1=t;
		t2=MakeCompound(A_I,dim);
		while(!is_empty_list(l1))
			{
			Term tt;
			j++;
			tt=ListFirst(l1);
			if(!is_parameter(tt))
				{
				ErrorInfo(220);
				WriteTerm(tt);
				puts(" undefined in OrthMatrix call");
				FreeAtomic(ll);
				FreeAtomic(t1);
				FreeAtomic(t2);
				return 0;
				}
			SetCompoundArg(t2,j,tt);
			l1=ListTail(l1);
			}
		SetCompoundArg(t1,i,t2);
		l=ListTail(l);
		}
	FreeAtomic(ll);	
	
setprp:
	for(i=1;i<=dim;i++)
	for(j=1;j<=dim;j++)
		{						
		Term tt, att;
		tt=CompoundArgN(CompoundArgN(t1,i),j);					
		t3=MakeCompound(A_I,4);
		SetCompoundArg(t3,1,lab);
		SetCompoundArg(t3,2,NewInteger(dim));
		SetCompoundArg(t3,3,NewInteger(i));
		SetCompoundArg(t3,4,NewInteger(j));
		SetAtomProperty(tt,A_ORTH_MATR,t3);
		att=GetAtomProperty(tt,A_ANTI);
		if(att)
			SetAtomProperty(att,A_ORTH_MATR,t3);
			
		}
		
	herm_matrices=AppendLast(herm_matrices,MakeCompound2(OPR_DIV,lab,t1));
/*	DumpList(herm_matrices);*/
/*	DisplayTerm(t1); puts("");
	WriteTerm(herm_matrices); puts("");*/
	return 0;
	
	}




	
static List rcos_mlt(List l, Atom par, int pw)
	{
	List l1;
	l1=l;
	while(!is_empty_list(l1))
		{
		Term pp;
		pp=ListFirst(l1);
		if(CompoundArg1(pp)==par)
			{
			int n;
			n=pw+(int)IntegerValue(CompoundArg2(pp));
			if(n)
				SetCompoundArg(pp,2,NewInteger(n));
			else
				l=CutFromList(l,l1);
			return l;
			}
		l1=ListTail(l1);
		}
	l=AppendLast(l,MakeCompound2(OPR_POW,par,NewInteger(pw)));
	l=SortedList(l,prtcmp);
	return l;
	} 		


static List rcos_add(List l, Term m2)
	{
	List ml=l;
/*	WriteTerm(l); printf("+"); WriteTerm(m2); puts("");*/
	while(!is_empty_list(l))
		{
		Term t;
		t=ListFirst(l);
		if(EqualTerms(CompoundArg2(t),CompoundArg2(m2)) &&
			EqualTerms(CompoundArgN(t,3),CompoundArgN(m2,3)))
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

void alg2_red_cos(Term a2)
	{
	List l1,l2,l3;
	int cfl;
beg:
	
	cfl=0;
	/*WriteTerm(a2); puts("");*/
	l1=l2=ConsumeCompoundArg(a2,5);
	l3=NewList();
	while(!is_empty_list(l1))
		{
		Term m1,m11;
		List lp;
		m1=ListFirst(l1);
		ChangeList(l1,0);
		lp=CompoundArg2(m1);
		while(!is_empty_list(lp))
			{
			Term ss,cs;
			int ipw;
			ipw=(int)IntegerValue(CompoundArg2(ListFirst(lp)));
			cs=CompoundArg1(ListFirst(lp));
			ss=GetAtomProperty(cs,A_COS);
			if(ss && ipw>1)
				{
				if(ipw==2)			
					{
					Term w1,w2;
					w1=ConsumeCompoundArg(m1,2);
					w1=CutFromList(w1,lp);
					m11=CopyTerm(m1);
					w2=CopyTerm(w1);
					SetCompoundArg(m1,2,w1);
					w1=w2;
					w1=rcos_mlt(w1,ss,2);
					SetCompoundArg(m11,2,w1);
					SetCompoundArg(m11,1,
					  NewInteger(-IntegerValue(CompoundArg1(m11))));
					l3=rcos_add(l3,m11);
					cfl=1;
					break;
					}
				if(ipw==3)		
					{
					Term w1,w2;
					w1=ConsumeCompoundArg(m1,2);
					w1=CutFromList(w1,lp);
					m11=CopyTerm(m1);
					w1=rcos_mlt(w1,cs,1);
					w2=CopyTerm(w1);
					SetCompoundArg(m1,2,w1);
					w1=w2;
					w1=rcos_mlt(w1,ss,2);
					SetCompoundArg(m11,2,w1);
					SetCompoundArg(m11,1,
					  NewInteger(-IntegerValue(CompoundArg1(m11))));
					l3=rcos_add(l3,m11);
					cfl=1;
					break;
					}
				if(ipw==4)
					{
					Term m12,w1,w2;
					w1=ConsumeCompoundArg(m1,2);
					w1=CutFromList(w1,lp);
					m11=CopyTerm(m1);
					m12=CopyTerm(m1);
					w2=CopyTerm(w1);
					SetCompoundArg(m1,2,w1);
					w1=CopyTerm(w2);
					w1=rcos_mlt(w1,ss,2);
					w2=rcos_mlt(w2,ss,4);
					SetCompoundArg(m11,2,w1);
					SetCompoundArg(m12,2,w2);
					SetCompoundArg(m11,1,
					  NewInteger(-2*IntegerValue(CompoundArg1(m11))));
					l3=rcos_add(l3,m11);
					l3=rcos_add(l3,m12);
					cfl=1;
					break;
					}
				if(ipw==5)
					{
					Term m12,w1,w2;
					w1=ConsumeCompoundArg(m1,2);
					w1=CutFromList(w1,lp);
					m11=CopyTerm(m1);
					m12=CopyTerm(m1);
					w1=rcos_mlt(w1,cs,1);
					w2=CopyTerm(w1);
					SetCompoundArg(m1,2,w1);
					w1=CopyTerm(w2);
					w1=rcos_mlt(w1,ss,2);
					w2=rcos_mlt(w2,ss,4);
					SetCompoundArg(m11,2,w1);
					SetCompoundArg(m12,2,w2);
					SetCompoundArg(m11,1,
					  NewInteger(-2*IntegerValue(CompoundArg1(m11))));
					l3=rcos_add(l3,m11);
					l3=rcos_add(l3,m12);
					cfl=1;
					break;
					}
				if(ipw==6)
					{
					Term m12,m13,w1,w2,w3;
					w1=ConsumeCompoundArg(m1,2);
					w1=CutFromList(w1,lp);
					m11=CopyTerm(m1);
					m12=CopyTerm(m1);
					m13=CopyTerm(m1);
					w2=CopyTerm(w1);
					SetCompoundArg(m1,2,w1);
					w1=CopyTerm(w2);
					w3=CopyTerm(w2);
					w1=rcos_mlt(w1,ss,2);
					w2=rcos_mlt(w2,ss,4);
					w3=rcos_mlt(w3,ss,6);
					SetCompoundArg(m11,2,w1);
					SetCompoundArg(m12,2,w2);
					SetCompoundArg(m13,2,w3);
					SetCompoundArg(m11,1,
					  NewInteger(-3*IntegerValue(CompoundArg1(m11))));
					SetCompoundArg(m12,1,
					  NewInteger(3*IntegerValue(CompoundArg1(m12))));
					SetCompoundArg(m13,1,
					  NewInteger(-IntegerValue(CompoundArg1(m13))));
					l3=rcos_add(l3,m11);
					l3=rcos_add(l3,m12);
					l3=rcos_add(l3,m13);
					cfl=1;
					break;
					}
				if(ipw==7)
					{
					Term m12,m13,w1,w2,w3;
					w1=ConsumeCompoundArg(m1,2);
					w1=CutFromList(w1,lp);
					m11=CopyTerm(m1);
					m12=CopyTerm(m1);
					m13=CopyTerm(m1);
					w1=rcos_mlt(w1,cs,1);
					w2=CopyTerm(w1);
					SetCompoundArg(m1,2,w1);
					w1=CopyTerm(w2);
					w3=CopyTerm(w2);
					w1=rcos_mlt(w1,ss,2);
					w2=rcos_mlt(w2,ss,4);
					w3=rcos_mlt(w3,ss,6);
					SetCompoundArg(m11,2,w1);
					SetCompoundArg(m12,2,w2);
					SetCompoundArg(m13,2,w3);
					SetCompoundArg(m11,1,
					  NewInteger(-3*IntegerValue(CompoundArg1(m11))));
					SetCompoundArg(m12,1,
					  NewInteger(3*IntegerValue(CompoundArg1(m12))));
					SetCompoundArg(m13,1,
					  NewInteger(-IntegerValue(CompoundArg1(m13))));
					l3=rcos_add(l3,m11);
					l3=rcos_add(l3,m12);
					l3=rcos_add(l3,m13);
					cfl=1;
					break;
					}
				if(ipw==8)
					{
					Term m12,m13,m14,w1,w2,w3,w4;
					w1=ConsumeCompoundArg(m1,2);
					w1=CutFromList(w1,lp);
					m11=CopyTerm(m1);
					m12=CopyTerm(m1);
					m13=CopyTerm(m1);
					m14=CopyTerm(m1);
					w2=CopyTerm(w1);
					SetCompoundArg(m1,2,w1);
					w1=CopyTerm(w2);
					w3=CopyTerm(w2);
					w4=CopyTerm(w2);
					w1=rcos_mlt(w1,ss,2);
					w2=rcos_mlt(w2,ss,4);
					w3=rcos_mlt(w3,ss,6);
					w4=rcos_mlt(w4,ss,8);
					SetCompoundArg(m11,2,w1);
					SetCompoundArg(m12,2,w2);
					SetCompoundArg(m13,2,w3);
					SetCompoundArg(m14,2,w4);
					SetCompoundArg(m11,1,
					  NewInteger(-4*IntegerValue(CompoundArg1(m11))));
					SetCompoundArg(m12,1,
					  NewInteger(6*IntegerValue(CompoundArg1(m12))));
					SetCompoundArg(m13,1,
					  NewInteger(-4*IntegerValue(CompoundArg1(m13))));
					l3=rcos_add(l3,m11);
					l3=rcos_add(l3,m12);
					l3=rcos_add(l3,m13);
					l3=rcos_add(l3,m14);
					cfl=1;
					break;
					}
				if(ipw>8)
					{
					
					Term m12,m13,w1,w2,w3;
					w1=ConsumeCompoundArg(m1,2);
					w1=CutFromList(w1,lp);
					m11=CopyTerm(m1);
					m12=CopyTerm(m1);
					m13=CopyTerm(m1);
					w1=rcos_mlt(w1,cs,ipw-6);
					w2=CopyTerm(w1);
					SetCompoundArg(m1,2,w1);
					w1=CopyTerm(w2);
					w3=CopyTerm(w2);
					w1=rcos_mlt(w1,ss,2);
					w2=rcos_mlt(w2,ss,4);
					w3=rcos_mlt(w3,ss,6);
					SetCompoundArg(m11,2,w1);
					SetCompoundArg(m12,2,w2);
					SetCompoundArg(m13,2,w3);
					SetCompoundArg(m11,1,
					  NewInteger(-3*IntegerValue(CompoundArg1(m11))));
					SetCompoundArg(m12,1,
					  NewInteger(3*IntegerValue(CompoundArg1(m12))));
					SetCompoundArg(m13,1,
					  NewInteger(-IntegerValue(CompoundArg1(m13))));
					l3=rcos_add(l3,m11);
					l3=rcos_add(l3,m12);
					l3=rcos_add(l3,m13);
					
					cfl=1;
					break;
					}
				printf("Warning: cosine in power %d\n",ipw);
				}
			lp=ListTail(lp);
			}
		l3=rcos_add(l3,m1);
		l1=ListTail(l1);
		}
	RemoveList(l2);
	SetCompoundArg(a2,5,l3);
	if(cfl)
		goto beg;
/*	WriteTerm(a2); puts(""); */
	}


