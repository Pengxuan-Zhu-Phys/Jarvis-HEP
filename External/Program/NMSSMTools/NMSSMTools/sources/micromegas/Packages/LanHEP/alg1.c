#include <setjmp.h>
#include "lanhep.h"
#include <time.h>


jmp_buf alg1_jmp_buf;

List DefaultIndex = 0;

extern int opSplitCol1;
int do_brst=0;

void alg1_derivp(Term a1);
void proc_hash(Term );
Term alg1_proc_brst(Term a1, Atom tp);
	
Term GetIndices(Term t, Term ownind)  /* Interface */
	{
	List l;
	Term t1;
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	t=t1;
	if(!is_compound(t) || (CompoundName(t)!=OPR_USCORE && CompoundName(t)!=OPR_CARET))
		return 0;		
	if(SplitIndices(t,&l)==0)
		return 0;
	else
		return l;
	}
	


int alg1_recurse_level = 0;

Term ExprTo1(Term t)
	{
	Term ret, oind,nind;
	int err;
	time_t ti[8];
	if(alg1_recurse_level==0)
		{
		err=setjmp(alg1_jmp_buf);
		if(err==1)
			{
			alg1_recurse_level=0;
			return 0;
			}
		}
	alg1_recurse_level++;
		
ti[0]=clock();

	ret=WheredTerm(t);   	/* alg1a.c */
	proc_hash(ret);

ti[1]=clock();

	ret=ExpandTerm(ret); 	/* alg1b.c */
	
ti[2]=clock();


	ret=SetInd1(ret,&oind,&nind); /* alg1c.c */
	
ti[3]=clock();

/*
if(GetAtomProperty(NewAtom("qq",0),NewAtom("qq",0)))
{
printf("to sl: %d terms\n",ListLength(ret));
DumpList(ret);
}
*/
		
	ret=SetIntAlgs(ret);

	ret=SetLets(ret);		/* alg1d.c */

	
/*
if(GetAtomProperty(NewAtom("qq",0),NewAtom("qq",0)))
{
printf("to wf: %d terms\n",ListLength(ret));
DumpList(ret);
AtomStat1();ListStat1();
exit(0);
}
*/
	
	/*printf("Free: "); WriteTerm(oind); puts(""); */

	ret=MakeCompound2(A_ALG1,ret,oind);

    alg1_fix_delta(ret);

ti[4]=clock();	

    alg1_rem_inf(ret);


	alg1_fix_wild(ret);
	
ti[5]=clock();	

	alg1_exp_wild(ret,nind);

	alg1_sum_wild(ret);

	
ti[6]=clock();	
	
	alg1_rem_inf(ret);

ti[7]=clock();	
	
	

/*	printf("e21 stages: %d %d %d %d %d %d %d\n",
			ti[1]-ti[0],ti[2]-ti[1],ti[3]-ti[2],ti[4]-ti[3],
			ti[5]-ti[4],ti[6]-ti[5],ti[7]-ti[6]);
*/
		
	alg1_recurse_level--;
	return ret;
	}

void alg1_kl_to_ia(List);
void alg1_rem_c4(Term);
void alg1_rem_sincos(Term), alg1_let5th(Term);

Term ExprTo1kl(Term t)
	{
	Term ret, oind, nind;
	List l;
	
	int err;

	if(alg1_recurse_level==0)
		{
		err=setjmp(alg1_jmp_buf);
		if(err==1)
			{
			alg1_recurse_level=0;
			return 0;
			}
		}
	alg1_recurse_level++;
	

	ret=WheredTerm(t);   	/* alg1a.c */
	proc_hash(ret);

	ret=ExpandTerm(ret); 	/* alg1b.c */
	

	ret=SetInd1(ret,&oind,&nind); /* alg1c.c */
	
	alg1_kl_to_ia(ret);
	
	ret=SetIntAlgs(ret);
	
	for(l=ret;l;l=ListTail(l))
	{
		List l1;
		for(l1=CompoundArgN(ListFirst(l),3);l1;l1=ListTail(l1))
			if(CompoundName(ListFirst(l1))==OPR_WILD)
		{
			ErrorInfo(345);
			printf("array object is not allowed in object in 'keep_lets'\n");
			return 0;
		}
	}
	
	ret=MakeCompound2(A_ALG1,ret,oind);
		
	alg1_recurse_level--;
	
	return ret;
	}

List fromdfdfc=0;
int allow_dfdfc=0;

Term ExprTo11(Term t, List *nind_out)
	{
	Term ret, oind,nind;
	int err;
	
	if(alg1_recurse_level==0)
		{
		err=setjmp(alg1_jmp_buf);
		if(err==1)
			{
			alg1_recurse_level=0;
			return 0;
			}
		}
	alg1_recurse_level++;
		
	ret=WheredTerm(t);   	/* alg1a.c */
	proc_hash(ret);

	fromdfdfc=0;
	allow_dfdfc=1;
	ret=ExpandTerm(ret); 	/* alg1b.c */
	allow_dfdfc=0;
	if(fromdfdfc)
		ret=ConcatList(ret,fromdfdfc);
	
	
	ret=SetInd1(ret,&oind,&nind); /* alg1c.c */
	
	if(!is_empty_list(oind))
		{
		ErrorInfo(324);
		printf("non-scalar lagrangian term\n");
		/*DumpList(ret);
		WriteTerm(oind);puts("");*/
		return 0;
		}
	
	ret=alg1_spl_col(ret);
		
	ret=SetIntAlgs(ret);
	
	alg1_recurse_level--;
	
	*nind_out=nind;
	return ret;
	}

extern time_t  tm_12, tm_sl, tm_sw, tm_rsc;
	
Term ExprTo12(Term ret, List nind)
{
	/*WriteTerm(ret);puts("");*/
	time_t st1, st2;
	st1=clock();
	
	alg1_let5th(ret);
	st2=clock();
	ret=SetLets(ret);		/* alg1d.c */
	/*WriteTerm(ret);puts("");*/
	tm_sl+=(clock()-st2);
	
	ret=MakeCompound2(A_ALG1,ret,NewList());
	
    alg1_fix_delta(ret);

	alg1_fix_wild(ret);
	
	alg1_exp_wild(ret,nind);


/*	printf("eto12: sw: %d -> ",ListLength(CompoundArg1(ret)));fflush(stdout);*/
	st2=clock();
	alg1_sum_wild(ret);
	tm_sw+=(clock()-st2);
/*    printf(" %d\n",ListLength(CompoundArg1(ret)));*/
	
	/*st2=clock();	*/
	alg1_rem_inf(ret);
	/*tm_rsc+=(clock()-st2);	*/
	alg1_rem_sincos(ret);
	
	if(do_brst)
	{
		ret=alg1_proc_brst(ret,A_BRST_TRANSF);
	}
	
	alg1_derivp(ret);

	
	if(opSplitCol1==-1)
		alg1_rem_c4(ret);
	
	tm_12+=(clock()-st1);

	return ret;
	}

extern int opDoSymmetrize, write_all_vertices;
		
Term ProcCheckBRST(Term t, Term ind)
{
	do_brst=1;
	opDoSymmetrize=0;
	write_all_vertices=1;
	return 0;
}
	
Term To_t1(Term t, Term ind)     /* Interface */
	{
	Term tt;
	tt=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	tt=ExprTo1(tt);
/*	WriteTerm(CompoundArg2(tt));
	puts(" : indices; alg = ");
	DumpList(CompoundArg1(tt));*/
	if(tt)
		alg1_dump(tt);
	return Alg1ToExpr(tt);
	}
	
Term GetDefIndex(Term t, Term ind)
	{
	if(is_compound(t))
		FreeAtomic(t);
	return CopyTerm(DefaultIndex);
	}
	
Term SetDefIndex(Term t, Term ind)
	{
	int ar,i;
	FreeAtomic(DefaultIndex);
	DefaultIndex=NewList();
	if(!is_compound(t))
		{ 
		return 0;	
		}	
	ar=CompoundArity(t);
	for(i=1;i<=ar;i++)
		{
		Term cc,g1;
		cc=ConsumeCompoundArg(t,i);
		g1=SpecToRepr(cc);
		if(g1!=0)
			DefaultIndex=AppendLast(DefaultIndex,g1);
		}
	FreeAtomic(t);
	return CopyTerm(DefaultIndex);
	}

			
Term ProcVEV(Term t, Term ind)
	{
	Term t1;
	if(!is_compound(t) || CompoundArity(t)!=1)
		{
		ErrorInfo(215);
		puts("Illegal argument in 'vev' call");
		return 0;
		}
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	return MakeCompound2(OPR_MLT,t1,A_VEV);
	}
	

void alg1_dump(Term t)
{
	List l1;
	
	printf("------------------------------------\n");
	printf("alg1: ");WriteTerm(CompoundArg2(t));puts("");
	for(l1=CompoundArg1(t);l1;l1=ListTail(l1))
	{
		Term m;
		m=ListFirst(l1);
		printf("mterm %ld/%ld\n",
				IntegerValue(CompoundArg1(m)),IntegerValue(CompoundArg2(m)));
		DumpList(CompoundArgN(m,3));
		printf("\t/");WriteTerm(CompoundArgN(m,4));puts("");
	}
	printf("------------------------------------\n");
}


