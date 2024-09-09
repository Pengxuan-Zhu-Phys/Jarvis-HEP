#include <stdio.h>
#include "lanhep.h"

void alg2_norm(Term);
List alg2_denorm(Term);
		
void SaveRules(char *fname)
	{
	List l,l1,li,lj;
	{
	int f;
	f=itrSetOut(fname);
	if(f==0)
		{
		printf("Can not open file '%s' for writing to save Feynman rules\n",
			fname);
		return;
		}
	}
	itrOut(NewAtom("parameters",0));
	l=all_param_list();
	while(!is_empty_list(l))
		{
		itrOut(ListFirst(l));
		l=ListTail(l);
		}
	itrOut(A_END);
	l=all_prtc_list();
	while(!is_empty_list(l))
		{
		Atom pp;
		List pl;
		pp=ListFirst(l);
		pl=GetProperties(MakeCompound1(A_I,pp),0);
		itrOut(MakeCompound1(OPR_PARTICLE,pp));
		itrOut(pl);
		FreeAtomic(pl);
		l=ListTail(l);
		}
	
	l1=l=all_vert_list();
	
	for(li=l;!is_empty_list(li);li=ListTail(li))
	{
		Term a2;
		List a2l;
		int sf=0;
		
		a2=ListFirst(li);
		if(CompoundArgN(a2,5)==0)
			continue;

/*		if(need_col_rdc(a2))
		{
			List l1,l2,l3=0;
			sf=1;
			
			l1=alg2_denorm(CopyTerm(a2));
			
			for(l2=l1;l2;l2=ListTail(l2))
				l3=ConcatList(l3,color_reduce(ListFirst(l2)));
			RemoveList(l1);
			a2l=l3;
			for(l1=a2l;l1;l1=ListTail(l1))
				alg2_norm(ListFirst(l1));
		}
		else
			a2l=MakeList1(a2);
		
		WriteTerm(a2l); puts("");	
		for(lj=a2l;lj;lj=ListTail(lj))
		{
			Term a2;
			a2=ListFirst(lj);
*/			if(is_atom(CompoundArg1(a2)))
				continue;

			alg2_symmetrize(a2);
			alg2_common_s(a2);
			alg2_common_n(a2);
			/*alg2_red_cos(a2);
			alg2_red_orth(a2);
			alg2_red_sico(a2);
			alg2_red_comsico(a2);*/
			alg2_red_1pm5(a2);
			alg2_recommon_n(a2);
			alg2_recommon_s(a2);
		
			if(CompoundArg2(a2)!=NewInteger(0) &&
			 !is_empty_list(CompoundArgN(a2,5)))
				{
				itrOut(a2);
				}
			
/*		if(sf)
			FreeAtomic(a2l);
		else
			RemoveList(a2l);

		}*/
	}
		
	FreeAtomic(l1);
	itrCloseOut();
	
	}
