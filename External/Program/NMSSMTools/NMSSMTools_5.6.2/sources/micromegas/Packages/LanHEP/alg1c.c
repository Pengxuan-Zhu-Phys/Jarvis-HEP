#include <setjmp.h>
#include "lanhep.h"

extern jmp_buf alg1_jmp_buf;

static void set_ind(Term t, List *ii1, List *ni)
	{
	List li,l1,ii;


	l1=CompoundArgN(t,3);
	ii=li=NewList();
	*ni=NewList();
/*	WriteTerm(l1);puts("");*/
	while(!is_empty_list(l1))
		{
		List l2;
		l2=CompoundArg1(ListFirst(l1));
		if(!is_empty_list(l2))
			{
			while(!is_empty_list(l2))
				{
				if(is_integer(CompoundArg2(ListFirst(l2))))
					{
					long int iv;
					Label lb;
					iv=IntegerValue(CompoundArg2(ListFirst(l2)));
					if(CompoundArity(ListFirst(l2))==2)
						{
						ErrorInfo(317);
						printf("index of '");
						WriteTerm(CompoundArg2(ListFirst(l1)));
						puts("' can not be number.");
						longjmp(alg1_jmp_buf,1);
						}
					if(iv<1 || iv>IntegerValue(CompoundArgN(ListFirst(l2),3)))
						{
						ErrorInfo(317);
						printf("index of '");
						WriteTerm(CompoundArg2(ListFirst(l1)));
						puts("' out of range.");
						longjmp(alg1_jmp_buf,1);
						}
					lb=NewLabel();
					SetCompoundArg(ListFirst(l2),2,lb);
					*ni=AppendLast(*ni,
							MakeCompound2(OPR_EQSIGN,lb,NewInteger(iv)));
					}
					else
						li=AppendLast(li,ListFirst(l2));
				l2=ListTail(l2);
				}
			li=AppendLast(li,0);
			}
		l1=ListTail(l1);
		}

/*	WriteTerm(li);WriteTerm(*ni);puts("");*/
				
	l1=li;
	while(!is_empty_list(l1))
		{
		Term ct;
		List l2;
		ct=ListFirst(l1);
		if(ct==0 || is_label(CompoundArg2(ct)))
			{
			l1=ListTail(l1);
			continue;
			}
		l2=ListTail(l1);
		if(CompoundArg2(ct)==0)
			while(!is_empty_list(l2) && ListFirst(l2)!=0)
				l2=ListTail(l2);
		while(!is_empty_list(l2))
			{
			Term ct1;
			ct1=ListFirst(l2);
			if(ct1==0)
				{
				l2=ListTail(l2);
				continue;
				}
			if(CompoundArg2(ct)==CompoundArg2(ct1) && equal_groups(ct,ct1))
				{
				Term lab;
				lab=NewLabel();
				SetCompoundArg(ct,2,lab);
				SetCompoundArg(ct1,2,lab);
				break;
				}
			if(CompoundArg2(ct)==CompoundArg2(ct1) && CompoundArg2(ct)!=0)
				{
				Term lab;
				ErrorInfo(311);
				printf(" index %s refers to different groups\n",
					AtomValue(CompoundArg2(ct)));
				longjmp(alg1_jmp_buf,1);
				lab=NewLabel();
				SetCompoundArg(ct,2,lab);
				SetCompoundArg(ct1,2,lab);
				break;
				}
			l2=ListTail(l2);
			}
			
		if(is_empty_list(l2))
			{
			if(CompoundArg2(ct)==0)
				SetCompoundArg(ct,2,NewLabel());
			ii=AppendLast(ii,ct);
			}
				
		l1=ListTail(l1);
		}			
	l1=li;
	while(!is_empty_list(l1))
		{
		ChangeList(l1,0);
		l1=ListTail(l1);
		}
	FreeAtomic(li);
	*ii1=ii;
	
	}


static List resort_ind(List ii, List oind)
	{
	List ret;
	ret=NewList();
	while(!is_empty_list(oind))
		{
		Term t;
		t=ListFirst(oind);
		if(is_atom(CompoundArg2(t)))
			{
			List l;
			l=ii;
			while(!is_empty_list(l))
				{
				Term t1;
				t1=ListFirst(l);
				if(CompoundArg2(t1)==CompoundArg2(t))
					{
					if(!EqualTerms(CompoundArg1(t),CompoundArg1(t1)))
						{
						ErrorInfo(312);
						printf(" misused index '%s'\n",AtomValue(CompoundArg2(t1)));
						longjmp(alg1_jmp_buf,1);
						}
					ChangeList(l,0);
					ret=AppendLast(ret,t1);
					ii=CutFromList(ii,l);
					goto fi;
					}
				l=ListTail(l);
				}
			return 0;
			}
		else
			{
			List l;
			l=ii;
			while(!is_empty_list(l))
				{
				Term t1;
				t1=ListFirst(l);
				if(is_label(CompoundArg2(t1)) &&
					EqualTerms(CompoundArg1(t),CompoundArg1(t1)))
						{
						ChangeList(l,0);
						ret=AppendLast(ret,t1);
						ii=CutFromList(ii,l);
						goto fi;
						}
				l=ListTail(l);
				}
			return 0;
			}
			
		fi:	
		oind=ListTail(oind);
		}
		
	return ret;
	}	
				
		



static void uni_ind(List ii, List oind)
	{
	while(!is_empty_list(ii))
		{
		SetCompoundArg(ListFirst(ii),2,CompoundArg2(ListFirst(oind)));
		ii=ListTail(ii);
		oind=ListTail(oind);
		}
	}

Term SetInd1(Term t, List *ir, List *ni)
	{
	List li,oind,nind;
	int first=1, fil=0;
	li=t;

		
	oind=NewList();
	nind=NewList();
	while(!is_empty_list(li))
		{
		Term ct,ii,jj;
		ct=ListFirst(li);
		ii=NewList();
		set_ind(ct,&ii,&jj);
		nind=ConcatList(nind,jj);
		if(first)
			{
/*			printf("Free: "); WriteTerm(ii); puts(""); */
			oind=CopyTerm(ii);
			fil=ListLength(oind);
			first=0;
			}
		else
			{
			if(ListLength(ii)!=fil)
				{
				ErrorInfo(313);
				printf("terms in expression have different indices.\n");
				longjmp(alg1_jmp_buf,1);
				}
			if(fil!=0)
				{
				ii=resort_ind(ii,oind);
				if(ii==0)
					{
					ErrorInfo(314);
					printf("terms in expression have different indices.\n");
					longjmp(alg1_jmp_buf,1);
					}
				uni_ind(ii,oind);
				}
			}
			
		li=ListTail(li);
		
		}
		
	*ir=oind;
	*ni=nind;
	return t;
	}
	

	
void alg1_cncl(Term m1)
{
	List l1,l2;
rpt:
	for(l1=CompoundArgN(m1,3);l1;l1=ListTail(l1))
	{
		Term sp=ListFirst(l1);
		if((CompoundName(sp)==OPR_PARAMETER || CompoundName(sp)==OPR_LET)
				&& CompoundArg1(sp)==0)
		{
			for(l2=CompoundArgN(m1,4);l2;l2=ListTail(l2))
			{
				Term sp2=ListFirst(l2);
				if(CompoundName(sp2)==CompoundName(sp) && CompoundArg1(sp2)==0
						&& CompoundArg2(sp)==CompoundArg2(sp2))
				{
					List tmp;
					tmp=ConsumeCompoundArg(m1,3);
					tmp=CutFromList(tmp,l1);
					SetCompoundArg(m1,3,tmp);
					tmp=ConsumeCompoundArg(m1,4);
					tmp=CutFromList(tmp,l2);
					SetCompoundArg(m1,4,tmp);
					goto rpt;
				}
			}
		}
	}
}

	
