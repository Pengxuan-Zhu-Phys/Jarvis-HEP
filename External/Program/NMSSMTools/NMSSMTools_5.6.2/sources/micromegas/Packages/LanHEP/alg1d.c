#include <setjmp.h>
#include "lanhep.h"

extern jmp_buf alg1_jmp_buf;

void alg1_rem_inf(Term a1);
List a1l_rem_inf(List);

static List labl=0, labm=0;

static void lablist(Term t)
{
	if(is_list(t))
	{
		List l;
		for(l=t;l;l=ListTail(l))
		{
			Term m=ListFirst(l);
			if(is_label(m) && !ListMember(labl,m))
				labl=AppendFirst(labl,m);
			else
				lablist(m);
		}
	}
	if(is_compound(t))
	{
		int i;
		for(i=1;i<=CompoundArity(t);i++)
		{
			Term m=CompoundArgN(t,i);
			if(is_label(m) && !ListMember(labl,m))
				labl=AppendFirst(labl,m);
			else
				lablist(m);
		}
	}
}

static void chlabs(Term t)
{
	if(is_list(t))
	{
		List l;
		for(l=t;l;l=ListTail(l))
		{
			Term m=ListFirst(l);
			if(is_label(m))
				ChangeList(l,ListNth(labm,ListMember(labl,m)));
			else
				chlabs(m);
		}
	}
	if(is_compound(t))
	{
		int i;
		for(i=1;i<=CompoundArity(t);i++)
		{
			Term m=CompoundArgN(t,i);
			if(is_label(m))
				SetCompoundArg(t,i,ListNth(labm,ListMember(labl,m)));
			else
				chlabs(m);
		}
	}
}


void renewlab(Term t)
{
	List l;
	labl=0;
	labm=0;
	lablist(t);
	for(l=labl;l;l=ListTail(l))
		labm=AppendFirst(labm,NewLabel());
	chlabs(t);
	if(labl) {FreeAtomic(labl);FreeAtomic(labm);}
}

static void repl_ind(Term a2, List oi, List ni)
	{
	List l1,l2,l3,oii=0,nii=0;
	
	a2=CompoundArg1(a2);
	
	for(l1=a2;l1;l1=ListTail(l1))
		for(l2=CompoundArgN(ListFirst(l1),3);l2;l2=ListTail(l2))
			for(l3=CompoundArg1(ListFirst(l2));l3;l3=ListTail(l3))
			{
				Label la;
				la=CompoundArg2(ListFirst(l3));
				if(!ListMember(oi,la) && !ListMember(oii,la))
					oii=AppendLast(oii,la);
			}
	
	for(l1=oii;l1;l1=ListTail(l1))
		nii=AppendLast(nii,NewLabel());
	
	while(!is_empty_list(a2))
		{
		List t;
		t=CompoundArgN(ListFirst(a2),3);
		while(!is_empty_list(t))
			{
			List t1;
			t1=CompoundArg1(ListFirst(t));
			while(!is_empty_list(t1))
				{
				Term a;
				List loi,lni;
				a=CompoundArg2(ListFirst(t1));
				loi=oi;
				lni=ni;
				while(!is_empty_list(loi))
					{
					if(a==ListFirst(loi))
						{
						SetCompoundArg(ListFirst(t1),2,CompoundArg2(ListFirst(lni)));
						break;
						}
					loi=ListTail(loi);
					lni=ListTail(lni);
					}
				if(!loi)
					{
					loi=oii;
					lni=nii;
					while(!is_empty_list(loi))
						{
						if(a==ListFirst(loi))
							{
							SetCompoundArg(ListFirst(t1),2,ListFirst(lni));
							break;
							}
						loi=ListTail(loi);
						lni=ListTail(lni);
						}
					}
				t1=ListTail(t1);
				}
			t=ListTail(t);
			}
		a2=ListTail(a2);
		}
	RemoveList(oii);
	RemoveList(nii);
	return ;
	}			
		
extern int infi_order;

static List mk_let(Term m1, List cut, Term a1)
	{
	List l,l2;
	List lb,le;
	List sd;
	long int num1,den1;
	long int io1=0,io2=0,rmn=0;
	
	for(l=CompoundArgN(m1,3);l;l=ListTail(l))
	{
		Term prp=0;
		if(CompoundName(ListFirst(l))==OPR_PARAMETER)
				prp=GetAtomProperty(CompoundArg2(ListFirst(l)),A_INFINITESIMAL);
		if(prp && IntegerValue(CompoundArg1(prp))>0)
				io1+=IntegerValue(CompoundArg1(prp));
		if(CompoundName(ListFirst(l))==A_INFINITESIMAL)
				io1+=IntegerValue(CompoundArg2(ListFirst(l)));
	}
	
	l=ConsumeCompoundArg(a1,1);
	FreeAtomic(a1);
	a1=l;
	
	for(l2=a1;l2;l2=ListTail(l2))
	{
		Term m11=ListFirst(l2);
		io2=0;
		for(l=CompoundArgN(m11,3);l;l=ListTail(l))
		{
			Term prp=0;
			if(CompoundName(ListFirst(l))==OPR_PARAMETER)
					prp=GetAtomProperty(CompoundArg2(ListFirst(l)),A_INFINITESIMAL);
			if(prp && IntegerValue(CompoundArg1(prp))>0)
					io2+=IntegerValue(CompoundArg1(prp));
			if(CompoundName(ListFirst(l))==A_INFINITESIMAL)
					io2+=IntegerValue(CompoundArg2(ListFirst(l)));
		}
		if(io1+io2>infi_order)
		{
			FreeAtomic(ListFirst(l2));
			ChangeList(l2,0);
			rmn++;
		}
	}
	
	if(rmn)
		do
		{
			for(l=a1;l;l=ListTail(l))
				if(ListFirst(l)==0)
				{
					a1=CutFromList(a1,l);
					break;
				}
		}
		while(l);
	
	num1=IntegerValue(ConsumeCompoundArg(m1,1));
	den1=IntegerValue(ConsumeCompoundArg(m1,2));
	l=ConsumeCompoundArg(m1,3);
	sd=ConsumeCompoundArg(m1,4);
	lb=ListSplit(l,cut,&le);
	FreeAtomic(cut);
	FreeAtomic(m1);
	
	l=a1;
	while(!is_empty_list(l))
		{
		long int n1,n2,d1,d2,num,den,cf;
		List lb1,lb1e,le1,le1e,lm;
		m1=ListFirst(l);
		lm=ConsumeCompoundArg(m1,3);
		lb1=CopyList(lb,&lb1e);
		le1=CopyList(le,&le1e);
		if(lm==0)
			lm=le1;
		else if(le1)
			ConcatList(lm,le1);
		if(lb1==0)
			lb1=lm;
		else if(lm)
			ListConcat(lb1e,lm);
		/*lb1=ConcatList(lb1,lm);
		lb1=ConcatList(lb1,le1);*/
		SetCompoundArg(m1,3,lb1);
		
		lm=ConsumeCompoundArg(m1,4);
		lb1=CopyList(sd,&lb1e);
		if(lb1==0)
			lb1=lm;
		else if(lm)
			ListConcat(lb1e,lm);
		/*lm=ConcatList(lm,CopyTerm(sd));*/
		SetCompoundArg(m1,4,lb1);
		
		n1=num1;
		d1=den1;
		n2=IntegerValue(CompoundArg1(m1));
		d2=IntegerValue(CompoundArg2(m1));
		num=n1*n2;
		den=d1*d2;
		if(den<0) den=-den;
		cf=gcf(num,den);
		num/=cf;
		den/=cf;
		if(d1<0 && d2<0)
			{
			num=-num;
			}
		else
			if((d1<0 && d2>0) || (d1>0 && d2<0))
				{
				den=-den;
				}
		SetCompoundArg(m1,1,NewInteger(num));
		SetCompoundArg(m1,2,NewInteger(den));
		l=ListTail(l);
		}
	FreeAtomic(lb);
	FreeAtomic(le);
	FreeAtomic(sd);
	return a1;
	}



static List mk_let_d(Term m1,  Term a1)
	{
	List l,l2;
	List lb;
	List sd;
	long int num1,den1;
	long int io1=0,io2=0,rmn=0;
	
	for(l=CompoundArgN(m1,3);l;l=ListTail(l))
	{
		Term prp=0;
		if(CompoundName(ListFirst(l))==OPR_PARAMETER)
				prp=GetAtomProperty(CompoundArg2(ListFirst(l)),A_INFINITESIMAL);
		if(prp && IntegerValue(CompoundArg1(prp))>0)
				io1+=IntegerValue(CompoundArg1(prp));
		if(CompoundName(ListFirst(l))==A_INFINITESIMAL)
				io1+=IntegerValue(CompoundArg2(ListFirst(l)));
	}


	l=ConsumeCompoundArg(a1,1);
	FreeAtomic(a1);
	a1=l;


	for(l2=a1;l2;l2=ListTail(l2))
	{
		Term m11=ListFirst(l2);
		io2=0;
		for(l=CompoundArgN(m11,3);l;l=ListTail(l))
		{
			Term prp=0;
			if(CompoundName(ListFirst(l))==OPR_PARAMETER)
					prp=GetAtomProperty(CompoundArg2(ListFirst(l)),A_INFINITESIMAL);
			if(prp && IntegerValue(CompoundArg1(prp))>0)
					io2+=IntegerValue(CompoundArg1(prp));
			if(CompoundName(ListFirst(l))==A_INFINITESIMAL)
					io2+=IntegerValue(CompoundArg2(ListFirst(l)));
		}
		if(io1+io2>infi_order)
		{
			FreeAtomic(ListFirst(l2));
			ChangeList(l2,0);
			rmn++;
		}
	}
	
	if(rmn)
		do
		{
			for(l=a1;l;l=ListTail(l))
				if(ListFirst(l)==0)
				{
					a1=CutFromList(a1,l);
					break;
				}
		}
		while(l);

	
	num1=IntegerValue(ConsumeCompoundArg(m1,1));
	den1=IntegerValue(ConsumeCompoundArg(m1,2));
	l=ConsumeCompoundArg(m1,3);
	sd=ConsumeCompoundArg(m1,4);
	lb=l;
	FreeAtomic(m1);
	
	l=a1;
	while(!is_empty_list(l))
		{
		long int n1,n2,d1,d2,num,den,cf;
		List lb1,lm;
		m1=ListFirst(l);
		lm=ConsumeCompoundArg(m1,3);
		lb1=CopyTerm(lb);
		lb1=ConcatList(lb1,lm);
		SetCompoundArg(m1,3,lb1);
		
		lm=ConsumeCompoundArg(m1,4);
		lm=ConcatList(lm,CopyTerm(sd));
		SetCompoundArg(m1,4,lm);
		
		n1=num1;
		d1=den1;
		n2=IntegerValue(CompoundArg1(m1));
		d2=IntegerValue(CompoundArg2(m1));
		num=n1*n2;
		den=d1*d2;
		if(den<0) den=-den;
		cf=gcf(num,den);
		num/=cf;
		den/=cf;
		if(d1<0 && d2<0)
			{
			num=-num;
			}
		else
			if((d1<0 && d2>0) || (d1>0 && d2<0))
				{
				den=-den;
				}
		SetCompoundArg(m1,1,NewInteger(num));
		SetCompoundArg(m1,2,NewInteger(den));
		l=ListTail(l);
		}
	FreeAtomic(lb);
	FreeAtomic(sd);
	return a1;
	}


static List s_l_1(Term m1)
	{
	List l,l1;
	
	l=CompoundArgN(m1,3);
	while(!is_empty_list(l))
		{
		Term t1;
		t1=ListFirst(l);
		if(CompoundName(t1)==OPR_LET)
			{
			Term sub,a1,ila,ill;
			a1=sub=ConsumeCompoundArg(t1,2);
			ill=ConsumeCompoundArg(t1,1);
			FreeAtomic(t1);
			ChangeList(l,0);
			
			sub=CopyTerm(GetAtomProperty(a1,OPR_LET));
			if(sub==0 || !is_compound(sub) || CompoundName(sub)!=OPR_LET)
				{
				sub=CopyTerm(GetAtomProperty(a1,PROP_TYPE));
				if(sub==0 || !is_compound(sub) || CompoundName(sub)!=OPR_LET)
					{
					printf("Internal error: inconsistent substitution '");
					WriteTerm(CompoundArg2(a1));
					puts("'");
					longjmp(alg1_jmp_buf,1);
					}
				}
			renewlab(sub);
				
			a1=ConsumeCompoundArg(sub,1);
			ila=ConsumeCompoundArg(sub,2);
			FreeAtomic(sub);
			repl_ind(a1,ila,ill);
			
			return SetLets(mk_let(m1,l,a1));
			
			}
		l=ListTail(l);
		}
	
	
	l=CompoundArgN(m1,4);
	while(!is_empty_list(l))
		{
		Term t1;
		t1=ListFirst(l);
		if(CompoundName(t1)==OPR_LET)
			{
			Term sub,a1,ila,ill;
			a1=sub=ConsumeCompoundArg(t1,2);
			ill=ConsumeCompoundArg(t1,1);
			if(ill)
				{
				puts("Integnal error (iinld)");
				longjmp(alg1_jmp_buf,1);
				}
			
			l1=ConsumeCompoundArg(m1,4);
			l1=CutFromList(l1,l);
			SetCompoundArg(m1,4,l1);
			sub=CopyTerm(GetAtomProperty(a1,OPR_LET));
			if(sub==0 || !is_compound(sub) || CompoundName(sub)!=OPR_LET)
				{
				sub=CopyTerm(GetAtomProperty(a1,PROP_TYPE));
				if(sub==0 || !is_compound(sub) || CompoundName(sub)!=OPR_LET)
					{
					printf("Internal error: inconsistent substitution '");
					WriteTerm(CompoundArg2(a1));
					puts("'");
					longjmp(alg1_jmp_buf,1);
					}
				}
			renewlab(sub);
				
			if(CompoundArgN(sub,3)==0)
				{
				ErrorInfo(378);
				printf("symbol '%s' can not be in denominator\n",AtomValue(a1));
				longjmp(alg1_jmp_buf,1);
				}

			a1=ConsumeCompoundArg(sub,3);
			ila=ConsumeCompoundArg(sub,2);
			FreeAtomic(sub);
			if(ila)
				{
				puts("Internal error (iidl1)");
				longjmp(alg1_jmp_buf,1);
				}
				
			
			return SetLets(mk_let_d(m1,a1));
			
			}
		l=ListTail(l);
		}
		
	

	
	return AppendFirst(NewList(),m1);
	}
			
			
			
static int slrl=0;

List SetLets(List l)
	{
	List l1,lr, lre=0;
	lr=NewList();
	slrl++;
	l1=l;
	while(!is_empty_list(l1))
		{
		List nl=s_l_1(ListFirst(l1));
		if(nl==0)
			{l1=ListTail(l1);continue;}
		if(lr==0)
			{lr=nl;lre=lr;}
		else
			ConcatList(lre,nl);
		while(ListTail(lre)) lre=ListTail(lre);
		/*lr=ConcatList(lr,);*/
		l1=ListTail(l1);
		}
	RemoveList(l);
	
	slrl--;
	if(slrl==0)
	{
		for(l1=lr;l1;l1=ListTail(l1))
		{
			List l2,sl=ConsumeCompoundArg(ListFirst(l1),3);
			int ch=0;
			for(l2=sl;l2;l2=ListTail(l2))
				if(CompoundName(ListFirst(l2))==A_INFINITESIMAL)
				{
					FreeAtomic(ListFirst(l2));
					ChangeList(l2,0);
					ch++;
				}
			if(ch)
			do
			{
				for(l2=sl;l2;l2=ListTail(l2))
				if(ListFirst(l2)==0)
				{
					sl=CutFromList(sl,l2);
					break;
				}
			} while(l2);
			SetCompoundArg(ListFirst(l1),3,sl);
		}
	}
						
		
	lr=a1l_rem_inf(lr);
	return lr;
	}
	
	
