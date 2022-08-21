#include <setjmp.h>
#include "lanhep.h"

extern jmp_buf alg1_jmp_buf;

void sub_ind_alg(Term a1, Atomic ind, int i); 
		
static List mk_sub1(Term m1, List cut, Term a1)
	{
	List l;
	List lb,le;
	List sd;
	long int num1,den1;
	MTHOPTDECL
	
	l=ConsumeCompoundArg(a1,1);
	FreeAtomic2(a1,mind);
	a1=l;
	
	num1=IntegerValue(ConsumeCompoundArg(m1,1));
	den1=IntegerValue(ConsumeCompoundArg(m1,2));
	l=ConsumeCompoundArg(m1,3);
	sd=ConsumeCompoundArg(m1,4);
	lb=ListSplit(l,cut,&le);
	FreeAtomic2(cut,mind);
	FreeAtomic2(m1,mind);
	
	
	
	l=a1;
	while(!is_empty_list(l))
		{
		long int n1,n2,d1,d2,num,den,cf;
		List lb1,le1,lm;
		m1=ListFirst(l);
		lm=ConsumeCompoundArg(m1,3);
		lb1=CopyTerm2(lb,mind);
		le1=CopyTerm2(le,mind);
		lb1=ConcatList(lb1,lm);
		lb1=ConcatList(lb1,le1);
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
	FreeAtomic2(lb,mind);
	FreeAtomic2(le,mind);
	FreeAtomic2(sd,mind);
	return a1;
	}

/*
static List mk_sub2(Term m1, List cut1, Term a11, List cut2, Term a12)
	{
	List l,l1,ret;
	List lb,le,lm;
	List sd;
	int num1,den1;
	
	
	num1=IntegerValue(ConsumeCompoundArg(m1,1));
	den1=IntegerValue(ConsumeCompoundArg(m1,2));
	l=ConsumeCompoundArg(m1,3);
	sd=ConsumeCompoundArg(m1,4);
	lb=ListSplit(l,cut1,&lm);
	lm=ListSplit(lm,cut2,&le);
	FreeAtomic(cut1);
	FreeAtomic(cut2);
	FreeAtomic(m1);
	
	ret=NewList();
	l=CompoundArg1(a11);	
	while(!is_empty_list(l))
		{
		l1=CompoundArg1(a12);
		while(!is_empty_list(l1))
			{
			int n1,n2,d1,d2,num,den,cf;
			List lb1,le1,lm1;
			Term m2,mr;
			m1=ListFirst(l);
			m2=ListFirst(l1);
			
			mr=MakeCompound(A_MTERM,4);
			lb1=CopyTerm(lb);
			le1=CopyTerm(le);
			lm1=CopyTerm(lm);
			lb1=ConcatList(lb1,CopyTerm(CompoundArgN(m1,3)));
			lb1=ConcatList(lb1,lm1);
			lb1=ConcatList(lb1,CopyTerm(CompoundArgN(m2,3)));
			lb1=ConcatList(lb1,le1);
			SetCompoundArg(mr,3,lb1);
			lb1=CopyTerm(sd);
			lb1=ConcatList(lb1,CopyTerm(CompoundArgN(m1,4)));
			lb1=ConcatList(lb1,CopyTerm(CompoundArgN(m2,4)));
			SetCompoundArg(mr,4,lb1);
		
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
			n1=num;
			d1=den;
			n2=IntegerValue(CompoundArg1(m2));
			d2=IntegerValue(CompoundArg2(m2));
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
			SetCompoundArg(mr,1,NewInteger(num));
			SetCompoundArg(mr,2,NewInteger(den));
			ret=AppendLast(ret,mr);
			
			l1=ListTail(l1);
			}
		l=ListTail(l);
		}
	FreeAtomic(lb);
	FreeAtomic(le);
	FreeAtomic(lm);
	FreeAtomic(sd);
	return ret;
	}
*/

		
static List mk_sub2(Term m1, List cut1, Term a11, List cut2, Term a12)
	{
	List l,l1,ret;
	List lb,le,lm;
	List sd;
	long int num1,den1;
	
	MTHOPTDECL
	
	num1=IntegerValue(ConsumeCompoundArg(m1,1));
	den1=IntegerValue(ConsumeCompoundArg(m1,2));
	l=ConsumeCompoundArg(m1,3);
	sd=ConsumeCompoundArg(m1,4);
	lb=ListSplit(l,cut1,&lm);
	lm=ListSplit(lm,cut2,&le);
	FreeAtomic2(cut1,mind);
	FreeAtomic2(cut2,mind);
	FreeAtomic2(m1,mind);
	
	ret=NewList();
	l=CompoundArg1(a11);	
	while(!is_empty_list(l))
		{
		l1=CompoundArg1(a12);
		while(!is_empty_list(l1))
			{
			long int n1,n2,d1,d2,num,den,cf;
			List lb1,lb1e,le1,le1e,lm1,lm1e,m1e,m2e;
			List base=0, end=0;
			Term m2,mr;
			m1=ListFirst(l);
			m2=ListFirst(l1);
			
			
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
			n1=num;
			d1=den;
			n2=IntegerValue(CompoundArg1(m2));
			d2=IntegerValue(CompoundArg2(m2));
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
			if(num==0)
			{
				l1=ListTail(l1);
				continue;
			}
			
			mr=MakeCompound(A_MTERM,4);
			SetCompoundArg(mr,1,NewInteger(num));
			SetCompoundArg(mr,2,NewInteger(den));
			
			
			lb1=CopyList2(lb,&lb1e,mind);
			le1=CopyList2(le,&le1e,mind);
			lm1=CopyList2(lm,&lm1e,mind);
			m1=CopyList2(CompoundArgN(m1,3),&m1e,mind);
			m2=CopyList2(CompoundArgN(m2,3),&m2e,mind);
			
			base=lb1; end=lb1e;
			if(m1)
			{
				if(base==0) base=m1;
				else ListConcat(end,m1); 
				end=m1e;
			}
			if(lm1)	
			{
				if(base==0) base=lm1;
				else ListConcat(end,lm1); 
				end=lm1e;
			}
			if(m2)
			{
				if(base==0) base=m2;
				else ListConcat(end,m2); 
				end=m2e;
			}
			if(le1)	
			{
				if(base==0) base=le1;
				else ListConcat(end,le1); 
				end=le1e;
			}
			
			SetCompoundArg(mr,3,base);
			lb1=CopyList2(sd,&lb1e,mind);
			lm1=CopyList2(CompoundArgN(ListFirst(l),4),&lm1e,mind);
			le1=CopyList2(CompoundArgN(ListFirst(l1),4),&le1e,mind);
			
			base=lb1; end=lb1e;
			if(lm1)	
			{
				if(base==0) base=lm1;
				else ListConcat(end,lm1); 
				end=lm1e;
			}
			if(le1)	
			{
				if(base==0) base=le1;
				else ListConcat(end,le1); 
				end=le1e;
			}
			
			SetCompoundArg(mr,4,base);
		
			ret=AppendLast(ret,mr);
			
			l1=ListTail(l1);
			}
		l=ListTail(l);
		}
	FreeAtomic2(lb,mind);
	FreeAtomic2(le,mind);
	FreeAtomic2(lm,mind);
	FreeAtomic2(sd,mind);
	return ret;
	}

static Term sub_ind_wild(Term w, int ino, int i)
	{
	List al,il,ilc;
	Term ret;
	MTHOPTDECL
	
	al=ConsumeCompoundArg(w,2);
	il=ConsumeCompoundArg(w,1);
	
	if(is_empty_list(al))
		{
		Atomic ind;	
		if(ListLength(il)==1)
			{
			FreeAtomic2(il,mind);
			FreeAtomic2(w,mind);
			return MakeCompound2(A_ALG1,0,0);
			}
			
		ilc=ListNthList(il,ino);
		ind=CompoundArg2(ListFirst(ilc));
		il=CutFromList(il,ilc);
		
		SetCompoundArg(w,1,il);
		
		ret=MakeCompound(A_MTERM,4);
		SetCompoundArg(ret,1,NewInteger(1));
		SetCompoundArg(ret,2,NewInteger(1));
		SetCompoundArg(ret,3,AppendFirst(NewList(),w));
		return MakeCompound2(A_ALG1,AppendFirst(NewList(),ret),
				CopyTerm(il));
		}
		
	if(ino!=ListLength(il))
		{
		List l;
		Atomic ind;
		
		ilc=ListNthList(il,ino);
		ind=CompoundArg2(ListFirst(ilc));
		il=CutFromList(il,ilc);
		l=al;
		while(!is_empty_list(l))
			{
			sub_ind_alg(ListFirst(l),ind,i);
			l=ListTail(l);
			}
		
		SetCompoundArg(w,1,il);
		SetCompoundArg(w,2,al);
		ret=MakeCompound(A_MTERM,4);
		SetCompoundArg(ret,1,NewInteger(1));
		SetCompoundArg(ret,2,NewInteger(1));
		SetCompoundArg(ret,3,AppendFirst(NewList(),w));
		return MakeCompound2(A_ALG1,AppendFirst(NewList(),ret),
				CopyTerm(il));
		}
	FreeAtomic2(w,mind);
	FreeAtomic2(il,mind);
	il=ListNthList(al,i);
	ret=ListFirst(il);
	ChangeList(il,0);
	FreeAtomic2(al,mind);
	return ret;
	}

	
void sub_ind_alg(Term a1, Atomic ind, int i)
	{
	List l,l1,l2,ret;
	/*WriteTerm(a1); printf("  "); WriteTerm(ind); printf("=%d\n",i);*/
	l=ConsumeCompoundArg(a1,2);
	l1=l;
	while(!is_empty_list(l1))
		{
		if(CompoundArg2(ListFirst(l1))==ind)
			{
			l=CutFromList(l,l1);
			SetCompoundArg(a1,2,l);
			goto p2;
			}
		l1=ListTail(l1);
		}
	puts("Internal error (iail failure)");
	longjmp(alg1_jmp_buf,1);
p2:
	ret=NewList();
	l=CompoundArg1(a1);
	while(!is_empty_list(l))
		{
		Term m1;
		m1=ListFirst(l);
		ChangeList(l,0);
		l1=CompoundArgN(m1,3);
		while(!is_empty_list(l1))
			{
			Term sp;
			int ino=1;
			Term wa1;
			sp=ListFirst(l1);
			if(CompoundName(sp)!=OPR_WILD)
				goto c1cnt;
			l2=CompoundArg1(sp);
			while(!is_empty_list(l2))
				{
				if(CompoundArg2(ListFirst(l2))==ind)
					{
					ChangeList(l1,0);
					wa1=sub_ind_wild(sp,ino,i);
					ret=ConcatList(ret,mk_sub1(m1,l1,wa1));
					goto ccnt;
					}
				ino++;
				l2=ListTail(l2);
				}		
		c1cnt:
			l1=ListTail(l1);
			}
		puts("Internal error (iawl failure)");
		longjmp(alg1_jmp_buf,1);
	ccnt:
		l=ListTail(l);
		}
	
	l=ConsumeCompoundArg(a1,1);
	RemoveList(l);		
	SetCompoundArg(a1,1,ret);
	/*WriteTerm(a1);
	puts("  sia finished");	*/			
	return;	
	}
	
	

static List sub_w_2(Term m1, Atomic ind, int ci1, int i1, int ci2, int i2, int i)
	{
	Term w1,w2;
	List wl1,wl2,ret;
	MTHOPTDECL
			
	/*puts("sw2 start");
	WriteTerm(m1);
	printf("\t where "); WriteTerm(ind); printf("=%d\n",i);
	*/
	wl1=ListNthList(CompoundArgN(m1,3),ci1);
	wl2=ListNthList(CompoundArgN(m1,3),ci2);
	w1=ListFirst(wl1);
	w2=ListFirst(wl2);
	ChangeList(wl1,0);
	ChangeList(wl2,0);
	w1=sub_ind_wild(w1,i1,i);
	w2=sub_ind_wild(w2,i2,i);
	/*
	printf("after sub:  ");
	WriteTerm(w1);
	printf("  ,  ");
	WriteTerm(w2);
	puts("\ngoing to mk_sub2");
	*/
	ret = mk_sub2(m1,wl1,w1,wl2,w2);
	FreeAtomic2(w1,mind);
	FreeAtomic2(w2,mind);
	return ret;
	}

	

static List sub_w_1(Term m1, Atomic ind, int c1, int i1, int i2, int i)
	{
	Term w;
	List wl,ret;
/*	puts("sw1 called");*/
	wl=ListNthList(CompoundArgN(m1,3),c1);
	w=ListFirst(wl);
	if(i1>i2)
		{
		int tmp;
		tmp=i1;
		i1=i2;
		i2=tmp;
		}
	w=sub_ind_wild(w,i2,i);
	sub_ind_alg(w,ind,i);
	ChangeList(wl,0);
	
/*	printf("sw1:  ");
	WriteTerm(w);
	puts("");
*/	
	ret=mk_sub1(m1,wl,w);
	return ret;
	}



List alg1_sub_w(Term m1, Atomic ind, int c1, int i1, int c2, int i2, int i)
	{
	MTHOPTDECL
		
	if(c1==c2)
		return sub_w_1(CopyTerm2(m1,mind),ind,c1,i1,i2,i);
	else
		return sub_w_2(CopyTerm2(m1,mind),ind,c1,i1,c2,i2,i);
	}
	
	
void alg1_exp_wild(Term a1, List nind)
{
	List ret, al, l1, l2, l3, l4;
	
	ret=0;
	al=ConsumeCompoundArg(a1,1);
	
	for(l1=nind;l1;l1=ListTail(l1))
	{
		Label il;
		long int   iv, c=0, i=0;
		il=CompoundArg1(ListFirst(l1));
		iv=IntegerValue(CompoundArg2(ListFirst(l1)));
		
		for(l2=al;l2;l2=ListTail(l2))
		{
			Term m1;
			m1=ListFirst(l2);
			for(c=1,l3=CompoundArgN(m1,3);l3;c++,l3=ListTail(l3))
			{
				for(i=1,l4=CompoundArg1(ListFirst(l3));l4;i++,l4=ListTail(l4))
					if(CompoundArg2(ListFirst(l4))==il)
						break;
				if(l4)
					break;
			}
			
			if(is_empty_list(l3))
				ret=AppendFirst(ret,m1);
			else
			{
				Term w,wl;
				wl=ListNthList(CompoundArgN(m1,3),(int)c);
				w=ListFirst(wl);
				w=sub_ind_wild(w,(int)i,(int)iv);
				ChangeList(wl,0);
				ret=ConcatList(mk_sub1(m1,wl,w),ret);
			}
		}
		
		RemoveList(al);
		al=ret;
		ret=0;
	}
	
	SetCompoundArg(a1,1,al);
}

		
