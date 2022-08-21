#include <setjmp.h>
#include <string.h>
#include "lanhep.h"

List alg2_denorm(Term);

		
extern jmp_buf alg2_jmp_buf;


static Integer p_num(Label l, List pl)
	{
	int nu=1;
	while(!is_empty_list(pl))
		{
		if(l==CompoundArg2(ListFirst(pl)))
			{ return NewInteger(nu);}
		nu++;
		pl=ListTail(pl);
		}
	printf("internal error (prt num lookup failed)\n");
	longjmp(alg2_jmp_buf,1);
	return 0;
	}
	
static void lbl_to_i(List l, Label from, Integer to)
	{
	List l1;
	while(!is_empty_list(l))
		{
		l1=CompoundArg2(ListFirst(l));
		while(!is_empty_list(l1))
			{
			if(ListFirst(l1)==from)
				ChangeList(l1,to);
			l1=ListTail(l1);
			}
		l=ListTail(l);
		}
	}	
		 	

static int sp_cmp(Term t1, Term t2)
	{
	Term f1,f2;
	int i1,i2;
	f1=CompoundName(t1);
	f2=CompoundName(t2);
	if(f1==A_MOMENT && f2==OPR_SPECIAL)
		{return -1;}
	if(f2==A_MOMENT && f1==OPR_SPECIAL)
		{return 1; }
	if(f1==A_MOMENT)
		{
		i1=(int)IntegerValue(CompoundArg1(t1));
		i2=(int)IntegerValue(CompoundArg1(t2));
		if(i1==i2) goto cmpl;
		if(i1<i2)
			return -1;
		else
			return 1;
		}
	
	i1=strcmp(AtomValue(CompoundArg1(t1)),AtomValue(CompoundArg1(t2)));
	if(i1!=0)
		return i1;
cmpl:

	f1=CompoundArg2(t1);
	f2=CompoundArg2(t2);
	
	i1=ListLength(f1)-ListLength(f2);
	if(i1!=0)
		return i1;
	
	while(!is_empty_list(f1))
		{
		Term l1,l2;
		l1=ListFirst(f1);
		l2=ListFirst(f2);
		if(is_integer(l1) && is_label(l2))
			return -1;
		if(is_integer(l2) && is_label(l1))
			return 1;
		if(is_integer(l1))
			{
			l1=IntegerValue(l1);
			l2=IntegerValue(l2);
			if(l1<l2)
				return -1;
			if(l2<l1)
				return 1;
			}
		else
			{
			if(l1<l2)
				return -1;
			if(l2<l1)
				return 1;
			}
		f1=ListTail(f1);
		f2=ListTail(f2);
		}
	return 0;
	}
	
static int inv_moment(List l)
	{
	int ch=1;
	while(!is_empty_list(l))
		{
		Term t;
		t=ListFirst(l);
		if(CompoundName(t)==A_MOMENT && IntegerValue(CompoundArg1(t))==2)
			{
			ch*=-1;
			SetCompoundArg(t,1,NewInteger(1));
			}
		l=ListTail(l);
		}
	return ch;
	}
/*
static void alg2_fix_eqv1(List pl, Term m2, int pos, int len)
{
	int lcarr[32];
	List l1, l2;
	int find;
	int i;
	if(len<2 || len>32)
		puts("Internal error (a2afev11)");
	
//	WriteVertex(pl);
//	printf(" pos=%d len=%d ",pos,len);
	
	find=IntegerValue(ListFirst(CompoundArg1(CompoundArg2(ListNth(pl,pos)))));
//	printf("find=%d carr=( ",find);
	
	for(i=find;i<find+len;i++)
	{
		int gno=0;
		for(l1=CompoundArgN(m2,3);l1;l1=ListTail(l1))
		{
			Term t;
			t=ListFirst(l1);
			if(CompoundName(t)==A_MOMENT && ListFirst(CompoundArg2(t))==NewInteger(i))
			{
				lcarr[i-find]=IntegerValue(CompoundArg1(t))+100;
				break;
			}
			
			if(CompoundArg1(t)==A_DELTA && 
					ListFirst(CompoundArg2(t))==NewInteger(i))
			{
				lcarr[i-find]=IntegerValue(ListFirst(ListTail(CompoundArg2(t))));
				break;
			}
			
			if(CompoundArg1(t)==A_DELTA && 
					ListFirst(ListTail(CompoundArg2(t)))==NewInteger(i))
			{
				lcarr[i-find]=IntegerValue(ListFirst(CompoundArg2(t)));
				break;
			}
			
			if(CompoundArg1(t)==A_GAMMA)
			{
				gno++;
				if(ListFirst(ListTail(ListTail(CompoundArg2(t))))==NewInteger(i))
				{
					lcarr[i-find]=gno+200;
					break;
				}
			}
		}
		
		if(is_empty_list(l1))
			puts("Internal error (a2afev12)");
	}
			
//	for(i=0;i<len;i++)
//		printf("%d ",lcarr[i]);
//	puts(")");
}

	
static void alg2_norm_eq1(List pl, Term m2)
{
	List l1;
	int i=0;
	
	for(l1=pl;l1&&ListTail(l1);l1=ListTail(l1))
	{
		Term t1,t2;
		int j;
		i++;
		t1=ListFirst(l1);
		t2=ListFirst(ListTail(l1));
		if(CompoundArg1(t1)!=CompoundArg1(t2))
			continue;
		if(CompoundName(CompoundArg2(t1))!=OPR_VECTOR)
			continue;
		if(ListLength(CompoundArg1(CompoundArg2(t1)))!=1)
			continue;
		j=1;
		while(ListTail(l1) && CompoundArg1(ListFirst(l1))==
						CompoundArg1(ListFirst(ListTail(l1))))
		{
			l1=ListTail(l1);
			j++;
		}
		
		alg2_fix_eqv1(pl, m2, i, j);
		i+=j-1;
	}
	
}
*/
	
void alg2_norm_delta(List pl, Term m2)
{
	List l1,l2,sl;
	Integer v1=0, v2=0, c1=0, c2=0;
	int cp1=0, cp2=0, cl1=0, cl2=0, pno;
	
	if(is_empty_list(pl))
		return;
/*	
	if(ListLength(pl)>2)
		alg2_norm_eq1(pl,m2);
*/
	sl=CompoundArgN(m2,3);
	
	for(l1=sl;l1;l1=ListTail(l1))
	{
		if(CompoundArg1(ListFirst(l1))==A_DELTA)
		{
			l2=CompoundArg2(ListFirst(l1));
			v1=ListFirst(l2);
			v2=ListFirst(ListTail(l2));
			if(IntegerValue(v1)>IntegerValue(v2))
			{
				ChangeList(l2,v2);
				ChangeList(ListTail(l2),v1);
			}
		}
	}
	
	if(ListLength(pl)<3)
		return;
		
	for(l1=pl,pno=1;ListTail(l1);l1=ListTail(l1),pno++)
	{
		Term t1,t2;
		t1=ListFirst(l1);
		t2=ListFirst(ListTail(l1));
		if(CompoundArg1(t1)!=CompoundArg1(t2))
			continue;
		if(CompoundName(CompoundArg2(t1))!=OPR_VECTOR)
			continue;
		if(ListLength(CompoundArg1(CompoundArg2(t1)))!=2)
			continue;
		if(ListTail(ListTail(l1)) && 
			CompoundArg1(ListFirst(ListTail(ListTail(l1))))==CompoundArg1(t1))
			return;
		v1=ListFirst(CompoundArg1(CompoundArg2(t1)));
		v2=ListFirst(CompoundArg1(CompoundArg2(t2)));
		c1=ListFirst(ListTail(CompoundArg1(CompoundArg2(t1))));
		c2=ListFirst(ListTail(CompoundArg1(CompoundArg2(t2))));
		break;
	}
	
	if(is_empty_list(ListTail(l1)))
		return;
			
	for(l1=sl;l1;l1=ListTail(l1))
		if(CompoundArg1(ListFirst(l1))==A_DELTA && 
			ListFirst(CompoundArg2(ListFirst(l1)))==v1 &&
			ListFirst(ListTail(CompoundArg2(ListFirst(l1))))==v2)
				break;
	
	if(is_empty_list(l1))
		return;
	
	for(l1=sl;l1;l1=ListTail(l1))
	{
		cl1++;
		if((cp1=ListMember(CompoundArg2(ListFirst(l1)),c1)))
		{
			l1=CompoundArg2(ListFirst(l1));
			break;
		}
	}
	
	for(l2=sl;l2;l2=ListTail(l2))
	{
		cl2++;
		if((cp2=ListMember(CompoundArg2(ListFirst(l2)),c2)))
		{
			l2=CompoundArg2(ListFirst(l2));
			break;
		}
	}
	
	if(is_empty_list(l1) || is_empty_list(l2))
	{
		puts("Internal error (a2nrdl)");
		return;
	}
	
	if(cl2>cl1 || (cl2==cl1 && cp2>cp1))
		return;
	
	ChangeList(ListNthList(l1,cp1),c2);
	ChangeList(ListNthList(l2,cp2),c1);
	
	for(l1=sl;l1;l1=ListTail(l1))
	{
		if(CompoundName(ListFirst(l1))==A_MOMENT)
		{
			if(CompoundArg1(ListFirst(l1))==NewInteger(pno))
				SetCompoundArg(ListFirst(l1),1,NewInteger(pno+1));
			else if(CompoundArg1(ListFirst(l1))==NewInteger(pno+1))
				SetCompoundArg(ListFirst(l1),1,NewInteger(pno));
		}
	}
}	
	
	
void alg2_norm(Term a2)
	{
	Term m2;
	List l,l1,l2,pl;
	int i;
	
/*	puts("\nvd in");
	WriteTerm(a2);
	puts("");
*/
		
	m2=ListFirst(CompoundArgN(a2,5));
	pl=CompoundArg1(a2);
	l1=ConsumeCompoundArg(m2,3);
	l=l1;

			/*   particle label -> number */
	
	while(!is_empty_list(l))
		{
		Term t;
		t=ListFirst(l);
		if(is_label(CompoundArg1(t)))
			SetCompoundArg(t,1, p_num(CompoundArg1(t),pl) );
		l=ListTail(l);
		}
		
	/*printf("vd1: "); WriteTerm(a2); puts("");	*/
	
	
		/*  p1 -> -p2 if two prtcls */
	
	if(ListLength(pl)==2 && !FAOutput)
		{
		long int num;
		num=IntegerValue(CompoundArg1(CompoundArg1(m2)));
		num*=inv_moment(l1);
		SetCompoundArg(CompoundArg1(m2),1,NewInteger(num));
		}
	
			/*	Rem prtcls from spec and move indices to plist */
	
	l=l1;
	while(!is_empty_list(l))
		{
		Term t,t1,t2;
		t=ListFirst(l);
		t1=CompoundName(t);
		if(t1==OPR_SCALAR || t1==OPR_SPINOR || t1==OPR_VECTOR || 
					t1==OPR_SPINOR3 || t1==OPR_TENSOR)
			{
			int no;
			no=(int)IntegerValue(CompoundArg1(t));
			t2=ListNth(pl,no);
			t1=MakeCompound1(t1,ConsumeCompoundArg(t,2));
			SetCompoundArg(t2,2,t1);
			l1=CutFromList(l1,l);
			l=l1;
			continue;
			}
		l=ListTail(l);
		}
		
		/*	Convert convols btw prtcls to deltas in specs, gen new ind */
	
	l2=NewList();
	l=pl;
	while(!is_empty_list(l))
		{
		List li;
		li=CompoundArg1(CompoundArg2(ListFirst(l)));
		while(!is_empty_list(li))
			{
			Label ll;
			ll=ListFirst(li);
			if(ListMember(l2,ll))
				{
				Label ll1;
				ll1=NewLabel();
				ChangeList(li,ll1);
				l1=AppendFirst(l1,MakeCompound2(OPR_SPECIAL,A_DELTA,
						AppendLast(AppendLast(NewList(),ll),ll1)));
				}
			else
				l2=AppendLast(l2,ll);
			li=ListTail(li);
			}
		l=ListTail(l);
		}
	RemoveList(l2);
	
	
		/*	Convert prtcls index labels to ordinal numders */
		
	l=pl;
	i=1;
	while(!is_empty_list(l))
		{
		List li;
		li=CompoundArg1(CompoundArg2(ListFirst(l)));
		while(!is_empty_list(li))
			{
			Label ll;
			ll=ListFirst(li);
			lbl_to_i(l1,ll,NewInteger(i));
			ChangeList(li,NewInteger(i));
			i++;
			li=ListTail(li);
			}
		l=ListTail(l);
		}
			
	/*printf("%d indices in prtcls\n",i-1);*/
	
	
		/*  Sort specs with labels 		*/
	
	{
	List gl=0;
rp:	for(l=l1;l;l=ListTail(l))
		{
		Term t=CompoundArg1(ListFirst(l));
		if(t==A_GAMMA || t==A_GAMMAP || t==A_GAMMAM || t==A_GAMMA5)
			{
			gl=AppendLast(gl,ListFirst(l));
			ChangeList(l,0);
			l1=CutFromList(l1,l);
			goto rp;
			}
		}
	l1=SortedList(l1,sp_cmp);
	l1=ConcatList(l1,gl);
	}
	
		/*  Convert to numbers trailing lab indices */
	
	l=l1;
	while(!is_empty_list(l))
		{
		List li;
		li=CompoundArg2(ListFirst(l));
		while(!is_empty_list(li))
			{
			Label ll;
			ll=ListFirst(li);
			if(is_label(ll))
				{
				lbl_to_i(l1,ll,NewInteger(i));
				i++;
				}
			li=ListTail(li);
			}
		l=ListTail(l);
		}
	
/*	WriteTerm(m2); puts("   ::"); DumpList(l1); puts(""); */
		
	/*	Here the place to apply symm rules to specials */
	/*  With subseq sorting   */
	
	
	SetCompoundArg(m2,3,l1);

	color_symm_f(pl,m2);
	epsv_symm(pl,m2);
	
	alg2_norm_delta(pl,m2);

	/*printf("vd2: "); WriteTerm(a2); puts("");*/
	
	
	}
	

List alg2_denorm(Term a2)
{
	List p1,q1,r1,p5,ret;
	List l1,l2,l3;
	int i,j;
	
	Label plabs[10], ilabs[100];
	
	p1=ConsumeCompoundArg(a2,1);
	p5=ConsumeCompoundArg(a2,5);
	FreeAtomic(a2);
	
	j=ListLength(p1);
	for(i=0;i<j;i++)
		plabs[i]=NewLabel();
	
	
	j=0;
	for(l1=p5;l1;l1=ListTail(l1))
		for(l2=CompoundArgN(ListFirst(l1),3);l2;l2=ListTail(l2))
			for(l3=CompoundArg2(ListFirst(l2));l3;l3=ListTail(l3))
			{
				int k;
				k=(int)IntegerValue(ListFirst(l3));
				if(k>j)
					j=k;
			}
	for(i=0;i<j;i++)
		ilabs[i]=NewLabel();
	
	i=0;
	q1=NewList();
	r1=NewList();
	for(l1=p1;l1;l1=ListTail(l1))
	{
		Term t;
		i++;
		q1=AppendLast(q1,MakeCompound2(OPR_DIV,CompoundArg1(ListFirst(l1)),
				plabs[i-1]));
		
		t=MakeCompound(CompoundName(CompoundArg2(ListFirst(l1))),2);
		SetCompoundArg(t,1,plabs[i-1]);
		SetCompoundArg(t,2,ConsumeCompoundArg(CompoundArg2(ListFirst(l1)),1));
		for(l2=CompoundArg2(t);l2;l2=ListTail(l2))
			ChangeList(l2,ilabs[IntegerValue(ListFirst(l2))-1]);
		r1=AppendLast(r1,t);
	}
	
	FreeAtomic(p1);
	
	for(l1=p5;l1;l1=ListTail(l1))
	for(l2=CompoundArgN(ListFirst(l1),3);l2;l2=ListTail(l2))
	{
		if(CompoundName(ListFirst(l2))==A_MOMENT)
			SetCompoundArg(ListFirst(l2),1,
					plabs[IntegerValue(CompoundArg1(ListFirst(l2)))-1]);
		for(l3=CompoundArg2(ListFirst(l2));l3;l3=ListTail(l3))
			ChangeList(l3,ilabs[IntegerValue(ListFirst(l3))-1]);
	}
	
	ret=NewList();
	
	for(l1=p5;l1;l1=ListTail(l1))
	{
		Term m2;
		List l4;
		
		m2=ListFirst(l1);
		ChangeList(l1,0);
		
		l2=ConsumeCompoundArg(m2,3);
		
		l2=ConcatList(CopyTerm(r1),l2);
	
	rep:	
		for(l3=l2;l3;l3=ListTail(l3))
			if(CompoundArg1(ListFirst(l3))==A_DELTA)
			{
				Label la1,la2;
				la1=ListFirst(CompoundArg2(ListFirst(l3)));
				la2=ListFirst(ListTail(CompoundArg2(ListFirst(l3))));
				l2=CutFromList(l2,l3);
				for(l3=l2;l3;l3=ListTail(l3))
					for(l4=CompoundArg2(ListFirst(l3));l4;l4=ListTail(l4))
						if(ListFirst(l4)==la1)
							ChangeList(l4,la2);
				goto rep;
			}
			
		SetCompoundArg(m2,3,l2);
		
		a2=MakeCompound(A_ALG2,5);
		SetCompoundArg(a2,1,CopyTerm(q1));
		SetCompoundArg(a2,5,MakeList1(m2));
		ret=AppendLast(ret,a2);
	}
	
	FreeAtomic(q1);
	FreeAtomic(r1);
	
	
	return ret;
}
