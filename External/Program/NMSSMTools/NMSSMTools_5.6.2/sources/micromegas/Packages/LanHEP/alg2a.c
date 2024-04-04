#include <setjmp.h>
#include <string.h>
#include "lanhep.h"


extern jmp_buf alg2_jmp_buf;
extern int NoColors;
int opOnlyMass=0;

static int ghost_factor(List *prt)
	{
	List gprt,l1,l2;
	int ret;
	int flag;
	ret=1;

	gprt=*prt;

	do
		{
		flag=0;
		l1=gprt;
		l2=ListTail(l1);
		while(!is_empty_list(l2))
			{
			Atom a1, a2;
			int spr;
			char *n1, *n2;
			a1=ListFirst(l1);
			a2=ListFirst(l2);
			n1=AtomValue(CompoundArg1(a1)); while(n1[1]) n1++;
			n2=AtomValue(CompoundArg1(a2)); while(n2[1]) n2++;
			if(n1[0]=='c' && n2[0]=='C')
				spr=1;
			else if(n2[0]=='c' && n1[0]=='C')
				spr=-1;
			else
				spr=strcmp(AtomValue(CompoundArg1(a1)),
						AtomValue(CompoundArg1(a2)));
			if(spr==0)
				return 0;
			if(spr>0)
				{
				ChangeList(l1,a2);
				ChangeList(l2,a1);
				flag=1;
				ret*=-1;
				}
			l1=l2;
			l2=ListTail(l2);
			}
		}
		while(flag);

	return ret;
	}


static List add_param(List l, Atom param, int p)
	{
	List li;
	li=l;
	while(!is_empty_list(li))
		{
		Term t;
		t=ListFirst(li);
		if(CompoundArg1(t)==param)
			{
			int pw;
			pw=(int)IntegerValue(CompoundArg2(t));
			pw+=p;
			if(pw==0)
				return CutFromList(l,li);
			SetCompoundArg(t,2,NewInteger(pw));
			return l;
			}
		li=ListTail(li);
		}
	return AppendLast(l,MakeCompound2(OPR_POW,param,NewInteger(p)));
	}


static int prmcmp(Term p1, Term p2)
	{
		if(CompoundArg1(p1)==A_I)
			return -1;
		if(CompoundArg1(p2)==A_I)
			return 1;
		return strcmp(AtomValue(CompoundArg1(p1)),
				AtomValue(CompoundArg1(p2)));
	}

extern void prt2cls(Atom *);
		
static int prtcmp(Term p1, Term p2)
	{
		Atom a1, a2;
		if(FAOutput)
		{
			Term prp;
			int s1=0,s2=0;
			p1=CompoundArg1(p1);
			p2=CompoundArg1(p2);
			prp=GetAtomProperty(p1,PROP_INDEX);
			if(prp && CompoundName(CompoundArg1(ListFirst(prp)))==A_LORENTZ)
				s1=2;
			prp=GetAtomProperty(p2,PROP_INDEX);
			if(prp && CompoundName(CompoundArg1(ListFirst(prp)))==A_LORENTZ)
				s2=2;
			if(s2>s1) return -1;
			if(s1>s2) return 1;
			a1=p1;
			a2=p2;
			prt2cls(&p1);
			prt2cls(&p2);
			if(p1!=p2) return strcmp(AtomValue(p1),AtomValue(p2));
			return strcmp(AtomValue(a1),AtomValue(a2));
		}
		a1=CompoundArg1(p1);
		a2=CompoundArg1(p2);
		prt2cls(&a1);
		prt2cls(&a2);
		if(a1!=a2) return strcmp(AtomValue(a1),AtomValue(a2));
		a1=CompoundArg1(p1);
		a2=CompoundArg1(p2);
		return strcmp(AtomValue(a1),AtomValue(a2));
	}

static int proc_moment(List l, List prtcls)
	{
	List pl;
	int retval=1;
/*	WriteTerm(l); puts(":pm");
	WriteTerm(prtcls); puts(":part");
*/	pl=l;
	while(!is_empty_list(l))
		{
		Term t;
		t=ListFirst(l);
		if(CompoundName(t)==OPR_SPECIAL && CompoundArg1(t)==A_MOMENT)
			{
			Term t1,t2;
			List l1;
			l1=ListTail(l);
			if(is_empty_list(l1))
				{
				printf("Error in vertex "); WriteVertex(prtcls); printf(":");
				printf(" moment must be followed with a field\n");
				return 0;
				}
			while(!is_empty_list(l1))
				{
				t1=ListFirst(l1);
				t2=CompoundName(t1);
				if(t2==OPR_SCALAR || t2==OPR_VECTOR || t2==OPR_SPINOR
						|| t2==OPR_SPINOR3 || t2==OPR_TENSOR)
					break;
				l1=ListTail(l1);
				}
			if(is_empty_list(l1))
				{
				printf("Error in vertex "); WriteVertex(prtcls); printf(":");
				printf("moment must be followed with a field\n");
				/*longjmp(alg2_jmp_buf,1);*/
				return 0;
				}
			t2=MakeCompound2(A_MOMENT,CompoundArg1(t1),
					ConsumeCompoundArg(t,2));
			if(CompoundArg1(t2)==A_VEV)
				retval=0;
			/*retval*=ch_s_moment(CompoundArg1(t2),prtcls);*/
			ChangeList(l,t2);
			}
		l=ListTail(l);
		}
	return retval;
	}

	
/*extern int verbose_mtermto2;*/
	
Term alg2_mterm_to_2(Term t)
	{
	List prt,params,spec,l,ls;
	long int num,den;
	Term a2, m2;

/*WriteTerm(t);puts("");*/

/*	
if(verbose_mtermto2){
DumpList(CompoundArgN(t,3));puts("m22 num");
DumpList(CompoundArgN(t,4));puts("m22 den");
}
*/	
	num=IntegerValue(CompoundArg1(t));
	den=IntegerValue(CompoundArg2(t));
	params=spec=prt=NewList();

	ls=l=ConsumeCompoundArg(t,3);


		/* field -> spin; divide spec & params with powers */

	while(!is_empty_list(l))
		{
		Term t1;
		t1=ListFirst(l);
		if(CompoundName(t1)==OPR_FIELD || CompoundName(t1)==OPR_SPECIAL)
			{
			Term l1,l2,l3,nm;
			l1=ConsumeCompoundArg(t1,1);
			nm=ConsumeCompoundArg(t1,2);
			if(CompoundName(t1)==OPR_FIELD)
				{
				Term lb,t3;
				lb=NewLabel();
				if(nm==A_VEV)
					lb=A_VEV;
				else
					prt=AppendLast(prt,MakeCompound2(OPR_DIV,nm,lb));
				nm=lb;
				FreeAtomic(t1);
				t1=OPR_SCALAR;
				if(l1!=0)
					{
					Term t2;
					t2=CompoundArg1(ListFirst(l1));
					if(CompoundName(t2)==A_LORENTZ)
						{
						if(CompoundArg1(t2)==NewInteger(2))
							t1=OPR_VECTOR;
						else
							t1=OPR_SPINOR;
						}
					t3=ListTail(l1);
					if(!is_empty_list(t3))
						{
						t3=CompoundArg1(ListFirst(t3));
						if(CompoundName(t3)==A_LORENTZ)
							{
							if(CompoundArg1(t2)!=NewInteger(2))
								t1=OPR_SPINOR3;
							else
								t1=OPR_TENSOR;
							}
						}

					}
				t1=MakeCompound(t1,2);
				}

			l2=l1;
			l3=NewList();
			while(!is_empty_list(l2))
				{
				l3=AppendLast(l3,CompoundArg2(ListFirst(l2)));
				l2=ListTail(l2);
				}
			FreeAtomic(l1);
			SetCompoundArg(t1,2,l3);
			SetCompoundArg(t1,1,nm);
			spec=AppendLast(spec,t1);
			goto ewhile;
			}
		if(CompoundName(t1)==OPR_PARAMETER)
			{
			params=add_param(params,CompoundArg2(t1),1);
			FreeAtomic(t1);
			goto ewhile;
			}

        printf("Error: fail to process %s in 1->2 (",
            AtomValue(CompoundName(t1)));
        WriteTerm(t1); puts("");
		FreeAtomic(t1);
		longjmp(alg2_jmp_buf,1);
	ewhile:
		l=ListTail(l);
		}
	RemoveList(ls);



			/*	move dens to params  */

	ls=l=ConsumeCompoundArg(t,4);
	while(!is_empty_list(l))
		{
		Term t1;
		t1=ListFirst(l);
		if(CompoundName(t1)==OPR_PARAMETER)
			{
			params=add_param(params,CompoundArg2(t1),-1);
			FreeAtomic(t1);
			}
		else
			{
			printf("Error: fail to process %s in 1->2 (den)\n",
									AtomValue(CompoundName(t1)));
			FreeAtomic(t1);
			longjmp(alg2_jmp_buf,1);
			}
		l=ListTail(l);
		}
	RemoveList(ls);
/*	if(is_empty_list(prt))
		{
		ErrorInfo();
		printf("constant term in lagrangian\n");
		return 0;
		}
*/
		/* sorting params list */


	params=SortedList(params,prmcmp);


	if(den<0)
		{
		den=-den;
		params=AppendFirst(params,MakeCompound2(OPR_POW,A_I,NewInteger(1)));
		}

		/* sorting particles */


	{
	List fr,pr,gr;
	fr=pr=gr=NewList();

    num*=alg2_updown(prt,&spec);

	ls=NewList();
aaa45:
	l=prt;
	while(!is_empty_list(l))
		{
		Atom p1;
		Term prop;
		p1=CompoundArg1(ListFirst(l));
		prop=GetAtomProperty(p1,PROP_INDEX);
		if(prop)
			{
			prop=CompoundArg1(ListFirst(prop));
			if(is_compound(prop) && CompoundName(prop)==A_LORENTZ &&
				CompoundArg1(prop)==NewInteger(1) &&
					CompoundArg2(prop)==NewInteger(0) )
				{
				fr=AppendLast(fr,ListFirst(l));
				ChangeList(l,0);
				pr=AppendLast(pr,NewInteger(1));
				prt=CutFromList(prt,l);
				goto aaa45;
				}
			if(is_compound(prop) && CompoundName(prop)==A_LORENTZ &&
				CompoundArg1(prop)==NewInteger(0) &&
					CompoundArg2(prop)==NewInteger(1) )
				{
				fr=AppendLast(fr,ListFirst(l));
				ChangeList(l,0);
				pr=AppendLast(pr,NewInteger(0));
				prt=CutFromList(prt,l);
				goto aaa45;
				}
			}
		l=ListTail(l);
		}

aaa46:
	l=prt;
	while(!is_empty_list(l))
		{
		Atom p1;
		Term prop;
		p1=CompoundArg1(ListFirst(l));
		prop=GetAtomProperty(p1,A_GRASS);
		if(prop)
			{
			gr=AppendLast(gr,ListFirst(l));
			ChangeList(l,0);
			prt=CutFromList(prt,l);
			goto aaa46;
			}
		l=ListTail(l);
		}

	prt=SortedList(prt,prtcmp);
	
	if(!is_empty_list(gr))
		{
		num*=ghost_factor(&gr);
		if(num==0)
			{
			return 0;
			}
		/*if(ghost_factor(&gr)==0)
			return 0;*/
		prt=ConcatList(gr,prt);
		}

	if(!is_empty_list(fr))
		{
		int cf;
		/*WriteTerm(fr);puts("");DumpList(spec);printf(" %d->",num);*/
		num*=alg2_refine_spinor(fr,&spec);
		/*printf("%d\n",num);DumpList(spec);puts("");*/
		if(num==0)
			{
			return 0;
			}
		cf=(int)gcf(num,den);
		num/=cf;
		den/=cf;
		prt=ConcatList(fr,prt);
		}

/*	if(!is_empty_list(fr))
		{
		num*=spin_factor(&fr,pr);
		if(num==0)
			{
			return 0;
			}
		prt=ConcatList(fr,prt);
		}*/



	}


	ls=NewList();

	FreeAtomic(t);


		/* spec(mom) -> mom(prt) */

	num*=proc_moment(spec,prt);
/*	printf("%d :",num);
	WriteTerm(spec);
	puts("");
*/
lab12:
	l=spec;
	while(!is_empty_list(l))
		{
		if(CompoundArg1(ListFirst(l))==A_VEV)
			{
			spec=CutFromList(spec,l);
			goto lab12;
			}
		l=ListTail(l);
		}


	if(num==0) return 0;
	if(opOnlyMass && ListLength(prt)>2)
		return 0;
		

	a2=MakeCompound(A_ALG2,5);
	m2=MakeCompound(A_MTERM,3);

	SetCompoundArg(m2,1,MakeCompound2(OPR_DIV,NewInteger(num),
				NewInteger(den)));
	SetCompoundArg(a2,1,prt);
	SetCompoundArg(m2,2,params);
	SetCompoundArg(m2,3,spec);
	SetCompoundArg(a2,5,AppendFirst(NewList(),m2));

/*
	WriteTerm(prt); puts("");
	DumpList(CompoundArgN(m2,3)); puts("");
*/

		/* symm by eq prtcls */ /*  Apply color reduction! */

/*
    if(ListLength(CompoundArg1(a2))>4)
        return 0;
*/
    
	return a2;
/*
	if((!TexOutput && !NoColors) || ForsedRedCol)
        {
		a2=color_reduce(a2);
        if(a2==0)
            return 0;
        }
	else
		a2=AppendLast(NewList(),a2);


	if(!TexOutput || ForsedRedCol)
		{
		ls=NewList();
		l=a2;
		while(!is_empty_list(l))
			{
			ls=ConcatList(ls,spinor_reduce(ListFirst(l)));
			l=ListTail(l);
			}
		RemoveList(a2);
		a2=ls;
		}


	ls=NewList();
	l=a2;
	while(!is_empty_list(l))
		{
		ls=ConcatList(ls,reduce_56(ListFirst(l)));
		l=ListTail(l);
		}
	RemoveList(a2);
	a2=ls;


	ls=NewList();

	l=a2;
	while(!is_empty_list(l))
		{
		ls=ConcatList(ls,symmetrize(ListFirst(l)));
		l=ListTail(l);
		}
	RemoveList(a2);
*/
/*
	puts("after symmetryze:");
	DumpList(ls);
*/
	return ls;

	}

	
int sp_cmp(Term t1, Term t2)
	{
	Term f1,f2;
	long int i1,i2;
	f1=CompoundName(t1);
	f2=CompoundName(t2);
	if(f1==A_MOMENT && f2==OPR_SPECIAL)
		{return -1;}
	if(f2==A_MOMENT && f1==OPR_SPECIAL)
		{return 1; }
	if(f1==A_MOMENT)
		{
		i1=IntegerValue(CompoundArg1(t1));
		i2=IntegerValue(CompoundArg1(t2));
		if(i1==i2) goto cmpl;
		if(i1<i2)
			return -1;
		else
			return 1;
		}
	i1=strcmp(AtomValue(CompoundArg1(t1)),AtomValue(CompoundArg1(t2)));
	if(i1!=0)
		return (int)i1;
cmpl:

	f1=CompoundArg2(t1);
	f2=CompoundArg2(t2);
	
	i1=ListLength(f1)-ListLength(f2);
	if(i1!=0)
		return (int)i1;
	
	while(!is_empty_list(f1))
		{
		Term l1,l2;
		l1=ListFirst(f1);
		l2=ListFirst(f2);
		
		l1=IntegerValue(l1);
		l2=IntegerValue(l2);
		
		if(l1<l2)
			return -1;
		if(l2<l1)
			return 1;
		
		f1=ListTail(f1);
		f2=ListTail(f2);
		}
	return 0;
	}

	
static int next_perm(int *arr, int size)
{
	int i,j,k,l;

	for(i=size-2;i>=0;i--)
	{
		
		if(arr[i]==size)
			continue;
		
		for(j=arr[i]+1;j<=size;j++)
		{
			int has=0;
			for(k=0;k<i;k++)
				if(arr[k]==j)
					has=1;
			if(has==0)
				break;
		}
		
		if(j>size)
			continue;
		
		arr[i]=j;
		break;
	}

	if(i<0)
		return 0;


	for(j=i+1;j<size;j++)
	{
		
		for(k=1;k<=size;k++)
		{
			int has=0;
			for(l=0;l<j;l++)
				if(arr[l]==k)
					has=1;
			if(has==0)
				break;
		}
		
		arr[j]=k;
	}
	
	return 1;
}
	
static void a2_symm_1(Term a2, int ppos, int pno)
{
	List alo;
	int rpl[10], i, j, pindno;
	List pind[10];
	
	for(i=ppos;i<ppos+pno;i++)
	{
		pind[i-ppos]=CompoundArg1(CompoundArg2(ListNth(CompoundArg1(a2),i)));
	}
	pindno=ListLength(pind[0]);
	
	alo=CopyTerm(CompoundArgN(a2,5));
	
	for(i=1;i<=pno;i++)
		rpl[i-1]=i;
	
	while(next_perm((int *)rpl,pno))
	{
		List ll;
		List rm=0, ri=0;
		List l1,l2,l3,l4,l5;
		
		ll=NewList();
		for(i=0;i<pno;i++)
			ll=AppendLast(ll,NewInteger(rpl[i]));
		
/*		WriteTerm(ll);puts("");*/
		
		for(i=1;i<=pno;i++)
			rm=AppendLast(rm,MakeCompound2(OPR_RARROW,NewInteger(ppos-1+i),
					NewInteger(IntegerValue(ListNth(ll,i))-1+ppos)));
		
		for(i=1;i<=pno;i++)
			for(j=0;j<pindno;j++)
				ri=AppendLast(ri,MakeCompound2(OPR_RARROW,
						ListNth(pind[i-1],j+1),
						ListNth(pind[IntegerValue(ListNth(ll,i))-1],j+1)));
		
		l1=CopyTerm(alo);
		
		for(l2=l1;l2;l2=ListTail(l2))
		for(l3=CompoundArgN(ListFirst(l2),3);l3;l3=ListTail(l3))
			if(CompoundName(ListFirst(l3))==A_MOMENT)
			{
				for(l4=rm;l4;l4=ListTail(l4))
					if(CompoundArg1(ListFirst(l4))==CompoundArg1(ListFirst(l3)))
					{
						SetCompoundArg(ListFirst(l3),1,CompoundArg2(ListFirst(l4)));
						break;
					}
			}
			
		for(l2=l1;l2;l2=ListTail(l2))
		for(l3=CompoundArgN(ListFirst(l2),3);l3;l3=ListTail(l3))
		for(l4=CompoundArg2(ListFirst(l3));l4;l4=ListTail(l4))
			{
				for(l5=ri;l5;l5=ListTail(l5))
					if(ListFirst(l4)==CompoundArg1(ListFirst(l5)))
					{
						ChangeList(l4,CompoundArg2(ListFirst(l5)));
						break;
					}
			}
			
		
		for(l2=l1;l2;l2=ListTail(l2))
			for(l3=CompoundArgN(ListFirst(l2),3);l3;l3=ListTail(l3))
				if(CompoundArg1(ListFirst(l3))==A_DELTA)
				{
					long int la1, la2;
					la1=IntegerValue(ListFirst(CompoundArg2(ListFirst(l3))));
					la2=IntegerValue(ListFirst(ListTail(
							CompoundArg2(ListFirst(l3)))));
					if(la1>la2)
					{
						ChangeList(CompoundArg2(ListFirst(l3)),NewInteger(la2));
						ChangeList(ListTail(CompoundArg2(ListFirst(l3))),
								NewInteger(la1));
					}
				}

		if(ListLength(CompoundArg1(a2))==2 && !FAOutput)
			for(l2=l1;l2;l2=ListTail(l2))
			{
				int s=1;
				for(l3=CompoundArgN(ListFirst(l2),3);l3;l3=ListTail(l3))
					if(CompoundName(ListFirst(l3))==A_MOMENT &&
							CompoundArg1(ListFirst(l3))==NewInteger(2))
					{
						SetCompoundArg(ListFirst(l3),1,NewInteger(1));
						s*=-1;
					}
				if(s==-1)
					SetCompoundArg(CompoundArg1(ListFirst(l2)),1,NewInteger(
						-IntegerValue(CompoundArg1(CompoundArg1(ListFirst(l2))))));
			}
			

								
		for(l2=l1;l2;l2=ListTail(l2))
		{
			List l,gl=0;
			l3=ConsumeCompoundArg(ListFirst(l2),3);
			
		rp:	for(l=l3;l;l=ListTail(l))
				{
				Term t=CompoundArg1(ListFirst(l));
				if(t==A_GAMMA || t==A_GAMMAP || t==A_GAMMAM || t==A_GAMMA5)
					{
					gl=AppendLast(gl,ListFirst(l));
					ChangeList(l,0);
					l3=CutFromList(l3,l);
					goto rp;
					}
				}
			
			l3=SortedList(l3,sp_cmp);
			l3=ConcatList(l3,gl);
			SetCompoundArg(ListFirst(l2),3,l3);
			
			epsv_symm(CompoundArg1(a2),ListFirst(l2));
			color_symm_f(CompoundArg1(a2),ListFirst(l2));
		}
		
		
		alg2_add_ml(a2,l1);
/*		
		puts("-------------------------------------------------------");
		WriteTerm(alo);puts("");
		
		WriteTerm(rm);puts("");
		WriteTerm(ri);puts("");
		
		WriteTerm(l1);puts("");
		puts("--------------------------------------------------------");
*/		
		FreeAtomic(ri);
		FreeAtomic(rm);
		FreeAtomic(ll);
	}
	
	FreeAtomic(alo);
/*
	WriteTerm(a2);
	puts("\n-------------------------------------------------------------");
*/
}
	
int opDoSymmetrize=1;

static void fix_im_in(List l3, int inomin)
{
	int inomax=inomin-1;
	List l1,l2;
	List ii=NewList();
	

	if(l3==0)
		return;
		
	for(l1=l3;l1;l1=ListTail(l1))
	{
		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
		{
			 int i;
			i=(int)IntegerValue(ListFirst(l2));
			if(i>=inomin && !ListMember(ii,NewInteger(i)))
			{
				ii=AppendLast(ii,NewInteger(i));
				if(i>inomax)
					inomax=i;
			}
		}
	}
	
	
	if(ListLength(ii)!=inomax-inomin+1)
	{
		puts("Internal error (twff)");
		return;
	}
	
	for(l1=l3;l1;l1=ListTail(l1))
	{
		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
		{
			Integer i;
			int j;
			i=ListFirst(l2);
			j=ListMember(ii,i);
			if(j)
			{
				ChangeList(l2,NewInteger(j+inomin-1));
			}
		}
	}
	
	FreeAtomic(ii);
	return;
}


void alg2_symmetrize(Term a2)
{
	List l, epl=0;
	int ps=0;
	int ino=1;
	
	if(!opDoSymmetrize || CompoundArgN(a2,5)==0)
		return;
	
	for(l=CompoundArg1(a2);l && ListTail(l);l=ListTail(l))
	{
		ps++;
		if(CompoundArg1(ListFirst(l))==CompoundArg1(ListFirst(ListTail(l))))
		{
			int ps1, sz;
			ps1=ps;
			while(ListTail(l) && 
						CompoundArg1(ListFirst(l))==
						CompoundArg1(ListFirst(ListTail(l))))
				l=ListTail(l),ps++;
			sz=ps-ps1;
			epl=AppendLast(epl,MakeCompound2(OPR_DIV,
					NewInteger(ps1),NewInteger(sz+1)));
		}
	}
	
	if(is_empty_list(epl))
		return;
/*	
	WriteTerm(a2);
	puts("");
	WriteTerm(epl);
	puts("");
*/	
	for(l=CompoundArg1(a2);l;l=ListTail(l))
		ino+=ListLength(CompoundArg1(CompoundArg2(ListFirst(l))));
	
	for(l=epl;l;l=ListTail(l))
		a2_symm_1(a2,(int)IntegerValue(CompoundArg1(ListFirst(l))),
					(int)IntegerValue(CompoundArg2(ListFirst(l))));

	for(l=CompoundArgN(a2,5);l;l=ListTail(l))
		fix_im_in(CompoundArgN(ListFirst(l),3),ino);

	FreeAtomic(epl);
}

	
