#include <stdio.h>
#include <string.h>
#include "lanhep.h"

extern int TexOutput;

static int prmcmp(Term p1, Term p2)
	{
	if(CompoundArg1(p1)==A_I)
		return -1;
	if(CompoundArg1(p2)==A_I)
		return 1;
	return strcmp(AtomValue(CompoundArg1(p1)),AtomValue(CompoundArg1(p2)));
	}

int verb_herm=0;
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
/*	return 1;*/
	}
	
	
static void ex_cc(Term ps)
{
	char buf[32];
	int len;
	strcpy(buf,AtomValue(CompoundArg1(ps)));
	if(is_integer(GetAtomProperty(CompoundArg1(ps),OPR_CLASS)))
		return;
	for(len=0;buf[len];len++);
	if(len>2 && buf[len-1]=='c' && buf[len-2]=='.')
		buf[len-2]=0;
	else
	{
		buf[len]='.';
		buf[len+1]='c';
		buf[len+2]=0;
	}
	SetCompoundArg(ps,1,NewAtom(buf,0));
}

static int fix_grass(List *prt)
{
	List l,gr=0;
    int f=1;
aaa46:

	for(l=(*prt);l;l=ListTail(l))
	{
		Atom p1;
		Term prop;
		p1=CompoundArg1(ListFirst(l));
		prop=GetAtomProperty(p1,A_GRASS);
		if(prop)
			{
			gr=AppendLast(gr,ListFirst(l));
			ChangeList(l,0);
			*prt=CutFromList(*prt,l);
			goto aaa46;
			}
	}
	*prt=SortedList(*prt,prtcmp);
	if(!is_empty_list(gr))
		{
		f=ghost_factor(&gr);
		*prt=ConcatList(gr,*prt);
		}
	return f;
}

static List inv_prt(Term a2c)
	{
		List l;
		int no=0,ghf;
		List ord;
		
/*WriteVertex(CompoundArg1(a2c));printf(" --> ");*/
		
		for(l=CompoundArg1(a2c);l;l=ListTail(l))
		{
			List l1;
			no++;
			l1=ConsumeCompoundArg(CompoundArg2(ListFirst(l)),1);
			l1=AppendFirst(l1,NewInteger(no));
			SetCompoundArg(CompoundArg2(ListFirst(l)),1,l1);
		}
			
		for(l=CompoundArg1(a2c);l;l=ListTail(l))
			{
				Atom p,ap;
				p=CompoundArg1(ListFirst(l));
				ap=GetAtomProperty(p,A_ANTI);
				if(ap==0)
				{
					printf("No antiparticle for %s.\n",AtomValue(p));
					return 0;
				}
				SetCompoundArg(ListFirst(l),1,ap);
			}
		l=CompoundArg1(a2c);
		
		if(CompoundName(CompoundArg2(ListFirst(l)))==OPR_SPINOR ||CompoundName(CompoundArg2(ListFirst(l)))==OPR_SPINOR3 )
		{
			Term p1, p2;
			int c1=0,c2=0,do_cc=0;
			Term prop;
			
			p1=ListFirst(l);
			p2=ListFirst(ListTail(l));
			l=ConsumeCompoundArg(a2c,1);
			ChangeList(l,0);
			ChangeList(ListTail(l),0);
			l=CutFromList(l,l);
			l=CutFromList(l,l);
			
			if(GetAtomProperty(CompoundArg1(p1),A_ANTI)==CompoundArg1(p1))
			{
				ex_cc(p1);
				do_cc=1;
			}
			if(GetAtomProperty(CompoundArg1(p2),A_ANTI)==CompoundArg1(p2))
			{
				ex_cc(p2);
				do_cc=1;
			}
            
            prop=GetAtomProperty(CompoundArg1(p1),PROP_TYPE);
			if(CompoundName(prop)==OPR_FIELD && CompoundArg2(prop)==NewInteger(4))
            c1=1;
			prop=GetAtomProperty(CompoundArg1(p2),PROP_TYPE);
			if(CompoundName(prop)==OPR_FIELD && CompoundArg2(prop)==NewInteger(4))
            c2=1;
			if((c1==0 && c2==1) || (c1==1 && c2==0))
            do_cc=1;
            if(c1==0 && c2==0)
                do_cc=0;
			ghf=fix_grass(&l);
			if((c1&&c2) || (do_cc && strcmp(AtomValue(CompoundArg1(p1)),
					AtomValue(CompoundArg1(p2)))<0))
			{
				ex_cc(p1);
				ex_cc(p2);
				l=AppendFirst(l,p2);
				l=AppendFirst(l,p1);
			}
			else
			{
				l=AppendFirst(l,p1);
				l=AppendFirst(l,p2);
			}
            
			SetCompoundArg(a2c,1,l);
		}
		else
		{
            l=ConsumeCompoundArg(a2c,1);
			ghf=fix_grass(&l);
			SetCompoundArg(a2c,1,l);
		}
		
		if(ghf<0)
		{
			List l;
			for(l=CompoundArgN(a2c,5);l;l=ListTail(l))
			{
				int num;
				num=(int)IntegerValue(CompoundArg1(CompoundArg1(ListFirst(l))));
				SetCompoundArg(CompoundArg1(ListFirst(l)),1,NewInteger(-num));
			}
		}
		
		ord=NewList();
		for(l=CompoundArg1(a2c);l;l=ListTail(l))
		{
			List l1;
			l1=ConsumeCompoundArg(CompoundArg2(ListFirst(l)),1);
			ord=AppendLast(ord,ListFirst(l1));
			l1=CutFromList(l1,l1);
			SetCompoundArg(CompoundArg2(ListFirst(l)),1,l1);
		}
		
/*WriteVertex(CompoundArg1(a2c));puts("");*/
		return ord;
	}

static void fix_mom(List ord, List ml)
{
	for(;ml;ml=ListTail(ml))
	{
		List l;
		for(l=CompoundArgN(ListFirst(ml),3);l;l=ListTail(l))
		{
			if(CompoundName(ListFirst(l))==A_MOMENT)
				SetCompoundArg(ListFirst(l),1,NewInteger(ListMember(ord,
						CompoundArg1(ListFirst(l)))));
		}
	}
}

static int frcch=0;

static void fix_mom_2(Term a2)
{
	List l1,l2;
	
	if((ListLength(CompoundArg1(a2))!=2 )/*||!frcch*/ )
		return;
	
	for(l1=CompoundArgN(a2,5);l1;l1=ListTail(l1))
	{
		Term m1;
		int s=1;
		
		m1=ListFirst(l1);

		for(l2=CompoundArgN(m1,3);l2;l2=ListTail(l2))
			if(CompoundName(ListFirst(l2))==A_MOMENT &&
					CompoundArg1(ListFirst(l2))==NewInteger(2))
			{
				SetCompoundArg(ListFirst(l2),1,NewInteger(1));
				s=-s;
			}
			
		if(s==-1)
			SetCompoundArg(CompoundArg1(m1),1,
					NewInteger(-IntegerValue(CompoundArg1(CompoundArg1(m1)))));
	

	}
	
	
}


static void fix_spin(List ord, List prt, List ml)
{
	int cc_done;
	if(CompoundName(CompoundArg2(ListFirst(prt)))!=OPR_SPINOR && CompoundName(CompoundArg2(ListFirst(prt)))!=OPR_SPINOR3)
		return;
	if(ListFirst(ord)==NewInteger(1))
		cc_done=1;
	else
		cc_done=0;
	if(!cc_done)
	{
		Label la1, la2;
		la1=ListFirst(CompoundArg1(CompoundArg2(ListFirst(prt))));
		la2=ListFirst(CompoundArg1(CompoundArg2(ListFirst(ListTail(prt)))));
		ChangeList(CompoundArg1(CompoundArg2(ListFirst(prt))),la2);
		ChangeList(CompoundArg1(CompoundArg2(ListFirst(ListTail(prt)))),la1);
	}
	
	if(!cc_done)
	{
		List l,l1,l2;
		for(l=ml;l;l=ListTail(l))
		{
			Term m2;
			List gl=0, gli=0;
			m2=ListFirst(l);
			for(l1=CompoundArgN(m2,3);l1;l1=ListTail(l1))
			{
				Term t;
				t=ListFirst(l1);
				if(CompoundArg1(t)==A_GAMMA)
				{
					gl=AppendLast(gl,t);
					gli=AppendFirst(gli,ListNth(CompoundArg2(t),3));
				}
			}

			for(l1=gl,l2=gli;l1;l1=ListTail(l1),l2=ListTail(l2))
			{
				Term t;
				t=CompoundArg2(ListFirst(l1));
				t=ListTail(ListTail(t));
				ChangeList(t,ListFirst(l2));
			}

			RemoveList(gl);
			RemoveList(gli);
		}
	}
			
		
	
	for(;ml;ml=ListTail(ml))
	{
		Term m2;
		Term gpm;
		int gno, g5no;
		List l;
		int chs=0;
		gno=0;
		g5no=0;
		gpm=0;
		m2=ListFirst(ml);
		for(l=CompoundArgN(m2,3);l;l=ListTail(l))
		{
			Term t;
			t=ListFirst(l);
			if(CompoundArg1(t)==A_GAMMA)
				gno++;
			if(CompoundArg1(t)==A_GAMMA5)
				g5no++;
			if(CompoundArg1(t)==A_GAMMAM || CompoundArg1(t)==A_GAMMAP)
				gpm=t;
		}
/*	printf("%d %d %d %d\n",cc_done,gno,g5no,gpm);*/
		if(!cc_done)
		{
			if((gno%2)==0 && g5no==1)
				chs=1;
			if((gno%2)==0 && gpm)
			{
				if(CompoundArg1(gpm)==A_GAMMAM)
					SetCompoundArg(gpm,1,A_GAMMAP);
				else
					SetCompoundArg(gpm,1,A_GAMMAM);
			}
		}
		else
		{
			if(((gno%2)==1 && g5no==0) || ((gno%2)==0 && g5no==1))
				chs=1;
			if(gpm)
			{
				if(CompoundArg1(gpm)==A_GAMMAM)
					SetCompoundArg(gpm,1,A_GAMMAP);
				else
					SetCompoundArg(gpm,1,A_GAMMAM);
				/*if((gno%2)==1)
					chs=1;*/
			}
		}
		if(chs)
		{
			int n;
			n=(int)IntegerValue(CompoundArg1(CompoundArg1(m2)));
			SetCompoundArg(CompoundArg1(m2),1,NewInteger(-n));
		}
		
	}
}


static void fix_color_lambda(List ml)
{
	List l;
	for(l=ml;l;l=ListTail(l))
	{
		List ll;
		for(ll=CompoundArgN(ListFirst(l),3);ll;ll=ListTail(ll))
		{
			
			Term sp;
			sp=ListFirst(ll);
			if(CompoundName(sp)==OPR_SPECIAL && 
					GetAtomProperty(CompoundArg1(sp),A_COLOR)==A_COLOR_LAMBDA)
			{
				List l1;
				Label la1,la2;
				l1=CompoundArg2(sp);
				la1=ListFirst(l1);
				la2=ListFirst(ListTail(l1));
				ChangeList(l1,la2);
				ChangeList(ListTail(l1),la1);
			}
		}
	}
}

static void fix_ind_order(List prt, List ml)
{
	List iord=NewList();
	List l;
	for(l=prt;l;l=ListTail(l))
	{
		List l1;
		for(l1=CompoundArg1(CompoundArg2(ListFirst(l)));l1;l1=ListTail(l1))
			iord=AppendLast(iord,ListFirst(l1));
	}
	for(l=ml;l;l=ListTail(l))
	{
		List l1;
		for(l1=CompoundArgN(ListFirst(l),3);l1;l1=ListTail(l1))
		{
			List l2;
			for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
				if(!ListMember(iord,ListFirst(l2)))
					iord=AppendLast(iord,ListFirst(l2));
		}
	}
	
	for(l=prt;l;l=ListTail(l))
	{
		List l1;
		for(l1=CompoundArg1(CompoundArg2(ListFirst(l)));l1;l1=ListTail(l1))
			ChangeList(l1,NewInteger(ListMember(iord,ListFirst(l1))));
	}
	for(l=ml;l;l=ListTail(l))
	{
		List l1;
		for(l1=CompoundArgN(ListFirst(l),3);l1;l1=ListTail(l1))
		{
			List l2;
			for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
				ChangeList(l2,NewInteger(ListMember(iord,ListFirst(l2))));
		}
	}
	
	
}

static void fix_i(List ml)
{
	List l;
	for(l=ml;l;l=ListTail(l))
	{
		List l1;
		Atom an;
		int s=1;
		int cprm=0;
		
		for(l1=CompoundArg2(ListFirst(l));l1;l1=ListTail(l1))
			if((an=GetAtomProperty(CompoundArg1(ListFirst(l1)),A_ANTI)))
			  {SetCompoundArg(ListFirst(l1),1,an);cprm=1;}
		
		if(cprm)
		{
		  List ll=ConsumeCompoundArg(ListFirst(l),2);
		  ll=SortedList(ll,prmcmp);
		  SetCompoundArg(ListFirst(l),2,ll);
		}
		
		for(l1=CompoundArg2(ListFirst(l));l1;l1=ListTail(l1))
			if(CompoundArg1(ListFirst(l1))==A_I)
				break;
		if(l1)
			s=(-s);
		
		/* looking for momenta */
		for(l1=CompoundArgN(ListFirst(l),3);l1;l1=ListTail(l1))
			if(CompoundName(ListFirst(l1))==A_MOMENT)
				s=(-s);
				
		for(l1=CompoundArgN(ListFirst(l),3);l1;l1=ListTail(l1))
			if(CompoundName(ListFirst(l1))==A_COLOR_F)
				s=(-s);
		
		if(s==-1)
		{
			int num;
			num=(int)IntegerValue(CompoundArg1(CompoundArg1(ListFirst(l))));
			SetCompoundArg(CompoundArg1(ListFirst(l)),1,NewInteger(-num));
		}
	}
	
}

static void fix_gh_C(List prt, List ml)
{
	List l;
	int s=1;
	for(l=prt;l;l=ListTail(l))
	{
		Atom p;
		Term pr;
		p=CompoundArg1(ListFirst(l));
		pr=GetAtomProperty(p,PROP_TYPE);
		if(is_compound(pr) && CompoundName(pr)==OPR_FIELD 
				&& CompoundArg2(pr)==NewInteger(3))
			s=-s;
	}
	
	if(s==-s)
	for(l=ml;l;l=ListTail(l))
	{
		int num;
		num=(int)IntegerValue(CompoundArg1(CompoundArg1(ListFirst(l))));
		SetCompoundArg(CompoundArg1(ListFirst(l)),1,NewInteger(-num));
	}
	
}

static void fix_delta(List pl, List ml)
{
	List l;
	for(l=ml;l;l=ListTail(l))
	{
		alg2_norm_delta(pl, ListFirst(l));
	}
}

Atom inv_color_eps(Atom);
Atom inv_color_k6(Atom);

static void fix_color_f(List ml, List prts)
{
	List l;
	int spcase=0;
	int hasceps=0;
	if(ListLength(prts)>2)
	{
	  Atom p1=CompoundArg1(ListFirst(prts));
	  Atom p2=CompoundArg1(ListFirst(ListTail(prts)));
	  Term prp=GetAtomProperty(p1,PROP_TYPE);
	  if(prp && is_compound(prp) && CompoundName(prp)==OPR_FIELD
	    && CompoundArg2(prp)==NewInteger(4) && CompoundArg1(prp)==p2)
	    spcase=1;
	  prp=GetAtomProperty(p2,PROP_TYPE);
	  if(prp && is_compound(prp) && CompoundName(prp)==OPR_FIELD
	    && CompoundArg2(prp)==NewInteger(4) && CompoundArg1(prp)==p1)
	    spcase=1;
	}
	
	for(l=ml;l;l=ListTail(l))
	{
		
		/*List l1;
		int s=1;
		for(l1=CompoundArgN(ListFirst(l),3);l1;l1=ListTail(l1))
			if(CompoundName(ListFirst(l1))==OPR_SPECIAL && (
				GetAtomProperty(CompoundArg1(ListFirst(l1)),A_COLOR)==A_COLOR_F ||
				GetAtomProperty(CompoundArg1(ListFirst(l1)),A_COLOR)==A_COLOR_EPS))
				{
					List il;
					int i1,i2,i3,tmp;
					il=CompoundArg2(ListFirst(l1));
					i1=IntegerValue(ListFirst(il));
					i2=IntegerValue(ListFirst(ListTail(il)));
					i3=IntegerValue(ListFirst(ListTail(ListTail(il))));
					if(i1>i2)
					{
						tmp=i2;	i2=i1;	i1=tmp; s=(-s);
					}
					if(i2>i3)
					{
						tmp=i3; i3=i2;	i2=tmp; s=(-s);
					}
					if(i1>i2)
					{
						tmp=i2;	i2=i1;	i1=tmp; s=(-s);
					}
					if(i2>i3)
					{
						tmp=i3; i3=i2;	i2=tmp; s=(-s);
					}
					ChangeList(il,NewInteger(i1));
					ChangeList(ListTail(il),NewInteger(i2));
					ChangeList(ListTail(ListTail(il)),NewInteger(i3));
				}
		if(s==-1)
		{
			int num;
			num=IntegerValue(CompoundArg1(CompoundArg1(ListFirst(l))));
			SetCompoundArg(CompoundArg1(ListFirst(l)),1,NewInteger(-num));
		}
		*/
		color_symm_f(0,ListFirst(l));
	}
	
	/* invert color_eps */
	
	for(l=ml;l;l=ListTail(l))
	{
		List l1;
		Term prp;
		for(l1=CompoundArgN(ListFirst(l),3);l1;l1=ListTail(l1))
			if(CompoundName(ListFirst(l1))==OPR_SPECIAL && 
				(prp=GetAtomProperty(CompoundArg1(ListFirst(l1)),A_COLOR)))
			{
			  if(prp==A_COLOR_EPS)
			  {
					SetCompoundArg(ListFirst(l1),1,
						inv_color_eps(CompoundArg1(ListFirst(l1))));
				hasceps=1;
			  }
			 if(prp==A_COLOR_K6)
					SetCompoundArg(ListFirst(l1),1,
						inv_color_k6(CompoundArg1(ListFirst(l1))));
			}
		if(hasceps && spcase) 
		{
		  Term rt=ConsumeCompoundArg(ListFirst(l),1);
		  Integer num=CompoundArg1(rt);
		  SetCompoundArg(rt,1,NewInteger(-IntegerValue(num)));
		  SetCompoundArg(ListFirst(l),1,rt);
		}
	}
}

static void fix_sp_ord(Term a2)
{
	List l1,l2;
	for(l1=CompoundArgN(a2,5);l1;l1=ListTail(l1))
	{
		Term m2;
		List l, gl=0;
		m2=ListFirst(l1);
		l2=ConsumeCompoundArg(m2,3);
	rp:	for(l=l2;l;l=ListTail(l))
			{
			Term t=CompoundArg1(ListFirst(l));
			if(t==A_GAMMA || t==A_GAMMAP || t==A_GAMMAM || t==A_GAMMA5)
				{
				gl=AppendLast(gl,ListFirst(l));
				ChangeList(l,0);
				l2=CutFromList(l2,l);
				goto rp;
				}
			}
		
		l2=SortedList(l2,sp_cmp);
		l2=ConcatList(l2,gl);
		SetCompoundArg(m2,3,l2);
	}	
}

static int wrt_mono(Term m2)
{
	int ret=0,n,d;
	List ms,l;
	ms=CopyTerm(CompoundArg2(m2));
	for(l=CompoundArgN(m2,3);l;l=ListTail(l))
		ms=AppendLast(ms,MakeCompound2(OPR_POW,ListFirst(l),NewInteger(1)));
	if(is_compound(CompoundArg1(m2)))
	{
		n=(int)IntegerValue(CompoundArg1(CompoundArg1(m2)));
		d=(int)IntegerValue(CompoundArg2(CompoundArg1(m2)));
		ret=printf("%d/%d",n,d);
	}
	else
	{
		n=(int)IntegerValue(CompoundArg1(m2));
		d=1;
		ret=printf("%d",n);
	}
	for(l=ms;l;l=ListTail(l))
	{
		Atom a;
		int p;
		a=CompoundArg1(ListFirst(l));
		p=(int)IntegerValue(CompoundArg2(ListFirst(l)));
		if(p>0)
			ret+=printf("*");
		else
			ret+=printf("/"),p=-p;
		ret+=printf("%s",AtomValue(a));
		if(p>1)
			ret+=printf("**%d",p);
	}
	FreeAtomic(ms);
	return ret;
}

static void compare_vertices(Term a2, Term a2cg, Term a2cf)
{
	List o, l1, l2;
	int mfound=0, pos;
	
	for(l1=CompoundArgN(a2cf,5);l1;l1=ListTail(l1))
	{
		o=ConsumeCompoundArg(ListFirst(l1),2);
		if(o) o=SortedList(o,prmcmp);
		SetCompoundArg(ListFirst(l1),2,o);
	}
	for(l1=CompoundArgN(a2cg,5);l1;l1=ListTail(l1))
	{
		o=ConsumeCompoundArg(ListFirst(l1),2);
		if(o) o=SortedList(o,prmcmp);
		SetCompoundArg(ListFirst(l1),2,o);
	}
	
	if(EqualTerms(a2cg,a2cf))
		return;
		
	a2=CopyTerm(a2);
	a2cf=CopyTerm(a2cf);
	a2cg=CopyTerm(a2cg);
	
	fix_mom_2(a2);
	fix_mom_2(a2cf);
	fix_mom_2(a2cg);
	
	alg2_symmetrize(a2);
	alg2_symmetrize(a2cf);
	alg2_symmetrize(a2cg);

	
	o=NewList();
	
	for(l1=CompoundArgN(a2cg,5);l1;l1=ListTail(l1))
	{
		pos=1;
		
		for(l2=CompoundArgN(a2cf,5);l2;l2=ListTail(l2))
			if(EqualTerms(ListFirst(l1),ListFirst(l2)))
				break;
			else
				pos++;
		if(l2)
			o=AppendLast(o,NewInteger(pos)),mfound++;
		else
			o=AppendLast(o,NewInteger(0));
	}
	
	if(mfound==ListLength(o) && mfound==ListLength(CompoundArgN(a2cf,5)))
	{
		FreeAtomic(o);
		FreeAtomic(a2);
		FreeAtomic(a2cf);
		FreeAtomic(a2cg);
		return;
	}


/*        puts("\n\n");
       WriteTerm(a2);puts(" : ori\n");
        WriteTerm(a2cg);puts(" : hc gen\n");
        WriteTerm(a2cf);puts(" : hc fou\n");
*/
	
	alg2_reduce(a2);
	alg2_reduce(a2cf);
	alg2_reduce(a2cg);

/*
	printf("\nori: ");WriteTerm(a2);puts("\n");
	printf("\ngen: ");WriteTerm(a2cg);puts("\n");
	printf("\nfou: ");WriteTerm(a2cf);puts("\n");
*/
	
	printf("CheckHerm: inconsistent conjugate vertices:\n\n");
	printf("    ");WriteVertex(CompoundArg1(a2));
	printf("                        ");WriteVertex(CompoundArg1(a2cf));puts("\n");

/*	DumpList(CompoundArgN(a2cf,5));
	DumpList(CompoundArgN(a2cg,5));
*/
	for(l1=o,l2=CompoundArgN(a2,5);l1;l1=ListTail(l1),l2=ListTail(l2))
	{
		int ll;
		if(ListFirst(l1)==NewInteger(0))
			continue;
		ll=printf("  ");
		ll+=wrt_mono(ListFirst(l2));
		WriteBlank(stdout,35-ll);
		printf(" <->   ");
		wrt_mono(ListNth(CompoundArgN(a2cf,5),(int)IntegerValue(ListFirst(l1))));
		puts("");
	}
	
	for(l1=o,l2=CompoundArgN(a2,5);l1;l1=ListTail(l1),l2=ListTail(l2))
	{
		int ll;
		if(ListFirst(l1)!=NewInteger(0))
			continue;
		ll=printf("  ");
		ll+=wrt_mono(ListFirst(l2));
		WriteBlank(stdout,35-ll);
		printf(" <->   (not found)\n");
	}
	
	pos=1;
	for(l1=CompoundArgN(a2cf,5);l1;l1=ListTail(l1),pos++)
	{
		int ll;
		if(ListMember(o,NewInteger(pos)))
			continue;
		ll=printf("   (not found)     ");
		WriteBlank(stdout,35-ll);
		printf(" <->   ");
		wrt_mono(ListFirst(l1));
		puts("");
	}
	
	if(verb_herm)
	{
		printf("Mononials of original vertex:\n");
		for(l1=CompoundArgN(a2,5);l1;l1=ListTail(l1))
		{
			wrt_mono(ListFirst(l1));puts("");
		}
		printf("\nMononials of generated conjugate vertex:\n");
		for(l1=CompoundArgN(a2cg,5);l1;l1=ListTail(l1))
		{
			wrt_mono(ListFirst(l1));puts("");
		}
		printf("\nMononials of found conjugate vertex:\n");
		for(l1=CompoundArgN(a2cf,5);l1;l1=ListTail(l1))
		{
			wrt_mono(ListFirst(l1));puts("");
		}
		puts("");
	}
	FreeAtomic(a2);
	FreeAtomic(a2cf);
	FreeAtomic(a2cg);
	FreeAtomic(o);
	puts("");

	
}

List cls_lagr_hook=0;
List sel_vrt(List, int,int);
	
Term ProcHermiticity(Term arg, Term ind)
	{
		List lagr, lagr1, ord;
		int ochk;
		List all_lagr;
		
		frcch=1;
		
		if(is_compound(arg))
		{
		if(CompoundArity(arg)!=1 || !is_list(CompoundArg1(arg)))
		{
			ErrorInfo(0);
			puts("Wring syntax in CheckHerm statement.");
			return 0;
		}
		all_lagr=sel_vrt(CompoundArg1(arg),1,0);
		}
		else
		if(cls_lagr_hook)
			{
			all_lagr=CopyTerm(cls_lagr_hook);
			for(lagr=all_lagr;lagr;lagr=ListTail(lagr))
				{
				List ml;
				SetCompoundArg(ListFirst(lagr),4,0);
				for(ml=CompoundArgN(ListFirst(lagr),5);ml;ml=ListTail(ml))
					{
					Integer ii=CompoundArg1(ListFirst(ml));
					SetCompoundArg(ListFirst(ml),1,MakeCompound2(
							OPR_DIV,ii,NewInteger(1)));
					}
				}
			}
		else
			all_lagr=all_vert_list();
			
		for(lagr=all_lagr;lagr;lagr=ListTail(lagr))
		{
			Term a2,a2c;
			List prt,l;
			int cf=0;
			a2=ListFirst(lagr);

			if(CompoundArgN(a2,4))
				continue;
			
			if(CompoundArg2(a2)==NewInteger(0) || 
					is_empty_list(CompoundArgN(a2,5)))
				continue;
			
			prt=CompoundArg1(a2);
			if(!TexOutput && ListLength(prt)>4)
				continue;
			
			for(l=prt;l;l=ListTail(l))
			{
				Atom p;
				Term pr;
				p=CompoundArg1(ListFirst(l));
				pr=GetAtomProperty(p,PROP_TYPE);
				if(!is_compound(pr))
				{
					cf=1;break;
					WriteVertex(prt);puts("internal error he01");
					exit(0);
				}
				if(CompoundName(pr)==OPR_PARTICLE)
                                { 
                                        Atom prp=CompoundArgN(pr,7);
                                        if(is_atom(prp) && (AtomValue(prp)[0]=='*'
                                            || AtomValue(prp)[1]=='*'))
					cf=1;
                                }
				if(CompoundName(pr)==OPR_FIELD &&
						CompoundArg2(pr)==NewInteger(5))
					cf=1;
			}
			
			if(cf)
				continue;

			if(!CompoundArg1(a2))
				continue;

						
			a2c=CopyTerm(a2);
			if((ord=inv_prt(a2c))==0)
			{
				FreeAtomic(all_lagr);
				frcch=0;
				return 0;
			}
			
			{
				List p1,p2;
				for(p1=CompoundArg1(a2),p2=CompoundArg1(a2c);
					p1&&p2;p1=ListTail(p1),p2=ListTail(p2))
						if(CompoundArg1(ListFirst(p1))
							!=CompoundArg1(ListFirst(p2)))
								break;
				if(p1==0)
				{
					FreeAtomic(a2c);
					FreeAtomic(ord);
					continue;
				}
		
				ochk=0;
				
				if(strcmp(AtomValue(CompoundArg1(ListFirst(p1))),
						AtomValue(CompoundArg1(ListFirst(p2))))<0)
				{
					ochk=1;
				}
		
			}

			fix_mom(ord,CompoundArgN(a2c,5));
			fix_spin(ord,CompoundArg1(a2c),CompoundArgN(a2c,5));
			fix_color_lambda(CompoundArgN(a2c,5));
			fix_ind_order(CompoundArg1(a2c),CompoundArgN(a2c,5));
			fix_color_f(CompoundArgN(a2c,5),CompoundArg1(a2c));
			fix_delta(CompoundArg1(a2c),CompoundArgN(a2c,5));
			fix_i(CompoundArgN(a2c,5));
			fix_gh_C(CompoundArg1(a2c),CompoundArgN(a2c,5));
			fix_mom_2(a2c);
			fix_sp_ord(a2c);
			FreeAtomic(ord);

			for(lagr1=all_lagr;lagr1;lagr1=ListTail(lagr1))
				if(EqualTerms(CompoundArg1(ListFirst(lagr1)),CompoundArg1(a2c)))
					break;

			if(is_empty_list(lagr1))
			{
				printf("CheckHerm: vertex ");
				WriteVertex(CompoundArg1(a2));
				printf(": conjugate ");
				WriteVertex(CompoundArg1(a2c));
				puts(" not found.");
				FreeAtomic(a2c);
				continue;
			}
	
			if(ochk)
			{
				FreeAtomic(a2c);
				continue;
			}
			compare_vertices(a2,a2c,ListFirst(lagr1));
			FreeAtomic(a2c);
		}
		
		FreeAtomic(all_lagr);
		frcch=0;
		return 0;
	}

List alg2_add_hermconj(List vert_list)
	{
	List lagr, ord, toadd;

	toadd=NewList();
	for(lagr=vert_list;lagr;lagr=ListTail(lagr))
		{
			Term a2,a2c;
			List prt,l;
			int cf=0;
			
			a2=ListFirst(lagr);
			if(CompoundArgN(a2,4))
				continue;
			
			if(is_empty_list(CompoundArg1(a2)))
				continue;
			
			if(CompoundArg2(a2)==NewInteger(0) || 
					is_empty_list(CompoundArgN(a2,5)))
				continue;
			
			prt=CompoundArg1(a2);
			for(l=prt;l;l=ListTail(l))
			{
				Atom p;
				Term pr;
				p=CompoundArg1(ListFirst(l));
				pr=GetAtomProperty(p,PROP_TYPE);
				if(!is_compound(pr))
				{
					puts("internal error he01");
					exit(0);
				}
				if(CompoundName(pr)==OPR_PARTICLE && 
						CompoundArgN(pr,7) == OPR_MLT)
					cf=1;
				if(CompoundName(pr)==OPR_FIELD &&
						CompoundArg2(pr)==NewInteger(5))
					cf=1;
			}
			
			if(cf)
			{
				ErrorInfo(445);
				printf("AddHermConj: can not process vertices with aux fields.\n");
				FreeAtomic(toadd);
				FreeAtomic(vert_list);
				return 0;
			}

			a2c=CopyTerm(a2);
			if((ord=inv_prt(a2c))==0)
				return 0;
			
			fix_mom(ord,CompoundArgN(a2c,5));
			fix_spin(ord,CompoundArg1(a2c),CompoundArgN(a2c,5));
			fix_color_lambda(CompoundArgN(a2c,5));
			fix_ind_order(CompoundArg1(a2c),CompoundArgN(a2c,5));
			fix_color_f(CompoundArgN(a2c,5),CompoundArg1(a2c));
			fix_delta(CompoundArg1(a2c),CompoundArgN(a2c,5));
			fix_i(CompoundArgN(a2c,5));
			fix_gh_C(CompoundArg1(a2c),CompoundArgN(a2c,5));
			fix_mom_2(a2c);
			fix_sp_ord(a2c);

			toadd=AppendLast(toadd,a2c);
		}
		
	vert_list=ConcatList(vert_list,toadd);
			
	return vert_list;
	}
