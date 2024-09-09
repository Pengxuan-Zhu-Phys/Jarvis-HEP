#include "lanhep.h"
#include <setjmp.h>
#include <string.h>
#include <ctype.h>

extern jmp_buf alg1_jmp_buf;
extern int alg1_recurse_level;

void alg1_derivp(Term);

/*
Term ProcAnti(Term p, Term i)
	{
	Term p1;
	p1=ConsumeCompoundArg(p,1);
	FreeAtomic(p);
	p=p1;
	
	if(!is_atom(p))
		{
		ErrorInfo(203);
		puts("wrong input in 'anti' call.");
		return 0;
		}
	
	p1=GetAtomProperty(p,PROP_TYPE);
	if(p1!=0 && is_compound(p1) && CompoundName(p1)==OPR_PARTICLE)
		{
		if(p==CompoundArg1(p1))
			return CompoundArg2(p1);
		else
			return CompoundArg1(p1);
		}
	if(p1!=0 && is_compound(p1) && CompoundName(p1)==OPR_LET)
		{
		p=gen_anti_name(p);
		p1=GetAtomProperty(p,PROP_TYPE);
		if(p1!=0 && is_compound(p1) && CompoundName(p1)==OPR_LET)
			return p;
		}
	ErrorInfo(203);
	puts("wrong input in 'ap' call.");
	return 0;
	}
	
Term ProcCC(Term p, Term ii)
	{
	Term p1;
	p1=ConsumeCompoundArg(p,1);
	FreeAtomic(p);
	p=p1;
	
	if(!is_atom(p))
		{
		ErrorInfo(204);
		puts("wrong input in 'cc' call.");
		return 0;
		}
	sprintf(mk_cc_buf,"%s.c",AtomValue(p));
	p=NewAtom(mk_cc_buf,0);
	if(!is_particle(p,NULL))
		{
		ErrorInfo(204);
		puts("wrong input in 'cc' call.");
		return 0;
		}
	return p;
	}	

	
*/
		
static int is_ferm(Atom p)
{
	Term prp;
	
	prp=GetAtomProperty(p,PROP_TYPE);
	if(prp==0)
		return 0;
	
        if(!is_compound(prp)) return 0;
        
	if(CompoundName(prp)==OPR_PARTICLE)
	{
		if(CompoundArgN(prp,4)==NewInteger(1))
			return 1;
		else
			return 0;
	}
	
	if(CompoundArg2(prp)==NewInteger(4))
		return 1;
	
	return 0;
}

static int is_maj(Atom p)
{
	Term prp;
	
	prp=GetAtomProperty(p,PROP_TYPE);
	if(prp==0)
		return 0;
	
         if(!is_compound(prp)) return 0;
         
	if(CompoundName(prp)==OPR_PARTICLE)
	{
		if((CompoundArgN(prp,4)==NewInteger(1)||CompoundArgN(prp,4)==NewInteger(3)) &&
				CompoundArg1(prp)==CompoundArg2(prp))
			return 1;
		else
			return 0;
	}
	
	if(CompoundArg2(prp)!=NewInteger(4))
		return 0;
	
	p=CompoundArg1(prp);
	prp=GetAtomProperty(p,PROP_TYPE);
	if(prp==0)
		return 0;
	if((CompoundArgN(prp,4)==NewInteger(1)||CompoundArgN(prp,4)==NewInteger(3)) &&
			CompoundArg1(prp)==CompoundArg2(prp))
		return 1;
		
	return 0;
}

static Atom anti_maj(Atom p)
{
	Term prp;
	char buf[32];
	
	prp=GetAtomProperty(p,PROP_TYPE);
	if(prp==0)
		return 0;
	
	if(CompoundName(prp)==OPR_PARTICLE)
	{
		sprintf(buf,"%s.c",AtomValue(p));
		return NewAtom(buf,0);
	}
	
	return CompoundArg1(prp);
}

Term inv_color_eps(Atom);	
int  alg1_anti_modf;
		
void alg1_anti(Term a1)
{
	List l1,l2,l3;
	
	for(l1=CompoundArg2(a1);l1;l1=ListTail(l1))
	{
		Term g, r1, r2;
		
		g=CompoundArg1(ListFirst(l1));
		r1=ConsumeCompoundArg(g,1);
		r2=ConsumeCompoundArg(g,2);
		if(r1!=r2)
			alg1_anti_modf=1;
		SetCompoundArg(g,1,r2);
		SetCompoundArg(g,2,r1);
	}
	
	for(l1=CompoundArg1(a1);l1;l1=ListTail(l1))
	{
		Term m1;
		m1=ListFirst(l1);
		
		if(IntegerValue(CompoundArg2(m1))<0)
		{
			alg1_anti_modf=1;
			SetCompoundArg(m1,1,NewInteger(-IntegerValue(CompoundArg1(m1))));
		}
		
		for(l2=CompoundArgN(m1,3);l2;l2=ListTail(l2))
		{
			if(CompoundName(ListFirst(l2))==OPR_PARAMETER)
			{
				Atom a=GetAtomProperty(CompoundArg2(ListFirst(l2)),A_ANTI);
				if(a)
					{
					SetCompoundArg(ListFirst(l2),2,a);
					alg1_anti_modf=1;
					}
				continue;
			}
			
			if(CompoundName(ListFirst(l2))!=OPR_SPECIAL)
			for(l3=CompoundArg1(ListFirst(l2));l3;l3=ListTail(l3))
			{
				Term g, r1, r2;

				g=CompoundArg1(ListFirst(l3));
				r1=ConsumeCompoundArg(g,1);
				r2=ConsumeCompoundArg(g,2);
				if(r1!=r2)
					alg1_anti_modf=1;
				SetCompoundArg(g,1,r2);
				SetCompoundArg(g,2,r1);
			}
			
			if(CompoundName(ListFirst(l2))==OPR_FIELD)
			{
				Atom p, ap;
				p=CompoundArg2(ListFirst(l2));
				ap=GetAtomProperty(p,A_ANTI);
				if(ap==0)
				{
					ErrorInfo(208);
					printf("anti: no antiparticle for '%s'.\n",AtomValue(p));
					if(alg1_recurse_level)
						longjmp(alg1_jmp_buf,1);
				}
				
				if(ap==p && is_maj(p))
					ap=anti_maj(p);
				
				if(p!=ap)
					alg1_anti_modf=1;
				SetCompoundArg(ListFirst(l2),2,ap);
				continue;
			}
			
			if(CompoundName(ListFirst(l2))==OPR_LET)
			{
				Atom p, ap;
				p=CompoundArg2(ListFirst(l2));
				ap=GetAtomProperty(p,A_ANTI);
				if(ap==0)
				{
					ErrorInfo(208);
					printf("anti: no anti-multiplet for '%s'.\n",AtomValue(p));
					if(alg1_recurse_level)
						longjmp(alg1_jmp_buf,1);
				}
				
				if(p!=ap)
					alg1_anti_modf=1;
				SetCompoundArg(ListFirst(l2),2,ap);
				continue;
			}
			
			if(CompoundName(ListFirst(l2))==OPR_WILD)
			{
				for(l3=CompoundArg2(ListFirst(l2));l3;l3=ListTail(l3))
					alg1_anti(ListFirst(l3));
				continue;
			}
			
			if(CompoundName(ListFirst(l2))==OPR_SPECIAL)
			{
				Atom s;
				s=CompoundArg2(ListFirst(l2));
				if(s==A_GAMMA || s==A_GAMMA5 || s==A_GAMMAP || s==A_GAMMAM)
				{
					Label la1,la2;
					l3=CompoundArg1(ListFirst(l2));
					la1=CompoundArg2(ListFirst(l3));
					la2=CompoundArg2(ListFirst(ListTail(l3)));
					SetCompoundArg(ListFirst(l3),2,la2);
					SetCompoundArg(ListFirst(ListTail(l3)),2,la1);
					if(s==A_GAMMAP)
						SetCompoundArg(ListFirst(l2),2,A_GAMMAM);
					if(s==A_GAMMAM)
						SetCompoundArg(ListFirst(l2),2,A_GAMMAP);
					if(s==A_GAMMA5)
						SetCompoundArg(m1,1,
								NewInteger(-IntegerValue(CompoundArg1(m1))));
					alg1_anti_modf=1;
					continue;
				}
				
				if(s==A_COLOR_LAMBDA)
				{
					Label la1,la2;
					l3=CompoundArg1(ListFirst(l2));
					la1=CompoundArg2(ListFirst(l3));
					la2=CompoundArg2(ListFirst(ListTail(l3)));
					SetCompoundArg(ListFirst(l3),2,la2);
					SetCompoundArg(ListFirst(ListTail(l3)),2,la1);
					alg1_anti_modf=1;
					continue;
				}
				
				if(GetAtomProperty(s,A_COLOR)==A_COLOR_EPS)
				{
					SetCompoundArg(ListFirst(l2),2,inv_color_eps(s));
					alg1_anti_modf=1;
					continue;
				}
				
				if(s==A_DELTA)
					continue;
				
				ErrorInfo(209);
				printf("anti: can not deal with '%s'.\n",AtomValue(s));
				if(alg1_recurse_level)
					longjmp(alg1_jmp_buf,1);
				return;
			}
			
			ErrorInfo(209);
			printf("anti: can not deal with '%s'.\n",
					AtomValue(CompoundName(ListFirst(l2))));
			if(alg1_recurse_level)
				longjmp(alg1_jmp_buf,1);
			return;
		}
		
	}
	
}
/*
void alg1_cc(Term a1)
{
	List l1,l2,l3;
	
	for(l1=CompoundArg2(a1);l1;l1=ListTail(l1))
	{
		Term g, r1, r2;
		
		g=CompoundArg1(ListFirst(l1));
		r1=ConsumeCompoundArg(g,1);
		r2=ConsumeCompoundArg(g,2);
		if(r1!=r2)
			alg1_anti_modf=1;
		SetCompoundArg(g,1,r2);
		SetCompoundArg(g,2,r1);
	}
	
	for(l1=CompoundArg1(a1);l1;l1=ListTail(l1))
	{
		Term m1;
		m1=ListFirst(l1);
		
		if(IntegerValue(CompoundArg2(m1))<0)
		{
			alg1_anti_modf=1;
			SetCompoundArg(m1,1,NewInteger(-IntegerValue(CompoundArg1(m1))));
		}
		
		for(l2=CompoundArgN(m1,3);l2;l2=ListTail(l2))
		{
			if(CompoundName(ListFirst(l2))==OPR_PARAMETER)
			{
				Atom a=GetAtomProperty(CompoundArg2(ListFirst(l2)),A_ANTI);
				if(a)
					{
					SetCompoundArg(ListFirst(l2),2,a);
					alg1_anti_modf=1;
					}
				continue;
			}
			
			if(CompoundName(ListFirst(l2))!=OPR_SPECIAL)
			for(l3=CompoundArg1(ListFirst(l2));l3;l3=ListTail(l3))
			{
				Term g, r1, r2;

				g=CompoundArg1(ListFirst(l3));
				r1=ConsumeCompoundArg(g,1);
				r2=ConsumeCompoundArg(g,2);
				if(r1!=r2)
					alg1_anti_modf=1;
				SetCompoundArg(g,1,r2);
				SetCompoundArg(g,2,r1);
			}
			
			if(CompoundName(ListFirst(l2))==OPR_FIELD || 
					(CompoundName(ListFirst(l2))==OPR_LET &&
					 is_particle(CompoundArg2(ListFirst(l2)),NULL)))
			{
				Atom p, ap;
				p=CompoundArg2(ListFirst(l2));
				ap=GetAtomProperty(p,A_ANTI);
				if(ap==0)
				{
					ErrorInfo(208);
					printf("anti: no antiparticle for '%s'.\n",AtomValue(p));
					if(alg1_recurse_level)
						longjmp(alg1_jmp_buf,1);
				}
				
				if(is_ferm(p))
					ap=anti_maj(p);
				
				if(p!=ap)
					alg1_anti_modf=1;
				SetCompoundArg(ListFirst(l2),2,ap);
				continue;
			}
			
			if(CompoundName(ListFirst(l2))==OPR_LET)
			{
				
				ErrorInfo(208);
				printf("cc: let-subst '%s' is forbidden.\n",
						AtomValue(CompoundArg2(ListFirst(l2))));
				if(alg1_recurse_level)
					longjmp(alg1_jmp_buf,1);
				continue;
			}
			
			if(CompoundName(ListFirst(l2))==OPR_WILD)
			{
				for(l3=CompoundArg2(ListFirst(l2));l3;l3=ListTail(l3))
					alg1_cc(ListFirst(l3));
				continue;
			}
			
			if(CompoundName(ListFirst(l2))==OPR_SPECIAL)
			{
				Atom s;
				s=CompoundArg2(ListFirst(l2));
				if(s==A_GAMMA || s==A_GAMMA5 || s==A_GAMMAP || s==A_GAMMAM)
				{
					Label la1,la2;
					l3=CompoundArg1(ListFirst(l2));
					la1=CompoundArg2(ListFirst(l3));
					la2=CompoundArg2(ListFirst(ListTail(l3)));
					SetCompoundArg(ListFirst(l3),2,la2);
					SetCompoundArg(ListFirst(ListTail(l3)),2,la1);
					if(s==A_GAMMA)
						SetCompoundArg(m1,1,
								NewInteger(-IntegerValue(CompoundArg1(m1))));
					alg1_anti_modf=1;
					continue;
				}
				
				if(s==A_COLOR_LAMBDA)
				{
					Label la1,la2;
					l3=CompoundArg1(ListFirst(l2));
					la1=CompoundArg2(ListFirst(l3));
					la2=CompoundArg2(ListFirst(ListTail(l3)));
					SetCompoundArg(ListFirst(l3),2,la2);
					SetCompoundArg(ListFirst(ListTail(l3)),2,la1);
					alg1_anti_modf=1;
					continue;
				}
				
				if(GetAtomProperty(s,A_COLOR)==A_COLOR_EPS)
				{
					SetCompoundArg(ListFirst(l2),2,inv_color_eps(s));
					alg1_anti_modf=1;
					continue;
				}
				
				if(s==A_DELTA)
					continue;
				
				ErrorInfo(209);
				printf("cc: can not deal with '%s'.\n",AtomValue(s));
				if(alg1_recurse_level)
					longjmp(alg1_jmp_buf,1);
				return;
			}
			
			ErrorInfo(209);
			printf("cc: can not deal with '%s'.\n",
					AtomValue(CompoundName(ListFirst(l2))));
			if(alg1_recurse_level)
				longjmp(alg1_jmp_buf,1);
			return;
		}
		
	}
	
}
*/

static Term app_ind(Term t, List ind)
{
	if(ind==0)
		return t;
	return il_to_caret(t,ind);
}
	
Term ProcAnti(Term p, Term i)
{
	Term t1, t2;
	
	Term tsv;
	
	if(!is_compound(p)||CompoundArity(p)!=1)
	{
		ErrorInfo(202);
		printf("'anti' should have one argument.\n");
		if(alg1_recurse_level)
			longjmp(alg1_jmp_buf,1);
		return 0;
	}
	
	tsv=CopyTerm(p);
	
	t1=ConsumeCompoundArg(p,1);
	FreeAtomic(p);
	
	if(is_atom(t1) && is_maj(t1))
		return app_ind(anti_maj(t1),i);
	
	if(is_atom(t1) && (t2=GetAtomProperty(t1,A_ANTI)))
	{
		return app_ind(t2,i);
	}
	
	t1=ExprTo1(t1);
	
	if(t1==0)
	{
		FreeAtomic(tsv);
		return 0;
	}
	
	if(CompoundArg1(t1)==0)
	{
		FreeAtomic(t1);
		return NewInteger(0);
	}
	
	alg1_anti_modf=0;
	
	alg1_anti(t1);
	
	if(CompoundArg1(t1)==0)
	{
		FreeAtomic(t1);
		FreeAtomic(tsv);
		return NewInteger(0);
	}
	
	if(!alg1_anti_modf)
	{
		WarningInfo(207);
		WriteTerm(tsv);
		printf(" is identical to its argument.\n");
	}
	
	FreeAtomic(tsv);
	
	if(i)
		t1=il_to_caret(t1,i);
		
	return t1;
	
}
	
Term ProcCC(Term p, Term ii)
	{
	Term t1, tsv;
	
	if(!is_compound(p)||CompoundArity(p)!=1)
	{
		ErrorInfo(202);
		printf("'cc' should have one argument.\n");
		if(alg1_recurse_level)
			longjmp(alg1_jmp_buf,1);
		return 0;
	}
	tsv=CopyTerm(p);
	
	t1=ConsumeCompoundArg(p,1);
	FreeAtomic(p);
	
	if(is_atom(t1) && is_particle(t1,NULL))
		return app_ind(anti_maj(t1),ii);
	
	ErrorInfo(0);
	puts("Argument of cc() function must be a fermion.\n");
	return 0;
	/*
	t1=ExprTo1(t1);
	
	if(t1==0)
	{
		FreeAtomic(tsv);
		return 0;
	}
	
	if(CompoundArg1(t1)==0)
	{
		FreeAtomic(t1);
		return NewInteger(0);
	}
	
	alg1_anti_modf=0;
	
	alg1_cc(t1);
	
	if(CompoundArg1(t1)==0)
	{
		FreeAtomic(t1);
		FreeAtomic(tsv);
		return NewInteger(0);
	}
	
	if(!alg1_anti_modf)
	{
		WarningInfo(207);
		WriteTerm(tsv);
		printf(" is identical to its argument.\n");
	}
	
	FreeAtomic(tsv);
	
	if(ii)
		t1=il_to_caret(t1,ii);
		
	return t1;
	*/
	}	

	
Term ProcDeriv(Term t, Term ind)
{
	Term il;
	Term t1;
	
	if(ListLength(ind)>1)
	{
		ErrorInfo(209);
		puts("deriv: too many indices");
		return 0;
	}
	
	if(ind)
		il=ListFirst(ind);
	else
		il=0;
	
	
	if(!is_compound(t) || CompoundArity(t)!=1)
	{
		ErrorInfo(210);
		puts("deriv: wrong syntax");
		return 0;
	}
	
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	
	if(il==0)
		t=A_DERIV_S;
	else
		t=MakeCompound2(OPR_CARET,A_DERIV_S,il);
	
	t1=MakeCompound2(OPR_MLT,t1,A_MOMENT_E);
	
	t=MakeCompound2(OPR_MLT,t,t1);
	return t;
	
}

static Term alg1_derivp11(Term m1, Term ml, int pos, List di)
{
	List ret=0;
	int i=0;
	List l1;
	Term t;
		
	t=MakeCompound2(OPR_SPECIAL,di,A_MOMENT);
	
	for(l1=ml;l1;l1=ListTail(l1))
	{
		i++;
		if(i==pos)
			ret=AppendLast(ret,t);
		ret=AppendLast(ret,ListFirst(l1));
	}
	
	RemoveList(ml);
	SetCompoundArg(m1,3,ret);
	return m1;
}
	

static void alg1_derivp1(Term m1, List *ml)
{
	List l1, l2, ml1, di;
	int pos=0;
	int fldpos[32];
	int fldno=0;
	int p1=0, p2=0;
	List mln;

	
	for(l1=CompoundArgN(m1,3);l1;l1=ListTail(l1))
	{
		pos++;
		if(CompoundName(ListFirst(l1))==OPR_SPECIAL && 
				CompoundArg2(ListFirst(l1))==A_MOMENT_S)
			p1=pos;
	}
			
	if(p1==0)
	{
		*ml=AppendFirst(*ml,m1);
		return;
	}
	
	pos=0;
	for(l1=CompoundArgN(m1,3);l1;l1=ListTail(l1))
	{
		Term t;
		
		pos++;
		if(pos<p1)
			continue;
		
		t=ListFirst(l1);
		
		if(p1 && CompoundName(t)==OPR_FIELD && CompoundArg2(t)!=A_VEV)
		{
			fldpos[fldno]=pos;
			fldno++;
		}
		
		
		if(CompoundName(t)==OPR_SPECIAL && CompoundArg2(t)==A_MOMENT_E)
		{
			if(p1==0)
			{
				puts("Internal error a1drv01");
				return;
			}
			p2=pos;
			break;
		}
	}
	
		
	if(fldno==0)
	{
		FreeAtomic(m1);
		return;
	}
	
	ml1=ConsumeCompoundArg(m1,3);
	
	if(p2)
	{
		l2=ListNthList(ml1,p2);
		ml1=CutFromList(ml1,l2);
	}
	
	l2=ListNthList(ml1,p1);
	di=ConsumeCompoundArg(ListFirst(l2),1);
	ml1=CutFromList(ml1,l2);
	
	mln=NewList();
	
	for(pos=0;pos<fldno;pos++)
		mln=AppendLast(mln,alg1_derivp11(
				CopyTerm(m1),CopyTerm(ml1),fldpos[pos]-1,CopyTerm(di)));
	
	FreeAtomic(m1);
	FreeAtomic(ml1);
	FreeAtomic(di);
	
/*	DumpList(mln);*/
	
	
	di=MakeCompound2(A_ALG1,mln,0);
	alg1_derivp(di);
	mln=ConsumeCompoundArg(di,1);
	FreeAtomic(di);
	*ml=ConcatList(mln,*ml);
	
}

void alg1_derivp(Term a1)
{
	
	List l1,l2,l3,l4;
	
	l1=ConsumeCompoundArg(a1,1);
	l2=NewList();
	
	for(l3=l1;l3;l3=ListTail(l3))
		alg1_derivp1(ListFirst(l3),&l2);
	
	RemoveList(l1);
	SetCompoundArg(a1,1,l2);
}

static Term xNewAtom(char *c, int n)
{
	int d=1;
	for(n=0;c[n];n++)
		if(!( (n==0&&c[n]=='-')||isdigit(c[n] )))
			d=0;
	if(d==0)
		return NewAtom(c,0);
	sscanf(c,"%d",&n);
	return NewInteger(n);
}


static void repla(Term t, char c, Term r)
{
	int i;
	char cbuf[1028];
	if(is_compound(t))
		for(i=1;i<=CompoundArity(t);i++)
		{
			rpt:
			if(is_atom(CompoundArgN(t,i)))
			{
				int j,lr;
				char *a=AtomValue(CompoundArgN(t,i));
				Atom na;
				if(!is_compound(r) && a[0]=='_'&&a[1]==c&&a[2]==0)
				{
					SetCompoundArg(t,i,r);
					continue;
				}
				for(j=0;a[j];j++)
				{
					cbuf[j]=a[j];
					if(a[j]=='_'&&a[j+1]==c)
						break;
				}
				if(a[j]==0)
				{
					if(isalpha(c) && is_atom(r) && isalpha(AtomValue(r)[0])
							&& AtomValue(r)[1]==0)
					{
						for(j=0;a[j];j++)
						{
							cbuf[j]=a[j];
							if(a[j]=='_' && a[j+1]==toupper(c))
								break;
						}
						if(a[j]==0)
							continue;
						cbuf[j]=toupper(AtomValue(r)[0]);
						strcpy(cbuf+j+1,a+j+2);
						na=xNewAtom(cbuf,0);
						SetCompoundArg(t,i,na);
					}
					continue;
				}
				lr=sWriteTerm(cbuf+j,r);
				strcpy(cbuf+j+lr,a+j+2);
				na=xNewAtom(cbuf,0);
				SetCompoundArg(t,i,na);
				goto rpt;
			}
			else
				repla(CompoundArgN(t,i),c,r);
		}
		
	if(is_list(t))
		for(;t;t=ListTail(t))
		{
			rpt2:
			if(is_atom(ListFirst(t)))
			{
				int j,lr;
				char *a=AtomValue(ListFirst(t));
				Atom na;
				if(!is_compound(r) && a[0]=='_'&&a[1]==c&&a[2]==0)
				{
					ChangeList(t,r);
					continue;
				}
				for(j=0;a[j];j++)
				{
					cbuf[j]=a[j];
					if(a[j]=='_'&&a[j+1]==c)
						break;
				}
				if(a[j]==0)
				{
					if(isalpha(c) && is_atom(r) && isalpha(AtomValue(r)[0])
							&& AtomValue(r)[1]==0)
					{
						for(j=0;a[j];j++)
						{
							cbuf[j]=a[j];
							if(a[j]=='_' && a[j+1]==toupper(c))
								break;
						}
						if(a[j]==0)
							continue;
						cbuf[j]=toupper(AtomValue(r)[0]);
						strcpy(cbuf+j+1,a+j+2);
						na=xNewAtom(cbuf,0);
						ChangeList(t,na);
					}
					continue;
				}
				lr=sWriteTerm(cbuf+j,r);
				strcpy(cbuf+j+lr,a+j+2);
				na=xNewAtom(cbuf,0);
				ChangeList(t,na);
				goto rpt2;
			}
			else
				repla(ListFirst(t),c,r);
		}
}
				 
				

Term ProcIn(Term t, List ind)
{
	Term sub, fun;
	Atom suba;
	char *subv;
	if(!is_compound(t)|| CompoundArity(t)!=2)
	{
		ErrorInfo(0);
		puts("Bad 'in' statement");
		return 0;
	}
	sub=ConsumeCompoundArg(t,1);
	fun=ConsumeCompoundArg(t,2);
	sub=ProcessAlias(sub);
	FreeAtomic(t);
	if(is_compound(sub) && CompoundArity(sub)==2 && CompoundName(sub)==OPR_COMMA)
	{
		Term s1,s2;
		s1=ConsumeCompoundArg(sub,1);
		s2=ConsumeCompoundArg(sub,2);
		FreeAtomic(sub);
		sub=s1;
		fun=MakeCompound2(OPR_IN,s2,fun);
	}
	if(!is_compound(sub) || CompoundArity(sub)!=2 ||
			CompoundName(sub)!=OPR_EQSIGN)
	{
		ErrorInfo(0);
		puts("in statement: bad subst rule");
		return 0;
	}
	suba=CompoundArg1(sub);
	if(!is_atom(suba) || (subv=AtomValue(suba))[0]!='_' || subv[2]!=0)
	{
		ErrorInfo(0);
		puts("in statement: bad subst name");
		return 0;
	}
	t=ConsumeCompoundArg(sub,2);
	FreeAtomic(sub);
	if(is_atom(t) || is_integer(t))
		t=MakeList1(t);
	if(is_compound(t)&&CompoundArity(t)==2 && CompoundName(t)==OPR_MINUS &&
			is_integer(CompoundArg1(t)) && is_integer(CompoundArg2(t)))
	{
		int i,i1,i2;
		i1=IntegerValue(CompoundArg1(t));
		i2=IntegerValue(CompoundArg2(t));
		for(i=i1;i<=i2;i++)
		{
			Term t1=CopyTerm(fun);
			repla(t1,subv[1],NewInteger(i));
			/*WriteTerm(t1);puts("");*/
			ProcessTerm(t1);
		}
	}
	else if(is_list(t))
	{
		List l;
		for(l=t;l;l=ListTail(l))
		{
			Term t1=CopyTerm(fun);
			repla(t1,subv[1],ListFirst(l));
			/*WriteTerm(t1);puts("");*/
			ProcessTerm(t1);
		}
	}
	else
	{
		ErrorInfo(0);
		puts("in statement: bad subst");
	}
	
	FreeAtomic(t);
	FreeAtomic(fun);
	return 0;

}

static List el_list=0;

Term ProcExtlib(Term t, Term ind)
{
	List l;
	if(!is_compound(t) || CompoundArity(t)!=1)
	{
		ErrorInfo(0);
		puts("wrong syntax in extlib statement");
		return 0;
	}
	t=ConsumeCompoundArg(t,1);
	t=CommaToList(t);
	for(l=t;l;l=ListTail(l))
	{
		if(!is_atom(ListFirst(l)))
		{
			ErrorInfo(9);
			puts("wrong syntax in extlib statement");
			return 0;
		}
	}
	el_list=t;
	return 0;
}

static char wpbuf[128];

void WriteExtlib(int fno, char *name)
{
	FILE *f;
	List l;
	if(el_list==0 || TexOutput || FAOutput)
		return;
	
	if(OutputDirectory!=NULL)
		sprintf(wpbuf,"%s/extlib%d.mdl",OutputDirectory,fno);
	else
		sprintf(wpbuf,"extlib%d.mdl",fno);
	f=fopen(wpbuf,"w");
	if(f==NULL)
		{
		printf("Can not open file \'%s\' for writing.\n",wpbuf);
		perror("");
		return;
		}
		
	fprintf(f,"%s\nLibraries \n",name);
	fprintf(f,
"External libraries  and citation                                      <|\n");
	for(l=el_list;l;l=ListTail(l))
		fprintf(f,"%s\n",AtomValue(ListFirst(l)));
	fprintf(f,
"========================================================================\n");
	fclose(f);
	
	
}

static List cp_list=0;

Term ProcCpart(Term t, Term ind)
{
	
	List l;
	if(!is_compound(t) || CompoundArity(t)!=1)
	{
		ErrorInfo(0);
		puts("wrong syntax in cpart statement");
		return 0;
	}
	t=ConsumeCompoundArg(t,1);
	t=CommaToList(t);
	
	for(l=t;l;l=ListTail(l))
	{
		Term u=ListFirst(l);
		List l2;
		if(!is_compound(u) || CompoundArity(u)!=2|| 
				!is_atom(CompoundArg1(u)) || !is_list(CompoundArg2(u)))
		{
			ErrorInfo(9);
			puts("wrong syntax in cpart statement");
			return 0;
		}
		for(l2=CompoundArg2(u);l2;l2=ListTail(l2))
			if(!is_atom(ListFirst(l2))||!is_particle(ListFirst(l2),NULL))
			{
				ErrorInfo(0);
				printf("cpart '%s': '",AtomValue(CompoundArg1(u)));
				WriteTerm(ListFirst(l2));
				puts("' is not a particle");
			}
		
	}
	
	cp_list=t;
    return 0;
}


void WriteCpart(int fno, char *name)
{
	FILE *f;
	List l;
	if(cp_list==0 || TexOutput || FAOutput)
		return;
	
	if(OutputDirectory!=NULL)
		sprintf(wpbuf,"%s/cpart%d.mdl",OutputDirectory,fno);
	else
		sprintf(wpbuf,"cpart%d.mdl",fno);
	f=fopen(wpbuf,"w");
	if(f==NULL)
		{
		printf("Can not open file \'%s\' for writing.\n",wpbuf);
		perror("");
		return;
		}
		
	fprintf(f,"%s\n Composite \n",name);
	fprintf(f,"Abr  |> elementary particles               <|> Comment                        <|\n");
	for(l=cp_list;l;l=ListTail(l))
	{
		int sp;
		List l2;
		sp=fprintf(f,"%s",AtomValue(CompoundArg1(ListFirst(l))));
		WriteBlank(f,5-sp);
		fprintf(f,"|");
		sp=0;
		for(l2=CompoundArg2(ListFirst(l));l2;l2=ListTail(l2))
			sp+=fprintf(f,"%s%c",AtomValue(ListFirst(l2)),ListTail(l2)?',':' ');
		WriteBlank(f,38-sp);
		fprintf(f,"|                                   \n");
	}
		
	
	fclose(f);

	
}


Term ProcPrm(Term t, Term ind)
{
	if(!is_compound(t) || CompoundArity(t)!=1)
	{
		ErrorInfo(0);
		puts("prm: bad syntax.");
		return 0;
	}
	if(!is_atom(CompoundArg1(t)))
	{
		ErrorInfo(0);
		puts("prm: argument must be inside quotation marks.");
		return 0;
	}
	
	SetAtomProperty(CompoundArg1(t),PROP_TYPE,OPR_PARAMETER);
	SetAtomProperty(CompoundArg1(t),A_DUMMY_PRM,NewInteger(0));
	return CompoundArg1(t);
	
}

