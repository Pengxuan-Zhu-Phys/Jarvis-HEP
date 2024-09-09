#include "lanhep.h"
#include <string.h>
#include <ctype.h>

static int cls_no=1;

static Atom cls_names[128];
static List cls_membs[128];
static char cls_ini[128];
static List cls_lens=0;

static List mtr_known=0;
static List mtr_no=0, dmtr_no=0, mtr_fit=0;
List mtr_visi=0;

int opclsreal=0, oprcmatr=1, opclsred=1;


Term ProcMtrOut(Term t, Term ind)
{
 char cbuf[16];
 List l;
 if(!is_compound(t) || CompoundArity(t)!=1 || !is_list(CompoundArg1(t)))
 {
   ErrorInfo(0);
   puts("wrong syntax.");
   return 0;
 }
 for(l=CompoundArg1(t);l;l=ListTail(l))
 {
   sprintf(cbuf,"MTR%03d",IntegerValue(ListFirst(l)));
   mtr_visi=AppendLast(mtr_visi,NewAtom(cbuf,0));
 }
 return 0; 
}

void prt2cls(Atom *p)
{
	Term prop=GetAtomProperty(*p,OPR_CLASS);
	if(prop && is_compound(prop))
	{
		*p=cls_names[IntegerValue(CompoundArg1(prop))];
	}
}

extern int lterm_used;

Term  ProcClass(Term cl, Term ind)
{
	Term t=ConsumeCompoundArg(cl,1);
	Atom clnm, anti=0;
	List clmem, l1;
	int i, ms=0,noms=0,spin=-1, ctype=0;
	char cbuf[128];
	FreeAtomic(cl);
	
	for(i=0;i<128;i++) cls_ini[i]=0;
	if(lterm_used)
	{
		ErrorInfo(0);
		puts("class: must be before the first 'lterm' statement.");
	}
	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==OPR_COMMA)
		{
		Term a1,a2;
		a1=ConsumeCompoundArg(t,1);
		a2=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		ProcClass(MakeCompound1(OPR_CLASS,a1),0);
		ProcClass(MakeCompound1(OPR_CLASS,a2),0);
		return 0;
		}

	if(!is_compound(t) || CompoundArity(t)!=2 || CompoundName(t)!=OPR_EQSIGN
			|| !is_atom(CompoundArg1(t)) || !is_list(CompoundArg2(t))
			|| CompoundArg2(t)==0 )
		{
		ErrorInfo(325);
		printf("bad argument in class statement\n");
		return 0;
		}

	clnm=CompoundArg1(t);
	clmem=ConsumeCompoundArg(t,2);
	FreeAtomic(t);
	if(GetAtomProperty(clnm,PROP_TYPE))
	{
		ErrorInfo(325);
		printf("class: name %s is already in use.\n",AtomValue(clnm));
		return 0;
	}
	sprintf(cbuf,"%s",AtomValue(clnm));
	if(strlen(cbuf)<3 || !isalpha(cbuf[0]) || !isalpha(cbuf[1]) ||
		!isalpha(cbuf[2]))
	{
		ErrorInfo(325);
		printf("class: name '%s': must be 3 letters.\n",AtomValue(clnm));
		return 0;
	}
	
	
	t=GetAtomProperty(ListFirst(clmem),PROP_INDEX);
	for(l1=clmem;l1;l1=ListTail(l1))
	{
		Term prp, bp;
		int ttype;
		if(!is_atom(ListFirst(l1)) || !is_particle(ListFirst(l1),0))
		{
			ErrorInfo(325);
			printf("class: '");WriteTerm(ListFirst(l1));
			printf("' is not a particle.\n");
			return 0;
		}
		if(!EqualTerms(t,GetAtomProperty(ListFirst(l1),PROP_INDEX)))
		{
			ErrorInfo(325);
			printf("class %s: ",AtomValue(clnm));
			printf("'%s' and '%s' have different indices.\n",
					AtomValue(ListFirst(clmem)),AtomValue(ListFirst(l1)));
			return 0;
		}
		if(GetAtomProperty(ListFirst(l1),OPR_CLASS))
		{
			ErrorInfo(325);
			printf("class %s: ",AtomValue(clnm));
			printf("'%s' is already a member of another class.\n",
					AtomValue(ListFirst(l1)));
			return 0;
		}
		bp=ListFirst(l1);
		prp=GetAtomProperty(bp,PROP_TYPE);
		if(CompoundName(prp)==OPR_FIELD && CompoundArg2(prp)==NewInteger(1))
			{
			bp=CompoundArg1(prp);
			prp=GetAtomProperty(bp,PROP_TYPE);
			}
		if(CompoundName(prp)!=OPR_PARTICLE)
		{
			ErrorInfo(325);
			printf("class: '");WriteTerm(ListFirst(l1));
			printf("' can not be class member.\n");
			return 0;		
		}
		if(CompoundArg1(prp)==CompoundArg2(prp))
			ttype=1;
		else if(CompoundArg1(prp)==bp)
			ttype=2;
		else
			ttype=3;
		if(ctype==0)
			ctype=ttype;
		else if(ctype!=ttype)
		{
			ErrorInfo(325);
			printf("class: '");WriteTerm(clnm);
			printf("': can not mix charged and neutral particles.\n");
			return 0;		
		}
		
			
	}
	
	t=AppendLast(CopyTerm(t),MakeCompound2(A_I,MakeCompound2(OPR_WILD,
			NewInteger(ListLength(clmem)),NewInteger(ListLength(clmem))),0));
	SetAtomProperty(clnm,PROP_INDEX,t);
	SetAtomProperty(clnm,OPR_CLASS,NewInteger(cls_no));
	if(!ListMember(cls_lens,NewInteger(ListLength(clmem))))
		cls_lens=AppendLast(cls_lens,NewInteger(ListLength(clmem)));
	SetAtomProperty(clnm,PROP_TYPE,MakeCompound1(OPR_CLASS,NewInteger(cls_no)));
	
	for(l1=clmem,i=1;l1;l1=ListTail(l1),i++)
	{
		int csp=-1;
		Term p2=GetAtomProperty(ListFirst(l1),PROP_TYPE);
		if(CompoundName(p2)==OPR_FIELD)
		{
			if(IntegerValue(CompoundArg2(p2))<4)
				csp=0;
			else csp=1;
			p2=GetAtomProperty(CompoundArg1(p2),PROP_TYPE);
		}
		else
			csp=(int)IntegerValue(CompoundArgN(p2,4));
		if(CompoundArgN(p2,5))
			ms++;
		else
			noms++;
		if(spin==-1)
			spin=csp;
		else if(spin!=csp)
		{
			ErrorInfo(144);
			printf("class %s: different spin particles.\n",
					AtomValue(clnm));
		}
		SetAtomProperty(ListFirst(l1),OPR_CLASS,MakeCompound2(OPR_CLASS,
				NewInteger(cls_no),NewInteger(i)));
	}
	
	t=GetAtomProperty(ListFirst(clmem),A_ANTI2);
	if(t) t=GetAtomProperty(t,A_ANTI);
	if(t==0) t=GetAtomProperty(ListFirst(clmem),A_ANTI);
	for(i=1;i<cls_no;i++)
		if(ListFirst(cls_membs[i])==t)
		{
			List l2;
			if(ListLength(clmem)!=ListLength(cls_membs[i]))
			{
				ErrorInfo(144);
				printf("classes %s/%s: different lengths.\n",
						AtomValue(cls_names[i]),AtomValue(clnm));
				return 0;
			}
			for(l1=clmem,l2=cls_membs[i];l1;l1=ListTail(l1),l2=ListTail(l2))
			{
				Term ta2;
				ta2=GetAtomProperty(ListFirst(l1),A_ANTI2);
				if(ta2) ta2=GetAtomProperty(ta2,A_ANTI);
				if(ta2==0)
					ta2=GetAtomProperty(ListFirst(l1),A_ANTI);
				if(ta2!=ListFirst(l2))
				{
					ErrorInfo(900);
					printf("classes %s/%s: '%s'/'%s' are not conjugated.\n",
							AtomValue(cls_names[i]),AtomValue(clnm),
							AtomValue(ListFirst(l1)),AtomValue(ListFirst(l2)));
					return 0;
				}
			}
			anti=cls_names[i];
			break;
		}
	
	if(ms&&noms)
	{
		WarningInfo(0);
		printf("class %s: mixing massive and massless particles.\n",
				AtomValue(clnm));
	}
	
	if(anti)
	{
		t=GetAtomProperty(anti,PROP_TYPE);
		SetCompoundArg(t,2,clnm);
		SetAtomProperty(clnm,PROP_TYPE,t);
	}
	else
	{	
		t=MakeCompound(OPR_PARTICLE,8);
		SetCompoundArg(t,1,clnm);
		SetCompoundArg(t,2,clnm);
		sWriteTerm(cbuf,clmem);
		SetCompoundArg(t,3,NewAtom(cbuf,0));
		SetCompoundArg(t,4,NewInteger(spin));		
		if(ms)
		{
			sprintf(cbuf,"%sMass",AtomValue(clnm));
			SetCompoundArg(t,5,NewAtom(cbuf,0));
		}
		SetAtomProperty(clnm,PROP_TYPE,t);
	}
	
	t=GetAtomProperty(ListFirst(clmem),A_COLOR);
	if(t)
		SetAtomProperty(clnm,A_COLOR,t);
	
	cls_names[cls_no]=clnm;
	cls_membs[cls_no]=clmem;
	cls_no++;
	if(cls_no==128)
	{
		puts("Internal error - too many classes");
		exit(0);
	}
	if(ctype==2)
	{
		List amem=0;
		Atom anm;
		for(l1=clmem;l1;l1=ListTail(l1))
			amem=AppendLast(amem,GetAtomProperty(ListFirst(l1),A_ANTI));
		sprintf(cbuf,"%s",AtomValue(clnm));
		if(islower(cbuf[0]))
			cbuf[0]=toupper(cbuf[0]);
		else
			cbuf[0]=tolower(cbuf[0]);
		anm=NewAtom(cbuf,0);
		anm=MakeCompound1(A_I,MakeCompound2(OPR_EQSIGN,anm,amem));
		ProcClass(anm,0);
	}
	
	return 0;	
}


void cls_set_em(void)
{
	int i;
	for(i=1;i<cls_no;i++)
	{
		Term ch;
		List l;
		ch=GetAtomProperty(ListFirst(cls_membs[i]),A_EM_CHARGE);
		for(l=ListTail(cls_membs[i]);l;l=ListTail(l))
			if(!EqualTerms(GetAtomProperty(ListFirst(cls_membs[i]),
						A_EM_CHARGE),ch))
		{
			printf("Error: SetEM: particles in class %s have different charge\n"
					, AtomValue(cls_names[i]));
		}
		SetAtomProperty(cls_names[i],A_EM_CHARGE,CopyTerm(ch));
		/*WriteTerm(cls_names[i]);
		printf(" -> ");
		WriteTerm(ch);
		puts("");*/
	}
}


void cls_wrt_nms(FILE *f, List cpml)
{
	List l1;
	int i;
	for(i=1;i<cls_no;i++)
	{
		Term fp=GetAtomProperty(cls_names[i],A_FANUM);
		int cl, no, j, m=0;
		if(fp==0)
			continue;
		cl=(int)IntegerValue(ListFirst(fp));
		no=(int)IntegerValue(ListFirst(ListTail(fp)));
		if(no<0) {m=1; no=-no;}
		for(l1=cls_membs[i],j=1;l1;l1=ListTail(l1),j++)
		{
			fprintf(f,"prt[\"");fWriteTerm(f,ListFirst(l1));fprintf(f,"\"] = ");
			if(m) fprintf(f,"-");
			switch(cl){
				case 0: fprintf(f,"S[%d, {%d}]\n",no,j); break;
				case 1: fprintf(f,"F[%d, {%d}]\n",no,j); break;
				case 2: fprintf(f,"V[%d, {%d}]\n",no,j); break;
					  }
		}
	}
	fprintf(f,"\n");
	
	for(i=1;i<cls_no;i++)
	{
		Term fp=GetAtomProperty(cls_names[i],PROP_TYPE);
		int j;
		List mmatr=0;
		
		Atom ms=CompoundArgN(fp,5);
		if(ms==0)
			continue;
		
		if(CompoundArg1(fp)!=CompoundArg2(fp) && CompoundArg2(fp)==cls_names[i])
			continue;
		
		for(j=1,l1=cls_membs[i];l1;l1=ListTail(l1),j++)
		{
			Term prp=GetAtomProperty(ListFirst(l1),PROP_TYPE);
			if(CompoundName(prp)==OPR_FIELD)
				prp=GetAtomProperty(CompoundArg1(prp),PROP_TYPE);
			fprintf(f,"%s[%d] = %s;\n",AtomValue(ms),j,CompoundArgN(prp,5)?
					AtomValue(CompoundArgN(prp,5)):"0");
			mmatr=AppendLast(mmatr,MakeCompound2(OPR_COLON,AppendFirst(0,NewInteger(j)),
					AppendFirst(0,MakeCompound3(A_MTERM,NewInteger(1),
					AppendFirst(0,MakeCompound2(OPR_POW,CompoundArgN(prp,5)?
					CompoundArgN(prp,5):NewAtom("0",0),NewInteger(1))),0))));
		}
		
		mtr_known=AppendLast(mtr_known,MakeCompound2(OPR_EQSIGN,ms,mmatr));
		SetAtomProperty(ms,OPR_CLASS,AppendFirst(0,NewInteger(ListLength(cls_membs[i]))));
		if(mtr_visi)
		  mtr_visi=AppendLast(mtr_visi,ms);
		for(l1=cpml;l1;l1=ListTail(l1))
		{
			if(ListFirst(l1)==ms)
			{
				fprintf(f,"%s[gen_, _] := %s[gen];\n",
						AtomValue(ms),AtomValue(ms));
				ChangeList(l1,0);
			}
		}
	}
		
	for(l1=cpml;l1;l1=ListTail(l1))
		{
			if(ListFirst(l1)==0)
				continue;
			fprintf(f,"%s[_] = %s\n",AtomValue(ListFirst(l1)),
					AtomValue(ListFirst(l1)));
		}
	fprintf(f,"\n");
}

/*
void cls_wrt_ind(FILE *f)
{
	int l;
	for(l=cls_lens;l;l=ListTail(l))
		fprintf(f,"IndexRange[ Index[TheyAre%ld] ] = Range[%ld]\n",
				IntegerValue(ListFirst(l)),IntegerValue(ListFirst(l)));
}
*/

void cls_wrt_ind(FILE *f)
{
	int i;
	for(i=1;i<cls_no;i++)
		{
		Term prop=GetAtomProperty(cls_names[i],PROP_TYPE);
		if(cls_names[i]!=CompoundArg1(prop))
			continue;
		fprintf(f,"IndexRange[ Index[%s] ] = Range[%d]\n",
			AtomValue(cls_names[i]),ListLength(cls_membs[i]));
		}
}

int cls_prt_info(Term *prp, Atom *a)
{
	Term cl=GetAtomProperty(*a,OPR_CLASS);
	int  n;
	if(cl==0)
		return 0;
	n=(int)IntegerValue(CompoundArg1(cl));
	if(cls_ini[n])
	{
		*prp=0;
		return 0;
	}
	cls_ini[n]=1;
	*prp=GetAtomProperty(cls_names[n],PROP_TYPE);
	*a=cls_names[n];
	return ListLength(cls_membs[n]);
}

extern int sp_cmp(Term,Term);


static void setzprm(List ml)
{
	for(;ml;ml=ListTail(ml))
	{
		List l1;
		List l2=ConsumeCompoundArg(ListFirst(ml),2);
		List l3=ConsumeCompoundArg(ListFirst(ml),3);
		l3=SortedList(l3,sp_cmp);
		for(l1=l2;l1;l1=ListTail(l1))
			if(GetAtomProperty(CompoundArg1(ListFirst(l1)),A_INFINITESIMAL) && 
                            !GetAtomProperty(CompoundArg1(ListFirst(l1)),PROP_TYPE)
                        )
			{
				if(oprcmatr)
					l3=AppendFirst(l3,MakeCompound2(OPR_CLASS,OPR_CLASS,0));
				else
				{
					l3=AppendFirst(l3,ListFirst(l1));
					ChangeList(l1,0);
					l2=CutFromList(l2,l1);
				}
				break;
			}
		if(opclsreal)
		for(l1=l2;l1;l1=ListTail(l1))
			if(CompoundArg1(ListFirst(l1))==A_I)
			{
				l3=AppendFirst(l3,ListFirst(l1));
				ChangeList(l1,0);
				l2=CutFromList(l2,l1);
				break;
			}
		SetCompoundArg(ListFirst(ml),2,l2);
		SetCompoundArg(ListFirst(ml),3,l3);
	}
}

static void unsetzprm(List ml)
{
	for(;ml;ml=ListTail(ml))
	{
		List l1;
		List l2=ConsumeCompoundArg(ListFirst(ml),2);
		List l3=ConsumeCompoundArg(ListFirst(ml),3);
		for(l1=l3;l1;l1=ListTail(l1))
			if(CompoundName(ListFirst(l1))==OPR_CLASS)
			{
				l3=CutFromList(l3,l1);
				break;
			}
		for(l1=l3;l1;l1=ListTail(l1))
			if(CompoundArg1(ListFirst(l1))==A_I)
			{
				l2=AppendFirst(l2,ListFirst(l1));
				ChangeList(l1,0);
				l3=CutFromList(l3,l1);
				break;
			}
		for(l1=l3;l1;l1=ListTail(l1))
			if(CompoundName(ListFirst(l1))==OPR_POW)
			{
				l2=AppendFirst(l2,ListFirst(l1));
				ChangeList(l1,0);
				l3=CutFromList(l3,l1);
				break;
			}
		SetCompoundArg(ListFirst(ml),2,l2);
		SetCompoundArg(ListFirst(ml),3,l3);
	}
}

static void add_elem(List *el, List prt, List ml)
{
	List l1;
	List pl=0;
	
	for(l1=prt;l1;l1=ListTail(l1))
	{
		Term prp=GetAtomProperty(CompoundArg1(ListFirst(l1)),OPR_CLASS);
		if(prp)
			pl=AppendLast(pl,CompoundArg2(prp));
	}
	
	for(l1=ml;l1;l1=ListTail(l1))
	{
		Term m2=ListFirst(l1);
		Term cme;
		List lt=ConsumeCompoundArg(m2,3);
		List l2;
		for(l2=(*el);l2;l2=ListTail(l2))
			if(EqualTerms(lt,CompoundArg1(ListFirst(l2))))
				break;
		if(l2==0)
		{
			cme=MakeCompound2(OPR_MINUS,lt,0);
			*el=AppendLast(*el,cme);
		}
		else
		{
			cme=ListFirst(l2);
			FreeAtomic(lt);
		}
		lt=ConsumeCompoundArg(cme,2);
		for(l2=lt;l2;l2=ListTail(l2))
			if(EqualTerms(pl,CompoundArg1(ListFirst(l2))))
				break;
		if(l2==0)
			lt=AppendLast(lt,MakeCompound2(OPR_COLON,CopyTerm(pl),MakeList1(m2)));
		else
			AppendLast(CompoundArg2(ListFirst(l2)),m2);
		SetCompoundArg(cme,2,lt);
	}
	
	FreeAtomic(pl);
	RemoveList(ml);
}

static int mtrcmp(Term t1, Term t2)
{
	List l1=CompoundArg1(t1),l2=CompoundArg1(t2);
	for(;l1;l1=ListTail(l1),l2=ListTail(l2))
	{
		int c=(int)IntegerValue(ListFirst(l1))-(int)IntegerValue(ListFirst(l2));
		if(c) return c;
	}
	puts("internal error: mtrcmp(cls)");
	return 0;
}

static int icmp(Term t1, Term t2)
{
	return (int)IntegerValue(t1)-(int)IntegerValue(t2);
}

static int spec_symm=0;

static Term mka(List bl, List el, List cl, List il)
{
	Term m;
	int p=0;
	char cbuf[32];
	Atom a;
	List l;
	p=sprintf(cbuf,"(");
	for(l=bl;l;l=ListTail(l))
		p+=sprintf(cbuf+p,"%ld,",IntegerValue(ListFirst(l)));
	for(l=il;l;l=ListTail(l))
		p+=sprintf(cbuf+p,(el||ListTail(l))?"%ld,":"%ld)",IntegerValue(ListFirst(l)));
	for(l=el;l;l=ListTail(l))
		p+=sprintf(cbuf+p,ListTail(l)?"%ld,":"%ld)",IntegerValue(ListFirst(l)));
	if(spec_symm==1)
		p+=sprintf(cbuf+p,")");
	l=CopyTerm(bl);
	l=ConcatList(l,CopyTerm(cl));
	l=ConcatList(l,CopyTerm(el));
	a=NewAtom(cbuf,0);
	if(spec_symm) 
		SetAtomProperty(a,A_I,NewInteger(spec_symm));
	m=MakeCompound3(A_MTERM,NewInteger(1),MakeList1(
		MakeCompound2(OPR_POW,a,NewInteger(1))),0);
	return MakeCompound2(OPR_COLON,l,MakeList1(m));
}

static List permute0(List l)
{
	List l1,ret=0;
	int i, n;
	if(ListLength(l)<2)
		return MakeList1(CopyTerm(l));
	n=ListLength(l);
	for(i=1;i<=n;i++)
		{
		List m=CopyTerm(l),m1;
		Term t=ListNth(l,i);
		m=CutFromList(m,ListNthList(m,i));
		m1=permute0(m);
		FreeAtomic(m);
		for(l1=m1;l1;l1=ListTail(l1))
			ChangeList(l1,AppendFirst(ListFirst(l1),t));
		ret=ConcatList(ret,m1);
		}
	return ret;
}
	
static int ListMember1(List l, Term m)
	{
	for(;l;l=ListTail(l))
		if(EqualTerms(ListFirst(l),m))
			return 1;
	return 0;
	}
	
static List permute(List l)
{
	List l1, l2;
	List ret=0;

	l1=permute0(l);
	for(l2=l1;l2;l2=ListTail(l2))
		{
		List l3=ListFirst(l2);
		if(!EqualTerms(l,l3) && !ListMember1(ret,l3))
			{
			ret=AppendLast(ret,l3);
			ChangeList(l2,0);
			}
		}
	FreeAtomic(l1);
	return ret;
}


static void symm_mtr(List *mtr, List cls)
{
	List l1,l2,l3;
	int  i,i1,n,i2;
	
	if((n=ListLength(cls))<2) return;
	if(spec_symm==-1) return;
	
	i2=0;
	for(i=1;i<=n-1;i++)
		if(ListNth(cls,i)==ListNth(cls,i+1)) i2++;
	if(i2==0)
		return;
	
/*	WriteTerm(cls);DumpList(*mtr);*/
	
	for(l1=cls,i1=1;ListTail(l1);l1=ListTail(l1),i1++)
	{
		if(ListFirst(l1)!=ListFirst(ListTail(l1)))
			continue;
		for(l2=ListTail(l1),i2=i1+1;l2&&ListFirst(l1)==ListFirst(l2);
			l2=ListTail(l2),i2++);
		i2--;
				
		for(l2=*mtr;l2;l2=ListTail(l2))
		{
			List il=0;
			l3=CompoundArg1(ListFirst(l2));
			for(i=i1;i<=i2;i++)
				il=AppendLast(il,ListNth(l3,i));
			if(spec_symm==2 && IntegerValue(ListFirst(il))>
					IntegerValue(ListNth(il,2)))
				{
				List l4;
				for(l4=CompoundArg2(ListFirst(l2));l4;l4=ListTail(l4))
					SetCompoundArg(ListFirst(l4),1,NewInteger(-IntegerValue(
							CompoundArg1(ListFirst(l4)))));
				}
			il=SortedList(il,icmp);
			for(i=1;l3;l3=ListTail(l3),i++)
				if(i>=i1 && i<=i2)
					ChangeList(l3,ListNth(il,i-i1+1));
			FreeAtomic(il);
		}
		for(i=i1;l1&&i<i2;i++) l1=ListTail(l1),i1++;
		if(l1==0 || ListTail(l1)==0) break;
	}
	
	for(l1=cls,i1=1;ListTail(l1);l1=ListTail(l1),i1++)
	{
		if(ListFirst(l1)!=ListFirst(ListTail(l1)))
			continue;
		for(l2=ListTail(l1),i2=i1+1;l2&&ListFirst(l1)==ListFirst(l2);
			l2=ListTail(l2),i2++);
		i2--;
		for(l2=*mtr;l2;l2=ListTail(l2))
		{
			List il=0,bl=0,el=0,ll;
			l3=CompoundArg1(ListFirst(l2));
			for(i=1;i<i1;i++)
				bl=AppendLast(bl,ListNth(l3,i));
			for(i=i1;i<=i2;i++)
				il=AppendLast(il,ListNth(l3,i));
			for(i=i2+1;i<=n;i++)
				el=AppendLast(el,ListNth(l3,i));
			ll=permute(il);
			for(l3=ll;l3;l3=ListTail(l3))
				*mtr=AppendFirst(*mtr,mka(bl,el,ListFirst(l3),il));
			FreeAtomic(ll);
			FreeAtomic(bl);
			FreeAtomic(il);
			FreeAtomic(el);
		}
		for(i=i1;l1&&i<i2;i++) l1=ListTail(l1),i1++;
		if(l1==0 || ListTail(l1)==0) break;
	}
	
	*mtr=SortedList(*mtr,mtrcmp);
/*	DumpList(*mtr);*/
	
	

}
					
static int mtr_infi=0;

static Atom mk_mtr(List mtr, Term *nf, List *sl, List mdim, List mcls)
{
	Term a2=MakeCompound(A_ALG2,5);
	List ml=0;
	List l1,l2;
	mtr=SortedList(mtr,mtrcmp);
	for(l1=mtr;l1;l1=ListTail(l1))
		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
			ml=AppendLast(ml,ListFirst(l2));
	SetCompoundArg(a2,5,ml);
	alg2_common_n(a2);
	alg2_common_s(a2);
	*nf=ConsumeCompoundArg(a2,2);
	*sl=ConsumeCompoundArg(a2,3);
	RemoveList(ConsumeCompoundArg(a2,5));
	if(ListLength(CompoundArg1(ListFirst(mtr)))==2)
	{
		int i, dlt=1, cnt=0;
		for(l1=mtr,i=1;l1;l1=ListTail(l1),i++)
		{
			Term vl=CompoundArg1(ListFirst(l1));
			cnt++;
			if(IntegerValue(ListFirst(vl))!=i || 
					ListFirst(vl)!=ListFirst(ListTail(vl)))
			{ dlt=0;break;}
			vl=CompoundArg2(ListFirst(l1));
			if(ListLength(vl)!=1) { dlt=0; break;}
			vl=ListFirst(vl);
			if(CompoundArg1(vl)!=NewInteger(1) || CompoundArg2(vl))
			{ dlt=0; break;}
		}
		if(IntegerValue(ListFirst(mdim))!=cnt ||
			ListFirst(mdim)!=ListFirst(ListTail(mdim)))
			dlt=0;
		if(dlt)
		{
			FreeAtomic(mtr);
			return A_DELTA;
		}
	}
	
	if(ListLength(CompoundArg1(ListFirst(mtr)))==1)
	{
		int i, dlt=1, cnt=0;
		for(l1=mtr,i=1;l1;l1=ListTail(l1),i++)
		{
			Term vl=CompoundArg1(ListFirst(l1));
			cnt++;
			if(IntegerValue(ListFirst(vl))!=i)
			{ dlt=0;break;}
			vl=CompoundArg2(ListFirst(l1));
			if(ListLength(vl)!=1) { dlt=0; break;}
			vl=ListFirst(vl);
			if(CompoundArg1(vl)!=NewInteger(1) || CompoundArg2(vl))
			{ dlt=0; break;}
		}
		if(IntegerValue(ListFirst(mdim))!=cnt)
			dlt=0;
		if(dlt)
		{
			FreeAtomic(mtr);
			return A_DELTA;
		}
	}

	FreeAtomic(a2);
	
	symm_mtr(&mtr, mcls);
	
/*	WriteTerm(mdim);WriteTerm(mcls);DumpList(mtr);*/
	
	for(l1=mtr_known;l1;l1=ListTail(l1))
		if(EqualTerms(CompoundArg2(ListFirst(l1)),mtr))
			{
			if(EqualTerms(
				GetAtomProperty(CompoundArg1(ListFirst(l1)),OPR_CLASS),mdim))
			break;
			}
	if( l1 )
	{
		FreeAtomic(mtr);
		mtr_fit++;
		return CompoundArg1(ListFirst(l1));
	}
	else
	{
		char cbuf[16];
		Atom nm;
		if(mtr_infi)
			{
			sprintf(cbuf,"dMTR%03d",++dmtr_no);
			nm=NewAtom(cbuf,0);
			SetAtomProperty(nm,A_INFINITESIMAL,NewInteger(1));
			}
		else
		{	
			sprintf(cbuf,"MTR%03d",++mtr_no);
			nm=NewAtom(cbuf,0);
		}
		SetAtomProperty(nm,OPR_CLASS,mdim);
		/*if(spec_symm)
			SetAtomProperty(nm,A_I,NewInteger(spec_symm));*/
		mtr_known=AppendLast(mtr_known,MakeCompound2(OPR_EQSIGN,nm,mtr));
		return nm;
	}
		
/*	WriteTerm(*nf);printf(" ");WriteTerm(*sl);puts("");
	DumpList(mtr);*/
	return A_DELTA;
}

extern int sp_cmp(Term, Term);

static void fixind(Term a2, int i, int j)
{
	List il1, il2, l;
	
	Atom p1=CompoundArg1(ListNth(CompoundArg1(a2),i));
	Atom p2=CompoundArg1(ListNth(CompoundArg1(a2),j));
	SetCompoundArg(ListNth(CompoundArg1(a2),i),1,p2);
	SetCompoundArg(ListNth(CompoundArg1(a2),j),1,p1);
	
	il1=CompoundArg1(CompoundArg2(ListNth(CompoundArg1(a2),i)));
	il2=CompoundArg1(CompoundArg2(ListNth(CompoundArg1(a2),j)));
	for(l=CompoundArgN(a2,5);l;l=ListTail(l))
		{
		List l1,l2,tl=ConsumeCompoundArg(ListFirst(l),3);
		for(l1=tl;l1;l1=ListTail(l1))
			{
			Term t=ListFirst(l1);
			List il;
			if(CompoundName(t)==A_MOMENT)
				{
				if(CompoundArg1(t)==NewInteger(i))
					SetCompoundArg(t,1,NewInteger(j));
				else if(CompoundArg1(t)==NewInteger(j))
					SetCompoundArg(t,1,NewInteger(i));
				}
			for(il=CompoundArg2(t);il;il=ListTail(il))
				{
				Term ii=ListFirst(il);
				int lm1=ListMember(il1,ii), lm2=ListMember(il2,ii);
				if(lm1)
					ChangeList(il,ListNth(il2,lm1));
				if(lm2)
					ChangeList(il,ListNth(il1,lm2));
				}
			}
		tl=SortedList(tl,sp_cmp);
		SetCompoundArg(ListFirst(l),3,tl);
		}
}
			
static void fixindf(Term a2, int i, int j)
{
	List il1, il2, l;
	
	Atom p1=CompoundArg1(ListNth(CompoundArg1(a2),i));
	Atom p2=CompoundArg1(ListNth(CompoundArg1(a2),j));
	SetCompoundArg(ListNth(CompoundArg1(a2),i),1,p2);
	SetCompoundArg(ListNth(CompoundArg1(a2),j),1,p1);
	
	il1=ListTail(CompoundArg1(CompoundArg2(ListNth(CompoundArg1(a2),i))));
	il2=ListTail(CompoundArg1(CompoundArg2(ListNth(CompoundArg1(a2),j))));

	for(l=CompoundArgN(a2,5);l;l=ListTail(l))
		{
		List l1,l2,tl; 
		/*WriteTerm(ListFirst(l));*/
		Term ga=0, g5=0;
		tl=ConsumeCompoundArg(ListFirst(l),3);
		for(l1=tl;l1;l1=ListTail(l1))
			{
			Term t=ListFirst(l1);
			List il;
			Atom s=CompoundArg1(t);
			if(CompoundName(t)==A_MOMENT)
				{
				if(CompoundArg1(t)==NewInteger(i))
					SetCompoundArg(t,1,NewInteger(j));
				else if(CompoundArg1(t)==NewInteger(j))
					SetCompoundArg(t,1,NewInteger(i));
				}
			if(s==A_GAMMA)
				ga=t;
			if(s==A_GAMMAP || s==A_GAMMAM || s==A_GAMMA5)
				g5=t;
			if(il1 || il2)
			for(il=CompoundArg2(t);il;il=ListTail(il))
				{
				Term ii=ListFirst(il);
				int lm1=ListMember(il1,ii), lm2=ListMember(il2,ii);
				if(lm1)
					ChangeList(il,ListNth(il2,lm1));
				if(lm2)
					ChangeList(il,ListNth(il1,lm2));
				}
			}
		if(ga)
			{
			Atom g5t=g5?CompoundArg1(g5):0;
			if(g5==0 || g5t!=A_GAMMA5)
				SetCompoundArg(CompoundArg1(ListFirst(l)),1,NewInteger(
					-IntegerValue(CompoundArg1(CompoundArg1(ListFirst(l))))));
			if(g5t==A_GAMMAP) 
				SetCompoundArg(g5,1,A_GAMMAM);
			if(g5t==A_GAMMAM)
				SetCompoundArg(g5,1,A_GAMMAP);
			}
		if(il1 || il2)
		tl=SortedList(tl,sp_cmp);
		SetCompoundArg(ListFirst(l),3,tl);
		/*printf("  ->  ");WriteTerm(ListFirst(l));puts("");*/
		}
}
		
void setcls_1(Term a2, List a2l)
{
	int i,j,ci,cf=0;
	int pno, is_cls[6], is_col[6];
	int ind_no[6]={1000,1000,1000,1000,1000,1000};
	int sym2=0, slp=0;
	Atom prts[6];
	List matr_dim=0;
	List matr_cls=0;
	List l1,l2,l3,elli=0,ml=0;
	
	if(CompoundArgN(a2,5)==0 || CompoundArg1(a2)==0 ||
			ListLength(CompoundArg1(a2))==1)
		return;
		
	/*WriteTerm(CompoundArg1(a2));puts("");*/
	pno=ListLength(CompoundArg1(a2));
	for(l1=CompoundArg1(a2),i=0;l1;l1=ListTail(l1),i++)
	{
		Atom p=CompoundArg1(ListFirst(l1));
		Term prp=GetAtomProperty(p,OPR_CLASS);
		if(prp)
		{
			is_cls[i]=1;
			prts[i]=cls_names[IntegerValue(CompoundArg1(prp))];
			matr_dim=AppendLast(matr_dim,NewInteger(
				ListLength(cls_membs[IntegerValue(CompoundArg1(prp))])));
			matr_cls=AppendLast(matr_cls,prts[i]);
			is_col[i]=(GetAtomProperty(prts[i],A_COLOR)!=0);
			cf++;
		}
		else
		{
			is_cls[i]=0;
			is_col[i]=0;
			prts[i]=p;
		}
	}
	
	if(cf==0)
		return;
	/*{WriteTerm(matr_cls);puts("");}*/

	/*WriteVertex(CompoundArg1(a2));WriteTerm(matr_dim);puts("");*/

	
	
	alg2_decommon_s(a2);
	alg2_decommon_n(a2);
	
	
	{
	List pl=CompoundArg1(a2);
	if(cf==4  &&
		   CompoundName(CompoundArg2(ListNth(pl,1)))==OPR_SCALAR &&
		   CompoundName(CompoundArg2(ListNth(pl,2)))==OPR_SCALAR &&
		   CompoundName(CompoundArg2(ListNth(pl,3)))==OPR_SCALAR &&
		   CompoundName(CompoundArg2(ListNth(pl,4)))==OPR_SCALAR)
		{
		if(is_col[0] && is_col[1] && is_col[2] && is_col[3])
			sym2=1;
		else if(!is_col[0] && !is_col[1] && !is_col[2] && !is_col[3] &&
		   prts[0]!=prts[2] && prts[0]!=prts[3] && prts[1]!=prts[2] && prts[1]!=prts[3])
		   { sym2=1; slp=1;}
		}
		
		if(sym2)
		{
		sym2=1;
		if(is_cls[0] && is_cls[1] && prts[0]==prts[1] &&
			CompoundArg1(ListNth(pl,1))!=CompoundArg1(ListNth(pl,2)))
			{
			Term a2c=CopyTerm(a2);
			fixind(a2c,1,2);
			
			if(is_cls[2] && is_cls[3] && prts[2]==prts[3] &&
			CompoundArg1(ListNth(pl,3))!=CompoundArg1(ListNth(pl,4)))
				{
				Term a2cc=CopyTerm(a2c);
				fixind(a2cc,3,4);
				setzprm(CompoundArgN(a2cc,5));
				add_elem(&elli,CompoundArg1(a2cc),ConsumeCompoundArg(a2cc,5));
				FreeAtomic(a2cc);
				}
				
			setzprm(CompoundArgN(a2c,5));
			add_elem(&elli,CompoundArg1(a2c),ConsumeCompoundArg(a2c,5));
			FreeAtomic(a2c);
			}
			
		if(is_cls[2] && is_cls[3] && prts[2]==prts[3] &&
			CompoundArg1(ListNth(pl,3))!=CompoundArg1(ListNth(pl,4)))
			{
			Term a2c=CopyTerm(a2);
			fixind(a2c,3,4);
			setzprm(CompoundArgN(a2c,5));
			add_elem(&elli,CompoundArg1(a2c),ConsumeCompoundArg(a2c,5));
			FreeAtomic(a2c);
			}
		}
	}
	
	{
	List pl=CompoundArg1(a2);
	if( (cf==2 || cf==3) && (pno==2 || pno==3) &&
		   CompoundName(CompoundArg2(ListNth(pl,1)))==OPR_SPINOR &&
		   CompoundName(CompoundArg2(ListNth(pl,2)))==OPR_SPINOR &&
		   is_cls[0] && is_cls[1] && prts[0]==prts[1] )
		{
		sym2=2;
		if( CompoundArg1(ListNth(pl,1))!=CompoundArg1(ListNth(pl,2)))
			{
			Term a2c=CopyTerm(a2);
			fixindf(a2c,1,2);				
			setzprm(CompoundArgN(a2c,5));
			add_elem(&elli,CompoundArg1(a2c),ConsumeCompoundArg(a2c,5));
			FreeAtomic(a2c);
			}			
		}
	}
			
			
/*	WriteTerm(a2);puts("");*/
	setzprm(CompoundArgN(a2,5));
	add_elem(&elli,CompoundArg1(a2),ConsumeCompoundArg(a2,5));
	
	for(l2=a2l;l2;l2=ListTail(l2))
	{
		Term a2c=ListFirst(l2);
		if(ListLength(CompoundArg1(a2c))!=pno)
			continue;
		for(l1=CompoundArg1(a2c),i=0;l1;l1=ListTail(l1),i++)
		{
			Atom p=CompoundArg1(ListFirst(l1));
			Term prp=GetAtomProperty(p,OPR_CLASS);
			if(prp)
				p=cls_names[IntegerValue(CompoundArg1(prp))];
			if(p!=prts[i])
				break;
		}
		if(l1)
			continue;

		alg2_decommon_s(a2c);
		alg2_decommon_n(a2c);

		{
		List pl=CompoundArg1(a2c);
		if(sym2==1)
		{   
		if(is_cls[0] && is_cls[1] && prts[0]==prts[1] &&
			CompoundArg1(ListNth(pl,1))!=CompoundArg1(ListNth(pl,2)))
			{
			Term a2cc=CopyTerm(a2c);
			fixind(a2cc,1,2);
			
			if(is_cls[2] && is_cls[3] && prts[2]==prts[3] &&
			CompoundArg1(ListNth(pl,3))!=CompoundArg1(ListNth(pl,4)))
				{
				Term a2ccc=CopyTerm(a2cc);
				fixind(a2ccc,3,4);
				setzprm(CompoundArgN(a2ccc,5));
				add_elem(&elli,CompoundArg1(a2ccc),ConsumeCompoundArg(a2ccc,5));
				FreeAtomic(a2ccc);
				}
				
			setzprm(CompoundArgN(a2cc,5));
			add_elem(&elli,CompoundArg1(a2cc),ConsumeCompoundArg(a2cc,5));
			FreeAtomic(a2cc);
			}
			
		if(is_cls[2] && is_cls[3] && prts[2]==prts[3] &&
			CompoundArg1(ListNth(pl,3))!=CompoundArg1(ListNth(pl,4)))
			{
			Term a2cc=CopyTerm(a2c);
			fixind(a2cc,3,4);
			setzprm(CompoundArgN(a2cc,5));
			add_elem(&elli,CompoundArg1(a2cc),ConsumeCompoundArg(a2cc,5));
			FreeAtomic(a2cc);
			}
		}


		}

		{
		List pl=CompoundArg1(a2c);
		if(sym2==2 )
			{
			if(CompoundArg1(ListNth(pl,1))!=CompoundArg1(ListNth(pl,2)))
				{
				Term a2cc=CopyTerm(a2c);
				fixindf(a2cc,1,2);				
				setzprm(CompoundArgN(a2cc,5));
				add_elem(&elli,CompoundArg1(a2cc),ConsumeCompoundArg(a2cc,5));
				FreeAtomic(a2cc);
				}			
			}
		}
		
		/*WriteVertex(CompoundArg1(a2c));puts("");*/
		setzprm(CompoundArgN(a2c,5));
		add_elem(&elli,CompoundArg1(a2c),ConsumeCompoundArg(a2c,5));
	}


	cf=0;
	ci=1;
	for(l1=CompoundArg1(a2),i=0;l1;l1=ListTail(l1),i++)
	{
		if(is_cls[i])
			SetCompoundArg(ListFirst(l1),1,prts[i]);
		for(l2=CompoundArg1(CompoundArg2(ListFirst(l1)));l2;l2=ListTail(l2))
		{
			ChangeList(l2,NewInteger(ci));
			ci++;
		}
		if(is_cls[i])
		{
			List ll=ConsumeCompoundArg(CompoundArg2(ListFirst(l1)),1);
			ll=AppendLast(ll,NewInteger(ci));
			SetCompoundArg(CompoundArg2(ListFirst(l1)),1,ll);
			ind_no[cf]=ci;
			ci++; cf++;
		}
	}
	
	for(l1=elli;l1;l1=ListTail(l1))
	for(l2=CompoundArg1(ListFirst(l1));l2;l2=ListTail(l2))
		if(CompoundName(ListFirst(l2))==OPR_SPECIAL || 
					CompoundName(ListFirst(l2))==A_MOMENT)
		{
			List l3;
			for(l3=CompoundArg2(ListFirst(l2));l3;l3=ListTail(l3))
			{
				int ov=(int)IntegerValue(ListFirst(l3));
				int dn=0;
				for(i=0;i<cf;i++)
					if(ov>=ind_no[i]-i) dn++;
				if(dn)
					ChangeList(l3,NewInteger(ov+dn));
			}
		}
	
	
	spec_symm=0;
	if(ListLength(CompoundArg1(a2))==3)
	{
		List pl=CompoundArg1(a2);
		
		if(CompoundName(CompoundArg2(ListNth(pl,1)))==OPR_SPINOR &&
		   CompoundName(CompoundArg2(ListNth(pl,2)))==OPR_SPINOR &&
		   CompoundName(CompoundArg2(ListNth(pl,3)))==OPR_VECTOR &&
		   CompoundArg1(ListNth(pl,1))==CompoundArg1(ListNth(pl,2)))
		   	spec_symm=1;
	}


	if(slp)
	for(l1=elli;l1;l1=ListTail(l1))
	{
		int c13=0, c14=0, fail=0;
		List ll, l13=0, l14=0;
		Term nn;
		
		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
		{
			Term k1,k2,k3,k4;
			l3=CompoundArg1(ListFirst(l2));
			k1=ListFirst(l3); l3=ListTail(l3);
			k2=ListFirst(l3); l3=ListTail(l3);
			k3=ListFirst(l3); l3=ListTail(l3);
			k4=ListFirst(l3);
			if(k1==k2 && k1==k3 && k1==k4)
				continue;
			if(k1==k3 && k2==k4) {c13++; continue;}
			if(k1==k4 && k2==k3) {c14++; continue;}
			fail++; break;
		}
		
		if(fail) continue;
		if(!(c13 && c14)) continue;
		
		/*WriteTerm(CompoundArg1(ListFirst(l1)));puts("");*/
		
		ll=ConsumeCompoundArg(ListFirst(l1),2);
		nn=CopyTerm(ListFirst(l1));
		/*DumpList(ll);*/
		for(l2=ll;l2;l2=ListTail(l2))
		{
			Term k1,k2,k3,k4;
			l3=CompoundArg1(ListFirst(l2));
			k1=ListFirst(l3); l3=ListTail(l3);
			k2=ListFirst(l3); l3=ListTail(l3);
			k3=ListFirst(l3); l3=ListTail(l3);
			k4=ListFirst(l3);
			if(k1==k4 && k2==k3)
				l14=AppendLast(l14,ListFirst(l2));
			else
				l13=AppendLast(l13,ListFirst(l2));
		}
		/*DumpList(l13);
		DumpList(l14);*/
		RemoveList(ll);
		SetCompoundArg(ListFirst(l1),2,l13);
		SetCompoundArg(nn,2,l14);
		AppendLast(elli,nn);
	}
			
	l2=0;
	for(i=0;i<cf;i++)
		l2=AppendLast(l2,NewInteger(ind_no[i]));

	for(l1=elli;l1;l1=ListTail(l1))
	{
		Term m2=MakeCompound(A_MTERM,3);
		Term nc, mtr;
		List tl,sl;
		List l2s, mcs, mds;
		int i1=0, j1=0, i2=0, j2=0;
		int j;
		tl=ConsumeCompoundArg(ListFirst(l1),1);
		
		
		if(ListLength(CompoundArg1(a2))==3)
		{
			List pl=CompoundArg1(a2);

			if(CompoundName(CompoundArg2(ListNth(pl,1)))==OPR_SCALAR &&
			   CompoundName(CompoundArg2(ListNth(pl,2)))==OPR_SCALAR &&
			   CompoundName(CompoundArg2(ListNth(pl,3)))==OPR_VECTOR &&
			   CompoundArg1(ListNth(pl,1))==CompoundArg1(ListNth(pl,2)))
		   		{
					spec_symm=0;
					List tl1=tl;
					if(ListLength(tl1)==2 && (CompoundName(ListFirst(tl1))==
							OPR_POW ||CompoundName(ListFirst(tl1))==OPR_CLASS ))
							tl1=ListTail(tl1);
					if(ListLength(tl1)!=1 || 
							CompoundName(ListFirst(tl1))!=A_MOMENT)
						{puts("Internal error (symmtr01)");
						/*WriteTerm(a2);puts("");WriteTerm(tl);puts("");*/}
						else if(CompoundArg1(ListFirst(tl1))!=NewInteger(3))
							spec_symm=2;
				}
		}
		
		if(sym2)
			spec_symm=-1;
			 
		/*{WriteTerm(l2);DumpList(CompoundArg2(ListFirst(l1)));}*/
		

		l2s=CopyTerm(l2);
		mcs=CopyTerm(matr_cls);
		mds=CopyTerm(matr_dim);

	if(opclsred)
	{
		
		for(i=1;i<ListLength(matr_cls);i++)
		{
			
			if( !sym2 && ((i>1 && ListNth(matr_cls,i-1)==ListNth(matr_cls,i)) ||
				ListNth(matr_cls,i)==ListNth(matr_cls,i+1)))
				continue;
			for(j=i+1;j<=ListLength(matr_cls);j++)
			{
				if(ListNth(matr_dim,i)!=ListNth(matr_dim,j))
					continue;
				if(!sym2 && (ListNth(matr_cls,j-1)==ListNth(matr_cls,j) ||
						(j<ListLength(matr_cls) && ListNth(matr_cls,j)==ListNth(matr_cls,j+1))))
					continue;
				for(l3=CompoundArg2(ListFirst(l1));l3;l3=ListTail(l3))
					if(ListNth(CompoundArg1(ListFirst(l3)),i)!=
							ListNth(CompoundArg1(ListFirst(l3)),j))
						break;
				if(l3)
					continue;
				if(i1==i||i1==j||j1==i||j1==j)
					continue;
				if(i1==0)
					{i1=i; j1=j;}
				else
					{i2=i;j2=j;break;}
			}
		}
	
		
		if(i2)
			{
			
			mtr=MakeCompound2(OPR_SPECIAL,A_DELTA,MakeList2(ListNth(l2,i2),ListNth(l2,j2)));
			tl=AppendLast(tl,mtr);
			for(l3=CompoundArg2(ListFirst(l1));l3;l3=ListTail(l3))
				{
				List il=ConsumeCompoundArg(ListFirst(l3),1);
				il=CutFromList(il,ListNthList(il,i2));
				SetCompoundArg(ListFirst(l3),1,il);
				}
			l2s=CutFromList(l2s,ListNthList(l2s,i2));
			mcs=CutFromList(mcs,ListNthList(mcs,i2));
			mds=CutFromList(mds,ListNthList(mds,i2));
			}
			
		if(i1)
			{
			mtr=MakeCompound2(OPR_SPECIAL,A_DELTA,MakeList2(ListNth(l2,i1),ListNth(l2,j1)));
			tl=AppendLast(tl,mtr);
			for(l3=CompoundArg2(ListFirst(l1));l3;l3=ListTail(l3))
				{
				List il=ConsumeCompoundArg(ListFirst(l3),1);
				il=CutFromList(il,ListNthList(il,i1));
				SetCompoundArg(ListFirst(l3),1,il);
				}
			l2s=CutFromList(l2s,ListNthList(l2s,i1));
			mcs=CutFromList(mcs,ListNthList(mcs,i1));
			mds=CutFromList(mds,ListNthList(mds,i1));
			}
	
		/*if(i1) {WriteTerm(l2s);DumpList(CompoundArg2(ListFirst(l1)));}*/
	
	}
		
		if(tl && CompoundName(ListFirst(tl))==OPR_CLASS)
			mtr_infi=1;
		else
			mtr_infi=0;
			
		mtr=mk_mtr(ConsumeCompoundArg(ListFirst(l1),2),&nc,&sl,
				CopyTerm(mds),mcs);
		
		if(mtr==A_DELTA && ListLength(l2s)==1)
		{
			if(!i1)
			{
				WarningInfo(0);puts("cls 1 dim delta");
				mtr=MakeCompound2(OPR_PARAMETER,mtr,MakeList2(ListFirst(l2s),ListFirst(l2s)));
			}
			else
				mtr=0;
		}
		else
			mtr=MakeCompound2(OPR_PARAMETER,mtr,CopyTerm(l2s));
			
		if(mtr) tl=AppendLast(tl,mtr);
		SetCompoundArg(m2,1,nc);
		SetCompoundArg(m2,2,sl);
		SetCompoundArg(m2,3,tl);
		ml=AppendLast(ml,m2);
		FreeAtomic(l2s);
		FreeAtomic(mds);
		FreeAtomic(mcs);
	}
	FreeAtomic(l2);
	unsetzprm(ml);
	SetCompoundArg(a2,5,ml);
	alg2_common_n(a2);
	alg2_common_s(a2);
/*	WriteTerm(CompoundArg1(a2));puts("");
	WriteTerm(CompoundArgN(a2,5));puts("");*/
	FreeAtomic(elli);
}

extern int opAbbrVrt;

void alg2_setcls(List a2l)
{
	List l1,l2;
	int i;
	
	if(cls_no==1)
		return;
		
	
	for(l1=a2l;l1;l1=ListTail(l1))
		for(l2=CompoundArg1(ListFirst(l1));l2;l2=ListTail(l2))
		{
			Term prp=GetAtomProperty(CompoundArg1(ListFirst(l2)),PROP_TYPE);
			if(CompoundName(prp)==OPR_FIELD && CompoundArg2(prp)==NewInteger(4))
				SetCompoundArg(ListFirst(l2),1,CompoundArg1(prp));
		}
	
	for(l1=a2l;l1;l1=ListTail(l1))
		setcls_1(ListFirst(l1),ListTail(l1));
		
	for(i=1;i<cls_no;i++)
		{
		Term prp=GetAtomProperty(cls_names[i],PROP_TYPE);
		Atom a1=CompoundArg1(prp), a2=CompoundArg2(prp);
		if(CompoundArgN(prp,4)==NewInteger(1))
			SetAtomProperty(cls_names[i],A_GRASS,NewInteger(1));
		if(a1==a2)
			{
			SetAtomProperty(a1,A_ANTI,a1);
			continue;
			}
		if(cls_names[i]==a2)
			continue;
		SetAtomProperty(a1,A_ANTI,a2);
		SetAtomProperty(a2,A_ANTI,a1);
		}
	
	if(mtr_visi && opAbbrVrt)
	{
	  for(l1=mtr_known;l1;l1=ListTail(l1))
	  {
	    if(!ListMember(mtr_visi,CompoundArg1(ListFirst(l1))))
	      continue;
	    for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
	    {
	      List l3;
	      for(l3=CompoundArg2(ListFirst(l2));l3;l3=ListTail(l3))
	      {
		List l4;
		for(l4=CompoundArg2(ListFirst(l3));l4;l4=ListTail(l4))
		{
		  Atom prm=CompoundArg1(ListFirst(l4));
		  if(is_atom(prm) && AtomValue(prm)[0]=='A')
		    SetAtomProperty(prm, OPR_LOCAL, A_I);
		}
	      }
	      
	    }
	  }
	}
	
/*	DumpList(mtr_known);*/

}


int cmplmtr(List l1)
{
	List ml,l2;
	int hasi=0;
		for(;l1;l1=ListTail(l1))
		for(ml=CompoundArg2(ListFirst(l1));ml;ml=ListTail(ml))
		for(l2=CompoundArg2(ListFirst(ml));l2;l2=ListTail(l2))
				if(CompoundArg1(ListFirst(l2))==A_I || 
					GetAtomProperty(CompoundArg1(ListFirst(l2)),A_ANTI))
				{
				hasi=1;
				break;
				}
	if(hasi) return 1; else return 0;
}

List cls_real_matr(void)
{
	List ret=0;
	List l;
	for(l=mtr_known;l;l=ListTail(l))
	{
		List l1=CompoundArg2(ListFirst(l));
		int hasi=cmplmtr(l1);
		if(GetAtomProperty(CompoundArg1(ListFirst(l)),A_INFINITESIMAL))
			continue;
		if(!hasi) ret=AppendLast(ret,CompoundArg1(ListFirst(l)));
	}
	return ret;
}

void cls_write_decl(FILE *f)
{
	int rmno=0, imno=0, i, p, first=1;
	int cmpno=0, decno=0;
	List lastr=0, lasti=0, lastm=0;
	List l;
	
	if(mtr_known==0)
		return;
	
	for(l=mtr_known;l;l=ListTail(l))
	{
		List l1=CompoundArg2(ListFirst(l));
		int hasi=cmplmtr(l1);
		if(GetAtomProperty(CompoundArg1(ListFirst(l)),A_INFINITESIMAL))
			continue;
		if(mtr_visi && !ListMember(mtr_visi, CompoundArg1(ListFirst(l))))
		      continue;
		if(hasi) {imno++; lasti=l;} else {rmno++; lastr=l;}
		lastm=l;
	}
	
	if(lastm==0) return;
	
	i=rmno;
	if(i) p=fprintf(f,"      double precision ");
	if(i)
	for(l=mtr_known;l;l=ListTail(l))
	{
		List l1=CompoundArg2(ListFirst(l));
		int hasi=cmplmtr(l1);
		if(GetAtomProperty(CompoundArg1(ListFirst(l)),A_INFINITESIMAL))
			continue;
		if(mtr_visi && !ListMember(mtr_visi, CompoundArg1(ListFirst(l))))
		      continue;
		if(hasi)
			continue;
		if(p>60)
		{
			fprintf(f,"\n");
			p=fprintf(f,"      double precision ");
			first=1;
		}
		if(!first)
			p+=fprintf(f,", ");
		else
			first=0;
		p+=fprintf(f,"%s",AtomValue(CompoundArg1(ListFirst(l))));
	}
	fprintf(f,"\n");
	
	first=1;
	i=imno;
	if(i) p=fprintf(f,"      double complex ");
	if(i)
	for(l=mtr_known;l;l=ListTail(l))
	{
		List l1=CompoundArg2(ListFirst(l));
		int hasi=cmplmtr(l1);
		if(GetAtomProperty(CompoundArg1(ListFirst(l)),A_INFINITESIMAL))
			continue;
		if(mtr_visi && !ListMember(mtr_visi, CompoundArg1(ListFirst(l))))
		      continue;
		if(!hasi)
			continue;
		if(p>60)
		{
			fprintf(f,"\n");
			p=fprintf(f,"      double complex ");
			first=1;
		}
		if(!first)
			p+=fprintf(f,", ");
		else
			first=0;
		p+=fprintf(f,"%s",AtomValue(CompoundArg1(ListFirst(l))));
	}
	fprintf(f,"\n\n");
	
	fprintf(f,"      dimension\n");
	p=fprintf(f,"     & ");
	for(l=mtr_known;l;l=ListTail(l))
	{
		List l1=GetAtomProperty(CompoundArg1(ListFirst(l)),OPR_CLASS);
		int comps=1;
		if(GetAtomProperty(CompoundArg1(ListFirst(l)),A_INFINITESIMAL))
			continue;
		if(mtr_visi && !ListMember(mtr_visi, CompoundArg1(ListFirst(l))))
		      continue;
		if(p>55)
		{
			fprintf(f,"\n");
			p=fprintf(f,"     & ");
			first=1;
		}
		p+=fprintf(f,"%s(",AtomValue(CompoundArg1(ListFirst(l))));
		for(;l1;l1=ListTail(l1))
		{
			p+=fprintf(f,"%ld%c",IntegerValue(ListFirst(l1)),
					ListTail(l1)?',':')');
			comps*=(int)IntegerValue(ListFirst(l1));
		}
		if(l!=lastm)
			p+=fprintf(f,", ");
		cmpno+=comps;
		decno+=ListLength(CompoundArg2(ListFirst(l)));
		if(comps==ListLength(CompoundArg2(ListFirst(l)))) 
			SetAtomProperty(CompoundArg1(ListFirst(l)),A_I,A_I);
	}
	fprintf(f,"\n\n");

	fprintf(f,"      common /mdl_mtrces/\n");
	p=fprintf(f,"     & ");
	for(l=mtr_known;l;l=ListTail(l))
	{
		if(GetAtomProperty(CompoundArg1(ListFirst(l)),A_INFINITESIMAL))
			continue;
		if(mtr_visi && !ListMember(mtr_visi, CompoundArg1(ListFirst(l))))
		      continue;
		if(p>55)
		{
			fprintf(f,"\n");
			p=fprintf(f,"     & ");
			first=1;
		}
		p+=fprintf(f,"%s",AtomValue(CompoundArg1(ListFirst(l))));
		if(l!=lastm)
			p+=fprintf(f,", ");
	}
	fprintf(f,"\n\n");
	
	printf("Matrices: %d real %d complex %4.1f%% zero components\n",
				rmno, imno, (1.0-(double)decno/cmpno)*100.0);



}

Term l_2_t(List l, int num, int den);
Atom fa_i, fa_cc=0;
extern void ccf_term(Term);
extern Atom toccf(Atom);

void cls_write_matr(FILE *f)
{
	
	List l1,l2,l3;
	int i,ldl;
	
	if(fa_i==0)
		fa_i=NewAtom("cI",0);
	if(fa_cc==0)
		fa_cc=NewAtom("Conjugate",0);
	
	for(l1=mtr_known;l1;l1=ListTail(l1))
	{
		Atom mtr=CompoundArg1(ListFirst(l1));
		List dl;
		if(GetAtomProperty(mtr,A_I))
			continue;
		if(GetAtomProperty(mtr,A_INFINITESIMAL))
			continue;
		
		dl=GetAtomProperty(mtr,OPR_CLASS);
		ldl=ListLength(dl);
		i=1;
		for(l2=dl;l2;l2=ListTail(l2))
			fprintf(f,"      do m%d = 1,%ld\n",i++,
					IntegerValue(ListFirst(l2)));
		for(l2=l1;l2;l2=ListTail(l2))
		{
			Atom mtr1=CompoundArg1(ListFirst(l2));
			if(GetAtomProperty(mtr1,A_INFINITESIMAL))
				continue;
			if(!EqualTerms(dl,GetAtomProperty(mtr1,OPR_CLASS)))
				continue;
			SetAtomProperty(mtr1,A_I,A_I);
			if(mtr_visi && !ListMember(mtr_visi, mtr1))
			  continue;
			
			fprintf(f,"      %s(",AtomValue(mtr1));
			for(i=1;i<=ldl;i++)
				fprintf(f,"m%d%c",i,i==ldl?')':',');
			fprintf(f,"=0D0\n");
		}
		for(i=0;i<ldl;i++)
			fprintf(f,"      enddo\n");
	}
	
	for(l1=mtr_known;l1;l1=ListTail(l1))
	if(GetAtomProperty(CompoundArg1(ListFirst(l1)),A_INFINITESIMAL)==0 &&
	  (mtr_visi==0 || ListMember(mtr_visi,CompoundArg1(ListFirst(l1)))))
	for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
	{
		Term t2=0;
		char mt[100024], *mtp;
		int p;
		
		p=fprintf(f,"      %s(",AtomValue(CompoundArg1(ListFirst(l1))));
		for(l3=CompoundArg1(ListFirst(l2));l3;l3=ListTail(l3))
			p+=fprintf(f,"%ld%c",IntegerValue(ListFirst(l3)),ListTail(l3)?',':')');
		p+=fprintf(f,"=");		
		
		for(l3=CopyTerm(CompoundArg2(ListFirst(l2)));l3;l3=ListTail(l3))
		{
			Term t3;
			List l4;
			int num=(int)IntegerValue(CompoundArg1(ListFirst(l3)));
			for(l4=CompoundArg2(ListFirst(l3));l4;l4=ListTail(l4))
			{
				Term pp=ListFirst(l4);
				if(CompoundArg1(pp)==A_I) SetCompoundArg(pp,1,fa_i);
				if(AtomValue(CompoundArg1(pp))[0]=='(')
					{
					if(GetAtomProperty(CompoundArg1(pp),A_I)==NewInteger(1))
						p+=fprintf(f,"dconjg(");
					if(GetAtomProperty(CompoundArg1(pp),A_I)==NewInteger(2))
						p+=fprintf(f,"-");
					p+=fprintf(f,"%s",AtomValue(CompoundArg1(ListFirst(l1))));
					}
			}
			t3=l_2_t(ConsumeCompoundArg(ListFirst(l3),2),num,1);
			if(t2==0)
			{
				t2=t3;
				if(num<0)
					t2=MakeCompound1(OPR_MINUS,t2);
			}
			else
				t2=MakeCompound2(num>0?OPR_PLUS:OPR_MINUS,t2,t3);
		}
		
		if(is_atom(t2))
			t2=toccf(t2);
		else
			ccf_term(t2);
		
		fortr_digi=1;
		sWriteTerm(mt,t2);
		fortr_digi=0;
		
		
		mtp=mt;
		while(p+strlen(mtp)>60)
		{
			char sv;
			char *sp=mtp+55-p;
			while(!isalnum(*sp)) sp++;
			while( isalnum(*sp) || *sp=='(') sp++;
			sv=(*sp);(*sp)=0;
			fprintf(f,"%s\n     &      ",mtp);
			(*sp)=sv;
			mtp=sp;
			p=12;
		}
		fprintf(f,"%s\n",mtp);
	}
	
}
	
static void wrt_expr(FILE *of, Term num, List sym, List ten)
{
	int f=1;
	List l;
	int sno=32;
	NoQuotes=1;
	if(is_integer(num) || (is_compound(num)&&IntegerValue(CompoundArg2(num))==1))
	{
		int n=(int)IntegerValue(is_integer(num)?num:CompoundArg1(num));
		if(n==-1 && ((sym&&IntegerValue(CompoundArg2(ListFirst(sym)))>0)||(ten&&(!sym))))
			sno+=fprintf(of,"- "),f=0;
		else if(n!=1 || (is_integer(num)&&(!sym)&&(!ten)) )
			sno+=fprintf(of,"%d ",n),f=0;
	}
	else
	{
		sno+=fprintf(of,"%ld/%ld ",IntegerValue(CompoundArg1(num)),
				IntegerValue(CompoundArg2(num)));
		f=0;
	}
	
	for(l=sym;l;l=ListTail(l))
	{
		Atom p=CompoundArg1(ListFirst(l));
		int  w=(int)IntegerValue(CompoundArg2(ListFirst(l)));
		if(w<0 && f==1)
		{
			sno+=fprintf(of,"1 ");
			f=0;
		}
		if(w<0)
		{
			sno+=fprintf(of,"/ ");
			w=-w;
			if(sno>75) {fprintf(of,"\n\t\t");sno=15;}
		}
		if(sno>75) {fprintf(of,"%c\n\t\t",w>0?'*':' ');sno=15;}
		if(p==A_I)
			sno+=fprintf(of,"I");
		else
			sno+=fWriteTerm(of,p);
		if(w==1)
			sno+=fprintf(of," ");
		else
			sno+=fprintf(of,"^%d ",w);
		f=0;
	}
	
	if(sno>55 && ten) {fprintf(of,"*\n\t\t");sno=15;}
	
	

	if(sno!=32) f=0;
	if(!is_integer(num) && !f)
		fprintf(of,"*");
	NoQuotes=0;
}

static List nexti(List c, List m)
{
	int i;
	if(c==0)
		{
		for(i=0;i<ListLength(m);i++)
			c=AppendFirst(c,NewInteger(1));
		return c;
		}
		
	for(i=ListLength(m);i;i--)
		{
		List p=ListNthList(c,i);
		if(ListFirst(p)==ListNth(m,i))
			{
			ChangeList(p,NewInteger(1));
			continue;
			}
		else
			{
			ChangeList(p,NewInteger(IntegerValue(ListFirst(p))+1));
			break;
			}
		}
	if(i)
		return c;
	else
		return 0;
}	
	

	
void cls_write_dmatr(FILE *f)
{
	
	List l1,l2,l3;
	List ci, mi;


	fprintf(f,"\n\n");
	
	for(l1=mtr_known;l1;l1=ListTail(l1))
	{
	
	if(GetAtomProperty(CompoundArg1(ListFirst(l1)),A_INFINITESIMAL)==0)
		continue;
	
		mi=GetAtomProperty(CompoundArg1(ListFirst(l1)),OPR_CLASS);


		for(ci=nexti(0,mi);ci;ci=nexti(ci,mi))
		{
			int p;	
			
			for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
				if(EqualTerms(ci,CompoundArg1(ListFirst(l2))))
					break;

			p=fprintf(f,"RenConst[%s[",AtomValue(CompoundArg1(ListFirst(l1))));
			for(l3=ci;l3;l3=ListTail(l3))
				p+=fprintf(f,"%ld%c",IntegerValue(ListFirst(l3)),ListTail(l3)?',':']');
			p+=fprintf(f,"] = ");		

			if(l2)
			{
			for(l3=CopyTerm(CompoundArg2(ListFirst(l2)));l3;l3=ListTail(l3))
			{
				Term m2=ListFirst(l3), pp;
				pp=ListFirst(CompoundArg2(m2));
				if(AtomValue(CompoundArg1(pp))[0]=='(')
					{
					char cbf[128];
					int x;
					if(GetAtomProperty(CompoundArg1(pp),A_I)==NewInteger(1))
						{puts("!");p+=fprintf(f,"Conjugate[");}
					if(GetAtomProperty(CompoundArg1(pp),A_I)==NewInteger(2))
						p+=fprintf(f,"-");
					p+=fprintf(f,"%s[",AtomValue(CompoundArg1(ListFirst(l1))));
					sprintf(cbf,"%s",AtomValue(CompoundArg1(pp)));
					for(x=0;cbf[x];x++)
						if(cbf[x]==')') cbf[x]=']';
					fprintf(f,"%s",cbf+1);
					}
				else
					wrt_expr(f,CompoundArg1(m2),CompoundArg2(m2),
							CompoundArgN(m2,3));
				if(ListTail(l3) &&
							IntegerValue(CompoundArg1(ListFirst(ListTail(l3))))>0)
						fprintf(f,"+ ");
			}
			}
			else
				fprintf(f,"0");
			
			fprintf(f,"\n");

		}
	
	}
	
}
	
