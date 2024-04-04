#include <setjmp.h>
#include "lanhep.h"

extern jmp_buf alg1_jmp_buf;

extern int remove_rc;
int infi_order=1;
int opMaxiLegs=0;
static List all_infi=0;

void clear_infi(Atom i)
{
		List l;
		SetAtomProperty(i,A_INFINITESIMAL,0);
		for(l=all_infi;l;l=ListTail(l))
			if(ListFirst(l)==i)
			{
				all_infi=CutFromList(all_infi,l);
				return;
			}
}


Term ProcInfsimal(Term t, Term ind)
	{
	Term t1;
	if(!is_compound(t) || CompoundArity(t)!=1)
		{
		ErrorInfo(701);
		puts("Wrong argument in 'infinitesimal' statement.");
		return 0;
		}
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	t=CommaToList(t1);
	

	for(t1=t;t1;t1=ListTail(t1))
		{
		int pw=1;
		Term val=0;
		Term aa=ListFirst(t1);
		Term at=0;
		
/*		if(is_compound(aa)&&CompoundName(aa)==OPR_DIV&&
				CompoundArg1(aa)==NewInteger(1))
		{
			aa=CompoundArg2(aa);
			pw=-1;
		}
*/
		if(is_compound(aa) && CompoundName(aa)==OPR_EQSIGN /*&&
			is_atom(CompoundArg1(aa))*/ )
		{
			val=ConsumeCompoundArg(aa,2);
			aa=CompoundArg1(aa);
		}

		
		if(is_compound(aa) && CompoundName(aa)==OPR_POW
				&& is_integer(CompoundArg2(aa)))
			{
				pw=(int)IntegerValue(CompoundArg2(aa));
				aa=CompoundArg1(aa);
			}			
		
		if(is_compound(aa) && CompoundName(aa)==OPR_DIV)
			{
			at=CompoundArg2(aa);
			aa=CompoundArg1(aa);
			}
		
		if(!is_atom(aa) || (at && !is_atom(at)))
			{
			ErrorInfo(702);
			puts("Wrong argument in 'infinitesimal' statement.");
			}
		else
			{
			if(!remove_rc)
				{
				SetAtomProperty(aa,A_INFINITESIMAL,
						MakeCompound2(OPR_COMMA,NewInteger(pw),val));
				if(at)
					{
					SetAtomProperty(at,A_INFINITESIMAL,
						MakeCompound2(OPR_COMMA,NewInteger(pw),0));
					SetAtomProperty(aa,A_ANTI,at);
					SetAtomProperty(at,A_ANTI,aa);
					SetAtomProperty(at,A_HERMC,aa);
					if(pw==1)
						SetAtomProperty(aa,OPR_COEFF,NewInteger(1));
					}
				all_infi=AppendLast(all_infi,aa);
				if(pw==1)
					SetAtomProperty(aa,OPR_COEFF,NewInteger(1));
				}
			else
				{
				Term f=MakeCompound1(OPR_LET,MakeCompound2(OPR_EQSIGN,
						aa,NewInteger(0)));
				CallFunction(f,0);
				if(at)
				{
				Term f=MakeCompound1(OPR_LET,MakeCompound2(OPR_EQSIGN,
						at,NewInteger(0)));
				CallFunction(f,0);
				}
				}
				
			}
		}
	FreeAtomic(t);		
	return 0;
	}

List a1l_rem_inf(List);

void alg1_rem_inf(Term a1)
	{
	List l;
	l=ConsumeCompoundArg(a1,1);
	l=a1l_rem_inf(l);
	SetCompoundArg(a1,1,l);	
	}

int inf_removed[10]={0,0,0,0,0,0,0,0,0,0};


			
List a1l_rem_inf(List l)
{
	List l1;
	int ch=0;
	
	for(l1=l;l1;l1=l1?ListTail(l1):0)
		{
		Term m1;
		List l2;
		int  ino;
		int fno;

		ino=0;
		fno=0;
		m1=ListFirst(l1);
		for(l2=CompoundArgN(m1,3);l2;l2=ListTail(l2))
			{
			Term prp=0;
			if(CompoundName(ListFirst(l2))==OPR_PARAMETER)
				prp=GetAtomProperty(CompoundArg2(ListFirst(l2)),A_INFINITESIMAL);
			if(prp && IntegerValue(CompoundArg1(prp))>0)
				ino+=(int)IntegerValue(CompoundArg1(prp));
			if(CompoundName(ListFirst(l2))==OPR_FIELD && 
					CompoundArg2(ListFirst(l2))!=A_VEV)
				fno++;
			if(CompoundName(ListFirst(l2))==OPR_WILD && 
					CompoundArg2(ListFirst(l2))==0)
				ino+=10000;
			}
		for(l2=CompoundArgN(m1,4);l2;l2=ListTail(l2))
			{
			Term prp=0;
			if(CompoundName(ListFirst(l2))==OPR_PARAMETER)
				prp=GetAtomProperty(CompoundArg2(ListFirst(l2)),A_INFINITESIMAL);
			if(prp && IntegerValue(CompoundArg1(prp))<0)
				ino++;
			}
			
		if(opMaxiLegs && fno>opMaxiLegs)
			{
			FreeAtomic(ListFirst(l1));
			ChangeList(l1,0);
			ch++;
			continue;
			}

		if(ino>infi_order)
			{
			if(ino>9) inf_removed[9]++;
			else
				inf_removed[ino]++;
			FreeAtomic(ListFirst(l1));
			ChangeList(l1,0);
			ch++;
			continue;
			}
		else
		for(l2=CompoundArgN(m1,4);l2;l2=ListTail(l2))
			{
			Term prp=0;
			if(CompoundName(ListFirst(l2))==OPR_PARAMETER)
				prp=GetAtomProperty(CompoundArg2(ListFirst(l2)),A_INFINITESIMAL);
			if(prp && IntegerValue(CompoundArg1(prp))>0)
				{
				ErrorInfo(703);
				printf("infinitesimal parameter '%s' in denominator\n",
					AtomValue(CompoundArg2(ListFirst(l2))));
				longjmp(alg1_jmp_buf,1);
				}
			}
			
		if(l1 && ino && ino==infi_order)
			{
			int insav=infi_order;
			infi_order=0;
			for(l2=CompoundArgN(m1,3);l2;l2=ListTail(l2))
				{
				List l3;
				if(CompoundName(ListFirst(l2))!=OPR_WILD)
					continue;
				for(l3=CompoundArg2(ListFirst(l2));l3;l3=ListTail(l3))
					alg1_rem_inf(ListFirst(l3));
				for(l3=CompoundArg2(ListFirst(l2));l3;l3=ListTail(l3))
					if(!is_empty_list(CompoundArg1(ListFirst(l3))))
						break;
				if(l3==0)
					ino+=10000;				
				}
			infi_order=insav;
			}
			
		if(ino>infi_order)
			{
			if(ino>1000) inf_removed[0]++;
			else if(ino>9) inf_removed[9]++;
			else
				inf_removed[ino]++;
			FreeAtomic(ListFirst(l1));
			ChangeList(l1,0);
			ch++;
			}
			
		}
		
	if(ch)
	{
		List ml2=0, mle=0;
		for(l1=l;l1;l1=ListTail(l1))
		if(ListFirst(l1))
		{
			if(ml2==0)
			{
			ml2=AppendLast(ml2,ListFirst(l1)); mle=ml2;
			}
		else
			{
			AppendLast(mle,ListFirst(l1));
			mle=ListTail(mle);
			}
		}
		RemoveList(l);
		l=ml2;
	}
		
	return l;
	}	

static Atom udprt(Atom p, char c)
{
  char cbuf[8];
  sprintf(cbuf,"%s.%c",AtomValue(p),c);
  return NewAtom(cbuf,0);
}

static void reg_spi_transf(Atom p, Term a1, int anti)
{
  Term tu=udprt(p,'u');
  Term td=udprt(p,'d');
  List l;
  int err=0;

  if(GetAtomProperty(tu,OPR_LET)||GetAtomProperty(td,OPR_LET))
	  return;
  
  for(l=CompoundArg1(a1);l;l=ListTail(l))
    {
      Atom prt=0, z=0, gpm=0;
      List l1;
      for(l1=CompoundArgN(ListFirst(l),3);l1;l1=ListTail(l1))
	{
	    Term sp=ListFirst(l1);
				if(CompoundName(sp)==OPR_FIELD)
				{
					if(prt) err=1;
					prt=CompoundArg2(sp);
					continue;
				}
				if(CompoundName(sp)==OPR_SPECIAL &&
					(CompoundArg2(sp)==A_GAMMAP || CompoundArg2(sp)==A_GAMMAM))
				{
					if(gpm) err=1;
					gpm=CompoundArg2(sp);
					continue;
				}
				if(CompoundName(sp)==OPR_PARAMETER /*&&
						GetAtomProperty(CompoundArg2(sp),A_INFINITESIMAL)*/ )
				{
					if(z) err=1;
					z=CompoundArg2(sp);
					continue;
				}
	}
	 if(prt && !z && !gpm)
		 continue;

      if(z==0||prt==0||gpm==0)
	  {
		  err=1;
		continue;
	  }
/*	WriteTerm(CompoundArg1(ListFirst(l)));puts("");
	WriteTerm(CompoundArg2(ListFirst(l)));puts("");
	DumpList(CompoundArgN(ListFirst(l),3));
	WriteTerm(CompoundArgN(ListFirst(l),4));puts("");
*/	
      if(gpm==A_GAMMAM)
	  {
	  Term mm=MakeCompound2(OPR_MLT,z,udprt(prt,'u'));
	  if(IntegerValue(CompoundArg2(ListFirst(l)))<0)
	  	mm=MakeCompound2(OPR_MLT,A_I,mm);
	  tu=MakeCompound2(
	  	IntegerValue(CompoundArg1(ListFirst(l)))>0?OPR_PLUS:OPR_MINUS,tu,
			 MakeCompound2(OPR_DIV,mm,NewInteger(2)));
	  }
      else
	  {
	  Term mm=MakeCompound2(OPR_MLT,z,udprt(prt,'d'));
	  if(IntegerValue(CompoundArg2(ListFirst(l)))<0)
	  	mm=MakeCompound2(OPR_MLT,A_I,mm);	  	
	  td=MakeCompound2(
	   IntegerValue(CompoundArg1(ListFirst(l)))>0?OPR_PLUS:OPR_MINUS,td,
			 MakeCompound2(OPR_DIV,mm,NewInteger(2)));
	  }
    }
	
	if(err)
	{
		ErrorInfo(0);
		printf("transformation for two-component up(%s) and down(%s) ",
				AtomValue(p), AtomValue(p));
		puts("are not set.");
		return;
	}
  if(!is_atom(tu))
    {
      tu=MakeCompound2(OPR_RARROW,udprt(p,'u'),tu);
      /*WriteTerm(tu);puts("");*/ 
      ProcLet(MakeCompound1(OPR_LET,tu),0);
    }
  if(!is_atom(td))
  {
    td=MakeCompound2(OPR_RARROW,udprt(p,'d'),td);
    /*WriteTerm(td);puts(""); */
    ProcLet(MakeCompound1(OPR_LET,td),0);
  }


}



Term ProcTransform(Term t, Term ind)
	{
	Term t1;
	
	if(remove_rc)
	{
		FreeAtomic(t);
		return 0;
	}
	
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	
	if(is_compound(t1) && CompoundArity(t1)==2 && CompoundName(t1)==OPR_COMMA)
		{
		Term a1,a2;
		a1=ConsumeCompoundArg(t1,1);
		a2=ConsumeCompoundArg(t1,2);
		FreeAtomic(t1);
		ProcTransform(MakeCompound1(OPR_LET,a1),0);
		ProcTransform(MakeCompound1(OPR_LET,a2),0);
		return 0;
		}
		
	if(!is_compound(t1) || CompoundArity(t1)!=2 || CompoundName(t1)!=OPR_RARROW)
		{
		ErrorInfo(725);
		printf("bad argument in 'transform' statement\n");
		return 0;
		}
	
	t=CompoundArg1(t1);
	ProcLet(MakeCompound1(OPR_LET,t1),0);
	if(!is_atom(t)) return 0;
	t1=GetAtomProperty(t,OPR_LET);
	if(t1==0) return 0;
	t1=CompoundArg1(t1);
	
/*	WriteTerm(t); printf(" -> "); WriteTerm(t1);puts("");*/ 
	
	if(is_parameter(t) && GetAtomProperty(t,OPR_MASS))
	{
		List inff=0;
		List l, l1;
		int slf=0;
		Atom prt=GetAtomProperty(t,OPR_MASS);
		for(l=CompoundArg1(t1);l;l=ListTail(l))
		{
			Term m1=ListFirst(l);
			if(ListLength(CompoundArgN(m1,3))==1 &&
				CompoundArg2(ListFirst(CompoundArgN(m1,3)))==t)
			{
				slf++;
				continue;
			}
			
			for(l1=CompoundArgN(m1,3);l1;l1=ListTail(l1))
			{
				Term sp=ListFirst(l1);
				if(CompoundName(sp)==OPR_PARAMETER && 
						GetAtomProperty(CompoundArg2(sp),A_INFINITESIMAL))
					inff=AppendLast(inff,CompoundArg2(sp));
			}
			
		}
			
		if(slf!=1 || ListLength(inff)!=1)
		{
			WarningInfo(0);
			printf("non-standard RC for %s (%s mass).\n",
					AtomValue(t),AtomValue(prt));
			return 0;
		}

		t1=GetAtomProperty(ListFirst(inff),A_INFINITESIMAL);
		if(CompoundArg2(t1)==0)
			SetCompoundArg(t1,2,MakeCompound1(OPR_MASS,
				MakeCompound2(OPR_DIV,prt,t)));
		return 0;
	}
			
	
	if(is_particle(t,0))
	{ 
		List l;
		int slf=0;
		Term fprp=0;
		int sp;
		Atom mass1;
		Term prp=GetAtomProperty(t,PROP_TYPE);
		if(CompoundName(prp)==OPR_PARTICLE &&
		   CompoundArgN(prp,4)==NewInteger(1))
		  reg_spi_transf(t,t1,CompoundArg1(prp)!=CompoundArg2(prp) &&
				t==CompoundArg2(prp));
		if(CompoundName(prp)==OPR_FIELD)
		{
			fprp=prp;
			prp=GetAtomProperty(CompoundArg1(fprp),PROP_TYPE);
			if(CompoundArg1(prp)!=CompoundArg2(prp) &&
					CompoundArg1(fprp)==CompoundArg2(prp))
				return 0;
		}
		if(CompoundName(prp)==OPR_PARTICLE &&
				CompoundArg1(prp)!=CompoundArg2(prp) &&
				t==CompoundArg2(prp))
			return 0;
		
		mass1=CompoundArgN(prp,5);
		
		l=CompoundArg2(t1);
		if(l==0)
			sp=0;
		else
		{
			Term isp=ListFirst(l);
			if(CompoundName(CompoundArg1(isp))!=A_LORENTZ)
				sp=0;
			else if(CompoundArg1(CompoundArg1(isp))==NewInteger(2))
				sp=2;
			else sp=1;
		}
		
		for(l=CompoundArg1(t1);l;l=ListTail(l))
		{
			Term m1=ListFirst(l);
			Term prp2, mass2, tp;
			List l1,fi=0, inff=0, gpm=0;
			
			if(ListLength(CompoundArgN(m1,3))==1 &&
				CompoundArg2(ListFirst(CompoundArgN(m1,3)))==t)
			{
				slf++;
				continue;
			}
			for(l1=CompoundArgN(m1,3);l1;l1=ListTail(l1))
			{
				Term sp=ListFirst(l1);
				if(CompoundName(sp)==OPR_FIELD)
				{
					fi=AppendLast(fi,CompoundArg2(sp));
					continue;
				}
				if(CompoundName(sp)==OPR_SPECIAL &&
					(CompoundArg2(sp)==A_GAMMAP || CompoundArg2(sp)==A_GAMMAM))
				{
					gpm=AppendLast(gpm,CompoundArg2(sp));
					continue;
				}
				if(CompoundName(sp)==OPR_PARAMETER &&
						GetAtomProperty(CompoundArg2(sp),A_INFINITESIMAL))
				{
					inff=AppendLast(inff,CompoundArg2(sp));
					continue;
				}
			}
			
			if(ListLength(inff)!=1 || ListLength(fi)!=1 ||
					(sp==1 && ListLength(gpm)!=1))
			{
				if(ListLength(inff))
				{
				WarningInfo(0);printf("non-standard rc for field %s\n",
						AtomValue(t));
				}
				continue;
			}
			
			prp2=GetAtomProperty(ListFirst(fi),PROP_TYPE);
			if(CompoundName(prp2)==OPR_FIELD)
				prp2=GetAtomProperty(CompoundArg1(prp2),PROP_TYPE);
			mass2=CompoundArgN(prp2,5);
			
			if(sp==0) tp=OPR_SCALAR;
			else if(sp==2) tp=OPR_VECTOR;
			else if(ListFirst(gpm)==A_GAMMAM) tp=A_LEFT;
			else tp=A_RIGHT;
			
			tp=MakeCompound1(tp,MakeCompound2(OPR_RARROW,
					MakeCompound2(OPR_DIV,t,mass1),
					MakeCompound2(OPR_DIV,ListFirst(fi),mass2)));
			RemoveList(fi);
			prp2=GetAtomProperty(ListFirst(inff),A_INFINITESIMAL);
			if(CompoundArg2(prp2)==0)
				SetCompoundArg(prp2,2,tp);
			else
				FreeAtomic(tp);
		}
	
	}
	
	return 0;
	}
	
static void alg1_dlt_to_wld(List ll, Label ind1, Label ind2, int d)
{
	List tind1, tind2;
	List l1;
	List r1,r2;
	Term ww;
	int i,j;
	
	l1=ConsumeCompoundArg(ListFirst(ll),1);
	FreeAtomic(ListFirst(ll));
	ChangeList(ll,0);
	
	tind1=ListFirst(l1);
	tind2=ListFirst(ListTail(l1));
	RemoveList(l1);
	
	
	r1=0;
	
	for(i=1;i<=d;i++)
	{
		Term mm, ww, aa;
		r2=0;
		
		for(j=1;j<=d;j++)
		{
			Term aa1;
			Term aa1ml=0;
			
			if(i==j)
			{	
				Term m1;
				m1=MakeCompound(A_MTERM,4);
				
				SetCompoundArg(m1,1,NewInteger(1));
				SetCompoundArg(m1,2,NewInteger(1));
				
				aa1ml=AppendLast(aa1ml,m1);
			}
			
			aa1=MakeCompound2(A_ALG1,aa1ml,0);
			r2=AppendLast(r2,aa1);
		}
		
		ww=MakeCompound2(OPR_WILD,MakeList1(CopyTerm(tind1)),r2);
		mm=MakeCompound(A_MTERM,4);
		SetCompoundArg(mm,1,NewInteger(1));
		SetCompoundArg(mm,2,NewInteger(1));
		SetCompoundArg(mm,3,MakeList1(ww));
		aa=MakeCompound2(A_ALG1,MakeList1(mm),MakeList1(CopyTerm(tind1)));
		r1=AppendLast(r1,aa);
		
	}
	
	ww=MakeCompound2(OPR_WILD,MakeList2(tind1,tind2),r1);
	ChangeList(ll,ww);
		
}
	


void alg1_fix_delta(Term a1)
{
	List l1,l2,il;
	
	il=NewList();
	for(l1=CompoundArg2(a1);l1;l1=ListTail(l1))
		il=AppendLast(il,CompoundArg2(ListFirst(l1)));
	
	
	for(l1=CompoundArg1(a1);l1;l1=ListTail(l1))
	{
		List ml;
		ml=ConsumeCompoundArg(ListFirst(l1),3);
		
	cnt:
		for(l2=ml;l2;l2=ListTail(l2))
		{
			if(CompoundArg2(ListFirst(l2))==A_DELTA)
			{
				Label in1, in2;
				int cno;
				
				in1=CompoundArg2(ListFirst(CompoundArg1(ListFirst(l2))));
				in2=CompoundArg2(ListNth(CompoundArg1(ListFirst(l2)),2));
				if(CompoundName(ListFirst(CompoundArg1(ListFirst(l2))))
																==OPR_WILD)
					cno=(int)IntegerValue(CompoundArgN(
							ListFirst(CompoundArg1(ListFirst(l2))),3));
				else
					cno=0;
				
				if(in1==in2)
				{
					Term r1,r2,mm;
					r1=CompoundArg1(ListFirst(CompoundArg1(ListFirst(l2))));
					r2=CompoundArg1(ListNth(CompoundArg1(ListFirst(l2)),2));
					mm=ListFirst(l1);
					if(CompoundName(r1)==A_LORENTZ 
						&& CompoundName(r2)==A_LORENTZ 
						&& CompoundArg1(r1)==NewInteger(2)
						&& CompoundArg1(r2)==NewInteger(2))
					{
						ml=CutFromList(ml,l2);
						SetCompoundArg(mm,1,
							NewInteger(4*IntegerValue(CompoundArg1(mm))));
						goto cnt;
					}
					puts("Internal error (rmdlt1)");
					WriteTerm(r1);puts("");
					WriteTerm(r2);puts("");
					exit(0);
				}
				
				if(!cno && ListMember(il,in1) && ListMember(il,in2))
					continue;
				if(ListMember(il,in1) && ListMember(il,in2))
				{
					alg1_dlt_to_wld(l2,in1,in2,cno);
					continue;
				}
				
				if(ListMember(il,in1))
				{
					Label tmp;
					tmp=in1;
					in1=in2;
					in2=tmp;
				}
				
				ml=CutFromList(ml,l2);
				
				for(l2=ml;l2;l2=ListTail(l2))
				{
					List l3;
					for(l3=CompoundArg1(ListFirst(l2));l3;l3=ListTail(l3))
						if(CompoundArg2(ListFirst(l3))==in1)
						{
							SetCompoundArg(ListFirst(l3),2,in2);
							break;
						}
					if(l3)
						break;
				}
				
				if(!l2)
				{
					puts("Internal error (rmdlt2)");
					exit(0);
				}
				
				goto cnt;
			}
		}
		SetCompoundArg(ListFirst(l1),3,ml);
	}
	
	
}

void inf_decl_hc(FILE *f)
{
	List l;
	for(l=all_infi;l;l=ListTail(l))
	{
		Atom z=ListFirst(l);
		Atom az=GetAtomProperty(z,A_ANTI);
		if(az)
			fprintf(f,"%s := Conjugate[%s]\n",AtomValue(az),AtomValue(z));
	}
}

void inf_write_rc(FILE *f)
{
	List l;
/*	f=stdout;
	for(l=all_infi;l;l=ListTail(l))
	{
		Term prp;
		printf("%s\t",AtomValue(ListFirst(l)));
		prp=GetAtomProperty(ListFirst(l),A_INFINITESIMAL);
		WriteTerm(prp);puts("");
	}
*/	
	for(l=all_infi;l;l=ListTail(l))
	{
		Atom z=ListFirst(l);
		Term t=CompoundArg2(GetAtomProperty(z,A_INFINITESIMAL));
		if(t==0)
		{
                        if(!GetAtomProperty(z,PROP_TYPE))
			printf("Undefinite RenConst[%s]\n",AtomValue(z));
			continue;
		}
		
		if(is_compound(t) && CompoundName(t)==OPR_MASS)
		{
			Atom prt=CompoundArg1(CompoundArg1(t));
			Atom mas=CompoundArg2(CompoundArg1(t));
			Term prp=GetAtomProperty(prt,PROP_TYPE);
			if(prp==0 || CompoundName(prp)!=OPR_PARTICLE)
			{
				puts("internal error a1hmsnp");
				continue;
			}
			fprintf(f,"RenConst[ %s ] := ",AtomValue(z));
			switch(IntegerValue(CompoundArgN(prp,4)))
			{
				case 0:
					fprintf(f,"ReTilde[SelfEnergy[prt[\"%s\"] -",AtomValue(prt));
					fprintf(f,"> prt[\"%s\"], %s]]\n",AtomValue(prt),AtomValue(mas));
					break;
				case 1:
					fprintf(f,"Block[ {sff}, sff = SelfEnergy[prt[\"%s\"] -> prt[\"%s\"]",
							AtomValue(prt),AtomValue(prt));
					fprintf(f,", %s];\n\tReTilde[ %s/2 (LVectorCoeff[sff]+RVectorCoeff[sff])",
							AtomValue(mas),AtomValue(mas));
					fprintf(f,"+LScalarCoeff[sff]]]\n");
					break;
				case 2:
					fprintf(f,"-VSESign ReTilde[SelfEnergy[prt[\"%s\"] -",AtomValue(prt));
					fprintf(f,"> prt[\"%s\"], %s]]\n",AtomValue(prt),AtomValue(mas));
					break;
			}
			continue;
		}
		
		if(is_compound(t) && CompoundName(t)==OPR_SCALAR)
		{
			Term t1=CompoundArg1(t);
			Atom prt1=CompoundArg1(CompoundArg1(t1));
			Atom mas1=CompoundArg2(CompoundArg1(t1));
			Atom prt2=CompoundArg1(CompoundArg2(t1));
			Atom mas2=CompoundArg2(CompoundArg2(t1));
			char *m1="0", *m2="0";
			/*if(prt1!=prt2)
			{
				printf("Can not generate RC[%s] for %s->%s\n",AtomValue(z),
						AtomValue(prt1),AtomValue(prt2));
				continue;
			}*/
			if(mas1) m1=AtomValue(mas1);
			if(mas2) m2=AtomValue(mas2);
			fprintf(f,"RenConst[ %s ] := ",AtomValue(z));
			if(prt1==prt2)
			fprintf(f,"-ReTilde[DSelfEnergy[prt[\"%s\"] -> prt[\"%s\"], %s]]\n",
					AtomValue(prt1),AtomValue(prt2),m2);
			else
			fprintf(f,"1/(%s^2-%s^2) ReTilde[SelfEnergy[prt[\"%s\"] -> prt[\"%s\"], %s]]\n",
				m1,m2,AtomValue(prt1),AtomValue(prt2),m2);
			continue;
		}
		if(is_compound(t) && CompoundName(t)==OPR_VECTOR)
		{
			Term t1=CompoundArg1(t);
			Atom prt1=CompoundArg1(CompoundArg1(t1));
			Atom mas1=CompoundArg2(CompoundArg1(t1));
			Atom prt2=CompoundArg1(CompoundArg2(t1));
			Atom mas2=CompoundArg2(CompoundArg2(t1));
			if(prt1!=prt2 && ((mas1&&mas2) || (!mas1 && !mas2)))
			{
				printf("Can not generate RC[%s] for %s->%s\n",AtomValue(z),
						AtomValue(prt1),AtomValue(prt2));
				continue;
			}
			fprintf(f,"RenConst[ %s ] := ",AtomValue(z));
			if(prt1==prt2)
			fprintf(f,"VSESign ReTilde[DSelfEnergy[prt[\"%s\"] -> prt[\"%s\"], %s]]\n",
					AtomValue(prt1),AtomValue(prt2),mas1?AtomValue(mas1):"0");
			else if(mas1==0)
			fprintf(f,"2 VSESign /%s^2 ReTilde[SelfEnergy[prt[\"%s\"] -> prt[\"%s\"], %s]]\n",
					AtomValue(mas2),AtomValue(prt1),AtomValue(prt2),AtomValue(mas2));
			else
			fprintf(f,"-2 VSESign /%s^2 ReTilde[SelfEnergy[prt[\"%s\"] -> prt[\"%s\"], %s]]\n",
					AtomValue(mas1),AtomValue(prt2),AtomValue(prt1),"0");
			continue;
		}
		if(is_compound(t) && (CompoundName(t)==A_LEFT || CompoundName(t)==A_RIGHT))
		{
			Term t1=CompoundArg1(t);
			Atom prt1=CompoundArg1(CompoundArg1(t1));
			Atom mas1=CompoundArg2(CompoundArg1(t1));
			Atom prt2=CompoundArg1(CompoundArg2(t1));
			Atom mas2=CompoundArg2(CompoundArg2(t1));
			char tp=(CompoundName(t)==A_LEFT)?'L':'R';
			
			if(prt1!=prt2 && (mas1==0 || mas2==0))
			{
				printf("Can not generate RC[%s] for %s->%s\n",AtomValue(z),
						AtomValue(prt1),AtomValue(prt2));
				continue;
			}
			
			fprintf(f,"RenConst[ %s ] := ",AtomValue(z));
			if(prt1==prt2)
			{
			fprintf(f,"Block[ {sff,dsff},\n\t sff=SelfEnergy[prt[\"%s\"] -> prt[\"%s\"] , %s];\n",
					AtomValue(prt1),AtomValue(prt1),mas1?AtomValue(mas1):"0");
			fprintf(f,"\tdsff=DSelfEnergy[prt[\"%s\"] -> prt[\"%s\"] , %s];\n",
					AtomValue(prt1),AtomValue(prt1),mas1?AtomValue(mas1):"0");
			fprintf(f,"\t-ReTilde[ ");
			if(mas1) 
			{
				fprintf(f,"%s^2 (LVectorCoeff[dsff] + RVectorCoeff[dsff]) +\n\t",AtomValue(mas1));
    			fprintf(f,"2 %s LScalarCoeff[dsff]+",AtomValue(mas1));
			}
			fprintf(f,"%cVectorCoeff[sff] ] ]\n",tp);
			}
			else
			{
			fprintf(f,"Block[ {sff},\n\t");
			fprintf(f,"sff = SelfEnergy[prt[\"%s\"] -> prt[\"%s\"], %s];\n\t",
					AtomValue(prt1),AtomValue(prt2),AtomValue(mas2));
			if(tp=='L')
			{
			fprintf(f,"2/(%s^2 - %s^2) ReTilde[ \n\t",AtomValue(mas1),AtomValue(mas2));
			fprintf(f,"%s^2 LVectorCoeff[sff] + %s %s RVectorCoeff[sff] +\n\t",
					AtomValue(mas2),AtomValue(mas1),AtomValue(mas2));
			fprintf(f,"(%s^2 + %s^2)/%s LScalarCoeff[sff] ] ]\n",
					AtomValue(mas1),AtomValue(mas2),AtomValue(mas1));
			}
			else
			{
			fprintf(f,"2 %s/(%s^2 - %s^2) ReTilde[\n\t",
					AtomValue(mas2),AtomValue(mas1),AtomValue(mas2));
			fprintf(f,"%s RVectorCoeff[sff] + %s LVectorCoeff[sff] + 2 LScalarCoeff[sff] ] ]\n",
					AtomValue(mas2),AtomValue(mas1));
			}
			}
			
			continue;
		}
		fprintf(f,"RenConst[ %s ] := ",AtomValue(z));
		NoQuotes=1;
		fWriteTerm(f,t);
		NoQuotes=0;
		fprintf(f,"\n");

	}
}

	
