#include "lanhep.h"

static List used_fields = 0;

static char ch00 = '~', ch1 = '0', ch2='0';

static char nnbuf[48];

static int intr_param = 0;

int verb_imprt = 0;


static Atom NewName(void)
	{
	Atom ret;
	do
		{
		sprintf(nnbuf,"~%c%c",ch1,ch2);
		ch2++;
		if(ch2-1=='9')
			ch2='A';
		if(ch2-1=='Z')
			ch2='a';
		if(ch2-1=='z')
			{
			ch2='0';
			ch1++;
			if(ch1-1=='9')
				ch1='A';
			if(ch1-1=='Z')
				ch1='a';
			if(ch1-1=='z')
				{
				if(ch00=='~')
					puts("Internal error: too many im particles, output incorrect...");
				ch1=ch2='0';
				ch00--;
				}
			}
		ret=NewAtom(nnbuf,3);
		} while(is_particle(ret,NULL));
	return ret;
	}

void check_hint(int col, int spin, int ch, List hint, Atom *n, Atom *an)
	{
	List l, l1;
	if(col==0)
		return;
	l=hint;
	while(!is_empty_list(l))
		{
		Atom a1,a2;
		Term cp, tp;
		a1=CompoundArg1(ListFirst(l));
		tp=GetAtomProperty(a1,PROP_TYPE);
		if(tp==0 || !is_compound(tp) || CompoundName(tp)!=OPR_PARTICLE)
			goto cnt;
		cp=GetAtomProperty(a1,A_COLOR);
		if(cp==0 || !is_compound(cp))
			goto cnt;
		cp=CompoundArg1(cp);
		if((col==3 && IntegerValue(cp)!=3) ||
			(col!=3 && IntegerValue(cp)==3))
			goto cnt;
		if(ch && CompoundArg1(tp)==CompoundArg2(tp))
			goto cnt;
		if(!ch && CompoundArg1(tp)!=CompoundArg2(tp))
			goto cnt;
		if(ch)
			{
			if(a1==CompoundArg1(tp))
				a2=CompoundArg2(tp);
			else
				a2=CompoundArg1(tp);
			}
		else
			a2=a1;
		if(ch)
			{
			l1=ListTail(l);
			while(!is_empty_list(l1))
				{
				if(a2==CompoundArg1(ListFirst(l1)))
					break;
				l1=ListTail(l1);
				}
			if(is_empty_list(l1))
				goto cnt;
			}
		l1=used_fields;
		while(!is_empty_list(l1))
			{
			Atom aa;
			aa=CompoundArg2(CompoundArg2(ListFirst(l1)));
			if(aa==a1 || aa==a2)
				goto cnt;
			l1=ListTail(l1);
			}
		a1=CompoundArg1(tp);
		a2=CompoundArg2(tp);
		if(spin==2)
			{
			*n=a1;
			*an=a2;
			return;
			}
		if(AtomValue(a1)[0]=='~' || AtomValue(a2)[0]=='~')
			goto cnt;
		sprintf(nnbuf,"~%s",AtomValue(a1));
		a1=NewAtom(nnbuf,0);
		sprintf(nnbuf,"~%s",AtomValue(a2));
		a2=NewAtom(nnbuf,0);
		if(is_particle(a1,NULL) || is_particle(a2,NULL))
			goto cnt;
		*n=a1;
		*an=a2;
		return;
	cnt:
		l=ListTail(l);
		}
	}


static int imppdg=100500;

List mk_im_field(int col, int spin, int ch, List hint)
	{
	Atom n, an, spp;
	int  checked_hint=0;
	Term prt;
	int hc_corr=0;
		
	if(ListLength(hint)==4)
	{
	  Atom nms[4];
	  Term prp;
	  int i;
	  nms[0]=CompoundArg1(ListFirst(hint));
	  nms[1]=CompoundArg1(ListNth(hint,2));
	  nms[2]=CompoundArg1(ListNth(hint,3));
	  nms[3]=CompoundArg1(ListNth(hint,4));
	  for(i=0;i<4;i++)
	  if((prp=GetAtomProperty(nms[i],PROP_TYPE)) && CompoundName(prp)==OPR_FIELD
	    && CompoundArg2(prp)==NewInteger(4))
	      nms[i]=CompoundArg1(prp);
	  if((nms[0]==GetAtomProperty(nms[2],A_ANTI) || nms[0]==GetAtomProperty(nms[3],A_ANTI)) &&
	      (nms[1]==GetAtomProperty(nms[2],A_ANTI) || nms[1]==GetAtomProperty(nms[3],A_ANTI)))
	    hc_corr=1;
	}
	
	n=0;
	check_hint(col, spin, ch, hint, &n, &an);
	
	if(n) checked_hint=1;

	if(n==0)
		{
		n=NewName();
		if(ch)
			an=NewName();
		else
			an=n;
		}

	sprintf(nnbuf,"%s%s%s%s",AtomValue(CompoundArg1(ListFirst(hint))),
			AtomValue(CompoundArg1(ListNth(hint,2))),
			AtomValue(CompoundArg1(ListNth(hint,3))),
			AtomValue(CompoundArg1(ListNth(hint,4))));
	spp=NewAtom(nnbuf,0);
	used_fields=AppendLast(used_fields, MakeCompound2(OPR_MINUS,
		spp, MakeCompound2(OPR_DIV, n, an)));

    
    
	if(verb_imprt)
		{
		printf("Intermediate field %s/%s for vertex ",AtomValue(n),
			AtomValue(an));
		WriteVertex(hint);
		printf(".\n");
		printf("hint= %d ",checked_hint);WriteTerm(hint);puts("");
		}

	if(!checked_hint || spin<2)
		{
		char bbb[64];
		prt=MakeCompound(OPR_PARTICLE,8);
		SetAtomProperty(n,NewAtom("_pdg",0),NewInteger(imppdg++));
		SetCompoundArg(prt,1,n);
		SetCompoundArg(prt,2,an);
		sprintf(bbb,"%s",AtomValue(spp));
		SetCompoundArg(prt,3,NewAtom(bbb,0));
		/*if(CalcOutput)
		  SetCompoundArg(prt,4,NewInteger(spin*2));
		else*/
		  SetCompoundArg(prt,4,NewInteger(spin?2:0));
		SetCompoundArg(prt,5,(spin==2&&!CalcOutput)?0:NewAtom("Maux",0));
		if(spin<2 || CalcOutput)
		{
		  if(CalcOutput && !hc_corr)
		    SetCompoundArg(prt,7,NewAtom("!*",0));
		  else
			SetCompoundArg(prt,7,OPR_MLT);
		}
	  else
		    SetCompoundArg(prt,7,A_GAUGE);
		if(col)
			{
			Term cl;
			if(col==6)
				cl=MakeCompound2(A_COLOR,NewAtom("c6",0),NewAtom("c6b",0));
			else if(col==3)
				cl=MakeCompound2(A_COLOR,NewAtom("c8",0),NewAtom("c8",0));
			else
				cl=MakeCompound2(A_COLOR,NewAtom("c3",0),NewAtom("c3b",0));
			SetCompoundArg(prt,8,AppendFirst(NewList(),cl));
			}

		AddIMParticle(prt);
		}

	if(spin<2)
		{

		if(!intr_param)
			{
			Term t1,t2,t3;

			t1=NewCompound(NewFunctor(OPR_PARAMETER,1));
			t2=NewCompound(NewFunctor(OPR_EQSIGN,2));
			t3=NewCompound(NewFunctor(OPR_COLON,2));
			SetCompoundArg(t1,1,t2);
			SetCompoundArg(t2,1,NewAtom("Maux",0));
			SetCompoundArg(t2,2,t3);
			SetCompoundArg(t3,1,NewInteger(1));
			SetCompoundArg(t3,2,NewAtom("mass of im particles",0));
			ProcessParameter(t1,0);
			intr_param=1;
			}
		return MakeList2(an,n);
		}

		
	/*if(CalcOutput && checked_hint==0)
	   return MakeList2(an,n);*/
	
	sprintf(nnbuf,"%s.t",AtomValue(n));
	n=NewAtom(nnbuf,0);

	sprintf(nnbuf,"%s.t",AtomValue(an));
	an=NewAtom(nnbuf,0);
	return MakeList2(an,n);
	}


Term ProcImPrt(Term t,  Term ind)
	{
	int spin, col, anti;
	Atom n, an;
	Term prt;
	Atom pcomm=0, ptname=0;
	
	if(is_compound(t) && CompoundArity(t)>2 && 
		is_compound(CompoundArg1(t))&&CompoundName(CompoundArg1(t))==OPR_DIV)
	{
		Term r;
		r=CompoundArg1(t);
		if(CompoundArg2(r)==NewInteger(2))
		{
			if(CompoundArg1(r)==NewInteger(1))
				SetCompoundArg(t,1,NewInteger(-1));
			else if(CompoundArg1(r)==NewInteger(3))
				SetCompoundArg(t,1,NewInteger(-3));
		}
	}
	
	if(!is_compound(t) || CompoundArity(t)<3 || CompoundArity(t)>5 ||
		!is_integer(CompoundArg1(t)) || !is_integer(CompoundArg2(t)) ||
		!is_integer(CompoundArgN(t,3)) ||
		(CompoundArity(t)>3 && !is_atom(CompoundArgN(t,4))) ||
		(CompoundArity(t)>4 && !is_atom(CompoundArgN(t,5))) )
		{
		ErrorInfo(502);
		printf("AuxPrt: bad arguments ");WriteTerm(t);puts("");
		return 0;
		}
	

	spin=(int)IntegerValue(CompoundArg1(t));
	col= (int)IntegerValue(CompoundArg2(t));
	anti=(int)IntegerValue(CompoundArgN(t,3));
	if(CompoundArity(t)>3)
		pcomm=CompoundArgN(t,4);
	if(CompoundArity(t)>4)
		ptname=CompoundArgN(t,5);
	FreeAtomic(t);

	n=NewName();
	if(anti)
		an=NewName();
	else
		an=n;


	if(1)
		{
		prt=MakeCompound(OPR_PARTICLE,8);
		SetAtomProperty(n,NewAtom("_pdg",0),NewInteger(imppdg++));
		SetCompoundArg(prt,1,n);
		SetCompoundArg(prt,2,an);
		SetCompoundArg(prt,3,pcomm?pcomm:NewAtom("im particle",0));
		/*if(CalcOutput)
		  SetCompoundArg(prt,4,NewInteger((spin>0)?(spin*2):(-spin)));
		else*/
		  SetCompoundArg(prt,4,NewInteger((spin>0)?(spin?2:0):(-spin)));
		SetCompoundArg(prt,5,(spin==2&&!CalcOutput)?0:NewAtom("Maux",0));
		if(spin<2 || CalcOutput || UFOutput)
			SetCompoundArg(prt,7,anti>0?OPR_MLT:NewAtom("!*",0));
    	else
		    SetCompoundArg(prt,7,A_GAUGE);

		if(col>1)
			{
			Term cl;
			if(col==8)
				cl=MakeCompound2(A_COLOR,NewAtom("c8",0),NewAtom("c8",0));
			else
				cl=MakeCompound2(A_COLOR,NewAtom("c3",0),NewAtom("c3b",0));
			SetCompoundArg(prt,8,AppendFirst(NewList(),cl));
			}

		AddIMParticle(prt);
		if(ptname)
			{
			char bbb[64];
			SetAtomProperty(n,A_TEXNAME,ptname);
			sprintf(bbb,"{\\bar{%s}}",AtomValue(ptname));
			SetAtomProperty(an,A_TEXNAME,NewAtom(bbb,0));
			}
		}

	if(spin<2)
		{
		if(!intr_param)
			{
			Term t1,t2,t3;
			t1=NewCompound(NewFunctor(OPR_PARAMETER,1));
			t2=NewCompound(NewFunctor(OPR_EQSIGN,2));
			t3=NewCompound(NewFunctor(OPR_COLON,2));
			SetCompoundArg(t1,1,t2);
			SetCompoundArg(t2,1,NewAtom("Maux",0));
			SetCompoundArg(t2,2,t3);
			SetCompoundArg(t3,1,NewInteger(1));
			SetCompoundArg(t3,2,NewAtom("mass of aux particles",0));
			ProcessParameter(t1,0);
			intr_param=1;
			}
		return n;
		}
	
	/*if(CalcOutput)
	    return n;*/
	sprintf(nnbuf,"%s.t",AtomValue(n));
	n=NewAtom(nnbuf,0);

	return n;
	}
	
