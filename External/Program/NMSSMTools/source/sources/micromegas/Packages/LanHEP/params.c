#include <math.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <dlfcn.h>

#include "lanhep.h"

int has_cmplx_params=0;

int check_funcs=1;

void proc_hash(Term);

typedef double (*extfunc)(double,...);

int P_D_WIDTH=85;
int longest_pdline=0;
Atom llparam=0;

static List params=0, cparams=0, rparams=0;
static Term rm_zero(Term);
static int legal_expr(Term t);
static List ExtFuncList=0;

extern int opTriHeu, FAOutput, doing_abbr;

static int chk_cmplx;

List all_param_list(void)
	{
	return params;
	}

List all_cparam_list(void)
	{
	return cparams;
	}

List all_rparam_list(void)
	{
	return rparams;
	}

static double cu_rt(double r, double i)
{
	double mod, ph;
/*	printf("cu_rt: %g %g; ",r,i);*/
	mod=sqrt(r*r+i*i);
	ph=atan2(i,r);
/*	printf("mod=%g, ph=%f\n",mod,ph);*/
	mod=pow(mod,1.0/3.0);
	ph/=3;
	return mod*cos(ph);
}

static void check_fortr(Atom p)
{
  char cbuf[64];
  int i;
  static int cnt=0;
  Atom fn, prp;
  strcpy(cbuf,AtomValue(p));
  
  for(i=0; cbuf[i]&&i<64; i++)
	cbuf[i]=(char)tolower(cbuf[i]);
  fn=NewAtom(cbuf,0);
  prp=GetAtomProperty(fn,A_FORT_NAME);
  if(prp==0)
  {
	SetAtomProperty(fn,A_FORT_NAME,p);
	return;
  }
  if(prp==p)
	return;
  cnt++;
  if(cnt==5) puts("More parameters which are equivalent in FORTRAN follow.");
  if(cnt>4) return;
  if(FAOutput)
	ErrorInfo(0);
  else
	WarningInfo(0);
  
  printf("parameters %s and %s are equivalent in FORTRAN\n",
		 AtomValue(p),AtomValue(prp));
}

int ufo_gg_injected=0;

static void proc_ufo_ext(Term t)
{
    static Atom sl, sl1=0, sl2=0, ec, eb, et, rc, rb, rt, aq;
    if(sl1==0) 
    {
        sl=NewAtom("slhaVal",0);
        sl1=NewAtom("slhaVal1",0);
        sl2=NewAtom("slhaVal2",0);
        ec=NewAtom("McEff",0);
        eb=NewAtom("MbEff",0);
        et=NewAtom("MtEff",0);
        rc=NewAtom("McRun",0);
        rb=NewAtom("MbRun",0);
        rt=NewAtom("MtRun",0);
        aq=NewAtom("alphaQCD",0);
    }
   // DisplayTerm(t);puts("");
    if(is_list(t))
    {
        List l;
        for(l=t;l;l=ListTail(l)) proc_ufo_ext(ListFirst(l));
        return;
    }
    if(is_compound(t) && CompoundName(t)==OPR_EQSIGN)
    {
       int hasc=is_compound(CompoundArg2(t)) &&
            CompoundName(CompoundArg2(t))==OPR_COLON;
        Term rp=hasc?CompoundArg1(CompoundArg2(t)):CompoundArg2(t);
        if(CompoundArg1(t)==A_GG)
        {
            Atom tn=NewAtom("aS",0);
            SetAtomProperty(tn,PROP_TYPE,OPR_PARAMETER);
            SetAtomProperty(tn,A_LHA,
                MakeCompound1(NewAtom("SMINPUTS",0),MakeList1(NewInteger(3))));
            Term t1=MakeCompound3(tn,NewFloat(0.118),0,0);
            params=AppendLast(params,t1);
            t1=MakeCompound2(OPR_MLT,MakeCompound1(NewAtom("cmath.acos",0),NewInteger(-1)),tn);
            t1=MakeCompound1(NewAtom("cmath.sqrt",0),t1);
            t1=MakeCompound2(OPR_MLT,NewInteger(2),t1);
            if(hasc) SetCompoundArg(CompoundArg2(t),1,t1);
            else SetCompoundArg(t,2,t1);
            ufo_gg_injected=1;
            printf("GG the strong coupling is automatically defined.\n");
            return;
        }
        if(is_compound(rp) && CompoundName(rp)==sl1)
        {
            Atom bl=CompoundArg1(rp);
            if(is_compound(bl)) bl=CompoundArg1(bl);
            Integer i1=CompoundArgN(rp,3);
            SetAtomProperty(CompoundArg1(t),A_LHA,
                MakeCompound1(bl,MakeList1(i1)));
            if(hasc) SetCompoundArg(CompoundArg2(t),1,NewInteger(0));
            else SetCompoundArg(t,2,NewInteger(0));
            return;
        }
        if(is_compound(rp) && CompoundName(rp)==sl && CompoundArity(rp)==4)
        {
            Atom bl=CompoundArg1(rp);
            if(is_compound(bl)) bl=CompoundArg1(bl);
            Integer i1=CompoundArgN(rp,4);
            SetAtomProperty(CompoundArg1(t),A_LHA,
                MakeCompound1(bl,MakeList1(i1)));
            if(hasc) SetCompoundArg(CompoundArg2(t),1,NewInteger(0));
            else SetCompoundArg(t,2,NewInteger(0));
            return;
        }
        if(is_compound(rp) && CompoundName(rp)==sl2)
        {
            Atom bl=CompoundArg1(rp);
            if(is_compound(bl)) bl=CompoundArg1(bl);
            Integer i1=CompoundArgN(rp,3);
            Integer i2=CompoundArgN(rp,4);
            SetAtomProperty(CompoundArg1(t),A_LHA,
                MakeCompound1(bl,MakeList2(i1,i2)));
            if(hasc) SetCompoundArg(CompoundArg2(t),1,NewInteger(0));
            else SetCompoundArg(t,2,NewInteger(0));
            return;
        }
        if(is_compound(rp) && CompoundName(rp)==sl && CompoundArity(rp)==5)
        {
            Atom bl=CompoundArg1(rp);
            if(is_compound(bl)) bl=CompoundArg1(bl);
            Integer i1=CompoundArgN(rp,4);
            Integer i2=CompoundArgN(rp,5);
            SetAtomProperty(CompoundArg1(t),A_LHA,
                MakeCompound1(bl,MakeList2(i1,i2)));
            if(hasc) SetCompoundArg(CompoundArg2(t),1,NewInteger(0));
            else SetCompoundArg(t,2,NewInteger(0));
            return;
        }
        
    }
    if(is_compound(t))
    {
        int i, ar=CompoundArity(t);
        for(i=1;i<=ar;i++)
        {
            Term ta, t1=CompoundArgN(t,i), tn;
            char cbuf[100];
            int cnt=1;
            if(!is_compound(t1)) {proc_ufo_ext(t1); continue;}
            ta=CompoundName(t1);
            if(!GetAtomProperty(ta,A_EXT_FUNC)) {proc_ufo_ext(t1);continue;}
            if(strncmp(AtomValue(ta),"cmath.",6)==0) {proc_ufo_ext(t1);continue;}
            if(strcmp(AtomValue(ta),"abs")==0) {proc_ufo_ext(t1);continue;}
            t1=ConsumeCompoundArg(t,i);
            do {
                sprintf(cbuf,"%s%d",AtomValue(ta),cnt++);
                tn=NewAtom(cbuf,0);
            } while(GetAtomProperty(tn,PROP_TYPE));
            SetCompoundArg(t,i,tn);
            SetAtomProperty(tn,PROP_TYPE,OPR_PARAMETER);
            printf("New parameter %s=",cbuf);WriteTerm(t1);puts("");
            Term nv=NewInteger(1);
            if(ta==ec) nv=NewFloat(0.66065);
            if(ta==eb) nv=NewFloat(3.2407);
            if(ta==et) nv=NewFloat(171.4);
            if(ta==rc) nv=NewFloat(0.59086);
            if(ta==rb) nv=NewFloat(2.8984);
            if(ta==rt) nv=NewFloat(161.91);
            if(ta==aq) nv=NewFloat(0.1156);
           // WriteTerm(ta); printf(" "); WriteTerm(nv);puts("");
            t1=MakeCompound3(tn,nv,0,0);
            params=AppendLast(params,t1);
        }
    }
}

static void proc_param(Term t)
	{
	Atom name, value, comment, aname=0;
	List li;
	
    if(UFOutput)
        proc_ufo_ext(t);
        
	if(is_compound(t) && CompoundName(t)==OPR_COMMA)
		{
		Term a,b;
		a=ConsumeCompoundArg(t,1);
		b=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		proc_param(a);
		proc_param(b);
		return;
		}

	if(is_atom(t) || (is_compound(t) && CompoundName(t)==OPR_DIV) )
		{ name=t; value=0; comment=0; }
	else 
		{
		Term t1;
		name=ConsumeCompoundArg(t,1);
		t1=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		t=t1;
		if(is_compound(t) && CompoundName(t)==OPR_COLON)
			{
			value=ConsumeCompoundArg(t,1);
			comment=ConsumeCompoundArg(t,2);
			FreeAtomic(t);
			}
		else
			{ value=t; comment=0; }
		}

	if(is_compound(name) && CompoundName(name)==OPR_DIV)
		{
/*		if(!FAOutput)
			{
			ErrorInfo(0);
			printf("complex parameter is not allowed.\n");
			}
*/
		aname=ConsumeCompoundArg(name,2);
		t=ConsumeCompoundArg(name,1);
		FreeAtomic(name);
		name=t;
		}

	if(!is_atom(name))
		{
		ErrorInfo(101);
		printf(" \'");
		WriteTerm(name);
		printf("\' is not appropriate name for parameter\n");
		FreeAtomic(name);
		FreeAtomic(value);
		FreeAtomic(comment);
		return;
		}

	if(aname && !is_atom(aname))
		{
		ErrorInfo(101);
		printf(" \'");
		WriteTerm(aname);
		printf("\' is not appropriate name for parameter\n");
		FreeAtomic(name);
		FreeAtomic(value);
		FreeAtomic(comment);
		return;
		}

	if(comment!=0 && !is_atom(comment))
		{
		ErrorInfo(102);
		printf(" \'");
		WriteTerm(comment);
		printf("\' is not appropriate comment for parameter \'%s\'\n",
				AtomValue(name));
		FreeAtomic(name);
		FreeAtomic(value);
		FreeAtomic(comment);
		return;
		}
	
	
	if(!doing_abbr)
	{
		value=WheredTerm(value);
		proc_hash(value);
		
		/*WriteTerm(value);printf(" -> ");*/
		value=rm_zero(value);
		/*WriteTerm(value);puts("");*/
	}

	

	chk_cmplx=0;

	if(is_compound(value)&&is_function(value,0)&&!
			GetAtomProperty(CompoundName(value),A_EXT_FUNC))
		value=CallFunction(value,0);
	
	if(/*doing_abbr==0 &&*/ value !=0 && !legal_expr(value))
		{
		ErrorInfo(103);
		printf(" \'");
		WriteTerm(value);
		printf("\' is not appropriate value for parameter \'%s\'\n",
				AtomValue(name));
		FreeAtomic(name);
		value=0;
                comment=0;
		}

	if(chk_cmplx) has_cmplx_params=1;
		
	if(chk_cmplx && aname==0)
		{
		char cbuf[80];
		if(!FAOutput && !CalcOutput)
			{
			ErrorInfo(0);
			printf("complex parameter is not allowed.\n");
			}

		sprintf(cbuf,"%scc",AtomValue(name));
		aname=NewAtom(cbuf,0);
		}

	if(AtomValue(name)[0]=='%' && value==0)
		value=A_I;

	if(!doing_abbr)
	{
	li=params;
	while(!is_empty_list(li))
		{
		Term t1;
		t1=ListFirst(li);
		if(FunctorName(CompoundFunctor(t1))==name)
			{
			Term ov=CopyTerm(CompoundArg1(t1));
			if(ov==0) ov=NewInteger(0);
			if(value!=0)
				SetCompoundArg(t1,1,value);
			if(comment!=0)
				SetCompoundArg(t1,2,comment);

			if(strcmp(AtomValue(name),"Maux"))
				{
				WarningInfo(104);
				printf("Changing parameter \'%s\':",AtomValue(name));
				printf(" it was ");	WriteTerm(ov);
				printf(", now "); WriteTerm(CompoundArg1(t1));
				puts("");
				}
			return;
			}
		li=ListTail(li);
		}
	 check_fortr(name);
	 if(aname) check_fortr(aname);
	}

/*	fu=NewFunctor(name,3);*/
	t=MakeCompound3(name,value,comment,0);
/*	SetCompoundArg(t,1,value);
	SetCompoundArg(t,2,comment);*/
	params=AppendLast(params,t);
	if(aname)
		cparams=AppendLast(cparams,name);
	else
	 	rparams=AppendLast(rparams,name);
	ReportRedefined(name,"parameter");
	SetAtomProperty(name,PROP_TYPE,OPR_PARAMETER);

	if(opTriHeu && !doing_abbr)
		tri_reg_prm(name,value);

	if(aname)
	{
	ReportRedefined(aname,"parameter");
	SetAtomProperty(aname,PROP_TYPE,OPR_PARAMETER);
	SetAtomProperty(name,A_ANTI,aname);
	SetAtomProperty(aname,A_ANTI,name);
	SetAtomProperty(aname,A_HERMC,name);
	}

/*	if(IsTermInput())
		{
		printf("Parameter \'%s\' ",AtomValue(name));
		if(value!=0)
			{ printf("   value= ");		WriteTerm(value);}
		if(comment!=0)
			{ printf("   comment= ");	WriteTerm(comment);}
		puts("");
		}*/
	}

static List params_tail=0, rparams_tail=0, cparams_tail=0;

void proc_abbr_param(Term name, Term value)
	{
	Atom comment=0, aname=0;
	Term t;
	if(params_tail==0)
		params_tail=params;
	while(ListTail(params_tail)) params_tail=ListTail(params_tail);

	if(rparams_tail==0)
		rparams_tail=rparams;
	while(ListTail(rparams_tail)) rparams_tail=ListTail(rparams_tail);

	if(cparams_tail==0)
		cparams_tail=cparams;
	while(cparams_tail && ListTail(cparams_tail)) cparams_tail=ListTail(cparams_tail);

	chk_cmplx=0;
	legal_expr(value);
	if(chk_cmplx)
	{
		char cbuf[80];
		Atom aname;
		sprintf(cbuf,"%scc",AtomValue(name));
		aname=NewAtom(cbuf,0);
		SetAtomProperty(aname,PROP_TYPE,OPR_PARAMETER);
		SetAtomProperty(name,A_ANTI,aname);
		SetAtomProperty(aname,A_ANTI,name);
		SetAtomProperty(aname,A_HERMC,name);
		if(cparams==0)
			cparams=AppendLast(cparams,name);
		else
			AppendLast(cparams_tail,name);
	}
	else
		AppendLast(rparams_tail,name);
	t=MakeCompound3(name,value,comment,0);
	AppendLast(params_tail,t);
	SetAtomProperty(name,PROP_TYPE,OPR_PARAMETER);

	}

Term ProcessParameter(Term t, Term ind)
	{
	Term arg;
	char *s;
	if(is_atom(t))
		arg=NewAtom("?",0);
	else
		arg=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	if(is_atom(arg))
		{
		s=AtomValue(arg);
		if(s[0]=='?' && s[1]==0)
			{
			List li;
			li=params;
			while(!is_empty_list(li))
				{
				Term t1,t2;
				t1=ListFirst(li);
				printf("%s",AtomValue(FunctorName(CompoundFunctor(t1))));
				t2=CompoundArg1(t1);
				if(t2!=0)
					{ printf(" = "); WriteTerm(t2); }
				t2=CompoundArg2(t1);
				if(t2!=0)
					{ printf(" : "); WriteTerm(t2); }
				if(CompoundArgN(t1,3))
					{ printf(" ( = %f )",FloatValue(CompoundArgN(t1,3))); }
				puts("");
				li=ListTail(li);
				}
			return 0;
			}
		}
	proc_param(arg);
	return 0;
	}

int is_parameter(Atom p)
	{
	if(p==A_SQRT2)
		return 1;
	if(GetAtomProperty(p,A_INFINITESIMAL))
		return 1;
    return (GetAtomProperty(p,PROP_TYPE)==OPR_PARAMETER);
    /*
	li=params;
	while(!is_empty_list(li))
		{
		if(FunctorName(CompoundFunctor(ListFirst(li)))==p)
			return 1;
		li=ListTail(li);
		}
	return 0;
    */
	}


static int legal_expr(Term t)
	{
	char c;
	Functor fu;
	int arity,i,ret;

	ret=1;
	c=AtomicType(t);
	switch(c)
		{
		case 'i':
		case 'f':
			return 1;
		case 'a':
			if(t==A_I)
				{
				chk_cmplx=1;
				return 1;
				}
			if(is_parameter(t))
				{
				if(GetAtomProperty(t,A_ANTI))
					chk_cmplx=1;
				return 1;
				}
			printf("Error: symbol \'%s\' was not declared.\n",AtomValue(t));
			return 0;
		case 'c':

			arity=CompoundArity(t);
			for(i=1;i<=arity;i++)
			{
				Term t1=CompoundArgN(t,i);
				if(is_compound(t1)&&is_function(t1,0)&&
						!GetAtomProperty(CompoundName(t1),A_EXT_FUNC))
				{
					t1=ConsumeCompoundArg(t,i);
					t1=CallFunction(t1,0);
					SetCompoundArg(t,i,t1);
				}
			}
			
			fu=CompoundFunctor(t);
			if( (CompoundArity(t)==2 && CompoundName(t)==OPR_PLUS)  ||
				(CompoundArity(t)==1 && CompoundName(t)==OPR_PLUS)  ||
				(CompoundArity(t)==2 && CompoundName(t)==OPR_MINUS) ||
				(CompoundArity(t)==1 && CompoundName(t)==OPR_MINUS) ||
				(CompoundArity(t)==2 && CompoundName(t)==OPR_DIV)   ||
				(CompoundArity(t)==2 && CompoundName(t)==OPR_MLT)   ||
				(CompoundArity(t)==2 && CompoundName(t)==OPR_POW)   ||
				(CompoundArity(t)==1 && CompoundName(t)==A_SIN)     ||
				(CompoundArity(t)==1 && CompoundName(t)==A_COS)     ||
				(CompoundArity(t)==1 && CompoundName(t)==A_ASIN)    ||
				(CompoundArity(t)==1 && CompoundName(t)==A_ACOS)    ||
				(CompoundArity(t)==1 && CompoundName(t)==A_TAN)     ||
				(CompoundArity(t)==1 && CompoundName(t)==A_ATAN)    ||
				(CompoundArity(t)==1 && CompoundName(t)==A_EXP)     ||
				(CompoundArity(t)==2 && CompoundName(t)==A_ATAN2)   ||
				(CompoundArity(t)==3 && CompoundName(t)==A_IF)      ||
				(CompoundArity(t)==1 && CompoundName(t)==A_LOG)     ||
				(CompoundArity(t)==1 && CompoundName(t)==A_SQRT)    ||
				(CompoundArity(t)==1 && CompoundName(t)==A_FABS)    ||
				(GetAtomProperty(CompoundName(t),A_EXT_FUNC))|| !check_funcs)
			{
			Integer efa=GetAtomProperty(CompoundName(t),A_EXT_FUNC);
			int cmplx_save=chk_cmplx;
			arity=FunctorArity(fu);

			if(efa && IntegerValue(efa)!=0  && IntegerValue(efa)!=arity)
				{
				printf("Wrong arguments number in function %s(...).\n",
					AtomValue(CompoundName(t)));
				return 0;
				}

			for(i=1;i<=arity;i++)
				{
				Term ca=CompoundArgN(t,i);
				if(is_compound(ca) && CompoundArity(ca)==1 &&
					is_atom(CompoundArg1(ca)) &&
					strcmp(AtomValue(CompoundName(ca)),"str")==0)
					{
					char sbuf[128];
					Atom va,qa;
					va=CompoundArg1(ca);
					if(FAOutput==0)
						sprintf(sbuf,"\"%s\"",AtomValue(va));
					else
						sprintf(sbuf,"'%s'",AtomValue(va));
					SetCompoundArg(t,i,(qa=NewAtom(sbuf,0)));
					SetAtomProperty(qa,A_I,va);
					continue;
					}
				if(!check_funcs && is_atom(ca) && !is_parameter(ca))
				{
					char sbuf[128];
					int j;
					Atom va,qa;
					va=ca;
					if(FAOutput==0)
						sprintf(sbuf,"\"%s\"",AtomValue(va));
					else
						sprintf(sbuf,"'%s'",AtomValue(va));
					for(j=0;sbuf[j];j++)
						if(sbuf[j]=='#') sbuf[j]='%';
					SetCompoundArg(t,i,(qa=NewAtom(sbuf,0)));
					SetAtomProperty(qa,A_I,va);
					continue;
				}
				if(!legal_expr(CompoundArgN(t,i)))
					ret=0;
				}
			if(efa && GetAtomProperty(CompoundName(t),A_ANTI)==0
					&& cmplx_save==0)
				chk_cmplx=0;
			return ret;
			}
			else
			{
                                ErrorInfo(0);
				printf("unknown function %s(). Use statement:\n",AtomValue(CompoundName(t)));
                                printf("external_func(%s,%d).\n",AtomValue(CompoundName(t)),CompoundArity(t));
				return 0;
			}
		default:
			printf("Internal error: bad expression '%c' ",c);
			DisplayTerm(t);
			puts("");
			return 0;
		}
	return 0;
	}


static char wpbuf[128];

static void texWriteParams(int fno, char *name)
	{
	FILE *f1;
	List li;
	if(OutputDirectory!=NULL)
		sprintf(wpbuf,"%s/vars%d.tex",OutputDirectory,fno);
	else
		sprintf(wpbuf,"vars%d.tex",fno);
	f1=fopen(wpbuf,"w");
	if(f1==NULL)
		{
		printf("Can not open file \'%s\' for writing.\n",wpbuf);
		perror("");
		return;
		}

	fprintf(f1,"\\begin{tabular}{|l|l|l|} \\hline\n");
	fprintf(f1,"Parameter & Value & Comment \\\\ \\hline\n");
	li=params;
	while(!is_empty_list(li))
		{
		Term ttt, val, comm;
		ttt=ListFirst(li);
		val=CompoundArg1(ttt);
		comm=CompoundArg2(ttt);
		if(val==0) val=NewInteger(0);
		if(is_integer(val) || is_float(val))
			{
			int sp;
			sp=fprintf(f1,"%s",AtomValue(FunctorName(CompoundFunctor(ttt))));
			WriteBlank(f1,6-sp);
			fprintf(f1,"&");
			sp=fWriteTerm(f1,val);
			WriteBlank(f1,20-sp);
			fprintf(f1,"&");
			if(comm!=0)
				{
				if(ListTail(li))
					fprintf(f1,"%s \\\\\n",AtomValue(comm));
				else
					fprintf(f1,"%s ",AtomValue(comm));
				}
			else
				{
				if(ListTail(li))
					fprintf(f1," \\\\\n");
				else
					fprintf(f1," ");
				}
			}
		else
			{
			int sp;
			sp=fprintf(f1,"%s",AtomValue(FunctorName(CompoundFunctor(ttt))));
			WriteBlank(f1,6-sp);
			fprintf(f1,"&");
			sp=fWriteTerm(f1,val);
			WriteBlank(f1,20-sp);
			fprintf(f1,"&");
			if(comm!=0)
				{
				if(ListTail(li))
					fprintf(f1,"%s \\\\\n",AtomValue(comm));
				else
					fprintf(f1,"%s ",AtomValue(comm));
				}
			else
				{
				if(ListTail(li))
					fprintf(f1," \\\\\n");
				else
					fprintf(f1," ");
				}
			}
		li=ListTail(li);
		}
	fprintf(f1,"\\\\ \\hline\n\\end{tabular}\n");
	fclose(f1);
	}

extern int EvalPrm;

static void rpl_pow(Term e)
{
	int a,i;
	if(!is_compound(e))
		return;
	a=CompoundArity(e);
	for(i=1;i<=a;i++)
		rpl_pow(CompoundArgN(e,i));
	if(CompoundName(e)==OPR_POW)
		SetCompoundName(e,OPR_CARET);
}

static int has_Q(Term v)
{
	int i;
	int r=0;
	if(is_atom(v) && AtomValue(v)[0]=='Q' && AtomValue(v)[1]==0)
		return 1;
	if(!is_compound(v))
		return 0;
	for(i=1;i<=CompoundArity(v);i++)
		r|=has_Q(CompoundArgN(v,i));
	return r;
}

static int micro_allpr(Atom p, Term val)
{
	char *nm=AtomValue(p);
	if(nm[0]=='B'&&nm[1]=='0')
		return 0;
	if(is_integer(val)||is_float(val))
		return 1;
	return !has_Q(val);
}

List cls_real_matr(void);
extern int opAbbArr, r_abbr_no, c_abbr_no;

void FADeclRealParam(FILE *f)
{
	List li,li2;
	int i;
	fprintf(f,"Scan[ (RealQ[#] = True)&,\n  { "); 
	for(li=rparams,i=1;li;li=ListTail(li),i++)
	{
		if(opAbbArr && GetAtomProperty(ListFirst(li),A_CHNAME))
		{
			fprintf(f,"AAABR");
			break;
		}
		fprintf(f,"%s%c ",AtomValue(ListFirst(li)),
				ListTail(li)?',':' ');
		if(i%10==0) fprintf(f,"\n    ");
	}
	fprintf(f,"} ]\n\n");
	li2=cls_real_matr();
	if(li2==0)
		return;

	fprintf(f,"Scan[ (RealQ[#] = True)&,\n  { "); 
	for(li=li2,i=1;li;li=ListTail(li),i++)
	{
		fprintf(f,"%s%c ",AtomValue(ListFirst(li)),
				ListTail(li)?',':' ');
		if(i%8==0) fprintf(f,"\n    ");
	}
	fprintf(f,"} ]\n\n");
	FreeAtomic(li2);
}

void cls_write_decl(FILE *);
void cls_write_matr(FILE *);
extern List fainclude;
extern Atom fa_i;
Atom fa_ccf=0, fa_rabb=0, fa_cabb=0;


Atom toccf(Atom a)
{
	Atom aa;
	Term  no;
	if(a==A_I)
		return fa_i;
	if(opAbbArr && (no=GetAtomProperty(a,A_CHNAME)))
	{
		no=IntegerValue(no);
		if(GetAtomProperty(a,A_ANTI))
		{
			aa=GetAtomProperty(a,A_HERMC);
			if(aa)
				return MakeCompound1(fa_ccf,
						MakeCompound1(fa_cabb,NewInteger(no)));
			else
				return MakeCompound1(fa_cabb,NewInteger(no));
		}
		else
			return MakeCompound1(fa_rabb,NewInteger(no));
	}
	if((aa=GetAtomProperty(a,A_HERMC)))
		return MakeCompound1(fa_ccf,aa);
	return a;
}

void ccf_term(Term t)
{
	if(is_compound(t))
	{
		int i, ar=CompoundArity(t);
		for(i=1;i<=ar;i++)
		{
			Term t1=CompoundArgN(t,i);
			if(is_atom(t1))
			{
				Atom a1=toccf(t1);
				if(a1!=t1)
					SetCompoundArg(t,i,a1);
				continue;
			}
			if(is_list(t1) || is_compound(t1))
				ccf_term(t1);
		}
		return;
	}
	if(is_list(t))
	{
		List l;
		for(l=t;l;l=ListTail(l))
		{
			Term t1=ListFirst(l);
			if(is_list(t1) || is_compound(t1))
				ccf_term(t1);
		}
	}
}

static void repl_pow(Term t)
{
	int i,ar;
	if(!is_compound(t))
		return;
	ar=CompoundArity(t);
	for(i=1;i<=ar;i++)
	{
		if(!is_compound(CompoundArgN(t,i)))
			continue;
		if(CompoundArity(CompoundArgN(t,i))==2 &&
			CompoundName(CompoundArgN(t,i))==OPR_POW)
			{
			Term t1=ConsumeCompoundArg(t,i);
			Atom n=CompoundArg1(t1);
            if(!is_integer(CompoundArg2(t1)))
            {
                Term n=MakeCompound1(NewAtom("math.log",0),ConsumeCompoundArg(t1,1));
                n=MakeCompound2(OPR_MLT,ConsumeCompoundArg(t1,2),n);
                n=MakeCompound1(NewAtom("math.exp",0),n);
                SetCompoundArg(t,i,n);
                continue;
            }
			int j, pw=(int)IntegerValue(CompoundArg2(t1));
			if(pw<2 || pw>15)
			{
			puts("Internal error (rplpowuf)");
			return;
			}
			t1=MakeCompound2(OPR_MLT,n,n);
			for(j=3;j<=pw;j++)
				t1=MakeCompound2(OPR_MLT,t1,n);
			SetCompoundArg(t,i,t1);
			}
		else
			repl_pow(CompoundArgN(t,i));
	}
}

extern int UFOutput, ufo_gg_injected;
extern List UFOparah,UFOcouph;
extern char *eff_infile;
List ufo_params_nolha=0;
static int ufo_params_nolha_no=1;

void UFWriteParameters(void)
{
	char cbuf[128];
	time_t tm;
	List li,lj;
	FILE *f;

	NoQuotes=1;
    
    if(!ufo_gg_injected)
        puts("Warning: strong coupling was not automatically defined.");
    
	if(OutputDirectory!=NULL)
		sprintf(cbuf,"%s/parameters.py",OutputDirectory);
	else
		sprintf(cbuf,"parameters.py");
	f=fopen(cbuf,"w");

	fprintf(f,"#     LanHEP output produced at ");
	time(&tm);
	fprintf(f,"%s",ctime(&tm));
	fprintf(f,"#     from the file '%s'\n",eff_infile);

	if(ModelName)
		fprintf(f,"#     Model named '%s'\n",ModelName);
	fprintf(f,"\n");

	
	for(li=UFOparah;li;li=ListTail(li))
		fprintf(f,"%s\n",AtomValue(ListFirst(li)));
	fprintf(f,"\n");
    
    fprintf(f,"ZERO = Parameter(name = 'ZERO',\n");
	fprintf(f,"                  nature = 'internal',\n");
	fprintf(f,"                  type = 'real',\n");
	fprintf(f,"                  value = '0.0',\n");
	fprintf(f,"                  texname = '0')\n\n");
    
	fprintf(f,"Sqrt2 = Parameter(name = 'Sqrt2',\n");
	fprintf(f,"                  nature = 'internal',\n");
	fprintf(f,"                  type = 'real',\n");
	fprintf(f,"                  value = 'cmath.sqrt(2.0)',\n");
	fprintf(f,"                  texname = '\\\\sqrt{2}')\n\n");
	

	for(li=params;li;li=ListTail(li))
	{
		Term t, val, tm;
		char mt[64];
		t=CompoundName(ListFirst(li));
/*		if(strcmp(AtomValue(t),"pi")==0)
		  continue;*/
		if(strcmp(AtomValue(t),"A00001")==0)
		  break;
		val=CopyTerm(CompoundArg1(ListFirst(li)));

		if(val==0) val=NewInteger(0);
		fprintf(f,"%s = Parameter(name = '%s',\n",AtomValue(t),AtomValue(t));
		fprintf(f,"               nature = '%s',\n",
				(is_integer(val)||is_float(val))?"external":"internal");
		fprintf(f,"               type = '%s',\n",
			GetAtomProperty(t,A_ANTI)?"complex":"real");
		fprintf(f,"               value = ");
		if(!is_integer(val) && !is_float(val))
			fprintf(f,"'");
		fortr_digi=1;
		repl_pow(val);
		fWriteTerm(f,val);
		if(!is_integer(val) && !is_float(val))
			fprintf(f,"'");
		/*FreeAtomic(val);*/
		fprintf(f,",\n");
		fortr_digi=0;
		tm=GetAtomProperty(t,A_TEXNAME);
		if(tm==0)
		{
			sprintf(mt,"\\\\text{%s}",AtomValue(t));
			tm=NewAtom(mt,0);
		}
		fprintf(f,"               texname = '%s'",AtomValue(tm));
		tm=GetAtomProperty(t,A_LHA);
		if(tm)
		{
			fprintf(f,",\n");
		fprintf(f,"               lhablock = '%s',\n",
						AtomValue(CompoundName(tm)));
		fprintf(f,"               lhacode = ");
		fWriteTerm(f,CompoundArg1(tm));
		}
		else if((is_integer(val)||is_float(val)))
		{
		  ufo_params_nolha=AppendLast(ufo_params_nolha,t);
		  fprintf(f,",\n               lhablock = 'USERDEF',\n");
		  fprintf(f,"               lhacode = [ %d ]",ufo_params_nolha_no++);
		}
		fprintf(f," )\n\n");

	}

	fclose(f);


	if(OutputDirectory!=NULL)
		sprintf(cbuf,"%s/couplings.py",OutputDirectory);
	else
		sprintf(cbuf,"couplings.py");
	f=fopen(cbuf,"w");


	fprintf(f,"#     LanHEP output produced at ");
	time(&tm);
	fprintf(f,"%s",ctime(&tm));
	fprintf(f,"#     from the file '%s'\n",eff_infile);

	if(ModelName)
		fprintf(f,"#     Model named '%s'\n",ModelName);
	fprintf(f,"\n");

	
	for(lj=UFOcouph;lj;lj=ListTail(lj))
		fprintf(f,"%s\n",AtomValue(ListFirst(lj)));
	fprintf(f,"\n");


	for(;li;li=ListTail(li))
	{
		Term t, val, tm, tm2;
		char mt[64];
		int no;
		t=CompoundName(ListFirst(li));
/*		if(strcmp(AtomValue(t),"pi")==0)
		  continue;*/

		sscanf(AtomValue(t)+1,"%d",&no);

		val=CopyTerm(CompoundArg1(ListFirst(li)));

		if(val==0) val=NewInteger(0);
		fprintf(f,"GC_%d = Coupling(name = 'GC_%d',\n",no,no);
		fprintf(f,"               value = ");
		fortr_digi=1;
		if(!is_integer(val) && !is_float(val))
			fprintf(f,"'");
		
		repl_pow(val);
		
		fWriteTerm(f,val);
		if(!is_integer(val) && !is_float(val))
			fprintf(f,"'");
		fortr_digi=1;
		fprintf(f,",\n");
		tm=GetAtomProperty(t,A_EE);
		tm2=GetAtomProperty(t,A_GG);

		fprintf(f,"               order = {");
		if(tm)
			fprintf(f,"'QED':%ld",IntegerValue(tm));
		if(tm&&tm2)
			fprintf(f,",");
		if(tm2)
			fprintf(f,"'QCD':%ld",IntegerValue(tm2));

		fprintf(f," })\n\n");

	}

	fclose(f);
	NoQuotes=0;

}

Term ProcLHA(Term t, Term ind)
{
	List l, l1;
	if(CompoundArity(t)!=1)
	{
		ErrorInfo(0);
		puts("lha: wrong syntax.");
		return 0;
	}

	l=CommaToList(CompoundArg1(t));
	for(l1=l;l1;l1=ListTail(l1))
	{
		Term t=ListFirst(l1), v, v1;
		if(!is_compound(t) || CompoundArity(t)!=2)
		{
			ErrorInfo(0);
			printf("lha: wrong construct '");
			WriteTerm(t);
			puts("'.");
			continue;
		}
		if(!is_parameter(CompoundArg1(t)))
		{
			WarningInfo(0);
			printf("lha:  '");
			WriteTerm(CompoundArg1(t));
			puts("' is not a parameter.");
			continue;
		}
		v=CompoundArg2(t);
		if(!is_compound(v) || CompoundArity(t)>2 || 
			!is_integer(CompoundArg1(v)) ||
			(CompoundArity(v)==2 && !is_integer(CompoundArg2(v))))
		{
			ErrorInfo(0);
			printf("lha: wrong block/code '");
			WriteTerm(v);
			puts("'.");
			continue;
		}
		if(CompoundArity(v)==1)
			v1=MakeCompound1(CompoundName(v),MakeList1(CompoundArg1(v)));
		else
			v1=MakeCompound1(CompoundName(v),MakeList2(CompoundArg1(v),
												CompoundArg2(v)));
		SetAtomProperty(CompoundArg1(t),A_LHA,v1);
	}
	FreeAtomic(l);
	return 0;
}

extern int FAver, FindVal;
extern List mtr_visi;

void FAWriteParameters(int fno)
{
	char cbuf[128];
	time_t tm;
	List li;
	FILE *f;
	int p;
	int has_GG=0;
	int mtrini=0;
	int first=1;
	int aano=-1, abno=0;

	if(UFOutput)
	{
		UFWriteParameters();
		return;
	}

	if(fa_i==0)
		fa_i=NewAtom("cI",0);
	if(fa_ccf==0)
	{
		fa_ccf=NewAtom("dconjg",0);
		fa_rabb=NewAtom("AAABR",0);
		fa_cabb=NewAtom("AAABC",0);
		SetAtomProperty(fa_rabb,A_CHNAME,A_I);
		SetAtomProperty(fa_cabb,A_CHNAME,A_I);
	}

	if(OutputDirectory!=NULL)
		sprintf(cbuf,"%s/model%d.h",OutputDirectory,fno);
	else
		sprintf(cbuf,"model%d.h",fno);
	f=fopen(cbuf,"w");


	fprintf(f,"*     LanHEP output produced at ");
	time(&tm);
	fprintf(f,"%s",ctime(&tm));
	if(ModelName)
		fprintf(f,"*     Model named '%s'\n",ModelName);
	fprintf(f,"\n");

	if(FAver<8)
	{
	fprintf(f,"      double precision Sqrt2, pi, degree, hbar_c2,bogus\n");
	fprintf(f,"      parameter (Sqrt2=1.41421356237309504880168872421D0)\n");
	fprintf(f,"      parameter (pi = 3.1415926535897932384626433832795029D0)\n");
	fprintf(f,"      parameter (degree = pi/180D0)\n");
	fprintf(f,"      parameter (hbar_c2 = 3.8937966D8)\n");
	fprintf(f,"      parameter (bogus = -1D123)\n");
	fprintf(f,"      double complex cI\n      parameter (cI = (0D0, 1D0))\n\n");

	fprintf(f,"      double precision Divergence\n      common /renorm/ Divergence\n\n");
	}

	p=fprintf(f,"      double precision ");

	for(li=rparams;li;li=ListTail(li))
	{
		Term ttt=ListFirst(li);
		if(strcmp(AtomValue(ttt),"pi")==0)
		  continue;
		if(ttt==A_GG) has_GG=1;
		if(opAbbArr && GetAtomProperty(ttt,A_CHNAME))
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
		p+=fprintf(f,"%s",AtomValue(ttt));
	}
	if(!has_GG)
	{
		if(!first) p+=fprintf(f,", ");
		fprintf(f,"GG");
	}
	if(cparams)
	{
		p=fprintf(f,"\n\n      double complex ")-2;
		first=1;
	}

	for(li=cparams;li;li=ListTail(li))
	{
		Term ttt=ListFirst(li);
		if(opAbbArr && GetAtomProperty(ttt,A_CHNAME))
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
		p+=fprintf(f,"%s",AtomValue(ttt));
	}
	fprintf(f,"\n\n");

	if(opAbbArr)
	{
		if(r_abbr_no)
			fprintf(f,"      double precision AAABR(%d)\n",r_abbr_no);
		if(c_abbr_no)
			fprintf(f,"      double complex AAABC(%d)\n",c_abbr_no);
	}
	fprintf(f,"\n");

	cls_write_decl(f);

	first=1;
	fprintf(f,"      common /mdl_para/\n     &    ");
	p=10;

	for(li=params;li;li=ListTail(li))
	{
		Term ttt=CompoundName(ListFirst(li));
		if(strcmp(AtomValue(ttt),"pi")==0)
		  continue;
		if(opAbbArr && GetAtomProperty(ttt,A_CHNAME))
			continue;
		if(p>60)
		{
			fprintf(f,",\n");
			p=fprintf(f,"     &    ");
			first=1;
		}
		if(!first)  p+=fprintf(f,", ");
		p+=fprintf(f,"%s",AtomValue(ttt));
		first=0;
	}
	if(!has_GG)
	{
		fprintf(f,", GG");
	}
	if(opAbbArr && r_abbr_no)
		p+=fprintf(f,", AAABR");
	if(opAbbArr && c_abbr_no)
		p+=fprintf(f,", AAABC");

	fprintf(f,"\n\n");
	fclose(f);

	NoQuotes=1;

	if(OutputDirectory!=NULL)
		sprintf(cbuf,"%s/mdl_ini%d.F",OutputDirectory,fno);
	else
		sprintf(cbuf,"mdl_ini%d.F",fno);
	f=fopen(cbuf,"w");


	if(EvalPrm)
	{
		for(li=params;li;li=ListTail(li))
		{
			Term t=CompoundArg1(ListFirst(li));
			if(t && !is_integer(t) && !is_float(t))
				SetCompoundArg(ListFirst(li),1,NewFloat(EvalParameter(t)));
		}
	}

	fprintf(f,"*     LanHEP output produced at ");
	time(&tm);
	fprintf(f,"%s",ctime(&tm));
	if(ModelName)
		fprintf(f,"*     Model named '%s'\n",ModelName);

	if(FAver<8)
		fprintf(f,"\n      subroutine ModelDefaults\n      implicit none\n\n");
	else
	{
		fprintf(f,"\n      subroutine ModelDefaults(argc,argv)\n");
		fprintf(f,"\n      implicit none\n");
		fprintf(f,"\n      integer argc\n");
		fprintf(f,"\n      character*128 argv(*)\n\n");
	}	
			
	if(FAver<8)
		fprintf(f,"#include \"model.h\"\n\n");
	else
		fprintf(f,"#include \"decl.h\"\n\n");
				
	for(li=params;li;li=ListTail(li))
	{
		Term t, val;
		char mt[64];
		t=CompoundName(ListFirst(li));
		if(strcmp(AtomValue(t),"pi")==0)
		  continue;
		val=CompoundArg1(ListFirst(li));
		if(val==0) val=NewInteger(0);
		if(is_float(val) || is_integer(val))
		{
		  if(FindVal)
		    fprintf(f,"      %s = findValW('%s')\n",AtomValue(t),AtomValue(t));
		  else
		  {
		  fortr_digi=1;
		  sWriteTerm(mt,val);
		  fortr_digi=0;
			fprintf(f,"      %s = %s\n",AtomValue(t),mt);
		  }
		}
	}

    fprintf(f,"\n      end\n\n");
    
    int firstll=1;
    fprintf(f,"      subroutine ModelVarFromFile(nFile)\n      implicit none\n");
    fprintf(f,"      double precision sqrtS\n      integer nFile\n");
    fprintf(f,"      double precision Alfas\n\n#include \"decl.h\"\n");
    fprintf(f,"      character*10 name\n       real*8  val\n\n");
    fprintf(f,"      character*128 argv(1)\n");
    fprintf(f,"      call ModelDefaults(0,argv)\n\n123   continue\n");
    fprintf(f,"      read(nFile,*,end=321) name, val\n");
    for(li=params;li;li=ListTail(li))
	{
		Term t, val;
		char mt[64];
		t=CompoundName(ListFirst(li));
		if(strcmp(AtomValue(t),"pi")==0)
		  continue;
        fprintf(f,"      %s",firstll?"":"else ");
        firstll=0;
        fprintf(f,"if(name.eq.\"%s\") then\n",AtomValue(t));
        fprintf(f,"         %s=val\n",AtomValue(t));
	}
    fprintf(f,"      endif\n      goto 123\n321   continue\n      end\n\n");
    
	if(FAver==4)
	fprintf(f,"\n      subroutine ModelConstIni(*)\n      implicit none\n\n");
	else
	{
		fprintf(f,"\n      subroutine ModelConstIni(fail)\n      implicit none\n");
		fprintf(f,"      integer fail\n\n");
	}
	if(FAver<8)
		fprintf(f,"#include \"model.h\"\n\n");
	else
		fprintf(f,"#include \"decl.h\"\n\n");

	if(ExtFuncList)
	{
		int g=0;
		fprintf(f,"      double precision ");
		for(li=ExtFuncList;li;li=ListTail(li),g++)
		{
			if(g%5==0) fprintf(f,"\n     > ");
			fprintf(f,"%s%c",AtomValue(ListFirst(li)),ListTail(li)?',':'\n');
		}
		fprintf(f,"\n");
	}

	if(FAver>4)
		fprintf(f,"      fail=0\n");

	if(FAver==4)
		fprintf(f,"      call ModelDefaults\n");

	for(li=params;li;li=ListTail(li))
	{
		Term t, val;
		char mt[102400], *mtp;
		t=CompoundName(ListFirst(li));
		if(strcmp(AtomValue(t),"pi")==0)
		  continue;
		val=CompoundArg1(ListFirst(li));
		if(val==0) val=NewInteger(0);
		if(is_float(val) || is_integer(val)) continue;

		if(strcmp(AtomValue(t),"A00001")==0)
		{
			int bno=ListLength(li)-1, j;
			fprintf(f,"\n");
			bno=(bno/1000)+1;
			for(j=1;j<=bno;j++)
				fprintf(f,"      call aaini%02d\n",j);
			fprintf(f,"      call mtrini\n      end\n\n");
			fprintf(f,"      subroutine aaini01\n      implicit none\n");
	if(FAver<8)
		fprintf(f,"#include \"model.h\"\n\n");
	else
		fprintf(f,"#include \"decl.h\"\n\n");
			abno=1;
			aano=0;
			mtrini=1;
		}

		if(aano==1000)
		{
			fprintf(f,"      end\n\n");
			fprintf(f,"      subroutine aaini%02d\n      implicit none\n",++abno);
	if(FAver<8)
		fprintf(f,"#include \"model.h\"\n\n");
	else
		fprintf(f,"#include \"decl.h\"\n\n");
			aano=0;
		}

		if(aano!=-1)
			aano++;

		if(opAbbArr && GetAtomProperty(t,A_CHNAME))
		{
			int no=(int)IntegerValue(GetAtomProperty(t,A_CHNAME));
			if(mtr_visi && GetAtomProperty(t,OPR_LOCAL)==0)
			  continue;
			p=fprintf(f,"      %s(%d) = ",
				GetAtomProperty(t,A_ANTI)?"AAABC":"AAABR",no);
		}
		else
		p=fprintf(f,"      %s = ",AtomValue(t));

		if(is_atom(val))
			val=toccf(val);
		else
			{
			val=CopyTerm(val);
			ccf_term(val);
			}
		fortr_digi=1;
		sWriteTerm(mt,val);
		/*WriteTerm(val);puts("\n");*/
		fortr_digi=0;
		FreeAtomic(val);
		mtp=mt;
		while(p+(int)strlen(mtp)>60)
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

	if(mtrini==0)
		fprintf(f,"      call mtrini\n");
	fprintf(f,"      end\n\n");
	fprintf(f,"      subroutine mtrini\n      implicit none\n");
	if(FAver<8)
		fprintf(f,"#include \"model.h\"\n\n");
	else
		fprintf(f,"#include \"decl.h\"\n\n");
	fprintf(f,"      integer m1,m2,m3,m4\n\n");

	cls_write_matr(f);

	fprintf(f,"\n      end\n\n*************************************");
	if(FAver==4)
	fprintf(f,"**********\n\n      subroutine ModelVarIni(sqrtS, *)\n\
      implicit none\n\
      double precision sqrtS\n\
      double precision Alfas\n\
\n\
#include \"model.h\"\n\
\n\
c      double precision ALPHAS2\n\
c      external ALPHAS2\n\
\n\
c      Alfas = ALPHAS2(sqrtS)\n\
c      GG = sqrt(4*pi*Alfas)\n\
      end\n\n");
	 else if(FAver==6)
	fprintf(f,"**********\n\n      subroutine ModelVarIni(fail, sqrtS)\n\
      implicit none\n\
      double precision sqrtS\n\
      integer fail\n\
      double precision Alfas\n\
\n\
#include \"model.h\"\n\
\n\
c      double precision ALPHAS2\n\
c      external ALPHAS2\n\
\n\
c      Alfas = ALPHAS2(sqrtS)\n\
c      GG = sqrt(4*pi*Alfas)\n\
      fail=0\n\
      end\n\n");
	  else
	fprintf(f,"**********\n\n      subroutine ModelVarIni(fail, sqrtS)\n\
      implicit none\n\
      double precision sqrtS\n\
      integer fail\n\
      double precision Alfas\n\
\n\
#include \"decl.h\"\n\
\n\
c      double precision ALPHAS2\n\
c      external ALPHAS2\n\
\n\
c      Alfas = ALPHAS2(sqrtS)\n\
c      GG = sqrt(4*pi*Alfas)\n\
      fail=0\n\
      end\n\n");
		
	 

	fprintf(f,"************************************************\n\n");
	fprintf(f,"\
      subroutine ModelDigest\n\
      implicit none\n\
\n");
	if(FAver<8)
		fprintf(f,"#include \"model.h\"\n\n");
	else
		fprintf(f,"#include \"decl.h\"\n\n");
				
	for(li=params;li;li=ListTail(li))
	  {
	    Term t=CompoundName(ListFirst(li));
		if(strcmp(AtomValue(t),"A00001")==0)
			break;
	    fprintf(f,"      write(16,*) '%s=',%s\n",AtomValue(t),AtomValue(t));
	  } 

      fprintf(f,"\n      end\n\n");


	for(li=fainclude;li;li=ListTail(li))
		fprintf(f,"#include \"%s\"\n\n",AtomValue(ListFirst(li)));

/*	fprintf(f,"#include \"alphas2.h\"\n\n");*/

	fclose(f);

	NoQuotes=0;
}

void prm_decl_hc(FILE *f, List a2l)
{
	List l,l1,l2;
	if(opAbbArr==-1)
	{
	for(l=a2l;l;l=ListTail(l))
		{
		for(l2=CompoundArgN(ListFirst(l),3);l2;l2=ListTail(l2))
			if(GetAtomProperty(CompoundArg1(ListFirst(l2)),A_CHNAME))
				SetAtomProperty(CompoundArg1(ListFirst(l2)),A_CHNAME,0);
		for(l1=CompoundArgN(ListFirst(l),5);l1;l1=ListTail(l1))
		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
			if(GetAtomProperty(CompoundArg1(ListFirst(l2)),A_CHNAME))
				SetAtomProperty(CompoundArg1(ListFirst(l2)),A_CHNAME,0);
		}
	}
	if(opAbbArr)
	for(l=params;l;l=ListTail(l))
	{
		Atom p=CompoundName(ListFirst(l));
		Integer no=GetAtomProperty(p,A_CHNAME);
		if(no==0)
			continue;
		fprintf(f,"%s := %s[%ld]\n",AtomValue(p),
				GetAtomProperty(p,A_ANTI)?"AAABC":"AAABR",IntegerValue(no));
	}
	for(l=cparams;l;l=ListTail(l))
	{
		Atom z=ListFirst(l);
		Atom az=GetAtomProperty(z,A_ANTI);
		if(az)
			fprintf(f,"%s := Conjugate[%s]\n",AtomValue(az),AtomValue(z));
	}
}

void FAsqparam(FILE *f)
{
	List l;
	int z=0;

	for(l=params;l;l=ListTail(l))
	{
		Term t=ListFirst(l);
		Term val=CompoundArg1(t);
		if(is_compound(val)&&CompoundName(val)==OPR_POW &&
				is_atom(CompoundArg1(val))&&CompoundArg2(val)==NewInteger(2) &&
				!GetAtomProperty(CompoundArg1(val),A_HERMC))
		{
			fprintf(f,"%s^(n_?EvenQ) ^:= %s^(n/2);\n",
					AtomValue(CompoundArg1(val)),AtomValue(CompoundName(t)));
			z++;
		}
	}
	fprintf(f,"\n");
/*	if(z)
	{
		fprintf(f,"FormSubst = \"\\\n");
		for(l=params;l;l=ListTail(l))
		{
			Term t=ListFirst(l);
			Term val=CompoundArg1(t);
			if(is_compound(val)&&CompoundName(val)==OPR_POW &&
					is_atom(CompoundArg1(val))&&CompoundArg2(val)==NewInteger(2))
			{
				z--;
				fprintf(f,"id %s^2 = %s;\\n\\\nid %s^-2 = %s^-1;\\n%c\n",
				AtomValue(CompoundArg1(val)),AtomValue(CompoundName(t)),
				AtomValue(CompoundArg1(val)),AtomValue(CompoundName(t)),
						z?'\\':'\"');
			}
		}
	}
	fprintf(f,"\n");
*/
}

int SecondVaFu=0;
extern int opAutoWidths, pformat_def(void);

void WriteParameters(int fno, char *name)
	{
	FILE *f1, *f2;
	List li;
	int pnamel=6, fnamel=6;
	if(TexOutput)
		{
		texWriteParams(fno, name);
		return;
		}
	if(FAOutput)
		return;
	if(OutputDirectory!=NULL)
		sprintf(wpbuf,"%s/vars%d.mdl",OutputDirectory,fno);
	else
		sprintf(wpbuf,"vars%d.mdl",fno);
	f1=fopen(wpbuf,"w");
	if(f1==NULL)
		{
		printf("Can not open file \'%s\' for writing.\n",wpbuf);
		perror("");
		return;
		}
	if(pformat_def())
	{
		pnamel=7;
		fnamel=7;
	}
	if(OutputDirectory!=NULL)
		sprintf(wpbuf,"%s/func%d.mdl",OutputDirectory,fno);
	else
		sprintf(wpbuf,"func%d.mdl",fno);
	f2=fopen(wpbuf,"w");
	if(f2==NULL)
		{
		printf("Can not open file \'%s\' for writing.\n",wpbuf);
		perror("");
		return;
		}
	fprintf(f1,"%s\n Variables \n",name);
	fprintf(f2,"%s\n Constraints \n",name);
	if(MicroOmega)
		fprintf(f1," "),fprintf(f2," ");
	for(li=params;li;li=ListTail(li))
	{
		Term ttt, val, attt;
		int tttl;
		ttt=ListFirst(li);
		val=CompoundArg1(ttt);
		tttl=(int)strlen(AtomValue(CompoundName(ttt)));
		attt=GetAtomProperty(CompoundName(ttt),A_ANTI);
		if(attt)
			tttl+=(int)strlen(AtomValue(attt))+1;
		if(val==0 || is_integer(val) || is_float(val))
		{
			if(tttl>pnamel) pnamel=tttl;
		}
		else
		{
			if(tttl>fnamel) fnamel=tttl;
		}
	}

	fprintf(f1," Name");
	WriteBlank(f1,pnamel-5);
	fprintf(f1,"| Value       |>  Comment                                   <|\n");
	fprintf(f2," Name");
	WriteBlank(f2,fnamel-5);
	fprintf(f2,"|> Expression");
	WriteBlank(f2,P_D_WIDTH-13);
	if(pformat_def()==0)
		fprintf(f2,"<|> Comment                      <|\n");
	else
		fprintf(f2,"<|\n");
	li=params;
	while(!is_empty_list(li))
		{
		Term ttt, val, comm, attt;
		ttt=ListFirst(li);
		val=CompoundArg1(ttt);
		comm=CompoundArg2(ttt);
		attt=GetAtomProperty(CompoundName(ttt),A_ANTI);
		if(val==0) val=NewInteger(0);
		if(SecondVaFu && !(is_integer(val)||is_float(val))
				&& micro_allpr(CompoundName(ttt),val))
			val=NewFloat(EvalParameter(CompoundName(ttt)));
		if(MicroOmega)
			fprintf((is_integer(val)||is_float(val))?f1:f2,
			"%c",(micro_allpr(CompoundName(ttt),val)&&SecondVaFu==0)?'*':' ');
		if(is_integer(val) || is_float(val) || EvalPrm)
			{
			int sp=0;
			if(!UFOutput && opAutoWidths && GetAtomProperty(CompoundName(ttt),OPR_WIDTH)
				&& val==NewInteger(0))
				{li=ListTail(li); continue;}

			if(CalcOutput && comm && AtomValue(comm)[0]=='*')
				sp+=fprintf((is_integer(val)||is_float(val))?f1:f2,"*");
			
			sp+=fprintf(f1,"%s",AtomValue(FunctorName(CompoundFunctor(ttt))));
			if(attt)
				sp+=fprintf(f1,"/%s",AtomValue(attt));
			WriteBlank(f1,pnamel-sp);
			fprintf(f1,"|");
			if(is_integer(val) || is_float(val))
				sp=fWriteTerm(f1,val);
			else
				sp=fprintf(f1,"%g",EvalParameter(CompoundName(ttt)));
			WriteBlank(f1,13-sp);
			fprintf(f1,"|");
			if(comm!=0)
				fprintf(f1,"%s\n",AtomValue(comm));
			else
				fprintf(f1,"\n");
/*			if(comm)
				sp=fprintf(f1,"%s",AtomValue(comm));
			else
				sp=0;
			WriteBlank(f1,23-sp);fprintf(f1,"|\n");*/
			}
		else
			{
			int sp=0;
			if(CalcOutput && comm && AtomValue(comm)[0]=='*')
				sp+=fprintf((is_integer(val)||is_float(val))?f1:f2,"*");

			sp+=fprintf(f2,"%s",AtomValue(FunctorName(CompoundFunctor(ttt))));
			if(attt)
				sp+=fprintf(f2,"/%s",AtomValue(attt));
			WriteBlank(f2,fnamel-sp);
			fprintf(f2,"|");
			NoQuotes=1;
			val=CopyTerm(val);
			if(ChepVersion>3)
				rpl_pow(val);
			if(AtomValue(CompoundName(ttt))[0]=='%'&&val==A_I)
				sp=0;
			else
				sp=fWriteTerm(f2,val);
			NoQuotes=0;
			if(sp>longest_pdline)
			{
				longest_pdline=sp;
				llparam=CompoundName(ttt);
			}
			FreeAtomic(val);
			if(pformat_def()==0)
			{
				WriteBlank(f2,P_D_WIDTH-sp);
				fprintf(f2,"|");
				if(comm!=0)
					fprintf(f2,"%s\n",AtomValue(comm));
				else
					fprintf(f2,"\n");
			}
			else
			{
				if(comm)
					fprintf(f2," %% %s",AtomValue(comm));
				fprintf(f2,"\n");
			}
/*			if(comm)
				sp=fprintf(f2,"%s",AtomValue(comm));
			else
				sp=0;
			WriteBlank(f2,47-sp);fprintf(f2,"|\n");*/
			}
		li=ListTail(li);
		}
	fclose(f1);
	fclose(f2);
	}

void ClearParameter(Atom p)
	{
	List l;
	l=params;
	while(!is_empty_list(l))
		{
		if(CompoundName(ListFirst(l))==p)
			{
			params=CutFromList(params,l);
			break;
			}
		l=ListTail(l);
		}
	l=rparams;
	while(!is_empty_list(l))
		{
		if(ListFirst(l)==p)
			{
			rparams=CutFromList(rparams,l);
			break;
			}
		l=ListTail(l);
		}
	l=cparams;
	while(!is_empty_list(l))
		{
		if(ListFirst(l)==p)
			{
			cparams=CutFromList(cparams,l);
			break;
			}
		l=ListTail(l);
		}
	RemoveAtomProperty(p,PROP_TYPE);
	}



void ChangeParameterValue(Atom t, double val)
	{
	List l;
	for(l=params;l;l=ListTail(l))
		{
		if(CompoundName(ListFirst(l))==t)
			{
			SetCompoundArg(ListFirst(l),1,NewFloat(val));
			return;
			}
		}
	puts("Internal error (prms01)");
	}


double sort4(double m1, double m2, double m3, double m4, double dn)
{
	int n;
	int i,f=0;
	double m[4];
	n=(int)floor(dn-0.5);
	m[0]=m1;
	m[1]=m2;
	m[2]=m3;
	m[3]=m4;

	do
	{
		f=0;
		for(i=0;i<3;i++)
			if(fabs(m[i])>fabs(m[i+1]))
			{
				double tmp;
				tmp=m[i];
				m[i]=m[i+1];
				m[i+1]=tmp;
				f=1;
			}

	} while(f);

	return m[n];
}

/*
extern int    slhaRead(char *fname,int mode);
extern double slhaVal(char * Block, double Q, int nKey, ...);
extern double rDiagonal(double nDim,...);
extern double rDiagonal2(double nDim,...);
extern double  MassArray(int id,  int i);
extern double  MixMatrix(int id,  int i, int j);
extern double  MixMatrixU(int id,  int i, int j);
*/

#include "SLHAplus/SLHAplus.h"

static double eval_ef(Atom, int, double *);

static cmplx cx_pow(cmplx c, int p)
{
    cmplx ret;
	int i;
	ret.r=1.0;
	ret.i=0.0;
	if(p<0)
	{
	  cmplx t=c;
	  c.r=t.r/(t.r*t.r+t.i*t.i);
	  c.i=-t.i/(t.r*t.r+t.i*t.i);
	  p=-p;
	}
	for(i=0;i<p;i++)
	{
	  cmplx t;
	  t.r=ret.r*c.r-ret.i*c.i;
	  t.i=ret.r*c.i+ret.i*c.r;
	  ret=t;
	}
	return ret;
}

cmplx cEvalParameter(Term t)
{
  cmplx ret;
  ret.r=ret.i=0.0;
  
  if(t==0)
		return ret;

	if(t==A_SQRT2)
	{
		ret.r=sqrt((double)2.0);
		return ret;
	}
	
	if(t==A_I)
	{
	  ret.i=1.0;
		return ret;
	}
  
  if(is_integer(t))
	{
		ret.r=(double)IntegerValue(t);
		return ret;
	}
	
	if(is_float(t))
		return ComplexValue(t);
	
	if(is_atom(t) && is_parameter(t))
		{
		List l;
		Atom at;

		l=GetAtomProperty(t,A_DUMMY_PRM);
		if(l)
		{
			
			List l1;
			for(l1=l;l1;l1=ListTail(l1))
			{
			  cmplx r1=cEvalParameter(CompoundArg2(ListFirst(l1)));
			  r1.r*=IntegerValue(CompoundArg1(ListFirst(l1)));
			  r1.i*=IntegerValue(CompoundArg1(ListFirst(l1)));
			  ret.r+=r1.r;
			  ret.i+=r1.i;
			}	
			return ret;
		}

		at=GetAtomProperty(t,A_ANTI);
		if(at==0)
		{
		  ret.r=EvalParameter(t);
		  return ret;
		}

		for(l=params;l;l=ListTail(l))
			{
			Atom curt=CompoundName(ListFirst(l));
			if(curt==t || curt==at)
				{
				Term t1;
				double dd;
				t1=CompoundArg1(ListFirst(l));
				if(is_integer(t1))
				{
					ret.r= (double)IntegerValue(t1);
					return ret;
				}
				
				if(is_float(t1))
				{
					ret= ComplexValue(t1);
					if(curt==at) ret.i=-ret.i;
					return ret;
				}
				if(t1==0)
					return ret;
				if(CompoundArgN(ListFirst(l),3)!=0)
				{
					ret=ComplexValue(CompoundArgN(ListFirst(l),3));
					if(curt==at) ret.i=-ret.i;
					return ret;
				}
				ret=cEvalParameter(t1);
				SetCompoundArg(ListFirst(l),3,NewComplex(ret));
				if(curt==at) ret.i=-ret.i;
				return ret;
				}
			}
		if(GetAtomProperty(t,A_EXT_FUNC))
		{
			if(strncmp(AtomValue(t),"initDiagonal",12)==0)
			{initDiagonal();
			  return ret;
			}
			else
				return ret;
		}

		printf("Internal error (prm01'%s')\n",AtomValue(t));
		return ret;
		}
	
	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==OPR_PLUS)
	{
	  cmplx c1=cEvalParameter(CompoundArg1(t)), c2=cEvalParameter(CompoundArg2(t));
	  ret.r=c1.r+c2.r;
	  ret.i=c1.i+c2.i;
	  return ret;
	}
	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==OPR_MINUS)
	{
	  cmplx c1=cEvalParameter(CompoundArg1(t)), c2=cEvalParameter(CompoundArg2(t));
	  ret.r=c1.r-c2.r;
	  ret.i=c1.i-c2.i;
	  return ret;
	}
	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==OPR_MLT)
	{
	  cmplx c1=cEvalParameter(CompoundArg1(t)), c2=cEvalParameter(CompoundArg2(t));
	  ret.r=c1.r*c2.r-c1.i*c2.i;
	  ret.i=c1.r*c2.i+c1.i*c2.r;
	  return ret;
	}
	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==OPR_DIV)
	{
	  cmplx c1=cEvalParameter(CompoundArg1(t)), c2=cEvalParameter(CompoundArg2(t));
	  ret.r=(c1.r*c2.r+c1.i*c2.i)/(c2.r*c2.r+c2.i*c2.i);
	  ret.i=(-c1.r*c2.i+c1.i*c2.r);
	  return ret;
	}
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==OPR_MINUS)
	{
	  ret=cEvalParameter(CompoundArg1(t));
	  ret.r=-ret.r;
	  ret.i=-ret.i;
	  return ret;
	}
	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==OPR_POW)
	{
	  cmplx c1=cEvalParameter(CompoundArg1(t));
	  Term pw=CompoundArg2(t);
	  if(!is_integer(pw))
	  {
		ErrorInfo(0);
		puts("Complex power is not implemented.");
		return ret;
	  }
	  ret=cx_pow(c1,(int)IntegerValue(pw));
	  return ret;
	}
	
	if(is_list(t))
		{
		List l;
		ret.r=1.0;
		for(l=t;l;l=ListTail(l))
			{
			cmplx tpv, tmp;
			int tpw;
			tpv=cEvalParameter(CompoundArg1(ListFirst(l)));
			tpw=(int)IntegerValue(CompoundArg2(ListFirst(l)));
			tpv=cx_pow(tpv,tpw);
			/*WriteTerm(ListFirst(l)); printf("=(%f,%f)\n",tpv.r,tpv.i);*/
			tmp.r=tpv.r*ret.r-tpv.i*ret.i;
			tmp.i=tpv.i*ret.r+tpv.r*ret.i;
			ret=tmp;
			}
		return ret;
		}
	
}

double EvalParameter(Term t)
	{

	if(t==0)
		return 0.0;

	if(t==A_SQRT2)
		return sqrt((double)2.0);
	if(t==A_I)
		return 1.0;

	if(is_integer(t))
	  return (double)IntegerValue(t);

	if(is_float(t))
		return FloatValue(t);
	
	if(is_atom(t) && is_parameter(t))
		{
		List l;


		l=GetAtomProperty(t,A_DUMMY_PRM);
		if(l)
		{
			double ret=0;
			List l1;
			for(l1=l;l1;l1=ListTail(l1))
				ret+=IntegerValue(CompoundArg1(ListFirst(l1)))*
						EvalParameter(CompoundArg2(ListFirst(l1)));
			return ret;
		}

		for(l=params;l;l=ListTail(l))
			{
			if(CompoundName(ListFirst(l))==t)
				{
				Term t1;
				double dd;
				t1=CompoundArg1(ListFirst(l));
				if(is_integer(t1))
					return (double)IntegerValue(t1);
				if(is_float(t1))
					return FloatValue(t1);
				if(t1==0)
					return 0.0;
				if(CompoundArgN(ListFirst(l),3)!=0)
					return FloatValue(CompoundArgN(ListFirst(l),3));
				dd=EvalParameter(t1);
				SetCompoundArg(ListFirst(l),3,NewFloat(dd));
				return dd;
				}
			}
		if(GetAtomProperty(t,A_EXT_FUNC))
		{
			if(strncmp(AtomValue(t),"initDiagonal",12)==0)
				return initDiagonal();
			else
				return 0.0;
		}

		if(GetAtomProperty(t,A_ANTI))
		{
			static int prm01w=0;
			if(prm01w==0)
			{
				prm01w=1;
				WarningInfo(0);
				puts("Numeric evaluation of complex mumbers is not implemented.");
			}
			return 0.0;
		}
		printf("Internal error (prm01'%s')\n",AtomValue(t));
		return 0.0;
		}

	if(is_compound(t) && GetAtomProperty(CompoundName(t),A_EXT_FUNC) &&
		GetAtomProperty(CompoundName(t),CompoundName(t)))
	{
		double arr[101];
		int i;
		for(i=1;i<=CompoundArity(t);i++)
			arr[i-1]=EvalParameter(CompoundArgN(t,i));
		return eval_ef(CompoundName(t),CompoundArity(t),(double *)arr);
	}

	if(is_compound(t) && strcmp(AtomValue(CompoundName(t)),"sort4")==0
			&& CompoundArity(t)==5)
	{
		double t1,t2,t3,t4,t5;
		t1=EvalParameter(CompoundArgN(t,1));
		t2=EvalParameter(CompoundArgN(t,2));
		t3=EvalParameter(CompoundArgN(t,3));
		t4=EvalParameter(CompoundArgN(t,4));
		t5=EvalParameter(CompoundArgN(t,5));
		t1=sort4(t1,t2,t3,t4,t5);
		return t1;

	}

	if(is_compound(t) && strcmp(AtomValue(CompoundName(t)),"MassArray")==0
			&& CompoundArity(t)==2)
	{
		double t1,t2;
		t1=EvalParameter(CompoundArgN(t,1));
		t2=EvalParameter(CompoundArgN(t,2));
		t2=MassArray((int)t1,(int)t2);
		return t2;
	}
	if(is_compound(t) && strcmp(AtomValue(CompoundName(t)),"MixMatrix")==0
			&& CompoundArity(t)==3)
	{
		double t1,t2,t3;
		t1=EvalParameter(CompoundArgN(t,1));
		t2=EvalParameter(CompoundArgN(t,2));
		t3=EvalParameter(CompoundArgN(t,3));
		t2=MixMatrix((int)t1,(int)t2,(int)t3);
		return t2;

	}
	if(is_compound(t) && strcmp(AtomValue(CompoundName(t)),"MixMatrixU")==0
			&& CompoundArity(t)==3)
	{
		double t1,t2,t3;
		t1=EvalParameter(CompoundArgN(t,1));
		t2=EvalParameter(CompoundArgN(t,2));
		t3=EvalParameter(CompoundArgN(t,3));
		t2=MixMatrixU((int)t1,(int)t2,(int)t3);
		return t2;

	}

	if(is_compound(t) && strcmp(AtomValue(CompoundName(t)),"slhaVal")==0
			&& CompoundArity(t)==4)
	{
		double t2,t3,t4;
		Atom va=GetAtomProperty(CompoundArg1(t),A_I);
		if(va==0)
		{
			printf("Error: first arg of slhaVal must be a string\n");
			return 0.0;
		}
		t2=EvalParameter(CompoundArgN(t,2));
		t3=EvalParameter(CompoundArgN(t,3));
		t4=EvalParameter(CompoundArgN(t,4));

		t2=slhaVal(AtomValue(va),t2,(int)t3,(int)t4);
		return t2;

	}
    
    if(is_compound(t) && strcmp(AtomValue(CompoundName(t)),"slhaVal1")==0
           && CompoundArity(t)==3)
        {
            double t2,t3,t4;
            Atom va=GetAtomProperty(CompoundArg1(t),A_I);
            if(va==0)
            {
                printf("Error: first arg of slhaVal must be a string\n");
                return 0.0;
            }
            t2=EvalParameter(CompoundArgN(t,2));
            t3=1.0;
            t4=EvalParameter(CompoundArgN(t,3));
            
            t2=slhaVal(AtomValue(va),t2,1,(int)t4);
            return t2;
            
        }
        if(is_compound(t) && strcmp(AtomValue(CompoundName(t)),"slhaVal2")==0
           && CompoundArity(t)==4)
        {
            double t2,t3,t4;
            Atom va=GetAtomProperty(CompoundArg1(t),A_I);
            if(va==0)
            {
                printf("Error: first arg of slhaVal must be a string\n");
                return 0.0;
            }
            t2=EvalParameter(CompoundArgN(t,2));
            t3=EvalParameter(CompoundArgN(t,3));
            t4=EvalParameter(CompoundArgN(t,4));
            
            t2=slhaVal(AtomValue(va),t2,2,(int)t3,(int)t4);
            return t2;
            
        }
    
        
	if(is_compound(t) && strcmp(AtomValue(CompoundName(t)),"slhaVal")==0
			&& CompoundArity(t)==5)
	{
		double t2,t3,t4,t5;
		Atom va=GetAtomProperty(CompoundArg1(t),A_I);
		if(va==0)
		{
			printf("Error: first arg of slhaVal must be a string\n");
			return 0.0;
		}
		t2=EvalParameter(CompoundArgN(t,2));
		t3=EvalParameter(CompoundArgN(t,3));
		t4=EvalParameter(CompoundArgN(t,4));
		t5=EvalParameter(CompoundArgN(t,5));
		t2=slhaVal(AtomValue(va),t2,(int)t3,(int)t4,(int)t5);
		return t2;

	}

	if(is_compound(t) && GetAtomProperty(CompoundName(t),A_EXT_FUNC))
	{
		return 0.0;
	}

	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==OPR_PLUS)
		return EvalParameter(CompoundArg1(t))+EvalParameter(CompoundArg2(t));
	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==OPR_MINUS)
		return EvalParameter(CompoundArg1(t))-EvalParameter(CompoundArg2(t));
	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==OPR_MLT)
		return EvalParameter(CompoundArg1(t))*EvalParameter(CompoundArg2(t));
	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==OPR_DIV)
		return EvalParameter(CompoundArg1(t))/EvalParameter(CompoundArg2(t));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==OPR_MINUS)
		return -EvalParameter(CompoundArg1(t));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==A_SQRT)
		return sqrt(EvalParameter(CompoundArg1(t)));
	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==A_CURT)
		return cu_rt(EvalParameter(CompoundArg1(t)),EvalParameter(CompoundArg2(t)));
	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==OPR_POW)
		return pow(EvalParameter(CompoundArg1(t)),EvalParameter(CompoundArg2(t)));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==A_SIN)
		return sin(EvalParameter(CompoundArg1(t)));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==A_COS)
		return cos(EvalParameter(CompoundArg1(t)));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==A_FABS)
		return fabs(EvalParameter(CompoundArg1(t)));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==A_TAN)
		return tan(EvalParameter(CompoundArg1(t)));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==A_ASIN)
		return asin(EvalParameter(CompoundArg1(t)));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==A_ACOS)
		return acos(EvalParameter(CompoundArg1(t)));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==A_ATAN)
		return atan(EvalParameter(CompoundArg1(t)));
	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==A_ATAN2)
		return atan2(EvalParameter(CompoundArg1(t)),EvalParameter(CompoundArg2(t)));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==A_LOG)
		return log(EvalParameter(CompoundArg1(t)));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==A_EXP)
		return exp(EvalParameter(CompoundArg1(t)));
	if(is_compound(t) && CompoundArity(t)==3 && CompoundName(t)==A_IF)
	{
		double t1;
		t1=EvalParameter(CompoundArg1(t));
		if(t1>0)
			t1=EvalParameter(CompoundArg2(t));
		else
			t1=EvalParameter(CompoundArgN(t,3));
		return t1;
	}

	if(is_list(t))
		{
		double ret=1.0;
		List l;
		for(l=t;l;l=ListTail(l))
			{
			double tpv;
			int tpw,i;
			tpv=EvalParameter(CompoundArg1(ListFirst(l)));
			tpw=(int)IntegerValue(CompoundArg2(ListFirst(l)));
			
			if(tpw>0)
				for(i=0;i<tpw;i++)
					ret*=tpv;
			else
				for(i=0;i>tpw;i--)
					ret/=tpv;
			}
		return ret;
		}

	printf("Error: can not evaluate ");
	WriteTerm(t);
	puts("");

	return 0.0;
	}

Term ProcslhaRead(Term t, Term ind)
	{
	if(!is_compound(t) || CompoundArity(t)!=1)
	{
		ErrorInfo(0);puts("wrong syntax in slhaRead statement.");
		return 0;
	}
	if(!is_atom(CompoundArg1(t)))
	{
		ErrorInfo(0);puts("slhaRead: file name is expected.");
		return 0;
	}
	slhaRead(AtomValue(CompoundArg1(t)),3);
	return 0;
	}

Term InterfEvalParam(Term t, Term ind)
	{
	cmplx ret;
	if(!is_compound(t) || CompoundArity(t)!=1)
		return 0;
	ret=cEvalParameter(CompoundArg1(t));
	WriteTerm(CompoundArg1(t));
	if(ret.i==0)
	  printf("=%f\n",ret.r);
	else
	  printf("=(%f,%f)\n",ret.r,ret.i);
	FreeAtomic(t);
	return 0;

	}

Term ProcTailPrm(Term t, Term ind)
{
	List l1,l2;

	if(!is_compound(t) || CompoundArity(t)!=1 || !is_list(CompoundArg1(t)))
	{
		ErrorInfo(112);
		puts("wrong syntax in 'tail_prm' call.");
		return 0;
	}

	l1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	t=l1;
	l2=NewList();
	for(l1=t;l1;l1=ListTail(l1))
	{
		if(!is_atom(ListFirst(l1)))
		{
			ErrorInfo(112);
			printf("tail_prm: '");
			WriteTerm(ListFirst(l1));
			puts("': only parameters expected.");
			continue;
		}
		if(!is_parameter(ListFirst(l1)))
		{
			ErrorInfo(112);
			printf("tail_prm: '");
			WriteTerm(ListFirst(l1));
			puts("': it is not a parameter.");
			continue;
		}
	}

	for(l1=params;l1;l1=ListTail(l1))
	{
		if(ListMember(t,CompoundName(ListFirst(l1))))
		{
			Term t1;
			t1=ListFirst(l1);
			l2=AppendLast(l2,t1);
			ChangeList(l1,0);
		}
	}

rpt:
	for(l1=params;l1;l1=ListTail(l1))
	{
		if(ListFirst(l1)==0)
		{
			params=CutFromList(params,l1);
			goto rpt;
		}
	}

	params=ConcatList(params,l2);
/*
	for(l1=params;l1;l1=ListTail(l1))
	{
		if(ListMember(t,CompoundName(ListFirst(l1))))
		{
			Term t1;
			double d;
			t1=ListFirst(l1);
			if(!is_integer(CompoundArg1(t1)) && !is_float(CompoundArg1(t1)))
			{
				d=EvalParameter(CompoundName(t1));
				SetCompoundArg(t1,1,NewFloat(d));
			}
		}
	}
*/
	return 0;
}

Term ProcCHEPPrm(Term t, Term ind)
{
	int mono;
	char fnbuf[512];
	int has_fn=0;
	FILE *fin, *fout;

	if(!is_compound(t) || CompoundArity(t)<1)
	{
		ErrorInfo(113);
		puts("wrong syntax in 'read_chep_prm' statement");
		return 0;
	}

	if(!is_integer(CompoundArg1(t)))
	{
		ErrorInfo(113);
		puts("'read_chep_prm': the first argument should be integer.");
		return 0;
	}

	mono=(int)IntegerValue(CompoundArg1(t));

	if(CompoundArity(t)>1)
		has_fn=1;

	if(has_fn && !is_atom(CompoundArg2(t)))
	{
		ErrorInfo(113);
		puts("'read_chep_prm': the second argument (directory) should be string const");
		return 0;
	}

	if(has_fn)
		sprintf(fnbuf,"%s/prm_from_chep.mdl",AtomValue(CompoundArg2(t)));
	else
		sprintf(fnbuf,"prm_from_chep.mdl");

	fout=fopen(fnbuf,"w");
	if(fout==NULL)
	{
		ErrorInfo(113);
		printf("'read_chep_prm': can not open output file '%s'\n",fnbuf);
		perror("'read_chep_prm':");
		return 0;
	}

	if(has_fn)
		sprintf(fnbuf,"%s/vars%d.mdl",AtomValue(CompoundArg2(t)),mono);
	else
		sprintf(fnbuf,"vars%d.mdl",mono);

	fin=fopen(fnbuf,"r");
	if(fin==NULL)
	{
		ErrorInfo(113);
		printf("'read_chep_prm': can not open input file '%s'\n",fnbuf);
		perror("'read_chep_prm':");
		return 0;
	}

	fprintf(fout,"%%\n%% %s\n%%\n\n",fnbuf);

	fgets(fnbuf,510,fin);
	fgets(fnbuf,510,fin);
	fgets(fnbuf,510,fin);

	while(fgets(fnbuf,510,fin)!=NULL)
	{
		int i, pos1=0, pos2=0;
		for(i=0;i<500;i++)
			if(fnbuf[i]=='|')
				break;

		if(i==500)
		{
			ErrorInfo(113);
			puts("'read_chep_prm': corrupted 'vars' file.");
			return 0;
		}

		pos1=i;

		for(i=pos1+1;i<500;i++)
			if(fnbuf[i]=='|')
				break;

		if(i==500)
		{
			ErrorInfo(113);
			puts("'read_chep_prm': corrupted 'vars' file.");
			return 0;
		}

		pos2=i;

		fnbuf[pos1]=0;
		fnbuf[pos2]=0;

		while(fnbuf[pos2-1]==' ')
		{
			pos2--;
			fnbuf[pos2]=0;
		}

		fprintf(fout,"parameter %s = %s.\n",fnbuf,fnbuf+pos1+1);
	}

	fclose(fin);

	if(has_fn)
		sprintf(fnbuf,"%s/func%d.mdl",AtomValue(CompoundArg2(t)),mono);
	else
		sprintf(fnbuf,"func%d.mdl",mono);

	fin=fopen(fnbuf,"r");
	if(fin==NULL)
	{
		ErrorInfo(113);
		printf("'read_chep_prm': can not open input file '%s'\n",fnbuf);
		perror("'read_chep_prm':");
		return 0;
	}

	fprintf(fout,"\n%%\n%% %s\n%%\n\n",fnbuf);

	fgets(fnbuf,510,fin);
	fgets(fnbuf,510,fin);
	fgets(fnbuf,510,fin);

	while(fgets(fnbuf,510,fin)!=NULL)
	{
		int i, pos1=0, pos2=0;
		for(i=0;i<500;i++)
			if(fnbuf[i]=='|')
				break;

		if(i==500)
		{
			ErrorInfo(113);
			puts("'read_chep_prm': corrupted 'vars' file.");
			return 0;
		}

		pos1=i;

		for(i=pos1+1;i<500;i++)
			if(fnbuf[i]=='|')
				break;

		if(i==500)
		{
			ErrorInfo(113);
			puts("'read_chep_prm': corrupted 'vars' file.");
			return 0;
		}

		pos2=i;

		fnbuf[pos1]=0;
		fnbuf[pos2]=0;

		while(fnbuf[pos2-1]==' ')
		{
			pos2--;
			fnbuf[pos2]=0;
		}

		fprintf(fout,"parameter %s = %s.\n",fnbuf,fnbuf+pos1+1);
	}

	fclose(fin);
	fprintf(fout,"\n");
	fclose(fout);

	if(has_fn)
		sprintf(fnbuf,"%s/prm_from_chep.mdl",AtomValue(CompoundArg2(t)));
	else
		sprintf(fnbuf,"prm_from_chep.mdl");

	ReadFile(fnbuf);

	FreeAtomic(t);

	return 0;
}

static struct efel
	{
	int no;
	Atom fname, ffile;
	void *handler;
	extfunc func;
	struct efel *next;
	} *flist=NULL;

static int extfuncno=0;
int extlib_problems=0;

Term ProcExtFunc(Term t, Term ind)
{
	Atom fl=0;


	if(!is_compound(t)|| CompoundArity(t)<2 || CompoundArity(t)>3 || 
		!is_atom(CompoundArg1(t)))
	{
		ErrorInfo(228);
		printf("Illegal syntax in external_func statement.\n");
		return 0;
	}

	if(is_integer(CompoundArg2(t))&&IntegerValue(CompoundArg2(t))==0)
	{
		char cbuf[128];
		Atom name=CompoundArg1(t),name2;
		Term t1;
		sprintf(cbuf,"%s()",AtomValue(name));
		name2=NewAtom(cbuf,0);
		SetAtomProperty(name2,PROP_TYPE,OPR_PARAMETER);
		t1=MakeCompound1(OPR_ALIAS,MakeCompound2(OPR_EQSIGN,name,name2));
		/*printf("extfu0: %s %s\n",AtomValue(name),AtomValue(name2));*/
		CallFunction(t1,0);
		if(CompoundArity(t)==3)
		{
			WarningInfo(0);
			printf(
			"external_func: evaluation of 0 arguments funcion is not supported\n");
		}
		SetAtomProperty(name2, A_EXT_FUNC, CompoundArg2(t));
		if(!ListMember(ExtFuncList,name))
		  ExtFuncList=AppendLast(ExtFuncList,name);
		return 0;
	}

	if(CompoundArg2(t)==OPR_MLT)
		SetCompoundArg(t,2,NewInteger(0));

	if(!is_integer(CompoundArg2(t)) || IntegerValue(CompoundArg2(t))<0 || 
			IntegerValue(CompoundArg2(t))>100)
	{
		ErrorInfo(229);
		puts("external_func: numer of arguments must be 0 to 100, or '*'.");
		return 0;
	}

	if(CompoundArity(t)==2 && 
			strcmp(AtomValue(CompoundArg1(t)),"initQCD")==0 &&
			CompoundArg2(t)==NewInteger(4))
	{
		struct efel *ep;
		void *ha;
		extfunc hf;
		hf=(extfunc)&initQCD;
		fl=A_I;
		extfuncno++;
		ep=(struct efel *)malloc(sizeof(struct efel));
		ep->no=extfuncno;
		ep->fname=CompoundArg1(t);
		ep->ffile=fl;
		ep->handler=ha;
		ep->func=hf;
		ep->next=flist;
		flist=ep;
	}


	if(CompoundArity(t)==2 && CompoundArg2(t)==NewInteger(1) &&
			(strcmp(AtomValue(CompoundArg1(t)),"McRun")==0 ||
			strcmp(AtomValue(CompoundArg1(t)),"MbRun")==0 ||
			strcmp(AtomValue(CompoundArg1(t)),"MtRun")==0 ||
			strcmp(AtomValue(CompoundArg1(t)),"alphaQCD")==0 ||
			strcmp(AtomValue(CompoundArg1(t)),"McEff")==0 ||
			strcmp(AtomValue(CompoundArg1(t)),"MbEff")==0 ||
			strcmp(AtomValue(CompoundArg1(t)),"MtEff")==0 ) )
	{
		struct efel *ep;
		void *ha;
		extfunc hf;
		if(strcmp(AtomValue(CompoundArg1(t)),"McRun")==0)
			hf=(extfunc)&McRun;
		else if(strcmp(AtomValue(CompoundArg1(t)),"MbRun")==0)
			hf=(extfunc)&MbRun;
		else if(strcmp(AtomValue(CompoundArg1(t)),"MtRun")==0)
			hf=(extfunc)&MtRun;
		else if(strcmp(AtomValue(CompoundArg1(t)),"McEff")==0)
			hf=(extfunc)&McEff;
		else if(strcmp(AtomValue(CompoundArg1(t)),"MbEff")==0)
			hf=(extfunc)&MbEff;
		else if(strcmp(AtomValue(CompoundArg1(t)),"MtEff")==0)
			hf=(extfunc)&MtEff;
		else
			hf=(extfunc)&alphaQCD;
		fl=A_I;
		extfuncno++;
		ep=(struct efel *)malloc(sizeof(struct efel));
		ep->no=extfuncno;
		ep->fname=CompoundArg1(t);
		ep->ffile=fl;
		ep->handler=ha;
		ep->func=hf;
		ep->next=flist;
		flist=ep;
	}

	if(CompoundArity(t)==2 && 
		(strcmp(AtomValue(CompoundArg1(t)),"rDiagonal")==0||
		strcmp(AtomValue(CompoundArg1(t)),"rDiagonalA")==0))
	{
		struct efel *ep;
		void *ha;
		extfunc hf;
		if(strcmp(AtomValue(CompoundArg1(t)),"rDiagonal")==0)
			hf=&rDiagonal;
		else
			hf=&rDiagonalA;
		fl=A_I;
		extfuncno++;
		ep=(struct efel *)malloc(sizeof(struct efel));
		ep->no=extfuncno;
		ep->fname=CompoundArg1(t);
		ep->ffile=fl;
		ep->handler=ha;
		ep->func=hf;
		ep->next=flist;
		flist=ep;
	}

	if(CompoundArity(t)==3 && is_atom(CompoundArgN(t,3)))
	{
		struct efel *ep;
		void *ha;
		extfunc hf;
		fl=CompoundArgN(t,3);
		for(ep=flist;ep;ep=ep->next)
			if(ep->ffile==fl)
				break;
		if(ep==0)
		{
			ha=dlopen(AtomValue(fl),RTLD_LAZY);
			if(ha==0)
				{
				char cbuf[128];
				sprintf(cbuf,"./%s",AtomValue(fl));
				ha=dlopen(cbuf,RTLD_LAZY);
				/*if(ha==0){dlerror();ha=dlopen(AtomValue(fl),RTLD_LAZY);}*/
				}
			if(ha==0)
				{
				WarningInfo(9);
				printf("external_func: %s\n",dlerror());
				extlib_problems++;
				fl=0;
				goto cnt;
				}
		}
		else
			ha=ep->handler;
		dlerror();
		hf=(extfunc)dlsym(ha,AtomValue(CompoundArg1(t)));
		if(hf==NULL)
		{
			WarningInfo(10);
			printf("external_func: %s\n",dlerror());
			extlib_problems++;
			if(ep==0)
				dlclose(ha);
			fl=0;
			goto cnt;
		}
		extfuncno++;
		ep=(struct efel *)malloc(sizeof(struct efel));
		ep->no=extfuncno;
		ep->fname=CompoundArg1(t);
		ep->ffile=fl;
		ep->handler=ha;
		ep->func=hf;
		ep->next=flist;
		flist=ep;
	}
	
	if(CompoundArity(t)==3 && is_compound(CompoundArgN(t,3)))
	{
		int i;
		Term t3=CompoundArgN(t,3);
		List ll=0;
		for(i=1;i<=CompoundArity(t3);i++)
			ll=AppendLast(ll,CompoundArgN(t3,i));
		SetAtomProperty(CompoundArg1(t),A_INTEGER,ll);
	}
		

cnt:
	SetAtomProperty(CompoundArg1(t), A_EXT_FUNC, CompoundArg2(t));
	if(fl)
		SetAtomProperty(CompoundArg1(t),CompoundArg1(t),NewInteger(extfuncno));
	if(!ListMember(ExtFuncList,CompoundArg1(t)))
	  ExtFuncList=AppendLast(ExtFuncList,CompoundArg1(t));
	return 0;
}

static double eval_ef(Atom fname, int ar, double *a)
	{
	struct efel *e;
	extfunc f;
	for(e=flist;e;e=e->next)
		if(e->fname==fname)
			break;
/*	printf("calling %s...\n",AtomValue(fname));*/
	if(e==0)
		{
		puts("Internal error (evlef01)");
		return 0.0;
		}
	f=e->func;
	switch(ar)
		{
		case 1: return f(a[0]);
		case 2: return f(a[0],a[1]);
		case 3: return f(a[0],a[1],a[2]);
		case 4: return f(a[0],a[1],a[2],a[3]);
		case 5: return f(a[0],a[1],a[2],a[3],a[4]);
		case 6: return f(a[0],a[1],a[2],a[3],a[4],a[5]);
		case 7: return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6]);
		case 8: return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7]);
		case 9: return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8]);
		case 10:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9]);

		case 11:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10]);
		case 12:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11]);
		case 13:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12]);
		case 14:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13]);
		case 15:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14]);
		case 16:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15]);
		case 17:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16]);
		case 18:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17]);
		case 19:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18]);
		case 20:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19]);
		case 21:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
				a[20]);
		case 22:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
				a[20],a[21]);
		case 23:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
				a[20],a[21],a[22]);
		case 24:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
				a[20],a[21],a[22],a[23]);
		case 25:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
				a[20],a[21],a[22],a[23],a[24]);
		case 26:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
				a[20],a[21],a[22],a[23],a[24],a[25]);
		case 27:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
				a[20],a[21],a[22],a[23],a[24],a[25],a[26]);
		case 28:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
				a[20],a[21],a[22],a[23],a[24],a[25],a[26],a[27]);
		case 29:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
				a[20],a[21],a[22],a[23],a[24],a[25],a[26],a[27],a[28]);
case  30:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30]);
case  31:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31]);
case  32:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32]);
case  33:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33]);
case  34:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34]);
case  35:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35]);
case  36:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36]);
case  37:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37]);
case  38:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38]);
case  39:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39]);
case  40:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40]);
case  41:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41]);
case  42:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42]);
case  43:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43]);
case  44:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44]);
case  45:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45]);
case  46:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46]);
case  47:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47]);
case  48:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48]);
case  49:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49]);
case  50:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50]);
case  51:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51]);
case  52:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52]);
case  53:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53]);
case  54:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54]);
case  55:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55]);
case  56:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56]);
case  57:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57]);
case  58:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58]);
case  59:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59]);
case  60:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60]);
case  61:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61]);
case  62:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62]);
case  63:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63]);
case  64:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64]);
case  65:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65]);
case  66:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66]);
case  67:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67]);
case  68:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68]);
case  69:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69]);
case  70:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70]);
case  71:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71]);
case  72:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72]);
case  73:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73]);
case  74:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74]);
case  75:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75]);
case  76:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76]);
case  77:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77]);
case  78:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78]);
case  79:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79]);
case  80:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80]);
case  81:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80],a[ 81]);
case  82:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80],a[ 81],a[ 82]);
case  83:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80],a[ 81],a[ 82],a[ 83]);
case  84:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80],a[ 81],a[ 82],a[ 83],a[ 84]);
case  85:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80],a[ 81],a[ 82],a[ 83],a[ 84],a[ 85]);
case  86:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80],a[ 81],a[ 82],a[ 83],a[ 84],a[ 85],a[ 86]);
case  87:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80],a[ 81],a[ 82],a[ 83],a[ 84],a[ 85],a[ 86],a[ 87]);
case  88:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80],a[ 81],a[ 82],a[ 83],a[ 84],a[ 85],a[ 86],a[ 87],a[ 88]);
case  89:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80],a[ 81],a[ 82],a[ 83],a[ 84],a[ 85],a[ 86],a[ 87],a[ 88],a[ 89]);
case  90:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80],a[ 81],a[ 82],a[ 83],a[ 84],a[ 85],a[ 86],a[ 87],a[ 88],a[ 89],a[ 90]);
case  91:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80],a[ 81],a[ 82],a[ 83],a[ 84],a[ 85],a[ 86],a[ 87],a[ 88],a[ 89],a[ 90],a[ 91]);
case  92:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80],a[ 81],a[ 82],a[ 83],a[ 84],a[ 85],a[ 86],a[ 87],a[ 88],a[ 89],a[ 90],a[ 91],a[ 92]);
case  93:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80],a[ 81],a[ 82],a[ 83],a[ 84],a[ 85],a[ 86],a[ 87],a[ 88],a[ 89],a[ 90],a[ 91],a[ 92],a[ 93]);
case  94:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80],a[ 81],a[ 82],a[ 83],a[ 84],a[ 85],a[ 86],a[ 87],a[ 88],a[ 89],a[ 90],a[ 91],a[ 92],a[ 93],a[ 94]);
case  95:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80],a[ 81],a[ 82],a[ 83],a[ 84],a[ 85],a[ 86],a[ 87],a[ 88],a[ 89],a[ 90],a[ 91],a[ 92],a[ 93],a[ 94],a[ 95]);
case  96:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80],a[ 81],a[ 82],a[ 83],a[ 84],a[ 85],a[ 86],a[ 87],a[ 88],a[ 89],a[ 90],a[ 91],a[ 92],a[ 93],a[ 94],a[ 95],a[ 96]);
case  97:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80],a[ 81],a[ 82],a[ 83],a[ 84],a[ 85],a[ 86],a[ 87],a[ 88],a[ 89],a[ 90],a[ 91],a[ 92],a[ 93],a[ 94],a[ 95],a[ 96],a[ 97]);
case  98:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80],a[ 81],a[ 82],a[ 83],a[ 84],a[ 85],a[ 86],a[ 87],a[ 88],a[ 89],a[ 90],a[ 91],a[ 92],a[ 93],a[ 94],a[ 95],a[ 96],a[ 97],a[ 98]);
case  99:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80],a[ 81],a[ 82],a[ 83],a[ 84],a[ 85],a[ 86],a[ 87],a[ 88],a[ 89],a[ 90],a[ 91],a[ 92],a[ 93],a[ 94],a[ 95],a[ 96],a[ 97],a[ 98],a[ 99]);
case 100:return f(a[  0],a[  1],a[  2],a[  3],a[  4],a[  5],a[  6],a[  7],a[  8],a[  9],a[ 10],a[ 11],a[ 12],a[ 13],a[ 14],a[ 15],a[ 16],a[ 17],a[ 18],a[ 19],a[ 20],a[ 21],a[ 22],a[ 23],a[ 24],a[ 25],a[ 26],a[ 27],a[ 28],a[ 29],a[ 30],a[ 31],a[ 32],a[ 33],a[ 34],a[ 35],a[ 36],a[ 37],a[ 38],a[ 39],a[ 40],a[ 41],a[ 42],a[ 43],a[ 44],a[ 45],a[ 46],a[ 47],a[ 48],a[ 49],a[ 50],a[ 51],a[ 52],a[ 53],a[ 54],a[ 55],a[ 56],a[ 57],a[ 58],a[ 59],a[ 60],a[ 61],a[ 62],a[ 63],a[ 64],a[ 65],a[ 66],a[ 67],a[ 68],a[ 69],a[ 70],a[ 71],a[ 72],a[ 73],a[ 74],a[ 75],a[ 76],a[ 77],a[ 78],a[ 79],a[ 80],a[ 81],a[ 82],a[ 83],a[ 84],a[ 85],a[ 86],a[ 87],a[ 88],a[ 89],a[ 90],a[ 91],a[ 92],a[ 93],a[ 94],a[ 95],a[ 96],a[ 97],a[ 98],a[ 99],a[100]);

		default:
				puts("Internal error (evlef02)");
				return 0.0;
		}
	}




Term rm_zero(Term t)
{
	int i;
	if(!is_compound(t))
		return t;

	

	for(i=1;i<=CompoundArity(t);i++)
	{
		SetCompoundArg(t,i,rm_zero(ConsumeCompoundArg(t,i)));
	}

	if(CompoundName(t)==OPR_POW && CompoundArity(t)==2 && 
			CompoundArg1(t)==NewInteger(0))
	{
		FreeAtomic(t); 
		return NewInteger(0);
	}
	
	if(CompoundName(t)==OPR_POW && CompoundArity(t)==2 && 
			is_integer(CompoundArg1(t)) && is_integer(CompoundArg2(t)))
	{
		int r=1;
		for(i=0;i<IntegerValue(CompoundArg2(t));i++)
			r*=(int)IntegerValue(CompoundArg1(t));
		return NewInteger(r);
	}
	
	if(CompoundName(t)==OPR_CARET && CompoundArity(t)==2 && 
			is_integer(CompoundArg1(t)) && is_integer(CompoundArg2(t)))
	{
		int r=1;
		for(i=0;i<IntegerValue(CompoundArg2(t));i++)
			r*=(int)IntegerValue(CompoundArg1(t));
		return NewInteger(r);
	}
	
	if(CompoundName(t)==OPR_MLT && CompoundArity(t)==2 && 
			(CompoundArg1(t)==NewInteger(0) || CompoundArg2(t)==NewInteger(0)))
	{
		FreeAtomic(t);
		return NewInteger(0);
	}

	if(CompoundName(t)==OPR_MLT && CompoundArity(t)==2 && 
			is_integer(CompoundArg1(t)) && is_integer(CompoundArg2(t)))
	{
		long int v=IntegerValue(CompoundArg1(t))*IntegerValue(CompoundArg2(t));
		FreeAtomic(t);
		return NewInteger(v);
	}
	if(CompoundName(t)==OPR_PLUS && CompoundArity(t)==2 &&CompoundArg2(t)==NewInteger(0))
	{
		Term t1=ConsumeCompoundArg(t,1);
		FreeAtomic(t);
		return rm_zero(t1);
	}
	if(CompoundName(t)==OPR_PLUS && CompoundArity(t)==2 &&CompoundArg1(t)==NewInteger(0))
	{
		Term t1=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		return rm_zero(t1);
	}
	if(CompoundName(t)==OPR_PLUS && CompoundArity(t)==2 &&
			is_integer(CompoundArg1(t)) && is_integer(CompoundArg2(t)))
	{
		int v=(int)IntegerValue(CompoundArg1(t))+(int)IntegerValue(CompoundArg2(t));
		FreeAtomic(t);
		return NewInteger(v);
	}
	if(CompoundName(t)==OPR_MINUS && CompoundArity(t)==2 &&CompoundArg2(t)==NewInteger(0))
	{
		Term t1=ConsumeCompoundArg(t,1);
		FreeAtomic(t);
		return rm_zero(t1);
	}
	if(CompoundName(t)==OPR_MINUS && CompoundArity(t)==2 &&CompoundArg1(t)==NewInteger(0))
	{
		Term t1=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		t1=rm_zero(t1);
		return t1==NewInteger(0)?NewInteger(0):MakeCompound1(OPR_MINUS,t1);
	}
	if(CompoundName(t)==OPR_MINUS && CompoundArity(t)==2 &&
			is_integer(CompoundArg1(t)) && is_integer(CompoundArg2(t)))
	{
		long int v=IntegerValue(CompoundArg1(t))-IntegerValue(CompoundArg2(t));
		FreeAtomic(t);
		return NewInteger(v);
	}
	return t;
}

