#include <string.h>
#include "lanhep.h"


Atomic 	A_RBRA,	A_RCET,	A_FBRA,	A_FCET,	A_QBRA,	A_QCET,
	A_POINT,A_COMMA,A_SECO,	A_QUOTE, A_RQUOTE, A_DQUOTE, A_QBRACET, A_FBRACET;

Atom OPR_MINUS, OPR_PLUS, OPR_MLT, OPR_DIV, OPR_CARET, OPR_POW, OPR_CLNMIN,
	OPR_COMMA, OPR_SECO, OPR_EQSIGN, OPR_COLON, OPR_USCORE, OPR_RARROW,
	OPR_EQEQSIGN, OPR_AND, OPR_OR, OPR_NOT, OPR_NOTEQ, OPR_DEFINED, OPR_DOLLAR, OPR_HASH;

Atom OPR_PARAMETER, OPR_PARTICLE, OPR_FIELD, OPR_GROUP, OPR_REPRES;
Atom OPR_SCALAR, OPR_VECTOR, OPR_SPINOR, OPR_TENSOR, OPR_SPINOR3, OPR_MASS, OPR_WIDTH,
      OPR_WILD, OPR_READ, OPR_EDIT, OPR_FUNCTION, OPR_CLASS, A_FORT_NAME;
Atom OPR_LET, OPR_LTERM, OPR_MODEL, OPR_CLEAR, OPR_SPECIAL;
Atom OP_PREFIX,OP_INFIX,OP_POSTFIX;
Atom A_COLOR, A_COLOR_LAMBDA, A_COLOR_F, A_COLOR_D, A_LEFT, A_RIGHT,
    A_GAUGE, A_GHOST, A_END, A_S_COLOR, A_COLOR_EPS, A_COLOR_K6, A_COLOR_L6, A_BRST_TRANSF,
	A_BRSTI_TRANSF, A_MOMENT_S, A_DERIV_S, A_DERIV, A_DERIV5, A_MOMENT_E, A_EPS_V;
Atom PROP_TYPE, PROP_INDEX;
Atom A_I, A_LORENTZ, A_MTERM, A_ALG1, A_ALG2, A_MOMENT, A_SQRT, A_COS, A_SIN,
	A_DELTA,A_GAMMA,A_GAMMA5, A_GAMMAP, A_GAMMAM, A_CURT, A_LHA,
	A_SQRT2, A_CC, A_GRASS, A_ANTI, A_ANTI2, A_GLUON, A_GG,
    A_VEV, A_ORTH_MATR, A_HERM_MATR, A_HERMC, A_CHNAME, A_TEXNAME, A_TEXLENGTH, 
	A_EM_CHARGE, A_CSPINOR, A_FANUM, A_SYMM, OPR_LOCAL,
    A_SET_ANGLE, A_TRIG_FU, A_INFINITESIMAL, A_DUMMY_PRM, A_KEEP_LETS;
Atom OPR_ALIAS, OPR_UNALIAS, OPR_WHERE, OPR_USE, OPR_COEFF, OPR_KEYS, OPR_DOIF, 
		OPR_DOELSEIF, OPR_DOELSE, OPR_ENDIF, OPR_OPTION, OPR_IN, OPR_EXTLIB,
		OPR_CPART, A_AUX;

Atom A_SIN, A_COS, A_TAN, A_ASIN, A_ACOS, A_ATAN, A_ATAN2, A_EXP, A_FABS,
		A_LOG, A_IF, A_EXT_FUNC, A_DFDFC, A_EE, A_GG, A_INTEGER;

/*
:- op(1150,fx,[parameter,model,particle,lterm,let,group,field]).
:- op(1050, fx,[vector,scalar,spinor]).
:- op(1140,xfx,where).
:- op(450,xfy,:).
:- op(800,fx,[mass,width,color,name]).
:- op(200,xfx,**).
*/

void InitAtoms1(void)
	{

	A_RBRA=NewAtom("(",0);
	A_RCET=NewAtom(")",0);
	A_FBRA=NewAtom("{",0);
	A_FCET=NewAtom("}",0);
	A_QBRA=NewAtom("[",0);
	A_QCET=NewAtom("]",0);
	A_FBRACET=NewAtom("{}",0);
	A_QBRACET=NewAtom("[]",0);
	A_POINT=NewAtom(".",0);
	OPR_COMMA=A_COMMA=NewAtom(",",0);
	OPR_SECO=A_SECO=NewAtom(";",0);
	A_QUOTE=NewAtom("'",0);
	A_RQUOTE=NewAtom("`",0);
	A_DQUOTE=NewAtom("\"",0);
	OPR_DOLLAR=NewAtom("$",0);
	OPR_HASH=NewAtom("#",0);
	A_I=NewAtom("i",0);
	A_LORENTZ=NewAtom("lorentz",0);
	A_MTERM=NewAtom("mterm",0);
	A_ALG1=NewAtom("alg1",0);
	A_ALG2=NewAtom("alg2",0);
	A_DELTA=NewAtom("delta",0);
	A_GAMMA=NewAtom("gamma",0);
	A_GAMMA5=NewAtom("gamma5",0);
	A_GAMMAP=NewAtom("(1+gamma5)/2",0);
	A_GAMMAM=NewAtom("(1-gamma5)/2",0);
	A_MOMENT=NewAtom("moment",0);
	OPR_WHERE=NewAtom("where",0);
	A_SQRT2=NewAtom("Sqrt2",0);
	A_SQRT=NewAtom("sqrt",0);
	A_CURT=NewAtom("Curt",0);
	A_COS=NewAtom("cos",0);
	A_SIN=NewAtom("sin",0);
	A_FABS=NewAtom("fabs",0);
	A_CC=NewAtom("cc",0);
	A_GRASS=NewAtom("grasman",0);
	A_ANTI=NewAtom("anti",0);
	A_ANTI2=NewAtom("anti2",0);
	A_GLUON=NewAtom("G",0);
	A_GG=NewAtom("GG",0);
	A_DFDFC=NewAtom("dfdfc",0);
	A_VEV=NewAtom("_vac_exp_val_",0);
	A_ORTH_MATR=NewAtom("OrthMatrix",0);
	A_HERM_MATR=NewAtom("HermMatrix",0);
	A_HERMC=NewAtom("hermconj",0);
	A_CHNAME=NewAtom("CompHEP name",0);
	A_TEXNAME=NewAtom("texname",0);
	A_TEXLENGTH=NewAtom("texlength",0);
	A_EM_CHARGE=NewAtom("em_charge",0);
	A_CSPINOR=NewAtom("cspinor",0);
	A_SET_ANGLE=NewAtom("_set_angle_",0);
	A_TRIG_FU=NewAtom("trig_func",0);
	A_INFINITESIMAL=NewAtom("infinitesimal",0);
	A_DUMMY_PRM=NewAtom("__dummy__parameter__",0);
	A_KEEP_LETS=NewAtom("keep_lets",0);
	A_BRST_TRANSF=NewAtom("brst_transform",0);
	A_BRSTI_TRANSF=NewAtom("brsti_transform",0);
	A_MOMENT_S=NewAtom("__moment__start__",0);
	A_DERIV_S=NewAtom("__deriv__start__",0);
	A_DERIV=NewAtom("deriv",0);
	A_DERIV5=NewAtom("deriv5",0);
	A_MOMENT_E=NewAtom("__moment__end__",0);
	A_EPS_V=NewAtom("epsv",0);
	OPR_LOCAL=NewAtom("local",0);

	OP_FX=NewAtom("fx",0);	OP_FY=NewAtom("fy",0);
	OP_XF=NewAtom("xf",0);	OP_YF=NewAtom("yf",0);
	OP_XFX=NewAtom("xfx",0);
	OP_YFX=NewAtom("yfx",0);OP_XFY=NewAtom("xfy",0);
	OP_INFIX=NewAtom("infix",0);
	OP_PREFIX=NewAtom("prefix",0);
	OP_POSTFIX=NewAtom("postfix",0);

	A_LEFT=NewAtom("left",0);
	A_RIGHT=NewAtom("right",0);
	A_COLOR=NewAtom("color",0);
	A_S_COLOR=NewAtom("split_color",0);
	A_COLOR_LAMBDA=NewAtom("__color_lambda__",0);
	A_COLOR_F=NewAtom("__color_stru_cont__",0);
	A_COLOR_D=NewAtom("__color_stru_symm__",0);
	A_COLOR_EPS=NewAtom("__color_eps__",0);
	A_COLOR_K6=NewAtom("__color_k6__",0);
	A_COLOR_L6=NewAtom("__color_l6__",0);
	A_GAUGE=NewAtom("gauge",0);
	A_GHOST=NewAtom("ghost",0);
	A_END=NewAtom("end",0);
	PROP_TYPE=NewAtom("type",0);
	PROP_INDEX=NewAtom("indices",0);
	OPR_PLUS=NewAtom("+",0);
	OPR_MINUS=NewAtom("-",0);
	OPR_MLT=NewAtom("*",0);
	OPR_DIV=NewAtom("/",0);
	OPR_POW=NewAtom("**",0);
	OPR_CARET=NewAtom("^",0);
	OPR_EQSIGN=NewAtom("=",0);
	OPR_EQEQSIGN=NewAtom("==",0);
	OPR_NOTEQ=NewAtom("!=",0);
	OPR_NOT=NewAtom("!",0);
	OPR_OR=NewAtom("||",0);
	OPR_AND=NewAtom("&&",0);
	OPR_DEFINED=NewAtom("defined",0);
	OPR_COLON=NewAtom(":",0);
	OPR_USCORE=NewAtom("_",0);
	OPR_RARROW=NewAtom("->",0);
	OPR_PARAMETER=NewAtom("parameter",0);
	OPR_PARTICLE=NewAtom("particle",0);
	OPR_FIELD=NewAtom("field",0);
	OPR_GROUP=NewAtom("group",0);
	OPR_REPRES=NewAtom("repres",0);
	OPR_READ=NewAtom("read",0);
	OPR_EDIT=NewAtom("edit",0);
	OPR_SCALAR=NewAtom("scalar",0);
	OPR_VECTOR=NewAtom("vector",0);
	OPR_SPINOR=NewAtom("spinor",0);
	OPR_SPINOR3=NewAtom("spinor3",0);
	OPR_TENSOR=NewAtom("tensor",0);
	OPR_LET=NewAtom("let",0);
	OPR_LTERM=NewAtom("lterm",0);
	OPR_MASS=NewAtom("mass",0);
	OPR_WIDTH=NewAtom("width",0);
	OPR_MODEL=NewAtom("model",0);
	OPR_CLEAR=NewAtom("clear",0);
	OPR_SPECIAL=NewAtom("special",0);
	OPR_WILD=NewAtom("wild",0);
	OPR_FUNCTION=NewAtom("function",0);
	OPR_ALIAS=NewAtom("alias",0);
	OPR_UNALIAS=NewAtom("unalias",0);
	OPR_USE=NewAtom("use",0);
	OPR_COEFF=NewAtom("coeff",0);
	OPR_KEYS=NewAtom("keys",0);
	OPR_DOIF=NewAtom("do_if",0);
	OPR_DOELSEIF=NewAtom("do_else_if",0);
	OPR_DOELSE=NewAtom("do_else",0);
	OPR_ENDIF=NewAtom("end_if",0);
	OPR_OPTION=NewAtom("option",0);
	OPR_CLASS=NewAtom("class",0);
	OPR_IN=NewAtom("in",0);
	OPR_EXTLIB=NewAtom("extlib",0);
	OPR_CPART=NewAtom("cpart",0);
	A_AUX=NewAtom("aux",0);

	SetAtomProperty(OPR_UNALIAS,OPR_UNALIAS,NewInteger(1));
	SetAtomProperty(NewAtom("GetProperties",0),OPR_UNALIAS,NewInteger(1));
	
	A_SIN=NewAtom("sin",0);
	A_COS=NewAtom("cos",0);
	A_TAN=NewAtom("tan",0);
	A_ASIN=NewAtom("asin",0);
	A_ACOS=NewAtom("acos",0);
	A_ATAN=NewAtom("atan",0);
	A_ATAN2=NewAtom("atan2",0);
	A_EXP=NewAtom("exp",0);
	A_LOG=NewAtom("log",0);
	A_IF=NewAtom("if",0);
	A_EXT_FUNC=NewAtom("external_func",0);
	A_FANUM=NewAtom("Feynart",0);
	A_SYMM=NewAtom("symmetry",0);
	A_LHA=NewAtom("lha",0);
	A_EE=NewAtom("EE",0);
	A_GG=NewAtom("GG",0);
	A_INTEGER=NewAtom("integer",0);
	A_FORT_NAME=NewAtom("fortran_name",0);

	SetOperator(OPR_DOLLAR,OP_FX,10);
	SetOperator(OPR_PLUS,OP_YFX,500);
	SetOperator(OPR_MINUS,OP_YFX,500);
	SetOperator(OPR_MLT,OP_YFX,400);
	SetOperator(OPR_DIV,OP_YFX,400);
	SetOperator(OPR_PLUS,OP_FX,10);
	SetOperator(OPR_MINUS,OP_FX,499);
	SetOperator(OPR_CARET,OP_XFY,200);
	SetOperator(OPR_USCORE,OP_XFY,200);
	SetOperator(OPR_POW,OP_XFX,100);
	SetOperator(OPR_COMMA,OP_XFY,1100);
	SetOperator(OPR_HASH,OP_XFY,11);
	SetOperator(OPR_SECO,OP_XFY,1200);
	SetOperator(OPR_EQSIGN,OP_XFY,700);
	SetOperator(OPR_EQEQSIGN,OP_XFY,700);
	SetOperator(OPR_NOT,OP_FX,600);
	SetOperator(OPR_NOTEQ,OP_XFY,700);
	SetOperator(OPR_AND,OP_XFY,750);
	SetOperator(OPR_OR,OP_XFY,800);
	SetOperator(OPR_RARROW,OP_XFY,700);
	SetOperator(OPR_COLON,OP_XFX,600);
	SetOperator(OPR_PARAMETER,OP_FX,1250);
	SetOperator(OPR_PARTICLE,OP_FX,1250);
	SetOperator(OPR_FIELD,OP_FX,1250);
	SetOperator(OPR_GROUP,OP_FX,1250);
	SetOperator(OPR_REPRES,OP_FX,1250);
	SetOperator(OPR_READ,OP_FX,1250);
	SetOperator(OPR_EDIT,OP_FX,1250);
	SetOperator(OPR_SCALAR,OP_FX,1150);
	SetOperator(OPR_VECTOR,OP_FX,1150);
	SetOperator(OPR_SPINOR,OP_FX,1150);
	SetOperator(OPR_TENSOR,OP_FX,1150);
	SetOperator(OPR_SPINOR3,OP_FX,1150);
	SetOperator(OPR_LET,OP_FX,1275);
	SetOperator(OPR_ALIAS,OP_FX,1275);
	SetOperator(OPR_UNALIAS,OP_FX,1275);
	SetOperator(OPR_LOCAL,OP_FX,1275);
	SetOperator(OPR_LTERM,OP_FX,1275);
	SetOperator(A_LHA,OP_FX,1275);
	SetOperator(OPR_WHERE,OP_XFX,/*1225*/ 2600);
	SetOperator(OPR_MASS,OP_FX,800);
	SetOperator(OPR_WIDTH,OP_FX,800);
	SetOperator(A_TEXNAME,OP_FX,800);
	SetOperator(NewAtom("atexname",0),OP_FX,800);
	SetOperator(NewAtom("echarge",0),OP_FX,800);
	SetOperator(NewAtom("fullname",0),OP_FX,800);
	SetOperator(A_AUX,OP_FX,800);
	SetOperator(OPR_MODEL,OP_FX,1250);
	SetOperator(OPR_CLEAR,OP_FX,1250);
	SetOperator(OPR_SPECIAL,OP_FX,1250);
	SetOperator(OPR_WILD,OP_FX,800);
	SetOperator(OPR_USE,OP_FX,1250);
	SetOperator(NewAtom("eval",0),OP_FX,1250);
	SetOperator(NewAtom("transform",0),OP_FX,1250);
	SetOperator(NewAtom("prtcformat",0),OP_FX,1250);
	SetOperator(NewAtom("prtcproperty",0),OP_FX,1250);
	SetOperator(NewAtom("prtcprop",0),OP_FX,1250);
	SetOperator(A_BRST_TRANSF,OP_FX,1250);
	SetOperator(A_BRSTI_TRANSF,OP_FX,1250);
	SetOperator(NewAtom("angle",0),OP_FX,1250);
	SetOperator(A_INFINITESIMAL,OP_FX,1250);
	SetOperator(OPR_COEFF,OP_FX,1250);
	SetOperator(OPR_KEYS,OP_FX,1250);
	SetOperator(OPR_DOIF,OP_FX,1250);
	SetOperator(OPR_DOELSEIF,OP_FX,1250);
	SetOperator(OPR_OPTION,OP_FX,1250);
	SetOperator(A_KEEP_LETS,OP_FX,1250);
	SetOperator(NewAtom("opt_lets",0),OP_FX,1250);
	SetOperator(NewAtom("ued_5th",0),OP_FX,1250);
	SetOperator(NewAtom("write",0),OP_FX,1250);
	SetOperator(NewAtom("display",0),OP_FX,1250);
	SetOperator(OPR_CLASS,OP_FX,1250);
	SetOperator(OPR_IN,OP_XFY,2500);
	SetOperator(OPR_EXTLIB,OP_FX,1250);
	SetOperator(OPR_CPART,OP_FX,1250);
	
	SetAtomProperty(A_COLOR_LAMBDA,A_COLOR,A_COLOR_LAMBDA);
	SetAtomProperty(A_COLOR_F,A_COLOR,A_COLOR_F);
	
	SetAtomProperty(A_VEV,A_ANTI,A_VEV);
	
		{
		Term t;
		t=MakeCompound(OPR_PARTICLE,8);
		SetCompoundArg(t,1,A_VEV);
		SetCompoundArg(t,2,A_VEV);
		SetCompoundArg(t,4,NewInteger(0));
		SetCompoundArg(t,8,OPR_MLT);
		SetAtomProperty(A_VEV,PROP_TYPE,t);
		}
	
	}


char *ATOM_buffers[256];
int  ATOM_buffill[256], ATOM_count=0, ATOM_stcount=0, ATOM_tscount=0;

void InitAtoms(void)
	{
	int i;
	for(i=0;i<256;i++)
		{
		ATOM_buffill[i]=0;
		ATOM_buffers[i]=NULL;
		}
	InitAtoms1();
	}
		

Atom NewAtom(char *s, int len)
	{
	int i,j,ret;
	if(len==0) len=(int)strlen(s);
	
	/*   searching..  */
	i=0;
	ATOM_stcount++;
	while(ATOM_buffers[i]!=NULL)
		{
		char *abuff;
		abuff=ATOM_buffers[i];
		j=0;
		while(j<ATOM_buffill[i])
			{
			
			int k=0;
			j+=(int)sizeof(List);
			while(s[k]==abuff[j+k] && k<len) k++;
			if(k==len && abuff[j+k]==0)
				return j-(int)sizeof(List)+i*0x10000+(1L<<56);
			j+=k;
			while(abuff[j]!=0) j++; j++;
			while(j%8)  j++;
			/*
			if(strcmp(s,ATOM_buffers[i]+j)==0)
				return j+i*0x10000+0x01000000;
			j+=strlen(ATOM_buffers[i]+j)+1;
			*/
			}
		i++;
		}
	/*  inserting ... */
	ATOM_count++;
	i=0;
	while(i<256 && ATOM_buffers[i]!=NULL && ATOM_buffill[i]+len+(int)sizeof(List)>0xffff)
		i++;
	if(i==256)
		{  puts("Internal error (too much atoms or memory lack)."); exit(0); }
	if(ATOM_buffers[i]==NULL)
		{
		ATOM_buffers[i]=(char *)malloc(0x10000);
		if(ATOM_buffers[i]==NULL)
			{  puts("Internal error (memory lack for atoms)."); exit(0); }
		}
	{ List *tmp; tmp=(List *)(ATOM_buffers[i]+ATOM_buffill[i]); *tmp=NewList(); }
	strncpy(ATOM_buffers[i]+ATOM_buffill[i]+(int)sizeof(List),s,(unsigned int)len);
	*(ATOM_buffers[i]+ATOM_buffill[i]+len+(int)sizeof(List))=0;
	ret=ATOM_buffill[i];
	ATOM_buffill[i]+=len+1+(int)sizeof(List);
	while(ATOM_buffill[i]%8) ATOM_buffill[i]++;
	return ret+i*0x10000+(1L<<56);
	
	}
	

static	char avcbuf[64];

		
char *AtomValue(Atom a)
	{
	int  bno,bpos;
	
	if(a==0)
		return "_";
	if(!is_atom(a))
	{
		if(is_float(a))
		{
			sprintf(avcbuf,"%e",FloatValue(a));
			return avcbuf;
		}
		return "?????";
	}
	
	bno=a&0xff0000;
	bno/=0x10000;
	bpos=a&0xffff;
	if(ATOM_buffers[bno]==NULL)
		{
		printf("Internal error (illegal atom '");
		if(!is_atom(a))
			WriteTerm(a);
		else
			printf("%x",a);
		printf("').\n");
		return "??";
		}
	ATOM_tscount++;
	return ATOM_buffers[bno]+bpos+sizeof(List);
	}

void SetAtomProperty(Atom a, Atom type, Term value)
     {
     int  bno,bpos;
     List *lp,ll;
     Term tt;
     bno=a&0xff0000;
     bno/=0x10000;
     bpos=a&0xffff;
     if(ATOM_buffers[bno]==NULL)
	     {  puts("Internal error (illegal atom)."); exit(0); }
     lp=(List *)(ATOM_buffers[bno]+bpos);
     ll=*lp;
     while(!is_empty_list(ll))
         {
         tt=ListFirst(ll);
         if(CompoundName(tt)==type)
             {
             Term kk;
             kk=ConsumeCompoundArg(tt,1);
             FreeAtomic(kk);
             SetCompoundArg(tt,1,value);
             return;
             }
         ll=ListTail(ll);
         }
     tt=MakeCompound(type,1);
     SetCompoundArg(tt,1,value);
     *lp=AppendLast(*lp,tt);
     }

Term GetAtomProperty(Atom a, Atom type)
     {
     int  bno,bpos;
     List *lp,ll;
     Term tt;
     if(!is_atom(a)) return 0;
     bno=a&0xff0000;
     bno/=0x10000;
     bpos=a&0xffff;
     if(ATOM_buffers[bno]==NULL)
	     {  puts("Internal error (illegal atom) (GetAtomProp)."); exit(0); }
     lp=(List *)(ATOM_buffers[bno]+bpos);
     ll=*lp;
     while(!is_empty_list(ll))
         {
         tt=ListFirst(ll);
         if(CompoundName(tt)==type)
             return CompoundArg1(tt);
         ll=ListTail(ll);
         }
     return (Atom)0;
     }

void RemoveAtomProperty(Atom a, Atom type)
     {
     int  bno,bpos;
     List *lp,ll;
     Term tt;
     bno=a&0xff0000;
     bno/=0x10000;
     bpos=a&0xffff;
     if(ATOM_buffers[bno]==NULL)
	     {  puts("Internal error (illegal atom)."); exit(0); }
     lp=(List *)(ATOM_buffers[bno]+bpos);
     ll=*lp;
     while(!is_empty_list(ll))
         {
         tt=ListFirst(ll);
         if(CompoundName(tt)==type)
             {
             *lp=CutFromList(*lp,ll);
             return;
             }
         ll=ListTail(ll);
         }
     return;
     }
     
Term SetProperty(Term t, Term ind)
	{
	if(!is_compound(t) || CompoundArity(t)!=3 ||
		!is_atom(CompoundArg1(t)) ||
		!is_atom(CompoundArg2(t)))
		{
		ErrorInfo(390);
		puts("bad arguments in SetProperty statement");
		}
	SetAtomProperty(CompoundArg1(t),CompoundArg2(t),
		ConsumeCompoundArg(t,3));
	FreeAtomic(t);
	return 0;
	}

Term GetProperties(Term t, Term ind)
	{
	int  bno,bpos;
	Atom a;
     List *lp;
	 WriteTerm(t);puts("");
     a=ConsumeCompoundArg(t,1);
     FreeAtomic(t);
     if(!is_atom(a))
     	return 0;
     bno=a&0xff0000;
     bno/=0x10000;
     bpos=a&0xffff;
     if(ATOM_buffers[bno]==NULL)
	     {  puts("Internal error (illegal atom)."); exit(0); }
     lp=(List *)(ATOM_buffers[bno]+bpos);
	 printf("%lx %lx\n",lp,*lp);
	return CopyTerm(*lp);
	}

List AtomPropertiesList(Atom a)
{
	int  bno,bpos;
	List *lp;
     if(!is_atom(a))
     	return 0;
     bno=a&0xff0000;
     bno/=0x10000;
     bpos=a&0xffff;
     if(ATOM_buffers[bno]==NULL)
	     {  puts("Internal error (illegal atom)."); exit(0); }
	lp=(List *)(ATOM_buffers[bno]+bpos);
	return *lp;
}
