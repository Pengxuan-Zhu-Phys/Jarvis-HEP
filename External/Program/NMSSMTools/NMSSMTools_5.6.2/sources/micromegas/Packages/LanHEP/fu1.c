#include "lanhep.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/stat.h>


extern FILE *log_file;
extern int VerbMode;
char *EditorName = "joe";
extern char *InputDirectory;
void clear_infi(Atom);
static char cbuf[80];

Term ProcessEdit(Term t, Term ind)
	{
	Term arg;
	arg=CompoundArg1(t);
	if(!is_atom(arg))
		{
		ErrorInfo(205);
		printf(" \'"); WriteTerm(arg);
		printf("\' is not appropriate file name.\n");
		FreeAtomic(t);
		return 0;
		}
	FreeAtomic(t);
	sprintf(cbuf,"%s %s",EditorName,AtomValue(arg));
	system(cbuf);
	return 0;
	}

static List current_keys = 0;

static void set_keys(Term k)
{
	Term k1;
	if(CompoundArity(k)!=1)
	{
		ErrorInfo(205);
		puts("keys: wrong syntax.");
		FreeAtomic(k);
		return ;
	}
	
	k1=ConsumeCompoundArg(k,1);
	FreeAtomic(k);
	if(is_atom(k1) && AtomValue(k1)[0]=='?')
	{
	  DumpList(current_keys);
	  return;
	}
	k=CommaToList(k1);
	for(k1=k;k1;k1=ListTail(k1))
	{
		Term k2,n,v;
		List l;
		k2=ListFirst(k1);
		if(!is_compound(k2)||CompoundArity(k2)!=2 || CompoundName(k2)!=OPR_EQSIGN)
		{
			ErrorInfo(206);
			printf("keys: incorrect syntax in '");
			WriteTerm(k2);
			puts("'.");
			continue;
		}
		n=CompoundArg1(k2);
		v=CompoundArg2(k2);
		if(!is_atom(n))
		{
			ErrorInfo(206);
			printf("keys: illegal key name '");
			WriteTerm(n);
			puts("'.");
			continue;
		}
		if(!is_atom(v) && !is_integer(v))
		{
			ErrorInfo(207);
			printf("keys: illegal key value '");
			WriteTerm(v);
			puts("'.");
			continue;
		}
		for(l=current_keys;l;l=ListTail(l))
			if(CompoundArg1(ListFirst(l))==n)
				break;
		if(l)
		{
			
			if(CompoundName(ListFirst(l))==A_I)
				continue;
			
			if(v!=CompoundArg2(ListFirst(l)))
			{
				WarningInfo(206);
				printf("keys: redefinition of '");
				WriteTerm(n);
				puts("'.");
			}
			SetCompoundArg(ListFirst(l),2,v);
			FreeAtomic(k2);
		}
		else
			current_keys=AppendLast(current_keys,k2);
	}
	RemoveList(k);
	
}

void SetKeyFromArg(char *arg)
{
	int pos=0;
	Atom key,val;
	while(arg[pos] && arg[pos]!='=')
		pos++;
	if(arg[pos]==0)
	{
		printf("Error: wrong key definition '%s' in command line.\n",arg);
		return;
	}
	key=NewAtom(arg,pos);
	if(isdigit(arg[pos+1]))
	{
		int n;
		sscanf(arg+pos+1,"%d",&n);
		val=NewInteger(n);
	}
	else
		val=NewAtom(arg+1+pos,0);
	
	current_keys=AppendLast(current_keys,MakeCompound2(A_I,key,val));
}

static Term read_to_else(void)
{
	int lev=1;
	Term a;
	while((a=ReadTerm()))
	{
		
		if(is_compound(a) && CompoundName(a)==OPR_DOIF)
		{
			lev++;
			continue;
		}
		
		if(a==OPR_ENDIF)
		{
			lev--;
			if(lev==0)
				return a;
			continue;
		}
		if(lev==1 && (a==OPR_DOELSE || 
				(is_compound(a) && CompoundName(a)==OPR_DOELSEIF)))
			return a;
	}
	
	ErrorInfo(209);
	printf("do_if statement not followed by 'end_if' one.\n");
	return 0;
	
}

static Term key_val(Term k)
{
	List l;
	for(l=current_keys;l;l=ListTail(l))
		if(CompoundArg1(ListFirst(l))==k)
			return CompoundArg2(ListFirst(l));
	return 0;
}
			
static int read_eval(Term t)
{
	List l;
	
	
	if(!is_compound(t))
	{
		ErrorInfo(208);
		printf("illegal condition '");
		WriteTerm(t);
		printf("' in do_if/do_else_if statement.\n");
		return 0;
	}
	
	if(CompoundName(t)==OPR_EQSIGN)
	{
		WarningInfo(211);
		printf("assuming '==' in '");
		WriteTerm(t);
		puts("' instead of '='.");
	}
	
	if(CompoundName(t)==OPR_EQSIGN || CompoundName(t)==OPR_EQEQSIGN
	  	|| CompoundName(t)==OPR_NOTEQ)
	{
		Term a1, a2, k1, k2;
		a1=CompoundArg1(t);
		a2=CompoundArg2(t);
		int inv=(CompoundName(t)==OPR_NOTEQ);
		k1=key_val(a1);
		k2=key_val(a2);
		
		if(k1==0)
		{
			ErrorInfo(0);
			printf("do_if: ");
			WriteTerm(a1);
			puts(" is not defined as key.");
			return 0;
		}
		/*if(k2==0)*/
			return inv?(k1!=a2):(k1==a2);
		/*if(k1==0)
			return inv?(k2!=a1):(k2==a1);
		return inv?(k1!=k2):(k1==k2);*/
	}
	
	if(CompoundName(t)==OPR_NOT)
	{
		int r=read_eval(CompoundArg1(t));
		return !r;
	}
	
	if(CompoundName(t)==OPR_OR)
	{
		int r1=read_eval(CompoundArg1(t)), r2=read_eval(CompoundArg2(t));
		return r1||r2;
	}
	
	if(CompoundName(t)==OPR_AND)
	{
		int r1=read_eval(CompoundArg1(t)), r2=read_eval(CompoundArg2(t));
		return r1&&r2;
	}
	
	if(CompoundName(t)==OPR_DEFINED)
	{
		Term r=key_val(CompoundArg1(t));
		return r!=0;
	}
	
	{
		ErrorInfo(208);
		printf("illegal condition '");
		WriteTerm(t);
		printf("' in do_if/do_else_if statement.\n");
		return 0;
	}
	
}

static Term read_to_endif(void)
{
	Term a;
	
	do
		a=read_to_else();
	while(a && a!=OPR_ENDIF);
	
	return a;
}

static Term read_file(int rec, char *file)
{
	
	Term a;
	int li=1;

	WritePrompt(li);
	while((a=ReadTerm())!=0)
		{
		if(a==A_END)
		{
			if(rec)
			{
				ErrorInfo(209);
				printf("'end' statement before 'end_if'\n");
				return 0;
			}
			else
				return 1;
		}
		
		if(a==OPR_ENDIF)
		{
			if(!rec)
			{
				ErrorInfo(209);
				printf("'end_if' statement without preceding 'do_if'\n");
				return 0;
			}
			else
				return 1;
		}
		
		if(a==OPR_DOELSE || (is_compound(a) && CompoundName(a)==OPR_DOELSEIF))
		{
			if(!rec)
			{
				ErrorInfo(209);
				printf("'do_else' statement without preceding 'do_if'\n");
				return 0;
			}
			else
				return read_to_endif();
		}
		
		if(a==OPR_KEYS || a==OPR_DOELSE || a==OPR_DOIF)
		{
			ErrorInfo(205);
			printf("'%s' statement needs arguments.\n",AtomValue(a));
			continue;
		}
		
		
		if(is_compound(a) && CompoundName(a)==OPR_KEYS)
		{
		  	set_keys(a);
			continue;
		}
		
		if(is_compound(a) && CompoundName(a)==OPR_DOIF)
		{
			int cimp, resc;
			cimp=CurrentInputLine();
			
			if(is_compound(CompoundArg1(a)) &&
				CompoundName(CompoundArg1(a))==OPR_COMMA)
			{
				List ll=CommaToList(ConsumeCompoundArg(a,1));
				if(ListLength(ll)>3)
				{
					ErrorInfo(0);
					puts("too many arguments in do_if.");
					continue;
				}
				if(read_eval(ListFirst(ll)))
					ProcessTerm(ListNth(ll,2));
				else
					if(ListLength(ll)==3)
						ProcessTerm(ListNth(ll,3));
				continue;
			}
			
			/*WarningInfo(0);printf("if "); WriteTerm(CompoundArg1(a));*/
			resc=read_eval(CompoundArg1(a));
			/*printf(" -> %d\n",resc);*/
			
			if(resc)
			{
				Term a1;
				a1=read_file(1,file);
				if(a1==0)
				{
					printf("\tdo_if block started at line %d.\n",cimp);
					return 0;
				}
				continue;
			}
			else
			{
				Term a;
				
				do
					a=read_to_else();
				while(a && a!=OPR_DOELSE && a!=OPR_ENDIF &&
						!(is_compound(a) && CompoundName(a)==OPR_DOELSEIF &&
						read_eval(CompoundArg1(a))));
				
				if(a==0)
				{
					printf("\tdo_if block started at line %d.\n",cimp);
					return 0;
				}
				if(a==OPR_ENDIF)
					continue;
				if(a==OPR_DOELSE || is_compound(a))
					{
						if(read_file(1,file)==0)
						{
							printf("\tdo_if block started at line %d.\n",cimp);
							return 0;
						}
						continue;
					}
				puts("Internal error (rfif)");
			}
		}
		
		
		if(a!=1)
			{
			if(log_file!=NULL && file==NULL)
				{
				fWriteTerm(log_file,a);
				fprintf(log_file,".\n");
				}
			ProcessTerm(a);
			}
/*		if(IsTermInput())
			puts("");*/
		li++;
		WritePrompt(li);
		}
	
	if(rec)
	{
		ErrorInfo(208);
		printf("file ended before 'end_if' statement.\n");
		return 0;
	}
	return 1;
}

char *eff_infile=0;
int doinitfile=0;

void ReadFile(char *file)
	{
	
	char ofbuf[1028];
	
	time_t  start_time;
	
	start_time=time(NULL);
	if(file==NULL)
		{
		SetInputFile(NULL);
		goto prc;
		}

	strcpy(ofbuf,file);
	if(SetInputFile(ofbuf)!=0)
		goto prc;


	strcat(ofbuf,".mdl");
	if(SetInputFile(ofbuf)!=0)
		goto prc;
	
	strcpy(ofbuf,file);
	strcat(ofbuf,".lh");
	if(SetInputFile(ofbuf)!=0)
		goto prc;
	
	strcpy(ofbuf,file);
	strcat(ofbuf,".lhep");
	if(SetInputFile(ofbuf)!=0)
		goto prc;

	strcpy(ofbuf,InputDirectory);
	strcat(ofbuf,file);
	if(SetInputFile(ofbuf)!=0)
		goto prc;

	strcat(ofbuf,".mdl");
	if(SetInputFile(ofbuf)!=0)
		goto prc;

	strcpy(ofbuf,InputDirectory);
	strcat(ofbuf,file);
	strcat(ofbuf,".lh");
	if(SetInputFile(ofbuf)!=0)
		goto prc;
	
	strcpy(ofbuf,InputDirectory);
	strcat(ofbuf,file);
	strcat(ofbuf,".lhep");
	if(SetInputFile(ofbuf)!=0)
		goto prc;

	printf("Error: file %s not found\n",file);
	return;


prc:
	
	if(eff_infile==0 && !doinitfile)
	{
		eff_infile=malloc(1038);
		if(eff_infile==NULL)
		{
			puts("internal error (fu1mef).");
			exit(0);
		}
		if(file==NULL)
			strcpy(eff_infile,"(stdin)");
		else
		{
		if(ofbuf[0]!='/')
		{
			getcwd(eff_infile,1028);
			if(!(ofbuf[0]=='.'&&ofbuf[1]=='/'))
				strcat(eff_infile,"/");
			strcat(eff_infile,(ofbuf[0]=='.'&&ofbuf[1]=='/')?ofbuf+1:ofbuf);
		}
		else
			strcpy(eff_infile,ofbuf);
		}
	}
		
	read_file(0,file);
		
	if(file!=NULL && strcmp(file,"lhep.rc") && !doinitfile)
		{
		long int min,sec;
		time_t ttt;
		ttt=time(NULL)-start_time;
		min=ttt/60;
		sec=ttt-min*60;
		printf("File %s processed, ",file);
		if(min!=0)
			printf("%ld min ",min);
		printf("%ld sec.\n",sec);
		}
	CloseInputFile();
	
	}



Term ProcessRead(Term t, Term ind)
	{
	Atom af;
	List l;
	char pubuf[64];
	Term t1;
	
	if(!is_compound(t) || CompoundArity(t)!=1)
	{
		ErrorInfo(201);
		puts("wrong arg in 'read' statement");
		return 0;
	}
	
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	t1=l=CommaToList(t1);
	for(;t1;t1=ListTail(t1))
		{
			af=ListFirst(t1);
			if(!is_atom(af))
			{
				ErrorInfo(202);
				printf("read: bad file name '");
				WriteTerm(af);
				printf("'\n");
				continue;
			}
			strcpy(pubuf,AtomValue(af));
			ReadFile(pubuf);
		}

	FreeAtomic(l);
	
	return 0;		
	}



Term ProcessModel(Term t, Term ind)
	{
	Term arg;
	arg=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	if(is_atom(arg))
		{
		ModelName=AtomValue(arg);
		if(ModelNumber==0)
			ModelNumber=1;
		goto outdir;
		}
	if(is_compound(arg) && FunctorName(CompoundFunctor(arg))==OPR_DIV
		&& FunctorArity(CompoundFunctor(arg))==2
		&& is_atom(CompoundArg1(arg))
		&& is_integer(CompoundArg2(arg)))
			{
			ModelName=AtomValue(CompoundArg1(arg));
			ModelNumber=(int)IntegerValue(CompoundArg2(arg));
			FreeAtomic(arg);
			goto outdir;
			}
	
    
	ErrorInfo(207);
	printf(" illegal arguments in 'model' statement.\n ");
	FreeAtomic(arg);
	return 0;
	outdir:
	if(UFOutput && OutputDirectory==0)
        {
            OutputDirectory=malloc(1000);
            strcpy(OutputDirectory,ModelName);
            for(int i=0;i<1000&&OutputDirectory[i];i++)
            {
                if(OutputDirectory[i]=='(' || OutputDirectory[i]==')' ||
                    OutputDirectory[i]==' ' || OutputDirectory[i]==',')
                    OutputDirectory[i]='_';
            }
            if(access(OutputDirectory,F_OK))
            {
                if(mkdir(OutputDirectory,0777))
                {
                    perror(OutputDirectory);
                    OutputDirectory=0;
                }
                else printf("Folder %s is created.\n",OutputDirectory);
            }
        }
    return 0;
	}

Term GetProp(Term t, Term ind)
	{
	Term at,prp,ll;
	if(CompoundArity(t)!=2)
		{
		ErrorInfo(208);
		printf("Wrong arguments number in 'getprop' call.\n");
		return 0;
		}
	at=ConsumeCompoundArg(t,1);
	prp=ConsumeCompoundArg(t,2);
	if(!is_atom(at) || !is_atom(prp))
		{
		ErrorInfo(208);
		printf("Arguments in 'getprop' call must be atoms.\n");
		return 0;
		}
	FreeAtomic(t);
	ll=GetAtomProperty(at,prp);
	return CopyTerm(ll);
	}


Term ProcessDisplay(Term t, Term ind)
	{
	Term t1;
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	DisplayTerm(t1);
	puts("");
	FreeAtomic(t1);
	return 0;
	}
	
	
Term ProcessWrite(Term t, Term ind)
	{
	Term t1;
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	if(is_atom(t1))
		puts(AtomValue(t1));
	else
	{
		WriteTerm(t1);
		puts("");
	}
	FreeAtomic(t1);
	return 0;
	}
		

Term ProcessDate(Term t, Term i)
	{
	time_t tttt;
	time(&tttt);
	return NewAtom(ctime(&tttt),0);
	}


static void p_clear(Term t1)
	{
	
	if(is_compound(t1) && CompoundName(t1)==OPR_COMMA)
		{
		Term a1,a2;
		a1=ConsumeCompoundArg(t1,1);
		a2=ConsumeCompoundArg(t1,2);
		FreeAtomic(t1);
		p_clear(a1);
		p_clear(a2);
		return;
		}
		
	if(is_compound(t1) && CompoundName(t1)==OPR_MINUS &&
			CompoundArity(t1)==2 && is_parameter(CompoundArg1(t1)) &&
			is_parameter(CompoundArg2(t1)))
	{
		List l, lc, lp=all_param_list();
		int f=0;
		lc=0;
		for(l=lp;l;l=ListTail(l))
		{
			Atom p=CompoundName(ListFirst(l));
			if(p==CompoundArg1(t1))
				f=1;
			if(f) lc=AppendLast(lc,p);
			if(p==CompoundArg2(t1))
				f=0;
		}
		if(lc==0)
		{
			WarningInfo(0);
			printf("clear ");WriteTerm(t1);puts(": no particles found.");
			return;
		}
		if(f)
		{
			WarningInfo(0);
			printf("clear ");WriteTerm(t1);puts(": second particle not found.");
			return;
		}
		for(l=lc;l;l=ListTail(l)) ClearParameter(ListFirst(l));
		FreeAtomic(lc);
		FreeAtomic(t1);
		return;
	}

	if(!is_atom(t1))
		{
		ErrorInfo(209);
		puts("wrong argument in call to 'clear'.");
		return ;
		}
	if(is_particle(t1,NULL))
		{
		ClearParticle(t1);
		return ;
		}
	if(is_let(t1,NULL))
		{
		ClearLet(t1);
		return ;
		}
	if(GetAtomProperty(t1,A_INFINITESIMAL))
	{
		clear_infi(t1);
		return;
	}
	if(is_parameter(t1))
		{
		ClearParameter(t1);
		return ;
		}
	if(is_special(t1,NULL))
		{
		ClearSpecial(t1);
		return ;
		}
	if(is_group(t1,NULL))
		{
		ClearGroup(t1);
		return ;
		}
	ErrorInfo(209);
	printf("don't know how to clear %s.\n",AtomValue(t1));
	return ;
	}
	
Term ProcClear(Term t,Term ind)
	{
	Term t1;
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	p_clear(t1);
	return 0;
	}
	
Term ProcStat(Term t, Term ind)
	{
    AtomStatistics();
	ListStatistics();
	return 0;
	}
	
Term ProcSetTex(Term t, Term ind)
	{
	List l;
	l=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	if(!is_list(l))
		{
		ErrorInfo(210);
		puts("illegal argument in SetTexName call.");
		return 0;
		}
	t=l;
	while(!is_empty_list(l))
		{
		Term t1;
		t1=ListFirst(l);
		if(!is_compound(t1) || CompoundName(t1)!=OPR_EQSIGN ||
			(!is_atom(CompoundArg1(t1))&& CompoundArg1(t1)!=NewInteger(0))|| 
			!is_atom(CompoundArg2(t1)))
				{
				ErrorInfo(210);
				puts("illegal argument in SetTexName call.");
				return 0;
				}
		if(CompoundArg1(t1)!=NewInteger(0))
			SetAtomProperty(CompoundArg1(t1),A_TEXNAME,CompoundArg2(t1));
		l=ListTail(l);
		}
	FreeAtomic(t);
	return 0;
	}
	
static List used_modules = 0;

Term ProcUse(Term t, Term i)
	{
	Atom af;
	List l;
	char pubuf[64];
	Term t1;
	
	if(!is_compound(t) || CompoundArity(t)!=1)
	{
		ErrorInfo(201);
		puts("wrong arg in 'use' statement");
		return 0;
	}
	
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	t1=l=CommaToList(t1);
	for(;t1;t1=ListTail(t1))
		{
			af=ListFirst(t1);
			if(!is_atom(af))
			{
				ErrorInfo(202);
				printf("use: bad file name '");
				WriteTerm(af);
				printf("'\n");
				continue;
			}
			if(ListMember(used_modules,af))
				continue;
			used_modules=AppendLast(used_modules,af);
			strcpy(pubuf,AtomValue(af));
			ReadFile(pubuf);
		}

	FreeAtomic(l);
	
	return 0;		

	}

extern int	opUndefAngleComb, opTriHeu, opBreakLines;
extern int C_F_WIDTH, L_P_WIDTH, P_D_WIDTH, MultByI, opOnlyMass;
extern int opReduceG5, opSetGpm, opRemDotWithFerm, opMaxiLegs;
extern int LagrHashSize, opSortLagr, opAutoWidths;
extern int opSplitCol1, opSplitCol2, opNeutrC8;
extern int opDoSymmetrize;
extern int infi_order, block_transf;
extern int WriteColors, write_all_vertices;
extern List opTexBOF, opTexEOF;
extern List opFAGS, opFAGE;
extern int opclsreal, oprcmatr, opclsred;
extern List abbr_coeffs, no_abbr_coeffs;

List UFOparah=0, UFOparth=0, UFOloreh=0, UFOverth=0, UFOcouph=0;

Term ProcOption(Term t, Term i)
{
	char *s;
	Term val;
	
	if(is_atom(t))
	{
		printf("  Defined options:\n");
		printf("UndefAngleComb   =   %d\n",opUndefAngleComb);
		printf("SmartAngleComb   =   %d\n",opTriHeu);
		printf("LagrHashSize     =   %d\n",LagrHashSize);
		printf("chepBreakLines   =   %d\n",opBreakLines);
		printf("SortLagr         =   %d\n",opSortLagr);
		printf("DoSymmetrize     =   %d\n",opDoSymmetrize);
		printf("OnlyMassTerms    =   %d\n",opOnlyMass);
		printf("ReduceGamma5     =   %d\n",opReduceG5);
		printf("SetGammaPM       =   %d\n",opSetGpm);
		printf("MultByI          =   %d\n",MultByI);
		printf("chepCFWidth      =   %d\n",C_F_WIDTH);
		printf("chepLPWidth      =   %d\n",L_P_WIDTH);
		printf("chepPDWidth      =   %d\n",P_D_WIDTH);
		printf("SplitCol1        =   %d\n",opSplitCol1);
		printf("SplitCol2        =   %d\n",opSplitCol2);
		
		puts("");
		return 0;
	}
	
	if(is_compound(t) && CompoundArity(t)==1)
	{
		Term t1;
		t1=ConsumeCompoundArg(t,1);
		FreeAtomic(t);
		t=CommaToList(t1);
		for(t1=t;t1;t1=ListTail(t1))
			ProcOption(ListFirst(t1),0);
		RemoveList(t);
		return 0;
	}
	
	if(!is_compound(t) || CompoundArity(t)!=2 || 
			CompoundName(t) != OPR_EQSIGN ||
			!is_atom(CompoundArg1(t)))
	{
		ErrorInfo(413);
		puts("wrong syntax in 'option' statement.");
		return 0;
	}
	
	 s=AtomValue(CompoundArg1(t));
	 val=ConsumeCompoundArg(t,2);
	 FreeAtomic(t);
	 
	 if(strcmp(s,"UndefAngleComb")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<0 ||
				 IntegerValue(val)>1)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'UndefAngleComb' option.");
			 return 0;
		 }
		 opUndefAngleComb=(int)IntegerValue(val);
		 return 0;
	 }
	 
	 if(strcmp(s,"SmartAngleComb")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<0 ||
				 IntegerValue(val)>4)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'SmartAngleComb' option.");
			 return 0;
		 }
		 opTriHeu=(int)IntegerValue(val);
		 return 0;
	 }
	 
	 if(strcmp(s,"LagrHashSize")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<1)
		 {
			 ErrorInfo(415);
			 puts("value for 'LagrHashSize' option must be >0.");
			 return 0;
		 }
		 LagrHashSize=(int)IntegerValue(val);
		 return 0;
	 }
	 
	 if(strcmp(s,"chepBreakLines")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<0 ||
				 IntegerValue(val)>1)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'chepBreakLines' option.");
			 return 0;
		 }
		 opBreakLines=(int)IntegerValue(val);
		 return 0;
	 }
	 
	 if(strcmp(s,"SortLagr")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<0 ||
				 IntegerValue(val)>1)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'SortLagr' option.");
			 return 0;
		 }
		 opSortLagr=(int)IntegerValue(val);
		 return 0;
	 }
	 
	 if(strcmp(s,"DoSymmetrize")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<0 ||
				 IntegerValue(val)>1)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'DoSymmetrize' option.");
			 return 0;
		 }
		 opDoSymmetrize=(int)IntegerValue(val);
		 return 0;
	 }
	 
	 if(strcmp(s,"AutoWidths")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<0 ||
				 IntegerValue(val)>1)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'AutoWidths' option.");
			 return 0;
		 }
		 opAutoWidths=(int)IntegerValue(val);
		 return 0;
	 }
	 
	 if(strcmp(s,"OnlyMassTerms")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<0 ||
				 IntegerValue(val)>1)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'OnlyMassTerms' option.");
			 return 0;
		 }
		 opOnlyMass=(int)IntegerValue(val);
		 return 0;
	 }

	 if(strcmp(s,"RemDotWithFerm")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<0 ||
				 IntegerValue(val)>1)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'RemDotWithFerm' option.");
			 return 0;
		 }
		 opRemDotWithFerm=(int)IntegerValue(val);
		 return 0;
	 }
	 
	 if(strcmp(s,"WriteAll")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<0 ||
				 IntegerValue(val)>1)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'WriteAll' option.");
			 return 0;
		 }
		 write_all_vertices=(int)IntegerValue(val);
		 return 0;
	 }
	 
	 if(strcmp(s,"MaxiLegs")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<0)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'MaxiLegs' option.");
			 return 0;
		 }
		 opMaxiLegs=(int)IntegerValue(val);
		 return 0;
	 }
	 
	 if(strcmp(s,"MultByI")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<0 ||
				 IntegerValue(val)>1)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'MultByI' option.");
			 return 0;
		 }
		 MultByI=(int)IntegerValue(val);
		 return 0;
	 }
	 
	 if(strcmp(s,"ReduceGamma5")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<0 ||
				 IntegerValue(val)>1)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'ReduceGamma5' option.");
			 return 0;
		 }
		 opReduceG5=(int)IntegerValue(val);
		 return 0;
	 }
	 
	 if(strcmp(s,"SetGammaPM")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<0 ||
				 IntegerValue(val)>1)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'SetGammaPM' option.");
			 return 0;
		 }
		 opSetGpm=(int)IntegerValue(val);
		 return 0;
	 }
	 
	if(strcmp(s,"chepCFWidth")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<=0)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'chepCFWidth' option.");
			 return 0;
		 }
		 C_F_WIDTH=(int)IntegerValue(val);
		 return 0;
	 }
	 
	 if(strcmp(s,"chepLPWidth")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<=0)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'chepLPWidth' option.");
			 return 0;
		 }
		 L_P_WIDTH=(int)IntegerValue(val);
		 return 0;
	 }
	 
	 if(strcmp(s,"chepPDWidth")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<0)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'chepPDWidth' option.");
			 return 0;
		 }
		 P_D_WIDTH=(int)IntegerValue(val);
		 return 0;
	 }
	 
	if(strcmp(s,"TexFileStart")==0)
	 {
		 if(!is_list(val))
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'TexFileStart' option.");
			 return 0;
		 }
		 opTexBOF=val;
		 return 0;
	 }
	 
	 if(strcmp(s,"TexFileEnd")==0)
	 {
		 if(!is_list(val))
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'TexFileEnd' option.");
			 return 0;
		 }
		 opTexEOF=val;
		 return 0;
	 }
	 
	if(strcmp(s,"FAGenStart")==0)
	 {
		 if(!is_list(val))
		 {
			 ErrorInfo(418);
			 puts("wrong value for 'FAGenStart' option.");
			 return 0;
		 }
		 opFAGS=val;
		 return 0;
	 }
	 
	 if(strcmp(s,"FAGenEnd")==0)
	 {
		 if(!is_list(val))
		 {
			 ErrorInfo(418);
			 puts("wrong value for 'GAGenEnd' option.");
			 return 0;
		 }
		 opFAGE=val;
		 return 0;
	 }
	 
	 if(strcmp(s,"SplitCol1")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<-1 || IntegerValue(val)>2)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'SplitCol1' option.");
			 return 0;
		 }
		 opSplitCol1=(int)IntegerValue(val);
		 return 0;
	 }
	 
	 if(strcmp(s,"SplitCol2")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<0 || IntegerValue(val)>1)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'SplitCol2' option.");
			 return 0;
		 }
		 opSplitCol2=(int)IntegerValue(val);
		 return 0;
	 }
	 
	 if(strcmp(s,"WriteColors")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<0 || IntegerValue(val)>1)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'WriteColors' option.");
			 return 0;
		 }
		 WriteColors=(int)IntegerValue(val);
		 return 0;
	 }
	 
	 if(strcmp(s,"InfiOrder")==0)
	   {
	     if(!is_integer(val))
	       {
			 ErrorInfo(415);
			 puts("wrong value for 'InfiOrder' option.");
			 return 0;
	       }
	     infi_order=(int)IntegerValue(val);
	     return 0;
	   }
	   
	 if(strcmp(s,"BlockTransf")==0)
	   {
	     if(!is_integer(val))
	       {
			 ErrorInfo(415);
			 puts("wrong value for 'BlockTransf' option.");
			 return 0;
	       }
	     block_transf=(int)IntegerValue(val);
	     return 0;
	   }

	 if(strcmp(s,"NeutralC8")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<0 || IntegerValue(val)>1)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'NeutralC8' option.");
			 return 0;
		 }
		 opNeutrC8=(int)IntegerValue(val);
		 return 0;
	 } 
	 if(strcmp(s,"clsRealMatr")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<0 || IntegerValue(val)>1)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'clsRealMatr' option.");
			 return 0;
		 }
		 opclsreal=(int)IntegerValue(val);
		 return 0;
	 } 
	 if(strcmp(s,"clsRedMatr")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<0 || IntegerValue(val)>1)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'clsRedMatr' option.");
			 return 0;
		 }
		 opclsred=(int)IntegerValue(val);
		 return 0;
	 } 
	 if(strcmp(s,"clsRCMatr")==0)
	 {
		 if(!is_integer(val) || IntegerValue(val)<0 || IntegerValue(val)>1)
		 {
			 ErrorInfo(415);
			 puts("wrong value for 'clsRCMatr' option.");
			 return 0;
		 }
		 oprcmatr=(int)IntegerValue(val);
		 return 0;
	 } 
	 if(strcmp(s,"AbbrCoeff")==0)
	 {
		 if(val==NewInteger(0))
		 {
			 abbr_coeffs=0;
			 return 0;
		 }
		 if(is_list(val))
		 {
			 abbr_coeffs=val;
			 return 0;
		 }
		 ErrorInfo(415);
		 puts("wrong value for 'AbbrCoeff' option.");
		 return 0;
	 }
	if(strcmp(s,"NoAbbrCoeff")==0)
	 {
		 if(val==NewInteger(0))
		 {
			 no_abbr_coeffs=0;
			 return 0;
		 }
		 if(is_list(val))
		 {
			 no_abbr_coeffs=val;
			 return 0;
		 }
		 ErrorInfo(415);
		 puts("wrong value for 'AbbrCoeff' option.");
		 return 0;
	 }  
	 if(strcmp(s,"UFOParaHdr")==0)
	 {
		 if(is_list(val))
		 {
			 List l;
			 for(l=val;l;l=ListTail(l))
				 if(!is_atom(ListFirst(l)))
					 break;
			 if(l==0)
			 {
			 UFOparah=val;
			 return 0;
		 	}
		 }
		 ErrorInfo(415);
		 puts("wrong value for 'UFOxxxHdr' option.");
		 return 0;
	 } 
	 if(strcmp(s,"UFOPartHdr")==0)
	 {
		 if(is_list(val))
		 {
			 List l;
			 for(l=val;l;l=ListTail(l))
				 if(!is_atom(ListFirst(l)))
					 break;
			 if(l==0)
			 {
			 UFOparth=val;
			 return 0;
		 	}
		 }
		 ErrorInfo(415);
		 puts("wrong value for 'UFOxxxHdr' option.");
		 return 0;
	 } 
	 if(strcmp(s,"UFOVertHdr")==0)
	 {
		 if(is_list(val))
		 {
			 List l;
			 for(l=val;l;l=ListTail(l))
				 if(!is_atom(ListFirst(l)))
					 break;
			 if(l==0)
			 {
			 UFOverth=val;
			 return 0;
		 	}
		 }
		 ErrorInfo(415);
		 puts("wrong value for 'UFOxxxHdr' option.");
		 return 0;
	 } 
	 if(strcmp(s,"UFOLoreHdr")==0)
	 {
		 if(is_list(val))
		 {
			 List l;
			 for(l=val;l;l=ListTail(l))
				 if(!is_atom(ListFirst(l)))
					 break;
			 if(l==0)
			 {
			 UFOloreh=val;
			 return 0;
		 	}
		 }
		 ErrorInfo(415);
		 puts("wrong value for 'UFOxxxHdr' option.");
		 return 0;
	 } 
	 if(strcmp(s,"UFOCoupHdr")==0)
	 {
		 if(is_list(val))
		 {
			 List l;
			 for(l=val;l;l=ListTail(l))
				 if(!is_atom(ListFirst(l)))
					 break;
			 if(l==0)
			 {
			 UFOcouph=val;
			 return 0;
		 	}
		 }
		 ErrorInfo(415);
		 puts("wrong value for 'UFOxxxHdr' option.");
		 return 0;
	 } 

	ErrorInfo(413);
	printf("unknown option '%s'.\n",s);	 
				
	return 0;
}
	

	
