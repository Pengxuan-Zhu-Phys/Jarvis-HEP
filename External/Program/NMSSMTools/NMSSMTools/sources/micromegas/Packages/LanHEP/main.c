#include <stdlib.h>
#include <time.h>
#include <signal.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/stat.h>
#include "lanhep.h"

#ifdef MTHREAD
#include <pthread.h>

pthread_key_t TermsKey;
#endif

char *ModelName=NULL;
int  ModelNumber=1;

FILE *log_file = NULL;

char *OutputDirectory = NULL;
char *InputDirectory  = /*"/home/semenov/lanhep/mdl/"*/ NULL;
char *InputFile = NULL;
int  VerbMode = 0;
int VerbT=0;
time_t tm_lt=0, tm_12=0, tm_sl=0, tm_sw=0, tm_a2=0, tm_add=0, tm_rl=0, tm_wr=0;
time_t tm_rsc=0;
char *InitFile = "lhep.rc";
int  TexOutput = 0, FAver=8;
int  ForsedRedCol = 0;
int	 NoColors = 0;
int WriteColors=0;
int  No4Color = 0;
int  EvalPrm  = 0, EvalVrt=0, FindVal=0;
int  ChepVersion=4;
int  opAbbrVrt=0, opAbbArr=0;
int remove_rc=0;
int opRemDotWithFerm=1;
int MicroOmega=0;
int th_no=0;

#ifdef CompHEP
int CompOutput=1;
#else
int CompOutput=0;
#endif

#ifdef CalcHEP
int CalcOutput=1;
#else
int CalcOutput=0;
#endif

#ifdef FeynArts
int FAOutput=1;
#else
int FAOutput=0;
#endif

#ifdef UFO
int UFOutput=1;
#else
int UFOutput=0;
#endif



int end_with_tty=0;

extern int opSplitCol1, opSplitCol2, opEvalVrt, opNo4Scal;
extern int verb_herm, err_cnt;
extern int write_all_vertices, opMaxiLegs;
extern void SaveRules(char *);
extern void SetKeyFromArg(char *); 
extern int off_srefine, ch_sign, longest_lpline, longest_cfline, longest_pdline;
extern int eval_vrt_len, eval_vrt_more, kill_gamma_pm, opNoDummies,opTriHeu;
extern int SecondVaFu;
extern int inf_removed[10];
extern List ufo_params_nolha;

extern int C_F_WIDTH, L_P_WIDTH, P_D_WIDTH, TEX_max_pno, doinitfile;

static char nmbuf[1280], nmbuf0[1280];

extern void abbr_stat(void);

#ifdef __DJGPP__
static char *find_path(char *nm, char *env[])
	{
	char *pp;
	strcpy(nmbuf,nm);
    pp=strrchr(nmbuf,'/');
	if(pp==NULL)
		{
		puts("Fail locating lhep.rc file");
		exit(0);
		}
	pp[1]=0;
    strcat(nmbuf,"mdl/");
	return nmbuf;
	}
#else
static char *find_path(char *nm, char *env[])
	{
	int i;
	char *pp;
	strcpy(nmbuf0,nm);
bgn:
	for(i=0;i<1280;nmbuf[i++]=0);
	if(readlink(nmbuf0,nmbuf,128)!=-1)
		{
		strcpy(nmbuf0,nmbuf);
		goto bgn;
		}
	if(access(nmbuf0,0)==0)
		{
		pp=strrchr(nmbuf0,'/');
		if(pp==NULL)
			return "./mdl/";
		else
			{
			pp[1]=0;
			strcat(nmbuf0,"mdl/");
			return nmbuf0;
			}
		}
	i=0;
	while(env[i] && strncmp(env[i],"PATH=",5))
		i++;
	if(env[i]==NULL)
		{
		puts("No PATH environment variable found.");
		exit(0);
		}
	pp=env[i];
	while(pp[0])
		{
		int ll;
		ll=0;
		while(pp[0] && pp[0]!=':')
			{
			nmbuf[ll]=pp[0];
			pp++;
			ll++;
			}
		if(pp[0]) pp++;
		nmbuf[ll]='/';
		ll++;
		strcpy(nmbuf+ll,nmbuf0);
		if(access(nmbuf,0)==0)
			{
			strcpy(nmbuf0,nmbuf);
			goto bgn;
			}
		}

	puts("Fail locating lhep.rc file");
	exit(0);
	return "???";
	}
#endif

static void wrt_time(void)
{
	double d=1.0/CLOCKS_PER_SEC;
	printf("Times: lterm %.2f (to1 %.2f (sl %.2f sw %.2f) to2 %.2f ha %.2f\n",
		(double)tm_lt*d, (double)tm_12*d, (double)tm_sl*d, (double)tm_sw*d, (double)tm_a2*d, (double)tm_add*d);
	/*printf("       rsc= %.2f\n",tm_rsc*d);*/
	printf("       RedL %.2f Wrt %.2f\n",(double)tm_rl*d, (double)tm_wr*d);
}


extern void proc_hash(Term t);

static List l_ali=0, l_pro=0;

void ProcessTerm(Term t)
	{
	Term res;
	Atom procc=0;
	char regbuf[800];

	if(is_compound(t)&&CompoundName(t)==OPR_LOCAL)
			{
			List l1,la=ConsumeCompoundArg(t,1);
			if(l_pro==0)
				{
				ErrorInfo(0);
				puts("local: must be within template.");
				return;
				}
			la=CommaToList(la);
			for(l1=la;l1;l1=ListTail(l1))
				l_ali=AppendLast(l_ali,
					InterfSetAlias(MakeCompound1(OPR_ALIAS,ListFirst(l1)),0));
			RemoveList(la);
			return;
			}

/*	if(is_compound(t) && GetAtomProperty(CompoundName(t),OPR_ALIAS))
	{*/
	
	if(!is_compound(t) || (CompoundName(t)!=OPR_ALIAS  && CompoundName(t)!=OPR_UNALIAS &&
			CompoundName(t)!=OPR_IN && CompoundName(t)!=OPR_WHERE  ) )
		t=ProcessAlias(t);

	
	if(!is_compound(t) || (CompoundName(t)!=OPR_ALIAS  && 
			CompoundName(t)!=OPR_IN && CompoundName(t)!=OPR_WHERE  ) )
			proc_hash(t);

			/*WriteTerm(t);puts("");*/
	/*}*/
	
	if(is_list(t))
	{
		List l, alisav;
		l_pro++;
		alisav=l_ali;
		l_ali=0;		
		for(l=t;l;l=ListTail(l))
			ProcessTerm(ListFirst(l));
		if(l_ali)
			RemoveAlias(l_ali);
		RemoveList(t);
		l_pro--;
		l_ali=alisav;
		return;
	}


    if(is_atom(t))
		procc=t;
	if(is_compound(t))
		{
		procc=CompoundName(t);
		}
    if(procc==0)
		{
        if(!IsTermInput())
			printf("File \"%s\", line %d: ",CurrentInputFile(),
				CurrentInputLine());
		printf("Not a valid operator \'");
		WriteTerm(t);
		printf("\'.\n");
		FreeAtomic(t);
		return;
		}

	if(!IsTermInput())
		sprintf(regbuf, "ProcTerm: File %s, line %d: %s statement.",
			CurrentInputFile(), CurrentInputLine(), AtomValue(procc));
	else
		sprintf(regbuf, "(tty): %s statement.",AtomValue(procc));

	RegisterLine(regbuf);

	if(!is_function(t,NULL))
		{
        if(!IsTermInput())
			printf("File \"%s\", line %d: ",CurrentInputFile(),
				CurrentInputLine());
		printf("Unknown operator \'%s\' in expression\n\t",AtomValue(procc));
		WriteTerm(t);
		puts("");
		FreeAtomic(t);
		UnregisterLine();
		return;
		}

    if(is_compound(t) && ((VerbMode && CompoundName(t)==OPR_LTERM) ||
			VerbMode==3))
    	{
    	char sment[100024];
    	sWriteTerm(sment,t);
    	sment[55]='.';
    	sment[56]='.';
    	sment[57]='.';
    	sment[58]=0;
    	printf("%s:%3d: %s\n",CurrentInputFile(),CurrentInputLine(),sment);
    	}
    	
	res=CallFunction(t,0);
	if(res!=0 && (IsTermInput()))
		{
		if(is_compound(res) && CompoundArity(t)==2 && 
				CompoundName(res)==A_ALG1)
			res=Alg1ToExpr(res);
		if(is_list(res))
			DumpList(res);
		else
			{
			WriteTerm(res);
			puts("");
			}
		FreeAtomic(res);
		}
	UnregisterLine();
	}


static void stop(int sig)
	{
	puts("\nInterrupt signal received, exiting...");
	DumpRegistered();
	AtomStat1();
	ListStat1();
	if(VerbT) wrt_time();
	exit(0);
	}

static void sigdmp(int sig)
	{
	puts("\nCurrent state:");
	DumpRegistered();
	AtomStat1();
	ListStat1();
	puts("");
	fflush(stdout);
	signal(SIGUSR1, sigdmp);
	}


static void sigsegv(int sig)
	{
	puts("Internal error (segmentation violation), aborted...");
	DumpRegistered();
	AtomStatistics();
	ListStatistics();
	abort();
	}


static void setoutput(int t)
{
	CompOutput=0;
	CalcOutput=0;
	FAOutput=0;
	UFOutput=0;
	eval_vrt_len=0;
	eval_vrt_more=0;
	switch(t)
	{
		case 1:
			CompOutput=1;
			ChepVersion=4;
			opSplitCol1=2;
			opSplitCol2=1;
			break;
		case 2:
			CalcOutput=1;
			eval_vrt_len=2;
			opNoDummies=1;
			eval_vrt_more=1;
			opSplitCol1=2;
			opSplitCol2=1;
			opRemDotWithFerm=0;
			break;
		case 3:
			FAOutput=1;
			FAver=8;
			opNoDummies=1;
			opSplitCol1=0;
			opSplitCol2=0;
			break;
		case 4:
			TexOutput=1;
			break;
		case 5:
			UFOutput=1;
			opNoDummies=1;
			opSplitCol1=0;
			opSplitCol2=0;
			break;
		default:
			puts("internal error (mainso)");
	}
}
	
extern Atom llparam;
static void exv(int, char **);


int main(int argc, char **argv, char **env)
	{
	int i;
	int save_flag=0;
	char *save_file=NULL;
#ifdef MTHREAD
	int mtermval=0;
	pthread_key_create(&TermsKey,NULL);
	pthread_setspecific(TermsKey,&mtermval);
#endif
	
	
	if(sizeof(Term)!=8)
		{
		puts("Compilation error");
		return -1;
		}

	/*{
		Atom a=NewAtom("test",0), b=NewAtom("hoho",0), c=NewAtom("oooh",0);
		printf("%lx\n",a);
		Term t=MakeCompound1(a,NewInteger(1));
		Term t2=MakeCompound2(b,CopyTerm(t),c);
//		printf("%lx %lx\n",a,t);
		WriteTerm(t2);puts("");
		List l=MakeList2(t,t2);
		DumpList(l);
		l=AppendLast(l,a);
		l=AppendLast(l,b);
		DumpList(l);
		SetAtomProperty(a,b,l);
		
		while(!is_empty_list(l))
		{
			printf("%lx %lx\n",l,ListTail(l));
			l=ListTail(l);
		}
	//	return 0;
	}*/
		
	if(argc>1 && strcmp(argv[1],"-exv")==0)
		{
		exv(argc-1,argv+1);
		return 0;
		}

	InitAtoms();
	InitFuncs();
	AlwaysBracets = 0;
	WideWriting = 0;
	signal(SIGINT, stop);
	signal(SIGSEGV, sigsegv);
	signal(SIGUSR1, sigdmp);
	
	if(CompOutput)
		setoutput(1);
	if(CalcOutput)
		setoutput(2);
	if(FAOutput)
		setoutput(3);
	if(UFOutput)
		setoutput(5);
		

	for(i=1; i<argc; i++)
		{

		if(strcmp(argv[i],"-mth")==0)
			{
			if(i<argc-1 && isdigit(argv[i+1][0]))
				{
				sscanf(argv[++i],"%d",&th_no);
				if(th_no<0 || th_no>16)
					th_no=2;
				}
			else
				th_no=2;
#ifndef MTHREAD
			puts("Error: multithread processing is switched off.");
			th_no=0;
#endif
			continue;
			}

		if(strcmp(argv[i],"-v")==0)
			{
			VerbMode=1;
			continue;
			}
		if(strcmp(argv[i],"-vv")==0)
			{
			VerbMode=2;
			continue;
			}
		if(strcmp(argv[i],"-vt")==0)
			{
			VerbT=1;
			continue;
			}
		if(strcmp(argv[i],"-vvv")==0)
			{
			VerbMode=3;
			continue;
			}
		if(strcmp(argv[i],"-findval")==0)
			{
			FindVal=1;
			continue;
			}
		if(strcmp(argv[i],"-c4")==0||
			strcmp(argv[i],"-co")==0||
			strcmp(argv[i],"-CompHEP")==0)
			{
			setoutput(1);
			continue;
			}	
		if(strcmp(argv[i],"-c3")==0)
			{
			setoutput(1);
			ChepVersion=3;
			continue;
			}	
		if(strcmp(argv[i],"-calchep")==0||
		   strcmp(argv[i],"-CalcHEP")==0||
		   strcmp(argv[i],"-ca")==0)
			{
			setoutput(2);
			continue;
			}
		if(strcmp(argv[i],"-mOmega")==0)
			{
			setoutput(2);
			ChepVersion=4;
			MicroOmega=1;
			continue;
			}	
		if(strcmp(argv[i],"-OutDir")==0)
			{
			OutputDirectory=argv[++i];
			if(access(OutputDirectory,F_OK))
            {
                if(mkdir(OutputDirectory,0777))
                {
                    perror(OutputDirectory);
                    OutputDirectory=0;
                }
                else printf("Folder %s is created.\n",OutputDirectory);
            }
			continue;
			}
		if(strcmp(argv[i],"-InDir")==0)
			{
			InputDirectory=argv[++i];
			continue;
			}
		if(strcmp(argv[i],"-allvrt")==0)
			{
			write_all_vertices=1;
			continue;
			}
		if(strcmp(argv[i],"-rc")==0)
			{
			InitFile=argv[++i];
			continue;
			}
		if(strcmp(argv[i],"-tex")==0)
			{
			setoutput(4);
			opSplitCol1=0;
			continue;
			}
		if(strcmp(argv[i],"-fa4")==0)
			{
			setoutput(3);
			FAver=4;
			continue;
			}
		if(strcmp(argv[i],"-fa6")==0)
			{
			setoutput(3);
			FAver=6;
			continue;
			}
		if(strcmp(argv[i],"-feynarts")==0 || strcmp(argv[i],"-FeynArts")==0 ||
				strcmp(argv[i],"-fa8")==0||strcmp(argv[i],"-fa")==0)
			{
			setoutput(3);
			continue;
			}
		if(strcmp(argv[i],"-ufo")==0 || strcmp(argv[i],"-UF")==0)
			{
			setoutput(5);
			continue;
			}
		if(strcmp(argv[i],"-save")==0)
			{
			save_flag=1;
			save_file=argv[++i];
			setoutput(1);
			ChepVersion=3;
			kill_gamma_pm=1;
			printf("Feynman rules will be saved in '%s' file.\n",save_file);
			continue;
			}
		if(strcmp(argv[i],"-eval-prm")==0)
			{
			EvalPrm=1;
			continue;
			}
		if(strcmp(argv[i],"-eval-vrt")==0)
			{
			EvalVrt=1;
			continue;
			}
		if(strcmp(argv[i],"-frc")==0)
			{
			ForsedRedCol=1;
			continue;
			}
		if(strcmp(argv[i],"-nocolor")==0)
			{
			NoColors=1;
			continue;
			}
		if(strcmp(argv[i],"-colors")==0)
			{
			WriteColors=1;
			continue;
			}
        if(strcmp(argv[i],"-no4color")==0)
			{
            No4Color=1;
			continue;
            }
		if(strcmp(argv[i],"-no4scal")==0)
			{
            opNo4Scal=1;
			continue;
            }
		if(strcmp(argv[i],"-nocdot")==0)
			{
			TEX_set_dot=0;
			continue;
			}
		if(strcmp(argv[i],"-v-charges")==0)
			{
			verb_charge=1;
			continue;
			}
		if(strcmp(argv[i],"-v-herm")==0)
			{
			verb_herm=1;
			continue;
			}
		if(strcmp(argv[i],"-v-imprt")==0)
			{
			verb_imprt=1;
			continue;
			}
		if(strcmp(argv[i],"-off-srefine")==0)
			{
			off_srefine=1;
			continue;
			}
		if(strcmp(argv[i],"-chep-srefine")==0)
			{
			ch_sign=1;
			continue;
			}
		if(strcmp(argv[i],"-texLines")==0)
			{
			sscanf(argv[++i],"%d",&TEX_lines);
			continue;
			}
		if(strcmp(argv[i],"-texLineLength")==0)
			{
			sscanf(argv[++i],"%d",&TEX_spec_in_line);
			continue;
			}
		if(strcmp(argv[i],"-texMaxPrtNo")==0)
			{
			sscanf(argv[++i],"%d",&TEX_max_pno);
			continue;
			}
		if(strncmp(argv[i],"-abbr",5)==0)
			{
			if(eval_vrt_len)
				{
				puts("Error: -evl and -abbr options are not compatible");
				continue;
				}
			opAbbrVrt=1;
			opAbbArr=1;
			opEvalVrt=0;
			opTriHeu=0;
			if(argv[i][5]>'1' && argv[i][5]<'9')
				opAbbrVrt=argv[i][5]-'1'+1;
			if(argv[i][5]=='A')
				opAbbArr=0;
			opNoDummies=1;
			continue;
			}
		if(strcmp(argv[i],"-evl")==0)
			{
			if(opAbbrVrt)
			{
				puts("Error: -evl and -abbr options are not compatible");
				i++;
				continue;
			}
			opNoDummies=0;
			eval_vrt_more=0;
			sscanf(argv[++i],"%d",&eval_vrt_len);
			if(eval_vrt_len==2)
			{
				opNoDummies=1;
				eval_vrt_more=1;
				/*kill_gamma_pm=1;*/
			}
			continue;
			}
			if(strcmp(argv[i],"-noabbr")==0)
			{
			eval_vrt_len=0;
			continue;
			}

		if(strcmp(argv[i],"-sleep")==0)
			{
			unsigned int sec;
			sscanf(argv[++i],"%d",&sec);
			sleep(sec*60);
			continue;
			}
		if(strcmp(argv[i],"-key")==0)
			{
			SetKeyFromArg(argv[++i]);
			continue;
			}
		if(strcmp(argv[i],"-edbg")==0)
			{
			end_with_tty=1;
			continue;
			}
		if(strcmp(argv[i],"-norc")==0)
			{
			remove_rc=1;
			continue;
			}
		if(argv[i][0]=='-')
			{
			printf("Error: unknown option %s.\n",argv[i]);
			continue;
			}
		if(InputFile)
			{
			ErrorInfo(0);
			printf(": unknown option %s.\n",argv[i]);
			continue;
			}
		InputFile=argv[i];
		}

	if(!write_all_vertices && !TexOutput)
	{
		if(UFOutput)
			opMaxiLegs=6;
		else
			opMaxiLegs=4;
	}

	if(CompOutput==0 && CalcOutput==0 && TexOutput==0 &&
			FAOutput==0 && UFOutput==0)
	{
		puts("Error: output format is not set. Use command line options:");
		puts("\t-CompHEP or -co for CompHEP format;");
		puts("\t-CalcHEP or -ca for CalcHEP format;");
		puts("\t-FeynArts or -fa for FeynArts format;");
		puts("\t-ufo for UFO format;");
		puts("\t-tex for LaTeX format.");
		return 0;
	}
	
	if(UFOutput)
	{
		FAOutput=1;
		opAbbrVrt=1;
		opAbbArr=0;
		opEvalVrt=0;
		opTriHeu=0;
		opNoDummies=1;
	}
		
	if(MicroOmega)
		SetKeyFromArg("MicrOmega=1");
	else
		SetKeyFromArg("MicrOmega=0");
	
	if(FAOutput)
		SetKeyFromArg("FeynArts=1");
	else
		SetKeyFromArg("FeynArts=0");

	if(UFOutput)
		SetKeyFromArg("UFO=1");
	else
		SetKeyFromArg("UFO=0");

	if(CalcOutput)
		SetKeyFromArg("CalcHEP=1");
	else
		SetKeyFromArg("CalcHEP=0");
		
	if(InputDirectory==NULL)
		{
		InputDirectory=find_path(argv[0],env);
		if(VerbMode)
			printf("Input directory is '%s'\n",InputDirectory);
		}
	doinitfile=1;	
   ReadFile(InitFile);
	doinitfile=0;
	if(InputFile!=NULL)
		ReadFile(InputFile);
	else
		{
		printf("Welcome to LanHEP                                Version 4.0.0  (Jan 16 2020)\n");
		/*
		log_file=fopen("lhep.log","w");
		if(log_file==NULL)
			printf("Warning: can not open file lhep.log for writing.\n");
		*/
		ReadFile(NULL);
		if(log_file!=NULL)
			{
			fclose(log_file);
			log_file=0;
			}
		}

	ProcReduceLagr(A_I,0);

puts("");

/*	AtomStatistics();
	ListStatistics();
*/
	if(save_flag)
		{
		Term t;
		SaveRules(save_file);
		if(itrSetIn("q.sav")==0)
			return 0;

		while((t=itrIn())!=0)
			{
			WriteTerm(t);
			puts("");
			}

		itrCloseIn();

		return 0;
		}

	if(ModelNumber!=0 || InputFile!=NULL)
		{
		time_t stt=clock();
		
		if(MicroOmega)
			ModelNumber=1;
		RegisterLine("MAIN: writing lagrangian.");
		WriteLagrFile(ModelNumber,ModelName);
		UnregisterLine();
		RegisterLine("MAIN: writing parameters and particles.");
		WriteParameters(ModelNumber,ModelName);
		WriteParticles(ModelNumber,ModelName);
		if(MicroOmega)
		{
			ModelNumber=2;
			SecondVaFu=1;
			WriteParameters(ModelNumber,ModelName);
		}
		UnregisterLine();
		WriteExtlib(ModelNumber,ModelName);
		WriteCpart(ModelNumber,ModelName);
		
		tm_wr=(clock()-stt);
		
		if(!TexOutput && !FAOutput && C_F_WIDTH<longest_cfline)
		{
			printf("Error: 'Common Factor' field longer than maximum of %d symbols.\n",
					C_F_WIDTH);
			printf("Use the statement 'option chepCFWidth=%d.'\n",longest_cfline+1);
		}
		if(!TexOutput && !FAOutput && L_P_WIDTH<longest_lpline)
		{
			printf("Error: 'Lorentz part' field longer than maximum of %d symbols.\n",
					L_P_WIDTH);
			printf("Use the statement 'option chepLPWidth=%d.'\n",longest_lpline+1);
		}
		if(!TexOutput && !FAOutput && P_D_WIDTH<longest_pdline)
		{
			if(P_D_WIDTH)
			{
			printf("Error: 'Parameter expression' field for '%s' is longer than maximum of %d symbols.\n",
					AtomValue(llparam),P_D_WIDTH);
			printf("Use the statement 'option chepPDWidth=%d.'\n",longest_pdline+1);
			}
			else if(longest_pdline>100)
				printf("Longest parameter expression is %d symbols for '%s'.\n",longest_pdline,
					AtomValue(llparam));
		}
		
		if(!TexOutput && !FAOutput && VerbMode)
			{
			printf("Longest Common factor line is %d symbols (max %d)\n",
											longest_cfline,C_F_WIDTH);
			printf("Longest Lorentz part  line is %d symbols (max %d)\n",
											longest_lpline,L_P_WIDTH);
			printf("Longest Parameter dep line is %d symbols (max %d)\n\n",
											longest_pdline,P_D_WIDTH);
			}
		}
/*
	else
		puts("Model # is zero");
*/


	if(UFOutput && ufo_params_nolha)
	{
	  printf("Warning: LHA block/index info is not provided for external parameters: ");
	  WriteTerm(ufo_params_nolha);
	  printf(".\nUse lha statements to set LHA data:\n\tlha parameter_name=block_name(index).\n\n");
	}
	if(VerbMode)
		{
		int i;
		printf("%5.1fMB of memory used.\n",
			(double)(ListMemory()+TermMemory())*0.001);
		abbr_stat();
		AtomStatistics();
		ListStatistics();
	/*
		for(i=0;i<10;i++)
			if(inf_removed[i]) printf("%d pwr of inf: %d\n",i,inf_removed[i]);*/
		}

	if(err_cnt)
		puts("!!!!!!!!!!!!! THERE WERE ERRORS DURING PROCESSING !!!!!!!!");
		
		if(end_with_tty)
		{
			ReadFile(NULL);
			puts("");
		}
	if(VerbT) wrt_time();	
	return 0;
	}



Term ProcessQuit(Term t, Term ind)
	{
	if(is_compound(t))
	{
	    AtomStatistics();
		ListStatistics();
	}
	if(log_file!=NULL)
		fclose(log_file);
	/*
	if(ModelNumber!=0)
		{
		WriteParameters(ModelNumber,ModelName);
		WriteLagrFile(ModelNumber,ModelName);
		WriteParticles(ModelNumber,ModelName);
		}
	*/
	exit(0);
	return 0;
	}


static char cbuf[5000];

static void exv(int ac, char *av[])
	{
	FILE *f;
	int pno;
	pno=ac-2;
	if(pno<1) return;
	sprintf(cbuf,"lgrng%s.mdl",av[1]);
	f=fopen(cbuf,"r");
	if(f==NULL)
		{
		printf("can not open file '%s'\n",cbuf);
		return;
		}
	
	fgets(cbuf,5000,f);
	fgets(cbuf,5000,f);
	fgets(cbuf,5000,f);
	
	while(fgets(cbuf,5000,f)==cbuf)
		{
		int ppos=0;
		int sflag=1;
		int i,j;
		for(i=0;i<4;i++)
			{
			int clen=0;
			char sav;
			int ffo=0;
			if(cbuf[ppos]==' ')
				break;
			while(cbuf[ppos+clen]!=' ' && cbuf[ppos+clen]!='|')
				clen++;
			sav=cbuf[ppos+clen];
			cbuf[ppos+clen]=0;
			for(j=0;j<pno;j++)
				if(strcmp(cbuf+ppos,av[j+2])==0)
					ffo=1;
			if(ffo==0)
				sflag=0;
			cbuf[ppos+clen]=sav;
			while(cbuf[ppos]!='|')
				ppos++;
			ppos++;
			}
		if(sflag)
			printf("%s",cbuf);
		}
	fclose(f);
	}

