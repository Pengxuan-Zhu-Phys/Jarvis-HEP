#include "lanhep.h"
#include <unistd.h>
#include <time.h>
#include <string.h>
#ifdef MTHREAD
#include <pthread.h>
#endif

int C_F_WIDTH = 50;
int L_P_WIDTH = 600;

extern int max_prt_lenl;

static int LagrPrtclsNo=4;
static int lagr_reduced=1;

void alg1_cncl(Term m1);
void alg2_abbr_vrt(Term a2);
void alg2_setcls(List a2l);

extern int VerbMode, opSetGpm, opAbbrVrt, UFOutput;
int  opEvalVrt = 1;

int write_all_vertices = 0;
int MultByI=0;
unsigned int LagrHashSize=1713;
int opSortLagr=1;

List *lagr_hash = NULL;

List sel_vrt(List, int, int);

int longest_lpline=0, longest_cfline=0;
int opBreakLines=0;

List alg2_denorm(Term);
void alg2_norm(Term);

void alg2_hash_add(List *, int , List);

List all_vert_list(void)
	{
		List l;
		int i;
		
		if(lagr_hash==NULL)
			return 0;

		ProcReduceLagr(A_I,0);

		l=NewList();
		for(i=0;i<LagrHashSize;i++)
			l=ConcatList(l,CopyTerm(lagr_hash[i]));
		
		return l;
	}

List all_vert_list2(void)
	{
		List l,l1;
		int i;
		
		if(lagr_hash==NULL)
			return 0;

		l=NewList();
		for(i=0;i<LagrHashSize;i++)
			for(l1=lagr_hash[i];l1;l1=ListTail(l1))
			{
			Term a2=ListFirst(l1);
			if(ListLength(CompoundArg1(a2))<3)
			l=AppendLast(l,CopyTerm(a2));
			}
		
		return l;
	}

static int need_herm_conj=0;

static int sort_lagr(Term, Term), sort_lagr_fa(Term, Term);

Term ProcAddHermConj(Term t, Term ind)
	{
	need_herm_conj=1;
	return NewInteger(0);
	}
	

static void plt_2(Term t, List nind);
	
int lterm_used=0;
extern time_t tm_lt, tm_add, tm_a2;
	
Term ProcLTerm(Term t, Term ind)
	{
	Term t1;
	char regbuf[80];
	time_t stt=clock();
	List nind=0;
	int e1_a, e1_i=0;
	
	lagr_reduced=0;
	lterm_used=1;
		
	if(lagr_hash==NULL)
	{
		int i;
		lagr_hash=malloc(sizeof(List)*LagrHashSize);
		if(lagr_hash==NULL)
		{
			puts("Error: can not allocate memory for Lagrangian"); 
			exit(0);
		}
		for(i=0;i<LagrHashSize;i++)
			lagr_hash[i]=NewList();
	}
	
	need_herm_conj=0;
	
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	if(is_atom(t1) && AtomValue(t1)[0]=='?')
		{
		WriteLagr(0,0);
		return 0;
		}

		
	RegisterLine("ProcLTerm: ExprTo11");
	
	t1 = ExprTo11(t1, &nind);
	
	UnregisterLine();
	if(t1==0)
		return 0;
	
	e1_a=ListLength(t1);
	
	for(t=t1;t;t=ListTail(t))
	{
		sprintf(regbuf,"ProcLTerm: %d monomials of %d\n",++e1_i,e1_a);
		alg1_cncl(ListFirst(t));
		RegisterLine(regbuf);
		
		if(VerbMode>1)
		{
			Term t3,t4;
			t3=MakeList1(CopyTerm(ListFirst(t)));
			t3=MakeCompound2(A_ALG1,t3,0);
			alg1_fix_delta(t3);
			t4=Alg1ToExpr(t3);
			WriteTerm(t4);
			puts("");
			FreeAtomic(t4);
		
/*			List l2,l3;
			int n,d;
			l2=CompoundArgN(ListFirst(t),3);
			l3=CompoundArgN(ListFirst(t),4);
			n=IntegerValue(CompoundArg1(ListFirst(t)));
			d=IntegerValue(CompoundArg2(ListFirst(t)));
			printf("\t\t%d/%d",n,d>0?d:-d);
			if(d<0)
				printf("*i");
			for(;l2;l2=ListTail(l2))
			{
				Term q;
				q=CompoundArg2(ListFirst(l2));
				if(is_atom(q))
					printf("*%s",AtomValue(q));
				else
					printf("*{...}");
			}
			for(;l3;l3=ListTail(l3))
			{
				Term q;
				q=CompoundArg2(ListFirst(l3));
				printf("/%s",AtomValue(q));
			}
			puts("");
*/
		}
		plt_2(MakeList1(ListFirst(t)),nind);
		UnregisterLine();
	}
	
	RemoveList(t1);
	tm_lt+=(clock()-stt);
	return 0;
}
	
static void plt_2(Term t, List nind)
{
	Term t1;
	char regbuf[80];
	time_t stt;
	
	t1=ExprTo12(t,nind);
	
	if(t1==0)
		return;
	
	
	sprintf(regbuf,"ProcLTerm: 1to2, %d mterms",
		(int)ListLength(CompoundArg1(t1)));
	RegisterLine(regbuf);
	
	stt=clock();
	t1=Alg1to2(t1);
	tm_a2+=(clock()-stt);
	
	UnregisterLine();

	if(need_herm_conj)
		{
		sprintf(regbuf,"ProcLTerm: alg2_add_herm_conj, %d terms",
			(int)ListLength(t1));
		RegisterLine(regbuf);
		t1=alg2_add_hermconj(t1);
		UnregisterLine();
		}

	sprintf(regbuf,"ProcLTerm: alg2_add, %d terms",
		(int)ListLength(t1));
	RegisterLine(regbuf);
	


/*	lgrng=alg2_add(lgrng,t1);*/
	
	stt=clock();
	alg2_hash_add(lagr_hash,(int)LagrHashSize,t1);
	tm_add+=(clock()-stt);
	
	UnregisterLine();

	return;
	}



static char wpbuf[128];


void texWriteLagr(int fno)
	{
	FILE *f;
	List l;

	if(OutputDirectory!=NULL)
		sprintf(wpbuf,"%s/lgrng%d.tex",OutputDirectory,fno);
	else
		sprintf(wpbuf,"lgrng%d.tex",fno);
	f=fopen(wpbuf,"w");
	if(f==NULL)
		{
		printf("Can not open file \'%s\' for writing.\n",wpbuf);
		perror("");
		return;
		}


	l=all_vert_list();
	if(opSortLagr)
		l=SortedList(l,sort_lagr);
	check_em_charge(l);
	tex_write_lagr(l,f);


	fclose(f);
	return;
	}

#include <dirent.h>
extern char *InputDirectory;

void FAWriteLagr(int fno)
	{
	FILE *f;
	List l;
	if(OutputDirectory)
		printf("Writing files to the folder %s\n",OutputDirectory); 
	if(UFOutput==0)
	{
	if(OutputDirectory!=NULL)
		sprintf(wpbuf,"%s/model%d.mod",OutputDirectory,fno);
	else
		sprintf(wpbuf,"model%d.mod",fno);
	f=fopen(wpbuf,"w");
	if(f==NULL)
		{
		printf("Can not open file \'%s\' for writing.\n",wpbuf);
		perror("");
		return;
		}
	else
		printf("Writing %s...\n",wpbuf);
	}
	else
		printf("Writing vertices.py...\n");
		
	l=all_vert_list();
                
	if(opSortLagr)
		l=SortedList(l,sort_lagr_fa);
	check_em_charge(l);
	if(UFOutput==0)
	{
		FA_write_lagr(l,f);
		fclose(f);
	}
	else
		UF_write_lagr(l);

	if(UFOutput==0)
		printf("Writing model%d.h, mdl_ini%d.F...\n",fno,fno);
	else
		printf("Writing parameters.py, couplings.py ...\n");
	FAWriteParameters(fno);
	
	if(UFOutput==0)
	{
	if(OutputDirectory!=NULL)
		sprintf(wpbuf,"%s/model%d.gen",OutputDirectory,fno);
	else
		sprintf(wpbuf,"model%d.gen",fno);
		f=fopen(wpbuf,"w");
	if(f==NULL)
		{
		printf("Can not open file \'%s\' for writing.\n",wpbuf);
		perror("");
		return;
		}
	printf("Writing %s...\n",wpbuf);
	FA_write_gen(f);
	fclose(f);
	}
	else
	{
		printf("Writing lorentz.py...\n");
		UF_write_gen();
	}
	if(UFOutput)
	{
	printf("\nCopying static files...\n");
	char cbuf[1000];
	sprintf(cbuf,"%s/ufo-static",InputDirectory);
	DIR *ufo=opendir(cbuf);
	if(ufo==0) printf("Error: can not find UFO static files\n");
	else
	{
		struct dirent *dp;
		while((dp=readdir(ufo))!=NULL)
		{
			if(dp->d_type==8)
			{
				char ibuf[1000], obuf[1000];
				char buf[100];
				unsigned long int brd;
				FILE *fi, *fo;
				sprintf(ibuf,"%s/ufo-static/%s",InputDirectory,dp->d_name);
				sprintf(obuf,"%s/%s",OutputDirectory?OutputDirectory:".",dp->d_name);
				fi=fopen(ibuf,"rb"); if(fi==0) continue;
				fo=fopen(obuf,"wb"); if(fo==0) {fclose(fi); continue;}
				while((brd=fread(buf,1,100,fi))!=0) fwrite(buf,1,brd,fo);
				fclose(fi);
				fclose(fo);
			}
		}
		closedir(ufo);
	}
	}
	FreeAtomic(l);
	
	return;
	}

void WriteLagrFile(int fno, char *name)
	{
	FILE *f;
	List l;
	int i;
	
	longest_lpline=0;
	longest_cfline=0;
	
	RegisterLine("WLF: counting particles");
	if(write_all_vertices && lagr_hash)
	{
		int mx=1;
		List ll;
	
		for(i=0;i<LagrHashSize;i++)
		for(ll=lagr_hash[i];ll;ll=ListTail(ll))
		{
			List pl;
			pl=CompoundArg1(ListFirst(ll));
			if(ListLength(pl)>mx)
				{
				/*WriteTerm(pl);puts("");*/
				mx=ListLength(pl);
				}
		}
		LagrPrtclsNo=mx;
	}
	UnregisterLine();
	
	if(TexOutput)
		{
		texWriteLagr(fno);
		return;
		}
	if(FAOutput)
		{
		FAWriteLagr(fno);
		return;
		}
	if(OutputDirectory!=NULL)
		sprintf(wpbuf,"%s/lgrng%d.mdl",OutputDirectory,fno);
	else
		sprintf(wpbuf,"lgrng%d.mdl",fno);
	f=fopen(wpbuf,"w");
	if(f==NULL)
		{
		printf("Can not open file \'%s\' for writing.\n",wpbuf);
		perror("");
		return;
		}

	fprintf(f,"%s\n Lagrangian \n",name);
	for(i=1;i<=LagrPrtclsNo;i++)
		{
		fprintf(f,"P%d",i);
		WriteBlank(f,max_prt_lenl-2);
		fprintf(f,"|");
		}
	fprintf(f,">   Factor ");
	WriteBlank(f,C_F_WIDTH-12);
	fprintf(f,"<|> dLagrangian/ dA(p1) dA(p2) dA(p3)");
	WriteBlank(f,L_P_WIDTH-33);
	fprintf(f,"<|\n");
	RegisterLine("WLF: Getting vertices list.");
	l=all_vert_list();
	UnregisterLine();
	RegisterLine("WLF: Sorting lagrangian.");
	if(opSortLagr)
		l=SortedList(l,sort_lagr);
	UnregisterLine();
	RegisterLine("WLF: checking EM charge.");
	check_em_charge(l);
	UnregisterLine();
	WriteLgrngn(l,f);
	return;
	}



Term WriteLagr(Term t, Term ind)
	{
	List l;
	l=all_vert_list();
	if(opSortLagr)
		l=SortedList(l,sort_lagr);
	WriteLgrngn(l,stdout);
	return 0;
	}

Term SelectVertices(Term t, Term ind)
	{
	List l;
	FILE *f;
	int sv;
	int  afl=0, delfl=0;
	int i;
	
	if(lagr_hash==NULL)
		return 0;
	
	if(!is_compound(t) || CompoundArity(t)<2
	   || !is_atom(CompoundArg1(t)) || !is_list(CompoundArg2(t)))
	   	{
	   	ErrorInfo(340);
	   	puts("wrong arguments in SelectVertices statement");
	   	return 0;
	   	}

	 for(i=3;i<=CompoundArity(t);i++)
	 	{
	 	if(CompoundArgN(t,i)==NewAtom("WithAnti",0))
	 		{
	 		afl=1;
	 		continue;
	 		}
	 	if(CompoundArgN(t,i)==NewAtom("Delete",0))
	 		{
	 		delfl=1;
	 		continue;
	 		}
	 	ErrorInfo(342);
	 	printf("SelectVertices: unknown option '");
	 	WriteTerm(CompoundArgN(t,i));
	 	puts("'");
	 	}
	 
	 l=sel_vrt(CompoundArg2(t),afl,delfl);
	 if(l==0)
	 	{
	 	if(access(AtomValue(CompoundArg1(t)),0)==0)
	 		unlink(AtomValue(CompoundArg1(t)));
	 	goto exi;
	 	}

	 f=fopen(AtomValue(CompoundArg1(t)),"w");

	 if(f==NULL)
	 	{
	 	ErrorInfo(341);
	 	printf("SelectVertices: can not open file %s for writing\n",
	 		AtomValue(CompoundArg1(t)));
	 	perror("");
		WriteTerm(t);puts("");
	 	return 0;
	 	}


	if(write_all_vertices)
	{
		int mx=1;
		List ll;
	
		for(ll=l;ll;ll=ListTail(ll))
		{
			List pl;
			pl=CompoundArg1(ListFirst(ll));
			if(ListLength(pl)>mx)
				mx=ListLength(pl);
		}
		LagrPrtclsNo=mx;
	}
	
	 sv=write_all_vertices;
	 write_all_vertices=1;

	 if(TexOutput)
	 	tex_write_lagr(l,f);
	 else
	 	WriteLgrngn(l,f);
	 write_all_vertices=sv;
	 fclose(f);
exi:

	 return 0;
	 }

static List add_fct(List m_l, Term p, Term m2)
	{
	List l1;
	Term t;
	l1=m_l;
	while(!is_empty_list(l1))
		{
		if(EqualTerms(CompoundArg1(ListFirst(l1)),p))
			{
			List ml;
			ml=ConsumeCompoundArg(ListFirst(l1),2);
			ml=AppendLast(ml,m2);
			SetCompoundArg(ListFirst(l1),2,ml);
			return m_l;
			}
		l1=ListTail(l1);
		}
	t=MakeCompound(OPR_DIV,4);
	SetCompoundArg(t,1,p);
	SetCompoundArg(t,2,AppendLast(NewList(),m2));
	return AppendLast(m_l,t);
	}

static int wrt_sln(List mm, FILE *f, int b_l, int is_0)
	{
	List s;
	int first=1;
	int cf;
	int w=0;
	s=mm;
	while(!is_empty_list(s))
		{
		Term m2;
		List lp;
		int ast=0;
		m2=ListFirst(s);
		cf=(int)IntegerValue(CompoundArg1(m2));
		lp=CompoundArg2(m2);
		if(!first && b_l)
			{
			int j;
			fprintf(f,"\n");
			for(j=0;j<LagrPrtclsNo;j++)
				fprintf(f,"     |");
			WriteBlank(f,C_F_WIDTH+1);
			if(is_0)
				fprintf(f,"        ");
			}
		if(cf<0)
			{
			w+=fprintf(f,"-");
			cf=-cf;
			}
		else
			if(!first)
				 w+=fprintf(f,"+");
		first=0;
		if(cf!=1 || (is_empty_list(lp)))
			{
			w+=fprintf(f,"%d",cf);
			ast=1;
			}
		while(!is_empty_list(lp))
			{
			int po;
			po=(int)IntegerValue(CompoundArg2(ListFirst(lp)));
			if(ast)
				w+=fprintf(f,"*");
			else
				ast=1;
			w+=fprintf(f,"%s",AtomValue(CompoundArg1(ListFirst(lp))));
			if(po!=1)
				w+=fprintf(f,"%s%d",ChepVersion==3?"**":"^",po);
			lp=ListTail(lp);
			}
		s=ListTail(s);
		}
	return w;
	}

static int exc_cnf(List ml)
{
	int cf;
	List l;
	cf=(int)IntegerValue(CompoundArg1(ListFirst(ml)));
	if(cf<0)
		cf=-cf;
	for(l=ListTail(ml);l;l=ListTail(l))
		cf=(int)gcf(cf,IntegerValue(CompoundArg1(ListFirst(l))));
	if(IntegerValue(CompoundArg1(ListFirst(ml)))<0)
		cf=-cf;
	for(l=ml;l;l=ListTail(l))
		SetCompoundArg(ListFirst(l),1,
				NewInteger(IntegerValue(CompoundArg1(ListFirst(l)))/cf));
	return cf;
}

static List exc_csf(List ml)
{
	List cf;
	List l;
	cf=NewList();
	for(l=CompoundArg2(ListFirst(ml));l;l=ListTail(l))
	{
		char *v;
		v=AtomValue(CompoundArg1(ListFirst(l)));
		if( (v[0]!='G' || (v[1]!='5' && v[1]!='(')) && 
			(v[0]!='(' && v[3]!='G' && v[4]!='5'))
				cf=AppendLast(cf,CopyTerm(ListFirst(l)));
	}
	
	for(l=ListTail(ml);l;l=ListTail(l))
	{
		List l1;
	xyz:
		for(l1=cf;l1;l1=ListTail(l1))
		{
			List l2;
			for(l2=CompoundArg2(ListFirst(l));l2;l2=ListTail(l2))
			{
				if(CompoundArg1(ListFirst(l1))==CompoundArg1(ListFirst(l2)))
				{
					if(IntegerValue(CompoundArg2(ListFirst(l1)))>
							IntegerValue(CompoundArg2(ListFirst(l2))))
						SetCompoundArg(ListFirst(l1),2,CompoundArg2(ListFirst(l2)));
					break;
				}
			}
			if(l2==0)
			{
				cf=CutFromList(cf,l1);
				goto xyz;
			}
		}
	}
	
	if(cf==0)
		return cf;

	for(l=ml;l;l=ListTail(l))
	{
		List l1;
		List pl;
		pl=ConsumeCompoundArg(ListFirst(l),2);
		for(l1=cf;l1;l1=ListTail(l1))
		{
			List l2;
			for(l2=pl;l2;l2=ListTail(l2))
			if(CompoundArg1(ListFirst(l2))==CompoundArg1(ListFirst(l1)))
			{
				if(CompoundArg2(ListFirst(l1))==CompoundArg2(ListFirst(l2)))
					pl=CutFromList(pl,l2);
				else
					SetCompoundArg(ListFirst(l2),2,NewInteger(
						IntegerValue(CompoundArg2(ListFirst(l2)))-
						IntegerValue(CompoundArg2(ListFirst(l1)))));
				break;
			}
		}
		SetCompoundArg(ListFirst(l),2,pl);
	}
	
	return cf;
}

static int sort_wsl(Term k1, Term k2)
{
	List l1,l2;
	int  i1, i2;
	Atom p1, p2;
	
	l1=CompoundArg2(k1);
	l2=CompoundArg2(k2);
	
	i1=ListLength(l1);
	i2=ListLength(l2);
	if(i1>i2)
		return 1;
	if(i1<i2)
		return -1;
	
	for(;l1;l1=ListTail(l1),l2=ListTail(l2))
	{
		p1=CompoundArg1(ListFirst(l1));
		p2=CompoundArg1(ListFirst(l2));
		i1=(int)IntegerValue(CompoundArg2(ListFirst(l1)));
		i2=(int)IntegerValue(CompoundArg2(ListFirst(l2)));
		if(p1!=p2)
			return strcmp(AtomValue(p1),AtomValue(p2));
		if(i1>i2)
			return 1;
		if(i1<i2)
			return -1;
	}
	
	puts("Internal error (lawss)");
	WriteTerm(k1);puts("");
	WriteTerm(k1);puts("");
	return 0;
}
	
static int wrt_scalar(List m2, FILE *ouf)
	{
	List l1,l0;
	List l_m;
	int w=0;
	int first=1;
	int break_lines=0;
	if(opBreakLines==1)
		break_lines=1;


	l1=m2;
	l_m=NewList();
	while(!is_empty_list(l1))
		{
		List lp,lp1=0;
		int prior=1000000;
		
		lp=CompoundArg2(ListFirst(l1));
		
		while(!is_empty_list(lp))
			{
			Term pp;
			pp=GetAtomProperty(CompoundArg1(ListFirst(lp)),OPR_COEFF);
			if(pp && is_integer(pp) && IntegerValue(pp)<prior)
				{
				lp1=lp;
				prior=(int)IntegerValue(pp);
				}
			lp=ListTail(lp);
			}
			
		lp=lp1;	
		if(!is_empty_list(lp))
			{
			Term mm2;
			Atom mp;
			mm2=ListFirst(l1);
			ChangeList(l1,0);
			mp=ListFirst(lp);
			ChangeList(lp,0);
			lp1=ConsumeCompoundArg(mm2,2);
			lp1=CutFromList(lp1,lp);
			SetCompoundArg(mm2,2,lp1);
			l_m=add_fct(l_m,mp,mm2);
			}
		l1=ListTail(l1);
		}

	l0=NewList();
	l1=m2;
	while(!is_empty_list(l1))
		{
		if(ListFirst(l1))
			l0=AppendLast(l0,ListFirst(l1));
		l1=ListTail(l1);
		}


	if(l0)
		{
		w=wrt_sln(l0,ouf,break_lines,0);
		first=0;
		RemoveList(l0);
		}


	for(l1=l_m;l1;l1=ListTail(l1))
	{
		int cn;
		Term cf;
		/*WriteTerm(ListFirst(l1));puts("");*/
		cf=exc_csf(CompoundArg2(ListFirst(l1)));
		SetCompoundArg(ListFirst(l1),4,cf);
		cf=ConsumeCompoundArg(ListFirst(l1),2);
		SetCompoundArg(ListFirst(l1),2,SortedList(cf,sort_wsl));
		cn=exc_cnf(CompoundArg2(ListFirst(l1)));
		SetCompoundArg(ListFirst(l1),3,NewInteger(cn));
		/*WriteTerm(ListFirst(l1));puts("\n");*/
	}
	
	
	l1=l_m;
	
	while(l1)
		{
		int c_n_f, has_sp;
		List l2,l3, c_s_f;
		
		List sim_te=NewList();
		
		if(ListFirst(l1)==0)
		{
			l1=ListTail(l1);
			continue;
		}
		
		l3=CompoundArg2(ListFirst(l1));
		
		for(l2=ListTail(l1);l2;l2=ListTail(l2))
		{
			if(EqualTerms(l3,CompoundArg2(ListFirst(l2))))
				sim_te=AppendLast(sim_te,l2);
		}
		
/*		printf("%d\n",ListLength(sim_te));*/
		
	/*	c_n_f=exc_cnf(l3);
		c_s_f=exc_csf(l3);
	*/
		c_n_f=(int)IntegerValue(CompoundArgN(ListFirst(l1),3));
		c_s_f=ConsumeCompoundArg(ListFirst(l1),4);
		
	/*	printf("%d*",c_n_f);WriteTerm(c_s_f);printf("*");WriteTerm(l3);puts("");*/
		
		has_sp=1;
		if(ListLength(l3)==1 && is_empty_list(CompoundArg2(ListFirst(l3))))
			has_sp=0;
		
		if(!first && break_lines)
			{
			int j;
			fprintf(ouf,"\n");
			for(j=0;j<LagrPrtclsNo;j++)
				fprintf(ouf,"     |");
			WriteBlank(ouf,C_F_WIDTH+1);
			fprintf(ouf,"|");
			}
			
		if(sim_te)
		{
			List i1,i2;
			List m2l;
			Term m2;
			m2=MakeCompound(A_MTERM,3);
			SetCompoundArg(m2,1,NewInteger(c_n_f));
			SetCompoundArg(m2,2,
					AppendFirst(c_s_f,ConsumeCompoundArg(ListFirst(l1),1)));
			m2l=AppendLast(NewList(),m2);
			for(i1=sim_te;i1;i1=ListTail(i1))
			{
				i2=ListFirst(i1);
				c_n_f=(int)IntegerValue(CompoundArgN(ListFirst(i2),3));
				c_s_f=ConsumeCompoundArg(ListFirst(i2),4);
				m2=MakeCompound(A_MTERM,3);
				SetCompoundArg(m2,1,NewInteger(c_n_f));
				SetCompoundArg(m2,2,
						AppendFirst(c_s_f,ConsumeCompoundArg(ListFirst(i2),1)));
				FreeAtomic(ListFirst(i2));
				m2l=AppendLast(m2l,m2);
				ChangeList(i2,0);
			}
			RemoveList(sim_te);
			if(!first)
				w+=fprintf(ouf,"+");
			w+=fprintf(ouf,"(");
			w+=wrt_sln(m2l,ouf,break_lines,1);
			/*wrt_sln(m2l,stdout,break_lines,1);puts("");*/
			w+=fprintf(ouf,")");
			if(has_sp)
				w+=fprintf(ouf,"*");
			first=0;
		}
		else
		{
		
			if(c_n_f<-1 || (c_n_f >1 && first))
				w+=fprintf(ouf,"%d*",c_n_f);
			if(c_n_f == -1)
				w+=fprintf(ouf,"-");
			if(c_n_f == 1 && !first)
				w+=fprintf(ouf,"+");
			if(c_n_f>1 && !first)
				w+=fprintf(ouf,"+%d*",c_n_f);
			first=0;

			c_s_f=AppendFirst(c_s_f,ConsumeCompoundArg(ListFirst(l1),1));
			for(l2=c_s_f;l2;l2=ListTail(l2))
			{
				w+=fprintf(ouf,"%s",AtomValue(CompoundArg1(ListFirst(l2))));
				if(IntegerValue(CompoundArg2(ListFirst(l2)))>1)
					w+=fprintf(ouf,"%s%ld",ChepVersion==3?"**":"^",
							IntegerValue(CompoundArg2(ListFirst(l2))));
				if(ListTail(l2) || has_sp)
				w+=fprintf(ouf,"*");
			}
			FreeAtomic(c_s_f);

		}
		
		if(break_lines && ListLength(l3)>1)
			{
			int j;
			fprintf(ouf,"\n");
			for(j=0;j<LagrPrtclsNo;j++)
				fprintf(ouf,"     |");
			WriteBlank(ouf,C_F_WIDTH+1);
			WriteBlank(ouf,7);
			}

		if(ListLength(l3)>1)
			w+=fprintf(ouf,"(");
		if(has_sp)
			w+=wrt_sln(l3,ouf,break_lines,1);
		if(ListLength(l3)>1)
			w+=fprintf(ouf,")");
		l1=ListTail(l1);
		}
	FreeAtomic(l_m);
	return w;
	}


static void write_l_line(FILE *f, Term a2)
	{

	
	if(CompoundArg2(a2)==NewInteger(0) ||
		 is_empty_list(CompoundArgN(a2,5)))
		return;

	
	{
	List pl;
	int i,ac,w;
	pl=CompoundArg1(a2);
	ac=ListLength(pl);
	if(ac<=2 && f!=stdout && !write_all_vertices)
		return;
	for(i=1;i<=LagrPrtclsNo;i++)
		{
		if(i>ac)
			{
			WriteBlank(f,max_prt_lenl);
			fprintf(f,"|");
			}
		else
			{
			Atom a,a1;
			a=CompoundArg1(ListNth(pl,i));
			a1=GetAtomProperty(a,A_CHNAME);
			if(a1) a=a1;
			w=fprintf(f,"%s",AtomValue(a));
			WriteBlank(f,max_prt_lenl-w);
			fprintf(f,"|");
			}
		}
	}

	{
	int num,den,w;
	List s,sn,sd;
	num=(int)IntegerValue(CompoundArg1(CompoundArg2(a2)));
	den=(int)IntegerValue(CompoundArg2(CompoundArg2(a2)));

	sn=sd=NewList();
	for(s=CompoundArgN(a2,3);s;s=ListTail(s))
		{
		if(IntegerValue(CompoundArg2(ListFirst(s)))>0)
			sn=AppendLast(sn,ListFirst(s));
		else
			sd=AppendLast(sd,ListFirst(s));
		}
	if(den!=1)
		sd=AppendFirst(sd,NewInteger(den));
	if(num!=1 && num!=-1)
		sn=AppendFirst(sn,NewInteger(num));

	w=0;
	if(num==-1)
		w=fprintf(f,"-");
	if(is_empty_list(sn))
		w+=fprintf(f,"1");
	for(s=sn;s;s=ListTail(s))
		{
		Atom p;
		int po;
		if(is_integer(ListFirst(s)))
			{
			w+=fprintf(f,"%ld",IntegerValue(ListFirst(s)));
			if(ListTail(s))
				w+=fprintf(f,"*");
			}
		else
			{
			p=CompoundArg1(ListFirst(s));
			po=(int)IntegerValue(CompoundArg2(ListFirst(s)));
			w+=fprintf(f,"%s",AtomValue(p));
			if(po!=1)
				w+=fprintf(f,"%s%d",ChepVersion==3?"**":"^",po);
			if(ListTail(s))
				w+=fprintf(f,"*");
			}
		}
	if(!is_empty_list(sd))
		{
		int lle;
		lle=ListLength(sd);
		w+=fprintf(f,"/");
		if(lle>1)
			w+=fprintf(f,"(");
		for(s=sd;s;s=ListTail(s))
			{
			Atom p;
			int po;
			if(is_integer(ListFirst(s)))
				{
				w+=fprintf(f,"%ld",IntegerValue(ListFirst(s)));
				if(ListTail(s))
					w+=fprintf(f,"*");
				}
			else
				{
				p=CompoundArg1(ListFirst(s));
				po=-(int)IntegerValue(CompoundArg2(ListFirst(s)));
				w+=fprintf(f,"%s",AtomValue(p));
				if(po!=1)
					w+=fprintf(f,"%s%d",ChepVersion==3?"**":"^",po);
				if(ListTail(s))
					w+=fprintf(f,"*");
				}
			}
		if(lle>1)
			w+=fprintf(f,")");
		}
	RemoveList(sd);
	RemoveList(sn);
	if(w>longest_cfline)
		longest_cfline=w;
	WriteBlank(f,C_F_WIDTH-w);
	fprintf(f,"|");

	}

	{
	List lm,l2,l3;
	for(lm=CompoundArgN(a2,5);lm;lm=ListTail(lm))
		{
		Term m2;
		m2=ListFirst(lm);
		l2=ConsumeCompoundArg(m2,2);
		for(l3=CompoundArgN(m2,3);l3;l3=ListTail(l3))
		{
			List l4;
			Atom ts=ListFirst(l3);
			for(l4=l2;l4;l4=ListTail(l4))
				if(CompoundArg1(ListFirst(l4))==ts)
					break;
			if(l4)
				SetCompoundArg(ListFirst(l4),2,	NewInteger(
						IntegerValue(CompoundArg2(ListFirst(l4)))+1));
			else
			l2=AppendLast(l2,
				MakeCompound2(OPR_POW,ListFirst(l3),NewInteger(1)));
		}
		l3=ConsumeCompoundArg(m2,3);
		FreeAtomic(l3);
		SetCompoundArg(m2,2,l2);
		}
	}

	/*DumpList(CompoundArgN(a2,5));*/


	{
	int w;
	w=wrt_scalar(CompoundArgN(a2,5),f);
	if(w>longest_lpline)
		longest_lpline=w;
	}

	fprintf(f,"\n");

	}


extern int UFOutput;

void WriteLgrngn(Term l, FILE *fout)
	{
	List li,lj,lg=0;
	Term prp;
	char cbuf[64];
	int pn, pni=0;
	pn=ListLength(l);
	
	prp=GetAtomProperty(NewAtom("chepOutput",0),NewAtom("cfWidth",0));
	if(is_integer(prp))
		C_F_WIDTH=(int)IntegerValue(prp);
/*
	if(MultByI)
	{
		for(li=l;li;li=ListTail(li))
			alg2_multbyi(ListFirst(li));
	}
*/
	li=l;


	for(li=l;!is_empty_list(li);li=ListTail(li))
		{
		Term a2;
		List a2l;
		int sf=0;
		
		pni++;
		
		a2=ListFirst(li);
		if(CompoundArgN(a2,5)==0)
			continue;
		sprintf(cbuf,"Writing lagrangian line %d of %d.\n",pni,pn);
		RegisterLine(cbuf);
		
		if(!write_all_vertices && ListLength(CompoundArg1(a2))>4)
		{
			static int repno=0;
			if(repno<10)
			{
				printf("Vertex ");
				WriteVertex(CompoundArg1(a2));
				printf(" with more than 4 particles.\n");
			}
			if(repno==10)
				puts("More vertices with more than 4 particles follow");
			repno++;
			UnregisterLine();
			continue;
		}

		if(need_col_rdc(a2))
		{
			List l1,l2,l3=0;
			sf=1;

			l1=alg2_denorm(CopyTerm(a2));
			
			for(l2=l1;l2;l2=ListTail(l2))
				l3=ConcatList(l3,color_reduce(ListFirst(l2)));
			RemoveList(l1);
			a2l=l3;
			for(l1=a2l;l1;l1=ListTail(l1))
				alg2_norm(ListFirst(l1));
		}
		else if(need_spin_rdc(a2))
		{
			List l1,l2,l3=0;
			sf=1;
			
			l1=alg2_denorm(CopyTerm(a2));
			
			for(l2=l1;l2;l2=ListTail(l2))
				l3=ConcatList(l3,spinor_reduce(ListFirst(l2)));
			RemoveList(l1);
			a2l=l3;
			for(l1=a2l;l1;l1=ListTail(l1))
				alg2_norm(ListFirst(l1));
			
		}
		else
			a2l=MakeList1(a2);
		
		
			
		for(lj=a2l;lj;lj=ListTail(lj))
		{
			Term a2;
			a2=ListFirst(lj);
			if(is_atom(CompoundArg1(a2)))
				continue;
						
			alg2_symmetrize(a2);
			
			
			alg2_common_s(a2);

			alg2_common_n(a2);
			
			/*alg2_red_cos(a2);

			alg2_red_orth(a2);

			if(CompoundArgN(a2,5)==0)
				continue;
			alg2_red_sico(a2);

			alg2_red_comsico(a2);*/

			alg2_recommon_n(a2);

			alg2_recommon_s(a2);
			
			if(MultByI) alg2_multbyi(a2);
			alg2_red_1pm5(a2);
			
			if(opAbbrVrt)
			{
				alg2_decommon_s(a2);
				alg2_abbr_vrt(a2);
				alg2_common_s(a2);
			}

			alg2_recommon_n(a2);

			if(opEvalVrt)
				alg2_eval_vrt(a2);

/*			alg2_reduce(a2);

						
			write_l_line(fout,a2);*/
			
			lg=AppendLast(lg,a2);

		}
		
		if(!sf)
			RemoveList(a2l);
        
		UnregisterLine();

		}

	RemoveList(l);
	
	alg2_setcls(lg);
	
	for(li=lg;li;li=ListTail(li))
	{
		
		alg2_reduce(ListFirst(li));
		write_l_line(fout,ListFirst(li));
	}
	
	FreeAtomic(lg);

	}

static int sort_lagr(Term a1, Term a2)
	{
	List l1,l2;
	l1=CompoundArg1(a1);
	l2=CompoundArg1(a2);
	/*WriteTerm(l1); WriteTerm(l2); puts("comp");*/
		{
		int e1,e2;
		e1=ListLength(l1);
		e2=ListLength(l2);
		if(e1!=e2)
			{
			if(e1>e2)
				return 1;
			else
				return -1;
			}
		}
	while(!is_empty_list(l1))
		{
		Term e1,e2;
		e1=CompoundArg1(ListFirst(l1));
		e2=CompoundArg1(ListFirst(l2));
		if(e1!=e2)
			return strcmp(AtomValue(e1), AtomValue(e2));
		l1=ListTail(l1);
		l2=ListTail(l2);
		}
	return 0;
	}
	
extern void prt2cls(Atom *);

static int sort_lagr_fa(Term a1, Term a2)
	{
	List l1,l2;
	l1=CompoundArg1(a1);
	l2=CompoundArg1(a2);
	/*WriteTerm(l1); WriteTerm(l2); puts("comp");*/
		{
		int e1,e2;
		e1=ListLength(l1);
		e2=ListLength(l2);
		if(e1!=e2)
			{
			if(e1>e2)
				return 1;
			else
				return -1;
			}
		}
	for(;l1;l1=ListTail(l1),l2=ListTail(l2))
	{
		Atom p1, p2, s1, s2;
		p1=CompoundArg1(ListFirst(l1));
		s1=CompoundName(CompoundArg2(ListFirst(l1)));
		if(s1==OPR_SCALAR && GetAtomProperty(p1,A_GRASS))
			s1=A_GRASS;
		p2=CompoundArg1(ListFirst(l2));
		s2=CompoundName(CompoundArg2(ListFirst(l2)));
		if(s2==OPR_SCALAR && GetAtomProperty(p2,A_GRASS))
			s2=A_GRASS;
		if(s1!=s2)
			return strcmp(AtomValue(s1),AtomValue(s2));
	}
		
	l1=CompoundArg1(a1);
	l2=CompoundArg1(a2);
	
	while(!is_empty_list(l1))
		{
		Term e1,e2;
		e1=CompoundArg1(ListFirst(l1));
		e2=CompoundArg1(ListFirst(l2));
		prt2cls(&e1);
		prt2cls(&e2);
		if(e1!=e2)
			return strcmp(AtomValue(e1), AtomValue(e2));
		l1=ListTail(l1);
		l2=ListTail(l2);
		}
	return 0;
	}

Term ProcSortL(Term t, Term ind)
	{
	FreeAtomic(t);
	return 0;
	}



void Write2Vertex(FILE *f, Term prt)
	{
	List l,l1;
	l1=l=all_vert_list();
	if(opSortLagr)
		l1=l=SortedList(l,sort_lagr);
	while(!is_empty_list(l))
		{
		List pl;
		pl=CompoundArg1(ListFirst(l));
		if(ListLength(pl)==2)
			{
			Atom a1, a2;
			a1=CompoundArg1(ListFirst(pl));
			a2=CompoundArg1(ListNth(pl,2));
			if((a1==CompoundArg1(prt) && a2==CompoundArg2(prt)) ||
				(a2==CompoundArg1(prt) && a1==CompoundArg2(prt)))
				{
				Term a2;
				a2=CopyTerm(ListFirst(l));
				alg2_common_s(a2);
				alg2_common_n(a2);
				alg2_common_t(a2);
				alg2_red_cos(a2);
				alg2_red_orth(a2);
				tex_wrt_2vrt(f,a2);
				return;
				}
			}
		l=ListTail(l);
		}
	FreeAtomic(l1);
	fprintf(f," 0 ");
	}

static int chk_sel(List sel, List prt, int typ, int afl)
	{
	List l;
	char c_used[4];
	Atom p,ap;
	Term prp;
	int i;
	c_used[0]=c_used[1]=c_used[2]=c_used[3]=0;

	l=prt;
	while(!is_empty_list(l))
		{
		int cprp=0;
		p=CompoundArg1(ListFirst(l));
		prp=GetAtomProperty(p,PROP_TYPE);
		if((CompoundName(prp)==OPR_PARTICLE
			&& CompoundArgN(prp,7)==OPR_MLT)
			|| (CompoundName(prp)==OPR_FIELD
				&& CompoundArg2(prp)==NewInteger(5)))
				cprp=1;
		if(CompoundName(prp)==OPR_FIELD
			&& CompoundArg2(prp) == NewInteger(4))
				p=CompoundArg1(prp);
		ap=GetAtomProperty(p,A_ANTI);
		if(typ==0)
			{
			int tf;
			tf=ListMember(sel,p);
			if(!(tf || (afl && ListMember(sel,ap))))
				return 0;
			}
		else
			{
			List l1;
			int cno=0;
			l1=sel;
			while(!is_empty_list(l1))
				{
				if(c_used[cno])
					{
					l1=ListTail(l1);
					cno++;
					continue;
					}
				if(ListFirst(l1)==OPR_MLT && cprp)
					{
					c_used[cno]=1;
					break;
					}
				if(is_list(ListFirst(l1)) &&
					(ListMember(ListFirst(l1),p) ||
					(afl && ListMember(ListFirst(l1),ap))))
					{
					c_used[cno]=1;
					break;
					}
				cno++;
				l1=ListTail(l1);
				}
			if(is_empty_list(l1))
				return 0;
			}

		l=ListTail(l);
		}
	if(typ==1)
		{
		for(i=0;i<ListLength(sel);i++)
			if(!c_used[i])
				return 0;
		}
	return 1;
	}

List sel_vrt(List sel, int afl, int rmfl)
	{
	List ret;
	List l1,l2;
	int typ, i;
	ret=NewList();
	typ=0;
	l1=sel;
	while(!is_empty_list(l1))
		{
		if(!is_atom(ListFirst(l1)) || ListFirst(l1)==OPR_MLT)
			{
			typ=1;
			break;
			}
		l1=ListTail(l1);
		}
	if(typ)
		{
		l1=sel;
		while(!is_empty_list(l1))
			{
			if(ListFirst(l1)==OPR_MLT)
				{
				l1=ListTail(l1);
				continue;
				}
			if(is_atom(ListFirst(l1)))
				{
				ChangeList(l1,AppendLast(NewList(),ListFirst(l1)));
				l1=ListTail(l1);
				continue;
				}
			if(!is_compound(ListFirst(l1)))
				{
				ErrorInfo(344);
				puts("SelectVertices: wrong format of particle list.");
				return 0;
				}
			l2=Oper1ToList(ListFirst(l1),CompoundName(ListFirst(l1)));
			ChangeList(l1,l2);
			while(!is_empty_list(l2))
				{
				if(!is_atom(ListFirst(l2)))
					{
					ErrorInfo(345);
					puts("SelectVertices: wrong format of particle list.");
					return 0;
					}
				l2=ListTail(l2);
				}
			l1=ListTail(l1);
			}
		}
	if(typ==1 && ListLength(sel)>4)
		{
		ErrorInfo(347);
		puts("SelectVertex: too many particle clusters.");
		return 0;
		}

	for(i=0;i<LagrHashSize;i++)
	{
		l1=lagr_hash[i];
		while(!is_empty_list(l1))
			{
			List l12;
			int trfl;
			trfl=chk_sel(sel,CompoundArg1(ListFirst(l1)),typ,afl);
			l12=l1;
			l1=ListTail(l1);
			if(trfl)
				{
				ret=AppendLast(ret,CopyTerm(ListFirst(l12)));
				if(rmfl)
					{

					lagr_hash[i]=CutFromList(lagr_hash[i],l12);
					}
				}
			}
	}
	/*
	if(is_empty_list(ret))
		{
		ErrorInfo(346);
		printf("SelectVertices: no vertices found for ");
		WriteTerm(sel);
		puts(" pattern");
		}
	*/

	return SortedList(ret,sort_lagr);
	}

Term ProcLongestLine(Term t, Term ind)
	{
	printf("Longest lagrangian line: %d symbols\n",longest_lpline);
	return 0;
	}

extern time_t tm_rl;
extern int a2_mono_no;

extern int th_no;

#ifdef MTHREAD
extern pthread_key_t TermsKey;
#endif

static void *frl(void *d)
{
	int sh, st, i, ind;
	List l;
	
	ind=sh= *(int *)d;
	
#ifdef MTHREAD
	if(ind>=0)
		pthread_setspecific(TermsKey,&ind);
#endif
	
	if(sh<0) sh=0;
	
	st=th_no;
	if(st<1)
		st=1;
	
	for(i=sh;i<LagrHashSize;i+=st)
	{
	for(l=lagr_hash[i];l;l=ListTail(l))
		{
		Term a2;
		
		a2=ListFirst(l);
		
		alg2_common_n(a2);
		alg2_common_s(a2);
		alg2_red_cos(a2);
		alg2_red_orth(a2);
		/*if(ListLength(CompoundArg1(a2))>2)*/
		{
		alg2_recommon_s(a2);
		alg2_red_sico(a2);
		alg2_red_comsico(a2);
		}
		alg2_decommon_n(a2);
		alg2_decommon_s(a2);

		}
	}
    return 0;
}


Term ProcReduceLagr(Term t, Term ind)
	{
	int vrb=0;
	List l;
	int i, cnt=0;
	int ll_1=0, ll_2=0, ll_m=0;
	char cbuf[64];
	time_t stt;
	if(lagr_hash==NULL)
		return 0;
/*	puts("RL");*/
	if(lagr_reduced)
		return 0;
	stt=clock();
	lagr_reduced=1;
	
	if(!is_atom(t))
		vrb=1;
	FreeAtomic(t);
			
	for(i=0;i<LagrHashSize;i++)
	{
		int ll;
		ll=ListLength(lagr_hash[i]);
		cnt+=ll;
		ll_1+=ll;
		ll_2+=ll*ll;
		if(ll>ll_m)
			ll_m=ll;
	}
	
	if(VerbMode)
	{
		printf("ReduceLagr: %d vertices in %d lists, ave %d, sigma %d, max %d; ",
				cnt, LagrHashSize, ll_1/LagrHashSize,
				(ll_2-ll_1*ll_1/(int)LagrHashSize)/(int)LagrHashSize, ll_m);
		fflush(stdout);
	}

#ifdef MTHREAD
	if(th_no>0 && a2_mono_no>1000)
	{
		pthread_t thid[16];
		int thsh[16];
		for(i=0;i<th_no;i++)
		{
			thsh[i]=i;
			pthread_create(&thid[i],NULL,frl,(void *)&thsh[i]);
		}
		for(i=0;i<th_no;i++)
			pthread_join(thid[i],NULL);
	}
	else
#endif
	{
	int sh=-1;
	frl((void *)&sh);
	}	
	
	/*for(i=0;i<LagrHashSize;i++)
	{
	for(l=lagr_hash[i];l;l=ListTail(l))
		{
		Term a2;
		a2=ListFirst(l);
		
		alg2_common_n(a2);
		alg2_common_s(a2);
		alg2_red_cos(a2);
		alg2_red_orth(a2);
		if(ListLength(CompoundArg1(a2))>2)
		{
		alg2_recommon_s(a2);
		alg2_red_sico(a2);
		alg2_red_comsico(a2);
		}
		alg2_decommon_n(a2);
		alg2_decommon_s(a2);

		}
	}*/
	if(vrb)
		{
		printf("PRLst: %d vrts, ",cnt);
		AtomStat1();
		ListStat1();
		}
	RegisterLine("ReduceLagr: cutting empty vertices.");
	
	for(i=0;i<LagrHashSize;i++)
	{
		l=lagr_hash[i];
		while(!is_empty_list(l))
		{
		List l1;
		if(CompoundArgN(ListFirst(l),5)==0)
			{
			l1=l;
			l=ListTail(l);
			FreeAtomic(ListFirst(l1));
			ChangeList(l1,0);
			lagr_hash[i]=CutFromList(lagr_hash[i],l1);
			cnt--;
			}
		else
			l=ListTail(l);
		}
	}

	UnregisterLine();
	
	if(VerbMode)
		printf("%d non-zero\n",cnt);
	
	if(vrb)
		{
		printf("PRLfi: %d vrts, ",cnt);
		AtomStat1();
		ListStat1();
		}
	tm_rl=(clock()-stt);
	return 0;
	}

Term ProcSaveLagr(Term t, Term ind)
	{
	List l;
	int i;
	
	if(!is_compound(t) || !is_atom(CompoundArg1(t)))
		{
		ErrorInfo(480);
		printf("SaveLagr: bad arguments\n");
		return 0;
		}
	if(!itrSetOut(AtomValue(CompoundArg1(t))))
		{
		ErrorInfo(481);
		printf("SaveLagr: can not open output file\n");
		perror("");
		return 0;
		}

	ProcReduceLagr(MakeCompound1(A_I,A_I),0);

	for(i=0;i<LagrHashSize;i++)
	{
		l=lagr_hash[i];
		while(!is_empty_list(l))
			{
			itrOut(ListFirst(l));
			l=ListTail(l);
			}
		FreeAtomic(lagr_hash[i]);
		lagr_hash[i]=NewList();
	}
	
	itrCloseOut();

	return 0;
	}

Term ProcLoadLagr(Term t, Term ind)
	{
	Term a2;

	if(!is_compound(t) || !is_atom(CompoundArg1(t)))
		{
		ErrorInfo(482);
		printf("LoadLagr: bad arguments\n");
		return 0;
		}

	if(!itrSetIn(AtomValue(CompoundArg1(t))))
		{
		ErrorInfo(483);
		printf("LoadLagr: can not open input file\n");
		perror("");
		return 0;
		}

	while((a2=itrIn()))
		{
		List l1,l2;
		l1=ConsumeCompoundArg(a2,5);
		l2=l1;
		while(!is_empty_list(l1))
			{
			Term a21;
			a21=CopyTerm(a2);
			SetCompoundArg(a21,5,AppendLast(NewList(),ListFirst(l1)));
			alg2_hash_add(lagr_hash,(int)LagrHashSize,AppendLast(NewList(),a21));
			l1=ListTail(l1);
			}
		RemoveList(l2);
		FreeAtomic(a2);
		}
	itrCloseIn();

	return 0;
	}


static int co_no=1;

Term ProcCoeff(Term t, Term ind)
	{
	Term t1;
	List l;
	if(!is_compound(t) || CompoundArity(t)!=1)
		{
		ErrorInfo(578);
		puts("Illegal argument in coeff statement.");
		return 0;
		}
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	t=CommaToList(t1);

	for(l=t;l;l=ListTail(l))
		{
		Atom a;
		a=ListFirst(l);
		if(!is_atom(a))
			{
			ErrorInfo(579);
			printf("wrong argument `");
			WriteTerm(a);
			puts("' in coeff statement");
			}
		else
			SetAtomProperty(a,OPR_COEFF,NewInteger(co_no++));
		}
	FreeAtomic(t);
	return 0;
	}



	
