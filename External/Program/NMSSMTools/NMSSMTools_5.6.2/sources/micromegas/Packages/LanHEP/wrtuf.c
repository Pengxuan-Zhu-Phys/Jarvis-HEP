#include "lanhep.h"
#include <time.h>
#include <string.h>

extern char *eff_infile;
extern char *ModelName;

extern int opEvalVrt, EvalVrt, NoColors, opAbbrVrt, allow_sym_div, write_all_vertices;
extern int FAver;

static char nv=0, nf=0, ns=0, nu=0;

static List gkl=0;
static Atom ili[26]={0,0,0,0,0,0}, imom[16];
static Integer pind[16], wind[16];

static Term conv_lor(List pl, Term m2, int eff);
void alg2_rem_col(Term a2);

static void wrt_expr(FILE *, Term n, List s, List t);
static List cv;

static List cpml=0;
extern List fainclude;

void inf_write_rc(FILE *), cls_write_dmatr(FILE *);
void alg2_fix_uuv(Term), alg2_fix_ff(Term), alg2_red_rc(Term);
void alg2_abbr_vrt(Term), alg2_eval_vrtn(Term);

void cls_wrt_ind(FILE *), cls_wrt_nms(FILE *, List), alg2_setcls(List);
int  cls_prt_info(Term *, Atom *);

extern List cls_lagr_hook;
extern void inf_decl_hc(FILE *), prm_decl_hc(FILE *, List);

static int pcmp(Term mm1, Term mm2)
{
	Term m1=CompoundArg2(mm1);
	Term m2=CompoundArg2(mm2);
	if(m1==0 && m2) return 1;
	if(m2==0 && m1) return -1;
	if(m1==0 && m2==0)
		{
		m1=CompoundArgN(mm1,3);
		m2=CompoundArgN(mm2,3);
		}
	while(is_float(CompoundArg1(ListFirst(m1)))||
			CompoundArg1(ListFirst(m1))==A_I) m1=ListTail(m1);
	while(is_float(CompoundArg1(ListFirst(m2)))||
			CompoundArg1(ListFirst(m2))==A_I) m2=ListTail(m2);
	return strcmp(AtomValue(CompoundArg1(ListFirst(m1))),
				AtomValue(CompoundArg1(ListFirst(m2))));
}
	

static char cnbuf[256];
static char *convname(Atom a)
{
	int i=0, j=0;
	char *s=AtomValue(a);
    int ccf=0;
    Term prop=GetAtomProperty(a,PROP_TYPE);
    if(prop && is_compound(prop) && CompoundName(prop)==OPR_FIELD
       && CompoundArg2(prop)==NewInteger(4)) ccf=1;
	for(i=0;s[i];i++)
	{
		if(s[i]=='~')
		{
			cnbuf[j++]='_';
			cnbuf[j++]='_';
			cnbuf[j++]='t';
			cnbuf[j++]='i';
			cnbuf[j++]='l';
			cnbuf[j++]='d';
			cnbuf[j++]='e';
			cnbuf[j++]='_';
			cnbuf[j++]='_';
			continue;
		}
		if(s[i]=='+')
		{
			cnbuf[j++]='_';
			cnbuf[j++]='_';
			cnbuf[j++]='p';
			cnbuf[j++]='l';
			cnbuf[j++]='u';
			cnbuf[j++]='s';
			cnbuf[j++]='_';
			cnbuf[j++]='_';
			continue;
		}
		if(s[i]=='-')
		{
			cnbuf[j++]='_';
			cnbuf[j++]='_';
			cnbuf[j++]='m';
			cnbuf[j++]='i';
			cnbuf[j++]='n';
			cnbuf[j++]='u';
			cnbuf[j++]='s';
			cnbuf[j++]='_';
			cnbuf[j++]='_';
			continue;
		}
		if(s[i]=='.' && s[i+1]=='c')
		{
            if(ccf) {cnbuf[j]=0; return cnbuf;}
			cnbuf[j++]='g';
			cnbuf[j++]='h';
			cnbuf[j++]=0;
			return cnbuf;
		}
		if(s[i]=='.' && s[i+1]=='f')
		{
			cnbuf[j++]='g';
			cnbuf[j++]='s';
			cnbuf[j++]='b';
			cnbuf[j++]=0;
			return cnbuf;
		}
		if(s[i]=='.' && s[i+1]=='C')
		{
			cnbuf[j++]='g';
			cnbuf[j++]='h';
			cnbuf[j++]='_';
			cnbuf[j++]='_';
			cnbuf[j++]='t';
			cnbuf[j++]='i';
			cnbuf[j++]='l';
			cnbuf[j++]='d';
			cnbuf[j++]='e';
			cnbuf[j++]='_';
			cnbuf[j++]='_';
			cnbuf[j++]=0;
			return cnbuf;
		}
		cnbuf[j++]=s[i];
	}
	cnbuf[j]=0;
	return cnbuf;
}


		
extern List UFOparth,UFOverth;
static int ufo_auto_pdg=1000501;

void UF_write_lagr(List l)
{
	time_t tm;
	FILE *f;
	Term gkl1;
	List l1;
	int firstprt=1;
	char cbuf[256];
	int vno=1;
	
	if(OutputDirectory!=NULL)
		sprintf(cbuf,"%s/particles.py",OutputDirectory);
	else
		sprintf(cbuf,"particles.py");
	f=fopen(cbuf,"w");
	
	fprintf(f,"#\tLanHEP output produced at ");
	time(&tm);
	fprintf(f,"%s",ctime(&tm));
	fprintf(f,"#\tfrom the file '%s'\n",eff_infile);
	if(ModelName)
		fprintf(f,"#\tModel named '%s'\n",ModelName);
	fprintf(f,"\n\n");
	
	
	for(l1=UFOparth;l1;l1=ListTail(l1))
		fprintf(f,"%s\n",AtomValue(ListFirst(l1)));
	fprintf(f,"\n");
	
	for(l1=all_prtc_list();l1;l1=ListTail(l1))
	{
		Term p,col,tnm;
		int dim=0;
		Atom a=ListFirst(l1);
		Atom pdg=NewAtom("_pdg",0);
		Term v;
		p=GetAtomProperty(a,PROP_TYPE);
		if(p==0) continue;
		dim=cls_prt_info(&p, &a);
		if(p==0) continue;
		if(!is_compound(p)) continue;
		col=GetAtomProperty(a,A_COLOR);
		if(col) col=IntegerValue(CompoundArg1(col));
		col++;
		if(col==2) col=3;
		if(col==4) col=8;
		
		/*tnm=GetAtomProperty(a,A_TEXNAME);if(tnm==0)*/ tnm=a;
		if(CompoundName(p)==OPR_PARTICLE)
		{
			if(CompoundArg1(p)!=CompoundArg2(p) && CompoundArg2(p)==a)
				continue;

			if(IntegerValue(CompoundArgN(p,4))==2 && GetAtomProperty(a,A_GHOST))
			{
				Term g =GetAtomProperty(a,A_GHOST);
				Term ag=GetAtomProperty(CompoundArg2(p),A_GHOST);
				nu++;
				
				fprintf(f,"%s =  Particle(  ",convname(CompoundArg1(g)));
				v=GetAtomProperty(CompoundArg1(g),pdg);
				fprintf(f,"pdg_code = %ld,\n",v?IntegerValue(v):0);
				fprintf(f,"              name = '%s',\n",
									AtomValue(CompoundArg1(g)));
				fprintf(f,"              antiname = '%s',\n",
									AtomValue(CompoundArg2(ag)));
				fprintf(f,"              spin = -1,\n");
				fprintf(f,"              color = %d,\n",col);
				v=CompoundArgN(p,5);
				fprintf(f,"              mass = '%s',\n",v?AtomValue(v):"ZERO");
				v=CompoundArgN(p,6);
				fprintf(f,"              width = '%s',\n",v?AtomValue(v):"ZERO");
				v=GetAtomProperty(CompoundArg1(g),A_TEXNAME);
				fprintf(f,"              texname = '%s',\n",
						v?AtomValue(v):AtomValue(CompoundArg1(g)));
				v=GetAtomProperty(CompoundArg2(ag),A_TEXNAME);
				fprintf(f,"              antitexname = '%s',\n",
						v?AtomValue(v):AtomValue(CompoundArg2(ag)));
				v=GetAtomProperty(CompoundArg1(p),A_EM_CHARGE);
				if(v==0 || CompoundArg2(v)==NewInteger(1))
				fprintf(f,"              charge = %ld,\n",
						v?IntegerValue(CompoundArg1(v)):0);
				else
				fprintf(f,"              charge = %ld/%ld,\n",
						IntegerValue(CompoundArg1(v)),
						IntegerValue(CompoundArg2(v)));
				
				fprintf(f,"              LeptonNumber = 0,\n");
				fprintf(f,"              GhostNumber = 1)\n\n");
					
				fprintf(f,"%s = ",convname(CompoundArg2(ag)));
				fprintf(f,"%s.anti()\n\n",convname(CompoundArg1(g)));
				
				if(CompoundArg1(p)!=CompoundArg2(p))
				{
				nu++;
				fprintf(f,"%s =  Particle(  ",convname(CompoundArg2(g)));
				v=GetAtomProperty(CompoundArg2(g),pdg);
				fprintf(f,"pdg_code = %ld,\n",v?IntegerValue(v):0);
				fprintf(f,"              name = '%s',\n",
									AtomValue(CompoundArg2(g)));
				fprintf(f,"              antiname = '%s',\n",
									AtomValue(CompoundArg1(ag)));
				fprintf(f,"              spin = -1,\n");
				fprintf(f,"              color = %d,\n",col);
				v=CompoundArgN(p,5);
				fprintf(f,"              mass = '%s',\n",v?AtomValue(v):"ZERO");
				v=CompoundArgN(p,6);
				fprintf(f,"              width = '%s',\n",v?AtomValue(v):"ZERO");
				v=GetAtomProperty(CompoundArg2(g),A_TEXNAME);
				fprintf(f,"              texname = '%s',\n",
						v?AtomValue(v):AtomValue(CompoundArg1(g)));
				v=GetAtomProperty(CompoundArg1(ag),A_TEXNAME);
				fprintf(f,"              antitexname = '%s',\n",
						v?AtomValue(v):AtomValue(CompoundArg2(ag)));
				v=GetAtomProperty(CompoundArg1(p),A_EM_CHARGE);
				if(v==0 || CompoundArg2(v)==NewInteger(1))
				fprintf(f,"              charge = %ld,\n",
						v?IntegerValue(CompoundArg1(v)):0);
				else
				fprintf(f,"              charge = %ld/%ld,\n",
						IntegerValue(CompoundArg1(v)),
						IntegerValue(CompoundArg2(v)));
				fprintf(f,"              LeptonNumber = 0,\n");
				fprintf(f,"              GhostNumber = 1)\n\n");
				
				fprintf(f,"%s = ",convname(CompoundArg1(ag)));
				fprintf(f,"%s.anti()\n\n",convname(CompoundArg2(g)));

				}
			}
			
			if(col && CompoundArgN(p,5))
				cpml=AppendLast(cpml,CompoundArgN(p,5));
				
			switch(IntegerValue(CompoundArgN(p,4)))
			{
				case 2:
                case 4:
				case 1:
				case 0:
					nv++;
                int aux=0;
                v=CompoundArgN(p,7);
				if(v && is_atom(v))
				{
					int i, len;
					char *av=AtomValue(v);
					len=strlen(av);
					for(i=0;i<len;i++) if(av[i]=='*') aux=1;
                }
				fprintf(f,"%s =  Particle(  ",convname(CompoundArg1(p)));
				v=GetAtomProperty(CompoundArg1(p),pdg);
				if(v==0 && !aux)
					printf("Error: PDG code is not provided for particle %s.\n",AtomValue(CompoundArg1(p)));
                if(v==0) v=NewInteger(ufo_auto_pdg++);
                
				fprintf(f,"pdg_code = %ld,\n",v?IntegerValue(v):0);
				fprintf(f,"              name = '%s',\n",
									AtomValue(CompoundArg1(p)));
				fprintf(f,"              antiname = '%s',\n",
									AtomValue(CompoundArg2(p)));
				fprintf(f,"              spin = %ld,\n",
						IntegerValue(CompoundArgN(p,4))+1);
				fprintf(f,"              color = %d,\n",col);
				v=CompoundArgN(p,5);
				fprintf(f,"              mass = '%s',\n",v?AtomValue(v):"ZERO");
				v=CompoundArgN(p,6);
				fprintf(f,"              width = '%s',\n",v?AtomValue(v):"ZERO");
				v=GetAtomProperty(CompoundArg1(p),A_TEXNAME);
				fprintf(f,"              texname = '%s',\n",
						v?AtomValue(v):AtomValue(CompoundArg1(p)));
				v=GetAtomProperty(CompoundArg2(p),A_TEXNAME);
				fprintf(f,"              antitexname = '%s',\n",
						v?AtomValue(v):AtomValue(CompoundArg2(p)));
				v=GetAtomProperty(CompoundArg1(p),A_EM_CHARGE);
				if(v==0 || CompoundArg2(v)==NewInteger(1))
				fprintf(f,"              charge = %ld,\n",
						v?IntegerValue(CompoundArg1(v)):0);
				else
				fprintf(f,"              charge = %ld/%ld,\n",
						IntegerValue(CompoundArg1(v)),
						IntegerValue(CompoundArg2(v)));
				
					if(aux)
					{
						int sp=(int)IntegerValue(CompoundArgN(p,4));
						fprintf(f,"              propagating = False,\n");
						if(sp==0) fprintf(f,"              propagator = Prop.ConstS,\n");
						else if(sp==1) fprintf(f,"              propagator = Prop.ConstF,\n");
						else if(sp==2) fprintf(f,"              propagator = Prop.ConstV,\n");
						else if(sp==4) fprintf(f,"              propagator = Prop.ConstT,\n");
						else
						{
							ErrorInfo(700);
							printf("Constant propagator in UFO undefined for spin %d/2\n",sp);
						}
					}
						
				
				fprintf(f,"              LeptonNumber = 0,\n");
				fprintf(f,"              GhostNumber = 0)\n\n");
				
				if(CompoundArg2(p)!=CompoundArg1(p))
				{
					fprintf(f,"%s = ",convname(CompoundArg2(p)));
					fprintf(f,"%s.anti()\n\n",convname(CompoundArg1(p)));
				}

					break;
			}
			continue;
		}
		if(CompoundName(p)==OPR_FIELD)
		{
			Atom bp=CompoundArg1(p);
			Term bpp=GetAtomProperty(bp,PROP_TYPE);
						
			if(CompoundArg1(bpp)!=CompoundArg2(bpp) && CompoundArg2(bpp)==bp)
				continue;

			if(CompoundArg2(p)==NewInteger(1) && CompoundArgN(bpp,5))
			{
				ns++;
				Atom aa=GetAtomProperty(a,A_ANTI);
				fprintf(f,"%s =  Particle(  ",convname(a));
				v=GetAtomProperty(a,pdg);
				fprintf(f,"pdg_code = %ld,\n",v?IntegerValue(v):0);
				fprintf(f,"              name = '%s',\n",
									AtomValue(a));
				fprintf(f,"              antiname = '%s',\n",
									AtomValue(aa));
				fprintf(f,"              spin = 1,\n");
				fprintf(f,"              color = %d,\n",col);
				v=CompoundArgN(bpp,5);
				fprintf(f,"              mass = '%s',\n",v?AtomValue(v):"ZERO");
				v=CompoundArgN(bpp,6);
				fprintf(f,"              width = 'ZERO',\n");
				v=GetAtomProperty(a,A_TEXNAME);
				fprintf(f,"              texname = '%s',\n",
						v?AtomValue(v):AtomValue(CompoundArg1(p)));
				v=GetAtomProperty(aa,A_TEXNAME);
				fprintf(f,"              antitexname = '%s',\n",
						v?AtomValue(v):AtomValue(CompoundArg2(p)));
				v=GetAtomProperty(CompoundArg1(bpp),A_EM_CHARGE);
				if(v==0 || CompoundArg2(v)==NewInteger(1))
				fprintf(f,"              charge = %ld,\n",
						v?IntegerValue(CompoundArg1(v)):0);
				else
				fprintf(f,"              charge = %ld/%ld,\n",
						IntegerValue(CompoundArg1(v)),
						IntegerValue(CompoundArg2(v)));
				fprintf(f,"              LeptonNumber = 0,\n");
				fprintf(f,"              GoldstoneBoson = True,\n");
				fprintf(f,"              GhostNumber = 0)\n\n");
				
				if(CompoundArg2(bpp)!=CompoundArg1(bpp))
				{
					fprintf(f,"%s = ",convname(aa));
					fprintf(f,"%s.anti()\n\n",convname(a));
				}
				continue;
			}
			
			
		}
	}

	fclose(f);

	for(l1=l;l1;l1=ListTail(l1))
	{
		Term a2=ListFirst(l1);
		List l2,lp;
		if(CompoundArgN(a2,5)==0 || ListLength(CompoundArg1(a2))<2)
			continue;
		if(NoColors)
		{
			alg2_rem_col(a2);
			if(CompoundArgN(a2,5)==0)
				continue;
		} 
		
		
		if(ListLength(CompoundArg1(a2))==2 && !write_all_vertices)
		{
			List ml1=ConsumeCompoundArg(a2,5);
			List ml2=NewList();
			for(l2=ml1;l2;l2=ListTail(l2))
			{
				List l3;
				for(l3=CompoundArg2(ListFirst(l2));l3;l3=ListTail(l3))
					if(GetAtomProperty(CompoundArg1(ListFirst(l3)),A_INFINITESIMAL))
						break;
				if(l3)
				{
					ml2=AppendLast(ml2,ListFirst(l2));
					ChangeList(l2,0);
				}
			}
			FreeAtomic(ml1);
			SetCompoundArg(a2,5,ml2);
			if(CompoundArgN(a2,5)==0)
				continue;
		}
		alg2_red_rc(a2);
		if(CompoundArgN(a2,5)==0) continue;
		alg2_symmetrize(a2);
		alg2_common_n(a2);
		alg2_common_s(a2);
		alg2_fix_uuv(a2);
/*		alg2_red_cos(a2);
		alg2_red_orth(a2);
		

		alg2_red_sico(a2);
		alg2_red_comsico(a2);*/

		if(CompoundArgN(a2,5)==0) continue;
		//WriteTerm(a2); puts(" --> ");
        //if(ListLength(CompoundArg1(a2))<4) alg2_red_1pm5(a2);
        //WriteTerm(a2); puts("\n -------------------");
        
		alg2_fix_ff(a2);
		alg2_multbyi(a2);
		if(opAbbrVrt)
		{
			alg2_decommon_s(a2);
			alg2_abbr_vrt(a2);
			/*alg2_common_s(a2);*/
		}
		alg2_recommon_n(a2);
		if(!opAbbrVrt)
			{
			alg2_decommon_s(a2);
			allow_sym_div=1;
			alg2_common_s(a2);
			allow_sym_div=0;
			}
		if(opEvalVrt)
			alg2_eval_vrt(a2);
		if(EvalVrt)
		{
			alg2_decommon_s(a2);
			alg2_decommon_n(a2);
			alg2_eval_vrtn(a2);
			alg2_common_n(a2);
			alg2_common_s(a2);
		}
			
	}
	
	alg2_setcls(l);

/*	cls_lagr_hook=l;
	ProcHermiticity(A_I,0);
	cls_lagr_hook=0;
*/	
	for(l1=l;l1;l1=ListTail(l1))
	{
		Term a2=ListFirst(l1);
		
		List l2,lp;
		if(CompoundArgN(a2,5)==0 || ListLength(CompoundArg1(a2))<2)
			continue;
		if(ListLength(CompoundArg1(a2))>6)
		{
			static int repno=0;
			if(repno<10)
			{
				printf("Vertex ");
				WriteVertex(CompoundArg1(a2));
				printf(" with more than 6 particles.\n");
			}
			if(repno==10)
				puts("More vertices with more than 6 particles follow");
			repno++;
			continue;
		}
		/*printf("%d\n",ListLength(CompoundArg1(a2)));*/

		/*alg2_common_t(a2);*/
		
			
/*		fWriteTerm(f,a2);
		fprintf(f,"\n");*/
		lp=NewList();
		for(l2=CompoundArg1(a2);l2;l2=ListTail(l2))
		{
			Atom p, s;
			p=CompoundArg1(ListFirst(l2));
			s=CompoundName(CompoundArg2(ListFirst(l2)));
			if(s==OPR_SCALAR && GetAtomProperty(p,A_GRASS))
				s=A_GRASS;
			lp=AppendLast(lp,s);
		}
		
		for(l2=gkl;l2;l2=ListTail(l2))
			if(EqualTerms(lp,CompoundArg1(ListFirst(l2))))
				break;
		if(l2==0)
			gkl=AppendLast(gkl,gkl1=MakeCompound2(OPR_MINUS,lp,NewList()));
		else
		{
			FreeAtomic(lp);gkl1=ListFirst(l2);
		}
		
		
		for(l2=CompoundArgN(a2,5);l2;l2=ListTail(l2))
		{
			List l3;
			Term ls=conv_lor(CompoundArg1(a2),ListFirst(l2),0);
			for(l3=CompoundArg2(gkl1);l3;l3=ListTail(l3))
				if(EqualTerms(ls,ListFirst(l3))) break;
			if(l3)
			{
				if(ls) FreeAtomic(ls);
				continue;
			}
			else
			{
				l3=ConsumeCompoundArg(gkl1,2);
				l3=AppendLast(l3,ls);
				SetCompoundArg(gkl1,2,l3);
			}
		}
		
		
		
		/*fWriteTerm(f,lp);fprintf(f,"\n");*/
	}
/*
DumpList(gkl);
*/	
	if(NoColors)
	for(l1=all_prtc_list();l1;l1=ListTail(l1))
		if(GetAtomProperty(ListFirst(l1),A_COLOR))
			SetAtomProperty(ListFirst(l1),A_COLOR,0);

	if(OutputDirectory!=NULL)
		sprintf(cbuf,"%s/vertices.py",OutputDirectory);
	else
		sprintf(cbuf,"vertices.py");
	f=fopen(cbuf,"w");
	
	fprintf(f,"#\tLanHEP output produced at ");
	time(&tm);
	fprintf(f,"%s",ctime(&tm));
	fprintf(f,"#\tfrom the file '%s'\n",eff_infile);
	if(ModelName)
		fprintf(f,"#\tModel named '%s'\n",ModelName);
	fprintf(f,"\n\n");
	
/*
	inf_decl_hc(f);
	prm_decl_hc(f,l);
*/
	for(l1=UFOverth;l1;l1=ListTail(l1))
	fprintf(f,"%s\n",AtomValue(ListFirst(l1)));
	fprintf(f,"\n");
		
	
	for(l1=l;l1;l1=ListTail(l1))
	{
		Term a2=ListFirst(l1);
		List l2,lp,lc;
		List *lv, *lv1,*lv2;
		int i,lvl,n,hasz=0;
		cv=CompoundArg1(a2);
		
		if(CompoundArgN(a2,5)==0 || ListLength(CompoundArg1(a2))<2 || 
				ListLength(CompoundArg1(a2))>6)
			continue;
		lp=NewList();
		fprintf(f,"V_%d = Vertex(name='V_%d',\n",vno,vno);vno++;
		fprintf(f,"              particles = [ ");
		for(l2=CompoundArg1(a2),i=0;l2;l2=ListTail(l2),i++)
		{
			Atom a;
			fprintf(f," P.%s",convname(CompoundArg1(ListFirst(l2))));
			if(ListTail(l2)) fprintf(f,", ");
			a=CompoundName(CompoundArg2(ListFirst(l2)));
			if(a==OPR_SCALAR && 
				GetAtomProperty(CompoundArg1(ListFirst(l2)),A_GRASS)) a=A_GRASS;
			lp=AppendLast(lp,a);
		}
		fprintf(f," ],\n");
						
		for(i=0;i<6;i++) pind[i]=wind[i]=0;
		for(l2=CompoundArg1(a2),i=0;l2;l2=ListTail(l2),i++)
		{
			Atom a=CompoundName(CompoundArg2(ListFirst(l2)));
			if(GetAtomProperty(CompoundArg1(ListFirst(l2)),A_COLOR))
				pind[i]=ListNth(CompoundArg1(CompoundArg2(ListFirst(l2))),
						(a==OPR_VECTOR||a==OPR_SPINOR)?2:1);
			else
				pind[i]=0;
			if(GetAtomProperty(CompoundArg1(ListFirst(l2)),OPR_CLASS))
			{
				List l4;
				for(l4=CompoundArg1(CompoundArg2(ListFirst(l2)));
					l4&&ListTail(l4);l4=ListTail(l4));
				wind[i]=l4?ListFirst(l4):0;
			}
			
		}

		/*WriteVertex(cv);printf(" ( ");
		for(i=0;i<4;i++){printf(" ");WriteTerm(pind[i]);}printf(" ) ( ");
		for(i=0;i<4;i++){printf(" ");WriteTerm(wind[i]);}puts(" )");*/
		

		

		/*		
		wrt_expr(f,CompoundArgN(a2,2),CompoundArgN(a2,3),CompoundArgN(a2,4));
		fprintf(f,"\n");
		*/
		if(CompoundArg1(CompoundArg2(a2))!=NewInteger(1) ||
			CompoundArg2(CompoundArg2(a2))!=NewInteger(1))
			puts("Internal error (wrtufcnf1)");
		if(CompoundArgN(a2,3) || CompoundArgN(a2,4))
			puts("Internal error (wrtufcsf1)");
		
		for(l2=gkl;l2;l2=ListTail(l2))
			if(EqualTerms(CompoundArg1(ListFirst(l2)),lp))
				break;
		if(l2==0)
		{
			WriteTerm(CompoundArg1(a2));
			WriteTerm(lp);puts("Internal error wrtfa05");
			continue;
		}
		
				
		FreeAtomic(lp);
		lp=CompoundArg2(ListFirst(l2));
		lv=calloc(lvl=ListLength(lp),sizeof(List));
		lv1=calloc(lvl,sizeof(List));
		lv2=calloc(lvl,sizeof(List));
		for(i=0;i<lvl;i++)
			lv[i]=lv1[i]=0;
			
			
		{
			int i;
			List l3;
			char cbuf[16];
			
			for(l3=CompoundArg1(ListFirst(l2)),i=0;l3;l3=ListTail(l3),i++)
			{
				Atom a=ListFirst(l3);
				if(a==OPR_SCALAR) 
					{cbuf[i]='S';}
				else if(a==OPR_SPINOR) 
					{cbuf[i]='F';}
				else if(a==OPR_VECTOR) 
					{cbuf[i]='V';}
				else if(a==OPR_TENSOR) 
					{cbuf[i]='T';}
				else if(a==A_GRASS) 
					{cbuf[i]='U';}
				else {puts("Internal error wrtuf01");cbuf[i]='?';}

			}
			n=i; cbuf[i]=0;
			fprintf(f,"              lorentz = [ ");
			for(i=1;i<=lvl;i++)
				fprintf(f,"L.%s_%d%c ",cbuf,i,i==lvl?' ':',');
			fprintf(f,"],\n");
		}

			
		lc=0;
		for(l2=CompoundArgN(a2,5);l2;l2=ListTail(l2))
		{
			List l3;
			Term pr;
			int io=0;
			Term ls=conv_lor(CompoundArg1(a2),ListFirst(l2),1);
			for(l3=lp,i=0;lp;l3=ListTail(l3),i++)
				if(EqualTerms(ls,ListFirst(l3)))
					break;
			if(l3==0 || i>=lvl)
			{puts("Internal error wrtfa06");continue;}
			for(l3=CompoundArg2(ListFirst(l2));l3;l3=ListTail(l3))
				if((pr=GetAtomProperty(CompoundArg1(ListFirst(l3)),A_INFINITESIMAL)))
					io+=(int)IntegerValue(CompoundArg1(pr))*
						(int)IntegerValue(CompoundArg2(ListFirst(l3)));
			for(l3=CompoundArgN(ListFirst(l2),3);l3;l3=ListTail(l3))
				if((is_atom(CompoundArg1(ListFirst(l3)))&& 
				  (pr=GetAtomProperty(CompoundArg1(ListFirst(l3)),A_INFINITESIMAL))))
					io+=1;
			if(io)
				puts("Error: Ren.Const. are not supported in this format.");
			if(io==1)
				lv1[i]=AppendLast(lv1[i],ListFirst(l2)),(hasz==0?hasz=1:0);
			else if(io==0)
				lv[i]=AppendLast(lv[i],ListFirst(l2));
			else
				lv2[i]=AppendLast(lv2[i],ListFirst(l2)),hasz=2;
			{
			for(l3=lc;l3;l3=ListTail(l3))
				if(EqualTerms(ListFirst(l3),CompoundArgN(ListFirst(l2),3)))
					break;
			if(!l3) lc=AppendLast(lc,CompoundArgN(ListFirst(l2),3));
			}
		}
		
		
		
		
		for(i=0;i<6;i++) pind[i]=0;
		for(l2=CompoundArg1(a2),i=0;l2;l2=ListTail(l2),i++)
		{
			Atom a=CompoundName(CompoundArg2(ListFirst(l2)));
			if(GetAtomProperty(CompoundArg1(ListFirst(l2)),A_COLOR))
				pind[i]=ListNth(CompoundArg1(CompoundArg2(ListFirst(l2))),
						a==OPR_TENSOR?3:((a==OPR_VECTOR||a==OPR_SPINOR)?2:1));
			else
				pind[i]=0;
		}

		fprintf(f,"              color = [");
		for(l2=lc;l2;l2=ListTail(l2))
		{
			fprintf(f," '");
			wrt_expr(f,NewInteger(1),0,CopyTerm(ListFirst(l2)));
			fprintf(f,"'");
			if(ListTail(l2))
				fprintf(f,",");
		}
		fprintf(f," ],\n");
		fprintf(f,"              couplings = {");

		firstprt=1;
		for(i=0;i<lvl;i++)
		for(l2=lv[i];l2;l2=ListTail(l2))
		{
			int j,no;
			List l3;
			if(firstprt)
				firstprt=0;
			else
				fprintf(f,",");
			for(l3=lc,j=0;l3;l3=ListTail(l3),j++)
				if(EqualTerms(ListFirst(l3),CompoundArgN(ListFirst(l2),3)))
					break;
			if(l3==0)
				puts("Internal error (wrtufcsk)");
			if(ListLength(CompoundArg2(ListFirst(l2)))!=1)
				puts("Internal error (wrtufmsf)");
			else
				sscanf(AtomValue(CompoundArg1(ListFirst(
						CompoundArg2(ListFirst(l2)))))+1,"%d",&no);
			fprintf(f,"(%d,%d):C.GC_%d",j,i,no);
		}
		fprintf(f,"})\n\n");
			
			

		/*
		for(i=0;i<lvl;i++)
		{
			fWriteTerm(f,lv[i]);fprintf(f,"\n");
			if(lv1[i]) lv1[i]=SortedList(lv1[i],pcmp);
			
			
			if(lv[i]==0 && lv1[i]==0 && lv2[i]==0)
			{
				if(hasz==1)
					fprintf(f," { 0, 0 }%c\n",i==lvl-1?' ':',');
				else if(hasz==2)
					fprintf(f," { 0, 0, 0 }%c\n",i==lvl-1?' ':',');
				else
					fprintf(f," { 0 }%c\n",i==lvl-1?' ':',');
				continue;
			}
			fprintf(f," { ");
			if(lv[i]==0)
				fprintf(f,"0 ");
			else
			for(l2=lv[i];l2;l2=ListTail(l2))
			{
				Term m2=ListFirst(l2);
				wrt_expr(f,CompoundArg1(m2),CompoundArg2(m2),
						CompoundArgN(m2,3));
				if(ListTail(l2) &&
						IntegerValue(CompoundArg1(ListFirst(ListTail(l2))))>0)
					fprintf(f,"+ ");
			}
			if(lv1[i])
			{
			fprintf(f,", ");
			for(l2=lv1[i];l2;l2=ListTail(l2))
			{
				Term m2=ListFirst(l2);
				wrt_expr(f,CompoundArg1(m2),CompoundArg2(m2),
						CompoundArgN(m2,3));
				if(ListTail(l2) &&
						IntegerValue(CompoundArg1(ListFirst(ListTail(l2))))>0)
					fprintf(f,"+ ");
			}
			}
			else if(hasz)
				fprintf(f,", 0");
			if(lv2[i])
			{
			fprintf(f,", ");
			for(l2=lv2[i];l2;l2=ListTail(l2))
			{
				Term m2=ListFirst(l2);
				wrt_expr(f,CompoundArg1(m2),CompoundArg2(m2),
						CompoundArgN(m2,3));
				if(ListTail(l2) &&
						IntegerValue(CompoundArg1(ListFirst(ListTail(l2))))>0)
					fprintf(f,"+ ");
			}
			}
			else if(hasz==2)
				fprintf(f,", 0");
			fprintf(f,"}%c\n",i==lvl-1?' ':',');
		}*/
		
		free(lv);free(lv1);free(lv2);
	}
	

	fclose(f);

/*
	for(l1=all_prtc_list();l1;l1=ListTail(l1))
	{
		List ll=AtomPropertiesList(ListFirst(l1));
		fWriteTerm(f,ListFirst(l1));
		fprintf(f,"\t");
		fWriteTerm(f,ll);fprintf(f,"\n");
	}
*/
}

static Atom fili(Integer ii, char *tns)
{
	int i;
	for(i=0;i<6;i++)
	{
		if(ii==pind[i]) return tns[i]?ili[i+10]:ili[i];
		if(ii==pind[i+10]) return ili[i+20];
	}
	return 0;
}


static List mkgammal(Term m2, List nc, Term curi, Term lasti, List pl, char *tens, int eff)
{
	List l,l2;
	static int gim=-11;
	if(nc==0) gim=-11;
		
		for(l=CompoundArgN(m2,3);l&&curi!=lasti;l=ListTail(l))
		{
			rpt:
			
			if(ListFirst(l)==0) continue;
			if(curi==lasti) break;
			if(CompoundName(ListFirst(l))==OPR_SPECIAL && 
					ListFirst(CompoundArg2(ListFirst(l)))==curi)
			{
				curi=ListNth(CompoundArg2(ListFirst(l)),2);
				if(CompoundArg1(ListFirst(l))==A_GAMMA)
				{
					Integer gi=ListNth(CompoundArg2(ListFirst(l)),3);
					List l2;
					Atom iii;
					if(is_atom(gi)) iii=gi;
					else iii=fili(gi,tens);
					if(iii)
					{
						nc=AppendLast(nc,MakeCompound1(A_GAMMA,iii));
						if(eff){FreeAtomic(ListFirst(l));ChangeList(l,0);}
						l=CompoundArgN(m2,3);
						goto rpt;
					}
					for(l2=CompoundArgN(m2,3);l2;l2=ListTail(l2))
						if(ListFirst(l2) &&
								CompoundName(ListFirst(l2))==A_MOMENT &&
								ListFirst(CompoundArg2(ListFirst(l2)))==gi)
					{
						nc=AppendLast(nc,MakeCompound1(A_GAMMA,
								imom[IntegerValue(CompoundArg1(ListFirst(l2)))-1]));
						if(eff)
						{
							FreeAtomic(ListFirst(l));ChangeList(l,0);
							FreeAtomic(ListFirst(l2));ChangeList(l2,0);
						}
						l=CompoundArgN(m2,3);
						goto rpt;
					}
					for(l2=CompoundArgN(m2,3);l2;l2=ListTail(l2))
					{
						if(l2==l) continue;
						if(ListFirst(l2) &&
								CompoundArg1(ListFirst(l2))==A_GAMMA &&
								ListNth(CompoundArg2(ListFirst(l2)),3)==gi)
						{
							char cbuf[80];
							Atom ni;
							sprintf(cbuf,"li%d",gim--);
							ni=NewAtom(cbuf,0);
							nc=AppendLast(nc,MakeCompound1(A_GAMMA,ni));
							ChangeList(ListNthList(CompoundArg2(ListFirst(l2)),3),ni);
							if(eff){FreeAtomic(ListFirst(l));ChangeList(l,0);}
							else ChangeList(ListNthList(CompoundArg2(ListFirst(l)),3),ni);
							l=CompoundArgN(m2,3);
							goto rpt;
						}
					}
					printf("Internal error: ");WriteVertex(pl);printf("lost gamma index:");
					WriteTerm(m2);puts("");
					if(eff){FreeAtomic(ListFirst(l));ChangeList(l,0);}
					l=CompoundArgN(m2,3);
					goto rpt;
				}
				if(CompoundArg1(ListFirst(l))!=A_DELTA)
					nc=AppendLast(nc,CompoundArg1(ListFirst(l)));
				else
				{
					if(pind[0]==ListFirst(CompoundArg2(ListFirst(l))) &&
							pind[1]==ListNth(CompoundArg2(ListFirst(l)),2))
						pind[0]=pind[1]=0;
					if(pind[2]==ListFirst(CompoundArg2(ListFirst(l))) &&
							pind[3]==ListNth(CompoundArg2(ListFirst(l)),2))
						pind[2]=pind[3]=0;
					nc=AppendLast(nc,A_DELTA);
				}
				if(eff){FreeAtomic(ListFirst(l));ChangeList(l,0);}
				l=CompoundArgN(m2,3);
				goto rpt;
			}
		}
	return nc;	
	}

static Term conv_lor(List pl, Term m2, int eff)
{
	int i;
	List l;
	List nc=0, cc=0;
	char tens[10];
  //  WriteTerm(m2);puts("");
	//WriteTerm(pl);puts("");
	if(ili[0]==0)
	{
		char cbuf[16];
		for(i=1;i<=6;i++)
		{
			sprintf(cbuf,"li%d",i); ili[i-1]=NewAtom(cbuf,0);
			sprintf(cbuf,"li%d",i+1000); ili[i-1+10]=NewAtom(cbuf,0);
			sprintf(cbuf,"li%d",i+2000); ili[i-1+20]=NewAtom(cbuf,0);
			sprintf(cbuf,"mom%d",i); imom[i-1]=NewAtom(cbuf,0);
		}
	}
    for(i=0;i<16;i++) pind[i]=0;
    
	for(l=pl,i=0;l;l=ListTail(l),i++)
	{
		Atom a;
		tens[i]=0;
		a=CompoundName(CompoundArg2(ListFirst(l)));
		if(a==OPR_VECTOR || a==OPR_SPINOR)
			pind[i]=ListFirst(CompoundArg1(CompoundArg2(ListFirst(l))));
		else if(a==OPR_TENSOR)
			pind[i]=ListFirst(CompoundArg1(CompoundArg2(ListFirst(l)))),
			pind[i+10]=ListFirst(ListTail(CompoundArg1(CompoundArg2(ListFirst(l))))),
			tens[i]=1;
		else
			pind[i]=0; 
	}
	
	
	if(CompoundName(CompoundArg2(ListFirst(pl)))==OPR_SPINOR)
		nc=mkgammal(m2,nc,pind[0],pind[1],pl,tens,eff);
	if(ListLength(pl)==4 && CompoundName(CompoundArg2(ListNth(pl,3)))==OPR_SPINOR)
	{nc=AppendLast(nc,0); nc=mkgammal(m2,nc,pind[2],pind[3],pl,tens,eff);}

		
	for(l=CompoundArgN(m2,3);l;l=ListTail(l))
	{
		if(ListFirst(l)==0)
			continue;
		if(CompoundArg1(ListFirst(l))==A_DELTA)
		{
			Atom i1,i2;
			Integer ci=ListFirst(CompoundArg2(ListFirst(l)));
			for(i=0;i<6;i++)
				if(ci==pind[i]||ci==pind[i+10])
					break;
			if(i==6)
				continue;
			if(tens[i]==0) i1=ili[i];
			else if(ci==pind[i]) i1=ili[i+10];
			else i1=ili[i+20];
			ci=ListFirst(ListTail(CompoundArg2(ListFirst(l))));
			for(;i<6;i++)
				if(ci==pind[i]||ci==pind[i+10])
					break;
			if(i==6)
			{
			printf("Internal error: ");WriteVertex(pl);printf("lost delta index:");
			WriteTerm(m2);
			puts("");
			continue;
			}
			if(tens[i]==0) i2=ili[i];
			else if(ci==pind[i]) i2=ili[i+10];
			else i2=ili[i+20];
			if(eff){FreeAtomic(ListFirst(l));ChangeList(l,0);}
			cc=AppendLast(cc,MakeCompound2(A_DELTA,i1,i2));
		}
	}
	
	for(l=CompoundArgN(m2,3);l;l=ListTail(l))
	{
		if(ListFirst(l)==0)
			continue;
		if(CompoundArg1(ListFirst(l))==A_EPS_V)
		{
			Atom in[6];
			int j;
			List lj,l2;
			for(j=0,lj=CompoundArg2(ListFirst(l));j<4;j++,lj=ListTail(lj))
			{
			in[j]=0;
			Atom iii=fili(ListFirst(lj),tens);
			if(iii)
				{
				in[j]=iii;
				continue;
				}
			for(l2=CompoundArgN(m2,3);l2;l2=ListTail(l2))
				if(CompoundName(ListFirst(l2))==A_MOMENT && ListFirst(lj)==
					ListFirst(CompoundArg2(ListFirst(l2))))
						break;
			if(l2)
				in[j]=imom[IntegerValue(CompoundArg1(ListFirst(l2)))-1];
			
			if(l2&& eff){FreeAtomic(ListFirst(l2));ChangeList(l2,0);}
			}
			if(eff){FreeAtomic(ListFirst(l));ChangeList(l,0);}
			lj=MakeCompound(A_EPS_V,4);
			SetCompoundArg(lj,1,in[0]);
			SetCompoundArg(lj,2,in[1]);
			SetCompoundArg(lj,3,in[2]);
			SetCompoundArg(lj,4,in[3]);
			cc=AppendLast(cc,lj);
		}
	}
	
	for(l=CompoundArgN(m2,3);l;l=ListTail(l))
	{
		if(ListFirst(l)==0)
			continue;
		if(CompoundName(ListFirst(l))==A_MOMENT)
		{
			int p1=(int)IntegerValue(CompoundArg1(ListFirst(l)));
			Term pi=ListFirst(CompoundArg2(ListFirst(l)));
			List l2;
			Atom iii=fili(pi,tens);
			if(iii)
			{
				cc=AppendLast(cc,MakeCompound2(A_DELTA,iii,imom[p1-1]));
				if(eff){FreeAtomic(ListFirst(l));ChangeList(l,0);}
				continue;
			}
			for(l2=ListTail(l);l2;l2=ListTail(l2))
			{
				if(ListFirst(l2)==0) continue;
				if(CompoundName(ListFirst(l2))==A_MOMENT && 
						pi==ListFirst(CompoundArg2(ListFirst(l2))))
				{
					int p2=(int)IntegerValue(CompoundArg1(ListFirst(l2)));
					if(p1<p2)
						cc=AppendLast(cc,MakeCompound2(A_DELTA,imom[p1-1],imom[p2-1]));
					else
						cc=AppendLast(cc,MakeCompound2(A_DELTA,imom[p2-1],imom[p1-1]));
					if(eff)
					{
						FreeAtomic(ListFirst(l));ChangeList(l,0);
						FreeAtomic(ListFirst(l2));ChangeList(l2,0);
					}
					break;
				}
			}
		}
	}
	

	if(eff)
	{
		List l2;
		l=ConsumeCompoundArg(m2,3);
		rpt2:
		for(l2=l;l2;l2=ListTail(l2))
			if(ListFirst(l2)==0)
			{
				l=CutFromList(l,l2);
				goto rpt2;
			}
		SetCompoundArg(m2,3,l);
	}
 //   WriteTerm(cc); printf(" "); WriteTerm(nc);puts("\n");
	return MakeCompound2(OPR_PLUS,cc,nc);
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
			sno+=fprintf(of,"%d",n),f=0;
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
	for(l=ten;l;l=ListTail(l))
	{
		
		if(CompoundArg1(ListFirst(l))==A_DELTA)
		{
			Integer in1,in2;
			int i,il1, il2;
			in1=ListFirst(CompoundArg2(ListFirst(l)));
			in2=ListFirst(ListTail(CompoundArg2(ListFirst(l))));
			for(i=0;i<6;i++) if(in1==pind[i]) break;
			if(i==6)
			{
				for(i=0;i<6;i++) if(in1==wind[i]) break;
				il1=i;
				if(i==6) puts("Internal error wrtfa07");
				for(i=0;i<6;i++) if(in2==wind[i]) break;
				il2=i;
				if(i==6) puts("Internal error wrtfa07");
				sno+=fprintf(of,"Identity(t%d, t%d)",il1+1,il2+1);
				if(sno>60 && ListTail(l)) {fprintf(of,"*\n\t\t");sno=15;}
				continue;
			}
			il1=i;
			for(i=0;i<6;i++) if(in2==pind[i]) break;
			if(i==6){puts("Internal error wrtfa07");}
			il2=i;
			
			sno+=fprintf(of,"Identity(%d,%d)",il1+1,il2+1);
			if(sno>60 && ListTail(l)) {fprintf(of,"*\n\t\t");sno=15;}
			continue;
		}
		
		if(CompoundName(ListFirst(l))==OPR_PARAMETER)
		{
			Integer in1;
			int i;
			List l1;
			sno+=fprintf(of,"%s[",AtomValue(CompoundArg1(ListFirst(l))));
			for(l1=CompoundArg2(ListFirst(l));l1;l1=ListTail(l1))
			{
				in1=ListFirst(l1);
				for(i=0;i<6;i++) if(in1==wind[i]) break;
				if(i==6) puts("Internal error wrtfa07");
				sno+=fprintf(of,"t%d",i+1);
				if(ListTail(l1)) sno+=fprintf(of,", ");
			}
			sno+=fprintf(of,"] ");
			if(sno>60 && ListTail(l)) {fprintf(of,"*\n\t\t");sno=15;}
			continue;
		}
		if(CompoundName(ListFirst(l))!=OPR_SPECIAL)
		{
			puts("Internal error wrtfaus");WriteTerm(ListFirst(l));puts("");
			fWriteTerm(of,ListFirst(l));
			continue;
		}
		
		if(GetAtomProperty(CompoundArg1(ListFirst(l)),A_COLOR)
			==A_COLOR_LAMBDA)
		{
			Integer in1,in2,in3,in4;
			int i,il1, il2, il3, il4;
			List l2;
			
			in1=ListFirst(CompoundArg2(ListFirst(l)));
			in2=ListFirst(ListTail(CompoundArg2(ListFirst(l))));
			in3=ListFirst(ListTail(ListTail(CompoundArg2(ListFirst(l)))));
			
			if(in1==0 || in2==0 || in3==0) continue;
			
			for(i=0;i<6;i++) if(in1==pind[i]) break;
			il1=(i==6?(int)IntegerValue(in1)+4:i);
			
			for(i=0;i<6;i++) if(in2==pind[i]) break;
			il2=(i==6?(int)IntegerValue(in2)+4:i);
			
			for(i=0;i<6;i++) if(in3==pind[i]) break; il3=i;
			il3=(i==6?(int)IntegerValue(in3)+4:i);
						
			if(il1<5 && il2<5 && il3<5)
			{
			sno+=fprintf(of,"T(%d,%d,%d)",il3+1,il2+1,il1+1);
			if(sno>60 && ListTail(l)) {fprintf(of,"*\n\t\t");sno=15;}
			continue;
			}
			for(l2=ListTail(l);l2;l2=ListTail(l2))
				if(GetAtomProperty(CompoundArg1(ListFirst(l2)),A_COLOR)
					==A_COLOR_LAMBDA )
					break;
			
			if(il2>=5 && l2 && ListFirst(CompoundArg2(ListFirst(l2)))==in2)
			{
			in3=ListFirst(ListTail(CompoundArg2(ListFirst(l2))));
			in4=ListFirst(ListTail(ListTail(CompoundArg2(ListFirst(l2)))));
			for(i=0;i<6;i++) if(in3==pind[i]) break; il2=i;
			if(i==6){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(2)");il3=(int)IntegerValue(in3);}
			for(i=0;i<6;i++) if(in4==pind[i]) break; il4=i;
			if(i==6){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(2)");il4=(int)IntegerValue(in4);}
			sno+=fprintf(of,"T(%d,%d,-1)*T(%d,-1,%d)",
					il3+1,il2+1,il4+1,il1+1);
			if(sno>60 && ListTail(l) && (ListTail(l)!=l2||ListTail(ListTail(l))))
				 {fprintf(of,"*\n\t\t");sno=15;}
			l2=CompoundArg2(ListFirst(l2));
			ChangeList(l2,0);
			continue;
			}
			if(il3>=5 && l2 && ListNth(CompoundArg2(ListFirst(l2)),3)==in3)
			{
			in3=ListFirst(CompoundArg2(ListFirst(l2)));
			in4=ListFirst(ListTail(CompoundArg2(ListFirst(l2))));
			for(i=0;i<6;i++) if(in3==pind[i]) break; il3=i;
			if(i==6){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(3)");il3=(int)IntegerValue(in3);}
			for(i=0;i<6;i++) if(in4==pind[i]) break; il4=i;
			if(i==6){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(3)");il4=(int)IntegerValue(in4);}
			sno+=fprintf(of,"T(-1,%d,%d)*T(-1,%d,%d)",
					il2+1,il1+1,il4+1,il3+1);
			if(sno>60 && ListTail(l) && (ListTail(l)!=l2||ListTail(ListTail(l)))) 
				{fprintf(of,"*\n\t\t");sno=15;}
			l2=ListTail(ListTail(CompoundArg2(ListFirst(l2))));
			ChangeList(l2,0);
			continue;
			}
			
			WriteVertex(cv);WriteTerm(ListFirst(l));
			puts(": color structure error");il2=(int)IntegerValue(in2);
			sno+=fprintf(of,"SUNT[c%d, c%d, c%d] ",il3+1,il1+1,il2+1);
			if(sno>60 && ListTail(l)) {fprintf(of,"*\n\t\t");sno=15;}
			continue;
		}
		if(GetAtomProperty(CompoundArg1(ListFirst(l)),A_COLOR)
			==A_COLOR_F)
		{
			Integer in1,in2,in3,in4;
			int i,il1, il2, il3,il4;
			List l2;
			in1=ListFirst(CompoundArg2(ListFirst(l)));
			in2=ListFirst(ListTail(CompoundArg2(ListFirst(l))));
			in3=ListFirst(ListTail(ListTail(CompoundArg2(ListFirst(l)))));
			
			/*WriteTerm(cv);puts("");
			WriteTerm(ListFirst(l));
			printf("; inds: ");WriteTerm(in1);printf(" ");WriteTerm(in2);
			printf(" ");WriteTerm(in3);printf(" pind[]: ");
			WriteTerm(pind[0]);printf(" ");WriteTerm(pind[1]);printf(" ");
			WriteTerm(pind[2]);printf(" ");WriteTerm(pind[3]);printf("\n");*/
			
			if(in3==0) continue;
			for(i=0;i<6;i++) if(in1==pind[i]) break; il1=i;
			if(i==6){WriteVertex(cv);WriteTerm(ListFirst(l));
			puts(": color structure error(1)");il1=(int)IntegerValue(in1);}
			for(i=0;i<6;i++) if(in2==pind[i]) break; il2=i;
			if(i==6){WriteVertex(cv);WriteTerm(ListFirst(l));
			puts(": color structure error(1)");il2=(int)IntegerValue(in2);}
			for(i=0;i<6;i++) if(in3==pind[i]) break; il3=i;
			if(i<6)
			{
				sno+=fprintf(of,"f(%d,%d,%d)",il1+1,il2+1,il3+1);
				if(sno>60 && ListTail(l)) {fprintf(of,"*\n\t\t");sno=15;}
				continue;
			}
			for(l2=ListTail(l);l2;l2=ListTail(l2))
				if(GetAtomProperty(CompoundArg1(ListFirst(l2)),A_COLOR)
				==A_COLOR_F && ListNth(CompoundArg2(ListFirst(l2)),3)==in3)
					break;
			if(l2)
			{
			in3=ListFirst(CompoundArg2(ListFirst(l2)));
			in4=ListFirst(ListTail(CompoundArg2(ListFirst(l2))));
			for(i=0;i<6;i++) if(in3==pind[i]) break; il3=i;
			if(i==6){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(2)");il3=(int)IntegerValue(in3);}
			for(i=0;i<6;i++) if(in4==pind[i]) break; il4=i;
			if(i==6){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(2)");il4=(int)IntegerValue(in4);}
			sno+=fprintf(of,"f(%d,%d,-1)*f(%d,%d,-1)",
					il1+1,il2+1,il3+1,il4+1);
			if(sno>60 && ListTail(l)) {fprintf(of,"*\n\t\t");sno=15;}
			l2=ListTail(ListTail(CompoundArg2(ListFirst(l2))));
			ChangeList(l2,0);
			continue;
			}
			for(l2=ListTail(l);l2;l2=ListTail(l2))
				if(GetAtomProperty(CompoundArg1(ListFirst(l2)),A_COLOR)
				==A_COLOR_D && ListNth(CompoundArg2(ListFirst(l2)),3)==in3)
					break;
			if(l2)
			{
			in3=ListFirst(CompoundArg2(ListFirst(l2)));
			in4=ListFirst(ListTail(CompoundArg2(ListFirst(l2))));
			for(i=0;i<6;i++) if(in3==pind[i]) break; il3=i;
			if(i==6){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(2)");il3=(int)IntegerValue(in3);}
			for(i=0;i<6;i++) if(in4==pind[i]) break; il4=i;
			if(i==6){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(2)");il4=(int)IntegerValue(in4);}
			sno+=fprintf(of,"f(%d,%d,-1)*d(%d,%d,-1)",
					il1+1,il2+1,il3+1,il4+1);
			if(sno>60 && ListTail(l)) {fprintf(of,"*\n\t\t");sno=15;}
			l2=ListTail(ListTail(CompoundArg2(ListFirst(l2))));
			ChangeList(l2,0);
			continue;
			}
			{WriteVertex(cv);WriteTerm(ListFirst(l));
			puts(": color structure error(3)");il3=(int)IntegerValue(in3);}
			sno+=fprintf(of,"SUNT[c%d, c%d, c%d] ",il1+1,il2+1,il3+1);
			if(sno>60 && ListTail(l)) {fprintf(of,"*\n\t\t");sno=15;}
			continue;
		}
		
		if(GetAtomProperty(CompoundArg1(ListFirst(l)),A_COLOR)
			==A_COLOR_D)
		{
			Integer in1,in2,in3,in4;
			int i,il1, il2, il3,il4;
			List l2;
			in1=ListFirst(CompoundArg2(ListFirst(l)));
			in2=ListFirst(ListTail(CompoundArg2(ListFirst(l))));
			in3=ListFirst(ListTail(ListTail(CompoundArg2(ListFirst(l)))));
			
			/*WriteTerm(cv);puts("");
			WriteTerm(ListFirst(l));
			printf("; inds: ");WriteTerm(in1);printf(" ");WriteTerm(in2);
			printf(" ");WriteTerm(in3);printf(" pind[]: ");
			WriteTerm(pind[0]);printf(" ");WriteTerm(pind[1]);printf(" ");
			WriteTerm(pind[2]);printf(" ");WriteTerm(pind[3]);printf("\n");*/
			
			if(in3==0) continue;
			for(i=0;i<6;i++) if(in1==pind[i]) break; il1=i;
			if(i==6){WriteVertex(cv);WriteTerm(ListFirst(l));
			puts(": color structure error(1)");il1=(int)IntegerValue(in1);}
			for(i=0;i<6;i++) if(in2==pind[i]) break; il2=i;
			if(i==6){WriteVertex(cv);WriteTerm(ListFirst(l));
			puts(": color structure error(1)");il2=(int)IntegerValue(in2);}
			for(i=0;i<6;i++) if(in3==pind[i]) break; il3=i;
			if(i<6)
			{
				sno+=fprintf(of,"d(%d,%d,%d)",il1+1,il2+1,il3+1);
				if(sno>60 && ListTail(l)) {fprintf(of,"*\n\t\t");sno=15;}
				continue;
			}
			for(l2=ListTail(l);l2;l2=ListTail(l2))
				if(GetAtomProperty(CompoundArg1(ListFirst(l2)),A_COLOR)
				==A_COLOR_F && ListNth(CompoundArg2(ListFirst(l2)),3)==in3)
					break;
			if(l2)
			{
			in3=ListFirst(CompoundArg2(ListFirst(l2)));
			in4=ListFirst(ListTail(CompoundArg2(ListFirst(l2))));
			for(i=0;i<6;i++) if(in3==pind[i]) break; il3=i;
			if(i==6){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(2)");il3=(int)IntegerValue(in3);}
			for(i=0;i<6;i++) if(in4==pind[i]) break; il4=i;
			if(i==6){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(2)");il4=(int)IntegerValue(in4);}
			sno+=fprintf(of,"d(%d,%d,-1)*f(%d,%d,-1)",il1+1,il2+1,il3+1,il4+1);
			if(sno>60 && ListTail(l)) {fprintf(of,"*\n\t\t");sno=15;}
			l2=ListTail(ListTail(CompoundArg2(ListFirst(l2))));
			ChangeList(l2,0);
			continue;
			}
			for(l2=ListTail(l);l2;l2=ListTail(l2))
				if(GetAtomProperty(CompoundArg1(ListFirst(l2)),A_COLOR)
				==A_COLOR_D && ListNth(CompoundArg2(ListFirst(l2)),3)==in3)
					break;
			if(l2)
			{
			in3=ListFirst(CompoundArg2(ListFirst(l2)));
			in4=ListFirst(ListTail(CompoundArg2(ListFirst(l2))));
			for(i=0;i<6;i++) if(in3==pind[i]) break; il3=i;
			if(i==6){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(2)");il3=(int)IntegerValue(in3);}
			for(i=0;i<6;i++) if(in4==pind[i]) break; il4=i;
			if(i==6){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(2)");il4=(int)IntegerValue(in4);}
			sno+=fprintf(of,"d(%d,%d,-1)*d(%d,%d,-1)",il1+1,il2+1,il3+1,il4+1);
			if(sno>60 && ListTail(l)) {fprintf(of,"*\n\t\t");sno=15;}
			l2=ListTail(ListTail(CompoundArg2(ListFirst(l2))));
			ChangeList(l2,0);
			continue;
			}
			{WriteVertex(cv);WriteTerm(ListFirst(l));
			puts(": color structure error(3)");il3=(int)IntegerValue(in3);}
			sno+=fprintf(of,"SUND[c%d, c%d, c%d] ",il1+1,il2+1,il3+1);
			if(sno>60 && ListTail(l)) {fprintf(of,"*\n\t\t");sno=15;}
			continue;
		}
		
		fWriteTerm(of,ListFirst(l));
		
	}

	if(sno!=32) f=0;
	if(!is_integer(num) && !f)
		fprintf(of,"*");
	NoQuotes=0;
}

extern List opFAGS, opFAGE, UFOloreh;

void UF_write_gen(void)
{
	List l1;
	time_t tm;
	FILE *f;
	char cbuf[256];
	
	if(OutputDirectory!=NULL)
		sprintf(cbuf,"%s/lorentz.py",OutputDirectory);
	else
		sprintf(cbuf,"lorentz.py");
	f=fopen(cbuf,"w");
	
	fprintf(f,"#\tLanHEP output produced at ");
	time(&tm);
	fprintf(f,"%s",ctime(&tm));
	fprintf(f,"#\tfrom the file '%s'\n",eff_infile);
	if(ModelName)
		fprintf(f,"#\tModel named '%s'\n",ModelName);
	fprintf(f,"\n\n");
	
	
	for(l1=UFOloreh;l1;l1=ListTail(l1))
		fprintf(f,"%s\n",AtomValue(ListFirst(l1)));
	fprintf(f,"\n");

		
	for(l1=gkl;l1;l1=ListTail(l1))
	{
		char cbuf[16];
		int i,j=1,n;
		List l2;
		List spl=NewList();
		
		for(l2=CompoundArg1(ListFirst(l1)),i=0;l2;l2=ListTail(l2),i++)
		{
			Atom a=ListFirst(l2);
			if(a==OPR_SCALAR) 
				{cbuf[i]='S';spl=AppendLast(spl,NewInteger(1));}
			else if(a==OPR_SPINOR) 
				{cbuf[i]='F';spl=AppendLast(spl,NewInteger(2));}
			else if(a==OPR_VECTOR) 
				{cbuf[i]='V';spl=AppendLast(spl,NewInteger(3));}
			else if(a==OPR_TENSOR) 
				{cbuf[i]='T';spl=AppendLast(spl,NewInteger(5));}
			else if(a==A_GRASS) 
				{cbuf[i]='U';spl=AppendLast(spl,NewInteger(-1));}
			else {puts("Internal error wrtuf01");cbuf[i]='?';}
			
		}
		n=i; cbuf[i]=0;
		
		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
		{
			List cc=CompoundArg1(ListFirst(l2)),nc=CompoundArg2(ListFirst(l2));
			List l3;
			int lin=1;
			fprintf(f,"%s_%d = Lorentz(name = '%s_%d',\n",cbuf,j,cbuf,j);
			fprintf(f,"                  spins = ");
			fWriteTerm(f,spl);
			fprintf(f,",\n");
			fprintf(f,"                  structure = '");
			if(cc==0 && nc==0)
			{
				fprintf(f,"1')\n\n");
				continue;
			}
			for(l3=cc;l3;l3=ListTail(l3))
			{
				Atom a1=CompoundArg1(ListFirst(l3)),
					 a2=CompoundArg2(ListFirst(l3));
				if(CompoundArity(ListFirst(l3))==4)
				{
					List lp=0,l6;
					int l4;
					fprintf(f," Epsilon(");
					for(l4=1;l4<=4;l4++)
					{
						if(AtomValue(CompoundArgN(ListFirst(l3),l4))[0]=='l')
							fprintf(f,"%s",AtomValue(CompoundArgN(ListFirst(l3),l4))+2);
						else
							{
							int l5;
							fprintf(f,"%d",-lin);
							sscanf(AtomValue(CompoundArgN(ListFirst(l3),l4))+3,"%d",&l5);
							lp=AppendLast(lp,MakeCompound2(OPR_DIV,NewInteger(-lin),NewInteger(l5)));
							lin++;
							}
						if(l4<4)
							fprintf(f,",");
					}
					fprintf(f,")");
					for(l6=lp;l6;l6=ListTail(l6))
						fprintf(f,"*P(%ld,%ld)",IntegerValue(CompoundArg1(ListFirst(l6))),
							IntegerValue(CompoundArg2(ListFirst(l6))));
					if(lp)
						FreeAtomic(lp);
				}		
				else
				if(AtomValue(a1)[0]=='l' && AtomValue(a2)[0]=='l')
					fprintf(f," Metric(%s,%s)",AtomValue(a1)+2,AtomValue(a2)+2);
				else if(AtomValue(a1)[0]=='l' && AtomValue(a2)[0]=='m')
					fprintf(f," P(%s,%s)",AtomValue(a1)+2,AtomValue(a2)+3);
				else if(AtomValue(a1)[0]=='m' && AtomValue(a2)[0]=='m')
					{fprintf(f," P(%d,%s)*P(%d,%s)",-lin,AtomValue(a1)+3,-lin,AtomValue(a2)+3);lin++;}
				else {printf("Internal error wrtfa02; ");WriteTerm(ListFirst(l3));puts("");}
				if(ListTail(l3) || nc) fprintf(f,"*");
			}
			
				
			if(nc)
			{
				int su=0;
				int twoli=0;
				for(l3=nc;l3;l3=ListTail(l3))
				{
                    int addin=0;
					if(ListFirst(l3)==0)
					{ twoli=2; su=0; continue;}
					if(ListFirst(l3)==A_GAMMAP) fprintf(f,"ProjP(");
					else if(ListFirst(l3)==A_GAMMAM) fprintf(f,"ProjM(");
					else if(ListFirst(l3)==A_GAMMA5) fprintf(f,"Gamma5(");
					else if(ListFirst(l3)==A_DELTA) fprintf(f,"Identity(");
					else if(is_compound(ListFirst(l3)) &&
						CompoundName(ListFirst(l3))==A_GAMMA)
					{
						Atom a = CompoundArg1(ListFirst(l3));
						if(AtomValue(a)[0]=='l')
							fprintf(f,"Gamma(%s,",AtomValue(a)+2);
						else if(AtomValue(a)[0]=='m')
                            {fprintf(f,"P(%d,%s)*Gamma(%d,",-lin-su,AtomValue(a)+3,-lin-su);lin++;addin=1;}
						else {puts("Internal error wrtuf03");}
					}
					else {puts("Internal error wrtuf04");}
					if(su==0)
						fprintf(f,"%d,",1+twoli);
					else
                    {fprintf(f,"%d,",-lin+addin);lin++;addin=0;}
					if(ListTail(l3) && ListFirst(ListTail(l3)))
					{
						fprintf(f,"%d)",-lin); su=1;
					}
					else
						fprintf(f,"%d)",2+twoli);
						
					if(ListTail(l3)) fprintf(f,"*");
				}
			}
			
			fprintf(f,"')\n\n");
			j++;
		}

										
	}
	
}
