#include <string.h>
#include "lanhep.h"
#include <time.h>

extern char *eff_infile;
extern char *ModelName;

extern int opEvalVrt, EvalVrt, NoColors, opAbbrVrt, allow_sym_div, write_all_vertices;
extern int FAver;

static char nv=0, nf=0, ns=0, nu=0;

static List gkl=0;
static Atom ili[4]={0,0,0,0}, imom[4];
static Integer pind[4], wind[4];

static Term conv_lor(List pl, Term m2, int eff);
void alg2_rem_col(Term a2);

static void wrt_expr(FILE *, Term n, List s, List t);
static List cv;

static List cpml=0;
List fainclude=0;

void inf_write_rc(FILE *), cls_write_dmatr(FILE *);
void alg2_fix_uuv(Term), alg2_fix_ff(Term), alg2_fix_mom2(Term), alg2_red_rc(Term);
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
	
static int mlsort(Term m1, Term m2)
{
	List l1,l2;
	l1=CompoundArg2(m1);
	l2=CompoundArg2(m2);
	if(ListLength(l1)!=ListLength(l2))
		return ListLength(l1)-ListLength(l2);
	for(;l1;l1=ListTail(l1),l2=ListTail(l2))
	{
		if(CompoundArg1(ListFirst(l1))!=CompoundArg1(ListFirst(l2)))
			return strcmp(AtomValue(CompoundArg1(ListFirst(l1))),
					AtomValue(CompoundArg1(ListFirst(l2))));
		if(CompoundArg2(ListFirst(l1))!=CompoundArg2(ListFirst(l2)))
			return (int)IntegerValue(CompoundArg2(ListFirst(l1)))-
					(int)IntegerValue(CompoundArg2(ListFirst(l2)));
	}
	l1=CompoundArgN(m1,3);
	l2=CompoundArgN(m2,3);
	if(ListLength(l1)!=ListLength(l2))
		return ListLength(l1)-ListLength(l2);
	for(;l1;l1=ListTail(l1),l2=ListTail(l2))
	{
		List l3,l4;
		if(CompoundName(ListFirst(l1))!=CompoundName(ListFirst(l2)))
			return (int)(CompoundName(ListFirst(l1))-CompoundName(ListFirst(l2)));
		if(CompoundArg1(ListFirst(l1))!=CompoundArg1(ListFirst(l2)))
			return (int)(CompoundArg1(ListFirst(l1))-CompoundArg1(ListFirst(l2)));
		l3=CompoundArg2(ListFirst(l1));
		l4=CompoundArg2(ListFirst(l2));
		for(;l3;l3=ListTail(l3),l4=ListTail(l4))
			if(ListFirst(l3)!=ListFirst(l4))
				return ListFirst(l3)!=ListFirst(l4);
	}
	WriteTerm(l1);WriteTerm(l2);puts("");
	return 0;
}

void FA_write_lagr(List l, FILE *f)
{
	time_t tm;
	Term gkl1;
	List l1;
	int firstprt=1;
	
	fprintf(f,"(*\n\tLanHEP output produced at ");
	time(&tm);
	fprintf(f,"%s",ctime(&tm));
	fprintf(f,"\tfrom the file '%s'\n",eff_infile);
	if(ModelName)
		fprintf(f,"\tModel named '%s'\n",ModelName);
	fprintf(f,"*)\n\n");
	
	
	if(!NoColors)
	{
	fprintf(f,"\nIndexRange[ Index[Colour] ] = NoUnfold[Range[3]]\n");
	fprintf(f,"IndexRange[ Index[Gluon] ] = NoUnfold[Range[8]]\n");
	}

	if(FAver==4)
		fprintf(f,"\nVSESign := 1\n");
	else
		fprintf(f,"\nVSESign := -1\n");
	
	cls_wrt_ind(f);
	
	fprintf(f,"\n\t\t(* Model particles  *)\n\n");
	fprintf(f,"M$ClassesDescription = {\n\n");
	
	for(l1=all_prtc_list();l1;l1=ListTail(l1))
	{
		Term p,col,tnm;
		int dim=0;
		char indsc[60];
		Atom a=ListFirst(l1);
		p=GetAtomProperty(a,PROP_TYPE);
		if(p==0) continue;
		dim=cls_prt_info(&p, &a);
		if(p==0) continue;
		if(!is_compound(p)) continue;
		col=GetAtomProperty(a,A_COLOR);
		if(col) col=IntegerValue(CompoundArg1(col));
		if(dim==0 && (col==0||NoColors))
			sprintf(indsc,"\tIndices -> {},\n");
		else if(dim && (col==0||NoColors))
			sprintf(indsc,"\tIndices -> {Index[%s]},\n",AtomValue(a));
		else if(dim==0 && !(col==0||NoColors))
			sprintf(indsc,"\tIndices -> {Index[%s]},\n",
						col==1?"Colour":"Gluon");
		else 
			sprintf(indsc,"\tIndices -> {Index[%s], Index[%s]},\n",
						AtomValue(a),col==1?"Colour":"Gluon");
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
				if(!firstprt) fprintf(f,",\n\n"); else firstprt=0;
				fprintf(f,"  U[%d] ==  { (* ",nu);
				fprintf(f,"%s/%s",AtomValue(CompoundArg1(g)),AtomValue(CompoundArg2(ag)));
				fprintf(f," *)\n\tSelfConjugate -> %s,\n","False");
				fprintf(f,"%s",indsc);
				if(col && NoColors)
					fprintf(f,"\tMatrixTraceFactor -> %d,\n",col==1?3:8);
				fprintf(f,"\tMass -> %s,\n",CompoundArgN(p,5)?
						AtomValue(CompoundArgN(p,5)):"0");
					if(defined_em_charge())
						{
						Term ch=GetAtomProperty(CompoundArg1(p),A_EM_CHARGE);
						fprintf(f,"\tCharge -> ");
						if(ch==0)
							fprintf(f,"0,\n");
						else if(CompoundArg2(ch)==NewInteger(1))
							fprintf(f,"%ld,\n",IntegerValue(CompoundArg1(ch)));
						else
							fprintf(f,"%ld/%ld,\n",IntegerValue(CompoundArg1(ch)),
											IntegerValue(CompoundArg2(ch)));
                                                if(ch)
                                                {
                                                    fprintf(f,"\tQuantumNumbers -> {");
                                                    if(CompoundArg2(ch)==NewInteger(1))
							fprintf(f," %ld Charge },\n",IntegerValue(CompoundArg1(ch)));
                                                    else
							fprintf(f," %ld/%ld Charge },\n",IntegerValue(CompoundArg1(ch)),
											IntegerValue(CompoundArg2(ch)));
                                                }
						}
				fprintf(f,"\tPropagatorLabel -> \"%s\",\n",AtomValue(CompoundArg1(g)));
				fprintf(f,"\tPropagatorType -> GhostDash,\n");
				fprintf(f,"\tPropagatorArrow -> %s }","Forward");
				SetAtomProperty(CompoundArg1(g),A_FANUM,
						MakeList2(NewInteger(3),NewInteger(nu)));
				SetAtomProperty(CompoundArg2(ag),A_FANUM,
						MakeList2(NewInteger(3),NewInteger(-nu)));
				if(CompoundArg1(p)!=CompoundArg2(p))
				{
				nu++;
				if(!firstprt) fprintf(f,",\n\n"); else firstprt=0;
				fprintf(f,"  U[%d] ==  { (* ",nu);
				fprintf(f,"%s/%s",AtomValue(CompoundArg1(ag)),AtomValue(CompoundArg2(g)));
				fprintf(f," *)\n\tSelfConjugate -> %s,\n","False");
				fprintf(f,"%s",indsc);
				if(col && NoColors)
					fprintf(f,"\tMatrixTraceFactor -> %d,\n",col==1?3:8);
				fprintf(f,"\tMass -> %s,\n",CompoundArgN(p,5)?
						AtomValue(CompoundArgN(p,5)):"0");
					if(defined_em_charge())
						{
						Term ch=GetAtomProperty(CompoundArg2(p),A_EM_CHARGE);
						fprintf(f,"\tCharge -> ");
						if(ch==0)
							fprintf(f,"0,\n");
						else if(CompoundArg2(ch)==NewInteger(1))
							fprintf(f,"%ld,\n",IntegerValue(CompoundArg1(ch)));
						else
							fprintf(f,"%ld/%ld,\n",IntegerValue(CompoundArg1(ch)),
											IntegerValue(CompoundArg2(ch)));
						}
				fprintf(f,"\tPropagatorLabel -> \"%s\",\n",AtomValue(CompoundArg1(g)));
				fprintf(f,"\tPropagatorType -> GhostDash,\n");
				fprintf(f,"\tPropagatorArrow -> %s }","Forward");
				SetAtomProperty(CompoundArg1(ag),A_FANUM,
						MakeList2(NewInteger(3),NewInteger(nu)));
				SetAtomProperty(CompoundArg2(g),A_FANUM,
						MakeList2(NewInteger(3),NewInteger(-nu)));		
				}
			}
			
			if(col && CompoundArgN(p,5))
				cpml=AppendLast(cpml,CompoundArgN(p,5));
				
			switch(IntegerValue(CompoundArgN(p,4)))
			{
				case 2:
					nv++;
					if(!firstprt) fprintf(f,",\n\n"); else firstprt=0;
					fprintf(f,"  V[%d] ==  { (* ",nv);
					fWriteTerm(f,CompoundArgN(p,3));
					fprintf(f," *)\n\tSelfConjugate -> %s,\n",
							CompoundArg1(p)==CompoundArg2(p)?"True":"False");
					fprintf(f,"%s",indsc);
					if(col && NoColors)
						fprintf(f,"\tMatrixTraceFactor -> %d,\n",col==1?3:8);
					fprintf(f,"\tMass -> %s,\n",CompoundArgN(p,5)?
							AtomValue(CompoundArgN(p,5)):"0");
					if(defined_em_charge())
						{
						Term ch=GetAtomProperty(a,A_EM_CHARGE);
						fprintf(f,"\tCharge -> ");
						if(ch==0)
							fprintf(f,"0,\n");
						else if(CompoundArg2(ch)==NewInteger(1))
							fprintf(f,"%ld,\n",IntegerValue(CompoundArg1(ch)));
						else
							fprintf(f,"%ld/%ld,\n",IntegerValue(CompoundArg1(ch)),
											IntegerValue(CompoundArg2(ch)));
						}
					fprintf(f,"\tPropagatorLabel -> \"%s\",\n",AtomValue(tnm));
					fprintf(f,"\tPropagatorType -> %s,\n",
							col==0?"Sine":"Cycles");
					fprintf(f,"\tPropagatorArrow -> %s }",
						CompoundArg1(p)==CompoundArg2(p)?"None":"Forward");
					SetAtomProperty(a,A_FANUM,MakeList2(NewInteger(2),NewInteger(nv)));
					if(CompoundArg1(p)!=CompoundArg2(p))
						SetAtomProperty(CompoundArg2(p),A_FANUM,
								MakeList2(NewInteger(2),NewInteger(-nv)));
					
					break;
				case 1:
					nf++;
					if(!firstprt) fprintf(f,",\n\n"); else firstprt=0;
					fprintf(f,"  F[%d] ==  { (* ",nf);
					fWriteTerm(f,CompoundArgN(p,3));
					fprintf(f," *)\n\tSelfConjugate -> %s,\n",
							CompoundArg1(p)==CompoundArg2(p)?"True":"False");
					fprintf(f,"%s",indsc);
					if(col && NoColors)
						fprintf(f,"\tMatrixTraceFactor -> %d,\n",col==1?3:8);
					fprintf(f,"\tMass -> %s,\n",CompoundArgN(p,5)?
							AtomValue(CompoundArgN(p,5)):"0");
					if(defined_em_charge())
						{
						Term ch=GetAtomProperty(a,A_EM_CHARGE);
						fprintf(f,"\tCharge -> ");
						if(ch==0)
							fprintf(f,"0,\n");
						else if(CompoundArg2(ch)==NewInteger(1))
							fprintf(f,"%ld,\n",IntegerValue(CompoundArg1(ch)));
						else
							fprintf(f,"%ld/%ld,\n",IntegerValue(CompoundArg1(ch)),
											IntegerValue(CompoundArg2(ch)));
						}
					fprintf(f,"\tPropagatorLabel -> \"%s\",\n",AtomValue(tnm));
					fprintf(f,"\tPropagatorType -> Straight,\n");
					fprintf(f,"\tPropagatorArrow -> %s }",
						CompoundArg1(p)==CompoundArg2(p)?"None":"Forward");
					SetAtomProperty(a,A_FANUM,MakeList2(NewInteger(1),NewInteger(nf)));
					if(CompoundArg1(p)!=CompoundArg2(p))
						SetAtomProperty(CompoundArg2(p),A_FANUM,
								MakeList2(NewInteger(1),NewInteger(-nf)));
					break;
				case 0:
					ns++;
					if(!firstprt) fprintf(f,",\n\n"); else firstprt=0;
					fprintf(f,"  S[%d] ==  { (* ",ns);
					fWriteTerm(f,CompoundArgN(p,3));
					fprintf(f," *)\n\tSelfConjugate -> %s,\n",
							CompoundArg1(p)==CompoundArg2(p)?"True":"False");
					fprintf(f,"%s",indsc);
					if(col && NoColors)
						fprintf(f,"\tMatrixTraceFactor -> %d,\n",col==1?3:8);
					fprintf(f,"\tMass -> %s,\n",CompoundArgN(p,5)?
							AtomValue(CompoundArgN(p,5)):"0");
					if(defined_em_charge())
						{
						Term ch=GetAtomProperty(a,A_EM_CHARGE);
						fprintf(f,"\tCharge -> ");
						if(ch==0)
							fprintf(f,"0,\n");
						else if(CompoundArg2(ch)==NewInteger(1))
							fprintf(f,"%ld,\n",IntegerValue(CompoundArg1(ch)));
						else
							fprintf(f,"%ld/%ld,\n",IntegerValue(CompoundArg1(ch)),
											IntegerValue(CompoundArg2(ch)));
						}
					fprintf(f,"\tPropagatorLabel -> \"%s\",\n",AtomValue(tnm));
					fprintf(f,"\tPropagatorType -> ScalarDash,\n");
					fprintf(f,"\tPropagatorArrow -> %s }",
						CompoundArg1(p)==CompoundArg2(p)?"None":"Forward");
					SetAtomProperty(a,A_FANUM,MakeList2(NewInteger(0),NewInteger(ns)));
					if(CompoundArg1(p)!=CompoundArg2(p))
						SetAtomProperty(CompoundArg2(p),A_FANUM,
								MakeList2(NewInteger(0),NewInteger(-ns)));
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
				if(!firstprt) fprintf(f,",\n\n"); else firstprt=0;
				fprintf(f,"  S[%d] ==  { (* ",ns);
				fWriteTerm(f,a);
				fprintf(f," *)\n\tSelfConjugate -> %s,\n",
						CompoundArg1(bpp)==CompoundArg2(bpp)?"True":"False");
					fprintf(f,"%s",indsc);
					if(col && NoColors)
						fprintf(f,"\tMatrixTraceFactor -> %d,\n",col==1?3:8);
				fprintf(f,"\tMass -> %s,\n",CompoundArgN(bpp,5)?
						AtomValue(CompoundArgN(bpp,5)):"0");
					if(defined_em_charge())
						{
						Term ch=GetAtomProperty(bp,A_EM_CHARGE);
						fprintf(f,"\tCharge -> ");
						if(ch==0)
							fprintf(f,"0,\n");
						else if(CompoundArg2(ch)==NewInteger(1))
							fprintf(f,"%ld,\n",IntegerValue(CompoundArg1(ch)));
						else
							fprintf(f,"%ld/%ld,\n",IntegerValue(CompoundArg1(ch)),
											IntegerValue(CompoundArg2(ch)));
						}
				fprintf(f,"\tPropagatorLabel -> \"%s\",\n",AtomValue(tnm));
				fprintf(f,"\tPropagatorType -> ScalarDash,\n");
				fprintf(f,"\tPropagatorArrow -> %s }",
					CompoundArg1(bpp)==CompoundArg2(bpp)?"None":"Forward");
				SetAtomProperty(a,A_FANUM,MakeList2(NewInteger(0),NewInteger(ns)));
				if(CompoundArg1(bpp)!=CompoundArg2(bpp))
					SetAtomProperty(GetAtomProperty(a,A_ANTI),A_FANUM,
							MakeList2(NewInteger(0),NewInteger(-ns)));
				continue;
			}
			
			
		}
	}

	fprintf(f,"}\n\n");
	NoQuotes=1;
	for(l1=all_prtc_list();l1;l1=ListTail(l1))
	{
		Term p;
		Atom a=ListFirst(l1);
		p=GetAtomProperty(a,PROP_TYPE);
		if(p==0) continue;
		if(!is_compound(p)) continue;
		
		if(CompoundName(p)==OPR_FIELD && GetAtomProperty(a,A_FANUM))
		{
			Term fp=GetAtomProperty(a,A_FANUM);
			int cl=(int)IntegerValue(ListFirst(fp));
			int no=(int)IntegerValue(ListFirst(ListTail(fp)));
			fprintf(f,"prt[\"");fWriteTerm(f,a);fprintf(f,"\"] = ");
			if(no<0) { fprintf(f,"-"); no=-no;}
			switch(cl){
				case 0: fprintf(f,"S[%d]\n",no); break;
				case 1: fprintf(f,"F[%d]\n",no); break;
				case 2: fprintf(f,"V[%d]\n",no); break;
				case 3: fprintf(f,"U[%d]\n",no); break;
					  }
		}
			
		if(CompoundName(p)==OPR_FIELD && CompoundArg2(p)==NewInteger(4))
			SetAtomProperty(a,A_FANUM,GetAtomProperty(CompoundArg1(p),A_FANUM));
		if(CompoundName(p)==OPR_PARTICLE)
		{
			Term fp=GetAtomProperty(a,A_FANUM);
			int cl, no;
			if(fp==0) continue;
			cl=(int)IntegerValue(ListFirst(fp));
			no=(int)IntegerValue(ListFirst(ListTail(fp)));
			fprintf(f,"prt[\"");fWriteTerm(f,a);fprintf(f,"\"] = ");
			if(no<0) { fprintf(f,"-"); no=-no;}
			switch(cl){
				case 0: fprintf(f,"S[%d]\n",no); break;
				case 1: fprintf(f,"F[%d]\n",no); break;
				case 2: fprintf(f,"V[%d]\n",no); break;
					  }
		}
	}
	
	cls_wrt_nms(f,NoColors?0:cpml);
	if(cpml) FreeAtomic(cpml);

			
	NoQuotes=0;

	fprintf(f,"\nGaugeXi[_] = 1\n\n");

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
					if(GetAtomProperty(CompoundArg1(ListFirst(l3)),A_INFINITESIMAL) &&
                                        !GetAtomProperty(CompoundArg1(ListFirst(l3)),PROP_TYPE)
                                        )
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
		alg2_fix_mom2(a2);
/*		alg2_red_cos(a2);
		alg2_red_orth(a2);
		
		{
			List ml=ConsumeCompoundArg(a2,5);
			ml=SortedList(ml,mlsort);
			SetCompoundArg(a2,5,ml);
		}

		alg2_red_sico(a2);
		alg2_red_comsico(a2);
*/
		if(CompoundArgN(a2,5)==0) continue;
       
		alg2_red_1pm5(a2);
       
		alg2_fix_ff(a2);
		if(opAbbrVrt)
		{
			alg2_decommon_s(a2);
			alg2_abbr_vrt(a2);
			alg2_common_s(a2);
		}
		alg2_recommon_n(a2);
		if(!opAbbrVrt)
			{
			alg2_decommon_s(a2);
			allow_sym_div=1;
			alg2_common_s(a2);
			allow_sym_div=0;
			}
		alg2_multbyi(a2);
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
		if(ListLength(CompoundArg1(a2))>4)
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
			continue;
		}
		//printf("%d\n",ListLength(CompoundArg1(a2)));

		alg2_common_t(a2);
		
			
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

	inf_decl_hc(f);
	prm_decl_hc(f,l);

	
	fprintf(f,"\nM$CouplingMatrices = {\n\n");
	
	
	for(l1=l;l1;l1=ListTail(l1))
	{
		Term a2=ListFirst(l1);
		List l2,lp;
		List *lv, *lv1,*lv2;
		int i,lvl,hasz=0;
		cv=CompoundArg1(a2);
		if(CompoundArgN(a2,5)==0 || ListLength(CompoundArg1(a2))<2 || 
				ListLength(CompoundArg1(a2))>4)
			continue;
		lp=NewList();
		fprintf(f,"  (*------   ");
		for(l2=CompoundArg1(a2),i=0;l2;l2=ListTail(l2),i++)
		{
			Atom a;
			fprintf(f," %s ",AtomValue(CompoundArg1(ListFirst(l2))));
			a=CompoundName(CompoundArg2(ListFirst(l2)));
			if(a==OPR_SCALAR && 
				GetAtomProperty(CompoundArg1(ListFirst(l2)),A_GRASS)) a=A_GRASS;
			lp=AppendLast(lp,a);
		}
		
						
		for(i=0;i<4;i++) pind[i]=wind[i]=0;
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
		
		fprintf(f," ------*)\n   C[ ");
		for(l2=CompoundArg1(a2),i=0;l2;l2=ListTail(l2),i++)
		{
			Atom a;
			Integer fai;
			char c;
			int ii=0;
			a=CompoundName(CompoundArg2(ListFirst(l2)));
			if(a==OPR_SCALAR && GetAtomProperty(CompoundArg1(ListFirst(l2)),A_GRASS))
				a=A_GRASS;
			if(a==OPR_SCALAR) c='S'; else if(a==OPR_SPINOR) c='F';
			else if(a==OPR_VECTOR) c='V'; else if(a==A_GRASS) c='U'; else c='?';
			fai=GetAtomProperty(CompoundArg1(ListFirst(l2)),A_FANUM);
			if(fai==0) printf("Internal error - miss no for %s\n",
				AtomValue(CompoundArg1(ListFirst(l2))));
			else
				ii=(int)IntegerValue(ListFirst(ListTail(fai)));
			if(ii<0) {fprintf(f,"-");ii=-ii;}
			fprintf(f,"%c[%d",c,ii);
			if(pind[i]||wind[i]) 
			{
				if(pind[i] && wind[i])
					fprintf(f,",{t%d, c%d}]",i+1,i+1);
				else
					fprintf(f,",{%c%d}]",pind[i]?'c':'t',i+1);
			}
			else        fprintf(f,"]");
			if(ListTail(l2)) fprintf(f,", ");
		}
		fprintf(f," ] == ");

		

				
		wrt_expr(f,CompoundArgN(a2,2),CompoundArgN(a2,3),CompoundArgN(a2,4));
		fprintf(f,"\n");
		
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
				if((pr=GetAtomProperty(CompoundArg1(ListFirst(l3)),A_INFINITESIMAL))
                                    && !GetAtomProperty(CompoundArg1(ListFirst(l3)),PROP_TYPE)
                                )
					io+=(int)IntegerValue(CompoundArg1(pr))*
						(int)IntegerValue(CompoundArg2(ListFirst(l3)));
			for(l3=CompoundArgN(ListFirst(l2),3);l3;l3=ListTail(l3))
				if((is_atom(CompoundArg1(ListFirst(l3)))&& 
				  (pr=GetAtomProperty(CompoundArg1(ListFirst(l3)),A_INFINITESIMAL))))
					io+=1;
			if(io==1)
				lv1[i]=AppendLast(lv1[i],ListFirst(l2)),(hasz==0?hasz=1:0);
			else if(io==0)
				lv[i]=AppendLast(lv[i],ListFirst(l2));
			else
				lv2[i]=AppendLast(lv2[i],ListFirst(l2)),hasz=2;
				
		}
		
		
		
		fprintf(f,"{ \n");
		
		for(i=0;i<4;i++) pind[i]=0;
		for(l2=CompoundArg1(a2),i=0;l2;l2=ListTail(l2),i++)
		{
			Atom a=CompoundName(CompoundArg2(ListFirst(l2)));
			if(GetAtomProperty(CompoundArg1(ListFirst(l2)),A_COLOR))
				pind[i]=ListNth(CompoundArg1(CompoundArg2(ListFirst(l2))),
						(a==OPR_VECTOR||a==OPR_SPINOR)?2:1);
			else
				pind[i]=0;
		}
		
		for(i=0;i<lvl;i++)
		{
			/*fWriteTerm(f,lv1[i]);fprintf(f,"\n");fflush(f);*/
			if(lv1[i]) lv1[i]=SortedList(lv1[i],pcmp);
			/*fWriteTerm(f,lv1[i]);fprintf(f,"\n");*/
			
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
		}
		fprintf(f,"}%c\n",ListTail(l1)?',':' ');
		free(lv);free(lv1);free(lv2);
	}
	
	fprintf(f,"}\n\n");
	
	fprintf(f,"M$LastModelRules = {}\n\n");

	FADeclRealParam(f);
		
	FAsqparam(f);
	
	inf_write_rc(f);
	
	cls_write_dmatr(f);
	
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


static Term conv_lor(List pl, Term m2, int eff)
{
	int i;
	List l;
	List nc=0, cc=0;
	if(ili[0]==0)
	{
		char cbuf[16];
		for(i=1;i<=4;i++)
		{
			sprintf(cbuf,"li%d",i); ili[i-1]=NewAtom(cbuf,0);
			sprintf(cbuf,"mom%d",i); imom[i-1]=NewAtom(cbuf,0);
		}
	}
	
	for(l=pl,i=0;l;l=ListTail(l),i++)
	{
		Atom a;
		a=CompoundName(CompoundArg2(ListFirst(l)));
		if(a==OPR_VECTOR || a==OPR_SPINOR)
			pind[i]=ListFirst(CompoundArg1(CompoundArg2(ListFirst(l))));
		else
			pind[i]=0; 
	}
	
	
	if(CompoundName(CompoundArg2(ListFirst(pl)))==OPR_SPINOR)
	{
		long int curi=pind[0];
		for(l=CompoundArgN(m2,3);l&&curi!=pind[1];l=ListTail(l))
		{
			rpt:
			if(ListFirst(l)==0) continue;
			if(curi==pind[1]) break;
			if(CompoundName(ListFirst(l))==OPR_SPECIAL && 
					ListFirst(CompoundArg2(ListFirst(l)))==curi)
			{
				curi=ListNth(CompoundArg2(ListFirst(l)),2);
				if(CompoundArg1(ListFirst(l))==A_GAMMA)
				{
					Integer gi=ListNth(CompoundArg2(ListFirst(l)),3);
					List l2;
					for(i=0;i<4;i++) if(pind[i]==gi) break;
					if(i<4)
					{
						nc=AppendLast(nc,MakeCompound1(A_GAMMA,ili[i]));
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
				}
				if(eff){FreeAtomic(ListFirst(l));ChangeList(l,0);}
				l=CompoundArgN(m2,3);
				goto rpt;
			}
		}
		
	}

		
	for(l=CompoundArgN(m2,3);l;l=ListTail(l))
	{
		if(ListFirst(l)==0)
			continue;
		if(CompoundArg1(ListFirst(l))==A_DELTA)
		{
			Atom i1,i2;
			for(i=0;i<4;i++)
				if(ListFirst(CompoundArg2(ListFirst(l)))==pind[i])
					break;
			if(i==4)
				continue;
			i1=ili[i];
			for(i++;i<4;i++)
				if(ListFirst(ListTail(CompoundArg2(ListFirst(l))))==pind[i])
					break;
			if(i==4)
			{
			printf("Internal error: ");WriteVertex(pl);printf("lost delta index:");
			WriteTerm(m2);
			puts("");
			continue;
			}
			i2=ili[i];
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
			Atom in[4];
			int j;
			List lj,l2;
			for(j=0,lj=CompoundArg2(ListFirst(l));j<4;j++,lj=ListTail(lj))
			{
			in[j]=0;
			for(i=0;i<4;i++)
				if(ListFirst(lj)==pind[i])
					break;
			if(i<4)
				{
				in[j]=ili[i];
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
			for(i=0;i<4;i++)
				if(pi==pind[i])
					break;
			if(i<4)
			{
				cc=AppendLast(cc,MakeCompound2(A_DELTA,ili[i],imom[p1-1]));
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
	
	for(l=ten;l;l=ListTail(l))
	{
		
		if(CompoundArg1(ListFirst(l))==A_DELTA)
		{
			Integer in1,in2;
			int i,il1, il2;
			in1=ListFirst(CompoundArg2(ListFirst(l)));
			in2=ListFirst(ListTail(CompoundArg2(ListFirst(l))));
			for(i=0;i<4;i++) if(in1==pind[i]) break;
			if(i==4)
			{
				for(i=0;i<4;i++) if(in1==wind[i]) break;
				il1=i;
				if(i==4) puts("Internal error wrtfa07");
				for(i=0;i<4;i++) if(in2==wind[i]) break;
				il2=i;
				if(i==4) puts("Internal error wrtfa07");
				sno+=fprintf(of,"IndexDelta[t%d, t%d] ",il1+1,il2+1);
				if(sno>60 && ListTail(l)) {fprintf(of,"*\n\t\t");sno=15;}
				continue;
			}
			il1=i;
			for(i=0;i<4;i++) if(in2==pind[i]) break;
			if(i==4){puts("Internal error wrtfa07");}
			il2=i;
			
			sno+=fprintf(of,"IndexDelta[c%d, c%d] ",il1+1,il2+1);
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
				for(i=0;i<4;i++) if(in1==wind[i]) break;
				if(i==4) puts("Internal error wrtfa07");
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
			
			for(i=0;i<4;i++) if(in1==pind[i]) break;
			il1=(i==4?(int)IntegerValue(in1)+4:i);
			
			for(i=0;i<4;i++) if(in2==pind[i]) break;
			il2=(i==4?(int)IntegerValue(in2)+4:i);
			
			for(i=0;i<4;i++) if(in3==pind[i]) break; il3=i;
			il3=(i==4?(int)IntegerValue(in3)+4:i);
						
			if(il1<5 && il2<5 && il3<5)
			{
			sno+=fprintf(of,"SUNT[c%d, c%d, c%d] ",il3+1,il1+1,il2+1);
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
			for(i=0;i<4;i++) if(in3==pind[i]) break; il2=i;
			if(i==4){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(2)");il3=(int)IntegerValue(in3);}
			for(i=0;i<4;i++) if(in4==pind[i]) break; il4=i;
			if(i==4){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(2)");il4=(int)IntegerValue(in4);}
			sno+=fprintf(of,"SUNT[c%d, c%d, c%d, c%d] ",
					il3+1,il4+1,il1+1,il2+1);
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
			for(i=0;i<4;i++) if(in3==pind[i]) break; il3=i;
			if(i==4){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(3)");il3=(int)IntegerValue(in3);}
			for(i=0;i<4;i++) if(in4==pind[i]) break; il4=i;
			if(i==4){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(3)");il4=(int)IntegerValue(in4);}
			sno+=fprintf(of,"SUNTSum[c%d, c%d, c%d, c%d] ",
					il3+1,il4+1,il1+1,il2+1);
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
			for(i=0;i<4;i++) if(in1==pind[i]) break; il1=i;
			if(i==4){WriteVertex(cv);WriteTerm(ListFirst(l));
			puts(": color structure error(1)");il1=(int)IntegerValue(in1);}
			for(i=0;i<4;i++) if(in2==pind[i]) break; il2=i;
			if(i==4){WriteVertex(cv);WriteTerm(ListFirst(l));
			puts(": color structure error(1)");il2=(int)IntegerValue(in2);}
			for(i=0;i<4;i++) if(in3==pind[i]) break; il3=i;
			if(i<4)
			{
				sno+=fprintf(of,"SUNF[c%d, c%d, c%d] ",il1+1,il2+1,il3+1);
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
			for(i=0;i<4;i++) if(in3==pind[i]) break; il3=i;
			if(i==4){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(2)");il3=(int)IntegerValue(in3);}
			for(i=0;i<4;i++) if(in4==pind[i]) break; il4=i;
			if(i==4){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(2)");il4=(int)IntegerValue(in4);}
			sno+=fprintf(of,"SUNF[c%d, c%d, c%d, c%d] ",il1+1,il2+1,il3+1,il4+1);
			if(sno>60 && ListTail(l) && (ListTail(l)!=l2 || ListTail(l2))) 
				{fprintf(of,"*\n\t\t");sno=15;}
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
			for(i=0;i<4;i++) if(in1==pind[i]) break; il1=i;
			if(i==4){WriteVertex(cv);WriteTerm(ListFirst(l));
			puts(": color structure error(1)");il1=(int)IntegerValue(in1);}
			for(i=0;i<4;i++) if(in2==pind[i]) break; il2=i;
			if(i==4){WriteVertex(cv);WriteTerm(ListFirst(l));
			puts(": color structure error(1)");il2=(int)IntegerValue(in2);}
			for(i=0;i<4;i++) if(in3==pind[i]) break; il3=i;
			if(i<4)
			{
				sno+=fprintf(of,"SUND[c%d, c%d, c%d] ",il1+1,il2+1,il3+1);
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
			for(i=0;i<4;i++) if(in3==pind[i]) break; il3=i;
			if(i==4){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(2)");il3=(int)IntegerValue(in3);}
			for(i=0;i<4;i++) if(in4==pind[i]) break; il4=i;
			if(i==4){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(2)");il4=(int)IntegerValue(in4);}
			sno+=fprintf(of,"SUNDF[c%d, c%d, c%d, c%d] ",il1+1,il2+1,il3+1,il4+1);
			if(sno>60 && ListTail(l)&& (ListTail(l)!=l2 || ListTail(l2))) {fprintf(of,"*\n\t\t");sno=15;}
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
			for(i=0;i<4;i++) if(in3==pind[i]) break; il3=i;
			if(i==4){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(2)");il3=(int)IntegerValue(in3);}
			for(i=0;i<4;i++) if(in4==pind[i]) break; il4=i;
			if(i==4){WriteVertex(cv);WriteTerm(ListFirst(l2));
			puts(": color structure error(2)");il4=(int)IntegerValue(in4);}
			sno+=fprintf(of,"SUND[c%d, c%d, c%d, c%d] ",il1+1,il2+1,il3+1,il4+1);
			if(sno>60 && ListTail(l)&& (ListTail(l)!=l2 || ListTail(l2))) {fprintf(of,"*\n\t\t");sno=15;}
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

List opFAGS=0, opFAGE=0;

void FA_write_gen(FILE *f)
{
	List l1;
	time_t tm;
	
	NoQuotes=1;
	fprintf(f,"(*\n\tLanHEP output produced at ");
	time(&tm);
	fprintf(f,"%s",ctime(&tm));
	fprintf(f,"\tfrom the file '%s'\n",eff_infile);
	if(ModelName)
		fprintf(f,"\tModel named '%s'\n",ModelName);
	fprintf(f,"*)\n\n");
	
	for(l1=opFAGS;l1;l1=ListTail(l1))
	{
		if(ListFirst(l1)!=NewInteger(0))
			fWriteTerm(f,ListFirst(l1));
		fprintf(f,"\n");
	}
	
	fprintf(f,"\n\t(* Generic analytical couplings for the model *)\n\n");
	fprintf(f,"M$GenericCouplings = {\n\n");
	
	for(l1=gkl;l1;l1=ListTail(l1))
	{
		char cbuf[16];
		int i,n;
		List l2;
		fprintf(f,"\t(* ");
		for(l2=CompoundArg1(ListFirst(l1)),i=0;l2;l2=ListTail(l2),i++)
		{
			Atom a=ListFirst(l2);
			if(a==OPR_SCALAR) cbuf[i]='S';
			else if(a==OPR_SPINOR) cbuf[i]='F';
			else if(a==OPR_VECTOR) cbuf[i]='V';
			else if(a==A_GRASS) cbuf[i]='U';
			else {puts("Internal error wrta01");cbuf[i]='?';}
			fprintf(f,"%c%s",cbuf[i],ListTail(l2)?"-":" *)\n");
		}
		n=i; cbuf[i]=0;
		fprintf(f,"    AnalyticalCoupling[ ");
		for(i=1;i<=n;i++)
		{
			fprintf(f,"s%d %c[j%d, mom%d",i,cbuf[i-1],i,i);
			if(cbuf[i-1]=='V') fprintf(f,", {li%d}]",i);
			else fprintf(f,"]");
			if(i<n) fprintf(f,", ");
			if(i<n && (i%2==0)) fprintf(f,"\n\t\t");
		}
		fprintf(f,"] ==\n   G[%d][",(cbuf[0]=='F'&&cbuf[1]=='F'&&
			(cbuf[2]==0||(cbuf[2]=='V'&&cbuf[3]==0)))?-1:1);
		for(i=1;i<=n;i++)
			fprintf(f," s%d %c[j%d]%c",i,cbuf[i-1],i,i==n?']':',');
		fprintf(f," .\n    {");
		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
		{
			List cc=CompoundArg1(ListFirst(l2)),nc=CompoundArg2(ListFirst(l2));
			List l3;
			for(l3=cc;l3;l3=ListTail(l3))
			{
				Atom a1=CompoundArg1(ListFirst(l3)),
					 a2=CompoundArg2(ListFirst(l3));
				if(CompoundArity(ListFirst(l3))==4)
					fprintf(f," Epsilon[%s, %s, %s, %s]",AtomValue(a1),AtomValue(a2),
						AtomValue(CompoundArgN(ListFirst(l3),3)),AtomValue(CompoundArgN(ListFirst(l3),4)));
				else
				if(AtomValue(a1)[0]=='l' && AtomValue(a2)[0]=='l')
					fprintf(f," MetricTensor[%s, %s]",AtomValue(a1),AtomValue(a2));
				else if(AtomValue(a1)[0]=='l' && AtomValue(a2)[0]=='m')
					fprintf(f," FourVector[%s, %s]",AtomValue(a2),AtomValue(a1));
				else if(AtomValue(a1)[0]=='m' && AtomValue(a2)[0]=='m')
					fprintf(f," ScalarProduct[%s, %s]",AtomValue(a1),AtomValue(a2));
				else {printf("Internal error wrtfa02; ");WriteTerm(ListFirst(l3));puts("");}
			}
			if(nc)
			{
				fprintf(f," NonCommutative[");
				for(l3=nc;l3;l3=ListTail(l3))
				{
					if(ListFirst(l3)==A_GAMMAP) fprintf(f,"ChiralityProjector[+1]");
					else if(ListFirst(l3)==A_GAMMAM) fprintf(f,"ChiralityProjector[-1]");
					else if(ListFirst(l3)==A_GAMMA5) 
					{
						puts("warn: gamma5");
						fprintf(f,"(ChiralityProjector[+1]-ChiralityProjector[-1])");
					}
					else if(is_compound(ListFirst(l3)) &&
						CompoundName(ListFirst(l3))==A_GAMMA)
					{
						Atom a = CompoundArg1(ListFirst(l3));
						if(AtomValue(a)[0]=='l')
							fprintf(f,"DiracMatrix[%s]",AtomValue(a));
						else if(AtomValue(a)[0]=='m')
							fprintf(f,"DiracSlash[%s]",AtomValue(a));
						else {puts("Internal error wrtfa03");}
					}
					else {puts("Internal error wrtfa04");}
					if(ListTail(l3)) fprintf(f,", ");
				}
				fprintf(f,"]");
			}
			if(cc==0 && nc==0) fprintf(f," 1");
			if(ListTail(l2)) fprintf(f,",\n     ");
			else fprintf(f," }");
		}
		if(ListTail(l1)) fprintf(f,",\n\n");
		else fprintf(f,"\n}\n\n");

										
	}
	
	for(l1=opFAGE;l1;l1=ListTail(l1))
	{
		if(ListFirst(l1)!=NewInteger(0))
				fWriteTerm(f,ListFirst(l1));
		fprintf(f,"\n");
	}
	
	NoQuotes=0;
}

Term ProcFAGC(Term t, Term ind)
{
	t=ConsumeCompoundArg(t,1);
	gkl=AppendLast(gkl,t);
	return 0;
}

void alg2_fix_ff(Term a2)
{
  List l;
  l=CompoundArg1(a2);
  if(ListLength(l)!=2)
    return;
  if(CompoundName(CompoundArg2(ListFirst(l)))!=OPR_SPINOR ||
     CompoundName(CompoundArg2(ListFirst(ListTail(l))))!=OPR_SPINOR)
    return;
  for(l=CompoundArgN(a2,5);l;l=ListTail(l))
    {
      Term m=ListFirst(l);
      List l1;
      Term mom=0, gpm=0;
      for(l1=CompoundArgN(m,3);l1;l1=ListTail(l1))
	{
	  Term s=ListFirst(l1);
	  if(CompoundName(s)==A_MOMENT)
	    mom=s;
	  if(CompoundName(s)==OPR_SPECIAL &&(CompoundArg1(s)==A_GAMMAP ||
					     CompoundArg1(s)==A_GAMMAM))
	    gpm=CompoundArg1(s);
	}
      if(mom && CompoundArg1(mom)==NewInteger(1) && gpm==A_GAMMAP)
	{
	  int cf;
	  SetCompoundArg(mom,1,NewInteger(2));
	  cf=(int)IntegerValue(CompoundArg1(m));
	  SetCompoundArg(m,1,NewInteger(-cf));
	}
      if(mom && CompoundArg1(mom)==NewInteger(2) && gpm==A_GAMMAM)
	{
	  int cf;
	  SetCompoundArg(mom,1,NewInteger(1));
	  cf=(int)IntegerValue(CompoundArg1(m));
	  SetCompoundArg(m,1,NewInteger(-cf));
	}
	
    }
	  
}

void alg2_fix_mom2(Term a2)
{
	List l=CompoundArg1(a2);
	return;
	if(ListLength(l)!=2)
		return;
	for(l=CompoundArgN(a2,5);l;l=ListTail(l))
	{
	Term d1=0, d2=0;
	int dno=0;
    List l2;
	for(l2=CompoundArgN(ListFirst(l),3);l2;l2=ListTail(l2))
		{
		Term m=ListFirst(l2);
	    if(is_compound(m) && CompoundName(m)==A_MOMENT)
		{
			dno++;
			if(dno==1) d1=m;
			if(dno==2) d2=m;
		}
		}
	if(dno==2)
	{
		int i1=(int)IntegerValue(ListFirst(CompoundArg2(d1)));
		int i2=(int)IntegerValue(ListFirst(CompoundArg2(d2)));
		int ch=0;
		if(i1!=i2 && i1<3 && i2<3 &&
			CompoundArg1(d1)==CompoundArg1(d2))
		{
			SetCompoundArg(d1,1,NewInteger(i2));
			SetCompoundArg(d2,1,NewInteger(i1));
			ch=1;
		}
		else
		if(CompoundArg1(d1)==NewInteger(1) && CompoundArg1(d2)==NewInteger(1))
		{
			SetCompoundArg(d2,1,NewInteger(2));
			ch=1;
		} else
		if(CompoundArg1(d1)==NewInteger(2) && CompoundArg1(d2)==NewInteger(2))
		{
			SetCompoundArg(d1,1,NewInteger(1));
			ch=1;
		}
		if(ch)
		{
			
			Term m=ListFirst(l);
		//	WriteTerm(m);puts("");
			SetCompoundArg(m,1,NewInteger(-IntegerValue(CompoundArg1(m))));
		}
		
	}
	}
	
}

void alg2_fix_uuv(Term a2)
{
	List l,ml,sl;
	int i;
	Term m,m2;
	l=CompoundArg1(a2);
	if(ListLength(l)!=3)
		return;
	for(i=1;l;l=ListTail(l),i++)
	{
		Atom s,p;
		p=CompoundArg1(ListFirst(l));
		s=CompoundName(CompoundArg2(ListFirst(l)));
		if(s==OPR_SCALAR && GetAtomProperty(p,A_GRASS))
			s=A_GRASS;
		if( ((i==1||i==2)&&s!=A_GRASS) || (i==3&&s!=OPR_VECTOR))
			return;
	}
	/*	WriteTerm(CompoundArg2(a2));
		DumpList(CompoundArgN(a2,5));*/
	ml=ConsumeCompoundArg(a2,5);
	for(l=ml;l;l=ListTail(l))
	{
		m2=ListFirst(l);
		sl=CompoundArgN(m2,3);
		if(CompoundArg1(m2)!=NewInteger(1) || CompoundArg2(m2)!=0 ||
				ListLength(sl)!=1)
			continue;
		m=ListFirst(sl);
		if(CompoundName(m)!=A_MOMENT || CompoundArg1(m)!=NewInteger(3))
			continue;
		break;
	}
	if(l==0)
	  {
	    SetCompoundArg(a2,5,ml);
	    return;
	  }
	ChangeList(l,0);
	ml=CutFromList(ml,l);
	
	SetCompoundArg(m,1,NewInteger(2));
	for(l=ml;l;l=ListTail(l))
		if(EqualTerms(ListFirst(l),m2))
			break;
	if(l)
		ml=CutFromList(ml,l);
	else
	{
		Term mc=CopyTerm(m2);
		SetCompoundArg(mc,1,NewInteger(-IntegerValue(CompoundArg1(mc))));
		ml=AppendFirst(ml,mc);
	}
	
	SetCompoundArg(m,1,NewInteger(1));
	for(l=ml;l;l=ListTail(l))
		if(EqualTerms(ListFirst(l),m2))
			break;
	if(l)
		ml=CutFromList(ml,l);
	else
	{
		Term mc=m2;
		SetCompoundArg(mc,1,NewInteger(-IntegerValue(CompoundArg1(mc))));
		ml=AppendFirst(ml,mc);
	}
	
	for(l=ml;l;l=ListTail(l))
	{
		m2=ListFirst(l);
		SetCompoundArg(m2,1,NewInteger(-IntegerValue(CompoundArg1(m2))));
	}
	
	m2=CompoundArg2(a2);
	SetCompoundArg(m2,1,NewInteger(-IntegerValue(CompoundArg1(m2))));
	
	/*	DumpList(ml);*/
	SetCompoundArg(a2,5,ml);
}

int opNo4Scal=0;

void alg2_red_rc(Term a2)
{
	List l1,l2,ml,mln;
	int rcf=1,n=0;
	
	if(opNo4Scal==0)
		return;
	/*WriteTerm(CompoundArg1(a2));puts("");*/
	if(ListLength(CompoundArg1(a2))!=4)
		return;

	for(l1=CompoundArg1(a2);l1;l1=ListTail(l1))
	{
		Atom p=CompoundName(CompoundArg2(ListFirst(l1)));
		if(p!=OPR_SCALAR)
			rcf=0;
	}
	if(rcf==0)
		return;
	{
		/*WriteTerm(CompoundArg1(a2));puts("");*/
		FreeAtomic(ConsumeCompoundArg(a2,5));
	}
	return;
	for(l1=CompoundArg1(a2);l1;l1=ListTail(l1))
	{
		Atom p=CompoundArg1(ListFirst(l1));
		Term prp=GetAtomProperty(p,PROP_TYPE);
		if(prp && CompoundName(prp)==OPR_FIELD &&
				CompoundArg2(prp)==NewInteger(4))
			p=CompoundArg1(prp);
		if(GetAtomProperty(p,OPR_LET)==0)
			rcf=0,n++;
	}
	
	if(rcf)
		return;
/*	if(n==4)
	{
		WriteTerm(CompoundArg1(a2));puts("");
		FreeAtomic(ConsumeCompoundArg(a2,5));
		return;
	}
*/	
	/*return;*/
	
	ml=ConsumeCompoundArg(a2,5);
	mln=NewList();
	
	
	for(l1=ml;l1;l1=ListTail(l1))
	{
		rcf=0;
		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
			if(GetAtomProperty(CompoundArg1(ListFirst(l2)),A_INFINITESIMAL) &&
                            !GetAtomProperty(CompoundArg1(ListFirst(l2)),PROP_TYPE)
                        )
				rcf++;
		if(rcf)
			FreeAtomic(ListFirst(l1));
		else
			mln=AppendLast(mln,ListFirst(l1));
	}
	
	RemoveList(ml);
	SetCompoundArg(a2,5,mln);
	
}


Term ProcFainclude(Term t, Term ind)
	{
	int i, n;
	n=CompoundArity(t);
	for(i=1;i<=n;i++)
		{
		Atom a=CompoundArgN(t,i);
		if(!is_atom(a))
			{
			ErrorInfo(2000);
			printf("fainclude: ");WriteTerm(a);
			puts("is not a string (use quotes).");
			return 0;
			}
		fainclude=AppendLast(fainclude,a);
		}
	return 0;
	}

