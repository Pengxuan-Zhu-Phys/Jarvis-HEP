#include <stdio.h>
#include <unistd.h>
#include "lanhep.h"
#include <math.h>
#include <string.h>

extern int extlib_problems;
		
static int acmp(Term a1, Term a2)
	{
	return strcmp(AtomValue(a1),AtomValue(a2));
	}

static double mysqrt(double x)
{
	if(x<0)
		return -sqrt(-x);
	else
		return sqrt(x);
}
		
static int mcmp(Term m1, Term m2)
	{
	List l1,l2;
	l1=CompoundArg1(m1);
	l2=CompoundArg1(m2);
	if(ListLength(l1)>ListLength(l2))
		return 1;
	if(ListLength(l1)<ListLength(l2))
		return -1;
	return strcmp(AtomValue(ListFirst(l1)),
					AtomValue(ListFirst(l2)));
	}

static int prt_spin1(Atom p)
{
	Term prop;
	int spi;
	
	prop=GetAtomProperty(p,PROP_TYPE);
	if(CompoundName(prop)==OPR_FIELD)
		prop=GetAtomProperty(CompoundArg1(prop),PROP_TYPE);
	spi=(int)IntegerValue(CompoundArgN(prop,4));
	if(spi==2 && CompoundArgN(prop,7)!=A_GAUGE)
		spi++;
	return spi;
}

	
static int prt_spin(Atom p)
	{
	Term prop;
	prop=GetAtomProperty(p,PROP_INDEX);
	if(is_empty_list(prop) || CompoundName(CompoundArg1(ListFirst(prop)))!=A_LORENTZ)
		return 0;
	prop=CompoundArg1(ListFirst(prop));
	if(CompoundArg1(prop)==NewInteger(2))
		return 2;
	return 1;
	}

	
static Atom remake_sp(List m3)
{
	Term mom1, mom2;
	char cbuf[64];
	List l;
	int ga=0;
	Atom ag=0, re=0;
	mom1=mom2=0;
	for(l=m3;l;l=ListTail(l))
	{
		if(CompoundName(ListFirst(l))==A_MOMENT)
		{
			if(mom1==0)
				mom1=ListFirst(l);
			else
				mom2=ListFirst(l);
		}
		if(CompoundName(ListFirst(l))==OPR_SPECIAL &&
				CompoundArg1(ListFirst(l))==A_GAMMA)
			ga++;
		if(CompoundName(ListFirst(l))==OPR_SPECIAL &&
				CompoundArg1(ListFirst(l))==A_GAMMA5)
			ag=A_GAMMA5;
		if(CompoundName(ListFirst(l))==OPR_SPECIAL &&
				CompoundArg1(ListFirst(l))==A_GAMMAM)
			ag=A_GAMMAM;
		if(CompoundName(ListFirst(l))==OPR_SPECIAL &&
				CompoundArg1(ListFirst(l))==A_GAMMAP)
			ag=A_GAMMAP;
	}
	
	if(mom2==0 && ga)
		re=NewAtom("G(p1)",0);
	else if(mom2==0)
		re= NewAtom("p1^mu",0);
	else if(ListFirst(CompoundArg2(mom1))==ListFirst(CompoundArg2(mom2)))
		re= NewAtom("p1**2",0);
	else re=NewAtom("p1^mu*p1^nu",0);
	if(ag==A_GAMMA5)
		{
		sprintf(cbuf,"%s*gamma5",AtomValue(re)); re=NewAtom(cbuf,0);
		}
	if(ag==A_GAMMAP)
		{
		sprintf(cbuf,"%s*(1+gamma5)/2",AtomValue(re)); re=NewAtom(cbuf,0);
		}
	if(ag==A_GAMMAM)
		{
		sprintf(cbuf,"%s*(1-gamma5)/2",AtomValue(re)); re=NewAtom(cbuf,0);
		}
	return re;
}


static void wrt_ml(FILE *f, List ml)
	{
	List gl;
	List l;
	if(ml==0)
		{
		fprintf(f,"null");
		return;
		}
	if(is_atom(ml))
		{
		fWriteTerm(f,ml);
		return;
		}
	gl=NewList();
	for(l=ml;l;l=ListTail(l))
		{
		Term sp;
		if(is_empty_list(CompoundArgN(ListFirst(l),3)))
			sp=0;
		else
			sp=CompoundArg1(ListFirst(CompoundArgN(ListFirst(l),3)));
		if(is_integer(sp))
			sp=remake_sp(CompoundArgN(ListFirst(l),3));
		if(!ListMember(gl,sp))
			gl=AppendLast(gl,sp);
		}
	
	for(l=gl;l;l=ListTail(l))
		{
		List l1;
		fprintf(f,"(");
		for(l1=ml;l1;l1=ListTail(l1))
			{
			Term sp;
			List l2;
			if(is_empty_list(CompoundArgN(ListFirst(l1),3)))
				sp=0;
			else
				sp=CompoundArg1(ListFirst(CompoundArgN(ListFirst(l1),3)));
			if(is_integer(sp))
				sp=remake_sp(CompoundArgN(ListFirst(l1),3));
			if(ListFirst(l)!=sp)
				continue;
			{
			long n,d;
			n=IntegerValue(CompoundArg1(CompoundArg1(ListFirst(l1))));
			d=IntegerValue(CompoundArg2(CompoundArg1(ListFirst(l1))));
			fprintf(f,"%+ld",n);
			if(d!=1)
				fprintf(f,"/%ld",d);
			}
			for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
				{
				Atom p;
				int pw;
				p=CompoundArg1(ListFirst(l2));
				pw=(int)IntegerValue(CompoundArg2(ListFirst(l2)));
				if(pw>0)
					fprintf(f,"*%s",AtomValue(p));
				else
					fprintf(f,"/%s",AtomValue(p));
				if(pw<0)
					pw=-pw;
				if(pw>1)
					fprintf(f,"**%d",pw);
				}
			}
			
		fprintf(f,")");
		if(ListFirst(l) && ListFirst(l)!=A_DELTA)
			{
			if(is_atom(ListFirst(l)))
				fprintf(f,"*%s",AtomValue(ListFirst(l)));
			else
				fprintf(f,"*p%ld",IntegerValue(ListFirst(l)));
			}
		if(ListTail(l))
			fprintf(f,"+");
		}
	FreeAtomic(gl);
	}

/*
static void wrt_num(FILE *f, List ml, int spi)
	{
	List gl;
	List l;
	gl=NewList();
	for(l=ml;l;l=ListTail(l))
		{
		Term sp;
		if(is_empty_list(CompoundArgN(ListFirst(l),3)))
			sp=0;
		else
			sp=CompoundArg1(ListFirst(CompoundArgN(ListFirst(l),3)));
		if(!ListMember(gl,sp))
			gl=AppendLast(gl,sp);
		}
	
	for(l=gl;l;l=ListTail(l))
		{
		List l1;
		double vi,vr;
		int ni,nr;
		vi=0.0;
		vr=0.0;
		ni=0;
		nr=0;
		
		for(l1=ml;l1;l1=ListTail(l1))
			{
			Term sp;
			double v1;
			List l2;
			int is_i=0;
			if(is_empty_list(CompoundArgN(ListFirst(l1),3)))
				sp=0;
			else
				sp=CompoundArg1(ListFirst(CompoundArgN(ListFirst(l1),3)));
			if(ListFirst(l)!=sp)
				continue;
			v1=(double)IntegerValue(CompoundArg1(CompoundArg1(ListFirst(l1))));
			v1/=(double)IntegerValue(CompoundArg2(CompoundArg1(ListFirst(l1))));
			if(CompoundArg2(ListFirst(l1)))
				{
				v1*=EvalParameter(CompoundArg2(ListFirst(l1)));
				for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
					{
					if(CompoundArg1(ListFirst(l2))==A_I)
						is_i=1;
					}
				}
			else
				nr++;
			if(is_i)
				{
				vi+=v1;
				ni++;
				}
			else
				{
				vr+=v1;
				nr++;
				}
			}
		if(nr && ni)
			{
			if(spi==1)
				fprintf(f,"(%+f%+f*i)",vr,vi);
			else
				fprintf(f,"(%+f**2%+f**2*i)",mysqrt(vr),mysqrt(vi));
			}
		else
			{
			if(nr)
			{
				if(spi==1)
					fprintf(f,"%+f",vr);
				else
					fprintf(f,"%+f**2",mysqrt(vr));
			}
			else
				fprintf(f,"%+f*i",vi);
			}
		if(ListFirst(l) && ListFirst(l)!=A_DELTA)
			{
			if(is_atom(ListFirst(l)))
				fprintf(f,"*%s",AtomValue(ListFirst(l)));
			else
				fprintf(f,"*p%ld",IntegerValue(ListFirst(l)));
			}
		}
	FreeAtomic(gl);
	}
*/

static void wrt_num(FILE *f, List ml, int spi)
	{
	List gl;
	List l;
	gl=NewList();
	for(l=ml;l;l=ListTail(l))
		{
		Term sp;
		if(is_empty_list(CompoundArgN(ListFirst(l),3)))
			sp=0;
		else
			sp=CompoundArg1(ListFirst(CompoundArgN(ListFirst(l),3)));
		if(!ListMember(gl,sp))
			gl=AppendLast(gl,sp);
		}
	
	for(l=gl;l;l=ListTail(l))
		{
		List l1;
		cmplx vc;
		vc.r=0.0;
		vc.i=0.0;
		
		for(l1=ml;l1;l1=ListTail(l1))
			{
			Term sp;
			double v1;
			List l2;
			cmplx vv;
			vv.r=vv.i=0.0;

			if(is_empty_list(CompoundArgN(ListFirst(l1),3)))
				sp=0;
			else
				sp=CompoundArg1(ListFirst(CompoundArgN(ListFirst(l1),3)));
			if(ListFirst(l)!=sp)
				continue;
			v1=(double)IntegerValue(CompoundArg1(CompoundArg1(ListFirst(l1))));
			v1/=(double)IntegerValue(CompoundArg2(CompoundArg1(ListFirst(l1))));
			if(CompoundArg2(ListFirst(l1)))
				{
				double vvr;
				vv=cEvalParameter(CompoundArg2(ListFirst(l1)));
				vvr=EvalParameter(CompoundArg2(ListFirst(l1)));
				vc.r+=v1*vv.r;
				vc.i+=v1*vv.i;
				}
			else
				vc.r+=v1;
			}
		if(vc.r!=0.0 && vc.i!=0.0)
			{
			if(spi==1)
				fprintf(f,"(%+f%+f*i)",vc.r,vc.i);
			else
				fprintf(f,"(%+f**2%+f**2*i)",mysqrt(vc.r),mysqrt(vc.i));
			}
		else
			{
			if(vc.r!=0.0)
			{
				if(spi==1)
					fprintf(f,"%+f",vc.r);
				else
					fprintf(f,"%+f**2",mysqrt(vc.r));
			}
			else
				fprintf(f,"%+f*i",vc.i);
			}
		if(ListFirst(l) && ListFirst(l)!=A_DELTA)
			{
			if(is_atom(ListFirst(l)))
				fprintf(f,"*%s",AtomValue(ListFirst(l)));
			else
				fprintf(f,"*p%ld",IntegerValue(ListFirst(l)));
			}
		}
	FreeAtomic(gl);
	}
	
static double chk_num(List ml, int spi, int *hc)
	{
	List gl;
	List l;
	gl=NewList();
	for(l=ml;l;l=ListTail(l))
		{
		Term sp;
		if(is_empty_list(CompoundArgN(ListFirst(l),3)))
			sp=0;
		else
			sp=CompoundArg1(ListFirst(CompoundArgN(ListFirst(l),3)));
		if(!ListMember(gl,sp))
			gl=AppendLast(gl,sp);
		}
	/*if(ListLength(gl)>1)
		{
		*hc=1;
		return 0.0;
		}
	*/
	for(l=gl;l;l=ListTail(l))
		{
		List l1;
		double vi,vr;
		int ni,nr;
		vi=0.0;
		vr=0.0;
		ni=0;
		nr=0;
		
		for(l1=ml;l1;l1=ListTail(l1))
			{
			Term sp;
			double v1;
			List l2;
			int is_i=0;
			if(is_empty_list(CompoundArgN(ListFirst(l1),3)))
				sp=0;
			else
				sp=CompoundArg1(ListFirst(CompoundArgN(ListFirst(l1),3)));
			if(ListFirst(l)!=sp)
				continue;
			v1=(double)IntegerValue(CompoundArg1(CompoundArg1(ListFirst(l1))));
			v1/=(double)IntegerValue(CompoundArg2(CompoundArg1(ListFirst(l1))));
			if(CompoundArg2(ListFirst(l1)))
				{
				v1*=EvalParameter(CompoundArg2(ListFirst(l1)));
				for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
					{
					if(CompoundArg1(ListFirst(l2))==A_I)
						is_i=1;
					}
				}
			else
				nr++;
			if(is_i)
				{
				vi+=v1;
				ni++;
				}
			else
				{
				vr+=v1;
				nr++;
				}
			}
		
		*hc=ni;
		if(spi==1)
			return vr;
		else
			return mysqrt(vr);
		}
	FreeAtomic(gl);
	return 0.0;
	}
	

static List prt_with_kt=0;
static List prt_with_mt=0;
static int correct_kt(Term m2, int spi);
static int correct_ktv(List m2);
static int correct_mt(Term t, Atom bn, Atom *m);

	
static void check_momenta(Term a2, Atom b1, Atom b2, FILE *f)
{
	List l1,l2,l3,lm;
	Term a2o;
/*	WriteTerm(a2);puts("");	*/
	l1=ConsumeCompoundArg(a2,4);
rpt1:
	for(l2=l1;l2;l2=ListTail(l2))
	{
		for(l3=CompoundArgN(ListFirst(l2),3);l3;l3=ListTail(l3))
		{
			if(CompoundName(ListFirst(l3))==A_MOMENT)
			{
				l1=CutFromList(l1,l2);
				goto rpt1;
			}
		}
	}
	SetCompoundArg(a2,4,l1);
	lm=NewList();
	l1=ConsumeCompoundArg(a2,5);
rpt2:
	for(l2=l1;l2;l2=ListTail(l2))
	{
		for(l3=CompoundArgN(ListFirst(l2),3);l3;l3=ListTail(l3))
		{
			if(CompoundName(ListFirst(l3))==A_MOMENT)
			{
				lm=AppendLast(lm,ListFirst(l2));
				ChangeList(l2,0);
				l1=CutFromList(l1,l2);
				goto rpt2;
			}
		}
	}
	SetCompoundArg(a2,5,l1);
	if(is_empty_list(lm))
		return;
	
	
/*	WriteTerm(lm); printf(" -> ");*/
	a2o=MakeCompound(A_ALG2,5);
	SetCompoundArg(a2o,1,MakeList1(A_I));
	SetCompoundArg(a2o,5,lm);
	alg2_red_orth(a2o);
	lm=ConsumeCompoundArg(a2o,5);
	FreeAtomic(a2o);
	
	{
		Term prp=GetAtomProperty(b2,PROP_TYPE);
		if(is_compound(prp) && CompoundName(prp)==OPR_PARTICLE &&
			(CompoundArgN(prp,7)==A_LEFT || CompoundArgN(prp,7)==A_RIGHT))
			{
			List ll1,ll2,ll3;
			for(ll1=lm;ll1;ll1=ListTail(ll1))
			{
				int n,d,g;
				Term m4=ListFirst(ll1);
				n=(int)IntegerValue(CompoundArg1(CompoundArg1(m4)));
				d=(int)IntegerValue(CompoundArg2(CompoundArg1(m4)));
				ll2=ConsumeCompoundArg(m4,3);
				for(ll3=ll2;ll3;ll3=ListTail(ll3))
					if(CompoundArg1(ListFirst(ll3))==A_GAMMAP ||
						CompoundArg1(ListFirst(ll3))==A_GAMMAM)
						{
						ll2=CutFromList(ll2,ll3);
						n*=2;
						break;
						}
				SetCompoundArg(m4,3,ll2);
				g=(int)gcf(n,d);
				SetCompoundArg(CompoundArg1(m4),1,NewInteger(n/g));
				SetCompoundArg(CompoundArg1(m4),2,NewInteger(d/g));
			}
			}
	}
	
	
        {
	Term prp=GetAtomProperty(b2,PROP_TYPE);
		if(is_compound(prp) && CompoundName(prp)==OPR_PARTICLE &&
			CompoundArgN(prp,4)==NewInteger(1) && ListLength(lm)==2)
		{
			List j1=CopyTerm(CompoundArgN(ListFirst(lm),3));
			List j2=CopyTerm(CompoundArgN(ListFirst(ListTail(lm)),3));
			for(l1=j1;l1;l1=ListTail(l1))
				if(CompoundArg1(ListFirst(l1))==A_GAMMAP ||
					CompoundArg1(ListFirst(l1))==A_GAMMAM)
					{
					j1=CutFromList(j1,l1);
					break;
					}
			for(l1=j2;l1;l1=ListTail(l1))
				if(CompoundArg1(ListFirst(l1))==A_GAMMAP ||
					CompoundArg1(ListFirst(l1))==A_GAMMAM)
					{
					j2=CutFromList(j2,l1);
					break;
					}		
			if(j1 && EqualTerms(j1,j2) && 
                            EqualTerms(CompoundArg2(ListFirst(lm)),
                                CompoundArg2(ListFirst(ListTail(lm)))))
			{
				int n,d,g;
				Term m4=ListFirst(lm);
				FreeAtomic(j2);
				CutFromList(lm,ListTail(lm));
				SetCompoundArg(ListFirst(lm),3,j1);
				n=(int)IntegerValue(CompoundArg1(CompoundArg1(m4)));
				d=(int)IntegerValue(CompoundArg2(CompoundArg1(m4)));
				n*=2;
				g=(int)gcf(n,d);
				SetCompoundArg(CompoundArg1(m4),1,NewInteger(n/g));
				SetCompoundArg(CompoundArg1(m4),2,NewInteger(d/g));
			}
		}
	}
				
	/*DumpList(lm); puts("");*/
        
	if(is_empty_list(lm))
		return;
	
/*	
	p1=CompoundArg1(ListFirst(CompoundArg1(a2)));
	p2=CompoundArg1(ListFirst(ListTail(CompoundArg1(a2))));
*/	
	
/*	WriteTerm(CompoundArg1(a2));printf(" base ");
	WriteTerm(b1);printf("/");WriteTerm(b2);printf("; sp=%d\n",prt_spin(b1));
	puts("");
*/

	
	
	if(b1!=b2 && GetAtomProperty(b1,A_ANTI2)!=b2)
	{
		fprintf(f,"\n----------------------------------------------------------------\n");
		fprintf(f,"Mixing of particles ");
		fWriteTerm(f,CompoundArg1(a2));
		fprintf(f," in kinetic term:\n");
		wrt_ml(f,lm);/*WriteTerm(lm);puts("");*/
		FreeAtomic(lm);
		return;
	}
	
	prt_with_kt=AppendLast(prt_with_kt,b1);
	
	if(ListLength(lm)==2 && prt_spin1(b1)==3 && correct_ktv(lm))
		return;
	
	if(ListLength(lm)>1 || !correct_kt(ListFirst(lm),prt_spin(b1)))
	{
		fprintf(f,"\n----------------------------------------------------------------\n");
		fprintf(f,"Incorrect kinetic term for particle ");
		fWriteTerm(f,CompoundArg1(a2));
		fprintf(f," :\n");
		wrt_ml(f,lm); /*DumpList(lm);puts("");*/
		FreeAtomic(lm);
		return;	
	}
	
}
	
static void look_for_orth(List prt, List al, int c, int s, FILE *f);
static void look_for_orth2(List prt, List al, int c, int s, FILE *f, Label,Label);

extern Label orth_del_offdiag;
extern Label orth_red_diag;
static int do_modify=0;
static int numcheck=0;
extern int FAOutput;
extern List all_vert_list2(void);

Term ProcCheckMasses(Term t, Term ind)
	{
	int fasav;
	List l, all_lagr;
	FILE *outf;
	List alone_prt=0, ds_prt=0;
	List masses=0;
	Atom p1, p2, b1, b2;

	if(FAOutput)
		{
		WarningInfo(11);
		puts("CheckMasses: not compatible with '-feynarts' switch.");
		unlink("masses.chk");
		return 0;
		}	
	do_modify=0;
	
	if(is_compound(t))
	{
		int i;
		for(i=1;i<=CompoundArity(t);i++)
		{
			Term a;
			a=CompoundArgN(t,i);
			if(is_atom(a) && strcmp(AtomValue(a),"modify")==0)
			{
				do_modify=1;
				continue;
			}
			if(is_atom(a) && strcmp(AtomValue(a),"numcheck")==0)
			{
				numcheck=1;
				continue;
			}
			if(is_compound(a) && CompoundName(a)==OPR_EQSIGN &&
					CompoundArity(a)==2 &&
					CompoundArg1(a)==NewAtom("DelOffDiag",0) &&
					is_atom(CompoundArg2(a)))
			{
				Atom v;
				Term prop;
				v=CompoundArg2(a);
				prop=GetAtomProperty(v,A_ORTH_MATR);
				if(!is_compound(prop))
				{
					WarningInfo(703);
					printf("CheckMasses: '%s' is not a matrix element.\n",
							AtomValue(v));
					continue;
				}
				orth_del_offdiag=CompoundArg1(prop);
				continue;
			}
			if(is_compound(a) && CompoundName(a)==OPR_EQSIGN &&
					CompoundArity(a)==2 &&
					CompoundArg1(a)==NewAtom("RedDiag",0) &&
					is_atom(CompoundArg2(a)))
			{
				Atom v;
				Term prop;
				v=CompoundArg2(a);
				prop=GetAtomProperty(v,A_ORTH_MATR);
				if(!is_compound(prop))
				{
					WarningInfo(703);
					printf("CheckMasses: '%s' is not a matrix element.\n",
							AtomValue(v));
					continue;
				}
				orth_red_diag=CompoundArg1(prop);
				continue;
			}
			ErrorInfo(709);
			printf("CheckMasses: invalid argument '");WriteTerm(a);puts("'.");
			return 0;
		}
	}
				
	
	outf=fopen("masses.chk","w");
	if(extlib_problems)
	{
		fprintf(outf,
		"WARNING: some functions from external libraries were not loaded.\n");
		fprintf(outf,"WARNING: Numerical values may be incorrect.\n");
	}
	fasav=FAOutput;
	/*FAOutput=0;*/
		
	all_lagr=all_vert_list2();
	
	for(l=all_lagr;l;l=ListTail(l))
		{
		Term a2;
		Term a2_so;
		Term prtl;
		List l1, l2;
		
		a2=ListFirst(l);
		prtl=CompoundArg1(a2);
		
				
		if(ListLength(prtl)>2 || is_empty_list(prtl))
			continue;

		{
			List ml1=ConsumeCompoundArg(a2,5);
			List ml2=NewList();
			for(l2=ml1;l2;l2=ListTail(l2))
			{
				List l3;
				for(l3=CompoundArg2(ListFirst(l2));l3;l3=ListTail(l3))
					if(GetAtomProperty(CompoundArg1(ListFirst(l3)),A_INFINITESIMAL) 
						&&GetAtomProperty(CompoundArg1(ListFirst(l3)),PROP_TYPE)==0 )
						break;
				if(!l3)
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

					
		a2=CopyTerm(a2);
		
		alg2_symmetrize(a2);
		alg2_common_s(a2);
		alg2_common_n(a2);
/*		alg2_common_t(a2);
		FreeAtomic(ConsumeCompoundArg(a2,4));
*/
		alg2_red_cos(a2);

		a2_so=CopyTerm(a2);
		alg2_red_orth(a2);

		alg2_red_sico(a2);
		alg2_red_sico(a2_so);

		alg2_red_comsico(a2);
		alg2_red_comsico(a2_so);
		
		alg2_recommon_n(a2);
		alg2_recommon_n(a2_so);
		
		alg2_recommon_s(a2);
		alg2_recommon_s(a2_so);
		
		alg2_red_1pm5(a2);
		alg2_red_1pm5(a2_so);
		
		if(is_empty_list(CompoundArgN(a2,5)))
			{
			FreeAtomic(a2);
			FreeAtomic(a2_so);
			continue;
			}
		
		alg2_decommon_s(a2);
		alg2_decommon_s(a2_so);
		
		alg2_decommon_n(a2);
		alg2_decommon_n(a2_so);

		
		if(ListLength(prtl)==1)
			{
			Atom p;
			prtl=ConsumeCompoundArg(a2,1);
			p=CompoundArg1(ListFirst(prtl));
			FreeAtomic(prtl);
			SetCompoundArg(a2,1,p);
			alone_prt=AppendFirst(alone_prt,a2);
			FreeAtomic(a2_so);
			continue;
			}
		
		SetCompoundArg(a2,4,ConsumeCompoundArg(a2_so,5));
		FreeAtomic(a2_so);
		
			
		prtl=ConsumeCompoundArg(a2,1);
		p1=CompoundArg1(ListFirst(prtl));
		p2=CompoundArg1(ListFirst(ListTail(prtl)));
		FreeAtomic(prtl);
		
		SetCompoundArg(a2,1,MakeCompound2(OPR_DIV,p1,p2));

		if(prt_spin(p1)!=prt_spin(p2))
			{
			ds_prt=AppendLast(ds_prt,a2);
			continue;
			}
		
		
		prtl=GetAtomProperty(p1,PROP_TYPE);
		if(CompoundName(prtl)==OPR_FIELD && CompoundArg2(prtl)==NewInteger(4))
			p1=CompoundArg1(prtl);
			
		prtl=GetAtomProperty(p2,PROP_TYPE);
		if(CompoundName(prtl)==OPR_FIELD && CompoundArg2(prtl)==NewInteger(4))
			p2=CompoundArg1(prtl);
			
		
		prtl=GetAtomProperty(p1,PROP_TYPE);
		if(CompoundName(prtl)==OPR_FIELD)
			{
			Atom bp;
			bp=CompoundArg1(prtl);
			prtl=GetAtomProperty(bp,PROP_TYPE);
			if(bp!=CompoundArg1(prtl))
				b1=GetAtomProperty(p1,A_ANTI);
			else
				b1=p1;
			}
		else
			b1=CompoundArg1(prtl);
		
		prtl=GetAtomProperty(p2,PROP_TYPE);
		if(CompoundName(prtl)==OPR_FIELD)
			{
			Atom bp;
			bp=CompoundArg1(prtl);
			prtl=GetAtomProperty(bp,PROP_TYPE);
			if(bp!=CompoundArg1(prtl))
				b2=GetAtomProperty(p2,A_ANTI);
			else
				b2=p2;
			}
		else
			b2=CompoundArg1(prtl);
		
		SetCompoundArg(a2,1,MakeCompound2(OPR_DIV,p1,p2));
		
		check_momenta(a2,b1,b2,outf);
		
		if(prt_spin(b1)<2)
		{
			for(l1=CompoundArgN(a2,5);l1;l1=ListTail(l1))
			{
				int n;
				n=(int)IntegerValue(CompoundArg1(CompoundArg1(ListFirst(l1))));
				n=-n;
				SetCompoundArg(CompoundArg1(ListFirst(l1)),1,NewInteger(n));
			}
			
			for(l1=CompoundArgN(a2,4);l1;l1=ListTail(l1))
			{
				int n;
				n=(int)IntegerValue(CompoundArg1(CompoundArg1(ListFirst(l1))));
				n=-n;
				SetCompoundArg(CompoundArg1(ListFirst(l1)),1,NewInteger(n));
			}
		}
		
		if(is_empty_list(CompoundArgN(a2,5)))
		{
			FreeAtomic(a2);
			continue;
		}
		
		if(!ListMember(prt_with_mt,b1))
			prt_with_mt=AppendLast(prt_with_mt,b1);
		if(!ListMember(prt_with_mt,b2))
			prt_with_mt=AppendLast(prt_with_mt,b2);
		
		
		for(l1=masses;l1;l1=ListTail(l1))
			{
			Term t;
			t=ListFirst(l1);

			if(ListMember(CompoundArg1(t),b1) || ListMember(CompoundArg1(t),b2))
				{
				if(!ListMember(CompoundArg1(t),b1))
					{
					l2=ConsumeCompoundArg(t,1);
					l2=AppendLast(l2,b1);
					SetCompoundArg(t,1,l2);
					}
				if(!ListMember(CompoundArg1(t),b2))
					{
					l2=ConsumeCompoundArg(t,1);
					l2=AppendLast(l2,b2);
					SetCompoundArg(t,1,l2);
					}

				l2=ConsumeCompoundArg(t,2);
				l2=AppendLast(l2,a2);
				SetCompoundArg(t,2,l2);
				break;
				}
			}
			
		if(is_empty_list(l1))
			masses=AppendLast(masses,MakeCompound2(OPR_DIV,
				b1==b2?MakeList1(b1):MakeList2(b1,b2),
				MakeList1(a2)));
		}
	
		
	FreeAtomic(all_lagr);
	FAOutput=fasav;

	/* end of all_vrt processing */
	
	{
		List l1,l2;
		for(l1=masses;l1;l1=ListTail(l1))
		{
		rpt:
			for(l2=ListTail(l1);l2;l2=ListTail(l2))
			{
				List l3;
				int ff=0;
				for(l3=CompoundArg1(ListFirst(l1));l3;l3=ListTail(l3))
					if(ListMember(CompoundArg1(ListFirst(l2)),ListFirst(l3)))
						ff=1;
				for(l3=CompoundArg1(ListFirst(l2));l3;l3=ListTail(l3))
					if(ListMember(CompoundArg1(ListFirst(l1)),ListFirst(l3)))
						ff=1;
				if(ff)
				{
					List lp,la;
					lp=ConsumeCompoundArg(ListFirst(l2),1);
					la=ConsumeCompoundArg(ListFirst(l2),2);
					CutFromList(masses,l2);
					l2=ConsumeCompoundArg(ListFirst(l1),1);
					for(l3=lp;l3;l3=ListTail(l3))
						if(!ListMember(l2,ListFirst(l3)))
							l2=AppendLast(l2,ListFirst(l3));
					SetCompoundArg(ListFirst(l1),1,l2);
					RemoveList(lp);
					l2=ConsumeCompoundArg(ListFirst(l1),2);
					l2=ConcatList(l2,la);
					SetCompoundArg(ListFirst(l1),2,l2);
					goto rpt;
				}
			}
		}
	}
	
	
	{
		List l1,l2;
		l2=NewList();
		for(l1=all_prtc_list();l1;l1=ListTail(l1))
		{
			Atom p;
			Term prop;
			p=ListFirst(l1);
			prop=GetAtomProperty(p,PROP_TYPE);
			if(CompoundName(prop)==OPR_FIELD)
				continue;
			if(p!=CompoundArg1(prop))
				continue;
			if(CompoundArgN(prop,7)==OPR_MLT)
				continue;
			if(!ListMember(prt_with_kt,ListFirst(l1)))
				l2=AppendLast(l2,ListFirst(l1));
		}
		if(l2)
		{
			fprintf(outf,"\n----------------------------------------------------------------\n");
			fprintf(outf,"Particles without kinetic term:\n");
			fWriteTerm(outf,l2);
			fprintf(outf,"\n");
			FreeAtomic(l2);
		}
	}
	
	{
		List l1,l2;
		l2=NewList();
		for(l1=all_prtc_list();l1;l1=ListTail(l1))
		{
			Atom p;
			Term prop;
			p=ListFirst(l1);
			prop=GetAtomProperty(p,PROP_TYPE);
			if(CompoundName(prop)==OPR_FIELD)
				continue;
			if(p!=CompoundArg1(prop))
				continue;
			if(CompoundArgN(prop,7)==OPR_MLT)
				continue;
			if(CompoundArgN(prop,5)==0)
				continue;
			if(!ListMember(prt_with_mt,ListFirst(l1)))
				l2=AppendLast(l2,ListFirst(l1));
		}
		if(l2)
		{
			fprintf(outf,"\n----------------------------------------------------------------\n");
			fprintf(outf,"Particles without mass term:\n");
			fWriteTerm(outf,l2);
			fprintf(outf,"\n");
			FreeAtomic(l2);
		}
	}
	
		
	for(l=masses;l;l=ListTail(l))
		if(ListLength(CompoundArg1(ListFirst(l)))>1)
			{
			List l1;
			l1=ConsumeCompoundArg(ListFirst(l),1);
			l1=SortedList(l1,acmp);
			SetCompoundArg(ListFirst(l),1,l1);
			}
			
	masses=SortedList(masses,mcmp);
		
	
	if(alone_prt)
		fprintf(outf,"\n----------------------------------------------------------------\n");
	for(l=alone_prt;l;l=ListTail(l))
		{
		fprintf(outf,"Linear term with particle %s\n",AtomValue(CompoundArg1(ListFirst(l))));
		fprintf(outf,"Expression: ");
		wrt_ml(outf,CompoundArgN(ListFirst(l),5));
		fprintf(outf,"\nNumeric value: ");
		wrt_num(outf,CompoundArgN(ListFirst(l),5),1);
		fprintf(outf,"\n\n");
		}
		
	if(ds_prt)
		fprintf(outf,"\n----------------------------------------------------------------\n");
	for(l=ds_prt;l;l=ListTail(l))
		{
		fprintf(outf,"Mixing particles with different spin: ");
		fWriteTerm(outf,CompoundArg1(ListFirst(l)));
		fprintf(outf,"\n");
		fprintf(outf,"Expression: ");
		wrt_ml(outf,CompoundArgN(ListFirst(l),5));
		/*WriteTerm(CompoundArgN(ListFirst(l),5));puts("");*/
		fprintf(outf,"\nNumeric value: ");
		wrt_num(outf,CompoundArgN(ListFirst(l),5),1);
		fprintf(outf,"\n\n");
		}
	
	for(l=masses;l;l=ListTail(l))
		{
		int is_c, spi;
		List l1,l2,al;
		int i1,i2; 
		
		l1=CompoundArg1(ListFirst(l));
		
		if(ListLength(l1)==1 || (ListLength(l1)==2 &&
			GetAtomProperty(ListFirst(l1),A_ANTI2)==ListNth(l1,2)) )
		{
			Atom tmass;
			Term a2;
			Atom bn;
			bn=ListFirst(l1);
			a2=ListFirst(CompoundArg2(ListFirst(l)));
			
			if(!correct_mt(a2,bn,&tmass))
			{
				double v;
				int hc=1;
				if(numcheck)
				{
					v=chk_num(CompoundArgN(a2,5),prt_spin(bn),&hc);
					if(hc==0 && is_atom(tmass) && 
						(fabs(EvalParameter(tmass)-v)/v)<1.0e-6)
						continue;
				}
				fprintf(outf,"\n---------------------------------------------------------------\n");
				fprintf(outf,"Mass term for particle ");
				fWriteTerm(outf,bn);
				fprintf(outf," :\n");
				wrt_ml(outf,CompoundArgN(a2,5));
				fprintf(outf," = ");
				wrt_num(outf,CompoundArgN(a2,5),prt_spin(bn));
				fprintf(outf,"\n");
				if(is_atom(tmass))
				{
					fprintf(outf,"This particles was declared with mass ");
					fWriteTerm(outf,tmass);
					fprintf(outf," = %f\n",EvalParameter(tmass));
				}
				else
				{
					if(tmass==0)
						fprintf(outf,"This particles was declared massless\n");
					else
						fprintf(outf,"This (derived) particle should not have mass term\n");
				}
			}
				
			continue;
		}
		
		if(numcheck)
		{
		int ok=1, hc=1;
		double v,v1;
		
		al=CompoundArg2(ListFirst(l));
		spi=prt_spin(ListFirst(l1));
		
		for(i1=1;i1<=ListLength(l1);i1++)
		for(i2=1;i2<=ListLength(l1);i2++)
			{
			if(i1>i2)
				continue;
				
			p1=GetAtomProperty(ListNth(l1,i1),A_ANTI);
			p2=ListNth(l1,i2);
			if(spi!=1 && (strcmp(AtomValue(p1),AtomValue(p2))>0))
				{
				Atom tmp;
				tmp=p1;
				p1=p2;
				p2=tmp;
				}
			v=0.0;
			for(l2=al;l2;l2=ListTail(l2))
				{
				if(CompoundArg1(CompoundArg1(ListFirst(l2)))==p1 &&
					CompoundArg2(CompoundArg1(ListFirst(l2)))==p2)
					{
					v=chk_num(CompoundArgN(ListFirst(l2),5),spi,&hc);
					break;
					}
				}
			if(hc) {ok=0; break;}
			v1=0.0;
			if(i1==i2)
				{
				Term ms=0,p;
				p=ListNth(l1,i1);
				p=GetAtomProperty(p,PROP_TYPE);
				if(p && CompoundName(p)==OPR_PARTICLE)
					ms=CompoundArgN(p,5);
				if(ms)
					v1=EvalParameter(ms);
				}
			if(v1==0.0 && fabs(v)<1.0e-4)
				continue;
			if(fabs(v1-v)/v1<1.0e-6)
				continue;
			ok=0;
			break;
			}
		
		if(ok)
			continue;
		
		}
		
		fprintf(outf,"\n---------------------------------------------------------------\n");
		fprintf(outf,"Mixing of particles ");
		fWriteTerm(outf,l1);
		fprintf(outf,"\n");
		
		is_c=(ListFirst(l1)!=GetAtomProperty(ListFirst(l1),A_ANTI));
		spi=prt_spin(ListFirst(l1));
		for(l2=l1;l2;l2=ListTail(l2))
			{
			int is_c1;
			is_c1=(ListFirst(l2)!=GetAtomProperty(ListFirst(l2),A_ANTI));
			if(is_c1 != is_c)
				{
				is_c=1;
				fprintf(outf,"Warning: charged and uncharged particles mixed\n");
				break;
				}
			}
			
		al=ConsumeCompoundArg(ListFirst(l),2);
		
		for(i1=1;i1<=ListLength(l1);i1++)
		for(i2=1;i2<=ListLength(l1);i2++)
			{
			if(is_c==0 && i1>i2)
				continue;
				
			p1=GetAtomProperty(ListNth(l1,i1),A_ANTI);
			p2=ListNth(l1,i2);
			if(spi!=1 && (strcmp(AtomValue(p1),AtomValue(p2))>0))
				{
				Atom tmp;
				tmp=p1;
				p1=p2;
				p2=tmp;
				}
			fprintf(outf,"\nM%d%d (%s/%s) = ",i1,i2,AtomValue(p1),AtomValue(p2));
			for(l2=al;l2;l2=ListTail(l2))
				{
				if(CompoundArg1(CompoundArg1(ListFirst(l2)))==p1 &&
					CompoundArg2(CompoundArg1(ListFirst(l2)))==p2)
					{
					wrt_ml(outf,CompoundArgN(ListFirst(l2),5));
/*					fprintf(outf,"\n    = ");
					wrt_ml(outf,CompoundArgN(ListFirst(l2),4));*/
					fprintf(outf,"\n    = ");
					wrt_num(outf,CompoundArgN(ListFirst(l2),5),spi);
					fprintf(outf,"\n");

					break;
					}
				}
			if(is_empty_list(l2))
				{
				fprintf(outf,"0\n");
				}
				
			if(i1==i2)
				{
				Term ms=0,p;
				p=ListNth(l1,i1);
				p=GetAtomProperty(p,PROP_TYPE);
				if(p && CompoundName(p)==OPR_PARTICLE)
					ms=CompoundArgN(p,5);
				if(ms)
					fprintf(outf,"    Declared mass %s=%f\n",AtomValue(ms),
							EvalParameter(ms));
				}
			fprintf(outf,"\n");	
			}
			
		fprintf(outf,"\n");
		
		look_for_orth(l1,al,is_c,spi,outf);
			
		}
		
        if(ftell(outf)==0)
            fprintf(outf,"No problems were found.\n");
	fclose(outf);
	puts("File 'masses.chk' is created"); 
	return 0;
	}

static List get_expr(List ml, int i, int j, int fixel, int is_c)
{
	List l1,l2,l3;
	List res=0;
	

	if(i>j)
	{
		int tmp;
		tmp=i;
		i=j;
		j=tmp;
	}
	
	for(l1=ml;l1;l1=ListTail(l1))
	{
		int ti,tj;
		List p1=0,p2=0;
		Term prop;
		for(l3=CompoundArg2(ListFirst(l1));l3;l3=ListTail(l3))
		{
			if((prop=GetAtomProperty(CompoundArg1(ListFirst(l3)),A_ORTH_MATR))&&
					CompoundArg1(prop)!=orth_red_diag)
			{
				if(p1==0)
					p1=l3;
				else
					p2=l3;
			}
		}
		if(p1==0)
			continue;
			
		prop=GetAtomProperty(CompoundArg1(ListFirst(p1)),A_ORTH_MATR);
		ti=(int)IntegerValue(CompoundArgN(prop,5-fixel));
		if(p2==0)
			tj=ti;
		else
		{
			prop=GetAtomProperty(CompoundArg1(ListFirst(p2)),A_ORTH_MATR);
			tj=(int)IntegerValue(CompoundArgN(prop,5-fixel));
		}
		
		if(tj<ti)
		{
			int tmp;
			tmp=ti;
			ti=tj;
			tj=tmp;
		}
		
		if(ti!=i || tj!=j)
		{
			/*printf("found %d %d; cont...\n",ti,tj);*/
			continue;
		}
		
		p1=CopyTerm(ListFirst(l1));
		
		l2=ConsumeCompoundArg(p1,2);

rpt:	for(l3=l2;l3;l3=ListTail(l3))
		{
			Term pp;
			if((pp=GetAtomProperty(CompoundArg1(ListFirst(l3)),A_ORTH_MATR))&&
					CompoundArg1(pp)!=orth_red_diag)
			{
				l2=CutFromList(l2,l3);
				goto rpt;
			}
		}
		SetCompoundArg(p1,2,l2);
		for(l2=res;l2;l2=ListTail(l2))
			if(EqualTerms(CompoundArg2(p1),CompoundArg2(ListFirst(l2))) &&
			   EqualTerms(CompoundArgN(p1,3),CompoundArgN(ListFirst(l2),3)))
				break;
		if(l2)
		{
			long int n,d,n1,d1,n2,d2,cf;
			n1=IntegerValue(CompoundArg1(CompoundArg1(p1)));
			d1=IntegerValue(CompoundArg2(CompoundArg1(p1)));
			n2=IntegerValue(CompoundArg1(CompoundArg1(ListFirst(l2))));
			d2=IntegerValue(CompoundArg2(CompoundArg1(ListFirst(l2))));
			n=n1*d2+n2*d1;
			d=d1*d2;
			cf=gcf(n,d);
			n/=cf;
			d/=cf;
			SetCompoundArg(CompoundArg1(ListFirst(l2)),1,NewInteger(n));
			SetCompoundArg(CompoundArg1(ListFirst(l2)),2,NewInteger(d));
			FreeAtomic(p1);
		}
		else
			res=AppendLast(res,p1);
	}
	
	for(l2=res;l2;l2=ListTail(l2))
	{
		long int n,d,cf;
		n=IntegerValue(CompoundArg1(CompoundArg1(ListFirst(l2))));
		d=IntegerValue(CompoundArg2(CompoundArg1(ListFirst(l2))));
/*		if(is_c==0)
			d*=2;*/
		if(i!=j)
			d*=2;
		cf=gcf(n,d);
		n/=cf;
		d/=cf;
		SetCompoundArg(CompoundArg1(ListFirst(l2)),1,NewInteger(n));
		SetCompoundArg(CompoundArg1(ListFirst(l2)),2,NewInteger(d));
	}
	/*printf("%d %d fixel=%d ",i,j,fixel);WriteTerm(res);puts("");*/

	return res;
}


static List get_expr2(List ml, int i, int j, int fixel, int is_c, Label ola, Label ola2)
{
	List l1,l2,l3;
	List res=0;
	
	
	for(l1=ml;l1;l1=ListTail(l1))
	{
		int ti,tj;
		List p1=0,p2=0;
		Term prop, prop2;
		for(l3=CompoundArg2(ListFirst(l1));l3;l3=ListTail(l3))
		{
			if((prop=GetAtomProperty(CompoundArg1(ListFirst(l3)),A_ORTH_MATR))&&
					CompoundArg1(prop)!=orth_red_diag)
			{
				if(p1==0)
					p1=l3;
				else
					p2=l3;
			}
		}
		if(p1==0 || p2==0)
		{
			printf("Error in determining chargino-like mass matrix: missing mixing element\n");
			return NewAtom("???",0);
		}
		
		prop =GetAtomProperty(CompoundArg1(ListFirst(p1)),A_ORTH_MATR);
		prop2=GetAtomProperty(CompoundArg1(ListFirst(p2)),A_ORTH_MATR);
		
		if(CompoundArg1(prop)==ola2)
			{
			Term tmp;
			tmp=prop;
			prop=prop2;
			prop2=tmp;
			}
		if(CompoundArg1(prop)!=ola || CompoundArg1(prop2)!=ola2)
		{
			printf("Error in determining chargino-like mass matrix: missing mixing element\n");
			return NewAtom("???",0);
		}
		
		ti=(int)IntegerValue(CompoundArgN(prop,5-fixel));
		tj=(int)IntegerValue(CompoundArgN(prop2,5-fixel));
		
		
		if(ti!=i || tj!=j)
		{
			/*printf("found %d %d; cont...\n",ti,tj);*/
			continue;
		}
		
		p1=CopyTerm(ListFirst(l1));
		
		l2=ConsumeCompoundArg(p1,2);

rpt:	for(l3=l2;l3;l3=ListTail(l3))
		{
			Term pp;
			if((pp=GetAtomProperty(CompoundArg1(ListFirst(l3)),A_ORTH_MATR))&&
					CompoundArg1(pp)!=orth_red_diag)
			{
				l2=CutFromList(l2,l3);
				goto rpt;
			}
		}
		SetCompoundArg(p1,2,l2);
		for(l2=res;l2;l2=ListTail(l2))
			if(EqualTerms(CompoundArg2(p1),CompoundArg2(ListFirst(l2))) &&
			   EqualTerms(CompoundArgN(p1,3),CompoundArgN(ListFirst(l2),3)))
				break;
		if(l2)
		{
			long int n,d,n1,d1,n2,d2,cf;
			n1=IntegerValue(CompoundArg1(CompoundArg1(p1)));
			d1=IntegerValue(CompoundArg2(CompoundArg1(p1)));
			n2=IntegerValue(CompoundArg1(CompoundArg1(ListFirst(l2))));
			d2=IntegerValue(CompoundArg2(CompoundArg1(ListFirst(l2))));
			n=n1*d2+n2*d1;
			d=d1*d2;
			cf=gcf(n,d);
			n/=cf;
			d/=cf;
			SetCompoundArg(CompoundArg1(ListFirst(l2)),1,NewInteger(n));
			SetCompoundArg(CompoundArg1(ListFirst(l2)),2,NewInteger(d));
			FreeAtomic(p1);
		}
		else
			res=AppendLast(res,p1);
	}
	
/*	for(l2=res;l2;l2=ListTail(l2))
	{
		int n,d,cf;
		n=IntegerValue(CompoundArg1(CompoundArg1(ListFirst(l2))));
		d=IntegerValue(CompoundArg2(CompoundArg1(ListFirst(l2))));
		if(i!=j)
			d*=2;
		cf=gcf(n,d);
		n/=cf;
		d/=cf;
		SetCompoundArg(CompoundArg1(ListFirst(l2)),1,NewInteger(n));
		SetCompoundArg(CompoundArg1(ListFirst(l2)),2,NewInteger(d));
	}*/
	/*printf("%d %d fixel=%d ",i,j,fixel);WriteTerm(res);puts("");*/
	return res;
}		
static void look_for_orth(List prt, List al, int is_c, int spi, FILE *f)
{
	List l1,l2;
	Label ola,ola2;
	int prtlen, fixel=0, i, j, nomix=0;
	
	Atom prta[32];
	int  prtp[32];
	Term mmatr[32][32];
	
	for(i=0;i<32;i++)
	for(j=0;j<32;j++)
		mmatr[i][j]=0;
	
/*
fprintf(f,"look for orth: ");fWriteTerm(f,prt);fprintf(f,"\n");
*/	
	fprintf(f,"Looking for the mixing matrix...\n");	
	ola=0;
	ola2=0;
	for(l1=al;l1;l1=ListTail(l1))
	{
		for(l2=CompoundArgN(ListFirst(l1),4);l2;l2=ListTail(l2))
		{
			int macoe=0,macoe2=0;
			List l3;
			for(l3=CompoundArg2(ListFirst(l2));l3;l3=ListTail(l3))
			{
				Term prop;
				Label th;
				prop=GetAtomProperty(CompoundArg1(ListFirst(l3)),A_ORTH_MATR);
				if(!prop)
					continue;
				if(CompoundArg1(prop)==orth_red_diag)
					continue;
				th=CompoundArg1(prop);
				if(th==ola)
					{
					macoe+=(int)IntegerValue(CompoundArg2(ListFirst(l3)));
					continue;
					}
				if(th==ola2)
					{
					macoe2+=(int)IntegerValue(CompoundArg2(ListFirst(l3)));
					continue;
					}
				if(ola==0)
					{
					ola=th;
					macoe+=(int)IntegerValue(CompoundArg2(ListFirst(l3)));
					continue;
					}
				if(ola2==0)
					{
					ola2=th;
					macoe2+=(int)IntegerValue(CompoundArg2(ListFirst(l3)));
					continue;
					}
				fprintf(f,"More than 2 mixing matrices found, vertex ");
				fWriteTerm(f,CompoundArg1(ListFirst(l1)));
				fprintf(f,"\n");
				return;
			}
			if(macoe+macoe2==0)
				{
				nomix++;
				continue;
				}
			if(macoe+macoe2!=2)
				{
				fprintf(f,"Expression with %d mixing elements, vertex ",
						macoe+macoe2);
				fWriteTerm(f,CompoundArg1(ListFirst(l1)));
				fprintf(f,"\n");
			    return;
				}
		}
	if(nomix)
		SetCompoundArg(ListFirst(l1),2,NewInteger(nomix));
	nomix=0;
	}
	
	if(ola==0)
		{
		fprintf(f,"No mixing matrices found.\n");
		return;
		}
	

	prtlen=ListLength(prt);
	if(is_c==0 && 2*ListLength(al)!=prtlen*(prtlen+1))
		fprintf(f,"%d particles, expecting %d vertices, found %d.\n",
			prtlen,prtlen*(prtlen+1)/2,ListLength(al));
			
	if(is_c!=0 &&   ListLength(al)!=prtlen*prtlen)
		fprintf(f,"%d particles, expecting %d vertices, found %d.\n",
			prtlen,prtlen*prtlen,ListLength(al));

	if(ola2)
	{
		look_for_orth2(prt, al, is_c, spi, f, ola, ola2);
		return;
	}


	for(i=1;i<=prtlen;i++)
	{
		prta[i]=ListNth(prt,i);
		prtp[i]=0;
	}
	
	for(i=1;i<=prtlen;i++)
	{
		Atom p1,p2;
		p2=prta[i];
		p1=GetAtomProperty(p2,A_ANTI);
		if(spi!=1 && (strcmp(AtomValue(p1),AtomValue(p2))>0))
				{
				Atom tmp;
				tmp=p1;
				p1=p2;
				p2=tmp;
				}
		for(l1=al;l1;l1=ListTail(l1))
		{
			if(CompoundArg1(CompoundArg1(ListFirst(l1)))==p1 &&
				CompoundArg2(CompoundArg1(ListFirst(l1)))==p2)
			{
				List oeused=0, e1used=0, e2used=0;
				List l3;
						
				for(l2=CompoundArgN(ListFirst(l1),4);l2;l2=ListTail(l2))
					for(l3=CompoundArg2(ListFirst(l2));l3;l3=ListTail(l3))
					{
						Atom p;
						Term prp1;
						p=CompoundArg1(ListFirst(l3));
						if((prp1=GetAtomProperty(p,A_ORTH_MATR)) && 
								CompoundArg1(prp1)!=orth_red_diag &&
								!ListMember(oeused,p))
							oeused=AppendLast(oeused,p);
					}
					
				for(l2=oeused;l2;l2=ListTail(l2))
				{
					Term prop;
					Integer i,j;
					prop=GetAtomProperty(ListFirst(l2),A_ORTH_MATR);
					i=CompoundArgN(prop,3);
					j=CompoundArgN(prop,4);
					if(!ListMember(e1used,i)) e1used=AppendLast(e1used,i);
					if(!ListMember(e2used,j)) e2used=AppendLast(e2used,j);
				}


				if((ListLength(e1used)==1&&ListLength(e2used)==1) ||
						(ListLength(e1used)!=1&&ListLength(e2used)!=1))
					return;
				
				if(ListLength(e1used)==1)
				{
					if(fixel==2) return;
					fixel=1;
					prtp[i]=(int)IntegerValue(ListFirst(e1used));
				}
				else
				{
					if(fixel==1) return;
					fixel=2;
					prtp[i]=(int)IntegerValue(ListFirst(e2used));
				}
				
				
				FreeAtomic(e1used);
				FreeAtomic(e2used);
				FreeAtomic(oeused);
				break;
			}
		}
		
	}
	
	
	for(l1=al;l1;l1=ListTail(l1))
	{
		if(is_integer(CompoundArg2(ListFirst(l1))))
				{
				fprintf(f,"Warning: terms without mixing in vertex ");
				fWriteTerm(f,CompoundArg1(ListFirst(l1)));
				fprintf(f,"\n");
				continue;
				}
		
		for(i=1;i<=prtlen;i++)
		for(j=i;j<=prtlen;j++)
		{
			Term thexpr;
			thexpr=get_expr(CompoundArgN(ListFirst(l1),4),
					prtp[i],prtp[j],fixel,is_c);
			if(mmatr[i][j]==0)
				mmatr[i][j]=thexpr;
			else
			{
				if(!EqualTerms(thexpr,mmatr[i][j]))
				{
					
					fprintf(f,"Warning: redefined M%d%d ",i,j);				
					fprintf(f," -> ");
					wrt_ml(f,thexpr);
					fprintf(f,"\n");
				}
				FreeAtomic(thexpr);
			}
		}
	}
	
	fprintf(f,"It is recognized that these fields are rotated by matrix:\n\n");
	
	for(i=1;i<=prtlen;i++)
	{
		fprintf(f,"  ( ");
		for(j=1;j<=prtlen;j++)
			fprintf(f," %s ",fixel==2?
			AtomValue(herm_matr_el(ola,prtp[i],prtp[j])):
			AtomValue(herm_matr_el(ola,prtp[j],prtp[i])));
		fprintf(f," )\n");
	}
	
	fprintf(f,"\nThe mass matrix before introducing this rotation:\n\n");
	
	for(i=1;i<=prtlen;i++)
	for(j=i;j<=prtlen;j++)
	{
		fprintf(f,"M%d%d = ",i,j);
		if(mmatr[i][j])
			wrt_ml(f,mmatr[i][j]);
		else
			fprintf(f,"0");
		fprintf(f,"\n");

	}
	
		
}

static void look_for_orth2(List prt, List al, int is_c, int spi, FILE *f,
				Label ola, Label ola2)
{
	List l1,l2;
	int prtlen, fixel=0, i, j;
	
	Atom prta[32];
	int  prtp[32];
	Term mmatr[32][32];
	
	for(i=0;i<32;i++)
	for(j=0;j<32;j++)
		mmatr[i][j]=0;
	
/*
fprintf(f,"look for orth: ");fWriteTerm(f,prt);fprintf(f,"\n");
*/	
	prtlen=ListLength(prt);

	for(i=1;i<=prtlen;i++)
	{
		prta[i]=ListNth(prt,i);
		prtp[i]=0;
	}
	
	for(i=1;i<=prtlen;i++)
	{
		Atom p1,p2;
		p2=prta[i];
		p1=GetAtomProperty(p2,A_ANTI);
		if(spi!=1 && (strcmp(AtomValue(p1),AtomValue(p2))>0))
				{
				Atom tmp;
				tmp=p1;
				p1=p2;
				p2=tmp;
				}
		for(l1=al;l1;l1=ListTail(l1))
		{
			if(CompoundArg1(CompoundArg1(ListFirst(l1)))==p1 &&
				CompoundArg2(CompoundArg1(ListFirst(l1)))==p2)
			{
				List oeused=0, e1used=0, e2used=0;
				List l3;
						
				for(l2=CompoundArgN(ListFirst(l1),4);l2;l2=ListTail(l2))
					for(l3=CompoundArg2(ListFirst(l2));l3;l3=ListTail(l3))
					{
						Atom p;
						Term prp1;
						p=CompoundArg1(ListFirst(l3));
						if((prp1=GetAtomProperty(p,A_ORTH_MATR)) && 
								CompoundArg1(prp1)!=orth_red_diag &&
								!ListMember(oeused,p))
							oeused=AppendLast(oeused,p);
					}
					
				for(l2=oeused;l2;l2=ListTail(l2))
				{
					Term prop;
					Integer i,j;
					prop=GetAtomProperty(ListFirst(l2),A_ORTH_MATR);
					i=CompoundArgN(prop,3);
					j=CompoundArgN(prop,4);
					if(!ListMember(e1used,i)) e1used=AppendLast(e1used,i);
					if(!ListMember(e2used,j)) e2used=AppendLast(e2used,j);
				}


				if((ListLength(e1used)==1&&ListLength(e2used)==1) ||
						(ListLength(e1used)!=1&&ListLength(e2used)!=1))
					goto fail1;
				
				if(ListLength(e1used)==1)
				{
					if(fixel==2) goto fail1;
					fixel=1;
					prtp[i]=(int)IntegerValue(ListFirst(e1used));
				}
				else
				{
					if(fixel==1) goto fail1;
					fixel=2;
					prtp[i]=(int)IntegerValue(ListFirst(e2used));
				}
				
				
				FreeAtomic(e1used);
				FreeAtomic(e2used);
				FreeAtomic(oeused);
				break;
			}
		}
		
	}
	
	
	for(l1=al;l1;l1=ListTail(l1))
	{
		Atom p1,p2;
		p1=CompoundArg1(CompoundArg1(ListFirst(l1)));
		p2=CompoundArg2(CompoundArg1(ListFirst(l1)));
		if(GetAtomProperty(p1,A_ANTI)!=p2)
			continue;
		for(i=1;i<=prtlen;i++)
		for(j=1;j<=prtlen;j++)
		{
			Term thexpr;
			thexpr=get_expr2(CompoundArgN(ListFirst(l1),4),
					prtp[i],prtp[j],fixel,is_c,ola,ola2);
			if(mmatr[i][j]==0)
				mmatr[i][j]=thexpr;
			else
			{
				if(!EqualTerms(thexpr,mmatr[i][j]))
				{
					
					fprintf(f,"Warning: redefined M%d%d ",i,j);				
					fWriteTerm(f,CompoundArg1(ListFirst(l1)));
					fprintf(f," -> ");
					wrt_ml(f,thexpr);
					fprintf(f,"\n");
					
					/*return;*/
				}
				FreeAtomic(thexpr);
			}
		}
	}
	
	fprintf(f,"Chargino-like mixing is recognized with rotation matrices:\n\n");
	
	for(i=1;i<=prtlen;i++)
	{
		fprintf(f,"  ( ");
		for(j=1;j<=prtlen;j++)
			fprintf(f," %s ",AtomValue(herm_matr_el(ola,prtp[i],prtp[j])));
		fprintf(f," )    ");
		fprintf(f,"  ( ");
		for(j=1;j<=prtlen;j++)
			fprintf(f," %s ",AtomValue(herm_matr_el(ola2,prtp[i],prtp[j])));
		fprintf(f," )\n");
	}
	
	fprintf(f,"\nThe mass matrix before introducing this rotation:\n\n");
	
	for(i=1;i<=prtlen;i++)
	for(j=1;j<=prtlen;j++)
	{
		fprintf(f,"M%d%d = ",i,j);
		if(mmatr[i][j])
			wrt_ml(f,mmatr[i][j]);
		else
			fprintf(f,"0");
		fprintf(f,"\n");

	}
	
	return;
	fail1:
	fprintf(f,"Mixing matrix can not be determined.\n");
	
	return;
		
}

static int correct_kt(Term m2, int spi)
{
	List l1,l2=0;
	Term ga=0;
	if(CompoundArg2(m2))
		return 0;
	if(CompoundArg2(CompoundArg1(m2))!=NewInteger(1))
		return 0;
	if((spi==0||spi==3) && CompoundArg1(CompoundArg1(m2))!=NewInteger(1))
		return 0;
	if((spi!=0&&spi!=3) && CompoundArg1(CompoundArg1(m2))!=NewInteger(-1))
		return 0;
	
	for(l1=CompoundArgN(m2,3);l1;l1=ListTail(l1))
	{
		Term t;
		t=ListFirst(l1);
		if(CompoundName(t)==OPR_SPECIAL && CompoundArg1(t)==A_DELTA)
			continue;
		if(CompoundName(t)==A_MOMENT)
		{
			l2=AppendLast(l2,t);
			continue;
		}
		if(CompoundName(t)==OPR_SPECIAL && CompoundArg1(t)==A_GAMMA && !ga)
		{
			ga=t;
			continue;
		}
		RemoveList(l2);
		return 0;
	}
	

	if(spi==1)
	{
		if(ListLength(l2)!=1 || ga==0 || 
				ListFirst(CompoundArg2(ListFirst(l2)))
				!= ListNth(CompoundArg2(ga),3))
		{
			RemoveList(l2);
			return 0;
		}
		RemoveList(l2);
		return 1;
	}
	
	if(spi==3)
	{
		if(ListLength(l2)==2)
		{
			RemoveList(l2);
			return 1;
		}
		RemoveList(l2);
		return 0;
	}
	
	if(ListLength(l2)!=2 || 
			ListFirst(CompoundArg2(ListFirst(l2)))!=
			ListFirst(CompoundArg2(ListFirst(ListTail(l2)))))
	{
		RemoveList(l2);
		return 0;
	}
	RemoveList(l2);
	return 1;
}


static int correct_ktv(List ml)
{
	if(correct_kt(ListFirst(ml),2) && correct_kt(ListFirst(ListTail(ml)),3))
		return 1;
	if(correct_kt(ListFirst(ml),3) && correct_kt(ListFirst(ListTail(ml)),2))
		return 1;
	return 0;
}
	
static int correct_mt(Term a2, Atom bn, Atom *m)
{
	int spi;
	Term prop;
	
	
	spi=prt_spin(bn);
	prop=GetAtomProperty(bn,PROP_TYPE);
	if(CompoundName(prop)==OPR_FIELD)
	{
		if(CompoundArg2(prop)!=NewInteger(1) &&
			CompoundArg2(prop)!=NewInteger(2) && 
			 CompoundArg2(prop)!=NewInteger(3))
		{
			*m=1;
			return 0;
		}
		prop=GetAtomProperty(CompoundArg1(prop),PROP_TYPE);
	}
	*m=CompoundArgN(prop,5);
	
	
	prop=CompoundArgN(a2,5);
	if(ListLength(prop)!=1)
		return 0;
	
	prop=ListFirst(prop);
	if(CompoundArg2(CompoundArg1(prop))!=NewInteger(1))
		return 0;
	if(CompoundArg1(CompoundArg1(prop))!=NewInteger(1))
		return 0;	
	if(ListLength(CompoundArg2(prop))!=1)
		return 0;
	if(CompoundArg1(ListFirst(CompoundArg2(prop)))!=(*m))
		return 0;
	if(spi==1 && CompoundArg2(ListFirst(CompoundArg2(prop)))!=NewInteger(1))
		return 0;
	if(spi!=1 && CompoundArg2(ListFirst(CompoundArg2(prop)))!=NewInteger(2))
		return 0;
	for(prop=CompoundArgN(prop,3);prop;prop=ListTail(prop))
		if(CompoundName(ListFirst(prop))!=OPR_SPECIAL ||
				CompoundArg1(ListFirst(prop))!=A_DELTA)
			return 0;
	return 1;
	
}


