#include <string.h>
#include "lanhep.h"

static List particles = 0;
static List pformat = 0, newprops=0;
extern int max_prt_lenp, max_prt_lenl, NoQuotes;

int pformat_def(void) { return pformat!=0;}

/*            1     2      3         4     5     6      7    8         */
/*   particle(Name, AName, FullName, Spin, Mass, Width, Aux, grlist ). */

int opAutoWidths=0;

static char *props[] = {
	"fullname",
	"name",
	"aname",
	"spin2",
	"mass",
	"width",
	"echarge",
	"echarge3",
	"color",
	"aux",
	"texname",
	"atexname"};
	
#define PNO (sizeof(props)/sizeof(char *))
	
	
static Atom mkprop(Atom a)
{
	char cbuf[128];
	sprintf(cbuf,"_%s",AtomValue(a));
	return NewAtom(cbuf,0);
}


extern int NoColors;

void register_particle(Term);

void AddIMParticle(Term prt)
	{
	register_particle(prt);
	particles=AppendLast(particles,prt);
	}

int is_particle(Atom p, List *ind)
	{
	Term prt;
	if(p==A_VEV)
		{
		if(ind!=NULL)
			*ind=NewList();
		return 1;
		}
    prt=GetAtomProperty(p,PROP_TYPE);
    if( !(is_compound(prt) && (CompoundName(prt)==OPR_PARTICLE ||
    	CompoundName(prt)==OPR_FIELD)))
        return 0;
    if(ind!=NULL)
        *ind=GetAtomProperty(p,PROP_INDEX);
    return prt;
	}


static char ganbuf[32];
Atom gen_anti_name(Atom name)
	{
	int i;
	char *s;
	i=0;
	s=AtomValue(name);
	do
		{
		ganbuf[i]=s[i];
		if(i==0)
			{
			if(s[0]>='a' && s[0]<='z')
				ganbuf[0]=s[0]-('a'-'A');
			else
				ganbuf[0]=s[0]+('a'-'A');
			}
		if(s[i]=='+') ganbuf[i]='-';
		if(s[i]=='-') ganbuf[i]='+';
		if(i!=0 && s[i]=='p') ganbuf[i]='m';
		if(i!=0 && s[i]=='m') ganbuf[i]='p';
		i++;
		} while(s[i]!=0);
    ganbuf[i]=0;
	return NewAtom(ganbuf,0);
	}			
	

int set_particle_name(Term t, Term sv, int pflag)
	{
	if(is_atom(t))
		{
		Atom aname;
		aname=gen_anti_name(t);
		SetCompoundArg(sv,1,t);	
		SetCompoundArg(sv,2,aname);
		return 1;
		}
	if(!is_compound(t))
		{
		ErrorInfo(110);
		printf(":  \'");
		WriteTerm(t);
		printf("\' illegal as particle name.\n");
		FreeAtomic(t);
		return 0;
		}
	if(FunctorName(CompoundFunctor(t))==OPR_DIV)
		{
		if(!is_atom(CompoundArg1(t)))
			{
			ErrorInfo(110);
			printf(" unexpected \'");
			WriteTerm(CompoundArg1(t));
			printf("\' illegal as particle name.\n");
			FreeAtomic(t);
			return 0;
			}
		if(!is_atom(CompoundArg2(t)))
			{
			ErrorInfo(110);
			printf(" \'");
			WriteTerm(CompoundArg2(t));
			printf("\' illegal as particle name.\n");
			FreeAtomic(t);
			return 0;
			}
		SetCompoundArg(sv,1,ConsumeCompoundArg(t,1));	
		SetCompoundArg(sv,2,ConsumeCompoundArg(t,2));
		FreeAtomic(t);
		return 1;
		}
	ErrorInfo(110);
	printf(" \'");
	WriteTerm(CompoundArg2(t));
	printf("\' illegal as particle name.\n");
	FreeAtomic(t);
	return 0;
	}
	
static char mnbuf[64];
	
static Atom make_param(Term p, char *type, Atom pname)
	{
	if(is_atom(p))
		{
		if(!is_parameter(p))
			{
			ErrorInfo(111);
				printf(" %s (%s of %s field) was not declared.\n",
					AtomValue(p),type,AtomValue(pname));
			}
		return p;
		}
	if(is_compound(p) && FunctorName(CompoundFunctor(p))==OPR_EQSIGN)
		{
		Term nm,val,comm,t1,t2,t3;
		nm=ConsumeCompoundArg(p,1);
		val=ConsumeCompoundArg(p,2);
		FreeAtomic(p);
		sprintf(mnbuf,"%s of %s",type,AtomValue(pname));
		comm=NewAtom(mnbuf,0);
		t1=NewCompound(NewFunctor(OPR_PARAMETER,1));
		t2=NewCompound(NewFunctor(OPR_EQSIGN,2));
		t3=NewCompound(NewFunctor(OPR_COLON,2));
		SetCompoundArg(t1,1,t2);
		SetCompoundArg(t2,1,nm);
		SetCompoundArg(t2,2,t3);
		SetCompoundArg(t3,1,val);
		SetCompoundArg(t3,2,comm);
		ProcessParameter(t1,0);
		return nm;
		}
	ErrorInfo(112);
	printf(" bad %s specification \'\n",type);
	WriteTerm(p);
	printf("\'.\n");
	return 0;
	}

static int set_particle_specs(Term t, Term sv, int pflag, int first)
	{
	Term repres;
	
	if(is_compound(t) && FunctorName(CompoundFunctor(t))==OPR_COMMA)
		{
		Term t1,t2;
		int ret=1;
		t1=ConsumeCompoundArg(t,1);
		t2=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		if(set_particle_specs(t1,sv,pflag,first)==0) ret=0;
		if(set_particle_specs(t2,sv,pflag,0)==0) ret=0;
		return ret;
		}
	
	if(t==A_LEFT || t==A_RIGHT || t==A_GAUGE || t==OPR_MLT)
		{
		SetCompoundArg(sv,7,t);
		return 1;
		}
		
	if(is_atom(t) && first==1)
		{
		SetCompoundArg(sv,3,t);
		return 1;
		}
		
	if(t==NewInteger(0) || 
		(is_atom(t) && AtomValue(t)[0]=='_' && AtomValue(t)[1]==0))
		return 1;
		
	if(is_compound(t) && strcmp(AtomValue(CompoundName(t)),"fullname")==0)
	{
		
		Atom n1;
		n1=CompoundArg1(t);
		if(!is_atom(n1))
			{
			ErrorInfo(113);
			printf(" bad specification \'");
			WriteTerm(t);
			printf("\'.\n");
			return 0;
			}
		SetCompoundArg(sv,3,n1);
		return 1;		
	}
	
	if(is_compound(t) && strcmp(AtomValue(CompoundName(t)),"atexname")==0)
	{
		Atom n1;
		n1=CompoundArg1(t);
		
		if(!is_atom(n1))
			{
			ErrorInfo(113);
			printf(" bad specification \'");
			WriteTerm(t);
			printf("\'.\n");
			return 0;
			}
		SetAtomProperty(CompoundArg2(sv),A_TEXNAME,n1);
		return 1;		
	}
	
	if(is_compound(t) && strcmp(AtomValue(CompoundName(t)),"aux")==0)
	{
		Atom n1;
		n1=CompoundArg1(t);
		
		if(!is_atom(n1) && !is_integer(n1))
			{
			ErrorInfo(113);
			printf(" bad specification \'");
			WriteTerm(t);
			printf("\'.\n");
			return 0;
			}
		SetAtomProperty(CompoundArg2(sv),A_AUX,n1);
		SetCompoundArg(sv,7,n1);
		return 1;		
	}
	
	if(is_compound(t) && strcmp(AtomValue(CompoundName(t)),"echarge")==0)
	{
		int n,d=1;
		Term t1, p1, p2;
		t1=CompoundArg1(t);
		if(is_integer(t1))
			n=(int)IntegerValue(t1);
		else
		{
			if(!is_compound(t1)||CompoundName(t1)!=OPR_DIV ||
					!is_integer(CompoundArg1(t1)) ||
					!is_integer(CompoundArg2(t1)))
				
			{
			ErrorInfo(113);
			printf(" bad specification \'");
			WriteTerm(t);
			printf("\'.\n");
			return 0;
			}
			n=(int)IntegerValue(CompoundArg1(t1));
			d=(int)IntegerValue(CompoundArg2(t1));
		}
		p1=CompoundArg1(sv);
		p2=CompoundArg2(sv);
		
		SetAtomProperty(p1,A_EM_CHARGE,
			MakeCompound2(OPR_DIV,NewInteger(n), NewInteger(d)));
		
		if(p1!=p2)
			SetAtomProperty(p2,A_EM_CHARGE,
				MakeCompound2(OPR_DIV,NewInteger(-n), NewInteger(d)));
		if(p1==p2 && n)
		{
			WarningInfo(900);
			printf(
				"particle %s identical to antiparticle with nonzero charge\n",
					AtomValue(p1));
		}
		return 1;
	}
		
	
	if(is_compound(t) && CompoundName(t)==A_TEXNAME)
		{
		Atom n1,n2;
		Term t1; 
		t1=CompoundArg1(t);
		if(is_compound(t1) && CompoundName(t1)==OPR_DIV)
			{
			n1=ConsumeCompoundArg(t1,1);
			n2=ConsumeCompoundArg(t1,2);
			FreeAtomic(t1);
			}
		else
			{
			n1=t1;
			n2=0;
			}
		if(!is_atom(n1) || (n2 && !is_atom(n2)))
			{
			ErrorInfo(113);
			printf(" bad specification \'");
			WriteTerm(t);
			printf("\'.\n");
			return 0;
			}
		SetAtomProperty(CompoundArg1(sv),A_TEXNAME,n1);
		if(n2 && CompoundArg1(sv)!=CompoundArg2(sv))
			SetAtomProperty(CompoundArg2(sv),A_TEXNAME,n2);
		return 1;
		}		
		
	if(is_compound(t) && FunctorName(CompoundFunctor(t))==OPR_MASS &&
		FunctorArity(CompoundFunctor(t))==1)
			{
			Term m,m1;
			m=ConsumeCompoundArg(t,1);
			
			if(m==NewInteger(0))
				return 1;
			FreeAtomic(t);
			if(CompoundArgN(sv,3)==0)
				m1=make_param(m,"mass",CompoundArg1(sv));
			else
				m1=make_param(m,"mass",CompoundArgN(sv,3));
			SetAtomProperty(m1,OPR_MASS,CompoundArg1(sv));
			SetCompoundArg(sv,5,m1);
			return 1;
			}
			
	if(is_compound(t) && FunctorName(CompoundFunctor(t))==OPR_WIDTH &&
		FunctorArity(CompoundFunctor(t))==1)
			{
			Term m,m1;
			m=ConsumeCompoundArg(t,1);
			FreeAtomic(t);
			if(is_compound(m)&&CompoundArity(m)==2&&is_atom(CompoundArg1(m))&&
				is_atom(CompoundArg2(m))&&
				strcmp(AtomValue(CompoundArg2(m)),"auto")==0)
			{
				if(UFOutput)
                {
                    SetCompoundArg(m,2,NewInteger(0));
                    if(CompoundArgN(sv,3)==0)
                    m1=make_param(m,"width",CompoundArg1(sv));
                    else
                    m1=make_param(m,"width",CompoundArgN(sv,3));
                }
                else {
                if(is_parameter(CompoundArg1(m)))
					{WarningInfo(0);
					 printf("%s is already declared parameter.\n",
					 	AtomValue(CompoundArg1(m)));
					}
				SetAtomProperty(CompoundArg1(m),OPR_WIDTH,A_I);
                    m1=CompoundArg1(m);}
			}
			else
			{
				if(CompoundArgN(sv,3)==0)
					m1=make_param(m,"width",CompoundArg1(sv));
				else
					m1=make_param(m,"width",CompoundArgN(sv,3));
			}
			SetCompoundArg(sv,6,m1);
			if(opAutoWidths && GetAtomProperty(m1,OPR_WIDTH)==0)
				SetAtomProperty(m1,OPR_WIDTH,CompoundArg1(sv));
			return 1;
			}
			
	if(is_compound(t) && ListMember(newprops,CompoundName(t)))
	{
		SetAtomProperty(CompoundArg1(sv),mkprop(CompoundName(t)),
				ConsumeCompoundArg(t,1));
		return 1;
	}
		
	
	repres=SpecToRepr(t);
	if(repres==0)
		{
		ErrorInfo(114);
		printf(" bad specification \'");
		WriteTerm(t);
		printf("\'.\n");
		return 0;
		}
		
	repres=AppendLast(ConsumeCompoundArg(sv,8),repres);
	SetCompoundArg(sv,8,repres);	
	return 1;
	}
	
	
static void proc_prtcl1(Term t, int pflag, int spin)
	{
	Term sv;
	
	if(is_compound(t) && FunctorName(CompoundFunctor(t))==OPR_COMMA)
		{
		Term t1,t2;
		t1=ConsumeCompoundArg(t,1);
		t2=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		proc_prtcl1(t1,pflag,spin);
		proc_prtcl1(t2,pflag,spin);
		return;
		}
		
	if(pflag)
		sv=NewCompound(NewFunctor(OPR_PARTICLE,8));
	else
		sv=NewCompound(NewFunctor(OPR_PARTICLE,8));

	if(is_atom(t))
		{
		set_particle_name(t,sv,pflag);
		SetCompoundArg(sv,4,NewInteger(spin));
		register_particle(sv);
		particles=AppendLast(particles,sv);
		return;
		}
		
	if(!is_compound(t))
		{
		ErrorInfo(115);
		printf(" unexpected \'");
		WriteTerm(t);
		printf("\' in particle declaration statement.\n");
		FreeAtomic(t);
		return;
		}
		
	if(FunctorName(CompoundFunctor(t))==OPR_DIV)
		{
		if(set_particle_name(t,sv,pflag)==0)
			return;
		SetCompoundArg(sv,4,NewInteger(spin));
        register_particle(sv);
		particles=AppendLast(particles,sv);
		return;
		}
		
	if(FunctorName(CompoundFunctor(t))!=OPR_COLON)
		{
		ErrorInfo(116);
		printf(" unexpected operator  \'%s\' in \'",
			AtomValue(FunctorName(CompoundFunctor(t))));
		WriteTerm(t);
		printf("\' in particle declaration statement.\n");
		FreeAtomic(t);
		return;
		}
		
	if(set_particle_name(ConsumeCompoundArg(t,1),sv,pflag)==0)
		return ;
	
	if(set_particle_specs(ConsumeCompoundArg(t,2),sv,pflag,1)==0)
		return ;
		
	FreeAtomic(t);
	SetCompoundArg(sv,4,NewInteger(spin));
    register_particle(sv);
	particles=AppendLast(particles,sv);
	return;
	} 
	
	
	
static void proc_prtcl(Term t, int pflag)
	{
	
	if(is_compound(t) && FunctorName(CompoundFunctor(t))==OPR_SECO)
		{
		Term t1,t2;
		t1=ConsumeCompoundArg(t,1);
		t2=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		proc_prtcl(t1,pflag);
		proc_prtcl(t2,pflag);
		return;
		}
	if(is_compound(t) && FunctorName(CompoundFunctor(t))==OPR_SCALAR &&
		FunctorArity(CompoundFunctor(t))==1)
			{
			Term t1;
			t1=ConsumeCompoundArg(t,1);
			FreeAtomic(t);
			proc_prtcl1(t1,pflag,0);
			return;
			}
	if(is_compound(t) && FunctorName(CompoundFunctor(t))==OPR_SPINOR &&
		FunctorArity(CompoundFunctor(t))==1)
			{
			Term t1;
			t1=ConsumeCompoundArg(t,1);
			FreeAtomic(t);
			proc_prtcl1(t1,pflag,1);
			return;
			}
	if(is_compound(t) && FunctorName(CompoundFunctor(t))==OPR_VECTOR &&
		FunctorArity(CompoundFunctor(t))==1)
			{
			Term t1;
			t1=ConsumeCompoundArg(t,1);
			FreeAtomic(t);
			proc_prtcl1(t1,pflag,2);
			return;
			}
	if(is_compound(t) && FunctorName(CompoundFunctor(t))==OPR_SPINOR3 &&
		FunctorArity(CompoundFunctor(t))==1)
			{
			Term t1;
			t1=ConsumeCompoundArg(t,1);
			FreeAtomic(t);
			proc_prtcl1(t1,pflag,3);
			return;
			}
	if(is_compound(t) && FunctorName(CompoundFunctor(t))==OPR_TENSOR &&
		FunctorArity(CompoundFunctor(t))==1)
			{
			Term t1;
			t1=ConsumeCompoundArg(t,1);
			FreeAtomic(t);
			proc_prtcl1(t1,pflag,4);
			return;
			}
	ErrorInfo(117);
	printf(" expected 'scalar', 'spinor', 'vector', got  '");
	WriteTerm(t);
	printf("\'.\n");
	FreeAtomic(t);

	}
	
/*   particle(Name, AName, FullName, Spin, Mass, Width, Aux, ..... ). */

static void ParticleList(int spin)
	{
	List li;
	li=particles;
	while(!is_empty_list(li))
		{
		Term t;
		t=ListFirst(li);
		if(spin==-1 || spin==IntegerValue(CompoundArgN(t,4)))
			{
			
			List l;
			if(spin==-1)
				switch(IntegerValue(CompoundArgN(t,4)))
					{
					case 0:
						printf("scalar\t"); break;
					case 1:
						printf("spinor\t"); break;
					case 2:
						printf("vector\t"); break;
					default:
						printf("\t");
					}
			printf("%s/%s\t",AtomValue(CompoundArg1(t)),
				AtomValue(CompoundArg2(t)));
			if(CompoundArgN(t,3)!=0)
				printf("\'%s\'   ",AtomValue(CompoundArgN(t,3)));
			printf("mass ");
			if(CompoundArgN(t,5)==0)
				printf("0");
			else
				WriteTerm(CompoundArgN(t,5));
			printf("   ");
			printf("  width "); 
			if(CompoundArgN(t,6)==0)
				printf("0");
			else
				WriteTerm(CompoundArgN(t,6));
			printf("   ");
			if(CompoundArgN(t,7)!=0)
				{printf("  "); WriteTerm(CompoundArgN(t,7)); printf("  "); }
			l=CompoundArgN(t,8);
			while(!is_empty_list(l))
				{
				WriteTerm(ListFirst(l));
				l=ListTail(l);
				}
			puts("");
			}
		li=ListTail(li);	
		}
	}
				
Term ProcessParticle(Term t, Term ind)
	{
	Term arg;
	char *s;
	int pflag;
	if(CompoundName(t)==OPR_PARTICLE)
		pflag=1;
	else
		pflag=0;
	arg=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	if(is_atom(arg))
		{
		s=AtomValue(arg);
		if(s[0]=='?' && s[1]==0)
			{
			ParticleList(-1);
			return 0;
			}
		}
	proc_prtcl(arg,pflag);
	return 0;
	}
	
	
	
Term ProcessScalar(Term t, Term ind)
	{
	Term arg;
	char *s;
	arg=CompoundArg1(t);
	if(is_atom(arg) && (s=AtomValue(arg)) && s[0]=='?')
		{
		ParticleList(0);
		FreeAtomic(t);
		return 0;
		}
	proc_prtcl(t,1);
	return 0;
	}
	
Term ProcessSpinor(Term t, Term ind)
	{
	Term arg;
	char *s;
	arg=CompoundArg1(t);
	if(is_atom(arg) && (s=AtomValue(arg)) && s[0]=='?')
		{
		ParticleList(1);
		FreeAtomic(t);
		return 0;
		}
	proc_prtcl(t,1);
	return 0;
	}
	
Term ProcessVector(Term t, Term ind)
	{
	Term arg;
	char *s;
	arg=CompoundArg1(t);
	if(is_atom(arg) && (s=AtomValue(arg)) && s[0]=='?')
		{
		ParticleList(2);
		FreeAtomic(t);
		return 0;
		}
	proc_prtcl(t,1);
	return 0;
	}

Term ProcessSpinor3(Term t, Term ind)
	{
	Term arg;
	char *s;
	arg=CompoundArg1(t);
	if(is_atom(arg) && (s=AtomValue(arg)) && s[0]=='?')
		{
		ParticleList(3);
		FreeAtomic(t);
		return 0;
		}
	proc_prtcl(t,1);
	return 0;
	}
	
Term ProcessTensor(Term t, Term ind)
	{
	Term arg;
	char *s;
	arg=CompoundArg1(t);
	if(is_atom(arg) && (s=AtomValue(arg)) && s[0]=='?')
		{
		ParticleList(4);
		FreeAtomic(t);
		return 0;
		}
	proc_prtcl(t,1);
	return 0;
	}


static char wpbuf[128];

static void texWritePrt(int fno)
	{
	FILE *f;
	List li;
	
	if(OutputDirectory!=NULL)
		sprintf(wpbuf,"%s/prtcls%d.tex",OutputDirectory,fno);
	else
		sprintf(wpbuf,"prtcls%d.tex",fno);
		
	f=fopen(wpbuf,"w");
	if(f==NULL)
		{
		printf("Can not open file \'%s\' for writing.\n",wpbuf);
		perror("");
		return;
		}
	
	if(defined_em_charge())
		{
		fprintf(f,"\\begin{tabular}{|cc|l|c|c|c|l|} \\hline\n");	
		fprintf(f,"P & aP & Name & Spin  & EM charge & Color & Comment \\\\ \\hline\n");	
		}
	else
		{	
		fprintf(f,"\\begin{tabular}{|cc|l|c|l|l|c|l|} \\hline\n");	
		fprintf(f,"P & aP & Name & Spin & Mass & Width & Color &  Comment \\\\ \\hline\n");
		}
	li=particles;
	while(!is_empty_list(li))
		{
		Term ttt;
		int sp;
		List l;
		int col;
		Term prp1, prp2;
		ttt=ListFirst(li);
		
		col=1;
		l=CompoundArgN(ttt,8);
		while(!is_empty_list(l))
			{
			Term t1;
			t1=ListFirst(l);
			if(CompoundArity(t1)==2 && CompoundName(t1)==A_COLOR)
				{
				if(EqualTerms(CompoundArg1(t1),CompoundArg2(t1)))
					col=8;
				else
					col=3;
				break;
				}
			l=ListTail(l);
			}
		
		prp1=GetAtomProperty(CompoundArg1(ttt),A_TEXNAME);
		prp2=GetAtomProperty(CompoundArg2(ttt),A_TEXNAME);
		
		if(prp1)
			sp=fprintf(f,"$%s{}_{",AtomValue(prp1));
		else
			sp=fprintf(f,"$%s_{",AtomValue(CompoundArg1(ttt)));
		if(CompoundArgN(ttt,4)==NewInteger(2))
			sp+=fprintf(f,"\\mu ");
		if(CompoundArgN(ttt,4)==NewInteger(1))
			sp+=fprintf(f,"a");
		if(CompoundArgN(ttt,4)==NewInteger(3))
			sp+=fprintf(f,"a\\mu ");
		if(CompoundArgN(ttt,4)==NewInteger(4))
			sp+=fprintf(f,"\\mu\\nu");
		if(col>1)
			sp+=fprintf(f,"p");
		sp+=fprintf(f,"}$");
			
		WriteBlank(f,10-sp);
		fprintf(f,"&");
		
		if(prp2)
			sp=fprintf(f,"$%s{}_{",AtomValue(prp2));
		else
			if(prp1)
				{
				if(CompoundArg1(ttt)==CompoundArg2(ttt))
					sp=fprintf(f,"$%s{}_{",AtomValue(prp1));
				else
					sp=fprintf(f,"$\\bar{%s}{}_{",AtomValue(prp1));
				}
			else
				sp=fprintf(f,"$%s_{",AtomValue(CompoundArg2(ttt)));
				
		if(CompoundArgN(ttt,4)==NewInteger(2))
			sp+=fprintf(f,"\\mu ");
		if(CompoundArgN(ttt,4)==NewInteger(1))
			sp+=fprintf(f,"a");
		if(CompoundArgN(ttt,4)==NewInteger(4))
			sp+=fprintf(f,"\\mu\\nu ");
		if(CompoundArgN(ttt,4)==NewInteger(3))
			sp+=fprintf(f,"a\\mu");
		if(col>1)
			sp+=fprintf(f,"q");
		sp+=fprintf(f,"}$");
			
		WriteBlank(f,10-sp);
		fprintf(f,"&");				
		
		if(CompoundArgN(ttt,3)!=0)
			sp=fprintf(f,"%s",AtomValue(CompoundArgN(ttt,3)));
		WriteBlank(f,13-sp);
		fprintf(f," &");
		
		if(IntegerValue(CompoundArgN(ttt,4))%2)
			sp=fprintf(f,"$%ld/2$",
					IntegerValue(CompoundArgN(ttt,4)));
		else
			sp=fprintf(f,"%ld",IntegerValue(CompoundArgN(ttt,4))/2);
		WriteBlank(f,12-sp);
		fprintf(f,"&");
		
		if(!defined_em_charge())
			{
			if(CompoundArgN(ttt,5)!=0)
				sp=fprintf(f,"%s",AtomValue(CompoundArgN(ttt,5)));
			else
				sp=fprintf(f,"0");
			WriteBlank(f,6-sp);
			fprintf(f,"&");
			if(CompoundArgN(ttt,6)!=0)
				sp=fprintf(f,"%s",AtomValue(CompoundArgN(ttt,6)));
			else
				sp=fprintf(f,"0");
			WriteBlank(f,6-sp);
			fprintf(f,"&");
			}
		
		if(defined_em_charge())
			{
			Term prp;
			prp=GetAtomProperty(CompoundArg1(ttt),A_EM_CHARGE);
			if(prp==0)
				fprintf(f," $\\phantom{-}0$ &");
			else				
				{
				int num, den;
				num=(int)IntegerValue(CompoundArg1(prp));
				den=(int)IntegerValue(CompoundArg2(prp));
				if(num<0)
					{
					fprintf(f,"$-");
					num*=-1;
					}
				else
					fprintf(f,"$\\phantom{-}");
				if(den==1)
					fprintf(f,"%d$ &",num);
				else
					fprintf(f,"\\frac{%d}{%d}$ &",num,den);
				}
			}
			
		sp=fprintf(f,"%d",col);
		WriteBlank(f,5-sp);
		
/*		fprintf(f,"&");
		Write2Vertex(f,ttt); */
		fprintf(f,"&");
		sp=0;
		if(CompoundArgN(ttt,7)==A_LEFT)
			sp=fprintf(f,"left");
		if(CompoundArgN(ttt,7)==A_RIGHT)
			sp=fprintf(f,"right");
		if(CompoundArgN(ttt,7)==A_GAUGE)
			sp=fprintf(f,"gauge");
		if(CompoundArgN(ttt,7)==OPR_MLT)
			sp=fprintf(f,"intermed. prt.");
		WriteBlank(f,3-sp);

		li=ListTail(li);
		if(li)
			fprintf(f,"\\\\\n");
		else
			fprintf(f,"\\\\ \\hline\n");
		}
	fprintf(f,"\\end{tabular}\n");
	fclose(f);
	}

	
void WriteParticles(int fno, char *name)
	{
	FILE *f;
	List li;
	if(TexOutput)
		{
		texWritePrt(fno);
		return;
		}
	if(FAOutput)
		return;
	if(OutputDirectory!=NULL)
		sprintf(wpbuf,"%s/prtcls%d.mdl",OutputDirectory,fno);
	else
		sprintf(wpbuf,"prtcls%d.mdl",fno);
	f=fopen(wpbuf,"w");
	if(f==NULL)
		{
		printf("Can not open file \'%s\' for writing.\n",wpbuf);
		perror("");
		return;
		}
		
/*printf("%d %d\n",max_prt_lenp, max_prt_lenl);		*/

	fprintf(f,"%s\n Particles \n",name);
	if(pformat==0)
	{
	fprintf(f," Full  name  | P | aP|2*spin| mass |width |color|aux|>      LaTeX(A)      <|>     LaTeX(A+)      <|\n");
	}
	else
	{
		for(li=pformat;li;li=ListTail(li))
			fprintf(f,"%s|",AtomValue(CompoundArg1(ListFirst(li))));
		fprintf(f,"\n");
	}
	
	if(pformat)
	{
		for(li=particles;li;li=ListTail(li))
		{
			List l;
			Term ttt=ListFirst(li);
			if(CompoundName(ttt)!=OPR_PARTICLE)
				continue;
			for(l=pformat;l;l=ListTail(l))
			{
				int w=strlen(AtomValue(CompoundArg1(ListFirst(l)))),sp;
				Atom tp=CompoundName(ListFirst(l));
				Term prp;
				int i;
				for(i=0;i<PNO;i++)
					if(strcmp(props[i],AtomValue(tp))==0)
						break;
				sp=0;
				switch(i)
				{
					case 0:
						if(CompoundArgN(ttt,3)!=0)
							sp=fprintf(f,"%s",AtomValue(CompoundArgN(ttt,3)));
						else
							sp=fprintf(f,"%s",AtomValue(CompoundArg1(ttt)));
						WriteBlank(f,w-sp);
						fprintf(f,"|");
						break;
						
					case 1:
						sp=fprintf(f,"%s",AtomValue(CompoundArg1(ttt)));
						WriteBlank(f,w-sp);
						fprintf(f,"|");
						break;
						
					case 2:
						sp=fprintf(f,"%s",AtomValue(CompoundArg2(ttt)));
						WriteBlank(f,w-sp);
						fprintf(f,"|");
						break;
					case 3:
						sp=fprintf(f,"%ld",IntegerValue(CompoundArgN(ttt,4)));
						WriteBlank(f,w-sp);
						fprintf(f,"|");
						break;
					case 4:
						if(CompoundArgN(ttt,5)!=0)
							sp=fprintf(f,"%s",AtomValue(CompoundArgN(ttt,5)));
						else
							sp=fprintf(f,"0");
						WriteBlank(f,w-sp);
						fprintf(f,"|");
						break;
					case 5:
						if(CompoundArgN(ttt,6)!=0)
							{
							Atom aw=CompoundArgN(ttt,6);
							if(opAutoWidths)
							{
								List ll;
								for(ll=all_param_list();ll;ll=ListTail(ll))
									if(CompoundName(ListFirst(ll))==aw)
										break;
								if(ll && (CompoundArg1(ListFirst(ll))==0 ||
									CompoundArg1(ListFirst(ll))==NewInteger(0)))
								sp=fprintf(f,"!%s",AtomValue(aw));
								else
								sp=fprintf(f,"%s",AtomValue(aw));
							}
							else
								{
								if(GetAtomProperty(aw,OPR_WIDTH)==A_I)
									sp=fprintf(f,"!%s",AtomValue(aw));
								else
									sp=fprintf(f,"%s",AtomValue(aw));
								}
							}
						else
							sp=fprintf(f,"0");
						WriteBlank(f,w-sp);
						fprintf(f,"|");
						break;
					case 6:
						prp=GetAtomProperty(CompoundArg1(ttt),A_EM_CHARGE);
						if(prp==0)
						{
							fprintf(f,"0");
							WriteBlank(f,w-1);
							fprintf(f,"|");
							break;
						}
						if(CompoundArg2(prp)==NewInteger(1))
						{
							sp=fprintf(f,"%ld",IntegerValue(CompoundArg1(prp)));
							WriteBlank(f,w-sp);
							fprintf(f,"|");
							break;
						}
						sp=fprintf(f,"%ld/%ld",IntegerValue(CompoundArg1(prp)),
								IntegerValue(CompoundArg2(prp)));
						WriteBlank(f,w-sp);
						fprintf(f,"|");
						break;
					case 7:
						prp=GetAtomProperty(CompoundArg1(ttt),A_EM_CHARGE);
						if(prp==0)
						{
							fprintf(f,"0");
							WriteBlank(f,w-1);
							fprintf(f,"|");
							break;
						}
						if(CompoundArg2(prp)==NewInteger(1))
						{
							sp=fprintf(f,"%ld",3*IntegerValue(CompoundArg1(prp)));
							WriteBlank(f,w-sp);
							fprintf(f,"|");
							break;
						}
						if(CompoundArg2(prp)==NewInteger(3))
						{
							sp=fprintf(f,"%ld",IntegerValue(CompoundArg1(prp)));
							WriteBlank(f,w-sp);
							fprintf(f,"|");
							break;
						}
						sp=fprintf(f,"%ld/%ld",3*IntegerValue(CompoundArg1(prp)),
								IntegerValue(CompoundArg2(prp)));
						WriteBlank(f,w-sp);
						fprintf(f,"|");
						printf("Warning: particle '%s' with charge not N/3\n",
								AtomValue(CompoundArg1(ttt)));
						break;
					case 8:
						
					{
					List l;
					int col=1;
					l=CompoundArgN(ttt,8);
					while(!is_empty_list(l))
						{
						Term t1;
						t1=ListFirst(l);
						if(CompoundArity(t1)==2 && CompoundName(t1)==A_COLOR)
							{
							char *v=AtomValue(CompoundArg1(t1));
							col=v[1]-'0';
							/*
							if(EqualTerms(CompoundArg1(t1),CompoundArg2(t1)))
								col=8;
							else
								col=3;
							*/
							break;
							}
						l=ListTail(l);
						}
					if(NoColors && !FAOutput)
						col=1;
					sp=fprintf(f,"%d",col);
					WriteBlank(f,w-sp);
					}
					fprintf(f,"|");
					break;
					
				case 9:				
					
					sp=0;
					if(CompoundArgN(ttt,7)==A_LEFT)
						sp=fprintf(f,"L");
					else
					if(CompoundArgN(ttt,7)==A_RIGHT)
						sp=fprintf(f,"R");
					else
					if(CompoundArgN(ttt,7)==A_GAUGE)
						sp=fprintf(f,"G");
					else
					if(CompoundArgN(ttt,7)==OPR_MLT)
						sp=fprintf(f,"*");
					else
					if(CompoundArgN(ttt,7)!=0)
					{
					  NoQuotes=1;
						sp=fWriteTerm(f,CompoundArgN(ttt,7));
					  NoQuotes=0;
					}

					WriteBlank(f,3-sp);

					fprintf(f,"|");
					break;
					
				case 10:
		
				{
					Atom a;
					a=GetAtomProperty(CompoundArg1(ttt),A_TEXNAME);
					if(a==0)
						a=CompoundArg1(ttt);
					if(AtomValue(a)[0]=='~')
						sp=fprintf(f,"(%s)",AtomValue(a)+1);
					else
						sp=fprintf(f,"%s",AtomValue(a));
					WriteBlank(f,w-sp);
					fprintf(f,"|");
				}
					break;
				case 11:
				{
					Atom a1,a2;
					a1=GetAtomProperty(CompoundArg2(ttt),A_TEXNAME);
					a2=GetAtomProperty(CompoundArg1(ttt),A_TEXNAME);
					if(a1==0 && a2==0)
						a1=CompoundArg2(ttt);
					if(a1)
					{
						if(AtomValue(a1)[0]=='~')
							sp=fprintf(f,"(%s)",AtomValue(a1)+1);
						else
							sp=fprintf(f,"%s",AtomValue(a1));
					}
					else sp=fprintf(f,"\\bar{%s}",AtomValue(a2));
					WriteBlank(f,w-sp);
					fprintf(f,"|");
				}
					break;
				case 12:
				{
					Term prop=GetAtomProperty(CompoundArg1(ttt),mkprop(tp));
					sp=0;
					if(prop==0)
						prop=CompoundArg2(ListFirst(l));
					if(prop)
						sp=fWriteTerm(f,prop);
					else
						if(CompoundArgN(ttt,7)!=OPR_MLT && 
						      AtomValue(CompoundArgN(ttt,7))[1]!='*'
                                                  && 
						      AtomValue(CompoundArgN(ttt,7))[0]!='*')
						printf(
					"Warning: property '%s' is not defined for particle '%s'\n",
						AtomValue(tp),AtomValue(CompoundArg1(ttt)));
					WriteBlank(f,w-sp);
					fprintf(f,"|");
				}
			}
			
			
			if(w-sp<0 )
			{
				printf("Warning: Particle table: Property %s for particle %s is too wide, make the title wider.\n",
						AtomValue(tp),AtomValue(CompoundArg1(ttt)));
			}
		}
		fprintf(f,"\n");
	}
	fclose(f);
	return;
	}
	
	
	
	
	li=particles;
	while(!is_empty_list(li))
		{
		Term ttt;
		int sp;
		ttt=ListFirst(li);
		if(FunctorName(CompoundFunctor(ttt))!=OPR_PARTICLE)
			goto eend;
		if(CompoundArgN(ttt,3)!=0)
			sp=fprintf(f,"%s",AtomValue(CompoundArgN(ttt,3)));
		else
			sp=fprintf(f,"%s",AtomValue(CompoundArg1(ttt)));
		WriteBlank(f,13-sp);
		fprintf(f,"|");
		sp=fprintf(f,"%s",AtomValue(CompoundArg1(ttt)));
		WriteBlank(f,3-sp);
		fprintf(f,"|");
		sp=fprintf(f,"%s",AtomValue(CompoundArg2(ttt)));
		WriteBlank(f,3-sp);
		fprintf(f,"|");
		sp=fprintf(f,"%ld",IntegerValue(CompoundArgN(ttt,4)));
		WriteBlank(f,6-sp);
		fprintf(f,"|");
		if(CompoundArgN(ttt,5)!=0)
			sp=fprintf(f,"%s",AtomValue(CompoundArgN(ttt,5)));
		else
			sp=fprintf(f,"0");
		WriteBlank(f,6-sp);
		fprintf(f,"|");
		if(CompoundArgN(ttt,6)!=0)
			{
			if(GetAtomProperty(CompoundArgN(ttt,6),OPR_WIDTH)==A_I)
				sp=fprintf(f,"!%s",AtomValue(CompoundArgN(ttt,6)));
			else
				sp=fprintf(f,"%s",AtomValue(CompoundArgN(ttt,6)));
			}
		else
			sp=fprintf(f,"0");
		WriteBlank(f,6-sp);
		fprintf(f,"|");
		{
		List l;
		int col=1;
		l=CompoundArgN(ttt,8);
		while(!is_empty_list(l))
			{
			Term t1;
			t1=ListFirst(l);
			if(CompoundArity(t1)==2 && CompoundName(t1)==A_COLOR)
				{
				char *v=AtomValue(CompoundArg1(t1));
				col=v[1]-'0';
				/*if(EqualTerms(CompoundArg1(t1),CompoundArg2(t1)))
					col=8;
				else
					col=3;*/
				break;
				}
			l=ListTail(l);
			}
		if(NoColors && !FAOutput)
			col=1;
		sp=fprintf(f,"%d",col);
		WriteBlank(f,5-sp);
		}
		fprintf(f,"|");
		sp=0;
		if(CompoundArgN(ttt,7)==A_LEFT)
			sp=fprintf(f,"L");
		else
		if(CompoundArgN(ttt,7)==A_RIGHT)
			sp=fprintf(f,"R");
		else
		if(CompoundArgN(ttt,7)==A_GAUGE)
			sp=fprintf(f,"G");
		else
		if(CompoundArgN(ttt,7)==OPR_MLT)
			sp=fprintf(f,"*");
		else
		if(CompoundArgN(ttt,7)!=0)
			sp=fprintf(f,"%s",AtomValue(CompoundArgN(ttt,7)));
		WriteBlank(f,3-sp);

		fprintf(f,"|");
		
		{
			Atom a;
			a=GetAtomProperty(CompoundArg1(ttt),A_TEXNAME);
			if(a==0)
				a=CompoundArg1(ttt);
			if(AtomValue(a)[0]=='~')
				sp=fprintf(f,"(%s)",AtomValue(a)+1);
			else
				sp=fprintf(f,"%s",AtomValue(a));
			WriteBlank(f,22-sp);
			fprintf(f,"|");
			a=GetAtomProperty(CompoundArg2(ttt),A_TEXNAME);
			if(a==0)
				a=CompoundArg2(ttt);
			if(AtomValue(a)[0]=='~')
				sp=fprintf(f,"(%s)",AtomValue(a)+1);
			else
				sp=fprintf(f,"%s",AtomValue(a));
			WriteBlank(f,22-sp);
			fprintf(f,"|\n");
		}

eend:	li=ListTail(li);
		}
	fclose(f);
	}


void ClearParticle(Atom p)
	{
	List l;
	Atom ap;
	l=particles;
	while(!is_empty_list(l))
		{
		if(CompoundArg1(ListFirst(l))==p)
			{
			ap=CompoundArg2(ListFirst(l));
			ChangeList(l,0);
			particles=CutFromList(particles,l);
			RemoveAtomProperty(p,PROP_TYPE);
			RemoveAtomProperty(p,PROP_INDEX);
			RemoveAtomProperty(p,A_GRASS);
			if(ap!=p)
				{
				RemoveAtomProperty(ap,PROP_TYPE);
				RemoveAtomProperty(ap,PROP_INDEX);
				RemoveAtomProperty(ap,A_ANTI);
				RemoveAtomProperty(ap,A_GRASS);
				}
			return;
			}
		if(CompoundArg2(ListFirst(l))==p)
			{
			ap=CompoundArg1(ListFirst(l));
			ChangeList(l,0);
			particles=CutFromList(particles,l);
			RemoveAtomProperty(p,PROP_TYPE);
			RemoveAtomProperty(p,PROP_INDEX);
			RemoveAtomProperty(p,A_ANTI);
			RemoveAtomProperty(p,A_GRASS);
			if(ap!=p)
				{
				RemoveAtomProperty(ap,PROP_TYPE);
				RemoveAtomProperty(ap,PROP_INDEX);
				RemoveAtomProperty(ap,A_GRASS);
				}
			return;
			}
		l=ListTail(l);
		}
	ErrorInfo(0);
	puts("internal error (cplu)"); 
	}
	
Term ProcGhost(Term t, Term ind)
	{
	
	Atom t1,pp;
	Term gp;
	if(!is_compound(t) || CompoundArity(t)!=1 || !is_atom(CompoundArg1(t)))
		{
		ErrorInfo(202);
		puts("wrong input in ghost (ccghost) call"); 
		return 0;
		}
	t1=CompoundName(t);
	pp=CompoundArg1(t);
	FreeAtomic(t);
	gp=GetAtomProperty(pp,A_GHOST);
	if(!is_compound(gp) || CompoundArity(gp)!=2 )
		{
		ErrorInfo(202);
		printf("can not find ghost for '%s'\n",AtomValue(pp));
		return 0;
		}
	if(AtomValue(t1)[0]=='g')
		pp=CompoundArg1(gp);
	else
		pp=CompoundArg2(gp);
	return pp;
	
	}

static char mk_cc_buf[32];

Term ProcUp(Term p, Term ii)
	{
	Term p1;
	Term prp;
	int has_cc=0;
	
	p1=ConsumeCompoundArg(p,1);
	FreeAtomic(p);
	p=p1;
	if(is_function(p,NULL))
		p=CallFunction(p,0);
	
	if(!is_atom(p))
		{
		ErrorInfo(204);
		puts("wrong input in 'up' call.");
		return 0;
		}
		
	prp=GetAtomProperty(p,PROP_TYPE);
	
	if(CompoundName(prp)==OPR_FIELD && CompoundArg2(prp)==NewInteger(4))
	{
		has_cc=1;
		p=CompoundArg1(prp);
		prp=GetAtomProperty(p,PROP_TYPE);
	}
	
	if(CompoundName(prp)!=OPR_PARTICLE || CompoundArgN(prp,4)!=NewInteger(1))
	{
		ErrorInfo(204);
		puts("wrong input in 'up' call.");
		return 0;
	}
	
	if(CompoundArg1(prp)==CompoundArg2(prp))
	{
		if(has_cc)
			sprintf(mk_cc_buf,"%s.d",AtomValue(p));
		else
			sprintf(mk_cc_buf,"%s.u",AtomValue(p));
	}
	else
	{
		if(has_cc)
			sprintf(mk_cc_buf,"%s.u",
					AtomValue(CompoundArg1(GetAtomProperty(p,A_ANTI))));
		else
			sprintf(mk_cc_buf,"%s.u",AtomValue(p));
	}	

	p=NewAtom(mk_cc_buf,0);
	if(!is_particle(p,NULL))
		{
		ErrorInfo(204);
		puts("wrong input in 'up' call.");
		return 0;
		}
	return p;
	}	
Term ProcDown(Term p, Term ii)
	{
	Term p1;
	Term prp;
	int has_cc=0;
	
	p1=ConsumeCompoundArg(p,1);
	FreeAtomic(p);
	p=p1;
	if(is_function(p,NULL))
		p=CallFunction(p,0);
	
	if(!is_atom(p))
		{
		ErrorInfo(204);
		puts("wrong input in 'down' call.");
		return 0;
		}
		
	prp=GetAtomProperty(p,PROP_TYPE);
	
	if(CompoundName(prp)==OPR_FIELD && CompoundArg2(prp)==NewInteger(4))
	{
		has_cc=1;
		p=CompoundArg1(prp);
		prp=GetAtomProperty(p,PROP_TYPE);
	}
	
	if(CompoundName(prp)!=OPR_PARTICLE || CompoundArgN(prp,4)!=NewInteger(1))
	{
		ErrorInfo(204);
		puts("wrong input in 'down' call.");
		return 0;
	}
	
	if(CompoundArg1(prp)==CompoundArg2(prp))
	{
		if(has_cc)
			sprintf(mk_cc_buf,"%s.u",AtomValue(p));
		else
			sprintf(mk_cc_buf,"%s.d",AtomValue(p));
	}
	else
	{
		if(has_cc)
			sprintf(mk_cc_buf,"%s.d",
					AtomValue(CompoundArg1(GetAtomProperty(p,A_ANTI))));
		else
			sprintf(mk_cc_buf,"%s.d",AtomValue(p));
	}	

	p=NewAtom(mk_cc_buf,0);
	if(!is_particle(p,NULL))
		{
		ErrorInfo(204);
		puts("wrong input in 'down' call.");
		return 0;
		}
	return p;
	}	

	

Term ProcPrtcFormat(Term t, Term ind)
{
	List l,fl;
	if(!is_compound(t))
	{
		ErrorInfo(205);
		puts("wrong syntax in 'PrtcFormat' statement");
		return 0;
	}
	
	if(pformat)
	{
		WarningInfo(89);
		printf("Particle table format is redefined\n");
		FreeAtomic(pformat);
		pformat=0;
	}
	
	fl=CommaToList(ConsumeCompoundArg(t,1));
	FreeAtomic(t);
	for(l=fl;l;l=ListTail(l))
	{
		Term t=ListFirst(l);
		Term n,tt,d=0;
		if(is_atom(t))
		{
			pformat=AppendLast(pformat,MakeCompound2(t,t,0));
			continue;
		}
		if(is_compound(t)&&CompoundArity(t)==2&&CompoundName(t)==OPR_EQSIGN
				&&is_atom(CompoundArg1(t)) && is_atomic(CompoundArg2(t)))
		{
			n=CompoundArg1(t);
			d=CompoundArg2(t);
			pformat=AppendLast(pformat,MakeCompound2(n,n,d));
			continue;
		}
		if(is_compound(t)&&CompoundArity(t)==2&&CompoundName(t)==OPR_COLON)
		{
			n=CompoundArg1(t);
			tt=CompoundArg2(t);
			if(!is_atom(n)||!is_atom(tt)) goto err;
			pformat=AppendLast(pformat,MakeCompound2(n,tt,0));
			continue;
		}
		if(is_compound(t)&&CompoundArity(t)==2&&CompoundName(t)==OPR_EQSIGN)
		{
			d=CompoundArg2(t);
			t=CompoundArg1(t);
			if(!is_compound(t)||CompoundArity(t)!=2||CompoundName(t)!=OPR_COLON)
				goto err;
			n=CompoundArg1(t);
			tt=CompoundArg2(t);
			if(!is_atom(n)||!is_atom(tt)) goto err;
			pformat=AppendLast(pformat,MakeCompound2(n,tt,d));
			continue;
		}
	err:
		ErrorInfo(999);
		printf("prtcformat: wrong specification ");
		WriteTerm(ListFirst(l));
		puts(".");
	}
	
	for(l=pformat;l;l=ListTail(l))
	{
		Atom n=CompoundName(ListFirst(l));
		int i;
		for(i=0;i<PNO;i++)
			if(strcmp(AtomValue(n),props[i])==0)
				break;
		if(i==PNO&&!ListMember(newprops,n))
		{
			newprops=AppendLast(newprops,n);
			SetAtomProperty(n,PROP_TYPE,GetAtomProperty(OPR_MASS,PROP_TYPE));
			SetOperator(n,OP_FX,800);
		}
	}
	
	return 0;
}

Term ProcPrtcProp(Term t, Term ind)
{
	
	List l; 
	if(!is_compound(t))
	{
		ErrorInfo(205);
		puts("wrong syntax in 'PrtcFormat' statement");
		return 0;
	}
	
	t=CommaToList(ConsumeCompoundArg(t,1));
	
	for(l=t;l;l=ListTail(l))
	{
		Term t1=ListFirst(l), l1;
		Atom a,aprp;
		int i;
		
		if(is_atom(t1))
		{
			if(ListMember(newprops,t1))
				continue;
			newprops=AppendLast(newprops,t1);
			SetOperator(t1,OP_FX,800);
			SetAtomProperty(t1,PROP_TYPE,GetAtomProperty(OPR_MASS,PROP_TYPE));
			continue;
		}
			
		if(!is_compound(t1)|| !is_atom(CompoundArg1(t1)))
		{
			ErrorInfo(999);
			printf("prtcproperty: wrong specification ");
			WriteTerm(ListFirst(l));
			puts(".");
			continue;
		}
		
		a=CompoundArg1(t1);
		t1=CommaToList(ConsumeCompoundArg(t1,2));
		for(i=0;i<PNO;i++)
			if(strcmp(AtomValue(a),props[i])==0)
				break;
		if(i!=PNO && i!=4 && i!=5 && i!=9)
		{
			ErrorInfo(999);
			printf("prtcproperty: %s can not be set in this way, use different name.\n",
					props[i]);
			continue;
		}
		
		if(!ListMember(newprops,a))
		{
			newprops=AppendLast(newprops,a);
			SetAtomProperty(a,PROP_TYPE,GetAtomProperty(OPR_MASS,PROP_TYPE));
			SetOperator(a,OP_FX,800);
		}
		
		aprp=mkprop(a);
		for(l1=t1;l1;l1=ListTail(l1))
		{
			Term t2=ListFirst(l1);
			if(!is_compound(t2) || CompoundName(t2)!=OPR_EQSIGN
					|| !is_atom(CompoundArg1(t2)))
			{
				ErrorInfo(89);
				printf("prtcproperty: wrong specification ");
				WriteTerm(t2);
				puts(".");
			}
			else
			{
			if(i==4||i==5||i==9)
				{
				Term prp=GetAtomProperty(CompoundArg1(t2),PROP_TYPE);
				if(!is_compound(prp)||CompoundName(prp)!=OPR_PARTICLE)
				{
					ErrorInfo(112);
					WriteTerm(CompoundArg1(t2));
					puts(" is not a particle.");
					continue;
				}
				if(i!=9 && !is_parameter(CompoundArg2(t2)))
				{
					ErrorInfo(113);
					WriteTerm(CompoundArg2(t2));
					puts(" is not a parameter.");
					continue;
				}
				if(i==4 || i==5)
					SetCompoundArg(prp,i+1,CompoundArg2(t2));
				else
					SetCompoundArg(prp,7,CompoundArg2(t2));
				}
			else
			SetAtomProperty(CompoundArg1(t2),aprp,
					ConsumeCompoundArg(t2,2));
			}
		}
		FreeAtomic(t1);
		
	}
	FreeAtomic(t);
	return 0;
}

Term ProcPrtcFunction(Term t, Term ind)
{
	Atom a,p;
	int i;
	Term ttt;
	if(!is_compound(t))
		return 0;
	a=CompoundName(t);
	p=CompoundArg1(t);
	if(!is_atom(p)||!is_particle(p,NULL))
			{
				ErrorInfo(89);
				WriteTerm(p);
				printf("is not a particle in ");
				WriteTerm(t);
				puts(".");
			}
	ttt=GetAtomProperty(p,PROP_TYPE);
	for(i=0;i<PNO;i++)
		if(strcmp(props[i],AtomValue(a))==0)
			break;
	switch(i)
	{
		
	case 4:
		if(CompoundName(ttt)==OPR_PARTICLE)
		{
			if(CompoundArgN(ttt,5)!=0)
				return CompoundArgN(ttt,5);
			else
				return NewInteger(0);
		}
		return 0;
		break;
	case 5:
		if(CompoundName(ttt)==OPR_PARTICLE)
		{
			if(CompoundArgN(ttt,6)!=0)
				return CompoundArgN(ttt,6);
			else
				return NewInteger(0);
		}
		return 0;
		break;
	case 6:
		ttt=GetAtomProperty(CompoundArg1(ttt),A_EM_CHARGE);
		if(ttt==0)
		{
			return 0;
		}
		return CopyTerm(ttt);
		break;
	case 12:
		ttt=GetAtomProperty(p,mkprop(a));
		if(ttt)
			return CopyTerm(ttt);
	default:
			
		ErrorInfo(26);
		WriteTerm(t);
		puts(" is not defined\n");
		
	}
	return 0;
}
