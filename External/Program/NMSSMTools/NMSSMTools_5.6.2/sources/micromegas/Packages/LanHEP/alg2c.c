#include <string.h>	
#include <math.h>
#include "lanhep.h"

extern int opTriHeu, opAbbrVrt, has_cmplx_params;

List abbr_coeffs=0, no_abbr_coeffs=0;

void alg2_abbr_coeff(Term);

static int mcmp_err=0;

static int prmcmp(Term p1, Term p2)
	{
	if(CompoundArg1(p1)==A_I)
		return -1;
	if(CompoundArg1(p2)==A_I)
		return 1;
	return strcmp(AtomValue(CompoundArg1(p1)),AtomValue(CompoundArg1(p2)));
	}

	static int mcmp(Term p1, Term p2)
{
	p1=CompoundArg2(p1);
	p2=CompoundArg2(p2);
	if(ListLength(p1)<ListLength(p2))
		return -1;
	if(ListLength(p1)>ListLength(p2))
		return 1;
	for(;p1;p1=ListTail(p1),p2=ListTail(p2))
	{
		int s=strcmp(AtomValue(CompoundArg1(ListFirst(p1))),
				AtomValue(CompoundArg1(ListFirst(p2))));
		if(s)
			return s;
		s=(int)IntegerValue(CompoundArg2(ListFirst(p1)))
				-(int)IntegerValue(CompoundArg2(ListFirst(p2)));
		if(s)
			return s;
	}
	puts("Internal error (abbr-mcmp0)");
	mcmp_err=1;
	return 0;
}

		
void alg2_add_ml(Term a1, List ml)
	{
	
	List ml1, l1,l2;
	
	ml1=ConsumeCompoundArg(a1,5);
	
	for(l1=ml;l1;l1=ListTail(l1))
	{
		Term m2;
		m2=ListFirst(l1);
		
		for(l2=ml1;l2;l2=ListTail(l2))
		{
		Term t;
		t=ListFirst(l2);
		if(EqualTerms(CompoundArg2(t),CompoundArg2(m2)) &&
			EqualTerms(CompoundArgN(t,3),CompoundArgN(m2,3)))
				{
				long int d1,n1,d2,n2,d,n,gc;
				n1=IntegerValue(CompoundArg1(CompoundArg1(t)));
				d1=IntegerValue(CompoundArg2(CompoundArg1(t)));
				n2=IntegerValue(CompoundArg1(CompoundArg1(m2)));
				d2=IntegerValue(CompoundArg2(CompoundArg1(m2)));
				d=d1*d2;
				n=n1*d2+n2*d1;
				if(n==0)
					{
					ml1=CutFromList(ml1,l2);
					}
				else
					{
					gc=gcf(n,d);
					d/=gc;
					n/=gc;
					SetCompoundArg(CompoundArg1(t),1,NewInteger(n));
					SetCompoundArg(CompoundArg1(t),2,NewInteger(d));
					}
				FreeAtomic(m2);
				break;
				}
		}
		
		if(is_empty_list(l2))
			ml1=AppendLast(ml1,m2);
	}
	
	RemoveList(ml);
	SetCompoundArg(a1,5,ml1);
}
			
int a2_mono_no=0;

static void alg2_add_11(Term a1, Term a2)
	{
	Term m2;
	List ml,l;
	m2=ConsumeCompoundArg(a2,5);
	FreeAtomic(a2);
	a2=ListFirst(m2);
	RemoveList(m2);
	m2=a2;
	ml=ConsumeCompoundArg(a1,5);
	l=ml;
	while(!is_empty_list(l))
		{
		Term t;
		t=ListFirst(l);
		if(EqualTerms(CompoundArg2(t),CompoundArg2(m2)) &&
			EqualTerms(CompoundArgN(t,3),CompoundArgN(m2,3)))
				{
				long int d1,n1,d2,n2,d,n,gc;
				n1=IntegerValue(CompoundArg1(CompoundArg1(t)));
				d1=IntegerValue(CompoundArg2(CompoundArg1(t)));
				n2=IntegerValue(CompoundArg1(CompoundArg1(m2)));
				d2=IntegerValue(CompoundArg2(CompoundArg1(m2)));
				d=d1*d2;
				n=n1*d2+n2*d1;
				if(n==0)
					{
					ml=CutFromList(ml,l);
					}
				else
					{
					gc=gcf(n,d);
					d/=gc;
					n/=gc;
					SetCompoundArg(CompoundArg1(t),1,NewInteger(n));
					SetCompoundArg(CompoundArg1(t),2,NewInteger(d));
					}
				FreeAtomic(m2);
				SetCompoundArg(a1,5,ml);
				return ;
				}
		l=ListTail(l);
		}
		
	ml=AppendLast(ml,m2);
	a2_mono_no++;
	SetCompoundArg(a1,5,ml);
	}
			
	

static List alg2_add_1(List l, Term a2)
	{
	List l1;
	Term pl;
	pl=CompoundArg1(a2);
	l1=l;
	while(!is_empty_list(l1))
		{
		Term t1;
		t1=ListFirst(l1);
		if(EqualTerms(CompoundArg1(t1),pl))
			{
			alg2_add_11(t1,a2);
			return l;
			}
		l1=ListTail(l1);
		}
	return AppendLast(l,a2);
	}


List alg2_add(List a1, List a2)
	{
	List a2s;
	int alen, acur;
	char buf[40];
	alen=ListLength(a2);
	acur=0;
	a2s=a2;
	while(!is_empty_list(a2))
		{
		acur++;
		sprintf(buf,"alg2_add: %d of %d",acur,alen);
		RegisterLine(buf);
		a1=alg2_add_1(a1,ListFirst(a2));
		a2=ListTail(a2);
		UnregisterLine();
		}
	RemoveList(a2s);
	return a1;	
	}	
	

static int h_value(List prt)
{
	List l;
	int ret=43;
	for(l=prt;l;l=ListTail(l))
	{
		char *p;
		int nn;
		p=AtomValue(CompoundArg1(ListFirst(l)));
		nn=p[0]*13+p[1]*107+p[2]*511;
		ret+=ret*nn;
	}
	return ret;
}
	
void alg2_hash_add(List *lagr, int hsz, List a2)
	{
	List a2s;
	int alen, acur;
	char buf[40];
		
	alen=ListLength(a2);
	acur=0;
	a2s=a2;
	while(!is_empty_list(a2))
		{
		int hv;
		acur++;
		sprintf(buf,"alg2_add: %d of %d",acur,alen);
		RegisterLine(buf);
		hv=h_value(CompoundArg1(ListFirst(a2)));
		hv%=hsz;
		if(hv<0)
			hv=-hv;
		
		lagr[hv]=alg2_add_1(lagr[hv],ListFirst(a2));
		a2=ListTail(a2);
		UnregisterLine();
		}
	RemoveList(a2s);
	return;
	}	
	

	
/************* alg2_eval_vrt *********************/
	
	
static Term ev_mlt(Term t, Term t1)
{
	Atom p;
	int pw;
	p=CompoundArg1(t1);
	pw=(int)IntegerValue(CompoundArg2(t1));
	if(pw>1 || pw<-1)
		p=MakeCompound2(OPR_POW,p,pw>0?NewInteger(pw):NewInteger(-pw));
	if(t==0 && pw<0)
		t=NewInteger(1);
	if(t==0 && pw>0)
		return p;
	if(pw>0)
		return MakeCompound2(OPR_MLT,t,p);
	else
		return MakeCompound2(OPR_DIV,t,p);
}

Term l_2_t(List l, int num, int den)
{
	List l1;
	Term res;
	
	res=0;
	if(num<0)
		num=-num;
	if(num!=1)
		res=NewInteger(num);
	
	if(l)
	{
		for(l1=l;l1;l1=ListTail(l1))
			if(IntegerValue(CompoundArg2(ListFirst(l1)))>0)
				break;

		if(l1==0)
			l1=l;

		res=ev_mlt(res,ListFirst(l1));
		l=CutFromList(l,l1);

		for(l1=l;l1;l1=ListTail(l1))
			res=ev_mlt(res,ListFirst(l1));

		FreeAtomic(l);
	}
	else
		if(res==0)
			res=NewInteger(1);
	
	if(den!=1)
		res=MakeCompound2(OPR_DIV,res,NewInteger(den));
	
	return res;
}

void proc_abbr_param(Term name, Term value);


static int var_no=0;

int eval_vrt_len = 0;
int eval_vrt_more = 0;
void alg2_eval_more(Term a2);
int doing_abbr=0;

void alg2_eval_vrt(Term a2)
{
	
	List l1,l2;
	
	Term ctf;
	
	
	opTriHeu=0;
	doing_abbr=1;
	
	alg2_abbr_coeff(a2);
	
	if(eval_vrt_len==0)
		return;
	
	l1=CompoundArg1(a2);
	
	/* Do not split 2-leg vertices */
	
	if(ListLength(l1)<3)
		return;
	
	/* Do not split non-scalar vertices (?) */
	/*
	for(l2=l1;l2;l2=ListTail(l2))
		if(CompoundName(CompoundArg2(ListFirst(l2)))!=OPR_SCALAR)
			return;
	*/
		
	l1=CompoundArgN(a2,5);
	if(ListLength(l1)<eval_vrt_len)
		return;
	

	/*
	for(l2=CompoundArgN(a2,3);l2;l2=ListTail(l2))
			if(CompoundArg1(ListFirst(l2))==A_GG)
				return;
	*/	
	
	for(l2=l1;l2;l2=ListTail(l2))
	{
		List l3;
		for(l3=CompoundArg2(ListFirst(l2));l3;l3=ListTail(l3))
			if(CompoundArg1(ListFirst(l3))==A_GG)
				return;
	}
		
	
	
	if(!CalcOutput || has_cmplx_params==0)
	{
	int hi1=0,hi2=0;
	List l3;
	/*WriteTerm(l1);puts("");*/
	for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
		if(CompoundArg1(ListFirst(l2))==A_I)
			hi1=1;
	for(l2=ListTail(l1);l2;l2=ListTail(l2))
	{
		hi2=0;
		for(l3=CompoundArg2(ListFirst(l2));l3;l3=ListTail(l3))
			if(CompoundArg1(ListFirst(l3))==A_I)
				hi2=1;
		if(hi2!=hi1)
		{
			for(l2=l1;l2;l2=ListTail(l2))
			{
				List l4;
				l3=ConsumeCompoundArg(ListFirst(l2),2);
				for(l4=l3;l4;l4=ListTail(l4))
					if(CompoundArg1(ListFirst(l4))==A_I)
					{
						List lt=ConsumeCompoundArg(ListFirst(l2),3);
						lt=AppendFirst(lt,ListFirst(l4));
						SetCompoundArg(ListFirst(l2),3,lt);
						ChangeList(l4,0);
						l3=CutFromList(l3,l4);
						break;
					}
				SetCompoundArg(ListFirst(l2),2,l3);
			}
			alg2_eval_more(a2);
			for(l2=CompoundArgN(a2,5);l2;l2=ListTail(l2))
			{
				List lt=ConsumeCompoundArg(ListFirst(l2),3);
				if(lt && CompoundArg1(ListFirst(lt))==A_I)
				{
					List lp=ConsumeCompoundArg(ListFirst(l2),2);
					lp=AppendFirst(lp,ListFirst(lt));
					SetCompoundArg(ListFirst(l2),2,lp);
					ChangeList(lt,0);
					lt=CutFromList(lt,lt);
				}
				SetCompoundArg(ListFirst(l2),3,lt);
			}
			
			return;
			
			
		}
	}
	}
	
	ctf=CopyTerm(CompoundArgN(ListFirst(l1),3));
		
	for(l2=l1;l2;l2=ListTail(l2))
	{
		if(!EqualTerms(CompoundArgN(ListFirst(l2),3),ctf))
		{
			FreeAtomic(ctf);
			if(eval_vrt_more)
				alg2_eval_more(a2);
			return;
		}
	}

	alg2_eval_more(a2);
	
	return;
	
}

/*static List evll=0;*/

static List evll[1713];
static int evll_inited=0;
static int b_value(List ml)
{
	List l1,l2;
	int ret=43;
	int no=1;
	for(l1=ml;l1;l1=ListTail(l1))
	for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
	{
		char *p;
		int nn;
		p=AtomValue(CompoundArg1(ListFirst(l2)));
		nn=(p[0]*511+p[1]*107+p[2]*59+p[3]*37)*no;
		no++;
		ret+=nn;
	}
	return ret%1713;
}

void alg2_eval_more(Term a2)
{
	
	List l0,l1,l2,newml=0;
	int num,den,cnt,i;
	Term t2,  ctf=0;
	char cbuf[32];
	Atom newprm;
	
	if(!evll_inited)
	{
		for(i=0;i<1713;i++)
			evll[i]=0;
		evll_inited=1;
	}
	
		
	l1=CompoundArgN(a2,5);
	if(ListLength(l1)<eval_vrt_more)
		return;
	
	for(l1=CompoundArgN(a2,5);l1;l1=ListTail(l1))
	{
		l2=ConsumeCompoundArg(ListFirst(l1),2);
		if(l2) l2=SortedList(l2,prmcmp);
		SetCompoundArg(ListFirst(l1),2,l2);
	}
	
	ctf=NewList();
		
	for(l1=CompoundArgN(a2,5);l1;l1=ListTail(l1))
	{
		Term tf;
		tf=CompoundArgN(ListFirst(l1),3);
		for(l2=ctf;l2;l2=ListTail(l2))
		{
			if(EqualTerms(CompoundArgN(ListFirst(ListFirst(l2)),3),tf))
			{
				AppendLast(ListFirst(l2),ListFirst(l1));
				break;
			}
		}
		if(l2==0)
			ctf=AppendLast(ctf,MakeList1(ListFirst(l1)));
	}
	
	num=0;
	for(l1=ctf;l1;l1=ListTail(l1))
		num+=ListLength(ListFirst(l1))-1;
	
	
	if(num<eval_vrt_more)
	{
		if(ctf)
			RemoveList(ctf);
		return;
	}
	
		
	RemoveList(ConsumeCompoundArg(a2,5));
	for(l1=ctf;l1;l1=ListTail(l1))
		for(l2=ListTail(ListFirst(l1));l2;l2=ListTail(l2))
			FreeAtomic(ConsumeCompoundArg(ListFirst(l2),3));
	
	for(l0=ctf;l0;l0=ListTail(l0))
	{
		Term tf;
		List l3,l1c;
		int numf=1;
		
		l1=ListFirst(l0);
		if(ListLength(l1)==1)
		{
			newml=AppendLast(newml,ListFirst(l1));
			RemoveList(l1);
			continue;
		}
		
		tf=ConsumeCompoundArg(ListFirst(l1),3);
		
		l1=SortedList(l1,mcmp);
		if(IntegerValue(CompoundArg1(ListFirst(l1)))<0)
		{
			numf*=-1;
			for(l3=l1;l3;l3=ListTail(l3))
			{
				int nn=(int)IntegerValue(CompoundArg1(ListFirst(l3)));
				SetCompoundArg(ListFirst(l3),1,NewInteger(-nn));
			}
		}
		
		
		num=den=1;
		
		i=b_value(l1);
		if(i<0) i=-i;
		for(l3=evll[i];l3;l3=ListTail(l3))
			if(EqualTerms(l1,CompoundArg2(ListFirst(l3))))
		{
			newprm=CompoundArg1(ListFirst(l3));
			if(is_compound(newprm))
			{
				numf*=-1;
				newprm=CompoundArg1(newprm);
			}
			FreeAtomic(l1);
			goto addprm;
		}
		
		l1c=CopyTerm(l1);
		
		t2=0;
		cnt=0;
		for(l2=l1;l2;l2=ListTail(l2))
		{
			Term t3;
			cnt++;

			num=(int)IntegerValue(CompoundArg1(ListFirst(l2)));
			t3=l_2_t(ConsumeCompoundArg(ListFirst(l2),2),num,1);
			if(t2==0)
			{
				t2=t3;
				if(num<0)
					t2=MakeCompound1(OPR_MINUS,t2);
			}
			else
			{
				t2=MakeCompound2(num>0?OPR_PLUS:OPR_MINUS,t2,t3);
			}

			if((cnt-5*(abbr_coeffs!=0))==5 || ListTail(l2)==0)
			{
				sprintf(cbuf,"B%05d",var_no++);
				newprm=NewAtom(cbuf,0);
				/*newprm=add_parameter(newprm,t2);*/
				proc_abbr_param(newprm,t2);

				t2=newprm;
				cnt=0;
			}

		}
		
		
		FreeAtomic(l1);
		evll[i]=AppendLast(evll[i],MakeCompound2(OPR_EQSIGN,newprm,l1c));
		if(GetAtomProperty(newprm,A_ANTI))
		{
			List l4,l5, l1cc=CopyTerm(l1c);
			int j;
			Term ccprm=GetAtomProperty(newprm,A_ANTI);
			
			for(l4=l1cc;l4;l4=ListTail(l4))
				for(l5=CompoundArg2(ListFirst(l4));l5;l5=ListTail(l5))
				{
					Atom p,ap;
					p=CompoundArg1(ListFirst(l5));
					if(p==A_I)
					{
						int f=(int)IntegerValue(CompoundArg1(ListFirst(l4)));
						SetCompoundArg(ListFirst(l4),1,NewInteger(-f));
						continue;
					}
					ap=GetAtomProperty(p,A_ANTI);
					if(ap)
						SetCompoundArg(ListFirst(l5),1,ap);
				}
				
			for(l4=l1cc;l4;l4=ListTail(l4))
			{
				l5=ConsumeCompoundArg(ListFirst(l4),2);
				if(l5) l5=SortedList(l5,prmcmp);
				SetCompoundArg(ListFirst(l4),2,l5);
			}
			
			l1cc=SortedList(l1cc,mcmp);	
			if(IntegerValue(CompoundArg1(ListFirst(l1cc)))<0)
			{
				for(l4=l1cc;l4;l4=ListTail(l4))
				{
					int nn=(int)IntegerValue(CompoundArg1(ListFirst(l4)));
					SetCompoundArg(ListFirst(l4),1,NewInteger(-nn));
				}
				ccprm=MakeCompound1(OPR_MINUS,ccprm);
			}
			j=b_value(l1cc);
			evll[j]=AppendLast(evll[j],MakeCompound2(OPR_EQSIGN,
					ccprm,l1cc));
			
		}

addprm:
		
		t2=MakeCompound(A_MTERM,3);
		SetCompoundArg(t2,1,NewInteger(numf));
		/*sprintf(cbuf,"B%05d",var_no++);
		newprm=NewAtom(cbuf,0);*/
		/*newprm=add_parameter(newprm,t1);*/
		SetCompoundArg(t2,2,
					MakeList1(MakeCompound2(OPR_POW,newprm,NewInteger(1))));
		SetCompoundArg(t2,3,tf);
		newml=AppendLast(newml,t2);

/*		t2=MakeCompound1(OPR_PARAMETER,MakeCompound2(OPR_EQSIGN,newprm,t1));
		CallFunction(t2,0);*/
	}
	
/*	DumpList(newml);*/
	RemoveList(ctf);
	SetCompoundArg(a2,5,newml);
}



static List prmli[1713];
static int prmli_inited =0;

static int abbr_no=0, abbr_no2=0;

void abbr_stat(void)
	{

	int i,no,max,ave;

	
	no=0;max=0;ave=0;
	if(evll_inited)
	{
	for(i=0;i<1713;i++)
	{
		int len;
		if(evll[i]==0)
			no++;
		else
		{
			len=ListLength(evll[i]);
			if(len>max) max=len;
			ave+=len;
		}
	}
	printf("evl_hash: %d zero, avelen %d, maxlen %d\n",no,ave/(1713-no),max);
	}
	
	no=0;max=0;ave=0;
	if(prmli_inited)
	{
	for(i=0;i<1713;i++)
	{
		int len;
		if(prmli[i]==0)
			no++;
		else
		{
			len=ListLength(prmli[i]);
			if(len>max) max=len;
			ave+=len;
		}
	}
	printf("prm_hash: %d zero, avelen %d, maxlen %d\n",no,ave/(1713-no),max);
	}
	
	
	}

int r_abbr_no=0, c_abbr_no=0;
extern int UFOutput;

static Atom abbr_find(List pl)
{
	List l, csf=0, plsv;
	Atom a,b;
	Term t;
	char cbuf[16];
	int i;
	
	if(prmli_inited==0)
	{
		for(i=0;i<1713;i++)
			prmli[i]=0;
		prmli_inited=1;
	}
		
	i=b_value(pl);
	
	for(l=prmli[i];l;l=ListTail(l))
	{
		if(EqualTerms(CompoundArg2(ListFirst(l)),pl))
		{
			int fit;
			t=GetAtomProperty(CompoundArg1(ListFirst(l)),OPR_PLUS);
			if(t)
				fit=(int)IntegerValue(t);
			else
				fit=0;
			SetAtomProperty(CompoundArg1(ListFirst(l)),OPR_PLUS,
				NewInteger(fit+1));
			FreeAtomic(pl);
			return CompoundArg1(ListFirst(l));
		}
	}
	sprintf(cbuf,"A%05d",++abbr_no);
	a=NewAtom(cbuf,0);
	prmli[i]=AppendLast(prmli[i],MakeCompound2(OPR_EQSIGN,a,(plsv=CopyTerm(pl))));
	
	if(UFOutput)
	{
	List l1,l2;
	int eo=0, co=0;
	Atom ni=NewAtom("complex(0,1)",0);
	SetAtomProperty(ni,PROP_TYPE,OPR_PARAMETER);
	for(l1=pl;l1;l1=ListTail(l1))
	for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
		{
		if(CompoundArg1(ListFirst(l2))==A_EE)
			{
			int d=(int)IntegerValue(CompoundArg2(ListFirst(l2)));
			if(d>eo) eo=d;
			}
		if(CompoundArg1(ListFirst(l2))==A_GG)
			{
			int d=(int)IntegerValue(CompoundArg2(ListFirst(l2)));
			if(d>co) co=d;
			}
		if(CompoundArg1(ListFirst(l2))==A_I)
			{
			SetCompoundArg(ListFirst(l2),1,ni);
			}
		}
	if(eo)
		SetAtomProperty(a,A_EE,NewInteger(eo));
		
	if(co)
		SetAtomProperty(a,A_GG,NewInteger(co));
	}	
			
	if(ListLength(pl)>1)
	{
		List l1,l2,csf1=CopyTerm(CompoundArg2(ListFirst(pl)));
		for(l=csf1;l;l=ListTail(l))
		{
			for(l1=pl;l1;l1=ListTail(l1))
			{
				for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
					if(EqualTerms(ListFirst(l),ListFirst(l2)))
						break;
				if(l2==0)
					break;
			}
			if(l1==0)
				csf=AppendLast(csf,CopyTerm(ListFirst(l)));  
		}
		FreeAtomic(csf1);
		if(csf)
		{
			for(l=pl;l;l=ListTail(l))
			{
				csf1=ConsumeCompoundArg(ListFirst(l),2);
				for(l1=csf;l1;l1=ListTail(l1))
					for(l2=csf1;l2;l2=ListTail(l2))
						if(EqualTerms(ListFirst(l1),ListFirst(l2)))
						{
							csf1=CutFromList(csf1,l2);
							break;
						}
				SetCompoundArg(ListFirst(l),2,csf1);
			}
			csf=l_2_t(csf,1,1);
		}
	}
		
	if(1 || opAbbrVrt==1 || ListLength(pl)<=opAbbrVrt)
	{
		Atom aa;
		if(UFOutput==0)
		{
			t=l_2_t(ConsumeCompoundArg(ListFirst(pl),2),(int)IntegerValue(CompoundArg1(ListFirst(pl))),1);
			if(IntegerValue(CompoundArg1(ListFirst(pl)))<0)
				t=MakeCompound1(OPR_MINUS,t);
			for(l=ListTail(pl);l;l=ListTail(l))
				t=MakeCompound2(IntegerValue(CompoundArg1(ListFirst(l)))>0?OPR_PLUS:OPR_MINUS,
						t,l_2_t(ConsumeCompoundArg(ListFirst(l),2),
						(int)IntegerValue(CompoundArg1(ListFirst(l))),1));
		}
		else
		{
			Term r;
			r=ConsumeCompoundArg(ListFirst(pl),1);
			t=l_2_t(ConsumeCompoundArg(ListFirst(pl),2),(int)IntegerValue(CompoundArg1(r)),
														(int)IntegerValue(CompoundArg2(r)));
			if(IntegerValue(CompoundArg1(r))<0)
				t=MakeCompound1(OPR_MINUS,t);
			for(l=ListTail(pl);l;l=ListTail(l))
			{
				r=ConsumeCompoundArg(ListFirst(l),1);
				t=MakeCompound2(IntegerValue(CompoundArg1(r))>0?OPR_PLUS:OPR_MINUS,
						t,l_2_t(ConsumeCompoundArg(ListFirst(l),2),
						(int)IntegerValue(CompoundArg1(r)),(int)IntegerValue(CompoundArg2(r))));
			}
		}
		FreeAtomic(pl);
		if(csf)
			t=MakeCompound2(OPR_MLT,CopyTerm(csf),t);
		proc_abbr_param(a,t);
		/*t=MakeCompound2(OPR_EQSIGN,a,t);
		t=MakeCompound1(OPR_PARAMETER,t);
		CallFunction(t,0);*/
		if(csf) FreeAtomic(csf);
		if((aa=GetAtomProperty(a,A_ANTI)))
			{
			List apl=CopyTerm(plsv), l1;
			for(l=apl;l;l=ListTail(l))
			for(l1=CompoundArg2(ListFirst(l));l1;l1=ListTail(l1))
				{
				Atom prm=CompoundArg1(ListFirst(l1));
				Atom apr=GetAtomProperty(prm,A_ANTI);
				if(apr)
					SetCompoundArg(ListFirst(l1),1,apr);
				}
			i=b_value(apl);
			prmli[i]=AppendLast(prmli[i],MakeCompound2(OPR_EQSIGN,aa,apl));
			c_abbr_no++;
			SetAtomProperty(a,A_CHNAME,NewInteger(c_abbr_no));
			SetAtomProperty(aa,A_CHNAME,NewInteger(c_abbr_no));
			}
		else
			{
			r_abbr_no++;
			SetAtomProperty(a,A_CHNAME,NewInteger(r_abbr_no));
			}
		return a;
	}
	
	t=l_2_t(ConsumeCompoundArg(ListFirst(pl),2),(int)IntegerValue(CompoundArg1(ListFirst(pl))),1);
	if(IntegerValue(CompoundArg1(ListFirst(pl)))<0)
		t=MakeCompound1(OPR_MINUS,t);
	pl=CutFromList(pl,pl);
	
	
	while(ListLength(pl)>opAbbrVrt)
	{
		int i;
		sprintf(cbuf,"B%05d",++abbr_no2);
		b=NewAtom(cbuf,0);
		for(i=0;i<opAbbrVrt;i++)
		{
			t=MakeCompound2(IntegerValue(CompoundArg1(ListFirst(pl)))>0?OPR_PLUS:OPR_MINUS,
					t,l_2_t(ConsumeCompoundArg(ListFirst(pl),2),
					(int)IntegerValue(CompoundArg1(ListFirst(pl))),1));
			pl=CutFromList(pl,pl);
		}
		
		t=MakeCompound2(OPR_EQSIGN,b,t);
		t=MakeCompound1(OPR_PARAMETER,t);
		CallFunction(t,0);
		t=b;
	}
	
	for(l=pl;l;l=ListTail(l))
		t=MakeCompound2(IntegerValue(CompoundArg1(ListFirst(l)))>0?OPR_PLUS:OPR_MINUS,
				t,l_2_t(ConsumeCompoundArg(ListFirst(l),2),
				(int)IntegerValue(CompoundArg1(ListFirst(l))),1));
	FreeAtomic(pl);
	if(csf)
		t=MakeCompound2(OPR_MLT,CopyTerm(csf),t);
	t=MakeCompound2(OPR_EQSIGN,a,t);
	t=MakeCompound1(OPR_PARAMETER,t);
	CallFunction(t,0);
	if(csf) FreeAtomic(csf); 
	return a;
	
}


void alg2_abbr_vrt(Term a2)
{
	List l1,l2;
	int has_i=0, has_inf=0;
	List nml;
	
	if(no_abbr_coeffs)
	for(l1=CompoundArgN(a2,5);l1;l1=ListTail(l1))
	{
	  Term m=ListFirst(l1);
	  for(l2=CompoundArg2(m);l2;l2=ListTail(l2))
	  {
	    if(ListMember(no_abbr_coeffs,CompoundArg1(ListFirst(l2))))
	      return;
	  }
	}
	
	if(doing_abbr==0)
	{
		doing_abbr=1;
	}
	
	alg2_abbr_coeff(a2);
	
	l1=CompoundArgN(a2,5);
	if(is_empty_list(l1))
		return;
	

	if(UFOutput==1)
	{
		alg2_decommon_n(a2);
		SetCompoundArg(a2,2,MakeCompound2(OPR_DIV,NewInteger(1),NewInteger(1)));
	}
	

	if(ListLength(l1)==1 && UFOutput!=1)
	{
		l2=CompoundArg2(ListFirst(l1));
		if(l2==0 || (ListLength(l2)==1 && CompoundArg2(ListFirst(l2))==
				NewInteger(1)) || (ListLength(l2)==1 && 
				CompoundArg1(ListFirst(l2))==A_I &&
				CompoundArg2(ListFirst(ListTail(l2)))==NewInteger(1)))
			return;
	}
	
	if(UFOutput==0)
	for(;l1;l1=ListTail(l1))
	{
		List l3;
		Term m2=ListFirst(l1);
		for(l3=CompoundArg2(m2);l3;l3=ListTail(l3))
		if(CompoundArg1(ListFirst(l3))==A_I)
		{
			l2=ConsumeCompoundArg(m2,2);
			l2=CutFromList(l2,l3);
			SetCompoundArg(m2,2,l2);
			l2=ConsumeCompoundArg(m2,3);
			l2=AppendFirst(l2,MakeCompound2(OPR_SPECIAL,A_I,0));
			SetCompoundArg(m2,3,l2);
			has_i=1;
			break;
		}
	}
	
	for(l1=CompoundArgN(a2,5);l1;l1=ListTail(l1))
	{
		Term m2=ListFirst(l1);
		List l2=CompoundArg2(m2);
		for(;l2;l2=ListTail(l2))
		if(GetAtomProperty(CompoundArg1(ListFirst(l2)),
				A_INFINITESIMAL) && !GetAtomProperty(CompoundArg1(ListFirst(l2)),
				PROP_TYPE))
		{
			Atom prm=CompoundArg1(ListFirst(l2));
			List lp=ConsumeCompoundArg(m2,2);
			lp=CutFromList(lp,l2);
			SetCompoundArg(m2,2,lp);
			l2=ConsumeCompoundArg(m2,3);
			l2=AppendFirst(l2,MakeCompound2(OPR_PARAMETER,prm,0));
			SetCompoundArg(m2,3,l2);
			has_inf=1;
			break;
		}
	}
	
	l1=ConsumeCompoundArg(a2,5);
	nml=NewList();
	while(!is_empty_list(l1))
	{
		Term m2=ListFirst(l1);
		int cnum=1;
		Atom np;
		List pl=0;
		ChangeList(l1,0);
		l1=CutFromList(l1,l1);
		pl=MakeList1(MakeCompound2(OPR_MLT,ConsumeCompoundArg(m2,1),
				ConsumeCompoundArg(m2,2)));
		l2=l1;
		while(l2)
		{
			if(EqualTerms(CompoundArgN(m2,3),CompoundArgN(ListFirst(l2),3)))
			{
				pl=AppendLast(pl,MakeCompound2(OPR_MLT,
					ConsumeCompoundArg(ListFirst(l2),1),
					ConsumeCompoundArg(ListFirst(l2),2)));
				l1=CutFromList(l1,l2);
				l2=l1;
				continue;
			}
			l2=ListTail(l2);
		}
		if(UFOutput!=1 && ListLength(pl)==1 && (CompoundArg2(ListFirst(pl))==0 || (
				ListLength(CompoundArg2(ListFirst(pl)))==1 &&
				CompoundArg2(ListFirst(CompoundArg2(ListFirst(pl))))==NewInteger(1))))
		{
			SetCompoundArg(m2,1,CompoundArg1(ListFirst(pl)));
			SetCompoundArg(m2,2,ConsumeCompoundArg(ListFirst(pl),2));
			FreeAtomic(pl);
			nml=AppendLast(nml,m2);
			continue;
		}
		
		mcmp_err=0;
		pl=SortedList(pl,mcmp);
		/*if(mcmp_err)
		{
		WriteTerm(a2);puts("");}*/
		if(UFOutput!=1)
		{
		cnum=(int)IntegerValue(CompoundArg1(ListFirst(pl)));
		if(cnum<0) cnum=-cnum;
		for(l2=ListTail(pl);l2;l2=ListTail(l2))
			cnum=(int)gcf(cnum,IntegerValue(CompoundArg1(ListFirst(l2))));
		if(IntegerValue(CompoundArg1(ListFirst(pl)))<0)
			cnum=-cnum;
		for(l2=pl;l2;l2=ListTail(l2))
			SetCompoundArg(ListFirst(l2),1,
					NewInteger(IntegerValue(CompoundArg1(ListFirst(l2)))/cnum));
		}
		np=abbr_find(pl);
		SetCompoundArg(m2,1,(UFOutput?NewInteger(1):NewInteger(cnum)));
		SetCompoundArg(m2,2,MakeList1(MakeCompound2(OPR_POW,np,NewInteger(1))));
		nml=AppendLast(nml,m2);
	}
	if(has_inf)
	for(l1=nml;l1;l1=ListTail(l1))
	{
		Term m2=ListFirst(l1);
		if(CompoundArgN(m2,3) && CompoundName(ListFirst(CompoundArgN(m2,3)))
				==OPR_PARAMETER)
		{
			Atom prm;
			prm=CompoundArg1(ListFirst(CompoundArgN(m2,3)));
			l2=ConsumeCompoundArg(m2,3);
			l2=CutFromList(l2,l2);
			SetCompoundArg(m2,3,l2);
			l2=ConsumeCompoundArg(m2,2);
			l2=AppendFirst(l2,MakeCompound2(OPR_POW,prm,NewInteger(1)));
			SetCompoundArg(m2,2,l2);
		}
	}
	
	if(has_i)
	for(l1=nml;l1;l1=ListTail(l1))
	{
		Term m2=ListFirst(l1);
		if(CompoundArgN(m2,3) && CompoundArg1(ListFirst(CompoundArgN(m2,3)))==A_I)
		{
			l2=ConsumeCompoundArg(m2,3);
			l2=CutFromList(l2,l2);
			SetCompoundArg(m2,3,l2);
			l2=ConsumeCompoundArg(m2,2);
			l2=AppendFirst(l2,MakeCompound2(OPR_POW,A_I,NewInteger(1)));
			SetCompoundArg(m2,2,l2);
			has_i=1;
		}
	}
	
	SetCompoundArg(a2,5,nml);

}

void alg2_abbr_coeff(Term a2)
{
	List l1,l2;
	int has_i=0, has_c=0;
	List nml;

	if(abbr_coeffs==0)
		return;
		
	if(doing_abbr==0)
		doing_abbr=1;
	
	l1=CompoundArgN(a2,5);
	if(is_empty_list(l1))
		return;
	
	

	if(ListLength(l1)==1)
	{
			return;
	}
	
	for(l1=CompoundArgN(a2,5);l1;l1=ListTail(l1))
	{
		Term m2=ListFirst(l1);
		List l2;
		rpt:
		l2=CompoundArg2(m2);
		for(;l2;l2=ListTail(l2))
		if(ListMember(abbr_coeffs,CompoundArg1(ListFirst(l2))))
		{
			Atom prm=CompoundArg1(ListFirst(l2));
			Integer pw=CompoundArg2(ListFirst(l2));
			List lp=ConsumeCompoundArg(m2,2);
			lp=CutFromList(lp,l2);
			SetCompoundArg(m2,2,lp);
			l2=ConsumeCompoundArg(m2,3);
			l2=AppendFirst(l2,MakeCompound2(OPR_PARAMETER,prm,pw));
			SetCompoundArg(m2,3,l2);
			has_c=1;
			goto rpt;
		}
	}


	
	if(has_c==0)
		return;
	
	if(!has_cmplx_params)				
	for(;l1;l1=ListTail(l1))
	{
		List l3;
		Term m2=ListFirst(l1);
		for(l3=CompoundArg2(m2);l3;l3=ListTail(l3))
		if(CompoundArg1(ListFirst(l3))==A_I)
		{
			l2=ConsumeCompoundArg(m2,2);
			l2=CutFromList(l2,l3);
			SetCompoundArg(m2,2,l2);
			l2=ConsumeCompoundArg(m2,3);
			l2=AppendFirst(l2,MakeCompound2(OPR_SPECIAL,A_I,0));
			SetCompoundArg(m2,3,l2);
			has_i=1;
			break;
		}
	}
	/*
	for(l1=CompoundArgN(a2,5);l1;l1=ListTail(l1))
	{
		Term m2=ListFirst(l1);
		List l2=CompoundArg2(m2);
		for(;l2;l2=ListTail(l2))
		if(GetAtomProperty(CompoundArg1(ListFirst(l2)),
				A_INFINITESIMAL))
		{
			Atom prm=CompoundArg1(ListFirst(l2));
			List lp=ConsumeCompoundArg(m2,2);
			lp=CutFromList(lp,l2);
			SetCompoundArg(m2,2,lp);
			l2=ConsumeCompoundArg(m2,3);
			l2=AppendFirst(l2,MakeCompound2(OPR_PARAMETER,prm,NewInteger(1)));
			SetCompoundArg(m2,3,l2);
			has_inf=1;
			break;
		}
	}
	*/
	l1=ConsumeCompoundArg(a2,5);
	nml=NewList();
	while(!is_empty_list(l1))
	{
		Term m2=ListFirst(l1);
		int cnum;
		Atom np;
		List pl=0;
		ChangeList(l1,0);
		l1=CutFromList(l1,l1);
		pl=MakeList1(MakeCompound2(OPR_MLT,ConsumeCompoundArg(m2,1),
				ConsumeCompoundArg(m2,2)));
		l2=l1;
		while(l2)
		{
			if(EqualTerms(CompoundArgN(m2,3),CompoundArgN(ListFirst(l2),3)))
			{
				pl=AppendLast(pl,MakeCompound2(OPR_MLT,
					ConsumeCompoundArg(ListFirst(l2),1),
					ConsumeCompoundArg(ListFirst(l2),2)));
				l1=CutFromList(l1,l2);
				l2=l1;
				continue;
			}
			l2=ListTail(l2);
		}
		if( ListLength(pl)==1 && (CompoundArg2(ListFirst(pl))==0 || (
				ListLength(CompoundArg2(ListFirst(pl)))==1 &&
				CompoundArg2(ListFirst(CompoundArg2(ListFirst(pl))))==NewInteger(1))))
		{
			SetCompoundArg(m2,1,CompoundArg1(ListFirst(pl)));
			SetCompoundArg(m2,2,ConsumeCompoundArg(ListFirst(pl),2));
			FreeAtomic(pl);
			nml=AppendLast(nml,m2);
			continue;
		}
		
		mcmp_err=0;
		pl=SortedList(pl,mcmp);
		
		/*if(mcmp_err)
		{
		WriteTerm(a2);puts("");}*/
		
		cnum=(int)IntegerValue(CompoundArg1(ListFirst(pl)));
		if(cnum<0) cnum=-cnum;
		for(l2=ListTail(pl);l2;l2=ListTail(l2))
			cnum=(int)gcf(cnum,IntegerValue(CompoundArg1(ListFirst(l2))));
		if(IntegerValue(CompoundArg1(ListFirst(pl)))<0)
			cnum=-cnum;
		for(l2=pl;l2;l2=ListTail(l2))
			SetCompoundArg(ListFirst(l2),1,
					NewInteger(IntegerValue(CompoundArg1(ListFirst(l2)))/cnum));
		
		
		np=abbr_find(pl);
		
		SetCompoundArg(m2,1,NewInteger(cnum));
		SetCompoundArg(m2,2,MakeList1(MakeCompound2(OPR_POW,np,NewInteger(1))));
		nml=AppendLast(nml,m2);
	}
	
		
	for(l1=nml;l1;l1=ListTail(l1))
	{
		Term m2=ListFirst(l1);
		while(CompoundArgN(m2,3) && CompoundName(ListFirst(CompoundArgN(m2,3)))
				==OPR_PARAMETER)
		{
			Atom prm;
			Integer pw;
			prm=CompoundArg1(ListFirst(CompoundArgN(m2,3)));
			pw=CompoundArg2(ListFirst(CompoundArgN(m2,3)));
			l2=ConsumeCompoundArg(m2,3);
			l2=CutFromList(l2,l2);
			SetCompoundArg(m2,3,l2);
			l2=ConsumeCompoundArg(m2,2);
			l2=AppendFirst(l2,MakeCompound2(OPR_POW,prm,pw));
			SetCompoundArg(m2,2,l2);
		}
	}

	
	if(has_i)
	for(l1=nml;l1;l1=ListTail(l1))
	{
		Term m2=ListFirst(l1);
		if(CompoundArgN(m2,3) && CompoundArg1(ListFirst(CompoundArgN(m2,3)))==A_I)
		{
			l2=ConsumeCompoundArg(m2,3);
			l2=CutFromList(l2,l2);
			SetCompoundArg(m2,3,l2);
			l2=ConsumeCompoundArg(m2,2);
			l2=AppendFirst(l2,MakeCompound2(OPR_POW,A_I,NewInteger(1)));
			SetCompoundArg(m2,2,l2);
			has_i=1;
		}
	}

	
	SetCompoundArg(a2,5,nml);

}

extern int sp_cmp(Term,Term);

static void setzprm(List ml)
{
	for(;ml;ml=ListTail(ml))
	{
		List l1;
		List l2=ConsumeCompoundArg(ListFirst(ml),2);
		List l3=ConsumeCompoundArg(ListFirst(ml),3);
		l3=SortedList(l3,sp_cmp);
		for(l1=l2;l1;l1=ListTail(l1))
			if(GetAtomProperty(CompoundArg1(ListFirst(l1)),A_INFINITESIMAL) &&
                            !GetAtomProperty(CompoundArg1(ListFirst(l1)),PROP_TYPE) )
			{
				l3=AppendFirst(l3,ListFirst(l1));
				ChangeList(l1,0);
				l2=CutFromList(l2,l1);
				break;
			}
		
		for(l1=l2;l1;l1=ListTail(l1))
			if(CompoundArg1(ListFirst(l1))==A_I)
			{
				l3=AppendFirst(l3,ListFirst(l1));
				ChangeList(l1,0);
				l2=CutFromList(l2,l1);
				break;
			}
		SetCompoundArg(ListFirst(ml),2,l2);
		SetCompoundArg(ListFirst(ml),3,l3);
	}
}

static void unsetzprm(List ml)
{
	for(;ml;ml=ListTail(ml))
	{
		List l1;
		List l2=ConsumeCompoundArg(ListFirst(ml),2);
		List l3=ConsumeCompoundArg(ListFirst(ml),3);
		for(l1=l3;l1;l1=ListTail(l1))
			if(CompoundArg1(ListFirst(l1))==A_I)
			{
				l2=AppendFirst(l2,ListFirst(l1));
				ChangeList(l1,0);
				l3=CutFromList(l3,l1);
				break;
			}
		for(l1=l3;l1;l1=ListTail(l1))
			if(CompoundName(ListFirst(l1))==OPR_POW)
			{
				l2=AppendFirst(l2,ListFirst(l1));
				ChangeList(l1,0);
				l3=CutFromList(l3,l1);
				break;
			}
		SetCompoundArg(ListFirst(ml),2,l2);
		SetCompoundArg(ListFirst(ml),3,l3);
	}
}

static double eval_mterm(Term m)
{
	double ret=1.0;
	List l;
	for(l=CompoundArg2(m);l;l=ListTail(l))
		ret*=pow(EvalParameter(CompoundArg1(ListFirst(l))),
				(double)IntegerValue(CompoundArg2(ListFirst(l))));
	ret*=(double)IntegerValue(CompoundArg1(CompoundArg1(m)));
	ret/=(double)IntegerValue(CompoundArg2(CompoundArg1(m)));
	return ret;
}

static int mlcmp(Term m1, Term m2)
{
	Term p1=CompoundArg1(ListFirst(CompoundArg2(m1)));
	Term p2=CompoundArg1(ListFirst(CompoundArg2(m2)));
	return FloatValue(p1)>FloatValue(p2);
}

	

void alg2_eval_vrtn(Term a2)
{
	List ml=CompoundArgN(a2,5),l1,l2;
	setzprm(ml);
/*	if(ml)
		WriteTerm(ListFirst(ml));puts("");*/
	for(l1=ml;l1;l1=ListTail(l1))
	{
		double v;
		if(ListFirst(l1)==0)
			continue;
		v=eval_mterm(ListFirst(l1));
		if(ListTail(l1))
			for(l2=ListTail(l1);l2;l2=ListTail(l2))
			{
				if(ListFirst(l2)==0)
					continue;
				if(EqualTerms(CompoundArgN(ListFirst(l1),3),
						CompoundArgN(ListFirst(l2),3)))
				{
					v+=eval_mterm(ListFirst(l2));
					ChangeList(l2,0);
				}
			}
		SetCompoundArg(ListFirst(l1),1,
				MakeCompound2(OPR_DIV,NewInteger(1),NewInteger(1)));
		SetCompoundArg(ListFirst(l1),2,MakeList1(
				MakeCompound2(OPR_POW,NewFloat(v),NewInteger(1))));
	}
rpt:
	for(l1=ml;l1;l1=ListTail(l1))
		if(ListFirst(l1)==0)
		{		
			ml=CutFromList(ml,l1);
			goto rpt;
		}
	unsetzprm(ml);
	ml=SortedList(ml,mlcmp);
	/*if(ml)
		WriteTerm(ListFirst(ml));puts("");*/
	SetCompoundArg(a2,5,ml);
	
	
}
