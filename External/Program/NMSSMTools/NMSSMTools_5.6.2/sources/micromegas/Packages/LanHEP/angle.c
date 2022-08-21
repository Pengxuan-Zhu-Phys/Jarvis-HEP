#include <string.h>
#include "lanhep.h"

static List s_l = 0;
static List a_l = 0;
static List u_l = 0;


static List cplist;

static int vrb=0, rdc=0;
extern int VerbMode;
int opUndefAngleComb=0;
void tri_wrt_sln(List, FILE *);
List tri_find_sub(List);

static int mcmp(Term m1, Term m2)
	{
	int rs;
	m1=CompoundArg2(m1);
	m2=CompoundArg2(m2);
	rs=ListLength(m1)-ListLength(m2);

	if(rs)
		return rs;

	while(!is_empty_list(m1))
		{
		int rs;
		rs=strcmp(AtomValue(CompoundArg1(ListFirst(m1))),
			AtomValue(CompoundArg1(ListFirst(m2))));
		if(rs)
			return rs;
		rs=(int)IntegerValue(CompoundArg2(ListFirst(m1))) -
			(int)IntegerValue(CompoundArg2(ListFirst(m2)));
		if(rs)
			return rs;
		m1=ListTail(m1);
		m2=ListTail(m2);
		}
	puts("eq mt");
	return 0;
	}

extern int infi_order;

Term ProcSetAngle(Term t, Term ind)
	{
	List l;
	Term t1,res;
	if(is_atom(t))
	{
		List l1,l2;
		for(l1=s_l,l2=u_l;l1;l1=ListTail(l1),l2=ListTail(l2))
			{
			printf(" --> ");
			tri_wrt_sln(CompoundArg2(ListFirst(l1)),stdout);
			printf(" : ");
			WriteTerm(ListFirst(l2));
			printf(" matches\n <-- ");
			tri_wrt_sln(CompoundArg1(ListFirst(l1)),stdout);
			puts("");
			}
		return 0;
        }
	
	if(!is_compound(t)||CompoundArity(t)!=1)
		{
		ErrorInfo(564);
		puts("SetAngle: wrong arguments");
		return 0;
		}
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);

	if(t1==NewInteger(1))
		{
		vrb=1;
		return 0;
		}

    if(t1==NewInteger(2))
		{
        WriteTerm(a_l);puts("");
        DumpList(s_l);
		return 0;
        }

	if(t1==NewInteger(3))
		{
		List l1,l2;
		for(l1=s_l,l2=u_l;l1;l1=ListTail(l1),l2=ListTail(l2))
			{
			WriteTerm(CompoundArg2(ListFirst(l1)));
			printf(" : ");
			WriteTerm(ListFirst(l2));
			puts("");
			}
		return 0;
        }

	if(t1==NewInteger(4))
		{
		rdc=1;
		return 0;
		}

	if(!is_compound(t1) || CompoundArity(t1)!=2)
		{
		ErrorInfo(565);
		puts("SetAngle: wrong arguments");
		return 0;
		}
	t=ConsumeCompoundArg(t1,1);
	res=ConsumeCompoundArg(t1,2);
	FreeAtomic(t1);

    if(res==OPR_USCORE)
        goto xyz;
	
	{
		int iosav=infi_order;
		infi_order=0;
	    res=ExprTo1(res);
		infi_order=iosav;
	}
	if(res==0)
		return 0;
	
	if(CompoundArg2(res))
		{
		ErrorInfo(567);
		puts("SetAngle: right expression is not scalar.");
		return 0;
		}
		
	res=Alg1to2(res);
	if(res==0)
		return 0;
	res=alg2_add(NewList(),res);
	
	for(l=res;l;l=ListTail(l))
		{
		Term a2,nn,ll;
		int ncf;
		a2=ListFirst(l);
		alg2_common_n(a2);
		nn=CompoundArg2(a2);
        if(CompoundArg2(nn)!=NewInteger(1))
            {
                ErrorInfo(569);
                puts("SetAngle: expression has to have no common numeric divider.");
                return 0;
            }
		ncf=(int)IntegerValue(CompoundArg1(nn));
        for(ll=CompoundArgN(a2,5);ll;ll=ListTail(ll))
           {
           int n;
           n=(int)IntegerValue(CompoundArg1(ListFirst(ll)));
           SetCompoundArg(ListFirst(ll),1,NewInteger(n*ncf));
           }
		}
	t1=ConsumeCompoundArg(ListFirst(res),5);
	FreeAtomic(res);
	res=t1;	
xyz:
	
	{
		int iosav=infi_order;
		infi_order=0;
		t=ExprTo1(t);
		infi_order=iosav;
	}
	if(t==0)
		return 0;
	if(CompoundArg2(t))
		{
		ErrorInfo(567);
		puts("SetAngle: left expression is not scalar.");
		return 0;
		}
	t=Alg1to2(t);
	if(t==0)
		return 0;
	t=alg2_add(NewList(),t);
	

	for(l=t;l;l=ListTail(l))
		{
        Term a2,nn,ll;
		int ncf;
		a2=ListFirst(l);
		alg2_common_n(a2);
        nn=CompoundArg2(a2);
        if(CompoundArg2(nn)!=NewInteger(1))
            {
                ErrorInfo(569);
                puts("SetAngle: expression has to have no common numeric divider.");
                return 0;
            }
		ncf=(int)IntegerValue(CompoundArg1(nn));
        for(ll=CompoundArgN(a2,5);ll;ll=ListTail(ll))
           {
           int n;
           n=(int)IntegerValue(CompoundArg1(ListFirst(ll)));
           SetCompoundArg(ListFirst(ll),1,NewInteger(n*ncf));
           }
		alg2_red_cos(a2);
		}
	t1=ConsumeCompoundArg(ListFirst(t),5);
	FreeAtomic(t);

	if(is_empty_list(t1))
		{
		ErrorInfo(568);
		puts("SetAngle: expression calcels");
		return 0;
		}
	t1=SortedList(t1,mcmp);
    for(l=t1;l;l=ListTail(l))
        {
        Term m2;
        List l1;
        m2=ListFirst(l);
        for(l1=CompoundArg2(m2);l1;l1=ListTail(l1))
            {
            Atom n;
            n=CompoundArg1(ListFirst(l1));
            if(!GetAtomProperty(n,A_SET_ANGLE))
                {
                SetAtomProperty(n,A_SET_ANGLE,NewInteger(1));
                a_l=AppendFirst(a_l,n);
                }
            }
        }

	l=tri_find_sub(t1);
	if(l)
	{
		if(!EqualTerms(l,res))
		{
			if(VerbMode)
			{
			printf("Warning: SetAngle: redefinition of substitution for ");
			tri_wrt_sln(t1,stdout);puts("");
			printf("\tpreviuos: ");
			tri_wrt_sln(l,stdout);
			printf("\nnew:      ");
			tri_wrt_sln(res,stdout);
			puts("");
			}
			s_l=AppendFirst(s_l,MakeCompound2(OPR_DIV,t1,res));
			u_l=AppendFirst(u_l,NewInteger(0));
		}
		else
		{
			FreeAtomic(t1);
			FreeAtomic(res);
		}
		FreeAtomic(l);
	}
	else
	{
		s_l=AppendFirst(s_l,MakeCompound2(OPR_DIV,t1,res));
		u_l=AppendFirst(u_l,NewInteger(0));
	}

	return 0;

	}



void tri_wrt_sln(List mm, FILE *f)
	{
	List s;
	int first=1;
	int cf;
	s=mm;
	while(!is_empty_list(s))
		{
		Term m2;
		List lp;
		int ast=0;
		m2=ListFirst(s);
		cf=(int)IntegerValue(CompoundArg1(m2));
		lp=CompoundArg2(m2);
		if(cf<0)
			{
			fprintf(f,"-");
			cf=-cf;
			}
		else
			if(!first)
				 fprintf(f,"+");
		first=0;
		if(cf!=1 || (is_empty_list(lp)))
			{
			fprintf(f,"%d",cf);
			ast=1;
			}
		while(!is_empty_list(lp))
			{
			int po;
			po=(int)IntegerValue(CompoundArg2(ListFirst(lp)));
			if(ast)
				fprintf(f,"*");
			else
				ast=1;
			fprintf(f,"%s",AtomValue(CompoundArg1(ListFirst(lp))));
			if(po!=1)
				fprintf(f,"**%d",po);
			lp=ListTail(lp);
			}
		s=ListTail(s);
		}
	}




static List add_fct(List m_l, Term p, Term m2)
	{
	List l1;
    List tf;
	l1=m_l;
    tf=CompoundArgN(m2,3);
	
	while(!is_empty_list(l1))
		{
        Term tt;
        tt=CompoundArg2(ListFirst(l1));
        tt=CompoundArgN(ListFirst(tt),3);
        if(EqualTerms(CompoundArg1(ListFirst(l1)),p) &&
           EqualTerms(tt,tf))
			{
			List ml;
			ml=ConsumeCompoundArg(ListFirst(l1),2);
			ml=AppendLast(ml,m2);
			SetCompoundArg(ListFirst(l1),2,ml);
			FreeAtomic(p);
			return m_l;
			}
		l1=ListTail(l1);
		}
	return AppendLast(m_l,MakeCompound2(OPR_DIV,p,
		AppendLast(NewList(),m2)));
	}

static void l_fit(List, List, int *);
static void l_set(List, int *, Term *);
int tri_heu(List);

List tri_find_sub(List ml)
{
	List l1;
	int nco=0, cl;
	List sco=0;
	
	cl=ListLength(ml);
	
	for(l1=s_l;l1;l1=ListTail(l1))
		{
		List l3;
		l3=CompoundArg1(ListFirst(l1));
		if(ListLength(l3)!=cl)
			continue;
		l_fit(l3,ml,&nco);
		if(nco)
			{
			sco=CopyTerm(CompoundArg2(ListFirst(l1)));
			break;
			}
		}
		
	if(nco==0)
		return 0;
	
	for(l1=sco;l1;l1=ListTail(l1))
		SetCompoundArg(ListFirst(l1),1,
				NewInteger(nco*IntegerValue(CompoundArg1(ListFirst(l1)))));
	
	return sco;
}

static int fflag=0;

static void hl2m2(List hhl)
	{
	Term hh;
	List ml,l1,l2;
	Term sco,m2,tf,m22;
	int cl,nco;
	fflag=0;
	
	hh=ListFirst(hhl);
    if(is_compound(hh) && CompoundName(hh)==A_MTERM)
        return;

	ml=CompoundArg2(hh);
	cl=ListLength(ml);
	if(cl==1 && ListLength(CompoundArg2(ListFirst(ml)))<2)
		{
		ml=ConsumeCompoundArg(hh,2);
		m2=ListFirst(ml);
		RemoveList(ml);
		if(CompoundArg1(hh))
			{
			l1=ConsumeCompoundArg(m2,2);
			l1=ConcatList(l1,ConsumeCompoundArg(hh,1));
			SetCompoundArg(m2,2,l1);
			}
		ChangeList(hhl,m2);
		FreeAtomic(hh);
		return;
		}
	nco=0;
	sco=0;
	for(l1=s_l,l2=u_l;l1;l1=ListTail(l1),l2=ListTail(l2))
		{
		List l3;
		l3=CompoundArg1(ListFirst(l1));
		if(ListLength(l3)!=cl)
			continue;
		l_fit(l3,ml,&nco);
		if(nco)
			{
			sco=CopyTerm(CompoundArg2(ListFirst(l1)));
			ChangeList(l2,
				NewInteger(IntegerValue(ListFirst(l2))+1));
			break;
			}
		}

	if(nco==0 && tri_heu(ml))
	{
	for(l1=s_l,l2=u_l;l1;l1=ListTail(l1),l2=ListTail(l2))
		{
		List l3;
		l3=CompoundArg1(ListFirst(l1));
		if(ListLength(l3)!=cl)
			continue;
		l_fit(l3,ml,&nco);
		if(nco)
			{
			sco=CopyTerm(CompoundArg2(ListFirst(l1)));
			ChangeList(l2,
				NewInteger(IntegerValue(ListFirst(l2))+1));
			break;
			}
		}
	/*if(nco==0)
		{
			fflag=1;
		puts("Internal error (triheu failed)");
		DumpList(ml);puts("");
		}*/
	}	
		
    if(sco==OPR_USCORE || (nco==0 && opUndefAngleComb==0))
        {
        List l;
        ml=ConsumeCompoundArg(hh,2);
        if(CompoundArg1(hh))
            {
            Term c1;
            c1=ConsumeCompoundArg(hh,1);
            for(l=ml;l;l=ListTail(l))
                {
                List cl;
                cl=ConsumeCompoundArg(ListFirst(l),2);
                cl=ConcatList(cl,CopyTerm(c1));
                SetCompoundArg(ListFirst(l),2,cl);
                }
            FreeAtomic(c1);
            }
        FreeAtomic(hh);
        l=ml;
        ChangeList(hhl,ListFirst(l));
        l=ListTail(l);
        for(;l;l=ListTail(l))
            {
            InsertList(hhl,ListFirst(l));
            hhl=ListTail(hhl);
            }
        return;
        }

	if(nco==0)
		l_set(ml,&nco,&sco);


	tf=CompoundArgN(ListFirst(ml),3);
	m22=ConsumeCompoundArg(hh,1);
	for(l1=sco;l1;l1=ListTail(l1))
	{
		m2=MakeCompound(A_MTERM,3);
		SetCompoundArg(m2,1,NewInteger(nco*IntegerValue(
				CompoundArg1(ListFirst(l1)))));
		if(m22)
			SetCompoundArg(m2,2,ConcatList(
					ConsumeCompoundArg(ListFirst(l1),2),
					CopyTerm(m22)));
		else
			SetCompoundArg(m2,2,
			    ConsumeCompoundArg(ListFirst(l1),2));
		
		SetCompoundArg(m2,3,CopyTerm(tf));
		if(l1==sco)
			ChangeList(hhl,m2);
		else
			InsertList(hhl,m2);
	}
	FreeAtomic(hh);
	FreeAtomic(sco);
	
	return;
		
    if(is_list(sco))
        {
        m2=MakeCompound(A_MTERM,3);
        SetCompoundArg(m2,1,NewInteger(nco));
        SetCompoundArg(m2,3,CopyTerm(CompoundArgN(ListFirst(ml),3)));
        l1=sco;
        if(CompoundArg1(hh))
            l1=ConcatList(l1,ConsumeCompoundArg(hh,1));
        SetCompoundArg(m2,2,l1);
        ChangeList(hhl,m2);
        FreeAtomic(hh);
        return;
        }
    else
        {
        Term tf;
        tf=CopyTerm(CompoundArgN(ListFirst(ml),3));
        m2 =MakeCompound(A_MTERM,3);
        SetCompoundArg(m2,1,
            NewInteger(nco*IntegerValue(CompoundArg2(sco))));
        l1=ConsumeCompoundArg(sco,1);
        if(CompoundArg1(hh))
            l1=ConcatList(l1,CopyTerm(CompoundArg1(hh)));
        SetCompoundArg(m2,2,l1);
        SetCompoundArg(m2,3,CopyTerm(tf));
        ChangeList(hhl,m2);
        m2=MakeCompound(A_MTERM,3);
        SetCompoundArg(m2,1,
            NewInteger(nco*IntegerValue(CompoundArgN(sco,4))));
        l1=ConsumeCompoundArg(sco,3);
        if(CompoundArg1(hh))
            l1=ConcatList(l1,ConsumeCompoundArg(hh,1));
        SetCompoundArg(m2,2,l1);
        SetCompoundArg(m2,3,tf);
        FreeAtomic(sco);
        FreeAtomic(hh);
        InsertList(hhl,m2);
        return;
        }

	}

	
void alg2_red_sico(Term a2)
	{
	List m2,l1,l2;
	List hl;
    int has_cs=0;
	int u;
	if(!s_l)
		return;

	m2=CompoundArgN(a2,5);
	cplist=CompoundArg1(a2);
	l1=m2;
	if(is_empty_list(l1) || is_empty_list(ListTail(l1)))
		return;
	while(!is_empty_list(l1))
		{
		int nocsc;
		nocsc=0;
		l2=CompoundArg2(ListFirst(l1));
		while(!is_empty_list(l2))
			{
			Atom p;
			p=CompoundArg1(ListFirst(l2));
            if(!GetAtomProperty(p,A_SET_ANGLE))
				nocsc++;
            else
                has_cs=1;
			l2=ListTail(l2);
			}
		l1=ListTail(l1);
		}

    if(!has_cs)
        return;

    if(vrb)
		{
		WriteVertex(cplist);
		puts("");
        }

	m2=ConsumeCompoundArg(a2,5);

	u=(int)IntegerValue(CompoundArg2(CompoundArg2(a2)));
	SetCompoundArg(CompoundArg2(a2),2,NewInteger(u*16384));
	for(l1=m2;l1;l1=ListTail(l1))
		{
		u=(int)IntegerValue(CompoundArg1(ListFirst(l1)));
		SetCompoundArg(ListFirst(l1),1,NewInteger(u*16384));
		}

	hl=NewList();
	l1=m2;

	while(!is_empty_list(l1))
		{
		Term mm,nct;
		List ll;
		mm=ListFirst(l1);
		ll=ConsumeCompoundArg(mm,2);
		nct=0;
	xyz1:
		for(l2=ll;l2;l2=ListTail(l2))
			{
            Atom p;
			p=CompoundArg1(ListFirst(l2));
            if(!GetAtomProperty(p,A_SET_ANGLE))
				{
				nct=AppendLast(nct,ListFirst(l2));
				ChangeList(l2,0);
				ll=CutFromList(ll,l2);
				goto xyz1;
				}
			}
		SetCompoundArg(mm,2,ll);
		
		hl=add_fct(hl,nct,mm);
		l1=ListTail(l1);
		}

    if(vrb)
		{
        DumpList(hl);
        }
	RemoveList(m2);

	for(l1=hl;l1;l1=ListTail(l1))
		{
		l2=ConsumeCompoundArg(ListFirst(l1),2);
		l2=SortedList(l2,mcmp);
		SetCompoundArg(ListFirst(l1),2,l2);
		}

    if(vrb)
		{
        DumpList(hl);
        }

	for(l1=hl;l1;l1=ListTail(l1))
	{
		hl2m2(l1);
		if(fflag)
		{
			WriteVertex(CompoundArg1(a2));
			puts("");
		}
	}

	if(vrb)
		DumpList(hl);

	SetCompoundArg(a2,5,hl);
	}


Term ProcEval(Term t, Term ind)
	{
	List l;
	t=ExprTo1(ConsumeCompoundArg(t,1));
	if(t==0)
		return 0;
	t=Alg1to2(t);
	if(t==0)
		return 0;
	t=alg2_add(NewList(),t);

	for(l=t;l;l=ListTail(l))
		{
		Term a2;
		a2=ListFirst(l);
		alg2_common_n(a2);
		alg2_common_s(a2);
		alg2_red_cos(a2);
		alg2_red_orth(a2);
		alg2_recommon_n(a2);
		alg2_recommon_s(a2);
		}
	printf("\t");
	WriteTerm(CompoundArg2(ListFirst(t)));
	printf("*");
	WriteTerm(CompoundArgN(ListFirst(t),3));
	printf("\n");
	DumpList(CompoundArgN(ListFirst(t),5));
	FreeAtomic(t);
	return 0;
	}



static int gcf_list(List l)
	{
	int ret;
	List sv;
	sv=l;
	ret=(int)IntegerValue(CompoundArg1(ListFirst(l)));
	l=ListTail(l);
	while(!is_empty_list(l))
		{
		ret=(int)gcf(ret,IntegerValue(CompoundArg1(ListFirst(l))));
		l=ListTail(l);
		}
	for(l=sv;l;l=ListTail(l))
		{
		int i;
		i=(int)IntegerValue(CompoundArg1(ListFirst(l)));
		SetCompoundArg(ListFirst(l),1,NewInteger(i/ret));
		}
	return ret;
	}


static int new_p=0;

static void to_sums(List);

void tri_set_sub(List ml)
	{
	int nco;
	Term sco;
	
	l_set(ml,&nco,&sco);
	FreeAtomic(sco);
	}
	

static void l_set(List l, int *nco, Term *sco)
	{
    char buf[8];
    Atom np;
	Term m2;
    l=CopyTerm(l);
    sprintf(buf,"aa%03d",new_p++);
    np=NewAtom(buf,5);
 /*   *sco=AppendFirst(NewList(),MakeCompound2(OPR_POW,np,NewInteger(1)));*/
	*nco=gcf_list(l);
	if(rdc==0)
		{
		puts("Warning: unknown angle combination; use:");
		printf("\tSetAngle(");
		tri_wrt_sln(l,stdout);
    	printf("=%5.5s).\n",buf);
		to_sums(l);
		}
	else
		{
		printf("%5.5s:=\n",buf);
		tri_wrt_sln(l,stdout);
		printf(";\n");
		}

	m2=MakeCompound(A_MTERM,3);
	SetCompoundArg(m2,1,NewInteger(1));
	SetCompoundArg(m2,2,AppendLast(NewList(),
			MakeCompound2(OPR_POW,np,NewInteger(1))));
	s_l=AppendFirst(s_l,MakeCompound2(OPR_DIV,l,
            AppendLast(NewList(),m2)));
	*sco=CopyTerm(CompoundArg2(ListFirst(s_l)));
	u_l=AppendFirst(u_l,NewInteger(1));
	SetAtomProperty(np,PROP_TYPE,OPR_PARAMETER);
	}

static void l_fit(List l1, List l2, int *nco)
	{
	int ic;
	ic=(int)IntegerValue(CompoundArg1(ListFirst(l2)))
		/ (int)IntegerValue(CompoundArg1(ListFirst(l1)));
	if(ic==0)
		return;
	while(!is_empty_list(l1))
		{
		int i1,i2;
		i1=(int)IntegerValue(CompoundArg1(ListFirst(l1)));
		i2=(int)IntegerValue(CompoundArg1(ListFirst(l2)));
		if(i2!=ic*i1)
			return;
		if(!EqualTerms(CompoundArg2(ListFirst(l1)),
				CompoundArg2(ListFirst(l2))))
			return;
		l1=ListTail(l1);
		l2=ListTail(l2);
		}
	*nco=ic;
	}

static int ccmp(Term a1, Term a2)
	{
	char *s1, *s2;
	s1=AtomValue(CompoundArg2(a1));
	s2=AtomValue(CompoundArg2(a2));
	return s1[0]-s2[0];
	}

static List tr_add(List l, Term a)
	{
	List l1;

	{
	List ll;
	ll=CompoundArg1(CompoundArg2(a));
	if(!is_empty_list(ll))
		{
		ll=ConsumeCompoundArg(CompoundArg2(a),1);
		ll=SortedList(ll,ccmp);
		if(IntegerValue(CompoundArg1(ListFirst(ll)))<0)
			{
			for(l1=ll;l1;l1=ListTail(l1))
				SetCompoundArg(ListFirst(l1),1,
				NewInteger(-IntegerValue(CompoundArg1(
				    ListFirst(l1)))));
			if(AtomValue(CompoundName(CompoundArg2(a)))[0]=='s')
				SetCompoundArg(CompoundArg1(a),1,
					NewInteger(-IntegerValue(
				    	CompoundArg1(CompoundArg1(a)))));
			}
		SetCompoundArg(CompoundArg2(a),1,ll);
		}
	}


	{
	Term a2;
	int n,d,f;
	a2=CompoundArg1(a);
	n=(int)IntegerValue(CompoundArg1(a2));
	d=(int)IntegerValue(CompoundArg2(a2));
	f=(int)gcf(n,d);
	if(f!=1)
		{
		SetCompoundArg(a2,1,NewInteger(n/f));
		SetCompoundArg(a2,2,NewInteger(d/f));
		}
	}

	for(l1=l;l1;l1=ListTail(l1))
		{
		Term a2;
		a2=ListFirst(l1);
		if(EqualTerms(CompoundArg2(a2),CompoundArg2(a)))
			{
			long int n,d,n1,d1,n2,d2,f;
			n1=IntegerValue(CompoundArg1(CompoundArg1(a)));
			d1=IntegerValue(CompoundArg2(CompoundArg1(a)));
			n2=IntegerValue(CompoundArg1(CompoundArg1(a2)));
			d2=IntegerValue(CompoundArg2(CompoundArg1(a2)));
			n=n1*d2+n2*d1;
			d=d1*d2;

			if(n==0)
				{
				FreeAtomic(a);
				return CutFromList(l,l1);
				}

			f=gcf(n,d);
			n/=f;
			d/=f;
			SetCompoundArg(CompoundArg1(a2),1,NewInteger(n));
			SetCompoundArg(CompoundArg1(a2),2,NewInteger(d));
			FreeAtomic(a);
			return l;
			}
		}

	return AppendLast(l,a);

	}

static List tr_mlt1(Term a1, Term a2)
	{
	long int n,d;
	List l1,l2,ls,ld,ll;
	char *f1, *f2;

	{
	long int n1,n2,d1,d2;
	n1=IntegerValue(CompoundArg1(CompoundArg1(a1)));
	d1=IntegerValue(CompoundArg2(CompoundArg1(a1)));
	n2=IntegerValue(CompoundArg1(CompoundArg1(a2)));
	d2=IntegerValue(CompoundArg2(CompoundArg1(a2)));
	n=n1*n2;
	d=d1*d2*2;
	}

	l1=CompoundArg1(CompoundArg2(a1));
	l2=CompoundArg1(CompoundArg2(a2));
	f1=AtomValue(CompoundName(CompoundArg2(a1)));
	f2=AtomValue(CompoundName(CompoundArg2(a2)));

	if(f1[0]!='c' && f1[0]!='s')
		{
		printf("SetAngle: unknown trig function %s\n",f1);
		return 0;
		}
	if(f2[0]!='c' && f2[0]!='s')
		{
		printf("SetAngle: unknown trig function %s\n",f2);
		return 0;
		}

	ls=CopyTerm(l1);
	ld=CopyTerm(l1);

	for(ll=l2;ll;ll=ListTail(ll))
		{
		Atom aa;
		List ll1;
		aa=CompoundArg2(ListFirst(ll));
		for(ll1=ls;ll1;ll1=ListTail(ll1))
			{
			if(CompoundArg2(ListFirst(ll1))==aa)
				{
				long int n;
				n=IntegerValue(CompoundArg1(ListFirst(ll1)))+
					IntegerValue(CompoundArg1(
						ListFirst(ll)));
				if(n==0)
					ls=CutFromList(ls,ll1);
				else
					SetCompoundArg(ListFirst(ll1),1,
					    NewInteger(n));
				break;
				}
			}
		if(is_empty_list(ll1))
			ls=AppendLast(ls,CopyTerm(ListFirst(ll)));

		for(ll1=ld;ll1;ll1=ListTail(ll1))
			{
			if(CompoundArg2(ListFirst(ll1))==aa)
				{
				long int n;
				n=IntegerValue(CompoundArg1(ListFirst(ll1)))-
					IntegerValue(CompoundArg1(ListFirst(ll)));
				if(n==0)
					ld=CutFromList(ld,ll1);
				else
					SetCompoundArg(ListFirst(ll1),1,
					    NewInteger(n));
				break;
				}
			}
		if(is_empty_list(ll1))
			{
			Term q;
			q=CopyTerm(ListFirst(ll));
			SetCompoundArg(q,1,NewInteger(-IntegerValue(
				CompoundArg1(q))));
			ld=AppendLast(ld,q);
			}
		}

	if(f1[0]=='c' && f2[0]=='c')
		{
		List res;
		res=AppendLast(NewList(),
			MakeCompound2(OPR_MLT,
				MakeCompound2(OPR_DIV,NewInteger(n),
					NewInteger(d)),
				MakeCompound1(A_COS,ls)));
		res=AppendLast(res,
			MakeCompound2(OPR_MLT,
				MakeCompound2(OPR_DIV,NewInteger(n),
					NewInteger(d)),
				MakeCompound1(A_COS,ld)));
		return res;
		}

	if(f1[0]=='s' && f2[0]=='s')
		{
		List res;
		res=AppendLast(NewList(),
			MakeCompound2(OPR_MLT,
				MakeCompound2(OPR_DIV,NewInteger(-n),
					NewInteger(d)),
				MakeCompound1(A_COS,ls)));
		res=AppendLast(res,
			MakeCompound2(OPR_MLT,
				MakeCompound2(OPR_DIV,NewInteger(n),
					NewInteger(d)),
				MakeCompound1(A_COS,ld)));
		return res;
		}

	if(f1[0]=='s' && f2[0]=='c')
		{
		List res;
		res=NewList();
		if(!is_empty_list(ls))
		res=AppendLast(res,
			MakeCompound2(OPR_MLT,
				MakeCompound2(OPR_DIV,NewInteger(n),
			    		NewInteger(d)),
				MakeCompound1(A_SIN,ls)));

		if(!is_empty_list(ld))
		res=AppendLast(res,
			MakeCompound2(OPR_MLT,
				MakeCompound2(OPR_DIV,NewInteger(n),
			   	 	NewInteger(d)),
				MakeCompound1(A_SIN,ld)));
		return res;
		}

	if(f1[0]=='c' && f2[0]=='s')
		{
		List res;
		res=NewList();
		if(!is_empty_list(ls))
		res=AppendLast(res,
			MakeCompound2(OPR_MLT,
				MakeCompound2(OPR_DIV,NewInteger(n),
			    		NewInteger(d)),
				MakeCompound1(A_SIN,ls)));

		if(!is_empty_list(ld))
		res=AppendLast(res,
			MakeCompound2(OPR_MLT,
				MakeCompound2(OPR_DIV,NewInteger(-n),
			    		NewInteger(d)),
				MakeCompound1(A_SIN,ld)));
		return res;
		}

	puts("Internal error (mktf)");
	return 0;

	}



static List tr_mlt(List l, Term a)
	{
	List l1;
	List res;
	res=NewList();
	for(l1=l;l1;l1=ListTail(l1))
		res=ConcatList(res,tr_mlt1(ListFirst(l1),a));
	FreeAtomic(l);
	FreeAtomic(a);
	return res;
	}


static void to_sums(List l)
	{
	List l1,l2;
	List res=0;
	
	for(l1=l;l1;l1=ListTail(l1))
		{
		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
			{
			if(!GetAtomProperty(CompoundArg1(ListFirst(l2)),A_TRIG_FU))
				{
		/*	printf("(symbol %s is not a trigonometric function)",
				AtomValue(CompoundArg1(ListFirst(l2)))); */
				return;
				}
			}
		}

	printf("  == ");
		
	for(l1=l;l1;l1=ListTail(l1))
		{
		List lcuf, res1;
		lcuf=NewList();
		res1=NewList();

		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
			{
			Term cuf;
			int pw;
			cuf=GetAtomProperty(CompoundArg1(ListFirst(l2)),
							A_TRIG_FU);
			pw=(int)IntegerValue(CompoundArg2(ListFirst(l2)));
			cuf=CopyTerm(cuf);
			cuf=MakeCompound2(OPR_MLT,
				MakeCompound2(OPR_DIV,NewInteger(1),
					NewInteger(1)),cuf);
			for(;pw>1;pw--)
				lcuf=AppendLast(lcuf,CopyTerm(cuf));
			lcuf=AppendLast(lcuf,cuf);
			}
/*
		WriteTerm(CompoundArg1(ListFirst(l1)));
		printf(" "); WriteTerm(lcuf);
		puts("");
 */

		if(is_empty_list(lcuf))
			{
			res=tr_add(res,MakeCompound2(OPR_MLT,
				MakeCompound2(OPR_DIV,
					CompoundArg1(ListFirst(l1)),
					NewInteger(1)),
					MakeCompound1(A_COS,0)));
			}
		else
			{
			res1=AppendFirst(NewList(),ListFirst(lcuf));
			ChangeList(lcuf,0);
			lcuf=CutFromList(lcuf,lcuf);
			for(l2=lcuf;l2;l2=ListTail(l2))
				res1=tr_mlt(res1,ListFirst(l2));
			for(l2=res1;l2;l2=ListTail(l2))
				{
				Term a1;
				a1=ListFirst(l2);
				SetCompoundArg(CompoundArg1(a1),1,NewInteger(
					IntegerValue(CompoundArg1(
						CompoundArg1(a1)))*
					IntegerValue(CompoundArg1(
						ListFirst(l1)))));
				res=tr_add(res,a1);
				}
			RemoveList(lcuf);
			RemoveList(res1);
			}
		}

	for(l1=res;l1;l1=ListTail(l1))
		{
		int n,d;
		Term ff;
		n=(int)IntegerValue(CompoundArg1(CompoundArg1(ListFirst(l1))));
		d=(int)IntegerValue(CompoundArg2(CompoundArg1(ListFirst(l1))));
		ff=CompoundArg2(ListFirst(l1));
		printf("%+d",n);
		if(d!=1)
			printf("/%d",d);
		printf("*");
		if(AtomValue(CompoundName(ff))[0]=='s')
			printf("sin(");
		else
			printf("cos(");
		if(is_empty_list(CompoundArg1(ff)))
			printf("0");
		else
			{
			int first = 1;
			for(l2=CompoundArg1(ff);l2;l2=ListTail(l2))
				{
				int n;
				char *v;
				n=(int)IntegerValue(CompoundArg1(ListFirst(l2)));
				v=AtomValue(CompoundArg2(ListFirst(l2)));
				if(first)
					{
					if(n!=1)
						printf("%d*%s",n,v);
					else
						printf("%s",v);
					}
				else
					{
					if(n==1 || n==-1)
						printf("%c%s",n==1?'+':'-',v);
					else
						printf("%+d%s",n,v);
					}
				first=0;
				}
			}
		printf(")");
		}

	puts("");

	}
