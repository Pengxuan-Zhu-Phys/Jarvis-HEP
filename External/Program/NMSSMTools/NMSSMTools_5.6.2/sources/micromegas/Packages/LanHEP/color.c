#include "lanhep.h"


Atom ColorLambda=0, ColorF=0, ColorD=0, ColorEps=0, ColorEpsB=0, 
		ColorK6=0, ColorK6b=0, ColorL6=0;

int is_color_eps(Atom s)
{
	if(s!=ColorEps && s!=ColorEpsB)
	{
		puts("Internal error (ceps)");
		return 0;
	}
	return s==ColorEps;
}

int is_color_d(Atom s)
{
	if(ColorD==0)
		return 0;
	return ColorD==s;
}

Term inv_color_eps(Atom s)
{
	if(s!=ColorEps && s!=ColorEpsB)
	{
		puts("Internal error (cepi)");
		return s;
	}
	
	return s==ColorEps?ColorEpsB:ColorEps;
}

Term inv_color_k6(Atom s)
{
	if(s!=ColorK6 && s!=ColorK6b)
	{
		puts("Internal error (cepi)");
		return s;
	}
	
	return s==ColorK6?ColorK6b:ColorK6;
}
int color_repres(Term ind)
	{
	Atom a1,a2;
	char *v1, *v2;
	ind=CompoundArg1(ind);
	if(CompoundName(ind)!=A_COLOR)
		return 0;
	a1=CompoundArg1(ind);
	a2=CompoundArg2(ind);
	v1=AtomValue(a1);
	if(v1[1]=='8')
	  return 3;
	if(v1[1]=='3' && v1[2]==0)
	  return 1;
	if(v1[1]=='3' && v1[2]=='b')
	  return 2;
	
	if(v1[1]=='6' && v1[2]==0)
	  return 6;
	if(v1[1]=='6' && v1[2]=='b')
	  return 7;
	if(a1==a2)
		return 3;
	if(a1<a2)
		return 1;
	return 2;
	}

static int f_symm_factor(List *ll)
	{
	List gprt,l1,l2;
	int ret;
	int flag;
	ret=1;
	gprt=*ll;

/*	WriteTerm(gprt); printf("  %d \n",ret);*/
		
	if(ListLength(gprt)<2)
		{
		
		return 1;
		}
	do
		{
		flag=0;
		l1=gprt;
		l2=ListTail(l1);
		while(!is_empty_list(l2))
			{
			Atom a1, a2;
			a1=ListFirst(l1);
			a2=ListFirst(l2);
			if(IntegerValue(a1)>IntegerValue(a2))
				{
				ChangeList(l1,a2);
				ChangeList(l2,a1);
				flag=1;
				ret*=-1;
				}
			l1=l2;
			l2=ListTail(l2);
			}
		}
		while(flag);
				

	*ll=gprt;
	return ret;
	}

void color_symm_f(List pl, Term m2)
	{
	
	int fact=1;
	List l;
	l=CompoundArgN(m2,3);
	while(!is_empty_list(l))
		{
		Term t;
		t=ListFirst(l);
		
		if(CompoundName(t)==OPR_SPECIAL && (
				GetAtomProperty(CompoundArg1(t),A_COLOR)==A_COLOR_F ||
				GetAtomProperty(CompoundArg1(t),A_COLOR)==A_COLOR_EPS ||
				GetAtomProperty(CompoundArg1(t),A_SYMM)==NewInteger(-1)))
			{
			List ll;
			ll=ConsumeCompoundArg(t,2);
			fact*=f_symm_factor(&ll);
			SetCompoundArg(t,2,ll);
			}
		if(CompoundName(t)==OPR_SPECIAL && (
				GetAtomProperty(CompoundArg1(t),A_COLOR)==A_COLOR_D ||
				GetAtomProperty(CompoundArg1(t),A_SYMM)==NewInteger(1)))
			{
			List ll;
			ll=ConsumeCompoundArg(t,2);
			f_symm_factor(&ll);
			SetCompoundArg(t,2,ll);
			}
		if(CompoundName(t)==OPR_SPECIAL && (
				GetAtomProperty(CompoundArg1(t),A_COLOR)==A_COLOR_K6))
		{
		  List ill=CompoundArg2(t);
		  Integer in1=ListFirst(ill);
		  Integer in2=ListFirst(ListTail(ill));
		  if(IntegerValue(in1)>IntegerValue(in2))
		  {
		    ChangeList(ill,in2);
		    ChangeList(ListTail(ill),in1);
		  }
		}
		
		l=ListTail(l);
		}

	if(fact-1)
		{
		Term t;
		t=CompoundArg1(m2);
		SetCompoundArg(t,1,NewInteger(-IntegerValue(CompoundArg1(t))));
		}		
/*	
	l=CompoundArgN(m2,3);
	while(!is_empty_list(l))
		{
		Term t;
		t=ListFirst(l);	
		if(CompoundName(t)==OPR_SPECIAL && 
			(GetAtomProperty(CompoundArg1(t),A_COLOR)==A_COLOR_F||
				GetAtomProperty(CompoundArg1(t),A_COLOR)==A_COLOR_D))
			{
			List ll;
			ll=ListTail(l);
			while(!is_empty_list(ll))
				{
				Term t1;
				t1=ListFirst(ll);
				if(CompoundName(t1)==OPR_SPECIAL && 
					(GetAtomProperty(CompoundArg1(t1),A_COLOR)==A_COLOR_F
					||GetAtomProperty(CompoundArg1(t1),A_COLOR)==A_COLOR_D) &&
					IntegerValue(ListFirst(CompoundArg2(t))) >
					IntegerValue(ListFirst(CompoundArg2(t1)))    )
					{
					List tl1,tl2;
					tl1=ConsumeCompoundArg(t,2);
					tl2=ConsumeCompoundArg(t1,2);
					SetCompoundArg(t,2,tl2);
					SetCompoundArg(t1,2,tl1);
					}
				ll=ListTail(ll);
				}
			}
		l=ListTail(l);
		}
		
	l=CompoundArgN(m2,3);
	while(!is_empty_list(l))
		{
		Term t;
		t=ListFirst(l);	
		if(CompoundName(t)==OPR_SPECIAL && 
			GetAtomProperty(CompoundArg1(t),A_COLOR)==A_COLOR_EPS)
			{
			List ll;
			ll=ListTail(l);
			while(!is_empty_list(ll))
				{
				Term t1;
				t1=ListFirst(ll);
				if(CompoundName(t1)==OPR_SPECIAL && 
					GetAtomProperty(CompoundArg1(t1),A_COLOR)==A_COLOR_EPS &&
					IntegerValue(ListFirst(CompoundArg2(t))) >
					IntegerValue(ListFirst(CompoundArg2(t1)))    )
					{
					List tl1,tl2;
					tl1=ConsumeCompoundArg(t,2);
					tl2=ConsumeCompoundArg(t1,2);
					SetCompoundArg(t,2,tl2);
					SetCompoundArg(t1,2,tl1);
					}
				ll=ListTail(ll);
				}
			}
		l=ListTail(l);
		}
*/		l=CompoundArgN(m2,3);

	while(!is_empty_list(l))
		{
		Term t;
		t=ListFirst(l);	
		if(CompoundName(t)==OPR_SPECIAL && 
			(GetAtomProperty(CompoundArg1(t),A_COLOR)==A_COLOR_F||
				GetAtomProperty(CompoundArg1(t),A_COLOR)==A_COLOR_D||
				GetAtomProperty(CompoundArg1(t),A_COLOR)==A_COLOR_EPS||
				GetAtomProperty(CompoundArg1(t),A_SYMM)))
			{
			List ll;
			ll=ListTail(l);
			while(!is_empty_list(ll))
				{
				Term t1;
				t1=ListFirst(ll);
				if(CompoundName(t1)==OPR_SPECIAL && 
					CompoundArg1(t1)==CompoundArg1(t) &&
					IntegerValue(ListFirst(CompoundArg2(t))) >
					IntegerValue(ListFirst(CompoundArg2(t1)))    )
					{
					List tl1,tl2;
					tl1=ConsumeCompoundArg(t,2);
					tl2=ConsumeCompoundArg(t1,2);
					SetCompoundArg(t,2,tl2);
					SetCompoundArg(t1,2,tl1);
					}
				ll=ListTail(ll);
				}
			}
		l=ListTail(l);
		}
	
	}
	
void epsv_symm(List pl, Term m2)
	{
	int fact=1;
	List l;
	l=CompoundArgN(m2,3);
	while(!is_empty_list(l))
		{
		Term t;
		t=ListFirst(l);
		
		if(CompoundName(t)==OPR_SPECIAL && CompoundArg1(t)==A_EPS_V)
			{
			List ll;
			ll=ConsumeCompoundArg(t,2);
			fact*=f_symm_factor(&ll);
			SetCompoundArg(t,2,ll);
			}
		l=ListTail(l);
		}

	if(fact-1)
		{
		Term t;
		t=CompoundArg1(m2);
		SetCompoundArg(t,1,NewInteger(-IntegerValue(CompoundArg1(t))));
		}		
	
	l=CompoundArgN(m2,3);
	while(!is_empty_list(l))
		{
		Term t;
		t=ListFirst(l);	
		if(CompoundName(t)==OPR_SPECIAL && CompoundArg1(t)==A_EPS_V)
			{
			List ll;
			ll=ListTail(l);
			while(!is_empty_list(ll))
				{
				Term t1;
				t1=ListFirst(ll);
				if(CompoundName(t1)==OPR_SPECIAL && CompoundArg1(t1)==A_EPS_V &&
					IntegerValue(ListFirst(CompoundArg2(t))) >
					IntegerValue(ListFirst(CompoundArg2(t1)))    )
					{
					List tl1,tl2;
					tl1=ConsumeCompoundArg(t,2);
					tl2=ConsumeCompoundArg(t1,2);
					SetCompoundArg(t,2,tl2);
					SetCompoundArg(t1,2,tl1);
					}
				ll=ListTail(ll);
				}
			}
		l=ListTail(l);
		}
		
	}

void color_check_spec(Atom name, List ind)
	{
	int t1, t2, t3, ito=0;
	List l;
	List itc=0;
	
	for(l=ind;l;l=ListTail(l))
	{
		int tp;
		tp=color_repres(ListFirst(l));
		if(tp)
			itc=AppendLast(itc,NewInteger(tp));
		else
			ito=1;
	}
	
	if(itc==0)
		return;
	
	/*printf("%s: ",AtomValue(name));
	WriteTerm(itc);
	puts("");
	*/
	if(ito)
	{
		ErrorInfo(71);
		printf("Warning: special %s breaks CompHEP color conventions.\n",
					AtomValue(name));
		return ;
	}
	
	if(ListLength(itc)!=3)
	{
		ErrorInfo(72);
		printf("Warning: special %s breaks CompHEP color conventions.\n",
					AtomValue(name));
		return ;
	}
	
	t1=(int)IntegerValue(ListFirst(itc));
	t2=(int)IntegerValue(ListFirst(ListTail(itc)));
	t3=(int)IntegerValue(ListFirst(ListTail(ListTail(itc))));
	
	if(t1==1 && t2==2 && t3==3)
		{
		if(ColorLambda!=0)
			{
			printf("Warning: special %s breaks CompHEP color conventions.\n",
				AtomValue(name));
			return ;
			}
		ColorLambda=name;
		SetAtomProperty(name,A_COLOR,A_COLOR_LAMBDA);
		}
		
	if(t1==3 && t2==3 && t3==3)
		{
		if(ColorF!=0)
			{
			/*printf("Warning: special %s breaks CompHEP color conventions.\n",
				AtomValue(name));*/
			if(ColorD!=0)
			{
				printf("Warning: (color) special %s unsupported.\n",
				AtomValue(name));
				return;
			}
			ColorD=name;
			SetAtomProperty(name,A_COLOR,A_COLOR_D);			
			return ;
			}
		ColorF=name;
		SetAtomProperty(name,A_COLOR,A_COLOR_F);
		}
		
	if(t1==1 && t2==1 && t3==1)
		{
		if(ColorEps!=0)
			{
			printf("Warning: special %s breaks CompHEP color conventions.\n",
				AtomValue(name));
			return ;
			}
		ColorEps=name;
		SetAtomProperty(name,A_COLOR,A_COLOR_EPS);
		}
		
	if(t1==2 && t2==2 && t3==2)
		{
		if(ColorEpsB!=0)
			{
			printf("Warning: special %s breaks CompHEP color conventions.\n",
				AtomValue(name));
			return ;
			}
		ColorEpsB=name;
		SetAtomProperty(name,A_COLOR,A_COLOR_EPS);
		}
		
	if(t1==1 && t2==1 &&  t3==7)
	{
	  ColorK6=name;
	  SetAtomProperty(name,A_COLOR,A_COLOR_K6);
	}
	if(t1==2 && t2==2 &&  t3==6)
	{
	  ColorK6b=name;
	  SetAtomProperty(name,A_COLOR,A_COLOR_K6);
	}
	if(t1==6 && t2==7 && t3==3)
	{
	  ColorL6=name;
	  SetAtomProperty(name,A_COLOR,A_COLOR_L6);
	}
	}

extern int check_splitting_c2;
		
Term alg2_rem_lambdaf(Term a2)
{
	Term cf=0, cl=0, cf2=0, cl2=0;
	Label la, la2, la3, i,j;
	Term a2c;
	List l=CompoundArgN(a2,5);
	List fisav=0;
	Atom cfsav;
	List l2;
	
	if(ListLength(CompoundArg1(a2))!=4)
		return 0;
	if(ListLength(l)!=1)
	{
		puts("Internal error (a2remlf1)");
		return 0;
	}
	
	for(l=CompoundArgN(ListFirst(l),3);l;l=ListTail(l))
	{
		Term sp=ListFirst(l);
		if(CompoundName(sp)==OPR_SPECIAL && GetAtomProperty(CompoundArg1(sp),
				A_COLOR)==A_COLOR_F)
		{
			if(cf)
				return 0;
			else
				cf=sp;
		}
		if(CompoundName(sp)==OPR_SPECIAL && GetAtomProperty(CompoundArg1(sp),
				A_COLOR)==A_COLOR_LAMBDA)
		{
			if(cl)
			{
				if(cf==0 && !FAOutput && !TexOutput)
				{
					cf=sp;
					goto llrdc;
				}
				return 0;
			}
			else
				cl=sp;
		}
	}
	if(cf==0 || cl==0)
		return 0;

	i =ListNth(CompoundArg2(cl),1);
	j =ListNth(CompoundArg2(cl),2);
	la=ListNth(CompoundArg2(cl),3);
	
	if(ListNth(CompoundArg2(cf),1)==la)
		la2=ListNth(CompoundArg2(cf),2),la3=ListNth(CompoundArg2(cf),3);
	else if(ListNth(CompoundArg2(cf),2)==la)
		la2=ListNth(CompoundArg2(cf),3),la3=ListNth(CompoundArg2(cf),1);
	else if(ListNth(CompoundArg2(cf),3)==la)
		la2=ListNth(CompoundArg2(cf),1),la3=ListNth(CompoundArg2(cf),2);
	else
		return 0;
	
	fisav=CopyTerm(CompoundArg2(cf));
	cfsav=CompoundArg1(cf);
	
	a2c=CopyTerm(a2);
	for(l=CompoundArgN(ListFirst(CompoundArgN(a2c,5)),3);l;l=ListTail(l))
	{
		Term sp=ListFirst(l);
		if(CompoundName(sp)==OPR_SPECIAL && GetAtomProperty(CompoundArg1(sp),
				A_COLOR)==A_COLOR_F)
			cf2=sp;
		if(CompoundName(sp)==OPR_SPECIAL && GetAtomProperty(CompoundArg1(sp),
				A_COLOR)==A_COLOR_LAMBDA)
			cl2=sp;
	}
	
	/*WriteTerm(a2);
	puts("->");*/
	
		
			
	SetCompoundArg(cf,1,CompoundArg1(cl));
	l=CompoundArg2(cl);
	ChangeList(ListTail(l),la);
	ChangeList(ListTail(ListTail(l)),la2);
	l=CompoundArg2(cf);
	ChangeList(l,la);
	ChangeList(ListTail(l),j);
	ChangeList(ListTail(ListTail(l)),la3);
	
	SetCompoundArg(cf2,1,CompoundArg1(cl2));
	l=CompoundArg2(cl2);
	ChangeList(ListTail(l),la);
	ChangeList(ListTail(ListTail(l)),la3);
	l=CompoundArg2(cf2);
	ChangeList(l,la);
	ChangeList(ListTail(l),j);
	ChangeList(ListTail(ListTail(l)),la2);
	
	if(!FAOutput && !TexOutput)
	{
		check_splitting_c2=1;	
		l=color_reduce(a2);
		if(l!=0 && !is_list(l) /*&& EqualTerms(l,a2)*/)
			l=color_reduce(a2c);
		check_splitting_c2=0;

		if(l==0 || is_list(l) /*|| !EqualTerms(l,a2c)*/) 
		{
			WriteTerm(l);puts("roll back");
			
			FreeAtomic(a2c);
			SetCompoundArg(cf,1,cfsav);
			SetCompoundArg(cf,2,fisav);
			l=CompoundArg2(cl);
			ChangeList(ListTail(l),j);
			ChangeList(ListTail(ListTail(l)),la);
			return 0;
		}
	}
	
	l=CompoundArg1(ListFirst(CompoundArgN(a2,5)));
	SetCompoundArg(l,1,NewInteger(-IntegerValue(CompoundArg1(l))));
	
	/*WriteTerm(a2);puts("");
	WriteTerm(a2c);puts("");*/
	return a2c;
	
llrdc:
	check_splitting_c2=1;	
	l=color_reduce(a2);
	check_splitting_c2=0;
	if(l!=0 && EqualTerms(l,a2))
		return 0;
	
	
	a2c=CopyTerm(a2);
	for(l=CompoundArgN(ListFirst(CompoundArgN(a2c,5)),3);l;l=ListTail(l))
	{
		Term sp=ListFirst(l);
		if(CompoundName(sp)==OPR_SPECIAL && GetAtomProperty(CompoundArg1(sp),
				A_COLOR)==A_COLOR_LAMBDA)
        {
			if(cl2==0)
				cl2=sp;
			else
				cf2=sp;
        }
	}
	
	if(ColorF==0) ColorF=NewAtom("f_SU3",0);
	l=CompoundArg2(cl);
	l2=CompoundArg2(cf);
	
	if(ListNth(l, 2)==ListNth(l2,1))
	{
		i=ListNth(l,1);
		j=ListNth(l2,2);
		la=ListFirst(l2);
		la2=ListNth(l,3);
		la3=ListNth(l2,3);
		ChangeList(ListTail(ListTail(l)),la3);
		ChangeList(ListTail(ListTail(l2)),la2);
	}
	else if(ListNth(l2, 2)==ListNth(l,1))
	{
		i=ListNth(l2,1);
		j=ListNth(l,2);
		la=ListFirst(l);
		la2=ListNth(l2,3);
		la3=ListNth(l,3);
		ChangeList(ListTail(ListTail(l)),la2);
		ChangeList(ListTail(ListTail(l2)),la3);
	}
	else return 0;
	
	l=CompoundArg2(cl2);
	ChangeList(l,i);
	ChangeList(ListTail(l),j);
	ChangeList(ListTail(ListTail(l)),la);
	SetCompoundArg(cf2,1,ColorF);
	l=CompoundArg2(cf2);
	ChangeList(l,la);
	ChangeList(ListTail(l),la2);
	ChangeList(ListTail(ListTail(l)),la3);
	
	l=CompoundArg1(ListFirst(CompoundArgN(a2c,5)));
	SetCompoundArg(l,1,NewInteger(-IntegerValue(CompoundArg1(l))));
	
	
	return a2c;
}

void alg2_rem_col(Term a2)
{
	List l;
	Integer ci[2];
	int cno=0;
	ci[0]=ci[1]=0;
	for(l=CompoundArg1(a2);l;l=ListTail(l))
	{
		Term a=ListFirst(l);
		if(GetAtomProperty(CompoundArg1(a),A_COLOR))
		{
			if(cno>1)
			{
				SetCompoundArg(a2,5,0);
				return;
			}
			if(CompoundName(CompoundArg2(a))==OPR_SCALAR)
				ci[cno]=ListFirst(CompoundArg1(CompoundArg2(a)));
			else
				ci[cno]=ListFirst(ListTail(CompoundArg1(CompoundArg2(a))));
			cno++;
		}
	}
	if(cno==1)
		puts("Internal error colremcol01");
	if(cno!=2)
		return;
	for(l=CompoundArgN(a2,5);l;l=ListTail(l))
	{
		Term m2=ListFirst(l);
		List mt=ConsumeCompoundArg(m2,3);
		List l2;
		for(l2=mt;l2;l2=ListTail(l2))
		{
			if(CompoundArg1(ListFirst(l2))==A_DELTA &&
					ListFirst(CompoundArg2(ListFirst(l2)))==ci[0])
			{
				mt=CutFromList(mt,l2);
				break;
			}
		}
		SetCompoundArg(m2,3,mt);
	}
	
}

