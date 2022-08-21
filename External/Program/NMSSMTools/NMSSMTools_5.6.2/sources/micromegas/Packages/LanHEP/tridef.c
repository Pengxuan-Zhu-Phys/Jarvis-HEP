#include <string.h>
#include <stdio.h>
#include "lanhep.h"

int tri_mcmp(Term, Term);
extern int tri_dbg_mode;
				
extern List tri_si_co_list, tri_sc_2_list;

/* tri_si_co list structure:
     1 sin
	 2 cos
	 3 sin 2a
	 4 bias #
	 5 cos 2a
	 6 tan 
	 7 ang tex name
*/

static Term check_sine(Atom nm, Term val)
	{
	Term t;
	if(is_compound(val) && CompoundArity(val)==1 
			&& CompoundName(val)==OPR_MINUS)
		val=CompoundArg1(val);
	
	if(!(is_compound(val) && CompoundName(val)==A_SQRT))
		return 0;
	t=CompoundArg1(val);
	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==OPR_MINUS &&
		CompoundArg1(t)==NewInteger(1))
			{
			Term t1;
			t1=CompoundArg2(t);
			if(is_compound(t1) && CompoundArity(t1)==2 && 
				CompoundName(t1)==OPR_POW && CompoundArg2(t1)==NewInteger(2) &&
				is_atom(CompoundArg1(t1)))
					{
					SetAtomProperty(nm,A_COS,CompoundArg1(t1));
					SetAtomProperty(CompoundArg1(t1),A_SIN,nm);
					return CompoundArg1(t1);
					}
			}
	return 0;
	}
	

static List gte_used = 0;

static int good_tr_expr(Term t)
{
	if(is_integer(t))
		return 1;
	if(is_atom(t))
	{
		List l;
		for(l=tri_si_co_list;l;l=ListTail(l))
			if(t==CompoundArg1(ListFirst(l)) || t==CompoundArg2(ListFirst(l)))
			{
				if(!ListMember(gte_used,t))
					gte_used=AppendLast(gte_used,t);
				return 1;
			}
		return 0;
	}
	
	if(!is_compound(t))
		return 0;
	if(CompoundArity(t)==1 && CompoundName(t)==OPR_MINUS)
		return good_tr_expr(CompoundArg1(t));
	if(CompoundArity(t)!=2)
		return 0;
	if(CompoundName(t)==OPR_POW)
		return good_tr_expr(CompoundArg1(t)) && is_integer(CompoundArg2(t));
	if(CompoundName(t)==OPR_PLUS)
		return good_tr_expr(CompoundArg1(t)) && good_tr_expr(CompoundArg2(t));
	if(CompoundName(t)==OPR_MINUS)
		return good_tr_expr(CompoundArg1(t)) && good_tr_expr(CompoundArg2(t));
	if(CompoundName(t)==OPR_MLT)
		return good_tr_expr(CompoundArg1(t)) && good_tr_expr(CompoundArg2(t));
	return 0;
}



void tri_reg_prm(Atom p, Term val)
{
	
	Atom c1, s1, c2, s2;
	Atom a_11, a_12, a_21, a_22;
	int g1, g2;
	
	List l;
	Term valsv, p1, p2;
	
	int type;
	
	p1=check_sine(p,val);
	if(p1)
	{
		Term ns;
		ns=MakeCompound(A_I,8);
		SetCompoundArg(ns,1,p1);
		SetCompoundArg(ns,2,p);
		tri_si_co_list=AppendLast(tri_si_co_list,ns);
		ns=MakeCompound2(OPR_MINUS,NewInteger(1),
				MakeCompound2(OPR_POW,p1,NewInteger(2)));
		ns=MakeCompound2(OPR_EQSIGN,ns,MakeCompound2(OPR_POW,p,NewInteger(2)));
		ProcSetAngle(MakeCompound1(A_I,ns),0);
		return;
	}
	
	gte_used=0;
	
	if(is_compound(val) && CompoundArity(val)==2 && CompoundName(val)==OPR_DIV
			&& is_atom(CompoundArg1(val)) && is_atom(CompoundArg2(val)))
	{
		Atom s,c;
		s=CompoundArg1(val);
		c=CompoundArg2(val);
		for(l=tri_si_co_list;l;l=ListTail(l))
			if(CompoundArg1(ListFirst(l))==s && CompoundArg2(ListFirst(l))==c)
			{
				char buf[32];
				Atom nn;
				SetCompoundArg(ListFirst(l),6,p);
				if(CompoundArgN(ListFirst(l),7) && !GetAtomProperty(p,A_TEXNAME))
				{
					sprintf(buf,"t_{%s}",AtomValue(CompoundArgN(ListFirst(l),7)));
					nn=NewAtom(buf,0);
					SetAtomProperty(p,A_TEXNAME,nn);
				}
				return;
			}
		return;
	}
	
	if(!good_tr_expr(val))
	{
		if(gte_used)
			FreeAtomic(gte_used);
		return;
	}
	
	if(gte_used==0 || (ListLength(gte_used)>2 && ListLength(gte_used)!=4))
		return;
	
	valsv=val;
	val=ExprTo1(CopyTerm(val));
	if(val==0)
		return ;
	
	val=Alg1to2(val);
	if(val==0)
		return ;
	
	val=alg2_add(NewList(),val);


	for(l=val;l;l=ListTail(l))
		{
        Term a2,nn,ll;
		int ncf;
		a2=ListFirst(l);
		alg2_common_n(a2);
        nn=CompoundArg2(a2);

		ncf=(int)IntegerValue(CompoundArg1(nn));

        for(ll=CompoundArgN(a2,5);ll;ll=ListTail(ll))
           {
           int n;
           n=(int)IntegerValue(CompoundArg1(ListFirst(ll)));
           SetCompoundArg(ListFirst(ll),1,NewInteger(n*ncf));
           }
		   
		alg2_red_cos(a2);
		}
	l=ConsumeCompoundArg(ListFirst(val),5);
	FreeAtomic(val);
	val=l;

	if(is_empty_list(val))
		{
		ErrorInfo(568);
		puts("empty expression in parameter statement");
		return;
		}
	val=SortedList(val,tri_mcmp);
	
	valsv=MakeCompound2(OPR_EQSIGN,CopyTerm(valsv),p);
	ProcSetAngle(MakeCompound1(A_I,valsv),0);
	
	
	if(ListLength(gte_used)<=2)
	{
		List l1;

		for(l=tri_si_co_list;l;l=ListTail(l))
			if(ListFirst(gte_used)==CompoundArg1(ListFirst(l)) || 
					ListFirst(gte_used)==CompoundArg2(ListFirst(l)))
				break;
		
		if(is_empty_list(l))
			goto xyz;

				
		if(ListLength(val)!=2 || CompoundArg2(ListFirst(val)) ||
				ListLength(CompoundArg2(ListFirst(ListTail(val))))!=1)
			goto xyz;
		
		l1=CompoundArg2(ListFirst(ListTail(val)));
		if(CompoundArg1(ListFirst(l1))!=CompoundArg1(ListFirst(l)) ||
				CompoundArg2(ListFirst(l1))!=NewInteger(2))
			goto xyz;
		
		
		
		if(CompoundArg1(ListFirst(val))==NewInteger(1) &&
				CompoundArg1(ListFirst(ListTail(val)))==NewInteger(-2))
			SetCompoundArg(ListFirst(l),5,MakeCompound1(OPR_PLUS,p));
		
		if(CompoundArg1(ListFirst(val))==NewInteger(-1) &&
				CompoundArg1(ListFirst(ListTail(val)))==NewInteger(2))
			SetCompoundArg(ListFirst(l),5,MakeCompound1(OPR_MINUS,p));
		
		if(CompoundArgN(ListFirst(l),5) && CompoundArgN(ListFirst(l),7) &&
				CompoundArg1(CompoundArgN(ListFirst(l),5))==p &&
				GetAtomProperty(p,A_TEXNAME)==0)
		{
			char buf[32];
			Atom nn;
			sprintf(buf,"c_{2%s}",AtomValue(CompoundArgN(ListFirst(l),7)));
			nn=NewAtom(buf,0);
			SetAtomProperty(p,A_TEXNAME,nn);
		}
		
	}
			
xyz:			
	
	if(ListLength(gte_used)==2)
	{
		List l1;
		c1=ListFirst(gte_used);
		s1=ListFirst(ListTail(gte_used));

		for(l=tri_si_co_list;l;l=ListTail(l))
			if(
			(c1==CompoundArg1(ListFirst(l)) && s1==CompoundArg2(ListFirst(l))) ||
			(s1==CompoundArg1(ListFirst(l)) && c1==CompoundArg2(ListFirst(l))))
				break;
		if(is_empty_list(l))
		{
			FreeAtomic(val);
			return;
		}

		if(ListLength(val)!=1 || CompoundArg1(ListFirst(val))!=NewInteger(2))
		{
			FreeAtomic(val);
			return;
		}

		l1=CompoundArg2(ListFirst(val));
		if(ListLength(l1)!=2)
		{
			FreeAtomic(val);
			return;
		}

		if(
				(CompoundArg1(ListFirst(l1))!=c1 || 
					CompoundArg1(ListFirst(ListTail(l1)))!=s1) 
				&&
				(CompoundArg1(ListFirst(l1))!=s1 || 
					CompoundArg1(ListFirst(ListTail(l1)))!=c1))
		{
			FreeAtomic(val);
			return;
		}

		if(CompoundArg2(ListFirst(l1))!=NewInteger(1) ||
				CompoundArg2(ListFirst(ListTail(l1)))!=NewInteger(1))
		{
			FreeAtomic(val);
			return;
		}

		SetCompoundArg(ListFirst(l),3,p);
		if(!GetAtomProperty(p,A_TEXNAME) && CompoundArgN(ListFirst(l),7))
		{
			char buf[32];
			Atom nn;
			sprintf(buf,"s_{2%s}",AtomValue(CompoundArgN(ListFirst(l),7)));
			nn=NewAtom(buf,0);
			SetAtomProperty(p,A_TEXNAME,nn);
		}
		FreeAtomic(val);
		return;
	}
	
	if(ListLength(val)!=2)
	{
		FreeAtomic(val);
		return;
	}
	
	if(ListLength(CompoundArg2(ListFirst(val)))!=2 ||
			ListLength(CompoundArg2(ListFirst(ListTail(val))))!=2)
	{
		FreeAtomic(val);
		return;
	}
	
	a_11=ListFirst(CompoundArg2(ListFirst(val)));
	a_12=ListFirst(ListTail(CompoundArg2(ListFirst(val))));
	a_21=ListFirst(CompoundArg2(ListFirst(ListTail(val))));
	a_22=ListFirst(ListTail(CompoundArg2(ListFirst(ListTail(val)))));
	
	if(CompoundArg2(a_11)!=NewInteger(1) || 
			CompoundArg2(a_12)!=NewInteger(1) ||
			CompoundArg2(a_21)!=NewInteger(1) ||
			CompoundArg2(a_22)!=NewInteger(1))
	{
		FreeAtomic(val);
		return;
	}
	
	a_11=CompoundArg1(a_11);
	a_12=CompoundArg1(a_12);
	a_21=CompoundArg1(a_21);
	a_22=CompoundArg1(a_22);
	
	g1=(int)IntegerValue(CompoundArg1(ListFirst(val)));
	g2=(int)IntegerValue(CompoundArg1(ListFirst(ListTail(val))));
	
	if((g1!=1 && g1!=-1) || (g2!=1 && g2!=-1))
	{
		FreeAtomic(val);
		return;
	}
	
	p1=p2=0;
	
	for(l=tri_si_co_list;l;l=ListTail(l))
	{
		if(ListMember(gte_used,CompoundArg1(ListFirst(l))) &&
				ListMember(gte_used,CompoundArg2(ListFirst(l))))
		{
			if(p1==0)
				p1=ListFirst(l);
			else
				p2=ListFirst(l);
		}
	}
	
	if(p1==0 || p2==0)
	{
		FreeAtomic(val);
		return;
	}
	
	s1=CompoundArg1(p1);
	s2=CompoundArg1(p2);
	c1=CompoundArg2(p1);
	c2=CompoundArg2(p2);
/*	
	printf("IIII ");WriteTerm(val);puts("");
	WriteTerm(s1);printf("/");WriteTerm(c1);printf(" & ");
	WriteTerm(s2);printf("/");WriteTerm(c2);printf(" ; %d",g1);
	WriteTerm(a_11);printf("*");WriteTerm(a_12);printf("  %+d ",g2);
	WriteTerm(a_21);printf("*");WriteTerm(a_22);printf("\n");
*/
			
	type=0;
	
	if( ((a_11==s1 && a_12==c2) || (a_12==s1 && a_11==c2)) &&
		((a_21==c1 && a_22==s2) || (a_22==c1 && a_21==s2)) )
	{
		if(g1==1 && g2==1)
			type=1;
		if(g1==1 && g2==-1)
			type=2;
		if(g1==-1 && g2==-1)
			type=-1;
		if(g1==-1 && g2==1)
			type=-2;
	}
	
	if( ((a_21==s1 && a_22==c2) || (a_22==s1 && a_21==c2)) &&
		((a_11==c1 && a_12==s2) || (a_12==c1 && a_11==s2)) )
	{
		if(g1==1 && g2==1)
			type=1;
		if(g1==1 && g2==-1)
			type=-2;
		if(g1==-1 && g2==-1)
			type=-1;
		if(g1==-1 && g2==1)
			type=2;
	}
	
	if( ((a_11==c1 && a_12==c2) || (a_12==c1 && a_11==c2)) &&
		((a_21==s1 && a_22==s2) || (a_22==s1 && a_21==s2)) )
	{
		if(g1==1 && g2==1)
			type=4;
		if(g1==1 && g2==-1)
			type=3;
		if(g1==-1 && g2==-1)
			type=-4;
		if(g1==-1 && g2==1)
			type=-3;
	}
	
	if( ((a_21==c1 && a_22==c2) || (a_22==c1 && a_21==c2)) &&
		((a_11==s1 && a_12==s2) || (a_12==s1 && a_11==s2)) )
	{
		if(g1==1 && g2==1)
			type=-4;
		if(g1==1 && g2==-1)
			type=-3;
		if(g1==-1 && g2==-1)
			type=4;
		if(g1==-1 && g2==1)
			type=3;
	}
	
	g1=ListMember(tri_si_co_list,p1);
	g2=ListMember(tri_si_co_list,p2);
	
	FreeAtomic(val);
	
	if(g1==0 || g2==0 || type==0)
		return;
	
	for(l=tri_sc_2_list;l;l=ListTail(l))
		if(CompoundArg1(ListFirst(l))==NewInteger(g1) ||
				CompoundArg2(ListFirst(l))==NewInteger(g1) ||
				CompoundArg1(ListFirst(l))==NewInteger(g2) ||
				CompoundArg2(ListFirst(l))==NewInteger(g2) )
			break;
	if(is_empty_list(l))
	{
		Term t;
		int s;
		t=MakeCompound(OPR_POW,7);
		SetCompoundArg(t,1,NewInteger(g1));
		SetCompoundArg(t,2,NewInteger(g2));
		SetCompoundArg(t,3,NewInteger(1));
		tri_sc_2_list=AppendLast(tri_sc_2_list,t);
		SetCompoundArg(p1,4,NewInteger(ListLength(tri_sc_2_list)));
		SetCompoundArg(p2,4,NewInteger(ListLength(tri_sc_2_list)));
		if(type>0)
			s=1;
		else
			s=-1,type=-type;
		SetCompoundArg(t,type+3,MakeCompound1(s==1?OPR_PLUS:OPR_MINUS,p));
		if(CompoundArgN(p1,7) && CompoundArgN(p2,7) && !GetAtomProperty(p,A_TEXNAME))
		{
			char buf[32];
			Atom nn=0;
			if(type==1 && s==1)
			{
				sprintf(buf,"s_{%s+%s}",AtomValue(CompoundArgN(p1,7)),
						AtomValue(CompoundArgN(p2,7)));
				nn=NewAtom(buf,0);
			}
			if(type==2 && s==1)
			{
				sprintf(buf,"s_{%s-%s}",AtomValue(CompoundArgN(p1,7)),
						AtomValue(CompoundArgN(p2,7)));
				nn=NewAtom(buf,0);
			}
			if(type==2 && s==-1)
			{
				sprintf(buf,"s_{%s-%s}",AtomValue(CompoundArgN(p2,7)),
						AtomValue(CompoundArgN(p1,7)));
				nn=NewAtom(buf,0);
			}
			if(type==3 && s==1)
			{
				sprintf(buf,"c_{%s+%s}",AtomValue(CompoundArgN(p1,7)),
						AtomValue(CompoundArgN(p2,7)));
				nn=NewAtom(buf,0);
			}
			if(type==4 && s==1)
			{
				sprintf(buf,"c_{%s-%s}",AtomValue(CompoundArgN(p1,7)),
						AtomValue(CompoundArgN(p2,7)));
				nn=NewAtom(buf,0);
			}
			if(nn)
			{
				SetAtomProperty(p,A_TEXNAME,nn);
				SetAtomProperty(p,A_TEXLENGTH,NewInteger(2));
			}
		}
		return;
	}
	
	if(CompoundArg1(ListFirst(l))!=NewInteger(g1) ||
			CompoundArg2(ListFirst(l))!=NewInteger(g2))
		return;
	
	{
		Term t;
		int s;
		t=ListFirst(l);
		if(type>0)
			s=1;
		else
			s=-1,type=-type;
		if(CompoundArgN(t,type+3))
			return;
		SetCompoundArg(t,type+3,MakeCompound1(s==1?OPR_PLUS:OPR_MINUS,p));
		SetCompoundArg(t,3,NewInteger(1+IntegerValue(CompoundArgN(t,3))));
		if(CompoundArgN(p1,7) && CompoundArgN(p2,7) && !GetAtomProperty(p,A_TEXNAME))
		{
			char buf[32];
			Atom nn=0;
			if(type==1 && s==1)
			{
				sprintf(buf,"s_{%s+%s}",AtomValue(CompoundArgN(p1,7)),
						AtomValue(CompoundArgN(p2,7)));
				nn=NewAtom(buf,0);
			}
			if(type==2 && s==1)
			{
				sprintf(buf,"s_{%s-%s}",AtomValue(CompoundArgN(p1,7)),
						AtomValue(CompoundArgN(p2,7)));
				nn=NewAtom(buf,0);
			}
			if(type==2 && s==-1)
			{
				sprintf(buf,"s_{%s-%s}",AtomValue(CompoundArgN(p2,7)),
						AtomValue(CompoundArgN(p1,7)));
				nn=NewAtom(buf,0);
			}
			if(type==3 && s==1)
			{
				sprintf(buf,"c_{%s+%s}",AtomValue(CompoundArgN(p1,7)),
						AtomValue(CompoundArgN(p2,7)));
				nn=NewAtom(buf,0);
			}
			if(type==4 && s==1)
			{
				sprintf(buf,"c_{%s-%s}",AtomValue(CompoundArgN(p1,7)),
						AtomValue(CompoundArgN(p2,7)));
				nn=NewAtom(buf,0);
			}
			if(nn)
				SetAtomProperty(p,A_TEXNAME,nn);
		}
	}
	
	
	/*printf("Recognized type is %d (%d/%d)\n",type,g1,g2);*/
	
	
}

Term ProcAngle(Term t, Term ind)
{
	Term t1;
	List l1;
	int i;
	
	Atom cos=0, sin=0, tan=0, cos2=0, sin2=0, tex=0;
	
	if(!is_compound(t) || CompoundArity(t)!=1)
	{
		ErrorInfo(591);
		puts("wrong syntax in angle statement");
		FreeAtomic(t);
		return 0;
	}
	
	t1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	t=t1;
	
	t=CommaToList(t);
	
	for(l1=t;l1;l1=ListTail(l1))
	{
		char *opt;
		t1=ListFirst(l1);
		if(!is_compound(t1) || CompoundArity(t1)!=2 || 
				CompoundName(t1)!=OPR_EQSIGN || !is_atom(CompoundArg1(t1)) ||
				!is_atom(CompoundArg2(t1)))
		{
			ErrorInfo(584);
			printf("angle: illegal construction '");WriteTerm(t1);puts("'");
			return 0;
		}
		opt=AtomValue(CompoundArg1(t1));
		
		if(strcmp(opt,"texname")==0)
		{
			tex=CompoundArg2(t1);
			continue;
		}
		
		if(!is_parameter(CompoundArg2(t1)))
		{
			ErrorInfo(588);
			printf("angle: unknown parameter '%s'.\n",
					AtomValue(CompoundArg2(t1)));
			return 0;
		}
		
		if(strcmp(opt,"cos")==0)
		{
			cos=CompoundArg2(t1);
			continue;
		}
		if(strcmp(opt,"sin")==0)
		{
			sin=CompoundArg2(t1);
			continue;
		}
		if(strcmp(opt,"cos2")==0)
		{
			cos2=CompoundArg2(t1);
			continue;
		}
		if(strcmp(opt,"sin2")==0)
		{
			sin2=CompoundArg2(t1);
			continue;
		}
		if(strcmp(opt,"tan")==0)
		{
			tan=CompoundArg2(t1);
			continue;
		}
		
		ErrorInfo(584);
		printf("angle: illegal construction '");WriteTerm(t1);puts("'");
		return 0;
	}
	
	if(cos==0 || sin==0)
	{
		ErrorInfo(584);
		printf("angle: definition for sinus and cosine must present");
		return 0;
	}
	
	for(l1=tri_si_co_list;l1;l1=ListTail(l1))
		if((CompoundArg1(ListFirst(l1))==sin && 
					CompoundArg2(ListFirst(l1))==cos) ||
			(CompoundArg1(ListFirst(l1))==cos && 
					CompoundArg2(ListFirst(l1))==sin))
			break;
	if(is_empty_list(l1))
	{
		Term tt;
		tt=MakeCompound(A_I,8);
		tri_si_co_list=AppendLast(tri_si_co_list,tt);
		for(l1=tri_si_co_list;ListTail(l1);l1=ListTail(l1));
	}

	for(i=3;i<=8;i++)
		if(CompoundArgN(ListFirst(l1),i))
		{
			ErrorInfo(585);
			printf("angle statement should be next to declaring parameters %s,%s\n",
					AtomValue(sin),AtomValue(cos));
			return 0;
		}
	
	t1=ListFirst(l1);
	SetCompoundArg(t1,1,sin);
	SetCompoundArg(t1,2,cos);
	SetCompoundArg(t1,3,sin2);
	if(cos2)
		SetCompoundArg(t1,5,MakeCompound1(OPR_PLUS,cos2));
	SetCompoundArg(t1,6,tan);
	SetCompoundArg(t1,7,tex);
	
	if(GetAtomProperty(cos,A_COS))
		RemoveAtomProperty(cos,A_COS);
	if(GetAtomProperty(cos,A_SIN))
		RemoveAtomProperty(cos,A_SIN);
	if(GetAtomProperty(sin,A_COS))
		RemoveAtomProperty(sin,A_COS);
	if(GetAtomProperty(sin,A_SIN))
		RemoveAtomProperty(sin,A_SIN);
	
	
	SetAtomProperty(sin,A_SIN,cos);
	SetAtomProperty(cos,A_COS,sin);

	if(sin2)
	{
		Term t;
		t=MakeCompound2(OPR_MLT,MakeCompound2(OPR_MLT,NewInteger(2),sin),cos);
		t=MakeCompound2(OPR_EQSIGN,t,sin2);
		ProcSetAngle(MakeCompound1(A_I,t),0);
	}
	
	if(cos2)
	{
		Term t;
		t=MakeCompound2(OPR_MINUS,MakeCompound2(OPR_POW,cos,NewInteger(2)),
				MakeCompound2(OPR_POW,sin,NewInteger(2)));
		t=MakeCompound2(OPR_EQSIGN,t,cos2);
		ProcSetAngle(MakeCompound1(A_I,t),0);
	}

		
	if(tex)
	{
		char buf[32];
		Atom nn;
		if(!GetAtomProperty(cos,A_TEXNAME))
		{
			sprintf(buf,"c_{%s}",AtomValue(tex));
			nn=NewAtom(buf,0);
			SetAtomProperty(cos,A_TEXNAME,nn);
		}
		if(!GetAtomProperty(sin,A_TEXNAME))
		{
			sprintf(buf,"s_{%s}",AtomValue(tex));
			nn=NewAtom(buf,0);
			SetAtomProperty(sin,A_TEXNAME,nn);
		}
		if(cos2 && !GetAtomProperty(cos2,A_TEXNAME))
		{
			sprintf(buf,"c_{2%s}",AtomValue(tex));
			nn=NewAtom(buf,0);
			SetAtomProperty(cos2,A_TEXNAME,nn);
		}
		if(sin2 && !GetAtomProperty(sin2,A_TEXNAME))
		{
			sprintf(buf,"s_{2%s}",AtomValue(tex));
			nn=NewAtom(buf,0);
			SetAtomProperty(sin2,A_TEXNAME,nn);
		}
		if(tan && !GetAtomProperty(tan,A_TEXNAME))
		{
			sprintf(buf,"t_{%s}",AtomValue(tex));
			nn=NewAtom(buf,0);
			SetAtomProperty(tan,A_TEXNAME,nn);
		}
	}
	
	FreeAtomic(t);
	return 0;
}

Term ProcDbgTrig(Term t, Term ind)
{
	if(is_compound(t) && is_integer(CompoundArg1(t)))
	{
		tri_dbg_mode=(int)IntegerValue(CompoundArg1(t));
		FreeAtomic(t);
		return 0;
	}
	DumpList(tri_si_co_list);
	DumpList(tri_sc_2_list);
	return 0;
}
