#include <setjmp.h>
#include "lanhep.h"


static jmp_buf parse_jmp_buf;
static char parse_err_mess[80];
int ExitOnError=0;

/*   list of fu(obj, oper, prior, assoc).  */


static Term list_to_term(List);

static Term read_term(Atomic *stop)
	{
	Atomic l_a, c_a, n_a, end;
	Atom funame;
	List list;
	Functor funct;
	int phase=0;
	list=NewList();
	end = *stop;
	l_a=0;
	funame=NewAtom("fu",0);
	funct=NewFunctor(funame,4);
	c_a=ReadAtomic();


	while(1)
		{
		int prior;
		Atom  ass;



		if(phase==0)  /* waiting for object/fx */
			{
			if(c_a==end || c_a==A_RCET || c_a==A_FCET || c_a==A_QCET)
				{
				sprintf(parse_err_mess,"unexpected '%s'.",AtomValue(c_a));
				longjmp(parse_jmp_buf,3);
				exit(-2);
				}
			if(c_a==0)
				{
				sprintf(parse_err_mess,"unexpected end of file.");
				longjmp(parse_jmp_buf,4);
				exit(-2);
				}
			n_a=ReadAtomic();
			UnreadAtomic(n_a);
			if(c_a==OPR_MINUS && is_integer(n_a))
				{
				long a;
				a=IntegerValue(n_a);
				c_a=NewInteger(-a);
				ReadAtomic();
				}
			if(c_a==OPR_MINUS && is_float(n_a))
				{
				double a;
				a=FloatValue(n_a);
				FreeAtomic(n_a);
				c_a=NewFloat(-a);
				ReadAtomic();
				}
			if(is_integer(c_a) || is_float(c_a))
				{
				Term tt;
				tt=NewCompound(funct);
				SetCompoundArg(tt,1,c_a);
				SetCompoundArg(tt,2,OP_NONE);
				SetCompoundArg(tt,3,NewInteger(0));
				list=AppendLast(list,tt);
				goto lab1;
				}
			if( GetOperator(c_a,OP_PREFIX,&ass,&prior)!=0 && /*n_a!=A_RBRA &&*/
				n_a!=end && n_a!=A_RCET && n_a!=A_FCET && n_a!=A_QCET
				&& (n_a==OPR_MINUS || n_a==OPR_PLUS || GetOperator(n_a,OP_INFIX,NULL,NULL)==0))
				{
				Term tt;
				tt=NewCompound(funct);
				SetCompoundArg(tt,1,c_a);
				SetCompoundArg(tt,2,OP_PREFIX);
				SetCompoundArg(tt,3,NewInteger(prior));
				SetCompoundArg(tt,4,ass);
				list=AppendLast(list,tt);
				if(phase==0)
					phase=1;
				else
					phase=0;
				goto lab1;
				}

			if(c_a==A_RBRA)
				{
				Atom stop, newa;
				Term tt;
				stop=A_RCET;
				newa=read_term(&stop);
				if(stop!=A_RCET)
					{
					sprintf(parse_err_mess,"unexpected '%s'.",AtomValue(stop));
					longjmp(parse_jmp_buf,5);
					exit(-2);
					}
				tt=NewCompound(funct);
				SetCompoundArg(tt,1,newa);
				SetCompoundArg(tt,2,OP_NONE);
				SetCompoundArg(tt,3,NewInteger(0));
				list=AppendLast(list,tt);
				goto lab1;
				}

			if(c_a==A_FBRA)
				{
				Atom stop, newa, fu1;
				Term tt,tt1;
				stop=A_FCET;
				newa=read_term(&stop);
				if(stop!=A_FCET)
					{
					sprintf(parse_err_mess,"unexpected '%s'.",AtomValue(stop));
					longjmp(parse_jmp_buf,6);
					exit(-2);
					}
				fu1=NewFunctor(A_FBRACET,1);
				tt1=NewCompound(fu1);
				SetCompoundArg(tt1,1,newa);
				newa=tt1;
				tt=NewCompound(funct);
				SetCompoundArg(tt,1,newa);
				SetCompoundArg(tt,2,OP_NONE);
				SetCompoundArg(tt,3,NewInteger(0));
				list=AppendLast(list,tt);
				goto lab1;
				}

			if(c_a==A_QBRA)
				{
				Atom stop, newa/*, fu1*/;
				Term tt/*,tt1*/;
				stop=A_QCET;
				newa=read_term(&stop);
				if(stop!=A_QCET)
					{
					sprintf(parse_err_mess,"unexpected '%s'.",AtomValue(stop));
					longjmp(parse_jmp_buf,7);
					exit(-2);
					}
				/*fu1=NewFunctor(A_QBRACET,1);
				tt1=NewCompound(fu1);
				SetCompoundArg(tt1,1,newa);
				newa=tt1;*/
				newa=CommaToList(newa);
				tt=NewCompound(funct);
				SetCompoundArg(tt,1,newa);
				SetCompoundArg(tt,2,OP_NONE);
				SetCompoundArg(tt,3,NewInteger(0));
				list=AppendLast(list,tt);
				goto lab1;
				}


			if(n_a==A_RBRA)
			{
				Atom stop;
				Term tt,tt1;
				Functor fu1;
				int arity=0;
				List al;
				al=0;

				l_a=c_a;
				c_a=ReadAtomic();
				stop=A_COMMA;
				while(stop==A_COMMA)
					{ al=AppendFirst(al,read_term(&stop));
					  arity++;
					  if(arity==100)
					  		{ 
					  		sprintf(parse_err_mess,"too many arguments for '%s'.",
								AtomValue(l_a));
							longjmp(parse_jmp_buf,8);
							exit(-2);
							}
					}
				if(stop!=A_RCET)
					{
					sprintf(parse_err_mess,"unexpected '%s'.",
						AtomValue(stop));
					longjmp(parse_jmp_buf,1);
					exit(-2);
					}
				fu1=NewFunctor(l_a,arity);
				tt1=NewCompound(fu1);

				while(al!=0)
					{
					SetCompoundArg(tt1,arity,ListFirst(al));
					al=ListTail(al);
					arity--;
					}

				tt=NewCompound(funct);
				SetCompoundArg(tt,1,tt1);
				SetCompoundArg(tt,2,OP_NONE);
				SetCompoundArg(tt,3,NewInteger(0));
				list=AppendLast(list,tt);
				goto lab1;
			}

			{
				Term tt;
				tt=NewCompound(funct);

				SetCompoundArg(tt,1,c_a);
				SetCompoundArg(tt,2,OP_NONE);
				SetCompoundArg(tt,3,NewInteger(0));
				list=AppendLast(list,tt);
				goto lab1;
			}


			}
		else		/* waiting for end/xfx/xf */
			{
			if(c_a==end || c_a==A_RCET || c_a==A_FCET || c_a==A_QCET)
				{
				*stop=c_a;
				return list_to_term(list);
				}
			if(c_a==0)
				{
				sprintf(parse_err_mess,"unexpected end of file.");
				longjmp(parse_jmp_buf,9);
				exit(-2);
				}
			if(GetOperator(c_a,OP_POSTFIX,&ass,&prior)!=0)
				{
				Term tt;
				tt=NewCompound(funct);
				SetCompoundArg(tt,1,c_a);
				SetCompoundArg(tt,2,OP_POSTFIX);
				SetCompoundArg(tt,3,NewInteger(prior));
				SetCompoundArg(tt,4,ass);
				list=AppendLast(list,tt);
				if(phase==0)
					phase=1;
				else
					phase=0;
				goto lab1;
				}
			if(GetOperator(c_a,OP_INFIX,&ass,&prior)!=0)
				{
				Term tt;
				tt=NewCompound(funct);
				SetCompoundArg(tt,1,c_a);
				SetCompoundArg(tt,2,OP_INFIX);
				SetCompoundArg(tt,3,NewInteger(prior));
				SetCompoundArg(tt,4,ass);
				list=AppendLast(list,tt);
				goto lab1;
				}

			sprintf(parse_err_mess,"operator expected %s <HERE> %s.",
					AtomValue(l_a),AtomValue(c_a));
			longjmp(parse_jmp_buf,10);
			exit(-2);
			}
		lab1:
			l_a=c_a;
			c_a=ReadAtomic();
			if(phase==0)
				phase=1;
			else
				phase=0;
		}
	}

static void rpl_dollar(Term t)
{
	if(is_compound(t))
	{
		int i,a;
		a=CompoundArity(t);
		for(i=1;i<=a;i++)
		{
			Term t1;
			t1=CompoundArgN(t,i);
			if(is_compound(t1) && CompoundArity(t1)==1 && 
					CompoundName(t1)==OPR_DOLLAR && CompoundArg1(t1)==OPR_USCORE)
				SetCompoundArg(t,i,0);
			else
				rpl_dollar(t1);
		}
		return;
	}
	if(is_list(t))
	{
		List l;
		for(l=t;l;l=ListTail(l))
		{
			Term t1;
			t1=ListFirst(l);
			if(is_compound(t1) && CompoundArity(t1)==1 && 
					CompoundName(t1)==OPR_DOLLAR && CompoundArg1(t1)==OPR_USCORE)
			{
				FreeAtomic(t1);
				ChangeList(l,0);
			}
			else
				rpl_dollar(t1);
		}
		return;
	}
}



Term ReadTerm(void)
	{
	int fline;
	int err;
	Term  a;
	Atom point;

	point=A_POINT;
	a=ReadAtomic();
	if(a==0)
		return a;
	UnreadAtomic(a);
	fline=CurrentInputLine();

	err=setjmp(parse_jmp_buf);

	if(err!=0)
		{
		SkipInputLine();
		ErrorInfo(err);
		printf(" %s\n",parse_err_mess);
		if(ExitOnError)
			exit(-1);
		else
			return 1;
		}

	a = read_term(&point);

	if(point!=A_POINT)
		{
		sprintf(parse_err_mess, "expression ended with '%s').",AtomValue(point));
		longjmp(parse_jmp_buf,11);
		}
	rpl_dollar(a);
	return a;
	}

Term ReadTermM(char *line)
{
	Term t;
	SetInputFileM(line);
	t=ReadTerm();
	CloseInputFile();
	return t;
}


static Term list_to_term(List list)
	{
	List l1, lop=0, l2;
	Term ct;
	int minprior=-1;
    Atom  oper_type;
	Atom oper_assoc;
 	Atom oper;
	
	if(ListTail(list)==0)
		{
		Term ret;
		ret=ConsumeCompoundArg(ListFirst(list),1);
		FreeAtomic(list);
		return ret;
		}
		
/*	l1=list;
	while(l1!=0)
		{
		DisplayTerm(ListFirst(l1));
		puts("");
		l1=ListTail(l1);
		}
	puts("------------------------------------------");
*/

	/*  Find min prior */

	l1=list;

	while(l1!=0)
		{
		int pri;
		ct=ListFirst(l1);
		pri=(int)IntegerValue(CompoundArgN(ct,3));
		if(pri>minprior) minprior=pri;
		l1=ListTail(l1);
		}

	l1=list;


	while(l1!=0)
		{
		if(IntegerValue(CompoundArgN(ListFirst(l1),3))==minprior)
			break;
		l1=ListTail(l1);
		}

	oper=CompoundArgN(ListFirst(l1),1);
	oper_type=CompoundArgN(ListFirst(l1),2);
	oper_assoc=CompoundArgN(ListFirst(l1),4);
	l1=ListTail(l1);


	while(l1!=0)
		{
		if(IntegerValue(CompoundArgN(ListFirst(l1),3))==minprior)
			{
			if(oper_assoc==OP_XF || oper_assoc==OP_FX || oper_assoc==OP_XFX)
				{
				sprintf(parse_err_mess,"operators conflict in expression.");
				longjmp(parse_jmp_buf,12);
				exit(-2);
				}
			if(CompoundArgN(ListFirst(l1),4)!=oper_assoc)
				{
				sprintf(parse_err_mess,"operators conflict in expression.");
				longjmp(parse_jmp_buf,12);
				exit(-2);
				}
			}
		l1=ListTail(l1);
		}


	l1=list;

	while(l1!=0)
		{
		if(IntegerValue(CompoundArgN(ListFirst(l1),3))==minprior)
			{
			oper=CompoundArgN(ListFirst(l1),1);
			oper_type=CompoundArgN(ListFirst(l1),2);
			oper_assoc=CompoundArgN(ListFirst(l1),4);
			lop=l1;
			if(oper_assoc==OP_XFY || oper_assoc==OP_XFX ||
				 oper_assoc==OP_FX || oper_assoc==OP_FY)
				break;
			}

		l1=ListTail(l1);
		}

	l1=list;

	l1=ListSplit(l1,lop,&l2);

	FreeAtomic(lop);
    
	if(oper_type==OP_PREFIX)
		{
		Functor ff;
		Term tt;
		if(l2==0)
			{
			sprintf(parse_err_mess,"prefix operator %s without argument.",
				 AtomValue(oper));
			longjmp(parse_jmp_buf,13);
			exit(-2);
			}
		ff=NewFunctor(oper,1);
		tt=NewCompound(ff);
		SetCompoundArg(tt,1,list_to_term(l2));
		if(l1==0)
			return tt;
		else
			{
			Term funame, funct, nterm;
			funame=NewAtom("fu",0);
			funct=NewFunctor(funame,4);
			nterm=NewCompound(funct);
			SetCompoundArg(nterm,1,tt);
			SetCompoundArg(nterm,2,OP_NONE);
			SetCompoundArg(nterm,3,NewInteger(0));
			l1=AppendLast(l1,nterm);
			return list_to_term(l1);
			}
		}

		if(oper_type==OP_POSTFIX)
		{
		Functor ff;
		Term tt;
		if(l1==0)
			{
			sprintf(parse_err_mess,"postfix operator %s without argument.",
				 AtomValue(oper));
			longjmp(parse_jmp_buf,14);
			exit(-2);
			}
		ff=NewFunctor(oper,1);
		tt=NewCompound(ff);
		SetCompoundArg(tt,1,list_to_term(l1));
		if(l2==0)
			return tt;
		else
			{
			Term funame, funct, nterm;
			funame=NewAtom("fu",0);
			funct=NewFunctor(funame,4);
			nterm=NewCompound(funct);
			SetCompoundArg(nterm,1,tt);
			SetCompoundArg(nterm,2,OP_NONE);
			SetCompoundArg(nterm,3,NewInteger(0));
			l2=AppendFirst(l2,nterm);
			return list_to_term(l2);
			}
		}

	if(oper_type==OP_INFIX)
		{
		Functor ff;
		Term tt;
		if(l1==0)
			{
			sprintf(parse_err_mess,"infix operator %s without left argument.",
				 AtomValue(oper));
			longjmp(parse_jmp_buf,15);
			exit(-2);
			}
		if(l2==0)
			{
			sprintf(parse_err_mess,"infix operator %s without right argument.",
				 AtomValue(oper));
			longjmp(parse_jmp_buf,16);
			exit(-2);
			}
		ff=NewFunctor(oper,2);
		tt=NewCompound(ff);
		SetCompoundArg(tt,1,list_to_term(l1));
		SetCompoundArg(tt,2,list_to_term(l2));
		return tt;
		}

	puts("???????????????????????");
	printf("\t\t:   ");
	DisplayTerm(ListFirst(lop));
	puts("");
	if(l1!=0)  list_to_term(l1);
	if(l2!=0)  list_to_term(l2);


	return 0;
	}

static Term read_displ_term(Atomic *stop)
	{
	Atomic l_a, c_a, n_a, end;
	Atom funame;
	List list;
	Functor funct;
	int phase=0;
	list=NewList();
	end = *stop;
	l_a=0;
	funame=NewAtom("fu",0);
	funct=NewFunctor(funame,4);
	c_a=ReadAtomic();


	while(1)
		{

		if(phase==0)  /* waiting for object/fx */
			{
			if(c_a==end || c_a==A_RCET || c_a==A_FCET || c_a==A_QCET)
				{
				sprintf(parse_err_mess,"unexpected '%s'.",AtomValue(c_a));
				longjmp(parse_jmp_buf,17);
				exit(-2);
				}
			if(c_a==0)
				{
				sprintf(parse_err_mess,"unexpected end of file.");
				longjmp(parse_jmp_buf,18);
				exit(-2);
				}
			n_a=ReadAtomic();
			UnreadAtomic(n_a);
			if(c_a==OPR_MINUS && is_integer(n_a))
				{
				long a;
				a=IntegerValue(n_a);
				c_a=NewInteger(-a);
				ReadAtomic();
				}
			if(c_a==OPR_MINUS && is_float(n_a))
				{
				double a;
				a=FloatValue(n_a);
				FreeAtomic(n_a);
				c_a=NewFloat(-a);
				ReadAtomic();
				}
			if(is_integer(c_a) || is_float(c_a))
				{
				list=c_a;
				goto lab1;
				}

			if(c_a==A_RBRA)
				{
				Atom stop, newa;
				WriteTerm(l_a); puts("!!!!!!"); 
				stop=A_RCET;
				newa=read_displ_term(&stop);
				if(stop!=A_RCET)
					{
					sprintf(parse_err_mess,"unexpected '%s'.",AtomValue(stop));
					longjmp(parse_jmp_buf,19);
					exit(-2);
					}
				list=newa;
				goto lab1;
				}

			if(c_a==A_FBRA)
				{
				sprintf(parse_err_mess,"unexpected '{'.");
				longjmp(parse_jmp_buf,20);
				exit(-2);
				}

			if(0)
				{
				Atom stop, newa/*, fu1*/;
				Term tt/*,tt1*/;
				stop=A_QCET;
				newa=read_displ_term(&stop);
				if(stop!=A_QCET)
					{
					sprintf(parse_err_mess,"unexpected '%s'.",AtomValue(stop));
					longjmp(parse_jmp_buf,21);
					exit(-2);
					}
				/*fu1=NewFunctor(A_QBRACET,1);
				tt1=NewCompound(fu1);
				SetCompoundArg(tt1,1,newa);
				newa=tt1;*/
				newa=CommaToList(newa);
				tt=NewCompound(funct);
				SetCompoundArg(tt,1,newa);
				SetCompoundArg(tt,2,OP_NONE);
				SetCompoundArg(tt,3,NewInteger(0));
				list=AppendLast(list,tt);
				goto lab1;
				}


			if(c_a==A_QBRA)
			{
				Atom stop;
				List al;
				al=0;

				stop=A_COMMA;
				while(stop==A_COMMA)
					{ 
					al=AppendLast(al,read_displ_term(&stop));
					}
				if(stop!=A_QCET)
					{
					sprintf(parse_err_mess,"unexpected '%s'.",
						AtomValue(stop));
					longjmp(parse_jmp_buf,22);
					exit(-2);
					}

				list=al;
				goto lab1;
			}

			if(n_a==A_RBRA)
			{
				Atom stop;
				Term tt1;
				Functor fu1;
				int arity=0;
				List al;
				al=0;

				l_a=c_a;
				c_a=ReadAtomic();
				stop=A_COMMA;
				while(stop==A_COMMA)
					{
					al=AppendFirst(al,read_displ_term(&stop));
					  arity++;
					  if(arity==100)
					  		{ 
					  		sprintf(parse_err_mess,"too many arguments for '%s'.",
								AtomValue(l_a));
							longjmp(parse_jmp_buf,23);
							exit(-2);
							}
					}
				if(stop!=A_RCET)
					{
					sprintf(parse_err_mess,"unexpected '%s'.",
						AtomValue(stop));
					longjmp(parse_jmp_buf,24);
					exit(-2);
					}
				fu1=NewFunctor(l_a,arity);
				tt1=NewCompound(fu1);

				while(al!=0)
					{
					SetCompoundArg(tt1,arity,ListFirst(al));
					al=ListTail(al);
					arity--;
					}

				list=tt1;
				goto lab1;
			}

			list=c_a;
			goto lab1;


			}
		else		/* waiting for end/xfx/xf */
			{
			if(c_a==end || c_a==A_RCET || c_a==A_FCET || c_a==A_QCET)
				{
				*stop=c_a;
				return list;
				}
			if(c_a==0)
				{
				sprintf(parse_err_mess,"unexpected end of file.");
				longjmp(parse_jmp_buf,25);
				exit(-2);
				}
			sprintf(parse_err_mess,"operator expected %s <HERE> %s.",
					AtomValue(l_a),AtomValue(c_a));
			longjmp(parse_jmp_buf,25);
			exit(-2);
			}
		lab1:
			l_a=c_a;
			c_a=ReadAtomic();
			if(phase==0)
				phase=1;
			else
				phase=0;
		}
	}



Term ReadDisplTerm(void)
	{
	int fline, lline;
	int err;
	Term  a;
	Atom point;

	point=A_POINT;
	a=ReadAtomic();
	if(a==0)
		return a;
	UnreadAtomic(a);
	fline=CurrentInputLine();

	err=setjmp(parse_jmp_buf);

	if(err!=0)
		{
		SkipInputLine();
		lline=CurrentInputLine();
		if(IsTermInput()==0)
        {
			if(fline==lline)
				printf("File \"%s\", line %d: ",CurrentInputFile(),fline);
			else
				printf("File \"%s\", lines %d - %d : ",CurrentInputFile(),
					fline,lline);
        }
		printf("Syntax error: %s\n",parse_err_mess);
		if(ExitOnError)
			exit(-1);
		else
			return 1;
		}

	a = read_displ_term(&point);

	if(point!=A_POINT)
		{
		sprintf(parse_err_mess, "expression ended with '%s').",AtomValue(point));
		longjmp(parse_jmp_buf,26);
		}
	ReplAtom(a,OPR_USCORE,0);
	return a;
	}

