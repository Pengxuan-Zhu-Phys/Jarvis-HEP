#include <ctype.h>
#include <string.h>
#include "lanhep.h"


int AlwaysBracets = 0, WideWriting =0, NoQuotes=0, fortr_digi=0 ;


static int correct_atom(Atom a)
	{
	char *s;
	
	if(NoQuotes)
		return 1;
	
	s=AtomValue(a);
	if(!isalpha(s[0]) && s[0]!='_') return 0;
	s++;
	while(s[0]!=0)
		{
		if(!isalnum(s[0]) && s[0]!='_')
			return 0;
		s++;
		}
	return 1;
	}

int DisplayTerm(Term T)
	{
	char c;
	int arg,arity;
	Functor fu;
	c=AtomicType(T);
	switch(c)
		{
		case 'i':
			return printf("%ld",IntegerValue(T));
		case 'a':
			if(correct_atom(T))
				return printf("%s",AtomValue(T));
			else
				return printf("\'%s\'",AtomValue(T));
		case 'f':
			if(ImaginaryValue(T)==0.0)
			  return printf("%g",FloatValue(T));
			else
			  return printf("(%g,%g)",FloatValue(T),ImaginaryValue(T));
		case 'u':
			return printf("%s/%d",AtomValue(FunctorName(T)),FunctorArity(T));
		case 'e':
			return printf("_");
		case '?':
			return printf("?");
		case 'L':
			return printf("L%d",LabelValue(T));
		case 'l':
			{
			List li;
			int ret=2;
			li=T;
			printf("[");
			while(!is_empty_list(li))
				{
				ret+=DisplayTerm(ListFirst(li));
				li=ListTail(li);
				if(!is_empty_list(li))
					ret+=printf(", ");
				}
			printf("]");
			return ret;
			}
		case 'c':
			{
			int ret;
			fu=CompoundFunctor(T);
			arity=FunctorArity(fu);
			if(correct_atom(FunctorName(fu)))
				ret=printf("%s(",AtomValue(FunctorName(fu)));
			else
				ret=printf("'%s'(",AtomValue(FunctorName(fu)));
			for(arg=1;arg<=arity;arg++)
				{
				if(arg!=1)
					printf(", ");
				ret+=DisplayTerm(CompoundArgN(T,arg));
				}
			return ret+printf(")");
			}			
		}
	return 0;
	}
				
				
				
int fDisplayTerm(FILE *file, Term T)
	{
	char c;
	int arg,arity;
	Functor fu;
	c=AtomicType(T);
	switch(c)
		{
		case 'i':
			return fprintf(file,"%ld",IntegerValue(T));
		case 'a':
			if(correct_atom(T))
				return fprintf(file,"%s",AtomValue(T));
			else
				return fprintf(file,"\'%s\'",AtomValue(T));
		case 'f':
		  if(ImaginaryValue(T)==0.0)
			return fprintf(file,"%g",FloatValue(T));
		  else
			  return fprintf(file,"(%g,%g)",FloatValue(T),ImaginaryValue(T));
		case 'u':
			return fprintf(file,"%s/%d",AtomValue(FunctorName(T)),FunctorArity(T));
		case 'e':
			return fprintf(file,"_");
		case '?':
			return fprintf(file,"?");
		case 'L':
			return fprintf(file,"L%d",LabelValue(T));
		case 'l':
			{
			List li;
			int ret=2;
			li=T;
			fprintf(file,"[");
			while(!is_empty_list(li))
				{
				ret+=fDisplayTerm(file,ListFirst(li));
				li=ListTail(li);
				if(!is_empty_list(li))
					ret+=fprintf(file,", ");
				}
			fprintf(file,"]");
			return ret;
			}
		case 'c':
			{
			int ret;
			fu=CompoundFunctor(T);
			arity=FunctorArity(fu);
			if(correct_atom(FunctorName(fu)))
				ret=fprintf(file,"%s(",AtomValue(FunctorName(fu)));
			else
				ret=fprintf(file,"'%s'(",AtomValue(FunctorName(fu)));
			for(arg=1;arg<=arity;arg++)
				{
				if(arg!=1)
					fprintf(file,", ");
				ret+=fDisplayTerm(file,CompoundArgN(T,arg));
				}
			return ret+fprintf(file,")");
			}			
		}
	return 0;
	}
					


static int arg_prior(Term tt)
	{
	Functor fu;
	int arity, prior;
	Atom name, assoc;
	
	if(AtomicType(tt)!='c')
		return 0;
		
	fu=CompoundFunctor(tt);
	arity=FunctorArity(fu);
	if(arity>2) return 0;
	name=FunctorName(fu);
	if(arity==2)
		{
		if(GetOperator(name,OP_INFIX,&assoc,&prior)==0)
			return 0;
		else
			return prior;
		}
	if(GetOperator(name,OP_PREFIX,&assoc,&prior)!=0)
		return prior;
	if(GetOperator(name,OP_POSTFIX,&assoc,&prior)!=0)
		return prior;
	return 0;
	}
	


static int pri_begin(Atom a)
	{
	char *s;
	if(WideWriting) return 1;
	s=AtomValue(a);
	return  isalpha(s[0]) || isdigit(s[0]) ;
	}
	
static int pri_end(Atom a)
	{
	char *s;
	int len;
	if(WideWriting) return 1;
	s=AtomValue(a);
	len=(int)strlen(s);
	return  isalpha(s[len-1]) || isdigit(s[len-1]) ;
	}
	


int WriteTerm(Term T)
	{
	char c;
	int arg,arity;
	Functor fu;
	c=AtomicType(T);
	switch(c)
		{
		case 'i':
			return printf("%ld",IntegerValue(T));
		case 'a':
			if(correct_atom(T))
				return printf("%s",AtomValue(T));
			else
				return printf("\'%s\'",AtomValue(T));
		case 'f':
		  if(ImaginaryValue(T)==0.0)
			return printf("%g",FloatValue(T));
		  else
			  return printf("(%g,%g)",FloatValue(T),ImaginaryValue(T));
		case 'u':
			return printf("%s/%d",AtomValue(FunctorName(T)),FunctorArity(T));
		case 'e':
			return printf("_");
		case '?':
			return printf("?");
		case 'L':
			return printf("L%d",LabelValue(T));
		case 'l':
			{
			List li;
			int ret=2;
			li=T;
			printf("[");
			while(!is_empty_list(li))
				{
				ret+=WriteTerm(ListFirst(li));
				li=ListTail(li);
				if(!is_empty_list(li))
					ret+=printf(", ");
				}
			printf("]");
			return ret;
			}
		case 'c':
			{
			Atom oper;
			int retval=0,i;
			fu=CompoundFunctor(T);
			arity=FunctorArity(fu);
			oper=FunctorName(fu);
			if(oper==A_FBRACET)
				{
				retval=printf("{");
				for(i=1;i<=arity;i++)
					{
					retval+=WriteTerm(CompoundArgN(T,i));
					if(i!=arity)
						retval+=printf(", ");
					}
				retval+=printf("}");
				return retval;
				}
			if(arity==1)
				{
				Atom assoc;
				int prior, prior1;
				if(oper==A_FBRACET)
					{
					retval=printf("{");
					retval+=WriteTerm(CompoundArg1(T));
					retval+=printf("}");
					return retval;
					}
				if(oper==A_QBRACET)
					{
					retval=printf("[");
					retval+=WriteTerm(CompoundArg1(T));
					retval+=printf("]");
					return retval;
					}
				if(GetOperator(oper,OP_PREFIX,&assoc,&prior)!=0)
					{
					prior1=arg_prior(CompoundArg1(T));
					if(assoc==OP_FX) prior1++;
					retval=printf("%s",AtomValue(oper));
					if(pri_end(oper)) retval+=printf(" ");
					if(prior1>prior || (AlwaysBracets && prior1>1) ||
						(!pri_end(oper)&&is_integer(CompoundArg1(T))&&
							IntegerValue(CompoundArg1(T))<0))
						retval+=printf("(");
					retval+=WriteTerm(CompoundArg1(T));
					if(prior1>prior || (AlwaysBracets && prior1>1)||
						(!pri_end(oper)&&is_integer(CompoundArg1(T))&&
							IntegerValue(CompoundArg1(T))<0))
						retval+=printf(")");
					return retval;
					}
				if(GetOperator(oper,OP_POSTFIX,&assoc,&prior)!=0)
					{
					prior1=arg_prior(CompoundArg1(T));
					if(assoc==OP_XF) prior1++;
					if(prior1>prior || (AlwaysBracets && prior1>1))
						retval+=printf("(");
					retval+=WriteTerm(CompoundArg1(T));
					if(prior1>prior || (AlwaysBracets && prior1>1))
						retval+=printf(")");
					if(pri_begin(oper)) retval+=printf(" ");
					retval+=printf("%s",AtomValue(oper));
					return retval;
					}
				}
			if(arity==2)
				{
				Atom assoc;
				int prior, prior1, prior2;
				if(GetOperator(oper,OP_INFIX,&assoc,&prior)!=0)
					{
					prior1=arg_prior(CompoundArg1(T));
					prior2=arg_prior(CompoundArg2(T));
					if(assoc==OP_XFX || assoc==OP_XFY) prior1++;
					if(assoc==OP_XFX || assoc==OP_YFX) prior2++;
					if(prior1>prior || (AlwaysBracets && prior1>1))
						retval+=printf("(");
					retval+=WriteTerm(CompoundArg1(T));
					if(prior1>prior || (AlwaysBracets && prior1>1))
						retval+=printf(")");
					if(pri_begin(oper)) retval+=printf(" ");
					retval+=printf("%s",AtomValue(oper));
					if(pri_end(oper)) retval+=printf(" ");
					if(prior2>prior || (AlwaysBracets && prior2>1)||
						(!pri_end(oper)&&is_integer(CompoundArg2(T))&&
							IntegerValue(CompoundArg2(T))<0))
						retval+=printf("(");
					retval+=WriteTerm(CompoundArg2(T));
					if(prior2>prior || (AlwaysBracets && prior2>1)||
						(!pri_end(oper)&&is_integer(CompoundArg2(T))&&
							IntegerValue(CompoundArg2(T))<0))
						retval+=printf(")");
					return retval;
					}
				}
			retval+=printf("%s(",AtomValue(FunctorName(fu)));
			for(arg=1;arg<=arity;arg++)
				{
				if(arg!=1)
					retval+=printf(", ");
				retval+=WriteTerm(CompoundArgN(T,arg));
				}
			retval+=printf(")");
			return retval;
			}
			
		}
	return 0;
	}


int fWriteTerm(FILE *file, Term T)
	{
	char c;
	int arg,arity,i;
	Functor fu;
	c=AtomicType(T);
	switch(c)
		{
		case 'i':
			return fprintf(file,"%ld",IntegerValue(T));
		case 'a':
			if(correct_atom(T))
				return fprintf(file,"%s",AtomValue(T));
			else
				return fprintf(file,"\'%s\'",AtomValue(T));
		case 'f':
		  if(ImaginaryValue(T)==0.0)
			return fprintf(file,"%g",FloatValue(T));
		  else
			  return fprintf(file,"(%g,%g)",FloatValue(T),ImaginaryValue(T));
		case 'u':
			return fprintf(file,"%s/%d",AtomValue(FunctorName(T)),FunctorArity(T));
		case 'e':
			return fprintf(file,"_");
		case '?':
			return fprintf(file,"?");
		case 'L':
			return fprintf(file,"L%d",LabelValue(T));
		case 'l':
			{
			List li;
			int ret=2;
			li=T;
			fprintf(file,"[");
			while(!is_empty_list(li))
				{
				if(AlwaysBracets)
					ret+=fprintf(file,"(");
				ret+=fWriteTerm(file,ListFirst(li));
				if(AlwaysBracets)
					ret+=fprintf(file,")");
				li=ListTail(li);
				if(!is_empty_list(li))
					ret+=fprintf(file,", ");
				}
			fprintf(file,"]");
			return ret;
			}
		case 'c':
			{
			Atom oper;
			int retval=0;
			fu=CompoundFunctor(T);
			arity=FunctorArity(fu);
			oper=FunctorName(fu);
			if(oper==A_FBRACET)
				{
				retval=fprintf(file,"{");
				for(i=1;i<=arity;i++)
					{
					retval+=fWriteTerm(file,CompoundArgN(T,i));
					if(i!=arity)
						retval+=fprintf(file,", ");
					}
				retval+=fprintf(file,"}");
				return retval;
				}
			if(arity==1)
				{
				Atom assoc;
				int prior, prior1;
				if(oper==A_FBRACET)
					{
					retval=fprintf(file,"{");
					retval+=fWriteTerm(file,CompoundArg1(T));
					retval+=fprintf(file,"}");
					return retval;
					}
				if(oper==A_QBRACET)
					{
					retval=fprintf(file,"[");
					retval+=fWriteTerm(file,CompoundArg1(T));
					retval+=fprintf(file,"]");
					return retval;
					}
				if(GetOperator(oper,OP_PREFIX,&assoc,&prior)!=0)
					{
					prior1=arg_prior(CompoundArg1(T));
					if(assoc==OP_FX) prior1++;
					retval=fprintf(file,"%s",AtomValue(oper));
					if(pri_end(oper)) retval+=fprintf(file," ");
					if(prior1>prior || (AlwaysBracets && prior1>1)||
						(!pri_end(oper)&&is_integer(CompoundArg1(T))&&
							IntegerValue(CompoundArg1(T))<0))
						retval+=fprintf(file,"(");
					retval+=fWriteTerm(file,CompoundArg1(T));
					if(prior1>prior || (AlwaysBracets && prior1>1)||
						(!pri_end(oper)&&is_integer(CompoundArg1(T))&&
							IntegerValue(CompoundArg1(T))<0))
						retval+=fprintf(file,")");
					return retval;
					}
				if(GetOperator(oper,OP_POSTFIX,&assoc,&prior)!=0)
					{
					prior1=arg_prior(CompoundArg1(T));
					if(assoc==OP_XF) prior1++;
					if(prior1>prior || (AlwaysBracets && prior1>1))
						retval+=fprintf(file,"(");
					retval+=fWriteTerm(file,CompoundArg1(T));
					if(prior1>prior || (AlwaysBracets && prior1>1))
						retval+=fprintf(file,")");
					if(pri_begin(oper)) retval+=fprintf(file," ");
					retval+=fprintf(file,"%s",AtomValue(oper));
					return retval;
					}
				}
			if(arity==2)
				{
				Atom assoc;
				int prior, prior1, prior2;
				if(GetOperator(oper,OP_INFIX,&assoc,&prior)!=0)
					{
					prior1=arg_prior(CompoundArg1(T));
					prior2=arg_prior(CompoundArg2(T));
					if(assoc==OP_XFX || assoc==OP_XFY) prior1++;
					if(assoc==OP_XFX || assoc==OP_YFX) prior2++;
					if(prior1>prior || (AlwaysBracets && prior1>1))
						retval+=fprintf(file,"(");
					retval+=fWriteTerm(file,CompoundArg1(T));
					if(prior1>prior || (AlwaysBracets && prior1>1))
						retval+=fprintf(file,")");
					if(pri_begin(oper)) retval+=fprintf(file," ");
					retval+=fprintf(file,"%s",AtomValue(oper));
					if(pri_end(oper)) retval+=fprintf(file," ");
					if(prior2>prior || (AlwaysBracets && prior2>1)||
						(!pri_end(oper)&&is_integer(CompoundArg2(T))&&
							IntegerValue(CompoundArg2(T))<0))
						retval+=fprintf(file,"(");
					retval+=fWriteTerm(file,CompoundArg2(T));
					if(prior2>prior || (AlwaysBracets && prior2>1)||
						(!pri_end(oper)&&is_integer(CompoundArg2(T))&&
							IntegerValue(CompoundArg2(T))<0))
						retval+=fprintf(file,")");
					return retval;
					}
				}
			retval+=fprintf(file,"%s(",AtomValue(FunctorName(fu)));
			for(arg=1;arg<=arity;arg++)
				{
				if(arg!=1)
					retval+=fprintf(file,", ");
				retval+=fWriteTerm(file,CompoundArgN(T,arg));
				}
			retval+=fprintf(file,")");
			return retval;
			}
			
		}
	return 0;
	}
	

			


int sWriteTerm(char *string, Term T)
	{
	char c;
	int arg,arity,i,fdisv;
	Functor fu;
	c=AtomicType(T);
	switch(c)
		{
		case 'i':
		  if(fortr_digi)
		    return sprintf(string,"%ldD0",IntegerValue(T));
		  else
		    return sprintf(string,"%ld",IntegerValue(T));
		case 'a':
			if(correct_atom(T))
				return sprintf(string,"%s",AtomValue(T));
			else
				return sprintf(string,"\'%s\'",AtomValue(T));
		case 'f':
		  {
		    int ret;
			if(ImaginaryValue(T)==0.0)
			  ret=sprintf(string,"%g",FloatValue(T));
			else
			  ret=sprintf(string,"(%g,%g)",FloatValue(T),ImaginaryValue(T));
		    int i;
		    for(i=0;i<ret;i++)
		      {
			if(string[i]=='e' || string[i]=='E')
			  {
			    string[i]='D';
			    return ret;
			  }
		      }
		    ret+=sprintf(string+ret,"D0");
		    return ret;
		  }
		case 'u':
			return sprintf(string,"%s/%d",AtomValue(FunctorName(T)),FunctorArity(T));
		case 'e':
			return sprintf(string,"_");
		case '?':
			return sprintf(string,"?");
		case 'L':
			return sprintf(string,"L%d",LabelValue(T));
		case 'l':
			{
			List li;
			int ret=1;
			li=T;
			sprintf(string,"[");
			while(!is_empty_list(li))
				{
				if(AlwaysBracets)
					ret+=sprintf(string+ret,"(");
				ret+=sWriteTerm(string+ret,ListFirst(li));
				if(AlwaysBracets)
					ret+=sprintf(string+ret,")");
				li=ListTail(li);
				if(!is_empty_list(li))
					ret+=sprintf(string+ret,", ");
				}
			sprintf(string+ret,"]");
			ret++;
			return ret;
			}
		case 'c':
			{
			Atom oper;
			int retval=0;
			List iarg;
			fu=CompoundFunctor(T);
			arity=FunctorArity(fu);
			oper=FunctorName(fu);
			iarg=GetAtomProperty(oper,A_INTEGER);
			if(oper==A_FBRACET)
				{
				retval=sprintf(string+retval,"{");
				for(i=1;i<=arity;i++)
					{
					retval+=sWriteTerm(string+retval,CompoundArgN(T,i));
					if(i!=arity)
						retval+=sprintf(string+retval,", ");
					}
				retval+=sprintf(string+retval,"}");
				return retval;
				}
			if(arity==1)
				{
				Atom assoc;
				int prior, prior1;
				if(oper==A_FBRACET)
					{
					retval=sprintf(string+retval,"{");
					retval+=sWriteTerm(string+retval,CompoundArg1(T));
					retval+=sprintf(string+retval,"}");
					return retval;
					}
				if(oper==A_QBRACET)
					{
					retval=sprintf(string+retval,"[");
					retval+=sWriteTerm(string+retval,CompoundArg1(T));
					retval+=sprintf(string+retval,"]");
					return retval;
					}
				if(GetOperator(oper,OP_PREFIX,&assoc,&prior)!=0)
					{
					prior1=arg_prior(CompoundArg1(T));
					if(assoc==OP_FX) prior1++;
					retval=sprintf(string+retval,"%s",AtomValue(oper));
					if(pri_end(oper)) retval+=sprintf(string+retval," ");
					if(prior1>prior || (AlwaysBracets && prior1>1)||
						(!pri_end(oper)&&is_integer(CompoundArg1(T))&&
							IntegerValue(CompoundArg1(T))<0))
						retval+=sprintf(string+retval,"(");
					retval+=sWriteTerm(string+retval,CompoundArg1(T));
					if(prior1>prior || (AlwaysBracets && prior1>1)||
						(!pri_end(oper)&&is_integer(CompoundArg1(T))&&
							IntegerValue(CompoundArg1(T))<0))
						retval+=sprintf(string+retval,")");
					return retval;
					}
				if(GetOperator(oper,OP_POSTFIX,&assoc,&prior)!=0)
					{
					prior1=arg_prior(CompoundArg1(T));
					if(assoc==OP_XF) prior1++;
					if(prior1>prior || (AlwaysBracets && prior1>1))
						retval+=sprintf(string+retval,"(");
					retval+=sWriteTerm(string+retval,CompoundArg1(T));
					if(prior1>prior || (AlwaysBracets && prior1>1))
						retval+=sprintf(string+retval,")");
					if(pri_begin(oper)) retval+=sprintf(string+retval," ");
					retval+=sprintf(string+retval,"%s",AtomValue(oper));
					return retval;
					}
				if(fortr_digi && GetAtomProperty(oper,A_CHNAME) && 
								is_integer(CompoundArg1(T)))
					  {
					    retval+=sprintf(string+retval,"%s(%ld)",AtomValue(oper),
							IntegerValue(CompoundArg1(T)));
						return retval;
					  }

				}
			if(arity==2)
				{
				Atom assoc;
				int prior, prior1, prior2;
				if(GetOperator(oper,OP_INFIX,&assoc,&prior)!=0)
					{
					prior1=arg_prior(CompoundArg1(T));
					prior2=arg_prior(CompoundArg2(T));
					if(assoc==OP_XFX || assoc==OP_XFY) prior1++;
					if(assoc==OP_XFX || assoc==OP_YFX) prior2++;
					if(prior1>prior || (AlwaysBracets && prior1>1))
						retval+=sprintf(string+retval,"(");
					retval+=sWriteTerm(string+retval,CompoundArg1(T));
					if(prior1>prior || (AlwaysBracets && prior1>1))
						retval+=sprintf(string+retval,")");
					if(pri_begin(oper)) retval+=sprintf(string+retval," ");
					retval+=sprintf(string+retval,"%s",AtomValue(oper));
					if(pri_end(oper)) retval+=sprintf(string+retval," ");
					if(prior2>prior || (AlwaysBracets && prior2>1)||
						(!pri_end(oper)&&is_integer(CompoundArg2(T))&&
							IntegerValue(CompoundArg2(T))<0))
						retval+=sprintf(string+retval,"(");
					if(fortr_digi && oper==OPR_POW && is_integer(CompoundArg2(T)))
					  {
					    fortr_digi=0;
					    retval+=sWriteTerm(string+retval,CompoundArg2(T));
					    fortr_digi=1;
					  }
					else
					  retval+=sWriteTerm(string+retval,CompoundArg2(T));
					if(prior2>prior || (AlwaysBracets && prior2>1)||
						(!pri_end(oper)&&is_integer(CompoundArg2(T))&&
							IntegerValue(CompoundArg2(T))<0))
						retval+=sprintf(string+retval,")");
					return retval;
					}
				}
			fdisv=fortr_digi;
			/*fortr_digi=0;*/
			retval+=sprintf(string+retval,"%s(",AtomValue(FunctorName(fu)));
			for(arg=1;arg<=arity;arg++)
				{
				if(arg!=1)
					retval+=sprintf(string+retval,", ");
				if(iarg && ListMember(iarg,NewInteger(arg)))
				{
					if(is_integer(CompoundArgN(T,arg)))
						retval+=sprintf(string+retval,"%ld",IntegerValue(CompoundArgN(T,arg)));
					else
					{
						retval+=sprintf(string+retval,"INT(");
						retval+=sWriteTerm(string+retval,CompoundArgN(T,arg));
						retval+=sprintf(string+retval,")");
					}
				}
				else
					retval+=sWriteTerm(string+retval,CompoundArgN(T,arg));
				}
			retval+=sprintf(string+retval,")");
			fortr_digi=fdisv;
			return retval;
			}
			
		}
	return 0;
	}
	

static void dump_term(Term t, int marg)
{
	char c;
	
	WriteBlank(stdout,marg);
	
	c=AtomicType(t);
	switch(c)
		{
		 
		case 'i':
			printf("%ld\n",IntegerValue(t));
			return;
		case 'a':
			printf("%s\n",AtomValue(t));
			return;
		case 'f':
		  if(ImaginaryValue(t)==0.0)
			printf("%g\n",FloatValue(t));
		  else
			  printf("(%g,%g)",FloatValue(t),ImaginaryValue(t));
			return;
		case 'u':
			printf("%s/%d\n",AtomValue(FunctorName(t)),
						FunctorArity(t));
			return;
		case 'e':
			printf("_\n");
			return;
		case '?':
			printf("(?)\n");
			return;
		case 'L':
			printf("L%d\n",LabelValue(t));
			return;
		case 'l':
			{
			List li;
			
			printf("[\n");
			for(li=t;li;li=ListTail(li))
				dump_term(ListFirst(li),marg+4);
			WriteBlank(stdout,marg);
			printf("]\n");
			return;
			}
		case 'c':
			{
			int i,a;
			a=CompoundArity(t);
			
			if(a==1 && is_atomic(CompoundArg1(t)))
			{
				printf("%s(",AtomValue(CompoundName(t)));
				WriteTerm(CompoundArg1(t));
				printf(")\n");
				return;
			}
			if(a==2 && is_atomic(CompoundArg1(t)) && is_atomic(CompoundArg2(t)))
			{
				WriteTerm(CompoundArg1(t));
				printf(" %s ",AtomValue(CompoundName(t)));
				WriteTerm(CompoundArg2(t));
				puts("");
				return;
			}
			
			printf("%s(\n",AtomValue(CompoundName(t)));
			for(i=1;i<=a;i++)
				dump_term(CompoundArgN(t,i),marg+4);
			WriteBlank(stdout,marg);
			printf(")\n");
			return;
			}
		default:
				printf("(unknown type)\n");
	}
			
}
	
void DumpTerm(Term t)
{
	dump_term(t,4);
}
			
