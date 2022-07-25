#include <stdio.h>
#include <string.h>
		
#include "lanhep.h"

extern Atom ColorLambda, ColorF, ColorEps, ColorEpsB;
	
extern int check_funcs;

static void rpl_caret(char *cbuf)
	{
	int i,j,len;
	len=strlen(cbuf);
	for(i=0;cbuf[i];i++)
		{
		if(cbuf[i]!='^')
			continue;
		for(j=len;j>i;j--)
			cbuf[j]=cbuf[j-1];
		cbuf[len+1]=0;
		len++;
		cbuf[i]='*';
		cbuf[i+1]='*';
		}
	}



static int EqualTermsF(Term t1, Term t2)
	{
	if(t1==t2)
		return 1;
	if(is_compound(t1) && is_compound(t2) && 
		CompoundFunctor(t1) == CompoundFunctor(t2))
		{
		int i,a;
		a=CompoundArity(t1);
		for(i=1;i<=a;i++)	
			if(!EqualTermsF(CompoundArgN(t1,i),CompoundArgN(t2,i)))
				return 0;
		return 1;
		}
	if(is_list(t1) && is_list(t2) && ListLength(t1)==ListLength(t2))
		{
		List l1,l2;
		l1=t1;
		l2=t2;
		while(!is_empty_list(l1))
			{
			if(!EqualTermsF(ListFirst(l1),ListFirst(l2)))
				return 0;
			l1=ListTail(l1);
			l2=ListTail(l2);
			}
		return 1;
		}
	if(is_float(t1) && is_float(t2) && FloatValue(t1)==FloatValue(t2))
		return 1;
	return 0;
	}
			

		
static int read_prm(int no, Atom dir)
{
	char cbuf[10250];
	FILE *f;
	int curl;
	
	if(dir)
		sprintf(cbuf,"%s/vars%d.mdl",AtomValue(dir),no);
	else
		sprintf(cbuf,"vars%d.mdl",no);
	
	f=fopen(cbuf,"r");
	if(f==NULL)
	{
		ErrorInfo(0);
		printf("ReadChep: can not open '%s'.\n",cbuf);
		return 0;
	}
	
	fgets(cbuf,10250,f);
	fgets(cbuf,10250,f);
	fgets(cbuf,10250,f);
	curl=0;
	
	while(fgets(cbuf,10250,f))
	{
		int pos1, pos2;
		Term p, v;
		List l;
		
		curl++;
		
		if(cbuf[0]=='=')
			continue;
		
		pos1=0;
		while(cbuf[pos1]==' ')
			pos1++;
		if(cbuf[pos1]=='%')
			continue;
		pos2=pos1;
		while(cbuf[pos2]!=' ' && cbuf[pos2]!='|')
			pos2++;
		p=NewAtom(cbuf+pos1, pos2-pos1);
		pos1=pos2;
		while(cbuf[pos1]!='|')
			pos1++;
		pos1++;
		pos2=pos1;
		while(cbuf[pos2]!='|')
			pos2++;
		cbuf[pos2]='.';
		v=ReadTermM(cbuf+pos1);
		if(v==1 || v==0 || (!is_integer(v) && !is_float(v)))
		{
			printf("ReadChep: error in file vars%d.mdl, expr for parameter %s\n",
					no, AtomValue(p));
			return 0;
		}
		
		 if(!is_parameter(p))
		 {
			 Term t;
			 t=MakeCompound1(OPR_PARAMETER,MakeCompound2(OPR_EQSIGN,p,v));
			 CallFunction(t,0);
		 }
		 else
		 {
			for(l=all_param_list();l;l=ListTail(l))
				 if(CompoundName(ListFirst(l))==p)
					 break;
			 if(l==0)
			 {
				 puts("Internal error (rcrprm01)");
				 return 1;
			 }
			 if(!EqualTermsF(CompoundArg1(ListFirst(l)),v))
			 {
				 printf("ReadChep: change of parameter %s: ",AtomValue(p));
				 WriteTerm(CompoundArg1(ListFirst(l)));
				 printf(" -> ");WriteTerm(v);puts("");
				 FreeAtomic(v);
			 } 
		 }
	 }
	 
	 printf("ReadChep: file vars%d.mdl, %d parameters read.\n",no,curl);
	 fclose(f);
	 
	 if(dir)
		sprintf(cbuf,"%s/func%d.mdl",AtomValue(dir),no);
	else
		sprintf(cbuf,"func%d.mdl",no);
	
	f=fopen(cbuf,"r");
	if(f==NULL)
	{
		ErrorInfo(0);
		printf("ReadChep: can not open '%s'.\n",cbuf);
		return 0;
	}
	
	fgets(cbuf,10250,f);
	fgets(cbuf,10250,f);
	fgets(cbuf,10250,f);
	curl=0;
	
	while(fgets(cbuf,10250,f))
	{
		int pos1, pos2;
		Term p, v;
		List l;
		
		curl++;
		
		if(cbuf[0]=='=')
			continue;
		
		rpl_caret(cbuf);
		
		pos1=0;
		while(cbuf[pos1]==' ')
			pos1++;
		if(cbuf[pos1]=='%')
			continue;
		pos2=pos1;
		while(cbuf[pos2]!=' ' && cbuf[pos2]!='|')
			pos2++;
		p=NewAtom(cbuf+pos1, pos2-pos1);
		pos1=pos2;
		while(cbuf[pos1]!='|')
			pos1++;
		pos1++;
		pos2=pos1;
		while(cbuf[pos2] && cbuf[pos2]!='|' && cbuf[pos2]!='%'
				&& cbuf[pos2]!='"')
			pos2++;
		if(cbuf[pos2]=='"')
		{
			pos2++;
			while(cbuf[pos2] && cbuf[pos2]!='"')
			{
				if(cbuf[pos2]=='%')
					cbuf[pos2]='#';
				pos2++;
			}
			pos2++;
			while(cbuf[pos2] && cbuf[pos2]!='|' && cbuf[pos2]!='%'
				&& cbuf[pos2]!='"') pos2++;
		}
		cbuf[pos2]='.';
		cbuf[pos2+1]=0;
		
		v=ReadTermM(cbuf+pos1);
		
		if(v==1 || v==0 )
		{
			printf("ReadChep: error in file func%d.mdl, expr for parameter %s\n",
					no, AtomValue(p));
			return 0;
		}
		
		
		 if(!is_parameter(p))
		 {
			 Term t;
			 t=MakeCompound1(OPR_PARAMETER,MakeCompound2(OPR_EQSIGN,p,v));
			 CallFunction(t,0);
		 }
		 else
		 {
			 for(l=all_param_list();l;l=ListTail(l))
				 if(CompoundName(ListFirst(l))==p)
					 break;
			 if(l==0)
			 {
				 puts("Internal error (rcrprm01)");
				 return 1;
			 }
			 if(!EqualTermsF(CompoundArg1(ListFirst(l)),v))
			 {
				 printf("ReadChep: change of parameter %s: ",AtomValue(p));
				 WriteTerm(CompoundArg1(ListFirst(l)));
				 printf(" -> ");WriteTerm(v);puts("");
				 FreeAtomic(v);
			 }
		 }
	 }
	 
	 printf("ReadChep: file func%d.mdl, %d parameters read.\n",no,curl);
	 fclose(f);
	 
	 return 1;
	 
	}
	
static Atom read_atom(char *cbuf, int *pos)
{
	int i,j,k;
	i=k=*pos;
	while(cbuf[i]==' ') i++;
	k=i;
	if(cbuf[i]=='%')
		return 0;
	while(cbuf[i]!='|')
		i++;
	i++;
	j=i-2;
	while(cbuf[j]==' ')
		j--;
	*pos=i;
	if(j<k)
		return 0;
	return NewAtom(cbuf+k,j-k+1);
}
	
static int read_prt(int no, Atom dir)
{
	char cbuf[1250];
	FILE *f;
	int curl=0, curcol=0, ca=0;
	
	if(dir)
		sprintf(cbuf,"%s/prtcls%d.mdl",AtomValue(dir),no);
	else
		sprintf(cbuf,"prtcls%d.mdl",no);
	
	f=fopen(cbuf,"r");
	if(f==NULL)
	{
		ErrorInfo(0);
		printf("ReadChep: can not open '%s'.\n",cbuf);
		return 0;
	}
	
	fgets(cbuf,1250,f);
	fgets(cbuf,1250,f);
	fgets(cbuf,1250,f);
	do
	{
		if(cbuf[curl]=='|')
			curcol++;
		curl++;
	} while(curcol!=3);
	
	while(cbuf[curl]==' ') curl++;
	if(cbuf[curl]=='n' || cbuf[curl]=='N'||cbuf[curl]=='p' || cbuf[curl]=='P')
		ca=1;
	curl=0;
	
	while(fgets(cbuf,1250,f))
	{
		Atom name, pname, apname, pdg, spa, mass, width, col, aux, lt, alt;
		int pos1;
		Term prt;
		
		if(cbuf[0]=='=')
			continue;
		
		
		pos1=0;
		name=read_atom(cbuf,&pos1);
		pname=read_atom(cbuf,&pos1);
		apname=read_atom(cbuf,&pos1);
		if(ca)
			pdg=read_atom(cbuf,&pos1);
		spa=read_atom(cbuf,&pos1);
		mass=read_atom(cbuf,&pos1);
		width=read_atom(cbuf,&pos1);
		col=read_atom(cbuf,&pos1);
		aux=read_atom(cbuf,&pos1);
		/*lt=read_atom(cbuf,&pos1);
		alt=read_atom(cbuf,&pos1);*/
		
		prt=MakeCompound(OPR_PARTICLE,8);
		SetCompoundArg(prt,1,pname);
		SetCompoundArg(prt,2,apname);
		SetCompoundArg(prt,3,name);
		SetCompoundArg(prt,4,NewInteger(AtomValue(spa)[0]-'0'));
		if(mass && AtomValue(mass)[0]!='0')
			SetCompoundArg(prt,5,mass);
		if(width && AtomValue(width)[0]!='0')
			SetCompoundArg(prt,6,width);
		
		if(aux)
			switch(AtomValue(aux)[0])
			{
				case 'G':
						SetCompoundArg(prt,7,A_GAUGE); break;
				case 'L':
						SetCompoundArg(prt,7,A_LEFT); break;
				case 'R':
						SetCompoundArg(prt,7,A_RIGHT); break;
				default:
						SetCompoundArg(prt,7,aux);
			}
		
		if(ca)
		{
			int val;
			sscanf(AtomValue(pdg),"%d",&val);
			SetAtomProperty(pname,NewAtom("_pdg",0),NewInteger(val));
		}
			
		/*if(lt)
			SetAtomProperty(pname,A_TEXNAME,lt);
		if(alt)
			SetAtomProperty(apname,A_TEXNAME,alt);
		*/
		if(col && AtomValue(col)[0]=='3')
			SetCompoundArg(prt,8,MakeList1(MakeCompound2(A_COLOR,
					NewAtom("c3",0),NewAtom("c3b",0))));
		
		if(col && AtomValue(col)[0]=='8')
			SetCompoundArg(prt,8,MakeList1(MakeCompound2(A_COLOR,
					NewAtom("c8",0),NewAtom("c8",0))));
		
		if(!is_particle(pname,NULL))
			AddIMParticle(prt);
		else
		{
			Term oprt;
			oprt=GetAtomProperty(pname,PROP_TYPE);
			
			if(CompoundName(oprt)!=OPR_PARTICLE)
			{
				printf("ReadChep: particle %s is previously ",AtomValue(pname));
				WriteTerm(oprt);
				puts(".");
				return 0;
			}
			if(CompoundArg2(oprt)!=apname)
			{
				printf("ReadChep: particle %s : antiparticle changed %s->%s\n",
						AtomValue(pname),AtomValue(CompoundArg2(oprt)),AtomValue(apname));
				return 0;
			}
			if(CompoundArgN(oprt,3)!=name)
			{
				printf("ReadChep: particle %s : full name changed '%s'->'%s'\n",
						AtomValue(pname),AtomValue(CompoundArgN(oprt,3)),AtomValue(name));
			}
			if(CompoundArgN(oprt,4)!=CompoundArgN(prt,4))
			{
				printf("ReadChep: particle %s : spin changed %ld->%s\n",
						AtomValue(pname),IntegerValue(CompoundArgN(oprt,3)),AtomValue(spa));
				return 0;
			}
			if(CompoundArgN(oprt,5)!=CompoundArgN(oprt,5))
			{
				printf("ReadChep: particle %s : mass changed %s->%s\n",
						AtomValue(pname),AtomValue(CompoundArgN(oprt,5)),AtomValue(mass));
			}
			
			if(CompoundArgN(oprt,6)!=CompoundArgN(oprt,6))
			{
				printf("ReadChep: particle %s : width changed %s->%s\n",
						AtomValue(pname),AtomValue(CompoundArgN(oprt,6)),AtomValue(width));
			}
			
			if(CompoundArgN(oprt,7)!=CompoundArgN(oprt,7))
			{
				printf("ReadChep: particle %s : aux field changed %s->%s\n",
						AtomValue(pname),AtomValue(CompoundArgN(oprt,7)),AtomValue(CompoundArgN(prt,7)));
			}
			if(!EqualTerms(CompoundArgN(oprt,8),CompoundArgN(prt,8)))
			{
				printf("ReadChep: particle %s : color changed ",AtomValue(pname));
				WriteTerm(CompoundArgN(oprt,8));printf("->");
				WriteTerm(CompoundArgN(prt,8));puts("");
			}
	
		}
			
/*		printf("%10s|%4s|%4s|%3s|%5s|%5s|%3s|%2s|%10s|%10s\n",AtomValue(name),AtomValue(pname),
				AtomValue(apname),
				AtomValue(spa), AtomValue(mass), AtomValue(width),
				AtomValue(col), 
				AtomValue(aux), AtomValue(lt), AtomValue(alt));
*/			
		

	}

		
	return 1;
	
}

static List l_mlt(List t1, List t2)
{
	List l1,l2,ret;
	
	ret=NewList();
	for(l1=t1;l1;l1=ListTail(l1))
	for(l2=t2;l2;l2=ListTail(l2))
	{
		int i1,i2;
		List q1,q2;
		i1=IntegerValue(CompoundArg1(ListFirst(l1)));
		i2=IntegerValue(CompoundArg1(ListFirst(l2)));
		q1=CopyTerm(CompoundArg2(ListFirst(l1)));
		q2=CopyTerm(CompoundArg2(ListFirst(l2)));
		ret=AppendLast(ret,MakeCompound2(OPR_MLT,
				NewInteger(i1*i2),ConcatList(q1,q2)));
	}
	
	FreeAtomic(t1);
	FreeAtomic(t2);
	return ret;
}

extern int opSetGpm;

static List e_2_l(Term t)
{
	Term t1, t2;
	
	if(is_integer(t))
		return MakeList1(MakeCompound2(OPR_MLT,t,0));
	
	if(is_atom(t) && strcmp(AtomValue(t),"G5")==0)
		return MakeList1(MakeCompound2(OPR_MLT,NewInteger(1),
				MakeList1(A_GAMMA5)));
	
	if(is_atom(t))
		return MakeList1(MakeCompound2(OPR_MLT,NewInteger(1),
				MakeList1(t)));
	
	if(!is_compound(t))
	{
		printf("ReadChep: e2l: bad expression '");
		WriteTerm(t);
		puts("'.");
		return NewList();
	}

	if(opSetGpm && CompoundArity(t)==2 && CompoundArg1(t)==NewInteger(1)
			&& is_atom(CompoundArg2(t))
			&& ( strcmp(AtomValue(CompoundArg2(t)),"G5")==0 )
			&& (CompoundName(t)==OPR_PLUS || CompoundName(t)==OPR_MINUS))
	{
		if(CompoundName(t)==OPR_PLUS)
			return MakeList1(MakeCompound2(OPR_MLT,NewInteger(2),
				MakeList1(A_GAMMAP)));
		else
			return MakeList1(MakeCompound2(OPR_MLT,NewInteger(2),
				MakeList1(A_GAMMAM)));
	}
	
	
	if(strcmp(AtomValue(CompoundName(t)),"G")==0 && CompoundArity(t)==1 &&
			is_atom(CompoundArg1(t)))
	{
		Atom v;
		v=ConsumeCompoundArg(t,1);
		FreeAtomic(t);
		return MakeList1(MakeCompound2(OPR_MLT,NewInteger(1),
				MakeList1(MakeCompound1(A_GAMMA,v))));
	}
	
	if(CompoundName(t)==OPR_PLUS && CompoundArity(t)==1)
	{
		t1=ConsumeCompoundArg(t,1);
		FreeAtomic(t);
		return e_2_l(t1);
	}
	
	if(CompoundName(t)==OPR_PLUS && CompoundArity(t)==2)
	{
		t1=ConsumeCompoundArg(t,1);
		t2=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		return ConcatList(e_2_l(t1),e_2_l(t2));
	}
	
	if(CompoundName(t)==OPR_MINUS && CompoundArity(t)==1)
	{
		List l;
		t1=ConsumeCompoundArg(t,1);
		FreeAtomic(t);
		t1 = e_2_l(t1);
		for(l=t1;l;l=ListTail(l))
		{
			int i;
			i=IntegerValue(CompoundArg1(ListFirst(l)));
			SetCompoundArg(ListFirst(l),1,NewInteger(-i));
		}
		return t1;
	}
	
	if(CompoundName(t)==OPR_MINUS && CompoundArity(t)==2)
	{
		List l;
		t1=ConsumeCompoundArg(t,1);
		t2=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		t1=e_2_l(t1);
		t2=e_2_l(t2);
		for(l=t2;l;l=ListTail(l))
		{
			int i;
			i=IntegerValue(CompoundArg1(ListFirst(l)));
			SetCompoundArg(ListFirst(l),1,NewInteger(-i));
		}
		return ConcatList(t1,t2);
	}
	
	if(CompoundName(t)==OPR_MLT && CompoundArity(t)==2)
	{
		t1=ConsumeCompoundArg(t,1);
		t2=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		t1=e_2_l(t1);
		t2=e_2_l(t2);
		return l_mlt(t1,t2);
	}
	
	if(CompoundName(t)==OPR_POW && CompoundArity(t)==2)
	{
		int po;
		t1=ConsumeCompoundArg(t,1);
		t2=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		
		if(!is_integer(t2) || IntegerValue(t2)<1)
		{
			printf("ReadChep: e2l: illegal power '");
			WriteTerm(t2);
			printf("'.\n");
			return 0;
		}
	
		po=IntegerValue(t2);
		t1=e_2_l(t1);
		t2=CopyTerm(t1);

		while(po>1)
		{
			t2=l_mlt(t2,CopyTerm(t1));
			po--;
		}
		FreeAtomic(t1);
		return t2;
	}
		
	
	if(CompoundName(t)==OPR_CARET && CompoundArity(t)==2)
	{
		List l1,l2,ret;
		t1=ConsumeCompoundArg(t,1);
		t2=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		t1=e_2_l(t1);
		t2=e_2_l(t2);
		ret=NewList();
		for(l1=t1;l1;l1=ListTail(l1))
		for(l2=t2;l2;l2=ListTail(l2))
		{
			int i1,i2;
			List q1,q2;
			i1=IntegerValue(CompoundArg1(ListFirst(l1)));
			i2=IntegerValue(CompoundArg1(ListFirst(l1)));
			q1=CompoundArg2(ListFirst(l1));
			q2=CompoundArg2(ListFirst(l2));
			if(ListLength(q1)!=1 || ListLength(q2)!=1)
			{
				puts("Internal error (rche2l2)");
				return 0;
			}
			q1=ListFirst(q1);
			q2=ListFirst(q2);
			if(!is_atom(q1) || !is_atom(q2))
			{
				puts("Internal error (rche2l2a)");
				return 0;
			}
			ret=AppendLast(ret,MakeCompound2(OPR_MLT,
					NewInteger(i1*i2),MakeList1(MakeCompound2(A_POINT,q1,q2))));
		}
		FreeAtomic(t1);
		FreeAtomic(t2);
		return ret;
	}
	
	printf("ReadChep: e2l: can not parse expr '");
	WriteTerm(t);
	puts("'.");
	
	return NewList();
	
}

static Term c_mlt(Term m1, Term m2)
{
	Term ret;
	
	ret=MakeCompound(A_MTERM,4);
	SetCompoundArg(ret,1,NewInteger(
			IntegerValue(CompoundArg1(m1))*IntegerValue(CompoundArg1(m2))));
	SetCompoundArg(ret,2,NewInteger(
			IntegerValue(CompoundArg2(m1))*IntegerValue(CompoundArg2(m2))));
	
	SetCompoundArg(ret,3,ConcatList(
			ConsumeCompoundArg(m1,3),ConsumeCompoundArg(m2,3)));
	SetCompoundArg(ret,4,ConcatList(
			ConsumeCompoundArg(m1,4),ConsumeCompoundArg(m2,4)));
	
	FreeAtomic(m1);
	FreeAtomic(m2);
	return ret;
}
	

static Term c_2_l1(Term t)
{
	Term t1, t2;
	Term ret;
		
	if(t==0)
		return 0;
	
	if(is_integer(t))
	{
		ret=MakeCompound(A_MTERM,4);
		SetCompoundArg(ret,2,NewInteger(1));
		SetCompoundArg(ret,1,t);
		return ret;
	}
		
	if(is_atom(t))
	{
		ret=MakeCompound(A_MTERM,4);
		SetCompoundArg(ret,1,NewInteger(1));
		SetCompoundArg(ret,2,NewInteger(1));
		SetCompoundArg(ret,3,MakeList1(t));
		return ret;
	}
	
	if(!is_compound(t))
	{
		printf("ReadChep: c2l: bad expression '");
		WriteTerm(t);
		puts("'.");
		return 0;
	}
	
	if(CompoundName(t)==OPR_PLUS && CompoundArity(t)==1)
	{
		t1=ConsumeCompoundArg(t,1);
		FreeAtomic(t);
		return c_2_l1(t1);
	}
	
	if(CompoundName(t)==OPR_MINUS && CompoundArity(t)==1)
	{
		t1=ConsumeCompoundArg(t,1);
		FreeAtomic(t);
		t1 = c_2_l1(t1);
		SetCompoundArg(t1,1,NewInteger(-IntegerValue(CompoundArg1(t1))));
		return t1;
	}
	
	
	if(CompoundName(t)==OPR_MLT && CompoundArity(t)==2)
	{
		t1=ConsumeCompoundArg(t,1);
		t2=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		t1=c_2_l1(t1);
		t2=c_2_l1(t2);
		return c_mlt(t1,t2);
	}
	
	if(CompoundName(t)==OPR_POW && CompoundArity(t)==2)
	{
		int po;
		t1=ConsumeCompoundArg(t,1);
		t2=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		
		if(!is_integer(t2) || IntegerValue(t2)<1)
		{
			printf("ReadChep: c2l: illegal power '");
			WriteTerm(t2);
			printf("'.\n");
			return 0;
		}
	
		po=IntegerValue(t2);
		t1=c_2_l1(t1);
		t2=CopyTerm(t1);

		while(po>1)
		{
			t2=c_mlt(t2,CopyTerm(t1));
			po--;
		}
		FreeAtomic(t1);
		return t2;
	}
		
	if(CompoundName(t)==OPR_DIV && CompoundArity(t)==2)
	{
		Term q1,q2;

		t1=ConsumeCompoundArg(t,1);
		t2=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		t1=c_2_l1(t1);
		t2=c_2_l1(t2);
		q1=CompoundArg1(t2);
		q2=CompoundArg2(t2);
		SetCompoundArg(t2,1,q2);
		SetCompoundArg(t2,2,q1);
		q1=ConsumeCompoundArg(t2,3);
		q2=ConsumeCompoundArg(t2,4);
		SetCompoundArg(t2,3,q2);
		SetCompoundArg(t2,4,q1);
		return c_mlt(t1,t2);
	}
	
	
	printf("ReadChep: e2l: can not parse expr '");
	WriteTerm(t);
	puts("'.");
	
	return NewList();
	
	
}

static Term c_2_l(Term cf, Atom *prt)
{
	Atom sprt[4];
	int i,j, n,d,c;
	List l1,l2;
	int pno, sf=1;

	sprt[0]=prt[1];
	sprt[1]=prt[2];
	sprt[2]=prt[3];
	sprt[3]=prt[4];
	if(prt[4]==0)
		pno=3;
	else
		pno=4;
	for(i=0;i<5;i++)
	for(j=0;j<pno-1;j++)
		if(sprt[j]>sprt[j+1])
		{
			Atom tmp;
			tmp=sprt[j];
			sprt[j]=sprt[j+1];
			sprt[j+1]=tmp;
		}
	
	if(pno==3)
	{
		if(sprt[0]==sprt[1] && sprt[1]==sprt[2])
			sf=6;
		else if(sprt[0]==sprt[1] || sprt[1]==sprt[2])
			sf=2;
	}
	
	if(pno==4)
	{
		if(sprt[0]==sprt[1] && sprt[1]==sprt[2] && sprt[2]==sprt[3])
			sf=24;
		else if((sprt[0]==sprt[1] && sprt[1]==sprt[2]) ||
			(sprt[1]==sprt[2] && sprt[2]==sprt[3]))
			sf=6;
		else if(sprt[0]==sprt[1] && sprt[2]==sprt[3])
			sf=4;
		else if(sprt[0]==sprt[1] || sprt[1]==sprt[2] || sprt[2]==sprt[3])
			sf=2;
	}

	cf=c_2_l1(cf);

	
	if(cf==0)
		return 0;
	n=IntegerValue(CompoundArg1(cf));
	d=IntegerValue(CompoundArg2(cf));
	if(d<0)
	{
		d=-d;
		n=-n;
	}
	d*=sf;
	c=gcf(n,d);
	n/=c;
	d/=c;
	
	l1=ConsumeCompoundArg(cf,3);
	for(l2=l1;l2;l2=ListTail(l2))
		if(ListFirst(l2)==A_I)
		{
			l1=CutFromList(l1,l2);
			d=-d;
			break;
		}
	for(l2=l1;l2;l2=ListTail(l2))
	{
		Term t;
		t=MakeCompound2(OPR_PARAMETER,0,ListFirst(l2));
		ChangeList(l2,t);
	}
	SetCompoundArg(cf,3,l1);
	SetCompoundArg(cf,1,NewInteger(n));
	SetCompoundArg(cf,2,NewInteger(d));
	
	l1=ConsumeCompoundArg(cf,4);
	for(l2=l1;l2;l2=ListTail(l2))
	{
		Term t;
		t=MakeCompound2(OPR_PARAMETER,0,ListFirst(l2));
		ChangeList(l2,t);
	}
	SetCompoundArg(cf,4,l1);
	
	return cf;
}

static Term mk_a1(Term cf, List lp, Term colobj, Label fei1, Label fei2,
		List prtl)
{
	List l1, l2, l3;
	
	List la1=0;
	
	for(l1=lp;l1;l1=ListTail(l1))
	{
		Term cf1;
		Term lp1;
		Term lf=0;
		Label fei11=fei1;
		Label pind[5];
		int pcn=0;
		
		pind[1]=pind[2]=pind[3]=pind[4]=0;
		cf1=CopyTerm(cf);
		l2=ConsumeCompoundArg(cf1,3);
		lp1=ListFirst(l1);
		
		SetCompoundArg(cf1,1,NewInteger(
				IntegerValue(CompoundArg1(cf1))*IntegerValue(CompoundArg1(lp1))));
		
		for(l3=CompoundArg2(lp1);l3;l3=ListTail(l3))
		{
			Term cuo;
			cuo=ListFirst(l3);
			if(cuo==A_GAMMA5 || cuo==A_GAMMAP || cuo==A_GAMMAM)
			{
				Term t1,t2;
				List li;
				t1=MakeCompound2(A_I,
					MakeCompound2(A_LORENTZ,NewInteger(1),NewInteger(0)),fei11);
				t2=MakeCompound2(A_I,
					MakeCompound2(A_LORENTZ,NewInteger(0),NewInteger(1)),0);
				if(lf)
					SetCompoundArg(lf,2,fei11);
				lf=t2;
				fei11=NewLabel();
				li=MakeList2(t1,t2);
				t1=MakeCompound2(OPR_SPECIAL,li,cuo);
				l2=AppendLast(l2,t1);
				continue;
			}
			if(is_compound(cuo) && CompoundName(cuo)==A_GAMMA)
			{
				Term t1,t2,t3;
				List li;
				Atom ind, lind;
				char *cind;
				ind=CompoundArg1(cuo);
				cind=AtomValue(ind);
				if(cind[0]=='m')
					lind=ind;
				else
				{
					lind=NewLabel();
					pind[cind[1]-'0']=lind;
				}
				t1=MakeCompound2(A_I,
					MakeCompound2(A_LORENTZ,NewInteger(1),NewInteger(0)),fei11);
				t2=MakeCompound2(A_I,
					MakeCompound2(A_LORENTZ,NewInteger(0),NewInteger(1)),0);
				t3=MakeCompound2(A_I,
					MakeCompound2(A_LORENTZ,NewInteger(2),NewInteger(2)),lind);
				if(lf)
					SetCompoundArg(lf,2,fei11);
				lf=t2;
				fei11=NewLabel();
				li=MakeList3(t1,t2,t3);
				t1=MakeCompound2(OPR_SPECIAL,li,A_GAMMA);
				l2=AppendLast(l2,t1);
				continue;
			}
			
			if(is_atom(cuo))
			{
				Term t1;
				t1=MakeCompound2(OPR_PARAMETER,0,cuo);
				l2=AppendLast(l2,t1);
				continue;
			}
			
			if(!is_compound(cuo) || CompoundArity(cuo)!=2)
			{
				puts("Internal error (rch05)");
				return 0;
			}
			
			{
				Atom i1, i2;
				char *ci1, *ci2;
				Label li1, li2;
				Term t1, t2;
				
				i1=CompoundArg1(cuo);
				i2=CompoundArg2(cuo);
				ci1=AtomValue(i1);
				ci2=AtomValue(i2);
				
				if(ci1[0]=='p')
				{
					li1=NewLabel();
					pind[ci1[1]-'0']=li1;
				}
				else
					li1=i1;
				
				if(ci2[0]=='p')
				{
					li2=NewLabel();
					pind[ci2[1]-'0']=li2;
				}
				else
					li2=i2;
				
				t1=MakeCompound2(A_I, MakeCompound2(
						A_LORENTZ,NewInteger(2),NewInteger(2)),li1);
				t2=MakeCompound2(A_I, MakeCompound2(
						A_LORENTZ,NewInteger(2),NewInteger(2)),li2);
				
				l2=AppendLast(l2,MakeCompound2(
						OPR_SPECIAL,MakeList2(t1,t2),A_DELTA));
			}
			
		}
		
		if(fei11)
		{
			if(lf)
				SetCompoundArg(lf,2,fei2);
			else
			{
				Term t1, t2;
				t1=MakeCompound2(A_I, MakeCompound2(
						A_LORENTZ,NewInteger(1),NewInteger(0)),fei11);
				t2=MakeCompound2(A_I, MakeCompound2(
						A_LORENTZ,NewInteger(0),NewInteger(1)),fei2);
				l2=AppendLast(l2,MakeCompound2(
						OPR_SPECIAL,MakeList2(t1,t2),A_DELTA));
			}
		}
		
		if(colobj)
			l2=AppendLast(l2,CopyTerm(colobj));
		
		for(l3=prtl;l3;l3=ListTail(l3))
		{
			Term ps;
			pcn++;
			ps=CopyTerm(ListFirst(l3));
			if(pind[pcn])
			{
				Term t1;
				t1=MakeCompound2(A_I, MakeCompound2(
						A_LORENTZ,NewInteger(2),NewInteger(2)),pind[pcn]);
				l2=AppendLast(l2,MakeCompound2(
						OPR_SPECIAL,MakeList1(t1),A_MOMENT));
			}
			l2=AppendLast(l2,ps);
		}
		
		SetCompoundArg(cf1,3,l2);
		la1=AppendLast(la1,cf1);
	}
	
	
	FreeAtomic(lp);
	FreeAtomic(cf);
	if(colobj)
		FreeAtomic(colobj);
	FreeAtomic(prtl);
	
	return MakeCompound2(A_ALG1,la1,0);
	
}

extern int LagrHashSize;
void alg2_hash_add(List *, int , List);
extern List *lagr_hash;	

static int subtr_lagr=0;

static int read_vrt(int no, Atom dir)
{
	char cbuf[12500];
	FILE *f;
	int curl;
	
	if(dir)
		sprintf(cbuf,"%s/lgrng%d.mdl",AtomValue(dir),no);
	else
		sprintf(cbuf,"lgrng%d.mdl",no);
	
	f=fopen(cbuf,"r");
	if(f==NULL)
	{
		ErrorInfo(0);
		printf("ReadChep: can not open '%s'.\n",cbuf);
		return 0;
	}
	
	SetOperator(OPR_PLUS,OP_FX,499);
	
	fgets(cbuf,12500,f);
	fgets(cbuf,12500,f);
	fgets(cbuf,12500,f);
	curl=0;
	
	while(fgets(cbuf,12500,f))
	{
		int pos1;
		int pos2;
		int prtno;
		int fcnt;
		
		Atom  prt[5];
		Label fei[5];
		Label coi[5];
		Term  cos[5];
		int   cot[5];
		List prtl;
		Label fei1=0, fei2=0;
		
		Term colobj=1;
		
		Term cf, lp;
		
		cot[1]=cot[2]=cot[3]=cot[4]=0;
		
		if(cbuf[0]=='=')
			continue;
		
		rpl_caret(cbuf);
		
		pos1=0;
		prt[1]=read_atom(cbuf,&pos1);
		prt[2]=read_atom(cbuf,&pos1);
		prt[3]=read_atom(cbuf,&pos1);
		prt[4]=read_atom(cbuf,&pos1);
		
		if(prt[1]==0)
			continue;
				
		
		pos2=pos1;
		while(cbuf[pos2]!='|')
		{
			if(cbuf[pos2]=='%') cbuf[pos2]='.';
			pos2++;
		}
		cbuf[pos2]='.';
		
		cf=ReadTermM(cbuf+pos1);
		if(cf==0 || cf==1)
			return 0;
		
		pos1=pos2+1;
		pos2=pos1;
		while(cbuf[pos2])
			pos2++;
		pos2--;
		while(cbuf[pos2]==' ' || cbuf[pos2]=='\n')
			cbuf[pos2--]=0;
		
		cbuf[pos2+1]='.';
		cbuf[pos2+2]=0;
		while(pos2>=pos1)
		{
			if(cbuf[pos2]=='.')
				cbuf[pos2]='^';
			pos2--;
		}
		
		lp=ReadTermM(cbuf+pos1);
		if(lp==0 || lp==1)
			return 0;
		
		if(prt[4]==0)
			prtno=3;
		else
			prtno=4;

		cf=c_2_l(cf,prt);
		if(subtr_lagr)
			SetCompoundArg(cf,1,NewInteger(-IntegerValue(CompoundArg1(cf))));

		
  /****************** Particle list *********************/
				
		prtl=NewList();
		fcnt=0;
		for(pos1=1;pos1<=prtno;pos1++)
		{
			Term ind;
			int ll;
			coi[pos1]=cot[pos1]=cos[pos1]=0;
			ind=CopyTerm(GetAtomProperty(prt[pos1],PROP_INDEX));
			ll=ListLength(ind);
			for(pos2=1;pos2<=ll;pos2++)
			{
				Term it;
				Label la=0;
				it=CompoundArg1(ListNth(ind,pos2));
				if(CompoundName(it)==A_LORENTZ && CompoundArg1(it)==NewInteger(2))
				{
					char cbuf[4];
					if(pos2==1)
						sprintf(cbuf,"m%d",pos1);
					else
						sprintf(cbuf,"M%d",pos1);
					la=NewAtom(cbuf,0);
				}
				if(CompoundName(it)==A_LORENTZ && CompoundArg1(it)==NewInteger(1))
				{
					la=NewLabel();
					fei[pos1]=la;
					if(fcnt==0)
						fei1=la;
					else
						fei2=la;
					if(fcnt==0)
					{
						char cbuf[8];
						strcpy(cbuf,AtomValue(prt[pos1]));
						strcat(cbuf,".c");
						prt[pos1]=NewAtom(cbuf,0);
					}
					fcnt++;
				}
				if(CompoundName(it)==A_LORENTZ && CompoundArg1(it)==NewInteger(0))
				{
					la=NewLabel();
					fei[pos1]=la;
					if(fcnt==0)
						fei1=la;
					else
						fei2=la;
					if(fcnt==1)
					{
						char cbuf[8];
						strcpy(cbuf,AtomValue(prt[pos1]));
						strcat(cbuf,".c");
						prt[pos1]=NewAtom(cbuf,0);
					}
					fcnt++;
				}
				if(CompoundName(it)==A_COLOR)
				{
					la=NewLabel();
					cos[pos1]=CopyTerm(it);
					coi[pos1]=la;
					cot[pos1]=IntegerValue(CompoundArg1(
							GetAtomProperty(prt[pos1],A_COLOR)));
				}
				if(la==0)
				{
					printf("internal error (rdvrt01) ");
					WriteTerm(it);
					puts("");
					return 0;
				}
				SetCompoundArg(ListNth(ind,pos2),2,la);
			}
			
			prtl=AppendLast(prtl,MakeCompound2(OPR_FIELD,ind,prt[pos1]));
				
		}
		
		/*********************************************************
		*               Determine color object                   *
		**********************************************************/
		
		{
			int ctypes[4];
			ctypes[1]=ctypes[2]=ctypes[3]=0;
			
			for(pos1=1;pos1<=prtno;pos1++)
			{
				ctypes[cot[pos1]]++;
			}
			
			if(ctypes[1]==0 && ctypes[2]==0 && ctypes[3]==0)
				colobj=0;
			else if(ctypes[1]==1 && ctypes[2]==1 && ctypes[3]==0)
			{
				Term c1=0, c2=0, i1=0, i2=0;
				for(pos1=1;pos1<=prtno;pos1++)
				{
					if(cot[pos1]==1)
						c1=cos[pos1],i2=coi[pos1];
					if(cot[pos1]==2)
						c2=cos[pos1],i1=coi[pos1];
				}
				colobj=MakeCompound2(OPR_SPECIAL,MakeList2(
						MakeCompound2(A_I,c1,i1),MakeCompound2(A_I,c2,i2)),
						A_DELTA);
				
			}
			else if(ctypes[1]==0 && ctypes[2]==0 && ctypes[3]==2)
			{
				Term c1=0, c2=0, i1=0, i2=0;
				for(pos1=1;pos1<=prtno;pos1++)
				{
					if(cot[pos1]==3)
					{
						if(c1==0)
						c1=cos[pos1],i2=coi[pos1];
							else
						c2=cos[pos1],i1=coi[pos1];
					}
				}
				colobj=MakeCompound2(OPR_SPECIAL,MakeList2(
						MakeCompound2(A_I,c1,i1),MakeCompound2(A_I,c2,i2)),
						A_DELTA);
				
			}
			else if(ctypes[1]==1 && ctypes[2]==1 && ctypes[3]==1)
			{
				Term c1=0, c2=0, c3=0, i1=0, i2=0, i3=0;
				for(pos1=1;pos1<=prtno;pos1++)
				{
					if(cot[pos1]==1)
						c1=cos[pos1],i2=coi[pos1];
					if(cot[pos1]==2)
						c2=cos[pos1],i1=coi[pos1];
					if(cot[pos1]==3)
						c3=cos[pos1],i3=coi[pos1];
				}
				colobj=MakeCompound2(OPR_SPECIAL,MakeList3(
						MakeCompound2(A_I,c1,i1),MakeCompound2(A_I,c2,i2),
						MakeCompound2(A_I,c3,i3)),
						ColorLambda);
			}
			else if(ctypes[1]==0 && ctypes[2]==0 && ctypes[3]==3)
			{
				Term c[3], i[3], no=0;
				for(pos1=1;pos1<=prtno;pos1++)
				{
					if(cot[pos1]==3)
					{
						c[no]=cos[pos1],i[no]=coi[pos1];
						no++;
					}
				}
				colobj=MakeCompound2(OPR_SPECIAL,MakeList3(
						MakeCompound2(A_I,c[0],i[0]),MakeCompound2(A_I,c[1],i[1]),
						MakeCompound2(A_I,c[2],i[2])),
						ColorF);
			}
			else
			{
				puts("internal error (rdvrt02)");
				colobj=0;
				WriteTerm(prtl);puts("");
				printf("%d %d %d %d\n",cot[1],cot[2],cot[3],cot[4]);
			}
		}

		
/*		WriteTerm(prtl);puts("");
		WriteTerm(colobj);puts("");
		
		WriteTerm(prt[1]);printf(" ");
		WriteTerm(prt[2]);printf(" ");
		WriteTerm(prt[3]);printf(" ");
		WriteTerm(prt[4]);printf(" ");
		WriteTerm(cf);printf(" ");
		WriteTerm(lp);printf("\n");
*/
		lp=e_2_l(lp);
/*		WriteTerm(cf);puts("");
		WriteTerm(lp);puts("\n");
*/
		
		lp=mk_a1(cf,lp,colobj,fei1,fei2,prtl);
		alg1_fix_delta(lp);
		
		if(lagr_hash==NULL)
		{
			int i;
			lagr_hash=malloc(sizeof(List)*LagrHashSize);
			if(lagr_hash==NULL)
			{
				puts("Error: can not allocate memory for Lagrangian"); 
				exit(0);
			}
			for(i=0;i<LagrHashSize;i++)
				lagr_hash[i]=NewList();
		}
		
		lp=Alg1to2(lp);
		alg2_hash_add(lagr_hash,LagrHashSize,lp);
		
	}
	
	return 1;
}
				
Term ProcReadChep(Term t, Term ind)
{
	Term t1;
	int modno;
	
	if(!is_compound(t) || CompoundArity(t)>2)
	{
		ErrorInfo(0);
		puts("ReadChep: need model number.");
		return 0;
	}
	
	
	if(strcmp(AtomValue(CompoundName(t)),"ReadChep")==0)
		subtr_lagr=0;
	else if(strcmp(AtomValue(CompoundName(t)),"ReadChepM")==0)
		subtr_lagr=1;
	else
	{
		ErrorInfo(0);
		printf("ReadChep: unknown mode\n");
		return 0;
	}
	
	t1=CompoundArg1(t);
	if(!is_integer(t1))
	{
		ErrorInfo(0);
		puts("ReadChep: incorrect model number.");
		return 0;
	}
	modno=IntegerValue(t1);
	
	t1=0;
	
	if(CompoundArity(t)==2)
	{
		t1=CompoundArg2(t);
		if(!is_atom(t1))
		{
			ErrorInfo(0);
			printf("ReadChep: incorrect directory name.\n");
			return 0;
		}
	}
	check_funcs=0;
	
	if(read_prm(modno,t1)==0)
	{
		ErrorInfo(0);
		puts("ReadChep: error reading vars&funcs");
		return 0;
	}

	if(read_prt(modno,t1)==0)
	{
		puts("ReadChep: error reading particles");
		return 0;
	}

	if(read_vrt(modno,t1)==0)
	{
		puts("ReadChep: error reading vertices");
		return 0;
	}
	check_funcs=1;
	
	return 0;
}
