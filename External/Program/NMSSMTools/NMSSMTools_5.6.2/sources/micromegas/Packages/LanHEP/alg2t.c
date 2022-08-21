#include "lanhep.h"


int TEX_lines = 25;
int TEX_spec_in_line = 30;
int TEX_set_dot=1;
int TEX_max_pno=10;

extern int MultByI;

static int cur_cnt;

static void fix_ino(List l3, int inomin)
{
	int inomax=inomin-1;
	List l1,l2;
	List ii=NewList();
	
	
	for(l1=l3;l1;l1=ListTail(l1))
	{
		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
		{
			int i;
			i=(int)IntegerValue(ListFirst(l2));
			if(i>=inomin && !ListMember(ii,NewInteger(i)))
			{
				ii=AppendLast(ii,NewInteger(i));
				if(i>inomax)
					inomax=i;
			}
		}
	}
	
/*	WriteTerm(l3); WriteTerm(ii); printf(" %d %d\n",inomin,inomax);*/
	
/*	if(ListLength(ii)!=inomax-inomin+1)
	{
		puts("Internal error (twff)");
	WriteTerm(l3); WriteTerm(ii); printf(" %d %d\n",inomin,inomax);
		return;
	}*/
	
	for(l1=l3;l1;l1=ListTail(l1))
	{
		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
		{
			Integer i;
			int j;
			i=ListFirst(l2);
			j=ListMember(ii,i);
			if(j)
			{
				ChangeList(l2,NewInteger(j+inomin-1));
			}
		}
	}
	
	FreeAtomic(ii);
	return;
}

static void fix_ctf(List prt, List ctf, List lp)
{
	List ctfi=0, ctfo=0;
	List l1,l2;
	int i,pino=0;
	
	if(is_empty_list(ctf))
		return;
	
	for(l1=prt;l1;l1=ListTail(l1))
		pino+=ListLength(CompoundArg1(CompoundArg2(ListFirst(l1))));
	
/*	WriteTerm(prt);printf(" : %d inds\n",pino);*/
	
	for(l1=ctf;l1;l1=ListTail(l1))
		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
			if(IntegerValue(ListFirst(l2))>pino && 
					!ListMember(ctfi,ListFirst(l2)))
				ctfi=AppendLast(ctfi,ListFirst(l2));
	
	
	for(i=pino+1,l1=ctfi;l1;l1=ListTail(l1),i++)
		ctfo=AppendLast(ctfo,NewInteger(i));
	
	if(EqualTerms(ctfi, ctfo))
	{
		FreeAtomic(ctfi);
		FreeAtomic(ctfo);
		return;
	}
	
	for(l1=ctfo;l1;l1=ListTail(l1))
	{
		if(ListMember(ctfi,ListFirst(l1)))
			continue;
		for(l2=ctfi;l2;l2=ListTail(l2))
			if(!ListMember(ctfo,ListFirst(l2)))
				break;
		if(is_empty_list(l2))
		{
			puts("Internal error (twl4fct)");
			return;
		}
		ctfi=AppendLast(ctfi,ListFirst(l1));
		ctfo=AppendLast(ctfo,ListFirst(l2));
	}
	
	
/*	WriteTerm(ctf);WriteTerm(ctfo);puts("");
	WriteTerm(ListFirst(lp));puts("");
*/	
	for(l1=ctf;l1;l1=ListTail(l1))
		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
			if((i=ListMember(ctfi,ListFirst(l2))))
				ChangeList(l2,ListNth(ctfo,i));
	
	for(;lp;lp=ListTail(lp))
		for(l1=CompoundArgN(ListFirst(lp),3);l1;l1=ListTail(l1))
		for(l2=CompoundArg2(ListFirst(l1));l2;l2=ListTail(l2))
			if((i=ListMember(ctfi,ListFirst(l2))))
				ChangeList(l2,ListNth(ctfo,i));	
	
	FreeAtomic(ctfi);
	FreeAtomic(ctfo);
}

static void wrt_prt_name(FILE *f, Atom prt)
	{
	Atom prop;
	prop=GetAtomProperty(prt, A_TEXNAME);
	if(prop)
		fprintf(f,"$%s{}_{",AtomValue(prop));
	else
		if(AtomValue(prt)[0]=='~')
			fprintf(f,"\\~{}${%s}_{",AtomValue(prt)+1);
		else
			fprintf(f,"${%s}_{",AtomValue(prt));
	}

static void wrt_cnf(FILE *f, Term dv)
	{
	long int num, den;
	num=IntegerValue(CompoundArg1(dv));
	den=IntegerValue(CompoundArg2(dv));

	if(num==1 && den==1)
		return;

	if(num==-1 && den==1)
		{
		fprintf(f,"-");
		return;
		}

	if(den==1)
		{
		fprintf(f,"%d",num);
		return;
		}
	if(num<0)
	 	{
	 	fprintf(f,"-");
	 	num*=-1;
	 	}
	 fprintf(f,"\\frac{%d}{%d}",num,den);
	 }


static void wrt_param(FILE *f, Atom p)
	{
	Atom prop;
	prop=GetAtomProperty(p,A_TEXNAME);
	if(prop)
		fprintf(f," %s",AtomValue(prop));
	else
		fprintf(f," %s",AtomValue(p));
	}


static void wrt_csf(FILE *f, List l)
	{
	List num,  den;
	num=den=NewList();

	if(is_empty_list(l))
		return;

	while(!is_empty_list(l))
		{
		Term t;
		t=ListFirst(l);
		if(IntegerValue(CompoundArg2(t))<0)
			den=AppendLast(den,t);
		else
			num=AppendLast(num,t);
		l=ListTail(l);
		}

	if(is_empty_list(den))
		{
		l=num;
		while(!is_empty_list(l))
			{
			if(CompoundArg1(ListFirst(l))==A_SQRT2)
				fprintf(f,"\\sqrt{2}");
			else
				{
				int po;
				po=(int)IntegerValue(CompoundArg2(ListFirst(l)));
				wrt_param(f,CompoundArg1(ListFirst(l)));
				if(po>1)
					fprintf(f,"{}^%d ",po);
				}
			l=ListTail(l);
			if(l && TEX_set_dot)
				fprintf(f," \\cdot");
			}
		RemoveList(num);
		return;
		}

	fprintf(f,"\\frac{");
	if(is_empty_list(num))
		fprintf(f,"1");
	else
		{
		l=num;
		while(!is_empty_list(l))
			{
			if(CompoundArg1(ListFirst(l))==A_SQRT2)
				fprintf(f,"\\sqrt{2}");
			else
				{
				int po;
				po=(int)IntegerValue(CompoundArg2(ListFirst(l)));
				wrt_param(f,CompoundArg1(ListFirst(l)));
				if(po>1)
					fprintf(f,"{}^%d ",po);
				}
			l=ListTail(l);
			if(l && TEX_set_dot)
				fprintf(f," \\cdot");
			}
		RemoveList(num);
		}
	fprintf(f,"}{");
	l=den;
	while(!is_empty_list(l))
		{
		if(CompoundArg1(ListFirst(l))==A_SQRT2)
			fprintf(f,"\\sqrt{2}");
		else
				{
				int po;
				po=-(int)IntegerValue(CompoundArg2(ListFirst(l)));
				wrt_param(f,CompoundArg1(ListFirst(l)));
				if(po>1)
					fprintf(f,"{}^%d ",po);
				}
		l=ListTail(l);
		if(l && TEX_set_dot)
			fprintf(f," \\cdot");
		}
	RemoveList(den);
	fprintf(f,"}");
	}


static char *vinarr[10] = {
	"\\mu", "\\nu",  "\\rho", "\\sigma", "\\lambda",
	"\\tau", "\\pi", "\\theta", "\\varepsilon", "\\eta" };

static char ibuf[4];

static int twl_1(FILE *f, Term a2)
	{
	List prt, ind, ind1, vind, vind1, l, l1, l2, m2l;
	int ino=1, vino=0, sino='a', oino='p', fflag, silcnt=0, lines=1;
	
	
	fix_ctf(CompoundArg1(a2),CompoundArgN(a2,4),CompoundArgN(a2,5));
	
	prt=CompoundArg1(a2);
	l=prt;
/*	WriteTerm(a2);
	puts("\n");*/
	ind=vind=NewList();
	while(!is_empty_list(l))
		{
		Term t;
		Atom n,n1;

		t=CompoundArg2(ListFirst(l));
		n=CompoundArg1(ListFirst(l));
		n1=GetAtomProperty(n,A_TEXNAME);
/*		if(n1)
			n=n1;*/
		l1=CompoundArg1(t);
		if(l1)
			{
			if(IntegerValue(ListFirst(l1))!=ino)
				{
				puts("Internal error (twl1).");
				return lines;
				}

			if(CompoundName(t)==OPR_VECTOR)
				{
				ind=AppendLast(ind, NewAtom(vinarr[vino++],0));
				vind=AppendLast(vind, NewInteger(ino));
				ino++;
				if(ListLength(l1)==2)
					{
					ibuf[0]=oino++;
					ind=AppendLast(ind,NewAtom(ibuf,1));
					ino++;
					}
				if(ListLength(l1)>2)
					{
					puts("Internal error (twl2v).");
					return lines;
					}
				goto cnt1;
				}

			if(CompoundName(t)==OPR_SPINOR)
				{
				ibuf[0]=sino++;
				ind=AppendLast(ind, NewAtom(ibuf,1));
				ino++;
				if(ListLength(l1)==2)
					{
					ibuf[0]=oino++;
					ind=AppendLast(ind,NewAtom(ibuf,1));
					ino++;
					}
				if(ListLength(l1)>2)
					{
					puts("Internal error (twl2sp).");
					return lines;
					}
				goto cnt1;
				}
				
			if(CompoundName(t)==OPR_SPINOR3)
				{
				ibuf[0]=sino++;
				ind=AppendLast(ind, NewAtom(ibuf,1));
				ino++;
				
				ind=AppendLast(ind, NewAtom(vinarr[vino++],0));
				vind=AppendLast(vind, NewInteger(ino));
				ino++;
				
				if(ListLength(l1)==3)
					{
					ibuf[0]=oino++;
					ind=AppendLast(ind,NewAtom(ibuf,1));
					ino++;
					}
				if(ListLength(l1)>3)
					{
					puts("Internal error (twl2sp).");
					return lines;
					}
				goto cnt1;
				}

			if(CompoundName(t)==OPR_TENSOR)
				{
				ind=AppendLast(ind, NewAtom(vinarr[vino++],0));
				vind=AppendLast(vind, NewInteger(ino));
				ino++;
				ind=AppendLast(ind, NewAtom(vinarr[vino++],0));
				vind=AppendLast(vind, NewInteger(ino));
				ino++;
				if(ListLength(l1)==3)
					{
					ibuf[0]=oino++;
					ind=AppendLast(ind,NewAtom(ibuf,1));
					ino++;
					}
				if(ListLength(l1)>3)
					{
					puts("Internal error (twl2t).");
					return lines;
					}
				goto cnt1;
				}

			if(CompoundName(t)==OPR_SCALAR)
				{
				if(ListLength(l1)==1)
					{
					ibuf[0]=oino++;
					ind=AppendLast(ind,NewAtom(ibuf,1));
					ino++;
					}
				if(ListLength(l1)>1)
					{
					puts("Internal error (twl2sc).");
					return lines;
					}
				goto cnt1;
				}
			puts("Internal error (twl2).");
			return lines;
			}
	cnt1:
		/*if(AtomValue(n)[0]=='~')
			fprintf(f,"\\~{}${%s}^{",AtomValue(n)+1);
		else
			fprintf(f,"${%s}_{",AtomValue(n));*/
		wrt_prt_name(f,n);
		while(!is_empty_list(l1))
			{
			fprintf(f,"%s ",AtomValue(ListNth(ind,(int)IntegerValue(ListFirst(l1)))));
			l1=ListTail(l1);
			}
		fprintf(f,"}$ \\phantom{-} ");
		l=ListTail(l);
		}

	fprintf(f," &\n\t");


	fprintf(f, "$");
	wrt_cnf(f, CompoundArg2(a2));
	wrt_csf(f, CompoundArgN(a2,3));
	silcnt=5;

	l2=CompoundArgN(a2,4);
		while(!is_empty_list(l2))
			{
			Term sp;
			sp=ListFirst(l2);
			if(CompoundName(sp)==A_MOMENT)
				{
				int ii;
				ii=(int)IntegerValue(ListFirst(CompoundArg2(sp)));
				if(ii>ino)
					{
					puts("Internal error twl4m.");
					goto aaa;
					}
				if(ii==ino)
					{
					ind=AppendLast(ind, NewAtom(vinarr[vino++],0));
					vind=AppendLast(vind, NewInteger(ino));
					ino++;
					}
				goto cnt7;
				}
			if(CompoundArg1(sp)==A_GAMMA5 || CompoundArg1(sp)==A_GAMMAP ||
					 CompoundArg1(sp)==A_GAMMAM || CompoundArg1(sp)==A_GAMMA)
				{
				int ii;
				ii=(int)IntegerValue(ListFirst(CompoundArg2(sp)));
				if(ii>ino)
					{
					puts("Internal error twl4g5.");
					goto aaa;
					}
				if(ii==ino)
					{
					ibuf[0]=sino++;
					ind=AppendLast(ind, NewAtom(ibuf,1));
					ino++;
					}
				ii=(int)IntegerValue(ListNth(CompoundArg2(sp),2));
				if(ii>ino)
					{
					puts("Internal error twl4g.");
					goto aaa;
					}
				if(ii==ino)
					{
					ibuf[0]=sino++;
					ind=AppendLast(ind, NewAtom(ibuf,1));
					ino++;
					}
				if(CompoundArg1(sp)==A_GAMMA)
					{
					ii=(int)IntegerValue(ListNth(CompoundArg2(sp),3));
					if(ii>ino)
						{
						puts("Internal error twl4gv.");
						goto aaa;
						}
					if(ii==ino)
						{
						ind=AppendLast(ind, NewAtom(vinarr[vino++],0));
						vind=AppendLast(vind, NewInteger(ino));
						ino++;
						}
					}
				goto cnt7;
				}
			l=CompoundArg2(sp);
			while(!is_empty_list(l))
				{
				int ii;
				ii=(int)IntegerValue(ListFirst(l));
				if(ii>ino)
					{
					puts("Internal error twl4s.");
					WriteTerm(CompoundArg1(a2));puts("");
					WriteTerm(CompoundArgN(a2,4));puts("");
					WriteTerm(ListFirst(CompoundArgN(a2,5)));puts("");
					WriteTerm(ListFirst(ListTail(CompoundArgN(a2,5))));puts("");
					goto aaa;
					}
				if(ii==ino)
					{
					if(CompoundArg1(sp)==A_DELTA)
						puts("Internal error (twl4d).");
					ibuf[0]=oino++;
					ind=AppendLast(ind, NewAtom(ibuf,1));
					ino++;
					}
				l=ListTail(l);
				}
		cnt7:
			l2=ListTail(l2);
			}



	l2=CompoundArgN(a2,4);
		while(!is_empty_list(l2))
			{
			Term sp;
			sp=ListFirst(l2);
			if(CompoundName(sp)==A_MOMENT)
				{
				fprintf(f,"p_%ld^%s ",IntegerValue(CompoundArg1(sp)),
					AtomValue(ListNth(ind,(int)IntegerValue(ListFirst(CompoundArg2(sp))))));
				silcnt++;
				goto cnt8;
				}
			if(CompoundArg1(sp)==A_GAMMA)
				{
				fprintf(f,"\\gamma_{%s %s}^%s ",
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),3)))));
				silcnt+=2;
				goto cnt8;
				}
			if(CompoundArg1(sp)==A_GAMMA5)
				{
				fprintf(f,"\\gamma_{%s %s}^5 ",
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),2)))));
				silcnt+=2;
				goto cnt8;
				}
			if(CompoundArg1(sp)==A_GAMMAP)
				{
				fprintf(f,"{(1+\\gamma^5)_{%s %s}\\over 2} ",
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),2)))));
				silcnt+=6;
				goto cnt8;
				}
			if(CompoundArg1(sp)==A_GAMMAM)
				{
				fprintf(f,"{(1-\\gamma^5)_{%s %s}\\over 2} ",
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),2)))));
				silcnt+=6;
				goto cnt8;
				}
			if(CompoundArg1(sp)==A_DELTA)
				{
				int n1,n2;
				n1=(int)IntegerValue(ListNth(CompoundArg2(sp),1));
				n2=(int)IntegerValue(ListNth(CompoundArg2(sp),2));
				if(ListMember(vind,NewInteger(n1)) && ListMember(vind,NewInteger(n2)))
					fprintf(f,"g^{%s %s} ",
						AtomValue(ListNth(ind,n1)),
						AtomValue(ListNth(ind,n2)));
				else
					fprintf(f,"\\delta_{%s %s} ",
						AtomValue(ListNth(ind,n1)),
						AtomValue(ListNth(ind,n2)));
				silcnt+=2;
				goto cnt8;
				}
			if(CompoundArg1(sp)==A_EPS_V)
			{
				fprintf(f,"\\varepsilon_{%s %s %s %s} ",
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),3)))),
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),4)))));
				silcnt+=3;
				goto cnt8;
			}

			if(GetAtomProperty(CompoundArg1(sp),A_COLOR)==A_COLOR_LAMBDA)
				{
				fprintf(f,"\\lambda_{%s %s}^%s ",
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),3)))));
				silcnt+=2;
				goto cnt8;
				}
			if(GetAtomProperty(CompoundArg1(sp),A_COLOR)==A_COLOR_F)
				{
				fprintf(f,"f_{%s %s %s} ",
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),3)))));
				silcnt+=3;
				goto cnt8;
				}
			if(GetAtomProperty(CompoundArg1(sp),A_COLOR)==A_COLOR_D)
				{
				fprintf(f,"d_{%s %s %s} ",
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),3)))));
				silcnt+=3;
				goto cnt8;
				}
			if(GetAtomProperty(CompoundArg1(sp),A_COLOR)==A_COLOR_EPS)
				{
				fprintf(f,"\\varepsilon_{%s %s %s} ",
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind,(int)IntegerValue(ListNth(CompoundArg2(sp),3)))));
				silcnt+=3;
				goto cnt8;
				}	
			puts("Internal error (twl5).");
		cnt8:
			l2=ListTail(l2);
			}





	m2l=CompoundArgN(a2,5);
	if(ListLength(m2l)==1 && CompoundArgN(ListFirst(m2l),3)==0)
		{
		fprintf(f,"$");
		return lines;
		}

	if(ListLength(m2l)==1 )
		{
		if(TEX_set_dot && CompoundArgN(a2,3))
			fprintf(f,"\\cdot ");
		}
	else
		fprintf(f,"\\big(");

	fflag=1;
	l1=m2l;
	while(!is_empty_list(l1))
		{
		int ino1, vino1, sino1, oino1;
		Term m2;
		m2=ListFirst(l1);
		{
		int num;
		num=(int)IntegerValue(CompoundArg1(m2));
		if(num==1)
			{
			if(fflag)
				fflag=0;
			else
				fprintf(f,"+");
			if(is_empty_list(CompoundArg2(m2)) &&
		       is_empty_list(CompoundArgN(m2,3)))
		    	fprintf(f,"1");
			goto s1;
			}
		if(num==-1)
			{
			fflag=0;
			fprintf(f,"-");
			if(
		    is_empty_list(CompoundArg2(m2)) &&
		    is_empty_list(CompoundArgN(m2,3)))
		    	fprintf(f,"1");
			goto s1;
			}

		if(fflag)
			{
			silcnt+=fprintf(f,"%ld",IntegerValue(CompoundArg1(m2)));
			fflag=0;
			}
		else
			silcnt+=fprintf(f,"%+ld",IntegerValue(CompoundArg1(m2)));
		}
	s1:
		l2=CompoundArg2(m2);
		while(!is_empty_list(l2))
			{
			if(CompoundArg1(ListFirst(l2))==A_SQRT2)
				fprintf(f,"\\sqrt{2}");
			else
				{
				int po;
				Integer sle;
				po=(int)IntegerValue(CompoundArg2(ListFirst(l2)));
				wrt_param(f,CompoundArg1(ListFirst(l2)));
				if(po>1)
					fprintf(f,"{}^%d ",po);
				sle=GetAtomProperty(CompoundArg1(ListFirst(l2)),A_TEXLENGTH);
				if(sle)
					silcnt+=2*(int)IntegerValue(sle);
				else
					silcnt+=2;
				}
			l2=ListTail(l2);
			if(TEX_set_dot && (l || CompoundArgN(m2,3)))
				fprintf(f, "\\cdot ");
			}


		ind1=CopyTerm(ind);
		vind1=CopyTerm(vind);
		ino1=ino;
		vino1=vino;
		oino1=oino;
		sino1=sino;


		l2=CompoundArgN(m2,3);
		fix_ino(l2,ino);
		while(!is_empty_list(l2))
			{
			Term sp;
			sp=ListFirst(l2);
			if(CompoundName(sp)==A_MOMENT)
				{
				int ii;
				ii=(int)IntegerValue(ListFirst(CompoundArg2(sp)));
				if(ii>ino1)
					{
					puts("Internal error twl3m.");
					WriteVertex(CompoundArg1(a2));puts("");
					WriteTerm(m2);puts("\n");
					goto aaa;
					}
				if(ii==ino1)
					{
					ind1=AppendLast(ind1, NewAtom(vinarr[vino1++],0));
					vind1=AppendLast(vind1, NewInteger(ino1));
					ino1++;
					}
				goto cnt5;
				}
			if(CompoundArg1(sp)==A_GAMMA5 || CompoundArg1(sp)==A_GAMMAP ||
					 CompoundArg1(sp)==A_GAMMAM || CompoundArg1(sp)==A_GAMMA)
				{
				int ii;
				ii=(int)IntegerValue(ListFirst(CompoundArg2(sp)));
				if(ii>ino1)
					{
					puts("Internal error twl3g5.");
					goto aaa;
					}
				if(ii==ino1)
					{
					ibuf[0]=sino1++;
					ind1=AppendLast(ind1, NewAtom(ibuf,1));
					ino1++;
					}
				ii=(int)IntegerValue(ListNth(CompoundArg2(sp),2));
				if(ii>ino1)
					{
					puts("Internal error twl3g.");
					goto aaa;
					}
				if(ii==ino1)
					{
					ibuf[0]=sino1++;
					ind1=AppendLast(ind1, NewAtom(ibuf,1));
					ino1++;
					}
				if(CompoundArg1(sp)==A_GAMMA)
					{
					ii=IntegerValue(ListNth(CompoundArg2(sp),3));
					if(ii>ino1)
						{
						puts("Internal error twl3gv.");
						goto aaa;
						}
					if(ii==ino1)
						{
						ind1=AppendLast(ind1, NewAtom(vinarr[vino1++],0));
						vind1=AppendLast(vind1, NewInteger(ino1));
						ino1++;
						}
					}
				goto cnt5;
				}
			l=CompoundArg2(sp);
			while(!is_empty_list(l))
				{
				int ii;
				ii=(int)IntegerValue(ListFirst(l));
				if(ii>ino1)
					{
					puts("Internal error twl3s.");
					goto aaa;
					}
				if(ii==ino1)
					{
					if(CompoundArg1(sp)==A_DELTA)
						puts("Internal error (twl3d).");
					ibuf[0]=oino1++;
					ind1=AppendLast(ind1, NewAtom(ibuf,1));
					ino1++;
					}
				l=ListTail(l);
				}
		cnt5:
			l2=ListTail(l2);
			}


		l2=CompoundArgN(m2,3);
		while(!is_empty_list(l2))
			{
			Term sp;
			sp=ListFirst(l2);
			if(CompoundName(sp)==A_MOMENT)
				{
				fprintf(f,"p_%ld^%s ",IntegerValue(CompoundArg1(sp)),
					AtomValue(ListNth(ind1,IntegerValue(ListFirst(CompoundArg2(sp))))));
				silcnt++;
				goto cnt6;
				}
			if(CompoundArg1(sp)==A_GAMMA)
				{
				fprintf(f,"\\gamma_{%s %s}^%s ",
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),3)))));
				silcnt+=2;
				goto cnt6;
				}
			if(CompoundArg1(sp)==A_GAMMA5)
				{
				fprintf(f,"\\gamma_{%s %s}^5 ",
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),2)))));
				silcnt+=2;
				goto cnt6;
				}
			if(CompoundArg1(sp)==A_GAMMAP)
				{
				fprintf(f,"{(1+\\gamma^5)_{%s %s}\\over 2} ",
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),2)))));
				silcnt+=6;
				goto cnt6;
				}
			if(CompoundArg1(sp)==A_GAMMAM)
				{
				fprintf(f,"{(1-\\gamma^5)_{%s %s}\\over 2} ",
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),2)))));
				silcnt+=6;
				goto cnt6;
				}
			if(CompoundArg1(sp)==A_DELTA)
				{
				int n1,n2;
				n1=IntegerValue(ListNth(CompoundArg2(sp),1));
				n2=IntegerValue(ListNth(CompoundArg2(sp),2));
				if(ListMember(vind1,NewInteger(n1)) && ListMember(vind1,NewInteger(n2)))
					fprintf(f,"g^{%s %s} ",
						AtomValue(ListNth(ind1,n1)),
						AtomValue(ListNth(ind1,n2)));
				else
					fprintf(f,"\\delta_{%s %s} ",
						AtomValue(ListNth(ind1,n1)),
						AtomValue(ListNth(ind1,n2)));
				silcnt+=2;
				goto cnt6;
				}
			if(CompoundArg1(sp)==A_EPS_V)
			{
				fprintf(f,"\\varepsilon_{%s %s %s %s} ",
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),3)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),4)))));
				silcnt+=3;
				goto cnt6;
			}

			if(GetAtomProperty(CompoundArg1(sp),A_COLOR)==A_COLOR_LAMBDA)
				{
				fprintf(f,"\\lambda_{%s %s}^%s ",
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),3)))));
//				silcnt+=2;
				goto cnt6;
				}
			if(GetAtomProperty(CompoundArg1(sp),A_COLOR)==A_COLOR_F)
				{
				fprintf(f,"f_{%s %s %s} ",
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),3)))));
				silcnt+=3;
				goto cnt6;
				}
			if(GetAtomProperty(CompoundArg1(sp),A_COLOR)==A_COLOR_D)
				{
				fprintf(f,"d_{%s %s %s} ",
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),3)))));
				silcnt+=3;
				goto cnt6;
				}
			if(GetAtomProperty(CompoundArg1(sp),A_COLOR)==A_COLOR_EPS)
				{
				fprintf(f,"\\varepsilon_{%s %s %s} ",
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),3)))));
				silcnt+=3;
				goto cnt6;
				}
			{
				int  iii;
				Atom NM;
				NM=GetAtomProperty(CompoundArg1(sp),A_TEXNAME);
				if(NM==0)
					NM=CompoundArg1(sp);
				fprintf(f,"%s_{",AtomValue(NM));
				for(iii=1;iii<=ListLength(CompoundArg2(sp));iii++)
					fprintf(f,"%s ",
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),iii)))));
				silcnt+=ListLength(CompoundArg2(sp));
				fprintf(f,"} ");
				goto cnt6;
			}
			puts("Internal error (twl6).");
		cnt6:
			l2=ListTail(l2);
			}

		FreeAtomic(ind1);
		FreeAtomic(vind1);
		if(silcnt>=TEX_spec_in_line && ListTail(l1) && lines+cur_cnt>=TEX_lines)
		{
			cur_cnt=0;
			lines=0;
			fprintf(f,"$");
			fprintf(f,"\\\\ \\hline\n\\end{tabular}\n\n");
			fprintf(f,"\\begin{tabular}{|l|l|} \\hline\n");
			fprintf(f,"Fields in the vertex & ");
			fprintf(f,"Variational derivative of Lagrangian by fields \\\\ \\hline\n");
			fprintf(f," & $");
		
		}
		if(silcnt>=TEX_spec_in_line && ListTail(l1))
			{
			silcnt=0;
			fprintf(f,"$ \\\\[2mm]\n  & $");
			lines++;
			}

		l1=ListTail(l1);
		}

	if(ListLength(m2l)>1)
		fprintf(f,"\\big)");

aaa:
	fprintf(f,"$");

	return lines;

	}

extern int write_all_vertices;

List opTexBOF=0, opTexEOF=0;

void tex_write_lagr(List l, FILE *f)
	{
	List l1;
	char cbuf[64];
	int pn, pni;
	
	pni=0;
	pn=ListLength(l);
		
	cur_cnt=0;
	
	for(l1=opTexBOF;l1;l1=ListTail(l1))
	{
		if(is_atom(ListFirst(l1)))
			fprintf(f,"%s\n",AtomValue(ListFirst(l1)));
		else
			fprintf(f,"\n");
	}
	
	
	l1=l;

	fprintf(f,"\\begin{tabular}{|l|l|} \\hline\n");
	fprintf(f,"Fields in the vertex & ");
	fprintf(f,"Variational derivative of Lagrangian by fields \\\\ \\hline\n");

	while(!is_empty_list(l1))
		{
		Term a2;
		
		pni++;
		
		a2=ListFirst(l1);
		if(
		  (ListLength(CompoundArg1(a2))>2 || write_all_vertices)
			  && ListLength(CompoundArg1(a2))<=TEX_max_pno
		      && CompoundArgN(a2,5))
			{
			sprintf(cbuf,"Writing lagrangian line %d of %d.\n",pni,pn);
			RegisterLine(cbuf);
			alg2_symmetrize(a2);
			alg2_common_n(a2);
			alg2_common_s(a2);
			alg2_common_t(a2);
			/*alg2_red_cos(a2);
			alg2_red_orth(a2);
			alg2_red_sico(a2);
			alg2_red_comsico(a2);*/
			alg2_red_1pm5(a2);
			alg2_recommon_n(a2);
			alg2_recommon_s(a2);
			if(MultByI) alg2_multbyi(a2);
			if(CompoundArgN(a2,5))
				{
				if(cur_cnt) fprintf(f,"\\\\[2mm]\n");
				cur_cnt+=twl_1(f, a2);
				if(cur_cnt>=TEX_lines && ListTail(l1))
					{
					cur_cnt=0;
					fprintf(f,"\\\\ \\hline\n\\end{tabular}\n\n");
					fprintf(f,"\\begin{tabular}{|l|l|} \\hline\n");
					fprintf(f,"Fields in the vertex & ");
					fprintf(f,"Variational derivative of Lagrangian by fields \\\\ \\hline\n");
					}
				/*else
					{
					if(ListTail(l1)) fprintf(f,"\\\\[2mm]\n");
					}*/
				}
			l1=ListTail(l1);
			UnregisterLine();
			}
		else
			{
			l1=ListTail(l1);
			}
		}
	fprintf(f,"\\\\ \\hline\n\\end{tabular}\n");
	FreeAtomic(l);
	
	for(l1=opTexEOF;l1;l1=ListTail(l1))
	{
		if(is_atom(ListFirst(l1)))
			fprintf(f,"%s\n",AtomValue(ListFirst(l1)));
		else
			fprintf(f,"\n");
	}
	
	
	}




void tex_wrt_2vrt(FILE *f, Term a2)
	{
	List prt, ind, ind1, vind, vind1, l, l1, l2, m2l;
	int ino=1, vino=0, sino='a', oino='p', fflag;
	prt=CompoundArg1(a2);
	l=prt;
/*	WriteTerm(a2);
	puts("\n");*/
	ind=vind=NewList();
	while(!is_empty_list(l))
		{
		Term t;
		Atom n,n1;

		t=CompoundArg2(ListFirst(l));
		n=CompoundArg1(ListFirst(l));
		n1=GetAtomProperty(n,A_TEXNAME);
/*		if(n1)
			n=n1;*/
		l1=CompoundArg1(t);
		if(l1)
			{
			if(IntegerValue(ListFirst(l1))!=ino)
				{
				puts("Internal error (twv1).");
				return ;
				}

			if(CompoundName(t)==OPR_VECTOR)
				{
				ind=AppendLast(ind, NewAtom(vinarr[vino++],0));
				vind=AppendLast(vind, NewInteger(ino));
				ino++;
				if(ListLength(l1)==2)
					{
					ibuf[0]=oino++;
					ind=AppendLast(ind,NewAtom(ibuf,1));
					ino++;
					}
				if(ListLength(l1)>2)
					{
					puts("Internal error (twv2v).");
					return ;
					}
				goto cnt1;
				}

			if(CompoundName(t)==OPR_SPINOR)
				{
				ibuf[0]=sino++;
				ind=AppendLast(ind, NewAtom(ibuf,1));
				ino++;
				if(ListLength(l1)==2)
					{
					ibuf[0]=oino++;
					ind=AppendLast(ind,NewAtom(ibuf,1));
					ino++;
					}
				if(ListLength(l1)>2)
					{
					puts("Internal error (twv2sp).");
					return;
					}
				goto cnt1;
				}

			if(CompoundName(t)==OPR_SPINOR3)
				{
				ibuf[0]=sino++;
				ind=AppendLast(ind, NewAtom(ibuf,1));
				ino++;
				
				ind=AppendLast(ind, NewAtom(vinarr[vino++],0));
				vind=AppendLast(vind, NewInteger(ino));
				ino++;
				
				if(ListLength(l1)==3)
					{
					ibuf[0]=oino++;
					ind=AppendLast(ind,NewAtom(ibuf,1));
					ino++;
					}
				if(ListLength(l1)>3)
					{
					puts("Internal error (twv2sp).");
					return;
					}
				goto cnt1;
				}

			if(CompoundName(t)==OPR_TENSOR)
				{
				ind=AppendLast(ind, NewAtom(vinarr[vino++],0));
				vind=AppendLast(vind, NewInteger(ino));
				ino++;
				ind=AppendLast(ind, NewAtom(vinarr[vino++],0));
				vind=AppendLast(vind, NewInteger(ino));
				ino++;
				if(ListLength(l1)==3)
					{
					ibuf[0]=oino++;
					ind=AppendLast(ind,NewAtom(ibuf,1));
					ino++;
					}
				if(ListLength(l1)>3)
					{
					puts("Internal error (twv2t).");
					return;
					}
				goto cnt1;
				}

			if(CompoundName(t)==OPR_SCALAR)
				{
				if(ListLength(l1)==1)
					{
					ibuf[0]=oino++;
					ind=AppendLast(ind,NewAtom(ibuf,1));
					ino++;
					}
				if(ListLength(l1)>1)
					{
					puts("Internal error (twv2sc).");
					return;
					}
				goto cnt1;
				}
			puts("Internal error (twv2).");
			return;
			}
	cnt1:
		/*if(AtomValue(n)[0]=='~')
			fprintf(f,"\\~{}${%s}^{",AtomValue(n)+1);
		else
			fprintf(f,"${%s}_{",AtomValue(n));*/
		/*wrt_prt_name(f,n);
		while(!is_empty_list(l1))
			{
			fprintf(f,"%s ",AtomValue(ListNth(ind,IntegerValue(ListFirst(l1)))));
			l1=ListTail(l1);
			}
		fprintf(f,"}$ & ");*/
		l=ListTail(l);
		}

/*
	if(ListLength(prt)==3)
		fprintf(f," & ");
*/

	fprintf(f, "$");
	wrt_cnf(f, CompoundArg2(a2));
	wrt_csf(f, CompoundArgN(a2,3));

	l2=CompoundArgN(a2,4);
		while(!is_empty_list(l2))
			{
			Term sp;
			sp=ListFirst(l2);
			if(CompoundName(sp)==A_MOMENT)
				{
				int ii;
				ii=IntegerValue(ListFirst(CompoundArg2(sp)));
				if(ii>ino)
					{
					puts("Internal error twv4m.");
					goto aaa;
					}
				if(ii==ino)
					{
					ind=AppendLast(ind, NewAtom(vinarr[vino++],0));
					vind=AppendLast(vind, NewInteger(ino));
					ino++;
					}
				goto cnt7;
				}
			if(CompoundArg1(sp)==A_GAMMA5 || CompoundArg1(sp)==A_GAMMAP ||
					 CompoundArg1(sp)==A_GAMMAM || CompoundArg1(sp)==A_GAMMA)
				{
				int ii;
				ii=IntegerValue(ListFirst(CompoundArg2(sp)));
				if(ii>ino)
					{
					puts("Internal error twv4g5.");
					goto aaa;
					}
				if(ii==ino)
					{
					ibuf[0]=sino++;
					ind=AppendLast(ind, NewAtom(ibuf,1));
					ino++;
					}
				ii=IntegerValue(ListNth(CompoundArg2(sp),2));
				if(ii>ino)
					{
					puts("Internal error twv4g.");
					goto aaa;
					}
				if(ii==ino)
					{
					ibuf[0]=sino++;
					ind=AppendLast(ind, NewAtom(ibuf,1));
					ino++;
					}
				if(CompoundArg1(sp)==A_GAMMA)
					{
					ii=IntegerValue(ListNth(CompoundArg2(sp),3));
					if(ii>ino)
						{
						puts("Internal error twv4gv.");
						goto aaa;
						}
					if(ii==ino)
						{
						ind=AppendLast(ind, NewAtom(vinarr[vino++],0));
						vind=AppendLast(vind, NewInteger(ino));
						ino++;
						}
					}
				goto cnt7;
				}
			l=CompoundArg2(sp);
			while(!is_empty_list(l))
				{
				int ii;
				ii=IntegerValue(ListFirst(l));
				if(ii>ino)
					{
					puts("Internal error twv4s.");
					goto aaa;
					}
				if(ii==ino)
					{
					if(CompoundArg1(sp)==A_DELTA)
						puts("Internal error (twv4d).");
					ibuf[0]=oino++;
					ind=AppendLast(ind, NewAtom(ibuf,1));
					ino++;
					}
				l=ListTail(l);
				}
		cnt7:
			l2=ListTail(l2);
			}



	l2=CompoundArgN(a2,4);
		while(!is_empty_list(l2))
			{
			Term sp;
			sp=ListFirst(l2);
			if(CompoundName(sp)==A_MOMENT)
				{
				fprintf(f,"p_%ld^%s ",IntegerValue(CompoundArg1(sp)),
					AtomValue(ListNth(ind,IntegerValue(ListFirst(CompoundArg2(sp))))));
				goto cnt8;
				}
			if(CompoundArg1(sp)==A_GAMMA)
				{
				fprintf(f,"\\gamma_{%s %s}^%s ",
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),3)))));
				goto cnt8;
				}
			if(CompoundArg1(sp)==A_GAMMA5)
				{
				fprintf(f,"\\gamma_{%s %s}^5 ",
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),2)))));
				goto cnt8;
				}
			if(CompoundArg1(sp)==A_GAMMAP)
				{
				fprintf(f,"{(1+\\gamma^5)_{%s %s}\\over 2} ",
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),2)))));
				goto cnt8;
				}
			if(CompoundArg1(sp)==A_GAMMAM)
				{
				fprintf(f,"{(1-\\gamma^5)_{%s %s}\\over 2} ",
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),2)))));
				goto cnt8;
				}
			if(CompoundArg1(sp)==A_DELTA)
				{
				int n1,n2;
				n1=IntegerValue(ListNth(CompoundArg2(sp),1));
				n2=IntegerValue(ListNth(CompoundArg2(sp),2));
				if(ListMember(vind,NewInteger(n1)) && ListMember(vind,NewInteger(n2)))
					fprintf(f,"g^{%s %s} ",
						AtomValue(ListNth(ind,n1)),
						AtomValue(ListNth(ind,n2)));
				else
					fprintf(f,"\\delta_{%s %s} ",
						AtomValue(ListNth(ind,n1)),
						AtomValue(ListNth(ind,n2)));
				goto cnt8;
				}
				
			if(CompoundArg1(sp)==A_EPS_V)
			{
				fprintf(f,"\\varepsilon_{%s %s %s %s} ",
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),3)))),
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),4)))));
				goto cnt8;
			}

			if(GetAtomProperty(CompoundArg1(sp),A_COLOR)==A_COLOR_LAMBDA)
				{
				fprintf(f,"\\lambda_{%s %s}^%s ",
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),3)))));
				goto cnt8;
				}
			if(GetAtomProperty(CompoundArg1(sp),A_COLOR)==A_COLOR_F)
				{
				fprintf(f,"f_{%s %s %s} ",
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),3)))));
				goto cnt8;
				}
			if(GetAtomProperty(CompoundArg1(sp),A_COLOR)==A_COLOR_D)
				{
				fprintf(f,"d_{%s %s %s} ",
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),3)))));
				goto cnt8;
				}
			if(GetAtomProperty(CompoundArg1(sp),A_COLOR)==A_COLOR_EPS)
				{
				fprintf(f,"\\varepsilon_{%s %s %s} ",
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),3)))));
				goto cnt8;
				}
			puts("Internal error (twv5).");
		cnt8:
			l2=ListTail(l2);
			}





	m2l=CompoundArgN(a2,5);
	if(ListLength(m2l)==1 && CompoundArgN(ListFirst(m2l),3)==0)
		{
		fprintf(f,"$");
		return;
		}

	if(ListLength(m2l)>1)
		fprintf(f,"\\left(");
	else
		{
		if(TEX_set_dot && CompoundArgN(a2,3))
			fprintf(f,"\\cdot ");
		}


	fflag=1;
	l1=m2l;
	while(!is_empty_list(l1))
		{
		int ino1, vino1, sino1, oino1;
		Term m2;
		m2=ListFirst(l1);
		{
		int num;
		num=IntegerValue(CompoundArg1(m2));
		if(num==1)
			{
			if(fflag)
				fflag=0;
			else
				fprintf(f,"+");
			goto s1;
			}
		if(num==-1)
			{
			fflag=0;
			fprintf(f,"-");
			goto s1;
			}
		if(fflag)
			{
			fprintf(f,"%ld",IntegerValue(CompoundArg1(m2)));
			fflag=0;
			}
		else
			fprintf(f,"%+ld",IntegerValue(CompoundArg1(m2)));
		}
	s1:
		l2=CompoundArg2(m2);
		while(!is_empty_list(l2))
			{
			if(CompoundArg1(ListFirst(l2))==A_SQRT2)
				fprintf(f,"\\sqrt{2}");
			else
				{
				int po;
				po=IntegerValue(CompoundArg2(ListFirst(l2)));
				wrt_param(f,CompoundArg1(ListFirst(l2)));
				if(po>1)
					fprintf(f,"{}^%d ",po);
				}
			l2=ListTail(l2);
			if(TEX_set_dot && (l || CompoundArgN(m2,3)))
				fprintf(f, "\\cdot ");
			}

		ind1=CopyTerm(ind);
		vind1=CopyTerm(vind);
		ino1=ino;
		vino1=vino;
		oino1=oino;
		sino1=sino;

		l2=CompoundArgN(m2,3);
		while(!is_empty_list(l2))
			{
			Term sp;
			sp=ListFirst(l2);
			if(CompoundName(sp)==A_MOMENT)
				{
				int ii;
				ii=IntegerValue(ListFirst(CompoundArg2(sp)));
				if(ii>ino1)
					{
					puts("Internal error twv3m.");
					goto aaa;
					}
				if(ii==ino1)
					{
					ind1=AppendLast(ind1, NewAtom(vinarr[vino1++],0));
					vind1=AppendLast(vind1, NewInteger(ino1));
					ino1++;
					}
				goto cnt5;
				}
			if(CompoundArg1(sp)==A_GAMMA5 || CompoundArg1(sp)==A_GAMMAP ||
					 CompoundArg1(sp)==A_GAMMAM || CompoundArg1(sp)==A_GAMMA)
				{
				int ii;
				ii=IntegerValue(ListFirst(CompoundArg2(sp)));
				if(ii>ino1)
					{
					puts("Internal error twv3g5.");
					goto aaa;
					}
				if(ii==ino1)
					{
					ibuf[0]=sino1++;
					ind1=AppendLast(ind1, NewAtom(ibuf,1));
					ino1++;
					}
				ii=IntegerValue(ListNth(CompoundArg2(sp),2));
				if(ii>ino1)
					{
					puts("Internal error twv3g.");
					goto aaa;
					}
				if(ii==ino1)
					{
					ibuf[0]=sino1++;
					ind1=AppendLast(ind1, NewAtom(ibuf,1));
					ino1++;
					}
				if(CompoundArg1(sp)==A_GAMMA)
					{
					ii=IntegerValue(ListNth(CompoundArg2(sp),3));
					if(ii>ino1)
						{
						puts("Internal error twv3gv.");
						goto aaa;
						}
					if(ii==ino1)
						{
						ind1=AppendLast(ind1, NewAtom(vinarr[vino1++],0));
						vind1=AppendLast(vind1, NewInteger(ino1));
						ino1++;
						}
					}
				goto cnt5;
				}
			l=CompoundArg2(sp);
			while(!is_empty_list(l))
				{
				int ii;
				ii=IntegerValue(ListFirst(l));
				if(ii>ino1)
					{
					puts("Internal error twv3s.");
					goto aaa;
					}
				if(ii==ino1)
					{
					if(CompoundArg1(sp)==A_DELTA)
						puts("Internal error (twv3d).");
					ibuf[0]=oino1++;
					ind1=AppendLast(ind1, NewAtom(ibuf,1));
					ino1++;
					}
				l=ListTail(l);
				}
		cnt5:
			l2=ListTail(l2);
			}


		l2=CompoundArgN(m2,3);
		while(!is_empty_list(l2))
			{
			Term sp;
			sp=ListFirst(l2);
			if(CompoundName(sp)==A_MOMENT)
				{
				fprintf(f,"p_%ld^%s ",IntegerValue(CompoundArg1(sp)),
					AtomValue(ListNth(ind1,IntegerValue(ListFirst(CompoundArg2(sp))))));
				goto cnt6;
				}
			if(CompoundArg1(sp)==A_GAMMA)
				{
				fprintf(f,"\\gamma_{%s %s}^%s ",
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),3)))));
				goto cnt6;
				}
			if(CompoundArg1(sp)==A_GAMMA5)
				{
				fprintf(f,"\\gamma_{%s %s}^5 ",
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),2)))));
				goto cnt6;
				}
			if(CompoundArg1(sp)==A_GAMMAP)
				{
				fprintf(f,"{(1+\\gamma^5)_{%s %s}\\over 2} ",
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),2)))));
				goto cnt6;
				}
			if(CompoundArg1(sp)==A_GAMMAM)
				{
				fprintf(f,"{(1-\\gamma^5)_{%s %s}\\over 2} ",
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),2)))));
				goto cnt6;
				}
			if(CompoundArg1(sp)==A_DELTA)
				{
				int n1,n2;
				n1=IntegerValue(ListNth(CompoundArg2(sp),1));
				n2=IntegerValue(ListNth(CompoundArg2(sp),2));
				if(ListMember(vind1,NewInteger(n1)) && ListMember(vind1,NewInteger(n2)))
					fprintf(f,"g^{%s %s} ",
						AtomValue(ListNth(ind1,n1)),
						AtomValue(ListNth(ind1,n2)));
				else
					fprintf(f,"\\delta_{%s %s} ",
						AtomValue(ListNth(ind1,n1)),
						AtomValue(ListNth(ind1,n2)));
				goto cnt6;
				}
				
			if(CompoundArg1(sp)==A_EPS_V)
			{
				fprintf(f,"\\varepsilon_{%s %s %s %s} ",
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),3)))),
					AtomValue(ListNth(ind,IntegerValue(ListNth(CompoundArg2(sp),4)))));
				goto cnt6;
			}

			if(GetAtomProperty(CompoundArg1(sp),A_COLOR)==A_COLOR_LAMBDA)
				{
				fprintf(f,"\\lambda_{%s %s}^%s ",
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),3)))));
				goto cnt6;
				}
			if(GetAtomProperty(CompoundArg1(sp),A_COLOR)==A_COLOR_F)
				{
				fprintf(f,"f_{%s %s %s} ",
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),3)))));
				goto cnt6;
				}
			if(GetAtomProperty(CompoundArg1(sp),A_COLOR)==A_COLOR_D)
				{
				fprintf(f,"d_{%s %s %s} ",
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),3)))));
				goto cnt6;
				}
			if(GetAtomProperty(CompoundArg1(sp),A_COLOR)==A_COLOR_EPS)
				{
				fprintf(f,"\\varepsilon_{%s %s %s} ",
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),1)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),2)))),
					AtomValue(ListNth(ind1,IntegerValue(ListNth(CompoundArg2(sp),3)))));
				goto cnt6;
				}
			puts("Internal error (twv6).");
		cnt6:
			l2=ListTail(l2);
			}

		FreeAtomic(ind1);
		FreeAtomic(vind1);

		l1=ListTail(l1);
		}

	if(ListLength(m2l)>1)
		fprintf(f,"\\right)");

aaa:
	fprintf(f,"$");



	}
	
