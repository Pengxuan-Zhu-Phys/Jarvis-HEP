#include <string.h>
#include <stdio.h>
#include "lanhep.h"

int off_srefine=0, ch_sign=0;
int opReduceG5=0;

static void inv_gamma(List);

int alg2_refine_spinor(List prt, List *spec)
	{
	int vtype;
	int retf=1;
	int i;
	List gl,l1;
	Label lab1,lab2,lab0;
	Label clab1=0, clab2=0;
	
/*	WriteTerm(prt);
	WriteTerm(*spec);
	puts("");
*/	
	vtype=ListLength(prt);
	if(vtype!=2 && vtype!=4)
		{
		printf("Error: vertex with anomalous fermion number:");
		WriteVertex(prt);
		puts("");
		return 0;
		}


	if(vtype==2)
		{
		int ity[2];
		for(i=0;i<2;i++)
		ity[i]=(int)IntegerValue(CompoundArg1(CompoundArg1(ListFirst(
		    GetAtomProperty(CompoundArg1(ListNth(prt,i+1)),PROP_INDEX) ))));
		if(ity[0]==0 && ity[1]==1)
			goto aaa1;
		if(ity[0]==1 && ity[1]==0)
			{
			Term p1,p2;
			p1=ListFirst(prt);
			p2=ListFirst(ListTail(prt));
			ChangeList(prt,p2);
			ChangeList(ListTail(prt),p1);
			retf=-1;
			goto aaa1;
			}
		printf("Internal error (a2s01)");
		return 0;
		}

	if(vtype==4)
		{
		List l1,l2;
		int flag;

		do
			{
			flag=0;
			l1=prt;
			l2=ListTail(l1);
			while(!is_empty_list(l2))
				{
				Atom a1, a2;
				a1=ListFirst(l1);
				a2=ListFirst(l2);
				if(strcmp(AtomValue(CompoundArg1(a1)),
							AtomValue(CompoundArg1(a2)))==1)
					{
					ChangeList(l1,a2);
					ChangeList(l2,a1);
					flag=1;
					retf*=-1;
					}
			l1=l2;
			l2=ListTail(l2);
			}
		}
		while(flag);
		return retf;
		}

aaa1:

/*
	WriteTerm(prt); puts("");
	DumpList(*spec);
	puts(" - refine start");
*/

	if(off_srefine)
		return retf;


	lab1=CompoundArg2(ListFirst(prt));
	l1=*spec;
	while(!is_empty_list(l1))
		{
		Term y1;
		y1=ListFirst(l1);
		if((CompoundName(y1)==OPR_SPINOR || CompoundName(y1)==OPR_SPINOR3)
					&& CompoundArg1(y1)==lab1)
			{
			List uu;
			uu=CompoundArg2(y1);
			lab1=ListFirst(uu);
			if(ListLength(uu)>1)
				clab1=ListNth(uu,2);
			break;
			}
		l1=ListTail(l1);
		}
	if(is_empty_list(l1))
		{
		puts("Internal error (a2s04)");
		return 0;
		}
	lab2=CompoundArg2(ListFirst(ListTail(prt)));
	l1=*spec;
	while(!is_empty_list(l1))
		{
		Term y1;
		y1=ListFirst(l1);
		if((CompoundName(y1)==OPR_SPINOR || CompoundName(y1)==OPR_SPINOR3) 
				&& CompoundArg1(y1)==lab2)
			{
			List uu;
			uu=CompoundArg2(y1);
			lab2=ListFirst(uu);
			if(ListLength(uu)>1)
				clab2=ListNth(uu,2);
			break;
			}
		l1=ListTail(l1);
		}
	if(is_empty_list(l1))
		{
		puts("Internal error (a2s04a)");
		return 0;
		}


/*WriteTerm(lab1); printf(" "); WriteTerm(lab2); puts("");*/

	gl=NewList();
	lab0=lab1;
	while(lab0!=lab2)
		{
		l1=*spec;
		while(!is_empty_list(l1))
			{
			Term y1;
			y1=ListFirst(l1);
			if(CompoundName(y1)==OPR_SPECIAL &&
			   CompoundArg1(y1)==A_GAMMA &&
			   ListFirst(CompoundArg2(y1))==lab0)
			   	{
			   	lab0=ListNth(CompoundArg2(y1),2);
			   	gl=AppendLast(gl,y1);
			   	ChangeList(l1,0);
			   	*spec=CutFromList(*spec,l1);
			   	goto aa154;
			   	}
			if(CompoundName(y1)==OPR_SPECIAL &&
			   (CompoundArg1(y1)==A_GAMMA5 || CompoundArg1(y1)==A_GAMMAP ||
					CompoundArg1(y1)==A_GAMMAM) &&
			   ListFirst(CompoundArg2(y1))==lab0)
			   	{
			   	lab0=ListNth(CompoundArg2(y1),2);
			   	gl=AppendLast(gl,y1);
			   	ChangeList(l1,0);
			   	*spec=CutFromList(*spec,l1);
			   	goto aa154;
			   	}
			l1=ListTail(l1);
			}
		puts("Internal error (a2s05)");
	WriteTerm(prt); puts("");
	DumpList(*spec);
	break;
    aa154: ;
		}

/*DumpList(gl);*/

	{
	Term prp,prp1;
	int dist=0;
	prp=GetAtomProperty(CompoundArg1(ListFirst(prt)),PROP_TYPE);
	if(!is_compound(prp) || CompoundName(prp)!=OPR_PARTICLE)
		goto lfg52;
	prp=CompoundArgN(prp,7);
	if(prp!=A_LEFT && prp!=A_RIGHT)
		goto lfg52;
	prp1=GetAtomProperty(CompoundArg1(ListFirst(ListTail(prt))),PROP_TYPE);
	if(!is_compound(prp1) || CompoundName(prp1)!=OPR_PARTICLE)
		goto lfg52;
	prp1=CompoundArgN(prp1,7);
	if(prp1!=A_LEFT && prp1!=A_RIGHT)
		goto lfg52;
	l1=gl;
	while(!is_empty_list(l1))
		{
		if(CompoundArg1(ListFirst(l1))==A_GAMMA)
			dist++;
	 	l1=ListTail(l1);
	 	}
	 if(prp1==prp && (dist%2)==0)
	 	return 0;
	 if(prp1!=prp && (dist%2)==1)
	 	return 0;
	}

lfg52:

	if(is_empty_list(gl))
		goto checkcc;


	l1=gl;
	while(!is_empty_list(l1))
		{
		List l2;
		int dist=0;
		Atom thisg,thatg;
		thisg=CompoundArg1(ListFirst(l1));
		if(thisg==A_GAMMA5 || thisg==A_GAMMAM || thisg==A_GAMMAP)
			{
			l2=ListTail(l1);
			while(!is_empty_list(l2))
				{
				thatg=CompoundArg1(ListFirst(l2));
				if(thatg==A_GAMMA5)
					{
					if(thisg==A_GAMMA5)
						{
						if(dist%2)
							retf*=-1;
						SetCompoundArg(ListFirst(l1),1,A_DELTA);
	 					SetCompoundArg(ListFirst(l2),1,A_DELTA);
	 					goto lfg52;
	 					}
					if(thisg==A_GAMMAP)
						{
						SetCompoundArg(ListFirst(l1),1,A_DELTA);
						SetCompoundArg(ListFirst(l2),1,A_GAMMAP);
						goto lfg52;
						}
					if(thisg==A_GAMMAM)
						{
						SetCompoundArg(ListFirst(l1),1,A_DELTA);
						SetCompoundArg(ListFirst(l2),1,A_GAMMAM);
						retf*=-1;
						goto lfg52;
						}
					puts("Internal error (a2s005)");
					}
				if(thatg==A_GAMMAM)
					{
					if(thisg==A_GAMMA5)
						{
						if((dist+1)%2)
							retf*=-1;
						SetCompoundArg(ListFirst(l1),1,A_DELTA);
						goto lfg52;
						}
					if(thisg==A_GAMMAM)
						{
						SetCompoundArg(ListFirst(l1),1,A_DELTA);
						goto lfg52;
						}
					if(thisg==A_GAMMAP)
						return 0;
					puts("Internal error (a2s005a)");
					}
				if(thatg==A_GAMMAP)
					{
					if(thisg==A_GAMMA5)
						{
						if(dist%2)
							retf*=-1;
						SetCompoundArg(ListFirst(l1),1,A_DELTA);
						goto lfg52;
						}
					if(thisg==A_GAMMAP)
						{
						SetCompoundArg(ListFirst(l1),1,A_DELTA);
						goto lfg52;
						}
					if(thisg==A_GAMMAM)
						return 0;
					puts("Internal error (a2s005b)");
					}			
					
	 			if(CompoundArg1(ListFirst(l2))==A_GAMMA)
					{
		 			dist++;
					if(thisg==A_GAMMAP)
						thisg=A_GAMMAM;
					else
						if(thisg==A_GAMMAM)
							thisg=A_GAMMAP;
					}
					
	 			l2=ListTail(l2);
	 			}
	 		break;
	 		}
	 	l1=ListTail(l1);
	 	}


		if(opReduceG5)
		{
		
			{
			Term prp;
			int dist=0;
			prp=GetAtomProperty(CompoundArg1(ListFirst(prt)),PROP_TYPE);
			if(!is_compound(prp))
			{
				puts("Internal error (a2s0a1)");
				return 0;
			}
			if(CompoundName(prp)==OPR_FIELD)
			{
				Term prp1,prp2;
				if(CompoundArg2(prp)!=NewInteger(4))
					goto lfpdc;
				prp1=CompoundArg1(prp);
				prp2=GetAtomProperty(prp1,PROP_TYPE);
				prp=CompoundArgN(prp2,7);
				if(prp!=A_LEFT && prp!=A_RIGHT)
					goto lfpdc;
				if(prp1==CompoundArg1(prp2))
				{
					if(prp==A_LEFT)
						prp=A_RIGHT;
					else
						prp=A_LEFT;
				}
			}
			else
				prp=CompoundArgN(prp,7);

			if(prp!=A_LEFT && prp!=A_RIGHT)
				goto lfpdc;
			l1=gl;
			while(!is_empty_list(l1))
				{
				Atom thisg;
				thisg=CompoundArg1(ListFirst(l1));
				if(thisg==A_GAMMA5)
					{
					if(dist%2)
						retf*=-1;
					SetCompoundArg(ListFirst(l1),1,A_DELTA);
					if(prp==A_RIGHT)
						retf*=-1;
	 				goto gmsh;
	 				}
				if(thisg==A_GAMMAP)
					{
					SetCompoundArg(ListFirst(l1),1,A_DELTA);
					if( ((dist%2)==0 && prp==A_RIGHT) ||
						((dist%2)==1 && prp==A_LEFT) )
						return 0;
	 				goto gmsh;
	 				}
				if(thisg==A_GAMMAM)
					{
					SetCompoundArg(ListFirst(l1),1,A_DELTA);
					if( ((dist%2)==1 && prp==A_RIGHT) ||
						((dist%2)==0 && prp==A_LEFT) )
						return 0;
	 				goto gmsh;
	 				}
				if(CompoundArg1(ListFirst(l1))==A_GAMMA)
	 				dist++;
				l1=ListTail(l1);
				}
			}

		lfpdc:

			{
			Term prp;
			prp=GetAtomProperty(CompoundArg1(ListFirst(ListTail(prt))),PROP_TYPE);
			if(!is_compound(prp))
			{
				puts("Internal error (a2s0a2)");
				return 0;
			}

			if(CompoundName(prp)==OPR_FIELD)
			{
				Term prp1,prp2;
				if(CompoundArg2(prp)!=NewInteger(4))
					goto gmsh;
				prp1=CompoundArg1(prp);
				prp2=GetAtomProperty(prp1,PROP_TYPE);
				prp=CompoundArgN(prp2,7);
				if(prp!=A_LEFT && prp!=A_RIGHT)
					goto gmsh;
				if(prp1!=CompoundArg1(prp2))
				{
					if(prp==A_LEFT)
						prp=A_RIGHT;
					else
						prp=A_LEFT;
				}
			}
			else
				prp=CompoundArgN(prp,7);

			if(prp!=A_LEFT && prp!=A_RIGHT)
				goto gmsh;
			l1=gl;
			while(!is_empty_list(l1))
				{
				Atom thisg;
				thisg=CompoundArg1(ListFirst(l1));
				if(thisg==A_GAMMA5 || thisg==A_GAMMAP || thisg==A_GAMMAM)
					{
					int dist=0;
					List l2;
					l2=ListTail(l1);
					while(!is_empty_list(l2))
						{
						if(CompoundArg1(ListFirst(l2))==A_GAMMA)
			 				dist++;
						l2=ListTail(l2);
						}

					if(thisg==A_GAMMA5)
						{
						if(dist%2)
							retf*=-1;
						SetCompoundArg(ListFirst(l1),1,A_DELTA);
						if(prp==A_LEFT)
							retf*=-1;
	 					goto gmsh;
						}
					if(thisg==A_GAMMAP)
						{
						SetCompoundArg(ListFirst(l1),1,A_DELTA);
						if( ((dist%2)==1 && prp==A_RIGHT) ||
							((dist%2)==0 && prp==A_LEFT) )
								return 0;
	 					goto gmsh;
	 					}
					if(thisg==A_GAMMAM)
						{
						SetCompoundArg(ListFirst(l1),1,A_DELTA);
						if( ((dist%2)==0 && prp==A_RIGHT) ||
							((dist%2)==1 && prp==A_LEFT) )
								return 0;
	 					goto gmsh;
	 					}
	 				}
				l1=ListTail(l1);
				}
			}
		
		}
gmsh:

	l1=gl;
	while(!is_empty_list(l1))
		{
		Atom thisg;
		thisg=CompoundArg1(ListFirst(l1));
		if(thisg==A_GAMMA5 || thisg==A_GAMMAP || thisg==A_GAMMAM)
			{
			int dist=0;
			List l2;
			l2=ListTail(l1);
			while(!is_empty_list(l2))
				{
				if(CompoundArg1(ListFirst(l2))==A_GAMMA)
						dist++;
				l2=ListTail(l2);
				}
			if((dist%2) && thisg==A_GAMMA5)
				retf*=-1;
			if(thisg!=A_GAMMA5 && (dist%2))
			{
				if(thisg==A_GAMMAP)
					thisg=A_GAMMAM;
				else
					thisg=A_GAMMAP;
			}
			if(dist)
				{
				SetCompoundArg(ListFirst(l1),1,thisg);
				gl=AppendLast(gl,CopyTerm(ListFirst(l1)));
				SetCompoundArg(ListFirst(l1),1,A_DELTA);
				}
			goto rmdlt;
			}
		l1=ListTail(l1);
		}

	
		
rmdlt:
	l1=gl;
	while(!is_empty_list(l1))
		{
		if(CompoundArg1(ListFirst(l1))==A_DELTA)
			{
			FreeAtomic(ListFirst(l1));
			ChangeList(l1,0);
			gl=CutFromList(gl,l1);
			goto rmdlt;
			}
		l1=ListTail(l1);
		}

	if(is_empty_list(gl))
		{
		if(lab1==lab2)
			{
			puts("Internal error (a2s08)");
			return 0;
			}
		lab2=CompoundArg2(ListFirst(ListTail(prt)));
		l1=*spec;
		while(!is_empty_list(l1))
			{
			Term y1;
			y1=ListFirst(l1);
			if((CompoundName(y1)==OPR_SPINOR || CompoundName(y1)==OPR_SPINOR3)
					&& CompoundArg1(y1)==lab2)
				{
				ChangeList(CompoundArg2(y1),lab1);
				goto tend;
				}
			l1=ListTail(l1);
			}
		if(is_empty_list(l1))
			{
			puts("Internal error (a2s07)");
			return 0;
			}
		}

	ChangeList(CompoundArg2(ListFirst(gl)),lab1);
	l1=gl;
	while(!is_empty_list(ListTail(l1)))
		{
		Label lab;
		lab=NewLabel();
		ChangeList(ListTail(CompoundArg2(ListFirst(l1))),lab);
		ChangeList(CompoundArg2(ListFirst(ListTail(l1))),lab);
		l1=ListTail(l1);
		}
	ChangeList(ListTail(CompoundArg2(ListFirst(l1))),lab2);

tend:

	*spec=ConcatList(*spec,gl);
/*
DumpList(*spec);printf("retf=%d\n",retf);
*/
		
checkcc:


		
	{
	int typ;
	Term pr1, pr2;
	int mfl=0;
	Atom bnm=0,cnm=0;
	int fl=0;
	Atom o1,o2;
		{

		pr1=GetAtomProperty(CompoundArg1(ListFirst(prt)),PROP_TYPE);
		pr2=GetAtomProperty(CompoundArg1(ListFirst(ListTail(prt))),PROP_TYPE);

		if(is_compound(pr1) && CompoundName(pr1)==OPR_FIELD &&
		   CompoundArg2(pr1)==NewInteger(4))
		   	pr1=CompoundArg1(pr1);
		else
			pr1=0;
		if(is_compound(pr2) && CompoundName(pr2)==OPR_FIELD &&
		   CompoundArg2(pr2)==NewInteger(4))
		   	pr2=CompoundArg1(pr2);
		else
			pr2=0;
		if((pr1==0 && pr2==0) || (pr1!=0 && pr2!=0))
			typ=0;
		else
			if(pr2==0)
				{
				typ=1;
				bnm=pr1;
				cnm=CallFunction(MakeCompound1(A_CC,
					CompoundArg1(ListFirst(ListTail(prt)))),0);
				if(bnm==CompoundArg1(ListFirst(ListTail(prt))))
					mfl=1;
				}
			else
				{
				typ=2;
				bnm=pr2;
				cnm=CallFunction(MakeCompound1(A_CC,
					CompoundArg1(ListFirst(prt))),0);
				if(bnm==CompoundArg1(ListFirst(prt)))
					mfl=1;
				}
		}
	if(pr1 && pr2)
		typ=3;

	if(typ==0)
		goto theend;

	if(mfl)
		{
		Term prp;
		int csf=1;
		int gsf=1;
		int half=0;
		
		

/*		prp=GetAtomProperty(bnm,PROP_TYPE);
		if(CompoundArg1(prp)!=CompoundArg2(prp))
			{
			printf("Warning: vertex with 2 equivalent Dirac spinors '");
			WriteTerm(bnm);
			puts("'");
			}*/
		prp=GetAtomProperty(bnm,A_COLOR);
		if(is_compound(prp) && CompoundArg1(prp)==NewInteger(3))
			{
			csf=0;
			if(!clab1 || !clab2)
				{
				puts("Internal error (a2s09)");
				return 0;
				}
			if(clab1==clab2)
				csf=1;
			else
				{
				l1=*spec;
				while(!is_empty_list(l1))
					{
					Term sp;
					sp=ListFirst(l1);
					if(is_compound(sp) && is_atom(CompoundArg1(sp))
						&& (A_COLOR_F==GetAtomProperty(CompoundArg1(sp),A_COLOR) ||
							A_COLOR_EPS==GetAtomProperty(CompoundArg1(sp),A_COLOR))
						&& ListMember(CompoundArg2(sp),clab1)
						&& ListMember(CompoundArg2(sp),clab2) )
							csf=-1;
					l1=ListTail(l1);
					}
				}
			}
			
		l1=*spec;
		fl=0;
		while(!is_empty_list(l1))
			{
			Term t;
			t=ListFirst(l1);
			if(CompoundName(t)==OPR_SPINOR || CompoundName(t)==OPR_SPINOR3)
				{
				fl=(!fl);
				goto ui;
				}
			if(CompoundName(t)==OPR_SPECIAL && CompoundArg1(t)==A_MOMENT)
				{
				if(fl)
					gsf*=-1;
				goto ui;
				}
			ui:
			l1=ListTail(l1);
			}

		if(ListLength(gl)==1)
			{
			if(CompoundArg1(ListFirst(gl))==A_GAMMA)
				gsf*=-1;
			}
			else
				if(ListLength(gl)==2)
					{
					if(CompoundArg1(ListFirst(ListTail(gl)))==A_GAMMA)
						gsf*=-1;
					else
						{
						Atom gpm=0;
						if(CompoundArg1(ListFirst(ListTail(gl)))!=A_GAMMA5)
							gpm=CompoundArg1(ListFirst(ListTail(gl)));
						if(gpm==0)
						{
							gsf*=1;
						}
						else
						if(gpm && gsf*csf==1)
							{
							half=1;
							if(gpm==A_GAMMAM)
								half=-1;
							SetCompoundArg(ListFirst(ListTail(gl)),1,A_GAMMA5);
							}
						else if(gpm && gsf*csf==-1)
							{
							Term lab1,lab2;
							gsf*=-1;
							half=1;
							lab1=ListFirst(CompoundArg2(ListFirst(ListTail(gl))));
							lab2=ListNth(CompoundArg2(ListFirst(ListTail(gl))),2);
							for(l1=(*spec);l1;l1=ListTail(l1))
								{
								Term t;
								t=ListFirst(l1);
								if((CompoundName(t)==OPR_SPINOR || CompoundName(t)==OPR_SPINOR3)
										&& ListFirst(CompoundArg2(t))==lab2)
									{
									ChangeList(CompoundArg2(t),lab1);
									break;
									}
								}
							if(!l1)
								puts("Internal error (a2s13)");
							CutFromList(gl,ListTail(gl));
							}
						
						}
					}
/*printf("%d %d %d\n",csf,gsf,ch_sign);*/
		if(csf*gsf==-1 && ch_sign==0)
			return 0;	
		return (half?half:2)*retf;
		}

	o1=CompoundArg1(ListFirst(prt));
	o2=CompoundArg1(ListFirst(ListTail(prt)));
	if(typ==3 || (typ==1 && strcmp(AtomValue(bnm),AtomValue(o2))>0)
		|| (typ==2 && strcmp(AtomValue(o1),AtomValue(bnm))>0) )
		{
		int gsf=0;
		Term p1,p2;

		if(is_empty_list(gl))
			gsf=1;
		else
			if(ListLength(gl)==1)
				{
				if(CompoundArg1(ListFirst(gl))==A_GAMMA)
					gsf=-1;
				else
					gsf=1;
				}
			else
				if(ListLength(gl)==2)
					{
					if(CompoundArg1(ListFirst(ListTail(gl)))==A_GAMMA)
						{
						gsf=1;
						inv_gamma(gl);
						}
					else
						if(CompoundArg1(ListFirst(ListTail(gl)))==A_GAMMA5)
							gsf=1;
						else
							{
							gsf=-1;
							if(CompoundArg1(ListFirst(ListTail(gl)))==A_GAMMAP)
								SetCompoundArg(ListFirst(ListTail(gl)),1,A_GAMMAM);
							else
								SetCompoundArg(ListFirst(ListTail(gl)),1,A_GAMMAP);
							}
					}
				else
					if(ListLength(gl)==3)
					{
					inv_gamma(gl);
					
					if(CompoundArg1(ListFirst(ListTail(ListTail(gl))))==A_GAMMA)
						gsf=-1;
					else
						gsf=1;
					}
						
					
				
		if(gsf==0)
			{
			puts("Internal error (a2s10), gamma matrices in cc vertex:");
			DumpList(gl);
			return retf;
			}
		retf*=gsf;
		if(ch_sign)
			retf*=-gsf;
		p1=ListFirst(prt);
		p2=ListFirst(ListTail(prt));
		ChangeList(prt,p2);
		ChangeList(ListTail(prt),p1);
		if(typ==1)
			{
			SetCompoundArg(p1,1,bnm);
			SetCompoundArg(p2,1,cnm);
			}
		else
			if(typ==2)
				{
				SetCompoundArg(p2,1,bnm);
				SetCompoundArg(p1,1,cnm);
				}
			else
				{
				SetCompoundArg(p1,1,pr1);
				SetCompoundArg(p2,1,pr2);
				}
		l1=gl;
		while(!is_empty_list(l1))
			{
			List l2;
			l2=CompoundArg2(ListFirst(l1));
			while(!is_empty_list(l2))
				{
				if(ListFirst(l2)==lab1)
					ChangeList(l2,lab2);
				else
					if(ListFirst(l2)==lab2)
						ChangeList(l2,lab1);
				l2=ListTail(l2);
				}
			l1=ListTail(l1);
			}
		}
	}

theend:
/*DumpList(*spec);	*/
	return retf;
	}

int alg2_updown(List prt, List *spec)
	{
	List l1;
	Label plab[2];
	int   ptyp[2];
	int   pbar[2];
	Atom  bname[2];
	Atom  baname[2];
	Term  pprt[2];
	Term  pspec[2];
	Label ilab[4];
	Label sind=0, sind1=0;
	
	int pno=0, sigma=0;

	
	for(l1=prt;l1;l1=ListTail(l1))
		{
		Term prop;
		prop=GetAtomProperty(CompoundArg1(ListFirst(l1)),PROP_TYPE);
		if(prop && is_compound(prop) && CompoundName(prop)==OPR_FIELD &&
				(CompoundArg2(prop)==NewInteger(6) ||
					CompoundArg2(prop)==NewInteger(7)))
			{
			if(pno>1)
				{
				pno++;
				continue;
				}
			ptyp[pno]=(int)IntegerValue(CompoundArg2(prop));
			bname[pno]=CompoundArg1(prop);
			plab[pno]=CompoundArg2(ListFirst(l1));
			pprt[pno]=ListFirst(l1);
			prop=GetAtomProperty(bname[pno],PROP_TYPE);
			if(bname[pno]==CompoundArg1(prop))
				{
				baname[pno]=CompoundArg2(prop);
				pbar[pno]=0;
				}
			else
				{
				baname[pno]=bname[pno];
				bname[pno]=CompoundArg1(prop);
				pbar[pno]=1;
				}
				
			pno++;
			}
		}
		

	if(pno && pno!=2)
		{
		ErrorInfo(918);
		printf(" %d 2-component fermions in vertex ",pno);
		WriteVertex(prt);
		puts("");
		return 0;
		}

	if(pno==0)
		return 1;

		
	pspec[0]=0;
	pspec[1]=0;
	
	for(l1=(*spec);l1;l1=ListTail(l1))
		{
		Term t;
		t=ListFirst(l1);
		if((CompoundName(t)==OPR_SPINOR || CompoundName(t)==OPR_SPINOR3) 
						&& CompoundArg1(t)==plab[0])
			pspec[0]=t;
		if((CompoundName(t)==OPR_SPINOR || CompoundName(t)==OPR_SPINOR3) 
						 && CompoundArg1(t)==plab[1])
			pspec[1]=t;
		}
	if(pspec[0]==0 || pspec[1]==0)
		{
		puts("Internal error (a2s18)");
		return 0;
		}

	ilab[0]=ListFirst(CompoundArg2(pspec[0]));
	ilab[1]=ListFirst(CompoundArg2(pspec[1]));

	if(ilab[0]!=ilab[1])
		{
		for(l1=*spec;l1;l1=ListTail(l1))
			{
			Term t;
			t=ListFirst(l1);
			if(CompoundName(t)==OPR_SPECIAL && ListLength(CompoundArg2(t))==3 && 
				ListFirst(CompoundArg2(t))==ilab[0] &&
				ListFirst(ListTail(CompoundArg2(t)))==ilab[1])
				{
				sind=ListFirst(ListTail(ListTail(CompoundArg2(t))));
				*spec=CutFromList(*spec,l1);
				sigma=1;
				break;
				}
			if(CompoundName(t)==OPR_SPECIAL && ListLength(CompoundArg2(t))==4 && 
				ListFirst(CompoundArg2(t))==ilab[0] &&
				ListFirst(ListTail(CompoundArg2(t)))==ilab[1])
				{
				sind=ListFirst(ListTail(ListTail(CompoundArg2(t))));
				sind1=ListFirst(ListTail(ListTail(ListTail(CompoundArg2(t)))));
				*spec=CutFromList(*spec,l1);
				sigma=2;
				break;
				}
			}
		if(sigma==0)
			{
			puts("Internal error (a2s19)");
			return 0;
			}
		}

	/*	printf("%s %s %s\n",AtomValue(CompoundArg1(pprt[0])),sigma?"sigma":"",
		AtomValue(CompoundArg1(pprt[1])));*/

	if(sigma==0)
		{
		Term newsp;
		
		if((ptyp[0]==6 && ptyp[1]==7) || (ptyp[0]==7 && ptyp[1]==6))
			{
			ErrorInfo(914);
			printf("wrong combination of 2-component spinors in vertex ");
			WriteVertex(prt);
			printf("\n");
			return 0;
			}
					
		ilab[0]=NewLabel();
		ilab[1]=NewLabel();
		newsp=MakeCompound(OPR_SPECIAL,2);
		SetCompoundArg(newsp,2,MakeList2(ilab[0],ilab[1]));
		if(ptyp[0]==6)
			SetCompoundArg(newsp,1,A_GAMMAM);
		else
			SetCompoundArg(newsp,1,A_GAMMAP);
		*spec=AppendLast(*spec,newsp);
		ChangeList(CompoundArg2(pspec[0]),ilab[0]);
		ChangeList(CompoundArg2(pspec[1]),ilab[1]);
		
		if(pbar[0]==1)
			SetCompoundArg(pprt[0],1,baname[0]);
		else
			SetCompoundArg(pprt[0],1,
				CallFunction(MakeCompound1(A_CC,bname[0]),0));
				
		if(pbar[1]==0)
			SetCompoundArg(pprt[1],1,bname[1]);
		else
			SetCompoundArg(pprt[1],1,
				CallFunction(MakeCompound1(A_CC,baname[1]),0));

		}
	else if(sigma==1)
		{
		Term newsp;
		
		if(ptyp[0]!=7 || ptyp[1]!=6)
			{
			ErrorInfo(915);
			printf("wrong combination of 2-component spinors in vertex ");
			WriteVertex(prt);
			puts("");
			return 0;
			}
/*		
		if(pbar[0]==0 && pbar[1]==0 &&
				bname[0]!=baname[0] && bname[1]==baname[1])
		{
			char bbuf[16];
			pbar[1]=1;
			sprintf(bbuf,"%s.c",AtomValue(bname[1]));
			bname[1]=NewAtom(bbuf,0);
		}
			
		if(pbar[0]==0 && pbar[1]==0 && 
				bname[0]==baname[0] && bname[1]==baname[1])
		{
			char bbuf[16];
			pbar[0]=1;
			sprintf(bbuf,"%s.c",AtomValue(bname[0]));
			baname[0]=NewAtom(bbuf,0);
		}
*/			
/*		if(pbar[0]+pbar[1]!=1 && bname[0]!=baname[0] && bname[1]!=baname[1])
			{
			ErrorInfo(916);
			printf("wrong combination of 2-component spinors in vertex ");
			WriteVertex(prt);
			puts("");
			return 0;
			}
*/			
		ilab[0]=NewLabel();
		ilab[1]=NewLabel();
		ilab[2]=NewLabel();
		
		newsp=MakeCompound(OPR_SPECIAL,2);
		SetCompoundArg(newsp,2,MakeList3(ilab[0],ilab[1],sind));
		SetCompoundArg(newsp,1,A_GAMMA);
		*spec=AppendLast(*spec,newsp);
		
		newsp=MakeCompound(OPR_SPECIAL,2);
		SetCompoundArg(newsp,2,MakeList2(ilab[1],ilab[2]));
		SetCompoundArg(newsp,1,A_GAMMAM);
		*spec=AppendLast(*spec,newsp);
		
		ChangeList(CompoundArg2(pspec[0]),ilab[0]);
		ChangeList(CompoundArg2(pspec[1]),ilab[2]);
		
		if(pbar[0]==1)
			SetCompoundArg(pprt[0],1,baname[0]);
		else
			SetCompoundArg(pprt[0],1,
				CallFunction(MakeCompound1(A_CC,bname[0]),0));
				
		if(pbar[1]==0)
			SetCompoundArg(pprt[1],1,bname[1]);
		else
			SetCompoundArg(pprt[1],1,
				CallFunction(MakeCompound1(A_CC,baname[1]),0));

		}
	else
		{
		Term newsp;
		
		if(ptyp[0]!=6 || ptyp[1]!=6)
			{
			ErrorInfo(915);
			printf("wrong combination of 2-component spinors in vertex ");
			WriteVertex(prt);
			puts("");
			return 0;
			}
/*		
		if(pbar[0]==0 && pbar[1]==0 &&
				bname[0]!=baname[0] && bname[1]==baname[1])
		{
			char bbuf[16];
			pbar[1]=1;
			sprintf(bbuf,"%s.c",AtomValue(bname[1]));
			bname[1]=NewAtom(bbuf,0);
		}
			
		if(pbar[0]==0 && pbar[1]==0 && 
				bname[0]==baname[0] && bname[1]==baname[1])
		{
			char bbuf[16];
			pbar[0]=1;
			sprintf(bbuf,"%s.c",AtomValue(bname[0]));
			baname[0]=NewAtom(bbuf,0);
		}
*/			
/*		if(pbar[0]+pbar[1]!=1 && bname[0]!=baname[0] && bname[1]!=baname[1])
			{
			ErrorInfo(916);
			printf("wrong combination of 2-component spinors in vertex ");
			WriteVertex(prt);
			puts("");
			return 0;
			}
*/			
		ilab[0]=NewLabel();
		ilab[1]=NewLabel();
		ilab[2]=NewLabel();
		ilab[3]=NewLabel();
		
		newsp=MakeCompound(OPR_SPECIAL,2);
		SetCompoundArg(newsp,2,MakeList3(ilab[0],ilab[1],sind));
		SetCompoundArg(newsp,1,A_GAMMA);
		*spec=AppendLast(*spec,newsp);

		newsp=MakeCompound(OPR_SPECIAL,2);
		SetCompoundArg(newsp,2,MakeList3(ilab[1],ilab[2],sind1));
		SetCompoundArg(newsp,1,A_GAMMA);
		*spec=AppendLast(*spec,newsp);
		
		newsp=MakeCompound(OPR_SPECIAL,2);
		SetCompoundArg(newsp,2,MakeList2(ilab[2],ilab[3]));
		SetCompoundArg(newsp,1,A_GAMMAM);
		*spec=AppendLast(*spec,newsp);
		
		ChangeList(CompoundArg2(pspec[0]),ilab[0]);
		ChangeList(CompoundArg2(pspec[1]),ilab[3]);
		
		if(pbar[0]==1)
			SetCompoundArg(pprt[0],1,baname[0]);
		else
			SetCompoundArg(pprt[0],1,
				CallFunction(MakeCompound1(A_CC,bname[0]),0));
				
		if(pbar[1]==0)
			SetCompoundArg(pprt[1],1,bname[1]);
		else
			SetCompoundArg(pprt[1],1,
				CallFunction(MakeCompound1(A_CC,baname[1]),0));

		}
				

	return 1;
	}

static void inv_gamma(List gamma)
{

	List l1,l2;
	List gl=0, gli=0;
	
	for(l1=gamma;l1;l1=ListTail(l1))
	{
		Term t;
		t=ListFirst(l1);
		if(CompoundArg1(t)==A_GAMMA)
		{
			gl=AppendLast(gl,t);
			gli=AppendFirst(gli,ListNth(CompoundArg2(t),3));
		}
	}

	for(l1=gl,l2=gli;l1;l1=ListTail(l1),l2=ListTail(l2))
	{
		Term t;
		t=CompoundArg2(ListFirst(l1));
		t=ListTail(ListTail(t));
		ChangeList(t,ListFirst(l2));
	}

	RemoveList(gl);
	RemoveList(gli);
}
