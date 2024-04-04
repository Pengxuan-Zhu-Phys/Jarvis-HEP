#include <string.h>
#include "lanhep.h"

int max_prt_lenp=3, max_prt_lenl=5;

static char mk_g_buf[32];

static List all_particles = 0;


List all_prtc_list(void)
	{
	return all_particles;
	}

static void trigg1_color(Atom name, List ind)
	{
	int i=1;
	List l=ind;
	int r=0;
	while(!is_empty_list(l))
		{
		r+=color_repres(ListFirst(l));
		if(r<0 && !CalcOutput)
			{
			printf("Warning: particle %s breaks CompHEP color conventions.\n",
				AtomValue(name));
			return;
			}
		if(r!=0)
			{
			SetAtomProperty(name,A_COLOR,MakeCompound2(OPR_DIV,
				NewInteger(r),NewInteger(i)));
			r=-100;
			}
		i++;
		l=ListTail(l);
		}
	}


static void run_trigg_1(Atom name, List ind)
	{
	int pl=strlen(AtomValue(name));
	if(pl>max_prt_lenl)
		max_prt_lenl=pl;
		
	all_particles=AppendLast(all_particles,name);
	trigg1_color(name, ind);
	}


static void trigg_gsb(Term prt, List ind, int anti)
	{
	Atom gh,agh,a,ap;
	if(CompoundArgN(prt,7)!=A_GAUGE || 
			CompoundArgN(prt,4)!=NewInteger(2)||
			CompoundArgN(prt,5)==0)
		return;
		
	ind=CopyTerm(ListTail(ind));
	if(anti)
		a=CompoundArg2(prt);
	else
		a=CompoundArg1(prt);
	if(anti)
		ap=CompoundArg1(prt);
	else
		ap=CompoundArg2(prt);
		
	sprintf(mk_g_buf,"%s.f",AtomValue(a));
	gh=NewAtom(mk_g_buf,0);
	sprintf(mk_g_buf,"%s.f",AtomValue(ap));
	agh=NewAtom(mk_g_buf,0);
	SetAtomProperty(gh,A_ANTI,agh);
	SetAtomProperty(gh,PROP_TYPE,MakeCompound2(OPR_FIELD,a,NewInteger(1)));
	if(!is_empty_list(ind))
		SetAtomProperty(gh,PROP_INDEX,ind);	
	run_trigg_1(gh,ind);
	}

/*
static void trigg_b(Term prt, List ind, int anti)
	{
	Atom gh,agh,a,ap;
	if(CompoundArgN(prt,7)!=A_GAUGE || 
			CompoundArgN(prt,4)!=NewInteger(2))
		return;
		
	ind=CopyTerm(ListTail(ind));
	if(anti)
		a=CompoundArg2(prt);
	else
		a=CompoundArg1(prt);
	if(anti)
		ap=CompoundArg1(prt);
	else
		ap=CompoundArg2(prt);
		
	sprintf(mk_g_buf,"%s.b",AtomValue(a));
	gh=NewAtom(mk_g_buf,0);
	sprintf(mk_g_buf,"%s.b",AtomValue(ap));
	agh=NewAtom(mk_g_buf,0);
	SetAtomProperty(gh,A_ANTI,agh);
	SetAtomProperty(gh,PROP_TYPE,MakeCompound2(OPR_FIELD,a,NewInteger(8)));
	if(!is_empty_list(ind))
		SetAtomProperty(gh,PROP_INDEX,ind);	
	run_trigg_1(gh,ind);
	}
*/


static void trigg_ghost(Term prt, List ind, int anti)
	{
	int i;
	Atom gh,agh,cgh,acgh,a,aa;
	
	if(CompoundArgN(prt,7)!=A_GAUGE) return;
	
	ind=CopyTerm(ListTail(ind));
	
	if(anti)
	{
		a=CompoundArg2(prt);
		aa=CompoundArg1(prt);
	}
	else
	{
		a=CompoundArg1(prt);
		aa=CompoundArg2(prt);
	}
	
	i=sprintf(mk_g_buf,"%s",AtomValue(a));
	
	sprintf(mk_g_buf+i,".c");
	gh=NewAtom(mk_g_buf,0);
	sprintf(mk_g_buf+i,".C");
	cgh=NewAtom(mk_g_buf,0);
	
	i=sprintf(mk_g_buf,"%s",AtomValue(aa));
	
	sprintf(mk_g_buf+i,".c");
	agh=NewAtom(mk_g_buf,0);
	sprintf(mk_g_buf+i,".C");
	acgh=NewAtom(mk_g_buf,0);
	
	SetAtomProperty(a,A_GHOST,MakeCompound2(A_GHOST,gh,cgh));
	SetAtomProperty(gh,A_GHOST,a);
	SetAtomProperty(gh,A_GRASS,NewInteger(1));
	SetAtomProperty(cgh,A_GHOST,a);
	SetAtomProperty(cgh,A_GRASS,NewInteger(1));
	
	SetAtomProperty(gh,A_ANTI,agh);
	SetAtomProperty(cgh,A_ANTI,acgh);
	
	SetAtomProperty(gh,A_ANTI2,cgh);
	SetAtomProperty(cgh,A_ANTI2,gh);
	
	SetAtomProperty(gh,PROP_TYPE,MakeCompound2(OPR_FIELD,a,NewInteger(2)));
	SetAtomProperty(cgh,PROP_TYPE,MakeCompound2(OPR_FIELD,a,NewInteger(3)));
	if(!is_empty_list(ind))
		{
		SetAtomProperty(gh,PROP_INDEX,ind);
		ind=CopyTerm(ind);
		SetAtomProperty(cgh,PROP_INDEX,ind);
		};
	run_trigg_1(gh,ind);
	run_trigg_1(cgh,ind);
	}

		
		
static void trigg_cc(Term prt, List ind, int anti)
	{
	List ii;
    Atom a, ca, aa, caa;
    Term pf;
	if(CompoundArgN(prt,4)!=NewInteger(1) &&
			CompoundArgN(prt,4)!=NewInteger(3))
		return;
    
    pf=GetAtomProperty(NewAtom("option",0),NewAtom("ccTranspAll",0));
    if(pf==NewInteger(0))
    	pf=0;
    	
    if(anti)
		a=CompoundArg2(prt);
	else
		a=CompoundArg1(prt);
	if(anti)
		aa=CompoundArg1(prt);
	else
		aa=CompoundArg2(prt);
        	
    sprintf(mk_g_buf,"%s.c",AtomValue(a));
    ca=NewAtom(mk_g_buf,0);
    sprintf(mk_g_buf,"%s.c",AtomValue(aa));
    caa=NewAtom(mk_g_buf,0);
    ii=CopyTerm(ind);
    if(pf)
    	{
    	List l1;
    	l1=ii;
    	while(!is_empty_list(l1))
    		{
    		Term tt, in1,in2;
			tt=CompoundArg1(ListFirst(l1));
			in1=ConsumeCompoundArg(tt,1);
			in2=ConsumeCompoundArg(tt,2);
			SetCompoundArg(tt,1,in2);
			SetCompoundArg(tt,2,in1);
			l1=ListTail(l1);
			}
		}
	else
    	if(!is_empty_list(ii) && 
        	CompoundName(CompoundArg1(ListFirst(ii)))==A_LORENTZ)
				{
				Term tt, in1,in2;
				tt=CompoundArg1(ListFirst(ii));
				in1=ConsumeCompoundArg(tt,1);
				in2=ConsumeCompoundArg(tt,2);
				SetCompoundArg(tt,1,in2);
				SetCompoundArg(tt,2,in1);
				}
	SetAtomProperty(ca,A_ANTI,caa);
	SetAtomProperty(ca,PROP_INDEX,ii);
	SetAtomProperty(ca,PROP_TYPE,MakeCompound2(OPR_FIELD,a,NewInteger(4)));
	SetAtomProperty(ca,A_CHNAME,a);
	SetAtomProperty(ca,A_GRASS,NewInteger(1));
	
	run_trigg_1(ca,ii);
    }
    
static void trigg_tens(Term prt, List ind, int anti)
	{
	List ii;
    Atom a, ca, aa, caa;
    
	if(CompoundArgN(prt,4)!=NewInteger(2))
		return;
    
    if(anti)
		a=CompoundArg2(prt);
	else
		a=CompoundArg1(prt);
	if(anti)
		aa=CompoundArg1(prt);
	else
		aa=CompoundArg2(prt);
    
    sprintf(mk_g_buf,"%s.t",AtomValue(a)); ca=NewAtom(mk_g_buf,0);
    sprintf(mk_g_buf,"%s.t",AtomValue(aa)); caa=NewAtom(mk_g_buf,0);
    ii=CopyTerm(ind);

    ii=AppendFirst(ii,MakeCompound2(A_I,SpecToRepr(OPR_VECTOR),0));
    SetAtomProperty(ca,A_ANTI,caa);
    SetAtomProperty(ca,PROP_INDEX,ii);
    SetAtomProperty(ca,PROP_TYPE,MakeCompound2(OPR_FIELD,a,NewInteger(5)));
    run_trigg_1(ca,ii);
    }
    
static void trigg_updown(Term prt, List ind, int anti)
	{
	List ii;
    Atom a, ca, ca1, aa, caa, caa1;

	if(CompoundArgN(prt,4)!=NewInteger(1) 
				&& CompoundArgN(prt,4)!=NewInteger(3))
		return;
    	
    if(anti)
		a=CompoundArg2(prt);
	else
		a=CompoundArg1(prt);
		
	if(anti)
		aa=CompoundArg1(prt);
	else
		aa=CompoundArg2(prt);
        	
    sprintf(mk_g_buf,"%s.u",AtomValue(a));
    ca=NewAtom(mk_g_buf,0);
    sprintf(mk_g_buf,"%s.u",AtomValue(aa));
    caa=NewAtom(mk_g_buf,0);
    ii=CopyTerm(ind);
    
    if(!is_empty_list(ii) && 
        CompoundName(CompoundArg1(ListFirst(ii)))==A_LORENTZ)
				{
				Term tt;
				tt=CompoundArg1(ListFirst(ii));
				SetCompoundArg(tt,1,NewInteger(0));
				SetCompoundArg(tt,2,NewInteger(0));
				}
	else
		puts("Internal error (rprt mkud)");
		
	
	sprintf(mk_g_buf,"%s.d",AtomValue(a));
    ca1=NewAtom(mk_g_buf,0);
    sprintf(mk_g_buf,"%s.d",AtomValue(aa));
    caa1=NewAtom(mk_g_buf,0);
	
	SetAtomProperty(ca,A_ANTI,caa1);
	SetAtomProperty(ca,PROP_INDEX,CopyTerm(ii));
	SetAtomProperty(ca,PROP_TYPE,MakeCompound2(OPR_FIELD,a,NewInteger(6)));
	SetAtomProperty(ca,A_GRASS,NewInteger(1));
    
	SetAtomProperty(ca1,A_ANTI,caa);
	SetAtomProperty(ca1,PROP_INDEX,ii);
	SetAtomProperty(ca1,PROP_TYPE,MakeCompound2(OPR_FIELD,a,NewInteger(7)));
	SetAtomProperty(ca1,A_GRASS,NewInteger(1));
	
	run_trigg_1(ca,ii);
	run_trigg_1(ca1,ii);
    }
 
	

static void run_triggers(Term prt, List ind, int anti)
	{

	int pl=strlen(AtomValue(CompoundArg1(prt)));
	if(pl>max_prt_lenp)
		max_prt_lenp=pl;
		
	if(anti)
		run_trigg_1(CompoundArg2(prt),ind);
	else
		run_trigg_1(CompoundArg1(prt),ind);

	trigg_gsb(prt,ind,anti);
	trigg_ghost(prt,ind,anti);
	trigg_cc(prt,ind,anti);
	trigg_tens(prt,ind,anti);
	trigg_updown(prt,ind,anti);
/*	trigg_b(prt,ind,anti);*/
	
	}


	

void register_particle(Term prt)
    {
    List l,li,li1,prt1;
	li=NewList();
	
	if(IntegerValue(CompoundArgN(prt,4))==1)
		li=AppendLast(li,MakeCompound2(A_I,SpecToRepr(OPR_SPINOR),0));
	if(IntegerValue(CompoundArgN(prt,4))==2)
		li=AppendLast(li,MakeCompound2(A_I,SpecToRepr(OPR_VECTOR),0));
	if(IntegerValue(CompoundArgN(prt,4))==3)
		{
		li=AppendLast(li,MakeCompound2(A_I,SpecToRepr(OPR_SPINOR),0));
		li=AppendLast(li,MakeCompound2(A_I,SpecToRepr(OPR_VECTOR),0));
		}
	if(IntegerValue(CompoundArgN(prt,4))==4)
		{
		li=AppendLast(li,MakeCompound2(A_I,SpecToRepr(OPR_VECTOR),0));
		li=AppendLast(li,MakeCompound2(A_I,SpecToRepr(OPR_VECTOR),0));
		}
		
	l=CompoundArgN(prt,8);
	while(!is_empty_list(l))
		{
		Term t1,t2;
		t1=ListFirst(l);
		if(CompoundArity(t1)==2)
			{
			t2=MakeCompound2(A_I,CopyTerm(t1),0);
			li=AppendLast(li,t2);
			}
		l=ListTail(l);
		}
		
	ReportRedefined(CompoundArg1(prt),"particle");
	ReportRedefined(CompoundArg2(prt),"particle");
		
    if(CompoundArg1(prt)!=CompoundArg2(prt))
        {
        prt1=CopyTerm(prt);
        li1=l=ConsumeCompoundArg(prt1,8);
        while(!is_empty_list(l))
        	{
        	Term t;
        	t=ListFirst(l);
        	if(CompoundArity(t)==2)
        		{
        		Term t1,t2;
        		t1=ConsumeCompoundArg(t,1);
        		t2=ConsumeCompoundArg(t,2);
        		SetCompoundArg(t,1,t2);
        		SetCompoundArg(t,2,t1);
        		}
        	else
        		{
        		Term t1;
        		t1=CompoundArg1(CompoundArg1(t));
        		t1=NewInteger(-IntegerValue(t1));
        		SetCompoundArg(CompoundArg1(t),1,t1);
        		}
        	l=ListTail(l);
        	}
        SetCompoundArg(prt1,8,li1);
        
        li1=l=CopyTerm(li);
        while(!is_empty_list(l))
        	{
        	Term t,t1,t2;
        	t=ListFirst(l);
        	t1=ConsumeCompoundArg(CompoundArg1(t),1);
        	t2=ConsumeCompoundArg(CompoundArg1(t),2);
        	SetCompoundArg(CompoundArg1(t),1,t2);
        	SetCompoundArg(CompoundArg1(t),2,t1);
        	l=ListTail(l);
        	}
        
        run_triggers(prt,li1,1);
        	
        SetAtomProperty(CompoundArg2(prt),PROP_TYPE,prt1);
		SetAtomProperty(CompoundArg2(prt),A_ANTI,CompoundArg1(prt));
        if(CompoundArgN(prt,4)==NewInteger(1))
        	SetAtomProperty(CompoundArg2(prt),A_GRASS,NewInteger(1));
        if(!is_empty_list(li1))
            SetAtomProperty(CompoundArg2(prt),PROP_INDEX,li1);
        
       }
       
	run_triggers(prt,li,0);
	        
    SetAtomProperty(CompoundArg1(prt),PROP_TYPE,prt);
    if(CompoundArgN(prt,4)==NewInteger(1))
        SetAtomProperty(CompoundArg1(prt),A_GRASS,NewInteger(1));
    if(!is_empty_list(li))
        SetAtomProperty(CompoundArg1(prt),PROP_INDEX,li);
    SetAtomProperty(CompoundArg1(prt),A_ANTI,CompoundArg2(prt));
    
    return;
    }
    
 
