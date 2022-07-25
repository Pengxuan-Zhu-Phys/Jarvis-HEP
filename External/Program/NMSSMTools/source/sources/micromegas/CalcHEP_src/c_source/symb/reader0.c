/*
 Copyright (C) 1997,2006, Alexander Pukhov 
*/

#include "physics.h"
#include "syst.h"
#include "syst2.h"
#include "reader0.h"

#include "prepdiag.h"
#include "parser.h"
#include "writeF.h"
#include "crt.h"
#include "process.h"


static char momsubst[4], indsubst[4];
static int  r_reading0 = 0;
static int maxMom=1000;
static int forR_code;

static void * bact0(char ch,void* mm1,void* mm2)
{
 char    *m1, *m2, *ans;
 int sgn;
 
  if ( r_reading0 && (ch == '*' ))
  {
     m1 = (char *) mm2; 
     m2 = (char *) mm1;
  }
  else
  {
     m1 = (char *) mm1;
     m2 = (char *) mm2;
  }
  if (ch == '+' || ch=='-')
  {
     lShift(m1,2);
     lShift(m2,2);
  }
  else
  {
     if (m1[0] == 'P' || ch == '^' )
     {  lShift(m1,1);
        m1[0]='(';
        strcat(m1,")");
     }  else lShift(m1,2);
     if (m2[0] == 'P' || ch== '^')
     {  lShift(m2,1);
        m2[0]='(';
        strcat(m2,")");
     } else lShift(m2,2);
  }

  ans=(char *) m_alloc(strlen(m1)+strlen(m2)+8);
  switch (ch)
  {
      case '+':
         if (m2[0] == '-') sprintf(ans,"P|%s%s",m1,m2);
         else              sprintf(ans,"P|%s+%s",m1,m2);
         break;

      case '-':
         if (m2[0] == '-') sprintf(ans,"P|%s+%s",m1,m2+1);
         else              sprintf(ans,"P|%s-%s",m1,m2);
         break;

      case '*':
	 sgn=1;
         if (m1[0]=='-') { lShift(m1,1); sgn= -sgn;}
         if (m2[0]=='-') { lShift(m2,1); sgn= -sgn;}
         if (sgn==1) sprintf(ans,"M|%s*%s",m1,m2);
	 else        sprintf(ans,"M|-%s*%s",m1,m2);
	 break;

      case '.':
               if (m2[0] != '-') sprintf(ans,"M|%s.%s",m1,m2);
         else  if (m1[0] != '-') sprintf(ans,"M|%s.%s",m2,m1);
         else
         {
            lShift(m1,1);
            lShift(m2,1);
            sprintf(ans,"M|%s.%s",m1,m2);
         }
         break;

      case '^':    sprintf(ans,"M|%s^%s",m1,m2);
  }  
  return (void*) ans;
}


static void *  uact0(char* ch,void* mm)
{char  *m, *ans;

   m = (char *) mm;
   ans=(char *) m_alloc(strlen(m)+10);

   if (strcmp(ch,"-") == 0)
   {  if (m[0] == 'M')
     {
        if (m[2] == '-') sprintf(ans,"M|%s",m+3);
        else  sprintf(ans,"M|-%s",m+2);
     }  else  sprintf(ans,"M|-(%s)",m+2);
   }
   
   if (strcmp(ch,"G") == 0)
   {  if(forR_code)
      {
        if (r_reading0) sprintf(ans,"M|-G(ln,%s)",m+2);
        else            sprintf(ans,"M|G(ln,%s)",m+2);
      }else
      {
        if (r_reading0) sprintf(ans,"M|-G(%s)",m+2);
        else            sprintf(ans,"M|G(%s)",m+2);
      }      
   }   
   return (void*) ans;
}

static void * act_rcode(char * ch, int n, void**args)
{  if(n==1) return uact0(ch,args[0]);
   if(n==2) return bact0(ch[0],args[0],args[1]);
   if(n==4 && !strcmp(ch,"eps"))
   { int l=15+strlen(args[0])+strlen(args[1])+strlen(args[2])+strlen(args[3]);
     char * ans=(char *) m_alloc(l);     
     sprintf(ans,"M|eps(%s,%s,%s,%s)",(char*)args[0]+2,(char*)args[1]+2,
     (char*)args[2]+2,(char*)args[3]+2);
     return ans;
   }
   return NULL;
}


static void*  rd_rcode(char* s)
{  char    * p;
   int  num;
   p = (char *) m_alloc(4+VAR_NAME_SIZE);
   p[0]=0;
   if (strlen(s) == 2 && s[1] > '0' && s[1] <= '9')
   {
      switch (s[0])
      {
        case 'p':
        case 'P':
               num = s[1] - '0';
               num = momsubst[num-1];
               if(num==maxMom) // for  parseOneVertex
               { p=realloc(p,3*maxMom+10);
                 strcpy(p,"M|(");
                 for(int i=1;i<maxMom;i++) sprintf(p+strlen(p),"-p%d",i);
                 strcat(p,")");
               } else 
               {  
                 if(num > 0) sprintf(p,"M|p%d",num); else    sprintf(p,"M|-p%d",-num);
               }  
            break;

        case 'm':
               num = s[1] - '0';
               num = indsubst[num-1];
	       sprintf(p,"M|m%d",num);
               break;           
        case 'M':
               num = s[1] - '0';
               num = indsubst[num-1]-1;
	       sprintf(p,"M|m%d",num);
      } 
      if (strcmp(s,"G5") == 0) 
      { if(forR_code) strcpy(p,"M|G(ln,A)"); else strcpy(p,"M|G5");}
   }
   if (!strlen(p)) sprintf(p,"M|%s",s);   
   return (void*) p;
}


char* parseVertex(int v, int forReduce)
{ 
  int i,j;
  int sgn=1;
  char *pstr, *pstr2;
  maxMom=1000; // To avoid substitutions P_maxMom=(-p1-p2 ...)
  for (j = 0; j < MAX(1,vcs.valence[v]); j++)
  {
      momsubst[j]=vcs.vertlist[v][vertexes[v].subst[j]-1].moment;
      indsubst[j]=vcs.vertlist[v][vertexes[v].subst[j]-1].lorentz;
  }
  r_reading0=vertexes[v].r_vert;

  forR_code=forReduce;                                        
  pstr = (char *) readExpression(vertexes[v].lgrnptr->description,
                                  rd_rcode, act_rcode, free);
  if(rderrcode)
  {  if(forReduce) outFileClose();
     finish();  
     sortie(60);
  }
  
  int nc=0;
  int cdi[MAXVALENCE],si[MAXVALENCE],rsubst[MAXVALENCE];
  for(i=0;i<vcs.valence[v];i++) rsubst[vertexes[v].subst[i]-1]=i+1;
  for(i=0;i<vcs.valence[v];i++)if(prtclbase1[vcs.vertlist[v][i].partcl].cdim!=1)
  {  si[nc]=rsubst[i];
     cdi[nc++]=prtclbase1[vcs.vertlist[v][i].partcl].cdim;
  }

  if(nc==3 && cdi[0]==cdi[1] && cdi[0]==cdi[2] && (abs(cdi[0])==3 || cdi[0]==8))
  {
         if (si[0] > si[1]) sgn*= -1;
         if (si[1] > si[2]) sgn*= -1;
         if (si[0] > si[2]) sgn*= -1;
  }
   
  if(sgn==1) 
  {  pstr2=m_alloc(strlen(pstr));
     strcpy(pstr2,pstr+2);
  } else
  { pstr2=m_alloc(strlen(pstr)+8);   
    sprintf(pstr2,"(-1)*(%s)",pstr+2);
  }                                                                                                                       
  free(pstr);
  return pstr2;                                                                                                                         
}

char* parseOneVertex(char*v, int * perm ,int * field)
{ 
  int i,j;
  int sgn=1;
  int val= field[3]==0? 3:4;
  char *pstr, *pstr2;

  maxMom=val;
  
  for(i = 0; i < val; i++)
  {
      momsubst[i]=perm[i];
      indsubst[i]=perm[i];
  }

  int spin[2]={0,0};
  for(i=0,j=0;i<val;i++) if(prtclbase1[field[i]].spin&1) spin[j++]=i;
  
  if(!spin[0]) r_reading0=0;
  else r_reading0= ( perm[spin[0]]> perm[spin[1]]);
  forR_code=0;
                                          
  pstr = (char *) readExpression(v,rd_rcode, act_rcode, free);
  if(rderrcode){ finish(); sortie(60);}
  
  int nc=0;
  int cdi[MAXVALENCE],si[MAXVALENCE],rsubst[MAXVALENCE];
  for(i=0;i<val;i++) rsubst[perm[i]-1]=i+1;
  for(i=0;i<val;i++)if(prtclbase1[field[i]].cdim!=1)
  {  si[nc]=rsubst[i];
     cdi[nc++]=prtclbase1[field[i]].cdim;
  }

  if(nc==3 && cdi[0]==cdi[1] && cdi[0]==cdi[2] && (abs(cdi[0])==3 || cdi[0]==8))
  {
         if (si[0] > si[1]) sgn*= -1;
         if (si[1] > si[2]) sgn*= -1;
         if (si[0] > si[2]) sgn*= -1;
  }
   
  if(sgn==1) 
  {  pstr2=m_alloc(strlen(pstr));
     strcpy(pstr2,pstr+2);
  } else
  { pstr2=m_alloc(strlen(pstr)+8);   
    sprintf(pstr2,"(-1)*(%s)",pstr+2);
  }                                                                                                                       
  free(pstr);
  return pstr2;                                                                                                                         
}


char * fermPropagTxt(int v,int l,int forReduce)
{  char * proptxt=(char*)malloc(100);
   int P=vcs.vertlist[v][l].partcl;
   int aux=prtclbase1[P].hlp;
   char * mass=prtclbase1[P].massidnt;
   char sgn[4]="+";
   int pa=vcs.vertlist[v][l].moment;
   char Ln[4]="";
   char G5[8]="G5";
   char dtwo[4]="";
      
   if(forReduce){strcpy(Ln,"ln,");strcpy(G5,"G(ln,A)");strcpy(dtwo,"/2");}
   
   if(pa<0){ pa=-pa; strcpy(sgn,"-");}
   int pb =  (pa==1 ? 2 : 1);
   if(aux=='*') strcpy(proptxt,mass);

   else if(strcmp(mass,"0"))
   {  sprintf(proptxt,"(%s%sG(%sp%d))",mass,sgn,Ln,pa);
      if(PLR_PRTCL&vcs.vertlist[v][l].prop)
      sprintf(proptxt+strlen(proptxt),"*(1-2*Helicity%d*%s*(N_pol_%d%d_*G(%sp%d)+N_pol_%d%d_*G(%sp%d)))",pa,G5,pa,pa,Ln,pa,pa,pb,Ln,pb);
//                                                   pa,G5,        pa,pa, Ln,pa,       pa,pb,  Ln,pb 
   }   
   else  if(!strchr("LR",aux) && !(PLR_PRTCL&vcs.vertlist[v][l].prop)  ) 
              sprintf(proptxt,"G(%s%sp%d)",Ln,sgn,pa);
   else if((fermionp(P)&&aux=='L')||(a_fermionp(P)&&aux=='R')) 
              sprintf(proptxt,"%sG(%sp%d)*(1-%s)%s",sgn,Ln,pa,G5,dtwo);
   else if((fermionp(P)&&aux=='R')||(a_fermionp(P)&&aux=='L'))                           
              sprintf(proptxt,"%sG(%sp%d)*(1+%s)%s",sgn,Ln,pa,G5,dtwo);
   else     sprintf(proptxt,"%sG(%sp%d)*(1%s2*Helicity%d*%s)",sgn,Ln,pa,sgn,pa,G5);
   return proptxt;
}   



char * tPropagator6mass4(int v,int l,set*indexs)
{ int vv=vcs.vertlist[v][l].link.vno;
  int ll=vcs.vertlist[v][l].link.edno;
  int m2=vcs.vertlist[v][l].lorentz,   m1=m2-1;
  int n2=vcs.vertlist[vv][ll].lorentz, n1=n2-1;
  int P=abs(vcs.vertlist[v][l].moment);
  char *mass=prtclbase1[vcs.vertlist[v][l].partcl].massidnt;
  char*txt=(char*)malloc(250); 

  if(indexs) *indexs=set_constr(m1,m2,n1,n2,_E);
  
  sprintf(txt,
     "3*(m%d.m%d*m%d.m%d+m%d.m%d*m%d.m%d-m%d.m%d*m%d.m%d)*%s^4",
     m1,n1, m2,n2,  m1,n2,m2,n1,  m1,m2,n1,n2, mass);

  sprintf(txt+strlen(txt),
     "-3*(m%d.m%d*p%d.m%d*p%d.m%d+m%d.m%d*p%d.m%d*p%d.m%d+"
     "m%d.m%d*p%d.m%d*p%d.m%d+m%d.m%d*p%d.m%d*p%d.m%d )*%s^2",
     m1,n1,P,m2,P,n2, m2,n2,P,m1,P,n1, m1,n2,P,m2,P,n1, m2,n1,P,m1,P,n2,mass);

  sprintf(txt+strlen(txt),           
     "+(%s^2*m%d.m%d+2*p%d.m%d*p%d.m%d)*(%s^2*m%d.m%d+2*p%d.m%d*p%d.m%d)",
     mass,m1,m2, P,m1,P,m2,mass, n1,n2, P,n1,P,n2);

  return txt; 
}

char * spin3_2_propagator(int v,int l,int forReduce)
{
  char*txt=(char*)malloc(250);
  char *Ln="",*sgn="";
  int vv=vcs.vertlist[v][l].link.vno;
  int ll=vcs.vertlist[v][l].link.edno;
  int m1=vcs.vertlist[v][l].lorentz;
  int m2=vcs.vertlist[vv][ll].lorentz;
  int P=vcs.vertlist[v][l].moment;
  char * mass=prtclbase1[vcs.vertlist[v][l].partcl].massidnt;

  if(P<0) {P=-P; sgn="-";} 
  if(forReduce) Ln="ln,";
  sprintf(txt, "-3*(%sG(%sp%d)+%s)"
               "*(%s^2*m%d.m%d - p%d.m%d*p%d.m%d)"
               "-(G(%sm%d)*%s +(%sp%d.m%d))"
               "*(%sG(%sp%d)-%s)"
               "*(G(%sm%d)*%s+(%sp%d.m%d))"

        ,sgn,Ln,P,mass
        , mass,m1,m2,P,m1,P,m2
       ,  Ln,m1,mass,sgn,P,m1
     ,    sgn,Ln,P,mass,
         Ln,m2,mass,sgn,P,m2
         
);   

return txt;
}
#ifdef XXX


/*4=============================================================*/
static symb_data tPropagator(int v,int l)
{ char txt[250];

  int vv=vcs.vertlist[v][l].link.vno;
  int ll=vcs.vertlist[v][l].link.edno;
  int m2=vcs.vertlist[v][l].lorentz,   m1=m2-1;
  int n2=vcs.vertlist[vv][ll].lorentz, n1=n2-1;
  int p=abs(vcs.vertlist[v][l].moment);
  char *mass=prtclbase1[vcs.vertlist[v][l].partcl].massidnt;

   
  indsubst[0]=m1;
  indsubst[1]=m2;
  indsubst[2]=n1;
  indsubst[3]=n2;
  momsubst[0]=p;

sprintf(txt,"(m1.m3*m2.m4+m1.m4*m2.m3-m1.m2*m3.m4)*3*%s^4"
            "-(m1.m3*p1.m2*p1.m4+m2.m4*p1.m1*p1.m3+m1.m4*p1.m2*p1.m3+m2.m3*p1.m1*p1.m4 )*3*%s^2"
            "+(%s^2*m1.m2+2*p1.m1*p1.m2)*(%s^2*m3.m4+2*p1.m3*p1.m4)",
            mass,mass,mass,mass); 

  return symb_read(txt);
}



#endif
