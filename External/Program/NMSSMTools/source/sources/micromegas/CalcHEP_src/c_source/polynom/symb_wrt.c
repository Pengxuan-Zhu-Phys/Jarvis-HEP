/*
 Copyright (C) 1997,2003  Alexander Pukhov 
*/

#include "polynom.h"
#include "syst2.h"
#include "pvars.h"
#include "tensor.h"
#include "symb_wrt.h"
#include "writeF.h"


void  writepoly(poly p)
{  char txt[STRSIZ], numtxt[STRSIZ];
   int  i, deg;
   int  plus, first;
   int  wpos;
   unsigned long   wpower;

   if(!p){wrt_0("0"); return;}
   if(p->next) plus=0; else plus=1;
 

   for(;p;p=p->next)
   {
      strcpy(txt,"");
      i = 1;
      first = 1;
      wpos = 1;
      wpower = p->power[wpos-1];
      for (i = 1; i <= vardef->nvar; i++)
      {
         deg = (p->power[vardef->vars[i-1].wordpos-1] /
                vardef->vars[i-1].zerodeg) %
                vardef->vars[i-1].maxdeg;

         if(deg > 0)
         {
            if(first)  first = 0; else        strcat(txt,"*");
            strcat(txt,vardef->vars[i-1].name);
            if(deg > 1)  sprintf(txt+strlen(txt),"^%d",deg);
         }
      }

      sprintf(numtxt,"%"NUM_STR, p->num);
       
      if(plus && numtxt[0]!='-') wrt_0("+"); 

      if(strlen(txt))
      {  if(strcmp(numtxt,"1"))
         {
            if(strcmp(numtxt,"-1")==0) wrt_0("-");
            else {wrt_0(numtxt);wrt_0("*");}
         }    
         wrt_0(txt); 
      }
      else wrt_0(numtxt);
      plus=1;
   }
}  


void  writetens(tensor p)
{  char  txt[STRSIZ];
   int   i,s;
   int   plus, first;
   int   isMono;

   if(!p) {wrt_0("0"); return;}

   if(p->next) plus=0; else plus=1;    

   for(; p; p=p->next)
   {
      strcpy(txt,"");
      first = 1;
      for (i = 1; i <= maxIndex; i++)
      {
         s = p->tens[i-1];
         if(s < i && s)
         {
            if(first) first = 0; else strcat(txt,"*"); 
            if(s > 0) sprintf(txt+strlen(txt),"m%d.m%d",i,s); 
            else      sprintf(txt+strlen(txt),"p%d.m%d",-s,i); 
         } 
      }
      
      isMono=0;
      if(p->re && p->re->next==NULL && p->im==NULL) isMono=1;
      if(p->im && p->im->next==NULL && p->re==NULL) isMono=1;

      if(!isMono){if(plus) wrt_0("+("); else wrt_0("(");}
      plus=1;

      if(p->re) writepoly(p->re); 
      if(p->im)
      { if(p->im->next) 
        {  wrt_0("+i*(");
           writepoly(p->im);
           wrt_0(")");
        } else {writepoly(p->im);wrt_0("*i");}
      }      

      if(!isMono) wrt_0(")");
       
      if(strlen(txt)) { wrt_0("*"); wrt_0(txt);} 
   } 
} 


void  writeEtens(Etens p)
{  char  txt[20]; 
   int   i,l; 
   int plus=0;
   if(!p) {wrt_0("0"); return;} 
   for(;p;p=p->next) 
   { 
      if((p->tcoef->next||((p->tcoef->re&&p->tcoef->im))) && p->eps[0]!=X_MARK)
      {  if(plus) wrt_0("+(");else wrt_0("(");   
         writetens(p->tcoef);   wrt_0(")");
      }
      else   writetens(p->tcoef);
      plus=1;
   
      if(p->eps[0]!=X_MARK)
      {  sprintf(txt,"*eps(");
         for(i=0;i<4;i++)
         { int  c=p->eps[i];
           l=strlen(txt);
           if(c>0)sprintf(txt+l,"m%d,",c);else sprintf(txt+l,"p%d,",-c);
         }
         txt[strlen(txt)-1]=')';
         wrt_0(txt);
      }     
   }
}

void  writespinor(SpinTensor p)
{  char  txt[200]; 
   int   i,l; 
   int plus=0;
   if(!p) {wrt_0("0"); return;} 
   for(;p;p=p->next) 
   { 
      if(p->tcoef->next && (p->g5 || p->l))
      { { if(plus) wrt_0("+(");else wrt_0("(");}  writetens(p->tcoef); wrt_0(")");}
      else   writetens(p->tcoef);
         
      if(p->g5 || p->l)
      { 
         sprintf(txt,"*g(ln");
         if(p->g5) strcat(txt,",a");
         for(i=0;i<p->l;i++)
         { int  c=p->g[i];
           l=strlen(txt);
           if(c>0)sprintf(txt+l,",m%d",c);else sprintf(txt+l,",p%d",-c);
         }
         strcat(txt,")");
         wrt_0(txt);
      }     
   }
}


void  symb_print(char* txt, symb_data  m)
{  char s[5];
   int n;
   tensor p;

   wrt_0(txt);
   switch(m.type)  
   { case numbertp:
     case polytp : writepoly(m.expr.p);  break;
     case tenstp : writetens(m.expr.t);  break;
     case spintp : writespinor(m.expr.s);break;
     case etenstp: writeEtens(m.expr.et);break;
     default:
        if(m.type == vectortp)				       
        {								       
           p=m.expr.t;
           if(p == NULL) wrt_0("0");					       
           else							       
           for(;p;p=p->next)							       
           {  n=-p->tens[0];
	      sprintf(s,"p%d*(",n);
              wrt_0(s);						       
              if(p->re) writepoly(p->re);	       
              if(p->im)					       
              { wrt_0("i*("); writepoly(p->re);wrt_0(")");}	       
              wrt_0(")");						       
              if(p->next) wrt_0("+");				       
           }   				       
        }								       
        else								       
        {	
           n = m.expr.t->tens[0];	 
           if(!n) n = 1;						      
           sprintf(s,"l%d",n);
           wrt_0(s);
           writetens(m.expr.t);
           
        }							       
   }
   writeF(";\n");
}


static int  expand_poly(poly p, char***txt, poly**exp)
{
   int Nterm=0;
   *txt=NULL;
   *exp=NULL;
   char buff[STRSIZ]={};
   poly p_=NULL;
   for(;p;)
   { 
     int  i, deg;
     char text[STRSIZ]={}; 
     poly q=p; 
     p=p->next; q->next=NULL; 
     for(i = 0; i < vardef->nvar; i++) if(strchr(vardef->vars[i].name,'.')) 
     {
         deg = (q->power[vardef->vars[i].wordpos-1]/vardef->vars[i].zerodeg) %
                vardef->vars[i].maxdeg;
         if(deg > 0)
         {  q->power[vardef->vars[i].wordpos-1]-=vardef->vars[i].zerodeg*deg; 
            if(strlen(text)) strcat(text,"*");
            strcat(text,vardef->vars[i].name);
            if(deg > 1)  sprintf(text+strlen(text),"^%d",deg);
         }
     }

     if(!p_){ p_=q; strcpy(buff,text);}  
     else  if(strcmp(text,buff)==0)  sewpoly(&p_,&q);
     else 
     {  
        *txt=(char**)realloc(*txt,sizeof(char*)*(Nterm+1));
        (*txt)[Nterm]=malloc(strlen(buff)+10);
        strcpy((*txt)[Nterm],buff);
        *exp=(poly*)realloc(*exp,sizeof(poly)*(Nterm+1)); 
        (*exp)[Nterm]=p_;
        Nterm++;
        p_=q;
        strcpy(buff,text); 
      }      
   }
   if(p_)
     { 
        *txt=(char**)realloc(*txt,sizeof(char*)*(Nterm+1));
        (*txt)[Nterm]=malloc(strlen(buff)+10);
        strcpy((*txt)[Nterm],buff);
        *exp=(poly*)realloc(*exp,sizeof(poly)*(Nterm+1));
        (*exp)[Nterm]=p_;
        Nterm++;
      }  
   return Nterm;  
}  

static int  expand_tensor(tensor p, char***txt, poly**exp)
{  char  text[STRSIZ],buff[STRSIZ];
   int   i,s;
   int   plus, first;
   int   isMono;

   *txt=NULL;
   *exp=NULL;   
    int Nterm=0; 

   for(; p; p=p->next)
   {
      strcpy(text,"");
      first = 1;
      for (i = 1; i <= maxIndex; i++)
      {
         s = p->tens[i-1];
         if(s < i && s)
         {
            if(first) first = 0; else strcat(text,"*"); 
            if(s > 0) sprintf(text+strlen(text),"m%d.m%d",i,s); 
            else      sprintf(text+strlen(text),"p%d.m%d",-s,i); 
         } 
      }

      char**txtRe=NULL,**txtIm=NULL;  
      poly*expRe=NULL,*expIm=NULL;
      int Nre=0,Nim=0;
      if(p->re)
      { Nre=expand_poly(p->re,&txtRe,&expRe);
        for(i=0;i<Nre;i++) 
        { txtRe[i]=realloc(txtRe[i], strlen(txtRe[i])+strlen(text)+2);
          if(strlen(txtRe[i])&& strlen(text)) { sprintf(buff,"%s*%s",text,txtRe[i]); strcpy(txtRe[i],buff); }
          else if(strlen(text)) strcpy(txtRe[i],text); 
        }
      }       

      if(p->im)
      { Nim=expand_poly(p->im,&txtIm,&expIm);
        for(i=0;i<Nim;i++) 
        { txtIm[i]=realloc(txtIm[i], strlen(txtIm[i])+strlen(text)+4);
          if(strlen(txtIm[i])&&strlen(text)) { sprintf(buff,"%s*i*%s",text,txtIm[i]);strcpy(txtIm[i],buff);}
          else if(strlen(txtIm[i])==0 && strlen(text)==0) strcpy(txtIm[i],"i");
          else if(strlen(text))      sprintf(txtIm[i],"%s*i",text);
          else if(strlen(txtIm[i])) { sprintf(buff, "i*%s",txtIm[i]); strcpy(txtIm[i],buff);}                   
        }
      }

      *txt=(char**)realloc(*txt,(Nterm+Nre+Nim)*sizeof(char*));
      *exp=(poly*) realloc(*exp,(Nterm+Nre+Nim)*sizeof(poly));
      for(i=0;i<Nre;i++) {(*txt)[Nterm+i]=txtRe[i]; (*exp)[Nterm+i]=expRe[i];}
      Nterm+=Nre;free(txtRe);free(expRe);  
      for(i=0;i<Nim;i++) {(*txt)[Nterm+i]=txtIm[i]; (*exp)[Nterm+i]=expIm[i];}
      Nterm+=Nim; free(txtIm); free(expIm);     
   }
   return Nterm;   
} 

static int  expand_spinor(SpinTensor p, char***txt, poly**exp)
{  char  text[STRSIZ],buff[STRSIZ];
   int   i,l;
   *txt=NULL;
   *exp=NULL;   
    int Nterm=0; 

   for(; p; p=p->next)
   {
      strcpy(text,"");
      if(p->g5 || p->l)
      { 
         if(p->g5) strcat(text,"G5");
         for(i=0;i<p->l;i++)
         { int  c=p->g[i];
           l=strlen(text);
           if(l) {strcat(text,"*");l++;}
           if(c>0)sprintf(text+l,"G(m%d)",c);else sprintf(text+l,"G(p%d)",-c);
         }
      }     
      char**txtT=NULL;  
      poly*expT=NULL;
      int NT=0;
      
      NT=expand_tensor(p->tcoef,&txtT,&expT);
        for(i=0;i<NT;i++) 
        {
          txtT[i]=realloc(txtT[i], strlen(txtT[i])+strlen(text)+2);
          if(strlen(txtT[i])) { sprintf(buff,"%s*%s",text,txtT[i]); strcpy(txtT[i],buff); }
          else strcpy(txtT[i],text); 
        }      

      *txt=(char**)realloc(*txt,(Nterm+NT)*sizeof(char*));
      *exp=(poly*) realloc(*exp,(Nterm+NT)*sizeof(poly));
      for(i=0;i<NT;i++) {(*txt)[Nterm+i]=txtT[i]; (*exp)[Nterm+i]=expT[i];}
      Nterm+=NT;free(txtT);free(expT);
   }
   return Nterm;   
} 

static int  expand_Etens(Etens p, char***txt, poly**exp)
{  char  text[STRSIZ],buff[STRSIZ];
   int   i,l;
   *txt=NULL;
   *exp=NULL;   
    int Nterm=0; 

   for(; p; p=p->next)
   {
      strcpy(text,"");
      if(p->eps[0]!=X_MARK)
      {  sprintf(text,"eps(");
         for(i=0;i<4;i++)
         { int  c=p->eps[i];
           l=strlen(text);
           if(c>0)sprintf(text+l,"m%d,",c);else sprintf(text+l,"p%d,",-c);
         }
         text[strlen(text)-1]=')';
      }     
      
      char**txtT=NULL;  
      poly*expT=NULL;
      int NT=0;
      NT=expand_tensor(p->tcoef,&txtT,&expT);
        for(i=0;i<NT;i++) 
        { txtT[i]=realloc(txtT[i], strlen(txtT[i])+strlen(text)+2);
          if(strlen(txtT[i])) { sprintf(buff,"%s*%s",text,txtT[i]); strcpy(txtT[i],buff); }
          else strcpy(txtT[i],text); 
        }      

      *txt=(char**)realloc(*txt,(Nterm+NT)*sizeof(char*));
      *exp=(poly*) realloc(*exp,(Nterm+NT)*sizeof(poly));
      for(i=0;i<NT;i++) {(*txt)[Nterm+i]=txtT[i]; (*exp)[Nterm+i]=expT[i];}
      Nterm+=NT;free(txtT);free(expT);
   }
   return Nterm;   
} 

int   symb_expand(symb_data  m, char***txt,poly**exp)
{  char s[5];
   int Nterm;

   switch(m.type)  
   { case numbertp:
     case polytp : Nterm= expand_poly(m.expr.p,  txt,exp);break;
     case tenstp : Nterm= expand_tensor(m.expr.t,  txt,exp);break;
     case spintp : Nterm= expand_spinor(m.expr.s,txt,exp);
     break;
     case etenstp: Nterm= expand_Etens(m.expr.et,txt,exp);break;
     default:      Nterm=0; *txt=NULL; *exp=NULL;
   } 
   return Nterm;     
}
