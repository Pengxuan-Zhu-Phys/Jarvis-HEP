#include<stdio.h>
#include"micromegas_aux.h"

void  pythonversion_(int *n1,int *n2)
{ 
    char version[40];     
    system("python --version 2> python_version");
    FILE *f=fopen("python_version","r");
    if(!f){ *n1=0; *n2=0; return;}
    fscanf(f,"%*s %d.%d",n1,n2);
    fclose(f);
    unlink("python_version");
}    