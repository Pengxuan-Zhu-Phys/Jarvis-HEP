#include <stdio.h>
#include <string.h>
#include "lanhep.h"

#define ITR_OA 50
#define ITR_NA 51
#define ITR_CO 52
#define ITR_LI 53
#define ITR_LE 54
#define ITR_LA 55
#define ITR_IN 56
#define ITR_FP 57
#define ITR_EM 58


#define AHASHSZ 2000

static Atom at_out[AHASHSZ], at_in[AHASHSZ];
static FILE *inf=NULL, *outf=NULL;
static int in_pos, out_pos;


int itrSetIn(char *file)
	{
	int i;
	for(i=0;i<AHASHSZ;i++)
		at_in[i]=0;
	in_pos=0;
	inf=fopen(file,"r");
	if(inf==NULL)
		{
		return 0;
		}
	return 1;
	}

Term itrIn(void)
	{
	int c;
	c=getc(inf);
	if(c==EOF)
		return 0;
	switch(c)
		{
		case ITR_LE:
			return ITR_LE;
		case ITR_EM:
			return 0;
		case ITR_IN:
			{
			int i;
			fread(&i,sizeof(int),1,inf);
			return NewInteger(i);
			}
		case ITR_FP:
			{
			double d;
			fread(&d,sizeof(double),1,inf);
			return NewFloat(d);
			}
		case ITR_OA:
			{
			unsigned int n;
			unsigned char b;
			b=getc(inf);
			n=b;
			b=getc(inf);
			n+=b*256;
			return at_in[n];
			}
		case ITR_NA:
			{
			unsigned int n;
			unsigned char b;
			char buf[1024];
			Atom a;
			b=getc(inf);
			n=b;
			b=getc(inf);
			n+=b*256;
			fread(&buf,n,1,inf);
            a=NewAtom(buf,n);
			at_in[in_pos]=a;
			in_pos++;
			if(in_pos==AHASHSZ)
				in_pos=0;
			return a;
			}
		case ITR_LA:
			{
			Label lb;
			fread(&lb,sizeof(Label),1,inf);
			return lb;
			}

		case ITR_LI:
			{
			List l;
			Term t;
			l=NewList();
			while((t=itrIn()))
				{
				if(t==ITR_LE)
					return l;
				l=AppendLast(l,t);
				}
			puts("Internal error (itril)");
			return 0;
			}
		case ITR_CO:
			{
			Term t;
			int ar,i;
			Atom nm;
			ar=getc(inf);
			nm=itrIn();
			if(!is_atom(nm))
				{
				puts("Internal error (itria)");
				return 0;
				}
			t=MakeCompound(nm,ar);
			for(i=1;i<=ar;i++)
				SetCompoundArg(t,i,itrIn());
			return t;
			}
		default:
			printf("Internal error (itr %d)\n",c);
			return 0;
		}
	}

int itrSetOut(char *file)
	{
	int i;
	for(i=0;i<AHASHSZ;i++)
		at_out[i]=0;
	out_pos=0;
	outf=fopen(file,"w");
	if(outf==NULL)
		{
		ErrorInfo(78);
		printf("can not open file %s for writing\n",file);
		perror("");
		return 0;
		}
	return 1;
	}

static void out_a(Atom a)
	{
	int i;
	unsigned char tp;
	for(i=0;i<AHASHSZ;i++)
		{
		if(at_out[i]==a)
			{
			tp=ITR_OA;
			putc(tp,outf);
			tp=(i&255);
			putc(tp,outf);
			tp=i/256;
			putc(tp,outf);
			return;
			}
		}

	tp=ITR_NA;
	putc(tp,outf);
	at_out[out_pos]=a;
	out_pos++;
	if(out_pos==AHASHSZ)
		out_pos=0;
	i=strlen(AtomValue(a));
	tp=(i&255);
	putc(tp,outf);
	tp=i/256;
	putc(tp,outf);
	fwrite(AtomValue(a),i,1,outf);
	}

static void out_c(Compound c)
	{
	unsigned char tp;
	int i;
	putc(ITR_CO,outf);
	tp=CompoundArity(c);
	putc(tp,outf);
	out_a(CompoundName(c));
	for(i=1;i<=tp;i++)
		itrOut(CompoundArgN(c,i));
	}

static void out_l(List l)
	{
	unsigned char tp;
	tp=ITR_LI;
	putc(tp,outf);
	while(!is_empty_list(l))
		{
		itrOut(ListFirst(l));
		l=ListTail(l);
		}
	tp=ITR_LE;
	putc(tp,outf);
	}
	
			
	
void itrOut(Term t)
	{
	char typ;
	typ=AtomicType(t);
	switch(typ)
		{
		case 'a':
			out_a(t);
			break;
		case 'c':
			out_c(t);
			break;
		case 'l':
			out_l(t);
			break;
		case 'i':
			{
			char tp;
			int oi;
			tp=ITR_IN;
			oi=IntegerValue(t);
			putc(tp,outf);
			fwrite(&oi,sizeof(int),1,outf);
			break;
			}
		case 'f':
			{
			char tp;
			double of;
			tp=ITR_FP;
			of=FloatValue(t);
			putc(tp,outf);
			fwrite(&of,sizeof(double),1,outf);
			break;
			}
		case 'L':
			{
			char tp;
			tp=ITR_LA;
			putc(tp,outf);
			fwrite(&t,sizeof(Label),1,outf);
			break;
			}		
		case 'e':
			{
			char tp;
			tp=ITR_EM;
			putc(tp,outf);
			break;
			}
		default:
			printf("Internal error (itrio t %c)\n",typ);
		}
	
	}
	
void itrCloseOut(void)
	{
	fclose(outf);
	}
	
void itrCloseIn(void)
	{
	fclose(inf);
	}
