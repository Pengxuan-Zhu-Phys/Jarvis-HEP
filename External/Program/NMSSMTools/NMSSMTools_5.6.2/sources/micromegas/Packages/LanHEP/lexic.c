#include <string.h>
#include <ctype.h>
#include <math.h>
#include <setjmp.h>
#include <errno.h>
#ifndef __ZTC__
#include <unistd.h>
#endif
#include "lanhep.h"


static jmp_buf lex_jmp_buf;

static struct rrec
	{
	int  term;
	int  diag, sec_prompt;
	
	char name[80];
	FILE *file;
	Atomic a[5024];
	int  line;
	int  eof,comment;
	int acnt, cura;
	Atomic unread;
	int was_unread;
	}   rrecs[16];

static int currec = -1;

static char cbuf[10240];

static void flushbuff(void)
	{
	if(currec==-1)
		{puts("Internal error (no open file)"); exit(0); }
	if(rrecs[currec].file==NULL)
	{
		if(rrecs[currec].line)
			rrecs[currec].eof=1;
		else
			rrecs[currec].line++;
		return;
	}
			
	if(rrecs[currec].diag && rrecs[currec].sec_prompt)
		 printf("   | "); /* printf("  | ");*/
rd:	errno=0;
	if(fgets(cbuf,10240,rrecs[currec].file)==NULL)
		{
		if(errno)
			{
			printf("Error reading file: "); perror("");
			}
		rrecs[currec].eof=1;
		}
	if(rrecs[currec].eof!=1 && cbuf[0]=='#' && cbuf[1]=='!')
		goto rd;
/*	if(rrecs[currec].eof!=1 && cbuf[0]=='C')
		goto rd;
*/
	rrecs[currec].sec_prompt=1;
	rrecs[currec].line++;

	}
	
	
void WritePrompt(int no)
	{
	if(rrecs[currec].diag)
		{
		printf("%-3d> ",no);  /*printf("- > ");*/
		rrecs[currec].sec_prompt=0;
		rrecs[currec].line=no-1;
		}
	}
	
int SetInputFile(char *name)
	{
	currec++;
	if(currec==16)
		{currec--; return 0; }
	rrecs[currec].eof=0;
	rrecs[currec].comment=0;
	rrecs[currec].line=0;
	rrecs[currec].acnt=0;
	rrecs[currec].cura=0;
	rrecs[currec].diag=rrecs[currec].sec_prompt=0;
	rrecs[currec].was_unread=0;
	if(name==NULL)
		{
		rrecs[currec].term=1;
		rrecs[currec].file=stdin; 
		}
	else
		{
		rrecs[currec].term=0;
		rrecs[currec].file=fopen(name,"r");
		if(rrecs[currec].file==NULL)
			{ currec--; return 0; }
		strncpy(rrecs[currec].name,name,70);
		}
	if(isatty(fileno(rrecs[currec].file))) { rrecs[currec].diag=1; }
	return 1;
	}

int SetInputFileM(char *line)
	{
	currec++;
	if(currec==16)
		{currec--; return 0; }
	rrecs[currec].eof=0;
	rrecs[currec].comment=0;
	rrecs[currec].line=0;
	rrecs[currec].acnt=0;
	rrecs[currec].cura=0;
	rrecs[currec].diag=rrecs[currec].sec_prompt=0;
	rrecs[currec].was_unread=0;
	rrecs[currec].file=NULL;
	rrecs[currec].term=0;
	strcpy(rrecs[currec].name,"(memory)");
	strcpy(cbuf,line);
	return 1;
	}
	
int IsTermInput(void)
	{
	return rrecs[currec].term;
	}
	
void CloseInputFile(void)
	{
	if(rrecs[currec].term==0 && rrecs[currec].file!=NULL )
		fclose(rrecs[currec].file);	
	currec--;
	}
	
void SkipInputLine(void)
	{
	rrecs[currec].was_unread=0;
	rrecs[currec].acnt=0;
	}
	
	
static Atomic read_a(char *, int *);

static int skip_quoted(char *s, char stop)
	{
	int alen=0;
	s++;
	while(s[alen]!=0)
		{
		if(s[alen]==stop)
			  return alen+2;
		if(s[alen]=='\\')
			{
			int i=alen;
			while(s[i]!=0)
				{ s[i]=s[i+1]; i++; }
			}
		alen++;
		}
	return alen+1;
	}

static int end_of_comment(int curpos)
	{
	curpos+=2;
	while(cbuf[curpos]!=0)
		{
		if(cbuf[curpos]=='*' && cbuf[curpos+1]=='/')
			return curpos+2;
		if(cbuf[curpos]=='\'')
			{ 	curpos+=skip_quoted(cbuf+curpos,'\'');
				continue; }
		if(cbuf[curpos]=='\"')
			{ 	curpos+=skip_quoted(cbuf+curpos,'\"');
				continue; }
		curpos++;
		}
	return curpos;
	}
				
	
	
static void readatomics(void)
	{
	int curpos,len;
	rrecs[currec].acnt=0;
	rrecs[currec].cura=0;
begin:
	flushbuff();
	if(rrecs[currec].eof==1 && rrecs[currec].comment==1)
		{
		longjmp(lex_jmp_buf,2);
		printf("\nEnd of file \"%s\" within comment /*...*/ \n",
		CurrentInputFile());
		exit(0);
		}
	
	if(rrecs[currec].eof==1) return;
	len=(int)strlen(cbuf);
	curpos=0;
	if(rrecs[currec].comment==1)
		{
		curpos=end_of_comment(-2);
		if(cbuf[curpos]!=0)
			rrecs[currec].comment=0;
		else
			goto begin;
		}
	
	while((isspace(cbuf[curpos])||iscntrl(cbuf[curpos])) && curpos<len)
			curpos++;
	if(curpos==len) goto begin;
	if(cbuf[curpos]=='%') goto begin;
	
	while(curpos<len)
		{
		int alen;
		if(cbuf[curpos]=='%')
			{
			if(rrecs[currec].acnt==0) goto begin;
			else return;
			}	
		if(cbuf[curpos]=='/' && cbuf[curpos+1]=='*')
			{
			curpos=end_of_comment(curpos);
			if(cbuf[curpos]==0)
				rrecs[currec].comment=1;
			goto lll;
			}
			
		rrecs[currec].a[rrecs[currec].acnt]=read_a(cbuf+curpos,&alen);
		rrecs[currec].acnt++;
		curpos+=alen;
lll:		while((isspace(cbuf[curpos])||iscntrl(cbuf[curpos])) 
				&& curpos<len)
			curpos++;
		}
		
	if(rrecs[currec].acnt==0) goto begin;
	}
				

Atomic ReadAtomic(void)
	{
	int err;
	if(rrecs[currec].was_unread)
		{
		rrecs[currec].was_unread=0;
		return rrecs[currec].unread;
		}
	err=setjmp(lex_jmp_buf);
	if(err!=0)
		{
		ErrorInfo(err);
		printf(" end of file \"%s\" within comment /*...*/ \n",
				CurrentInputFile());
		
		return 0;	
		}
			
	if(rrecs[currec].cura>=rrecs[currec].acnt)
		readatomics();
	if(rrecs[currec].eof==1 || rrecs[currec].acnt==0)
		return 0;
	rrecs[currec].cura++;
	return rrecs[currec].a[rrecs[currec].cura-1];
	}
	
	
static int alone(char c)
	{
	if(	c=='.' ||
		c==',' ||
		c==';' ||
		c=='(' ||
		c==')' ||
		c=='[' ||
		c==']' ||
		c=='{' ||
		c=='}' ||
		c=='\"'||
		c=='\'')
			return 1;
	return 0;
	}
	
static Atomic read_quoted(char *s, int *len, char stop)
	{
	int alen=0;
	s++;
	while(s[alen]!=0)
		{
		if(s[alen]==stop)
			{ *len=alen+2;
			  return NewAtom(s,alen);
			  }
		if(s[alen]=='\\')
			{
			int i=alen;
			while(s[i]!=0)
				{ s[i]=s[i+1]; i++; }
			}
			
		alen++;
		}
	*len=alen+1;
	ErrorInfo(1);
	printf(" unbalanced %c \n",stop);
	return NewAtom(s,alen); 
	}
	

		
	

static Atomic read_a(char *s, int *len)
	{
	if(s[0]=='(')
		{ *len=1; return A_RBRA; }
	if(s[0]==')')
		{ *len=1; return A_RCET; }
	if(s[0]=='{')
		{ *len=1; return A_FBRA; }
	if(s[0]=='}')
		{ *len=1; return A_FCET; }
	if(s[0]=='[')
		{ *len=1; return A_QBRA; }
	if(s[0]==']')
		{ *len=1; return A_QCET; }
	if(s[0]=='.')
		{ *len=1; return A_POINT; }
	if(s[0]==';')
		{ *len=1; return A_SECO; }
	if(s[0]==',')
		{ *len=1; return A_COMMA; }
	if(s[0]=='\'')
		return read_quoted(s,len,'\'');
	if(s[0]=='\"')
		return read_quoted(s,len,'\"');
	if(s[0]=='`')
		{ *len=1; return A_RQUOTE; }
	
	if(isalpha(s[0]) || s[0]=='_' || s[0]=='~')
		{
		int alen;
		alen=0;
		while(isalnum(s[alen]) || s[alen]=='_' || s[alen]=='~' )
			alen++;
		*len=alen;
		return NewAtom(s,alen);
		}
	/* test for numbers */
	if(isdigit(s[0]) /*|| (s[0]=='-' && isdigit(s[1]))*/ )
		{
		int negflag=0;
		int alen=0;
		int lll;
		if(s[0]=='-')
			{negflag=1; s++; alen++; }
		lll=0;
		while(isdigit(s[lll])) lll++;
		if((s[lll]=='e' || s[lll]=='E') && (isdigit(s[lll+1]) || s[lll+1]=='-'))
			{
			double f=0.0;
			int eneg=0,exp=0;
			while(s[0]!='e' && s[0]!='E')
				{
				f*=10.0;
				f+=s[0]-'0';
				s++;
				alen++;
				}
			if(negflag) f=-f;
			s++; alen++;
			if(s[0]=='-')
				{
				s++; alen++; eneg=1;
				}
			if(s[0]=='+')
				{
				s++; alen++;
				}
			while(isdigit(s[0]))
				{
				exp*=10;
				exp+=s[0]-'0';
				s++; alen++;
				}
			if(eneg) exp=-exp;
			f*=pow(10.0,exp);
			*len=alen;
			return NewFloat(f);
			}
				
		if(s[lll]=='.' && isdigit(s[lll+1]))
			{
			double f=0.0;
			double fact=1.0;
			while(s[0]!='.')
				{
				f*=10.0;
				f+=s[0]-'0';
				s++;
				alen++;
				}
			s++; alen++;
			while(isdigit(s[0]))
				{
				fact/=10.0;
				f+=fact*(s[0]-'0');
				s++; alen++;
				}
			if(negflag) f=-f;
			if(s[0]=='e' || s[0]=='E')
				{
				int eneg=0;
				int exp=0;
				s++; alen++;
				if(s[0]=='-')
					{
					s++; alen++; eneg=1;
					}
				if(s[0]=='+')
					{
					s++; alen++;
					}
				while(isdigit(s[0]))
					{
					exp*=10;
					exp+=s[0]-'0';
					s++; alen++;
					}
				if(eneg) exp=-exp;
				f*=pow(10.0,exp);
				}
			*len=alen;
			
			return NewFloat(f);
			}
		lll=0;
		while(isdigit(s[0]))
			{
			lll*=10;
			lll+=s[0]-'0';
			s++;
			alen++;
			}
		if(negflag) lll=-lll;
		*len=alen;
		return NewInteger(lll);
		}
				
			
		
	if(ispunct(s[0]))
		{
		int alen;
		alen=0;
		while(ispunct(s[alen]) && !alone(s[alen]) && 
			!(alen>0 && s[alen-1] == '=' && s[alen]!= '=') )
			alen++;
		*len=alen;
		return NewAtom(s,alen);
		}
	
			
	*len=1;
	return NewAtom("???",0);
	}
	
int CurrentInputLine(void)
	{
	return rrecs[currec].line;
	}
	
char *CurrentInputFile(void)
	{
	if(rrecs[currec].term==1)
		return "tty";
	else
		return rrecs[currec].name;
	}
	
void UnreadAtomic(Atomic a)
	{
	rrecs[currec].was_unread=1;
	rrecs[currec].unread=a;
	}
