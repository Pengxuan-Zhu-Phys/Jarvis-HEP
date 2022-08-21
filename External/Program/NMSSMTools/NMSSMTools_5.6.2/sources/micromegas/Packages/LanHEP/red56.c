#include "lanhep.h"
#include "math.h"


static int tot_cnt=0;
		
List reduce_56(Term a2)
	{
	List prt;
	prt=CompoundArg1(a2);
	if(ListLength(prt)<5)
		return AppendFirst(NewList(),a2);
	if(tot_cnt<10)
	{
		printf("Error: Vertex with more than 4 particles: ");
		WriteVertex(prt);
		puts("");
	}
	if(tot_cnt==10)
		puts("More vertices with more than 4 particles follow...");
	tot_cnt++;
	/*
	DumpList(CompoundArgN(ListFirst(CompoundArgN(a2,5)),3));
	*/
	return NewList();
	}
	
