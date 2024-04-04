#include "lanhep.h"


Term wdbg(Term,Term), ProcslhaRead(Term, Term), ProcLHA(Term,Term),
		ProcExtlib(Term, Term), ProcCpart(Term, Term), ProcPrm(Term, Term),
		ProcGenDiagr(Term, Term), ProcCoefVrt(Term, Term), 
		ProcDelVertex(Term, Term), ProcMtrOut(Term, Term);

static struct
	{
	char *cname;
	Atom aname;
	Term (*fu)(Term,Term);
	}  funcs[] = {
		{ "parameter",0,ProcessParameter },
		{ "particle",0,	ProcessParticle },
		{ "field",0,    ProcessParticle },
		{ "scalar",0,   ProcessScalar   },
		{ "spinor",0,   ProcessSpinor   },
		{ "vector",0,   ProcessVector   },
		{ "spinor3",0,   ProcessSpinor3 },
		{ "tensor",0,   ProcessTensor   },
		{ "prtcformat",0,ProcPrtcFormat },
		{ "prtcproperty",0,ProcPrtcProp },
		{ "prtcprop", 0,ProcPrtcProp },
		{ "mass",0,     ProcPrtcFunction},
		{ "width",0,    ProcPrtcFunction},
		{ "echarge",0,  ProcPrtcFunction},
		{ "group",0, 	ProcessGroup    },
		{ "repres",0,   ProcessRepres   },
		{ "special",0,  ProcessSpecial  },
		{ "model",0,    ProcessModel    },
		{ "edit",0,     ProcessEdit     },
		{ "read",0,     ProcessRead     },
		{ "quit",0,     ProcessQuit    },
		{ "SetDefIndex",0,SetDefIndex   },
		{ "GetDefIndex",0,GetDefIndex   },
		{ "getind",0,   GetIndices	  },
		{ "GetProperty",0,  GetProp   },
		{ "display",0,  ProcessDisplay },
		{ "write",0,    ProcessWrite },
		{ "AtomsInTerm",0,AtomsInTerm },
		{ "alias",      0,InterfSetAlias },
		{ "GetProperties",0,GetProperties},
		{ "SetProperty",0,SetProperty},
		{ "SpecToGroup",0,InterfSpecToGroup},
		{ "SpecToRepr",0,InterfSpecToRepr},
		{ "to1",    0,  To_t1 },
		{ "lterm",  0,  ProcLTerm },
		{ "let",   0,   ProcLet },
		{ "ProcAlias",0,InterfProcessAlias },
		{ "unalias",0,InterfRemAlias },
		{ "WriteL",  0, WriteLagr     },
		{ "SortLagr", 0,  ProcSortL  },
		{ "ReduceLagr",0,ProcReduceLagr },
		{ "SaveLagr", 0, ProcSaveLagr },
		{ "LoadLagr", 0,ProcLoadLagr},
		{ "clear", 0,   ProcClear },
		{ "ghost" , 0,  ProcGhost },
		{ "ccghost", 0, ProcGhost  },
		{ "anti",     0,  ProcAnti  },
		{ "cc",      0, ProcCC },
		{ "up",      0, ProcUp },
		{ "down",      0, ProcDown },
		{ "AuxPrt", 0,  ProcImPrt },
		{ "vev",     0, ProcVEV },
		{ "delta", 0,   ProcDelta},
		{ "gauge",  0,  ProcGauge  },
		{ "GenGauge", 0,ProcGenGauge  },
		{ "OrthMatrix", 0, ProcRegMatrix },
		{ "HermMatrix", 0, ProcRegMatrix },
		{ "SetEM", 0, ProcSetEM },
		{ "wdbg", 0,    wdbg	},
		{ "use",  0,    ProcUse 	},
		{ "Statistics",0,ProcStat },
		{ "SetTexName", 0, ProcSetTex },
		{ "SelectVertices", 0, SelectVertices},
		{ "RepLongestLine", 0, ProcLongestLine},
		{ "eval", 0, ProcEval},
		{ "SetAngle", 0, ProcSetAngle},
		{ "dbg_trig", 0, ProcDbgTrig},
		{ "angle", 0, ProcAngle},
		{ "coeff", 0, ProcCoeff},
		{ "transform", 0, ProcTransform},
		{ "infinitesimal", 0, ProcInfsimal},
		{ "option", 0, ProcOption},
		{ "VarVer",0, VarVer },
		{ "CheckHerm",0,ProcHermiticity},
		{ "AddHermConj",0,ProcAddHermConj},
		{ "CheckMasses",0,ProcCheckMasses},
		{ "EvalParameter",0,InterfEvalParam },
		{ "keep_lets", 0, ProcKeepLets},
		{ "opt_lets",0,ProcOptLets},
		{ "df",       0, ProcDF},
		{ "dfdfc",    0, ProcDFDFC},
		{ "tail_prm", 0, ProcTailPrm},
		{ "read_chep_prm", 0, ProcCHEPPrm},
		{ "brst_transform", 0, ProcBRSTTransf},
		{ "brsti_transform", 0, ProcBRSTTransf},
		{ "brst", 0, ProcBRST},
		{ "brsti", 0, ProcBRST},
		{ "CheckBRST", 0, ProcCheckBRST},
		{ "derivp", 0, ProcDeriv},
		{ "ReadChep", 0,ProcReadChep},
		{ "ReadChepM",0,ProcReadChep},
		{ "external_func", 0, ProcExtFunc},
		{ "ChVertex", 0, ProcChVertex},
		{ "DelVertex", 0, ProcDelVertex},
		{ "fa_gencpl",0,ProcFAGC},
		{ "class", 0, ProcClass},
		{ "in", 0, ProcIn},
		{ "mkProc", 0, ProcMkProc},
		{ "fainclude",0, ProcFainclude},
		{ "ued_5th",  0, ProcUED5th},
		{ "where", 0, ProcWhere},
		{ "slhaRead", 0, ProcslhaRead},
		{ "lha", 0, ProcLHA},
		{ "extlib",0,ProcExtlib},
		{ "cpart",0,ProcCpart},
		{ "prm", 0, ProcPrm},
		{ "CoefVrt",0, ProcCoefVrt},
		{ "OutMatrices",0,ProcMtrOut},
		{ "date", 0,    ProcessDate }
		};


void InitFuncs(void)
	{
	int i;
	for(i=0;i<sizeof(funcs)/sizeof(funcs[0]);i++)
		{
		Term tt;
		funcs[i].aname=NewAtom(funcs[i].cname,0);
		tt=MakeCompound1(OPR_FUNCTION,NewInteger(i));
		SetAtomProperty(funcs[i].aname,PROP_TYPE,tt);
		
		}

	}

int is_function(Term t, List *indx)
	{
	Atom name,typ;
	int i;
	if(is_atom(t))
		name=t;
	else
		if(is_compound(t))
			name=CompoundName(t);
		else
			return 0;
	typ=GetAtomProperty(name,PROP_TYPE);
	if(typ==0)
		return 0;
	if(!is_compound(typ) || CompoundName(typ)!=OPR_FUNCTION)
		return 0;
	if(indx!=NULL)
		*indx=GetAtomProperty(name,PROP_INDEX);
	return 1;
	}

int fnotfound(Term t)
	{
	if(!IsTermInput())
		printf("File \"%s\", line %d: ",CurrentInputFile(),
			CurrentInputLine());
		printf("Error: unknown function \'"); WriteTerm(t);
		printf("\'.\n");
		FreeAtomic(t);
		return 0;
	}



Term CallFunction(Term t, Term ind)
	{
	Atom name;
	Term prop;
	int i;
	Term (*fu)(Term,Term);
    if(is_atom(t))
		name=t;
	else
		if(is_compound(t))
			name=CompoundName(t);
		else
			return fnotfound(t);
	prop=GetAtomProperty(name,PROP_TYPE);
	if(!is_compound(prop) ||
		CompoundName(prop)!=OPR_FUNCTION ||
		!is_integer(CompoundArg1(prop)))
			return fnotfound(t);
	i=(int)IntegerValue(CompoundArg1(prop));
	fu=funcs[i].fu;
	return fu(t,ind);
	}

