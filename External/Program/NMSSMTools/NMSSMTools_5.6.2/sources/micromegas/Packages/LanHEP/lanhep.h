#include <stdio.h>
#include <stdlib.h>

#include "terms.h"



	/*	File main.c   */
void ReadFile(char *);
Term ProcessQuit(Term,Term);
Term ProcUse(Term, Term);
void ProcessTerm(Term);


	/* File  params.c */

Term ProcessParameter(Term , Term);
int  is_parameter(Atom);
void WriteParameters(int fno, char *name);
void WriteExtlib(int fno, char *name);
void WriteCpart(int fno, char *name);
void FAWriteParameters(int fno);
void FADeclRealParam(FILE *);
void FAsqparam(FILE *);
void ClearParameter(Atom);
List all_param_list(void);
Term InterfEvalParam(Term, Term);
double EvalParameter(Term);
cmplx cEvalParameter(Term);
void ChangeParameterValue(Atom, double);
Term ProcTailPrm(Term, Term);
Term ProcCHEPPrm(Term, Term);


	/* File field.c */
Term ProcessParticle(Term, Term);
Term ProcessScalar(Term, Term);
Term ProcessSpinor(Term, Term);
Term ProcessVector(Term, Term);
Term ProcessSpinor3(Term, Term);
Term ProcessTensor(Term, Term);
int is_particle(Atom, Term *);
void WriteParticles(int fno, char *name);
void FAWriteParticles(FILE *f);
Term ProcCC(Term, Term);
Term ProcUp(Term, Term);
Term ProcDown(Term, Term);
Term ProcGhost(Term,Term);
Term ProcAnti(Term, Term);
void ClearParticle(Atom);
void AddIMParticle(Term prt);
List all_prtc_list(void);
Term ProcExtFunc(Term, Term);
List all_prt_list(void);
Term ProcPrtcFormat(Term, Term);
Term ProcPrtcProp(Term, Term);
Term ProcPrtcFunction(Term, Term);

	/* File imprt.c */
List mk_im_field(int col, int spin, int ch, List hint);
Term ProcImPrt(Term, Term);

	/* File group.c */
Term InterfSpecToGroup(Term,Term);
Term SpecToGroup(Term);
Term InterfSpecToRepr(Term,Term);
Term SpecToRepr(Term);
Term ProcessGroup(Term,Term);
Term ProcessRepres(Term, Term);
int is_group(Atom, Term *);
int equal_groups(Term, Term);
void ClearGroup(Atom);


	/* File spec.c */
Term ProcessSpecial(Term,Term);
int is_special(Atom, Term *);
void ClearSpecial(Atom);
Term ProcDelta(Term, Term);

	/* File func.c */
void InitFuncs(void);
int  is_function(Term, Term *);
Term CallFunction(Term,Term);
void alg1_anti(Term);
Term ProcMkProc(Term, Term);

	/* File fu1.c */

Term ProcessModel(Term,Term);
Term ProcessEdit(Term,Term);
Term ProcessRead(Term,Term);
Term GetProp(Term,Term);
Term ProcessDisplay(Term,Term);
Term ProcessWrite(Term,Term);
Term ProcessDate(Term,Term);
Term ProcClear(Term,Term);
Term ProcStat(Term, Term);
Term ProcSetTex(Term, Term);

Term ProcDeriv(Term, Term);
Term ProcReadChep(Term, Term);

	/* File alg1.c */
Term GetIndices(Term,Term);
Term To_t1(Term,Term);
Term SetDefIndex(Term, Term);
Term GetDefIndex(Term, Term);
Term ExprTo1(Term);
Term ExprTo11(Term, List *);
Term ExprTo12(Term, List);
Term ProcVEV(Term, Term);
Term ProcDF(Term, Term);
Term ProcDFDFC(Term, Term);
void alg1_dump(Term);

	/* File alg1x.c */
Term WheredTerm(Term t);        /* a */
Term ProcWhere(Term, Term);
Term SplitIndices(Term t, List *ilist);  /* b */
Term ExpandTerm(Term);          /* b */
Term SetInd1(Term, Term *, Term *);/* c */
Term SetLets(Term);             /* d */
Term SetIntAlgs(Term);          /* i */
Term Alg1ToExpr(Term);          /* w */
Term alg1_mk_wild(Term, List *, List *); /* e */
void alg1_fix_wild(Term);       /* f */
void alg1_sum_wild(Term a1);    /* f */
void alg1_exp_wild(Term a1, List);/* f */
Term ProcTransform(Term, Term); /* h */
Term ProcInfsimal(Term, Term);  /* h */
void alg1_rem_inf(Term);        /* h */
void alg1_fix_delta(Term);      /* h */
Term alg1_guess_mpl(Term);      /* p */
Term alg1_guess_mtr(Term, int *);      /* p */
List alg1_spl_col(List);        /* s */
void alg1_let_cw(Atom);         /* s */
void alg1_let_sw(Atom);         /* s */
Term ProcUED5th(Term, Term);    /* i */

	/* File alg2.c  */
Term Alg1to2(Term);
List alg2_add(List, List);
void alg2_add_ml(Term, List);
void alg2_symmetrize(Term);
void alg2_reduce(List);
void alg2_red_cos(Term);
void alg2_red_orth(Term);
void alg2_red_1pm5(Term);
void alg2_kill_gpm(List);
void alg2_common_t(Term);
void alg2_common_s(Term);
void alg2_recommon_s(Term);
void alg2_common_n(Term);
void alg2_recommon_n(Term);
void alg2_decommon_n(Term);
void alg2_decommon_s(Term);
void alg2_multbyi(Term);
Term ProcRegMatrix(Term, Term);
int herm_matr_dim(Label);
Atom herm_matr_el(Label, int, int);
void tex_write_lagr(List, FILE *);
void tex_wrt_2vrt(FILE *, Term);
int alg2_refine_spinor(List, List *);
int alg2_updown(List prt, List *spec);
List reduce_56(Term);
void alg2_norm_delta(List, Term);
void alg2_eval_vrt(Term);
void FA_write_lagr(List,FILE *);
void FA_write_gen(FILE *);
void UF_write_lagr(List);
void UF_write_gen(void);
Term ProcFainclude(Term, Term);


	/* File util.c */
long int gcf(long int, long int);
List CommaToList(Term t);
List OperToList(Term t, Atom opr);
List Oper1ToList(Term t, Atom opr);
void ErrorInfo(int);
void WarningInfo(int);
List AtomsInTerm(Term, Term);
Term l2plus(List);
Term l2mult(List);
void WriteBlank(FILE *, int);
void RegisterLine(char *);
void UnregisterLine(void);
void DumpRegistered(void);
void ReportRedefined(Atom, char *);
void WriteVertex(List);
void ReplAtom(Term, Atom, Atom);
List XorList(List, Atomic);
Term ProcOption(Term, Term);
Term ProcHermiticity(Term, Term);
List alg2_add_hermconj(List lagr);
Term ProcCheckMasses(Term, Term);
Term ProcClass(Term, Term);
Term il_to_caret(Term, List);
Term ProcIn(Term, List);

	/* File alias.c */
Term SetAlias(Term, Term, int);
Term InterfSetAlias(Term, Term);
Term ProcessAlias(Term);
Term InterfProcessAlias(Term, Term);
void RemoveAlias(Term);
Term InterfRemAlias(Term, Term);

	/* File lagr.c */
Term ProcLTerm(Term,Term);
Term ProcAlg2(Term, Term);
Term WriteLagr(Term,Term);
void WriteLagrFile(int, char *);
void WriteLgrngn(Term, FILE *);
Term ProcSortL(Term, Term);
void Write2Vertex(FILE *, Term prt);
List all_vert_list(void);
Term SelectVertices(Term, Term);
Term ProcLongestLine(Term, Term);
Term ProcReduceLagr(Term, Term);
Term ProcSaveLagr(Term, Term);
Term ProcLoadLagr(Term, Term);
Term ProcCoeff(Term, Term);
Term ProcAddHermConj(Term, Term);
/*Term ProcGenDiagr(Term t, Term ind);*/

	/* File let.c */
Term ProcLet(Term,Term);
int  is_let(Term, List *);
void ClearLet(Atom);
Term ProcKeepLets(Term, Term);
Term ProcOptLets(Term, Term);
Term alg1_inv_alg(Term);

	/* File color.c */

int  color_repres(Term);
void color_symm_f(List, Term);
void epsv_symm(List,Term);
void color_check_spec(Atom, List);

	/* File redcol.c */

List color_reduce(Term);
int  need_col_rdc(Term);
void sort_spin_prt(Term);

	/* File redspin.c */

List spinor_reduce(Term);
int  need_spin_rdc(Term);

	/* File photon.c */

Term ProcSetEM(Term, Term);
void check_em_charge(List);
Atom defined_em_charge(void);

	/* File gauge.c */
Term ProcGauge(Term, Term);
Term ProcGenGauge(Term, Term);
Term ProcBRSTTransf(Term, Term);
Term ProcBRST(Term, Term);
Term ProcCheckBRST(Term, Term);

    /* varver.c */
    
Term VarVer(Term, Term);
Term ProcChVertex(Term, List);
Term ProcFAGC(Term, Term);

	/* angle.c */

Term ProcEval(Term, Term);
Term ProcSetAngle(Term, Term);
Term ProcAngle(Term, Term);
void alg2_red_sico(Term);
void alg2_red_comsico(Term);
void tri_reg_prm(Atom, Term);
Term ProcDbgTrig(Term, Term);

       /* postscript functions */
 void psOpen(const char *fname, int horiz, int pages_no);
 void psClose(void);
 void psNextPage(int horiz);
 void psThick(int thick);
 void psRing(int x, int y, int r1, int r2, int shade);
 void psLineStyle(int st);
 void psLine(int x1, int y1, int x2, int y2, int thick);
 void psRect(int x, int y, int w, int h, int th, int fill);
 void psOutText(int x, int y, char *s, int font, int size);
 int  psTextWidth(char *s, int font, int size);


	/*   Names of frequently used atoms */


extern char *ModelName;
extern int ModelNumber;
extern List DefaultIndex;
extern int  TexOutput, FAOutput, UFOutput, CalcOutput, CompOutput;
extern int  ForsedRedCol;
extern int MicroOmega;

extern int TEX_lines, TEX_spec_in_line, TEX_set_dot;
extern int verb_charge, verb_imprt;
extern int ChepVersion;
