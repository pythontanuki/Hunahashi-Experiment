#include "sbml/SBMLTypes.h"
#include <stdbool.h>
#include <stdio.h>
#define uint unsigned int
#define rep(i, n) for (int i = 0; i < n; i++)
#define rep2(i,a,b) for (int i = a; i < b; i++)
#define M 10000
#define getnewline() puts("")
#define Ndbg(X) printf("%.10lf\n", X)
#define Sdbg(X) printf("%s\n", X)

Species_t *Sp[M];
Reaction_t *Reac[M];
ListOf_t *React_List[M];
ListOf_t *Prod_List[M];
KineticLaw_t *KineticLaw[M];
const ASTNode_t *AST_NumericalFormula[M];
double Initial_Amount[M], val[M], ReacVals[M], dv[M];
const char *Ids[M];
uint ReactEachSize[M], ProdEachSize[M];

bool checkType(ASTNodeType_t type) {
	if (type == AST_INTEGER || type == AST_REAL)
		return true;
	return false;
}

uint getSize(const char *s) {
	uint n = 0;
	while (*s++ != '\0')
		n++;
	return (n);
}

bool compareId(const char *s1, const char *s2) {
	uint l1 = getSize(s1);
	uint l2 = getSize(s2);
	if (l1 != l2) {
		return false;
	}
	for (uint i = 0; i < l1; i++) {
		if (s1[i] != s2[i])
			return false;
	}
	return true;
}

double applyParameterAndVariables(const ASTNode_t *v,KineticLaw_t *kl,int Splist_size) {
	double res = 0;
	// variables or parameters
	// if (ASTNode_isName(v)) {
	const char *id = ASTNode_getName(v);
	bool find_Id = false;
	rep(i, Splist_size) if (compareId(id, Species_getId(Sp[i]))) {
		find_Id = true, res = Initial_Amount[i];
		break;
	}
	if (!find_Id) {
		Parameter_t *p = KineticLaw_getParameterById(kl, id);
		res = Parameter_getValue(p);
	}
	// } else
	// 	res = (ASTNode_isInteger(v) ? ASTNode_getInteger(v) :
	// ASTNode_getReal(v));
	return res;
}

double solve(ASTNodeType_t op, double left_val, double right_val) {
	if (op == AST_PLUS)
		return left_val + right_val;
	else if (op == AST_MINUS)
		return left_val - right_val;
	else if (op == AST_TIMES)
		return left_val * right_val;
	else if (op == AST_DIVIDE)
		return left_val / right_val;
	return pow(left_val, right_val);
}

double dfsOnAST_BinaryTree(const ASTNode_t *v, KineticLaw_t *kl, int Splist_size) {
	double left_val = 0, right_val = 0;
	// bool left_isNot_operator = false, right_isNot_operator = false;
	// If we get to the leaf, we should end this function.
	ASTNodeType_t op = ASTNode_getType(v);
	if (ASTNode_getNumChildren(v) == 0) { // leaf
		if (op == AST_INTEGER)
			return ASTNode_getInteger(v);
		if (op == AST_REAL)
			return ASTNode_getReal(v);
		if (ASTNode_isName(v))
			return applyParameterAndVariables(v, kl, Splist_size);
	}
	ASTNode_t *left = ASTNode_getLeftChild(v);
	ASTNode_t *right = ASTNode_getRightChild(v);
	// ASTNodeType_t left_type = ASTNode_getType(left);
	// ASTNodeType_t right_type = ASTNode_getType(right);
	// checkType(left_type) ? (left_isNot_operator = true) : (left_val =
	// dfs(left, kl, Splist_size)); checkType(right_type) ?
	// (right_isNot_operator = true) : (right_val = dfs(right, kl,
	// Splist_size));
	left_val = dfsOnAST_BinaryTree(left, kl, Splist_size);
	right_val = dfsOnAST_BinaryTree(right, kl, Splist_size);
	// In this state, We should confirm types and determine l's and r's
	// parameter. This will be acheived by 2 steps, Numbers State/Otherwise.
	// We do left part here, but right will be done as well.
	// if (left_isNot_operator)
	// 	left_val = applyParameterAndVariables(left, kl, Splist_size);
	// if (right_isNot_operator)
	// 	right_val = applyParameterAndVariables(right, kl, Splist_size);
	return solve(op, left_val, right_val);
}


int main() {
	SBMLDocument_t *d;
	d = readSBML("mapk.xml");
	Model_t *m = SBMLDocument_getModel(d);
	ListOf_t *Sp_list = Model_getListOfSpecies(m);
	ListOf_t *Reac_list = Model_getListOfReactions(m);

	uint Splist_size = ListOf_size(Sp_list), Reac_size = ListOf_size(Reac_list);

	// get each sp's information and puts them
	rep(i, Splist_size) {
		Sp[i] = (Species_t *)ListOf_get(Sp_list, i);
		Initial_Amount[i] = Species_getInitialAmount(Sp[i]);
		Ids[i] = Species_getId(Sp[i]);
	}
    rep(i,Splist_size) {
        Ndbg(Initial_Amount[i]);
        Sdbg(Species_getName(Sp[i]));
        Sdbg(Ids[i]);
    }
	// prepare for operating information with AST_Node type
	rep(i, Reac_size) {
		Reac[i] = (Reaction_t *)ListOf_get(Reac_list, i);
		React_List[i] = Reaction_getListOfReactants(Reac[i]);
		ReactEachSize[i] = ListOf_size(React_List[i]);
		Prod_List[i] = Reaction_getListOfProducts(Reac[i]);
		ProdEachSize[i] = ListOf_size(Prod_List[i]);
		KineticLaw[i] = Reaction_getKineticLaw(Reac[i]);
		AST_NumericalFormula[i] = KineticLaw_getMath(KineticLaw[i]);
	}
  double t = 0, t_max = 4000, dt = 1.0;
	double results[Splist_size][M];
	int cnt = 0;
	rep(i, Splist_size) results[i][0] = Initial_Amount[i];
	// Simulation State
	while (t <= t_max) {
		rep(i, Reac_size) {
			ReacVals[i] = dfsOnAST_BinaryTree(AST_NumericalFormula[i], KineticLaw[i], Splist_size);
			// Ndbg(ReacVals[i]);
		}
		t += dt;
		cnt++;
		rep(i, Splist_size) {
			dv[i] = 0;
			rep(j, Reac_size) {
				rep(k, ReactEachSize[j]) {
					SpeciesReference_t *R =
						(SpeciesReference_t *)ListOf_get(React_List[j], k);
					if (compareId(Ids[i], SpeciesReference_getSpecies(R)))
						dv[i] -= ReacVals[j];
				}
				rep(k, ProdEachSize[j]) {
					SpeciesReference_t *P =
						(SpeciesReference_t *)ListOf_get(Prod_List[j], k);
					if (compareId(Ids[i], SpeciesReference_getSpecies(P)))
						dv[i] += ReacVals[j];
				}
			}
			Initial_Amount[i] += dv[i] * dt;
			results[i][cnt] = Initial_Amount[i];
		}
	}
	// Confirm Results
	FILE *fp = fopen("res.csv", "w");
	rep(i, 4001) {
		double time = (double)i;
		fprintf(fp, "%lf,", time);
		rep(j, Splist_size) fprintf(fp, "%.10lf,", results[j][i]);
		fprintf(fp, "\n");
	}
	SBMLDocument_free(d);
	return 0;
}
