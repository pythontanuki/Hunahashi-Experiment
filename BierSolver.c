#include "sbml/SBMLTypes.h"
#include <stdbool.h>
#include <stdio.h>
#define uint unsigned int
#define rep(i, n) for (int i = 0; i < n; i++)
#define rep2(i, a, b) for (int i = a; i < b; i++)
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

double applyParameterAndVariables(const ASTNode_t *v,KineticLaw_t *kl,int Splist_size, Model_t* m) {
	double res = 0;
	const char *id = ASTNode_getName(v);
	bool find_Id = false;
	rep(i, Splist_size) if (compareId(id, Species_getId(Sp[i]))) {
		find_Id = true, res = Initial_Amount[i];
		break;
	}
	if (!find_Id) {
        ListOf_t* lop = Model_getListOfParameters(m);
        Parameter_t *p = ListOfParameters_getById(lop,id);
		res = Parameter_getValue(p);
	}
	return res;
}

double solve(ASTNodeType_t op, double left_val, double right_val) {
	if (op == AST_PLUS) return left_val + right_val;
	else if (op == AST_MINUS) return left_val - right_val;
	else if (op == AST_TIMES) return left_val * right_val;
	else if (op == AST_DIVIDE) return left_val / right_val;
	return pow(left_val, right_val);
}

ASTNode_t* makeBinarySubtree(const ASTNode_t *v) {
  ASTNode_t* SubTree = ASTNode_deepCopy(v);
  ASTNodeType_t operateType = ASTNode_getType(v);
  int childCnt = ASTNode_getNumChildren(v);
  ASTNode_t* tmpVertexVector[childCnt-2];
  rep(i,childCnt-2) {
    ASTNode_t* operateNode = ASTNode_create();
    ASTNode_setCharacter(operateNode, operateType);
    ASTNode_t* child = ASTNode_getChild(v, i+2);
    ASTNodeType_t childType = ASTNode_getType(child);
    ASTNode_t* newChild = ASTNode_createWithType(childType);
    int err = ASTNode_prependChild(operateNode, newChild);
    if(err == 0) {
      printf("Error:ASTNode_prependChild\n");
      exit(0);
    }
    tmpVertexVector[i] = operateNode;
    ASTNode_removeChild(SubTree, i+2);
  }
  rep(i,childCnt-3) ASTNode_addChild(tmpVertexVector[i], tmpVertexVector[i+1]);
  ASTNode_addChild(tmpVertexVector[childCnt-2], SubTree);
  return tmpVertexVector[0];
}

double dfsOnASTBinaryTree(const ASTNode_t *v, KineticLaw_t *kl, int Splist_size,Model_t* m) {
	double left_val = 0, right_val = 0;
	ASTNodeType_t op = ASTNode_getType(v);
	if (ASTNode_getNumChildren(v) == 0) { // leaf
		if (op == AST_INTEGER)
			return ASTNode_getInteger(v);
		if (op == AST_REAL)
			return ASTNode_getReal(v);
		if (ASTNode_isName(v))
			return applyParameterAndVariables(v, kl, Splist_size, m);
	}
	ASTNode_t *left = ASTNode_getLeftChild(v);
	ASTNode_t *right = ASTNode_getRightChild(v);
  if((ASTNode_getNumChildren(left) > 2) && (ASTNode_getNumChildren(right) > 2)) {
    left_val = dfsOnASTBinaryTree(makeBinarySubtree(left),kl,Splist_size,m);
    right_val = dfsOnASTBinaryTree(makeBinarySubtree(right),kl,Splist_size,m);
    return solve(op,left_val,right_val);
  }
	else if (ASTNode_getNumChildren(left) > 2) {
		left_val = dfsOnASTBinaryTree(makeBinarySubtree(left),kl,Splist_size,m);
    right_val = dfsOnASTBinaryTree(right, kl, Splist_size,m);
    return solve(op,left_val,right_val);
	}
	else if (ASTNode_getNumChildren(right) > 2) {
		right_val = dfsOnASTBinaryTree(makeBinarySubtree(right),kl,Splist_size,m);
    dfsOnASTBinaryTree(left, kl, Splist_size,m);
    return solve(op,left_val,right_val);
	}
  left_val = dfsOnASTBinaryTree(left, kl, Splist_size,m);
	right_val = dfsOnASTBinaryTree(right, kl, Splist_size,m);
	return solve(op, left_val, right_val);
}

int main() {
	SBMLDocument_t *d;
	d = readSBML("Bier2000_s.xml");
	Model_t *m = SBMLDocument_getModel(d);
	ListOf_t *Sp_list = Model_getListOfSpecies(m);
	ListOf_t *Reac_list = Model_getListOfReactions(m);

	uint Splist_size = ListOf_size(Sp_list), Reac_size = ListOf_size(Reac_list);

	rep(i, Splist_size) {
		Sp[i] = (Species_t *)ListOf_get(Sp_list, i);
		Initial_Amount[i] = Species_getInitialAmount(Sp[i]);
		Ids[i] = Species_getId(Sp[i]);
	}
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
			ReacVals[i] = dfsOnASTBinaryTree(AST_NumericalFormula[i], KineticLaw[i], Splist_size, m);
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
