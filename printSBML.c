#include "sbml/SBMLTypes.h"
#include <stdbool.h>
#include <stdio.h>
#define uint unsigned int
const char *MAPK = "MAPK";
const char *MAPK_PP = "MAPK_PP";
const char *J6 = "J6";
const char *KK7 = "KK7";

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

void printNameAndInitialAmount(Species_t *s,const char *id,const char *targetId) {
	if (compareId(id, targetId)) {
		const char *SpeciesName = Species_getName(s);
		double InitialAmount = Species_getInitialAmount(s);
		printf("%s Species Name: %s  Initial Amount %.10lf\n",targetId,SpeciesName,InitialAmount);
	}
}

void printParameterOfReaction(Model_t *m, const char *id, const char *targetId, const char *parameterId) {
	if (compareId(id, targetId)) {
		Reaction_t *r = Model_getReactionById(m, targetId);
		KineticLaw_t *kl = Reaction_getKineticLaw(r);
		Parameter_t *p = KineticLaw_getParameterById(kl, parameterId);
		double v = Parameter_getValue(p);
		printf("%s %.10lf\n", parameterId, v);
	}
}
int main(void) {
	//手順2
	SBMLDocument_t *d = readSBML("mapk.xml");
	uint SBML_Level = SBMLDocument_getLevel(d);
	uint SBML_version = SBMLDocument_getVersion(d);
	Model_t *m = SBMLDocument_getModel(d);
	uint SpeciesNum = Model_getNumSpecies(m);
	uint ReactionsNum = Model_getNumReactions(m);
	uint CompartmentsNum = Model_getNumCompartments(m);
	printf("Level = %d\n", SBML_Level);
	printf("version = %d\n", SBML_version);
	printf("num. of species = %d\n", SpeciesNum);
	printf("num. of reactions = %d\n", ReactionsNum);
	printf("num. of compartments = %d\n", CompartmentsNum);

	//手順3
	ListOf_t *ListOfSpecies = Model_getListOfSpecies(m);
	uint ListSize1 = ListOf_size(ListOfSpecies);
	for (uint i = 0; i < ListSize1; i++) {
		Species_t *s = (Species_t *)ListOf_get(ListOfSpecies, i);
		const char *id = Species_getId(s);
		printNameAndInitialAmount(s, id, MAPK);
		printNameAndInitialAmount(s, id, MAPK_PP);
	}

	ListOf_t *ListOfReactions = Model_getListOfReactions(m);
	uint ListSize2 = ListOf_size(ListOfReactions);
	for (uint i = 0; i < ListSize2; i++) {
		Reaction_t *r = (Reaction_t *)ListOf_get(ListOfReactions, i);
		const char *id = Reaction_getId(r);
		printParameterOfReaction(m, id, J6, KK7);
	}
  
	return 0;
}
