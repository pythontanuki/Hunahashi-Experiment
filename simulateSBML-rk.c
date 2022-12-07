#include "sbml/SBMLTypes.h"
#include <stdio.h>
#define uint unsigned int

double f(double x, double y);

int main() {
	//(1)state, read SBML file
	// get species names, initial amount, kineticlaw formula
	puts("(1)state");
	SBMLDocument_t *d = readSBML("simple.xml");
	Model_t *m = SBMLDocument_getModel(d);
	ListOf_t *ListOfSpecies = Model_getListOfSpecies(m);
	uint SpeciesListOfSize = ListOf_size(ListOfSpecies);

	const char *SpeciesNameVector[SpeciesListOfSize];
	double SpeciesInitialAmountVector[SpeciesListOfSize];

	for (uint i = 0; i < SpeciesListOfSize; i++) {
		Species_t *s = (Species_t *)ListOf_get(ListOfSpecies, i);
		SpeciesNameVector[i] = Species_getName(s);
		SpeciesInitialAmountVector[i] = Species_getInitialAmount(s);
	}

	for (uint i = 0; i < SpeciesListOfSize; i++) {
		//(s1 1.0), (s2, 0.0)
		printf("Species Name : %s\n", SpeciesNameVector[i]);
		printf("Initial Amount : %.10lf\n", SpeciesInitialAmountVector[i]);
	}

	ListOf_t *ListOfReactions = Model_getListOfReactions(m);
	uint ReactionsListOfSize = ListOf_size(ListOfReactions);

	double k;
	const char *KineticLaw_Formulas[ReactionsListOfSize];

	for (uint i = 0; i < ReactionsListOfSize; i++) {
		Reaction_t *r = (Reaction_t *)ListOf_get(ListOfReactions, i);
		KineticLaw_t *kl = Reaction_getKineticLaw(r);
		ListOf_t *ListOfParameters = KineticLaw_getListOfParameters(kl);
		uint ParametersListOfSize = ListOf_size(ListOfParameters);
		for (uint j = 0; j < ParametersListOfSize; j++) {
			Parameter_t *p = KineticLaw_getParameter(kl, j);
			double ParameterValue = Parameter_getValue(p);
			k = ParameterValue;
			printf("Parameter Value : %.10lf\n", ParameterValue);
		}
		const char *kf = KineticLaw_getFormula(kl);
		KineticLaw_Formulas[i] = kf;
	}

	for (uint i = 0; i < ReactionsListOfSize; i++) {
		printf("KineticLaw_Fomulas : %s\n", KineticLaw_Formulas[i]);
	}
	puts("");
	puts("(2)state");
	//(2)state, integrate numbers
	//
	double s1 = SpeciesInitialAmountVector[0];
	// s2
	double s2 = SpeciesInitialAmountVector[1];
	// parameter k
	// parameter h --> δt
	double h = 0.1;
	// count of simulate
        // printf("%.10lf\n", k);
	double N = 100;
        FILE *fp = fopen("res.csv", "w");
	for (double i = 0; i <= N; i++) {
		// printf("%.10lf %.10lf\n", s1, s2);
            double time = i*h;  
            fprintf(fp, "%.10lf,%.10lf,%.10lf\n", time, s1, s2);
		double d1 = h * f(k, s1);
		double d2 = h * f(k, s1 - d1 / 2.0);
		double d3 = h * f(k, s1 - d2 / 2.0);
		double d4 = h * f(k, s1 - d3);
		double dif = (d1 + 2.0 * d2 + 2.0 * d3 + d4) / 6.0;
		s1 -= dif;
		s2 += dif;
	}
	SBMLDocument_free(d);
	return 0;
}

double f(double k, double x) {
	return k * x;
}
