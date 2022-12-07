#include <stdio.h>
#include "sbml/SBMLTypes.h"

int main (void)
{
  SBMLDocument_t *d;
  int level, version;

  d = readSBML("mapk.xml");
  level = SBMLDocument_getLevel(d);
  version = SBMLDocument_getVersion(d);

  printf("Level %d, Version %d\n", level, version);
  SBMLDocument_free(d);
  return 0;
}
