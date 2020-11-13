#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _ichanR859C1_reg();
extern void _ichanWT2005_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," ichanR859C1.mod");
fprintf(stderr," ichanWT2005.mod");
fprintf(stderr, "\n");
    }
_ichanR859C1_reg();
_ichanWT2005_reg();
}
