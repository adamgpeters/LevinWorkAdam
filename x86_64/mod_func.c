#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _ichanR859C1_reg(void);
extern void _ichanWT2005_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," ichanR859C1.mod");
    fprintf(stderr," ichanWT2005.mod");
    fprintf(stderr, "\n");
  }
  _ichanR859C1_reg();
  _ichanWT2005_reg();
}
