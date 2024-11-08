#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _Traub_reg();
extern void _cad_reg();
extern void _cagk_reg();
extern void _cal_reg();
extern void _calH_reg();
extern void _car_reg();
extern void _cat_reg();
extern void _dCaAP_reg();
extern void _h_reg();
extern void _kca_reg();
extern void _myAMPA_reg();
extern void _myGABA_reg();
extern void _myNMDA_reg();
extern void _mySyn_reg();
extern void _na_reg();
extern void _na3_reg();
extern void _namir_reg();
extern void _nax_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," Traub.mod");
fprintf(stderr," cad.mod");
fprintf(stderr," cagk.mod");
fprintf(stderr," cal.mod");
fprintf(stderr," calH.mod");
fprintf(stderr," car.mod");
fprintf(stderr," cat.mod");
fprintf(stderr," dCaAP.mod");
fprintf(stderr," h.mod");
fprintf(stderr," kca.mod");
fprintf(stderr," myAMPA.mod");
fprintf(stderr," myGABA.mod");
fprintf(stderr," myNMDA.mod");
fprintf(stderr," mySyn.mod");
fprintf(stderr," na.mod");
fprintf(stderr," na3.mod");
fprintf(stderr," namir.mod");
fprintf(stderr," nax.mod");
fprintf(stderr, "\n");
    }
_Traub_reg();
_cad_reg();
_cagk_reg();
_cal_reg();
_calH_reg();
_car_reg();
_cat_reg();
_dCaAP_reg();
_h_reg();
_kca_reg();
_myAMPA_reg();
_myGABA_reg();
_myNMDA_reg();
_mySyn_reg();
_na_reg();
_na3_reg();
_namir_reg();
_nax_reg();
}
