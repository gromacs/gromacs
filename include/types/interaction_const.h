#ifndef _INTERACTION_CONST_
#define _INTERACTION_CONST_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    /* VdW */
    real rvdw;

    /* Cut-off */
    real rlist;

    /* PME/Ewald */
    real ewaldcoeff;
  
    /* Dielectric constant resp. multiplication factor for charges */
    real epsilon_r,epsilon_rf,epsfac;  
  
    /* Constants for reaction fields */
    real k_rf,c_rf;

    /* type of electrostatics (defined in enums.h) */
    int  eeltype;

} interaction_const_t;

#ifdef __cplusplus
}
#endif

#endif /* _INTERACTION_CONST_ */
