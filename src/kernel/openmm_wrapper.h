#ifndef _OPENMM_WRAPPER_H_
#define _OPENMM_WRAPPER_H_

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef GMX_OPENMM
void* openmm_init(FILE *fplog, const char *platformOptStr,
                    t_commrec *cr,t_inputrec *ir,
                    gmx_mtop_t *top_global, gmx_localtop_t *top,
                    t_mdatoms *mdatoms, t_forcerec *fr,t_state *state);

void openmm_take_one_step(void* data);

void openmm_copy_state(void *data,
                        t_state *state, double *time,
                        rvec f[], gmx_enerdata_t *enerd,
                        bool includePos, bool includeVel, bool includeForce, bool includeEnergy);

void openmm_cleanup(FILE *fplog, void* data);
#else 
/* dummy versions of the openmm wrapper functions to enable compilation of 
   do_md_openmm even when OpenMM is not used */ 
void* openmm_init(FILE *fplog, const char *platformOptStr,
                    t_commrec *cr,t_inputrec *ir,
                    gmx_mtop_t *top_global, gmx_localtop_t *top,
                    t_mdatoms *mdatoms, t_forcerec *fr,t_state *state){return NULL;}

void openmm_take_one_step(void* data){}

void openmm_copy_state(void *data,
                        t_state *state, double *time,
                        rvec f[], gmx_enerdata_t *enerd,
                        bool includePos, bool includeVel, bool includeForce, bool includeEnergy){}

void openmm_cleanup(FILE *fplog, void* data){}

#endif //GMX_OPENMM


#ifdef __cplusplus
} // extern "C"
#endif

#endif /* _OPENMM_WRAPPER_H_ */

