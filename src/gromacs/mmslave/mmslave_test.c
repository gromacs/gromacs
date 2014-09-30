#include <stdio.h>
#include <stdlib.h>
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/mmslave.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/types/hw_info.h"
#include "gromacs/legacyheaders/gmx_detect_hardware.h"
#include "gromacs/legacyheaders/gmx_omp_nthreads.h"

int main(int argc, char *argv[])
{
    struct t_commrec *cr;
    int natoms, ngroups;
    rvec *x, *v, *f, *A;
    real *phi;
    gmx_mmslave_t gms;
    double e0, e1, qq;
    int i;
    gmx_hw_info_t            *hwinfo       = NULL;

#ifdef GMX_LIB_MPI
    (void) MPI_Init(&argc, &argv);
#endif
    cr = init_commrec();
    
    /* Detect hardware, gather information. This is an operation that is
     * global for this process (MPI rank). */
    hwinfo = gmx_detect_hardware(stdout, cr, FALSE);
    
    gmx_omp_nthreads_init(stdout, cr,
                          hwinfo->nthreads_hw_avail,
                          0, /* hw_opt->nthreads_omp */
                          0, /* hw_opt->nthreads_omp_pme */
                          FALSE, /* (cr->duty & DUTY_PP) == 0 */
                          FALSE /* inputrec->cutoff_scheme == ecutsVERLET*/
                          );
    
    gms = mmslave_init(cr);
    if (argc > 1) 
    {
        gmx_bool bOK;
        
        bOK = mmslave_read_tpr(argv[1], gms);
        if (bOK)
        {
            ngroups = mmslave_ngroups(gms);
            natoms = mmslave_natoms(gms);
            printf("There are %d atoms in %d groups\n", natoms, ngroups);
            x = (rvec *)calloc(natoms, sizeof(rvec));
            v = (rvec *)calloc(natoms, sizeof(rvec));
            f = (rvec *)calloc(natoms, sizeof(rvec));
            A = (rvec *)calloc(natoms, sizeof(rvec));
            phi = (real *)calloc(natoms, sizeof(real));

            bOK = mmslave_copyX(gms, natoms, x);
            if ( 0 ) 
            {
                for(i=0; (i<natoms); i++)
                {
                    printf("Atom %5d Atomnumber %3d QM group %d\n",
                           i, mmslave_get_atomnumber(gms, i),
                           mmslave_get_group_id(gms, i));
                }
            }
        }
        else
        {
            printf("Could not read tpr file %s\n", argv[1]);
        }
        
        if (bOK)
        {
            bOK = mmslave_calc_energy(gms, stdout, (const rvec *)x, f, A, phi, &e0);
        }
        else
        {
            printf("Could not calculate energy\n");
        }
        
        if (bOK)
        {
            printf("The energy is %lf\n", e0);
            for(i=0; (i<3); i++)
            {
                printf("f[%d] = %10g  %10g  %10g\n", i, f[i][XX], f[i][YY], f[i][ZZ]);
                printf("A[%d] = %10g  %10g  %10g\n", i, A[i][XX], A[i][YY], A[i][ZZ]);
                printf("phi[%d] = %10g\n", i, phi[i]);
            }
            printf("Changing a coordinate and re-evaluating energy:\n");
            x[0][0] += 0.01;
            bOK = mmslave_calc_energy(gms, stdout, (const rvec *)x, f, A, phi, &e1);
        }
        else
        {
            printf("Could not calculate energy\n");
        }
        
        if (bOK)
        {
            printf("The energy is %lf\n", e1);
            for(i=0; (i<3); i++)
            {
                printf("f[%d] = %10g  %10g  %10g\n", i, f[i][XX], f[i][YY], f[i][ZZ]);
                printf("A[%d] = %10g  %10g  %10g\n", i, A[i][XX], A[i][YY], A[i][ZZ]);
                printf("phi[%d] = %10g\n", i, phi[i]);
            }
            bOK = mmslave_get_q(gms, 1, &qq);
        }
        
        if (bOK) 
        {
            printf("Changing an atom charge and re-evaluating energy:\n");
            bOK = mmslave_set_q(gms, 1, qq*2);
        }
        if (bOK)
        {
            bOK = mmslave_calc_energy(gms, stdout, (const rvec *)x, f, A, phi, &e0);
        }
        else
        {
            printf("Could not calculate energy\n");
        }
        
        if (bOK)
        {
            printf("The energy is %lf\n", e0);
        }
    }
    else
    {
        printf("Usage: %s file.tpr\n", argv[0]);
    }
        
    mmslave_done(gms);

#ifdef GMX_LIB_MPI
    MPI_Finalize();
#endif    
    return 0;
}
