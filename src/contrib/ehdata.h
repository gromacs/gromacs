extern real get_omega(real ekin,int *seed,FILE *fp);

extern real get_q_inel(real ekin,real omega,int *seed,FILE *fp);

extern real get_theta_el(real ekin,int *seed,FILE *fp);

extern real cross_inel(real ekin);

extern real cross_el(real ekin);

extern real band_ener(int *seed,FILE *fp);

extern void test_tables(int *seed,char *fn);

