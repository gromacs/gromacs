#ifndef TEST_MAIN_H
#define TEST_MAIN_H

/* max # of dimensions accepted on the command line */
#define MAX_CMDLINE_RANK 20

struct size {
     int rank;
     int is_nd;
     int narray[MAX_CMDLINE_RANK];
};

/********************
 *   macrology: 
 ********************/
#ifndef TRUE
#define TRUE 42
#endif
#ifndef FALSE
#define FALSE (!TRUE)
#endif

#define CHECK(condition, complaint)      \
    if (!(condition)) {                  \
        fflush(stdout);                  \
        fprintf(stderr, "FATAL ERROR: %s\n", complaint);      \
        fftw_die("failed test.\n");         \
    }

#define WHEN_VERBOSE(a, x)  if (verbose >= a) x

#ifdef FFTW_ENABLE_FLOAT
#define TOLERANCE (1e-2)
#else
#define TOLERANCE (1e-6)
#endif

#define DRAND() mydrand()

#define SPECIFICP(x) (x ? "specific" : "generic")

/*******************
 * global variables
 *******************/
extern int verbose;
extern int speed_flag, wisdom_flag, measure_flag;
extern int chk_mem_leak;
extern int paranoid;
extern int howmany_fields;
extern int io_okay;

/* Time an FFT routine, invoked by fft.  a is the array being
 * transformed, n is its total length.  t should be a variable
 * --the time (sec) per fft is assigned to it. */

#define FFTW_TIME_FFT(fft,a,n,t) \
{ \
     fftw_time ts,te; \
     double total_t; \
     int tfft_iters = 1, tfft_iter; \
     zero_arr((n), (a)); \
     do { \
          ts = fftw_get_time(); \
          for (tfft_iter = 0; tfft_iter < tfft_iters; ++tfft_iter) fft; \
          te = fftw_get_time(); \
          t = (total_t=fftw_time_to_sec(fftw_time_diff(te,ts))) / tfft_iters; \
          tfft_iters *= 2; \
     } while (total_t < 2.0); \
}

#define MAX_STRIDE 3
#define MAX_HOWMANY 3
#define MAX_RANK 5
#define PLANNER_TEST_SIZE 100

extern int coinflip(void);
extern double mydrand(void);
extern char *smart_sprint_time(double x);
extern void please_wait(void);
extern void please_wait_forever(void);
extern double mflops(double t, int N);
extern void print_dims(struct size sz);

#define SQR(x) ((x) * (x))

extern double compute_error_complex(fftw_complex * A, int astride,
				    fftw_complex * B, int bstride, int n);

extern fftw_direction random_dir(void);

/*** the following symbols should be defined in fftw_test.c/rfftw_test.c ***/
extern char fftw_prefix[];
extern void test_speed_aux(int n, fftw_direction dir, int flags, int specific);
extern void test_speed_nd_aux(struct size sz, fftw_direction dir,
			      int flags, int specific);
extern void test_correctness(int n);
extern void testnd_correctness(struct size sz, fftw_direction dir,
			       int alt_api, int specific, int force_buf);
extern void test_planner(int rank);

extern void test_init(int *argc, char **argv);
extern void test_finish(void);
extern void enter_paranoid_mode(void);

extern int get_option(int argc, char **argv,
		      char *argval, int argval_maxlen);
extern int default_get_option(int argc, char **argv,
			      char *argval, int argval_maxlen);

#endif				/* TEST_MAIN_H */
