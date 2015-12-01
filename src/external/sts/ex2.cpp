#include <cmath>
#include "sts.h"

const int niters = 10000000;
float A[niters];
float B[niters/3];
float C[niters/3];
float D[niters/3];

void do_something_A(const char* s, int i, int step) {
    // fprintf(stderr, "%s: i=%d step=%d tid=%d\n", s, i, step, Thread::getId());
    A[i] = sinf(i);
}

void do_something_B(const char* s, int i, int step) {
    // fprintf(stderr, "%s: i=%d step=%d tid=%d\n", s, i, step, Thread::getId());
    B[i] = sinf(i);
}

void do_something_C(const char* s, int i, int step) {
    // fprintf(stderr, "%s: i=%d step=%d tid=%d\n", s, i, step, Thread::getId());
    C[i] = sinf(i);
}

void do_something_D(const char* s, int i, int step) {
    // fprintf(stderr, "%s: i=%d step=%d tid=%d\n", s, i, step, Thread::getId());
    D[i] = sinf(i);
}

void f(int step) {
    // fprintf(stderr, "F: step=%d tid=%d\n", step, Thread::getId());

    parallel_for("TASK_F_0", 0, niters, [=](size_t i) {do_something_A("F0", i, step);});
}

void g(int step) {
    // fprintf(stderr, "G: step=%d tid=%d\n", step, Thread::getId());

    parallel_for("TASK_G_0", 0, niters/3, [=](size_t i) {do_something_B("G0", i, step);});

    for(int i=0; i<niters/3; i++) {do_something_C("G1", i, step);}

    parallel_for("TASK_G_1", 0, niters/3, [=](size_t i) {do_something_D("G2", i, step);});
}

int main(int argc, char **argv)
{
  const int nthreads = 3;
  const int nsteps = 3;

  setNumThreads(nthreads);

  for (int step=0; step<nsteps; step++)
  {
      /*
      if(step==2) 
          reschedule(); //can be done every step if desired
      if(step==3) 
          nextStep();
      */
      reschedule();
      run("TASK_F", [=]{f(step);});
      run("TASK_G", [=]{g(step);});
      wait();
      printf("%f\n", A[niters/4] + B[niters/4] + C[niters/4] + D[niters/4]);
  }
}
