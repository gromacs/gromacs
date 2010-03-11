#ifndef _powerspect_h
#define _powerspect_h

#include <gmx_fft.h>
#include "interf.h"

extern void powerspectavg(real ***interface, int t, int xbins, int ybins, char **outfiles);
extern void powerspectavg_intf(t_interf ***if1,t_interf ***if2,int t, int xbins, int ybins, char **outfiles);

#endif
