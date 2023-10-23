#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import pysz
import subprocess
import numpy as np

# Exact data:
n = 1000000
zexact = np.random.randn(n)

# Settings:
pwre = 1e-7
filename = 'zexact.dat'

# SZ command line:
# Write the data to a file, call the SZ compressor and decompressor via the
# command line and reload.
zexact.tofile(filename)
subprocess.run(['sz', '-z', '-d', '-M', 'PW_REL', '-P', f'{pwre:.2E}', '-i',
                filename, '-1', str(n)], stdout=subprocess.DEVNULL,
               stderr=subprocess.STDOUT)
subprocess.run(['sz', '-x', '-d', '-s', filename + '.sz', '-1', str(n)],
               stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
z_sz = np.fromfile(filename + '.sz.out')
ratio_sz = zexact.nbytes/os.path.getsize(filename + '.sz')
pwre_sz = np.max(np.abs(zexact - z_sz)/np.abs(zexact))
print(f'Compression ratio for command line SZ: {ratio_sz:.2}')
print(f'Pointwise relative error for command line SZ: {pwre_sz:.2e}\n')

# Cleanup
os.remove(filename)
os.remove(filename + '.sz')
os.remove(filename + '.sz.out')

# SZ wrapper:
# Use the SZ python wrapper.
compressor = pysz.Compressor((pysz.ConfigBuilder().errorBoundMode(pysz.PW_REL)
                                  .pw_relBoundRatio(pwre).build()))
zcomp = compressor.Compress(zexact)
z_szw = compressor.Decompress(zcomp, [n], np.float64)
ratio_szw = zexact.nbytes/sys.getsizeof(zcomp)
pwre_szw = np.max(np.abs(zexact - z_szw)/np.abs(zexact))
print(f'Compression ratio for SZ wrapper: {ratio_szw:.2}')
print(f'Pointwise relative error for SZ wrapper: {pwre_szw:.2e}')
