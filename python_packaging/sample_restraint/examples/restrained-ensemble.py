#!/usr/bin/env python
"""Run restrained-ensemble sampling and biasing workflow.

Irrgang, M. E., Hays, J. M., & Kasson, P. M.
gmxapi: a high-level interface for advanced control and extension of molecular dynamics simulations.
*Bioinformatics* 2018.
DOI: `10.1093/bioinformatics/bty484 <https://doi.org/10.1093/bioinformatics/bty484>`_
"""

# Restrained-ensemble formalism is a variant of that defined by Roux et al., 2013

import os
import sys

# Of course, the location of the Python plugin module is user-specific and could be
# passed by PYTHONPATH instead of programatically here.
sys.path.append("/home/mei2n/sample_restraint/build/src/pythonmodule")

import gmx

import logging

logging.getLogger().setLevel(logging.DEBUG)
# create console handler
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
# create formatter and add it to the handler
formatter = logging.Formatter("%(asctime)s:%(name)s:%(levelname)s: %(message)s")
ch.setFormatter(formatter)
# add the handlers to the logger
logging.getLogger().addHandler(ch)
logger = logging.getLogger()

import myplugin

logger.info("myplugin is {}".format(myplugin.__file__))

if len(sys.argv) > 1:
    size = int(sys.argv[1])
else:
    size = 20
input_dir_list = ["aa_{:02d}".format(i) for i in range(size)]
print("Input directory list: {}".format(input_dir_list))

tpr_list = [
    os.path.abspath(os.path.join(directory, "mRMR.tpr")) for directory in input_dir_list
]

# dt = 0.002
# First restraint applied between atoms 387 and 2569
# Second restraint applied between atom 1330 and 2520
# Restraint site coordinates relative to atom 1735
# Gathers 50 distance samples over 10ps, then averages the histogram across the ensemble to
# get a smooth histogram for the sample window. At each update (10 ps), updates the bias
# potential with the average statistics from the last 20 windows.
params = {
    "sites": [387, 1735, 2569],
    "k": 100.0,
    "sigma": 0.2,
    "nbins": 70,
    "binWidth": 0.1,
    "max_dist": 6.0,
    "min_dist": 1.9,
    "experimental": [
        1.799741371805743e-21,
        1.394386099050501e-19,
        8.502972718446353e-18,
        4.085581973134053e-16,
        1.548764419813057e-14,
        4.638792711370235e-13,
        1.0996581066864261e-11,
        2.0673577567003664e-10,
        3.0896196375481035e-09,
        3.680765701683053e-08,
        3.5071251461730994e-07,
        2.683161946307409e-06,
        1.6559233749685584e-05,
        8.289071953350906e-05,
        0.00033870482482321125,
        0.0011381605345541928,
        0.0031720369601255603,
        0.007403098031042665,
        0.014627984136430199,
        0.024781412546281113,
        0.036538520471987135,
        0.04776926450937005,
        0.056703917249967935,
        0.06290983130284482,
        0.06723680313442071,
        0.07080726976949929,
        0.07402160052652465,
        0.0765402296381977,
        0.07828089499810763,
        0.08047585508402359,
        0.08570745927698609,
        0.09674399700081715,
        0.11510583738691207,
        0.14051216332953817,
        0.17140325711911122,
        0.20554101371952832,
        0.23970384865306096,
        0.2690817283268939,
        0.2883208679737597,
        0.2947929622276882,
        0.2912919631495654,
        0.28445334990710786,
        0.27916526960634136,
        0.2740440567694397,
        0.2627885645671059,
        0.24084583312911054,
        0.2119920710658402,
        0.18961266641719554,
        0.19160148241368294,
        0.23183520488057596,
        0.31258390292452404,
        0.4212625964074282,
        0.5329315503397933,
        0.6177862653849041,
        0.6510130306580447,
        0.621354685844679,
        0.5350130330039692,
        0.4131558700737114,
        0.28383821525616004,
        0.17174126508493523,
        0.0904869618596048,
        0.04102083542122555,
        0.0158113507543527,
        0.00512381828028717,
        0.0013817258262933331,
        0.00030725601604590235,
        5.5898129564429734e-05,
        8.270145724798172e-06,
        1.066958972950409e-06,
        8.525674649177577e-07,
    ],
    "nsamples": 5,  # window size: 100 ps
    "sample_period": 10000 * 0.002,  # 20 ps
    "nwindows": 100,  # averaging period: 10 ns
}

potential1 = gmx.workflow.WorkElement(
    namespace="myplugin", operation="ensemble_restraint", depends=[], params=params
)
potential1.name = "ensemble_restraint_1"

params["sites"] = [1330, 1735, 2520]
params["experimental"] = [
    8.750538172089207e-20,
    4.963054541010076e-18,
    2.2585602895138136e-16,
    8.296295971141421e-15,
    2.4727072025999945e-13,
    6.001704592322284e-12,
    1.188164346191529e-10,
    1.917922517069427e-09,
    2.52033705254017e-08,
    2.691327568605142e-07,
    2.3326113123212396e-06,
    1.641120092881496e-05,
    9.38976910432654e-05,
    0.0004385427322814972,
    0.0016815454179692623,
    0.005334700709331302,
    0.014137960662648396,
    0.03164990069467528,
    0.06058416858044986,
    0.10044829089636065,
    0.1462620702334438,
    0.19008247987025195,
    0.22511473247247643,
    0.2496090827957139,
    0.2672559820151655,
    0.28375718473756745,
    0.3027960501795795,
    0.3245600144118122,
    0.346907132615181,
    0.3674499825469857,
    0.3851742255922683,
    0.40052661472987944,
    0.4135408960076755,
    0.4218337765359014,
    0.42147736029942173,
    0.4107213639593705,
    0.39238589395514606,
    0.37160127657310765,
    0.3512844290946452,
    0.3307655537130543,
    0.30883452118638083,
    0.2868146457416455,
    0.2677676912313197,
    0.2532530719313838,
    0.2417481938021803,
    0.23018156335307405,
    0.21644933221299598,
    0.20050997978961216,
    0.18377185607646654,
    0.16786393271557773,
    0.15362016010974153,
    0.1405238761751657,
    0.12680145743115132,
    0.11034786248372368,
    0.09017088591908493,
    0.06741463593670398,
    0.045069844741355614,
    0.026422842565654883,
    0.013359579997979203,
    0.005742341678209477,
    0.0020723966306519146,
    0.0006212767709345187,
    0.00015329363473072974,
    3.088693304526221e-05,
    5.048219168357052e-06,
    6.655228160961857e-07,
    7.04490968489331e-08,
    6.246463184699169e-09,
    4.5050655561489735e-09,
    4.676670950356501e-08,
]

potential2 = gmx.workflow.WorkElement(
    namespace="myplugin", operation="ensemble_restraint", depends=[], params=params
)
potential2.name = "ensemble_restraint_2"


# Settings for a 20 core HPC node. Use 18 threads for domain decomposition for pair potentials
# and the remaining 2 threads for PME electrostatics.
md = gmx.workflow.from_tpr(
    tpr_list, tmpi=20, grid=[3, 3, 2], ntomp_pme=1, npme=2, ntomp=1
)
md.add_dependency(potential1)
md.add_dependency(potential2)

context = gmx.context.ParallelArrayContext(md)

with context as session:
    session.run()
