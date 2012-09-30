/*
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
 *
 * Gromacs is a library for molecular simulation and trajectory analysis,
 * written by Erik Lindahl, David van der Spoel, Berk Hess, and others - for
 * a full list of developers and information, check out http://www.gromacs.org
 *
 * This program is free software; you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the Free 
 * Software Foundation; either version 2 of the License, or (at your option) any 
 * later version.
 * As a special exception, you may use this file as part of a free software
 * library without restriction.  Specifically, if other files instantiate
 * templates or use macros or inline functions from this file, or you compile
 * this file and link it with other files to produce an executable, this
 * file does not by itself cause the resulting executable to be covered by
 * the GNU Lesser General Public License.  
 *
 * In plain-speak: do not worry about classes/macros/templates either - only
 * changes to the library have to be LGPL, not an application linking with it.
 *
 * To help fund GROMACS development, we humbly ask that you cite
 * the papers people have written on it - you can find them on the website!
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>

#include "types/nrnb.h"
#include "nb_kernel_c.h"
#include "../nb_kernel.h"


/* Include standard kernel headers in local directory */
#include "nb_kernel010.h"
#include "nb_kernel020.h"
#include "nb_kernel030.h"
#include "nb_kernel100.h"
#include "nb_kernel101.h"
#include "nb_kernel102.h"
#include "nb_kernel103.h"
#include "nb_kernel104.h"
#include "nb_kernel110.h"
#include "nb_kernel111.h"
#include "nb_kernel112.h"
#include "nb_kernel113.h"
#include "nb_kernel114.h"
#include "nb_kernel120.h"
#include "nb_kernel121.h"
#include "nb_kernel122.h"
#include "nb_kernel123.h"
#include "nb_kernel124.h"
#include "nb_kernel130.h"
#include "nb_kernel131.h"
#include "nb_kernel132.h"
#include "nb_kernel133.h"
#include "nb_kernel134.h"
#include "nb_kernel200.h"
#include "nb_kernel201.h"
#include "nb_kernel202.h"
#include "nb_kernel203.h"
#include "nb_kernel204.h"
#include "nb_kernel210.h"
#include "nb_kernel211.h"
#include "nb_kernel212.h"
#include "nb_kernel213.h"
#include "nb_kernel214.h"
#include "nb_kernel220.h"
#include "nb_kernel221.h"
#include "nb_kernel222.h"
#include "nb_kernel223.h"
#include "nb_kernel224.h"
#include "nb_kernel230.h"
#include "nb_kernel231.h"
#include "nb_kernel232.h"
#include "nb_kernel233.h"
#include "nb_kernel234.h"
#include "nb_kernel300.h"
#include "nb_kernel301.h"
#include "nb_kernel302.h"
#include "nb_kernel303.h"
#include "nb_kernel304.h"
#include "nb_kernel310.h"
#include "nb_kernel311.h"
#include "nb_kernel312.h"
#include "nb_kernel313.h"
#include "nb_kernel314.h"
#include "nb_kernel320.h"
#include "nb_kernel321.h"
#include "nb_kernel322.h"
#include "nb_kernel323.h"
#include "nb_kernel324.h"
#include "nb_kernel330.h"
#include "nb_kernel331.h"
#include "nb_kernel332.h"
#include "nb_kernel333.h"
#include "nb_kernel334.h"
#include "nb_kernel400.h"
#include "nb_kernel410.h"
#include "nb_kernel420.h"
#include "nb_kernel430.h"


nb_kernel_info_t
kernellist_c[] =
{
    {nb_kernel010, "nb_kernel010", "C", "None","LennardJones","ParticleParticle","", "VF" },
    {nb_kernel020, "nb_kernel020", "C", "None","Buckingham","ParticleParticle","", "VF" },
    {nb_kernel030, "nb_kernel030", "C", "None","CubicSplineTable","ParticleParticle","", "VF" },
    {nb_kernel100, "nb_kernel100", "C", "Coulomb","None","ParticleParticle","", "VF" },
    {nb_kernel101, "nb_kernel101", "C", "Coulomb","None","Water3Particle","", "VF" },
    {nb_kernel102, "nb_kernel102", "C", "Coulomb","None","Water3Water3","", "VF" },
    {nb_kernel103, "nb_kernel103", "C", "Coulomb","None","Water4Particle","", "VF" },
    {nb_kernel104, "nb_kernel104", "C", "Coulomb","None","Water4Water4","", "VF" },
    {nb_kernel110, "nb_kernel110", "C", "Coulomb","LennardJones","ParticleParticle","", "VF" },
    {nb_kernel111, "nb_kernel111", "C", "Coulomb","LennardJones","Water3Particle","", "VF" },
    {nb_kernel112, "nb_kernel112", "C", "Coulomb","LennardJones","Water3Water3","", "VF" },
    {nb_kernel113, "nb_kernel113", "C", "Coulomb","LennardJones","Water4Particle","", "VF" },
    {nb_kernel114, "nb_kernel114", "C", "Coulomb","LennardJones","Water4Water4","", "VF" },
    {nb_kernel120, "nb_kernel120", "C", "Coulomb","Buckingham","ParticleParticle","", "VF" },
    {nb_kernel121, "nb_kernel121", "C", "Coulomb","Buckingham","Water3Particle","", "VF" },
    {nb_kernel122, "nb_kernel122", "C", "Coulomb","Buckingham","Water3Water3","", "VF" },
    {nb_kernel123, "nb_kernel123", "C", "Coulomb","Buckingham","Water4Particle","", "VF" },
    {nb_kernel124, "nb_kernel124", "C", "Coulomb","Buckingham","Water4Water4","", "VF" },
    {nb_kernel130, "nb_kernel130", "C", "Coulomb","CubicSplineTable","ParticleParticle","", "VF" },
    {nb_kernel131, "nb_kernel131", "C", "Coulomb","CubicSplineTable","Water3Particle","", "VF" },
    {nb_kernel132, "nb_kernel132", "C", "Coulomb","CubicSplineTable","Water3Water3","", "VF" },
    {nb_kernel133, "nb_kernel133", "C", "Coulomb","CubicSplineTable","Water4Particle","", "VF" },
    {nb_kernel134, "nb_kernel134", "C", "Coulomb","CubicSplineTable","Water4Water4","", "VF" },
    {nb_kernel200, "nb_kernel200", "C", "ReactionField","None","ParticleParticle","", "VF" },
    {nb_kernel201, "nb_kernel201", "C", "ReactionField","None","Water3Particle","", "VF" },
    {nb_kernel202, "nb_kernel202", "C", "ReactionField","None","Water3Water3","", "VF" },
    {nb_kernel203, "nb_kernel203", "C", "ReactionField","None","Water4Particle","", "VF" },
    {nb_kernel204, "nb_kernel204", "C", "ReactionField","None","Water4Water4","", "VF" },
    {nb_kernel210, "nb_kernel210", "C", "ReactionField","LennardJones","ParticleParticle","", "VF" },
    {nb_kernel211, "nb_kernel211", "C", "ReactionField","LennardJones","Water3Particle","", "VF" },
    {nb_kernel212, "nb_kernel212", "C", "ReactionField","LennardJones","Water3Water3","", "VF" },
    {nb_kernel213, "nb_kernel213", "C", "ReactionField","LennardJones","Water4Particle","", "VF" },
    {nb_kernel214, "nb_kernel214", "C", "ReactionField","LennardJones","Water4Water4","", "VF" },
    {nb_kernel220, "nb_kernel220", "C", "ReactionField","Buckingham","ParticleParticle","", "VF" },
    {nb_kernel221, "nb_kernel221", "C", "ReactionField","Buckingham","Water3Particle","", "VF" },
    {nb_kernel222, "nb_kernel222", "C", "ReactionField","Buckingham","Water3Water3","", "VF" },
    {nb_kernel223, "nb_kernel223", "C", "ReactionField","Buckingham","Water4Particle","", "VF" },
    {nb_kernel224, "nb_kernel224", "C", "ReactionField","Buckingham","Water4Water4","", "VF" },
    {nb_kernel230, "nb_kernel230", "C", "ReactionField","CubicSplineTable","ParticleParticle","", "VF" },
    {nb_kernel231, "nb_kernel231", "C", "ReactionField","CubicSplineTable","Water3Particle","", "VF" },
    {nb_kernel232, "nb_kernel232", "C", "ReactionField","CubicSplineTable","Water3Water3","", "VF" },
    {nb_kernel233, "nb_kernel233", "C", "ReactionField","CubicSplineTable","GEOM","", "VF" },
    {nb_kernel234, "nb_kernel234", "C", "ReactionField","CubicSplineTable","Water4Water4","", "VF" },
    {nb_kernel300, "nb_kernel300", "C", "CubicSplineTable","None","ParticleParticle","", "VF" },
    {nb_kernel301, "nb_kernel301", "C", "CubicSplineTable","None","Water3Particle","", "VF" },
    {nb_kernel302, "nb_kernel302", "C", "CubicSplineTable","None","Water3Water3","", "VF" },
    {nb_kernel303, "nb_kernel303", "C", "CubicSplineTable","None","Water4Particle","", "VF" },
    {nb_kernel304, "nb_kernel304", "C", "CubicSplineTable","None","Water4Water4","", "VF" },
    {nb_kernel310, "nb_kernel310", "C", "CubicSplineTable","LennardJones","ParticleParticle","", "VF" },
    {nb_kernel311, "nb_kernel311", "C", "CubicSplineTable","LennardJones","Water3Particle","", "VF" },
    {nb_kernel312, "nb_kernel312", "C", "CubicSplineTable","LennardJones","Water3Water3","", "VF" },
    {nb_kernel313, "nb_kernel313", "C", "CubicSplineTable","LennardJones","Water4Particle","", "VF" },
    {nb_kernel314, "nb_kernel314", "C", "CubicSplineTable","LennardJones","Water4Water4","", "VF" },
    {nb_kernel320, "nb_kernel320", "C", "CubicSplineTable","Buckingham","ParticleParticle","", "VF" },
    {nb_kernel321, "nb_kernel321", "C", "CubicSplineTable","Buckingham","Water3Particle","", "VF" },
    {nb_kernel322, "nb_kernel322", "C", "CubicSplineTable","Buckingham","Water3Water3","", "VF" },
    {nb_kernel323, "nb_kernel323", "C", "CubicSplineTable","Buckingham","Water4Particle","", "VF" },
    {nb_kernel324, "nb_kernel324", "C", "CubicSplineTable","Buckingham","Water4Water4","", "VF" },
    {nb_kernel330, "nb_kernel330", "C", "CubicSplineTable","CubicSplineTable","ParticleParticle","", "VF" },
    {nb_kernel331, "nb_kernel331", "C", "CubicSplineTable","CubicSplineTable","Water3Particle","", "VF" },
    {nb_kernel332, "nb_kernel332", "C", "CubicSplineTable","CubicSplineTable","Water3Water3","", "VF" },
    {nb_kernel333, "nb_kernel333", "C", "CubicSplineTable","CubicSplineTable","Water4Particle","", "VF" },
    {nb_kernel334, "nb_kernel334", "C", "CubicSplineTable","CubicSplineTable","Water4Water4","", "VF" },

    {nb_kernel300, "nb_kernel300", "C", "Ewald","None","ParticleParticle","", "VF" },
    {nb_kernel301, "nb_kernel301", "C", "Ewald","None","Water3Particle","", "VF" },
    {nb_kernel302, "nb_kernel302", "C", "Ewald","None","Water3Water3","", "VF" },
    {nb_kernel303, "nb_kernel303", "C", "Ewald","None","Water4Particle","", "VF" },
    {nb_kernel304, "nb_kernel304", "C", "Ewald","None","Water4Water4","", "VF" },
    {nb_kernel310, "nb_kernel310", "C", "Ewald","LennardJones","ParticleParticle","", "VF" },
    {nb_kernel311, "nb_kernel311", "C", "Ewald","LennardJones","Water3Particle","", "VF" },
    {nb_kernel312, "nb_kernel312", "C", "Ewald","LennardJones","Water3Water3","", "VF" },
    {nb_kernel313, "nb_kernel313", "C", "Ewald","LennardJones","Water4Particle","", "VF" },
    {nb_kernel314, "nb_kernel314", "C", "Ewald","LennardJones","Water4Water4","", "VF" },
    {nb_kernel320, "nb_kernel320", "C", "Ewald","Buckingham","ParticleParticle","", "VF" },
    {nb_kernel321, "nb_kernel321", "C", "Ewald","Buckingham","Water3Particle","", "VF" },
    {nb_kernel322, "nb_kernel322", "C", "Ewald","Buckingham","Water3Water3","", "VF" },
    {nb_kernel323, "nb_kernel323", "C", "Ewald","Buckingham","Water4Particle","", "VF" },
    {nb_kernel324, "nb_kernel324", "C", "Ewald","Buckingham","Water4Water4","", "VF" },
    {nb_kernel330, "nb_kernel330", "C", "Ewald","CubicSplineTable","ParticleParticle","", "VF" },
    {nb_kernel331, "nb_kernel331", "C", "Ewald","CubicSplineTable","Water3Particle","", "VF" },
    {nb_kernel332, "nb_kernel332", "C", "Ewald","CubicSplineTable","Water3Water3","", "VF" },
    {nb_kernel333, "nb_kernel333", "C", "Ewald","CubicSplineTable","Water4Particle","", "VF" },
    {nb_kernel334, "nb_kernel334", "C", "Ewald","CubicSplineTable","Water4Water4","", "VF" },

    {nb_kernel400, "nb_kernel400", "C", "GeneralizedBorn","None","ParticleParticle","", "VF" },
    {nb_kernel410, "nb_kernel410", "C", "GeneralizedBorn","LennardJones","ParticleParticle","", "VF" },
    {nb_kernel420, "nb_kernel420", "C", "GeneralizedBorn","Buckingham","ParticleParticle","", "VF" },
    {nb_kernel430, "nb_kernel430", "C", "GeneralizedBorn","CubicSplineTable","ParticleParticle","", "VF" }
};

int
kernellist_c_size = sizeof(kernellist_c)/sizeof(kernellist_c[0]);

