/*
 * This source file is part of the Alexandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef MOLSELECT_HPP
#define MOLSELECT_HPP

typedef struct gmx_molselect *gmx_molselect_t;

enum iMolSelect { imsTrain, imsTest, imsIgnore, imsUnknown, imsNR };

extern const char *ims_names[imsNR];

gmx_molselect_t gmx_molselect_init(const char *fn);

void gmx_molselect_done(gmx_molselect_t gms);

iMolSelect gmx_molselect_status(gmx_molselect_t gms,const char *iupac);

int gmx_molselect_index(gmx_molselect_t gms,const char *iupac);

#endif
