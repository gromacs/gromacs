/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef DEBUG
#define DEBUG_BOND(log, x, ai, aj, dx, dr2, dr, bondparams, delta_r, vbond, fbond, fij)
#define DEBUG_ANGLE_UPD(log, f_i, f_j, f_k)
#define DEBUG_ANGLES(log, x, ai, aj, ak, theta, bondparams, dtheta, fijk, va)
#define DEBUG_ANGLE4S(log, x, ai, aj, ak, theta, bondparams, dth1, dth2, dth12, fijk, va)
#define DEBUG_DIH_ANGLE(log, xi, xj, xk, xl, r_ij, r_kj, r_kl, m, n, cos_phi, phi, ipr, \
                        sign)
#define DEBUG_DIH_FUP(log, i, j, k, l, ddphi, ipr, p, q, r_ij, r_kj, r_kl, m, n, f_i, f_j, \
                      f_k, f_l, uvec, vvec, svec)
#define DEBUG_PDIHS(log, ai, aj, ak, al, x, r_ij, r_kj, r_kl, m, n, ph0, cp, mult, mdphi, \
                    ddphi, vpd)
#define DEBUG_IDIHS(log, ai, aj, ak, al, x, ph0, cp, dphi, ddphi, vid)
#define DEBUG_RBDIHS(log, rbc, ddphi, v)
#define DEBUG_POSRES(log, ai, x, pos0, fc, dx, dr2, fi, v)
#else
#define vlist(v)    (v)[0], (v)[1], (v)[2]
#define ivlist(i, v) i, (v)[i][0], (v)[i][1], (v)[i][2]
#define VECFMT(v)   #v "=(%e,%e,%e) "
#define IVECFMT(i, v)    #i "=%d (%e,%e,%e) "
#define INTFMT(i)   #i "=%d "
#define REALFMT(r)  #r "=%e "
#define EOLNFMT     "\n",

#define DEBUG_BOND(log, x, ai, aj, dx, dr2, dr, bondparams, delta_r, vbond, fbond, fij) \
    fprintf(log, "bonds: " \
            IVECFMT(ai, x) \
            IVECFMT(aj, x) \
            VECFMT (dx) \
            REALFMT(dr2) \
            REALFMT(dr) \
            REALFMT(b0) \
            REALFMT(cb) \
            REALFMT(delta_r) \
            REALFMT(vbond) \
            REALFMT(fbond) \
            VECFMT (fij) \
            EOLNFMT \
            ivlist(i, x), ivlist(aj, x), \
            vlist(dx), dr2, dr, bondparams->bonds.b0, bondparams->bonds.cb, \
            delta_r, vbond, fbond, vlist(fij))
#define DEBUG_ANGLE_UPD(log, f_i, f_j, f_k) \
    fprintf(log, "do_fang_upd: " \
            VECFMT (f_i) \
            VECFMT (f_j) \
            VECFMT (f_k) \
            EOLNFMT \
            vlist(f_i), vlist(f_j), vlist(f_k))
#define DEBUG_ANGLES(log, x, ai, aj, ak, theta, bondparams, dtheta, fijk, va) \
    fprintf(log, "angles: " \
            IVECFMT(ai, x) \
            IVECFMT(aj, x) \
            IVECFMT(ak, x) \
            REALFMT(theta) \
            REALFMT(th0) \
            REALFMT(ct) \
            REALFMT(theta) \
            REALFMT(fijk) \
            REALFMT(va) \
            EOLNFMT \
            ivlist(ai, x), ivlist(aj, x), ivlist(ak, x), theta, \
            bondparams->angles.th0, bondparams->angles.ct, dtheta, fijk, va)
#define DEBUG_ANGLE4S(log, x, ai, aj, ak, theta, bondparams, dth1, dth2, dth12, fijk, v) \
    fprintf(log, "angle4s: " \
            IVECFMT(ai, x) \
            IVECFMT(aj, x) \
            IVECFMT(ak, x) \
            REALFMT(theta) \
            REALFMT(th0) \
            REALFMT(ct) \
            REALFMT(dth1) \
            REALFMT(dth2) \
            REALFMT(dth12) \
            REALFMT(fijk) \
            REALFMT(va) \
            EOLNFMT \
            ai, vlist(x[ai]), aj, vlist(x[aj]), ak, vlist(x[ak]), theta, \
            bondparams->angles.th0, bondparams->angles.ct, dth1, dth2, dth12, fijk, va)
#define DEBUG_DIH_ANGLE(log, xi, xj, xk, xl, r_ij, r_kj, r_kl, m, n, cos_phi, phi, ipr, \
                        sign) \
    fprintf(log, "dih_angle: " \
            VECFMT (xi) \
            VECFMT (xj) \
            VECFMT (xk) \
            VECFMT (xl) \
            VECFMT (r_ij) \
            VECFMT (r_kj) \
            VECFMT (r_kl) \
            VECFMT (m) \
            VECFMT (n) \
            REALFMT(cos_phi) \
            REALFMT(phi) \
            REALFMT(ipr) \
            EOLNFMT \
            vlist(xi), vlist(xj), vlist(xk), vlist(xl), \
            vlist(r_ij), vlist(r_kj), vlist(r_kl), vlist(m), vlist(n), \
            *(cos_phi), phi, ipr, *(sign))
#define DEBUG_DIH_FUP(log, i, j, k, l, ddphi, ipr, p, q, r_ij, r_kj, r_kl, m, n, f_i, f_j, \
                      f_k, f_l, uvec, vvec, svec) \
    fprintf(log, "dih_fup: " \
            INTFMT (i) \
            INTFMT (j) \
            INTFMT (k) \
            INTFMT (l) \
            REALFMT(ddphi) \
            REALFMT(ipr) \
            REALFMT(p) \
            REALFMT(q) \
            VECFMT (r_ij) \
            VECFMT (r_kj) \
            VECFMT (r_kl) \
            VECFMT (m) \
            VECFMT (n) \
            VECFMT (f_i) \
            VECFMT (f_j) \
            VECFMT (f_k) \
            VECFMT (f_l) \
            VECFMT (uvec) \
            VECFMT (vvec) \
            VECFMT (svec) \
            EOLNFMT \
            i, j, k, l, ddphi, ipr, p, q, vlist(r_ij), vlist(r_kj), vlist(r_kl), \
            vlist(m), vlist(n), \
            vlist(f_i), vlist(f_j), vlist(f_k), vlist(f_l), \
            vlist(uvec), vlist(vvec), vlist(svec))
#define DEBUG_PDIHS(log, ai, aj, ak, al, x, r_ij, r_kj, r_kl, m, n, ph0, cp, mult, mdphi, \
                    ddphi, vpd) \
    fprintf(log, "pdihs: " \
            IVECFMT(ai, x) \
            IVECFMT(aj, x) \
            IVECFMT(ak, x) \
            IVECFMT(al, x) \
            VECFMT(r_ij) \
            VECFMT(r_kj) \
            VECFMT(r_kl) \
            VECFMT(m) \
            VECFMT(n) \
            REALFMT(ph0) \
            REALFMT(cp) \
            INTFMT(mult) \
            REALFMT(mdphi) \
            REALFMT(ddphi) \
            REALFMT(vpd) \
            EOLNFMT \
            ivlist(ai, x), ivlist(aj, x), ivlist(ak, x), ivlist(al, x), vlist(r_ij), \
            vlist(r_kj), vlist(r_kl), vlist(m), vlist(n), ph0, cp, mult, mdphi, ddphi, \
            vpd)
#define DEBUG_IDIHS(log, ai, aj, ak, al, x, ph0, cp, dphi, ddphi, vid) \
    fprintf(log, "idihs: " \
            IVECFMT(ai, x) \
            IVECFMT(aj, x) \
            IVECFMT(ak, x) \
            IVECFMT(al, x) \
            REALFMT(ph0) \
            REALFMT(cp) \
            REALFMT(dphi) \
            REALFMT(ddphi) \
            REALFMT(vid) \
            EOLNFMT \
            ivlist(ai, x), ivlist(aj, x), ivlist(ak, x), ivlist(al, x), ph0, cp, dphi, \
            ddphi, vid)
#define DEBUG_RBDIHS(log, rbc, ddphi, v) \
    fprintf(log, "rbdihs: " \
            REALFMT(rbc0) \
            REALFMT(rbc1) \
            REALFMT(rbc2) \
            REALFMT(rbc3) \
            REALFMT(rbc4) \
            REALFMT(rbc5) \
            REALFMT(ddphi) \
            REALFMT(v) \
            EOLNFMT \
                (rbc)[0], (rbc)[1], (rbc)[2], (rbc)[3], (rbc)[4], (rbc)[5], ddphi, v)
#define DEBUG_POSRES(log, ai, x, pos0, fc, dx, dr2, fi, v) \
    fprintf(log, "rbdihs: " \
            IVECFMT(ai, x) \
            VECFMT(pos0) \
            REALFMT(fc) \
            VECFMT(dx) \
            REALFMT(dr2) \
            VECFMT(fi) \
            REALFMT(v) \
            EOLNFMT \
            ivlist(ai, x), vlist(pos0), fc, vlist(dx), dr2, vlist(fi), v)
#endif
