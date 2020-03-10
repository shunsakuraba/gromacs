/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
#include "gmxpre.h"

#include "zmm_util.h"

#include <cmath>
#include <cstdio>

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/pleasecite.h"

void calc_zmmfac_double(FILE *fplog, int eel, int zmm_degree, double zmm_alpha,
                        double Rc,
                        double *c2, double *c4, double *c6, double *c)
{
    /* Calculate the constants for Zero-multipole methods.
     * epsfac q_i q_j (erf(alpha r)/r + c2 r^2 + c4 r^4 + c6 r^6 + c)
     * note that the constant term c0 in paper and herein have a different sign.
     */
    if(eel == eelZMM)
    {
        real arc, d1, d2, d3;
        arc = zmm_alpha * Rc;

        switch(zmm_degree)
        {
        case 0:
            /* Zero-monopole (Wolf) method */
            *c = - std::erfc(arc) / Rc;
            *c2 = 0;
            *c4 = 0;
            *c6 = 0;
            if(fplog)
            {
                please_cite(fplog, "Wolf1999");
                fprintf(fplog, "Wolf method:\n"
                        "alpha = %g, rc = %g, u(Rc) = %g\n",
                        zmm_alpha, Rc, *c);
            }
            break;
        case 1:
            /* Zero-dipole method */
            d1 = 1. * pow(Rc, -2.0) * (std::erfc(arc) + 2. / sqrt(M_PI) * exp(-arc * arc) * arc);
            *c2 = 1. / 2. * d1 * pow(Rc, -1.0);
            *c =  - (std::erfc(arc) / Rc + *c2 * pow(Rc, 2.0));
            *c4 = 0;
            *c6 = 0;
            if(fplog)
            {
                please_cite(fplog, "Fukuda2011");
                fprintf(fplog, "%s (up to dipole):\n"
                        "alpha = %g, rc = %g, d1 = %g, c2 = %g, u(Rc) = %g\n",
                        eel_names[eel], zmm_alpha, Rc, d1, *c2, *c);
            }
            break;
        case 2:
            /* Zero-quadrupole method */
            d1 = 1. * pow(Rc, -2.0) * (std::erfc(arc) + 2. / sqrt(M_PI) * exp(-arc * arc) * arc);
            d2 = 2. * pow(Rc, -3.0) * (std::erfc(arc) + 2. / sqrt(M_PI) * exp(-arc * arc) * (arc + pow(arc, 3.0)));

            *c2 =   3. / 4. * d1 * pow(Rc, -1.0) + 1. / 4. * d2;
            *c4 = - 1. / 8. * d1 * pow(Rc, -3.0) - 1. / 8. * d2 * pow(Rc, -2.0);
            *c6 = 0.;

            *c = - (std::erfc(arc) / Rc + *c2 * pow(Rc, 2.0) + *c4 * pow(Rc, 4.0));

            if(fplog)
            {
                please_cite(fplog, "Fukuda2014");
                please_cite(fplog, "Sakuraba2018");
                fprintf(fplog, "%s (up to quadrupole):\n"
                        "alpha = %g, rc = %g, d1 = %g, d2 = %g, c2 = %g, c4 = %g, u(Rc) = %g\n",
                        eel_names[eel], zmm_alpha, Rc, d1, d2, *c2, *c4, *c);
            }

            break;
        case 3:
            /* Zero-octupole method */
            d1 = 1. * pow(Rc, -2.0) * (std::erfc(arc) + 2. / sqrt(M_PI) * exp(-arc * arc) * arc);
            d2 = 2. * pow(Rc, -3.0) * (std::erfc(arc) + 2. / sqrt(M_PI) * exp(-arc * arc) * (arc + pow(arc, 3.0)));
            d3 = 6. * pow(Rc, -4.0) * (std::erfc(arc) + 2. / sqrt(M_PI) * exp(-arc * arc) * (arc + 2. / 3. * pow(arc, 3.0) + 2. / 3. * pow(arc, 5.0)));

            *c2 =  15. / 16. * d1 * pow(Rc, -1.0) + 7. / 16. * d2 * pow(Rc,  0.0) + 1. / 16. * d3 * pow(Rc,  1.0);
            *c4 = - 5. / 16. * d1 * pow(Rc, -3.0) - 5. / 16. * d2 * pow(Rc, -2.0) - 1. / 16. * d3 * pow(Rc, -1.0);
            *c6 =   1. / 16. * d1 * pow(Rc, -5.0) + 1. / 16. * d2 * pow(Rc, -4.0) + 1. / 48. * d3 * pow(Rc, -3.0);
            *c = - (std::erfc(arc) / Rc + *c2 * pow(Rc, 2.0) + *c4 * pow(Rc, 4.0) + *c6 * pow(Rc, 6.0));
            if(fplog)
            {
                please_cite(fplog, "Fukuda2014");
                please_cite(fplog, "Sakuraba2018");
                fprintf(fplog, "%s (up to octupole):\n"
                        "alpha = %g, rc = %g, d1 = %g, d2 = %g, d3 = %g, c2 = %g, c4 = %g, c6 = %g, u(Rc) = %g\n",
                        eel_names[eel], zmm_alpha, Rc, d1, d2, d3, *c2, *c4, *c6, *c);
            }
            break;
        default:
            gmx_fatal(FARGS, "Zero-Multipole method is only supported with 0 <= degree <= 3");
        }
    }
}

void calc_zmmfac(FILE *fplog, int eel, int zmm_degree, real zmm_alpha,
                 real Rc,
                 real *c2, real *c4, real *c6, real *c)
{
    double c2d, c4d, c6d, cd;
    calc_zmmfac_double(fplog, eel, zmm_degree, zmm_alpha,
                       Rc,
                       &c2d, &c4d, &c6d, &cd);
    *c2 = (real)c2d;
    *c4 = (real)c4d;
    *c6 = (real)c6d;
    *c = (real)cd;
}
