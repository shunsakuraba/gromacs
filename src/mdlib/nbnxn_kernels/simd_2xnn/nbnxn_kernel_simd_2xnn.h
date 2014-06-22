/*
 * This file is part of the GROMACS molecular simulation package.
 *
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
/*
 * Note: this file was generated by the Verlet kernel generator for
 * kernel type 2xnn.
 */

#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/*! \brief Run-time dispatcher for nbnxn kernel functions. */
void
nbnxn_kernel_simd_2xnn(nbnxn_pairlist_set_t       *nbl_list,
                       const nbnxn_atomdata_t     *nbat,
                       const interaction_const_t  *ic,
                       int                         ewald_excl,
                       rvec                       *shift_vec,
                       int                         force_flags,
                       int                         clearF,
                       real                       *fshift,
                       real                       *Vc,
                       real                       *Vvdw);

/* Need an #include guard so that sim_util.c can include all
 * such files. */
#ifndef _nbnxn_kernel_simd_include_h
#define _nbnxn_kernel_simd_include_h
/*! \brief Typedefs for declaring kernel functions. */
typedef void (nbk_func_ener)(const nbnxn_pairlist_t     *nbl,
                             const nbnxn_atomdata_t     *nbat,
                             const interaction_const_t  *ic,
                             rvec                       *shift_vec,
                             real                       *f,
                             real                       *fshift,
                             real                       *Vvdw,
                             real                       *Vc);
typedef nbk_func_ener *p_nbk_func_ener;

typedef void (nbk_func_noener)(const nbnxn_pairlist_t     *nbl,
                               const nbnxn_atomdata_t     *nbat,
                               const interaction_const_t  *ic,
                               rvec                       *shift_vec,
                               real                       *f,
                               real                       *fshift);
typedef nbk_func_noener *p_nbk_func_noener;
#endif

nbk_func_ener         nbnxn_kernel_simd_2xnn_rf_comb_geom_energrp;
nbk_func_ener         nbnxn_kernel_simd_2xnn_rf_comb_lb_energrp;
nbk_func_ener         nbnxn_kernel_simd_2xnn_rf_comb_none_energrp;
nbk_func_ener         nbnxn_kernel_simd_2xnn_zq_comb_geom_energrp;
nbk_func_ener         nbnxn_kernel_simd_2xnn_zq_comb_lb_energrp;
nbk_func_ener         nbnxn_kernel_simd_2xnn_zq_comb_none_energrp;
nbk_func_ener         nbnxn_kernel_simd_2xnn_tab_comb_geom_energrp;
nbk_func_ener         nbnxn_kernel_simd_2xnn_tab_comb_lb_energrp;
nbk_func_ener         nbnxn_kernel_simd_2xnn_tab_comb_none_energrp;
nbk_func_ener         nbnxn_kernel_simd_2xnn_tab_twin_comb_geom_energrp;
nbk_func_ener         nbnxn_kernel_simd_2xnn_tab_twin_comb_lb_energrp;
nbk_func_ener         nbnxn_kernel_simd_2xnn_tab_twin_comb_none_energrp;
nbk_func_ener         nbnxn_kernel_simd_2xnn_ewald_comb_geom_energrp;
nbk_func_ener         nbnxn_kernel_simd_2xnn_ewald_comb_lb_energrp;
nbk_func_ener         nbnxn_kernel_simd_2xnn_ewald_comb_none_energrp;
nbk_func_ener         nbnxn_kernel_simd_2xnn_ewald_twin_comb_geom_energrp;
nbk_func_ener         nbnxn_kernel_simd_2xnn_ewald_twin_comb_lb_energrp;
nbk_func_ener         nbnxn_kernel_simd_2xnn_ewald_twin_comb_none_energrp;

nbk_func_noener       nbnxn_kernel_simd_2xnn_rf_comb_geom_noener;
nbk_func_noener       nbnxn_kernel_simd_2xnn_rf_comb_lb_noener;
nbk_func_noener       nbnxn_kernel_simd_2xnn_rf_comb_none_noener;
nbk_func_noener       nbnxn_kernel_simd_2xnn_zq_comb_geom_noener;
nbk_func_noener       nbnxn_kernel_simd_2xnn_zq_comb_lb_noener;
nbk_func_noener       nbnxn_kernel_simd_2xnn_zq_comb_none_noener;
nbk_func_noener       nbnxn_kernel_simd_2xnn_tab_comb_geom_noener;
nbk_func_noener       nbnxn_kernel_simd_2xnn_tab_comb_lb_noener;
nbk_func_noener       nbnxn_kernel_simd_2xnn_tab_comb_none_noener;
nbk_func_noener       nbnxn_kernel_simd_2xnn_tab_twin_comb_geom_noener;
nbk_func_noener       nbnxn_kernel_simd_2xnn_tab_twin_comb_lb_noener;
nbk_func_noener       nbnxn_kernel_simd_2xnn_tab_twin_comb_none_noener;
nbk_func_noener       nbnxn_kernel_simd_2xnn_ewald_comb_geom_noener;
nbk_func_noener       nbnxn_kernel_simd_2xnn_ewald_comb_lb_noener;
nbk_func_noener       nbnxn_kernel_simd_2xnn_ewald_comb_none_noener;
nbk_func_noener       nbnxn_kernel_simd_2xnn_ewald_twin_comb_geom_noener;
nbk_func_noener       nbnxn_kernel_simd_2xnn_ewald_twin_comb_lb_noener;
nbk_func_noener       nbnxn_kernel_simd_2xnn_ewald_twin_comb_none_noener;

nbk_func_ener         nbnxn_kernel_simd_2xnn_rf_comb_geom_ener;
nbk_func_ener         nbnxn_kernel_simd_2xnn_rf_comb_lb_ener;
nbk_func_ener         nbnxn_kernel_simd_2xnn_rf_comb_none_ener;
nbk_func_ener         nbnxn_kernel_simd_2xnn_zq_comb_geom_ener;
nbk_func_ener         nbnxn_kernel_simd_2xnn_zq_comb_lb_ener;
nbk_func_ener         nbnxn_kernel_simd_2xnn_zq_comb_none_ener;
nbk_func_ener         nbnxn_kernel_simd_2xnn_tab_comb_geom_ener;
nbk_func_ener         nbnxn_kernel_simd_2xnn_tab_comb_lb_ener;
nbk_func_ener         nbnxn_kernel_simd_2xnn_tab_comb_none_ener;
nbk_func_ener         nbnxn_kernel_simd_2xnn_tab_twin_comb_geom_ener;
nbk_func_ener         nbnxn_kernel_simd_2xnn_tab_twin_comb_lb_ener;
nbk_func_ener         nbnxn_kernel_simd_2xnn_tab_twin_comb_none_ener;
nbk_func_ener         nbnxn_kernel_simd_2xnn_ewald_comb_geom_ener;
nbk_func_ener         nbnxn_kernel_simd_2xnn_ewald_comb_lb_ener;
nbk_func_ener         nbnxn_kernel_simd_2xnn_ewald_comb_none_ener;
nbk_func_ener         nbnxn_kernel_simd_2xnn_ewald_twin_comb_geom_ener;
nbk_func_ener         nbnxn_kernel_simd_2xnn_ewald_twin_comb_lb_ener;
nbk_func_ener         nbnxn_kernel_simd_2xnn_ewald_twin_comb_none_ener;



#if 0
{
#endif
#ifdef __cplusplus
}
#endif
