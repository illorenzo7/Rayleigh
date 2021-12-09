!
!  Copyright (C) 2018 by the authors of the RAYLEIGH code.
!
!  This file is part of RAYLEIGH.
!
!  RAYLEIGH is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 3, or (at your option)
!  any later version.
!
!  RAYLEIGH is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with RAYLEIGH; see the file LICENSE.  If not see
!  <http://www.gnu.org/licenses/>.
!


#include "indices.F"
Module Diagnostics_Base
    !//////////////////////////////////////////////////////////
    ! This module holds common variables that may be accessed
    ! by any diagnostics routine.  These variables are primarily
    ! temporary, allocatable arrays and the output menu codes.
    Use ProblemSize
    Use Spherical_IO
    Use Fields
    Use Math_Constants
    Use PDE_Coefficients

    Implicit None

    !/////////////////////////////////////////////////////////
    !  Diagnostic Quantity Codes
    !  Reynolds decompositions are often used in the outputs.
    !  As a result, some shorthand is used as follows:
    !   "m" and "< >" denote the azimuthal OR spherical mean.
    !   "p" and " ' " denote perturbations about that mean





    !////////////////////////////////////////////////////////////
    ! First, we include all system variables and their gradients

    include "velocity_field_codes.F"
    include "mass_flux_codes.F"
    include "vorticity_field_codes.F"
    include "kinetic_energy_codes.F"

    include "thermal_field_codes.F"
    include "thermal_energy_codes.F"

    include "magnetic_field_codes.F"
    include "current_density_codes.F"
    include "magnetic_energy_codes.F"

    include "momentum_equation_codes.F"
    include "thermal_equation_codes.F"
    include "induction_equation_codes.F"

    include "amom_equation_codes.F"
    include "ke_equation_codes.F"
    include "me_equation_codes.F"

    

    ! We have some "known" outputs as well that allow us to verify that
    ! the spherical_io interface is functional
    Integer, Parameter :: dcheck_off = meq_off+ 100  !2100
    Integer, Parameter :: diagnostic1 = dcheck_off+1
    Integer, Parameter :: diagnostic2 = dcheck_off+2
    Integer, Parameter :: test_y11 = dcheck_off+3
    Integer, Parameter :: test_y22 = dcheck_off+4
    Integer, Parameter :: test_y22_sq = dcheck_off+5


    !//////////////////////////////////////////////////////////
    !  Custom Outputs:  range from ...
    Integer, Parameter :: custom_offset = dcheck_off+100 !2200
    Integer, Parameter :: set_offset = 15 
    Integer, Parameter :: set_offset2 = 12
    Integer, Parameter :: set_offset3 = 9

    ! total
    Integer, Parameter :: induct_work_r     = custom_offset + 1
    Integer, Parameter :: induct_work_t     = custom_offset + 2
    Integer, Parameter :: induct_work_p     = custom_offset + 3

    Integer, Parameter :: ishear_work_r     = custom_offset + 4
    Integer, Parameter :: ishear_work_t     = custom_offset + 5
    Integer, Parameter :: ishear_work_p     = custom_offset + 6

    Integer, Parameter :: iadvec_work_r     = custom_offset + 7
    Integer, Parameter :: iadvec_work_t     = custom_offset + 8
    Integer, Parameter :: iadvec_work_p     = custom_offset + 9

    Integer, Parameter :: icomp_work_r     = custom_offset + 10
    Integer, Parameter :: icomp_work_t     = custom_offset + 11
    Integer, Parameter :: icomp_work_p     = custom_offset + 12

    Integer, Parameter :: idiff_work_r     = custom_offset + 13
    Integer, Parameter :: idiff_work_t     = custom_offset + 14
    Integer, Parameter :: idiff_work_p     = custom_offset + 15

    !pmp
    Integer, Parameter :: induct_work_pmp_r     = custom_offset + 1 + set_offset
    Integer, Parameter :: induct_work_pmp_t     = custom_offset + 2 + set_offset
    Integer, Parameter :: induct_work_pmp_p     = custom_offset + 3 + set_offset

    Integer, Parameter :: ishear_work_pmp_r     = custom_offset + 4 + set_offset
    Integer, Parameter :: ishear_work_pmp_t     = custom_offset + 5 + set_offset
    Integer, Parameter :: ishear_work_pmp_p     = custom_offset + 6 + set_offset

    Integer, Parameter :: iadvec_work_pmp_r     = custom_offset + 7 + set_offset
    Integer, Parameter :: iadvec_work_pmp_t     = custom_offset + 8 + set_offset
    Integer, Parameter :: iadvec_work_pmp_p     = custom_offset + 9 + set_offset

    Integer, Parameter :: icomp_work_pmp_r     = custom_offset + 10 + set_offset
    Integer, Parameter :: icomp_work_pmp_t     = custom_offset + 11 + set_offset
    Integer, Parameter :: icomp_work_pmp_p     = custom_offset + 12 + set_offset

    !ppm
    Integer, Parameter :: induct_work_ppm_r     = custom_offset + 1 + set_offset + set_offset2
    Integer, Parameter :: induct_work_ppm_t     = custom_offset + 2 + set_offset + set_offset2
    Integer, Parameter :: induct_work_ppm_p     = custom_offset + 3 + set_offset + set_offset2

    Integer, Parameter :: ishear_work_ppm_r     = custom_offset + 4 + set_offset + set_offset2
    Integer, Parameter :: ishear_work_ppm_t     = custom_offset + 5 + set_offset + set_offset2
    Integer, Parameter :: ishear_work_ppm_p     = custom_offset + 6 + set_offset + set_offset2

    Integer, Parameter :: iadvec_work_ppm_r     = custom_offset + 7 + set_offset + set_offset2
    Integer, Parameter :: iadvec_work_ppm_t     = custom_offset + 8 + set_offset + set_offset2
    Integer, Parameter :: iadvec_work_ppm_p     = custom_offset + 9 + set_offset + set_offset2

    Integer, Parameter :: icomp_work_ppm_r     = custom_offset + 10 + set_offset + set_offset2
    Integer, Parameter :: icomp_work_ppm_t     = custom_offset + 11 + set_offset + set_offset2
    Integer, Parameter :: icomp_work_ppm_p     = custom_offset + 12 + set_offset + set_offset2
!mmm
    Integer, Parameter :: induct_work_mmm_r     = custom_offset + 1 + set_offset + 2*set_offset2
    Integer, Parameter :: induct_work_mmm_t     = custom_offset + 2 + set_offset + 2*set_offset2
    Integer, Parameter :: induct_work_mmm_p     = custom_offset + 3 + set_offset + 2*set_offset2

    Integer, Parameter :: ishear_work_mmm_r     = custom_offset + 4 + set_offset + 2*set_offset2
    Integer, Parameter :: ishear_work_mmm_t     = custom_offset + 5 + set_offset + 2*set_offset2
    Integer, Parameter :: ishear_work_mmm_p     = custom_offset + 6 + set_offset + 2*set_offset2

    Integer, Parameter :: iadvec_work_mmm_r     = custom_offset + 7 + set_offset + 2*set_offset2
    Integer, Parameter :: iadvec_work_mmm_t     = custom_offset + 8 + set_offset + 2*set_offset2
    Integer, Parameter :: iadvec_work_mmm_p     = custom_offset + 9 + set_offset + 2*set_offset2

    Integer, Parameter :: icomp_work_mmm_r     = custom_offset + 10 + set_offset + 2*set_offset2
    Integer, Parameter :: icomp_work_mmm_t     = custom_offset + 11 + set_offset + 2*set_offset2
    Integer, Parameter :: icomp_work_mmm_p     = custom_offset + 12 + set_offset + 2*set_offset2

    Integer, Parameter :: idiff_work_mm_r     = custom_offset + 13 + set_offset + 2*set_offset2
    Integer, Parameter :: idiff_work_mm_t     = custom_offset + 14 + set_offset + 2*set_offset2
    Integer, Parameter :: idiff_work_mm_p     = custom_offset + 15 + set_offset + 2*set_offset2

    !mpp
    Integer, Parameter :: induct_work_mpp_r     = custom_offset + 1 + 2*set_offset + 2*set_offset2
    Integer, Parameter :: induct_work_mpp_t     = custom_offset + 2 + 2*set_offset + 2*set_offset2
    Integer, Parameter :: induct_work_mpp_p     = custom_offset + 3 + 2*set_offset + 2*set_offset2

    Integer, Parameter :: ishear_work_mpp_r     = custom_offset + 4 + 2*set_offset + 2*set_offset2
    Integer, Parameter :: ishear_work_mpp_t     = custom_offset + 5 + 2*set_offset + 2*set_offset2
    Integer, Parameter :: ishear_work_mpp_p     = custom_offset + 6 + 2*set_offset + 2*set_offset2

    Integer, Parameter :: iadvec_work_mpp_r     = custom_offset + 7 + 2*set_offset + 2*set_offset2
    Integer, Parameter :: iadvec_work_mpp_t     = custom_offset + 8 + 2*set_offset + 2*set_offset2
    Integer, Parameter :: iadvec_work_mpp_p     = custom_offset + 9 + 2*set_offset + 2*set_offset2

    !ppp
    Integer, Parameter :: induct_work_ppp_r     = custom_offset + 1 + 2*set_offset + 2*set_offset2 + set_offset3
    Integer, Parameter :: induct_work_ppp_t     = custom_offset + 2 + 2*set_offset + 2*set_offset2 + set_offset3
    Integer, Parameter :: induct_work_ppp_p     = custom_offset + 3 + 2*set_offset + 2*set_offset2 + set_offset3

    Integer, Parameter :: ishear_work_ppp_r     = custom_offset + 4 + 2*set_offset + 2*set_offset2 + set_offset3
    Integer, Parameter :: ishear_work_ppp_t     = custom_offset + 5 + 2*set_offset + 2*set_offset2 + set_offset3
    Integer, Parameter :: ishear_work_ppp_p     = custom_offset + 6 + 2*set_offset + 2*set_offset2 + set_offset3

    Integer, Parameter :: iadvec_work_ppp_r     = custom_offset + 7 + 2*set_offset + 2*set_offset2 + set_offset3
    Integer, Parameter :: iadvec_work_ppp_t     = custom_offset + 8 + 2*set_offset + 2*set_offset2 + set_offset3
    Integer, Parameter :: iadvec_work_ppp_p     = custom_offset + 9 + 2*set_offset + 2*set_offset2 + set_offset3

    Integer, Parameter :: icomp_work_ppp_r     = custom_offset + 10 + 2*set_offset + 2*set_offset2 + set_offset3
    Integer, Parameter :: icomp_work_ppp_t     = custom_offset + 11 + 2*set_offset + 2*set_offset2 + set_offset3
    Integer, Parameter :: icomp_work_ppp_p     = custom_offset + 12 + 2*set_offset + 2*set_offset2 + set_offset3

    Integer, Parameter :: idiff_work_pp_r     = custom_offset + 13 + 2*set_offset + 2*set_offset2 + set_offset3
    Integer, Parameter :: idiff_work_pp_t     = custom_offset + 14 + 2*set_offset + 2*set_offset2 + set_offset3
    Integer, Parameter :: idiff_work_pp_p     = custom_offset + 15 + 2*set_offset + 2*set_offset2 + set_offset3

    ! ALTERNATE induction
    Integer, Parameter :: inductalt_work_r     = custom_offset + 100 + 1
    Integer, Parameter :: inductalt_work_t     = custom_offset + 100 + 2
    Integer, Parameter :: inductalt_work_p     = custom_offset + 100 + 3

    Integer, Parameter :: ialtshear_work_r     = custom_offset + 100 + 4
    Integer, Parameter :: ialtshear_work_t     = custom_offset + 100 + 5
    Integer, Parameter :: ialtshear_work_p     = custom_offset + 100 + 6

    Integer, Parameter :: ialtadvec_work_r     = custom_offset + 100 + 7
    Integer, Parameter :: ialtadvec_work_t     = custom_offset + 100 + 8
    Integer, Parameter :: ialtadvec_work_p     = custom_offset + 100 + 9

    Integer, Parameter :: ialtcomp_work_r     = custom_offset + 100 + 10
    Integer, Parameter :: ialtcomp_work_t     = custom_offset + 100 + 11
    Integer, Parameter :: ialtcomp_work_p     = custom_offset + 100 + 12

    ! Alternate induction, each term broken up more
    Integer, Parameter :: ialtshear_work_r1     = custom_offset + 200 + 1
    Integer, Parameter :: ialtshear_work_r2     = custom_offset + 200 + 2
    Integer, Parameter :: ialtshear_work_t1     = custom_offset + 200 + 3
    Integer, Parameter :: ialtshear_work_t2     = custom_offset + 200 + 4
    Integer, Parameter :: ialtshear_work_p1     = custom_offset + 200 + 5
    Integer, Parameter :: ialtshear_work_p2     = custom_offset + 200 + 6

    Integer, Parameter :: ialtadvec_work_r1     = custom_offset + 200 + 7
    Integer, Parameter :: ialtadvec_work_r2     = custom_offset + 200 + 8
    Integer, Parameter :: ialtadvec_work_r3     = custom_offset + 200 + 9
    Integer, Parameter :: ialtadvec_work_r4     = custom_offset + 200 + 10
    Integer, Parameter :: ialtadvec_work_t1     = custom_offset + 200 + 11
    Integer, Parameter :: ialtadvec_work_t2     = custom_offset + 200 + 12
    Integer, Parameter :: ialtadvec_work_t3     = custom_offset + 200 + 13
    Integer, Parameter :: ialtadvec_work_t4     = custom_offset + 200 + 14
    Integer, Parameter :: ialtadvec_work_p1     = custom_offset + 200 + 15
    Integer, Parameter :: ialtadvec_work_p2     = custom_offset + 200 + 16
    Integer, Parameter :: ialtadvec_work_p3     = custom_offset + 200 + 17
    Integer, Parameter :: ialtadvec_work_p4     = custom_offset + 200 + 18

    Integer, Parameter :: ialtcomp_work_r1     = custom_offset + 200 + 19
    Integer, Parameter :: ialtcomp_work_r2     = custom_offset + 200 + 20
    Integer, Parameter :: ialtcomp_work_t1     = custom_offset + 200 + 21
    Integer, Parameter :: ialtcomp_work_t2     = custom_offset + 200 + 22
    Integer, Parameter :: ialtcomp_work_p1     = custom_offset + 200 + 23
    Integer, Parameter :: ialtcomp_work_p2     = custom_offset + 200 + 24

    include "turbKE_codes.F"
    include "axial_field_codes.F"
    include "linear_system_codes.F"
    

    !///////////////////////////////////
    Real*8, Allocatable :: qty(:,:,:)   ! This variable holds each quantity that we output
    Real*8, Allocatable :: tmp1(:,:,:), tmp4(:,:,:)  ! Work arrays  -- so sorry about the 4
    Real*8, Allocatable :: rweights(:), tweights(:), tmp1d(:)

    !//////////////////////////////////
    ! The ell0 and m0 _ values arrays contain, yes, the ell = 0 and m = 0 values of
    ! everything in buffer at output time.
    Real*8, Allocatable :: ell0_values(:,:), m0_values(:,:,:)

    ! This array will hold fluctuating quantities from the buffer { q - <q>}
    Real*8, Allocatable :: fbuffer(:,:,:,:)

    Logical :: azimuthal_mean = .true. ! when false, the m0_values are overwritten with the ell0_values

    !///////////////////////////////////////////////////////////////////////////
    ! A special buffer used for holding second derivatives at output time
    Type(SphericalBuffer) :: d2buffer
    ! ell0 and m0 values of those variables stored in d2buffer
    Real*8, Allocatable :: d2_ell0(:,:), d2_m0(:,:,:)

    ! This array will hold fluctuating quantities from the d2buffer { q - <q>}
    Real*8, Allocatable :: d2_fbuffer(:,:,:,:)


    ! Indices within the d2buffer
    Integer :: dvrdrdr, dvrdtdt, dvrdpdp, dvrdrdt, dvrdrdp, dvrdtdp
    Integer :: dvtdrdr, dvtdtdt, dvtdpdp, dvtdrdt, dvtdrdp, dvtdtdp
    Integer :: dvpdrdr, dvpdtdt, dvpdpdp, dvpdrdt, dvpdrdp, dvpdtdp

    Integer :: dtdrdr, dtdtdt, dtdpdp, dtdrdt, dtdrdp, dtdtdp
    Integer :: dpdrdr, dpdtdt, dpdpdp, dpdrdt, dpdrdp, dpdtdp

    Integer :: dbrdrdr, dbrdtdt, dbrdpdp, dbrdrdt, dbrdrdp, dbrdtdp
    Integer :: dbtdrdr, dbtdtdt, dbtdpdp, dbtdrdt, dbtdrdp, dbtdtdp
    Integer :: dbpdrdr, dbpdtdt, dbpdpdp, dbpdrdt, dbpdrdp, dbpdtdp

    Logical :: need_second_derivatives = .false.


    !////////////////////////////////////////////////////////////////////////////
    ! Variables related to mean-correction
    ! (we only correct radial terms, but retain logic for horizontal terms)
    Integer :: ncorrect = 0  ! Number of fields whose ell=0 mean we need to substract
    Real*8, Allocatable :: mean_3dbuffer(:,:,:,:)
    Real*8, Allocatable :: mean_ell0buffer(:,:)
    Integer :: cforce_r, cforce_theta, cforce_phi
    Integer :: aforce_r, aforce_theta, aforce_phi
    Integer :: aforcepp_r, aforcepp_theta, aforcepp_phi
    Integer :: aforcemm_r, aforcemm_theta, aforcemm_phi
    Integer :: vforce_r
    Integer :: lforce_r, lforcepp_r, lforcemm_r
Contains

    Subroutine Generate_Diagnostic_Labels()
        ! Define labels for our quantity codes
        Write(6,*)'A line of code.'
        !Call Load_Label(v_r,'V_r')
        !Call Load_Label(v_theta,'V_theta')
        !Call Load_Label(v_phi, 'V_phi')
    End Subroutine Generate_Diagnostic_Labels

    Subroutine Initialize_Diagnostics_Buffer()
        Logical :: dbtrans, dbconfig
        Logical :: test_reduce

        dbtrans = .false.
        dbconfig = .false.
        test_reduce = .false.

        Call cobuffer%init(field_count = cbfcount, config = 'p1a', &
            dynamic_transpose =dbtrans, dynamic_config = dbconfig, &
            hold_cargo = test_reduce, padding = pad_alltoall)
    End Subroutine Initialize_Diagnostics_Buffer


    Subroutine Compute_Fluctuations(buffer)
        Implicit None
        Integer :: r,k, t, j,jmax
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        jmax = size(buffer,4)
        Allocate(fbuffer(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:jmax))


        Do j = 1, jmax
            DO_PSI
                fbuffer(PSI,j) = buffer(PSI,j) - m0_values(PSI2,j)
            END_DO
        Enddo

    End Subroutine Compute_Fluctuations

    Subroutine DeAllocate_Fluctuations()
        Implicit None
        DeAllocate(fbuffer)
    End Subroutine DeAllocate_Fluctuations


End Module Diagnostics_Base
