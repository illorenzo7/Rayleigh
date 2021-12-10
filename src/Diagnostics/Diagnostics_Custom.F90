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
!///////////////////////////////////////////////////////////////////
!               DIAGNOSTICS_CUSTOM
!   This is where users can implement their own diagnostics
!   Custom diagnostics are most easily formed using the contents of the buffer
!
!   That array is dimensioned as:
!   buffer(1:n_phi+2, my_r%min:my_r%max, my_theta%min:my_theta%max,1:nvariables)
!
!   The extra 2 in the first index is needed for the in-place FFTs.  Care should be taken
!   to only loop over 1 to n_phi.
!
!   Each index along the 4th dimension of buffer corresponds to a different variable.
!   These indices (left) and the variables they correspond to (right) are given below.
!
! Field variables:
!   vr      -- radial velocity
!   vtheta  -- theta velocity
!   vphi    -- phi velocity
!   tvar    -- temperature or entropy
!   pvar    -- pressure
!   zvar    -- l(l+1)*Z/r^2  where Z is the toroidal streamfunction
!
! Radial Derivatives:
!   dvrdr   -- d(v_r)/dr
!   dvtdr   -- d(v_theta)/dr
!   dvpdr   -- d(v_phi)/dr
!   dtdr    -- d(temperature or entropy)/dr
!
!
! Theta Derivatives:
!   dvrdt   -- d(v_r)/dtheta
!   dvtdt   -- d(v_theta)/dtheta
!   dvpdt   -- d(v_phi)/dtheta
!   dtdt    -- (1/r)*d(temperature or entropy)/dtheta (<--- Note 1/r)
!
!
! Phi Derivatives:
!   dvrdp   --  d(v_r)/dphi
!   dvtdp   --  d(v_theta)/dphi
!   dvpdp   --  d(v_phi)/dphi
!   dtdp    --  (1/r)*d(temperature or entropy)/dphi   (<--- Note 1/r)
!
!
! If Magnetism is On, six additional variables are present:
!   br      -- radial magnetic field
!   btheta  -- theta magnetic field
!   bphi    -- phi magnetic field
!   curlbr      -- [Del x B]_r
!   curlbtheta  -- [Del x B]_theta
!   curlbphi    -- [Del x B]_phi
!
!
! If Induction Output is needed for this iteration,
!   the buffer also holds the derivatives of each
!   component of B.
!
! Radial Derivatives:
!   dbrdr   -- d(b_r)/dr
!   dbtdr   -- d(b_theta)/dr
!   dbpdr   -- d(b_phi)/dr
!
!
! Theta Derivatives:
!   dbrdt   -- d(b_r)/dtheta
!   dbtdt   -- d(b_theta)/dtheta
!   dbpdt   -- d(b_phi)/dtheta


! Phi Derivatives:
!   dbrdp   --  d(b_r)/dphi
!   dbtdp   --  d(b_theta)/dphi
!   dbpdp   --  d(b_phi)/dphi
!///////////////////////////////////////////////////////////////////

Module Diagnostics_Custom
    Use Diagnostics_Base
    Use Diagnostics_ADotGradB
    Implicit None
Contains

    Subroutine Custom_MHD_Diagnostics(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Integer :: j, jmax, offset, offset2, indextmp

        !=============================================
        ! Edit Below This Line (you may define your own variables below)

        Real*8, Allocatable :: ind_work_r(:,:,:), ind_work_t(:,:,:), ind_work_p(:,:,:)
        Real*8, Allocatable :: ind_work_r2(:,:,:), ind_work_t2(:,:,:), ind_work_p2(:,:,:)
        Real*8, Allocatable :: cbuffer(:,:,:,:)
        Real*8, Allocatable :: buff1(:,:,:,:),buff2(:,:,:,:)
        Real*8, Allocatable :: inducttmp(:,:,:,:), sheartmp(:,:,:,:), advtmp(:,:,:,:), comptmp(:,:,:,:), bfieldtmp(:,:,:,:)
        Real*8, Allocatable :: sheartmp1(:,:,:,:), sheartmp2(:,:,:,:)
        Real*8, Allocatable :: advtmp1(:,:,:,:), advtmp2(:,:,:,:), advtmp3(:,:,:,:), advtmp4(:,:,:,:), advtmp5(:,:,:,:)
        Real*8, Allocatable :: comptmp1(:,:,:,:), comptmp2(:,:,:,:)


        Real*8 :: del2b
        Real*8, Allocatable :: ovstheta(:), ovs2theta(:)

        jmax = size(buffer,4)
        Allocate(cbuffer(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))
        Allocate(ind_work_r(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Allocate(ind_work_t(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Allocate(ind_work_p(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Allocate(ind_work_r2(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Allocate(ind_work_t2(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Allocate(ind_work_p2(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Allocate(ovstheta(1:N_theta), ovs2theta(1:N_theta)) ! 1/sin; 1/sin^2

        Allocate(buff1(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:jmax))
        Allocate(buff2(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:jmax))

        Allocate(inducttmp(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))
        Allocate(sheartmp(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))
        Allocate(advtmp(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))
        Allocate(comptmp(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))
        Allocate(bfieldtmp(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))

        Allocate(sheartmp1(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))
        Allocate(sheartmp2(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))

        Allocate(advtmp1(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))
        Allocate(advtmp2(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))
        Allocate(advtmp3(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))
        Allocate(advtmp4(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))
        Allocate(advtmp5(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))

        Allocate(comptmp1(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))
        Allocate(comptmp2(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))

        ovstheta = 1.0d0/sintheta
        ovs2theta = 1.0d0/sin2theta

        !//////////////////////////////////////////////////
        !   Part 1.    Terms resulting full v cross full B.
        ! /////////////////////////////////////////////////

        !shear: B dot grad v
        Call ADotGradB(buffer,buffer,cbuffer,aindices = bindex, bindices=vindex)
        ! Remove canceling curvature terms
        DO_PSI
            cbuffer(PSI,1) = cbuffer(PSI,1) + (buffer(PSI,vtheta)*buffer(PSI,btheta) +&
               &buffer(PSI,vphi)*buffer(PSI,bphi))/radius(r)
            cbuffer(PSI,2) = cbuffer(PSI,2) + (buffer(PSI,vphi)*buffer(PSI,bphi))*cottheta(t)/radius(r)
        END_DO


        If (compute_quantity(ishear_work_r)) Then
            DO_PSI
                qty(PSI) = cbuffer(PSI,1)*buffer(PSI,br)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(ishear_work_t)) Then
            DO_PSI
                qty(PSI) = cbuffer(PSI,2)*buffer(PSI,btheta)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(ishear_work_p)) Then
            DO_PSI
                qty(PSI) = cbuffer(PSI,3)*buffer(PSI,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(induct_work_r)) Then
            DO_PSI
                ind_work_r(PSI) = cbuffer(PSI,1)*buffer(PSI,br)
            END_DO
        Endif

        If (compute_quantity(induct_work_t)) Then
            DO_PSI
                ind_work_t(PSI) = cbuffer(PSI,2)*buffer(PSI,btheta)
            END_DO
        Endif

        If (compute_quantity(induct_work_p)) Then
            DO_PSI
                ind_work_p(PSI) = cbuffer(PSI,3)*buffer(PSI,bphi)
            END_DO
        Endif

        ! advection: -v dot grad B
        Call ADotGradB(buffer,buffer,cbuffer,aindices = vindex, bindices=bindex)
        ! Remove canceling curvature terms
        DO_PSI
            cbuffer(PSI,1) = cbuffer(PSI,1) + (buffer(PSI,vtheta)*buffer(PSI,btheta) +&
               &buffer(PSI,vphi)*buffer(PSI,bphi))/radius(r)
            cbuffer(PSI,2) = cbuffer(PSI,2) + (buffer(PSI,vphi)*buffer(PSI,bphi))*cottheta(t)/radius(r)
        END_DO

        If (compute_quantity(iadvec_work_r)) Then
            DO_PSI
                qty(PSI) = -cbuffer(PSI,1)*buffer(PSI,br)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(iadvec_work_t)) Then
            DO_PSI
                qty(PSI) = -cbuffer(PSI,2)*buffer(PSI,btheta)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(iadvec_work_p)) Then
            DO_PSI
                qty(PSI) = -cbuffer(PSI,3)*buffer(PSI,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(induct_work_r)) Then
            DO_PSI
                ind_work_r(PSI) = ind_work_r(PSI)-&
                                cbuffer(PSI,1)*buffer(PSI,br)
            END_DO
        Endif
        If (compute_quantity(induct_work_t)) Then
            DO_PSI
                ind_work_t(PSI) = ind_work_t(PSI)-&
                                cbuffer(PSI,2)*buffer(PSI,btheta)
            END_DO
        Endif
        If (compute_quantity(induct_work_p)) Then
            DO_PSI
                ind_work_p(PSI) = ind_work_p(PSI)-&
                                cbuffer(PSI,3)*buffer(PSI,bphi)
            END_DO
        Endif

        ! compression: -B (div v)
        If (compute_quantity(induct_work_r)) Then
            DO_PSI
                ind_work_r(PSI) = ind_work_r(PSI)+&
                    buffer(PSI,vr)*ref%dlnrho(r)*buffer(PSI,br)**2
            END_DO
            Call Add_Quantity(ind_work_r)
        Endif

        If (compute_quantity(induct_work_t)) Then
            DO_PSI
                ind_work_t(PSI) = ind_work_t(PSI)+&
                    buffer(PSI,vr)*ref%dlnrho(r)*buffer(PSI,btheta)**2
            END_DO
            Call Add_Quantity(ind_work_t)
        Endif
        If (compute_quantity(induct_work_p)) Then
            DO_PSI
                ind_work_p(PSI) = ind_work_p(PSI)+&
                    buffer(PSI,vr)*ref%dlnrho(r)*buffer(PSI,bphi)**2
            END_DO
            Call Add_Quantity(ind_work_p)
        Endif

        If (compute_quantity(icomp_work_r)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,vr)*ref%dlnrho(r)*&
                    buffer(PSI,br)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(icomp_work_t)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,vr)*ref%dlnrho(r)*&
                    buffer(PSI,btheta)**2
            END_DO
            Call Add_Quantity(qty)
        Endif
        
        If (compute_quantity(icomp_work_p)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,vr)*ref%dlnrho(r)*&
                    buffer(PSI,bphi)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! diffusion: del X (eta del X B)

        If (compute_quantity(idiff_work_r)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,bdiff_rs)*buffer(PSI,br)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(idiff_work_t)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,bdiff_ts)*buffer(PSI,btheta)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(idiff_work_p)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,bdiff_ps)*buffer(PSI,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !//////////////////////////////////////////////////
        !   Part 2.    Terms resulting <v> cross B'
        ! /////////////////////////////////////////////////

        !shear: B dot grad v
        Call ADotGradB(fbuffer,m0_values,cbuffer,aindices = bindex, bindices=vindex)
        ! Remove canceling curvature terms
        DO_PSI
            cbuffer(PSI,1) = cbuffer(PSI,1) + (m0_values(PSI2,vtheta)*fbuffer(PSI,btheta) +&
               &m0_values(PSI2,vphi)*fbuffer(PSI,bphi))/radius(r)
            cbuffer(PSI,2) = cbuffer(PSI,2) + (m0_values(PSI2,vphi)*fbuffer(PSI,bphi))*cottheta(t)/radius(r)
        END_DO

        If (compute_quantity(ishear_work_pmp_r)) Then
            DO_PSI
                qty(PSI) = cbuffer(PSI,1)*fbuffer(PSI,br)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(ishear_work_pmp_t)) Then
            DO_PSI
                qty(PSI) = cbuffer(PSI,2)*fbuffer(PSI,btheta)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(ishear_work_pmp_p)) Then
            DO_PSI
                qty(PSI) = cbuffer(PSI,3)*fbuffer(PSI,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(induct_work_pmp_r)) Then
            DO_PSI
                ind_work_r(PSI) = cbuffer(PSI,1)*fbuffer(PSI,br)
            END_DO
        Endif

        If (compute_quantity(induct_work_pmp_t)) Then
            DO_PSI
                ind_work_t(PSI) = cbuffer(PSI,2)*fbuffer(PSI,btheta)
            END_DO
        Endif

        If (compute_quantity(induct_work_pmp_p)) Then
            DO_PSI
                ind_work_p(PSI) = cbuffer(PSI,3)*fbuffer(PSI,bphi)
            END_DO
        Endif

        ! advection: -v dot grad B
        Call ADotGradB(m0_values,fbuffer,cbuffer,aindices = vindex, bindices=bindex)
        ! Remove canceling curvature terms
        DO_PSI
            cbuffer(PSI,1) = cbuffer(PSI,1) + (m0_values(PSI2,vtheta)*fbuffer(PSI,btheta) +&
               &m0_values(PSI2,vphi)*fbuffer(PSI,bphi))/radius(r)
            cbuffer(PSI,2) = cbuffer(PSI,2) + (m0_values(PSI2,vphi)*fbuffer(PSI,bphi))*cottheta(t)/radius(r)
        END_DO

        If (compute_quantity(iadvec_work_pmp_r)) Then
            DO_PSI
                qty(PSI) = -cbuffer(PSI,1)*fbuffer(PSI,br)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(iadvec_work_pmp_t)) Then
            DO_PSI
                qty(PSI) = -cbuffer(PSI,2)*fbuffer(PSI,btheta)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(iadvec_work_pmp_p)) Then
            DO_PSI
                qty(PSI) = -cbuffer(PSI,3)*fbuffer(PSI,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(induct_work_pmp_r)) Then
            DO_PSI
                ind_work_r(PSI) = ind_work_r(PSI)-&
                                cbuffer(PSI,1)*fbuffer(PSI,br)
            END_DO
        Endif
        If (compute_quantity(induct_work_pmp_t)) Then
            DO_PSI
                ind_work_t(PSI) = ind_work_t(PSI)-&
                                cbuffer(PSI,2)*fbuffer(PSI,btheta)
            END_DO
        Endif
        If (compute_quantity(induct_work_pmp_p)) Then
            DO_PSI
                ind_work_p(PSI) = ind_work_p(PSI)-&
                                cbuffer(PSI,3)*fbuffer(PSI,bphi)
            END_DO
        Endif

        ! compression: -B (div v)
        If (compute_quantity(induct_work_pmp_r)) Then
            DO_PSI
                ind_work_r(PSI) = ind_work_r(PSI)+&
                    m0_values(PSI2,vr)*ref%dlnrho(r)*fbuffer(PSI,br)**2
            END_DO
            Call Add_Quantity(ind_work_r)
        Endif

        If (compute_quantity(induct_work_pmp_t)) Then
            DO_PSI
                ind_work_t(PSI) = ind_work_t(PSI)+&
                    m0_values(PSI2,vr)*ref%dlnrho(r)*fbuffer(PSI,btheta)**2
            END_DO
            Call Add_Quantity(ind_work_t)
        Endif
        If (compute_quantity(induct_work_pmp_p)) Then
            DO_PSI
                ind_work_p(PSI) = ind_work_p(PSI)+&
                    m0_values(PSI2,vr)*ref%dlnrho(r)*fbuffer(PSI,bphi)**2
            END_DO
            Call Add_Quantity(ind_work_p)
        Endif

        If (compute_quantity(icomp_work_pmp_r)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vr)*ref%dlnrho(r)*&
                    fbuffer(PSI,br)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(icomp_work_pmp_t)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vr)*ref%dlnrho(r)*&
                    fbuffer(PSI,btheta)**2
            END_DO
            Call Add_Quantity(qty)
        Endif
        
        If (compute_quantity(icomp_work_pmp_p)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vr)*ref%dlnrho(r)*&
                    fbuffer(PSI,bphi)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! diffusion: calculate in ppp and mmm

        !//////////////////////////////////////////////////
        !   Part 3.    Terms resulting v' cross <B>
        ! /////////////////////////////////////////////////

        !shear: B dot grad v
        Call ADotGradB(m0_values,fbuffer,cbuffer,aindices = bindex, bindices=vindex)
        ! Remove canceling curvature terms
        DO_PSI
            cbuffer(PSI,1) = cbuffer(PSI,1) + (m0_values(PSI2,btheta)*fbuffer(PSI,vtheta) +&
               &m0_values(PSI2,bphi)*fbuffer(PSI,vphi))/radius(r)
            cbuffer(PSI,2) = cbuffer(PSI,2) + (m0_values(PSI2,bphi)*fbuffer(PSI,vphi))*cottheta(t)/radius(r)
        END_DO

        If (compute_quantity(ishear_work_ppm_r)) Then
            DO_PSI
                qty(PSI) = cbuffer(PSI,1)*fbuffer(PSI,br)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(ishear_work_ppm_t)) Then
            DO_PSI
                qty(PSI) = cbuffer(PSI,2)*fbuffer(PSI,btheta)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(ishear_work_ppm_p)) Then
            DO_PSI
                qty(PSI) = cbuffer(PSI,3)*fbuffer(PSI,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(induct_work_ppm_r)) Then
            DO_PSI
                ind_work_r(PSI) = cbuffer(PSI,1)*fbuffer(PSI,br)
            END_DO
        Endif

        If (compute_quantity(induct_work_ppm_t)) Then
            DO_PSI
                ind_work_t(PSI) = cbuffer(PSI,2)*fbuffer(PSI,btheta)
            END_DO
        Endif

        If (compute_quantity(induct_work_ppm_p)) Then
            DO_PSI
                ind_work_p(PSI) = cbuffer(PSI,3)*fbuffer(PSI,bphi)
            END_DO
        Endif

        ! advection: -v' dot grad <B>
        Call ADotGradB(fbuffer,m0_values,cbuffer,aindices = vindex, bindices=bindex)
        ! Remove canceling curvature terms
        DO_PSI
            cbuffer(PSI,1) = cbuffer(PSI,1) + (m0_values(PSI2,btheta)*fbuffer(PSI,vtheta) +&
               &m0_values(PSI2,bphi)*fbuffer(PSI,vphi))/radius(r)
            cbuffer(PSI,2) = cbuffer(PSI,2) + (m0_values(PSI2,bphi)*fbuffer(PSI,vphi))*cottheta(t)/radius(r)
        END_DO

        If (compute_quantity(iadvec_work_ppm_r)) Then
            DO_PSI
                qty(PSI) = -cbuffer(PSI,1)*fbuffer(PSI,br)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(iadvec_work_ppm_t)) Then
            DO_PSI
                qty(PSI) = -cbuffer(PSI,2)*fbuffer(PSI,btheta)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(iadvec_work_ppm_p)) Then
            DO_PSI
                qty(PSI) = -cbuffer(PSI,3)*fbuffer(PSI,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(induct_work_ppm_r)) Then
            DO_PSI
                ind_work_r(PSI) = ind_work_r(PSI)-&
                                cbuffer(PSI,1)*fbuffer(PSI,br)
            END_DO
        Endif
        If (compute_quantity(induct_work_ppm_t)) Then
            DO_PSI
                ind_work_t(PSI) = ind_work_t(PSI)-&
                                cbuffer(PSI,2)*fbuffer(PSI,btheta)
            END_DO
        Endif
        If (compute_quantity(induct_work_ppm_p)) Then
            DO_PSI
                ind_work_p(PSI) = ind_work_p(PSI)-&
                                cbuffer(PSI,3)*fbuffer(PSI,bphi)
            END_DO
        Endif

        ! compression: -B (div v)
        If (compute_quantity(induct_work_ppm_r)) Then
            DO_PSI
                ind_work_r(PSI) = ind_work_r(PSI)+&
                    fbuffer(PSI,vr)*ref%dlnrho(r)*fbuffer(PSI,br)*m0_values(PSI2,br)
            END_DO
            Call Add_Quantity(ind_work_r)
        Endif

        If (compute_quantity(induct_work_ppm_t)) Then
            DO_PSI
                ind_work_t(PSI) = ind_work_t(PSI)+&
                    fbuffer(PSI,vr)*ref%dlnrho(r)*fbuffer(PSI,btheta)*m0_values(PSI2,btheta)
            END_DO
            Call Add_Quantity(ind_work_t)
        Endif
        If (compute_quantity(induct_work_ppm_p)) Then
            DO_PSI
                ind_work_p(PSI) = ind_work_p(PSI)+&
                    fbuffer(PSI,vr)*ref%dlnrho(r)*fbuffer(PSI,bphi)*m0_values(PSI2,bphi)
            END_DO
            Call Add_Quantity(ind_work_p)
        Endif

        If (compute_quantity(icomp_work_ppm_r)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,vr)*ref%dlnrho(r)*fbuffer(PSI,br)*m0_values(PSI2,br)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(icomp_work_ppm_t)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,vr)*ref%dlnrho(r)*fbuffer(PSI,btheta)*m0_values(PSI2,btheta)
            END_DO
            Call Add_Quantity(qty)
        Endif
        
        If (compute_quantity(icomp_work_ppm_p)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,vr)*ref%dlnrho(r)*fbuffer(PSI,bphi)*m0_values(PSI2,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! diffusion: calculate in ppp and mmm

        !//////////////////////////////////////////////////
        !   Part 4.    Terms resulting from <v> X <B>
        ! /////////////////////////////////////////////////

        !shear: <B> dot grad <v>
        Call ADotGradB(m0_values,m0_values,cbuffer,aindices = bindex, bindices=vindex)
        ! Remove canceling curvature terms
        DO_PSI
            cbuffer(PSI,1) = cbuffer(PSI,1) + (m0_values(PSI2,btheta)*m0_values(PSI2,vtheta) +&
               &m0_values(PSI2,bphi)*m0_values(PSI2,vphi))/radius(r)
            cbuffer(PSI,2) = cbuffer(PSI,2) + (m0_values(PSI2,bphi)*m0_values(PSI2,vphi))*cottheta(t)/radius(r)
        END_DO

        If (compute_quantity(ishear_work_mmm_r)) Then
            DO_PSI
                qty(PSI) = cbuffer(PSI,1)*m0_values(PSI2,br)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(ishear_work_mmm_t)) Then
            DO_PSI
                qty(PSI) = cbuffer(PSI,2)*m0_values(PSI2,btheta)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(ishear_work_mmm_p)) Then
            DO_PSI
                qty(PSI) = cbuffer(PSI,3)*m0_values(PSI2,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(induct_work_mmm_r)) Then
            DO_PSI
                ind_work_r(PSI) = cbuffer(PSI,1)*m0_values(PSI2,br)
            END_DO
        Endif

        If (compute_quantity(induct_work_mmm_t)) Then
            DO_PSI
                ind_work_t(PSI) = cbuffer(PSI,2)*m0_values(PSI2,btheta)
            END_DO
        Endif

        If (compute_quantity(induct_work_mmm_p)) Then
            DO_PSI
                ind_work_p(PSI) = cbuffer(PSI,3)*m0_values(PSI2,bphi)
            END_DO
        Endif

        ! advection: -v dot grad B
        Call ADotGradB(m0_values,m0_values,cbuffer,aindices = vindex, bindices=bindex)
        ! Remove canceling curvature terms
        DO_PSI
            cbuffer(PSI,1) = cbuffer(PSI,1) + (m0_values(PSI2,btheta)*m0_values(PSI2,vtheta) +&
               &m0_values(PSI2,bphi)*m0_values(PSI2,vphi))/radius(r)
            cbuffer(PSI,2) = cbuffer(PSI,2) + (m0_values(PSI2,bphi)*m0_values(PSI2,vphi))*cottheta(t)/radius(r)
        END_DO

        If (compute_quantity(iadvec_work_mmm_r)) Then
            DO_PSI
                qty(PSI) = -cbuffer(PSI,1)*m0_values(PSI2,br)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(iadvec_work_mmm_t)) Then
            DO_PSI
                qty(PSI) = -cbuffer(PSI,2)*m0_values(PSI2,btheta)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(iadvec_work_mmm_p)) Then
            DO_PSI
                qty(PSI) = -cbuffer(PSI,3)*m0_values(PSI2,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(induct_work_mmm_r)) Then
            DO_PSI
                ind_work_r(PSI) = ind_work_r(PSI)-&
                                cbuffer(PSI,1)*m0_values(PSI2,br)
            END_DO
        Endif
        If (compute_quantity(induct_work_t)) Then
            DO_PSI
                ind_work_t(PSI) = ind_work_t(PSI)-&
                                cbuffer(PSI,2)*m0_values(PSI2,btheta)
            END_DO
        Endif
        If (compute_quantity(induct_work_mmm_p)) Then
            DO_PSI
                ind_work_p(PSI) = ind_work_p(PSI)-&
                                cbuffer(PSI,3)*m0_values(PSI2,bphi)
            END_DO
        Endif

        ! compression: -B (div v)
        If (compute_quantity(induct_work_mmm_r)) Then
            DO_PSI
                ind_work_r(PSI) = ind_work_r(PSI)+&
                    m0_values(PSI2,vr)*ref%dlnrho(r)*m0_values(PSI2,br)**2
            END_DO
            Call Add_Quantity(ind_work_r)
        Endif

        If (compute_quantity(induct_work_mmm_t)) Then
            DO_PSI
                ind_work_t(PSI) = ind_work_t(PSI)+&
                    m0_values(PSI2,vr)*ref%dlnrho(r)*m0_values(PSI2,btheta)**2
            END_DO
            Call Add_Quantity(ind_work_t)
        Endif
        If (compute_quantity(induct_work_mmm_p)) Then
            DO_PSI
                ind_work_p(PSI) = ind_work_p(PSI)+&
                    m0_values(PSI2,vr)*ref%dlnrho(r)*m0_values(PSI2,bphi)**2
            END_DO
            Call Add_Quantity(ind_work_p)
        Endif

        If (compute_quantity(icomp_work_mmm_r)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vr)*ref%dlnrho(r)*&
                    m0_values(PSI2,br)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(icomp_work_mmm_t)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vr)*ref%dlnrho(r)*&
                    m0_values(PSI2,btheta)**2
            END_DO
            Call Add_Quantity(qty)
        Endif
        
        If (compute_quantity(icomp_work_mmm_p)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vr)*ref%dlnrho(r)*&
                    m0_values(PSI2,bphi)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! diffusion: del X (eta del X B)

        If (compute_quantity(idiff_work_mm_r)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,bdiff_rs)*m0_values(PSI2,br)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(idiff_work_mm_t)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,bdiff_ts)*m0_values(PSI2,btheta)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(idiff_work_mm_p)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,bdiff_ps)*m0_values(PSI2,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !//////////////////////////////////////////////////
        !   Part 5.    Terms resulting v' cross B': ppp and mpp 
        ! /////////////////////////////////////////////////

        !shear: B' dot grad v'
        Call ADotGradB(fbuffer,fbuffer,cbuffer,aindices = bindex, bindices=vindex)
        ! Remove canceling curvature terms
        DO_PSI
            cbuffer(PSI,1) = cbuffer(PSI,1) + (fbuffer(PSI,vtheta)*fbuffer(PSI,btheta) +&
               &fbuffer(PSI,vphi)*fbuffer(PSI,bphi))/radius(r)
            cbuffer(PSI,2) = cbuffer(PSI,2) + (fbuffer(PSI,vphi)*fbuffer(PSI,bphi))*cottheta(t)/radius(r)
        END_DO

        If (compute_quantity(ishear_work_ppp_r)) Then
            DO_PSI
                qty(PSI) = cbuffer(PSI,1)*fbuffer(PSI,br)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(ishear_work_mpp_r)) Then
            DO_PSI
                qty(PSI) = cbuffer(PSI,1)*m0_values(PSI2,br)
            END_DO
            Call Add_Quantity(qty)
        Endif        

        If (compute_quantity(ishear_work_ppp_t)) Then
            DO_PSI
                qty(PSI) = cbuffer(PSI,2)*fbuffer(PSI,btheta)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(ishear_work_mpp_t)) Then
            DO_PSI
                qty(PSI) = cbuffer(PSI,2)*m0_values(PSI2,btheta)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(ishear_work_ppp_p)) Then
            DO_PSI
                qty(PSI) = cbuffer(PSI,3)*fbuffer(PSI,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(ishear_work_mpp_p)) Then
            DO_PSI
                qty(PSI) = cbuffer(PSI,3)*m0_values(PSI2,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(induct_work_ppp_r)) Then
            DO_PSI
                ind_work_r(PSI) = cbuffer(PSI,1)*fbuffer(PSI,br)
            END_DO
        Endif

        If (compute_quantity(induct_work_mpp_r)) Then
            DO_PSI
                ind_work_r2(PSI) = cbuffer(PSI,1)*m0_values(PSI2,br)
            END_DO
        Endif

        If (compute_quantity(induct_work_ppp_t)) Then
            DO_PSI
                ind_work_t(PSI) = cbuffer(PSI,2)*fbuffer(PSI,btheta)
            END_DO
        Endif

        If (compute_quantity(induct_work_mpp_t)) Then
            DO_PSI
                ind_work_t2(PSI) = cbuffer(PSI,2)*m0_values(PSI2,btheta)
            END_DO
        Endif

        If (compute_quantity(induct_work_ppp_p)) Then
            DO_PSI
                ind_work_p(PSI) = cbuffer(PSI,3)*fbuffer(PSI,bphi)
            END_DO
        Endif

        If (compute_quantity(induct_work_mpp_p)) Then
            DO_PSI
                ind_work_p2(PSI) = cbuffer(PSI,3)*m0_values(PSI2,bphi)
            END_DO
        Endif

        ! advection: -v dot grad B
        Call ADotGradB(fbuffer,fbuffer,cbuffer,aindices = vindex, bindices=bindex)
        ! Remove canceling curvature terms
        DO_PSI
            cbuffer(PSI,1) = cbuffer(PSI,1) + (fbuffer(PSI,vtheta)*fbuffer(PSI,btheta) +&
               &fbuffer(PSI,vphi)*fbuffer(PSI,bphi))/radius(r)
            cbuffer(PSI,2) = cbuffer(PSI,2) + (fbuffer(PSI,vphi)*fbuffer(PSI,bphi))*cottheta(t)/radius(r)
        END_DO

        If (compute_quantity(iadvec_work_ppp_r)) Then
            DO_PSI
                qty(PSI) = -cbuffer(PSI,1)*fbuffer(PSI,br)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(iadvec_work_mpp_r)) Then
            DO_PSI
                qty(PSI) = -cbuffer(PSI,1)*m0_values(PSI2,br)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(iadvec_work_ppp_t)) Then
            DO_PSI
                qty(PSI) = -cbuffer(PSI,2)*fbuffer(PSI,btheta)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(iadvec_work_mpp_t)) Then
            DO_PSI
                qty(PSI) = -cbuffer(PSI,2)*m0_values(PSI2,btheta)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(iadvec_work_ppp_p)) Then
            DO_PSI
                qty(PSI) = -cbuffer(PSI,3)*fbuffer(PSI,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(iadvec_work_mpp_p)) Then
            DO_PSI
                qty(PSI) = -cbuffer(PSI,3)*m0_values(PSI2,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(induct_work_ppp_r)) Then
            DO_PSI
                ind_work_r(PSI) = ind_work_r(PSI)-&
                                cbuffer(PSI,1)*fbuffer(PSI,br)
            END_DO
        Endif
        If (compute_quantity(induct_work_ppp_t)) Then
            DO_PSI
                ind_work_t(PSI) = ind_work_t(PSI)-&
                                cbuffer(PSI,2)*fbuffer(PSI,btheta)
            END_DO
        Endif
        If (compute_quantity(induct_work_ppp_p)) Then
            DO_PSI
                ind_work_p(PSI) = ind_work_p(PSI)-&
                                cbuffer(PSI,3)*fbuffer(PSI,bphi)
            END_DO
        Endif

        ! mpp
        If (compute_quantity(induct_work_mpp_r)) Then
            DO_PSI
                ind_work_r2(PSI) = ind_work_r2(PSI)-&
                                cbuffer(PSI,1)*m0_values(PSI2,br)
            END_DO
        Endif
        If (compute_quantity(induct_work_mpp_t)) Then
            DO_PSI
                ind_work_t2(PSI) = ind_work_t2(PSI)-&
                                cbuffer(PSI,2)*m0_values(PSI2,btheta)
            END_DO
        Endif
        If (compute_quantity(induct_work_mpp_p)) Then
            DO_PSI
                ind_work_p2(PSI) = ind_work_p2(PSI)-&
                                cbuffer(PSI,3)*m0_values(PSI2,bphi)
            END_DO
        Endif

        ! compression: -B (div v): don't need mpp (except for inductive bit),
        ! since that equals ppm (already computed)

        !ppp, induct
        If (compute_quantity(induct_work_ppp_r)) Then
            DO_PSI
                ind_work_r(PSI) = ind_work_r(PSI)+&
                    fbuffer(PSI,vr)*ref%dlnrho(r)*fbuffer(PSI,br)**2
            END_DO
            Call Add_Quantity(ind_work_r)
        Endif
        If (compute_quantity(induct_work_ppp_t)) Then
            DO_PSI
                ind_work_t(PSI) = ind_work_t(PSI)+&
                    fbuffer(PSI,vr)*ref%dlnrho(r)*fbuffer(PSI,btheta)**2
            END_DO
            Call Add_Quantity(ind_work_t)
        Endif
        If (compute_quantity(induct_work_ppp_p)) Then
            DO_PSI
                ind_work_p(PSI) = ind_work_p(PSI)+&
                    fbuffer(PSI,vr)*ref%dlnrho(r)*fbuffer(PSI,bphi)**2
            END_DO
            Call Add_Quantity(ind_work_p)
        Endif
        
        !mpp,induct
        If (compute_quantity(induct_work_mpp_r)) Then
            DO_PSI
                ind_work_r2(PSI) = ind_work_r2(PSI)+&
                    fbuffer(PSI,vr)*ref%dlnrho(r)*fbuffer(PSI,br)*m0_values(PSI2,br)
            END_DO
            Call Add_Quantity(ind_work_r2)
        Endif
        If (compute_quantity(induct_work_mpp_t)) Then
            DO_PSI
                ind_work_t2(PSI) = ind_work_t2(PSI)+&
                    fbuffer(PSI,vr)*ref%dlnrho(r)*fbuffer(PSI,btheta)*m0_values(PSI2,btheta)
            END_DO
            Call Add_Quantity(ind_work_t2)
        Endif
        If (compute_quantity(induct_work_mpp_p)) Then
            DO_PSI
                ind_work_p2(PSI) = ind_work_p2(PSI)+&
                    fbuffer(PSI,vr)*ref%dlnrho(r)*fbuffer(PSI,bphi)*m0_values(PSI2,bphi)
            END_DO
            Call Add_Quantity(ind_work_p2)
        Endif
   
        !ppp, compression
        If (compute_quantity(icomp_work_ppp_r)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,vr)*ref%dlnrho(r)*&
                    fbuffer(PSI,br)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(icomp_work_ppp_t)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,vr)*ref%dlnrho(r)*&
                    fbuffer(PSI,btheta)**2
            END_DO
            Call Add_Quantity(qty)
        Endif
        
        If (compute_quantity(icomp_work_ppp_p)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,vr)*ref%dlnrho(r)*&
                    fbuffer(PSI,bphi)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! diffusion: del X (eta del X B')

        If (compute_quantity(idiff_work_pp_r)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,bdiff_rs)*fbuffer(PSI,br)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(idiff_work_pp_t)) Then

            DO_PSI
                qty(PSI) = fbuffer(PSI,bdiff_ts)*fbuffer(PSI,btheta)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(idiff_work_pp_p)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,bdiff_ps)*fbuffer(PSI,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !//////////////////////////////////////////////////
        !   ALTERNATE INDUCTION
        !   Part 1.    Terms resulting full v cross full B.
        ! /////////////////////////////////////////////////

        ! total production (triple product)
        offset = 0
        offset2 = 0
        DO_PSI
            Do j=1,jmax
                buff1(PSI,j) = buffer(PSI,j) ! v in del x (v x B)
                buff2(PSI,j) = buffer(PSI,j) ! B in del x (v x B)
            Enddo
            Do j=1,3
                If (j .eq. 1) Then
                    indextmp = br
                Elseif (j .eq. 2) Then
                    indextmp = btheta
                Elseif (j .eq. 3) Then
                    indextmp = bphi
                Endif
                bfieldtmp(PSI,j) = buffer(PSI,indextmp) ! bfieldtmp -- the one that gets dotted into induction
            Enddo
        END_DO

        ! Get Shear
        DO_PSI
            sheartmp1(PSI,1) = buff2(PSI,btheta)*buff1(PSI,dvrdt)/radius(r)
            sheartmp2(PSI,1) = buff2(PSI,bphi)*buff1(PSI,dvrdp)/radius(r)/sintheta(t)
            sheartmp1(PSI,2) = buff2(PSI,br)*(buff1(PSI,dvtdr) - buff1(PSI,vtheta)/radius(r))
            sheartmp2(PSI,2) = buff2(PSI,bphi)*buff1(PSI,dvtdp)/radius(r)/sintheta(t)
            sheartmp1(PSI,3) = buff2(PSI,br)*(buff1(PSI,dvpdr) - buff1(PSI,vphi)/radius(r))
            sheartmp2(PSI,3) = buff2(PSI,btheta)*(buff1(PSI,dvpdt) - buff1(PSI,vphi)*cottheta(t))/radius(r)
            Do j = 1,3
                sheartmp(PSI,j) = sheartmp1(PSI,j) + sheartmp2(PSI,j)
            Enddo
        END_DO    

        ! Get Advection
        DO_PSI
            advtmp1(PSI,1) = -buff1(PSI,vr)*buff2(PSI,dbrdr)
            advtmp2(PSI,1) = -buff1(PSI,vtheta)*buff2(PSI,dbrdt)/radius(r)
            advtmp3(PSI,1) = -buff1(PSI,vphi)*buff2(PSI,dbrdp)/radius(r)/sintheta(t)
            advtmp4(PSI,1) = -2.0d0*buff1(PSI,vr)*buff2(PSI,br)/radius(r)
            advtmp5(PSI,1) = 0.0d0

            advtmp1(PSI,2) = -buff1(PSI,vr)*buff2(PSI,dbtdr)
            advtmp2(PSI,2) = -buff1(PSI,vtheta)*buff2(PSI,dbtdt)/radius(r)
            advtmp3(PSI,2) = -buff1(PSI,vphi)*buff2(PSI,dbtdp)/radius(r)/sintheta(t)
            advtmp4(PSI,2) = buff1(PSI,vr)*buff2(PSI,btheta)/radius(r)
            advtmp5(PSI,2) = -buff1(PSI,vtheta)*buff2(PSI,btheta)*cottheta(t)/radius(r)

            advtmp1(PSI,3) = -buff1(PSI,vr)*buff2(PSI,dbpdr)
            advtmp2(PSI,3) = -buff1(PSI,vtheta)*buff2(PSI,dbpdt)/radius(r)
            advtmp3(PSI,3) = -buff1(PSI,vphi)*buff2(PSI,dbpdp)/radius(r)/sintheta(t)
            advtmp4(PSI,3) = buff1(PSI,vr)*buff2(PSI,bphi)/radius(r)
            advtmp5(PSI,3) = buff1(PSI,vtheta)*buff2(PSI,bphi)*cottheta(t)/radius(r)

            Do j = 1,3
                advtmp(PSI,j) = advtmp1(PSI,j) + advtmp2(PSI,j) + advtmp3(PSI,j) +&
                    &advtmp4(PSI,j) + advtmp5(PSI,j)
            Enddo
        END_DO

        ! get compression
        ! ALTERNATE compression: -B (div v) (but TRANSVERSE compression only)
        DO_PSI
            comptmp1(PSI,1) = -buff2(PSI,br)*(buff1(PSI,dvtdt)+buff1(PSI,vtheta)*cottheta(t))/radius(r)
            comptmp2(PSI,1) = -buff2(PSI,br)*buff1(PSI,dvpdp)/radius(r)/sintheta(t)
            comptmp1(PSI,2) = -buff2(PSI,btheta)*(buff1(PSI,dvrdr)+2.0d0*buff1(PSI,vr)/radius(r))
            comptmp2(PSI,2) = -buff2(PSI,btheta)*buff1(PSI,dvpdp)/radius(r)/sintheta(t)
            comptmp1(PSI,3) = -buff2(PSI,bphi)*(buff1(PSI,dvrdr)+2.0d0*buff1(PSI,vr)/radius(r))
            comptmp2(PSI,3) = -buff2(PSI,bphi)*(buff1(PSI,dvtdt)+buff1(PSI,vtheta)*cottheta(t))/radius(r)
            Do j = 1,3
                comptmp(PSI,j) = comptmp1(PSI,j) + comptmp2(PSI,j)
            Enddo
        END_DO    

        ! get total induction
        ! ALTERNATE total induction
        DO_PSI
            Do j=1,3
                inducttmp(PSI,j) = sheartmp(PSI,j) +  advtmp(PSI,j) + comptmp(PSI,j)
            Enddo
        END_DO

        ! Now add quantities we need

        ! total production (shear, comp, advec)
        Do j=1,3
            ! shear terms
            If (compute_quantity(ialtshear_work_r+j-1+offset)) Then
                DO_PSI
                    qty(PSI) = sheartmp(PSI,j)*bfieldtmp(PSI,j)
                END_DO
                Call Add_Quantity(qty)
            Endif

            ! advect. terms
            If (compute_quantity(ialtadvec_work_r+j-1+offset)) Then
                DO_PSI
                    qty(PSI) = advtmp(PSI,j)*bfieldtmp(PSI,j)
                END_DO
                Call Add_Quantity(qty)
            Endif
            
            ! Comp. terms
            If (compute_quantity(ialtcomp_work_r+j-1+offset)) Then
                DO_PSI
                    qty(PSI) = comptmp(PSI,j)*bfieldtmp(PSI,j)
                END_DO
                Call Add_Quantity(qty)
            Endif
            
            ! tot ind. terms
            If (compute_quantity(inductalt_work_r+j-1+offset)) Then
                DO_PSI
                    qty(PSI) = inducttmp(PSI,j)*bfieldtmp(PSI,j)
                END_DO
                Call Add_Quantity(qty)
            Endif

        Enddo

        ! production, broken up more
        Do j=1,3
            ! shear terms
            If (compute_quantity(ialtshear_work_r1+j-1+offset2)) Then
                DO_PSI
                    qty(PSI) = sheartmp1(PSI,j)*bfieldtmp(PSI,j)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(ialtshear_work_r2+j-1+offset2)) Then
                DO_PSI
                    qty(PSI) = sheartmp2(PSI,j)*bfieldtmp(PSI,j)
                END_DO
                Call Add_Quantity(qty)
            Endif

            ! advect. terms
            If (compute_quantity(ialtadvec_work_r1+j-1+offset2)) Then
                DO_PSI
                    qty(PSI) = advtmp1(PSI,j)*bfieldtmp(PSI,j)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(ialtadvec_work_r2+j-1+offset2)) Then
                DO_PSI
                    qty(PSI) = advtmp2(PSI,j)*bfieldtmp(PSI,j)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(ialtadvec_work_r3+j-1+offset2)) Then
                DO_PSI
                    qty(PSI) = advtmp3(PSI,j)*bfieldtmp(PSI,j)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(ialtadvec_work_r4+j-1+offset2)) Then
                DO_PSI
                    qty(PSI) = advtmp4(PSI,j)*bfieldtmp(PSI,j)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(ialtadvec_work_r5+j-1+offset2)) Then
                DO_PSI
                    qty(PSI) = advtmp5(PSI,j)*bfieldtmp(PSI,j)
                END_DO
                Call Add_Quantity(qty)
            Endif
            
            ! Comp. terms
            If (compute_quantity(ialtcomp_work_r1+j-1+offset2)) Then
                DO_PSI
                    qty(PSI) = comptmp1(PSI,j)*bfieldtmp(PSI,j)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(ialtcomp_work_r2+j-1+offset2)) Then
                DO_PSI
                    qty(PSI) = comptmp2(PSI,j)*bfieldtmp(PSI,j)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Enddo

        DeAllocate(ind_work_r)
        DeAllocate(ind_work_t)
        DeAllocate(ind_work_p)
        DeAllocate(ind_work_r2)
        DeAllocate(ind_work_t2)
        DeAllocate(ind_work_p2)
        DeAllocate(cbuffer)

        DeAllocate(buff1)
        DeAllocate(buff2)

        DeAllocate(inducttmp)
        DeAllocate(sheartmp)
        DeAllocate(advtmp)
        DeAllocate(comptmp)
        DeAllocate(bfieldtmp)

        DeAllocate(sheartmp1)
        DeAllocate(sheartmp2)

        DeAllocate(advtmp1)
        DeAllocate(advtmp2)
        DeAllocate(advtmp3)
        DeAllocate(advtmp4)
        DeAllocate(advtmp5)

        DeAllocate(comptmp1)
        DeAllocate(comptmp2)

        ! Edit Above This Line
        !=============================================


    End Subroutine Custom_MHD_Diagnostics

    Subroutine Custom_Hydro_Diagnostics(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        !=============================================
        ! Edit Below This Line (you may define your own variables below)



        ! Edit Above This Line
        !=============================================

    End Subroutine Custom_Hydro_Diagnostics

End Module Diagnostics_Custom

