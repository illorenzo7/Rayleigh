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

        !=============================================
        ! Edit Below This Line (you may define your own variables below)

        Real*8, Allocatable :: ind_work_r(:,:,:), ind_work_t(:,:,:), ind_work_p(:,:,:)
        Real*8, Allocatable :: ind_work_r2(:,:,:), ind_work_t2(:,:,:), ind_work_p2(:,:,:)
        Real*8, Allocatable :: cbuffer(:,:,:,:)
        Real*8 :: del2b
        Real*8, Allocatable :: ovstheta(:), ovs2theta(:)
        Allocate(cbuffer(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))
        Allocate(ind_work_r(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Allocate(ind_work_t(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Allocate(ind_work_p(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Allocate(ind_work_r2(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Allocate(ind_work_t2(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Allocate(ind_work_p2(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Allocate(ovstheta(1:N_theta), ovs2theta(1:N_theta)) ! 1/sin; 1/sin^2
        ovstheta = 1.0d0/sintheta
        ovs2theta = 1.0d0/sin2theta

        !//////////////////////////////////////////////////
        !   Part 1.    Terms resulting full v cross full B.
        ! /////////////////////////////////////////////////

        !shear: B dot grad v
        Call ADotGradB(buffer,buffer,cbuffer,aindices = bindex, bindices=vindex)

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
                ! Del^2 {B_r}
                del2b = DDBUFF(PSI,dbrdrdr)+Two_Over_R(r)*buffer(PSI,dbrdr)
                del2b = del2b+OneOverRSquared(r)*(DDBUFF(PSI,dbrdtdt)+cottheta(t)*buffer(PSI,dbrdt))
                del2b = del2b+OneOverRSquared(r)*DDBUFF(PSI,dbrdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{B} }_r
                del2b = del2b-2.0d0*OneOverRsquared(r)*( &
                        buffer(PSI,br) + &
                        buffer(PSI,dbtdt)+buffer(PSI,btheta)*cottheta(t) + &
                        ovstheta(t)*buffer(PSI,dbpdp) )

                qty(PSI) = eta(r)*del2b
                qty(PSI) = qty(PSI)*buffer(PSI,br)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(idiff_work_t)) Then

            DO_PSI
                ! Del^2 {B_theta}
                del2b = DDBUFF(PSI,dbtdrdr)+Two_Over_R(r)*buffer(PSI,dbtdr)
                del2b = del2b+OneOverRSquared(r)*(DDBUFF(PSI,dbtdtdt)+cottheta(t)*buffer(PSI,dbtdt))
                del2b = del2b+OneOverRSquared(r)*DDBUFF(PSI,dbtdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{B} }_theta
                del2b = del2b +OneOverRSquared(r)*( 2.0d0*buffer(PSI,dbrdt) - &
                        ovs2theta(t)*(   buffer(PSI,btheta) + &
                        2.0d0*costheta(t)*buffer(PSI,dbpdp) ) )

                ! Add the contribution from a gradient in eta
                qty(PSI) = eta(r)*(del2b+buffer(PSI,curlbphi)*dlneta(r))
                qty(PSI) = qty(PSI)*buffer(PSI,btheta)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(idiff_work_p)) Then
            DO_PSI
                ! build Del^2{B_phi}
                del2b = DDBUFF(PSI,dbpdrdr)+Two_Over_R(r)*buffer(PSI,dbpdr)
                del2b = del2b+OneOverRSquared(r)*(DDBUFF(PSI,dbpdtdt)+cottheta(t)*buffer(PSI,dbpdt))
                del2b = del2b+OneOverRSquared(r)*DDBUFF(PSI,dbpdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{u} }_phi
                del2b = del2b +OneOverRSquared(r)*( 2.0d0*buffer(PSI,dbrdp)*ovstheta(t) - &
                        ovs2theta(t)*(   buffer(PSI,bphi) - &
                        2.0d0*costheta(t)*buffer(PSI,dbtdp) ) )

                qty(PSI) = eta(r)*(del2b-buffer(PSI,curlbtheta)*dlneta(r))
                qty(PSI) = qty(PSI)*buffer(PSI,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !//////////////////////////////////////////////////
        !   Part 2.    Terms resulting <v> cross B'
        ! /////////////////////////////////////////////////

        !shear: B dot grad v
        Call ADotGradB(fbuffer,m0_values,cbuffer,aindices = bindex, bindices=vindex)

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
                ! Del^2 {B_r}
                del2b = d2_m0(PSI2,dbrdrdr)+Two_Over_R(r)*m0_values(PSI2,dbrdr)
                del2b = del2b+OneOverRSquared(r)*(d2_m0(PSI2,dbrdtdt)+cottheta(t)*m0_values(PSI2,dbrdt))
                del2b = del2b+OneOverRSquared(r)*d2_m0(PSI2,dbrdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{B} }_r
                del2b = del2b-2.0d0*OneOverRsquared(r)*( &
                        m0_values(PSI2,br) + &
                        m0_values(PSI2,dbtdt)+m0_values(PSI2,btheta)*cottheta(t) + &
                        ovstheta(t)*m0_values(PSI2,dbpdp) )

                qty(PSI) = eta(r)*del2b
                qty(PSI) = qty(PSI)*m0_values(PSI2,br)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(idiff_work_mm_t)) Then

            DO_PSI
                ! Del^2 {B_theta}
                del2b = d2_m0(PSI2,dbtdrdr)+Two_Over_R(r)*m0_values(PSI2,dbtdr)
                del2b = del2b+OneOverRSquared(r)*(d2_m0(PSI2,dbtdtdt)+cottheta(t)*m0_values(PSI2,dbtdt))
                del2b = del2b+OneOverRSquared(r)*d2_m0(PSI2,dbtdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{B} }_theta
                del2b = del2b +OneOverRSquared(r)*( 2.0d0*m0_values(PSI2,dbrdt) - &
                        ovs2theta(t)*(   m0_values(PSI2,btheta) + &
                        2.0d0*costheta(t)*m0_values(PSI2,dbpdp) ) )

                ! Add the contribution from a gradient in eta
                qty(PSI) = eta(r)*(del2b+m0_values(PSI2,curlbphi)*dlneta(r))
                qty(PSI) = qty(PSI)*m0_values(PSI2,btheta)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(idiff_work_mm_p)) Then
            DO_PSI
                ! build Del^2{B_phi}
                del2b = d2_m0(PSI2,dbpdrdr)+Two_Over_R(r)*m0_values(PSI2,dbpdr)
                del2b = del2b+OneOverRSquared(r)*(d2_m0(PSI2,dbpdtdt)+cottheta(t)*m0_values(PSI2,dbpdt))
                del2b = del2b+OneOverRSquared(r)*d2_m0(PSI2,dbpdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{u} }_phi
                del2b = del2b +OneOverRSquared(r)*( 2.0d0*m0_values(PSI2,dbrdp)*ovstheta(t) - &
                        ovs2theta(t)*(   m0_values(PSI2,bphi) - &
                        2.0d0*costheta(t)*m0_values(PSI2,dbtdp) ) )

                qty(PSI) = eta(r)*(del2b-m0_values(PSI2,curlbtheta)*dlneta(r))
                qty(PSI) = qty(PSI)*m0_values(PSI2,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !//////////////////////////////////////////////////
        !   Part 5.    Terms resulting v' cross B': ppp and mpp 
        ! /////////////////////////////////////////////////

        !shear: B' dot grad v'
        Call ADotGradB(fbuffer,fbuffer,cbuffer,aindices = bindex, bindices=vindex)

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
            Call Add_Quantity(ind_work_r)
        Endif
        If (compute_quantity(induct_work_mpp_t)) Then
            DO_PSI
                ind_work_t2(PSI) = ind_work_t2(PSI)+&
                    fbuffer(PSI,vr)*ref%dlnrho(r)*fbuffer(PSI,btheta)*m0_values(PSI2,btheta)
            END_DO
            Call Add_Quantity(ind_work_t)
        Endif
        If (compute_quantity(induct_work_mpp_p)) Then
            DO_PSI
                ind_work_p2(PSI) = ind_work_p2(PSI)+&
                    fbuffer(PSI,vr)*ref%dlnrho(r)*fbuffer(PSI,bphi)*m0_values(PSI2,bphi)
            END_DO
            Call Add_Quantity(ind_work_p)
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
                ! Del^2 {B_r}
                del2b = d2_fbuffer(PSI,dbrdrdr)+Two_Over_R(r)*fbuffer(PSI,dbrdr)
                del2b = del2b+OneOverRSquared(r)*(d2_fbuffer(PSI,dbrdtdt)+cottheta(t)*fbuffer(PSI,dbrdt))
                del2b = del2b+OneOverRSquared(r)*d2_fbuffer(PSI,dbrdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{B} }_r
                del2b = del2b-2.0d0*OneOverRsquared(r)*( &
                        fbuffer(PSI,br) + &
                        fbuffer(PSI,dbtdt)+fbuffer(PSI,btheta)*cottheta(t) + &
                        ovstheta(t)*fbuffer(PSI,dbpdp) )

                qty(PSI) = eta(r)*del2b
                qty(PSI) = qty(PSI)*fbuffer(PSI,br)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(idiff_work_pp_t)) Then

            DO_PSI
                ! Del^2 {B_theta}
                del2b = d2_fbuffer(PSI,dbtdrdr)+Two_Over_R(r)*fbuffer(PSI,dbtdr)
                del2b = del2b+OneOverRSquared(r)*(d2_fbuffer(PSI,dbtdtdt)+cottheta(t)*fbuffer(PSI,dbtdt))
                del2b = del2b+OneOverRSquared(r)*d2_fbuffer(PSI,dbtdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{B} }_theta
                del2b = del2b +OneOverRSquared(r)*( 2.0d0*fbuffer(PSI,dbrdt) - &
                        ovs2theta(t)*(   fbuffer(PSI,btheta) + &
                        2.0d0*costheta(t)*fbuffer(PSI,dbpdp) ) )

                ! Add the contribution from a gradient in eta
                qty(PSI) = eta(r)*(del2b+fbuffer(PSI,curlbphi)*dlneta(r))
                qty(PSI) = qty(PSI)*fbuffer(PSI,btheta)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(idiff_work_pp_p)) Then
            DO_PSI
                ! build Del^2{B_phi}
                del2b = d2_fbuffer(PSI,dbpdrdr)+Two_Over_R(r)*fbuffer(PSI,dbpdr)
                del2b = del2b+OneOverRSquared(r)*(d2_fbuffer(PSI,dbpdtdt)+cottheta(t)*fbuffer(PSI,dbpdt))
                del2b = del2b+OneOverRSquared(r)*d2_fbuffer(PSI,dbpdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{u} }_phi
                del2b = del2b +OneOverRSquared(r)*( 2.0d0*fbuffer(PSI,dbrdp)*ovstheta(t) - &
                        ovs2theta(t)*(   fbuffer(PSI,bphi) - &
                        2.0d0*costheta(t)*fbuffer(PSI,dbtdp) ) )

                qty(PSI) = eta(r)*(del2b-fbuffer(PSI,curlbtheta)*dlneta(r))
                qty(PSI) = qty(PSI)*fbuffer(PSI,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        DeAllocate(ind_work_r)
        DeAllocate(ind_work_t)
        DeAllocate(ind_work_p)
        DeAllocate(ind_work_r2)
        DeAllocate(ind_work_t2)
        DeAllocate(ind_work_p2)
        DeAllocate(cbuffer)



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

