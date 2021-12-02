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
Module Diagnostics_Linear_System
    Use Diagnostics_Base
    Implicit None

Contains

    Subroutine Compute_Linear_System_Mag(buffer)
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        
        
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! A-related quantities
        If (compute_quantity(avar_out)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,avar_str)
            END_DO
            Call Add_Quantity(qty)            
        Endif
        
        If (compute_quantity(dadr_out)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dadr_str)
            END_DO
            Call Add_Quantity(qty)            
        Endif
        
        If (compute_quantity(d2adr2_out)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,d2adr2_str)
            END_DO
            Call Add_Quantity(qty)            
        Endif        

        If (compute_quantity(a_HLaplace)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,avar_lap_str)
            END_DO
            Call Add_Quantity(qty)            
        Endif        
        
        If (compute_quantity(avar_diff)) Then
            DO_PSI
                qty(PSI) = (buffer(PSI,avar_lap_str)+buffer(PSI,d2adr2_str))*eta(r)
            END_DO
            DO_PSI
                qty(PSI) = qty(PSI)+buffer(PSI,dadr_str)*A_Diffusion_Coefs_1(r)
            END_DO
            Call Add_Quantity(qty)            
        Endif
        
        
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! C-related variables
        If (compute_quantity(cvar_out)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,cvar_str)
            END_DO
            Call Add_Quantity(qty)            
        Endif
        
        If (compute_quantity(dcdr_out)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dcdr_str)
            END_DO
            Call Add_Quantity(qty)            
        Endif
        
        If (compute_quantity(d2cdr2_out)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,d2cdr2_str)
            END_DO
            Call Add_Quantity(qty)            
        Endif        

        If (compute_quantity(c_HLaplace)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,cvar_lap_str)
            END_DO
            Call Add_Quantity(qty)            
        Endif        
        
        If (compute_quantity(cvar_diff)) Then
            DO_PSI
                qty(PSI) = (buffer(PSI,cvar_lap_str)+buffer(PSI,d2cdr2_str))*eta(r)
            END_DO
            Call Add_Quantity(qty)            
        Endif        
        
        If (compute_quantity(bdiffr_out)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,bdiff_rs)
            END_DO
            Call Add_Quantity(qty)
        Endif
        
         If (compute_quantity(bdifft_out)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,bdiff_ts)
            END_DO
            Call Add_Quantity(qty)
         Endif       
         
         If (compute_quantity(bdiffp_out)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,bdiff_ps)
            END_DO
            Call Add_Quantity(qty)
         Endif        
    End Subroutine Compute_Linear_System_Mag



End Module Diagnostics_Linear_System
