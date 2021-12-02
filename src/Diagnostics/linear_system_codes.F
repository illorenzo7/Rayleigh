    !//////////////////////////////////////////////////////////////////////////
    !///////////////////////////////////////////////////
    !       Linear System Outputs  (2900)
    Integer, Parameter :: linear_off = axial_off+100 ! :OFFSET CODE: 
    
    Integer, Parameter :: avar_out    = linear_off+1 ! :tex: $A$
    Integer, Parameter :: dadr_out    = linear_off+2 ! :tex: $\frac{\partial A}{\partial r}$
    Integer, Parameter :: d2adr2_out  = linear_off+3 ! :tex: $\frac{\partial^2 A}{\partial r^2}$        
    Integer, Parameter :: a_HLaplace  = linear_off+4 ! :tex: $-\frac{\ell(\ell+1)}{r^2}A$ 
    Integer, Parameter :: avar_diff   = linear_off+5 ! :tex: $\eta\left( \frac{\partial^2 A}{\partial r^2} -\frac{\ell(\ell+1)}{r^2}A  \right) +\frac{\partial\eta}{\partial r}\frac{\partial A}{\partial r} $
    
    Integer, Parameter :: cvar_out    = linear_off+6 ! :tex: $C$
    Integer, Parameter :: dcdr_out    = linear_off+7 ! :tex: $\frac{\partial C}{\partial r}$
    Integer, Parameter :: d2cdr2_out  = linear_off+8 ! :tex: $\frac{\partial^2 C}{\partial r^2}$        
    Integer, Parameter :: c_HLaplace  = linear_off+9 ! :tex: $-\frac{\ell(\ell+1)}{r^2}C$     
    Integer, Parameter :: cvar_diff   = linear_off+10 ! :tex: $\eta\left( \frac{\partial^2 C}{\partial r^2} -\frac{\ell(\ell+1)}{r^2}C  \right)$   
    
    Integer, Parameter :: bdiffr_out   = linear_off+11 ! :tex: $B_r\,\,Diffusion$
    Integer, Parameter :: bdifft_out   = linear_off+12 ! :tex: $B_\theta\,\,Diffusion$    
    Integer, Parameter :: bdiffp_out   = linear_off+13 ! :tex: $B_\phi\,\,Diffusion$    

