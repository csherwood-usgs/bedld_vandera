      PROGRAM VANDERA 

      USE utils_constants
      implicit none 

!      real(kind=8) :: smgd, theta_mean, theta_wav, w_asym
!      real(kind=8) :: theta_max, theta_max1, theta_max2
!      real(kind=8) :: phi_x1, phi_x2, phi_x, phi_y
!      real(kind=8) :: bedld_x, bedld_y, tau_cur, waven, wavec

      integer      :: wavecycle    
      real(kind=8) :: Hs, Td, depth, d50, d90
      real(kind=8) :: umag_curr, phi_curwave 
      real(kind=8) :: y, kh, uhat, ahat 
      real(kind=8) :: k, c_w
      real(kind=8) :: smgd, osmgd
      real(kind=8) :: bed_frac 
!
      real(kind=8) :: r, phi, Su, Au
      real(kind=8) :: Sk, Ak
      real(kind=8) :: T_cu, T_tu, T_t, T_c, umax, umin 
      real(kind=8) :: theta_c, theta_t, tau_wRe 
      real(kind=8) :: theta_cx, theta_cy, theta_tx, theta_ty
      real(kind=8) :: mag_theta_c, mag_theta_t

      real(kind=8) :: eta, dsf                
      real(kind=8) :: uc_r, ut_r, uhat_c, uhat_t
      real(kind=8) :: uc_rx, uc_ry, alpha
      real(kind=8) :: ut_rx, ut_ry, Rcheck
      real(kind=8) :: mag_uc, mag_ut 

      real(kind=8) :: om_cc, om_tt, om_ct, om_tc
        
      real(kind=8) :: cff, cff1, cff2   
      real(kind=8) :: smgd_3, bedld_c, bedld_t, bedld
!            
!     variables that may be not needed by COAWST but declared
! to be deleted when in  COAWST 
        real ::  dtr,deg2rad,pi
      real(kind=8) ::  g, vk, rho0, rhos, nu
        pi=3.1416_r8                           
        dtr = pi/180.0_r8! % degrees to radians
        deg2rad=dtr   
!      DO ised=NCS+1,NST
!
!  Dimensionalizing bedload to get volumetric bed load transport per
!  unit width equation 2, page VA-2013
!
!        sd50_3=Sd50(ised,ng)**3
!      Input to this 
!      Hs=significant wave height
!      Td=dominant wave period 
!     depth=water depth(m)
!  
        Hs = 2.0_r8
        Td = 10.0_r8
        depth =  10.0_r8
! 
! change it in coawst for general sed type
       ! dsf=Sd50(ised,ng)
        d50 = 0.15e-3
!
!TSK Check d90 definitiaon
!
        d90 = 1.5*d50 !% check this...van Rijn has an eqn.

        ! constants
         g = 9.81;
        vk = 0.41;
        rho0 = 1027.
        rhos= 2650; 
         nu = 1.36e-6;

! umag_curr is the current velocity magnitude (MAKE Sure its the magnitude)
! direction of currents is counter-clockwise from wave direction
        umag_curr=ABS(0.5_r8)
        phi_curwave=45.0_r8*dtr

        y=kh(Td,depth)
!
! uhat is the wave orbital velocity defiend for the entire wave cycle
! TSK --. HARDWIRED for now
        uhat=depth*0.5_r8*((2.0_r8*g*y)/(depth*sinh(2.0_r8*y)))**0.5_r8 
        uhat=0.5472_r8           
! TSK check the origin of the equaiton
!
        ahat=uhat*Td/(2.0_r8*pi)

        k=kh(Td,depth)/depth     ! Wave number 
        c_w=2*pi/(k*Td)          ! Wave speed

!        smgdr=SQRT(smgd)


!        DO j=Jstrm1,Jendp1
!          DO i=Istrm1,Iendp1

! VA-2013 equation 1 is solved in 3 sub-steps
!
!-----------------------------------------------------------------------
!
! Ruessink et al. provides equations for calculating skewness parameters
! Uses Malarkey and Davies equations to get "bb" and "r"
! COMMON TO both crest and trough
!
            CALL skewness_params(Hs, Td, depth, r, phi, Su, Au)
!        
! Abreu et al. use skewness params to get representative critical orbital
! velocity for crest and trough cycles , get r and phi from above
! 
            CALL abreu_points(r, phi, uhat, Sk, Ak, Td,  T_c, T_t,      &
     &                                          T_cu, T_tu, umax, umin)

!
!------------          CREST HALF CYCLE --------------------------------
!-----------------------------------------------------------------------
! Get the "representative crest half cycle water particle velocity
!    as well as full cycle orbital velocity and excursion 
! CRS COMMENTS CHECK 
!% Crest and trough "representative" velocities:
!% (coordinate system is aligned with wave direction; A13, Fig. 2)
! % TODO - check to make sure these orbital velocities are correct and don't
! % have to be made into "representative" by multiplying by, say, sqrt(2)
!-----------------------------------------------------------------------
!
            ! from Abreu points 
            uhat_c=umax
!TSK ->
            uhat_t=-umin 
!
!           VanderA equation 10, 11
!
            uc_r=0.5_r8*sqrt(2.0_r8)*uhat_c
            ut_r=0.5_r8*sqrt(2.0_r8)*uhat_t

            uc_rx=uc_r+umag_curr*cos(phi_curwave)
            uc_ry=umag_curr*sin(phi_curwave)
            mag_uc=sqrt(uc_rx*uc_rx+uc_ry*uc_ry)
            print*, "crest vel.",uc_rx, uc_ry
!            alpha=umag_wave/(umag_wave+uhat)
! 
!           TROUGH HALF CYCLE 
!-----------------------------------------------------------------------
! 1. Get the representative trough half cycle water particle velocity 
!    as well as full cycle orbital velocity and excursion 
!-----------------------------------------------------------------------
!
            ut_rx=-ut_r+umag_curr*cos(phi_curwave)
            ut_ry=umag_curr*sin(phi_curwave)
            mag_ut=sqrt(ut_rx*ut_rx+ut_ry*ut_ry)
            print*, "trough vel.", ut_rx, ut_ry
!
! Need to chekc velocity skewnes parameter
!
            Rcheck = uhat_c/(uhat_c+uhat_t)
!            print*, "Rcheck",Rcheck
!            alpha=umag_wave/(umag_wave+uhat)
! 
!
!-----------------------------------------------------------------------
! 2. Bed shear stress (Shields parameter) for Crest half cycle 
!    alpha VA2013 Eqn. 19  
!-----------------------------------------------------------------------
! TSK       smgd=(Srho(ised,ng)/rho0-1.0_r8)*g*sd50(ised,ng)
        smgd=(rhos/rho0-1.0_r8)*g*d50
        osmgd=1.0_r8/smgd
! 
! TSK what happens when there are multiple sediments 
! 
!      DO ised=NCS+1,NST
         CALL stress_progressive_surface_waves(d50, d90, osmgd,         &
     &                                        Td, depth, rho0,          &
     &                                   umag_curr, uhat, uhat_c, ahat, &
     &                                        T_cu, T_c, mag_uc,        &
     &                                        eta, dsf,                 &
     &                                        mag_theta_c, tau_wRe) 

!
! Shields parameter in trough cycle
!
         theta_cx=ABS(mag_theta_c)*uc_rx/(mag_uc)+tau_wRe*osmgd/rho0
         theta_cy=ABS(mag_theta_c)*uc_ry/(mag_uc)
         mag_theta_c=sqrt(theta_cx*theta_cx+theta_cy*theta_cy)
         print*, "computed", mag_theta_c
         print*, "theta_cx",theta_cx
         print*, "theta_cy",theta_cy
        
!    
!-----------------------------------------------------------------------
! 2. Bed shear stress (Shields parameter) for Trough half cycle 
!    alpha VA2013 Eqn. 19  
!-----------------------------------------------------------------------
!
         CALL stress_progressive_surface_waves(d50, d90, osmgd,         &
     &                                        Td, depth, rho0,          &
     &                                   umag_curr, uhat, uhat_t, ahat, &
     &                                        T_tu, T_t, mag_ut,        &
     &                                        eta, dsf,                 &
     &                                        theta_t, tau_wRe) 
!
! Shields parameter in crest cycle
! rho0 is required for non-dimensionalizing 
!
         theta_tx=ABS(theta_t)*ut_rx/(mag_ut)+tau_wRe*osmgd/rho0
         theta_ty=ABS(theta_t)*ut_ry/(mag_ut)
         mag_theta_t=sqrt(theta_tx*theta_tx+theta_ty*theta_ty)
         print*, "theta_tx",theta_tx
         print*, "theta_ty",theta_tx
!
!-----------------------------------------------------------------------
! 3. Compute sediment load entrained during each crest half cycle
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!      Crest half cycle
!-----------------------------------------------------------------------
!
        wavecycle=1
        CALL sandload_vandera(wavecycle, d50, c_w,                      &
     &                              eta, dsf,                           &
     &                              T_c, T_cu, uhat_c, theta_c,         &
     &                              om_cc, om_ct)
!
!        omega_i=theta_crit
!
!-----------------------------------------------------------------------
!       Trough half cycle 
!-----------------------------------------------------------------------
!
        wavecycle=-1
        CALL sandload_vandera(wavecycle, d50, c_w,                      &
     &                              eta, dsf,                           &
     &                              T_t, T_tu, uhat_t, theta_t,         &
     &                              om_tt, om_tc)
!
!-----------------------------------------------------------------------
! VA2013  Use the velocity-load equation 1. 
! Non-dimensional net transport rate
!-----------------------------------------------------------------------
!
        cff1=T_c/(2.0_r8*T_cu)
        cff2=theta_c/mag_theta_c
        bedld_c=sqrt(mag_theta_c)*T_c*(om_cc+cff*om_tc)*cff2 
       
        cff1=T_t/(2.0_r8*T_tu)
        cff2=theta_t/mag_theta_t      
        bedld_t=sqrt(mag_theta_t)*T_t*(om_tt+cff*om_ct)*cff2 

        smgd_3=sqrt((rhos/rho0-1.0_r8)*g*d50**3.0_r8)
!
! The units of these are m2 sec-1
!
        bedld=smgd_3*(bedld_c+bedld_t)/Td 
!
! TSK-Check the units of this rhos
!
        bedld=bed_frac*rhos*bedld


!          END DO
!        END DO
!     
        END PROGRAM VANDERA
 
        SUBROUTINE sandload_vandera(wavecycle, d50, c_w,                &
     &                              eta, dsf,                           &
     &                              T_i, T_iu, uhat_i, theta_i,         &
     &                              om_ii, om_iy)
!
        USE utils_constants
        implicit none 

        integer,  intent(in) :: wavecycle 
        real(r8), intent(in) :: d50, c_w
        real(r8), intent(in) :: eta, dsf
        real(r8), intent(in) :: T_i, T_iu
        real(r8), intent(in) :: uhat_i, theta_i
        real(r8), intent(out):: om_ii, om_iy

        ! local variables
! 
! VA2012 Text under equation 37
! 
        real(r8), parameter :: m=11.0_r8, n=1.2_r8, alpha=8.2_r8
        real(r8) :: eps_eff
        real(r8) :: om_i
        real(r8) :: theta_ieff, theta_cr, theta_cr_calc 
        real(r8) :: w_sc_eta, w_s, w_min_eta, w_min_dsf 
        real(r8) :: w_sc_dsf, w_max_eta, w_max_dsf
        real(r8) :: w_s_calc
        real(r8) :: cff, cff1_eta, cff1_dsf
        real(r8) :: xi
        real(r8) :: P 
!
! TSK delete it in COAWST 
!      
        real(r8) :: g, rho0, rhos, nu
        g = 9.81_r8
        rho0 = 1027.0_r8
        rhos= 2650.0_r8
        nu = 1.36e-6_r8;
!
        eps_eff=(dsf/d50)**0.25_r8 
        theta_ieff=eps_eff*theta_i
! 
! Find critical Shields parameters based on Soulsby (1997).
!
        theta_cr=theta_cr_calc(0.8_r8*d50, rhos, rho0, nu)
        
        IF(theta_ieff.le.theta_cr) THEN 
          om_i=0.0_r8
        ELSE
          om_i=m*(theta_ieff-theta_cr)
        ENDIF 
! 
! Find settling velocity based on Soulsby (1997). 
! VA2012 Use 0.8*d50 for settling velocity (text under equation 28).
!
        w_s=w_s_calc(0.8_r8*d50, rhos, rho0, nu)       
!
! VA2012 Equation 29, for crest cycle
!
        IF(wavecycle.eq.1) THEN 
          w_sc_eta=w_s-w_min_eta
          w_sc_dsf=w_s-w_min_dsf
        ENDIF 
!
! VA2012 Equation 30, for trough cycle
!
        IF(wavecycle.eq.-1) THEN 
          w_sc_eta=MAX(w_s-w_max_eta,0.0_r8)
          w_sc_dsf=MAX(w_s-w_max_dsf,0.0_r8)
        ENDIF 
!
! VA2012 Equation 33, Phase lag parameter
!
        cff=(1.0_r8-wavecycle*xi*uhat_i)/c_w
        cff1_eta=(1.0_r8/(2.0_r8*(T_i-T_iu)*w_sc_eta))
        cff1_dsf=(1.0_r8/(2.0_r8*(T_i-T_iu)*w_sc_dsf))
!
! For ripple regime 
!
        IF(eta.gt.0.0_r8) THEN 
          P=alpha*eta*cff*cff1_eta
        ELSEIF(eta.eq.0.0_r8)THEN
!
! For sheet flow regime 
!
          P=alpha*dsf*cff*cff1_dsf
        ENDIF 
!
! VA2012 Equation 23-26, Sandload entrained during half cycle 
! 
        IF(P.le.1.0_r8) THEN 
          om_ii=om_i
        ELSE
          om_ii=om_i/P 
        ENDIF 
!
        IF(P.le.1.0_r8) THEN 
          om_iy=0.0_r8
        ELSE 
          cff=1.0_r8/P
          om_iy=om_i*(1.0_r8-cff) 
        ENDIF
!     
        END SUBROUTINE sandload_vandera
!
        SUBROUTINE stress_progressive_surface_waves(d50, d90, osmgd,    &
     &                                        Td, depth, rho0,          &
     &                           umag_curr, uhat, uhat_i, ahat,         &
     &                                        T_iu, T_i, mag_ui,        &
     &                                        eta, dsf,                 &
     &                                     mag_theta_i, tau_wRe) 
!
        USE utils_constants
        implicit none 
!  
!
!
! Input the crest or trough half cycle velocity
! d50 -- grain size in meters
! TSK ---> eta, lambda are ripple height and length (should we compute them)
! Different for crest and trough half cycles 
!       
        real(r8), intent(in) :: d50, d90, osmgd
        real(r8), intent(in) :: rho0, Td, depth
        real(r8), intent(in) :: umag_curr, uhat, uhat_i, ahat
        real(r8), intent(in) :: T_iu, T_i, mag_ui
        real(r8), intent(inout) :: eta, dsf ! ripple height and sheet flow thickness
        real(r8), intent(out):: mag_theta_i, tau_wRe
!
        integer  :: iter, total_iters
        real(r8) :: d50_mm, d90_mm
        real(r8) :: mu, mu_calc, theta_timeavg, theta_hat_i

        real(r8) :: psi
        real(r8) :: lambda  ! ripple dimensions
        real(r8) :: ksd, ksd_calc
        real(r8) :: ksw, ksw_calc
        real(r8) :: fw, fw_calc
        real(r8) :: fd, fd_calc
        real(r8) :: fw_i, fwi_calc, alpha, fwd_i 
        real(r8) ::  dsf_calc
        real(r8) :: alpha_w, fwd, kh, k, c_w
!  delete this is in COAWST 
        real(r8) ::  dtr,deg2rad,pi
        real(r8) ::  g, vk
        pi=3.1416_r8                           
        dtr = pi/180.0_r8! % degrees to radians
!
! Iterative solution to obtain current and wave related bed roughness
! VA2013 Apendix A, Shields parameter (Stress) depends on bed roughness 
! Bed roughness computed from converged Shields parameter
!
!
! Bed roughness factor based on grain size
!
        d50_mm=d50*0.001_r8  ! convert m to mm
        d90_mm=d90*0.001_r8  ! convert m to mm
!        
! Maximum mobility number at crest and trough 
! For irregular waves, use Rayleigh distributed maximum value
! VA, text under equation Appendix B.4 
!
        psi=(1.27_r8*uhat)**2*osmgd
!
! Use Appendix B eqn B.1 and B.2 to get ripple height and length 
!   
        CALL ripple_dim(psi, d50_mm, eta, lambda)
        eta=eta*ahat
        lambda=lambda*ahat
!
! Initiliaze with theta_timeavg=0
! 
        theta_timeavg=0.0_r8
        DO iter=1,total_iters
!
! Calculate current-related bed roughness from VA2013 Appendix A.1
!
          ksd=ksd_calc(d50, d90, mu_calc(d50_mm), theta_timeavg,        &
     &                                                  eta, lambda)    
! 
! Calculate wave related bed roughness from VA2013 A.5     
! 
          ksw=ksw_calc(d50, mu, theta_timeavg, eta, lambda)
!
! Calculate full-cycle wave friction factor VA2013 Appendix Eqn. A.4 
!
          fw=fw_calc(ahat, ksw)
!
! Calculate full-cycle current friction factor from VA2013 Eqn. 20
!
! TSK confusion of units d50 is in meters, ksd in mm 
! 
          fd=fd_calc(d50, ksd)
!
! Time-averaged absolute Shields stress VA2013 Appendix Eq. A.3
! 
          theta_timeavg=osmgd*(0.5_r8*fd*umag_curr**2.0_r8+             &
     &                         0.25_r8*fw*uhat**2.0_r8)
             
          tlast=theta_timeavg
          IF((theta_timeavg-tlast).lt.0.001_r8) THEN 
             print*, "iterations finished"
          ENDIF 
         END DO 
!
! Recompute to use converged value of time averaged Shields Stress
!
! Calculate current-related bed roughness from VA2013 Appendix A.1
!
          ksd=ksd_calc(d50, d90, mu_calc(d50_mm), theta_timeavg,        &
     &                                                  eta, lambda)    
!
! Calculate full-cycle current friction factor from VA2013 Eqn. 20
!
! TSK confusion of units d50 is in meters, ksd in mm 
! 
          fd=fd_calc(d50, ksd)
! 
! Calculate wave related bed roughness from VA2013 A.5     
! 
          ksw=ksw_calc(d50, mu, theta_timeavg, eta, lambda)
!
! Calculate full-cycle wave friction factor VA2013 Appendix Eqn. A.4 
!
          fw=fw_calc(ahat, ksw)
!        
! Wave friction factor for wave and crest half cycle VA2013 Eqn. 21
! 
          fw_i=fwi_calc(T_iu, T_i, ahat, ksw)
!
! Wave current friction factor (Madsen and Grant) VA2013 Eqn. 18
! Different for crest and trough 
!
          alpha=umag_curr/(umag_curr+uhat)
          fwd_i=alpha*fd+(1.0_r8-alpha)*fw_i
!
! theta_hat_i required to get sheet flow layer thickness Appendix C
!
          theta_hat_i=0.5_r8*fwd_i*uhat_i**2*osmgd
!
! Sheet flow thickness VA2013 Appendix C C.1 
! Update from initial value 
!
          dsf=dsf_calc(d50, theta_hat_i) !this dsf is in m 
!
! Calculate wave Reynolds stress from full cycle wave and friction factor
! that were formed from the iterative cycle, VA2013, Eqn.22
!
        alpha_w=0.424_r8
        fwd=alpha*fd+(1.0_r8-alpha)*fw
!
        k=kh(Td,depth)/depth     ! Wave number 
        c_w=2*pi/(k*Td)          ! Wave speed
!
        tau_wRe=rho0*fwd*alpha_w*uhat**3.0_r8/(2.0_r8*c_w)
        print*, "------------------------------"
         print*, "tau_wRe",tau_wRe
!
! VA2013- Recalculate fd, fwd_i to get magnitude of Shields parameter
!
        mag_theta_i=0.5_r8*fwd_i*mag_ui**2*osmgd
        
        END SUBROUTINE stress_progressive_surface_waves 
!
        SUBROUTINE ripple_dim(psi, d50_mm, eta, lambda)
!
! Calculate ripple dimensions of O'Donoghue et al. 2006
! based on VA2013 Appendix B
! Returns eta-Ripple length and lambda-Ripple length 
!        
        USE utils_constants
!
        real(r8), intent(in)  :: psi, d50_mm
        real(r8), intent(out) :: eta, lambda
!
        real(r8) :: m_eta, m_lambda, n_eta, n_lambda 
! to be deleted when in  COAWST 
        real ::  dtr,deg2rad,pi
        pi=3.1416_r8                           
        dtr = pi/180.0_r8! % degrees to radians
!
        IF(d50_mm.lt.0.22_r8) THEN
          m_eta=0.55_r8
          m_lambda=0.73_r8
        ELSEIF(d50_mm.ge.0.22_r8.and.d50_mm.lt.0.30_r8) THEN
          m_eta=0.55_r8+(0.45_r8*(d50_mm-0.22_r8)/(0.30_r8-0.22_r8))
          m_lambda=0.73_r8+(0.27_r8*(d50_mm-0.22_r8)/(0.30_r8-0.22_r8))
        ELSE
          m_eta=1.0_r8
          m_lambda=1.0_r8
        ENDIF
! 
! Smooth transition between ripple regime and bed sheet flow regime 
!
        IF(psi.le.190.0_r8) THEN
          n_eta=1.0_r8
        ELSEIF(psi.gt.190.0_r8.and.psi.lt.240.0_r8) THEN
          n_eta=0.5_r8*(1.0_r8+cos(pi*(psi-190.0_r8)/(50.0_r8)))
        ELSEIF(psi.ge.240.0_r8) THEN
          n_eta=0.0_r8
        ENDIF
        n_lambda=n_eta
!
        eta=MAX(0.0_r8,m_eta*n_eta*(0.275_r8-0.022*psi**0.42_r8))
        lambda=MAX(0.0_r8,m_lambda*n_lambda*                            &
     &                             (1.97_r8-0.44_r8*psi**0.21_r8))
!
        END SUBROUTINE ripple_dim
!
        REAL FUNCTION theta_cr_calc(d50, rhos, rho0, nu)
!
! Critical Shields parameter from Soulsby (1997).
!
        USE utils_constants
! 
        real(r8) :: d50, rhos, rho0, nu
        real(r8) :: s, dstar 
        real(r8) :: cff1, cff2
!
        s=rhos/rho0
        dstar=(g*(s-1)/(nu*nu))**(1.0_r8/3.0_r8)*d50
        cff1=0.30_r8/(1.0_r8+1.2_r8*dstar)
        cff2=0.055_r8*(1.0_r8-EXP(-0.020_r8*dstar))
        theta_cr_calc=cff1+cff2
!
        END FUNCTION theta_cr_calc
!
        REAL FUNCTION w_s_calc(d50, rhos, rho0, nu)
!
! Critical Shields parameter from Soulsby (1997).
!
        USE utils_constants
! 
        real(r8) :: d50, rhos, rho0, nu
        real(r8) :: s, dstar 
        real(r8) :: cff, cff1
!
        s=rhos/rho0
        dstar=(g*(s-1)/(nu*nu))**(1.0_r8/3.0_r8)*d50
        cff=nu/d50
        cff1=10.36_r8
        w_s_calc=cff*(sqrt(cff1*cff1+1.049_r8*dstar**3.0_r8)+    &
     &                   - cff1)
!
        END FUNCTION 
!
        REAL FUNCTION mu_calc(d50_mm)
!
! Calculate bed roughness factor based on grain size
! VA2013 Appendix A., required for current related bed roughness
! and wave related bed roughness. 
!
        USE utils_constants
!
        real(r8) :: d50_mm
!
        IF(d50_mm.le.0.15_r8) THEN
          mu_calc=6.0_r8
        ELSEIF(0.15_r8.gt.d50_mm.lt.0.20_r8) THEN
          mu_calc=6.0_r8-5.0_r8*((d50_mm-0.15_r8)/(0.2_r8-0.15_r8))
        ELSEIF(d50_mm.gt.0.20_r8) THEN
          mu_calc=1.0_r8
        ENDIF
!
        RETURN
        END FUNCTION
!
        REAL FUNCTION ksd_calc(d50, d90, mu, theta_timeavg,             &
     &                         eta, lambda)
        USE utils_constants
!
! Calculate current-related bed roughness from VA2013 Appendix A.1
!
        real(r8) :: d50, d90, mu, theta_timeavg, eta, lambda
        real(r8) :: ripple_fac
!
        eta=MAX(eta,d50)
        lambda=MAX(lambda,d50)
        ripple_fac=0.4_r8*eta**2.0_r8/lambda
        ksd_calc=MAX( 3.0_r8*d90,                                       &
     &             d50*(mu+6.0_r8*(theta_timeavg-1.0_r8)) )+            &
     &             ripple_fac
!
        RETURN 
        END FUNCTION ksd_calc
!
        REAL FUNCTION ksw_calc(d50, mu, theta_timeavg, eta, lambda)
        USE utils_constants
! 
! Calculate wave related bed roughness from VA2013 Eqn. A.5     
! 
        real(r8) :: d50, mu, theta_timeavg, eta, lambda
        real(r8) :: ripple_fac, ksw
!
        eta=MAX(eta,d50)
        lambda=MAX(lambda,d50)
        ripple_fac=0.4_r8*eta**2.0_r8/lambda
        ksw_calc=MAX( d50,                                              &
     &                d50*(mu+6.0_r8*(theta_timeavg-1.0_r8)) )          &
     &               +ripple_fac
!
        END FUNCTION ksw_calc
!
        REAL FUNCTION fw_calc(ahat, ksw)
        USE utils_constants
!
! Calculate full-cycle wave friction factor from VA2013 Eqn. A.4 
!
        real(r8) :: ahat, ksw, ratio, fw
!
        ratio=ahat/ksw
        IF(ratio.gt.1.587_r8) THEN
          fw_calc=0.00251_r8*EXP(5.21_r8*(ratio)**(-0.19_r8))
        ELSE
          fw_calc=0.3_r8
        ENDIF
!
        END FUNCTION fw_calc 
!
        REAL FUNCTION fd_calc(dsf, ksd)
        USE utils_constants
!
! Calculate current related friction factor VA2013 Eqn. 20
! Assuming logarithmic velocity profile 
! 
        real(r8) :: dsf, ksd
        real(r8), parameter :: von_k=0.41_r8
! 
        fd_calc=2.0_r8*(von_k/LOG(30.0_r8*dsf/ksd))**2.0_r8
!       
        END FUNCTION fd_calc 
!
        REAL FUNCTION fwi_calc(T_iu, T_i, ahat, ksw)
        USE utils_constants
!
! Wave friction factor for wave and crest half cycle VA2013 Eqn. 21
!       
        real(r8) :: T_iu, T_i, ahat, ksw
        real(r8) :: c1, ratio, fwi
!
        c1=2.6_r8
        ratio=ahat/ksw
        IF(ratio.gt.1.587_r8) THEN
          cff=(2.0_r8*T_iu/T_i)**c1
          fwi_calc=0.00251_r8*EXP(5.21_r8*(cff*ratio)**(-0.19_r8))
        ELSE
          fwi_calc=0.3_r8
        ENDIF
!       
        END FUNCTION fwi_calc
!
        REAL FUNCTION dsf_calc(d50, theta_i)
        USE utils_constants
!
! Sheet flow thickness VA2013 Appendix C.1
!       
        real(r8) :: d50, theta_i
        real(r8) :: d50_mm
!
        d50_mm=d50*0.001_r8
        IF(d50_mm.le.0.15_r8)THEN
          cff=25.0_r8*theta_i
        ELSEIF(d50_mm.gt.0.15_r8.and.d50_mm.lt.0.20_r8)THEN 
          cff=25.0_r8-(12.0_r8*(d50_mm-0.15_r8)/0.05_r8)
        ELSEIF(d50_mm.ge.0.20_r8)THEN 
          cff=13.0_r8*theta_i
        ENDIF 
        dsf_calc=MAX(d50*cff,d50)
!
        END FUNCTION dsf_calc
!
!  End of functions for step 2 for Shear stress formulation
!
        SUBROUTINE skewness_params(H_s, T, depth, r, phi, Su, Au)
!        
! Ruessink et al. provides equations for calculating skewness parameters
! Uses Malarkey and Davies equations to get "bb" and "r"
! Given input of H_s, T and depth 
! r     - skewness/asymmetry parameter r=2b/(1+b^2)            [value]
! phi   - skewness/asymmetry parameter                         [value]
! Su     - umax/(umax-umin)                                    [value]
! Au   - amax/(amax-amin)                                      [value]
! alpha - tmax/pi                                              [value]

        use utils_constants
        implicit none 
!
        real, intent(in)  :: H_s, T, depth
        real, intent(out) :: r, phi, Su, Au
!
! Local variables 
! 
        real, parameter :: p1=0.0_r8
        real, parameter :: p2=0.857_r8
        real, parameter :: p3=-0.471_r8
        real, parameter :: p4=0.297_r8
        real, parameter :: p5=0.815_r8
        real, parameter :: p6=0.672_r8
        real :: a_w, Ur
        real :: B, psi, bb 
        real :: k, kh, cff
        real :: kh_calc 
! to be deleted when in  COAWST 
        real ::  dtr,deg2rad,pi
        pi=3.1416_r8                           
        dtr = pi/180.0_r8! % degrees to radians
        deg2rad=dtr   
!
! Ruessink et al., 2012, Coastal Engineering 65:56-63.
!
! k is the local wave number computed with linear wave theory.
!
          k=kh(T,depth)/depth       
!
!          H_s=sqrt(2.0_r8)*H_rms
          a_w=0.5_r8*H_s 
          Ur=0.75_r8*a_w*k/((k*depth)**3.0_r8)
!
! Ruessink et al., 2012 Equation 9.
!
          cff=EXP( (p3-log10(Ur)) /p4)
          B=p1+((p2-p1)/(1.0_r8+cff))
          psi=-90.0_r8*deg2rad*(1.0_r8-TANH(p5/Ur**p6))
! 
! Markaley and Davies, equation provides bb which is "b" in paper
! Check from where CRS found these equations
! 
          bb=sqrt(2.0_r8)*B/(sqrt(2.0_r8*B**2.0_r8+9.0_r8))
          r=2.0_r8*bb/(bb**2.0_r8+1.0_r8)
!
! Ruessink et al., 2012 under Equation 12.
!
          phi=-psi-0.5_r8*pi
!
! Where are these asymmetry Su, Au utilized 
! recreate the asymetry 
!          
          Su=B*cos(psi)
          Au=B*sin(psi)
     
        END SUBROUTINE skewness_params

        SUBROUTINE abreu_points(r, phi, Uw, Sk, Ak, T, DTc, DTt,        &
     &                                      DTcu, DTtu, umax, umin)
! 
!  Calculate umax, umin, and phases of asymmetrical wave orbital velocity 
!
! Why calling it abreu if all equations from Malarkey
!
!  Use the asymmetry parameters from Ruessink et al, 2012
!  to get the umax, umin and phases of asymettrical wave 
!  orbital velocity to be used by Van Der A. 
!  T_c is duration of crest
!  T_cu Duration of accerating flow within crest half cycle
        use utils_constants
!  
        implicit none 
!
        real, intent(in)  :: r, phi, Uw, T
        real, intent(out) :: Sk, Ak
        real, intent(out) :: Dtc, DTt, DTcu, DTtu 
        real, intent(out) :: umax, umin 
!
! Local variables 
! 
        real :: b, c, tmt, tmc, tzd, tzu
        real :: omega, w, phi_new 
        real :: RR, P, F0, betar_0, beta
        real :: T_tu, T_cu, T_c, T_t 
        real :: cff1, cff2, cff 
!
! to be deleted when in  COAWST 
        real ::  dtr,deg2rad,pi
        pi=3.1416_r8                           
        dtr = pi/180.0_r8! % degrees to radians
        deg2rad=dtr   

          omega=2.0_r8*pi/T
!
!
          phi_new=-phi

! Malarkey and Davies (Under equation 16b) 
          P=sqrt(1.0_r8-r*r)
!
! Malarkey and Davies (Under equation 16b) 
!
          b=r/(1.0_r8+P)
!
! Appendix E of Malarkey and Davies 
!
          c=b*sin(phi_new)
!
          cff1=4.0_r8*c*(b*b-c*c)+(1.0_r8-b*b)*(1.0_r8+b*b-2.0_r8*c*c)
          cff2=(1.0_r8+b*b)**2.0_r8-4.0_r8*c*c
          tmc=ASIN(cff1/cff2)
!
          cff1=4.0_r8*c*(b*b-c*c)-(1.0_r8-b*b)*(1.0_r8+b*b-2.0_r8*c*c)
          cff2=(1.0_r8+b*b)**2.0_r8-4.0_r8*c*c
          tmt=ASIN(cff1/cff2)
!           
          IF(tmt.lt.0.0_r8) THEN 
            tmt=tmt+2.0_r8*pi
          ENDIF 
          IF(tmt.lt.0.0_r8) THEN 
            tmc=tmc+2.0_r8*pi
          ENDIF
! 
! Non dimensional umax and umin, under E5 in Malarkey and Davies 
! 
          umax=1.0_r8+c
          umin=umax-2.0_r8
!
!         Dimensionalize
!
          umax=umax*Uw
          umin=umin*Uw
!
! phase of zero upcrossing and downcrossing (radians)
!
          tzu=ASIN(b*SIN(phi_new))
          tzd=2.0_r8*ACOS(c)+tzu 
! 
! MD, equation 17
!
          RR=0.5_r8*(1.0_r8+b*SIN(phi_new))
!
! 
! MD, under equation 18
! 
          IF(r.le.0.5_r8) THEN
            F0=1.0_r8-0.27_r8*(2.0_r8*r)**(2.1_r8)
          ELSE
            F0=0.59_r8+0.14_r8*(2.0_r8*r)**(-6.2_r8)
          ENDIF
!
! MD, Equation 15a,b 
!
          IF(0.0_r8.ge.r.lt.0.5)THEN
            betar_0=0.5_r8*(1.0_r8+r)
          ELSEIF(0.5_r8.le.r.gt.0.5_r8)THEN
            cff1=4.0_r8*r*(1.0_r8+r)
            cff2=cff1+1.0_r8
            betar_0=cff1/cff2
          ENDIF
!
! MD, Equation 18
!
          cff=sin((0.5_r8*pi-ABS(phi_new))*F0)/sin(0.5_r8*pi*F0)
          beta=0.5_r8+(betar_0-0.5_r8)*cff
!
! MD, Table 1, get asymmetry parameterization
! using GSSO (10a,b)
!
          cff=sqrt(2.0_r8*(1.0_r8+b*b)**3.0_r8)
          Sk=3.0_r8*sin(phi_new)/cff
          Ak=-3.0_r8*cos(phi_new)/cff
!
! These are the dimensional fractions of wave periods needed by Van der A eqn.
! TSK - Check source of these equations 
!
          w=1.0_r8/omega
          DTc=(tzd-tzu)*w
          DTt=T-DTc;
          DTcu=(tmc-tzu)*w
          DTtu=(tmt-tzd)*w
!
          T_tu=tzd*w
          T_cu=tzu*w
          T_c=tmc*w
          T_t=tmt*w
!
        END SUBROUTINE abreu_points

        REAL FUNCTION kh(Td,depth)
!
!  Calculate wave number from Wave period and depth 
!
! RL Soulsby (2006) "Simplified calculation of wave orbital velocities"
! HR Wallingford Report TR 155, February 2006
        use utils_constants
! 
        real(r8), parameter :: g=9.81_r8, pi=3.14_r8 
        real(r8) :: Td, depth
        real(r8) :: cff
        real(r8) :: x, y, t, omega 
!
        omega=2.0_r8*pi/Td
        x=omega**2.0_r8*depth/g
!
!TSK Double check the conditional
! 
        IF(x.lt.1.0_r8) THEN
          y=SQRT(x)
        ELSE
          y=x
        ENDIF
!
! Iteratively solving 3 times for eqn.7 of Soulsby 1997 by using 
! eqns. (12a-14)
!      
        t=TANH(y)
        cff=(y*t-x)/(t+y*(1.0_r8-t*t))
        y=y-cff
!
        t=TANH(y)
        cff=(y*t-x)/(t+y*(1.0_r8-t*t))
        y=y-cff

        t=TANH(y)
        cff=(y*t-x)/(t+y*(1.0_r8-t*t))
        y=y-cff
        kh=y
!       
        RETURN 
        END FUNCTION 
!
