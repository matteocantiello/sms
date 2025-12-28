! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use auto_diff
      
      implicit none
      
      include "test_suite_extras_def.inc"

      
      contains

      include "test_suite_extras.inc"
      include "other_cgrav.inc"
      
      
      subroutine extras_controls(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            
            s% extras_startup => extras_startup
            s% extras_check_model => extras_check_model
            s% extras_finish_step => extras_finish_step
            s% extras_after_evolve => extras_after_evolve
            s% how_many_extra_history_columns => how_many_extra_history_columns
            s% data_for_extra_history_columns => data_for_extra_history_columns
            s% how_many_extra_profile_columns => how_many_extra_profile_columns
            s% data_for_extra_profile_columns => data_for_extra_profile_columns  
            s% other_cgrav => my_other_cgrav
            s% other_accreting_state => sms_accretion
         
         end subroutine extras_controls


         subroutine sms_accretion(id, total_specific_energy, accretion_pressure, accretion_density, ierr)
         ! ------------------------------------------------------------------------------
         ! SMS COLD ACCRETION ROUTINE (v1.0 // M.Cantiello // 12.2025)
         !
         ! Purpose: 
         !   Overrides MESA's default accretion boundary conditions to model "cold" 
         !   accretion onto Supermassive Stars (SMS).
         !
         ! Physics:
         !   1. Calculates the density and pressure of a free-falling accretion stream
         !      coming from infinity (Ram Pressure boundary condition).
         !      - v_ff = sqrt(2GM/R)
         !      - P_ram = rho * v_ff^2
         !   2. Forces the accreted material to have a fixed, low temperature (T_target),
         !      reducing the specific entropy of the infalling mass.
         !   3. Uses MESA's EOS to determine the specific internal energy of the accreted material
         !         based on the calculated density and target temperature.
         !  
         ! Inlist Inputs :
         !   - s% x_ctrl(1): Accretion Rate in Msun/
         !   - s% x_ctrl(2): Target Temperature of Accreted Material in K        
         ! ------------------------------------------------------------------------------
 
            use star_def
            use const_def
            use eos_def, only: i_lnE, num_eos_basic_results
            use eos_lib, only: eosDT_get
            
            integer, intent(in) :: id
            real(dp), intent(out) :: total_specific_energy, accretion_pressure, accretion_density
            !real(dp) :: total_specific_energy, accretion_pressure, accretion_density

            integer, intent(out) :: ierr
            
            type(star_info), pointer :: s
            
            ! Physics variables
            real(dp) :: M_g, R_cm, Mdot_g_s, v_ff
            real(dp) :: T_target, rho_calc, P_ram
            
            ! EOS helper variables
            real(dp) :: log_rho, log_T
            
            ! EOS Output Arrays
            real(dp) :: res(num_eos_basic_results)
            real(dp) :: dres_dlnd(num_eos_basic_results)
            real(dp) :: dres_dlnT(num_eos_basic_results)
            real(dp), allocatable :: dres_dxa(:,:) 
            

            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return

            ! Initialize Outputs
            total_specific_energy = 0d0
            accretion_pressure = 0d0
            accretion_density = 0d0

            ! Gather Star Parameters (CGS)
            M_g = s% star_mass * Msun
            R_cm = s% r(1)          
            Mdot_g_s = s% x_ctrl(1) * Msun / secyer ! Input Accretion Rate in Msun/yr
            
            

            ! --- Safety / Fallback Block ---
            if (R_cm < 1.0d5 .or. Mdot_g_s <= 0.0_dp) then
               allocate(dres_dxa(num_eos_basic_results, s% species))
               
            call eosDT_get( &
                  s% eos_handle, &
                  s% species, s% chem_id, s% net_iso, s% xa(:,1), & ! Chemistry First
                  s% rho(1), s% lnd(1)/ln10, &                      ! Density
                  s% T(1), s% lnT(1)/ln10, &                        ! Temperature
                  res, dres_dlnd, dres_dlnT, dres_dxa, ierr)
                           
               total_specific_energy = exp(res(i_lnE))
               accretion_pressure = s% Peos(1)
               accretion_density = s% rho(1)

               deallocate(dres_dxa)
               return
            end if
            
            ! --- Physics Calculations (Cold Stream) ---
            
            ! 1. Free-fall velocity
            v_ff = sqrt(2.0_dp * standard_cgrav * M_g / R_cm)

            ! 2. Continuity Density
            rho_calc = Mdot_g_s / (4.0_dp * pi * R_cm**2 * v_ff)
            
            ! 3. Ram Pressure
            P_ram = rho_calc * v_ff**2
            
            ! 4. Target Temperature (Cold)
            T_target =  s% x_ctrl(2) !1.0d2 Temperature of Accreted Material (in K)

            ! 5. EOS Lookup
            log_rho = log10(rho_calc)
            log_T   = log10(T_target)
            
            allocate(dres_dxa(num_eos_basic_results, s% species))

         ! Call the public EOS routine 
            call eosDT_get( &
               s% eos_handle, &              
               s% species, s% chem_id, s% net_iso, s% xa(:,1), & ! Chemistry Inputs First
               rho_calc, log_rho, &                              ! Density Inputs
               T_target, log_T, &                                ! Temperature Inputs
               res, dres_dlnd, dres_dlnT, &                      ! Outputs
               dres_dxa, ierr)

            if (ierr /= 0) then
               write(*,*) 'Failed in sms_accretion eosDT_get'
               deallocate(dres_dxa) 
               return
            end if

            ! 6. Final Assignments
            accretion_density = rho_calc
            accretion_pressure = P_ram
            
            ! Exponentiate i_lnE to get specific energy in erg/g
            total_specific_energy = exp(res(i_lnE)) 
            
            deallocate(dres_dxa)

            write(*,*) 'Accretion properties (e, p_ram, rho_acc): ', total_specific_energy, accretion_pressure,accretion_density


         end subroutine sms_accretion




      subroutine check_tov_correction(id, s, ierr)
      ! ------------------------------------------------------------------------------
      ! TOV CORRECTION DIAGNOSTIC (v1.0 // M.Cantiello // 12.2025)
      !
      ! Purpose:
      !   Verifies that the General Relativistic gravity corrections (Tolman-Oppenheimer-Volkoff)
      !   are being correctly applied to the stellar model. It compares the theoretical 
      !   TOV enhancement factor against the effective gravity currently used by MESA.
      !
      ! Physics:
      !   Calculates the theoretical TOV factor (Ratio of GR Gravity / Newtonian Gravity):
      !   f_TOV = (1 + P/rho*c^2) * (1 + 4*pi*r^3*P / m*c^2) * (1 - 2Gm/rc^2)^-1
      !
      !   It then compares this to the ratio (s% c_grav / G_standard). Since MESA's 
      !   'other_gravity' hook modifies the effective gravitational constant G_eff, 
      !   these two values must match for the physics to be consistent.
      ! ------------------------------------------------------------------------------

         use star_def
         use const_def, only: clight, standard_cgrav
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         
         integer :: k
         real(dp) :: r, m, P, rho, tov_factor
         real(dp) :: c2, G_const, G_mesa

         ierr = 0
         c2 = clight**2
         G_const = standard_cgrav

         ! Loop over zones (MESA indexes surface=1 to center=s% nz)
         ! We skip the center (k=s% nz) to avoid division by zero if r=0
         do k = 1, s% nz - 1
            
            ! Variables at cell face k (between cell k and k-1? No, face k is outer boundary of cell k)
            ! Note: MESA variables are often cell-centered. c_grav is face-defined.
            ! We need to interpolate P and rho to the face k or use the values carefully.
            ! For a quick check, using cell k values for P/rho is usually sufficient 
            ! or simple averaging: P_face = 0.5*(s% P(k) + s% P(k-1))
            
            r = s% r(k)      ! Radius at face k
            m = s% m(k)      ! Mass enclosed within face k
            
            ! Approximate P and rho at face k (simple average)
            
            if (k > 1) then
               P = 0.5_dp * (s% Peos(k-1) + s% Peos(k))
               rho = 0.5_dp * (s% rho(k-1) + s% rho(k))
            else
               P = s% Peos(k)
               rho = s% rho(k)
            end if
            
            ! 1. Calculate TOV Factor terms
            ! Term 1: (1 + P / rho c^2) - Enthalpy correction
            ! Term 2: (1 + 4pi r^3 P / m c^2) - Pressure as source of gravity
            ! Term 3: (1 - 2Gm / r c^2)^-1 - Metric correction (Schwarzschild geometry)
            
            tov_factor = (1.0_dp + P / (rho * c2)) * &
                        (1.0_dp + (4.0_dp * 3.1415926535_dp * r**3 * P) / (m * c2)) * &
                        (1.0_dp - (2.0_dp * G_const * m) / (r * c2))**(-1)

           
            G_mesa = s% cgrav(k)  ! This is G_eff (~10^-8)

            ! 2. Print comparison for a few zones (e.g., deep inside where GR matters)
            if (mod(k, 100) == 0 .and. m > 0.1_dp*s% mstar) then
               write(*, '(A, I4, 4ES14.6)') 'TOV Check k=', k, &
                      tov_factor, G_mesa/G_const
               
               ! If ratio is ~1.0, TOV is NOT active or negligible.
               ! For SMS, tov_factor should be > 1.0 (e.g., 1.001 - 1.1 depending on compactness)
            end if
            
         end do

      end subroutine check_tov_correction

      
      ! subroutine check_gr_instability(id, s, ierr)
      ! ! ------------------------------------------------------------------------------
      ! ! GR INSTABILITY MONITOR (v1.0 // M.Cantiello // 12.2025)
      ! !
      ! ! Purpose:
      ! !   Monitors the stability of Supermassive Stars against General Relativistic 
      ! !   radial instability. It computes the Chandrasekhar (1964) stability integral
      ! !   to detect the onset of collapse.
      ! !
      ! ! Physics:
      ! !   1. Calculates the stability integral I, which includes:
      ! !      - Newtonian Stiffness Term: Integral[ (Gamma1 - 4/3) * (P/rho) dm ]
      ! !      - GR Destabilizing Term: Integral[ K * (2Gm/rc^2) * (P/rho) dm ]
      ! !      If I < 0, the fundamental radial mode frequency is imaginary (unstable).
      ! !
      ! !   2. Computes global diagnostics for comparison with literature (e.g. Nagele+2022):
      ! !      - Pressure-weighted average <Gamma1>
      ! !      - Global Compactness parameter (2GM / Rc^2)
      ! !      - Critical Gamma threshold (Gamma_crit = 4/3 + 1.12 * Compactness)
      ! ! ------------------------------------------------------------------------------   

      !    use star_def
      !    use const_def
      !    integer, intent(in) :: id
      !    type (star_info), pointer :: s
      !    integer, intent(out) :: ierr
         
      !    integer :: k
      !    real(dp) :: dm, r_outer, r_inner, r_mid, m_mid, P, rho, gam1
      !    real(dp) :: term_newton, term_gr, integral_val
      !    real(dp) :: c2, G_const, compactness, local_comp
         
      !    ! Variables for Global Averages
      !    real(dp) :: sum_weights, sum_weighted_gam, avg_gamma, crit_gamma, global_compactness
      !    real(dp) :: current_R
         
      !    ierr = 0
      !    c2 = clight**2
      !    G_const = standard_cgrav

      !    term_newton = 0.0_dp
      !    term_gr = 0.0_dp
      !    sum_weights = 0.0_dp
      !    sum_weighted_gam = 0.0_dp
         
      !    do k = 1, s% nz
            
      !       ! 1. Load Cell Properties
      !       dm = s% dm(k)
      !       P = s% Peos(k)      
      !       rho = s% rho(k)     
      !       gam1 = s% gamma1(k) 
            
      !       ! 2. Geometry Centering
      !       ! 's% dr' does not exist in some versions. We calculate r_mid manually.
      !       ! MESA indexing: k=1 (Surface) -> k=nz (Center)
      !       r_outer = s% r(k)
            
      !       if (k < s% nz) then
      !          r_inner = s% r(k+1)
      !       else
      !          r_inner = 0.0_dp
      !       end if
            
      !       r_mid = 0.5_dp * (r_outer + r_inner)
      !       m_mid = s% m(k) - 0.5_dp * dm 
            
      !       if (r_mid < 1.0d5) cycle 
            
      !       ! 3. Instability Integral Terms
      !       ! Newtonian Stiffness
      !       term_newton = term_newton + (gam1 - 4.0_dp/3.0_dp) * (P/rho) * dm
            
      !       ! GR Correction (local compactness at this shell)
      !       local_comp = (2.0_dp * G_const * m_mid) / (r_mid * c2)
      !       term_gr = term_gr + 1.12_dp * local_comp * (P/rho) * dm

      !       ! 4. Accumulate Averages for Plotting
      !       ! Weighting factor w = P/rho * dm (proportional to internal energy)
      !       sum_weights = sum_weights + (P/rho) * dm
      !       sum_weighted_gam = sum_weighted_gam + gam1 * (P/rho) * dm

      !    end do
         
      !    ! --- Instability Integral Check ---
      !    integral_val = term_newton - term_gr
         
      !    ! Use a safe extra column (ensure s% xtra1 is defined in your defaults, or just print)
      !    s% xtra(1) = integral_val 

      !    ! --- Global Parameters for Comparison with Papers ---
      !    if (sum_weights > 0.0_dp) then
      !       avg_gamma = sum_weighted_gam / sum_weights
      !    else
      !       avg_gamma = 0.0_dp
      !    end if

      !    ! Global Compactness (2GM / Rc^2)
      !    ! Use s% r(1) for the surface radius (R_phot is not always available)
      !    current_R = s% r(1)
         
      !    if (current_R > 0.0_dp) then
      !       global_compactness = (2.0_dp * G_const * s% star_mass * Msun) / (current_R * Rsun * c2)
      !    else
      !       global_compactness = 0.0_dp
      !    end if

      !    ! Approximate Critical Gamma (Gamma_crit = 4/3 + 1.12 * Compactness)
      !    crit_gamma = (4.0_dp/3.0_dp) + 1.12_dp * global_compactness

      !    ! Print Status
      !    ! Only print every 10 models or if unstable
      !    if (mod(s% model_number, 10) == 0 .or. integral_val < 0.0_dp) then
      !       write(*, '(A, I6)') '--- GR Monitor Model: ', s% model_number
      !       write(*, '(A, ES12.5, A, ES12.5)') '   Integral Val: ', integral_val, &
      !                                           ' (Neg = Unstable)'
      !       write(*, '(A, F10.6, A, F10.6)')   '   <Gamma_1>:    ', avg_gamma, &
      !                                           ' vs Crit: ', crit_gamma
      !       write(*, '(A, ES12.5)')            '   Compactness:  ', global_compactness
            
      !       if (avg_gamma < crit_gamma) then
      !          write(*,*) '   >>> WARNING: Gamma < Gamma_crit (Instability approaching)'
      !       end if
      !       write(*,*) '--------------------------------'
      !    end if

      ! end subroutine check_gr_instability


      subroutine check_gr_instability(id, s, ierr)
         
         ! GR INSTABILITY MONITOR
         ! ... (Header description we wrote earlier) ...

         use star_def
         use const_def
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         
         integer :: k
         real(dp) :: dm, P, rho, gam1
         real(dp) :: term_newton, term_gr, integral_val
         real(dp) :: c2, G_const
         
         ! Variables for geometry in CGS
         real(dp) :: r_outer_cm, r_inner_cm, r_mid_cm
         real(dp) :: m_mid_grams
         real(dp) :: local_comp
         
         ! Variables for Global Averages
         real(dp) :: sum_weights, sum_weighted_gam, avg_gamma, crit_gamma, global_compactness
         real(dp) :: current_R_cm, current_M_g
         
         ierr = 0
         c2 = clight**2
         G_const = standard_cgrav

         term_newton = 0.0_dp
         term_gr = 0.0_dp
         sum_weights = 0.0_dp
         sum_weighted_gam = 0.0_dp
         
         do k = 1, s% nz
            
            ! 1. Load Cell Properties (CGS standard in MESA)
            dm   = s% dm(k)     ! Grams
            P    = s% Peos(k)   ! dyn/cm^2
            rho  = s% rho(k)    ! g/cm^3
            gam1 = s% gamma1(k) ! dimensionless
            
            ! 2. Geometry Conversion (Solar Units -> CGS)
            ! s% r is in Rsun. Convert to cm.
            r_outer_cm = s% r(k) * Rsun
            
            if (k < s% nz) then
               r_inner_cm = s% r(k+1) * Rsun
            else
               r_inner_cm = 0.0_dp
            end if
            
            r_mid_cm = 0.5_dp * (r_outer_cm + r_inner_cm)
            
            ! s% m is in Msun. Convert to grams.
            ! Note: dm is ALREADY in grams.
            ! m(k) is mass enclosed by outer boundary of cell k.
            m_mid_grams = (s% m(k) * Msun) - 0.5_dp * dm 
            
            ! Safety check for center
            if (r_mid_cm < 1.0d5) cycle 
            
            ! 3. Instability Integral Terms
            ! weighting factor dV_P = (P/rho) * dm has units of Energy (ergs)
            
            ! Newtonian Stiffness (Dimensionless coefficient * ergs)
            term_newton = term_newton + (gam1 - 4.0_dp/3.0_dp) * (P/rho) * dm
            
            ! GR Correction
            ! Compactness = (2 G M / R c^2) is dimensionless.
            ! We MUST use CGS for G, M, R, c here.
            local_comp = (2.0_dp * G_const * m_mid_grams) / (r_mid_cm * c2)
            
            term_gr = term_gr + 1.12_dp * local_comp * (P/rho) * dm

            ! 4. Accumulate Averages
            sum_weights = sum_weights + (P/rho) * dm
            sum_weighted_gam = sum_weighted_gam + gam1 * (P/rho) * dm

         end do
         
         ! --- Instability Integral Check ---
         integral_val = term_newton - term_gr
         
         ! Store in extra column (ensure you use xtra1 or xtra(1) depending on MESA ver)
         s% xtra(1) = integral_val 

         ! --- Global Parameters for Comparison ---
         if (sum_weights > 0.0_dp) then
            avg_gamma = sum_weighted_gam / sum_weights
         else
            avg_gamma = 0.0_dp
         end if

         ! Global Compactness (2GM / Rc^2)
         current_R_cm = s% r(1) * Rsun
         current_M_g  = s% star_mass * Msun
         
         if (current_R_cm > 0.0_dp) then
            global_compactness = (2.0_dp * G_const * current_M_g) / (current_R_cm * c2)
         else
            global_compactness = 0.0_dp
         end if

         ! Approximate Critical Gamma
         crit_gamma = (4.0_dp/3.0_dp) + 1.12_dp * global_compactness

         ! Print Status
         if (mod(s% model_number, 10) == 0 .or. integral_val < 0.0_dp) then
            write(*, '(A, I6)') '--- GR Monitor Model: ', s% model_number
            write(*, '(A, ES12.5, A, ES12.5)') '   Integral Val: ', integral_val, &
                                                ' (Neg = Unstable)'
            write(*, '(A, F10.6, A, F10.6)')   '   <Gamma_1>:    ', avg_gamma, &
                                                ' vs Crit: ', crit_gamma
            write(*, '(A, ES12.5)')            '   Compactness:  ', global_compactness
            
            if (avg_gamma < crit_gamma) then
               write(*,*) '   >>> WARNING: Gamma < Gamma_crit (Instability approaching)'
            end if
            write(*,*) '--------------------------------'
         end if

      end subroutine check_gr_instability


      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
      end subroutine extras_startup
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve
      

      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !call check_tov_correction(id, s, ierr)  ! Uncomment to test TOV correction diagnostics
         call check_gr_instability(id, s, ierr)   ! Monitor GR Instability


         extras_check_model = keep_going         
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 1
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         names(1) = 'GI_integral_val'   ! Save in history the GR Instability Integral value
         vals(1)  =  s% xtra(1)

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
      end function extras_finish_step
      
      

      end module run_star_extras
      
