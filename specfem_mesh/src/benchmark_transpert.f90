! Program to benchmark a transversely isotropic perturbation
! to the inner core

program benchmark_trans_perturb

    use params, only: rho_spl, Wmat, Tell, Vell, Vcen, mu_spl, kappa_spl, nmodes, Vani
    use Integrate, only: integrate_r_traps
    use modes, only: Mode, get_mode
    use mineos_model, only: mineos, mineos_ptr
    use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
    use V_ani, only: integrate_GNIr2, save_Vani_matrix
    use w3j, only: thrj
    use splitting_function, only: get_Ssum_bounds, Hcomplex_to_cst_8, write_cst_complex_to_file

    implicit none
    include "constants.h"

    integer :: i, j, k, l, l2,  m1, m2, npoints, knot_lower, knot_upper , nrad, ier, iline
    character(len=80)  :: out_name
    real(kind=CUSTOM_REAL), allocatable :: r_lower, r_upper, radius(:), & 
                    Aspl(:), Cspl(:), Lspl(:), Nspl(:), Fspl(:), K_C(:), K_A(:), K_L(:), K_N(:), K_F(:)
    real(SPLINE_REAL), allocatable :: integrand(:)
    complex(SPLINE_REAL), allocatable :: W_s(:), W_a(:)

    real(kind=CUSTOM_REAL) :: om2, kf, dw, mf, df_from_c00
    real(SPLINE_REAL) ::  sum
    complex(kind=SPLINE_REAL), allocatable :: cst_imag(:,:)
    complex(kind=SPLINE_REAL) :: c00
    character(len=2)      :: nstr ,lstr
    type(Mode)            :: mode_1
    type(InterpPiecewise) :: interp

    character(len=1) :: t1
    integer :: n1, l1, m, s, n, tl1, smin, smax, num_s, ncols, is, it

    integer, dimension(5), parameter :: I_n = (/ 5, 3, 3, 1, 1/)

    ! 19 rn
    integer, dimension(27), parameter :: modeNs = (/2, 5, 6, 7, 8, 21, 7, 9, 3, 9, 9, 11, 11, 13, 13, 13, 13, 15, 15, 18, 18, 20, 21, 25, 27, 21, 16/)
    integer, dimension(27), parameter :: modeLs = (/3, 3, 3, 4, 5,  7, 5, 2, 2, 3, 4,  4,  5,  1,  2,  3,  6,  3,  4,  3,  4,  1,  6,  2,  2,  8,  7/)
    integer :: imode
    
    logical, parameter :: benchmark_isotropic_perturb = .false.

    ! Read mineos model 
    call mineos%process_mineos_model(.false.)
    mineos_ptr => mineos


    ! READ BENCHMARK PROFILES OF ACLNF and RADIUS: 
    if(benchmark_isotropic_perturb)then 
        open(unit=1,file='benchmarks/isotropic_perturb_ACLNF_profiles.txt', &
        status='old', form='formatted', iostat=ier)
        read(1,*)nrad
        allocate(radius(nrad), Aspl(nrad), Cspl(nrad), Lspl(nrad), Nspl(nrad), Fspl(nrad)) 
        do iline = 1, nrad
            read(1,*) radius(iline), Aspl(iline), Cspl(nrad), Lspl(iline), Nspl(nrad), Fspl(iline)
        enddo 
        close(1)
        ! Non-dimensionalise the ACLNF: 
        Aspl = Aspl/(SCALE_V*SCALE_V*RHOAV)
        Cspl = Cspl/(SCALE_V*SCALE_V*RHOAV)
        Lspl = Lspl/(SCALE_V*SCALE_V*RHOAV)
        Nspl = Nspl/(SCALE_V*SCALE_V*RHOAV)
        Fspl = Fspl/(SCALE_V*SCALE_V*RHOAV)
    else 
        open(unit=1,file='benchmarks/isoradial_perturb_ACLNF_profiles.txt', &
        status='old', form='formatted', iostat=ier)
        read(1,*)nrad
        allocate(radius(nrad), Aspl(nrad), Cspl(nrad), Lspl(nrad), Nspl(nrad), Fspl(nrad)) 

        do iline = 1, nrad
            read(1,*) radius(iline), Aspl(iline), Cspl(iline), Lspl(iline), Nspl(iline), Fspl(iline)
        enddo 
        close(1)
        ! Non-dimensionalise the ACLNF: 
        Aspl = Aspl/(SCALE_V*SCALE_V*RHOAV)
        Cspl = Cspl/(SCALE_V*SCALE_V*RHOAV)
        Lspl = Lspl/(SCALE_V*SCALE_V*RHOAV)
        Nspl = Nspl/(SCALE_V*SCALE_V*RHOAV)
        Fspl = Fspl/(SCALE_V*SCALE_V*RHOAV)
    endif 
    

    knot_lower = 1
    r_lower    = zero        
    knot_upper = mineos%disc(mineos%ndisc)
    r_upper    = mineos%rdisc(mineos%ndisc)
    npoints    = nrad

    ! Create interpolator with evenly spaced points in IC 
    interp = create_PieceInterp(npoints)
    interp%radial = radius
    call interp%setup()
    call interp%create_interpolation_radial_map()


    allocate(K_A(nrad))
    allocate(K_C(nrad))
    allocate(K_L(nrad))
    allocate(K_N(nrad))
    allocate(K_F(nrad))
    allocate(integrand(nrad))

    do imode = 1, nmodes
        ! METHOD FROM CHAPTER 9 of DT 98 
        n1 = modeNs(imode)
        l1 = modeLs(imode)

        tl1 = 2*l1 + 1
        allocate(Vani(tl1, tl1))
        Vani = SPLINE_iZERO

        mode_1 = get_mode(n1, 'S', l1 , mineos_ptr)
        call interp%interpolate_mode_eigenfunctions(mode_1)
        
        ! Non dimensional eigenfreq x 2
        om2 = mode_1%wcom * SCALE_T * two 
        kf  = mode_1%kf
        ! 9.31
        K_C = ((radius * mode_1%du_spl)**two)/om2
        
        ! 9.31 - note we are not using V/k as v, its kV (capital)
        K_A =  ((two*mode_1%u_spl - kf*mode_1%v_spl)**two)/om2

        ! 9.32 
        K_L = ((radius*mode_1%dv_spl -  mode_1%v_spl + kf*mode_1%u_spl)**two)/om2

        ! 9.33 
        K_N = (-(two*mode_1%u_spl - kf*mode_1%v_spl)**two + & 
                (kf**two - two)*(mode_1%v_spl**two))/om2
        ! 9.34
        K_F = (two*radius * mode_1%du_spl * (two*mode_1%u_spl - kf*mode_1%v_spl))/om2

        ! form the integrand, ignoring perturbations in d and rho:
        integrand = Aspl*K_A + Cspl*K_C + Lspl*K_L + Nspl*K_N + Fspl*K_F

        ! Integrate and re-dimensionalise 
        dw =  integrate_r_traps(interp%radial, integrand, nrad)/SCALE_T

        ! convert to microHz:
        dw = 1.0e6 * dw/(two*PI)

        if(benchmark_isotropic_perturb)then 
            ! NOW LETS MATCH IT WITH TROMP 1993 METHOD
            do m = -l1, l1
                mf = real(m, kind=CUSTOM_REAL)
        
                do s = 0, 4, 2 ! Loop with step of 2 from 0-4
                    ! Compute sum_N sum_I int \Gamma_{NI} r^2 dr
                    sum = zero
                    do N = 0, 4
                        do I = 1, I_n(N+1)
                            sum = sum + integrate_GNIr2(s, l1, N, I, mode_1%u_spl, mode_1%du_spl, mode_1%v_spl/kf, & 
                                                        mode_1%dv_spl/kf, nrad, radius, Aspl, Cspl, Lspl, Nspl, Fspl, 'S')
                        enddo
                    enddo 
        
                    ! Computing D.208 but not including the (2s + 1 / 4pi)^1/2 term 
                    ! since that is already added in to the Gamma_NI via the gammaD1_coeff 
                    ! function 
                    Vani(m+l1+1,m+l1+1) = Vani(m+l1+1,m+l1+1) + (-SPLINE_ONE)**mf * (two*mode_1%lf + one) * thrj(l1, s, l1, -m, 0, m) * sum 
                enddo 
            enddo 

            ! Next we need to do 1/2omega Vani (non dim)
            Vani = Vani/(two * mode_1%wcom * SCALE_T) * 1.0e6

            ! Compute the Csts
            call get_Ssum_bounds(l1, l1, smin, smax, num_s, ncols)
            allocate(cst_imag(num_s, ncols))
            call Hcomplex_to_cst_8(Vani, l1, l1, cst_imag, ncols, num_s, 'S', 'S', 1)

            ! To get the centre frequency shift from C00 we have 
            ! df = Re(C00)/(4*pi)**0.5
            ! first we re-dimensionalise and convert to Hz
            df_from_c00 = real(cst_imag(1,1)/(two*PI*SCALE_T)) ! 
            df_from_c00 = df_from_c00 / (four*PI)**half
            deallocate(cst_imag)

            write(*,'(i2, a, i3, a, f15.10, f15.10)') n1, ' S ', l1, ' ', dw, df_from_c00

        else
            ! CANT USE TROMP 93 for radial TI model so print Ch 9 method only
            write(*,'(i2, a, i2, f15.10, f15.10, f15.10, f15.10)') n1, ' S ', l1,  dw
        endif 

        deallocate(Vani)

    enddo ! imode

end program