
program tromp_1995
    use params, only: Vani, rho_spl, vp_spl, Arad, Crad, Lrad, Nrad, Frad, nprocs, nmodes
    use Integrate, only: integrate_r_traps
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use mesh_utils, only:  delta_spline
    use V_ani, only: integrate_GNIr2, save_Vani_matrix
    use w3j, only: thrj
    use mineos_model, only: mineos, mineos_ptr
    use modes, only: Mode, get_mode
    use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
    implicit none
    include "constants.h"

    integer :: imode, i, j, k, l, m, n, n1, s, q, ispec, myrank, lentrim, iproc, region, m1, m2, knot_lower, knot_upper, is
    character(len=1)   :: t1
    character(len=250)  :: out_name, outfmt
    real(SPLINE_REAL)  :: sum
    integer :: npoints
    real(kind=CUSTOM_REAL) r_lower, r_upper

    real(SPLINE_REAL) :: val, lf, mf, kf
    integer, dimension(5), parameter :: I_n = (/ 5, 3, 3, 1, 1/)
    real(kind=CUSTOM_REAL) :: dA, dC, dL, dN, dF, thirty, twone
    integer :: nl, ll 
    type(Mode)            :: mode_1
    type(InterpPiecewise) :: interp

    ! If we want to test lots: 
    integer, dimension(30), parameter :: modeNs = (/9, 2, 5, 6, 7, 8, 21, 7, 9, 2, 3, 9, 9, 11, 11, 13, 13, 13, 13, 15, 15, 18, 18, 20, 21, 25, 27, 21, 21, 16/)
    integer, dimension(30), parameter :: modeLs = (/3, 3, 3, 3, 4, 5,  7, 5, 2, 3, 2, 3, 4,  4,  5,  1,  2,  3,  6,  3,  4,  3,  4,  1,  6,  2,  2,  8,  6,  7/)
    

    logical, parameter :: tromp_93_model = .true.

    thirty = three * ten 
    twone  = three * seven 


    if(.not.tromp_93_model)then 
        ! constant benchmark
        dA =  0.04d0
        dC = -0.02d0
        dL =  0.03d0
        dN = -0.05d0
        dF =  0.01d0
        write(*,*)"Using constant..."
    endif

        
    ! Read mineos model 
    call mineos%process_mineos_model(.false.)
    mineos_ptr => mineos

    ! Choose a mode:
    do imode = 1, nmodes 
        
        nl = modeNs(imode) !9
        ll = modeLs(imode) !3

        mode_1 = get_mode(nl, 'S', ll, mineos_ptr)

   

        allocate(Vani(mode_1%tl1, mode_1%tl1))
        Vani = SPLINE_iZERO

        ! Values for the inner core
        knot_lower = 1
        r_lower    = zero        
        knot_upper = mineos%disc(2)
        r_upper    = mineos%rdisc(2)
        npoints    = 50*(knot_upper-knot_lower)

        write(*,*)"number of points ", npoints
    
        interp = create_PieceInterp(npoints)
        interp%radial = [((r_lower +  (real(j-1)/real(npoints-1))*(r_upper-r_lower))/scale_R, j = 1, npoints)] 

        call interp%setup()
        call interp%create_interpolation_radial_map()

        ! Interpolate the mode splines
        call interp%interpolate_mode_eigenfunctions(mode_1)

        
        ! Constant value over the radius
        allocate(Arad(npoints))
        allocate(Crad(npoints))
        allocate(Lrad(npoints))
        allocate(Nrad(npoints))
        allocate(Frad(npoints))


        if(tromp_93_model)then 

                ! Need to write the radial spline values to disk so python
                ! can determine the values: 
                write(*,*)'Max radius: ', maxval(interp%radial)
                open(1, file='Vani/Tromp_1993_model/ACLNF/radial/r_for_radmethod', form='formatted')
                do i = 1, npoints
                write(1,*)interp%radial(i)
                enddo 
                close(1)


                write(*,*)"Loading ACLNF from files"
                call load_ACLNF_from_files('/scratch/gpfs/TROMP/we3822/NMSplit90/specfem_mesh/tests/Vani/Tromp_1993_model/ACLNF/radial', & 
                                            npoints, '_0', 0)

                ! Now we load in actual PREM perturbations so we need to non-dimensionalise 
                Arad = Arad/(RHOAV*SCALE_V*SCALE_V)
                Crad = Crad/(RHOAV*SCALE_V*SCALE_V)
                Lrad = Lrad/(RHOAV*SCALE_V*SCALE_V)
                Nrad = Nrad/(RHOAV*SCALE_V*SCALE_V)
                Frad = Frad/(RHOAV*SCALE_V*SCALE_V)

                ! Interpolate the vp and rho to the relevant points
                ! since these need to be abs perturbations not % perturbs
                allocate(rho_spl(npoints))
                allocate(vp_spl(npoints))
    
                call interp%interpolate_mineos_variable(real(mineos%rho_mineos, kind=SPLINE_REAL), rho_spl)
                call interp%interpolate_mineos_variable(real(mineos%vp_mineos,  kind=SPLINE_REAL), vp_spl)

                ! Compute A0 for rel scaling -- kappa + 4*mu 
                !Arad = Arad * (vp_spl*vp_spl)*rho_spl
                !Crad = Crad * (vp_spl*vp_spl)*rho_spl
                !Lrad = Lrad * (vp_spl*vp_spl)*rho_spl
                !Nrad = Nrad * (vp_spl*vp_spl)*rho_spl
                !Frad = Frad * (vp_spl*vp_spl)*rho_spl


        else 
            Arad(:) = dA
            Crad(:) = dC
            Lrad(:) = dL
            Nrad(:) = dN
            Frad(:) = dF
        endif 

        write(*,*)" Min and max of params: "
        write(*,*) minval(Arad), maxval(Arad)
        write(*,*) minval(Crad), maxval(Crad)
        write(*,*) minval(Lrad), maxval(Lrad)
        write(*,*) minval(Nrad), maxval(Nrad)
        write(*,*) minval(Frad), maxval(Frad)


        do m = -mode_1%l, mode_1%l
            mf = real(m, kind=CUSTOM_REAL)

            do s = 0, 4, 2 ! Loop with step of 2 from 0-4
                ! Compute sum_N sum_I int \Gamma_{NI} r^2 dr
                sum = zero
                do N = 0, 4
                    do I = 1, I_n(N+1)
                        sum = sum + integrate_GNIr2(s, mode_1%l, N, I, mode_1%u_spl, & 
                                                    mode_1%du_spl, & 
                                                    mode_1%v_spl/mode_1%kf, mode_1%dv_spl/mode_1%kf, npoints, & 
                                                    interp%radial, Arad, Crad, Lrad, &
                                                    Nrad, Frad, mode_1%t)
                    enddo
                enddo 

                ! Computing D.208 but not including the (2s + 1 / 4pi)^1/2 term 
                ! since that is already added in to the Gamma_NI via the gammaD1_coeff 
                ! function 
                Vani(m+mode_1%l+1,m+mode_1%l+1) = Vani(m+mode_1%l+1,m+mode_1%l+1) + & 
                                                (-SPLINE_ONE)**mf * (two*mode_1%lf + one) * thrj(mode_1%l, s, mode_1%l, -m, 0, m) * sum 
            enddo 
        enddo 

        outfmt = '(a,' !i1,a,i1,a)'
        if(nl.ge.10)then 
            outfmt = trim(outfmt)//'i2'
        else 
            outfmt = trim(outfmt)//'i1'
        endif 
        outfmt = trim(outfmt)//',a,'
        if(ll.ge.10)then 
            outfmt = trim(outfmt)//'i2'
        else 
            outfmt = trim(outfmt)//'i1'
        endif 
        outfmt = trim(outfmt)//',a)'
        write(*,*)"WARNING: WE UPDATED save_VANI - need division by 1/(2omega * SCALE_T) here"


        ! Save the non-dim versoin
        write(out_name, trim(outfmt))'./v_ani_matrix/radial_', mode_1%n, mode_1%t, mode_1%l, '.txt'
        call save_Vani_matrix(mode_1%l, mode_1%l, out_name, .false.)

        write(*,*)"Real freq; ", mode_1%wcom, mode_1%wcom*SCALE_T / (two*PI)


             ! Write the eigenfunctions to disk for reference
        write(out_name, trim(outfmt))'./store_eigen/', mode_1%n, mode_1%t, mode_1%l,"_eigens.txt"
        open(unit=1, file=trim(out_name))
            do i = 1, mineos%NR
                write(1, "(E14.5,E14.5,E14.5,E14.5,E14.5)")mineos%rad_mineos(i), mode_1%u(i), mode_1%v(i), mode_1%du(i), mode_1%dv(i)
            enddo 
        close(1)

        stop 

        deallocate(Vani)
        deallocate(Arad)
        deallocate(Crad)
        deallocate(Lrad)
        deallocate(Nrad)
        deallocate(Frad)
        call deallocate_if_allocated(rho_spl)
        call deallocate_if_allocated(vp_spl)
    enddo

end program tromp_1995



