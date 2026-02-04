
program radial_perturb
    ! Computes splitting for a radial anisotropy perturbation 
    ! a la Lythgoe and Deuss 
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

    integer :: imode, i, j, k, l, m, n, n1, s, q, ispec, lentrim, iproc, region, m1, m2, knot_lower, knot_upper, is, nradpts, ier, iline
    character(len=1)   :: t1
    character(len=250)  :: out_name, outfmt
    real(SPLINE_REAL)  :: sum
    real(kind=CUSTOM_REAL) r_lower, r_upper

    real(SPLINE_REAL) :: val, lf, mf, kf
    integer, dimension(5), parameter :: I_n = (/ 5, 3, 3, 1, 1/)
    real(kind=CUSTOM_REAL) :: dA, dC, dL, dN, dF, thirty, twone
    integer :: nl, ll 
    type(Mode)            :: mode_1
    type(InterpPiecewise) :: interp

    real(kind=CUSTOM_REAL), allocatable :: radius(:)

    ! If we want to test lots: 
    integer, dimension(29), parameter :: modeNs = (/2, 5, 6, 7, 8, 21, 7, 9, 2, 3, 9, 9, 11, 11, 13, 13, 13, 13, 15, 15, 18, 18, 20, 21, 25, 27, 21, 21, 16/)
    integer, dimension(29), parameter :: modeLs = (/3, 3, 3, 4, 5,  7, 5, 2, 3, 2, 3, 4,  4,  5,  1,  2,  3,  6,  3,  4,  3,  4,  1,  6,  2,  2,  8,  6,  7/)
    

    ! WARNING THIS WAS WRITTEN WRONGLY FOR THE CASE OF TRANSVERSE ISOTROPY BUT THIS IS STILL a VTI calculation
    ! SEE BENCHMARK TRANSPERT 

    thirty = three * ten 
    twone  = three * seven 


    ! Load the radial ACLNF parameters: 
    open(unit=1,file='./Vani/radial_perturb_model.txt', &
        status='old', form='formatted', iostat=ier)
    read(1,*)nradpts

    ! Constant value over the radius
    allocate(radius(nradpts), Arad(nradpts), Crad(nradpts), Lrad(nradpts), Nrad(nradpts), Frad(nradpts)) 
    do iline = 1, nradpts
        read(1,*) radius(iline), Arad(iline), Crad(iline), Lrad(iline), Nrad(iline), Frad(iline)
    enddo 


    ! We need to non-dim the ACLNF: 

    Arad = Arad/(RHOAV*SCALE_V*SCALE_V)
    Crad = Crad/(RHOAV*SCALE_V*SCALE_V)
    Lrad = Lrad/(RHOAV*SCALE_V*SCALE_V)
    Nrad = Nrad/(RHOAV*SCALE_V*SCALE_V)
    Frad = Frad/(RHOAV*SCALE_V*SCALE_V)



    allocate(rho_spl(nradpts))
    allocate(vp_spl(nradpts))
        
    ! Read mineos model 
    call mineos%process_mineos_model(.false.)
    mineos_ptr => mineos


   interp = create_PieceInterp(nradpts)
    interp%radial = radius



    call interp%setup()
    call interp%create_interpolation_radial_map()




    ! Choose a mode:
    do imode = 1, nmodes 
        
        nl = modeNs(imode) !9
        ll = modeLs(imode) !3

        mode_1 = get_mode(nl, 'S', ll, mineos_ptr)

        allocate(Vani(mode_1%tl1, mode_1%tl1))
        Vani = SPLINE_iZERO

        ! Values for the inner core


     
        ! Interpolate the mode splines
        call interp%interpolate_mode_eigenfunctions(mode_1)


   
        call interp%interpolate_mineos_variable(real(mineos%rho_mineos, kind=SPLINE_REAL), rho_spl)
        call interp%interpolate_mineos_variable(real(mineos%vp_mineos,  kind=SPLINE_REAL), vp_spl)
    


        do m = -mode_1%l, mode_1%l
            mf = real(m, kind=CUSTOM_REAL)

            do s = 0, 4, 2 ! Loop with step of 2 from 0-4
                ! Compute sum_N sum_I int \Gamma_{NI} r^2 dr
                sum = zero
                do N = 0, 4
                    do I = 1, I_n(N+1)
                        sum = sum + integrate_GNIr2(s, mode_1%l, N, I, mode_1%u_spl, & 
                                                    mode_1%du_spl, & 
                                                    mode_1%v_spl/mode_1%kf, mode_1%dv_spl/mode_1%kf, nradpts, & 
                                                    interp%radial, Arad, Crad, Lrad, Nrad, Frad, mode_1%t)
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

        write(out_name, trim(outfmt))'./v_ani_matrix/radial_', mode_1%n, mode_1%t, mode_1%l, '.txt'
        call save_Vani_matrix(mode_1%l, mode_1%l, out_name)



        deallocate(Vani)

    enddo

end program radial_perturb



