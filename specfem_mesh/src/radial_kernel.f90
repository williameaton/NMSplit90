! Compute the radial kernels for a mode

program radial_kernel
        use params, only: Vani, nprocs, nmodes, rho_spl, vp_spl, vs_spl
        use v_ani, only: load_vani_from_file, convert_imag_to_real, save_Vani_real_matrix
        use splitting_function, only: get_Ssum_bounds, Hreal_to_cst, write_cst_to_file, Hcomplex_to_cst_8, write_cst_complex_to_file
        use specfem_mesh,       only: SetMesh, create_SetMesh
        use modes,              only: get_mode, Mode 
        use mineos_model,       only: mineos, mineos_ptr
        use ylm_plm,            only: ylm_real    
        use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
        use woodhouse_kernels, only: WK_Vkappa, WK_Vmu
        implicit none 
        include "constants.h"
        character(len=3) :: n1str, l1str, n2str, l2str
        integer   :: l1, n1, l2, n2, tl1, tl2,  & 
                     iproc, i, j, k, ispec,  & 
                     i_mode, npoints, knot_lower, knot_upper 
        character :: t1, t2
        real(kind=CUSTOM_REAL), allocatable :: r_lower, r_upper, omnondim

        character(len=250) :: modefile , out_name

        type(InterpPiecewise) :: interp
        type(Mode)            :: mode_1    
  

        complex(SPLINE_REAL), allocatable :: Vk(:), Vmu(:), As(:), Bs(:)

        !integer, dimension(34), parameter :: modeN1s = (/7, 27, 9, 5, 17, 16, 3, 23, 3, 8, 11, 18, 21, 3, 16, 13, 6, 13, 21, 2,  8,  7, 23, 11, 13, 18, 21,5, 27, 9, 22, 15, 14, 5/)
        !integer, dimension(34), parameter :: modeL1s = (/4, 1,  3, 3, 1,  7,  2, 4,  8, 5, 5,   3, 7,  1,  5,  3, 3,  2,  6,  3, 1,  5,  5,  4,  1,  4,  8,2,  2, 2,  1,  3,  4, 10/)
    
        integer, dimension(4), parameter :: modeN1s = (/ 5, 8, 11, 14/)
        integer, dimension(4), parameter :: modeL1s = (/10, 5, 7, 4/)
    

        ! Read mineos model 
        call mineos%process_mineos_model(.false.)
        mineos_ptr => mineos
    
        ! Values for the whole Earth
        knot_lower = 1
        r_lower    = zero        
        knot_upper = mineos%disc(mineos%ndisc)
        r_upper    = mineos%rdisc(mineos%ndisc)
        npoints    = 10*(knot_upper-knot_lower)
    
        ! Create interpolator with evenly spaced points in IC 
        interp = create_PieceInterp(npoints)
        interp%radial = [((r_lower +  (real(j-1)/real(npoints-1))*(r_upper-r_lower))/scale_R, j = 1, npoints)] 
        call interp%setup()
        call interp%create_interpolation_radial_map()
    

        allocate(rho_spl(npoints))
        allocate(vs_spl(npoints))
        allocate(vp_spl(npoints))

        ! interpolate relevant parameters 
        call interp%interpolate_mineos_variable(real(mineos%rho_mineos, kind=SPLINE_REAL), rho_spl)
        call interp%interpolate_mineos_variable(real(mineos%vp_mineos,  kind=SPLINE_REAL), vp_spl)
        call interp%interpolate_mineos_variable(real(mineos%vs_mineos,  kind=SPLINE_REAL), vs_spl)

        allocate(Vmu(npoints))
        allocate( Vk(npoints))
        allocate( As(npoints))
        allocate( Bs(npoints))


        ! Load from file
    do i_mode = 1,  4
        n1      = modeN1s(i_mode)
        t1      = 'S'
        l1      = modeL1s(i_mode)
        tl1     = l1*2 + 1
        omnondim= SCALE_T * mode_1%wcom

        mode_1 = get_mode(n1, 'S', l1, mineos_ptr)

        call interp%interpolate_mode_eigenfunctions(mode_1)


        call buffer_int(n1str, n1)
        call buffer_int(l1str, l1)


        ! Eqn 100 of Woodhouse & Dahlen 1978 if divided by 2 omega
        call WK_Vkappa(mode_1, mode_1, l1, interp%radial, Vk)
        Vk = Vk/(omnondim*two)

        ! Eqn 102 of Woodhouse & Dahlen 1978 if divided by 2 omega
        call WK_Vmu(mode_1, mode_1, l1, interp%radial, Vmu)
        Vk = Vmu/(omnondim*two)


        As = two * interp%radial * interp%radial * vp_spl * vp_spl * rho_spl * Vk / omnondim

        Bs = two * interp%radial * interp%radial * vs_spl * vs_spl * rho_spl * &
             (Vmu - (four/three)*Vk ) / omnondim

    
        write(out_name, '(a)')'./output/radial_kernels/'//trim(n1str)//mode_1%t//trim(l1str)//'.txt'
        open(1,file=trim(out_name), form='formatted')
        do i = 1, npoints
            write(1,*)interp%radial(i), real(As(i)), real(Bs(i))
        enddo 
        close(1)

        write(out_name, '(a)')'./output/eigenfunctions/'//trim(n1str)//mode_1%t//trim(l1str)//'.txt'
        open(1,file=trim(out_name), form='formatted')
        do i = 1, npoints
            write(1,*)interp%radial(i), real(mode_1%u_spl(i)), real(mode_1%v_spl(i))
        enddo 
        close(1)

    enddo ! Imode 
    
    end program radial_kernel