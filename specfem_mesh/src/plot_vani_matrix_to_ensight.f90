program plot_vani_to_ensight
! program loads in a Vani matrix (complex) and outputs the splitting functions
    use params, only: Vani, nprocs, nmodes
    use v_ani, only: load_vani_from_file, convert_imag_to_real, save_Vani_real_matrix
    use splitting_function, only: get_Ssum_bounds, Hreal_to_cst, write_cst_to_file, Hcomplex_to_cst_8, write_cst_complex_to_file
    use specfem_mesh,       only: SetMesh, create_SetMesh
    use modes,              only: get_mode, Mode 
    use mineos_model,       only: mineos, mineos_ptr
    use ylm_plm,            only: ylm_real    

    implicit none 
    include "constants.h"
    character(len=3) :: n1str, l1str, n2str, l2str
    integer   :: l1, n1, l2, n2, tl1, tl2, smin, smax, num_s, ncols, nrows, is, it, & 
                 iproc, i, j, k, ispec, s, t, ilat, ilon, nlat, nlon, i_mode
    character :: t1, t2
    character(len=20) :: model_ti
    character(len=250) :: modefile , out_name
    real(kind=CUSTOM_REAL)  :: dlat, dlon, theta, phi
    real(kind=SPLINE_REAL)  :: sum
    real(kind=SPLINE_REAL), allocatable :: Vani_real(:,:)
    real(kind=SPLINE_REAL), allocatable :: sig_st(:,:)
    complex(kind=SPLINE_REAL), allocatable :: cst_imag(:,:)
    real(kind=CUSTOM_REAL), allocatable :: sigma(:), ylm_global(:)
    type(SetMesh) :: sm
    logical, parameter :: vti_model = .false.
    logical, parameter :: output_to_ensight = .false.
    logical, parameter :: output_to_evengrid = .false.

    logical, parameter :: cross_coupled = .true.


    ! Modes: 
    ! integer, dimension(20), parameter :: modeNs = (/2, 3, 9, 9, 9, 11, 11, 13,13,13,13,15,15,18,18,20,21,25,27, 6/)
    ! integer, dimension(20), parameter :: modeLs = (/3, 2, 2, 3, 4,  4,  5,  1, 2, 3, 6, 3, 4, 3, 4, 1, 6,2,2, 10/)
    
    !integer, dimension(40), parameter :: modeNs = (/2, 3, 3, 5, 6, 8, 8, 9, 9, 9, 11, 11, 11, 11, 13, 13, 13, 13, 14, 15, 15, 16, 16, 16, 17, 17, 18, 18, 18, 20, 20, 21, 21, 21, 22, 23, 23, 25, 25, 27/)
    !integer, dimension(40), parameter :: modeLs = (/3, 1, 2, 2, 3, 1, 5, 2, 3, 4,  1,  4,  5,  6, 1,  2,  3,  6,   4,  3,  4,  5,  6,  7,  1,  8,  3,  4,  6,  1,  5,  6,  7,  8,  1,  4,  5,  1,  2,  2/)
    !integer, dimension(27), parameter :: modeNs = (/7, 27, 9, 5, 17, 16, 3, 23, 3, 8, 11, 18, 21, 3, 16, 13, 6, 13, 21, 2,  8,  7, 23, 11, 13, 18, 21/)
    !integer, dimension(27), parameter :: modeLs = (/4, 1,  3, 3, 1,  7,  2, 4,  8, 5, 5,   3, 7,  1,  5,  3, 3,  2,  6,  3, 1,  5,  5,  4,  1,  4,  8/)

    integer, dimension(1), parameter :: modeN1s = (/16/)
    integer, dimension(1), parameter :: modeL1s = (/5/)

    integer, dimension(1), parameter :: modeN2s = (/17/)
    integer, dimension(1), parameter :: modeL2s = (/4/)


    ! integer, dimension(33), parameter :: modeN1s = (/7, 27, 9, 5, 17, 16, 3, 23, 3, 8, 11, 18, 21, 3, 16, 13, 6, 13, 21, 2,  8,  7, 23, 11, 13, 18, 21,5, 27, 9, 22, 15, 14/)
    ! integer, dimension(33), parameter :: modeL1s = (/4, 1,  3, 3, 1,  7,  2, 4,  8, 5, 5,   3, 7,  1,  5,  3, 3,  2,  6,  3, 1,  5,  5,  4,  1,  4,  8,2,  2, 2,  1,  3,  4/)

    ! integer, dimension(1), parameter :: modeN2s = (/17/)
    ! integer, dimension(1), parameter :: modeL2s = (/4/)


    ! Load from file
do i_mode = 1,  nmodes
    n1      = modeN1s(i_mode)
    t1      = 'S'
    l1      = modeL1s(i_mode)
    tl1     = l1*2 + 1


    if(cross_coupled)then
        n2      = modeN2s(i_mode)
        t2      = 'S'
        l2      = modeL2s(i_mode)
        tl2     = l2*2 + 1

    else
        n2      = n1
        t2      = t1
        l2      = l1
        tl2     = tl1
    endif 

    call buffer_int(n1str, n1)
    call buffer_int(l1str, l1)
    call buffer_int(n2str, n2)
    call buffer_int(l2str, l2)

    if(vti_model)then 
        model_ti = '_VTI'
    else 
        model_ti = ''
    endif 
    
    if(cross_coupled)then 
        modefile = './output/vani'//trim(n1str)//trim(t1)//trim(l1str)//'_'//trim(n2str)//trim(t2)//trim(l2str)//trim(model_ti)//'.txt'
    else
        modefile = './output/N-0.01/sem_fast_'//trim(n1str)//trim(t1)//trim(l1str)//trim(model_ti)//'.txt'
    endif 
    call load_vani_from_file(l1, l2, modefile)
    
    

    ! Convert to real matrix 
    !allocate(Vani_real(tl1, tl2))
    !call convert_imag_to_real(l1, l1, Vani, Vani_real)
    !call get_Ssum_bounds(l1, l1, smin, smax, num_s, ncols)
    !allocate(sig_st(num_s, ncols))
    !call Hreal_to_cst(Vani_real, l1, l1, sig_st, ncols, num_s, t1, t1, 1)
    !out_name = 'output/sig_st_'//trim(nstr)//trim(t1)//trim(lstr)//trim(model_ti)//'.txt'
    !call write_cst_to_file(out_name, sig_st, ncols, num_s, smin, 2)

    call get_Ssum_bounds(l1, l2, smin, smax, num_s, ncols)
    allocate(cst_imag(num_s, ncols))
    call Hcomplex_to_cst_8(Vani, l1, l2, cst_imag, ncols, num_s, t1, t2, 1)

    if(cross_coupled)then 
        out_name =   'output/cst_'//trim(n1str)//trim(t1)//trim(l1str)//'_'//trim(n2str)//trim(t2)//trim(l2str)//trim(model_ti)//'.txt'
    else
        out_name = 'output/N-0.01/cst_'//trim(n1str)//trim(t1)//trim(l1str)//trim(model_ti)//'.txt'
    endif 

    call write_cst_complex_to_file(out_name, cst_imag, ncols, num_s, smin, 2)



    ! Output to Ensight: 
    if(output_to_ensight)then 
        write(*,*)'output_to_ensight', output_to_ensight
        call mineos%process_mineos_model(.true.)

        do iproc = 0, nprocs -1 
            sm = create_SetMesh(iproc, 3)
            ! Read the mesh info and coordinates
            call sm%read_proc_coordinates()
            call sm%load_ibool()
            ! needed for ensight geo file
            call sm%setup_global_coordinate_arrays(.false.)
            call sm%compute_rtp_from_xyz(.false.)

            
            allocate(sigma(sm%nglob))
            allocate(ylm_global(sm%nglob))

            ! Loop through the s 
            do is = 3, num_s, 2
                s = smin+is-1
                ! Loop through the t values: 
                do it = 1, 2*s +1
                    t = it - s - 1

                    ! Get the ylm as a global array: 
                    do ispec = 1, sm%nspec
                        do i = 1, sm%ngllx 
                            do j = 1, sm%nglly
                                do k = 1, sm%ngllz 
                                    ylm_global(sm%ibool(i,j,k,ispec)) = ylm_real(s, t, sm%thetastore(i,j,k,ispec),&
                                                                                    sm%phistore(i,j,k,ispec))
                                enddo 
                            enddo 
                        enddo 
                    enddo 
                    write(*,*) s, t
                    sigma = sigma + ylm_global * real(cst_imag(is,it))

                enddo 
            enddo 

    

            call create_ensight_file_prefix(iproc, 3)
            call create_proc_case_file()

            call create_proc_geo_file(sm, 1)

            call write_real_scalar_to_ensight(sm, sigma, 'sigma', 1)

            call sm%cleanup()

            deallocate(sigma)
            deallocate(ylm_global)


            write(*,*)'Finished processor ', iproc
        enddo 

    endif ! output_to_ensight



    ! Write to a surface grid for plotting in MPL 
    ! grid spacing in degrees: 
    if(output_to_evengrid)then
        dlat = 0.25 
        dlon = 0.25

        nlat = int(180.0d0/dlat) 
        nlon = int(360.0d0/dlon) 

        out_name = 'output/mpl_cst_'//trim(n1str)//trim(t1)//trim(l1str)//trim(model_ti)//'.txt'
        open(1,file=trim(out_name), form='formatted')
        write(*,*)'writing to '//trim(out_name)

        write(1,'(E15.6)', advance='yes')dlat 
        write(1,'(E15.6)', advance='yes')dlon 


        ! for each latitude we will write one row 
        ! output matrix is in format rows = lat, cols = lon 
        do ilat = 1, nlat 
            do ilon = 1, nlon 
                !write(*,*)ilat, ilon 
                ! Compute colatitude and longitude in radians
                theta = -90.0d0 + ilat*dlat ! latitude 
                !write(*,*)'theta: ', theta
                theta = 90.0d0 - theta 
                !write(*,*)'       ', theta
                theta = PI * theta/180.d0 
                !write(*,*)'       ', theta
                phi   = -180.0d0 + ilon*dlon
                !write(*,*)'phi: ', phi
                if (phi.lt.zero) phi = phi + 360.d0
                !write(*,*)'     ', phi
                phi = TWO_PI * phi/360.d0 
                !write(*,*)'     ', phi


            
                ! Loop through the s 
                sum = zero 
                do is = 3, num_s, 2
                    s = smin+is-1
                    ! Loop through the t values: 
                    do it = 1, 2*s +1
                        t = it - s - 1
                        sum = sum + ylm_real(s, t, theta, phi) * real(cst_imag(is,it))
                    enddo 
                enddo 

                if (ilon.eq.nlon)then 
                    write(1,'(E15.6)', advance='yes')sum 
                else 
                    write(1,'(E15.6)', advance='no')sum 
                endif 
            enddo 
        enddo 
    endif 

    deallocate(Vani)
    deallocate(cst_imag)
    !deallocate(sig_st)


enddo ! Imode 

end program