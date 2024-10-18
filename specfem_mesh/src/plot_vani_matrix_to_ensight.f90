program plot_vani_to_ensight
! program loads in a Vani matrix (complex) and outputs the splitting functions
    use params, only: Vani, nprocs, nmodes
    use v_ani, only: load_vani_from_file, convert_imag_to_real, save_Vani_real_matrix
    use splitting_function, only: get_Ssum_bounds, cst_to_H, H_to_cst, write_cst_to_file
    use specfem_mesh,       only: SetMesh, create_SetMesh
    use modes,              only: get_mode, Mode 
    use mineos_model,       only: mineos, mineos_ptr
    use ylm_plm,            only: ylm_real    

    implicit none 
    include "constants.h"
    character(len=3) :: nstr, lstr
    integer   :: l1, n1 , tl1, smin, smax, num_s, ncols, nrows, is, it, & 
                 iproc, i, j, k, ispec, s, t, ilat, ilon, nlat, nlon, i_mode
    character :: t1
    character(len=20) :: model_ti
    character(len=250) :: modefile , out_name
    real(kind=CUSTOM_REAL)  :: dlat, dlon, theta, phi
    real(kind=SPLINE_REAL)  :: sum
    real(kind=SPLINE_REAL), allocatable :: Vani_real(:,:)
    real(kind=SPLINE_REAL), allocatable :: cst(:,:)
    real(kind=CUSTOM_REAL), allocatable :: sigma(:), ylm_global(:)
    type(SetMesh) :: sm
    logical, parameter :: vti_model = .false.
    logical, parameter :: output_to_ensight = .false.

    ! Modes: 
    integer, dimension(20), parameter :: modeNs = (/2, 3, 9, 9, 9, 11, 11, 13,13,13,13,15,15,18,18,20,21,25,27, 6/)
    integer, dimension(20), parameter :: modeLs = (/3, 2, 2, 3, 4,  4,  5,  1, 2, 3, 6, 3, 4, 3, 4, 1, 6,2,2, 10/)
    


    ! Load from file
do i_mode = 1, 1! nmodes
    n1      = 13!modeNs(i_mode)
    t1      = 'S'
    l1      = 6!modeLs(i_mode)
    tl1     = l1*2 + 1

    call buffer_int(nstr, n1)
    call buffer_int(lstr, l1)
    
    if(vti_model)then 
        model_ti = '_VTI'
    else 
        model_ti = ''
    endif 
    modefile = './output/sem_fast_'//trim(nstr)//trim(t1)//trim(lstr)//trim(model_ti)//'.txt'
    call load_vani_from_file(l1, modefile)
    

    ! Convert to real matrix 
    allocate(Vani_real(tl1, tl1))
    call convert_imag_to_real(l1, l1, Vani, Vani_real)

    out_name =  './output/real_'//trim(nstr)//trim(t1)//trim(lstr)//trim(model_ti)//'.txt'
    call save_Vani_real_matrix(l1, Vani_real, out_name)

    call get_Ssum_bounds(l1, l1, smin, smax, num_s, ncols)
    allocate(cst(num_s, ncols))

    call H_to_cst(Vani_real, l1, l1, cst, ncols, num_s, t1, t1)

    out_name = 'output/cst_'//trim(nstr)//trim(t1)//trim(lstr)//trim(model_ti)//'.txt'
    call write_cst_to_file(out_name, cst, ncols, num_s, smin, 2)



    write(*,*)'output_to_ensight', output_to_ensight
    ! Output to Ensight: 
    if(output_to_ensight)then 
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
                    sigma = sigma + ylm_global * real(cst(is,it))

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
    dlat = 0.25 
    dlon = 0.25

    nlat = int(180.0d0/dlat) 
    nlon = int(360.0d0/dlon) 

    out_name = 'output/mpl_cst_'//trim(nstr)//trim(t1)//trim(lstr)//trim(model_ti)//'.txt'
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
                    sum = sum + ylm_real(s, t, theta, phi) * real(cst(is,it))
                enddo 
            enddo 

            if (ilon.eq.nlon)then 
                write(1,'(E15.6)', advance='yes')sum 
            else 
                write(1,'(E15.6)', advance='no')sum 
            endif 
        enddo 
    enddo 


    deallocate(Vani_real)
    deallocate(cst)


enddo ! Imode 

end program