program test_ellipticity_gll
    use params,             only: nprocs, nmodes, Tmat,Tell
    use specfem_mesh,       only: SetMesh, create_SetMesh
    use allocation_module,  only: allocate_if_unallocated, deallocate_if_allocated
    use modes,              only: get_mode, Mode 
    use mineos_model,       only: mineos, mineos_ptr
    use ylm_plm,            only: ylm_real    
    use model3d, only: M3D
    use ylm_plm, only: legendre
    use PREMModel, only: get_radial_density_derivative_at_radius
    implicit none 
    include "constants.h"
    character(len=3) :: nstr, lstr, procstr
    integer   :: l1, n1, n2, l2, tl1, ispec, i, j, iproc, k, ii, jj
    character :: t1, t2
    character(len=250) :: out_name
    real(kind=CUSTOM_REAL)  :: integral, drhodr
    real(kind=SPLINE_REAL), allocatable  :: ell_spline(:)
    real(kind=CUSTOM_REAL), allocatable  :: delta_rho(:,:,:,:)

    type(SetMesh) :: sm
    !type(Mode)    :: mode_1
    type(M3D)     :: model3D

    integer, parameter :: region = 1


    n1 = 16
    t1 = 'S'
    l1 = 5

    n2 = 16
    t2 = 'S'
    l2 = 5


    call mineos%process_mineos_model(.false.)
    mineos_ptr => mineos
    call mineos%load_ellipticity_from_file('../ellipticity/mineos_epsilon.txt')

    allocate(Tmat(2*l1+1, 2*l2+1))
    allocate(Tell(2*l1+1, 2*l2+1))


    do iproc = 0, nprocs -1 
        sm = create_SetMesh(iproc, region)
        ! Read the mesh info and coordinates
        call sm%read_proc_coordinates()
        call sm%load_ibool()
        call sm%setup_gll()
        call sm%load_original_boundaries()

        ! needed for ensight geo file
        call sm%setup_global_coordinate_arrays(.false.)
        call sm%compute_rtp_from_xyz(.false.)

        call sm%get_unique_radii(.false.)
        call sm%compute_jacobian(.false.)
        call sm%compute_surface_jacobian()

        call sm%compute_elliptical_boundary_perturbation('../ellipticity/epsi_discontinuities.txt', mineos%ndisc)


        ! Load epsilon and interpolate it to the unique radii for this set: 
        allocate(ell_spline(sm%interp%n_radial))
        call sm%interp%interpolate_mineos_variable(real(mineos%ell_mineos,kind=SPLINE_REAL), ell_spline)

        ! Fix some issues at the discontinuities
        do ii=1, sm%interp%n_radial
            do jj = 1, mineos%NR
                if(mineos%rdisc(jj)/SCALE_R .eq. sm%unique_r(ii))then 
                    ell_spline(ii) = mineos%ell_mineos(jj)
                endif 
            enddo 
        enddo 

    
        ! Compute delta rho as a function only of ellipticity: 
        allocate(delta_rho(sm%ngllx, sm%nglly, sm%ngllz, sm%nspec))


        do ispec = 1, sm%nspec
            do i = 1, sm%ngllx
                do j = 1, sm%nglly
                    do k = 1, sm%ngllz
                        call get_radial_density_derivative_at_radius(sm%rstore(i,j,k,ispec), drhodr, k)
                        delta_rho(i,j,k,ispec) = (two/three) * sm%rstore(i,j,k,ispec) * ell_spline(sm%rad_id(i,j,k,ispec)) * drhodr * legendre(2, dcos(sm%thetastore(i,j,k,ispec)))
                        
                        if(abs(ell_spline(sm%rad_id(i,j,k,ispec)).gt.100.0))then 
                            write(*,*)'Bad ell_spline(sm%rad_id(i,j,k,ispec)) value', ell_spline(sm%rad_id(i,j,k,ispec))
                            stop
                        endif
                        if(abs(drhodr.gt.100.0))then 
                            write(*,*)'Bad drhodr value', drhodr
                            stop
                        endif
                    
                        if(abs(legendre(2, dcos(sm%thetastore(i,j,k,ispec))).gt.100.0))then 
                            write(*,*)'Bad drhodr value', legendre(2, dcos(sm%thetastore(i,j,k,ispec)))
                            stop
                        endif
                    
                    
                    enddo 
                enddo 
            enddo 
        enddo

        write(*,*)minval(delta_rho), maxval(delta_rho)

        ! Compute the volumetric component of Tmat:
        ! First boolean says to compute the ellipticity integral over the surface 
        call compute_T_matrix(sm, sm%interp, delta_rho, n1, t1, l1, n2, t2, l2, .true., .false.)
    

        call buffer_int(nstr, n1)
        call buffer_int(lstr, l1)
        call buffer_int(procstr, iproc)
        write(out_name, '(a)')'./ellipticity/Tell_'//trim(nstr)//t1//trim(lstr)//'_'//trim(nstr)//t2//trim(lstr)//'proc'//trim(procstr)//'.txt'
        Tell = Tmat
        call save_Tell_matrix(l1, l1, out_name)
    

        deallocate(ell_spline)
        deallocate(delta_rho)
        deallocate(sm%jacobian2D_top)
        deallocate(sm%normal_top)
        deallocate(sm%jacobian2D_bottom)
        deallocate(sm%normal_bottom)

        deallocate(sm%delta_surf_top)
        deallocate(sm%delta_surf_bottom)
        call sm%cleanup()
   
    enddo 


    write(out_name, '(a)')'./ellipticity/Tell_'//trim(nstr)//t1//trim(lstr)//'_'//trim(nstr)//t2//trim(lstr)//'done'//'.txt'
    Tell = Tmat
    call save_Tell_matrix(l1, l1, out_name)

        

end program test_ellipticity_gll