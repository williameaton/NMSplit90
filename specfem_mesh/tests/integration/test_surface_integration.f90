program test_surface_integration
        use params,             only: nprocs, nmodes
        use specfem_mesh,       only: SetMesh, create_SetMesh
        use allocation_module,  only: allocate_if_unallocated, deallocate_if_allocated
        use modes,              only: get_mode, Mode 
        use mineos_model,       only: mineos, mineos_ptr
        use ylm_plm,            only: ylm_real    
        use model3d, only: M3D

        implicit none 
        include "constants.h"
        character(len=3) :: nstr, lstr
        integer   :: l1, n1, tl1, ispec, i, j, iproc
        character :: t1
        character(len=250) :: out_name
        real(kind=CUSTOM_REAL)  :: integral

        type(SetMesh) :: sm
        type(Mode)    :: mode_1
        type(M3D)     :: model3D

        integer, parameter :: region = 1


        call mineos%process_mineos_model(.false.)
        mineos_ptr => mineos


        integral = zero 

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
 
            call sm%compute_surface_jacobian()
            call sm%get_unique_radii(.false.)

            do ispec = 1, sm%nspec2D_top
                do i = 1, sm%ngllx
                    do j = 1, sm%nglly
                        integral = integral + sm%jacobian2D_top(i,j,ispec)*sm%wgll(i)*sm%wgll(j)
                    enddo 
                enddo 
            enddo 

            call sm%cleanup()
        enddo 

        if(abs(integral - four*PI) .gt. 1e-5)then 
            write(*,*)'ERROR in surface integral'
            stop 
        endif 
        
            
    
end program test_surface_integration