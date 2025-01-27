program integrate_mesh
    use params, only: nspec, ngllx, nglly, ngllz, rstore, thetastore, phistore, nprocs
    use allocation_module, only: allocate_if_unallocated
    use mesh_utils, only: read_proc_coordinates, load_ibool, & 
                          compute_jacobian, cleanup_for_mode, &
                          compute_rtp_from_xyz
    use gll, only: setup_gll
    use ylm_plm, only: ylm_complex
    implicit none 
    include "constants.h"
    
    integer :: iproc, i, j, k, ispec                  
    integer :: region                ! Region code
    real(kind=CUSTOM_REAL)    :: rreal, phreal, threal
    complex(kind=CUSTOM_REAL) :: totalint, ylm1, ylm2
    complex(kind=CUSTOM_REAL) :: integral
    complex(kind=CUSTOM_REAL), allocatable    :: integrand(:,:,:,:)
    integer :: n, l1, m1, l2, m2
        
    ! Setup parameters: 
    region = 3      ! Inner core
    
    n = 10
    l1 = 5
    m1 = 2
    
    l2 = 5
    m2 = 2

    
    ! Read mineos model 
    call process_mineos_model(.false.)
    

    totalint = zero 

    do iproc = 0, nprocs-1
        
        ! Read the mesh info and coordinates
        call read_proc_coordinates(iproc, region)
        
        ! Load ibool variable: 
        call load_ibool(iproc, region)
        call setup_gll()
        call compute_jacobian(.false.)
        
        ! Get unique mesh radii that are present
        call get_mesh_radii(iproc, .false.)
        
        ! We will need the r theta phi coordinates: 
        call compute_rtp_from_xyz(iproc, .false.)
        
        call allocate_if_unallocated(ngllx, nglly, ngllz, nspec, integrand)
        
        do ispec = 1, nspec 
            do i = 1, ngllx
                do j = 1, nglly
                    do k = 1, ngllz
                        rreal  = real(rstore(i,j,k,ispec),     CUSTOM_REAL)
                        threal = real(thetastore(i,j,k,ispec), CUSTOM_REAL)
                        phreal = real(phistore(i,j,k,ispec),   CUSTOM_REAL)

                        ylm1   = conjg(ylm_complex(l1, m1, threal, phreal))

                        ylm2   = ylm_complex(l2, m2, threal, phreal)
                        integrand(i,j,k,ispec) =  ylm1 * ylm2
                    enddo 
                enddo 
            enddo 
        enddo 

        integral = integrate_complex_mesh_scalar(integrand)

        write(*,*)'integral ', integral

        totalint = totalint + integral
        call cleanup_for_mode()
    enddo 


    ! Analytical value: 
    write(*,*)'Total value       : ', totalint 
    if (l1.eq.l2 .and. m1.eq.m2)then
        write(*,*)'Analytical value  : ', ((1221500.0d0/6371000.d0)**THREE)/THREE
    else
        write(*,*)'Analytical value  : ', 0.0d0
    endif
    
    
end program integrate_mesh
    
    
    
    