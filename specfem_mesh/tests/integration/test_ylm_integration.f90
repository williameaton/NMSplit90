program test_ylm_integration
    use params, only: nspec, ngllx, nglly, ngllz, rstore, thetastore, phistore
    use allocation_module, only: allocate_if_unallocated
    use mesh_utils, only: read_proc_coordinates, load_ibool, & 
                          compute_jacobian, cleanup_for_mode, &
                          compute_rtp_from_xyz
    use gll, only: setup_gll
    use integrate, only: integrate_over_mesh, integrate_complex_mesh_scalar
    use ylm_plm, only: ylm_complex
    implicit none 
    include "constants.h"
    
    integer :: iproc, i, j, k, ispec                  
    integer :: nprocs                ! Number of processors used by mesher 
    integer :: region                ! Region code
    real(kind=CUSTOM_REAL)    :: rreal, phreal, threal
    complex(kind=CUSTOM_REAL) :: totalint, ylm1, ylm2
    complex(kind=CUSTOM_REAL) :: integral
    complex(kind=CUSTOM_REAL), allocatable    :: integrand(:,:,:,:)
    integer :: l1, m1, l2, m2
    character(len=80)::outfname
        
    ! Setup parameters: 
    region = 3      ! Inner core
    nprocs = 6
    
    l1 = 5
    m1 = 2
    
    l2 = 5

    ! Read mineos model 
    call process_mineos_model()

    ! Loop over two cases: m=2 and m=3 -- should be 0 for m=3 and 
    ! non zero for m=2
    do m2 = 2,3
        ! Initialise the integral
        totalint = zero 

        do iproc = 0, nprocs-1
            call read_proc_coordinates(iproc, region)
            call load_ibool(iproc, region)
            call setup_gll()
            call compute_jacobian()
            call get_mesh_radii()
            call compute_rtp_from_xyz()
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
            totalint = totalint + integral
            call cleanup_for_mode()
        enddo 

        write(outfname, '(a,i1,a)')'test_ylm_int_', m2, '.txt'
        open(unit=1,file=trim(outfname), &
            status='unknown',form='formatted',action='write')

        ! Code value: 
        write(1,*)real(totalint), aimag(totalint)

        ! Analytical value
        if (l1.eq.l2 .and. m1.eq.m2)then

            ! Because of orthogonality B.35 the integral over the volume 
            ! collapses to 1/3 b^3 where b is the IC radius
            ! if l=l' and m=m'
            write(1,*)((1221500.0d0/6371000.d0)**THREE)/THREE
        else
            write(1,*)0.0d0
        endif
    enddo 
    
end program test_ylm_integration
    
    
    
    
