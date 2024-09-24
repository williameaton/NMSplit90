program integrate_mesh_1
    ! Code tests integration over the IC sphere for a function 
    ! f = 1 and f = r


    use params, only: nspec, ngllx, nglly, ngllz, rstore, nprocs
    use allocation_module, only: allocate_if_unallocated
    use mesh_utils, only: read_proc_coordinates, load_ibool, & 
                          compute_jacobian, cleanup_for_mode, &
                          compute_rtp_from_xyz
    use gll, only: setup_gll
    use integrate, only: integrate_over_mesh
    use ylm_plm, only: ylm_complex
    implicit none 
    include "constants.h"
    
    integer :: iproc, i, j, k, ispec                  
    integer :: region                ! Region code
    real(kind=CUSTOM_REAL) :: totalint_1, totalint_r
    real(kind=CUSTOM_REAL) :: integral_1, integral_r
    real(kind=CUSTOM_REAL), allocatable    :: integrand(:,:,:,:)
    character(len=80)::outfname
        
    ! Setup parameters: 
    region = 3      ! Inner core
    
    ! Read mineos model 
    call process_mineos_model()

    ! Loop over two cases: m=2 and m=3 -- should be 0 for m=3 and 
    ! non zero for m=2
    ! Initialise the integral
    totalint_1 = zero 
    totalint_r = zero 

    do iproc = 0, nprocs-1
        call read_proc_coordinates(iproc, region)
        call load_ibool(iproc, region)
        call setup_gll()
        call compute_jacobian(iproc, .false.)
        call get_mesh_radii(iproc, .false.)
        call compute_rtp_from_xyz(iproc, .false.)
        call allocate_if_unallocated(ngllx, nglly, ngllz, nspec, integrand)
        
        ! Just integration 1 over the whole IC 
        integrand(:,:,:,:) = ONE
        integral_1 = integrate_over_mesh(integrand)
        totalint_1 = totalint_1 + integral_1

        ! Set integral to = r for the IC 
        do ispec = 1, nspec 
            do i = 1, ngllx
                do j = 1, nglly
                    do k = 1, ngllz
                        integrand(i,j,k,ispec) =  real(rstore(i,j,k,ispec),CUSTOM_REAL)
                    enddo 
                enddo 
            enddo 
        enddo 
        integral_r = integrate_over_mesh(integrand)
        totalint_r = totalint_r + integral_r


        call cleanup_for_mode()
    enddo 

    open(unit=1,file='integration/test_1_int.txt', &
        status='unknown',form='formatted',action='write')
    ! Code value: 
    write(1,*)totalint_1
    ! Analytical value
    write(1,*)(FOUR/THREE)*PI*((1221.5d0/6371.d0)**THREE)


    open(unit=2,file='integration/test_r_int.txt', &
    status='unknown',form='formatted',action='write')
    ! Code value: 
    write(2,*)totalint_r
    ! Analytical value
    write(2,*) PI*((1221.5d0/6371.d0)**FOUR)


    close(1)
    close(2)    
end program integrate_mesh_1
    
    
    
    