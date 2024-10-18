program integrate_mesh_1
    ! Code tests integration over the IC sphere for a function 
    ! f = 1 and f = r
    use params, only: nprocs
    use allocation_module, only: allocate_if_unallocated
    use integrate, only: integrate_over_mesh
    use ylm_plm, only: ylm_complex
    use mineos_model, only: mineos
    use specfem_mesh, only: SetMesh, create_SetMesh
    implicit none 
    include "constants.h"
    
    integer :: iproc, i, j, k, ispec                  
    real(kind=CUSTOM_REAL) :: totalint_1, totalint_r
    real(kind=CUSTOM_REAL) :: integral_1, integral_r
    real(kind=CUSTOM_REAL), allocatable    :: integrand(:,:,:,:)
    character(len=80)::outfname
        
    type(SetMesh) :: sm 


    ! Initialise the integral
    totalint_1 = zero 
    totalint_r = zero 

    do iproc = 0, nprocs-1

        sm = create_SetMesh(iproc, 3) ! region is 3 for inner core
        call sm%read_proc_coordinates()
        call sm%load_ibool()
        call sm%setup_gll()
        call sm%compute_jacobian(.false.)
        call sm%compute_wglljac(.true.)

        call allocate_if_unallocated(sm%ngllx, sm%nglly, sm%ngllz, sm%nspec, integrand)


        ! Just integration 1 over the whole IC 
        integral_1 = zero 
        do ispec = 1, sm%nspec 
            do i = 1, sm%ngllx 
                do j = 1, sm%nglly 
                    do k = 1, sm%ngllz 
                        integral_1 = integral_1 + 1.0d0 * sm%wglljac(i,j,k,ispec)
                    enddo 
                enddo 
            enddo 
        enddo

        totalint_1 = totalint_1 + integral_1


        call sm%compute_rtp_from_xyz(.false.)
        ! Set integral to = r for the IC 
        do ispec = 1, sm%nspec 
            do i = 1, sm%ngllx
                do j = 1, sm%nglly
                    do k = 1, sm%ngllz
                        integrand(i,j,k,ispec) =  sm%rstore(i,j,k,ispec)
                    enddo 
                enddo 
            enddo 
        enddo 
        integral_r = integrate_over_mesh(sm, integrand)
        totalint_r = totalint_r + integral_r

        call sm%cleanup()
    enddo 

    open(unit=1,file='integration/test_1_int.txt', &
        status='unknown',form='formatted',action='write')
    ! Code value: 
    write(*,*)totalint_1
    write(1,*)totalint_1
    ! Analytical value
    write(*,*)(FOUR/THREE)*PI*((1221.5d0/6371.d0)**THREE)
    write(1,*)(FOUR/THREE)*PI*((1221.5d0/6371.d0)**THREE)


    open(unit=2,file='integration/test_r_int.txt', &
    status='unknown',form='formatted',action='write')
    ! Code value: 
    write(*,*)totalint_r
    write(2,*)totalint_r
    ! Analytical value
    write(*,*) PI*((1221.5d0/6371.d0)**FOUR)
    write(2,*) PI*((1221.5d0/6371.d0)**FOUR)


    close(1)
    close(2)    
end program integrate_mesh_1
    
    
    
    