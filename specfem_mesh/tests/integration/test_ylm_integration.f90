program test_ylm_integration
    use params, only: nprocs
    use allocation_module, only: allocate_if_unallocated
    use integrate, only: integrate_over_mesh, integrate_complex_mesh_scalar
    use ylm_plm, only: ylm_complex
    use mineos_model, only: mineos
    use specfem_mesh, only: SetMesh, create_SetMesh
    implicit none 
    include "constants.h"
    
    integer :: iproc, i, j, k, ispec                  
    real(kind=CUSTOM_REAL)    :: rreal, phreal, threal
    complex(kind=CUSTOM_REAL) :: totalint, ylm1, ylm2
    complex(kind=CUSTOM_REAL) :: integral
    complex(kind=CUSTOM_REAL), allocatable    :: integrand(:,:,:,:)
    integer :: l1, m1, l2, m2
    character(len=80)::outfname
        
    type(SetMesh) :: sm 


    ! Setup parameters:     
    l1 = 5
    m1 = 2
    
    l2 = 5

    ! Read mineos model 
    call mineos%process_mineos_model(.false.)

    ! Loop over two cases: m=2 and m=3 -- should be 0 for m=3 and 
    ! non zero for m=2
    do m2 = 2,3
        ! Initialise the integral
        totalint = zero 

        do iproc = 0, nprocs-1
            sm = create_SetMesh(iproc, 3)
            call sm%read_proc_coordinates()
            call sm%load_ibool()
            call sm%setup_gll()
            call sm%compute_jacobian(.false.)
            call sm%get_unique_radii(.false.)
            call sm%compute_rtp_from_xyz(.false.)

            call allocate_if_unallocated(sm%ngllx, sm%nglly, sm%ngllz, sm%nspec, integrand)
            
            do ispec = 1, sm%nspec 
                do i = 1, sm%ngllx
                    do j = 1, sm%nglly
                        do k = 1, sm%ngllz
                            rreal  = real(sm%rstore(i,j,k,ispec),     CUSTOM_REAL)
                            threal = real(sm%thetastore(i,j,k,ispec), CUSTOM_REAL)
                            phreal = real(sm%phistore(i,j,k,ispec),   CUSTOM_REAL)
                            ylm1   = conjg(ylm_complex(l1, m1, threal, phreal))
                            ylm2   = ylm_complex(l2, m2, threal, phreal)
                            integrand(i,j,k,ispec) =  ylm1 * ylm2
                        enddo 
                    enddo 
                enddo 
            enddo 

            integral = integrate_complex_mesh_scalar(sm, integrand)
            totalint = totalint + integral
            call sm%cleanup()
        enddo 

        write(outfname, '(a,i1,a)')'integration/test_ylm_int_', m2, '.txt'
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
    
    
    
    
