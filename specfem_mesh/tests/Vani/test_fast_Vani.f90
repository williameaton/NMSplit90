program test_constant_Vani_matrix
    use params, only: Vani, nprocs
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use v_ani, only: compute_Vani_matrix_stored, save_Vani_matrix, compute_Cxyz_at_gll_constantACLNF
    use mineos_model, only: mineos, mineos_ptr
    use specfem_mesh, only: SetMesh, create_SetMesh
    use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
    implicit none
    include "constants.h"

    real(kind=CUSTOM_REAL) :: A, C, L, N, F
    integer                :: iproc, i, j, k, ispec, l1, l2, n1, n2, region, h, b
    character              ::  t1, t2
    character(len=250)     :: out_name
    real(kind=CUSTOM_REAL), allocatable :: eta1(:), eta2(:)

    type(SetMesh)          :: sm 
    type(InterpPiecewise)  :: interp 

    ! Read mineos model 
    call mineos%process_mineos_model(.true.)
    mineos_ptr => mineos

    ! Choose modes: 
    n1      = 6
    t1      = 'S'
    l1      = 10

    n2      = 6
    t2      = 'S'
    l2      = 10

    A =  0.4d0
    C = -0.2d0
    L =  0.3d0
    N = -0.5d0
    F =  0.1d0

    region = 3

    ! The matrix should be 2l + 1 from -m to m 
    call allocate_if_unallocated(2*l1+1, 2*l2+1, Vani)
    Vani = SPLINE_iZERO

    do iproc = 0, nprocs-1
        sm = create_SetMesh(iproc, region)
        write(*,*)"Processor: ", iproc
        call sm%read_proc_coordinates()
        call sm%load_ibool()
        call sm%setup_gll()
        call sm%compute_jacobian(.true.)
        call sm%compute_rtp_from_xyz(.true.)
        call sm%setup_global_coordinate_arrays(.true.)
        call sm%get_unique_radii(.true.)
        call sm%compute_rotation_matrix()
        call sm%compute_wglljac(.true.)

        allocate(eta1(sm%nglob), eta2(sm%nglob))
        eta1 = zero 
        eta2 = zero 

        call compute_Cxyz_at_gll_constantACLNF(sm, A, C, L, N, F, eta1, eta2)
        call compute_Vani_matrix_stored(sm, t1, l1, n1, t2, l2, n2)

        call sm%cleanup()
        deallocate(eta1, eta2)
    enddo 

    write(out_name, '(a,i1,a,i1,a)')'./v_ani_matrix/stored_sem_', n1, t1, l1, '.txt'
    call save_Vani_matrix(l1, out_name)


    
end program test_constant_Vani_matrix