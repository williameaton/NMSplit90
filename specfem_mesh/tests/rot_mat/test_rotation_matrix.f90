
program test_W_matrix
    use params, only:  Wmat, rho_spl, nprocs
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use mineos_model, only: mineos, mineos_ptr
    use modes, only: Mode, get_mode
    use specfem_mesh, only: create_SetMesh, SetMesh
    use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
    implicit none
    include "constants.h"

    integer :: iproc, region, knot_lower, knot_upper, n1, n2, l1, l2
    character :: t1, t2
    type(SetMesh)         :: sm
    type(InterpPiecewise) :: interp
    character(len=30)     :: out_name


    ! Read mineos model 
    call mineos%process_mineos_model(.false.)
    mineos_ptr => mineos

 
    n1 = 0
    n2 = 0

    t1 = 'S'
    t2 = 'S'

    l1 = 2
    l2 = 2

    region = 3

    ! The matrix should be 2l + 1 from -m to m 
    allocate(Wmat(2*l1+1, 2*l2+1))
    Wmat = SPLINE_iZERO


    do iproc = 0, nprocs-1
        ! Things that need to be done for each processor
        sm = create_SetMesh(iproc, region)
        call sm%read_proc_coordinates()
        call sm%load_ibool()
        call sm%setup_gll()
        call sm%compute_jacobian(.true.)
        call sm%compute_rtp_from_xyz(.true.)

        call sm%setup_global_coordinate_arrays(.true.)
        call sm%compute_rtp_from_xyz(.true.)
        call sm%get_unique_radii(.true.)
        call sm%compute_rotation_matrix()

        call compute_W_matrix(sm, sm%interp, n1, t1, l1, & 
                                             n2, t2, l2, & 
                                             .true.)
        call sm%cleanup()
    enddo 


    ! Constants so multiply after
    Wmat = Wmat * OMEGA * iONE

    write(out_name, '(a,i1,a,i1,a,i1,a,i1,a)')'rot_mat/Wmat_', n1, t1, l1, '_', n2, t2, l2, '.txt'
    call save_W_matrix(l1, l2, trim(out_name))



end program test_W_matrix