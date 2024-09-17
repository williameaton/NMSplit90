
program Tromp93_SEM
    use params, only: Vani, datadir, n_unique_rad, rad_id, rho_spl, vp_spl, A0, & 
                      Arad, Crad, Frad, Lrad, Nrad, IC_ID, vp, unique_r, interp_id_r,& 
                       n_unique_rad, rho_mineos
    use spline, only: interpolate_mineos_variable
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use mesh_utils, only: cleanup_for_mode, compute_jacobian, compute_rotation_matrix, & 
                          compute_rtp_from_xyz, load_ibool, read_proc_coordinates, & 
                          setup_global_coordinate_arrays
    use gll, only: setup_gll
    use v_ani, only: compute_Vani_matrix, save_Vani_matrix, compute_Cxyz_at_gll_radialACLNF

    implicit none
    include "constants.h"

    real(kind=CUSTOM_REAL) :: A, C, L, N, F
    integer :: iproc, i,j,k,ispec, l1, l2, n1, n2, nproc, region, tl1, tl2, h, b
    character ::  type_1, type_2
    character(len=250) :: out_name, aclnf_dir
    real(kind=SPLINE_REAL) :: min_r, min_i, thirty, twone

    aclnf_dir = '/Users/eaton/Documents/Software/NMSplit90/& 
                specfem_mesh/tests/v_ani_matrix/Tromp_1993_model/ACLNF'


    ! Read mineos model 
    call process_mineos_model()

    ! Choose modes: 
    n1      = 9
    type_1 = 'S'
    l1      = 3

    region = 3
    nproc  = 6

    ! Setup Vani matrix
    tl1 = 2*l1 + 1

    ! The matrix should be 2l + 1 from -m to m 
    call allocate_if_unallocated(tl1, tl1, Vani)
    Vani = SPLINE_iZERO

    do iproc = 0, nproc-1
        write(*,*)"Processor: ", iproc
        ! Things that need to be done for each processor
        call read_proc_coordinates(iproc, region)

        call load_ibool(iproc, region)
        call setup_gll()
        call compute_jacobian()

        call setup_global_coordinate_arrays()
        call compute_rtp_from_xyz()
        call get_mesh_radii()
        call compute_rotation_matrix()

        call load_ACLNF_from_files(aclnf_dir, iproc, n_unique_rad)


        ! Interpolate the vp and rho to the relevant points 
        call deallocate_if_allocated(rho_spl)
        allocate(rho_spl(n_unique_rad))
        call interpolate_mineos_variable(real(rho_mineos, kind=SPLINE_REAL), 1, IC_ID, & 
                                        unique_r, n_unique_rad, rho_spl, interp_id_r)
        call deallocate_if_allocated(vp_spl)
        allocate(vp_spl(n_unique_rad))
        call interpolate_mineos_variable(real(vp, kind=SPLINE_REAL), 1, IC_ID, & 
                                         unique_r, n_unique_rad, vp_spl, interp_id_r)
        call deallocate_if_allocated(A0)
        allocate(A0(n_unique_rad))

        ! kappa + 4*mu 
        A0 = (vp_spl*vp_spl)*rho_spl

        Arad = Arad * A0 
        Crad = Crad * A0 
        Lrad = Lrad * A0 
        Nrad = Nrad * A0 
        Frad = Frad * A0 





        call compute_Cxyz_at_gll_radialACLNF(zero, zero, rad_id)

        call compute_Vani_matrix(type_1, l1, n1, type_1, l1, n1, .false., iproc)

        call cleanup_for_mode()
    enddo 

    write(out_name, '(a,i1,a,i1,a)')'./matrices/sem_', n1, type_1, l1, '.txt'
    call save_Vani_matrix(l1, out_name)

end program Tromp93_SEM