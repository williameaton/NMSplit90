
program test_rtp_xyz
    use params, only: nglob, x_glob, y_glob, z_glob, Rmat
    use math, only: sinp, cosp, atan2p, acosp, sqrtp
    use mesh_utils, only: compute_rotation_matrix
    implicit none
    include "constants.h"

    reaL(kind=CUSTOM_REAL) :: rt2, i_rt2, vec(3)

    rt2   = two**half
    i_rt2 = one/rt2

    nglob = 3

    allocate(x_glob(nglob))
    allocate(y_glob(nglob))
    allocate(z_glob(nglob))
    x_glob = zero
    y_glob = zero
    z_glob = zero

    ! Set up some xyz values: 

    ! Should be r: (1,  0,  0)
    !           θ: (0,  0, -1)
    !           φ: (0,  1,  0)
    x_glob(1) = one
    y_glob(1) = zero
    z_glob(1) = zero 


    x_glob(2) = zero
    y_glob(2) = one
    z_glob(2) = zero 


    call compute_rotation_matrix()

    ! Check rotation matrix effect on vector in RTP (1, 1, 0): 
    ! in XYZ this is (1, 0, -1)
    vec(1) = one
    vec(2) = one
    vec(3) = zero
    vec(:) = matmul(Rmat(:,:,1), vec)

    write(*,*)vec(1), one
    write(*,*)vec(2), zero
    write(*,*)vec(3), -one


    ! Check rotation matrix effect on vector in RTP (-1, -1, 0): 
    ! in XYZ this is (0, -1, 1)
    vec(1) = -one
    vec(2) = -one
    vec(3) = zero
    vec(:) = matmul(Rmat(:,:,2), vec)

    write(*,*)vec(1), zero
    write(*,*)vec(2), -one
    write(*,*)vec(3), one









end program test_rtp_xyz