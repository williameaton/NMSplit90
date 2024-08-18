
subroutine get_mesh_radii()
    ! Determines a list of the unique radii of the GLL point
    ! Points each GLL point to its radial value in the list 
    ! 

    use params, only: ngllx, nglly, ngllz, nspec, RA, &
                      xstore, ystore, zstore, unique_r, n_unique_rad, rad_id,&
                      interp_id_r, rad_mineos, IC_ID
    implicit none 
    include "constants.h"

    ! Local variables: 
    double precision :: rr(ngllx, nglly, ngllz, nspec) , tol, rgll
    double precision, allocatable :: unique_r_tmp(:)
    integer :: i, j, k, ispec, size_r, unique_id, ur, i_unq, i_knot, min_knot
    
    logical :: match

    ! Compute radii
    rr  = (xstore**TWO + ystore**TWO + zstore**TWO)**HALF

    ! Maximum number of unique radii is this
    size_r =  size(rr)
    allocate(unique_r_tmp(size_r))
    unique_r_tmp = ZERO

    ! Set tolerance for 'same radii' 
    tol = 1.0e-4


    ! Initial value
    unique_id = 1

    ! Allocate the array that will store the unique radius IDs: 
    allocate(rad_id(ngllx, nglly, ngllz, nspec))

    do ispec = 1, nspec
        do i = 1, ngllx
            do j = 1, nglly
                do k = 1, ngllz
                    ! Radius of this gll 
                    rgll = rr(i, j, k, ispec)
                    match = .false. 

                    ! Loop over current unique radii: 
                    do ur = 1, unique_id 
                        if (abs(unique_r_tmp(ur) - rgll) < tol)then 
                            ! already a 'unique' radius
                            match = .true.
                            rad_id(i,j,k,ispec) = ur
                        endif 
                    enddo ! ur 

                    ! No matches found
                    if (.not. match) then 
                        ! a new 'unique' radius 
                        rad_id(i,j,k,ispec)     = unique_id
                        unique_r_tmp(unique_id) = rgll
                        unique_id               = unique_id + 1
                    endif
                enddo 
            enddo 
        enddo 
    enddo 

    
    ! Allocate new unique radius array of the correct size
    ! It might be that the very last node tested in unique and so 
    ! unique_id is actually 1 too large -- test if last radius is 0 and if so cute 
    if (unique_r_tmp(unique_id).eq.ZERO)then 
        n_unique_rad = unique_id - 1
    else
        n_unique_rad = unique_id
    endif

    ! Copy over the unique values to the smaller array and deallocate long.
    allocate(unique_r(n_unique_rad)) 
    unique_r(1:n_unique_rad) = unique_r_tmp(1:n_unique_rad)
    deallocate(unique_r_tmp)

    ! Print some details
    write(*,'(/,a)')'• Determine unique mesh radii'
    write(*,'(a,es8.1)')  '  -- tolerance value       :', tol
    write(*,'(a,f8.2,a)') '  -- matches radii within  :', tol*RA, ' metres'
    write(*,'(a,i8,a,i8)')'  -- number of unique radii:', n_unique_rad, ' out of ', size_r
    
    ! Prints the difference between the original radius and the 'unique' radius it is pointed to 
    !do ispec = 1, nspec
    !    do i = 1, ngllx
    !        do j = 1, nglly
    !            do k = 1, ngllz
    !                if (abs(rr(i, j, k, ispec)- unique_r(rad_id(i,j,k,ispec))) > tol)&
    !                write(*,*) rr(i, j, k, ispec)- unique_r(rad_id(i,j,k,ispec)) 
    !            enddo 
    !        enddo 
    !    enddo 
    !enddo


    ! Find the knot id's of the unique radii 
    ! we will want to interpolate to 
    allocate(interp_id_r(n_unique_rad))
    do i_unq = 1, n_unique_rad
        ! For each radius we want to find the last knot that it is larger than
        do i_knot = 1, IC_ID+1
            if (unique_r(i_unq) .lt. rad_mineos(i_knot))then 
                ! This should be the first time that the knot is above the
                ! radius in question and so we want the i_knot - 1 to be stored
                interp_id_r(i_unq) = i_knot - 1 
                exit 
            else
            endif
        enddo 
    enddo 



end subroutine get_mesh_radii




subroutine compute_global_mode_displacement(mode_type, nord, l, m, disp)
    use params, only: rad_mineos, NL, IC_ID, n_unique_rad, nspec, ngllx, & 
    nglly, ngllz, thetastore, phistore, rad_id, rstore,& 
    unique_r, u_spl, v_spl, udot_spl, vdot_spl, xx, zz
    use spline, only: interpolate_mode_eigenfunctions
    use ylm_plm, only: ylm_complex, ylm_deriv
    implicit none
    include "constants.h"

    ! IO variables: 
    character(len=1) :: mode_type     ! mode type; S or T
    integer          :: nord          ! n value of mode
    integer          :: l             ! l value (degree) of mode
    integer          :: m             ! m value (order)  of mode
    complex(kind=CUSTOM_REAL) :: disp(3,ngllx, nglly, ngllz, nspec)
    
    ! Local variables: 
    real(kind=4)  :: omega            ! normalized anglar frequency   
    real(kind=4)  :: qval             ! normalized Q factor            
    real(kind=4)  :: u(NL), du(NL)    ! must be 4 bytes
    real(kind=4)  :: v(NL), dv(NL)    ! must be 4 bytes
    complex(kind=CUSTOM_REAL) :: ylm, dylm_theta, dylm_phi
    real(kind=CUSTOM_REAL)    ::  mf, lf, ll1, theta, phi, u_r, v_r,w_r
    integer :: i,j,k,ispec


    ! Load the mode from the database
    call get_mode(mode_type,nord,l,omega,qval,u,du,v,dv, .true.)

    ! Spline interpolation: 
    call interpolate_mode_eigenfunctions(mode_type, u, v, du, dv)
    
    ! DT98 below D.1: k = sqrt(l(l+1))
    lf  = real(l, kind=CUSTOM_REAL)
    mf  = real(m, kind=CUSTOM_REAL)
    ll1 = sqrt(lf*(lf+ONE))
    v_spl    = v_spl    / ll1
    vdot_spl = vdot_spl / ll1
    if(mode_type.eq.'T')then 
        ! Code uses u and du to store w and dw
        u_spl    = u_spl    / ll1
        udot_spl = udot_spl / ll1
    endif 


    if (mode_type.eq.'S')then 
        do ispec = 1, nspec
            do i = 1, ngllx
                do j = 1, nglly
                    do k = 1, ngllz
                    
                        ! Get ylm and the partial derivatives
                        theta = thetastore(i,j,k,ispec)
                        phi   = phistore(i,j,k,ispec)
                        ylm   = ylm_complex(l,m, theta, phi)
                        call ylm_deriv(l, m, theta, phi, dylm_theta, dylm_phi)
                        
                        ! u, v at at the radius of this node
                        u_r   =    u_spl(rad_id(i,j,k,ispec))
                        v_r   =    v_spl(rad_id(i,j,k,ispec))

                        ! S_r     (DT98 D.4)
                        disp(1,i,j,k,ispec) = ylm*u_r  

                        ! S_theta (DT98 D.5)
                        disp(2,i,j,k,ispec) = v_r*dylm_theta

                        ! S_phi   (DT98 D.6)
                        disp(3,i,j,k,ispec) = iONE * mf * v_r * ylm / sin(theta)
                    enddo 
                enddo 
            enddo 
        enddo 
    elseif (mode_type.eq.'T')then 
        do ispec = 1, nspec
            do i = 1, ngllx
                do j = 1, nglly
                    do k = 1, ngllz
                    
                        ! Get ylm and the partial derivatives, and 
                        ! w at the radius of this node
                        w_r   = u_spl(rad_id(i,j,k,ispec))
                        theta = thetastore(i,j,k,ispec)
                        phi   = phistore(i,j,k,ispec)
                        ylm   = ylm_complex(l,m, theta, phi)
                        call ylm_deriv(l, m, theta, phi, dylm_theta, dylm_phi)

                        ! S_r     (DT98 D.4)
                        disp(1,i,j,k,ispec) = (ZERO, ZERO) 

                        ! S_theta (DT98 D.5)
                        disp(2,i,j,k,ispec) = iONE * mf * w_r * ylm / sin(theta)

                        ! S_phi   (DT98 D.6)
                        disp(3,i,j,k,ispec) = -w_r * dylm_theta 
                    enddo 
                enddo 
            enddo 
        enddo 

    else
        write(*,*)'ERROR: Mode type must be S or T but was ', mode_type
        stop
    endif 

end subroutine compute_global_mode_displacement





subroutine compute_gll_mode_strain(mode_type, nord, l, m, strain)
    use params, only: rad_mineos, NL, IC_ID, n_unique_rad, nspec, ngllx, & 
                      nglly, ngllz, thetastore, phistore, rad_id, rstore,& 
                      unique_r, u_spl, v_spl, udot_spl, vdot_spl, xx, zz
    use ylm_plm, only: ylm_complex, ylm_deriv
    use spline, only: interpolate_mode_eigenfunctions
    implicit none
    include "constants.h"

    ! IO variables: 
    character(len=1) :: mode_type     ! mode type; S or T
    integer          :: nord          ! n value of mode
    integer          :: l             ! l value (degree) of mode
    integer          :: m             ! m value (order)  of mode
    complex(kind=CUSTOM_REAL) :: strain(6,ngllx, nglly, ngllz, nspec)
    ! Local variables: 
    real(kind=4)  :: omega  ! normalized anglar frequency   
    real(kind=4)  :: qval   ! normalized Q factor        

    complex(kind=CUSTOM_REAL) :: ylm, dylm_theta, dylm_phi
    real(kind=CUSTOM_REAL)    :: theta, phi, du_r, u_r, & 
                                 unq_r, dv_r, v_r, mf, lf, ll1, xx_r, & 
                                 zz_r, sinth, tanth, w_r, dw_r
    integer :: ispec, i,j,k
    ! Mineos eigenfunction values:
    ! Note that these NEED to be 4-byte because that is the mineos
    ! output format 
    real(kind=4)            :: u(NL), du(NL)
    real(kind=4)            :: v(NL), dv(NL)

    ! Load the mode from the database
    call get_mode(mode_type,nord,l,omega,qval,u,du,v,dv, .true.)

    ! Spline interpolation: 
    !   - Computes the value of u, v, du, dv at each of the unique radial values
    call interpolate_mode_eigenfunctions(mode_type, u, v, du, dv)
    
    ! DT98 below D.1: k = sqrt(l(l+1))
    lf  = real(l, kind=CUSTOM_REAL)
    mf  = real(m, kind=CUSTOM_REAL)
    ll1 = sqrt(lf*(lf+ONE))

    ! Convert eigenfunctions to auxillary form
    ! DT98 D.1: 
    !    u = U     v = V/k     w = W/k     p = P
    v_spl    = v_spl    / ll1
    vdot_spl = vdot_spl / ll1
    if(mode_type.eq.'T')then 
        ! Code uses u and du to store w and dw
        !TODO: Do we want to move this inside the if statement below? 
        u_spl    = u_spl    / ll1
        udot_spl = udot_spl / ll1
    endif 


    ! Compute x and z auxillary variables - DT98 D.20
    ! If spheroidal mode then u,v are non zero but w, wdot are 0
    !   --> z should be 0 
    ! If toroidal mode then u,v are zero but w, wdot are not
    !   --> x should be 0
    if (mode_type.eq.'S')then 
        allocate(xx(n_unique_rad))
        xx = udot_spl + (u_spl - v_spl)/unique_r

        ! Loop over each GLL point: 
        do ispec = 1, nspec
            do i = 1, ngllx
                do j = 1, nglly
                    do k = 1, ngllz
                        
                        ! Get ylm and the partial derivatives
                        theta = thetastore(i,j,k,ispec)
                        phi   = phistore(i,j,k,ispec)
                        ylm   = ylm_complex(l,m, theta, phi)
                        call ylm_deriv(l, m, theta, phi, dylm_theta, dylm_phi)
                        

                        ! u, v, udot, vdo, x, at at the radius of this node
                        u_r   =    u_spl(rad_id(i,j,k,ispec))
                        v_r   =    v_spl(rad_id(i,j,k,ispec))
                        du_r  = udot_spl(rad_id(i,j,k,ispec))
                        dv_r  = vdot_spl(rad_id(i,j,k,ispec))
                        xx_r  =       xx(rad_id(i,j,k,ispec))

                        ! r and theta at the node  
                        unq_r = rstore(i,j,k,ispec)
                        sinth = sin(theta)
                        tanth = tan(theta)

                        ! E_rr: DT98 D.14
                        strain(1,i,j,k,ispec) = ylm*du_r  
                        
                        ! E_tt: DT98 D.15
                        strain(2,i,j,k,ispec) = (ylm*u_r - v_r*(dylm_theta/tanth -  ylm*(mf/sinth)**TWO + ll1*ll1*ylm))/unq_r

                        ! E_pp: DT98 D.16
                        strain(3,i,j,k,ispec) = (ylm*u_r + v_r*(dylm_theta/tanth -  ylm*(mf/sinth)**TWO))/unq_r

                        ! E_rt: DT98 D.17
                        strain(4,i,j,k,ispec) = HALF * xx_r * dylm_theta

                        ! E_rp: DT98 D.18
                        strain(5,i,j,k,ispec) = HALF * iONE * mf * xx_r * ylm/sinth

                        ! E_tp: DT98 D.19
                        strain(6,i,j,k,ispec) = iONE * mf * v_r * (dylm_theta - ylm/tanth) / (unq_r * sinth)

                    enddo 
                enddo 
            enddo 
        enddo


    elseif (mode_type.eq.'T')then 
        ! I believe that get_mode stores w, wdot in the u and du terms
        allocate(zz(n_unique_rad))
        zz = udot_spl - u_spl/unique_r

        ! Loop over each GLL point: 
        do ispec = 1, nspec
            do i = 1, ngllx
                do j = 1, nglly
                    do k = 1, ngllz
                        
                        ! Get ylm and the partial derivatives
                        theta = thetastore(i,j,k,ispec)
                        phi   = phistore(i,j,k,ispec)
                        ylm   = ylm_complex(l,m, theta, phi)
                        call ylm_deriv(l, m, theta, phi, dylm_theta, dylm_phi)
                        

                        ! w, wdot, z, at at the radius of this node
                        w_r   =    u_spl(rad_id(i,j,k,ispec))
                        dw_r  = udot_spl(rad_id(i,j,k,ispec))
                        zz_r  =       zz(rad_id(i,j,k,ispec))

                        ! r and theta at the node  
                        unq_r = rstore(i,j,k,ispec)
                        sinth = sin(theta)
                        tanth = tan(theta)

                        ! E_rr: DT98 D.14
                        strain(1,i,j,k,ispec) = iZERO
                        
                        ! E_tt: DT98 D.15
                        strain(2,i,j,k,ispec) = iONE * mf * w_r * (dylm_theta - ylm/tanth) / (unq_r*sinth)

                        ! E_pp: DT98 D.16
                        strain(3,i,j,k,ispec) = - strain(2,i,j,k,ispec)

                        ! E_rt: DT98 D.17
                        strain(4,i,j,k,ispec) = HALF * iONE * mf * zz_r * ylm / sinth

                        ! E_rp: DT98 D.18
                        strain(5,i,j,k,ispec) =  - HALF * zz_r * dylm_theta

                        ! E_tp: DT98 D.19
                        strain(6,i,j,k,ispec) = (dylm_theta/tanth + & 
                                                 (HALF*ll1*ll1 - (mf/sinth)**TWO)*ylm & 
                                                )* w_r / unq_r

                    enddo 
                enddo 
            enddo 
        enddo




    else
        write(*,*)'Error: mode_type must be S or T'
        stop
    endif 


end subroutine


