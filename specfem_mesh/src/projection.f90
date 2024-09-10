
subroutine get_mesh_radii()
    ! Determines a list of the unique radii of the GLL point
    ! Points each GLL point to its radial value in the list 
    ! 
    use params, only: ngllx, nglly, ngllz, nspec, RA, &
                      xstore, ystore, zstore, unique_r, n_unique_rad, rad_id,&
                      interp_id_r, IC_ID, verbose
    use allocation_module, only: allocate_if_unallocated
    use spline, only: create_interpolation_radial_map
    implicit none 
    include "constants.h"

    ! Local variables: 
    real(kind=CUSTOM_REAL) :: rr(ngllx, nglly, ngllz, nspec) , tol, rgll
    real(kind=CUSTOM_REAL), allocatable :: unique_r_tmp(:)
    integer :: i, j, k, ispec, size_r, unique_id, ur
    
    logical :: match

    ! Compute radii
    rr  = (xstore**TWO + ystore**TWO + zstore**TWO)**HALF

    ! Maximum number of unique radii is this
    size_r =  size(rr)
    allocate(unique_r_tmp(size_r))
    unique_r_tmp = -ONE


    ! Set tolerance for 'same radii' 
    tol = 1.0e-4

    ! Initial value
    unique_id = 1

    ! Allocate the array that will store the unique radius IDs: 
    call allocate_if_unallocated(ngllx, nglly, ngllz, nspec, rad_id)

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
                        !write(*,*)'unique_id ', unique_id
                    endif
                enddo 
            enddo 
        enddo 
    enddo 




    ! Allocate new unique radius array of the correct size
    ! It might be that the very last node tested in unique and so 
    ! unique_id is actually 1 too large -- test if last radius is -1 and if so cut it 
    if (unique_r_tmp(unique_id).eq.-ONE)then 
        n_unique_rad = unique_id - 1
    else
        n_unique_rad = unique_id
    endif

    ! Copy over the unique values to the smaller array and deallocate long.
    allocate(unique_r(n_unique_rad)) 
    unique_r(1:n_unique_rad) = unique_r_tmp(1:n_unique_rad)
    deallocate(unique_r_tmp)

    ! Print some details
    if(verbose.ge.2) then 
        write(*,'(/,a)')'• Determine unique mesh radii'
        write(*,'(a,es8.1)')  '  -- tolerance value       :', tol
        write(*,'(a,f8.2,a)') '  -- matches radii within  :', tol*RA, ' metres'
        write(*,'(a,i8,a,i8)')'  -- number of unique radii:', n_unique_rad, ' out of ', size_r
    endif 

    if(minval(rad_id).le.0 .or. maxval(rad_id).gt.n_unique_rad)then 
        write(*,*)'Error with rad_id:' 
        write(*,'(a,i8)')  '  -- min id       :', minval(rad_id)
        write(*,'(a,i8)')  '  -- max id       :', maxval(rad_id)
        write(*,'(a,i8)')  '  -- n_unique_rad :', n_unique_rad
        stop
    endif        




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


    allocate(interp_id_r(n_unique_rad))
    call create_interpolation_radial_map(unique_r, interp_id_r, n_unique_rad, 1, IC_ID)



end subroutine get_mesh_radii





subroutine compute_global_mode_displacement(mode_type, l, m, u_spline, v_spline, n_spl, disp)
    use params, only:  nspec, ngllx, & 
    nglly, ngllz, thetastore, phistore, rad_id
    use spline, only: interpolate_mode_eigenfunctions
    use ylm_plm, only: ylm_complex, ylm_deriv
    use mesh_utils, only: delta_spline
    use math, only: sinp, sqrtp

    implicit none
    include "constants.h"

    ! IO variables: 
    character(len=1) :: mode_type     ! mode type; S or T
    integer          :: l             ! l value (degree) of mode
    integer          :: m             ! m value (order)  of mode
    integer          :: n_spl         ! length of spline
    complex(kind=SPLINE_REAL) :: disp(3,ngllx, nglly, ngllz, nspec)
    
    real(kind=SPLINE_REAL) :: u_spline(n_spl),  v_spline(n_spl)

    ! Local variables: 
    complex(kind=CUSTOM_REAL) :: ylm, dylm_theta, dylm_phi
    real(kind=CUSTOM_REAL)    ::  theta, phi 
    integer :: i,j,k,ispec

    real(kind=SPLINE_REAL)  ::  mf, lf, ll1, tl14p, mone_l, pref, u_r, v_r, w_r, dm0, dd1, sinth
    complex(kind=SPLINE_REAL) :: sp_ylm, spl_dylm_theta

    disp = SPLINE_iZERO

    ! DT98 below D.1: k = sqrt(l(l+1))
    lf  = real(l, kind=SPLINE_REAL)
    mf  = real(m, kind=SPLINE_REAL)


    ll1    = sqrtp(lf*(lf+SPLINE_ONE))                                              ! (l(l+1))^1/2
    tl14p  = ((SPLINE_TWO*lf + SPLINE_ONE)/(SPLINE_FOUR*SPLINE_PI))**SPLINE_HALF    ! (2l+1/4π)^1/2
    dm0    = delta_spline(m, 0)                                                     ! δ_{m0}
    dd1    = delta_spline(m, -1) - delta_spline(m, 1)                               ! δ_{m -1} - δ_{m 1}
    mone_l = (-SPLINE_ONE)**lf                                                      ! -1 ^l 
    pref   = SPLINE_HALF * tl14p * ll1                                              ! 0.5 (2l+1/4π)^1/2  (l(l+1))^1/2


    if (mode_type.eq.'S')then 
        do ispec = 1, nspec
            do i = 1, ngllx
                do j = 1, nglly
                    do k = 1, ngllz
                    
                        ! Get ylm and the partial derivatives
                        theta = real(thetastore(i,j,k,ispec), kind=CUSTOM_REAL)
                        phi   = real(phistore(i,j,k,ispec),   kind=CUSTOM_REAL)
                        ylm   = ylm_complex(l, m, theta, phi)
                        call ylm_deriv(l, m, theta, phi, dylm_theta, dylm_phi)

                        ! u, v at at the radius of this node
                        u_r   =    u_spline(rad_id(i,j,k,ispec))
                        v_r   =    v_spline(rad_id(i,j,k,ispec)) / ll1

                        sinth = real(sinp(theta), kind=SPLINE_REAL)

                        if (theta.ge.zero .and. theta.le.pole_tolerance) then 
                            ! North pole 
                            ! S_r     (DT98 D.8)
                            disp(1,i,j,k,ispec) = tl14p * u_r * dm0
                            ! S_theta (DT98 D.9)
                            disp(2,i,j,k,ispec) =  pref * v_r * dd1
                            ! S_phi   (DT98 D.10)
                            disp(3,i,j,k,ispec) = pref * SPLINE_iONE * mf * v_r * dd1


                        elseif (abs(theta-PI).le.pole_tolerance) then
                            ! South pole
                            ! S_r     (DT98 D.11)
                            disp(1,i,j,k,ispec) =  mone_l * tl14p * u_r * dm0
                            ! S_theta (DT98 D.12)
                            disp(2,i,j,k,ispec) =  mone_l * pref * v_r * dd1
                            ! S_phi   (DT98 D.13)
                            disp(3,i,j,k,ispec) = -mone_l * pref * SPLINE_iONE * mf * v_r * dd1
                        else 
                            ! Convert to SPLINE_REAL precision
                            sp_ylm         = cmplx(ylm,        kind=SPLINE_REAL)
                            spl_dylm_theta = cmplx(dylm_theta, kind=SPLINE_REAL)

                            ! S_r     (DT98 D.4)
                            disp(1,i,j,k,ispec) = sp_ylm*u_r  
                            ! S_theta (DT98 D.5)
                            disp(2,i,j,k,ispec) = v_r * spl_dylm_theta
                            ! S_phi   (DT98 D.6)
                            disp(3,i,j,k,ispec) = SPLINE_iONE * mf * v_r * sp_ylm / sinth
                        endif 

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
                        
                        w_r   = u_spline(rad_id(i,j,k,ispec))/ll1

                        theta = real(thetastore(i,j,k,ispec), kind=CUSTOM_REAL)
                        phi   = real(phistore(i,j,k,ispec), kind=CUSTOM_REAL)
                        ylm   = ylm_complex(l,m, theta, phi)
                        call ylm_deriv(l, m, theta, phi, dylm_theta, dylm_phi)

                        if (theta.ge.zero .and. theta.le.pole_tolerance) then 
                            ! North pole 
                            ! S_theta (DT98 D.9)
                            disp(2,i,j,k,ispec) = pref * SPLINE_iONE * mf * w_r * dd1
                            ! S_phi   (DT98 D.10)
                            disp(3,i,j,k,ispec) = - pref * w_r * dd1
                            
                        elseif (abs(theta-PI).le.pole_tolerance) then
                            ! South pole
                            ! S_theta (DT98 D.12)
                            disp(2,i,j,k,ispec) = - mone_l * pref * SPLINE_iONE * mf * w_r * dd1

                            ! S_phi   (DT98 D.13)
                            disp(3,i,j,k,ispec) = - mone_l * pref * w_r * dd1
                        else 

                            ! Convert to SPLINE_REAL precision
                            sp_ylm         = cmplx(ylm,        kind=SPLINE_REAL)
                            spl_dylm_theta = cmplx(dylm_theta, kind=SPLINE_REAL)

                            ! S_theta (DT98 D.5)
                            disp(2,i,j,k,ispec) = SPLINE_iONE * mf * w_r * sp_ylm / sinth
                            ! S_phi   (DT98 D.6)
                            disp(3,i,j,k,ispec) = -w_r * spl_dylm_theta 
                        endif 

                        ! No radial displacement for toroidal modes
                        disp(1,i,j,k,ispec) = SPLINE_iZERO
                    enddo 
                enddo 
            enddo 
        enddo 

    else
        write(*,*)'ERROR: Mode type must be S or T but was ', mode_type
        stop
    endif 

end subroutine compute_global_mode_displacement





subroutine compute_gll_mode_strain(mode_type, nord, l, m, usp, vsp, udotsp, vdotsp, spl_len, strain)
    ! Strain ordering is as follows: 
    ! 1  =  E_rr     4  =  E_tp
    ! 2  =  E_tt     5  =  E_rp
    ! 3  =  E_pp     6  =  E_rt
    ! Note that this is differnet from what I would usually do as an order (e.g. that of a moment tensor)
    ! This is voigt notation convention

    use params, only: nspec, ngllx, unique_r, & 
                      nglly, ngllz, thetastore, phistore, & 
                      rad_id, rstore, all_warnings
    use ylm_plm, only: ylm_complex, ylm_deriv
    use spline, only: interpolate_mode_eigenfunctions
    use mesh_utils, only: delta_spline
    use math, only: sinp, tanp, sqrtp
    use allocation_module, only: allocate_if_unallocated

    implicit none
    include "constants.h"

    ! IO variables: 
    character(len=1) :: mode_type     ! mode type; S or T
    integer          :: nord          ! n value of mode
    integer          :: l             ! l value (degree) of mode
    integer          :: m             ! m value (order)  of mode
    integer          :: spl_len       ! length of splined eigenfunctions
    
    ! Spline-interpolated eigenfunctions
    real(kind=SPLINE_REAL)    :: usp(spl_len)
    real(kind=SPLINE_REAL)    :: vsp(spl_len)
    real(kind=SPLINE_REAL)    :: udotsp(spl_len)
    real(kind=SPLINE_REAL)    :: vdotsp(spl_len)
    complex(kind=SPLINE_REAL) :: strain(6,ngllx, nglly, ngllz, nspec)

    ! Local variables: 
    complex(kind=CUSTOM_REAL) :: ylm, dylm_theta, dylm_phi
    real(kind=CUSTOM_REAL)    :: theta, phi
    complex(kind=SPLINE_REAL) :: sp_ylm, spl_dylm_theta

    real(kind=SPLINE_REAL)    :: sinth, tanth, unq_r, du_r, u_r, & 
                                 dv_r, v_r, w_r, dw_r, ff, &
                                 dm0, dd1, dd2

    ! Auxillary eigenfunction variables if needed
    real(kind=SPLINE_REAL), allocatable :: xx(:), zz(:)
    
    ! GLL level variables
    real(kind=SPLINE_REAL)  :: xx_r, zz_r, mf, lf, ll1, tl14p, kr2
    integer :: rid 
    ! Loop variables 
    integer :: ispec, i,j,k, h

    lf  = real(l, kind=SPLINE_REAL)           ! float l 
    mf  = real(m, kind=SPLINE_REAL)           ! float m
    ll1 = (lf*(lf+SPLINE_ONE))**SPLINE_HALF   ! float of k = √(l(l+1))
    tl14p = ((SPLINE_TWO*lf + SPLINE_ONE)/(SPLINE_FOUR*SPLINE_PI))**SPLINE_HALF    ! (2l+1/4π)^1/2
    kr2 = ll1*((ll1*ll1 - SPLINE_TWO)**SPLINE_HALF)      

    dm0 = delta_spline(m, 0)                         ! δ_{m0}
    dd1 = delta_spline(m, -1) - delta_spline(m, 1)   ! δ_{m -1} - δ_{m 1}
    dd2 = delta_spline(m, -2) + delta_spline(m, 2)   ! δ_{m -2} + δ_{m 2}

    ! Convert eigenfunctions to auxillary form
    ! DT98 D.1: 
    !    u = U     v = V/k     w = W/k     p = P

    ! Compute x and z auxillary variables - DT98 D.20
    ! If spheroidal mode then u,v are non zero but w, wdot are 0
    !   --> z should be 0 
    ! If toroidal mode then u,v are zero but w, wdot are not
    !   --> x should be 0
    if (mode_type.eq.'S')then 
        allocate(xx(spl_len))
        xx = vdotsp/ll1 + (usp - vsp/ll1)/real(unique_r, kind=SPLINE_REAL)

        ! Loop over each GLL point: 
        do ispec = 1, nspec
            do i = 1, ngllx
                do j = 1, nglly
                    do k = 1, ngllz
                        ! Get ylm and the partial derivatives
                        theta = thetastore(i,j,k,ispec)
                        phi   = phistore(i,j,k,ispec)
                        ylm   = ylm_complex(l, m, theta, phi)
                        call ylm_deriv(l, m, theta, phi, dylm_theta, dylm_phi)
                        
                        ! u, v, udot, vdo, x, at at the radius of this node
                        rid   = rad_id(i,j,k,ispec)
                        u_r   = usp(rid)
                        du_r  = udotsp(rid)
                        
                        v_r   = vsp(rid)/ll1
                        xx_r  = xx(rid)

                        ! r and theta at the node  
                        unq_r = real(rstore(i,j,k,ispec), kind=SPLINE_REAL)
                        sinth = real(sinp(theta),         kind=SPLINE_REAL)
                        tanth = real(tanp(theta),         kind=SPLINE_REAL)


                        if(unq_r.eq.0)then 
                            if(all_warnings) write(*,*)'Warning: setting strain at centre to 0 artificially'

                            strain(:,i,j,k,ispec) = SPLINE_iZERO
                        else

                            if (theta.ge.zero .and. theta.le.pole_tolerance) then 
                                ! Values at the pole 
                                ! D.28
                                ff = (SPLINE_TWO*u_r - ll1*ll1*v_r)/unq_r    
                                
                                ! E_rr: DT98 D.22
                                strain(1,i,j,k,ispec) = tl14p * du_r * dm0

                                ! E_tt: DT98 D.23
                                strain(2,i,j,k,ispec) = SPLINE_HALF * tl14p * &
                                                        (ff*dm0 + & 
                                                        kr2*SPLINE_TWO*v_r * & 
                                                        dd2/(SPLINE_FOUR*unq_r))

                                ! E_pp: DT98 D.24
                                strain(3,i,j,k,ispec) = SPLINE_HALF * tl14p * &
                                                        (ff*dm0 - & 
                                                        kr2*SPLINE_TWO*v_r * & 
                                                        dd2/(SPLINE_FOUR*unq_r))

                                ! E_rt: DT98 D.25
                                strain(6,i,j,k,ispec) = tl14p * ll1 * xx_r * dd1/SPLINE_FOUR
                                ! E_rp: DT98 D.26
                                strain(5,i,j,k,ispec) = tl14p * ll1 * SPLINE_iONE * mf * xx_r/SPLINE_FOUR
                                ! E_tp: DT98 D.27
                                strain(4,i,j,k,ispec) = tl14p * kr2 * SPLINE_iONE * mf * v_r * dd2/(SPLINE_EIGHT*unq_r)
                            else 

                                ! Convert to SPLINE_REAL precision
                                sp_ylm         = cmplx(ylm, kind=SPLINE_REAL)
                                spl_dylm_theta = cmplx(dylm_theta, kind=SPLINE_REAL)

                                ! Values away from the pole
                                ! E_rr: DT98 D.14
                                strain(1,i,j,k,ispec) = sp_ylm*du_r  
                                ! E_tt: DT98 D.15
                                strain(2,i,j,k,ispec) = (sp_ylm*u_r - v_r*(spl_dylm_theta/tanth -  sp_ylm*(mf/sinth)**SPLINE_TWO + ll1*ll1*sp_ylm))/unq_r
                                ! E_pp: DT98 D.16
                                strain(3,i,j,k,ispec) = (sp_ylm*u_r + v_r*(spl_dylm_theta/tanth -  sp_ylm*(mf/sinth)**SPLINE_TWO))/unq_r
                                ! E_rt: DT98 D.17
                                strain(6,i,j,k,ispec) = SPLINE_HALF * xx_r * spl_dylm_theta
                                ! E_rp: DT98 D.18
                                strain(5,i,j,k,ispec) = SPLINE_HALF * SPLINE_iONE * mf * xx_r * sp_ylm/sinth
                                ! E_tp: DT98 D.19
                                strain(4,i,j,k,ispec) = SPLINE_iONE * mf * v_r * (spl_dylm_theta - sp_ylm/tanth) / (unq_r * sinth)
                            endif        
                        endif !unique r



                        do h = 1, 6
                            if (abs(strain(h,i,j,k,ispec)).gt. 100.0)then
                                strain(h,i,j,k,ispec) = SPLINE_ZERO
                                write(*,*)'Setting strain values at pi to 0,'
                                !stop 
                            endif 
                        enddo 

                    enddo 
                enddo 
            enddo 
        enddo

    elseif (mode_type.eq.'T')then 
        ! I believe that get_mode stores w, wdot in the u and du terms
        allocate(zz(spl_len))
        zz = (udotsp - usp/real(unique_r, kind=SPLINE_REAL))/ll1

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
                        rid   = rad_id(i,j,k,ispec)
                        w_r   = usp(rid)/ll1
                        dw_r  = udotsp(rid)/ll1
                        zz_r  = zz(rid) ! No division as divided above

                        ! r and theta at the node  
                        unq_r = real(rstore(i,j,k,ispec), kind=SPLINE_REAL)
                        sinth = real(sinp(theta),         kind=SPLINE_REAL)
                        tanth = real(tanp(theta),         kind=SPLINE_REAL)

                        if (theta.ge.zero .and. theta.le.pole_tolerance) then 
                            ! Values at the pole 
                            ! E_rr: DT98 D.22
                            strain(1,i,j,k,ispec) = SPLINE_iZERO
                            ! E_tt: DT98 D.23
                            strain(2,i,j,k,ispec) = SPLINE_HALF * tl14p * SPLINE_iONE * mf * w_r * dd2 * kr2/(SPLINE_FOUR * unq_r)
                            ! E_pp: DT98 D.24
                            strain(3,i,j,k,ispec) = -SPLINE_HALF * tl14p * SPLINE_iONE * mf * w_r * dd2 * kr2/(SPLINE_FOUR * unq_r)
                            ! E_rt: DT98 D.25
                            strain(6,i,j,k,ispec) = tl14p * ll1 * SPLINE_iONE * mf * zz_r * dd1 / SPLINE_FOUR
                            ! E_rp: DT98 D.26
                            strain(5,i,j,k,ispec) = - tl14p * ll1 * zz_r * dd1 / SPLINE_FOUR
                            ! E_tp: DT98 D.27
                            strain(4,i,j,k,ispec) = - tl14p * kr2 * w_r * dd2 / (unq_r * SPLINE_FOUR)
                        else

                            ! Convert to SPLINE_REAL precision
                            sp_ylm         = cmplx(ylm, kind=SPLINE_REAL)
                            spl_dylm_theta = cmplx(dylm_theta, kind=SPLINE_REAL)

                            ! E_rr: DT98 D.14
                            strain(1,i,j,k,ispec) = SPLINE_iZERO
                            ! E_tt: DT98 D.15
                            strain(2,i,j,k,ispec) = SPLINE_iONE * mf * w_r * (spl_dylm_theta - sp_ylm/tanth) / (unq_r*sinth)
                            ! E_pp: DT98 D.16
                            strain(3,i,j,k,ispec) = - strain(2,i,j,k,ispec)
                            ! E_rt: DT98 D.17
                            strain(6,i,j,k,ispec) = SPLINE_HALF * SPLINE_iONE * mf * zz_r * sp_ylm / sinth
                            ! E_rp: DT98 D.18
                            strain(5,i,j,k,ispec) =  - SPLINE_HALF * zz_r * spl_dylm_theta
                            ! E_tp: DT98 D.19
                            strain(4,i,j,k,ispec) = (spl_dylm_theta/tanth + & 
                                                    (SPLINE_HALF*ll1*ll1 - (mf/sinth)**SPLINE_TWO)*sp_ylm & 
                                                    )* w_r / unq_r
                        endif 
                    enddo 
                enddo 
            enddo 
        enddo
    else
        write(*,*)'Error: mode_type must be S or T'
        stop
    endif 


end subroutine compute_gll_mode_strain




subroutine gll_strain_from_disp(disp, strain)
    ! Uses xyz displacement (d) to compute strain 
    use params, only: ngllx, nglly, ngllz, nspec, dgll, jacinv
    implicit none 
    include "constants.h"

    complex(SPLINE_REAL) :: disp(3, ngllx, nglly, ngllz, nspec)
    complex(SPLINE_REAL) :: strain(6, ngllx, nglly, ngllz, nspec)

    ! Local: 
    integer :: ispec, i, j, k, p, q, g, h
    complex(SPLINE_REAL) :: sum_c(3), strn_tmp(3,3)

    strain = SPLINE_iZERO


        ! Compute gradient of the displacement
    do ispec = 1, nspec
        do i = 1, ngllx 
            do j = 1, nglly
                do k = 1, ngllz

                    ! Loop over strain elementsa
                    do p = 1, 3
                        do q = 1, 3

                            sum_c(1) = SPLINE_iZERO
                            sum_c(2) = SPLINE_iZERO
                            sum_c(3) = SPLINE_iZERO

                            do g = 1, ngllx
                                sum_c(1) = sum_c(1) + disp(q, g, j, k, ispec) * real(dgll(g, i), kind=SPLINE_REAL)
                                sum_c(2) = sum_c(2) + disp(q, i, g, k, ispec) * real(dgll(g, j), kind=SPLINE_REAL)
                                sum_c(3) = sum_c(3) + disp(q, i, j, g, ispec) * real(dgll(g, k), kind=SPLINE_REAL)
                            enddo 

                            strn_tmp(q,p) = sum_c(1) * jacinv(1,p,i,j,k,ispec) +   &  
                                            sum_c(2) * jacinv(2,p,i,j,k,ispec) +   & 
                                            sum_c(3) * jacinv(3,p,i,j,k,ispec) 
                        enddo !q
                    enddo !p

                    ! Symmetric part of tensor
                    strn_tmp(:,:) =  SPLINE_HALF * (strn_tmp + transpose(strn_tmp) )

                    ! Save in voigt notation 
                    do h = 1, 3
                        strain(h, i, j, k, ispec) = strn_tmp(h,h)
                    enddo 
                    strain(4, i, j, k, ispec) = strn_tmp(2,3)
                    strain(5, i, j, k, ispec) = strn_tmp(1,3)
                    strain(6, i, j, k, ispec) = strn_tmp(1,2)


                enddo 
            enddo 
        enddo 
    enddo




end subroutine gll_strain_from_disp
