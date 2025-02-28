module Vanibindings

use iso_c_binding
  
interface
   ! The interface to the CPP code
   integer(c_int) function cudatestfunc(isize) bind(C, name="cudatestfunc")
     use iso_c_binding
     implicit none

     integer(c_int),value :: isize

   end function cudatestfunc
   !



   integer function copythisarraytodevice(hloc, size, varid) bind(C, name="copythisarraytodevice")
     use iso_c_binding
     implicit none
     type(C_PTR), value :: hloc
     integer(c_int),value :: size, varid
    end function 

    integer function copyfromdevice(hloc, size, varid) bind(C, name="copyfromdevice")
    use iso_c_binding
    implicit none
    type(C_PTR), value :: hloc
    integer(c_int),value :: size, varid
   end function 



    integer function allocate_eta_arrays(size) bind(C, name="allocate_eta_arrays")
      use iso_c_binding
      implicit none
      integer(c_int),value :: size
    end function 

    integer function allocate_Cxyz_array(size) bind(C, name="allocate_Cxyz_array")
      use iso_c_binding
      implicit none
      integer(c_int),value :: size
    end function 

    integer function cpp_project_eta_to_gll(model_npoints, nspec, ngll, modelv) bind(C, name="cpp_project_eta_to_gll")
      use iso_c_binding
      implicit none 
      integer(c_int),value :: model_npoints, ngll, nspec
      real(c_double), dimension(5) :: modelv

    end function 


end interface

end module Vanibindings
