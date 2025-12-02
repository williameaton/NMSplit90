module Vanibindings

use iso_c_binding
  
interface


    integer function fh_cuda_preinit() bind(C, name="fh_cuda_preinit")
    use iso_c_binding
      implicit none
    end function 

    integer function get_cpp_precision() bind(C, name="get_cpp_precision")
    use iso_c_binding
      implicit none
    end function 
    


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



    integer function allocate_Vani_arrays(size) bind(C, name="allocate_Vani_arrays")
      use iso_c_binding
      implicit none
      integer(c_int),value :: size
    end function 


    integer function allocate_eta_arrays(size) bind(C, name="allocate_eta_arrays")
      use iso_c_binding
      implicit none
      integer(c_int),value :: size
    end function 

    integer function copy_LUT_array(hloc, size) bind(C, name="copy_LUT_array")
      use iso_c_binding
      implicit none
      type(C_PTR),   value :: hloc
      integer(c_int),value :: size
    end function 
    


    integer function allocate_M3D_array(size) bind(C, name="allocate_M3D_array")
      use iso_c_binding
      implicit none
      integer(c_int),value :: size
    end function 


    integer function copy_M3D_array(hloc, size) bind(C, name="copy_M3D_array")
      use iso_c_binding
      implicit none
      type(C_PTR),   value :: hloc
      integer(c_int),value :: size
    end function 
  

    integer function assign_proc_to_device(nprocs, myrank) bind(C, name="assign_proc_to_device")
      use iso_c_binding
      implicit none
      integer(c_int),value :: nprocs, myrank
    end function 




    integer function copy_wgll_array(hloc, size) bind(C, name="copy_wgll_array")
      use iso_c_binding
      implicit none
      type(C_PTR),   value :: hloc
      integer(c_int),value :: size
    end function 

    integer function copy_allstrains(hloc_r, hloc_i, size) bind(C, name="copy_allstrains")
      use iso_c_binding
      implicit none
      type(C_PTR),   value :: hloc_r, hloc_i
      integer(c_int64_t),value :: size
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



    integer function launch_vanikernel(ngll, nspec, nn1_total, maxnn1,  maxtl1, nmodes, myrank) bind(C, name="launch_vanikernel")
      use iso_c_binding
      implicit none 
      integer(c_int),value :: ngll, nspec, nn1_total, maxtl1, maxnn1, nmodes, myrank
    end function 



end interface

end module Vanibindings
