VPATH = ../src

# cpp flag is needed for preprocessing #ifdef etc 
AM_FCFLAGS = -cpp -Wall -O2 -llapack -lblas  -I../src -I$(CORETRAN_MODDIR) -L$(CORETRAN_LIBDIR)


LDADD = -lcoretran


# List ALL test programs: 
bin_PROGRAMS = test_plm_function    \
               test_plm_values      \
               test_vani_fast       \
               test_vani_diag       \
               test_tromp1995       \
               test_1_integration   \
               test_ylm_integration \
               test_mode_norm       \
               rotation_semi        \
               semi_whole           \
               test_wrot            \
               test_wstored         \
               test_spline          \
               test_3j              \
               test_xlm             \
               test_ylm_function    \
               test_ylm_values      \
               test_cst_to_mat      \
               tmat_semi            \
               test_tmat            \
               viso_semi            \
               test_viso            \
               test_gmag            \
               vcen_semi            \
               test_amp             \
               reproduce_tromp93    \
               test_surf_int        \
               ell_gll              \
               radial_perturb       \
               test_hemisphere


d_cst  = cst
d_plm  = plm
d_vani = Vani
d_int  = integration
d_norm = mode_normalisation
d_wrm  = rot_mat
d_spl  = spline
d_wig  = wj3
d_xlm  = xlm
d_ylm  = ylm
d_tmat = Tmat
d_viso = Viso
d_grav = gravity
d_vcen = Vcen
d_amp  = amplitude
d_t93  = reproduce_Tromp93
d_ell  = ellipticity

# Define basic sources groups 
SETUP_SRCS = params.f90 allocation.f90 math.f90 mesh_utils.f90  \
             ylm_plm.f90 thrj.f90 

MED_SRCS =  $(SETUP_SRCS) gll.f90 spline.f90 mineos_model.f90 \
                          mode.f90 integrate.f90 gravity.f90 woodhouse.f90 piecewise_interpolation.f90  \
                          lagrange_any.f90 mesh.f90 prem.f90  \
                          output.f90 model3D_reader.f90


VANI_SRCS = $(MED_SRCS) elastic_tensor.f90 cuda_proxies.f90 V_ani_matrix.f90
VISO_SRCS = $(MED_SRCS) cst.f90 $(d_tmat)/rho_st.f90 V_iso_matrix.f90
TMAT_SRCS = $(MED_SRCS) cst.f90 $(d_tmat)/rho_st.f90 T_matrix.f90
INT_SRCS  = $(MED_SRCS)
WMAT_SRCS = $(MED_SRCS) W_matrix.f90 centrifugal.f90 V_cen_matrix.f90
SPL_SRCS  = $(MED_SRCS)
WIG_SRCS  = $(SETUP_SRCS)
SHRM_SRCS = $(SETUP_SRCS)
XLM_SRCS  = $(MED_SRCS)
CST_SRCS  = $(MED_SRCS) cuda_proxies.f90 V_ani_matrix.f90 cst.f90
GRAV_SRCS = $(MED_SRCS) 
VCEN_SRCS = $(MED_SRCS) V_cen_matrix.f90
AMP_SRCS  = $(MED_SRCS) 
T93_SRCS  = $(MED_SRCS) voronoi.f90 cuda_proxies.f90 V_ani_matrix.f90 cst.f90
ELL_SRCS  = $(MED_SRCS) T_matrix.f90 V_cen_matrix.f90

test_plm_function_SOURCES     = $(SETUP_SRCS) $(d_plm)/test_plm.f90
test_plm_values_SOURCES       = $(SETUP_SRCS) $(d_plm)/test_plm_values.f90
test_vani_fast_SOURCES        = $(VANI_SRCS)  $(d_vani)/test_fast_Vani.f90
test_vani_diag_SOURCES        = $(VANI_SRCS)  $(d_vani)/test_Vani.f90
test_tromp1995_SOURCES        = $(VANI_SRCS)  $(d_vani)/tromp_1995.f90
radial_perturb_SOURCES        = $(VANI_SRCS)  $(d_vani)/radial_perturb.f90
test_hemisphere_SOURCES        = $(VANI_SRCS) $(d_vani)/hemisphere.f90
test_1_integration_SOURCES    = $(INT_SRCS)   $(d_int)/test_1_integration.f90
test_surf_int_SOURCES         = $(INT_SRCS)   $(d_int)/test_surface_integration.f90
test_ylm_integration_SOURCES  = $(INT_SRCS)   $(d_int)/test_ylm_integration.f90
test_mode_norm_SOURCES        = $(MED_SRCS)   $(d_norm)/test_mode_normalisation.f90
test_wrot_SOURCES             = $(WMAT_SRCS)  $(d_wrm)/test_rotation_matrix.f90
rotation_semi_SOURCES         = $(WMAT_SRCS)  $(d_wrm)/rotation_semi_analytical.f90
semi_whole_SOURCES            = $(WMAT_SRCS)  $(d_wrm)/test_semi_whole.f90
test_wstored_SOURCES          = $(WMAT_SRCS)  $(d_wrm)/test_stored_Wmat.f90
test_spline_SOURCES           = $(SPL_SRCS)   $(d_spl)/test_spline.f90
test_3j_SOURCES               = $(SETUP_SRCS) $(d_wig)/test_wj3.f90
test_xlm_SOURCES              = $(XLM_SRCS)   $(d_xlm)/test_xlm.f90
test_ylm_function_SOURCES     = $(SHRM_SRCS)  $(d_ylm)/test_ylm.f90
test_ylm_values_SOURCES       = $(SHRM_SRCS)  $(d_ylm)/test_values.f90
test_cst_to_mat_SOURCES       = $(CST_SRCS)   $(d_cst)/test_cst_to_mat.f90
tmat_semi_SOURCES             = $(TMAT_SRCS)  $(d_tmat)/tmat_radial.f90
test_tmat_SOURCES             = $(TMAT_SRCS)  $(d_tmat)/tmat_sem.f90
test_viso_SOURCES             = $(VISO_SRCS)  $(d_viso)/viso_sem.f90
viso_semi_SOURCES             = $(VISO_SRCS)  $(d_viso)/viso_radial.f90
test_gmag_SOURCES             = $(GRAV_SRCS)  $(d_grav)/test_gmag.f90
vcen_semi_SOURCES             = $(VCEN_SRCS)  $(d_vcen)/vcen_woodhouse.f90
test_amp_SOURCES              = $(AMP_SRCS)   $(d_amp)/rotation_amplitude.f90
reproduce_tromp93_SOURCES     = $(T93_SRCS)   $(d_t93)/reproduce_vani.f90
ell_gll_SOURCES               = $(ELL_SRCS)   $(d_ell)/ell_gll.f90


#reproduce_tromp93_LDADD = -lcoretran 
#test_tromp1995_LDADD    = -lcoretran 

# Test groups: 
plm_tests:  test_ylm_function test_ylm_values
ylm_tests:  test_ylm_function test_plm_values
vani_tests: test_vani_fast test_vani_diag test_tromp1995
integration_tests: test_1_integration test_ylm_integration test_surf_int
wmat_tests: test_wstored rotation_semi test_wrot
cst_tests: test_cst_to_mat
viso_tests: test_viso viso_semi
grav_tests: test_gmag

.PHONY: plm_tests vani_tests integration_tests ylm_tests wmat_tests \
        cst_tests viso_tests grav_tests



