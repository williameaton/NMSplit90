VPATH = ../src

# cpp flag is needed for preprocessing #ifdef etc 
AM_FCFLAGS = -cpp -Wall -O2 -llapack -lblas  -I../src


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
               test_wrot            \
               test_wstored         \
               test_spline          \
               test_3j              \
               test_xlm             \
               test_ylm_function    \
               test_ylm_values      \
               test_cst_to_mat

d_cst  = cst
d_plm  = plm
d_vani = v_ani_matrix
d_int  = integration
d_norm = mode_normalisation
d_wrm  = rot_mat
d_spl  = spline
d_wig  = wj3
d_xlm  = xlm
d_ylm  = ylm


# Define basic sources groups 
SETUP_SRCS = params.f90 allocation.f90 math.f90 mesh_utils.f90  \
             ylm_plm.f90 thrj.f90

MED_SRCS =  $(SETUP_SRCS) gll.f90 spline.f90 mineos_model.f90 \
                          mode.f90 piecewise_interpolation.f90  \
                          lagrange_any.f90 mesh.f90 integrate.f90  output.f90 


VANI_SRCS = $(MED_SRCS) cuda_proxies.f90 V_ani_matrix.f90
INT_SRCS  = $(MED_SRCS)
WMAT_SRCS = $(MED_SRCS) W_matrix.f90
SPL_SRCS  = $(MED_SRCS)
WIG_SRCS  = $(SETUP_SRCS)
SHRM_SRCS = $(SETUP_SRCS)
XLM_SRCS  = $(MED_SRCS)
CST_SRCS  = $(MED_SRCS) cuda_proxies.f90 V_ani_matrix.f90 cst.f90

test_plm_function_SOURCES     = $(SETUP_SRCS) $(d_plm)/test_plm.f90
test_plm_values_SOURCES       = $(SETUP_SRCS) $(d_plm)/test_plm_values.f90
test_vani_fast_SOURCES        = $(VANI_SRCS)  $(d_vani)/test_fast_Vani.f90
test_vani_diag_SOURCES        = $(VANI_SRCS)  $(d_vani)/test_Vani.f90
test_tromp1995_SOURCES        = $(VANI_SRCS)  $(d_vani)/tromp_1995.f90
test_1_integration_SOURCES    = $(INT_SRCS)   $(d_int)/test_1_integration.f90
test_ylm_integration_SOURCES  = $(INT_SRCS)   $(d_int)/test_ylm_integration.f90
test_mode_norm_SOURCES        = $(MED_SRCS)   $(d_norm)/test_mode_normalisation.f90
test_wrot_SOURCES             = $(WMAT_SRCS)  $(d_wrm)/test_rotation_matrix.f90
rotation_semi_SOURCES         = $(WMAT_SRCS)  $(d_wrm)/rotation_semi_analytical.f90
test_wstored_SOURCES          = $(WMAT_SRCS)  $(d_wrm)/test_stored_Wmat.f90
test_spline_SOURCES           = $(SPL_SRCS)   $(d_spl)/test_spline.f90
test_3j_SOURCES               = $(SETUP_SRCS) $(d_wig)/test_wj3.f90
test_xlm_SOURCES              = $(XLM_SRCS)   $(d_xlm)/test_xlm.f90
test_ylm_function_SOURCES     = $(SHRM_SRCS)  $(d_ylm)/test_ylm.f90
test_ylm_values_SOURCES       = $(SHRM_SRCS)  $(d_ylm)/test_values.f90
test_cst_to_mat_SOURCES       = $(CST_SRCS)   $(d_cst)/test_cst_to_mat.f90

# Test groups: 
plm_tests:  test_ylm_function test_ylm_values
ylm_tests:  test_ylm_function test_plm_values
vani_tests: test_vani_fast test_vani_diag test_tromp1995
integration_tests: test_1_integration test_ylm_integration
wmat_tests: test_wstored rotation_semi test_wrot
cst_tests: test_cst_to_mat

.PHONY: plm_tests vani_tests integration_tests ylm_tests wmat_tests cst_tests



