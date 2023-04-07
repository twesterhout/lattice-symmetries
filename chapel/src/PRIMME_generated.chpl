// Generated with c2chapel version 0.1.0

// Header given to c2chapel:
require "primme_headers/primme.h";

// Note: Generated with fake std headers

use CPtr;
use SysCTypes;
use SysBasic;
// #define'd integer literals:
// Note: some of these may have been defined with an ifdef
extern const PRIMME_VERSION_MAJOR : int;
extern const PRIMME_VERSION_MINOR : int;

// End of #define'd integer literals

extern "struct _primme_half" record _primme_half {
  // var a : int short;
}

extern "struct _primme_complex_half" record _primme_complex_half {
  var r : _primme_half;
  var i : _primme_half;
}

extern proc hprimme(ref evals : _primme_half, ref evecs : _primme_half, ref resNorms : _primme_half, ref primme : primme_params) : c_int;

extern proc hprimme(evals : c_ptr(_primme_half), evecs : c_ptr(_primme_half), resNorms : c_ptr(_primme_half), primme : c_ptr(primme_params)) : c_int;

extern proc kprimme(ref evals : _primme_half, ref evecs : _primme_complex_half, ref resNorms : _primme_half, ref primme : primme_params) : c_int;

extern proc kprimme(evals : c_ptr(_primme_half), evecs : c_ptr(_primme_complex_half), resNorms : c_ptr(_primme_half), primme : c_ptr(primme_params)) : c_int;

extern proc sprimme(ref evals : c_float, ref evecs : c_float, ref resNorms : c_float, ref primme : primme_params) : c_int;

extern proc sprimme(evals : c_ptr(c_float), evecs : c_ptr(c_float), resNorms : c_ptr(c_float), primme : c_ptr(primme_params)) : c_int;

extern proc cprimme(ref evals : c_float, ref evecs : _Complex float, ref resNorms : c_float, ref primme : primme_params) : c_int;

extern proc cprimme(evals : c_ptr(c_float), evecs : c_ptr(_Complex float), resNorms : c_ptr(c_float), primme : c_ptr(primme_params)) : c_int;

extern proc dprimme(ref evals : c_double, ref evecs : c_double, ref resNorms : c_double, ref primme : primme_params) : c_int;

extern proc dprimme(evals : c_ptr(c_double), evecs : c_ptr(c_double), resNorms : c_ptr(c_double), primme : c_ptr(primme_params)) : c_int;

extern proc zprimme(ref evals : c_double, ref evecs : _Complex double, ref resNorms : c_double, ref primme : primme_params) : c_int;

extern proc zprimme(evals : c_ptr(c_double), evecs : c_ptr(_Complex double), resNorms : c_ptr(c_double), primme : c_ptr(primme_params)) : c_int;

extern proc magma_hprimme(ref evals : _primme_half, ref evecs : _primme_half, ref resNorms : _primme_half, ref primme : primme_params) : c_int;

extern proc magma_hprimme(evals : c_ptr(_primme_half), evecs : c_ptr(_primme_half), resNorms : c_ptr(_primme_half), primme : c_ptr(primme_params)) : c_int;

extern proc magma_kprimme(ref evals : _primme_half, ref evecs : _primme_complex_half, ref resNorms : _primme_half, ref primme : primme_params) : c_int;

extern proc magma_kprimme(evals : c_ptr(_primme_half), evecs : c_ptr(_primme_complex_half), resNorms : c_ptr(_primme_half), primme : c_ptr(primme_params)) : c_int;

extern proc magma_sprimme(ref evals : c_float, ref evecs : c_float, ref resNorms : c_float, ref primme : primme_params) : c_int;

extern proc magma_sprimme(evals : c_ptr(c_float), evecs : c_ptr(c_float), resNorms : c_ptr(c_float), primme : c_ptr(primme_params)) : c_int;

extern proc magma_cprimme(ref evals : c_float, ref evecs : _Complex float, ref resNorms : c_float, ref primme : primme_params) : c_int;

extern proc magma_cprimme(evals : c_ptr(c_float), evecs : c_ptr(_Complex float), resNorms : c_ptr(c_float), primme : c_ptr(primme_params)) : c_int;

extern proc magma_dprimme(ref evals : c_double, ref evecs : c_double, ref resNorms : c_double, ref primme : primme_params) : c_int;

extern proc magma_dprimme(evals : c_ptr(c_double), evecs : c_ptr(c_double), resNorms : c_ptr(c_double), primme : c_ptr(primme_params)) : c_int;

extern proc magma_zprimme(ref evals : c_double, ref evecs : _Complex double, ref resNorms : c_double, ref primme : primme_params) : c_int;

extern proc magma_zprimme(evals : c_ptr(c_double), evecs : c_ptr(_Complex double), resNorms : c_ptr(c_double), primme : c_ptr(primme_params)) : c_int;

extern proc hsprimme(ref evals : c_float, ref evecs : _primme_half, ref resNorms : c_float, ref primme : primme_params) : c_int;

extern proc hsprimme(evals : c_ptr(c_float), evecs : c_ptr(_primme_half), resNorms : c_ptr(c_float), primme : c_ptr(primme_params)) : c_int;

extern proc ksprimme(ref evals : c_float, ref evecs : _primme_complex_half, ref resNorms : c_float, ref primme : primme_params) : c_int;

extern proc ksprimme(evals : c_ptr(c_float), evecs : c_ptr(_primme_complex_half), resNorms : c_ptr(c_float), primme : c_ptr(primme_params)) : c_int;

extern proc magma_hsprimme(ref evals : c_float, ref evecs : _primme_half, ref resNorms : c_float, ref primme : primme_params) : c_int;

extern proc magma_hsprimme(evals : c_ptr(c_float), evecs : c_ptr(_primme_half), resNorms : c_ptr(c_float), primme : c_ptr(primme_params)) : c_int;

extern proc magma_ksprimme(ref evals : c_float, ref evecs : _primme_complex_half, ref resNorms : c_float, ref primme : primme_params) : c_int;

extern proc magma_ksprimme(evals : c_ptr(c_float), evecs : c_ptr(_primme_complex_half), resNorms : c_ptr(c_float), primme : c_ptr(primme_params)) : c_int;

extern proc kprimme_normal(ref evals : _primme_complex_half, ref evecs : _primme_complex_half, ref resNorms : _primme_half, ref primme : primme_params) : c_int;

extern proc kprimme_normal(evals : c_ptr(_primme_complex_half), evecs : c_ptr(_primme_complex_half), resNorms : c_ptr(_primme_half), primme : c_ptr(primme_params)) : c_int;

extern proc cprimme_normal(ref evals : _Complex float, ref evecs : _Complex float, ref resNorms : c_float, ref primme : primme_params) : c_int;

extern proc cprimme_normal(evals : c_ptr(_Complex float), evecs : c_ptr(_Complex float), resNorms : c_ptr(c_float), primme : c_ptr(primme_params)) : c_int;

extern proc zprimme_normal(ref evals : _Complex double, ref evecs : _Complex double, ref resNorms : c_double, ref primme : primme_params) : c_int;

extern proc zprimme_normal(evals : c_ptr(_Complex double), evecs : c_ptr(_Complex double), resNorms : c_ptr(c_double), primme : c_ptr(primme_params)) : c_int;

extern proc magma_kprimme_normal(ref evals : _primme_complex_half, ref evecs : _primme_complex_half, ref resNorms : _primme_half, ref primme : primme_params) : c_int;

extern proc magma_kprimme_normal(evals : c_ptr(_primme_complex_half), evecs : c_ptr(_primme_complex_half), resNorms : c_ptr(_primme_half), primme : c_ptr(primme_params)) : c_int;

extern proc magma_cprimme_normal(ref evals : _Complex float, ref evecs : _Complex float, ref resNorms : c_float, ref primme : primme_params) : c_int;

extern proc magma_cprimme_normal(evals : c_ptr(_Complex float), evecs : c_ptr(_Complex float), resNorms : c_ptr(c_float), primme : c_ptr(primme_params)) : c_int;

extern proc magma_zprimme_normal(ref evals : _Complex double, ref evecs : _Complex double, ref resNorms : c_double, ref primme : primme_params) : c_int;

extern proc magma_zprimme_normal(evals : c_ptr(_Complex double), evecs : c_ptr(_Complex double), resNorms : c_ptr(c_double), primme : c_ptr(primme_params)) : c_int;

extern proc kcprimme_normal(ref evals : _Complex float, ref evecs : _primme_complex_half, ref resNorms : c_float, ref primme : primme_params) : c_int;

extern proc kcprimme_normal(evals : c_ptr(_Complex float), evecs : c_ptr(_primme_complex_half), resNorms : c_ptr(c_float), primme : c_ptr(primme_params)) : c_int;

extern proc magma_kcprimme_normal(ref evals : _Complex float, ref evecs : _primme_complex_half, ref resNorms : c_float, ref primme : primme_params) : c_int;

extern proc magma_kcprimme_normal(evals : c_ptr(_Complex float), evecs : c_ptr(_primme_complex_half), resNorms : c_ptr(c_float), primme : c_ptr(primme_params)) : c_int;

extern proc primme_params_create() : c_ptr(primme_params);

extern proc primme_params_destroy(ref primme : primme_params) : c_int;

extern proc primme_params_destroy(primme : c_ptr(primme_params)) : c_int;

extern proc primme_initialize(ref primme : primme_params) : void;

extern proc primme_initialize(primme : c_ptr(primme_params)) : void;

extern proc primme_set_method(method : primme_preset_method, ref params : primme_params) : c_int;

extern proc primme_set_method(method : primme_preset_method, params : c_ptr(primme_params)) : c_int;

extern proc primme_display_params(in primme : primme_params) : void;

extern proc primme_free(ref primme : primme_params) : void;

extern proc primme_free(primme : c_ptr(primme_params)) : void;

extern proc primme_get_member(ref primme : primme_params, label_arg : primme_params_label, value : c_void_ptr) : c_int;

extern proc primme_get_member(primme : c_ptr(primme_params), label_arg : primme_params_label, value : c_void_ptr) : c_int;

extern proc primme_set_member(ref primme : primme_params, label_arg : primme_params_label, value : c_void_ptr) : c_int;

extern proc primme_set_member(primme : c_ptr(primme_params), label_arg : primme_params_label, value : c_void_ptr) : c_int;

extern proc primme_member_info(ref label_arg : primme_params_label, ref label_name : c_string, ref type_arg : primme_type, ref arity : c_int) : c_int;

extern proc primme_member_info(label_arg : c_ptr(primme_params_label), label_name : c_ptr(c_string), type_arg : c_ptr(primme_type), arity : c_ptr(c_int)) : c_int;

extern proc primme_constant_info(label_name : c_string, ref value : c_int) : c_int;

extern proc primme_constant_info(label_name : c_string, value : c_ptr(c_int)) : c_int;

extern proc hprimme_svds(ref svals : _primme_half, ref svecs : _primme_half, ref resNorms : _primme_half, ref primme_svds : primme_svds_params) : c_int;

extern proc hprimme_svds(svals : c_ptr(_primme_half), svecs : c_ptr(_primme_half), resNorms : c_ptr(_primme_half), primme_svds : c_ptr(primme_svds_params)) : c_int;

extern proc kprimme_svds(ref svals : _primme_half, ref svecs : _primme_complex_half, ref resNorms : _primme_half, ref primme_svds : primme_svds_params) : c_int;

extern proc kprimme_svds(svals : c_ptr(_primme_half), svecs : c_ptr(_primme_complex_half), resNorms : c_ptr(_primme_half), primme_svds : c_ptr(primme_svds_params)) : c_int;

extern proc sprimme_svds(ref svals : c_float, ref svecs : c_float, ref resNorms : c_float, ref primme_svds : primme_svds_params) : c_int;

extern proc sprimme_svds(svals : c_ptr(c_float), svecs : c_ptr(c_float), resNorms : c_ptr(c_float), primme_svds : c_ptr(primme_svds_params)) : c_int;

extern proc cprimme_svds(ref svals : c_float, ref svecs : _Complex float, ref resNorms : c_float, ref primme_svds : primme_svds_params) : c_int;

extern proc cprimme_svds(svals : c_ptr(c_float), svecs : c_ptr(_Complex float), resNorms : c_ptr(c_float), primme_svds : c_ptr(primme_svds_params)) : c_int;

extern proc dprimme_svds(ref svals : c_double, ref svecs : c_double, ref resNorms : c_double, ref primme_svds : primme_svds_params) : c_int;

extern proc dprimme_svds(svals : c_ptr(c_double), svecs : c_ptr(c_double), resNorms : c_ptr(c_double), primme_svds : c_ptr(primme_svds_params)) : c_int;

extern proc zprimme_svds(ref svals : c_double, ref svecs : _Complex double, ref resNorms : c_double, ref primme_svds : primme_svds_params) : c_int;

extern proc zprimme_svds(svals : c_ptr(c_double), svecs : c_ptr(_Complex double), resNorms : c_ptr(c_double), primme_svds : c_ptr(primme_svds_params)) : c_int;

extern proc magma_hprimme_svds(ref svals : _primme_half, ref svecs : _primme_half, ref resNorms : _primme_half, ref primme_svds : primme_svds_params) : c_int;

extern proc magma_hprimme_svds(svals : c_ptr(_primme_half), svecs : c_ptr(_primme_half), resNorms : c_ptr(_primme_half), primme_svds : c_ptr(primme_svds_params)) : c_int;

extern proc magma_kprimme_svds(ref svals : _primme_half, ref svecs : _primme_complex_half, ref resNorms : _primme_half, ref primme_svds : primme_svds_params) : c_int;

extern proc magma_kprimme_svds(svals : c_ptr(_primme_half), svecs : c_ptr(_primme_complex_half), resNorms : c_ptr(_primme_half), primme_svds : c_ptr(primme_svds_params)) : c_int;

extern proc magma_sprimme_svds(ref svals : c_float, ref svecs : c_float, ref resNorms : c_float, ref primme_svds : primme_svds_params) : c_int;

extern proc magma_sprimme_svds(svals : c_ptr(c_float), svecs : c_ptr(c_float), resNorms : c_ptr(c_float), primme_svds : c_ptr(primme_svds_params)) : c_int;

extern proc magma_cprimme_svds(ref svals : c_float, ref svecs : _Complex float, ref resNorms : c_float, ref primme_svds : primme_svds_params) : c_int;

extern proc magma_cprimme_svds(svals : c_ptr(c_float), svecs : c_ptr(_Complex float), resNorms : c_ptr(c_float), primme_svds : c_ptr(primme_svds_params)) : c_int;

extern proc magma_dprimme_svds(ref svals : c_double, ref svecs : c_double, ref resNorms : c_double, ref primme_svds : primme_svds_params) : c_int;

extern proc magma_dprimme_svds(svals : c_ptr(c_double), svecs : c_ptr(c_double), resNorms : c_ptr(c_double), primme_svds : c_ptr(primme_svds_params)) : c_int;

extern proc magma_zprimme_svds(ref svals : c_double, ref svecs : _Complex double, ref resNorms : c_double, ref primme_svds : primme_svds_params) : c_int;

extern proc magma_zprimme_svds(svals : c_ptr(c_double), svecs : c_ptr(_Complex double), resNorms : c_ptr(c_double), primme_svds : c_ptr(primme_svds_params)) : c_int;

extern proc hsprimme_svds(ref svals : c_float, ref svecs : _primme_half, ref resNorms : c_float, ref primme_svds : primme_svds_params) : c_int;

extern proc hsprimme_svds(svals : c_ptr(c_float), svecs : c_ptr(_primme_half), resNorms : c_ptr(c_float), primme_svds : c_ptr(primme_svds_params)) : c_int;

extern proc ksprimme_svds(ref svals : c_float, ref svecs : _primme_complex_half, ref resNorms : c_float, ref primme_svds : primme_svds_params) : c_int;

extern proc ksprimme_svds(svals : c_ptr(c_float), svecs : c_ptr(_primme_complex_half), resNorms : c_ptr(c_float), primme_svds : c_ptr(primme_svds_params)) : c_int;

extern proc magma_hsprimme_svds(ref svals : c_float, ref svecs : _primme_half, ref resNorms : c_float, ref primme_svds : primme_svds_params) : c_int;

extern proc magma_hsprimme_svds(svals : c_ptr(c_float), svecs : c_ptr(_primme_half), resNorms : c_ptr(c_float), primme_svds : c_ptr(primme_svds_params)) : c_int;

extern proc magma_ksprimme_svds(ref svals : c_float, ref svecs : _primme_complex_half, ref resNorms : c_float, ref primme_svds : primme_svds_params) : c_int;

extern proc magma_ksprimme_svds(svals : c_ptr(c_float), svecs : c_ptr(_primme_complex_half), resNorms : c_ptr(c_float), primme_svds : c_ptr(primme_svds_params)) : c_int;

extern proc primme_svds_params_create() : c_ptr(primme_svds_params);

extern proc primme_svds_params_destroy(ref primme_svds : primme_svds_params) : c_int;

extern proc primme_svds_params_destroy(primme_svds : c_ptr(primme_svds_params)) : c_int;

extern proc primme_svds_initialize(ref primme_svds : primme_svds_params) : void;

extern proc primme_svds_initialize(primme_svds : c_ptr(primme_svds_params)) : void;

extern proc primme_svds_set_method(method : primme_svds_preset_method, methodStage1 : primme_preset_method, methodStage2 : primme_preset_method, ref primme_svds : primme_svds_params) : c_int;

extern proc primme_svds_set_method(method : primme_svds_preset_method, methodStage1 : primme_preset_method, methodStage2 : primme_preset_method, primme_svds : c_ptr(primme_svds_params)) : c_int;

extern proc primme_svds_display_params(in primme_svds : primme_svds_params) : void;

extern proc primme_svds_free(ref primme_svds : primme_svds_params) : void;

extern proc primme_svds_free(primme_svds : c_ptr(primme_svds_params)) : void;

extern proc primme_svds_get_member(ref primme_svds : primme_svds_params, label_arg : primme_svds_params_label, value : c_void_ptr) : c_int;

extern proc primme_svds_get_member(primme_svds : c_ptr(primme_svds_params), label_arg : primme_svds_params_label, value : c_void_ptr) : c_int;

extern proc primme_svds_set_member(ref primme_svds : primme_svds_params, label_arg : primme_svds_params_label, value : c_void_ptr) : c_int;

extern proc primme_svds_set_member(primme_svds : c_ptr(primme_svds_params), label_arg : primme_svds_params_label, value : c_void_ptr) : c_int;

extern proc primme_svds_member_info(ref label_arg : primme_svds_params_label, ref label_name : c_string, ref type_arg : primme_type, ref arity : c_int) : c_int;

extern proc primme_svds_member_info(label_arg : c_ptr(primme_svds_params_label), label_name : c_ptr(c_string), type_arg : c_ptr(primme_type), arity : c_ptr(c_int)) : c_int;

extern proc primme_svds_constant_info(label_name : c_string, ref value : c_int) : c_int;

extern proc primme_svds_constant_info(label_name : c_string, value : c_ptr(c_int)) : c_int;

// ==== c2chapel typedefs ====

extern record JD_projectors {
  var LeftQ : c_int;
  var LeftX : c_int;
  var RightQ : c_int;
  var RightX : c_int;
  var SkewQ : c_int;
  var SkewX : c_int;
}

extern record correction_params {
  var precondition : c_int;
  var robustShifts : c_int;
  var maxInnerIterations : c_int;
  var projectors : JD_projectors;
  var convTest : primme_convergencetest;
  var relTolBase : c_double;
}

// primme_convergencetest enum
extern type primme_convergencetest = c_int;
extern const primme_full_LTolerance :primme_convergencetest;
extern const primme_decreasing_LTolerance :primme_convergencetest;
extern const primme_adaptive_ETolerance :primme_convergencetest;
extern const primme_adaptive :primme_convergencetest;


// primme_event enum
extern type primme_event = c_int;
extern const primme_event_outer_iteration :primme_event;
extern const primme_event_inner_iteration :primme_event;
extern const primme_event_restart :primme_event;
extern const primme_event_reset :primme_event;
extern const primme_event_converged :primme_event;
extern const primme_event_locked :primme_event;
extern const primme_event_message :primme_event;
extern const primme_event_profile :primme_event;


// primme_init enum
extern type primme_init = c_int;
extern const primme_init_default :primme_init;
extern const primme_init_krylov :primme_init;
extern const primme_init_random :primme_init;
extern const primme_init_user :primme_init;


// primme_op_datatype enum
extern type primme_op_datatype = c_int;
extern const primme_op_default :primme_op_datatype;
extern const primme_op_half :primme_op_datatype;
extern const primme_op_float :primme_op_datatype;
extern const primme_op_double :primme_op_datatype;
extern const primme_op_quad :primme_op_datatype;
extern const primme_op_int :primme_op_datatype;


// primme_orth enum
extern type primme_orth = c_int;
extern const primme_orth_default :primme_orth;
extern const primme_orth_implicit_I :primme_orth;
extern const primme_orth_explicit_I :primme_orth;


extern record primme_params {
  var n : int(64);
  var matrixMatvec : c_fn_ptr;
  var matrixMatvec_type : primme_op_datatype;
  var applyPreconditioner : c_fn_ptr;
  var applyPreconditioner_type : primme_op_datatype;
  var massMatrixMatvec : c_fn_ptr;
  var massMatrixMatvec_type : primme_op_datatype;
  var numProcs : c_int;
  var procID : c_int;
  var nLocal : int(64);
  var commInfo : c_void_ptr;
  var globalSumReal : c_fn_ptr;
  var globalSumReal_type : primme_op_datatype;
  var broadcastReal : c_fn_ptr;
  var broadcastReal_type : primme_op_datatype;
  var numEvals : c_int;
  var target : primme_target;
  var numTargetShifts : c_int;
  var targetShifts : c_ptr(c_double);
  var dynamicMethodSwitch : c_int;
  var locking : c_int;
  var initSize : c_int;
  var numOrthoConst : c_int;
  var maxBasisSize : c_int;
  var minRestartSize : c_int;
  var maxBlockSize : c_int;
  var maxMatvecs : int(64);
  var maxOuterIterations : int(64);
  var iseed : c_ptr(int(64));
  var aNorm : c_double;
  var BNorm : c_double;
  var invBNorm : c_double;
  var eps : c_double;
  var orth : primme_orth;
  var internalPrecision : primme_op_datatype;
  var printLevel : c_int;
  var outputFile : c_ptr(_file);
  var matrix : c_void_ptr;
  var preconditioner : c_void_ptr;
  var massMatrix : c_void_ptr;
  var ShiftsForPreconditioner : c_ptr(c_double);
  var initBasisMode : primme_init;
  var ldevecs : int(64);
  var ldOPs : int(64);
  var projectionParams : projection_params;
  var restartingParams : restarting_params;
  var correctionParams : correction_params;
  var stats : primme_stats;
  var convTestFun : c_fn_ptr;
  var convTestFun_type : primme_op_datatype;
  var convtest : c_void_ptr;
  var monitorFun : c_fn_ptr;
  var monitorFun_type : primme_op_datatype;
  var monitor : c_void_ptr;
  var queue : c_void_ptr;
  var profile : c_string;
}

// primme_params_label enum
extern type primme_params_label = c_int;
extern const PRIMME_n :primme_params_label;
extern const PRIMME_matrixMatvec :primme_params_label;
extern const PRIMME_matrixMatvec_type :primme_params_label;
extern const PRIMME_applyPreconditioner :primme_params_label;
extern const PRIMME_applyPreconditioner_type :primme_params_label;
extern const PRIMME_massMatrixMatvec :primme_params_label;
extern const PRIMME_massMatrixMatvec_type :primme_params_label;
extern const PRIMME_numProcs :primme_params_label;
extern const PRIMME_procID :primme_params_label;
extern const PRIMME_commInfo :primme_params_label;
extern const PRIMME_nLocal :primme_params_label;
extern const PRIMME_globalSumReal :primme_params_label;
extern const PRIMME_globalSumReal_type :primme_params_label;
extern const PRIMME_broadcastReal :primme_params_label;
extern const PRIMME_broadcastReal_type :primme_params_label;
extern const PRIMME_numEvals :primme_params_label;
extern const PRIMME_target :primme_params_label;
extern const PRIMME_numTargetShifts :primme_params_label;
extern const PRIMME_targetShifts :primme_params_label;
extern const PRIMME_locking :primme_params_label;
extern const PRIMME_initSize :primme_params_label;
extern const PRIMME_numOrthoConst :primme_params_label;
extern const PRIMME_maxBasisSize :primme_params_label;
extern const PRIMME_minRestartSize :primme_params_label;
extern const PRIMME_maxBlockSize :primme_params_label;
extern const PRIMME_maxMatvecs :primme_params_label;
extern const PRIMME_maxOuterIterations :primme_params_label;
extern const PRIMME_iseed :primme_params_label;
extern const PRIMME_aNorm :primme_params_label;
extern const PRIMME_BNorm :primme_params_label;
extern const PRIMME_invBNorm :primme_params_label;
extern const PRIMME_eps :primme_params_label;
extern const PRIMME_orth :primme_params_label;
extern const PRIMME_internalPrecision :primme_params_label;
extern const PRIMME_printLevel :primme_params_label;
extern const PRIMME_outputFile :primme_params_label;
extern const PRIMME_matrix :primme_params_label;
extern const PRIMME_massMatrix :primme_params_label;
extern const PRIMME_preconditioner :primme_params_label;
extern const PRIMME_ShiftsForPreconditioner :primme_params_label;
extern const PRIMME_initBasisMode :primme_params_label;
extern const PRIMME_projectionParams_projection :primme_params_label;
extern const PRIMME_restartingParams_maxPrevRetain :primme_params_label;
extern const PRIMME_correctionParams_precondition :primme_params_label;
extern const PRIMME_correctionParams_robustShifts :primme_params_label;
extern const PRIMME_correctionParams_maxInnerIterations :primme_params_label;
extern const PRIMME_correctionParams_projectors_LeftQ :primme_params_label;
extern const PRIMME_correctionParams_projectors_LeftX :primme_params_label;
extern const PRIMME_correctionParams_projectors_RightQ :primme_params_label;
extern const PRIMME_correctionParams_projectors_RightX :primme_params_label;
extern const PRIMME_correctionParams_projectors_SkewQ :primme_params_label;
extern const PRIMME_correctionParams_projectors_SkewX :primme_params_label;
extern const PRIMME_correctionParams_convTest :primme_params_label;
extern const PRIMME_correctionParams_relTolBase :primme_params_label;
extern const PRIMME_stats_numOuterIterations :primme_params_label;
extern const PRIMME_stats_numRestarts :primme_params_label;
extern const PRIMME_stats_numMatvecs :primme_params_label;
extern const PRIMME_stats_numPreconds :primme_params_label;
extern const PRIMME_stats_numGlobalSum :primme_params_label;
extern const PRIMME_stats_volumeGlobalSum :primme_params_label;
extern const PRIMME_stats_numBroadcast :primme_params_label;
extern const PRIMME_stats_volumeBroadcast :primme_params_label;
extern const PRIMME_stats_flopsDense :primme_params_label;
extern const PRIMME_stats_numOrthoInnerProds :primme_params_label;
extern const PRIMME_stats_elapsedTime :primme_params_label;
extern const PRIMME_stats_timeMatvec :primme_params_label;
extern const PRIMME_stats_timePrecond :primme_params_label;
extern const PRIMME_stats_timeOrtho :primme_params_label;
extern const PRIMME_stats_timeGlobalSum :primme_params_label;
extern const PRIMME_stats_timeBroadcast :primme_params_label;
extern const PRIMME_stats_timeDense :primme_params_label;
extern const PRIMME_stats_estimateMinEVal :primme_params_label;
extern const PRIMME_stats_estimateMaxEVal :primme_params_label;
extern const PRIMME_stats_estimateLargestSVal :primme_params_label;
extern const PRIMME_stats_estimateBNorm :primme_params_label;
extern const PRIMME_stats_estimateInvBNorm :primme_params_label;
extern const PRIMME_stats_maxConvTol :primme_params_label;
extern const PRIMME_stats_lockingIssue :primme_params_label;
extern const PRIMME_dynamicMethodSwitch :primme_params_label;
extern const PRIMME_convTestFun :primme_params_label;
extern const PRIMME_convTestFun_type :primme_params_label;
extern const PRIMME_convtest :primme_params_label;
extern const PRIMME_ldevecs :primme_params_label;
extern const PRIMME_ldOPs :primme_params_label;
extern const PRIMME_monitorFun :primme_params_label;
extern const PRIMME_monitorFun_type :primme_params_label;
extern const PRIMME_monitor :primme_params_label;
extern const PRIMME_queue :primme_params_label;
extern const PRIMME_profile :primme_params_label;


// primme_preset_method enum
extern type primme_preset_method = c_int;
extern const PRIMME_DEFAULT_METHOD :primme_preset_method;
extern const PRIMME_DYNAMIC :primme_preset_method;
extern const PRIMME_DEFAULT_MIN_TIME :primme_preset_method;
extern const PRIMME_DEFAULT_MIN_MATVECS :primme_preset_method;
extern const PRIMME_Arnoldi :primme_preset_method;
extern const PRIMME_GD :primme_preset_method;
extern const PRIMME_GD_plusK :primme_preset_method;
extern const PRIMME_GD_Olsen_plusK :primme_preset_method;
extern const PRIMME_JD_Olsen_plusK :primme_preset_method;
extern const PRIMME_RQI :primme_preset_method;
extern const PRIMME_JDQR :primme_preset_method;
extern const PRIMME_JDQMR :primme_preset_method;
extern const PRIMME_JDQMR_ETol :primme_preset_method;
extern const PRIMME_STEEPEST_DESCENT :primme_preset_method;
extern const PRIMME_LOBPCG_OrthoBasis :primme_preset_method;
extern const PRIMME_LOBPCG_OrthoBasis_Window :primme_preset_method;


// primme_projection enum
extern type primme_projection = c_int;
extern const primme_proj_default :primme_projection;
extern const primme_proj_RR :primme_projection;
extern const primme_proj_harmonic :primme_projection;
extern const primme_proj_refined :primme_projection;


extern record primme_stats {
  var numOuterIterations : int(64);
  var numRestarts : int(64);
  var numMatvecs : int(64);
  var numPreconds : int(64);
  var numGlobalSum : int(64);
  var numBroadcast : int(64);
  var volumeGlobalSum : int(64);
  var volumeBroadcast : int(64);
  var flopsDense : c_double;
  var numOrthoInnerProds : c_double;
  var elapsedTime : c_double;
  var timeMatvec : c_double;
  var timePrecond : c_double;
  var timeOrtho : c_double;
  var timeGlobalSum : c_double;
  var timeBroadcast : c_double;
  var timeDense : c_double;
  var estimateMinEVal : c_double;
  var estimateMaxEVal : c_double;
  var estimateLargestSVal : c_double;
  var estimateBNorm : c_double;
  var estimateInvBNorm : c_double;
  var maxConvTol : c_double;
  var estimateResidualError : c_double;
  var lockingIssue : int(64);
}

// primme_svds_operator enum
extern type primme_svds_operator = c_int;
extern const primme_svds_op_none :primme_svds_operator;
extern const primme_svds_op_AtA :primme_svds_operator;
extern const primme_svds_op_AAt :primme_svds_operator;
extern const primme_svds_op_augmented :primme_svds_operator;


extern record primme_svds_params {
  var primme : primme_params;
  var primmeStage2 : primme_params;
  var m : int(64);
  var n : int(64);
  var matrixMatvec : c_fn_ptr;
  var matrixMatvec_type : primme_op_datatype;
  var applyPreconditioner : c_fn_ptr;
  var applyPreconditioner_type : primme_op_datatype;
  var numProcs : c_int;
  var procID : c_int;
  var mLocal : int(64);
  var nLocal : int(64);
  var commInfo : c_void_ptr;
  var globalSumReal : c_fn_ptr;
  var globalSumReal_type : primme_op_datatype;
  var broadcastReal : c_fn_ptr;
  var broadcastReal_type : primme_op_datatype;
  var numSvals : c_int;
  var target : primme_svds_target;
  var numTargetShifts : c_int;
  var targetShifts : c_ptr(c_double);
  var method : primme_svds_operator;
  var methodStage2 : primme_svds_operator;
  var matrix : c_void_ptr;
  var preconditioner : c_void_ptr;
  var locking : c_int;
  var numOrthoConst : c_int;
  var aNorm : c_double;
  var eps : c_double;
  var precondition : c_int;
  var initSize : c_int;
  var maxBasisSize : c_int;
  var maxBlockSize : c_int;
  var maxMatvecs : int(64);
  var iseed : c_ptr(int(64));
  var printLevel : c_int;
  var internalPrecision : primme_op_datatype;
  var outputFile : c_ptr(_file);
  var stats : primme_svds_stats;
  var convTestFun : c_fn_ptr;
  var convTestFun_type : primme_op_datatype;
  var convtest : c_void_ptr;
  var monitorFun : c_fn_ptr;
  var monitorFun_type : primme_op_datatype;
  var monitor : c_void_ptr;
  var queue : c_void_ptr;
  var profile : c_string;
}

// primme_svds_params_label enum
extern type primme_svds_params_label = c_int;
extern const PRIMME_SVDS_primme :primme_svds_params_label;
extern const PRIMME_SVDS_primmeStage2 :primme_svds_params_label;
extern const PRIMME_SVDS_m :primme_svds_params_label;
extern const PRIMME_SVDS_n :primme_svds_params_label;
extern const PRIMME_SVDS_matrixMatvec :primme_svds_params_label;
extern const PRIMME_SVDS_matrixMatvec_type :primme_svds_params_label;
extern const PRIMME_SVDS_applyPreconditioner :primme_svds_params_label;
extern const PRIMME_SVDS_applyPreconditioner_type :primme_svds_params_label;
extern const PRIMME_SVDS_numProcs :primme_svds_params_label;
extern const PRIMME_SVDS_procID :primme_svds_params_label;
extern const PRIMME_SVDS_mLocal :primme_svds_params_label;
extern const PRIMME_SVDS_nLocal :primme_svds_params_label;
extern const PRIMME_SVDS_commInfo :primme_svds_params_label;
extern const PRIMME_SVDS_globalSumReal :primme_svds_params_label;
extern const PRIMME_SVDS_globalSumReal_type :primme_svds_params_label;
extern const PRIMME_SVDS_broadcastReal :primme_svds_params_label;
extern const PRIMME_SVDS_broadcastReal_type :primme_svds_params_label;
extern const PRIMME_SVDS_numSvals :primme_svds_params_label;
extern const PRIMME_SVDS_target :primme_svds_params_label;
extern const PRIMME_SVDS_numTargetShifts :primme_svds_params_label;
extern const PRIMME_SVDS_targetShifts :primme_svds_params_label;
extern const PRIMME_SVDS_method :primme_svds_params_label;
extern const PRIMME_SVDS_methodStage2 :primme_svds_params_label;
extern const PRIMME_SVDS_matrix :primme_svds_params_label;
extern const PRIMME_SVDS_preconditioner :primme_svds_params_label;
extern const PRIMME_SVDS_locking :primme_svds_params_label;
extern const PRIMME_SVDS_numOrthoConst :primme_svds_params_label;
extern const PRIMME_SVDS_aNorm :primme_svds_params_label;
extern const PRIMME_SVDS_eps :primme_svds_params_label;
extern const PRIMME_SVDS_precondition :primme_svds_params_label;
extern const PRIMME_SVDS_initSize :primme_svds_params_label;
extern const PRIMME_SVDS_maxBasisSize :primme_svds_params_label;
extern const PRIMME_SVDS_maxBlockSize :primme_svds_params_label;
extern const PRIMME_SVDS_maxMatvecs :primme_svds_params_label;
extern const PRIMME_SVDS_iseed :primme_svds_params_label;
extern const PRIMME_SVDS_printLevel :primme_svds_params_label;
extern const PRIMME_SVDS_internalPrecision :primme_svds_params_label;
extern const PRIMME_SVDS_outputFile :primme_svds_params_label;
extern const PRIMME_SVDS_stats_numOuterIterations :primme_svds_params_label;
extern const PRIMME_SVDS_stats_numRestarts :primme_svds_params_label;
extern const PRIMME_SVDS_stats_numMatvecs :primme_svds_params_label;
extern const PRIMME_SVDS_stats_numPreconds :primme_svds_params_label;
extern const PRIMME_SVDS_stats_numGlobalSum :primme_svds_params_label;
extern const PRIMME_SVDS_stats_volumeGlobalSum :primme_svds_params_label;
extern const PRIMME_SVDS_stats_numBroadcast :primme_svds_params_label;
extern const PRIMME_SVDS_stats_volumeBroadcast :primme_svds_params_label;
extern const PRIMME_SVDS_stats_numOrthoInnerProds :primme_svds_params_label;
extern const PRIMME_SVDS_stats_elapsedTime :primme_svds_params_label;
extern const PRIMME_SVDS_stats_timeMatvec :primme_svds_params_label;
extern const PRIMME_SVDS_stats_timePrecond :primme_svds_params_label;
extern const PRIMME_SVDS_stats_timeOrtho :primme_svds_params_label;
extern const PRIMME_SVDS_stats_timeGlobalSum :primme_svds_params_label;
extern const PRIMME_SVDS_stats_timeBroadcast :primme_svds_params_label;
extern const PRIMME_SVDS_stats_lockingIssue :primme_svds_params_label;
extern const PRIMME_SVDS_convTestFun :primme_svds_params_label;
extern const PRIMME_SVDS_convTestFun_type :primme_svds_params_label;
extern const PRIMME_SVDS_convtest :primme_svds_params_label;
extern const PRIMME_SVDS_monitorFun :primme_svds_params_label;
extern const PRIMME_SVDS_monitorFun_type :primme_svds_params_label;
extern const PRIMME_SVDS_monitor :primme_svds_params_label;
extern const PRIMME_SVDS_queue :primme_svds_params_label;
extern const PRIMME_SVDS_profile :primme_svds_params_label;


// primme_svds_preset_method enum
extern type primme_svds_preset_method = c_int;
extern const primme_svds_default :primme_svds_preset_method;
extern const primme_svds_hybrid :primme_svds_preset_method;
extern const primme_svds_normalequations :primme_svds_preset_method;
extern const primme_svds_augmented :primme_svds_preset_method;


extern record primme_svds_stats {
  var numOuterIterations : int(64);
  var numRestarts : int(64);
  var numMatvecs : int(64);
  var numPreconds : int(64);
  var numGlobalSum : int(64);
  var numBroadcast : int(64);
  var volumeGlobalSum : int(64);
  var volumeBroadcast : int(64);
  var numOrthoInnerProds : c_double;
  var elapsedTime : c_double;
  var timeMatvec : c_double;
  var timePrecond : c_double;
  var timeOrtho : c_double;
  var timeGlobalSum : c_double;
  var timeBroadcast : c_double;
  var lockingIssue : int(64);
}

// primme_svds_target enum
extern type primme_svds_target = c_int;
extern const primme_svds_largest :primme_svds_target;
extern const primme_svds_smallest :primme_svds_target;
extern const primme_svds_closest_abs :primme_svds_target;


// primme_target enum
extern type primme_target = c_int;
extern const primme_smallest :primme_target;
extern const primme_largest :primme_target;
extern const primme_closest_geq :primme_target;
extern const primme_closest_leq :primme_target;
extern const primme_closest_abs :primme_target;
extern const primme_largest_abs :primme_target;


// primme_type enum
extern type primme_type = c_int;
extern const primme_int :primme_type;
extern const primme_double :primme_type;
extern const primme_pointer :primme_type;
extern const primme_string :primme_type;


extern record projection_params {
  var projection : primme_projection;
}

extern record restarting_params {
  var maxPrevRetain : c_int;
}

// c2chapel thinks these typedefs are from the fake headers:
/*
extern type FILE = c_int;

// Opaque struct?
extern record MirBlob {};

// Opaque struct?
extern record MirBufferStream {};

// Opaque struct?
extern record MirConnection {};

// Opaque struct?
extern record MirDisplayConfig {};

extern type MirEGLNativeDisplayType = c_void_ptr;

extern type MirEGLNativeWindowType = c_void_ptr;

// Opaque struct?
extern record MirPersistentId {};

// Opaque struct?
extern record MirPromptSession {};

// Opaque struct?
extern record MirScreencast {};

// Opaque struct?
extern record MirSurface {};

// Opaque struct?
extern record MirSurfaceSpec {};

extern type _LOCK_RECURSIVE_T = c_int;

extern type _LOCK_T = c_int;

extern type __FILE = c_int;

extern type __ULong = c_int;

extern type __builtin_va_list = c_int;

extern type __dev_t = c_int;

extern type __gid_t = c_int;

extern type __gnuc_va_list = c_int;

extern type __int16_t = c_int;

extern type __int32_t = c_int;

extern type __int64_t = c_int;

extern type __int8_t = c_int;

extern type __int_least16_t = c_int;

extern type __int_least32_t = c_int;

extern type __loff_t = c_int;

extern type __off_t = c_int;

extern type __pid_t = c_int;

extern type __s16 = c_int;

extern type __s32 = c_int;

extern type __s64 = c_int;

extern type __s8 = c_int;

extern type __sigset_t = c_int;

extern type __tzinfo_type = c_int;

extern type __tzrule_type = c_int;

extern type __u16 = c_int;

extern type __u32 = c_int;

extern type __u64 = c_int;

extern type __u8 = c_int;

extern type __uid_t = c_int;

extern type __uint16_t = c_int;

extern type __uint32_t = c_int;

extern type __uint64_t = c_int;

extern type __uint8_t = c_int;

extern type __uint_least16_t = c_int;

extern type __uint_least32_t = c_int;

extern type _flock_t = c_int;

extern type _fpos_t = c_int;

extern type _iconv_t = c_int;

extern type _mbstate_t = c_int;

extern type _off64_t = c_int;

extern type _off_t = c_int;

extern type _sig_func_ptr = c_int;

extern type _ssize_t = c_int;

extern type _types_fd_set = c_int;

extern type bool = _Bool;

extern type caddr_t = c_int;

extern type clock_t = c_int;

extern type clockid_t = c_int;

extern type cookie_close_function_t = c_int;

extern type cookie_io_functions_t = c_int;

extern type cookie_read_function_t = c_int;

extern type cookie_seek_function_t = c_int;

extern type cookie_write_function_t = c_int;

extern type daddr_t = c_int;

extern type dev_t = c_int;

extern type div_t = c_int;

extern type fd_mask = c_int;

extern type fpos_t = c_int;

extern type gid_t = c_int;

extern type ino_t = c_int;

extern type int16_t = c_int;

extern type int32_t = c_int;

extern type int64_t = c_int;

extern type int8_t = c_int;

extern type int_fast16_t = c_int;

extern type int_fast32_t = c_int;

extern type int_fast64_t = c_int;

extern type int_fast8_t = c_int;

extern type int_least16_t = c_int;

extern type int_least32_t = c_int;

extern type int_least64_t = c_int;

extern type int_least8_t = c_int;

extern type intmax_t = c_int;

extern type intptr_t = c_int;

extern type jmp_buf = c_int;

extern type key_t = c_int;

extern type ldiv_t = c_int;

extern type lldiv_t = c_int;

extern type mbstate_t = c_int;

extern type mode_t = c_int;

extern type nlink_t = c_int;

extern type off_t = c_int;

extern type pid_t = c_int;

extern type pthread_attr_t = c_int;

extern type pthread_barrier_t = c_int;

extern type pthread_barrierattr_t = c_int;

extern type pthread_cond_t = c_int;

extern type pthread_condattr_t = c_int;

extern type pthread_key_t = c_int;

extern type pthread_mutex_t = c_int;

extern type pthread_mutexattr_t = c_int;

extern type pthread_once_t = c_int;

extern type pthread_rwlock_t = c_int;

extern type pthread_rwlockattr_t = c_int;

extern type pthread_spinlock_t = c_int;

extern type pthread_t = c_int;

extern type ptrdiff_t = c_int;

extern type rlim_t = c_int;

extern type sa_family_t = c_int;

extern type sem_t = c_int;

extern type sig_atomic_t = c_int;

extern type siginfo_t = c_int;

extern type sigjmp_buf = c_int;

extern type sigset_t = c_int;

extern type size_t = c_int;

extern type ssize_t = c_int;

extern type stack_t = c_int;

extern type suseconds_t = c_int;

extern type time_t = c_int;

extern type timer_t = c_int;

extern type u_char = c_int;

extern type u_int = c_int;

extern type u_long = c_int;

extern type u_short = c_int;

extern type uid_t = c_int;

extern type uint = c_int;

extern type uint16_t = c_int;

extern type uint32_t = c_int;

extern type uint64_t = c_int;

extern type uint8_t = c_int;

extern type uint_fast16_t = c_int;

extern type uint_fast32_t = c_int;

extern type uint_fast64_t = c_int;

extern type uint_fast8_t = c_int;

extern type uint_least16_t = c_int;

extern type uint_least32_t = c_int;

extern type uint_least64_t = c_int;

extern type uint_least8_t = c_int;

extern type uintmax_t = c_int;

extern type uintptr_t = c_int;

extern type useconds_t = c_int;

extern type ushort = c_int;

extern type va_list = c_int;

extern type wchar_t = c_int;

extern type wint_t = c_int;

// Opaque struct?
extern record xcb_connection_t {};

extern type xcb_visualid_t = uint(32);

extern type xcb_window_t = uint(32);

extern type z_stream = c_int;

*/
