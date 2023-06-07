module FFI {
  use CTypes;
  use IO;

  pragma "no doc"
  pragma "fn synchronization free"
  private extern proc c_pointer_return(const ref x : ?t) : c_ptr(t);

  inline proc c_const_ptrTo(const ref arr: []) {
    if (!arr.isRectangular() || !arr.domain.dist._value.dsiIsLayout()) then
      compilerError("Only single-locale rectangular arrays support c_ptrTo() at present");

    if (arr._value.locale != here) then
      halt("c_ptrTo() can only be applied to an array from the locale on which it lives (array is on locale "
           + arr._value.locale.id:string + ", call was made on locale " + here.id:string + ")");
    return c_pointer_return(arr[arr.domain.low]);
  }
  inline proc c_const_ptrTo(const ref x) {
    return c_pointer_return(x);
  }

  inline proc GET(addr, node, rAddr, size) {
    __primitive("chpl_comm_get", addr, node, rAddr, size);
  }

  inline proc PUT(addr, node, rAddr, size) {
    __primitive("chpl_comm_put", addr, node, rAddr, size);
  }

  proc unsafeViewAsExternalArray(const ref arr: []): chpl_external_array {
    if !isIntegralType(arr.domain.idxType) {
      // Probably not reachable any more, but may become reachable again
      // once support for interoperability with array types expands.
      compilerError("cannot return an array with indices that are not " +
                    "integrals");
    }
    if arr.domain.stridable {
      compilerError("cannot return a strided array");
    }
    if arr.domain.rank != 1 {
      compilerError("cannot return an array with rank != 1");
    }

    var externalArr = chpl_make_external_array_ptr(
      c_const_ptrTo(arr[arr.domain.low]), arr.size: uint);
    return externalArr;
  }

  inline proc _makeInds(shape: int ...?n) {
    var inds : n * range;
    foreach i in 0 ..# n {
      inds[i] = 0 ..# shape[i];
    }
    return inds;
  }

  pragma "no copy return"
  proc makeArrayFromPtr(ptr : c_ptr, shape)
      where isTuple(shape) && isHomogeneousTuple(shape) && shape[0].type == int {
    var dom = defaultDist.dsiNewRectangularDom(rank=shape.size,
                                               idxType=shape[0].type,
                                               stridable=false,
                                               inds=_makeInds((...shape)));
    dom._free_when_no_arrs = true;
    var arr = new unmanaged DefaultRectangularArr(eltType=ptr.eltType,
                                                  rank=dom.rank,
                                                  idxType=dom.idxType,
                                                  stridable=dom.stridable,
                                                  dom=dom,
                                                  data=ptr:_ddata(ptr.eltType),
                                                  externFreeFunc=nil,
                                                  externArr=true,
                                                  _borrowed=true);
    dom.add_arr(arr, locking = false);
    return _newArray(arr);
  }

  proc logDebug(msg...) {
    try! stderr.writeln("[Debug]   [", here, "]   ", (...msg));
  }


  require "lattice_symmetries_haskell.h";

  extern type ls_hs_particle_type = c_int;
  extern const LS_HS_SPIN : ls_hs_particle_type;
  extern const LS_HS_SPINFUL_FERMION : ls_hs_particle_type;
  extern const LS_HS_SPINLESS_FERMION : ls_hs_particle_type;

  extern record ls_hs_basis_kernels {
    var state_index_kernel : c_fn_ptr;
    var state_index_data : c_void_ptr;
  }
  extern record ls_hs_basis {
    var number_sites : c_int;
    var number_particles : c_int;
    var number_up : c_int;
    var particle_type : ls_hs_particle_type;
    var spin_inversion : c_int;
    var state_index_is_identity : bool;
    var requires_projection : bool;
    var kernels : c_ptr(ls_hs_basis_kernels);
    var representatives : chpl_external_array;
    // ... other stuff ...
  }
  extern record ls_hs_expr {
    // ... other stuff ...
  }
  extern record ls_hs_nonbranching_terms {
    var number_terms : c_int;
    var number_bits : c_int;
    // ... other stuff ...
  }
  extern record ls_hs_operator {
    var basis : c_ptr(ls_hs_basis);
    var off_diag_terms : c_ptr(ls_hs_nonbranching_terms);
    var diag_terms : c_ptr(ls_hs_nonbranching_terms);
    // ... other stuff ...
  }

  extern record ls_hs_yaml_config {
    var basis : c_ptr(ls_hs_basis);
    var hamiltonian : c_ptr(ls_hs_operator);
    var number_observables : c_int;
    var observables : c_ptr(c_ptr(ls_hs_operator));
  }

  extern proc ls_hs_init();
  extern proc ls_hs_exit();

  proc initRuntime() {
    coforall loc in Locales do on loc {
      ls_hs_init();
    }
  }

  proc deinitRuntime() {
    coforall loc in Locales do on loc {
      ls_hs_exit();
    }
  }

  // extern proc ls_hs_create_basis(particleType : ls_hs_particle_type, numberSites : c_int,
  //                                numberParticles : c_int, numberUp : c_int) : c_ptr(ls_hs_basis);
  extern proc ls_hs_clone_basis(basis : c_ptr(ls_hs_basis)) : c_ptr(ls_hs_basis);
  extern proc ls_hs_destroy_basis(basis : c_ptr(ls_hs_basis));
  extern proc ls_hs_min_state_estimate(basis : c_ptr(ls_hs_basis)) : uint(64);
  extern proc ls_hs_max_state_estimate(basis : c_ptr(ls_hs_basis)) : uint(64);
  extern proc ls_hs_basis_number_bits(basis : c_ptr(ls_hs_basis)) : c_int;
  extern proc ls_hs_basis_number_words(basis : c_ptr(ls_hs_basis)) : c_int;
  extern proc ls_hs_basis_has_fixed_hamming_weight(basis : c_ptr(ls_hs_basis)) : bool;
  extern proc ls_hs_basis_has_spin_inversion_symmetry(basis : c_ptr(ls_hs_basis)) : bool;
  extern proc ls_hs_basis_has_permutation_symmetries(basis : c_ptr(ls_hs_basis)) : bool;
  extern proc ls_hs_basis_requires_projection(basis : c_ptr(ls_hs_basis)) : bool;

  extern proc ls_hs_basis_from_json(json_string : c_string) : c_ptr(ls_hs_basis);
  extern proc ls_hs_basis_to_json(basis : c_ptr(ls_hs_basis)) : c_string;
  extern proc ls_hs_destroy_string(str : c_string);

  // extern proc ls_hs_create_spin_basis_from_json(json_string : c_string) : c_ptr(ls_hs_basis);
  // extern proc ls_hs_create_spin_basis_from_yaml(yaml_filename : c_string) : c_ptr(ls_hs_basis);
  // extern proc ls_hs_create_spinful_fermion_basis_from_json(json_string : c_string) : c_ptr(ls_hs_basis);
  // extern proc ls_hs_create_spinless_fermion_basis_from_json(json_string : c_string) : c_ptr(ls_hs_basis);

  extern proc ls_hs_fixed_hamming_state_to_index(basis_state : uint(64)) : c_ptrdiff;
  extern proc ls_hs_fixed_hamming_index_to_state(state_index : c_ptrdiff, hamming_weight : c_int) : uint(64);

  extern proc ls_hs_basis_build(basis : c_ptr(ls_hs_basis));

  extern proc ls_hs_unchecked_set_representatives(basis : c_ptr(ls_hs_basis),
                                                  states : c_ptr(chpl_external_array));

  extern proc ls_hs_state_index(basis : c_ptr(ls_hs_basis), batch_size : c_ptrdiff,
    spins : c_ptr(uint(64)), spins_stride : c_ptrdiff, indices : c_ptr(c_ptrdiff),
    indices_stride : c_ptrdiff);

  extern proc ls_hs_is_representative(basis : c_ptr(ls_hs_basis), batch_size : c_ptrdiff,
    alphas : c_ptr(uint(64)), alphas_stride : c_ptrdiff, are_representatives : c_ptr(uint(8)),
    norms : c_ptr(real(64)));

  extern proc ls_hs_state_info(basis : c_ptr(ls_hs_basis), batch_size : c_ptrdiff,
                               alphas : c_ptr(uint(64)), alphas_stride : c_ptrdiff,
                               betas : c_ptr(uint(64)), betas_stride : c_ptrdiff,
                               characters : c_ptr(complex(128)), norms : c_ptr(real(64)));

  // extern proc ls_hs_create_basis_kernels(basis : c_ptr(ls_hs_basis)) : c_ptr(ls_hs_basis_kernels);
  // extern proc ls_hs_destroy_basis_kernels(kernels : c_ptr(ls_hs_basis_kernels));
  extern proc ls_hs_create_expr(expr : c_string);
  extern proc ls_hs_destroy_expr(expr : c_ptr(ls_hs_expr));
  extern proc ls_hs_expr_to_json(expr : c_ptr(ls_hs_expr)) : c_string;
  extern proc ls_hs_expr_from_json(json_string : c_string) : c_ptr(ls_hs_expr);

  extern proc ls_hs_create_operator(basis : c_ptr(ls_hs_basis),
                                    expr : c_ptr(ls_hs_expr)) : c_ptr(ls_hs_operator);
  extern proc ls_hs_clone_operator(op : c_ptr(ls_hs_operator)) : c_ptr(ls_hs_operator);
  extern proc ls_hs_operator_plus(a : c_ptr(ls_hs_operator), b : c_ptr(ls_hs_operator)) : c_ptr(ls_hs_operator);
  extern proc ls_hs_print_terms(op : c_ptr(ls_hs_operator));
  extern proc ls_hs_destroy_operator(op : c_ptr(ls_hs_operator));

  extern proc ls_hs_operator_max_number_off_diag(op : c_ptr(ls_hs_operator)) : c_int;
  extern proc ls_hs_operator_is_hermitian(op : c_ptr(ls_hs_operator)) : bool;
  extern proc ls_hs_operator_is_real(op : c_ptr(ls_hs_operator)) : bool;

  extern proc ls_hs_operator_get_expr(op : c_ptr(ls_hs_operator)) : c_ptr(ls_hs_expr);

  extern proc ls_hs_load_hamiltonian_from_yaml(filename : c_string) : c_ptr(ls_hs_operator);

  extern proc ls_hs_load_yaml_config(filename : c_string) : c_ptr(ls_hs_yaml_config);
  extern proc ls_hs_destroy_yaml_config(p : c_ptr(ls_hs_yaml_config));

  extern proc ls_hs_operator_apply_diag_kernel(op : c_ptr(ls_hs_operator),
    batchSize : c_ptrdiff, alphas : c_ptr(uint(64)), alphas_stride : c_ptrdiff,
    coeffs : c_ptr(complex(128)));

  extern proc ls_hs_operator_apply_off_diag_kernel(op : c_ptr(ls_hs_operator),
    batchSize : c_ptrdiff, alphas : c_ptr(uint(64)), alphas_stride : c_ptrdiff,
    betas : c_ptr(uint(64)), betas_stride : c_ptrdiff, coeffs : c_ptr(complex(128)));

  extern proc ls_internal_operator_apply_diag_x1(
    op : c_ptr(ls_hs_operator), batch_size : c_ptrdiff, alphas : c_ptr(uint(64)),
    ys : c_ptr(real(64)), xs : c_ptr(real(64)));
  extern proc ls_internal_operator_apply_off_diag_x1(
    op : c_ptr(ls_hs_operator), batch_size : c_ptrdiff, alphas : c_ptr(uint(64)),
    betas : c_ptr(uint(64)), coeffs : c_ptr(complex(128)), offsets : c_ptr(c_ptrdiff),
    xs : c_ptr(real));

  extern proc ls_hs_evaluate_wavefunction_via_statevector(
    kernels : c_ptr(ls_hs_basis_kernels), batch_size : c_ptrdiff,
    alphas : c_ptr(uint(64)), alphas_stride : c_ptrdiff, state_vector : c_void_ptr,
    element_size : uint(64), coeffs : c_void_ptr);


  extern record ls_chpl_kernels {
    var enumerate_states : c_fn_ptr;
    var operator_apply_off_diag : c_fn_ptr;
    var operator_apply_diag : c_fn_ptr;
    var matrix_vector_product : c_fn_ptr;
  }
  extern proc ls_hs_internal_set_chpl_kernels(kernels : c_ptr(ls_chpl_kernels));

  // extern proc print_external_string(s : c_string);
  // extern proc ls_hs_hdf5_get_dataset_rank(path: c_string, dataset: c_string):c_uint;
  // extern proc ls_hs_hdf5_get_dataset_shape(path: c_string, dataset: c_string,
  //                                          shape: c_ptr(uint(64)));
  // extern proc ls_hs_hdf5_create_dataset_u64(path: c_string,
  //     dataset: c_string, dim: c_uint, shape: c_ptr(uint(64)));
  // extern proc ls_hs_hdf5_create_dataset_f64(path: c_string,
  //     dataset: c_string, dim: c_uint, shape: c_ptr(uint(64)));
  // extern proc ls_hs_hdf5_write_chunk_u64(path: c_string, dataset: c_string,
  //     dim: c_uint, offset: c_ptr(uint(64)), shape: c_ptr(uint(64)), data: c_ptr(uint(64)));
  // extern proc ls_hs_hdf5_write_chunk_f64(path: c_string, dataset: c_string,
  //     dim: c_uint, offset: c_ptr(uint(64)), shape: c_ptr(uint(64)), data: c_ptr(real(64)));
  // extern proc ls_hs_hdf5_read_chunk_u64(path: c_string, dataset: c_string,
  //     dim: c_uint, offset: c_ptr(uint(64)), shape: c_ptr(uint(64)), data: c_ptr(uint(64)));
  // extern proc ls_hs_hdf5_read_chunk_f64(path: c_string, dataset: c_string,
  //     dim: c_uint, offset: c_ptr(uint(64)), shape: c_ptr(uint(64)), data: c_ptr(real(64)));
}
