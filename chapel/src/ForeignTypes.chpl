module ForeignTypes {
  use FFI;

  use CTypes;
  use ByteBufferHelpers;
  use Time;

  record Basis {
    var payload : c_ptr(ls_hs_basis);
    var owning : bool;
    var _origin : locale;
    var _json_repr : string;
    var _hasPermutationSymmetries : bool;

    proc init(p : c_ptr(ls_hs_basis), owning : bool = true) {
      assert(here == p.locale);
      this.payload = p;
      this.owning = owning;
      this._origin = here;
      complete();
      this._json_repr = _toJSON();
      this._hasPermutationSymmetries =
        ls_hs_basis_has_permutation_symmetries(this.payload);
    }
    proc init(jsonString : string) {
      const s = jsonString.localize();
      this.payload = ls_hs_basis_from_json(s.c_str());
      this.owning = true;
      this._origin = here;
      this._json_repr = s;
      this._hasPermutationSymmetries =
        ls_hs_basis_has_permutation_symmetries(this.payload);
    }

    proc init=(const ref from : Basis) {
      if here == from._origin {
        this.payload = ls_hs_clone_basis(from.payload);
        this.owning = true;
        this._origin = here;
        this._json_repr = from.json_repr;
        this._hasPermutationSymmetries =
          ls_hs_basis_has_permutation_symmetries(this.payload);
      }
      else {
        const s = from.json_repr.localize();
        this.payload = ls_hs_basis_from_json(s.c_str());
        this.owning = true;
        this._origin = here;
        this._json_repr = s;
        this._hasPermutationSymmetries =
          ls_hs_basis_has_permutation_symmetries(this.payload);
      }
    }

    proc _toJSON() : string {
      const c_str = ls_hs_basis_to_json(payload);
      defer ls_hs_destroy_string(c_str);
      return c_str:string;
    }

    proc json_repr const ref : string { return _json_repr; }

    proc _destroy() {
      if owning then
        ls_hs_destroy_basis(payload);
    }

    proc deinit() {
      _destroy();
    }

    proc build() { ls_hs_basis_build(payload); }

    proc uncheckedSetRepresentatives(representatives : [] uint(64)) {
      var arr = unsafeViewAsExternalArray(representatives);
      ls_hs_unchecked_set_representatives(payload, c_ptrTo(arr));
    }

    proc isSpinBasis() { return payload.deref().particle_type == LS_HS_SPIN; }
    proc isSpinfulFermionicBasis() { return payload.deref().particle_type == LS_HS_SPINFUL_FERMION; }
    proc isSpinlessFermionicBasis() { return payload.deref().particle_type == LS_HS_SPINLESS_FERMION; }
    proc isStateIndexIdentity() { return payload.deref().state_index_is_identity; }
    proc requiresProjection() { return payload.deref().requires_projection; }
    proc isHammingWeightFixed() {
      // logDebug("ls_hs_basis_has_fixed_hamming_weight");
      return ls_hs_basis_has_fixed_hamming_weight(payload);
    }
    proc hasSpinInversionSymmetry() {
      return ls_hs_basis_has_spin_inversion_symmetry(payload);
    }
    inline proc hasPermutationSymmetries() {
      return _hasPermutationSymmetries;
    }

    proc numberBits : int { return ls_hs_basis_number_bits(payload):int; }
    proc numberWords : int { return ls_hs_basis_number_words(payload):int; }
    proc numberSites() : int { return payload.deref().number_sites; }
    proc numberParticles() : int { return payload.deref().number_particles; }
    proc numberUp() : int { return payload.deref().number_up; }
    proc spinInversion : int { return payload.deref().spin_inversion:int; }

    proc minStateEstimate() : uint(64) {
      // logDebug("ls_hs_min_state_estimate");
      return ls_hs_min_state_estimate(payload);
    }
    proc maxStateEstimate() : uint(64) {
      // logDebug("ls_hs_max_state_estimate");
      return ls_hs_max_state_estimate(payload);
    }

    proc representatives() {
      ref rs = payload.deref().representatives;
      if rs.elts == nil then
        halt("basis is not built");
      return makeArrayFromExternArray(rs, uint(64));
    }
  }

  operator Basis.= (ref lhs : Basis, const ref rhs : Basis) {
    assert(lhs.locale == rhs.locale);
    assert(lhs._origin == lhs.locale);
    lhs._destroy();
    lhs.payload = ls_hs_clone_basis(rhs.payload);
    lhs.owning = true;
    lhs._origin = rhs._origin;
    lhs._json_repr = rhs._json_repr;
  }

  proc isRepresentative(const ref basis : Basis, const ref alphas : [?D] uint(64),
                        ref are_representatives : [?D2] uint(8),
                        ref norms : [?D3] real(64))
      where D.rank == 1 && D2.rank == 1 && D3.rank == 1 {
    const batchSize = alphas.size;
    assert(are_representatives.size == batchSize);
    assert(norms.size == batchSize);

    ls_hs_is_representative(
      basis.payload,
      batchSize,
      c_const_ptrTo(alphas), 1,
      c_ptrTo(are_representatives),
      c_ptrTo(norms)
    );
  }
  proc isRepresentative(const ref basis : Basis, const ref alphas : [?D] uint(64))
      where D.rank == 1 {
    const batchSize = alphas.size;
    var areRepresentatives : [0 ..# batchSize] uint(8);
    var norms : [0 ..# batchSize] real(64);
    isRepresentative(basis, alphas, areRepresentatives, norms);
    return (areRepresentatives, norms);
  }

  record Operator {
    var payload : c_ptr(ls_hs_operator);
    var basis : Basis;
    var owning : bool;
    var _origin : locale;

    // proc init(const ref basis : Basis, expression : string, const ref indices : [?D] ?i)
    //     where D.rank == 2 {
    //   assert(basis._origin == here);

    //   const c_indices = indices:c_int;
    //   const c_expr = expression.localize().c_str();
    //   print_external_string(c_expr);
    //   this.payload = ls_hs_create_operator(basis.payload, c_expr,
    //     indices.dim(0).size:c_int, indices.dim(1).size:c_int, c_const_ptrTo(c_indices));
    //   this.basis = new Basis(this.payload.deref().basis, owning=false);
    // }
    proc init(const ref basis : Basis, expression : c_ptr(ls_hs_expr)) {
      assert(basis._origin == here);
      assert(expression.locale == here);

      this.payload = ls_hs_create_operator(basis.payload, expression);
      this.basis = new Basis(this.payload.deref().basis, owning=false);
      this.owning = true;
      this._origin = here;
    }
    proc init(raw : c_ptr(ls_hs_operator), owning : bool = true) {
      assert(raw.locale == here);
      // logDebug("Creating Operator from ", raw, ", raw.deref().basis=", raw.deref().basis);
      this.payload = raw;
      this.basis = new Basis(raw.deref().basis, owning=false);
      this.owning = owning;
      this._origin = here;
    }

    proc init=(const ref from : Operator) {
      if from._origin == here {
        init(ls_hs_clone_operator(from.payload), owning=true);
      }
      else {
        var s : string;
        on from._origin do
          s = from.exprToJSON();

        var expr = ls_hs_expr_from_json(s.localize().c_str());
        defer ls_hs_destroy_expr(expr);

        var basis = from.basis;
        init(basis, expr);
      }
    }

    proc exprToJSON() {
      var expr = ls_hs_operator_get_expr(payload);
      if expr == nil then halt("failed to get expr");
      defer ls_hs_destroy_expr(expr);

      var c_str = ls_hs_expr_to_json(expr);
      defer ls_hs_destroy_string(c_str);

      return c_str:string;
    }

    proc deinit() {
      if owning then
        ls_hs_destroy_operator(payload);
    }

    inline proc numberDiagTerms() : int {
      const p = payload.deref().diag_terms;
      if (p == nil) { return 0; }
      return p.deref().number_terms:int;
    }

    inline proc numberOffDiagTerms() : int {
      return ls_hs_operator_max_number_off_diag(payload):int;
      // const p = payload.deref().off_diag_terms;
      // if (p == nil) { return 0; }
      // return p.deref().number_terms:int;
    }

    proc writeTerms() {
      ls_hs_print_terms(payload);
    }

    proc _localApplyDiag(const ref alphas : [?D] uint(64), ref coeffs : [?D2] complex(128))
        where D.rank == 1 && D2.rank == 1 {
      assert(alphas.size <= coeffs.size);
      const batchSize = alphas.size;
      ls_hs_operator_apply_diag_kernel(payload, batchSize, c_const_ptrTo(alphas), 1, c_ptrTo(coeffs));
    }

    proc _localApplyOffDiag(const ref alphas : [?D] uint(64), ref betas : [?D2] uint(64),
                            ref coeffs : [?D3] complex(128))
        where D.rank == 1 && D2.rank == 1 && D3.rank == 1 {
      const batchSize = alphas.size;
      const numberTerms = numberOffDiagTerms();
      assert(betas.size >= batchSize * numberTerms);
      assert(coeffs.size >= batchSize * numberTerms);
      ls_hs_operator_apply_off_diag_kernel(payload, batchSize,
          c_const_ptrTo(alphas), 1, c_ptrTo(betas), 1, c_ptrTo(coeffs));
    }

    proc isHermitian : bool { return ls_hs_operator_is_hermitian(payload); }
    proc isReal : bool { return ls_hs_operator_is_real(payload); }
  }

  proc loadConfigFromYaml(filename : string, param hamiltonian = false,
                                             param observables = false) {
    const configPtr = ls_hs_load_yaml_config(filename.localize().c_str());
    if configPtr == nil then
      halt("failed to load Config from '" + filename + "'");
    defer ls_hs_destroy_yaml_config(configPtr);
  
    // logDebug("built config; creating basis...");
    ref conf = configPtr.deref();
    const basis = new Basis(ls_hs_clone_basis(conf.basis), owning=true);

    // logDebug("built basis; creating hamiltonian...");
    if hamiltonian && conf.hamiltonian == nil then
      halt("'" + filename + "' does not contain a Hamiltonian");
    const h = if hamiltonian
                then new Operator(ls_hs_clone_operator(conf.hamiltonian), owning=true)
                else nil;
    // logDebug("built hamiltonian; creating observables...");

    var os : [0 ..# conf.number_observables:int] Operator =
      [i in 0 ..# conf.number_observables:int]
        new Operator(ls_hs_clone_operator(conf.observables[i]), owning=true);

    if !hamiltonian && !observables then return basis;
    if hamiltonian && !observables then return (basis, h);
    if !hamiltonian && observables then return (basis, os);
    if hamiltonian && observables then return (basis, h, os);
  }

  // proc loadHamiltonianFromYaml(filename : string) {
  //   var ptr = ls_hs_load_hamiltonian_from_yaml(filename.localize().c_str());
  //   if ptr == nil then
  //     halt("failed to load Hamiltonian from " + filename);
  //   return new Operator(ptr);
  // }

  // operator +(const ref a : Operator, const ref b : Operator) {
  //   return new Operator(ls_hs_operator_plus(a.payload, b.payload));
  // }
}
