module ForeignTypes {
  use FFI;
  // use CommonParameters;
  import Communication;

  use CTypes;
  use ChplConfig;
  // use ByteBufferHelpers;
  use IO;
  use JSON;
  // use Time;

  require "ForeignTypes.h";
  extern proc ls_chpl_get_basis_info(p : c_ptrConst(ls_hs_basis)) : c_ptrConst(ls_hs_basis_info);
  extern proc ls_chpl_get_is_representative_kernel(p : c_ptrConst(ls_hs_basis)) : c_fn_ptr;
  extern proc ls_chpl_invoke_is_representative_kernel(kernel : c_fn_ptr, count : int(64),
                                                      basis_states : c_ptrConst(uint(64)), norms : c_ptr(real(64)));

  // pragma "fn synchronization free"
  // private extern proc c_pointer_return(const ref x : ?t) : c_ptr(t);

  // inline proc c_const_ptrTo(const ref arr: []) {
  //   if (!arr.isRectangular() || !arr.domain.dist._value.dsiIsLayout()) then
  //     compilerError("Only single-locale rectangular arrays support c_ptrTo() at present");

  //   if (arr._value.locale != here) then
  //     halt("c_ptrTo() can only be applied to an array from the locale on which it lives (array is on locale "
  //          + arr._value.locale.id:string + ", call was made on locale " + here.id:string + ")");
  //   return c_pointer_return(arr[arr.domain.low]);
  // }
  // inline proc c_const_ptrTo(const ref x) {
  //   return c_pointer_return(x);
  // }

  /*
  inline proc GET(addr, node, rAddr, size) {
    // __primitive("chpl_comm_get", addr, node, rAddr, size);
    Communication.get(dest = addr,
                      src = rAddr,
                      srcLocID = node,
                      numBytes = size);
  }

  inline proc PUT(addr, node, rAddr, size) {
    // __primitive("chpl_comm_put", addr, node, rAddr, size);
    Communication.put(dest = rAddr,
                      src = addr,
                      destLocID = node,
                      numBytes = size);
  }
  */

  /*
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
      c_ptrToConst(arr[arr.domain.low]), arr.size: uint);
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
  proc makeArrayFromPtr(ptr : c_ptr(?), shape)
      where isTuple(shape) && isHomogeneousTuple(shape) && shape[0].type == int {
    var dom = defaultDist.dsiNewRectangularDom(rank=shape.size,
                                               idxType=shape[0].type,
                                               strides=strideKind.one,
                                               inds=_makeInds((...shape)));
    dom._free_when_no_arrs = true;
    var arr = new unmanaged DefaultRectangularArr(eltType=ptr.eltType,
                                                  rank=dom.rank,
                                                  idxType=dom.idxType,
                                                  strides=strideKind.one,
                                                  dom=dom,
                                                  data=ptr:_ddata(ptr.eltType),
                                                  externFreeFunc=nil,
                                                  externArr=true,
                                                  _borrowed=true);
    dom.add_arr(arr, locking = false);
    return _newArray(arr);
  }
  */

  proc logDebug(msg...) {
    try! stderr.writeln("[Debug]   [", here, "]   ", (...msg));
  }

  /*
  class _NilableClass {}

  record optional : serializable {
    type eltType;
    var _value : eltType;
    var _hasValue : bool;

    proc init(type t)
        where isDefaultInitializable(t) {
      this.eltType = t;
      this._hasValue = false;
    }
    proc init(type t, x : nothing) do this.init(t);
    proc init(type t, x : nil) do this.init(t);

    proc init(x : ?t)
        where isCopyable(t) {
      this.eltType = t;
      this._value = x;
      this._hasValue = true;
    }

    proc init(const ref x : optional(?t)) {
      this.eltType = t;
      this._value = x._value;
      this._hasValue = x._hasValue;
    }

    proc init=(x : optional(?t)) do this.init(x);
    proc init=(x : ?t) do this.init(x);

    proc init(type eltType, reader : fileReader(?), ref deserializer : ?dt) throws {
      if reader.matchLiteral("null") then
        init(eltType);
      else
        init(deserializer.deserializeType(reader, eltType));
    }

    operator :(x : ?eltType, type t : optional(eltType)) : optional(eltType) { return new optional(x); }
    operator :(x : nothing, type t : optional(eltType)) : optional(eltType) { return new optional(eltType); }
    operator :(x : nil, type t : optional(eltType)) : optional(eltType) { return new optional(eltType); }

    proc value() throws {
      if !hasValue then
        throw new Error("trying to dereference a nil optional");
      return _value;
    }

    inline proc valueOr(def : eltType) { return if hasValue then _value else def; }

    inline proc hasValue { return _hasValue; }

    proc serialize(writer : fileWriter(?), ref serializer : ?st) throws {
      if hasValue then
        serializer.serializeValue(writer, _value);
      else {
        // jsonSerializer misbehaves when we pass nil by itself :/ So we
        // pretend as if that nil represents some nilable class.
        serializer.serializeValue(writer, nil:unmanaged _NilableClass?);
      }
    }
  }
  */

  record HsResult : initDeserializable {
    type eltType;
    var value : eltType;

    proc init(value : ?eltType) {
      this.eltType = eltType;
      this.value = value;
    }

    proc init(type eltType, reader : fileReader(?), ref deserializer : ?dt) throws {
      this.eltType = eltType;
      init this;

      var m = deserializer.startMap(reader);
      const key = m.readKey(string);
      select key {
        when "Left" do halt(m.readValue(string));
        when "Right" do m.readValue(this.value);
        otherwise do halt("unexpected key when parsing HsResult: " + key);
      }
      m.endMap();
    }
  }

  record BasisForeignInterface {
    proc clone(p : c_ptrConst(ls_hs_basis)) {
      assert(p != nil);
      const mp = p:c_ptr(ls_hs_basis);
      ls_hs_internal_object_inc_ref_count(c_ptrTo(mp.deref().base));
      return mp;
    }

    proc destroy(p : c_ptr(ls_hs_basis)) {
      assert(p != nil);
      ls_hs_destroy_basis(p);
    }

    proc to_json(p : c_ptrConst(ls_hs_basis)) {
      const c_str = ls_hs_basis_to_json(p);
      defer ls_hs_destroy_string(c_str);
      return try! string.createCopyingBuffer(c_str);
    }

    proc from_json(s : c_ptrConst(c_char)) throws {
      const c_str = ls_hs_basis_from_json(s);
      defer ls_hs_destroy_string(c_str);

      var file = openMemFile();
      var writer = file.writer();
      writer.write(string.createBorrowingBuffer(c_str));
      writer.close();
      var r = file.reader(deserializer=new jsonDeserializer());
      var p = r.read(HsResult(uint(64))).value;
      return p:c_ptr(void):c_ptr(ls_hs_basis);
    }
  }

  /*
  record OperatorForeignInterface {
    proc clone(p : c_ptrConst(ls_hs_operator)) {
      assert(p != nil);
      const mp = p:c_ptr(ls_hs_operator);
      ls_hs_internal_object_inc_ref_count(c_ptrTo(mp.deref().base));
      return mp;
    }

    proc destroy(p : c_ptr(ls_hs_basis)) {
      assert(p != nil);
      ls_hs_destroy_basis(p);
    }

    proc to_json(p : c_ptrConst(ls_hs_basis)) {
      const c_str = ls_hs_basis_to_json(p);
      defer ls_hs_destroy_string(c_str);
      return try! string.createCopyingBuffer(c_str);
    }

    proc from_json(s : c_ptrConst(c_char)) throws {
      const c_str = ls_hs_basis_from_json(s);
      defer ls_hs_destroy_string(c_str);

      var file = openMemFile();
      var writer = file.writer();
      writer.write(string.createBorrowingBuffer(c_str));
      writer.close();
      var r = file.reader(deserializer=new jsonDeserializer());
      var p = r.read(HsResult(uint(64))).value;
      return p:c_ptr(void):c_ptr(ls_hs_basis);
    }
  }
  */

  proc foreignInterfaceType(type t) type {
    if t == ls_hs_basis then return BasisForeignInterface;
    compilerError("unknown type");
  }

  record HsWrapper {
    type eltType;
    var payload : c_ptr(eltType);
    var localeId : chpl_nodeID_t;

    inline proc foreignInterface { return new foreignInterfaceType(eltType)(); }

    proc init(ptr : c_ptr(?eltType), owning : bool = true) {
      this.eltType = eltType;
      init this;
      this.payload = if owning then ptr else foreignInterface.clone(ptr);
      this.localeId = chpl_nodeID;
    }

    proc init(const ref x : HsWrapper(?eltType)) {
      this.eltType = eltType;
      init this;

      if compiledForSingleLocale() || x.localeId == chpl_nodeID {
        this.payload = foreignInterface.clone(x.payload);
        this.localeId = chpl_nodeID;
      }
      else {
        // We serialize the value to JSON on x.localeId, copy the JSON representation to chpl_nodeID,
        // and then reconstruct the eltType.
        //
        // the following like is taken from the BigInteger module
        const xLoc = chpl_buildLocaleID(x.localeId, c_sublocid_any);
        const p = x.payload;
        var jsonString : string;
        // the following like is taken from the BigInteger module
        on __primitive("chpl_on_locale_num", xLoc) {
          jsonString = foreignInterface.to_json(p);
        }
        this.payload = try! foreignInterface.from_json(jsonString.localize().c_str());
        this.localeId = chpl_nodeID;
      }
    }

    proc init(type eltType, jsonString : string) {
      this.eltType = eltType;
      init this;
      const s = jsonString.localize();
      this.payload = try! foreignInterface.from_json(s.c_str());
      this.localeId = chpl_nodeID;
    }

    proc init=(const ref x : HsWrapper(?eltType)) do this.init(x);

    proc deinit() {
      assert(localeId == chpl_nodeID);
      foreignInterface.destroy(payload);
    }

    proc deref() ref : eltType {
      assert(compiledForSingleLocale() || localeId == chpl_nodeID);
      return payload.deref();
    }
  }

  record Basis {
    var wrapper : HsWrapper(ls_hs_basis);

    inline proc payload { return wrapper.payload; }
    inline proc info const ref {
      const info = ls_chpl_get_basis_info(payload);
      assert(info != nil);
      return info.deref();
    }

    proc init(p : c_ptr(ls_hs_basis), owning : bool = true) {
      this.wrapper = new HsWrapper(p, owning);
    }
    proc init(p : c_ptrConst(ls_hs_basis), owning : bool) {
      assert(owning == false);
      init(p:c_ptr(ls_hs_basis), owning);
    }
    proc init(jsonString : string) {
      this.wrapper = new HsWrapper(ls_hs_basis, jsonString);
    }
    proc init(const ref x : Basis) {
      this.wrapper = x.wrapper;
    }
    proc init=(const ref x : Basis) do this.init(x);

    // proc build() { ls_hs_basis_build(wrapper.payload); }

    // proc uncheckedSetRepresentatives(representatives : [] uint(64)) {
    //   var arr = unsafeViewAsExternalArray(representatives);
    //   ls_hs_unchecked_set_representatives(wrapper.payload, c_ptrTo(arr),
    //                                       kCacheNumberBits:c_int);
    // }

    // proc isSpinBasis() { return wrapper.deref().particle_type == LS_HS_SPIN; }
    // proc isSpinfulFermionicBasis() { return wrapper.deref().particle_type == LS_HS_SPINFUL_FERMION; }
    // proc isSpinlessFermionicBasis() { return wrapper.deref().particle_type == LS_HS_SPINLESS_FERMION; }

    // proc isHammingWeightFixed() {
    //   // logDebug("ls_hs_basis_has_fixed_hamming_weight");
    //   return info.hamming_weight.hasValue;
    //   // ls_hs_basis_has_fixed_hamming_weight(wrapper.payload);
    // }
    // proc hasSpinInversionSymmetry() {
    //   return info.spin_inversion.hasValue();
    //   // ls_hs_basis_has_spin_inversion_symmetry(wrapper.payload);
    // }
    // inline proc hasPermutationSymmetries() { return info.has_permutation_symmetries; }

    // proc numberBits : int { return info.number_bits; }
    // proc numberWords : int { return info.number_words; }
    // proc numberSites() : int { return info.number_sites; }
    // proc numberParticles() : int {
    //   const r = try! info.number_particles.value();
    //   return r;
    // }
    // proc numberUp() : int { return wrapper.deref().number_up; }
    // proc spinInversion : int { return wrapper.deref().spin_inversion:int; }

    // proc minStateEstimate() : uint(64) {
    //   // logDebug("ls_hs_min_state_estimate");
    //   return ls_hs_min_state_estimate(wrapper.payload);
    // }
    // proc maxStateEstimate() : uint(64) {
    //   // logDebug("ls_hs_max_state_estimate");
    //   return ls_hs_max_state_estimate(wrapper.payload);
    // }

    // proc representatives() {
    //   ref rs = wrapper.deref().representatives;
    //   if rs.elts == nil then
    //     halt("basis is not built");
    //   return makeArrayFromExternArray(rs, uint(64));
    // }
  }

  // operator Basis.= (ref lhs : Basis, const ref rhs : Basis) {
  //   assert(lhs.locale == rhs.locale);
  //   assert(lhs._origin == lhs.locale);
  //   lhs._destroy();
  //   lhs.payload = ls_hs_clone_basis(rhs.payload);
  //   lhs.owning = true;
  //   lhs._origin = rhs._origin;
  //   lhs._json_repr = rhs._json_repr;
  // }
  proc Basis.getLocaleIdxFn() {
    record LocaleIndexFn {
      inline proc hash64_01(in x: uint(64)): uint(64) {
        x = (x ^ (x >> 30)) * (0xbf58476d1ce4e5b9:uint(64));
        x = (x ^ (x >> 27)) * (0x94d049bb133111eb:uint(64));
        x = x ^ (x >> 31);
        return x;
      }

      inline proc this(basisState : uint(64)) : int where CHPL_COMM == "" { return 0; }
      inline proc this(basisState : uint(64)) : int where CHPL_COMM != "" { return (hash64_01(basisState) % numLocales:uint):int; }

      proc this(count : int, basisStates : c_ptrConst(uint(64)), targetLocales : c_ptr(uint(8))) {
        if numLocales > 1 then
          foreach i in 0 ..# count do
            targetLocales[i] = this(basisStates[i]):uint(8);
        else
          POSIX.memset(targetLocales, 0, count:c_size_t);
      }

    }

    return new LocaleIndexFn();
  }

  proc isRepresentative(const ref basis : Basis,
                        const ref alphas : [?D1] uint(64),
                        ref norms : [?D2] real(64)) where D1.rank == 1 && D2.rank == 1 {
    const count = alphas.size;
    assert(norms.size == count,
           "alphas.size (" + alphas.size:string + ") != norms.size ("
           + norms.size:string + ")");

    assert(basis.info.has_permutation_symmetries,
           "do not call isRepresentative on bases that do not have permutation symmetries");
    const kernel = ls_chpl_get_is_representative_kernel(basis.payload);
    assert(kernel != nil, "is_representative_kernel was not initialized");

    ls_chpl_invoke_is_representative_kernel(kernel, count, c_ptrToConst(alphas), c_ptrTo(norms));
  }
  proc isRepresentative(const ref basis : Basis, const ref alphas : [?D] uint(64)) where D.rank == 1 {
    var norms : [0 ..# alphas.size] real(64) = noinit;
    isRepresentative(basis, alphas, norms);
    return norms;
  }

  /*
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
      // assert(basis._origin == here);
      assert(expression.locale == here);

      this.payload = ls_hs_create_operator(basis.wrapper.payload, expression);
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
      if owning {
        logDebug("Destroying operator ", payload);
        ls_hs_destroy_operator(payload);
      }
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
  */

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
