module ForeignTypes {
  use FFI;
  // use CommonParameters;
  import Communication;
  use Timing;
  use Utils;

  import Reflection.getRoutineName;
  use CTypes;
  use ChplConfig;
  use IO;
  use JSON;
  // use Time;

  require "lattice_symmetries_chapel.h";

  extern proc ls_chpl_get_basis_info(p : c_ptrConst(ls_hs_basis)) : c_ptrConst(ls_hs_basis_info);
  extern proc ls_chpl_get_is_representative_kernel(p : c_ptrConst(ls_hs_basis)) : c_fn_ptr;
  extern proc ls_chpl_invoke_is_representative_kernel(kernel : c_fn_ptr, count : int(64),
                                                      basis_states : c_ptrConst(uint(64)), norms : c_ptr(uint(16)));
  extern proc ls_chpl_get_state_to_index_kernel(p : c_ptrConst(ls_hs_basis)) : c_fn_ptr;
  extern proc ls_chpl_invoke_state_to_index_kernel(kernel : c_fn_ptr, count : int(64),
                                                   basis_states : c_ptrConst(uint(64)), indices : c_ptr(int(64)));
  extern proc ls_chpl_get_state_info_kernel(p : c_ptrConst(ls_hs_basis)) : c_fn_ptr;
  extern proc ls_chpl_invoke_state_info_kernel(kernel : c_fn_ptr, count : int(64),
                                               basis_states : c_ptrConst(uint(64)),
                                               representatives : c_ptr(uint(64)), indices : c_ptr(int(32)));

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

  record OperatorForeignInterface {
    proc clone(p : c_ptrConst(ls_hs_operator)) {
      assert(p != nil);
      const mp = p:c_ptr(ls_hs_operator);
      ls_hs_internal_object_inc_ref_count(c_ptrTo(mp.deref().base));
      return mp;
    }

    proc destroy(p : c_ptr(ls_hs_operator)) {
      assert(p != nil);
      ls_hs_destroy_operator(p);
    }

    proc to_json(p : c_ptrConst(ls_hs_operator)) {
      const c_str = ls_hs_operator_to_json(p);
      defer ls_hs_destroy_string(c_str);
      return try! string.createCopyingBuffer(c_str);
    }

    proc from_json(s : c_ptrConst(c_char)) throws {
      const c_str = ls_hs_operator_from_json(s);
      defer ls_hs_destroy_string(c_str);
      return from_foreign_result(c_str);
    }

    proc from_foreign_result(c_str : c_ptrConst(c_char)) throws {
      var file = openMemFile();
      var writer = file.writer();
      writer.write(string.createBorrowingBuffer(c_str));
      writer.close();
      var r = file.reader(deserializer=new jsonDeserializer());
      var p = r.read(HsResult(uint(64))).value;
      return p:c_ptr(void):c_ptr(ls_hs_operator);
    }
  }

  proc foreignInterfaceType(type t) type {
    if t == ls_hs_basis then return BasisForeignInterface;
    if t == ls_hs_operator then return OperatorForeignInterface;
    compilerError("unknown type");
  }

  record HsWrapper {
    type eltType;
    var payload : c_ptr(eltType);
    var localeId : chpl_nodeID_t;

    inline proc foreignInterface { return new foreignInterfaceType(eltType)(); }

    proc init(type eltType) {
      this.eltType = eltType;
      this.payload = nil;
      this.localeId = chpl_nodeID;
    }

    proc init(ptr : c_ptr(?eltType), owning : bool = true) {
      this.eltType = eltType;
      init this;
      if owning && ptr == nil then
        halt("HsWrapper.init: if owning=true, then ptr must not be nil");
      this.payload = if owning then ptr else foreignInterface.clone(ptr);
      this.localeId = chpl_nodeID;
    }

    proc init(const ref x : HsWrapper(?eltType)) {
      this.eltType = eltType;
      this.localeId = chpl_nodeID;
      init this;

      if x.payload == nil { // Trivial case, nothing to initialize
        this.payload = nil;
      }
      else if compiledForSingleLocale() || x.localeId == chpl_nodeID { // Local copy constructor
        this.payload = foreignInterface.clone(x.payload);
      }
      else {
        // We serialize the value to JSON on x.localeId, copy the JSON representation to chpl_nodeID,
        // and then reconstruct the eltType.
        const xLoc = chpl_buildLocaleID(x.localeId, c_sublocid_any);
        const p = x.payload;
        var jsonString : string;
        on __primitive("chpl_on_locale_num", xLoc) do
          jsonString = foreignInterface.to_json(p);
        this.payload = try! foreignInterface.from_json(jsonString.localize().c_str());
      }
    }
    proc init=(const ref x : HsWrapper(?eltType)) do this.init(x);

    proc init(type eltType, jsonString : string) {
      this.eltType = eltType;
      this.localeId = chpl_nodeID;
      init this;

      const s = jsonString.localize();
      this.payload = try! foreignInterface.from_json(s.c_str());
    }

    proc deinit() {
      assert(localeId == chpl_nodeID);
      if payload != nil then
        foreignInterface.destroy(payload);
    }

    proc deref() ref : eltType {
      assert(compiledForSingleLocale() || localeId == chpl_nodeID);
      return payload.deref();
    }
  }

  operator HsWrapper.=(ref lhs : HsWrapper(?eltType), const ref rhs : HsWrapper(eltType)) {
    // logDebug("HsWrapper.=");
    var t = rhs;
    lhs.payload <=> t.payload;
    lhs.localeId <=> t.localeId;
  }

  proc getBasisInfo(ptr : c_ptr(ls_hs_basis)) const ref : ls_hs_basis_info {
      const info = ls_chpl_get_basis_info(ptr);
      if info == nil then
        halt("failed to initialize basis_info");
      return info.deref();
  }

  record Basis {
    var wrapper : HsWrapper(ls_hs_basis);

    inline proc payload {
      assert(wrapper.payload != nil);
      return wrapper.payload;
    }
    inline proc info const ref { return getBasisInfo(payload); }

    // TODO: figure out a way to remote this
    proc init() { }

    proc init(p : c_ptr(ls_hs_basis), owning : bool = true) {
      // logDebug("Basis.init(p, owning)");
      this.wrapper = new HsWrapper(p, owning);
    }
    proc init(p : c_ptrConst(ls_hs_basis), owning : bool) {
      assert(owning == false);
      init(p:c_ptr(ls_hs_basis), owning);
    }
    proc init(jsonString : string) {
      // logDebug("Basis.init(string)");
      this.wrapper = new HsWrapper(ls_hs_basis, jsonString);
    }
    proc init(const ref x : Basis) {
      // logDebug("Basis.init(Basis)");
      this.wrapper = x.wrapper;
    }
    proc init=(const ref x : Basis) do this.init(x);

    proc deinit() {
      // logDebug("Basis.deinit");
    }
  }

  operator Basis.=(ref lhs : Basis, const ref rhs : Basis) {
    // logDebug("Basis.=");
    var t = rhs;
    lhs.wrapper <=> t.wrapper;
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
      inline proc this(basisState : uint(64)) : int where CHPL_COMM != "" {
        return (hash64_01(basisState) % numLocales:uint):int;
      }

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

  proc Basis.getSpinInversionMask() : uint(64) {
    const n = info.number_sites;
    if n == 64 then return max(uint(64));
    if n == 0 then return 0;
    return (1:uint(64) << n) - 1;
  }

  proc isRepresentative(const ref basis : Basis,
                        count : int,
                        alphas : c_ptrConst(uint(64)),
                        norms : c_ptr(uint(16))) {
    const _timer = recordTime(getRoutineName());
    if count <= 0 then return;
    if basis.info.has_permutation_symmetries {

      const kernel = ls_chpl_get_is_representative_kernel(basis.payload);
      if kernel == nil then halt("is_representative_kernel was not initialized");
      const countWithPadding = roundUpToMaxBlockSize(count);
      if count != countWithPadding {
        // halt(try! "is_representative_kernel expects 'count' to be a multiple of %n".format(
        //           LS_HS_MAX_BLOCK_SIZE));

        // The slow path with an extra memory allocation
        var tempAlphas : [0 ..# countWithPadding] uint(64) = noinit;
        POSIX.memcpy(c_ptrTo(tempAlphas), alphas, count:c_size_t * c_sizeof(uint(64)));
        var tempNorms : [0 ..# countWithPadding] uint(16);

        // logDebug(try! "countWithPadding=%n".format(countWithPadding));
        ls_chpl_invoke_is_representative_kernel(kernel, countWithPadding, c_ptrToConst(tempAlphas), c_ptrTo(tempNorms));
        POSIX.memcpy(norms, c_ptrToConst(tempNorms), count:c_size_t * c_sizeof(uint(16)));
      }
      else {
        // Normal optimized path
        ls_chpl_invoke_is_representative_kernel(kernel, count, alphas, norms);
      }

    }
    else if basis.info.requires_projection {

      // We just have the spin_inversion symmetry
      assert(basis.info.spin_inversion != 0);
      const inversionMask = basis.getSpinInversionMask();
      foreach k in 0 ..# count do
        if alphas[k] ^ inversionMask < alphas[k] then
          norms[k] = 0;
        else
          norms[k] = 1;

    }
    else {

      // No projection required: all norms are 1
      foreach k in 0 ..# count do
        norms[k] = 1;

    }
  }
  proc isRepresentative(const ref basis : Basis,
                        const ref alphas : [?D1] uint(64),
                        ref norms : [?D2] uint(16)) where D1.rank == 1 && D2.rank == 1 {
    if alphas.size != norms.size then
      halt(try! "alphas.size (%i) != norms.size (%i)".format(alphas.size, norms.size));
    isRepresentative(basis, alphas.size, c_ptrToConst(alphas), c_ptrTo(norms));
  }
  proc isRepresentative(const ref basis : Basis, const ref alphas : [?D] uint(64)) where D.rank == 1 {
    var norms : [0 ..# alphas.size] uint(16) = noinit;
    isRepresentative(basis, alphas, norms);
    return norms;
  }

  proc basisStatesToIndices(const ref basis : Basis, count : int, alphas : c_ptrConst(uint(64)), indices : c_ptr(int(64))) {
    const _timer = recordTime(getRoutineName());
    const kernel = ls_chpl_get_state_to_index_kernel(basis.payload);
    if kernel == nil then halt("state_to_index kernel was not initialized");
    const countWithPadding = roundUpToMaxBlockSize(count);
    // logDebug(try! "%n -> %n".format(count, countWithPadding));
    if count != countWithPadding {

      // The slow path with an extra memory allocation
      var tempAlphas : [0 ..# countWithPadding] uint(64) = noinit;
      POSIX.memcpy(c_ptrTo(tempAlphas), alphas, count:c_size_t * c_sizeof(uint(64)));
      var tempIndices : [0 ..# countWithPadding] int(64);

      // logDebug("init tempAlphas:  ", tempAlphas);
      // logDebug("init tempIndices: ", tempIndices);
      ls_chpl_invoke_state_to_index_kernel(kernel, countWithPadding, c_ptrToConst(tempAlphas), c_ptrTo(tempIndices));
      // logDebug("then tempAlphas:  ", tempAlphas);
      // logDebug("then tempIndices: ", tempIndices);

      POSIX.memcpy(indices, c_ptrToConst(tempIndices), count:c_size_t * c_sizeof(int(64)));
    }
    else {
      // Normal optimized path
      ls_chpl_invoke_state_to_index_kernel(kernel, count, alphas, indices);
    }
  }

  export proc ls_chpl_is_representative(basisPtr : c_ptrConst(ls_hs_basis),
                                        count : int(64),
                                        alphas : c_ptrConst(uint(64)),
                                        norms : c_ptr(uint(16))) {
    isRepresentative(new Basis(basisPtr, owning=false), count, alphas, norms);
  }

  export proc ls_chpl_basis_state_to_index(basisPtr : c_ptrConst(ls_hs_basis),
                                           count : int(64),
                                           alphas : c_ptrConst(uint(64)),
                                           indices : c_ptr(int(64))) {
    basisStatesToIndices(new Basis(basisPtr, owning=false), count, alphas, indices);
  }

  record Operator {
    var wrapper : HsWrapper(ls_hs_operator);
    var basis : Basis;

    inline proc payload { return wrapper.payload; }

    proc _getReference() ref {
      assert(payload != nil);
      return payload.deref();
    }
    forwarding _getReference();

    proc init(p : c_ptr(ls_hs_operator), owning : bool = true) {
      this.wrapper = new HsWrapper(p, owning);
      init this;

      this.basis = new Basis(this.wrapper.payload.deref().basis, owning=false);
    }
    proc init(p : c_ptrConst(ls_hs_operator), owning : bool) {
      assert(owning == false);
      init(p:c_ptr(ls_hs_operator), owning);
    }
    proc init(basis : c_ptrConst(ls_hs_basis), expr : c_ptrConst(ls_hs_expr), owning : bool) {
      assert(owning == false);
      const c_str = ls_hs_create_operator(basis, expr);
      defer ls_hs_destroy_string(c_str);
      init(wrapper.foreignInterface.from_foreign_result(c_str), owning=true);
    }
    proc init(jsonString : string) {
      this.wrapper = new HsWrapper(ls_hs_operator, jsonString);
      init this;
      // logDebug("Operator.init(string): refcount=", this.payload.deref().base.refcount);

      this.basis = new Basis(this.wrapper.payload.deref().basis, owning=false);
      // logDebug("Operator.init(string): refcount=", this.payload.deref().base.refcount,
      //          ", basis refcount=", this.basis.payload.deref().base.refcount);
    }
    proc init(const ref x : Operator) {
      this.wrapper = x.wrapper;
      this.basis = x.basis;
    }
    proc init=(const ref x : Operator) do this.init(x);

    proc deinit() {
      // logDebug("Operator.deinit");
    }
  }

  operator Operator.=(ref lhs : Operator, const ref rhs : Operator) {
    // logDebug("Operator.=");
    var t = rhs;
    lhs.wrapper <=> t.wrapper;
    lhs.basis <=> t.basis;
  }

  proc loadOperatorFromFile(jsonFile : string) throws {
    var f = open(jsonFile, ioMode.r);
    var r = f.reader();
    return new Operator(r.readAll(string));
  }





}
