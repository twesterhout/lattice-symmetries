module LatticeSymmetries {
  public use FFI;
  // public use MyHDF5;
  public use ForeignTypes;
  public use StatesEnumeration;
  // public use MatrixVectorProduct;
  public use DistributedMatrixVector;
  public use BatchedOperator;
  public use ConcurrentAccessor;
  // public use CommunicationQueue;
  // public use MultiwayMerge;
  public use Vector;

  private use CTypes;

  proc initExportedKernels() {
    var kernels = new ls_chpl_kernels(
      enumerate_states=c_ptrTo(ls_chpl_enumerate_representatives),
      operator_apply_off_diag=c_ptrTo(ls_chpl_operator_apply_off_diag),
      operator_apply_diag=c_ptrTo(ls_chpl_operator_apply_diag),
      matrix_vector_product=c_ptrTo(ls_chpl_matrix_vector_product)
    );
    logDebug("Initializing chpl_kernels ...");
    ls_hs_internal_set_chpl_kernels(c_ptrTo(kernels));
  }

  export proc ls_chpl_init_kernels() {
    initExportedKernels();
  }

  initExportedKernels();
}
