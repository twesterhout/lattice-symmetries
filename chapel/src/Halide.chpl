module Halide {

// halide_type_code_t enum
extern type halide_type_code_t = c_int;
extern const halide_type_int :halide_type_code_t;
extern const halide_type_uint :halide_type_code_t;
extern const halide_type_float :halide_type_code_t;
extern const halide_type_handle :halide_type_code_t;
extern const halide_type_bfloat :halide_type_code_t;

extern "struct halide_type_t" record halide_type_t {
  var code : uint(8);
  var bits : uint(8);
  var lanes : uint(16);
}

extern "struct halide_device_interface_t" record halide_device_interface_t {
  var device_malloc : c_fn_ptr;
  var device_free : c_fn_ptr;
  var device_sync : c_fn_ptr;
  var device_release : c_fn_ptr;
  var copy_to_host : c_fn_ptr;
  var copy_to_device : c_fn_ptr;
  var device_and_host_malloc : c_fn_ptr;
  var device_and_host_free : c_fn_ptr;
  var buffer_copy : c_fn_ptr;
  var device_crop : c_fn_ptr;
  var device_slice : c_fn_ptr;
  var device_release_crop : c_fn_ptr;
  var wrap_native : c_fn_ptr;
  var detach_native : c_fn_ptr;
  var compute_capability : c_fn_ptr;
  var impl : c_ptrConst(halide_device_interface_impl_t);
}

extern record halide_dimension_t {
  var min : int(32);
  var extent : int(32);
  var stride : int(32);
  var flags : uint(32);
}

}
