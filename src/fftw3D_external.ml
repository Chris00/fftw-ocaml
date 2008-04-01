(* External functions for the fftw3 stubs (in double precision for the
   original file). *)

(* The types for Array1,... can be converted to this at no cost. *)
type 'l complex_array = (Complex.t, complex_elt, 'l) Genarray.t
type 'l float_array   = (float, float_elt, 'l) Genarray.t


(** {2 Execution of plans}
 ***********************************************************************)

external normalize :
  (* array *) (_,_,_) Genarray.t ->
  (* offset *) int ->
  (* strides *) int array ->
  (* dimensions *) int array ->
  (* multiply all entries by this number *) float -> unit
  = "fftw_ocaml_normalize" "noalloc"

external fftw_exec : 'a fftw_plan -> unit = "fftw_ocaml_execute" "noalloc"

external exec_dft : c2c fftw_plan -> 'l complex_array -> 'l complex_array -> unit
  = "fftw_ocaml_execute_dft" "noalloc"
external exec_split_dft : c2c fftw_plan -> 'l float_array -> 'l float_array ->
  'l float_array -> 'l float_array -> unit
  = "fftw_ocaml_execute_split_dft" "noalloc"


(** {2 Creating plans}
 ***********************************************************************)

(* BEWARE: wrapper functions are just thin wrappers around their C
   counterpart.  In particular, their arguments must be thought for
   the C layout. *)
external guru_dft :
  (* in *)  'l complex_array ->
  (* out *) 'l complex_array ->
  (* sign (forward/backward) *) int ->
  (* flags (GOOD: they do not use the 32th bit) *) int ->
  (* input offset (as 1D array) *) int ->
  (* output offset (as 1D array) *) int ->
  (* n (transform dimensions; its length = rank) *) int array ->
  (* istride (same length as [n]) *) int array ->
  (* ostride (same length as [n]) *) int array ->
  (* howmany (multiplicity dimensions; its length=howmany_rank) *) int array ->
  (* howmany input strides (same length as [howmany]) *) int array ->
  (* howmany output strides (same length as [howmany]) *) int array
  -> c2c fftw_plan
  = "fftw_ocaml_guru_dft_bc" "fftw_ocaml_guru_dft"
  (* Wrapper of fftw_plan_guru_dft.  No coherence check is done in the
     C code.  @raise Failure if the plan cannot be created. *)

external guru_r2c :
  (* in *) 'l float_array ->
  (* out *) 'l complex_array ->
  (* flags *) int ->
  (* input offset *) int ->
  (* output offset *) int ->
  (* n (transform dimensions) *) int array ->
  (* istride (same length as [n]) *) int array ->
  (* ostride (same length as [n]) *) int array ->
  (* howmany (multiplicity dimensions) *) int array ->
  (* howmany input strides (same length as [howmany]) *) int array ->
  (* howmany output strides (same length as [howmany]) *) int array
  -> r2c fftw_plan
  = "fftw_ocaml_guru_r2c_bc" "fftw_ocaml_guru_r2c"

external guru_c2r :
  (* in *) 'l complex_array ->
  (* out *) 'l float_array ->
  (* flags *) int ->
  (* input offset *) int ->
  (* output offset *) int ->
  (* n (transform dimensions) *) int array ->
  (* istride (same length as [n]) *) int array ->
  (* ostride (same length as [n]) *) int array ->
  (* howmany (multiplicity dimensions) *) int array ->
  (* howmany input strides (same length as [howmany]) *) int array ->
  (* howmany output strides (same length as [howmany]) *) int array
  -> c2r fftw_plan
  = "fftw_ocaml_guru_c2r_bc" "fftw_ocaml_guru_c2r"

external guru_r2r :
  (* in *) 'l float_array ->
  (* out *) 'l float_array ->
  (* kind (same length as [n]) *) int array ->
  (* flags *) int ->
  (* input offset *) int ->
  (* output offset *) int ->
  (* n (transform dimensions) *) int array ->
  (* istride (same length as [n]) *) int array ->
  (* ostride (same length as [n]) *) int array ->
  (* howmany (multiplicity dimensions) *) int array ->
  (* howmany input strides (same length as [howmany]) *) int array ->
  (* howmany output strides (same length as [howmany]) *) int array
  -> r2r fftw_plan
  = "fftw_ocaml_guru_r2r_bc" "fftw_ocaml_guru_r2r"
