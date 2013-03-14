(* Generate precision dependent files.  Avoid external dependencies
   such as "sed". *)

#load "str.cma";;
open Printf

let debug = false

let string_of_file fname =
  let buf = Buffer.create 8192 in
  let fh = open_in fname in
  let s = String.create 4096 in
  let n = ref 1 (* enter the loop *) in
  while !n > 0 do
    n := input fh s 0 4096;
    Buffer.add_substring buf s 0 !n
  done;
  close_in fh;
  Buffer.contents buf

let transform fin fout tr =
  let s = string_of_file fin in
  let replace s (re, sub) =
    Str.global_replace (Str.regexp re) sub s in
  let s = List.fold_left replace s tr in
  (* Output *)
  let fh = open_out fout in
  fprintf fh "(* AUTOMATICALLY GENERATED from %S. *)\n" fin;
  output_string fh s;
  close_out fh


let () =
  transform "fftw3SD.ml" "fftw3D.ml"
            ["floatXX_elt", "Bigarray.float64_elt";
             "floatXX", "Bigarray.float64";
             "complexXX_elt", "Bigarray.complex64_elt";
             "complexXX", "Bigarray.complex64";
             "$FFTW", "Fftw3.D"];
  transform "fftw3SD.ml" "fftw3S.ml"
            ["fftw_ocaml", "fftwf_ocaml"; (* C stubs *)
             "floatXX_elt", "Bigarray.float32_elt";
             "floatXX", "Bigarray.float32";
             "complexXX_elt", "Bigarray.complex32_elt";
             "complexXX", "Bigarray.complex32";
             "$FFTW", "Fftw3.S"]

let () =
  let for_dim = "FOR_DIM(\\([a-z0-9_]+\\), *\\([a-z0-9_]+\\))" in
  let for_hm =  "FOR_HM(\\([a-z0-9_]+\\), *\\([a-z0-9_]+\\))" in
  let debug = "DEBUG{\\([^}]+\\)};", if debug then "\\1" else "" in
  transform "fftw3_geomCF.ml" "fftw3_geomC.ml"
            ["$LAYOUT", "c_layout";
             "FIRST_INDEX", "0";
             "LAST_INDEX(\\([a-z0-9]+\\))", "\\1 - 1";
             for_dim, "\\1 = \\2 - 1 downto 0 "; (* do ... done *)
             for_hm,  "\\1 = 0 to \\2 - 1 ";
             "$LT", "<";
             "$GE", ">=";
             debug];
  transform "fftw3_geomCF.ml" "fftw3_geomF.ml"
            ["$LAYOUT", "fortran_layout";
             "FIRST_INDEX", "1";
             "LAST_INDEX(\\([a-z0-9]+\\))", "\\1";
             for_dim, "\\1 = 0 to \\2 - 1 ";
             for_hm, "\\1 = \\2 - 1 downto 0 ";
             "$LT", "<=";
             "$GE", ">";
             debug]
