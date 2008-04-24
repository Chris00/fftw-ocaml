(* File: fftw3SD.ml

   Copyright (C) 2006-2008

     Christophe Troestler <chris_77@users.sourceforge.net>
     WWW: http://math.umh.ac.be/an/software/

   This library is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License version 2.1 or
   later as published by the Free Software Foundation, with the special
   exception on linking described in the file LICENSE.

   This library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the file
   LICENSE for more details. *)

(* FFTW3 interface for Single/Double precision *)

open Bigarray
open Printf

type 'a fftw_plan (* single and double precision plans are different *)

(* Types of plans *)
type c2c
type r2c
type c2r
type r2r

IFDEF SINGLE_PREC THEN
INCLUDE "fftw3S_external.ml"
ELSE
INCLUDE "fftw3D_external.ml"
ENDIF
;;

type 'a plan = {
  plan: 'a fftw_plan;
  i : 'b 'c 'd. ('b,'c,'d) Genarray.t; (* hold input array => not
                                          freed by GC before the plan *)
  offseto : int; (* output offset; C-stubs *)
  strideo : int array; (* strides; C-stubs *)
  no : int array; (* dimensions *)
  o : 'b 'c 'd. ('b,'c,'d) Genarray.t; (* output array *)
  normalize : bool; (* whether to normalize the output *)
  normalize_factor : float; (* multiplication factor to normalize *)
}

type dir = Forward | Backward
type measure = Estimate | Measure | Patient | Exhaustive
type r2r_kind =
    | R2HC | HC2R | DHT | REDFT00
    | REDFT10 | REDFT01 | REDFT11 | RODFT00 | RODFT10 | RODFT01 | RODFT11


let sign_of_dir = function
  | Forward -> -1
  | Backward -> 1

(* WARNING: keep in sync with fftw3.h *)
let flags meas unaligned preserve_input : int =
  let f = match meas with
    | Measure -> 0 (* 0U *)
    | Exhaustive -> 8 (* 1U lsl 3 *)
    | Patient -> 32 (* 1U lsl 5 *)
    | Estimate -> 64 (* 1U lsl 6 *) in
  let f = if unaligned then f lor 2 (* 1U lsl 1 *) else f in
  if preserve_input then f lor 16 (* 1U lsl 4 *) else f lor 1 (* 1U lsl 0 *)

(* WARNING: keep in sync with fftw3.h *)
let int_of_r2r_kind = function
  | R2HC    -> 0
  | HC2R    -> 1
  | DHT     -> 2
  | REDFT00 -> 3
  | REDFT01 -> 4
  | REDFT10 -> 5
  | REDFT11 -> 6
  | RODFT00 -> 7
  | RODFT01 -> 8
  | RODFT10 -> 9
  | RODFT11 -> 10


(** {2 Execution of plans}
 ***********************************************************************)

let exec p =
  fftw_exec p.plan;
  if p.normalize then
    normalize p.o p.offseto p.strideo p.no p.normalize_factor

let exec_dft plan i o =
  (* how to check that the arrays conform to the plan specification? *)
  exec_dft plan i o

let exec_split_dft plan ri ii ro io =
  (* again, how to check conformance with the plan? *)
  exec_split_dft plan ri ii ro io



(** {2 Creating plans}
 ***********************************************************************)

module Genarray = struct
  external create: ('a, 'b) Bigarray.kind -> 'c Bigarray.layout ->
    int array -> ('a, 'b, 'c) Bigarray.Genarray.t
    = "fftw3_ocaml_ba_create"

  type 'l complex_array = (Complex.t, complex_elt, 'l) Genarray.t
  type 'l float_array   = (float, float_elt, 'l) Genarray.t

  let is_c_layout m = (Genarray.layout m = (Obj.magic c_layout : 'a layout))

  (** [get_rank default m] returns the length of by the first array in
      the list of options [m]. *)
  let rec get_rank default = function
    | [] -> default
    | None :: t -> get_rank default t
    | Some m :: _ -> Array.length m

  let get_mat_rank name rank default = function
    | None -> Array.make rank default (* Create matrix with default value *)
    | Some m ->
        if Array.length m <> rank then
          invalid_arg(sprintf "%s: expected length=%i, got=%i"
                         name rank (Array.length m));
        m

  (** Return a string showing the content of the array *)
  let string_of_array a =
    assert(Array.length a >= 1);
    let b = Buffer.create 80 in
    Buffer.add_string b "[|";
    Buffer.add_string b (string_of_int a.(0));
    for i = 1 to Array.length a - 1 do
      Buffer.add_string b "; ";
      Buffer.add_string b (string_of_int a.(i));
    done;
    Buffer.add_string b "|]";
    Buffer.contents b

  (* C layout *)
  module C =
  struct
    DEFINE LAYOUT = "c_layout";;
    DEFINE FIRST_INDEX = 0;;
    DEFINE LAST_INDEX(dim) = dim - 1;;
    DEFINE FOR_DIM(k, rank, expr) = for k = rank - 1 downto 0 do expr done;;
    DEFINE FOR_HM(k, rank, expr) = for k = 0 to rank - 1 do expr done;;
    DEFINE LT_DIM_SYM = "<";;
    DEFINE GE_DIM_SYM = ">=";;
    INCLUDE "fftw3SD_genarray.ml"
  end

  (* FORTRAN layout *)
  module F =
  struct
    DEFINE FORTRAN (* BEWARE it is still defined after this module! *)
    DEFINE LAYOUT = "fortran_layout";;
    DEFINE FIRST_INDEX = 1;;
    DEFINE LAST_INDEX(dim) = dim;;
    DEFINE FOR_DIM(k, rank, expr) = for k = 0 to rank - 1 do expr done;;
    DEFINE FOR_HM(k, rank, expr) = for k = rank - 1 downto 0 do expr done;;
    DEFINE LT_DIM_SYM = "<=";;
    DEFINE GE_DIM_SYM = ">";;
    INCLUDE "fftw3SD_genarray.ml"
  end

  (* Layout independent function *)
  let apply name wrapper n hm_n  hmi ofsi inci i  hmo ofso inco o =
    (if is_c_layout i then C.apply else F.apply)
      name wrapper n hm_n  hmi ofsi inci i  hmo ofso inco o

  let dft dir ?(meas=Measure) ?(normalize=false)
      ?(preserve_input=false) ?(unaligned=false) ?n ?(howmany_n=[| |])
      ?(howmanyi=[]) ?ofsi ?inci (i: 'l complex_array)
      ?(howmanyo=[]) ?ofso ?inco (o: 'l complex_array) =
    apply "Fftw3.D.Genarray.dft"
      (guru_dft i o (sign_of_dir dir) (flags meas unaligned preserve_input))
      n howmany_n  howmanyi ofsi inci i  howmanyo ofso inco o

  (* At the moment, in place transforms are not possible but they may
     be if OCaml bug 0004333 is resolved. *)
  let r2c ?(meas=Measure) ?(normalize=false)
      ?(preserve_input=false) ?(unaligned=false) ?n
      ?howmany_ranki ?howmanyi  ?ofsi ?inci (i: 'l float_array)
      ?howmany_ranko ?howmanyo  ?ofso ?inco (o: 'l complex_array) =
    apply "Fftw3.D.Genarray.r2c"
      (guru_r2c i o (flags meas unaligned preserve_input))
      howmany_rank howmany howmany_ofsi howmany_inci howmany_ofso howmany_inco
      n ofsi inci i ofso inco o

  let c2r ?(meas=Measure) ?(normalize=false)
      ?(preserve_input=false) ?(unaligned=false) ?n
      ?howmany_ranki ?howmanyi  ?ofsi ?inci (i: 'l float_array)
      ?howmany_ranko ?howmanyo  ?ofso ?inco (o: 'l complex_array) =
    (* FIXME:  Are the checks of apply appropriate here? *)
    apply "Fftw3.D.Genarray.c2r"
      (guru_c2r i o (flags meas unaligned preserve_input))
      howmany_rank howmany howmany_ofsi howmany_inci howmany_ofso howmany_inco
      n ofsi inci i ofso inco o

  let r2r kind ?(meas=Measure) ?(normalize=false)
      ?(preserve_input=true) ?(unaligned=false) ?n
      ?howmany_ranki ?howmanyi  ?ofsi ?inci (i: 'l float_array)
      ?howmany_ranko ?howmanyo  ?ofso ?inco (o: 'l complex_array) =
    (* FIXME: must check [kind] has the right length?? *)
    let kind = Array.map int_of_r2r_kind kind in
    apply "Fftw3.D.Genarray.r2r"
      (guru_r2r i o kind (flags meas unaligned preserve_input))
      (* howmany_rank: *)(Some(Genarray.num_dims i - Array.length kind))
      howmany howmany_ofsi howmany_inci howmany_ofso howmany_inco
      n ofsi inci i ofso inco o

end


module Array1 = struct
  external array1_of_ba : ('a,'b,'c) Bigarray.Genarray.t -> ('a,'b,'c) Array1.t
    = "%identity"
    (* We know that the bigarray will have only 1D, convert without check *)

  let create kind layout dim =
    array1_of_ba(Genarray.create kind layout [|dim|])

  let of_array kind layout data =
    let ba = create kind layout (Array.length data) in
    let ofs = if layout = (Obj.magic c_layout : 'a layout) then 0 else 1 in
    for i = 0 to Array.length data - 1 do ba.{i + ofs} <- data.(i) done;
    ba


  type 'l complex_array = (Complex.t, complex_elt, 'l) Array1.t
  type 'l float_array   = (float, float_elt, 'l) Array1.t

  let is_c_layout m = (Array1.layout m = (Obj.magic c_layout : 'a layout))


  let check_geom_c name ofsname ofs incname inc nname n dim =
    if inc = 0 then invalid_arg(sprintf "%s: %s = 0" name incname);
    if ofs < 0 then invalid_arg(sprintf "%s: %s < 0" name ofsname);
    if ofs >= dim then invalid_arg(sprintf "%s: %s >= %i" name ofsname dim);
    if inc > 0 && ofs + (n - 1) * inc >= dim then
      invalid_arg(sprintf "%s: %s + (%s-1)*%s = %i >= %i"
                     name ofsname nname incname (ofs + (n - 1) * inc) dim);
    if inc < 0 && ofs + (n - 1) * inc < 0 then
      invalid_arg(sprintf "%s: %s + (%s-1)*%s = %i < 0"
                     name ofsname nname incname (ofs + (n - 1) * inc))


  let apply_c name wrapper  hm hm_inci hm_inco
      ?n ?(ofsi=0) inci i ?(ofso=0) inco o =
    let dim = Array1.dim i in
    let n = match n with
      | None ->
          if inci = 0 then invalid_arg(sprintf "%s: inci = 0" name);
          if inci > 0 then 1 + (dim - 1 - ofsi) / inci else 1 - ofsi / inci
      | Some n -> n in
    check_geom_c name "ofsi" ofsi "inci" inci "n" n dim;
    check_geom_c name "ofso" ofso "inco" inco "n" n dim;
    (* If howmany > 1, check that matrices do not overlap. *)
    let howmany, hm_istride, hm_ostride =
      if hm <= 1 then [| |], [| |], [| |]
      else
        let hm_inci = match hm_inci with None -> ofsi + n*inci | Some i -> i
        and hm_inco = match hm_inco with None -> ofso + n*inco | Some i -> i in

        [|hm|], [|hm_inci|], [|hm_inco|] in
    wrapper ofsi ofso [|n|] [|inci|] [|inco|] howmany hm_istride hm_ostride


  let check_geom_f name ofsname ofs incname inc nname n dim =
    if inc = 0 then invalid_arg(sprintf "%s: %s = 0" name incname);
    if ofs < 0 then invalid_arg(sprintf "%s: %s < 0" name ofsname);
    if ofs > dim then invalid_arg(sprintf "%s: %s > %i" name ofsname dim);
    if inc > 0 && ofs + (n - 1) * inc > dim then
      invalid_arg(sprintf "%s: %s + (%s-1)*%s = %i > %i"
                     name ofsname nname incname (ofs + (n - 1) * inc) dim);
    if inc < 0 && ofs + (n - 1) * inc < 1 then
      invalid_arg(sprintf "%s: %s + (%s-1)*%s = %i < 1"
                     name ofsname nname incname (ofs + (n - 1) * inc))

  let apply_f name wrapper  hm hm_inci hm_inco
      ?n ?(ofsi=1) inci i ?(ofso=1) inco o =
    (* FIXME *)
    let dim = Array1.dim i in
    let n = match n with
      | None ->
          if inci = 0 then invalid_arg(sprintf "%s: inci = 0" name);
          if inci > 0 then 1 + (dim - ofsi) / inci
          else 1 - (ofsi - 1) / inci
      | Some n -> n in
    check_geom_f name "ofsi" ofsi "inci" inci "n" n dim;
    check_geom_f name "ofso" ofso "inco" inco "n" n dim;
    (* If howmany > 1, check that matrices do not overlap. *)
    let howmany, hm_istride, hm_ostride =
      if hm <= 1 then [| |], [| |], [| |]
      else
        let hm_inci = match hm_inci with None -> ofsi + n*inci | Some i -> i
        and hm_inco = match hm_inco with None -> ofso + n*inco | Some i -> i in

        [|hm|], [|hm_inci|], [|hm_inco|] in
    wrapper (ofsi - 1) (ofso - 1) [|n|] [|inci|] [|inco|]
      howmany hm_istride hm_ostride


  let apply name wrapper  hm hm_inci hm_inco ?n ?ofsi inci i ?ofso inco o =
    if is_c_layout i then
      apply_c name wrapper  hm hm_inci hm_inco ?n ?ofsi inci i ?ofso inco o
    else
      apply_f name wrapper  hm hm_inci hm_inco ?n ?ofsi inci i ?ofso inco o


  let dft dir ?(meas=Measure) ?(normalize=false)
      ?(preserve_input=true) ?(unaligned=false)
      ?n ?(howmany=1)
      ?howmany_inci ?ofsi ?(inci=1) (i: 'l complex_array)
      ?howmany_inco ?ofso ?(inco=1) (o: 'l complex_array) =
    let gi = genarray_of_array1 i
    and go = genarray_of_array1 o in
    apply "Fftw3.D.Array1.dft"
      (guru_dft gi go (sign_of_dir dir) (flags meas unaligned preserve_input))
      howmany howmany_inci howmany_inco ?n ?ofsi inci i ?ofso inco o

  let r2r kind ?(meas=Measure) ?(normalize=false)
      ?(preserve_input=true) ?(unaligned=false)
      ?n ?(howmany=1)
      ?howmany_inci ?ofsi ?(inci=1) (i: 'l float_array)
      ?howmany_inco ?ofso ?(inco=1) (o: 'l float_array) =
    let gi = genarray_of_array1 i
    and go = genarray_of_array1 o in
    let kind = int_of_r2r_kind kind in
    apply "Fftw3.D.Array1.r2r"
      (guru_r2r gi go [|kind|] (flags meas unaligned preserve_input))
      howmany howmany_inci howmany_inco ?n ?ofsi inci i ?ofso inco o
end


module Array2 = struct
  external array2_of_ba : ('a,'b,'c) Bigarray.Genarray.t -> ('a,'b,'c) Array2.t
    = "%identity"
    (* BEWARE: only for bigarray with 2D, convert without check *)

  let create kind layout dim1 dim2 =
    array2_of_ba(Genarray.create kind layout [|dim1; dim2|])

  type 'l complex_array = (Complex.t, complex_elt, 'l) Array2.t
  type 'l float_array   = (float, float_elt, 'l) Array2.t

  let is_c_layout m = (Array2.layout m = (Obj.magic c_layout : 'a layout))
end


module Array3 = struct
  external array3_of_ba : ('a,'b,'c) Bigarray.Genarray.t -> ('a,'b,'c) Array3.t
    = "%identity"
    (* BEWARE: only for bigarray with 3D, convert without check *)

  let create kind layout dim1 dim2 dim3 =
    array3_of_ba(Genarray.create kind layout [|dim1; dim2; dim3|])

  type 'l complex_array = (Complex.t, complex_elt, 'l) Array3.t
  type 'l float_array   = (float, float_elt, 'l) Array3.t

  let is_c_layout m = (Array3.layout m = (Obj.magic c_layout : 'a layout))
end
