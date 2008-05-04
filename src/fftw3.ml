(*pp camlp4o pa_macro.cmo $FFTW3F_EXIST *)
(* File: fftw3.ml

   Objective Caml interface for FFTW.

   Copyright (C) 2005-2008

     Christophe Troestler <chris_77@users.sourceforge.net>
     WWW: http://www.umh.ac.be/math/an/software/

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public License
   version 2.1 as published by the Free Software Foundation, with the
   special exception on linking described in file LICENSE.

   This library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the file
   LICENSE for more details.
*)

DEFINE DEBUG(expr) = expr;;
DEFINE DEBUG(expr) = ();;

open Bigarray
open Printf

module type Sig = sig
  type float_elt
    (** Precision of float numbers. *)
  type complex_elt
    (** Precision of complex numbers. *)

  val float : (float, float_elt) Bigarray.kind
  val complex : (Complex.t, complex_elt) Bigarray.kind

  type 'a plan (** Immutable FFTW plan. *)
  type c2c
  type r2c
  type c2r
  type r2r

  type dir = Forward | Backward

  type measure =
    | Estimate
    | Measure
    | Patient
    | Exhaustive


  type r2r_kind =
    | R2HC
    | HC2R
    | DHT
    | REDFT00
    | REDFT10
    | REDFT01
    | REDFT11
    | RODFT00
    | RODFT10
    | RODFT01
    | RODFT11

  val exec : 'a plan -> unit

  module Guru : sig

  end


  module Genarray :
  sig
    external create: ('a, 'b) kind -> 'c layout -> int array
      -> ('a, 'b, 'c) Bigarray.Genarray.t
      = "fftw3_ocaml_ba_create"

    type 'l complex_array = (Complex.t, complex_elt, 'l) Bigarray.Genarray.t
    type 'l float_array   = (float, float_elt, 'l) Bigarray.Genarray.t
    type coord = int array

    val dft : dir ->
      ?meas:measure -> ?normalize:bool ->
      ?preserve_input:bool -> ?unaligned:bool ->
      ?howmany_n:int array ->
      ?howmanyi: coord list ->
      ?ni: coord -> ?ofsi: coord -> ?inci: coord -> 'l complex_array ->
      ?howmanyo: coord list ->
      ?no: coord -> ?ofso: coord -> ?inco: coord -> 'l complex_array
      -> c2c plan

    val r2c : ?meas:measure -> ?normalize:bool ->
      ?preserve_input:bool -> ?unaligned:bool ->
      ?howmany_n:int array ->
      ?howmanyi: coord list ->
      ?ni: coord -> ?ofsi: coord -> ?inci: coord -> 'l float_array ->
      ?howmanyo: coord list ->
      ?no: coord -> ?ofso: coord -> ?inco: coord -> 'l complex_array
      -> r2c plan

    val c2r : ?meas:measure -> ?normalize:bool ->
      ?preserve_input:bool -> ?unaligned:bool ->
      ?howmany_n:int array ->
      ?howmanyi: coord list ->
      ?ni: coord -> ?ofsi: coord -> ?inci: coord -> 'l complex_array ->
      ?howmanyo: coord list ->
      ?no: coord -> ?ofso: coord -> ?inco: coord -> 'l float_array
      -> c2r plan

    val r2r : r2r_kind array ->
      ?meas:measure -> ?normalize:bool ->
      ?preserve_input:bool -> ?unaligned:bool ->
      ?howmany_n:int array ->
      ?howmanyi: coord list ->
      ?ni: coord -> ?ofsi: coord -> ?inci: coord -> 'l float_array ->
      ?howmanyo: coord list ->
      ?no: coord -> ?ofso: coord -> ?inco: coord -> 'l float_array
      -> r2r plan
  end


  module Array1 :
  sig
    val create: ('a, 'b) kind -> 'c layout -> int -> ('a, 'b, 'c) Array1.t
    val of_array : ('a, 'b) kind -> 'c layout -> 'a array -> ('a, 'b, 'c) Array1.t
    type 'l complex_array = (Complex.t, complex_elt, 'l) Array1.t
    type 'l float_array   = (float, float_elt, 'l) Array1.t

    val dft : dir -> ?meas:measure -> ?normalize:bool ->
      ?preserve_input:bool -> ?unaligned:bool ->
      ?howmany_n:int array ->
      ?howmanyi:int list ->
      ?ni:int -> ?ofsi:int -> ?inci:int -> 'l complex_array ->
      ?howmanyo:int list ->
      ?no:int -> ?ofso:int -> ?inco:int -> 'l complex_array
      -> c2c plan

    val r2r : r2r_kind -> ?meas:measure -> ?normalize:bool ->
      ?preserve_input:bool -> ?unaligned:bool ->
      ?howmany_n:int array ->
      ?howmanyi:int list ->
      ?ni:int -> ?ofsi:int -> ?inci:int -> 'l float_array ->
      ?howmanyo:int list ->
      ?no:int -> ?ofso:int -> ?inco:int -> 'l float_array
      -> r2r plan
  end

  module Array2 :
  sig
    val create: ('a, 'b) kind -> 'c layout -> int -> int -> ('a, 'b, 'c) Array2.t
    type 'l complex_array = (Complex.t, complex_elt, 'l) Array2.t
    type 'l float_array   = (float, float_elt, 'l) Array2.t
    type coord = int * int

    val dft : dir ->
      ?meas:measure -> ?normalize:bool ->
      ?preserve_input:bool -> ?unaligned:bool ->
      ?howmany_n:int array ->
      ?howmanyi: coord list ->
      ?ni: coord -> ?ofsi: coord -> ?inci: coord -> 'l complex_array ->
      ?howmanyo: coord list ->
      ?no: coord -> ?ofso: coord -> ?inco: coord -> 'l complex_array
      -> c2c plan

    val r2r : r2r_kind * r2r_kind -> ?meas:measure -> ?normalize:bool ->
      ?preserve_input:bool -> ?unaligned:bool ->
      ?howmany_n:int array ->
      ?howmanyi: coord list ->
      ?ni: coord -> ?ofsi: coord -> ?inci: coord -> 'l float_array ->
      ?howmanyo: coord list ->
      ?no: coord -> ?ofso: coord -> ?inco: coord -> 'l float_array
      -> r2r plan
  end

  module Array3 :
  sig
    val create: ('a, 'b) kind -> 'c layout -> int -> int -> int
      -> ('a, 'b, 'c) Array3.t

    type 'l complex_array = (Complex.t, complex_elt, 'l) Array3.t
    type 'l float_array   = (float, float_elt, 'l) Array3.t
    type coord = int * int * int

    val dft : dir ->
      ?meas:measure -> ?normalize:bool ->
      ?preserve_input:bool -> ?unaligned:bool ->
      ?howmany_n: int array ->
      ?howmanyi: coord list ->
      ?ni: coord -> ?ofsi: coord -> ?inci: coord -> 'l complex_array ->
      ?howmanyo: coord list ->
      ?no: coord -> ?ofso: coord -> ?inco: coord -> 'l complex_array
      -> c2c plan
  end
end


(** {2 Helper funs}
 ***********************************************************************)

(* specialized for speed *)
let min i j = if (i:int) < j then i else j

module List =
struct
  include List

  let rec list_iteri_loop f i = function
    | [] -> ()
    | a :: tl -> f i a; list_iteri_loop f (succ i) tl

  let iteri ~(f: int -> _ -> unit) l = list_iteri_loop f 0 l
end

let option_map f = function Some v -> Some(f v) | None -> None

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


let is_c_layout m =
  (Genarray.layout m = (Obj.magic c_layout : 'a layout))

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


(** {2 Geometry checks}
 ***********************************************************************)

(* This module perform some checks on the dimensions and howmany
   specifications that depend on the layout but not on the
   precision.  *)
module Geom =
struct
  (* C layout *)
  module C =
  struct
    DEFINE LAYOUT = "c_layout";;
    DEFINE FIRST_INDEX = 0;;
    DEFINE LAST_INDEX(dim) = dim - 1;;
    DEFINE FOR_DIM(k, ndims, expr) = for k = ndims - 1 downto 0 do expr done;;
    DEFINE FOR_HM(k, ndims, expr) = for k = 0 to ndims - 1 do expr done;;
    DEFINE LT_DIM_SYM = "<";;
    DEFINE GE_DIM_SYM = ">=";;
    INCLUDE "fftw3_geom.ml"
  end

  (* FORTRAN layout *)
  module F =
  struct
    DEFINE FORTRAN (* BEWARE it is still defined after this module! *)
    DEFINE LAYOUT = "fortran_layout";;
    DEFINE FIRST_INDEX = 1;;
    DEFINE LAST_INDEX(dim) = dim;;
    DEFINE FOR_DIM(k, ndims, expr) = for k = 0 to ndims - 1 do expr done;;
    DEFINE FOR_HM(k, ndims, expr) = for k = ndims - 1 downto 0 do expr done;;
    DEFINE LT_DIM_SYM = "<=";;
    DEFINE GE_DIM_SYM = ">";;
    INCLUDE "fftw3_geom.ml"
  end

end

(** {2 Precision dependent modules}
 ***********************************************************************)

module D =
struct
  type float_elt = Bigarray.float64_elt
  type complex_elt = Bigarray.complex64_elt
  let float = Bigarray.float64
  let complex = Bigarray.complex64
  ;;
  DEFINE FFTW = "Fftw3.D.";;
  INCLUDE "fftw3SD.ml"
end

IFDEF FFTW3F_EXIST THEN
module S =
struct
  type float_elt = Bigarray.float32_elt
  type complex_elt = Bigarray.complex32_elt
  let float = Bigarray.float32
  let complex = Bigarray.complex32
  ;;
  DEFINE SINGLE_PREC;;
  DEFINE FFTW = "Fftw3.S.";;
  INCLUDE "fftw3SD.ml"
end
ENDIF


module Wisdom =
struct
  external to_file : string -> unit = "fftw3_ocaml_export_wisdom_to_file"

  external to_string : unit -> string = "fftw3_ocaml_export_wisdom_to_string"

  external export : unit -> unit = "fftw3_ocaml_export_wisdom"

  let export write =
    Callback.register "fftw3_wisdom_export" (write: char -> unit);
    export() (* FIXME: is the callback for each char way really efficient??
                It is probably better to export strings to ocaml and have
                an output function -- like Unix.output *)


  external from_system : unit -> unit = "fftw3_ocaml_import_system_wisdom"

  external from_file : string -> unit = "fftw3_ocaml_import_wisdom_from_file"

  external from_string : string -> unit
    = "fftw3_ocaml_import_wisdom_from_string"

  external import : unit -> unit = "fftw3_ocaml_import_wisdom"

  let import read =
    let read_char () = try Char.code(read()) with End_of_file -> -1 (* EOF *) in
    Callback.register "fftw3_wisdom_import" (read_char: unit -> int);
    import()


  external forget : unit -> unit = "fftw3_ocaml_forget_wisdom"

end
