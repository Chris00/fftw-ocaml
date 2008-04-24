(*pp camlp4o pa_macro.cmo *)
(* File: fftw3.ml

   Objective Caml interface for FFTW.

   Copyright (C) 2005

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

  let iteri (f: int -> _ -> unit) l = list_iteri_loop f 0 l
end

(** positive part *)
let pos i = if i > 0 then i else 0

(** negative part *)
let neg i = if i < 0 then i else 0

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


(** {2 Precision dependent modules}
 ***********************************************************************)

module D =
struct
  type float_elt = Bigarray.float64_elt
  type complex_elt = Bigarray.complex64_elt
  let float = Bigarray.float64
  let complex = Bigarray.complex64

  INCLUDE "fftw3SD.ml"
end

module S =
struct
  type float_elt = Bigarray.float32_elt
  type complex_elt = Bigarray.complex32_elt
  let float = Bigarray.float32
  let complex = Bigarray.complex32

  DEFINE SINGLE_PREC
  INCLUDE "fftw3SD.ml"
end



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
