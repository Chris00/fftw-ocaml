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


module D =
struct
  type float_elt = Bigarray.float64_elt
  type complex_elt = Bigarray.complex64_elt
  let float = float64
  let complex = complex64

  INCLUDE "fftw3SD.ml"
end

module S =
struct
  type float_elt = Bigarray.float32_elt
  type complex_elt = Bigarray.complex32_elt
  let float = float32
  let complex = complex32

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
