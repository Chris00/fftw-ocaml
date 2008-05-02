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

type genarray
external genarray : (_,_,_) Genarray.t -> genarray = "%identity"
    (* Since we want the FFT functions to be polymorphic in the layout
       of the arrays, some back magic is unavoidable.  This one way
       conversion is actually safe, it is the use of [genarray] by C
       functions that must be taken care of. *)

type 'a plan = {
  plan: 'a fftw_plan;
  i : genarray; (* hold input array => not freed by GC before the plan *)
  offseto : int; (* output offset; C-stubs *)
  strideo : int array; (* strides; C-stubs *)
  no : int array; (* dimensions *)
  o : genarray; (* output array *)
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
(*   if p.normalize then *)
(*     normalize p.o p.offseto p.strideo p.no p.normalize_factor *)
;;

module Guru = struct

  let dft plan i o =
  (* how to check that the arrays conform to the plan specification? *)
  exec_dft plan i o

  let split_dft plan ri ii ro io =
    (* again, how to check conformance with the plan? *)
    exec_split_dft plan ri ii ro io

end

(** {2 Creating plans}
 ***********************************************************************)

module Genarray = struct
  external create: ('a, 'b) Bigarray.kind -> 'c Bigarray.layout ->
    int array -> ('a, 'b, 'c) Bigarray.Genarray.t
    = "fftw3_ocaml_ba_create"

  type 'l complex_array = (Complex.t, complex_elt, 'l) Genarray.t
  type 'l float_array   = (float, float_elt, 'l) Genarray.t


  (* Layout independent function *)
  let apply name make_plan n hm_n  hmi ofsi inci i  hmo ofso inco o  normalize =
    let make offseti offseto ni stridei strideo hm_ni hm_stridei hm_strideo =
      let p = make_plan offseti offseto ni stridei strideo
        hm_ni hm_stridei hm_strideo in
      { plan = p;
        i = genarray i;
        offseto = offseto;
        strideo = strideo;
        no = ni;
        o = genarray o;
        normalize = normalize;
        normalize_factor = 1.;          (* FIXME:  *)
      } in
    (if is_c_layout i then Geom.C.apply else Geom.F.apply) name make
      n hm_n  hmi ofsi inci i  hmo ofso inco o

  let dft_name = FFTW ^ "Genarray.dft"
  let dft dir ?(meas=Measure) ?(normalize=false)
      ?(preserve_input=false) ?(unaligned=false) ?n ?(howmany_n=[| |])
      ?(howmanyi=[]) ?ofsi ?inci (i: 'l complex_array)
      ?(howmanyo=[]) ?ofso ?inco (o: 'l complex_array) =
    apply dft_name
      (guru_dft i o (sign_of_dir dir) (flags meas unaligned preserve_input))
      n howmany_n  howmanyi ofsi inci i  howmanyo ofso inco o  normalize

  (* At the moment, in place transforms are not possible but they may
     be if OCaml bug 0004333 is resolved. *)
  let r2c_name = FFTW ^ "Genarray.r2c"
  let r2c ?(meas=Measure) ?(normalize=false)
      ?(preserve_input=false) ?(unaligned=false) ?n ?(howmany_n=[| |])
      ?(howmanyi=[]) ?ofsi ?inci (i: 'l float_array)
      ?(howmanyo=[]) ?ofso ?inco (o: 'l complex_array) =
    apply r2c_name
      (guru_r2c i o (flags meas unaligned preserve_input))
      n howmany_n  howmanyi ofsi inci i  howmanyo ofso inco o  normalize

  let c2r_name = FFTW ^ "Genarray.c2r"
  let c2r ?(meas=Measure) ?(normalize=false)
      ?(preserve_input=false) ?(unaligned=false) ?n ?(howmany_n=[| |])
      ?(howmanyi=[]) ?ofsi ?inci (i: 'l complex_array)
      ?(howmanyo=[]) ?ofso ?inco (o: 'l float_array) =
    (* FIXME:  Are the checks of apply appropriate here? *)
    apply c2r_name
      (guru_c2r i o (flags meas unaligned preserve_input))
      n howmany_n  howmanyi ofsi inci i  howmanyo ofso inco o  normalize

  let r2r_name = FFTW ^ "Genarray.r2r"
  let r2r kind ?(meas=Measure) ?(normalize=false)
      ?(preserve_input=true) ?(unaligned=false) ?n ?(howmany_n=[| |])
      ?(howmanyi=[]) ?ofsi ?inci (i: 'l float_array)
      ?(howmanyo=[]) ?ofso ?inco (o: 'l float_array) =
    (* FIXME: must check [kind] has the right length?? *)
    let kind = Array.map int_of_r2r_kind kind in
    apply r2r_name
      (guru_r2r i o kind (flags meas unaligned preserve_input))
      n howmany_n  howmanyi ofsi inci i  howmanyo ofso inco o  normalize
(*       (\* howmany_rank: *\)(Some(Genarray.num_dims i - Array.length kind)) *)


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


  let apply name make_plan n hm_n  hmi ofsi inci i  hmo ofso inco o  normalize =
    let n = option_map (fun n -> [| n |]) n in
    let hmi = List.map (fun v -> [| v |]) hmi in
    let ofsi = option_map (fun n -> [| n |]) ofsi in
    let inci = Some [| inci |] in
    let hmo = List.map (fun v -> [| v |]) hmo in
    let ofso = option_map (fun n -> [| n |]) ofso in
    let inco = Some [| inco |] in
    Genarray.apply name make_plan
      n hm_n  hmi ofsi inci i  hmo ofso inco o  normalize

  let dft_name = FFTW ^ "Array1.dft"
  let dft dir ?(meas=Measure) ?(normalize=false)
      ?(preserve_input=true) ?(unaligned=false) ?n ?(howmany_n=[| |])
      ?(howmanyi=[]) ?ofsi ?(inci=1) (i: 'l complex_array)
      ?(howmanyo=[]) ?ofso ?(inco=1) (o: 'l complex_array) =
    let gi = genarray_of_array1 i
    and go = genarray_of_array1 o in
    apply dft_name
      (guru_dft gi go (sign_of_dir dir) (flags meas unaligned preserve_input))
      n howmany_n  howmanyi ofsi inci gi  howmanyo ofso inco go  normalize

  let r2r_name = FFTW ^ "Array1.r2r"
  let r2r kind ?(meas=Measure) ?(normalize=false)
      ?(preserve_input=true) ?(unaligned=false) ?n ?(howmany_n=[| |])
      ?(howmanyi=[]) ?ofsi ?(inci=1) (i: 'l float_array)
      ?(howmanyo=[]) ?ofso ?(inco=1) (o: 'l float_array) =
    let gi = genarray_of_array1 i
    and go = genarray_of_array1 o in
    let kind = [| int_of_r2r_kind kind |] in
    apply r2r_name
      (guru_r2r gi go kind (flags meas unaligned preserve_input))
      n howmany_n  howmanyi ofsi inci gi howmanyo ofso inco go  normalize
end


module Array2 = struct
  external array2_of_ba : ('a,'b,'c) Bigarray.Genarray.t -> ('a,'b,'c) Array2.t
    = "%identity"
    (* BEWARE: only for bigarray with 2D, convert without check *)

  let create kind layout dim1 dim2 =
    array2_of_ba(Genarray.create kind layout [|dim1; dim2|])

  type 'l complex_array = (Complex.t, complex_elt, 'l) Array2.t
  type 'l float_array   = (float, float_elt, 'l) Array2.t

  let apply name make_plan n hm_n  hmi ofsi (inci1,inci2) i
      hmo ofso (inco1,inco2) o  normalize =
    let n = option_map (fun (n1,n2) -> [| n1; n2 |]) n in
    let hmi = List.map (fun (d1,d2) -> [| d1; d2 |]) hmi in
    let ofsi = option_map (fun (n1,n2) -> [| n1; n2 |]) ofsi in
    let inci = Some [| inci1; inci2 |] in
    let hmo = List.map (fun (d1,d2) -> [| d1; d2 |]) hmo in
    let ofso = option_map (fun (n1,n2) -> [| n1; n2 |]) ofso in
    let inco = Some [| inco1; inco2 |] in
    Genarray.apply name make_plan
      n hm_n  hmi ofsi inci i  hmo ofso inco o  normalize

  let dft_name = FFTW ^ "Array2.dft"
  let dft dir ?(meas=Measure) ?(normalize=false)
      ?(preserve_input=true) ?(unaligned=false)
      ?n ?(howmany_n=[| |])
      ?(howmanyi=[]) ?ofsi ?(inci=(1,1)) (i: 'l complex_array)
      ?(howmanyo=[]) ?ofso ?(inco=(1,1)) (o: 'l complex_array) =
    let gi = genarray_of_array2 i
    and go = genarray_of_array2 o in
    apply dft_name
      (guru_dft gi go (sign_of_dir dir) (flags meas unaligned preserve_input))
      n howmany_n  howmanyi ofsi inci gi  howmanyo ofso inco go  normalize

  let r2r_name = FFTW ^ "Array2.r2r"
  let r2r (kind1,kind2) ?(meas=Measure) ?(normalize=false)
      ?(preserve_input=true) ?(unaligned=false) ?n ?(howmany_n=[| |])
      ?(howmanyi=[]) ?ofsi ?(inci=(1,1)) (i: 'l float_array)
      ?(howmanyo=[]) ?ofso ?(inco=(1,1)) (o: 'l float_array) =
    let gi = genarray_of_array2 i
    and go = genarray_of_array2 o in
    let kind = [| int_of_r2r_kind kind1; int_of_r2r_kind kind2 |] in
    apply r2r_name
      (guru_r2r gi go kind (flags meas unaligned preserve_input))
      n howmany_n  howmanyi ofsi inci gi howmanyo ofso inco go  normalize
end


module Array3 = struct
  external array3_of_ba : ('a,'b,'c) Bigarray.Genarray.t -> ('a,'b,'c) Array3.t
    = "%identity"
    (* BEWARE: only for bigarray with 3D, convert without check *)

  let create kind layout dim1 dim2 dim3 =
    array3_of_ba(Genarray.create kind layout [|dim1; dim2; dim3|])

  type 'l complex_array = (Complex.t, complex_elt, 'l) Array3.t
  type 'l float_array   = (float, float_elt, 'l) Array3.t


  let apply name make_plan n hm_n  hmi ofsi (inci1,inci2,inci3) i
      hmo ofso (inco1,inco2,inco3) o  normalize =
    let n = option_map (fun (n1,n2,n3) -> [| n1; n2; n3 |]) n in
    let hmi = List.map (fun (d1,d2,d3) -> [| d1; d2; d3 |]) hmi in
    let ofsi = option_map (fun (n1,n2,n3) -> [| n1; n2; n3 |]) ofsi in
    let inci = Some [| inci1; inci2; inci3 |] in
    let hmo = List.map (fun (d1,d2,d3) -> [| d1; d2; d3 |]) hmo in
    let ofso = option_map (fun (n1,n2,n3) -> [| n1; n2; n3 |]) ofso in
    let inco = Some [| inco1; inco2; inco3 |] in
    Genarray.apply name make_plan
      n hm_n  hmi ofsi inci i  hmo ofso inco o  normalize

  let dft_name = FFTW ^ "Array3.dft"
  let dft dir ?(meas=Measure) ?(normalize=false)
      ?(preserve_input=true) ?(unaligned=false)
      ?n ?(howmany_n=[| |])
      ?(howmanyi=[]) ?ofsi ?(inci=(1,1,1)) (i: 'l complex_array)
      ?(howmanyo=[]) ?ofso ?(inco=(1,1,1)) (o: 'l complex_array) =
    let gi = genarray_of_array3 i
    and go = genarray_of_array3 o in
    apply dft_name
      (guru_dft gi go (sign_of_dir dir) (flags meas unaligned preserve_input))
      n howmany_n  howmanyi ofsi inci gi  howmanyo ofso inco go  normalize

end
