(* File: fftw3SD.ml

   Copyright (C) 2006

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
  i : 'b 'c 'd. ('b,'c,'d) Genarray.t; (* input array; => not freed by GC *)
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
let flags meas unaligned preserve_input =
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


(** {2 Helper funs}
 ***********************************************************************)

let min i j = if (i:int) < j then i else j

module List =
struct
  include List

  let rec list_iteri_loop f i = function
    | [] -> ()
    | a :: tl -> f i a; list_iteri_loop f (succ i) tl

  let iteri f l = list_iteri_loop f 0 l
end

(* positive part *)
let pos i = if i > 0 then i else 0

(* negative part *)
let neg i = if i < 0 then i else 0


(** {2 Creating plans}
 ***********************************************************************)

module Genarray = struct
  external create: ('a, 'b) Bigarray.kind -> 'c Bigarray.layout ->
    int array -> ('a, 'b, 'c) Bigarray.Genarray.t
    = "fftw3_ocaml_ba_create"

  type 'l complex_array = (Complex.t, complex_elt, 'l) Genarray.t
  type 'l float_array   = (float, float_elt, 'l) Genarray.t

  let is_c_layout m = (Genarray.layout m = (Obj.magic c_layout : 'a layout))

  (* [get_rank default m] returns the rank provided by the first
     matrix in the list of options [m]. *)
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

  (* Assume that [a] has at least one element -- which does not matter
     for here because bigarrays have at least one dim. *)
  let string_of_array a =
    let b = Buffer.create 80 in
    Buffer.add_string b "[|";
    Buffer.add_string b (string_of_int a.(0));
    for i = 1 to Array.length a - 1 do
      Buffer.add_string b "; ";
      Buffer.add_string b (string_of_int a.(i));
    done;
    Buffer.add_string b "|]";
    Buffer.contents b

  (* Check whether the matrix given by [ofs], [inc], [n] is a valid
     submatrix of the one whose dimensions are given by [dim].  Return
     the offset, stride and dimensions needed by the C wrappers. *)
  let get_geom_c name ofsname ofs incname inc nname n mat =
    let rank = Genarray.num_dims mat in
    let ofs = match ofs with
      | Some o ->
          if Array.length o = rank then o
          else invalid_arg(sprintf "%s: length %s <> %i" name ofsname rank)
      | None -> Array.make rank 0
    and inc = match inc with
      | Some i ->
          if Array.length i = rank then i
          else invalid_arg(sprintf "%s: length %s <> %i" name incname rank)
      | None -> Array.make rank 1
    and n = match n with
      | Some n ->
          if Array.length n = rank then
            Array.copy n (* will be modified when a value is 0 *)
          else invalid_arg(sprintf "%s: length %s <> %i" name nname rank)
      | None -> Array.make rank 0 in
    let offset = ref 0
    and stride = Array.make rank 0
    and pdim = ref 1 (* product of dims > k *) in
    for k = rank - 1 downto 0 do
      let dimk = Genarray.nth_dim mat k in
      stride.(k) <- inc.(k) * !pdim;
      offset := !offset + ofs.(k) * !pdim;
      if n.(k) < 0 then invalid_arg(sprintf "%s: %s.(%i) < 0" name nname k);
      if inc.(k) > 0 then (
        if ofs.(k) < 0 then
          invalid_arg(sprintf "%s: %s.(%i) < 0" name ofsname k);
        if n.(k) = 0 then (
          (* n.(k) = max { n : ofs.(k) + (n-1) * inc.(k) < dimk } *)
          n.(k) <- 1 + (dimk - ofs.(k) - 1) / inc.(k);
          if n.(k) < 1 then
            invalid_arg(sprintf "%s: dim %i empty; no n >= 1 s.t. \
		%i%+i*(n-1) < %i" name k ofs.(k) inc.(k) dimk);
        )
        else
          (* n.(k) > 0 is given; bound check *)
          let last = ofs.(k) + (n.(k) - 1) * inc.(k) in
          if last >= dimk then
            invalid_arg(sprintf "%s: %s.(%i) + (%s.(%i) - 1) * %s.(%i) = %i \
		>= %i where %s.(%i) = %i is the physical dim"
                           name ofsname k nname k incname k last dimk
                           nname k n.(k));
      )
      else if inc.(k) < 0 then (
        if ofs.(k) >= dimk then
          invalid_arg(sprintf "%s: %s.(%i) >= %i" name ofsname k dimk);
        if n.(k) = 0 then (
          (* n.(k) = max { n : 0 <= ofs.(k) + (n-1) * inc.(k) } *)
          n.(k) <- 1 - ofs.(k) / inc.(k);
          if n.(k) < 1 then
            invalid_arg(sprintf "%s: dim %i empty; no n >= 1 s.t. %i%+i*(n-1) \
		>= 0" name k ofs.(k) inc.(k));
        )
        else
          (* n.(k) > 0 is given; bound check *)
          let last = ofs.(k) + (n.(k) - 1) * inc.(k) in
          if last < 0 then
            invalid_arg(sprintf "%s: %s.(%i) + (%s.(%i) - 1) * %s.(%i) = %i \
		< 0 where %s.(%i) = %i is the physical dim"
                           name ofsname k nname k incname k last nname k n.(k));
      )
      else (* inc.(k) = 0 *)
        invalid_arg(sprintf "%s: %s.(%i) = 0" name incname k);
      pdim := !pdim * dimk;
    done;
    !offset, stride, ofs, inc, n

  let get_geom_hm_c name hm_nname hm_n hmname hm stride ofs inc nname n mat =
    let rank = Genarray.num_dims mat in
    if hm = [] then
      if hm_n = [| |] then
        [| |], [| |], stride, n (* only one transform *)
      else (
        (* Transforms indexed by the first dimensions of the submatrix
           of [mat] defined by [ofs] and [inc]. *)
        let hm_rank = Array.length hm_n in
        if hm_rank >= rank then
          invalid_arg(sprintf "%s: length %s >= %i = number of dimensions"
                         name hm_nname rank);
        let mat_rank = rank - hm_rank in
        let hm_stride = Array.sub stride 0 hm_rank in
        let stride = Array.sub stride hm_rank mat_rank in
        let hm_n = Array.mapi (fun i ni ->
          if ni < 0 then invalid_arg(sprintf "%s: %s.(%i) < 0" name hm_nname i)
          else if ni = 0 then n.(i)
          else if ni <= n.(i) then ni
          else invalid_arg(sprintf "%s: %s.(%i) = %i > %s.(%i) = %i"
                              name hm_nname i ni nname i n.(i))
        ) hm_n in
        let n = Array.sub n hm_rank mat_rank in
        hm_stride, hm_n, stride, n
      )
    else (
      (* Transforms indices = vectors of [hm] with desired dims [hm_n]. *)
      let hm_rank = List.length hm in
      let hm_n =
        if hm_n = [| |] then Array.make hm_rank 0
        else if Array.length hm_n = hm_rank then (
          Array.mapi (fun i ni ->
            if ni < 0 then
              invalid_arg(sprintf "%s: %s.(%i) < 0" name hm_nname i);
            ni
          ) hm_n; (* copy because 0 entries will be modified *)
        )
        else invalid_arg(sprintf "%s: length %s = %i <> length %s = %i"
                            name hm_nname (Array.length hm_n) hmname hm_rank) in
      let hm_stride = Array.make hm_rank 0 in
      List.iteri begin fun i v ->
        (* [i]th translation vector [v] *)
        if Array.length v <> rank then
          invalid_arg(sprintf "%s: length %ith element of %s <> %i \
		= number of dimensions" name i hmname rank);
        let hm_s = ref 0 in
        if hm_n.(i) = 0 then (
          (* dimension for [i]th "howmany vector" [v] to determine *)
          let ni = ref max_int in
          for k = 0 to rank - 1 do
            let dimk = Genarray.nth_dim mat k in
            if v.(k) > 0 then
              (* max{j : ofs.(k) + inc.(k)^+ *(n.(k)-1) + v.(k)*(j-1) < dimk} *)
              ni := min !ni
                (1 + (dimk - ofs.(k) - pos inc.(k) * (n.(k)-1) - 1) / v.(k))
            else if v.(k) < 0 then
              (* max{j : ofs.(k) + inc.(k)^- *(n.(k)-1) + v.(k)*(j-1) >= 0} *)
              ni := min !ni (1 - (ofs.(k) + neg inc.(k) * (n.(k)-1)) / v.(k));
            hm_s := !hm_s * dimk + v.(k); (* Horner *)
          done;
          hm_n.(i) <- !ni
        )
        else (
          (* dimension [hm_n.(k)] provided; bound check *)
          for k = 0 to rank - 1 do
            let dimk = Genarray.nth_dim mat k in
            if (v.(k) > 0 &&
                   ofs.(k) + pos inc.(k)*(n.(k)-1) + v.(k)*(hm_n.(i)-1) >= dimk)
              || (v.(k) < 0 &&
                     ofs.(k) + neg inc.(k) * (n.(k)-1) + v.(k)*(hm_n.(i)-1) < 0)
            then
              invalid_arg(sprintf "%s: translating %i times the %ith \
			vector of %s exceeds the %ith dim array bounds"
                             name hm_n.(i) i hmname k);
            hm_s := !hm_s * dimk + v.(k); (* Horner *)
          done;
        );
        if !hm_s = 0 then
          invalid_arg(sprintf "%s: %ith element of %s = [|0.;...;0.|]"
                         name i hmname);
        hm_stride.(i) <- !hm_s;
      end hm;
      hm_stride, hm_n, stride, n
    )

  (* Take the [wrapper], the dimensions, offsets and increments of
     input/output arrays and compute the informations needed by the
     wrapper.  Check the coherence of the data at the same time.
     There may be more input (resp. output) arrays than [i]
     (resp. [o]) but these must have the same dimensions. *)
  let apply_c name wrapper n hm_n  hmi ofsi inci i  hmo ofso inco o =
    let rank = Genarray.num_dims i in
    if rank <> Genarray.num_dims o then
      invalid_arg(name ^ ": input and output arrays do not have the same \
	NUMBER of dimensions");
    let offseti, stridei, ofsi, inci, ni =
      get_geom_c name "ofsi" ofsi "inci" inci "n" n i
    and offseto, strideo, ofso, inco, no =
      get_geom_c name "ofso" ofso "inco" inco "n" n o in
    if ni <> no then invalid_arg("dim input = " ^ string_of_array ni
                                  ^ " <> dim ouput = " ^ string_of_array no);
    let hm_stridei, hm_ni, stridei, ni =
      get_geom_hm_c name "howmany_n" hm_n "howmanyi" hmi stridei
        ofsi inci "n" ni i
    and hm_strideo, hm_no, strideo, no =
      get_geom_hm_c name "howmany_n" hm_n "howmanyo" hmo strideo
        ofso inco "n" no o in
    if hm_ni <> hm_no then
      invalid_arg("howmany dim input = " ^ string_of_array hm_ni
                   ^ " <> howmany dim output = " ^ string_of_array hm_no);
    wrapper offseti offseto ni stridei strideo hm_ni hm_stridei hm_strideo


  (* Like get_geom_c but for fortran layout *)
  let get_geom_f name ofsname ofs incname inc nname n mat =
    let rank = Genarray.num_dims mat in
    let ofs = match ofs with
      | Some o ->
          if Array.length o = rank then o
          else invalid_arg(sprintf "%s: length %s <> %i" name ofsname rank)
      | None -> Array.make rank 1 (* FORTRAN *)
    and inc = match inc with
      | Some i ->
          if Array.length i = rank then i
          else invalid_arg(sprintf "%s: length %s <> %i" name incname rank)
      | None -> Array.make rank 1
    and n = match n with
      | Some n ->
          if Array.length n = rank then
            Array.copy n (* will be modified when a value is 0 *)
          else invalid_arg(sprintf "%s: length %s <> %i" name nname rank)
      | None -> Array.make rank 0 in
    let offset = ref 0 (* FORTRAN: external functions: C layout *)
    and stride = Array.make rank 0
    and pdim = ref 1 (* FORTRAN: product of dims < k *) in
    for k = 0 to rank - 1 do
      let dimk = Genarray.nth_dim mat k in
      stride.(k) <- inc.(k) * !pdim;
      offset := !offset + (ofs.(k) - 1) * !pdim;
      if n.(k) < 0 then invalid_arg(sprintf "%s: %s.(%i) < 0" name nname k);
      if inc.(k) > 0 then (
        if ofs.(k) < 1 (* FORTRAN *) then
          invalid_arg(sprintf "%s: %s.(%i) < 1 (fortran_layout)"
                         name ofsname k);
        if n.(k) = 0 then (
          (* n.(k) = max { n : ofs.(k) + (n-1) * inc.(k) <= dimk } *)
          n.(k) <- 1 + (dimk - ofs.(k)) / inc.(k);
          if n.(k) < 1 then
            invalid_arg(sprintf "%s: dim %i empty; no n >= 1 s.t. \
		%i%+i*(n-1) <= %i" name k ofs.(k) inc.(k) dimk);
        )
        else
          (* n.(k) > 0 is given; bound check *)
          let last = ofs.(k) + (n.(k) - 1) * inc.(k) in
          if last > dimk then
            invalid_arg(sprintf "%s: %s.(%i) + (%s.(%i) - 1) * %s.(%i) = %i \
		> %i" name ofsname k nname k incname k last dimk);
      )
      else if inc.(k) < 0 then (
        if ofs.(k) > dimk then
          invalid_arg(sprintf "%s: %s.(%i) > %i" name ofsname k dimk);
        if n.(k) = 0 then (
          (* n.(k) = max { n : 1 <= ofs.(k) + (n-1) * inc.(k) } *)
          n.(k) <- 1 + (1 - ofs.(k)) / inc.(k);
          if n.(k) < 1 then
            invalid_arg(sprintf "%s: dim %i empty; no n >= 1 s.t. %i%+i*(n-1) \
		>= 1 (fortran_layout)" name k ofs.(k) inc.(k));
        )
        else
          (* n.(k) > 0 is given; bound check *)
          let last = ofs.(k) + (n.(k) - 1) * inc.(k) in
          if last < 1 then
            invalid_arg(sprintf "%s: %s.(%i) + (%s.(%i) - 1) * %s.(%i) = %i \
		< 1 (fortran_layout)" name ofsname k nname k incname k last);
      )
      else (* inc.(k) = 0 *)
        invalid_arg(sprintf "%s: %s.(%i) = 0" name incname k);
      pdim := !pdim * dimk;
    done;
    !offset, stride, ofs, inc, n

  (* Like get_geom_hm_c but for fortran layout *)
  let get_geom_hm_f name hm_nname hm_n hmname hm stride ofs inc nname n mat =
    let rank = Genarray.num_dims mat in
    if hm = [] then
      if hm_n = [| |] then
        [| |], [| |], stride, n (* only one transform *)
      else (
        (* Transforms indexed by the last (FORTRAN) dimensions of the
           submatrix of [mat] defined by [ofs] and [inc]. *)
        let hm_rank = Array.length hm_n in
        if hm_rank >= rank then
          invalid_arg(sprintf "%s: length %s >= %i = number of dimensions"
                         name hm_nname rank);
        let mat_rank = rank - hm_rank in
        let hm_stride = Array.sub stride mat_rank hm_rank in
        let stride = Array.sub stride 0 mat_rank in
        let hm_n = Array.mapi (fun i ni ->
          let j = mat_rank + i in
          if ni < 0 then invalid_arg(sprintf "%s: %s.(%i) < 0" name hm_nname i)
          else if ni = 0 then n.(j)
          else if ni <= n.(j) then ni
          else invalid_arg(sprintf "%s: %s.(%i) = %i > %s.(%i) = %i"
                              name hm_nname i ni nname j n.(j))
        ) hm_n in
        let n = Array.sub n 0 mat_rank in
        hm_stride, hm_n, stride, n
      )
    else (
      (* Transforms indices = vectors of [hm] with desired dims [hm_n]. *)
      let hm_rank = List.length hm in
      let hm_n =
        if hm_n = [| |] then Array.make hm_rank 0
        else if Array.length hm_n = hm_rank then (
          Array.mapi (fun i ni ->
            if ni < 0 then
              invalid_arg(sprintf "%s: %s.(%i) < 0" name hm_nname i);
            ni
          ) hm_n; (* copy because 0 entries will be modified *)
        )
        else invalid_arg(sprintf "%s: length %s = %i <> length %s = %i"
                            name hm_nname (Array.length hm_n) hmname hm_rank) in
      let hm_stride = Array.make hm_rank 0 in
      List.iteri begin fun i v ->
        (* [i]th translation vector [v] *)
        if Array.length v <> rank then
          invalid_arg(sprintf "%s: length %ith element of %s <> %i \
		= number of dimensions" name i hmname rank);
        let hm_s = ref 0 in
        if hm_n.(i) = 0 then (
          (* dimension for [i]th "howmany vector" [v] to determine *)
          let ni = ref max_int in
          for k = rank - 1 downto 0 do (* FORTRAN *)
            let dimk = Genarray.nth_dim mat k in
            if v.(k) > 0 then
              (* max{j: ofs.(k) + inc.(k)^+ *(n.(k)-1) + v.(k)*(j-1) <= dimk} *)
              ni := min !ni
                (1 + (dimk - ofs.(k) - pos inc.(k) * (n.(k)-1)) / v.(k))
            else if v.(k) < 0 then
              (* max{j: ofs.(k) + inc.(k)^- *(n.(k)-1) + v.(k)*(j-1) >= 1} *)
              ni := min !ni
                (1 + (1 - ofs.(k) - neg inc.(k) * (n.(k)-1)) / v.(k));
            hm_s := !hm_s * dimk + v.(k); (* Horner *)
          done;
          hm_n.(i) <- !ni
        )
        else (
          (* dimension [hm_n.(k)] provided; bound check *)
          for k = rank - 1 downto 0 do (* FORTRAN *)
            let dimk = Genarray.nth_dim mat k in
            if (v.(k) > 0 &&
                   ofs.(k) + pos inc.(k)*(n.(k)-1) + v.(k)*(hm_n.(i)-1) > dimk)
              || (v.(k) < 0 &&
                     ofs.(k) + neg inc.(k) * (n.(k)-1) + v.(k)*(hm_n.(i)-1) < 1)
            then
              invalid_arg(sprintf "%s: translating %i times the %ith \
			vector of %s exceeds the %ith dim array bounds"
                             name hm_n.(i) i hmname k);
            hm_s := !hm_s * dimk + v.(k); (* Horner *)
          done;
        );
        if !hm_s = 0 then
          invalid_arg(sprintf "%s: %ith element of %s = [|0.;...;0.|]"
                         name i hmname);
        hm_stride.(i) <- !hm_s;
      end hm;
      hm_stride, hm_n, stride, n
    )

  let apply_f name wrapper n hm_n  hmi ofsi inci i  hmo ofso inco o =
    let rank = Genarray.num_dims i in
    if rank <> Genarray.num_dims o then
      invalid_arg(name ^ ": input and output arrays do not have the same \
	NUMBER of dimensions");
    let offseti, stridei, ofsi, inci, ni =
      get_geom_f name "ofsi" ofsi "inci" inci "n" n i
    and offseto, strideo, ofso, inco, no =
      get_geom_f name "ofso" ofso "inco" inco "n" n o in
    if ni <> no then invalid_arg("dim input = " ^ string_of_array ni
                                  ^ " <> dim ouput = " ^ string_of_array no);
    let hm_stridei, hm_ni, stridei, ni =
      get_geom_hm_f name "howmany_n" hm_n "howmanyi" hmi stridei
        ofsi inci "n" ni i
    and hm_strideo, hm_no, strideo, no =
      get_geom_hm_f name "howmany_n" hm_n "howmanyo" hmo strideo
        ofso inco "n" no o in
    if hm_ni <> hm_no then
      invalid_arg("howmany dim input = " ^ string_of_array hm_ni
                   ^ " <> howmany dim output = " ^ string_of_array hm_no);
    wrapper offseti offseto ni stridei strideo hm_ni hm_stridei hm_strideo


  let apply name wrapper  n hm_ranki hmi ofsi inci i hm_ranko hmo ofso inco o =
    (if is_c_layout i then apply_c else apply_f)
      name wrapper n hm_ranki hmi ofsi inci i hm_ranko hmo ofso inco o


  let dft dir ?(meas=Measure) ?(normalize=false)
      ?(preserve_input=false) ?(unaligned=false) ?n
      ?howmany_ranki ?howmanyi  ?ofsi ?inci (i: 'l complex_array)
      ?howmany_ranko ?howmanyo  ?ofso ?inco (o: 'l complex_array) =
    apply "Fftw3.D.Genarray.dft"
      (guru_dft i o (sign_of_dir dir) (flags meas unaligned preserve_input))
      howmany_rank howmany n ofsi inci i ofso inco o

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
