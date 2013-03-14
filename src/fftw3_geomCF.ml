(* File: fftw3_geomCF.ml

   Copyright (C) 2008-

     Christophe Troestler <Christophe.Troestler@umons.ac.be>
     WWW: https://math.umons.ac.be/anum/software/

   This library is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License version 2.1 or
   later as published by the Free Software Foundation, with the special
   exception on linking described in the file LICENSE.

   This library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the file
   LICENSE for more details. *)

(** Geometry checks.  Uniform treatement of the C and FORTRAN layouts
    through macros (for the bigarray interface).  Does not depend on
    the precision of the arrays. *)

open Printf
open Bigarray
open Fftw3_utils

(* Check whether the matrix given by [ofs], [inc], [n] is a valid
   submatrix of [mat].  Return the (C) offset, stride array and
   (logical) dimensions needed by the C wrappers. *)
let get_geom name ofsname ofs incname inc nname n mat =
  let num_dims = Genarray.num_dims mat in
  let ofs = match ofs with
    | Some o ->
        if Array.length o = num_dims then o
        else invalid_arg(sprintf "%s: length %s <> %i" name ofsname num_dims)
    | None -> Array.make num_dims FIRST_INDEX
  and inc = match inc with
    | Some i ->
        if Array.length i = num_dims then i
        else invalid_arg(sprintf "%s: length %s <> %i" name incname num_dims)
    | None -> Array.make num_dims 1
  and n = match n with
    | Some n ->
        if Array.length n = num_dims then n
        else invalid_arg(sprintf "%s: length %s <> %i" name nname num_dims)
    | None -> Array.make num_dims 0 in
  let rank = ref 0 (* Number of dims <> 1 in the transform *)
  and n_sub = Array.make num_dims 1 (* dimensions of the sub-array *)
  and up = Array.make num_dims 0 (* greater submatrix corner *) in
  let set_n_sub v = n_sub.(!rank) <- v; incr rank in
  (* Return the upper bound for the dimension [k].  Set [n_sub] and
     increment [rank] if appropriate. *)
  let upper_bound k abs_inck dimk =
    if n.(k) = 0 then begin
      (* nk = max {n| ofs.(k) + (n-1) * abs inc.(k) <= LAST_INDEX(dimk)} *)
      let nk = 1 + (LAST_INDEX(dimk) - ofs.(k)) / abs_inck in
      if nk > 1 then set_n_sub nk
      else if nk < 1 then
        invalid_arg(sprintf "%s: dim %i empty; no n >= 1 s.t. %i\
	            %+i*(n-1) $LT %i" name k ofs.(k) inc.(k) dimk);
      ofs.(k) + (nk - 1) * abs_inck;
    end
    else if n.(k) > 1 then begin
      (* n.(k) > 1 is given; bound check.   *)
      let last = ofs.(k) + (n.(k) - 1) * abs_inck in
      if last > LAST_INDEX(dimk) then
        invalid_arg(sprintf "%s: %s.(%i) + (%s.(%i) - 1) * abs %s.(%i) \
	                     = %i $GE %i ($LAYOUT) where %s.(%i) = %i \
                             is the physical dim"
                            name ofsname k nname k incname k last dimk
                            nname k n.(k));
      set_n_sub n.(k);
      last;
    end
    else (* n.(k) = 1 => ignore this dimension (do not incr rank) *)
      ofs.(k)
  in
  let pdim = ref 1 (* product of physical dims > k; FORTRAN: < k *)
  and offset = ref 0 (* external functions use the C layout *)
  and stride = Array.make num_dims 0 in
  (* C: decreasing order of k; FORTRAN: increasing order of k *)
  for FOR_DIM(k, num_dims) do
          let dimk = Genarray.nth_dim mat k in
          if n.(k) < 0 then
            invalid_arg(sprintf "%s: %s.(%i) < 0" name nname k);
          if ofs.(k) < FIRST_INDEX then
            invalid_arg(sprintf "%s: %s.(%i) < FIRST_INDEX ($LAYOUT)"
                                name ofsname k);
          stride.(!rank) <- inc.(k) * !pdim;
          if inc.(k) > 0 then (
            up.(k) <- upper_bound k inc.(k) dimk;
            offset := !offset + (ofs.(k) - FIRST_INDEX) * !pdim;
          )
          else if inc.(k) < 0 then (
            up.(k) <- upper_bound k (-inc.(k)) dimk;
            offset := !offset + (up.(k) - FIRST_INDEX) * !pdim;
          )
          else ( (* inc.(k) = 0 => dimension ignored for the transform. *)
            if ofs.(k) > LAST_INDEX(dimk) then
              invalid_arg(sprintf "%s: %s.(%i) = %i $GE %i ($LAYOUT)"
                                  name ofsname k ofs.(k) dimk);
            up.(k) <- ofs.(k);
            offset := !offset + (ofs.(k) - FIRST_INDEX) * !pdim;
          );
          pdim := !pdim * dimk;
  done;
  DEBUG{eprintf "DEBUG: %s: n_sub=%s (rank=%i); offset=%i stride=%s\n%!" name
                (string_of_array n_sub) !rank !offset (string_of_array stride)};
  !offset, (Array.sub n_sub 0 !rank), stride, ofs, up
;;

(** [only_ones d] tells whether [d] entries are all [1] -- which is in
    particular the case if [d] is empty. *)
let only_ones d =
  try
    for i = 0 to Array.length d - 1 do if d.(i) <> 1 then raise Exit done;
    true
  with Exit -> false

(** Check whether the matrix given by [hm_n] (howmany matrix) and [hm]
    is a valid submatrix of the hermitian matrix [mat].  @return the
    [hm_stride], [hm_n] (howmany matrix), [stride] and [n] (logical
    dimensions). *)
let get_geom_hm name hm_nname hm_n hmname hm  nname n low up  mat =
  let num_dims = Genarray.num_dims mat in
  if hm = [] then
    if only_ones hm_n then
      [| |], [| |] (* only one transform *)
    else
      (* Dimensions but no the corresponding vectors *)
      invalid_arg(sprintf "%s: %s = %s but %s = []"
                    name hm_nname (string_of_array hm_n) hmname)
  else begin
    (* Transforms indices = vectors of [hm] with desired dims [hm_n]. *)
    let hm_rank = List.length hm in
    let hm_n =
      if hm_n = [| |] then Array.make hm_rank 0
      else if Array.length hm_n = hm_rank then (
        (* copy [hm_n] because 0 entries will be modified: *)
        let copy_hm i ni =
          if ni >= 0 then ni
          else invalid_arg(sprintf "%s: %s.(%i) < 0" name hm_nname i) in
        Array.mapi copy_hm hm_n
      )
      else invalid_arg(sprintf "%s: length %s = %i <> length %s = %i"
                         name hm_nname (Array.length hm_n) hmname hm_rank) in
    let hm_stride = Array.make hm_rank 0 in
    List.iteri hm ~f:begin fun i v ->
      (* [i]th translation vector [v] *)
      if Array.length v <> num_dims then
        invalid_arg(sprintf "%s: length %ith array of %s <> %i \
			= number of dimensions" name i hmname num_dims);
      let hm_s = ref 0 (* stride corresponding to [v] *) in
      if hm_n.(i) = 0 then (
        (* Dimension for [i]th "howmany vector" [v] to determine *)
        let ni = ref max_int in
        for FOR_HM(k, num_dims) do
               let dimk = Genarray.nth_dim mat k in
               if v.(k) > 0 then
                 (* max{j | up.(k) + v.(k)*(j-1) <= LAST_INDEX(dimk)} *)
                 ni := min !ni (1 + (LAST_INDEX(dimk) - up.(k)) / v.(k))
               else if v.(k) < 0 then
                 (* max{j | low.(k) + v.(k)*(j-1) >= FIRST_INDEX} *)
                 ni := min !ni (1 + (FIRST_INDEX - low.(k)) / v.(k));
               hm_s := !hm_s * dimk + v.(k); (* Horner *)
        done;
        hm_n.(i) <- !ni
      )
      else (
        (* dimension [hm_n.(i)] provided; bound check *)
        for FOR_HM(k, num_dims) do
               let dimk = Genarray.nth_dim mat k in
               if (v.(k) > 0 && up.(k) + v.(k)*(hm_n.(i)-1) > LAST_INDEX(dimk))
                 || (v.(k) < 0 && low.(k) + v.(k)*(hm_n.(i)-1) < FIRST_INDEX)
               then
                 invalid_arg(sprintf "%s: translating %i times by the %ith \
			       vector %s of %s exceeds the %ith dim bounds"
                               name hm_n.(i) i (string_of_array v) hmname k);
               hm_s := !hm_s * dimk + v.(k); (* Horner *)
        done;
      );
      if !hm_s = 0 then
        invalid_arg(sprintf "%s: %ith element of %s = [|0.;...;0.|]"
                      name i hmname);
      hm_stride.(i) <- !hm_s;
    end;
    DEBUG{eprintf "DEBUG: %s: hm_n=%s; hm_stride=%s\n%!" name
                  (string_of_array hm_n) (string_of_array hm_stride)};
    hm_n, hm_stride
  end


(* Take the [make_plan] function creating plans, the dimensions, offsets
   and increments of input/output arrays and compute the informations
   needed by [make_plan].  Check the coherence of the data at the same
   time.  There may be more input (resp. output) arrays than [i]
   (resp. [o]) but these must have the same dimensions. *)
let apply name make_plan hm_n  hmi ?ni ofsi inci i  hmo ?no ofso inco o
    ~logical_dims =
  let num_dims = Genarray.num_dims i in
  if num_dims <> Genarray.num_dims o then
    invalid_arg(name ^ ": input and output arrays do not have the same \
	NUMBER of dimensions");
  let offseti, ni, stridei, lowi, upi =
    get_geom name "ofsi" ofsi "inci" inci "n" ni i
  and offseto, no, strideo, lowo, upo =
    get_geom name "ofso" ofso "inco" inco "n" no o in
  let n =                               (* or raise invalid_arg *)
    logical_dims ni no
      (sprintf "%s: dim input = %s incompatible with dim ouput = %s"
         name (string_of_array ni) (string_of_array no)) in
  let hm_ni, hm_stridei =
    get_geom_hm name "howmany_n" hm_n "howmanyi" hmi  "n" ni lowi upi i
  and hm_no, hm_strideo =
    get_geom_hm name "howmany_n" hm_n "howmanyo" hmo  "n" no lowo upo o  in
  if hm_ni <> hm_no then
    invalid_arg(sprintf "%s: howmany dim input = %s <> howmany dim output = %s"
                  name (string_of_array hm_ni) (string_of_array hm_no));
  make_plan offseti offseto n stridei strideo  hm_ni hm_stridei hm_strideo
