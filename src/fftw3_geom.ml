(* File: fftw3_geom.ml

   Copyright (C) 2008

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

(** Geometry checks.  Uniform treatement of the C and FORTRAN layouts
    through macros (for the bigarray interface).  Does not depend on
    the precision of the arrays. *)


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
  and offset = ref 0 (* external functions use the C layout *)
  and stride = Array.make num_dims 0
  and pdim = ref 1 (* product of physical dims > k; FORTRAN: < k *)
  and low = Array.make num_dims 0 (* lower submatrix corner *)
  and up = Array.make num_dims 0 (* greater submatrix corner *) in
  let set_n_sub v = n_sub.(!rank) <- v; incr rank in
  (* C: decreasing order of k; FORTRAN: increasing order of k *)
  FOR_DIM(k, num_dims,
          let dimk = Genarray.nth_dim mat k in
          stride.(!rank) <- inc.(k) * !pdim;
          offset := !offset + (ofs.(k) - FIRST_INDEX) * !pdim;
          if n.(k) < 0 then
            invalid_arg(sprintf "%s: %s.(%i) < 0" name nname k);
          if inc.(k) > 0 then begin
            if ofs.(k) < FIRST_INDEX then
              invalid_arg(sprintf "%s: %s.(%i) < %i (%s)"
                            name ofsname k FIRST_INDEX LAYOUT);
            low.(k) <- ofs.(k);
            if n.(k) = 0 then (
              (* nk = max {n | ofs.(k) + (n-1) * inc.(k) < LAST_INDEX(dimk) }
                 with inc.(k) > 0 *)
              let nk = 1 + (LAST_INDEX(dimk) - ofs.(k)) / inc.(k) in
              if nk > 1 then set_n_sub nk
              else if nk < 1 then
                invalid_arg(sprintf "%s: dim %i empty; no n >= 1 s.t. %i\
	            %+i*(n-1) %s %i" name k ofs.(k) inc.(k) LT_DIM_SYM dimk);
              up.(k) <- ofs.(k) + (nk - 1) * inc.(k);
            )
            else if n.(k) > 1 then (
              (* n.(k) > 1 is given; bound check.   *)
              let last = ofs.(k) + (n.(k) - 1) * inc.(k) in
              if last > LAST_INDEX(dimk) then
                invalid_arg
                  (sprintf "%s: %s.(%i) + (%s.(%i) - 1) * %s.(%i) \
		     = %i %s %i (%s) where %s.(%i) = %i is the physical dim"
                     name ofsname k nname k incname k last GE_DIM_SYM dimk
                     LAYOUT nname k n.(k));
              set_n_sub n.(k);
              up.(k) <- last;
            )
            else (* n.(k) = 1 means: ignore this dimension *)
              up.(k) <- ofs.(k)
          end
          else if inc.(k) < 0 then begin
            if ofs.(k) > LAST_INDEX(dimk) then
              invalid_arg(sprintf "%s: %s.(%i) %s %i (%s)"
                            name ofsname k GE_DIM_SYM dimk LAYOUT);
            up.(k) <- ofs.(k);
            if n.(k) = 0 then (
              (* nk = max {n | FIRST_INDEX <= ofs.(k) + (n-1) * inc.(k) }
                 with inc.(k) < 0 *)
              let nk = 1 + (FIRST_INDEX - ofs.(k)) / inc.(k) in
              if nk > 1 then set_n_sub nk
              else if nk < 1 then
                invalid_arg(sprintf "%s: dim %i empty; no n >= 1 s.t. %i\
		  %+i*(n-1) >= %i" name k ofs.(k) inc.(k) FIRST_INDEX);
              low.(k) <- ofs.(k) + (nk - 1) * inc.(k);
            )
            else if n.(k) > 1 then (
              (* n.(k) > 1 is given; bound check *)
              let last = ofs.(k) + (n.(k) - 1) * inc.(k) in
              if last < FIRST_INDEX then
                invalid_arg
                  (sprintf "%s: %s.(%i) + (%s.(%i) - 1) * %s.(%i) \
		     = %i < %i (%s) where %s.(%i) = %i is the physical dim"
                     name ofsname k nname k incname k last FIRST_INDEX
                     LAYOUT nname k n.(k));
              set_n_sub n.(k);
              low.(k) <- last;
            )
            else (* n.(k) = 1 => ignore this dimension *)
              low.(k) <- ofs.(k)
          end
          else begin (* inc.(k) = 0 => dimension ignored for the transform. *)
            if ofs.(k) < FIRST_INDEX || ofs.(k) > LAST_INDEX(dimk) then
              invalid_arg(sprintf "%s: %s.(%i) = %i not in [%i, %i] (%s)"
                            name ofsname k ofs.(k) FIRST_INDEX
                            (LAST_INDEX(dimk)) LAYOUT);
            low.(k) <- ofs.(k);
            up.(k) <- ofs.(k);
          end;
          pdim := !pdim * dimk;
         );
  DEBUG(eprintf "DEBUG: %s: n_sub=%s (rank=%i); offset=%i stride=%s\n%!" name
          (string_of_array n_sub) !rank !offset (string_of_array stride));
  !offset, (Array.sub n_sub 0 !rank), stride, low, up
;;

(* Check whether the matrix given by [hm_n] (howmany matrix) and [hm]
   is a valid submatrix of the hermitian matrix [mat].  @return the
   [hm_stride], [hm_n] (howmany matrix), [stride] and [n] (logical
   dimensions). *)
let get_geom_hm name hm_nname hm_n hmname hm  nname n low up  mat =
  let num_dims = Genarray.num_dims mat in
  if hm = [] then
    if hm_n = [| |] then
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
        FOR_HM(k, num_dims,
               let dimk = Genarray.nth_dim mat k in
               if v.(k) > 0 then
                 (* max{j | up.(k) + v.(k)*(j-1) <= LAST_INDEX(dimk)} *)
                 ni := min !ni (1 + (LAST_INDEX(dimk) - up.(k)) / v.(k))
               else if v.(k) < 0 then
                 (* max{j | low.(k) + v.(k)*(j-1) >= FIRST_INDEX} *)
                 ni := min !ni (1 + (FIRST_INDEX - low.(k)) / v.(k));
               hm_s := !hm_s * dimk + v.(k); (* Horner *)
              );
        hm_n.(i) <- !ni
      )
      else (
        (* dimension [hm_n.(i)] provided; bound check *)
        FOR_HM(k, num_dims,
               let dimk = Genarray.nth_dim mat k in
               if (v.(k) > 0 && up.(k) + v.(k)*(hm_n.(i)-1) > LAST_INDEX(dimk))
                 || (v.(k) < 0 && low.(k) + v.(k)*(hm_n.(i)-1) < FIRST_INDEX)
               then
                 invalid_arg(sprintf "%s: translating %i times by the %ith \
			       vector %s of %s exceeds the %ith dim bounds"
                               name hm_n.(i) i (string_of_array v) hmname k);
               hm_s := !hm_s * dimk + v.(k); (* Horner *)
              );
      );
      if !hm_s = 0 then
        invalid_arg(sprintf "%s: %ith element of %s = [|0.;...;0.|]"
                      name i hmname);
      hm_stride.(i) <- !hm_s;
    end;
    DEBUG(eprintf "DEBUG: %s: hm_n=%s; hm_stride=%s\n%!" name
            (string_of_array hm_n) (string_of_array hm_stride));
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
  let hm_stridei, hm_ni =
    get_geom_hm name "howmany_n" hm_n "howmanyi" hmi  "n" ni lowi upi i
  and hm_strideo, hm_no =
    get_geom_hm name "howmany_n" hm_n "howmanyo" hmo  "n" no lowo upo o  in
  if hm_ni <> hm_no then
    invalid_arg(sprintf "%s: howmany dim input = %s <> howmany dim output = %s"
                  name (string_of_array hm_ni) (string_of_array hm_no));
  make_plan offseti offseto n stridei strideo  hm_ni hm_stridei hm_strideo
