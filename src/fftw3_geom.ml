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
  and n_log = Array.make num_dims 1 (* logical dimensions *)
  and offset = ref 0 (* external functions use the C layout *)
  and stride = Array.make num_dims 0
  and pdim = ref 1 (* product of physical dims > k; FORTRAN: < k *)
  and low = Array.make num_dims 0 (* lower submatrix corner *)
  and up = Array.make num_dims 0 (* greater submatrix corner *) in
  let set_n_log v = n_log.(!rank) <- v; incr rank in
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
              if nk > 1 then set_n_log nk
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
              set_n_log n.(k);
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
              if nk > 1 then set_n_log nk
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
              set_n_log n.(k);
              low.(k) <- last;
            )
            else (* n.(k) = 1 means: ignore this dimension *)
              low.(k) <- ofs.(k)
          end
          else begin (* inc.(k) = 0 =>  dimension ignored for the transform. *)
            if ofs.(k) < FIRST_INDEX || ofs.(k) > LAST_INDEX(dimk) then
              invalid_arg(sprintf "%s: %s.(%i) = %i not in [%i, %i] (%s)"
                            name ofsname k ofs.(k) FIRST_INDEX
                            (LAST_INDEX(dimk)) LAYOUT);
            low.(k) <- ofs.(k);
            up.(k) <- ofs.(k);
          end;
          pdim := !pdim * dimk;
         );
  DEBUG(eprintf "%s: n_log=%s; stride=%s\n%!" name (string_of_array n_log)
          (string_of_array stride));
  !offset, (Array.sub n_log 0 !rank), stride, low, up
;;

(* Check whether the matrix given by [hm_n] (howmany matrix) and [hm]
   is a valid submatrix of the hermitian matrix [mat].  @return the
   [hm_stride], [hm_n] (howmany matrix), [stride] and [n] (logical
   dimensions). *)
let get_geom_hm name hm_nname hm_n hmname hm  nname n low up  mat =
  let rank = Genarray.num_dims mat in
  if hm = [] then
    if hm_n = [| |] then
      [| |], [| |], stride, n (* only one transform *)
    else begin
      (* Transforms indexed by the first (FORTRAN: last) dimensions of
         the submatrix of [mat] defined by [ofs] and [inc]. *)
      let hm_rank = Array.length hm_n in
      if hm_rank >= rank then
        invalid_arg(sprintf "%s: length %s >= %i = number of dimensions"
                      name hm_nname rank);
      let mat_rank = rank - hm_rank in
      let hm_index = IFDEF FORTRAN THEN mat_rank ELSE 0 ENDIF in
      let hm_stride = Array.sub stride hm_index hm_rank in
      let stride_index = IFDEF FORTRAN THEN 0 ELSE hm_rank ENDIF in
      let stride = Array.sub stride stride_index mat_rank in
      let normalize_hm i ni =
        let j = IFDEF FORTRAN THEN mat_rank + i ELSE i ENDIF in
        if ni < 0 then invalid_arg(sprintf "%s: %s.(%i) < 0" name hm_nname i)
        else if ni = 0 then n.(j)
        else if ni <= n.(j) then ni
        else invalid_arg(sprintf "%s: %s.(%i) = %i > %s.(%i) = %i"
                           name hm_nname i ni nname j n.(j)) in
      let hm_n = Array.mapi normalize_hm hm_n in
      let n = Array.sub n stride_index mat_rank in
      hm_stride, hm_n, stride, n
    end
  else (
    (* Transforms indices = vectors of [hm] with desired dims [hm_n]. *)
    let hm_rank = List.length hm in
    let hm_n =
      if hm_n = [| |] then Array.make hm_rank 0
      else if Array.length hm_n = hm_rank then (
        let normalize_hm i ni =
          if ni < 0 then invalid_arg(sprintf "%s: %s.(%i) < 0" name hm_nname i);
          ni in
        (* copy [hm_n] because 0 entries will be modified: *)
        Array.mapi normalize_hm hm_n;
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
        FOR_HM(k, rank,
               let dimk = Genarray.nth_dim mat k in
               if v.(k) > 0 then
                 (* max{j : ofs.(k) + inc.(k)^+ * (n.(k)-1) + v.(k)*(j-1)
                    < dimk}.   FORTRAN: <= *)
                 let d = LAST_INDEX(dimk) - ofs.(k) - pos inc.(k) * (n.(k)-1) in
                 ni := min !ni (1 + d / v.(k))
               else if v.(k) < 0 then
                 (* max{j : ofs.(k) + inc.(k)^- * (n.(k)-1) + v.(k)*(j-1) >= 0}
                    FORTRAN: >= 1 *)
                 let d = FIRST_INDEX - ofs.(k) - neg inc.(k) * (n.(k)-1) in
                 ni := min !ni (1 + d / v.(k));
                 hm_s := !hm_s * dimk + v.(k); (* Horner *)
              );
        hm_n.(i) <- !ni
      )
      else (
        (* dimension [hm_n.(k)] provided; bound check *)
        FOR_HM(k, rank,
               let dimk = Genarray.nth_dim mat k in
               if (v.(k) > 0 &&
                     ofs.(k) + pos inc.(k)*(n.(k)-1) + v.(k)*(hm_n.(i)-1)
                   > LAST_INDEX(dimk))
                 || (v.(k) < 0 &&
                       ofs.(k) + neg inc.(k) * (n.(k)-1) + v.(k)*(hm_n.(i)-1)
                     < FIRST_INDEX)
               then
                 invalid_arg(sprintf "%s: translating %i times the %ith \
			       vector of %s exceeds the %ith dim array bounds"
                               name hm_n.(i) i hmname k);
               hm_s := !hm_s * dimk + v.(k); (* Horner *)
              );
      );
      if !hm_s = 0 then
        invalid_arg(sprintf "%s: %ith element of %s = [|0.;...;0.|]"
                      name i hmname);
      hm_stride.(i) <- !hm_s;
    end hm;
    hm_stride, hm_n, stride, n
  )


(* Take the [make_plan] function creating plans, the dimensions, offsets
   and increments of input/output arrays and compute the informations
   needed by [make_plan].  Check the coherence of the data at the same
   time.  There may be more input (resp. output) arrays than [i]
   (resp. [o]) but these must have the same dimensions. *)
let apply name make_plan n hm_n  hmi ofsi inci i  hmo ofso inco o =
  let rank = Genarray.num_dims i in
  if rank <> Genarray.num_dims o then
    invalid_arg(name ^ ": input and output arrays do not have the same \
	NUMBER of dimensions");
  let offseti, ni, stridei, lowi, upi =
    get_geom name "ofsi" ofsi "inci" inci "n" n i
  and offseto, no, strideo, lowo, upo =
    get_geom name "ofso" ofso "inco" inco "n" n o in
  if ni <> no then invalid_arg("dim input = " ^ string_of_array ni
                               ^ " <> dim ouput = " ^ string_of_array no);
  let hm_stridei, hm_ni, stridei, ni =
    get_geom_hm name "howmany_n" hm_n "howmanyi" hmi stridei ofsi inci "n" ni i
  and hm_strideo, hm_no, strideo, no =
    get_geom_hm name "howmany_n" hm_n "howmanyo" hmo strideo ofso inco "n" no o
  in
  if hm_ni <> hm_no then
    invalid_arg("howmany dim input = " ^ string_of_array hm_ni
                ^ " <> howmany dim output = " ^ string_of_array hm_no);
  make_plan offseti offseto ni stridei strideo  hm_ni hm_stridei hm_strideo
