(*pp camlp4o pa_macro.cmo $GNUPLOT_EXISTS *)
(** Maximum entropy method.

    Example demonstrating the "Maximum Entropy Method" as explained in
    the book "OCaml for Scientists" by Dr Jon D Harrop.  Although the
    ideas are those of the book, the code has been completely rewritten.

    Required modules: Lacaml
*)

open Printf
open Scanf
open Bigarray
open Lacaml.D
module FFT = Fftw3.D

IFDEF GNUPLOT_EXISTS THEN
module G = Gnuplot.Bigarray
ENDIF;;

let delta = sqrt epsilon_float

(* creates a vector (1D bigarray) guaranteed to be aligned in memory. *)
let create n : vec = FFT.Array1.create FFT.float fortran_layout n


(** [grad f x] approximate the gradient of [f] at [x] by means of
    divided differences (of order 1). *)
let grad f (x:vec) =
  let g = create (Vec.dim x) in
  let fx = f x in
  for i = 1 to Vec.dim x do
    let xi = x.{i} in
    x.{i} <- x.{i} +. delta;
    g.{i} <- (f x -. fx) /. delta;
    x.{i} <- xi
  done;
  g

(** [gradient_ascent f f' x] modifies the vector [x] until it reaches
    a local maximum of [f] using a simple gradient ascent algorithm. *)
let gradient_ascent ?(print=fun _ _ _ _ _ -> ()) (f: vec -> float) f' (x:vec) =
  let new_x = create (Vec.dim x) in
  (* The code is imperative so floats are unboxed by the compiler *)
  let lambda = ref delta
  and fx = ref (f x)
  and step = ref 1 in
  while !lambda >= delta do
    let f'x = f' x in
    print !step !lambda x !fx (nrm2 f'x);
    incr step;
    Array1.blit x new_x;
    axpy ~alpha:!lambda ~x:f'x new_x;
    let f_new_x = f new_x in
    if f_new_x <= !fx then lambda := 0.5 *. !lambda
    else (
      lambda := 1.1 *. !lambda;
      Array1.blit new_x x;
      fx := f_new_x;
    )
  done


(** [mem const n] returns a vector [x] of dimension [n] whose first
    components are [const] (the measured data) and the last
    components are obtained by maximising the Shannon entropy
    {% $$
    H(x) = \ln\Bigl(\sum p_k\Bigr) - \frac{\sum p_k \ln p_k}{\sum p_k}
    $$ %}
    where {% $p_k$ %} is the square root of [y.{k}] and [y] is the
    discrete sine transform of [x]. *)
let mem (consts:vec) n =
  let d = Array1.dim consts in
  if n <= d then invalid_arg(sprintf "mem: 'n' is too small, pick n > %i" d);
  let x = create n
  and y = create n (* Fourier variables *) in
  (* We will spend a lot of time doing FFT so it is worth making sure
     a good plan is generated: *)
  let ffst = FFT.Array1.r2r FFT.RODFT10 x y ~meas:FFT.Patient in
  ignore(copy consts ~y:x); (* first components of [x] are [const] *)
  let vars = Array1.sub x (d + 1) (n - d) (* last components of [x] to set *) in
  let entropy (v:vec) =
    Array1.blit v vars;
    FFT.exec ffst;
    let s = ref 0.
    and h = ref 0. in
    for k = 1 to n do
      let yk = abs_float y.{k} in
      if yk > 0. then begin
	s := !s +. yk;
	h := !h +. yk *. log yk; (* y=0 => y *. log y = 0 *)
      end
    done;
    (* Given the discussion in Jon D Harrop book, the term [!h] of the
       entropy should be divided by [!s] (this is a mistake in the
       book).  However, this gives non-physical results which is the
       reason why I left it commented here.  *)
    log !s -. !h (* /. !s *)
  in
  let v = Vec.make0 (n-d) (* starting point *) in
  gradient_ascent entropy (grad entropy) v
    ~print:(fun i l x h h' ->
              Format.eprintf
                "%3i: l = %#-6.g  x = %#-2g  H(x) = %1.14f  H'(x) = %.5g@\n%!"
                i l (nrm2 x) h h');
  Array1.blit v vars;
  x (* return [const] followed by the optimal [vars]. *)

(* Rename this alternative to [mem] function if you want to recover
   the very same results as given by the downloadable code from the
   "OCaml for Scientists" book. *)
let mem_harrop (consts: vec) n =
  let d = Array1.dim consts in
  if n <= d then invalid_arg(sprintf "mem: 'n' is too small, pick n > %i" d);
  let x = FFT.Array1.create FFT.complex fortran_layout (2*n)
  and y = FFT.Array1.create FFT.complex fortran_layout (2*n) in
  let ffst = FFT.Array1.dft FFT.Forward x y in
  x.{1} <- Complex.zero; (* ouch: [consts.{1}] treated as 0. *)
  for k = 2 to d do
    let z = {Complex.re = consts.{k}; im = 0.} in
    x.{k} <- z;
    x.{2*n + 2 - k} <- Complex.neg z
  done;
  x.{n+1} <- Complex.zero;
  let entropy (v:vec) =
    for k = d + 1 to n do
      let z = {Complex.re = v.{k - d}; im = 0.} in
      x.{k} <- z;
      x.{2*n + 2 - k} <- Complex.neg z
    done;
    FFT.exec ffst;
    let s = ref 0.
    and h = ref 0. in
    for k = 1 to n do
      let yk = abs_float y.{k}.Complex.im in
      if yk > 0. then begin
	s := !s +. yk;
	h := !h +. yk *. log yk;
      end
    done;
    log !s -. !h (* /. !s *)
  in
  let v = Vec.make0 (n-d) (* starting point *) in
  gradient_ascent entropy (grad entropy) v
    ~print:(fun i l x h h' ->
              Format.eprintf
                "%3i: l = %#-6.g  x = %#-2g  H(x) = %1.14f  H'(x) = %.5g@\n%!"
                i l (nrm2 x) h h');
  let x = copy consts ~y:(create n) (* first components of [x] are [const] *) in
  ignore(copy v ~y:x ~ofsy:(d + 1));
  x


(** Read data of the form {% $\\{ f_1, \ldots, f_p \\}$ %} from [fh]
    (very tolerant to errors). *)
let read_consts bh =
  let d = ref 0 and l = ref [] in
  try
    bscanf bh " { " ();
    while true do
      bscanf bh "%f%[, \n]" (fun f _ -> l := f :: !l; incr d)
    done;
    assert false
  with Scanf.Scan_failure _ ->
    let c = create !d in
    ignore(List.fold_left (fun k f -> c.{k} <- f; k-1) !d !l);
    c

(* Print a vector in Mathematica format (automatically splitted on
   several lines if needed). *)
let print_vec fh (x:vec) =
  Format.fprintf fh "{ @[";
  for i = 1 to Vec.dim x - 1 do Format.fprintf fh "%g,@ " x.{i} done;
  Format.fprintf fh "%f @]}@\n" x.{Vec.dim x}


let () =
  let input = ref "" in
  let n = ref 1200 in
  let output = ref "" in
  let usage = sprintf "Usage: %s" Sys.argv.(0) in
  let args = Arg.align [
    ("-n", Arg.Set_int n,
     Printf.sprintf "i extend the data to i elements (default: %i)" !n);
    ("--input", Arg.Set_string input, "fname input file.  Expected to contain \
	a string of the form {f1,...,fN} (N < n).");
    ("--output", Arg.Set_string output,
     "fname output file (default: standard output or graphical display \
      if gnuplot is found).");
  ] in
  Arg.parse args (fun s -> raise(Arg.Bad "no anonymous argument")) usage;
  if !input = "" then (Arg.usage args usage; exit 1);

  let consts = read_consts (Scanning.from_file !input) in
  let x = mem consts !n in
  if !output = "" then begin
    IFDEF GNUPLOT_EXISTS THEN
      (* Plot the graph of the result, indicating what is extended *)
      let g = G.init G.X in
      G.box g;
      G.pen g 1;
      G.x g ~n0:(Vec.dim consts) x ~ofsx:(Vec.dim consts);
      G.pen g 2;
      G.x g consts;
      G.close g;
    ELSE
      (* Spit out the result on the standard output *)
      print_vec Format.std_formatter x
    ENDIF
  end
  else begin
    let fh = open_out !output in
    fprintf fh "# Maximum Entropy Method (#const=%i, n=%i)\n"
      (Array1.dim consts) !n;
    for i = 1 to Vec.dim x - 1 do fprintf fh "%g\n" x.{i} done;
    close_out fh
  end
