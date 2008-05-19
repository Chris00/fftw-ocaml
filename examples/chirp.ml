(* FFT analysis of a "chirp" signal.

   Inspired by
   http://ccrma.stanford.edu/~jos/sasp/Computational_Examples_Matlab.html

   Remark: The could could have been shortened by using Lacaml.  We
   did not however to avoid this additional dependency.
*)

open Bigarray
module FFT = Fftw3.D
module G = Gnuplot.Bigarray

type vec = fortran_layout FFT.Array1.complex_array

let pi = 4. *. atan 1.

(** Chirp signal
 ***********************************************************************)

type shape = Linear | Quadratic | Logarithmic

let option_default f0 = function Some f -> f | None -> f0
let pi = 4. *. atan 1.

(** [chirp ?f0 ?t1 ?f1 shape] generates a linear swept-frequency
    cosine signal function, where [f0] is the instantaneous frequency
    at time 0, and [f1] is the instantaneous frequency at time [t1].
    [f0] and [f1] are both in hertz.  If unspecified, [f0] is e-6 for
    logarithmic chirp and 0 for all other methods, [t1] is 1, and [f1]
    is 100. *)
let chirp ?f0 ?(t1=1.) ?(f1=100.) ?(phase=0.) = function
  | Linear ->
      let f0 = option_default 0. f0 in
      let a = pi *. (f1 -. f0) /. t1
      and b = 2. *. pi *. f0  in
      (fun t -> cos(a *. t *. t +. b *. t +. phase))
  | Quadratic ->
      let f0 = option_default 0. f0 in
      let a = 2. /. 3. *. pi *. (f1 -. f0) /. (t1 *. t1)
      and b = 2. *. pi *. f0 in
      (fun t -> cos(a *. t *. t *. t +. b *. t +. phase))
  | Logarithmic ->
      let f0 = option_default (exp(-. 6.)) f0 in
      let df = f1 -. f0 in
      let a = 2. *. pi *. t1 /. log df
      and b = 2. *. pi *. f0
      and x = df**(1. /. t1) in
      (fun t -> cos(a *. x**t +. b *. t +. phase))

(** Filter function (akin matlab one)
 ***********************************************************************)

let create kind n = FFT.Array1.create kind fortran_layout n

let copy0 (x:vec) (y:vec) =
  for i = 1 to Array1.dim x do y.{i} <- x.{i} done;
  for i = Array1.dim x + 1 to Array1.dim y do y.{i} <- Complex.zero done

let scal s (x:vec) =
  for i = 1 to Array1.dim x do
    let xi = x.{i} in
    x.{i} <- { Complex.re = s *. xi.Complex.re; im = s *. xi.Complex.im }
  done

(** [filter b ?a x] returns [y] the data in vector [x] filtered with
    the filter described by vectors [a] and [b].  The filter
    implements of the standard difference equation:
    {v
    a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
                          - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
    v}
    for n = 1,..., [Array1.dim x].  *)
let filter (b:vec) ?a (x:vec) =
  let a, a_provided = match a with
    | None -> create FFT.complex 0, false
    | Some a -> a, true in
  let n = max (Array1.dim a) (Array1.dim b) + Array1.dim x - 1 in
  let y = create FFT.complex n
  and b' = create FFT.complex n in
  let fftb = FFT.Array1.dft FFT.Forward y b' in
  let x' = create FFT.complex n in
  let fftx = FFT.Array1.dft FFT.Forward y x' in
  let y' = create FFT.complex n in
  let iffty = FFT.Array1.dft FFT.Backward y' y in
  copy0 b y;  FFT.exec fftb;
  copy0 x y;  FFT.exec fftx;
  if a_provided then (
    let a' = create FFT.complex n in
    let ffta = FFT.Array1.dft FFT.Forward y a' in
    copy0 a y;  FFT.exec ffta;
    (* FIXME: not correct: *)
    for i = 1 to n do
      y'.{i} <- Complex.div (Complex.mul b'.{i} x'.{i}) a'.{i}
    done;
  )
  else
    (* Consider a.{1} = 1 and a.{i} = 0 for i > 1 *)
    for i = 1 to n do y'.{i} <- Complex.mul b'.{i} x'.{i} done;
  FFT.exec iffty;
  scal (1. /. float n) y;              (* normalize ifft *)
  Array1.sub y 1 (Array1.dim x)
;;

let () =
  let cpl x = {Complex.re = x; Complex.im = 0.} in
  let of_array a =
    Array1.of_array FFT.complex fortran_layout (Array.map cpl a) in
  let a = of_array [| 1.;2. |] in
  let b = of_array [| 1.; 1.; 1. |] in
  let x = create FFT.complex 20 in
  for i = 1 to Array1.dim x do x.{i} <- cpl(float i) done;
  let y  = filter ~a b x in
  for i = 1 to Array1.dim y do
    Printf.printf "%g " y.{i}.Complex.re
  done;
  Printf.printf "\n%!";
;;

let () =
  let n = 10          (* number of filters = DFT length *)
  and fs = 1000.      (* sampling frequency (arbitrary) *)
  and d = 1. in       (* duration in seconds *)
  let n_graphs = 5 in (* technically =n but not space for all graphs *)

  let len = truncate(ceil(fs *. d)) + 1 in (* signal duration (samples) *)
  (* sine sweep from 0 Hz to fs/2 Hz: *)
  let ch = chirp ~t1:d ~f1:(0.5 *. fs) Linear in
  let t = create FFT.float len
  and x = create FFT.float len in
  for i = 1 to len do
    t.{i} <- float(i-1) /. fs;
    x.{i} <- ch t.{i}
  done;
  let h = create FFT.complex n in
  Array1.fill h Complex.one;            (* h = [| 1.; ...; 1. |] *)

  let g = G.init G.X ~xsize:1000. ~ysize:800. ~nysub:(n_graphs + 1) in
  G.box g;
  G.pen g 1;
  G.xy g t x;
  let xk = create FFT.complex (Array1.dim x) in
  let y_re = create FFT.float (Array1.dim x) in
  for k = 1 to n_graphs do
    G.adv g;
    G.pen g 0;
    G.box g;
    (* Modulation by the complex exponential i -> exp(-I*wk*i) *)
    let wk = 2. *. pi *. float(k-1) /. float n in
    for i = 1 to Array1.dim x do
      let theta = -. wk *. float i in
      xk.{i} <- { Complex.re = x.{i} *. cos theta; im = x.{i} *. sin theta }
    done;
    (* Filter and display the real part *)
    let y = filter h xk in
    for i = 1 to Array1.dim y do y_re.{i} <- y.{i}.Complex.re done;
    G.pen g 3;
    G.xy g t y_re
  done;
  G.close g
