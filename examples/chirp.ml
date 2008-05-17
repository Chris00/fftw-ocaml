(* FFT analysis of a "chirp" signal.  *)

open Bigarray
module FFT = Fftw3.D
type vec = c_layout FFT.Array1.float_array

module G = Gnuplot.Bigarray

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



let () =
  let n = 10                            (* number of filters = DFT length *)
  and fs = 1000.                        (* sampling frequency (arbitrary) *)
  and d = 1. in                         (* duration in seconds *)

  let len = truncate(ceil(fs *. d)) + 1 in (* signal duration (samples) *)
  (* sine sweep from 0 Hz to fs/2 Hz: *)
  let ch = chirp ~t1:d ~f1:(0.5 *. fs) Linear in

  let t = FFT.Array1.create FFT.float fortran_layout len
  and x = FFT.Array1.create FFT.float fortran_layout len in
  let x' = FFT.Array1.create FFT.complex fortran_layout (len/2+1) in
  let plan = FFT.Array1.r2c x x' in
  for i = 1 to len do
    t.{i} <- float(i-1) /. fs;
    x.{i} <- ch t.{i}
  done;
  FFT.exec plan;

  let g = G.init G.X ~xsize:1000. ~ysize:300. ~nysub:2 in
  G.box g;
  G.pen g 1;
  G.xy g t x;
  G.adv g;
  G.pen g 0;
  G.box g;
  let y = FFT.Array1.create FFT.float fortran_layout (len/2+1) in
  for i = 1 to Array1.dim y do y.{i} <- Complex.norm x'.{i} done;
  G.pen g 3;
  G.xy g t y;
  G.close g
