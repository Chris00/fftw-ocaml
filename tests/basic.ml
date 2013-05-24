open Printf
open Bigarray
module FFT = Fftw3.D


let () =
  let d = 16 in
  let u = FFT.Array2.create float64 fortran_layout d d
  and fu = FFT.Array2.create complex64 fortran_layout (d/2 + 1) d in
  let forward = FFT.Array2.r2c u fu ~meas:FFT.Measure in
  let backward = FFT.Array2.c2r fu u ~meas:FFT.Measure in
  for i = 1 to d do
    for j = 1 to d do
      u.{i, j} <- float(i+j)
    done
  done;
  FFT.exec forward;
  FFT.exec backward;
  (* [u] must be multiplied by d² *)
  for i = 1 to d do
    for j = 1 to d do
      let r = float((i + j) * d * d) in
      if abs_float(u.{i,j} -. r) > 1e-10 then
        printf "u.{%i,%i} = %g (expected %g)\n" i j u.{i,j} r
    done
  done


(* Local Variables: *)
(* compile-command: "make -k -C .. tests" *)
(* End: *)
