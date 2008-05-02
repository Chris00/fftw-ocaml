(*pp camlp4o pa_macro.cmo $FFTW3F_EXIST *)
(* File: fftw3.mli

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



(** Interface for FFTW version 3.

    @author Christophe Troestler <chris_77\@users.sourceforge.net>
    @version 0.5.1
*)


(** Precision independent signature for FFTW3 submodules.

    We advise against opening this module as it contains submodules with
    the same names as the [Bigarray] ones.  Instead, declare
    {[
      module FFT = Fftw3.D 						]}
    or
    {[
      module FFT = Fftw3.S 						]}
    depending to the precision you need (this way of proceeding makes it
    easy to change the precision of the FFT sould it be necessary) and
    then use it as
    {[
      let x = FFT.Array1.create FFT.complex c_layout dim
      let y = FFT.Array1.create FFT.complex c_layout dim
      let plan = FFT.Array1.dft FFT.Forward x y
      (* fill x and y *)
      FFT.exec plan (* perform the DFT *) 				]}
    The last line can be repeated as many times as needed to compute the
    FFT of [x] into [y].  {b Beware} that creating the plan usually
    destroys the content of [x] and [y], so only fill them afterwards.

    HINT: Plan creation functions like {!Fftw3.Sig.Array1.dft} have
    many optional arguments for maximum flexibility.  The two
    important ones are [~meas] and [~normalize].  The other ones can
    be ignored at first.
*)
module type Sig = sig
  open Bigarray

  (** {2 Precision} *)

  type float_elt
    (** Precision of float numbers. *)
  type complex_elt
    (** Precision of complex numbers. *)

  val float : (float, float_elt) Bigarray.kind
    (** Float of the precision of this module.  Use this to create
        precision independent code. *)
  val complex : (Complex.t, complex_elt) Bigarray.kind
    (** Complex of the precision of this module.  Use this to create
        precision independent code. *)


  (** {2 Specifying plans} *)

  type 'a plan (** Immutable FFTW plan. *)
  type c2c     (** [c2c plan] usual discrete Fourier transform,
                   from complex to complex *)
  type r2c     (** [r2c plan] real to complex transform *)
  type c2r     (** [c2r plan] complex to real transform *)
  type r2r     (** [r2r plan] real to real transform *)

  (** Direction of the transform -- see the FFTW manual. *)
  type dir = Forward | Backward

  (** Planning-rigor flags. *)
  type measure =
    | Estimate (** No measurements are made, use a simple heuristic to
                   pick a (probably sub-optimal) plan quickly. *)
    | Measure (** Find an optimized plan by actually computing several
                  FFTs and measuring their execution time. *)
    | Patient (** Like [Measure], but considers a wider range of
                  algorithms and often produces a "more optimal" plan
                  at the expense of several times longer planning
                  time. *)
    | Exhaustive (** Like [Patient], but considers an even wider range
                     of algorithms, including many that we think are
                     unlikely to be fast, to produce the most optimal
                     plan but with a substantially increased planning
                     time. *)

  (** Real-to-Real transform kinds.  The real-even (resp. real-odd) DFT
      are somtimes called Discrete Cosine Transform (DCT)
      (resp. Discrete Sine Transform (DST)).  Note that the explanations
      of the various transforms are for an {i input} array of dimension
      [n] and C layout (i.e. the input array is [input[0..n-1]]).  The
      logical size [N] is [N=2(n-1)] for [REDFT00], [N=2(n+1)] for
      [RODFT00], and [N=2n] otherwise.  See the FFTW manual for more
      details. *)
  type r2r_kind =
    | R2HC (** real to halfcomplex *)
    | HC2R (** halfcomplex to real *)
    | DHT  (** discrete Hartley Transform *)
    | REDFT00 (** real-even DFT: even around j=0 and even around j=n-1 *)
    | REDFT10 (** real-even DFT: even around j=-0.5 and even around j=n-0.5 *)
    | REDFT01 (** real-even DFT: even around j=0 and odd around j=n *)
    | REDFT11 (** real-even DFT: even around j=-0.5 and odd around j=n-0.5 *)
    | RODFT00 (** real-odd DFT; odd around j=-1 and odd around j=n *)
    | RODFT10 (** real-odd DFT; odd around j=-0.5 and odd around j=n-0.5 *)
    | RODFT01 (** real-odd DFT; odd around j=-1 and even around j=n-1 *)
    | RODFT11 (** real-odd DFT; odd around j=-0.5 and even around j=n-0.5 *)


  (** {2 Executing plans} *)

  val exec : 'a plan -> unit
    (** [exec plan] executes the [plan] on the arrays given at the
        creation of this plan.  This is the normal way to execute any
        kind of plan.

        This function is thread safe.  *)

  (*
  (** {3 Guru execution of plans.}

    If you want to transform other arrays than those specified in the
    plan, you are advised to create a new plan -- it won't be too
    expensive if the wisdom can be reused.  To transform a known bunch
    of arrays of the same size, you should {b not} use the following
    functions but instead create a plan with [?howmany] set
    appropriately.

    These functions are thread safe.  You can even execute the {i same
    plan} in parallel by multiple threads by providing different
    arrays than the ones with which the plan was created.
  *)

    val exec_dft : c2c plan -> 'l complex_array -> 'l complex_array -> unit

    val exec_split_dft : c2c plan ->
    'l float_array -> 'l float_array ->
    'l float_array -> 'l float_array -> unit

    val exec_r2c : r2c plan -> 'l float_array -> 'l complex_array -> unit

    val exec_split_r2c : r2c plan ->
    'l float_array -> 'l float_array -> 'l float_array -> unit

    val exec_c2r : c2r plan -> 'l complex_array -> 'l float_array -> unit

    val exec_split_c2r : c2r plan ->
    'l float_array -> 'l float_array -> 'l float_array -> unit

    val exec_r2r : r2r plan -> 'l float_array -> 'l float_array -> unit
  *)


  (** {2 Creating plans} *)

  module Genarray :
  sig
    external create: ('a, 'b) kind -> 'c layout -> int array
      -> ('a, 'b, 'c) Bigarray.Genarray.t
      = "fftw3_ocaml_ba_create"
      (** Creates a new array, just as [Bigarray.Genarray.create] does,
	  but guarantees that it is aligned so one gets the better
	  performance from FFTW.

          Remark: In order to deserialize such a bigarray, this module
          must be linked to the program as the deserialization
          function also aligns the data. *)

    type 'l complex_array = (Complex.t, complex_elt, 'l) Bigarray.Genarray.t
        (** Double precision complex array. *)
    type 'l float_array   = (float, float_elt, 'l) Bigarray.Genarray.t
        (** Double precision float array. *)

    val dft : dir ->
      ?meas:measure -> ?normalize:bool ->
      ?preserve_input:bool -> ?unaligned:bool ->
      ?n:int array -> ?howmany_n:int array ->
      ?howmanyi:int array list ->
      ?ofsi:int array -> ?inci:int array -> 'l complex_array ->
      ?howmanyo:int array list ->
      ?ofso:int array -> ?inco:int array -> 'l complex_array
      -> c2c plan
      (** [dft dir in out] returns a plan for computing the FFT in the
	  direction [dir] from [in] to [out].  [in] and [out] must have
	  the same number of dimensions and may be equal (in which case
	  the transform is done in-place).  If they are not equal, they
	  should not overlap.
          @raise Failure if the plan cannot be created.

          {b Beware} that, unless [~meas] is [Estimate], creating a plan
          requires some trials which will destroy the content of the
          arrays.

	  - [meas] controls how much time is dedicated to the
	  creation of the plan.  Default: [Measure]

	  - [normalize] if [true], divide the result by sqrt(n), where n
	  is the size of the FFT vector -- the product of the logical
	  dimensions.  Forward and Backward normalized transforms are
	  inverse of each other.  Default: [false].

	  - [preserve_input] specifies that an out-of-place transform
	  must {i not change its input} array.  Default: [true] except for
	  c2r and HC2R. FIXME: [false] means that... (??)

          - [unaligned] specifies that the algorithm may not impose any
          unusual alignment requirements.  You normally do not need this
          flag unless you want to use the plan with {i other unaligned
          arrays} (using the guru interface).  Default: [false] meaning
          that alignment may be used to speed up the computations (when
          [in] and [out] are aligned of course).

          {9 Subarrays}

          Fftw3 allows you to perform the FFT transform on subarrays
          defined by offset, strides and dimensions.

          - [ofsi] the initial element in the input array.  Default:
          [[|0;...;0|]] for c_layout and [[|1;...;1|]] for fortran_layout.

          - [inci] an array of increments for each dimension.
          Default: [[|1;...;1|]].

          - [ofso] same as [ofsi] but for output.

          - [inco] same as [inci] but for output.

          {9 Multiple transforms}

	  - [howmany_n] allow to compute several transforms at once
	  by specifying the number of dimensions ([< Genarray.num_dims
	  i]) used to index the arrays.  Default: 0.  *)

    val r2c : ?meas:measure -> ?normalize:bool ->
      ?preserve_input:bool -> ?unaligned:bool ->
      ?n:int array -> ?howmany_n:int array ->
      ?howmanyi:int array list ->
      ?ofsi:int array -> ?inci:int array -> 'l float_array ->
      ?howmanyo:int array list ->
      ?ofso:int array -> ?inco:int array -> 'l complex_array
      -> r2c plan
      (** [r2c in out] returns a plan for computing the {i forward}
          transform

          - [n] is the array of {i logical} sizes of the transform.

	  See {!Fftw3.D.Genarray.dft} for the meaning of the other
	  optional parameters. *)

    val c2r : ?meas:measure -> ?normalize:bool ->
      ?preserve_input:bool -> ?unaligned:bool ->
      ?n:int array -> ?howmany_n:int array ->
      ?howmanyi:int array list ->
      ?ofsi:int array -> ?inci:int array -> 'l complex_array ->
      ?howmanyo:int array list ->
      ?ofso:int array -> ?inco:int array -> 'l float_array
      -> c2r plan
      (** [c2r in out] returns a plan for computing the {i backward}
          transform

          destroy its input

          - [n] is the array of {i logical} sizes of the transform.

	  See {!Fftw3.D.Genarray.dft} for the meaning of the other
	  optional parameters. *)

    val r2r : r2r_kind array ->
      ?meas:measure -> ?normalize:bool ->
      ?preserve_input:bool -> ?unaligned:bool ->
      ?n:int array -> ?howmany_n:int array ->
      ?howmanyi:int array list ->
      ?ofsi:int array -> ?inci:int array -> 'l float_array ->
      ?howmanyo:int array list ->
      ?ofso:int array -> ?inco:int array -> 'l float_array
      -> r2r plan
      (** [r2r kind in out]

          See {!Fftw3.D.Genarray.dft} for the meaning of optional parameters. *)
  end


  module Array1 :
  sig
    val create: ('a, 'b) kind -> 'c layout -> int -> ('a, 'b, 'c) Array1.t
      (** See {!Fftw3.D.Genarray.create}. *)

    val of_array : ('a, 'b) kind -> 'c layout -> 'a array -> ('a, 'b, 'c) Array1.t
      (** [of_array kind layout a] build a one-dimensional aligned big
          array initialized from the given array. *)

    type 'l complex_array = (Complex.t, complex_elt, 'l) Array1.t
        (** Double precision complex 1D array. *)
    type 'l float_array   = (float, float_elt, 'l) Array1.t
        (** Double precision float 1D array. *)


    val dft : dir -> ?meas:measure -> ?normalize:bool ->
      ?preserve_input:bool -> ?unaligned:bool ->
      ?n:int -> ?howmany_n:int array ->
      ?howmanyi:int list -> ?ofsi:int -> ?inci:int -> 'l complex_array ->
      ?howmanyo:int list -> ?ofso:int -> ?inco:int -> 'l complex_array
      -> c2c plan
      (** [dft dir x y] returns a plan to compute the DFT of [x] and
          store it in [y].

          The parameters [meas], [preserve_input], [unaligned],
          [normalize] are as for {!Fftw3.D.Genarray.dft}.

          @param n the logical length of the array.  If not provided, it
          is automatically computed from [ofsi], [inci] and [Array1.dim
          x].

          Remark: If you want to transform several 1D arrays at once,
          use {!Fftw3.D.Array2.dft} with [~many:true]. *)

    val r2r : r2r_kind -> ?meas:measure -> ?normalize:bool ->
      ?preserve_input:bool -> ?unaligned:bool ->
      ?n:int -> ?howmany_n:int array ->
      ?howmanyi:int list -> ?ofsi:int -> ?inci:int -> 'l float_array ->
      ?howmanyo:int list -> ?ofso:int -> ?inco:int -> 'l float_array
      -> r2r plan
  end


  module Array2 :
  sig
    val create: ('a, 'b) kind -> 'c layout -> int -> int -> ('a, 'b, 'c) Array2.t
      (** See {!Fftw3.D.Genarray.create}. *)

    type 'l complex_array = (Complex.t, complex_elt, 'l) Array2.t
        (** Double precision complex 2D array. *)
    type 'l float_array   = (float, float_elt, 'l) Array2.t
        (** Double precision float 2D array. *)

    val dft : dir ->
      ?meas:measure -> ?normalize:bool ->
      ?preserve_input:bool -> ?unaligned:bool ->
      ?n: int * int -> ?howmany_n:int array ->
      ?howmanyi:(int * int) list ->
      ?ofsi:int * int -> ?inci:int * int -> 'l complex_array ->
      ?howmanyo:(int * int) list ->
      ?ofso:int * int -> ?inco:int * int -> 'l complex_array
      -> c2c plan
      (** See {!Fftw3.D.Genarray.dft}. *)

    val r2r : r2r_kind * r2r_kind -> ?meas:measure -> ?normalize:bool ->
      ?preserve_input:bool -> ?unaligned:bool ->
      ?n:int * int -> ?howmany_n:int array ->
      ?howmanyi:(int * int) list ->
      ?ofsi:int * int -> ?inci:int * int -> 'l float_array ->
      ?howmanyo:(int * int) list ->
      ?ofso:int * int -> ?inco:int * int -> 'l float_array
      -> r2r plan
      (** See {!Fftw3.D.Genarray.r2r}. *)
  end


  module Array3 :
  sig
    val create: ('a, 'b) kind -> 'c layout -> int -> int -> int
      -> ('a, 'b, 'c) Array3.t
      (** See {!Fftw3.D.Genarray.create}. *)

    type 'l complex_array = (Complex.t, complex_elt, 'l) Array3.t
        (** Double precision complex 3D array. *)
    type 'l float_array   = (float, float_elt, 'l) Array3.t
        (** Double precision float 3D array. *)

    val dft : dir ->
      ?meas:measure -> ?normalize:bool ->
      ?preserve_input:bool -> ?unaligned:bool ->
      ?n: int * int * int -> ?howmany_n: int array ->
      ?howmanyi:(int * int * int) list ->
      ?ofsi:int * int * int -> ?inci:int * int * int -> 'l complex_array ->
      ?howmanyo:(int * int * int) list ->
      ?ofso:int * int * int -> ?inco:int * int * int -> 'l complex_array
      -> c2c plan
      (** See {!Fftw3.D.Genarray.dft}. *)
  end
end


(** Double precision FFTW. *)
module D : Sig
  with type float_elt = Bigarray.float64_elt
  and type complex_elt = Bigarray.complex64_elt


IFDEF FFTW3F_EXIST THEN
(** Single precision FFTW.  This is only available if the single
    precision FFTW3 library was available when this module was
    compiled. *)
module S : Sig
  with type float_elt = Bigarray.float32_elt
  and type complex_elt = Bigarray.complex32_elt
ENDIF

(** Managing wisdom.  Save and restore plans to/from disk or other
    media. *)
module Wisdom :
sig

  val export : (char -> unit) -> unit
    (** [Wisdom.export write] exports the current wisdom to any medium,
        as specified by the callback function [write].

        This function is not thread safe.
    *)

  val to_file : string -> unit
    (** [Wisdom.to_file fname] writes the current wisdom to the file
        [fname].
         *)
    (* FIXME: @param mode The file is created with permissions [mode].
        Default: [0o644].  *)

  val to_string : unit -> string
    (** [Wisdom.to_string()] exports the current wisdom as a string. *)


  val import : (unit -> char) -> unit
    (** [Wisdom.import read] imports wisdom from any input medium, as
        specified by the callback function [read].  If the end of the
        input data is reached (which should never happen for valid
        data), [read] should raise [End_of_file].  The imported wisdom
        {i replaces} any wisdom accumulated by the running program.

        This function is not thread safe.

        @raise Failure if the wisdom was not successfully read. *)
  val from_file : string -> unit
    (** [Widsom.from_file fname] replace the current wisdom with the
        one read from the file [fname].

        @raise Failure if the wisdom was not successfully read. *)
  val from_string : string -> unit
    (** [Wisdom.from_string s] replace the current wisdom whith the
        one read from [s].

        @raise Failure if the wisdom was not successfully read. *)
  val from_system : unit -> unit
    (** [Wisdom.from_system()] replace the current wisdom with one read
        from an implementation-defined standard file (e.g. /etc/fftw/wisdom).

        @raise Failure if the wisdom was not successfully read. *)


  val forget : unit -> unit
    (** [Wisdom.forget()] causes all accumulated wisdom to be
        discarded.  (New wisdom can be gathered subsequently
        though.) *)
end
