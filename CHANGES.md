0.8.5
-----

- Update to `dune-configurator`.

0.8.4 2019-12-30
----------------

- Compatibility with OCaml 4.10.

0.8.3 2019-03-19
----------------

- Compatibility with OCaml 4.08 and the newer `Bigarray` code.
- Upgrade to OPAM 2.
- Remove conditionals to support older versions of OCaml.

0.8.2 2018-09-05
----------------

- Fix enabling the single precision module `Fftw3.S`.
- Enable robustness tests at compile time (thanks to Dune.configurator).
- Switch to Dune.
- Rely on a `conf-fftw3` package for the presence of the C FFTW3 library.

0.8.1 2018-05-22
----------------

- Fix documentation (issue #11).
- Add bigarray dependency.

0.8 2017-11-23
--------------

- Conditionally define functions exported in `bigarray.h` in OCaml
  4.06.0.
- Port to `jbuilder` and `topkg`.
- Use [@@noalloc] when possible.
