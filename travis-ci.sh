# Hacking the build into Travis-CI "C" environment
# See http://anil.recoil.org/2013/09/30/travis-and-ocaml.html

OPAM_PACKAGES='ocamlfind lacaml archimedes'

export OPAMYES=1
opam init
eval `opam config env`
opam install -q -y ${OPAM_PACKAGES}


# compile & run tests
./bootstrap && ./configure && make
