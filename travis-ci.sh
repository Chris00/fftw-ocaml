# Hacking the build into Travis-CI "C" environment
# See http://anil.recoil.org/2013/09/30/travis-and-ocaml.html

DEBS="libfftw3-dev liblapack-dev"
OPAM_PACKAGES='ocamlfind lacaml archimedes'

# Install OCaml
echo "yes" | sudo add-apt-repository ppa:avsm/ppa
sudo apt-get update -qq
sudo apt-get install --force-yes ocaml ocaml-native-compilers opam ${DEBS}

export OPAMYES=1
opam init
eval `opam config env`
opam install -q -y ${OPAM_PACKAGES}


# compile & run tests
./bootstrap && ./configure && make
