# Hacking the build into Travis-CI "C" environment
# Inspired by
#     http://blog.mlin.net/2013/02/testing-ocaml-projects-on-travis-ci.html

DEBS=libfftw3-dev
# OPAM version to install:
export OPAM_VERSION=1.0.0
# OPAM packages needed to build tests:
export OPAM_PACKAGES='ocamlfind'

# Install OCaml
sudo apt-get update -q -y
sudo apt-get install -q -y ocaml ${DEBS}

# install opam
curl -L https://github.com/OCamlPro/opam/archive/${OPAM_VERSION}.tar.gz \
    | tar xz -C /tmp
pushd /tmp/opam-${OPAM_VERSION}
./configure
make
sudo make install
opam init --auto-setup
eval `opam config -env`
popd

# install packages from opam
opam install -q -y ${OPAM_PACKAGES}


# compile & run tests (here assuming OASIS DevFiles)
./bootstrap && ./configure && make
