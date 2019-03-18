PKGVERSION = $(shell git describe --always --dirty)

all build byte native:
	dune build @install @examples

test runtest:
	dune runtest --force

install uninstall:
	dune $@

doc: build
	sed -e 's/%%VERSION%%/$(PKGVERSION)/' src/fftw3.mli \
	  > _build/default/src/fftw3.mli
	dune build @doc

# Profiling
pnc:
#	$(MAKE) -C src -B pnc
	$(MAKE) -C examples -B pnc

lint:
	opam lint fftw3.opam

clean::
	dune clean

.PHONY: all byte native test runtest install uninstall doc pnc clean lint
