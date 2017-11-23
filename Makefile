
PKGVERSION = $(shell git describe --always --dirty)

all build byte native:
	jbuilder build @install --dev
	jbuilder build @examples

test runtest:
# Force the tests to be run
	$(RM) -rf _build/default/tests/
	jbuilder runtest

install uninstall:
	jbuilder $@

doc: build
	sed -e 's/%%VERSION%%/$(PKGVERSION)/' src/fftw3.mli \
	  > _build/default/src/fftw3.mli
	jbuilder build @doc

# Profiling
pnc:
#	$(MAKE) -C src -B pnc
	$(MAKE) -C examples -B pnc

lint:
	opam lint fftw3.opam

clean::
	jbuilder clean

.PHONY: all byte native test runtest install uninstall doc pnc clean lint
