
PKGVERSION = $(shell git describe --always --dirty)

all build byte native:
	jbuilder build @install #--dev
	jbuilder build @examples

test runtest:
# Force the tests to be run
	$(RM) -rf _build/default/tests/
	jbuilder runtest

install uninstall:
	jbuilder $@

doc:
	sed -e 's/%%VERSION%%/$(PKGVERSION)/' src/Root1D.mli \
	  > _build/default/src/Root1D.mli
	jbuilder build @doc
	echo '.def { background: #f9f9de; }' >> _build/default/_doc/odoc.css

# Profiling
pnc:
#	$(MAKE) -C src -B pnc
	$(MAKE) -C examples -B pnc

lint:
	opam lint fftw3.opam

clean::
	jbuilder clean

.PHONY: all byte native test runtest install uninstall doc pnc clean lint
