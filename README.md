Interface to FFTW, version 3
============================

Perquisites
-----------

The FFTW, version 3, with its development files (fftw3.h) must be
installed on your platform.


Compilation
-----------

    ./bootstrap (optional, requires autoconf)
    ./configure
    make

If your fftw3 header files or libraries are not where you C compiler
expects to find them, you need to tell configure by using the
following syntax:

    CPPFLAGS="-I/opt/local/include" LDFLAGS="-L/opt/local/lib" ./configure

Of course, replace /opt/local/include (resp. /opt/local/lib) by the
actual paths of your header files (resp. library).


Usage
-----

Fftw3 contains two submodules differing only by precision Fftw3.D for
double precision and Fftw3.S for single precision.  Note that the
single precision module will only exist if the corresponding FFTW3
library was discovered by the configure script.  See the mli file for
more details.


Examples
--------

Some examples need additional libraries:
- [Lacaml](https://bitbucket.org/mmottl/lacaml)
- [Gnuplot](http://sourceforge.net/projects/ocaml-gnuplot)

If these are not detected by the configure script, the corresponding
examples will not be compiled or will have less features.


Bugs
----

To report bugs, use the Github
[issues page](https://github.com/Chris00/fftw-ocaml/issues).


Development
-----------

If you want to participate to the development of Fftw3, make a
fork of the [project](https://github.com/Chris00/fftw-ocaml) and
submit a pull request.


Licence
-------

FFtw3 is released under the LGPL with the special axception of the
standard library.  See the file LICENSE for more details.

Although not required, it will be appreciated if your product mentions
it uses these bindings.
