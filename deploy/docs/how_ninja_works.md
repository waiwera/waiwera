In Waiwera directory, type `python config.py` to configure the Meson
build and config/build PETSc.

Then `ninja -C build` to build Waiwera (and FSON and Zofu as subprojects
if these are not already installed on the system).

Unit tests can be run with `python unit_tests.py`.

Install Waiwera with `ninja -C build install`.

By default this will do a release build and install to bin/ in the
user's home dir. These and other things can be overridden with
parameters to the `config.py` script.
