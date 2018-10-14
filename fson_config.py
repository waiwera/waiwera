# Additional configuration of FSON JSON library

import os

def write_fson_pkgconfig_file(base_path):
    """Write pkg-config file for FSON shared library, for specified FSON base path.
    FSON does not have version numbers so a version of 1.0.0 is assumed."""

    config_dir = os.path.join(base_path, "dist", "pkgconfig")
    if not os.path.isdir(config_dir): os.mkdir(config_dir)

    version = "1.0.0"

    filename = os.path.join(config_dir, "FSON.pc")
    f = open(filename, 'w')
    f.write("prefix=%s\n" % base_path)
    f.write("includedir=${prefix}/build\n")
    f.write("libdir=${prefix}/dist\n\n")

    f.write("Name: FSON\n")
    f.write("Description: Fortran 95 JSON parser library\n")
    f.write("Version: " + version + "\n")
    f.write("Cflags: -I${includedir}\n")
    f.write("Libs: -L${libdir} -lfson\n")

    f.close()
    
