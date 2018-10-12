# Additional configuration of FRUIT Fortran unit testing library

import os

def get_fruit_version(base_path):
    """Returns version string for FRUIT library"""

    changes_filename = os.path.join(base_path, "CHANGES.txt")
    changes = open(changes_filename)
    version = None
    for line in changes:
        if line.startswith("FRUIT"):
            version = line.split()[-1]
            break
    changes.close()
    return version

def write_fruit_pkgconfig_file(base_path):
    """Write pkg-config file for FRUIT shared library, for specified FRUIT base path"""

    filename = os.path.join(base_path, "pkgconfig", "FRUIT.pc")

    f = open(filename, 'w')
    f.write("prefix=%s\n" % base_path)
    f.write("includedir=${prefix}\n")
    f.write("libdir=${prefix}\n\n")

    f.write("Name: FRUIT\n")
    f.write("Description: Fortran unit test framework\n")
    version = get_fruit_version(base_path)
    if version is not None: f.write("Version: " + version + "\n")
    f.write("Cflags: -I${includedir}\n")
    f.write("Libs: -L${libdir} -lfruit\n")

    f.close()
    
