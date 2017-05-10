"""Analyses Fortan 90 source code module dependencies."""
from builtins import str
from builtins import zip
from builtins import object

def comment_line(line, i):
    """Returns True if line has a comment marker before the specified index. """
    return '!' in line[:i]

def line_unit_name(line, unit = 'module'):
    """Returns unit (module or program) name in line, or None if there
    isn't one."""
    name = None
    if line:
        lowerline = line.lower()
        imod = lowerline.find(unit)
        if imod >= 0:
            if not comment_line(lowerline, imod):
                if 'end' not in lowerline[:imod]:
                    items = line[imod:].strip().split()
                    if items[1].lower() != 'procedure':
                        name = items[1]
    return name

def change_ext(filename, ext):
    """Changes filename extension.."""
    from os.path import splitext
    base, oldext = splitext(filename)
    return base + ext

class module(object):
    """Fortran code module."""
    def __init__(self, f):
        self.name = None
        self.program = False
        self.depends = []
        self.used = []
        self.parse(f)

    def __repr__(self): return self.name + ':' + str(self.depends)

    def parse(self, f):
        """Parses module."""
        self.parse_module_name(f)
        self.parse_depends(f)

    def parse_module_name(self, f):
        """Parses module name from file f."""
        line, name, prog = f.readline(), None, False
        while line and name is None:
            name = line_unit_name(line, 'module')
            if name is None:
                name = line_unit_name(line, 'program')
                if name is None: line = f.readline()
                else: prog = True
        self.name = name
        self.filename = f.name
        self.program = prog

    def parse_depends(self, f):
        """Parses dependencies in module."""
        endstr = 'end'
        if self.program: endstr += ' program'
        else: endstr += ' module'
        line = True
        depends = set([])
        while line:
            line = f.readline()
            if line:
                lowerline = line.lower()
                if endstr in lowerline: line = None
                else:
                    iuse = lowerline.find('use ')
                    if iuse >= 0:
                        if not comment_line(lowerline, iuse):
                            name = line[iuse:].strip().split()[1]
                            name = name.split(',')[0] # for 'only'
                            depends.add(name)
        self.depends = list(depends)
        self.depends.sort()

class sourcefile(object):
    """Source code file."""
    def __init__(self, filename):
        from os.path import split
        path, name = split(filename)
        self.filename = name
        self.path = path
        self.modules = []
        self._file = open(filename)
        self.parse()
        self._file.close()

    def __repr__(self):
        return self.filename + ':' + str(self.modules)

    def parse(self):
        """Parses modules in source file."""
        done = False
        while not done:
            m = module(self._file)
            if m.name is None: done = True
            else: self.modules.append(m)

    def get_depends(self):
        """Returns list of all dependencies in the source file,
        from all modules."""
        deps = set([])
        for m in self.modules: deps = deps.union((set(m.depends)))
        deplist = list(deps)
        deplist.sort()
        return deplist
    depends = property(get_depends)

class dependencies(object):
    """Source code dependencies object."""
    def __init__(self, srcdir = './', ext = '.F90'):
        if isinstance(srcdir, str): self.srcdirs = [srcdir]
        else: self.srcdirs = srcdir
        self.ext = ext
        self.parse_sources()

    def parse_sources(self):
        """Parses source files."""
        from glob import glob
        from os.path import sep
        self.sourcefiles = {}
        self.modules = {}
        for srcdir in self.srcdirs:
            sourcenames = glob(srcdir + sep + '*' + self.ext)
            for filename in sourcenames:
                s = sourcefile(filename)
                self.sourcefiles[s.filename] = s
                for m in s.modules: self.modules[m.name] = m
        for mname in self.modules:
            m = self.modules[mname]
            for d in m.depends:
                if d in self.modules: self.modules[d].used.append(mname)

    def make_depends(self, objdirs = [], obj = '$(OBJ)', variables = [],
                     filename = 'depends.in', omit = []):
        """Writes file with all source file dependencies written
        one per line, suitable for an include file in a makefile.
        objdirs specifies the object directories corresponding to the
        source directories. obj is the string to be used for object file
        extensions. variables contains pairs of strings with variable
        definitions declared in the makefile, to be substituted in the
        dependency rules. filename is the name of the output include file.
        omit is a list of targets to omit from the output."""
        from os.path import split, sep
        if objdirs == []: objdirs = self.srcdirs
        objdirdict = dict(list(zip(self.srcdirs, objdirs)))
        if len(variables) > 0:
            def subst_vars(f):
                newf = f
                for (k,r) in variables:
                    newf = newf.replace(r, k)
                return newf
        else:
            def subst_vars(f): return f
        outfile = open(filename, 'w')
        lines = []
        for sourcename in self.sourcefiles:
            s = self.sourcefiles[sourcename]
            if len(s.depends) > 0:
                path = objdirdict[s.path] + sep
                fname = subst_vars(s.filename)
                objname = path + change_ext(fname, obj)
                if objname not in omit:
                    depends = []
                    for module_name in s.depends:
                        if module_name in self.modules:
                            dependfile = self.modules[module_name].filename
                            dependpath, fname = split(dependfile)
                            fname = subst_vars(fname)
                            dependobjdir = objdirdict[dependpath]
                            depend_objname = dependobjdir + sep + \
                                change_ext(fname, obj)
                            if depend_objname != objname:
                                depends.append(depend_objname)
                    if depends:
                        lines.append(objname + ":" + " ".join(depends) + '\n')
        for line in sorted(lines): outfile.write(line)
        outfile.close()

    def write_module_dot(self, path = './', filename = 'depends.dot',
                         exclude_modules = []):
        """Writes DOT file describing the module dependency graph.
        This can be used to produce a module dependency diagram
        using e.g. Graphviz."""
        outfile = open(filename, 'w')
        lines = []
        top = 'depends'
        for sourcename in self.sourcefiles:
            s = self.sourcefiles[sourcename]
            if s.path == path:
                for m in s.modules:
                    if len(m.used) == 0 and len(m.depends) > 0:
                        top = m.name
        outfile.write('digraph ' + top + '{\n')

        for sourcename in self.sourcefiles:
            s = self.sourcefiles[sourcename]
            if s.path == path:
                for m in s.modules:
                    for d in m.depends:
                        if d not in exclude_modules:
                            lines.append('        ' + m.name + ' -> ' + d + '\n')
        for line in sorted(lines): outfile.write(line)
        outfile.write('}')
        outfile.close()

if __name__ == '__main__':

    srcdirs = ['src', 'test/unit/src']
    objdirs = ['$(BUILD)', '$(TEST)/$(BUILD)']
    variables = [('$(TESTSUF)', '_test'), ('$(PROG)', 'waiwera'),
             ('$(TESTPROG)', 'test_all')]
    omit = ['$(BUILD)/$(PROG)$(OBJ)', '$(TEST)/$(BUILD)/$(TESTPROG)$(OBJ)']
    exclude_modules = ['petsc', 'petscvec', 'petscsnes', 'petscis', 'petscsys',
                       'kinds_module', 'fson_value_m', 'fson',
                       'iso_fortran_env']

    deps = dependencies(srcdirs)
    deps.make_depends(objdirs, variables = variables, omit = omit)
    deps.write_module_dot(srcdirs[0], 'doc/depends.dot', exclude_modules)
