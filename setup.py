import psyco
psyco.full()
from distutils.core import setup
import sys, os, py2exe

# data_files:   List of datafiles that need to be explicitly included
#               e.g., [('target/relative/path', ['file','list'])]
# excludes:     List of modules to exclude (takes precedence over includes)
# dll_excludes: List of dlls to exclude
# includes:     List of modules to include
# packages:     List of packages to include
data_files = []
excludes = [ 'numarray', 'Numeric', 'Tkinter', 'email','wxPython','OpenGL',"enthought","vtk",'vpython']
dll_excludes = [] # dlls to exclude
includes = [ 'matplotlib', 'wx', 'numpy', ]
packages = []  # May want 'wx' in packages if creating a generic environment

# Explicit package rules
if 'matplotlib' in includes:
    import matplotlib
    #matplotlib.use( 'WXAgg',warn=False )
    data_files += matplotlib.get_py2exe_datafiles()
    includes += ['matplotlib.numerix.'+pkg for pkg in
                 ['ma','mlab','fft','linear_algebra','npyma','random_array']]
    includes += ['matplotlib.backends.backend_'+pkg for pkg in
                 ['wx','wxagg','pdf','ps','svg','agg','emf']]

# Don't include anything explicitly excluded
includes = [x for x in includes if x not in excludes]

# Include the C/C++ runtime (there should be a way of putting these beside
# the python dll in the executable, but for now I put them in the executable
# directory. Note: the C++ runtime might be in the wx package instead of the
# python directory.
pythonpath = os.path.dirname(sys.executable)
c_runtime = pythonpath+"/MSVCR71.DLL"
cpp_runtime = pythonpath+"/MSVCP71.DLL"
data_files += [('.',[c_runtime,cpp_runtime])]
               #C:\Python25\code\huckel\data
root = 'C:\git/orbis'
atomic_data = root+'/data/atomic_data.yaml'
bond_data = root+'/data/bond_data.yaml'
def_atomic_data = root+'/data/default/atomic_data.yaml'
def_bond_data = root+'/data/default/bond_data.yaml'
data_files += [('./data',[atomic_data,bond_data]),('./data/default',[def_atomic_data,def_bond_data])]

doc_dir = root+'/docs'
doc_img = doc_dir+'/images'
images = [doc_img+'/'+x for x in os.listdir(doc_img) if x.find('.png')>=0]
docs = [doc_dir+'/'+x for x in os.listdir(doc_dir) if x.find('.css')>=0 or x.find('.html')>=0]
data_files+=[('./docs',docs),('./docs/images',images)]

image_dir = root+'/images'
image_files = [image_dir + x for x in ['/new.png','/open.png','/save.png','/icon.bmp']]
data_files += [('./images',image_files)]
data_files += [('.',[root+'/license.txt'])]
# End of configuration


py2exe_opts = dict(includes=includes, excludes=excludes,
                   skip_archive=0, compressed=1,
                   dll_excludes=dll_excludes, packages=packages,
                   bundle_files=1, optimize=2)

manifest = """
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<assembly xmlns="urn:schemas-microsoft-com:asm.v1"
manifestVersion="1.0">
<assemblyIdentity
    version="0.64.1.0"
    processorArchitecture="x86"
    name="Controls"
    type="win32"
/>
<description>Orbis - Simple Huckel Solver: Beta v0.1.5</description>
<dependency>
    <dependentAssembly>
        <assemblyIdentity
            type="win32"
            name="Microsoft.Windows.Common-Controls"
            version="6.0.0.0"
            processorArchitecture="X86"
            publicKeyToken="6595b64144ccf1df"
            language="*"
        />
    </dependentAssembly>
</dependency>
</assembly>
"""

setup(
    options = dict(py2exe=py2exe_opts),
    windows = [
        {
            "script": "orbis.py",
            "icon_resources": [(1, "icon.ico")],
            "other_resources": [(24,1,manifest)]
        }
    ],

    data_files = data_files,
    #zipfile = None,  # Bundle library.zip with the executable
    )

