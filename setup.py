import os
from numpy.distutils.core import setup, Extension

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

humidity = Extension(name = 'interimR.humidity_toolbox',
                sources = ['interimR/f90/humidity_toolbox.f90'])

setup(
    name = "interimR",
    version = "2.0",
    author = "Raphael Dussin",
    author_email = "raphael.dussin@gmail.com",
    description = ("A package to create forcing for ocean models " ),
    license = "DTC",
    keywords = "ocean forcing",
    url = "",
    packages=['interimR'],
    ext_modules = [humidity],
    scripts = ['interimR/scripts/process_ERAinterim']
)
#    long_description=read('README'),
#    classifiers=[
#        "Development Status :: 3 - Alpha",
#        "Topic :: Utilities",
#        "License :: OSI Approved :: BSD License",
#    ],
#)

os.system('./interimR/post/postinstall.bash')

