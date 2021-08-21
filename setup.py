import setuptools
import os, sys
sys.path.append(os.getcwd())
import mpdb

with open("README.md", "r") as fh:
    long_description = fh.read()

try:
    from Cython.Build import cythonize
    cem = cythonize("mpdb/eq/*.pyx")
except:
    cem = []
    
setuptools.setup(
    author="Stuart Nolan",
    author_email="61199416+stuart-nolan@users.noreply.github.com",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: POSIX :: Linux",
    ],
    description="material property data base",
    include_package_data=True,
    packages=setuptools.find_packages(),
    #packages = ['.'],
    install_requires=[
        'markdown',
        'appdirs',
        'tabulate',],
    extras_require={
        "cython": ["Cython",],
    },
    ext_modules=cem,
    long_description=long_description,
    long_description_content_type="text/markdown",
    name="mpdb", 
    url="https://github.com/stuart-nolan/mpdb.git",
    version=mpdb.__version__,
    zip_safe = False
)
"""
https://acp.copernicus.org/articles/15/4399/2015/acp-15-4399-2015-supplement.zip
"""
