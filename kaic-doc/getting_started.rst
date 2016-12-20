===============
Getting started
===============

.. contents::
   :depth: 2

Kai-C offers access to its Hi-C processing pipeline on multiple levels, including
a high-level executable and a low-level Python 2.7.x / 3.5.x API. Often the ``kaic``
executable will cover the analysis needs of basic users. For more advanced use cases, the API
is available, which gives maximum flexibility with regards to
integration of other data sets and custom analysis methods. The output of the command line tool
is fully compatible with the API.

Installation
~~~~~~~~~~~~

Before installing Kai-C, make sure you have all the prerequisites installed in your system.
Specifically, Kai-C uses the HDF5 (via PyTables), a file format that simplifies the storage and access to huge
amounts of data. The minimum required version is 1.8.4, but we recommend installing the latest version.

Prerequisite: HDF5
__________________

If you are on Linux, download the source code of the latest version from
the `HDF5 website <https://www.hdfgroup.org/HDF5/>`_ and unpack it.
There are two different HDF5 branches - note that the 1.10.x branch is not supported!

.. code:: bash

   # make a new directory
   mkdir hdf5-build
   cd hdf5-build
   # replace xx with current version number
   wget https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.xx.tar.gz
   # unpack
   tar xzf hdf5-1.8.xx.tar.gz
   cd hdf5-1.8.xx/
   # use --prefix to set the folder in which HDF5 should be installed
   # alternatively, you can omit --prefix=... here and run
   # sudo make install to install globally (requires admin rights)
   ./configure --prefix=/path/to/hdf5/dir
   make
   make install

If you are on OS X or macOS, we highly recommend using the fantastic `Homebrew <http://brew.sh/>`_.
Then you can simply:

.. code:: bash

   brew tap homebrew/science
   brew install hdf5

To ensure that PyTables, the Python library that uses HDF5, finds the correct HDF5 version, we
need to set an environment variable pointing to the installation directory:

.. code:: bash

   # on Linux, this is the same /path/to/hdf5/dir that you used above
   # on macOS, this is wherever brew installs its recipies (typically /usr/local/Cellar)
   #   - you can find this by running 'brew info hdf5'
   export HDF5_DIR=/path/to/hdf5/dir


Kai-C
_____

The simplest way to install Kai-C is via pip:

.. code:: bash

   pip install kaic

and that should be all you need! If you are not the owner of the Python installation,
try:

.. code:: bash

   pip install --user kaic

You can also directly download the Kai-C source code from Github by cloning its repository.
The installation is then done via setup.py:

.. code:: bash

   git clone http://www.github.com/vaquerizaslab/kaic
   cd kaic
   python setup.py install


Overview
~~~~~~~~

Kai-C can be accessed via command line (``kaic`` for analysis, ``klot`` for plotting) or as a Python 2.7.x / 3.5.x
module (``import kaic``).