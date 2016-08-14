===============
Getting started
===============

.. contents::
   :depth: 2

Kai-C offers access to its Hi-C processing pipeline on multiple levels, including
a high-level executable and a low-level Python 2.7 API. Often the ``kaic``
executable will cover the analysis needs of basic users, but its output is
compatible with the API, which gives maximum flexibility with regards to
integration of other data sets and custom analysis methods.

Installation
------------

You can download the Kai-C source code from GitLab by cloning its repository.

.. code:: bash

  git clone http://www.github.com/vaquerizaslab/kaic

kaic depends on several python packages which should mostly be installed
automatically during the installation process.

However, one core dependency is not installed automatically:
`PyTables <https://github.com/PyTables/PyTables>`_ is currently on version 3.2.2
and there is a critical bug affecting kaic performance that will be fixed in
version 3.3. To circumvent the bug, it is highly recommended to install the
PyTables development version. (Note: The kaic installer will now attempt to install
this directly).

PyTables depends on the `HDF5 <https://www.hdfgroup.org/HDF5/>`_ library. To
install the latest version it is easiest to use Homebrew (on OS X) or Linuxbrew.

.. code:: bash

  brew install hdf5

Or you can install
`from source <https://www.hdfgroup.org/HDF5/release/obtain5.html>`_.

To ensure this version of HDF5 is used during the installation of PyTables, set
the HDF5_DIR environment variable to the path of your HDF5 installation.

Then install PyTables:

.. code:: bash

  pip install git+https://github.com/PyTables/PyTables.git

or if you do not have root access to the Python installation:

.. code:: bash

  pip install --user git+https://github.com/PyTables/PyTables.git

Finally, from the root directory of kaic do

.. code:: bash

  pip install .

or, if you don't have permissions to install packages globally:

.. code:: bash

  pip install --user .



Overview
--------

Kai-C can be accessed via command line (``kaic`` for analysis, ``klot`` for plotting) or as a Python 2.7.x module (
``import kaic``.