===========
Basic usage
===========

``fanc`` uses subcommands to run all of its analyses. The ``fanc`` command itself can
be used to get an overview of available subcommands, to print the current FAN-C version,
or to set logging and notification parameters that affect all subcommands.

********
Overview
********

.. argparse::
   :module: fanc.commands.fanc_commands
   :func: fanc_parser
   :prog: fanc
   :nodescription:
   :nodefault:


*******
Logging
*******

You can set the verbosity level of any ``fanc`` subcommand with the ``-v`` option. Use
more or less v's for more or less logging output. The default is ``-vv``, which
corresponds to error, warning, and info messages. ``-vvv`` also displays debug messages,
which might be helpful to identify issues with an analysis. ``-v`` only displays error
and warning messages. To disable logging completely, use the ``-s`` option.

By default, logging output is sent to ``stderr``. You can redirect log messages to a file
using ``-l <file name>``.


*******************
Email notifications
*******************

Sometimes it is convenient to be notified by email if a ``fanc`` command finishes,
especially for long-running commands such as ``fanc pairs`` or ``fanc map``. You can
instruct ``fanc`` to send an email when a command finished using ``-m <email address>``.
You must also specify the SMTP server settings using the options
``--smtp-<server|username|password|sender-address>``, or you can pre-configure these using
the :ref:`fanc-config-files`.


***************
Temporary files
***************

Many of the more computationally intensive FAN-C commands support the :code:`-tmp`
argument. This instructs the command to copy all input files to a temporary directory
before processing. Similarly, output files will intitially generated in the temporary
directory and only copied to their intended output locations once the command completes.

This can be very effective when your data is located, for example, on a network file
system or a slow external HDD, while your local machine or computing node has access
to a fast SSD. Using :code:`-tmp`, and assuming the local machine's default temporary
directory resides on an SSD, files are copied from their original location to the SSD
at the start of the command, thus avoiding the slow file system access throughout the
remainder of the processing steps.

You can change the default temporary directory by setting the :code:`TMPDIR` environment
variable on your local system to a folder of your choice.


.. _fanc-config-files:

******************
FAN-C config files
******************

FAN-C supports configuration files, which can be located (in descending order of priority)

- in the current directory, named ``fanc.conf``
- in a path specified by the Unix environment variable ``FANC_CONF``
- in the user's home folder (named ``fanc.conf`` or ``.fanc.conf``)
- in the ``.config`` folder in a user's home directory, called ``fanc.conf``
- in ``/etc/fanc/fanc.conf``

Settings made in one config file are overridden by settings in a file with higher priority.

You can write the default config file to a location of your choice using ``fanc write-config``.

.. argparse::
   :module: fanc.commands.fanc_commands
   :func: write_config_parser
   :prog: fanc write-config
   :nodescription:
   :nodefault:

An explanation of the different settings can be found as comments in the default config file.
The file is written in `YAML <https://yaml.org/>`_


.. _fanc-numexpr:

********************************
NumExpr ThreadPool configuration
********************************

FAN-C uses PyTables for fast querying of most of its storage classes. `Condition-based queries
in PyTables <https://www.pytables.org/usersguide/libref/structured_storage.html#tables.Table.where>`_,
which are used, for example, to find regions and pixels in certain matrix subsets,
rely on the NumExpr package. NumExpr can be multi-threaded, and FAN-C uses the default NumExpr
ThreadPool configuration (typically 8 threads). There is generally no need to change this
preset, but if you want to optimise every single aspect of your pipeline, you may want to
take a look at
`this NumExpr help page
<https://numexpr.readthedocs.io/projects/NumExpr3/en/latest/user_guide.html#threadpool-configuration>`_
.
