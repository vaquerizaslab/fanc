===========
Basic usage
===========

``kaic`` uses subcommands to run all of its analyses. The ``kaic`` command itself can
be used to get an overview of available subcommands, to print the current Kai-C version,
or to set logging and notification parameters that affect all subcommands.

********
Overview
********

.. argparse::
   :module: kaic.commands.kaic_commands
   :func: kaic_parser
   :prog: kaic
   :nodescription:


*******
Logging
*******

You can set the verbosity level of any ``kaic`` subcommand with the ``-v`` option. Use
more or less ``v``s for more or less logging output. The default is ``-vv``, which
corresponds to error, warning, and info messages. ``-vvv`` also displays debug messages,
which might be helpful to identify issues with an analysis. ``-v`` only displays error
and warning messages. To disable logging completely, use the ``-s`` option.

By default, logging outout is sent to ``stderr``. You can redirect log messages to a file
using ``-l <file name>``.


*******************
Email notifications
*******************

Sometimes it is convenient to be notified by email if a ``kaic`` command finishes,
especially for long-running commands such as ``kaic pairs`` or ``kaic map``. You can
instruct ``kaic`` to send an email when a command finished using ``-m <email address>``.
You must also specify the SMTP server settings using the options
``--smtp-<server|username|password|sender-address>``, or you can pre-configure these using
the :ref:`kaic-config-files`.


.. _kaic-config-files:

******************
Kai-C config files
******************

Kai-C supports configuration files, which can be located (in descending order of priority)

- in the current directory, named ``kaic.conf``
- in a path specified by the Unix environment variable ``KAIC_CONF``
- in the user's home folder (named ``kaic.conf`` or ``.kaic.conf``)
- in the ``.config`` folder in a user's home directory, called ``kaic.conf``
- in ``/etc/kaic/kaic.conf``

Settings made in one config file are overridden by settings in a file with higher priority.

You can write the default config file to a location of your choice using ``kaic write-config``.

.. argparse::
   :module: kaic.commands.kaic_commands
   :func: write_config_parser
   :prog: kaic write-config
   :nodescription:

An explanation of the different settings can be found as comments in the default config file.
The file is written in `YAML <https://yaml.org/>`_
