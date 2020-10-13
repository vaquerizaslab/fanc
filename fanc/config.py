"""
Module that handles FAN-C configuration files.

Generally, some FAN-C behaviour can be controlled using a "fanc.conf"
file in one of the following locations, in order of preference:

1. fanc.conf in current directory
2. File specified by environment variable fanc_CONF
3. fanc.conf in home folder
4. .fanc.conf in home folder
5. ~/.config/fanc.conf
6. /etc/fanc/fanc.conf

All configuration files will be read, but the options of lower-ranking
locations will be overwritten.

To write a default configuration file use  :code:`fanc write-config`
from the command line or

.. code::

    from fanc.config import write_default_config

    write_default_config('~/.fanc.conf')  # or any other location

"""

from fanc.tools.general import Map
import yaml
import os
import copy
import warnings
import logging
logger = logging.getLogger(__name__)

_default_config_content = """\
#
# MEMORY
#
edge_buffer_size: 3G

#
# HDF5
#
hdf5_compression_level: 1
hdf5_compression_algorithm: blosc

#
# Mapping
#
bowtie2_path: bowtie2
bwa_path: bwa

#
# PLOTTING
#
colormap_hic: germany
adjustment_slider_height: .2
pad_empty_axis: .1
pad_with_label: .2
pad_with_ticks: .1
pad_next_title: .2
pad_with_tick_legend: .1
pdf_font_as_text: False

#
# EMAIL
#
email_from_address:
email_smtp_server:
email_smtp_port:
email_smtp_username:
email_smtp_password:

#
# ERROR HANDLING
#
trap_signals:
    - SIGUSR1
    - SIGUSR2
raise_exception_on_trapped_signal: True

#
# Progressbars
#
hide_progressbars: False

#
# Sun/Oracle Grid Engine
#
gridmap_tmpdir:

sge_qsub_path: qsub
sge_qdel_path: qdel
sge_default_queue: all.q
sge_log_dir:
sge_parallel_environment: smp
sge_nodes:
sge_qsub_options: "-V -notify"
sge_shell: /bin/bash

#
# SLURM
#

slurm_sbatch_path: sbatch
slurm_sbatch_options: ""
slurm_scancel_path: scancel
slurm_log_dir:
slurm_shell: /bin/bash

#
# Juicer
#
juicer_tools_jar_path: 
"""
default_config = yaml.load(_default_config_content, Loader=yaml.FullLoader)

_config_file_name = 'fanc.conf'
_environment_variable_name = 'FANC_CONF'

# check common places for config file
# modified from: http://stackoverflow.com/questions/7567642/where-to-put-a-configuration-file-in-python
# the order of files in this tuple determines config file priority
# options in files with higher priority override options in files with lower priority
# default options have lowest priority
_config_file_locations = (
    # current directory
    os.path.join(os.curdir, _config_file_name),
    # environment variable
    os.environ.get(_environment_variable_name),
    # home folder
    os.path.join(os.path.expanduser("~"), _config_file_name),
    # home folder (hidden file)
    os.path.join(os.path.expanduser("~"), '.' + _config_file_name),
    # home folder (.config directory)
    os.path.join(os.path.expanduser("~"), '.config', 'fanc', _config_file_name),
    # /etc/fanc (MySQL style)
    os.path.join("/etc/fanc/", _config_file_name)
)


# read the first config file that is found, checking in the order above
def read_file_configs(config_file_locations=None):
    """
    Read all fanc.conf configuration files found in
    1. current directory
    2. specified by environment variable fanc_CONF
    3. fanc.conf in home folder
    4. .fanc.conf in home folder
    5. ~/.config/fanc.conf
    6. /etc/fanc/fanc.conf

    :param config_file_locations: optional custom config file paths in a list
    :return: dictionaries of parsed configs
    """

    if config_file_locations is None:
        config_file_locations = _config_file_locations

    configs = []
    for location in config_file_locations:
        if location is None:
            continue

        config_path = os.path.expanduser(location)
        if os.path.exists(config_path):
            logger.debug("Loading config from {}".format(config_path))
            try:
                with open(config_path, 'r') as config_file:
                    config_file_content = config_file.read()
                    file_config = yaml.safe_load(config_file_content)
                    configs.append((config_path, file_config))
            except IOError as e:
                logger.error("Could not read config file {}".format(config_path), e)
                pass
    return configs


file_configs = read_file_configs()

_config_dict = copy.deepcopy(default_config)
for current_config_path, current_file_config in reversed(file_configs):
    for key, value in current_file_config.items():
        if key not in default_config:
            warnings.warn("Config option '{}' is not known".format(key))
        else:
            _config_dict[key] = value

config = Map(_config_dict)


def write_default_config(file_name=None, overwrite=False):
    """
    Write the default configuration file.

    :param file_name: Output path, e.g. ~/.fanc.conf
    :param overwrite: Overwrite existing file if True. Default: False
    """
    if file_name is None:
        file_name = os.path.expanduser('~/.fanc.conf')
    else:
        file_name = os.path.expanduser(file_name)
    if not overwrite and os.path.exists(file_name):
        raise IOError("Config file exists in this location ({}). "
                      "Use overwrite=True or remove file.".format(file_name))

    with open(file_name, 'w') as output_config_file:
        output_config_file.write(_default_config_content)
