from kaic.tools.general import Map
import yaml
import os
import copy
import warnings
import logging
logger = logging.getLogger(__name__)

_default_config_content = """\
#
# PLOTTING
#
colormap_hic: germany

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
"""
default_config = yaml.load(_default_config_content)

_config_file_name = 'kaic.conf'
_environment_variable_name = 'KAIC_CONF'

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
    os.path.join(os.path.expanduser("~"), '.config', 'kaic', _config_file_name),
    # /etc/kaic (MySQL style)
    os.path.join("/etc/kaic/", _config_file_name)
)


# read the first config file that is found, checking in the order above
def read_file_configs(config_file_locations):
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

file_configs = read_file_configs(_config_file_locations)

_config_dict = copy.deepcopy(default_config)
for current_config_path, current_file_config in reversed(file_configs):
    for key, value in current_file_config.items():
        if key not in default_config:
            warnings.warn("Config option '{}' is not known".format(key))
        else:
            _config_dict[key] = value

config = Map(_config_dict)


def write_default_config(file_name, overwrite=False):
    file_name = os.path.expanduser(file_name)
    if not overwrite and os.path.exists(file_name):
        raise IOError("Config file exists in this location ({}). Use overwrite=True or remove file.".format(file_name))

    with open(file_name, 'w') as output_config_file:
        output_config_file.write(_default_config_content)
