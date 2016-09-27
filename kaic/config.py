from kaic.tools.general import Map
import yaml
import os
import copy
import warnings
import logging
logger = logging.getLogger(__name__)

default_config_content = """\
colormap_hic: germany
"""
default_config = yaml.load(default_config_content)

_config_file_name = 'kaic.conf'
_environment_variable_name = 'KAIC_CONF'

# check common places for config file
# modified from: http://stackoverflow.com/questions/7567642/where-to-put-a-configuration-file-in-python
locations = (
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
file_config = None
for location in locations:
    if location is None:
        continue

    config_path = os.path.expanduser(location)
    if file_config is None and os.path.exists(config_path):
        logging.info("Loading config from {}".format(config_path))
        try:
            with open(config_path, 'r') as config_file:
                config_file_content = config_file.read()
                file_config = yaml.safe_load(config_file_content)
        except IOError, e:
            logging.error("Could not read config file {}".format(config_path), e)
            pass

config_dict = copy.deepcopy(default_config)
if file_config is not None:
    for key, value in file_config.iteritems():
        if key not in default_config:
            warnings.warn("Config option '{}' is not known".format(key))
        else:
            config_dict[key] = value

config = Map(config_dict)


def write_default_config(file_name, overwrite=False):
    file_name = os.path.expanduser(file_name)
    if not overwrite and os.path.exists(file_name):
        raise IOError("Config file exists in this location ({}). Use overwrite=True or remove file.".format(file_name))

    with open(file_name, 'w') as output_config_file:
        output_config_file.write(default_config_content)
