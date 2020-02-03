import signal
from fanc.config import config


class TrappedSignalException(Exception):
    pass


# trap signals requested by user and throw custom exception
def _trapped_signal_handler(signum, frame):
    try:
        message = "Trapped signal {} ({})".format(signum, signal.Signals(signum).name)
    except AttributeError:
        message = "Trapped signal {}".format(signum)

    raise TrappedSignalException(message)


def enable_signal_traps(signal_names=config.trap_signals):
    for signal_name in signal_names:
        s = getattr(signal, signal_name)
        signal.signal(s, _trapped_signal_handler)


def disable_signal_traps(signal_names=config.trap_signals):
    for signal_name in signal_names:
        s = getattr(signal, signal_name)
        signal.signal(s, signal.SIG_DFL)
