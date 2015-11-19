"""
Provides custom per-object instance LRU (last recently used) cache

from: http://stackoverflow.com/a/18723434
"""

import time
import functools
import collections


def lru_cache(maxsize=255, timeout=None):
    """
    lru_cache(maxsize = 255, timeout = None) --> returns a decorator which returns an instance (a descriptor).

    This decorator factory will wrap a function / instance method and will
    supply a caching mechanism to the function. For every given input params
    it will store the result in a queue of maxsize size, and will return a cached
    ret_val if the same parameters are passed.

    - If an instance method is wrapped, each instance will have it's own cache
      and it's own timeout.
    - The wrapped function will have a cache_clear variable inserted into
      it and may be called to clear it's specific cache.
    - The wrapped function will maintain the original function's docstring and name (wraps)
    - The type of the wrapped function will no longer be that of a function
      but either an instance of _LRU_Cache_class or a functool.partial type.

    No error handling is done, in case an exception is raised - it will permeate up.

    :param maxsize: (int) The cache size limit, anything added above
                    that will delete the first values enterred (FIFO).
                    This size is per instance, thus 1000 instances with
                    maxsize of 255, will contain at max 255K elements.
    :param timeout: (int / float / None) Every n seconds the cache is deleted,
                    regardless of usage. If None - cache will never be refreshed.
    """

    class _LRU_Cache_class(object):
        def __init__(self, input_func, max_size, timeout):
            self._input_func = input_func
            self._max_size = max_size
            self._timeout = timeout

            # This will store the cache for this function, format - {caller1 : [OrderedDict1, last_refresh_time1], caller2 : [OrderedDict2, last_refresh_time2]}.
            # In case of an instance method - the caller is the instance, in case called from a regular function - the caller is None.
            self._caches_dict = {}

        def cache_clear(self, caller=None):
            # Remove the cache for the caller, only if exists:
            if caller in self._caches_dict:
                del self._caches_dict[caller]
                self._caches_dict[caller] = [collections.OrderedDict(), time.time()]

        def __get__(self, obj, objtype):
            """ Called for instance methods """
            return_func = functools.partial(self._cache_wrapper, obj)
            return_func.cache_clear = functools.partial(self.cache_clear, obj)
            # Return the wrapped function and wraps it to maintain the docstring and the name of the original function
            return functools.wraps(self._input_func)(return_func)

        def __call__(self, *args, **kwargs):
            """ Called for regular functions """
            return self._cache_wrapper(None, *args, **kwargs)
        # Set the cache_clear function in the __call__ operator
        __call__.cache_clear = cache_clear

        def _cache_wrapper(self, caller, *args, **kwargs):
            # Create a unique key including the types (in order to differentiate between 1 and '1')
            kwargs_key = "".join(map(
                lambda x: str(x) + str(type(kwargs[x])) + str(kwargs[x]),
                sorted(kwargs)
            ))
            key = "".join(map(
                lambda x: str(type(x)) + str(x),
                args
            )) + kwargs_key

            # Check if caller exists, if not create one
            if caller not in self._caches_dict:
                self._caches_dict[caller] = [collections.OrderedDict(), time.time()]
            else:
                # Validate in case the refresh time has passed
                if self._timeout != None:
                    if time.time() - self._caches_dict[caller][1] > self._timeout:
                        self.cache_clear(caller)

            # Check if the key exists, if so - return it
            cur_caller_cache_dict = self._caches_dict[caller][0]
            if key in cur_caller_cache_dict:
                return cur_caller_cache_dict[key]

            # Validate we didn't exceed the max_size:
            if len(cur_caller_cache_dict) >= self._max_size:
                # Delete the first item in the dict:
                cur_caller_cache_dict.popitem(False)

            # Call the function and store the data in the cache (call it with the caller in case it's an instance function - Ternary condition):
            if caller is not None:
                cur_caller_cache_dict[key] = self._input_func(caller, *args, **kwargs)
            else:
                cur_caller_cache_dict[key] = self._input_func(*args, **kwargs)

            return cur_caller_cache_dict[key]

    # Return the decorator wrapping the class (also wraps the instance to maintain the docstring and the name of the original function):
    return (lambda input_func: functools.wraps(input_func)(_LRU_Cache_class(input_func, maxsize, timeout)))
