"""Decorator function for tracking memory usage
    Typical usage example:

    from utils import track

    @track
    def list_create(n):
        print("inside list create")
        return [1] * n

    Produces this output:
    >>> inside list create
    >>> list_create: memory before: 45,928,448, after: 46,211,072, consumed: 282,624; exec time: 00:00:00

    https://stackoverflow.com/questions/938733/total-memory-used-by-python-process

"""
import os
import time

import psutil


def elapsed_since(start):
    return time.strftime("%H:%M:%S", time.gmtime(time.time() - start))


def get_process_memory():
    process = psutil.Process(os.getpid())
    return process.memory_info().rss


def track(func):
    """

    :rtype: object
    """

    def wrapper(*args, **kwargs):
        mem_before = get_process_memory()
        start = time.time()
        result = func(*args, **kwargs)
        elapsed_time = elapsed_since(start)
        mem_after = get_process_memory()
        print("{}: memory before: {:,}, after: {:,}, consumed: {:,}; exec time: {}".format(
            func.__name__,
            mem_before, mem_after, mem_after - mem_before,
            elapsed_time))
        return result

    return wrapper
