#!/usr/bin/env python
""" 
Timing decorator.
"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015"
import time

def clockit(func):
    """
    A decorator to compute running time of a function.
    Usage::

        import clock
        @clock.clockit
        def func(self):
           pass 

    :param func: A function.
    """
    def wrapper(*arg, **kw):
        # Updated to use perf_counter (Python 3.3+) instead of deprecated clock()
        t0 = time.perf_counter()
        res = func(*arg, **kw)
        t1 = time.perf_counter()
        total_time = t1 - t0
        print("Time %s : %2.4f " % (func.__name__, total_time))
        return res
    
    return wrapper

