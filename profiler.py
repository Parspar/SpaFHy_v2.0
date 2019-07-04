# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 20:29:23 2019

@author: L1656
"""

import cProfile, pstats, io
from pstats import SortKey

def profile(fnc):

    def inner(*args, **kwargs):

        pr = cProfile.Profile()
        pr.enable()
        retval = fnc(*args, **kwargs)
        pr.disable()
        s = io.StringIO()
        sortby = SortKey.CUMULATIVE
        ps = pstats.Stats(pr, stream=s).strip_dirs().sort_stats(sortby)
        ps.print_stats(0.05)
        print(s.getvalue())
        return retval

    return inner