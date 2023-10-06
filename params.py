
__author__ = "Xiaohong Yuan"

PPM_MULTI = 'ppm_multi'
MASS_MIN = 'mass_min'
MASS_MAX = 'mass_max'
RT_MIN = 'rt_min'
RT_MAX = 'rt_max'

SCATTER_ALL = 'scatter_all'
USE_4_NUCLE = 'use4nucle'

DEFAULT_PPM_MULTI = 10 
DEFAULT_MASS_MIN = 600
DEFAULT_MASS_MAX = 8000
DEFAULT_RT_MIN = 0
DEFAULT_RT_MAX = 20
DEFAULT_SCATTER_ALL = False
DEFAULT_USE_4_NUCLE = False

class CmdParams:

    def __init__(self, kvals):
        self.ppm_multi = None
        self.mass_min = None
        self.mass_max = None
        self.rt_min = None
        self.rt_max = None
        self.scatter_all = None

        self.parse_argvs(kvals)
        
    def parse_argvs(self, kvals):
        ppm_multi = kvals.get(PPM_MULTI)
        if not ppm_multi:
            self.ppm_multi = DEFAULT_PPM_MULTI 
        else:
            self.ppm_multi = float(ppm_multi)
        if self.ppm_multi < 1 or self.ppm_multi > 150:
            self.ppm_multi = DEFAULT_PPM_MULTI 

        mass_min = kvals.get(MASS_MIN) 
        if not mass_min:
            self.mass_min = DEFAULT_MASS_MIN 
        else:
            self.mass_min = float(mass_min)
        if self.mass_min < 0:
            self.mass_min = DEFAULT_MASS_MIN 

        mass_max = kvals.get(MASS_MAX) 
        if not mass_max:
            self.mass_max = DEFAULT_MASS_MAX
        else:
            self.mass_max = float(mass_max)
        if self.mass_max < self.mass_min:
            self.mass_max = DEFAULT_MASS_MAX

        rt_min = kvals.get(RT_MIN)
        if not rt_min:
            self.rt_min = DEFAULT_RT_MIN
        else:
            self.rt_min = float(rt_min)
        if self.rt_min < DEFAULT_RT_MIN:
            self.rt_min = DEFAULT_RT_MIN 

        rt_max = kvals.get(RT_MAX)
        if not rt_max:
            self.rt_max = DEFAULT_RT_MAX
        else:
            self.rt_max = float(rt_max)
        if self.rt_max < self.rt_min:
            self.rt_max = DEFAULT_RT_MAX

        scatter_all = kvals.get(SCATTER_ALL)
        if not scatter_all:
            self.scatter_all = DEFAULT_SCATTER_ALL 
        else:
            self.scatter_all = not DEFAULT_SCATTER_ALL 

        use4nucle = kvals.get(USE_4_NUCLE)
        if not use4nucle:
            self.use4nucle = DEFAULT_USE_4_NUCLE
        else:
            self.use4nucle = not DEFAULT_USE_4_NUCLE

    def __str__(self):
        str = 'ppm_multi {}\nmass {}-{}\nRT {}-{}\nScatterAll {}\n4Nucle {}'.format(
                self.ppm_multi,
                self.mass_min, self.mass_max,
                self.rt_min, self.rt_max,
                self.scatter_all,
                self.use4nucle
                )
        return str
