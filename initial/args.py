
import sys
from os.path import expanduser

homefolder = expanduser('~')
sys.path.insert(0, homefolder+'/prog/four')

import numpy as np
import matplotlib.pyplot as plt

from FOUR.Transforms import *
