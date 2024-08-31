
import sys
from os.path import expanduser

homefolder = expanduser('~')
sys.path.insert(0, homefolder+'/prog/four')

import numpy as np
import matplotlib.pyplot as plt

from FOUR.Transforms import *

# Set global number print setting.
np.set_printoptions(precision=3, suppress=True, linewidth=np.inf)
