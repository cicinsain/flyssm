import sys
import glob
import matplotlib
matplotlib.use('Agg')
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Sans serif']})
matplotlib.rc('text', usetex=True)

import prettyplotlib as ppl
from prettyplotlib import plt
from prettyplotlib import mpl
import numpy as np

path = sys.argv[1]
# path = ""

files = glob.glob(path + "*.network")
# basefile =  "/home/amabdol/master/fly_nsga2/input/tllg58c13-v4.2/input_001.network"
# basedata = np.loadtxt(basefile).flatten()


basedata = [0.023445221,-0.020599646,0.009530062,0.017762582 
			,  0.029434593,0.003952001,0.016962347,0.022714978 
			,  0.004830379,0.007907694,0.030765266,-0.003459609 
			, -0.014494242,-0.028390379,-0.029581705,0.024965998 
			,0.03882552,0.00620087,0.01676870,-0.03275255 
			,0.07235224,0.01497864,-0.15010887,-0.00683359 
			,0.02629917,0.02338136,-0.07710051,-0.00214066 
			, -0.01333367,0.01843562,-0.07863500,-0.01182617]

fig = plt.figure()
ax = fig.add_subplot()




plt.plot(range(len(basedata)), basedata, color='red', linewidth=2)

for file in files:
	filedata = np.loadtxt(file).flatten()
	print file
	plt.plot(filedata, '-', color='blue', linewidth=.5)

fig.savefig(path + "perturb.png")