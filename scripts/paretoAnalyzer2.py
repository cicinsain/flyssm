import matplotlib
matplotlib.use('Agg')
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Sans serif']})
matplotlib.rc('text', usetex=True)


import prettyplotlib as ppl
from prettyplotlib import plt
from prettyplotlib import mpl


from mpl_toolkits.mplot3d import axes3d, Axes3D

import sys
import numpy as np

popSize = 0
nGen    = 0
gndp    = 0		# gene data point
tndp    = 0		# total datapoint

def pareto_plotter_2d(ngenes, objective_rms, mean_of_objective_rms = None):
	fig, ax = plt.subplots(1)
	cmhot = plt.get_cmap("jet")

	# print objective_rms[1]

	for i in xrange(ngenes):
		for j in range(i):
			ppl.scatter(ax, objective_rms[i], objective_rms[j], c = mean_of_objective_rms, cmap = cmhot, edgecolors='none')
			ax.set_xlim(0,100)
			ax.set_xlabel(r'$f_%s$ RMS' % i)
			ax.set_ylim(0,100)
			ax.set_ylabel(r'$f_%s$ RMS' % j)
			fig.savefig(filename + str(i) + str(j) + ".png", bbox_inches='tight')
			ax.clear()
			

def pareto_plotter_3d(ngenes, objective_rms, mean_of_objective_rms = None):
	fig   = plt.figure() 
	ax    = Axes3D(fig)
	cmhot = plt.get_cmap("jet")
	for i in xrange(ngenes):
		for j in range(i):
			for k in range(j):
				ax.scatter( objective_rms[i], objective_rms[j], objective_rms[k], c=mean_of_objective_rms, cmap=cmhot , edgecolors='none')
				ax.set_xlim(0,)
				ax.set_xlabel(r'$f_%s$ RMS' % i)
				ax.set_ylim(0,)
				ax.set_ylabel(r'$f_%s$ RMS' % j)
				ax.set_zlim(0,)
				ax.set_zlabel(r'$f_%s$ RMS' % k)
				fig.savefig(filename + str(i) + str(j) + str(k) + ".png", bbox_inches='tight')
				ax.clear()


def pareto_plotter_4d(ngenes, objective_rms, mean_of_objective_rms=None):
	fig = plt.figure() 
	ax = Axes3D(fig)
	cmhot = plt.get_cmap("jet")
	ax.scatter( objective_rms[0], objective_rms[1], objective_rms[2], c = objective_rms[3], s= mean_of_objective_rms.max() - mean_of_objective_rms,  cmap=cmhot, edgecolors='none', label=r'$f_4$ RMS')
	# mean_of_objective_rms.max() - mean_of_objective_rms	--> big blues are better (easier to see...)
	# mean_of_objective_rms - mean_of_objective_rms.min()	--> small blues are better
	ax.set_xlim(0,)
	ax.set_xlabel(r'$f_%s$ RMS' % 0)
	ax.set_ylim(0,)
	ax.set_ylabel(r'$f_%s$ RMS' % 1)
	ax.set_zlim(0,)
	ax.set_zlabel(r'$f_%s$ RMS' % 2)
	# ax.legend(loc='upper right', fontsize = 'small', frameon = False)
	fig.savefig(filename + "4d" + ".png", bbox_inches='tight')
	ax.clear()


def retrieve_best_result(objective_rms, total_rms, mean_of_objective_rms, fname):
	"""
	"""
	f = open("sorted_pop_" + fname, 'w')

	f.write("\n# Sorted solutions based on single_objective_rms\n")
	total_rms_sorted_index = total_rms.argsort()
	for i in total_rms_sorted_index:
		f.write(str(i) +" ")
	f.write("\n")

	for i in total_rms_sorted_index:
		f.write(str(total_rms[i]) + "\n")

	fig, ax = plt.subplots(1)
	ax.set_ylabel("Frequencies")
	ax.set_xlabel("RMS")
	ppl.hist(ax, total_rms)
	fig.savefig(filename + "_total_rms_histogram.png")


	#------------------------------ Sorted 

	f.write("\n# Sorted solutions based on mean rms of objectives\n")
	# print mean_of_objective_rms
	objective_rms_mean_sorted_index = mean_of_objective_rms.argsort()
	for i in objective_rms_mean_sorted_index:
		f.write(str(i) +" ")
	f.write("\n")

	for i in objective_rms_mean_sorted_index:
		f.write(str(mean_of_objective_rms[i]) + "\t" +  str(objective_rms[i]) +"\n")

	fig, ax = plt.subplots(1)
	ax.set_ylabel("Frequencies")
	ax.set_xlabel("RMS")
	ppl.hist(ax, objective_rms.mean(axis=1))
	fig.savefig(filename + "_objectives_mean_rms_histogram.png")

	# fig.clear()
	# plt.hist(objective_rms, stacked=True, normed= True)
	# fig.savefig(filename + "_objective_rms_stacked_histogram.png")

	#------------------------- Lexsort

	f.write("\n# Sorted solutions based on lexsort of objectives rms\n")
	objective_rms_lexsort_index = np.lexsort(objective_rms.T)
	for i in objective_rms_lexsort_index:
		f.write(str(i) +" ")
	f.write("\n")

	for i in objective_rms_lexsort_index:
		f.write(str(objective_rms[i]) +"\n")



if __name__ == '__main__':

	ngenes   = int(sys.argv[1])
	gndp     = (58*8 + 30)
	tndp     = ngenes*(58*8 + 30)
	popSize  = int(sys.argv[2])
	nGen     = int(sys.argv[3])
	filename = sys.argv[4]		

	data = np.genfromtxt(filename, delimiter = '\t')
	objective_values = data[:, :ngenes]

	objective_rms = np.sqrt( objective_values / gndp)
	total_rms = np.sqrt(objective_values.sum(axis=1)/tndp)
	mean_of_objective_rms = objective_rms.mean(axis=1)

	# objective_rms = objective_rms.T
	# print objective_rms
	# print mean_of_objective_rms
	# print total_rms

	retrieve_best_result(objective_rms, total_rms, mean_of_objective_rms, filename)

	pareto_plotter_2d(ngenes, objective_rms.T, mean_of_objective_rms)
	pareto_plotter_3d(ngenes, objective_rms.T, mean_of_objective_rms)
	pareto_plotter_4d(ngenes, objective_rms.T, mean_of_objective_rms)