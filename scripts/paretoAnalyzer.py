import matplotlib
matplotlib.use('Agg')
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Sans serif']})
matplotlib.rc('text', usetex=True)


import prettyplotlib as ppl
from prettyplotlib import plt
from prettyplotlib import mpl

# import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import axes3d, Axes3D

import sys
import numpy as np

popSize = 0
nGen    = 0
gndp    = 0
tndp    = 0

def pareto_plotter_2d(n, data, last_pop_obj_rms_mean):
	fig, ax = plt.subplots(1)
	cmhot = plt.get_cmap("jet")

	for i in xrange(n):
		for j in range(i):
			ppl.scatter(ax, data[i], data[j], c = last_pop_obj_rms_mean, cmap = cmhot, edgecolors='none')
			ax.set_xlim(0,100)
			ax.set_xlabel(r'$f_%s$ RMS' % i)
			ax.set_ylim(0,100)
			ax.set_ylabel(r'$f_%s$ RMS' % j)
			fig.savefig(filename + str(i) + str(j) + ".png", bbox_inches='tight')
			ax.clear()
			

def pareto_plotter_3d(n, data, last_pop_obj_rms_mean):
	fig   = plt.figure() 
	ax    = Axes3D(fig)
	cmhot = plt.get_cmap("jet")
	for i in xrange(n):
		for j in range(i):
			for k in range(j):
				ax.scatter( data[i], data[j], data[k], c=last_pop_obj_rms_mean, cmap=cmhot , edgecolors='none')
				ax.set_xlim(0,)
				ax.set_xlabel(r'$f_%s$ RMS' % i)
				ax.set_ylim(0,)
				ax.set_ylabel(r'$f_%s$ RMS' % j)
				ax.set_zlim(0,)
				ax.set_zlabel(r'$f_%s$ RMS' % k)
				fig.savefig(filename + str(i) + str(j) + str(k) + ".png", bbox_inches='tight')
				ax.clear()


def pareto_plotter_4d(n, data, last_pop_obj_rms_mean=None):
	fig = plt.figure() 
	ax = Axes3D(fig)
	cmhot = plt.get_cmap("jet")
	ax.scatter( data[0], data[1], data[2], c = data[3], s= last_pop_obj_rms_mean.max() - last_pop_obj_rms_mean,  cmap=cmhot, edgecolors='none', label=r'$f_4$ RMS')
	# last_pop_obj_rms_mean.max() - last_pop_obj_rms_mean	--> big blues are better (easier to see...)
	# last_pop_obj_rms_mean - last_pop_obj_rms_mean.min()	--> small blues are better
	ax.set_xlim(0,)
	ax.set_xlabel(r'$f_%s$ RMS' % 0)
	ax.set_ylim(0,)
	ax.set_ylabel(r'$f_%s$ RMS' % 1)
	ax.set_zlim(0,)
	ax.set_zlabel(r'$f_%s$ RMS' % 2)
	# ax.legend(loc='upper right', fontsize = 'small', frameon = False)
	fig.savefig(filename + "4d" + ".png", bbox_inches='tight')
	ax.clear()

def compute_stats(n, data):
	deap_generation_objectives_stats = []
	deap_generation_total_stats      = []

	dd = np.array_split(data, nGen)
	dd = np.array(dd)

	objective_rms      = np.sqrt(dd/gndp)					# RMS of all objective
	# print objective_rms.shape
	mean_objective_rms = objective_rms.mean(axis=1)		# mean RMS of each population for each objective
	total_rms          = np.sqrt(dd.sum(axis=2)/tndp)		# 
	
	max_generation_objective_rms  = objective_rms.max(axis=1)
	min_generation_objective_rms  = objective_rms.min(axis=1)
	mean_generation_objective_rms = objective_rms.mean(axis=1)
	std_generation_objective_rms  = objective_rms.std(axis=1)
	
	mean_of_mean_generation_objecitve_rms = mean_generation_objective_rms.mean(axis=1)
	# print mean_of_mean_generation_objecitve_rms
	
	max_generation_total_rms  = total_rms.max(axis=1)
	min_generation_total_rms  = total_rms.min(axis=1)
	mean_generation_total_rms = total_rms.mean(axis=1)
	std_generation_total_rms  = total_rms.std(axis=1)

	# Can be plotted using flyPlotter Stat plotter
	deap_generation_objectives_stats.append(np.concatenate([std_generation_objective_rms, \
												max_generation_objective_rms, \
												mean_generation_objective_rms, \
												min_generation_objective_rms]))

	# Should write another script for import
	deap_generation_total_stats.append([std_generation_total_rms, \
											max_generation_total_rms, \
											mean_generation_total_rms, \
											min_generation_total_rms, \
											mean_of_mean_generation_objecitve_rms])

	return deap_generation_objectives_stats, deap_generation_total_stats, objective_rms[nGen-1], objective_rms[nGen-1].mean(axis=1)


def plot_total_stats(data):
	l = nGen
	fig = plt.figure()
	ax = fig.add_subplot(111)
	# plt.gca().set_ylim(bottom=0)
	# ax.set_ylim(top=1)
	# ax.set_title('\normalsize Single objective convergence for %s genes (based on MOOP results)' % ngenes)
	ax.set_xlabel('Generations')
	ax.set_ylabel('RMS')
	ax.plot(data[0][1], label='max', linewidth = 2)#, color = 'blue')
	ax.plot(data[0][2], label='mean', linewidth = 2)#, color = 'green')
	ax.plot(data[0][3], label='min', linewidth = 2)#, color = 'red')
	ax.plot(data[0][4], label='mean rms \nof objectives', linewidth = 2)#, color = 'cyan')
	ax.errorbar(range(l), data[0][2], yerr = data[0][0], fmt='--', linewidth = 0, elinewidth=.3)#color='green', label = 'std' )
	ax.legend(loc='upper right', fontsize = 'small', frameon = False)
	fig.savefig(filename + "_Single_objective_convergence_for_%s_genes.png" % ngenes, bbox_inches='tight')
	# fig.show()
	ax.clear()

def plot_statistics(data):
	l = nGen
	st = []; mx = []; av = []; mn = []

	for d in data:
		spl = np.array_split(d, 4)
		st.append(spl[0])
		mx.append(spl[1])
		av.append(spl[2])
		mn.append(spl[3])

	st = np.array(st)
	mx  = np.array(mx)
	av = np.array(av)
	mn  = np.array(mn)

	fig = plt.figure()
	for i in xrange(ngenes):
		ax = fig.add_subplot(111)
		# plt.gca().set_ylim(bottom=0)
		# ax.set_ylim(0, 250)
		# ax.set_title('\huge $f_%s$ \normalsize convergence for %s optimization' % ((i+1), ngenes))
		# ax.set_title('f%s' %  i)
		ax.set_xlabel('Generations')
		ax.set_ylabel('RMS')
		ax.plot(range(l), mx[0].T[i], label='max', linewidth = 2, )# label='max' , color = 'blue')
		ax.plot(range(l), av[0].T[i], label='mean', linewidth = 2)#, label='mean' , color = 'green')
		ax.plot(range(l), mn[0].T[i], label='max', linewidth = 2)#, label='min' , color = 'red')
		ax.errorbar(range(l), av[0].T[i], yerr = st[0].T[i], fmt='--', linewidth = 0, elinewidth=.3)#color = 'green', label = 'st', )
		ax.legend(loc='upper right', fontsize = 'small', frameon = False)
		fig.savefig(filename + "_f_%s_convergence_for_%s_optimization.png" % (i, ngenes), bbox_inches='tight')
		# fig.show()
		ax.clear()







def retrieve_best_result(data, fname):
	"""
		gStats: generation_stat
		tStats: total_stat
	"""	
	bins = np.arange(0, 100, 10)

	dd = np.array_split(data, nGen)
	dd = np.array(dd)

	objective_rms      = np.sqrt(dd/gndp)					# RMS of all objective
	total_rms          = np.sqrt(dd.sum(axis=2)/tndp)		# 

	objective_rms = objective_rms[nGen-1]
	total_rms = total_rms[nGen-1]


	fname = fname.replace("all_pop_", "")
	f = open("sorted_pop_" + fname, 'w')

	f.write("\n# Sorted solutions based on single_objective_rms\n")
	total_rms_sorted_index = total_rms.argsort()
	for i in total_rms_sorted_index:
		f.write(str(i) +" ")
	f.write("\n")

	for i in total_rms_sorted_index:
		f.write(str(total_rms[i]) +"\n")

	fig, ax = plt.subplots(1)
	ax.set_ylabel("Frequencies")
	ax.set_xlabel("RMS")
	ppl.hist(ax, total_rms)
	fig.savefig(filename + "_total_rms_histogram.png")



	f.write("\n# Sorted solutions based on mean rms of objectives\n")
	objective_rms_mean_sorted_index = objective_rms.mean(axis=1).argsort()
	for i in objective_rms_mean_sorted_index:
		f.write(str(i) +" ")
	f.write("\n")

	for i in objective_rms_mean_sorted_index:
		f.write(str(objective_rms[i]) +"\n")

	fig, ax = plt.subplots(1)
	ax.set_ylabel("Frequencies")
	ax.set_xlabel("RMS")
	ppl.hist(ax, objective_rms.mean(axis=1))
	fig.savefig(filename + "_objectives_mean_rms_histogram.png")

	# fig.clear()
	# plt.hist(objective_rms, stacked=True, normed= True)
	# fig.savefig(filename + "_objective_rms_stacked_histogram.png")


	f.write("\n# Sorted solutions based on lexsort of objectives rms\n")
	objective_rms_lexsort_index = np.lexsort(objective_rms.T)
	for i in objective_rms_lexsort_index:
		f.write(str(i) +" ")
	f.write("\n")

	for i in objective_rms_lexsort_index:
		f.write(str(objective_rms[i]) +"\n")






def write_stats(data):
	pass


if __name__ == '__main__':
	
	ngenes   = int(sys.argv[1])
	gndp     = (58*8 + 30)
	tndp     = ngenes*(58*8 + 30)
	popSize  = int(sys.argv[2])
	nGen     = int(sys.argv[3])
	filename = sys.argv[4]		
	# directory = sys.argv[3] 	# 0 or 1

	data = np.genfromtxt(filename, delimiter='\t')
	odata = data[:,:ngenes].T
	generation_stat, total_stat, pdata, last_pop_obj_rms_mean = compute_stats(ngenes, odata.T)
	plot_statistics(generation_stat)
	plot_total_stats(total_stat)
	retrieve_best_result(odata.T, filename)

	pareto_plotter_2d(ngenes, pdata.T, last_pop_obj_rms_mean)
	pareto_plotter_3d(ngenes, pdata.T, last_pop_obj_rms_mean)
	pareto_plotter_4d(ngenes, pdata.T, last_pop_obj_rms_mean)
