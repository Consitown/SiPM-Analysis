import glob, os, sys
import numpy as np
# import ROOT as r
# from multiprocessing import Pool
import pandas as pd
import matplotlib
matplotlib.use('agg')
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import matplotlib.pyplot as plt
import uproot
import shutil
import math
from os import listdir
from os.path import isfile, join
#from IPython.display import display
#start from Visual Studio Code (F1 with RunExtension) or with F5 (Debugging Sessions)



#Function to calculate the dc probability for each slide
#Input: df-> dataframe with raw measurement data, looking like:
#			 Integral_ch0
#entry
#0          -0.075697
#1           0.237665
#2          -0.386032
#pe_val ->  Array with cut steps, looking like (0,1,2,3,4,5,6...)^T
#
def pe_slice_prob(df,pe_val, granularity):
	slice_names = []
	for i in pe_val:
		slice_names.append("N_pe{0}".format(i))

	sc = 0 # slice counter
	d_pe_collection = {}

	for k in pe_val:
		if k == 0:
			d_pe_collection[slice_names[sc]] = df.iloc[:,0][(df.iloc[:,0]<k+granularity/2) ] #All Values inside [-inf;0.5]
			#df.iloc[:,0] choose all values inside the first column (where actually the data is), if they are smaller than 0.5
		else:
			d_pe_collection[slice_names[sc]] = df.iloc[:,0][ (df.iloc[:,0]>=k-granularity/2) & (df.iloc[:,0]<k+granularity/2) ] #All Values in [k-0.5;k+0.5]
		sc = sc+1

	df_pe_slices = pd.DataFrame.from_dict(d_pe_collection) #create DF from slices
	slice_prob = np.zeros( (1,len(pe_val)) )

	#Calculate the probability from the entries
	sc = 0
	for i in range(0,len(pe_val)):
		slice_prob[0,i] = len( df_pe_slices[slice_names[sc]].dropna() ) / len(df) # Probability = events in range / all events
		sc = sc+1

	df_slice_sum_prob = pd.DataFrame(data=slice_prob.T,columns=df.columns) # new Dataframe with probabilities and column names

	return df_slice_sum_prob

#Function to calculate the cummulativ probabilities
#Input: Inidividual slice probabilities, pe_val array
#It just sums all the probabilities for all following slices starting from a specific slice. So cummulative Prob for slice 4 would be the prob for 4+5+6..30
def pe_cum_prob(df_pe_slices_prob,pe_val,):
	slice_cprob = np.zeros( (1,len(pe_val)) ) #1 Row, 30 Columns
	# sum
	for i in range(0,len(pe_val)):
		slice_cprob[0,i] = df_pe_slices_prob.iloc[i:len(pe_val),0].sum() #Sum over all slices starting from a a npe value-> first value should be near 1

	df_slice_cprob = pd.DataFrame(data=slice_cprob.T,columns=df_pe_slices_prob.columns)
	#		Integral_ch0
	#	0       0.015894
	#	1       0.015445


	return df_slice_cprob


################
## INITIALIZE ##
################

directory = os.path.join( os.getcwd())
print(directory)

if 'mnt' in directory:
	root_dir="/mnt/d/Programme/RootAnalysis/RootAnalysis/rootfiles"
	out_dir= "/mnt/d/Programme/RootAnalysis/RootAnalysis/DarkCounts/DCRPlots/"
else:
	root_dir="D:\\Programme\\RootAnalysis\\RootAnalysis\\rootfiles"
	out_dir= "D:\\Programme\\RootAnalysis\\RootAnalysis\\DarkCounts\\DCRPlots\\"

if not os.path.isdir(out_dir):
	os.system("mkdir -p %s"%out_dir)

run_list = [f for f in listdir(root_dir) if isfile(join(root_dir, f))]
run_list = [k for k in run_list if 'dc' in k]



#############
## ANALYZE ##
#############

# create N_pe threshold array
npe_max = 30

npe_granularity=0.1
npe_number=math.ceil(npe_max/npe_granularity)
channelXMax=15
sumXMax=30
channelXMin=-1.5
sumXMin=-5.5

threshold=0.05 #Percentage
#threshold=0.5 #0.5NPE 1 NPE ->
tresholdMode=0 #0 -> Find NPE, 1 -> Find Percentage



pe_val = np.zeros(npe_number+1)
for i in range(npe_number+1):
	pe_val[i] = i*npe_granularity

# store DC spectra
d_data = {}

# store prob results
d_prob = {}
d_cprob = {}
v_mean = []
triggerChannel=8
#__ loop over ALL RUNS ____

logEnabled=True
run_count = 0

for branch in range(0,3):
	if branch==0:
		ch_variable = "Integral"
		sum_variable = "chargeChannelSumWOM" #Keeps all WOM sums, for wavecatcher 16Ch with Channel 0-7 -> WOMID=3=D -> [3] keeps the results, proof with print(df_sum)
	if branch==1:
		ch_variable = "IntegralErrorP"
		sum_variable = "chargeChannelSumWOMErrorP"
	if branch==2:
		ch_variable = "IntegralErrorM"
		sum_variable = "chargeChannelSumWOMErrorM"



	for runName in run_list:
		runNameWithoutExtension=os.path.splitext(runName)[0]
		print("DOING: "+runNameWithoutExtension+" With Branch: "+ch_variable)
		path="{:s}{:s}\\".format(out_dir,runNameWithoutExtension,runNameWithoutExtension)

		if not os.path.isdir(path):
			#shutil.rmtree(path)
			os.system("mkdir %s"%path)

		tree = uproot.open(root_dir+"/{0}".format(runName))["T"]





		# get single channel charge branch
		df_temp = tree.pandas.df([ch_variable],flatten=False)
		df_ch = pd.DataFrame( df_temp[ch_variable].values.tolist(), index=df_temp.index )
		del df_temp

		clmns = []
		for i in range(0,8):
				clmns.append("{1}_ch{0}".format(i,ch_variable))

		df_ch = df_ch.drop([triggerChannel],axis=1) #drop first column -> Row Numbers
		df_ch.columns = clmns #New Column Names
		df_ch = df_ch.replace([np.inf, -np.inf], np.nan)

		# get sum branch
		df_sum = tree.pandas.df([sum_variable])
		df_sum = df_sum.drop(df_sum.columns[0:3],axis=1) #NEEDED cause df_sum= array with size 4, and 0-2 are holding bullshit
		df_sum.columns.values[0] = sum_variable #rename column back from sumXXX[3] to sumXXX

		df_sum = pd.DataFrame(df_sum.iloc[:,0],columns=[sum_variable])
		df_sum = df_sum.replace([np.inf, -np.inf], np.nan)

		# loop over datasets, 8 channels + 1 sum
		for i in range(0,9):

			df_charge = pd.DataFrame() #Create Empty Dataframe
			data_id = ""

			if i!=triggerChannel: #Channel
				df_charge = pd.DataFrame(df_ch.iloc[:,i],columns=[clmns[i]]) #Create Channel DF with 1 Column on ith place and column name from clmns
				data_id = "{:s}_ch{:d}".format(runNameWithoutExtension,i)
			else:
				df_charge = df_sum
				data_id = "{:s}_sum".format(runNameWithoutExtension)

			# mean of distribution
			mean_charge = df_charge.mean()

			# probability of PE slices
			df_pe_slices_prob = pe_slice_prob(df_charge,pe_val,npe_granularity)
			# cumulative probability of PE slices
			df_pe_slice_cprob = pe_cum_prob(df_pe_slices_prob,pe_val)


			#___ SAVE RESULTS ____
			d_data["{:s}".format(data_id)] = df_charge.iloc[:,0]
			d_prob["prob_{:s}".format(data_id)] = df_pe_slices_prob.iloc[:,0]
			d_cprob["cprob_{:s}".format(data_id)] = df_pe_slice_cprob.iloc[:,0]
			v_mean.append( mean_charge )

		run_count = run_count + 1


	df_data = pd.DataFrame.from_dict(d_data)
	df_pe_val = pd.DataFrame(data=pe_val,columns=["Thr"])
	df_prob = pd.DataFrame.from_dict(d_prob)
	df_prob = pd.concat(([df_pe_val,df_prob]),axis=1)
	df_cprob = pd.DataFrame.from_dict(d_cprob)
	df_cprob = pd.concat(([df_pe_val,df_cprob]),axis=1)

	df_mean = pd.DataFrame(np.array(v_mean))

	###########
	## PLOTS ##
	###########
	x_min=-6.5
	# create set of bin-defining histogram edges of variable size
	xbins = np.zeros(npe_number+2)
	xbins[0] = x_min
	xbins[1] = npe_granularity/2

	for i in range(npe_number+1):
		xbins[i+1] = i*npe_granularity+ npe_granularity/2




	ch_prob_color = "black"
	ch_hist_color = "dodgerblue"
	sum_prob_color = "black"
	sum_hist_color = "red"
	y_labelsize = 8
	x_labelsize = 8
	grid_alpha = 0.15
	leg_fsize = 9

	print("Calculations done!")

	file = open("{:s}Thresholds_{:s}.txt".format(out_dir,ch_variable),"a")
	# loop over runs
	for i in range(0,len(run_list)):
		try:
			runNameWithoutExtension=os.path.splitext(run_list[i])[0]
			print("Printing: {:s}".format(runNameWithoutExtension))
			#_____ SINGLE PLOTS ____
			# loop over datasets, 8 channels + 1 sum
			fig_title = "dark count pulse-charge spectrum"
			file.write("{:s}={{".format(runNameWithoutExtension))

			for k in range(0,9):
				data_id = ""
				single_fig_title = ""
				pdf_filename = ""
				if k!=triggerChannel:
					data_id = "{:s}_ch{:d}".format(runNameWithoutExtension,k)
					single_fig_title = str(fig_title+", ch{:d}".format(k))
					pdf_filename = "{:s}{:s}\\{:s}_ch{:d}_{:s}.pdf".format(out_dir,runNameWithoutExtension,runNameWithoutExtension,k,ch_variable)

					x_min = channelXMin
					x_max =channelXMax
					x_ticks = np.arange(0, x_max+1,1)
				else:
					data_id = "{:s}_sum".format(runNameWithoutExtension)
					single_fig_title = str(fig_title+", sum")
					pdf_filename = "{:s}{:s}\\{:s}_sum_{:s}.pdf".format(out_dir,runNameWithoutExtension,runNameWithoutExtension,ch_variable)
					x_min = sumXMin
					x_max = sumXMax
					x_ticks = np.arange(0, x_max+1,5)

				fig0, ax0 = plt.subplots( nrows=1, ncols=1, figsize=(8,5) )
				fig0.suptitle(single_fig_title,ha="left", x = 0.05, y = 0.96, fontsize=9)

				ax = ax0

				## weights normalize by number of histogram entries
				weights = np.ones_like(df_data[data_id].dropna().values)/float(len(df_data[data_id].dropna().values))


				## dc pulse-heigt, 10 times more bins -> quasi continuous
			#	nBins = int(((df_data[data_id].max() - df_data[data_id].min()) *10).round(0))
			#	ax.hist(df_data[data_id].dropna().values, weights=weights, bins=nBins, histtype="stepfilled",color=sum_hist_color,alpha=0.8, log=logEnabled,label="PCS raw-data".format(df_mean.iloc[(k+i*9),0]),zorder=1,density=False)

				## dc pulse-heigt, corrected binning
				binsPerNpe=int(1/npe_granularity)
				nBins = int(((df_data[data_id].max() - df_data[data_id].min()) *10).round(0))

				ax.hist(df_data[data_id].dropna().values, weights=weights, bins=xbins, histtype="stepfilled",color=ch_hist_color,alpha=0.8, log=logEnabled,zorder=2,density=False)

				if tresholdMode:
					matchingEntry=df_cprob[df_cprob.iloc[:,0] >= threshold]
				else:
					matchingEntry=df_cprob[df_cprob.iloc[:, (k+i*9)+1] <= threshold]





				rowOfMatching=matchingEntry.iloc[0,0]
				percentageOfMatching=matchingEntry.iloc[0,(k+i*9)+1]


				file.write("{:1.2f}/{:1.4f},".format(rowOfMatching,percentageOfMatching))

				#display(probMatchingPercent)
				#display(rowOfMatching)
				ax.plot(pe_val,df_cprob.iloc[:,(k+i*9)+1].values,color=ch_prob_color,linestyle=":",linewidth=1,ms=2,marker="o",zorder=2)


				legend_elements = [Line2D([0], [0], color='black', marker='o' , label=r"$P_{{DC}}(x)$"), Line2D([0], [0], color='black', lw=0, label=r"$P_{{thr}} \leq {:1.2f}\%$".format(percentageOfMatching*100)),Line2D([0], [0], color='black', lw=0, label=r"$N^{{thr}}_{{pe}} \geq {:1.2f}N_{{pe}}$".format(rowOfMatching))
				,Patch(facecolor=ch_hist_color, label=r"PCS rebinned with {:d} bins per $N_{{pe}}$".format(binsPerNpe,df_mean.iloc[(k+i*9),0]))]

				ax.legend(handles=legend_elements, loc="best",fontsize=10)


				# ax.set_title(single_fig_title,fontsize=10)
				ax.set_xlabel(r"$N_{{pe}}$",fontsize=10)
				ax.set_ylabel("$P_{DC}$",fontsize=10)
				ax.set_xticks(x_ticks)
				ax.set_xlim(x_min,x_max)
				# ax.tick_params(axis='y', labelcolor=ch_hist_color,labelsize=y_labelsize)
				ax.grid(True,"both","both",alpha=grid_alpha,color=ch_prob_color)

				plt.subplots_adjust(left=0.12, right=0.99, top=0.9, bottom=0.1, hspace=0.15, wspace=0.2)
				if os.path.isfile(pdf_filename):
  					os.remove(pdf_filename)
				plt.savefig(pdf_filename)
				plt.close()
		except KeyError as err:
			print("KeyError: {0}".format(err))
			print("Skipped: KeyError -> Missing Branches in Root File most likely!")
		file.write("}\n")

		#_____ COMBINED PLOT ____
		print("Printing Combined Plot...")
		fig1, ax1 = plt.subplots( nrows=3, ncols=3, figsize=(11,10) )
		fig1.suptitle(fig_title)
		ax_comb = []
		for j in range(0,3):
			for m in range(0,3):
				ax_comb.append(ax1[j,m])

		comb_pdf_filename = "{:s}{:s}\\Combined_{:s}_{:s}.pdf".format(out_dir,runNameWithoutExtension,runNameWithoutExtension,ch_variable)

		# loop over datasets, 8 channels + 1 sum
		for k in range(0,9):

			data_id = ""
			single_fig_title = ""
			if k!=8:
				data_id = "{:s}_ch{:d}".format(runNameWithoutExtension,k)
				single_fig_title = "ch{:d}".format(k)
				x_min = channelXMin
				x_max = channelXMax
				x_ticks = np.arange(0, x_max+1,2)
			else:
				data_id = "{:s}_sum".format(runNameWithoutExtension)
				single_fig_title = "sum"
				x_min = sumXMin
				x_max = sumXMax
				x_ticks = np.arange(0, x_max+1,5)

			ax = ax_comb[k]

			## weights normalize by number of histogram entries

			try:
				weights = np.ones_like(df_data[data_id].dropna().values)/float(len(df_data[data_id].dropna().values))
			except KeyError:
				continue

			## dc pulse-heigt, 1 bin per photoelctron
			#ax.hist(df_data[data_id].dropna().values, weights=weights, bins=xbins, histtype="step",color=ch_hist_color,alpha=0.8, log=logEnabled,label="1 bin per photoelectron\nmean: $N_{{pe}}$ = {:1.2f}\n ".format(df_mean.iloc[(k+i*9),0]),zorder=2,density=False)

			## dc pulse-heigt, 10 bin per photoelctron
			nBins = int(((df_data[data_id].max() - df_data[data_id].min()) *10).round(0))
			ax.hist(df_data[data_id].dropna().values, weights=weights, bins=nBins, histtype="stepfilled",color=sum_hist_color,alpha=0.8, log=logEnabled,label=r"PCS rebinned with {:d} bins per $N_{{pe}}$".format(binsPerNpe,df_mean.iloc[(k+i*9),0]),zorder=1,density=False)



			if tresholdMode:
				probMatchingPercent=df_cprob[df_cprob.iloc[:,0] >= threshold]
			else:
				probMatchingPercent=df_cprob[df_cprob.iloc[:, (k+i*9)+1] <= threshold]

			rowOfMatching=probMatchingPercent.iloc[0,0]
			percentageOfMatching=probMatchingPercent.iloc[0,(k+i*9)+1]




			ax.plot(pe_val,df_cprob.iloc[:,(k+i*9)+1].values,color=ch_prob_color,linestyle=":",linewidth=1,ms=1,marker="o",label=( r"$P_{{DC}}(\Lambda)$"+"\n"+r"$P_{{thr}} \leq {:1.2f}\%$".format(percentageOfMatching*100)+ "\n" +"$\Lambda_{{thr}} \geq {:1.2f}N_{{pe}}$".format(rowOfMatching)),zorder=2)

			ax.set_title(single_fig_title,fontsize=10)
			ax.set_xlabel(r"$N_{{pe}}$",fontsize=9)
			ax.set_ylabel(r"$P_{{DC}}$",fontsize = 9)

			legend_elements = [Line2D([0], [0], color='black', marker='o' , label=r"$P_{{DC}}(x)$"), Line2D([0], [0], color='black', lw=0, label=r"$P_{{thr}} \leq {:1.2f}\%$".format(percentageOfMatching*100)),Line2D([0], [0], color='black', lw=0, label=r"$N^{{thr}}_{{pe}} \geq {:1.2f}N_{{pe}}$".format(rowOfMatching))
			,Patch(facecolor=sum_hist_color, label=r"PCS rebinned with {:d} bins per $N_{{pe}}$".format(binsPerNpe,df_mean.iloc[(k+i*9),0]))]

			ax.legend(handles=legend_elements, loc="best",fontsize=6.5)


			ax.set_xticks(x_ticks)
			ax.set_xlim(x_min,x_max)
			# ax.tick_params(axis='y', labelcolor=ch_hist_color,labelsize=y_labelsize)
			ax.grid(True,"both","both",alpha=grid_alpha,color=ch_prob_color)

		plt.subplots_adjust(left=0.06, right=0.99, top=0.93, bottom=0.06, hspace=0.35, wspace=0.30)
		if os.path.isfile(comb_pdf_filename):
  			 os.remove(comb_pdf_filename)
		plt.savefig(comb_pdf_filename)
		plt.close()
	print("Done Branch:{:d}".format(branch))



