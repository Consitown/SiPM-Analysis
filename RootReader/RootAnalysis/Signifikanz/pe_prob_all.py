import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import uproot
import os

from os import listdir
from os.path import isfile, join
from PyPDF2 import PdfFileMerger
import glob
################
## INITIALIZE ##
################


#tree = uproot.open("../a_root_files/0_constBL_Calib_ampW50_intW25_scale0/{0}.root".format(runName))["T"]

runs = [
	#"dc_vb54",
	# "dc_vb55",
 	#"dc_vb56",
	 "dc_vb57",
	 "dc_vb58",
	 "dc_vb591",
 	"dc_vb60",
	 "dc_vb61",
	#"dc_vb62"
	] 


dirpath = os.getcwd()
print("current directory is : " + dirpath)
path = "E:\\Messungen\\Skripte\\DarkCounts\\"
os.chdir( path )

ampdirectory="./Plots/PE_PROB_AMP/"
chargedirectory="./Plots/PE_PROB_CH/"

files = glob.glob(ampdirectory+"*")
for f in files:
    os.remove(f)

files = glob.glob(chargedirectory+"*")
for f in files:
    os.remove(f)	

showHistogram=0
apo=0
api=0
fig0, ax0 = plt.subplots( nrows=3, ncols=3, figsize=(50,50) )

for m in [1]:
	for runID in runs:
		runName="./rootfiles/"+runID+".root"
		tree = uproot.open(runName)["T"]
		print(runID)
		is_amp = m
		if is_amp:
			variable = "chPE_amp"	
			variable2 = "amp_array"
			pe_def = "amplitude method"
			npe_idx = "amp"
			ch_npe_lim = 5
			sum_npe_lim = 10
		
		else:
			variable = "chPE_charge"
			variable2 = "charge_array"
			pe_def = "charge method"
			npe_idx = "charge"
			ch_npe_lim = 5
			sum_npe_lim = 15

		# indiv ch
		df_temp = tree.pandas.df([variable],flatten=False)
		df_ch = pd.DataFrame( df_temp[variable].values.tolist(), index=df_temp.index )
		del df_temp
		# sum
		df_sum = tree.pandas.df([variable2])

		clmns = []
		for i in range(0,8):
				clmns.append("{1}_ch{0}".format(i,variable))

		df_ch = df_ch.drop([8],axis=1)
		df_ch.columns = clmns

		#############################
		## PE SLICES of amp histos ##
		#############################

		pe_val = ([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25])
		pe_val_sum = ([0,1,2,3,4,5,6,7,8,9,10])
		slice_names = []
		slice_names_sum = []
		for j in range(0,8):
			for i in pe_val:
				slice_names.append("ch{0}_pe{1}".format(j,i))
				if j == 0:
					slice_names_sum.append("sum_pe{0}".format(i))	

		#__ DATAFRAME OF PE SLICES ___
		sc = 0 # slice counter
		d_ch_pe_collection = {}
		d_sum_pe_collection = {}
		# indiv ch
		for i in range(0,8):
			for k in pe_val:
				if k == 0:
					#ALLE WERTE die kleiner als 0.5 sind für jeden Channel
					d_ch_pe_collection[slice_names[sc]] = df_ch.iloc[:,i][(df_ch.iloc[:,i]<k+0.5) ]
				else:
					d_ch_pe_collection[slice_names[sc]] = df_ch.iloc[:,i][ (df_ch.iloc[:,i]>k-0.5) & (df_ch.iloc[:,i]<k+0.5) ]
				sc = sc+1
		sc = 0		
		# sum
		for k in pe_val:
			if k == 0:
				d_sum_pe_collection[slice_names_sum[sc]] = df_sum.iloc[:,0][(df_sum.iloc[:,0]<k+0.5) ]
			else:
				d_sum_pe_collection[slice_names_sum[sc]] = df_sum.iloc[:,0][ (df_sum.iloc[:,0]>k-0.5) & (df_sum.iloc[:,0]<k+0.5) ]
			sc = sc+1

		df_ch_pe_slices = pd.DataFrame.from_dict(d_ch_pe_collection)
		df_sum_pe_slices = pd.DataFrame.from_dict(d_sum_pe_collection)


		#__ DATAFRAME OF PE PROBABILITIES ____
		slice_ch_prob = np.zeros( (8,len(pe_val)) )
		slice_sum_prob = np.zeros( (1,len(pe_val)) )
		# indiv ch
		sc = 0
		for k in range(0,8):
			for i in range(0,len(pe_val)):
				slice_ch_prob[k,i] = len( df_ch_pe_slices[slice_names[sc]].dropna() ) / len(df_ch)
				sc = sc+1
		# sum
		sc = 0		
		for i in range(0,len(pe_val)):
			slice_sum_prob[0,i] = len( df_sum_pe_slices[slice_names_sum[sc]].dropna() ) / len(df_sum)
			sc = sc+1


		df_pe_val = pd.DataFrame(data=pe_val,columns=["Thr"])

		df_slice_ch_prob = pd.DataFrame(data=slice_ch_prob.T,columns=df_ch.columns)
		df_slice_ch_prob = pd.concat([df_pe_val,df_slice_ch_prob],axis=1)

		df_slice_sum_prob = pd.DataFrame(data=slice_sum_prob.T,columns=df_sum.columns)
		df_slice_sum_prob = pd.concat([df_pe_val,df_slice_sum_prob],axis=1)

		#__ DATAFRAME OF CUMMULATIVE PE PROBABILITIES ____
		slice_ch_cprob = np.zeros( (8,len(pe_val)) )
		slice_sum_cprob = np.zeros( (1,len(pe_val)) )
		# indiv ch
		sc = 0
		for k in range(0,8):
			for i in range(0,len(pe_val)):
				slice_ch_cprob[k,i] = df_slice_ch_prob.iloc[i:len(pe_val),k+1].sum()
				sc = sc+1
		# sum
		sc = 0		
		for i in range(0,len(pe_val)):
			slice_sum_cprob[0,i] = df_slice_sum_prob.iloc[i:len(pe_val),1].sum()
			# print(df_slice_sum_prob.iloc[i:-1,1])
			sc = sc+1

		df_slice_ch_cprob = pd.DataFrame(data=slice_ch_cprob.T,columns=df_ch.columns)
		df_slice_ch_cprob = pd.concat([df_pe_val,df_slice_ch_cprob],axis=1)
		df_slice_sum_cprob = pd.DataFrame(data=slice_sum_cprob.T,columns=df_sum.columns)
		df_slice_sum_cprob = pd.concat([df_pe_val,df_slice_sum_cprob],axis=1)

		###########
		## PLOTS ##
		###########

		nBins = 150
		x_min = 0
		if is_amp:
			sum_lim = 40
			ch_lim = 25
		else:
			sum_lim = 45
			ch_lim = 30


		#fig0.suptitle(r"$\bf{{Dark\ count\ probability\ and\ distribution\ of\ indiv.\ SiPM\ channels\ and\ SiPM\ array\ sum\ }} p.e. definition: {0}".format(pe_def), ha="left", x = 0.08, y = 0.98)
		if is_amp:
			savename="AMP"
		else:
			savename="CH"
		fig0.suptitle("PE PROB: " +savename)
		ax_ch = []
		axx_ch = []
		for i in range(0,3):
			for k in range(0,3):
				ax_ch.append(ax0[i,k])
				axx_ch.append(ax0[i,k])

		ch_prob_color = "black"
		ch_hist_color = "dodgerblue"
		sum_prob_color = "black"
		sum_hist_color = "red"
		y_labelsize = 7
		x_labelsize = 7
		grid_alpha = 0.15
		leg_fsize = 5

		# indiv ch
		for i in range(0,8):
			
			axx = axx_ch[i]
			ax = ax_ch[i]
			
			if showHistogram: 
				ax.hist(df_ch[clmns[i]][df_ch[clmns[i]]<ch_lim].dropna().values, bins=nBins, histtype="stepfilled",color=ch_hist_color,alpha=0.8, log=False,label="\nDark count dist.\nmean = {:1.2f} $N_{{pe}}$".format(df_ch[clmns[i]].mean()),zorder=2)

				ax.set_title("SiPM Channel {0}".format(i),fontsize=9)
			
				ax.legend(loc="lower left",fontsize=leg_fsize)
		#		ax.set_xticks((np.arange(x_min, ch_lim+1,1)))
				ax.set_xlim(x_min,ch_lim)
				ax.set_yscale("log")
				ax.tick_params(axis='y', labelsize=y_labelsize)
				ax.set_ylabel("Dark count distribution",fontsize = 8,color=ch_hist_color)


			axx.plot(pe_val,df_slice_ch_prob.iloc[:,i+1].values,linestyle="-",linewidth=1,ms=3,marker="s",label="\nDark count prob.\n$P(N_{{pe}}={:d})$ = {:1.1f} % for {:s}".format(ch_npe_lim,df_slice_ch_prob.iloc[ch_npe_lim,i+1]*100,runID),zorder=1)
			axx.plot(pe_val,df_slice_ch_cprob.iloc[:,i+1].values,linestyle=":",linewidth=1,ms=3,marker="o",label=("Dark count cumm. prob.\n"+r"$P^{{\ast}}(N_{{pe}}={:d})$ = {:1.1f} %".format(ch_npe_lim,df_slice_ch_cprob.iloc[ch_npe_lim,i+1]*100)),zorder=1)
			axx.legend(loc="upper right",fontsize=leg_fsize)
			axx.grid(True,"both","both",alpha=grid_alpha)
			axx.tick_params(axis='y',labelsize=y_labelsize)
			axx.tick_params(axis='x',labelsize=y_labelsize)
			axx.set_yscale("log")

			axx.set_xlabel("Threshold [$N_{pe}$]",fontsize = 8)
			axx.set_ylabel(r"Dark count probability $P^{{(\ast)}}(N_{pe})$",fontsize = 8,color=ch_prob_color,va="baseline")

		# sum
		axx = axx_ch[8]
		axx.plot(pe_val,df_slice_sum_prob.iloc[:,1].values,linestyle="-",linewidth=1,ms=3,marker="s",label="\nDark count prob.\n$P(N_{{pe}}={:d})$ = {:1.1f} %".format(sum_npe_lim,df_slice_sum_prob.iloc[sum_npe_lim,1]*100))
		axx.plot(pe_val,df_slice_sum_cprob.iloc[:,1].values,linestyle=":",linewidth=1,ms=3,marker="o",label=("Dark count cumm. prob.\n"+r"$P^{{\ast}}(N_{{pe}}={:d})$ = {:1.1f} %".format(sum_npe_lim,df_slice_sum_cprob.iloc[sum_npe_lim,1]*100)),zorder=1)
		axx.grid(True,"both","both",alpha=grid_alpha)
		axx.tick_params(axis='y', labelsize=y_labelsize)
		axx.tick_params(axis='x',labelsize=y_labelsize)
		axx.set_yscale("log")
		axx.legend(loc="upper right",fontsize=leg_fsize)

		if showHistogram: 
			ax = ax_ch[8]
			
			ax.hist(df_sum[variable2][df_sum[variable2]<sum_lim].dropna().values, bins=nBins, histtype="stepfilled",alpha=0.8, log=False,label="\nDark count dist.\nmean = {:1.2f} $N_{{pe}}$".format(df_sum[variable2].mean()))

			ax.set_title("SiPM Array Sum",fontsize=9)
			#ax.set_xticks((np.arange(x_min, sum_lim+1,2)))
			ax.set_xlim(x_min,sum_lim)
			# ax.set_yscale("log")
			ax.tick_params(axis='y', labelcolor=sum_hist_color,labelsize=y_labelsize)
			ax.legend(loc="lower left",fontsize=leg_fsize)
			ax.set_ylabel("Dark count distribution",fontsize = 8,color=sum_hist_color)

		axx.set_xlabel("Threshold [$N_{pe}$]",fontsize = 8)
		
		axx.set_ylabel(r"Dark count probability $P^{{(\ast)}}(N_{pe})$",fontsize = 8,va="baseline")
		

		
	if is_amp:
		savedir=ampdirectory
	else:
		savedir=chargedirectory

	#plt.savefig(savedir+"ALL.pdf")
	#plt.close()
	plt.show()

	"""onlyfiles = [f for f in listdir(savedir) if isfile(join(savedir, f))]
	pdf_merger = PdfFileMerger()
	for path in onlyfiles:
			pdf_merger.append(savedir+path)
	if is_amp:
		savename="AMP"
	else:
		savename="CH"
	with open(savedir+"Combined_"+savename+".pdf", 'wb') as fileobj:
			pdf_merger.write(fileobj)"""