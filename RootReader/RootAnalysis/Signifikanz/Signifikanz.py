import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import uproot
import os
from os import listdir
from os.path import isfile, join
import re
import math
"""Get into correct Dir"""
dirpath = os.getcwd()
path = "E:\\Master\\Skripte\\Auswertung\\Signifikanz\\"   
os.chdir( path )
runs = [
	#"dc_vb54",
	 #"dc_vb55",
	# "dc_vb56",
	 "dc_vb57",
	 "dc_vb58",
	"dc_vb59",
 	"dc_vb60",
	 "dc_vb61"
#	 "dc_vb62"
	] 

#Fit the order of runs
pdes=[0.56,0.6,0.64,0.65,0.655]

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


"""Pandas Dataframe-> Datenstruktur"""

is_amp = 1
trigger= 9 #Channel 9
drop= 8 #Unbenutzer Channel
runCount=len(runs)
channelCount=8
dcFractionArray=[]
dcFractionSumArray=[]
npe_cut=0.5
#tgate=100*10^-9
#tgate=0.1*math.pow(10, 0)
photonNrArray=[1] #just a scale factor -> useless
pe_val = ([1,2,3,4,5,6])
pe_val_sum = ([7,8,9,10,11]) #Number of PE Prob is calculated for/Thresholds

#For each run, store the dataframes with the probabilities to print them at the end
dataFrameArrayProbChannel=[]
dataFrameArrayProbSum=[]
dataFrameArrayCProbChannel=[]
dataFrameArrayCProbSum=[]



###########	
## CREATE ARRAY ##
########### 
for i in range(0, channelCount):
	dcFractionArray.append([])					#Channel
	for k in range(0,len(photonNrArray)):
		dcFractionArray[i].append([])			#PhotonNr
		for j in range(0,len(pe_val)):
			dcFractionArray[i][k].append([])	#PE VALS
			for l in range(0,runCount):      	#Runs
				dcFractionArray[i][k][j].append([])


for k in range(0,len(photonNrArray)):
		dcFractionSumArray.append([])			#PhotonNr
		for l in range(0,len(pe_val_sum)):      #PE VALS
			dcFractionSumArray[k].append([])
			for j in range(0,runCount):
				dcFractionSumArray[k][l].append([])	#Run
			


#print(dcFractionArray)
#Calculate
for runID in runs:
	for photonNr in photonNrArray: 
		runName="./rootfiles/"+runID+".root"
		tree = uproot.open(runName)["T"]
		indexInRuns=runs.index(runID)
		print(indexInRuns)
		if is_amp:
			# variable = "amp"
			# variable = "amp_inRange"
			variable = "chPE_amp"	#Alle Daten aus allen Channeln vermischt, siehe TBrowser b, 20K*9 Events im Histogram
			variable2 = "amp_array" #Summen
		else:
			variable = "Integral_inRange"
			# variable = "Integral"
			# variable = "chPE_charge"
			variable2 = "charge_array"

		"""Dataframe, bestehend aus Spalten, Reihen und Daten"""
		df_temp = tree.pandas.df([variable],flatten=False)

		"""ALLE Channel als Spalten, Nur Daten aus Dataframe"""
		df_ch = pd.DataFrame( df_temp[variable].to_numpy().tolist(), index=df_temp.index )
		del df_temp
		"""Summen, hier 20K Events"""
		df_sum = tree.pandas.df([variable2])


		"""Neue Spaltennamen zuweisen"""
		clmns = []
		for i in range(0,channelCount):
			clmns.append("{1}_ch{0}".format(i,variable))

		df_ch = df_ch.drop([drop],axis=1)
		df_ch.columns = clmns

		



		""" Erstellt Nummer von Subplots+ Figure, ax0=Array Axes=Darstellungsobjekt von Daten, fig0=Figure"""  
		"""  ax0 [Reihe, Spalte] ist der Subplot dann"""  

		"""           0         1         2         3         4          5         6         7    8
		entry
		0      1.124611  0.157786  0.719731 -0.124557 -0.029694   0.241216  0.285401  0.109142  0.0
		1      1.047989 -0.026520  0.053147  3.883717 -0.071948  -0.210920 -0.116301 -0.203468  0.0
		2     -0.047245 -0.229729 -0.142632 -0.069397 -0.160680  -0.030066 -0.073567  6.547208  0.0
		3      0.002334  0.134157  0.132391  0.279948  0.718192   1.283281  3.058850  0.041551  0.0
		4      0.024869  0.257027  0.412077  0.004149  1.605514   0.008689  1.199912  0.054224  0.0
		5      0.096984  0.550027  0.104423  0.974041  4.106917   4.685063  1.652895  0.295018  0.0
		6      1.174189  0.105802 -0.086695  0.073099  1.068895   1.235914 -0.398347  0.109142  0.0
		7     -0.299645 -0.031246  0.234943  0.100679 -0.278989   0.099116  0.986241  4.950365  0.0
		8     -0.259080  2.397813 -0.040081 -0.290036  0.929459  -0.150635  1.020429  0.320365  0.0
		9      0.146562  0.077447  1.255794 -0.032624  2.175935   0.107728  0.033269  0.185182  0.0
		10     1.958431  0.039641 -0.170601 -0.051011  0.464671   0.245522 -0.150488  0.020428  0.0
		11     0.813618  0.082173 -0.110002  0.082292 -0.215609   0.396234  1.020429  1.131459  0.0 """


		slice_names = []
		slice_names_sum = []
		for j in range(0,channelCount):
			for i in pe_val:
				slice_names.append("ch{0}_pe{1}".format(j,i))
				
		for i in pe_val_sum:
				slice_names_sum.append("sum_pe{0}".format(i))	

		#__ DATAFRAME OF PE SLICES ___
		sc = 0 # slice counter
		d_ch_pe_collection = {}
		d_sum_pe_collection = {}
		# indiv ch
		for i in range(0,channelCount):
			for k in pe_val:
				if k == 0:
					#ALLE WERTE die kleiner als 0.5 sind für jeden Channel
					d_ch_pe_collection[slice_names[sc]] = df_ch.iloc[:,i][(df_ch.iloc[:,i]<k+0.5) ]
				else:
					d_ch_pe_collection[slice_names[sc]] = df_ch.iloc[:,i][ (df_ch.iloc[:,i]>k-0.5) & (df_ch.iloc[:,i]<k+0.5) ]
				sc = sc+1
		sc = 0		
		# sum
		for k in pe_val_sum:
			if k == 0:
				d_sum_pe_collection[slice_names_sum[sc]] = df_sum.iloc[:,0][(df_sum.iloc[:,0]<k+0.5) ]
			else:
				d_sum_pe_collection[slice_names_sum[sc]] = df_sum.iloc[:,0][ (df_sum.iloc[:,0]>k-0.5) & (df_sum.iloc[:,0]<k+0.5) ]
			sc = sc+1

		df_ch_pe_slices = pd.DataFrame.from_dict(d_ch_pe_collection)
		df_sum_pe_slices = pd.DataFrame.from_dict(d_sum_pe_collection)

		print(df_sum_pe_slices)
		#ALL ROWS-> of i'th colummn
		#print(df_ch.iloc[:,i][(df_ch.iloc[:,i]<0+0.0) ])

		#__ DATAFRAME OF PE PROBABILITIES ____
		slice_ch_prob = np.zeros( (channelCount,len(pe_val)) )
		slice_sum_prob = np.zeros( (1,len(pe_val_sum)) )
		# indiv ch
		sc = 0
		for k in range(0,channelCount):
			for i in range(0,len(pe_val)):
				slice_ch_prob[k,i] = len( df_ch_pe_slices[slice_names[sc]].dropna() ) / len(df_ch)
				sc = sc+1

		# sum
		sc = 0		
		for i in range(0,len(pe_val_sum)):
			slice_sum_prob[0,i] = len( df_sum_pe_slices[slice_names_sum[sc]].dropna() ) / len(df_sum)
			sc = sc+1

		df_pe_val = pd.DataFrame(data=pe_val,columns=["Thr"])

		df_slice_ch_prob = pd.DataFrame(data=slice_ch_prob.T,columns=df_ch.columns)
		df_slice_ch_prob = pd.concat([df_pe_val,df_slice_ch_prob],axis=1)

		df_slice_sum_prob = pd.DataFrame(data=slice_sum_prob.T,columns=df_sum.columns)
		df_slice_sum_prob = pd.concat([df_pe_val,df_slice_sum_prob],axis=1)

		
	#	df_pe_slices_prob = pe_slice_prob(pd.DataFrame(),pe_val,1)	
	#	df_pe_slice_cprob = pe_cum_prob(df_pe_slices_prob,pe_val)


		dataFrameArrayProbChannel.append(df_slice_ch_prob)
		dataFrameArrayProbSum.append(df_slice_sum_prob)



		# print("PROB")	
		# print(df_slice_ch_prob)
		# print(df_slice_sum_prob)

		#__ DATAFRAME OF CUMMULATIVE PE PROBABILITIES ____
		slice_ch_cprob = np.zeros( (channelCount,len(pe_val)) )
		slice_sum_cprob = np.zeros( (1,len(pe_val_sum)) )
		# indiv ch
		sc = 0
		for k in range(0,channelCount):
			for i in range(0,len(pe_val)):
				slice_ch_cprob[k,i] = df_slice_ch_prob.iloc[i:len(pe_val),k+1].sum() # Spalte 0 ist Threshold
				sc = sc+1
		# sum
		sc = 0		
		for i in range(0,len(pe_val_sum)):
			slice_sum_cprob[0,i] = df_slice_sum_prob.iloc[i:len(pe_val_sum),1].sum()
			# print(df_slice_sum_prob.iloc[i:-1,1])
			sc = sc+1

		df_slice_ch_cprob = pd.DataFrame(data=slice_ch_cprob.T,columns=df_ch.columns)
		df_slice_ch_cprob = pd.concat([df_pe_val,df_slice_ch_cprob],axis=1)
		df_slice_sum_cprob = pd.DataFrame(data=slice_sum_cprob.T,columns=df_sum.columns)
		df_slice_sum_cprob = pd.concat([df_pe_val,df_slice_sum_cprob],axis=1)
		#print(df_slice_ch_cprob)

		dataFrameArrayCProbChannel.append(slice_ch_cprob) # [Channel,PeVal]
		dataFrameArrayCProbSum.append(slice_sum_cprob)

print(dataFrameArrayCProbSum)



'''Calculations'''

numberAllEvents=df_ch.shape[0]
#print(dataFrameArrayCProbChannel)
#print(dataFrameArrayCProbChannel[0])

"""Über Liste iterieren, dropna löscht fehlende Daten"""
for runID in runs:
	for photonNr in photonNrArray: 
		for dataframe in dataFrameArrayCProbChannel:
			for i in range(0,channelCount):
				for j in range(0,len(pe_val)):  # Für jeden Treshold Wert
					#numberGreaterThanCut=df_ch[clmns[i]][df_ch[clmns[i]]>npe_cut].dropna().values.size
					#dcr=math.log(numberGreaterThanCut/numberAllEvents)/tgate
				#	dcr=math.log(dataframe[i][j])/tgate
					#g=(photonNrArray[photonNrArray.index(photonNr)]*pdes[runs.index(runID)])/math.sqrt(dcr)
					dcFractionArray[i][photonNrArray.index(photonNr)][j][runs.index(runID)]=0


counter=0
for dataframe in dataFrameArrayCProbSum:
	for p in range(0,len(pe_val_sum)):  # Für jeden Treshold Wert
				#numberGreaterThanCut=df_ch[clmns[i]][df_ch[clmns[i]]>npe_cut].dropna().values.size
					#dcr=math.log(numberGreaterThanCut/numberAllEvents)/tgate
					#print(dataframe[0][j])
					v=dataframe[0][p]
					#dcr=math.log(v)/tgate
					dcr=v
					if dcr==0:
						g=1
					else:
						eff=pdes[counter]
						g=(eff)/math.sqrt(dcr)
					dcFractionSumArray[photonNrArray.index(photonNr)][p][counter]=g
	counter=counter+1






	'''	#numberGreaterThanCut=df_ch[clmns[i]][df_ch[clmns[i]]>npe_cut].dropna().values.size
					#dcr=math.log(numberGreaterThanCut/numberAllEvents)/tgate

					dcr=math.log(dataframe[j][i+1])/tgate
					g=(photonNrArray[photonNrArray.index(photonNr)]*pdes[runs.index(runID)])/math.sqrt(dcr)
					dcFractionArray[i][photonNrArray.index(photonNr)][runs.index(runID)]=g


numberGreaterThanCut=df_sum[:][df_sum[:]>npe_cut].dropna().values.size
dcr=math.log(numberGreaterThanCut/numberAllEvents)/tgate
g=(photonNrArray[photonNrArray.index(photonNr)]*pdes[runs.index(runID)])/math.sqrt(dcr)
dcFractionSumArray[photonNrArray.index(photonNr)][runs.index(runID)]=g'''

#print(dcFractionArray)







###########
## PLOTS ##
###########

ax0 = plt.subplot( )
#ax_chB = []
#for i in range(0,3):
	#for k in range(0,3):
	#	ax_chB.append(ax0[i,k])

voltage_array=[]
for j in range(0,runCount):
	str=runs[j]
	arr=re.findall(r'\d+', str)
	number=int(arr[0])
	voltage_array.append(number)


#dcFractionArray[Channel,PhotonNr,PE_VAL,run]


#[Channel, Photonnr, ]
#Create array for pe_val 0,1,2,3,...






lines2=[]
labels2=[]
for a in range(0,len(pe_val_sum)):
	lines2.append(0)
	labels2.append(0)
for j in range(0, len(photonNrArray)):
	for o in range(0,len(dcFractionSumArray[j])):
		ax = ax0
		x=ax.plot(voltage_array,dcFractionSumArray[j][o],marker='o',markersize=3)
		lines2[o]=x
		labels2[o]="N_pe={}".format(pe_val_sum[o])
		for xy in zip(voltage_array, dcFractionSumArray[j][o]):                                       # <--
			label=round(dataFrameArrayCProbSum[voltage_array.index(xy[0])][0][o]*100,2)
			ax.annotate('%s' % label, xy=xy, textcoords='data',fontsize=17) 
		ax.set_ylim(0,np.amax(dcFractionSumArray)+0.5)
		ax.set_ylabel("$\epsilon (V) / \sqrt{P_{DC}} $", fontsize=20)
		ax.set_xlabel("$U_{set} [mV]$", fontsize=20)
		ax.tick_params(axis='y', which='major', labelsize=17)
		ax.tick_params(axis='y', which='minor', labelsize=17)
		ax.tick_params(axis='x', which='major', labelsize=17)
		ax.tick_params(axis='x', which='minor', labelsize=17)
	#	ax.set_title("Sum Channel",fontsize=13)






#plt.figlegend( lines,labels=labels,loc = 'lower left', labelspacing=2.,prop={'size': 13})
plt.figlegend( lines2,labels=labels2,loc = 'lower left', labelspacing=0.1,prop={'size': 15})


##fig2, (ax2, ax3) = plt.subplots(nrows=2, ncols=1) # two axes on figure
#ax2.plot(x, z)
#ax3.plot(x, -z)



#print(dcFractionArray)
#plt.subplots_adjust(left=0.08, right=0.97, top=0.85, bottom=0.06, hspace=0.48, wspace=0.23)
#fig0.suptitle("Significance")
plt.show()



#ax0[1,0].plot([1, 2, 3, 4])



