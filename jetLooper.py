import argparse
from plothelper import *
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(usage=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s", "--sampleName", help="sample name", default="higgs_test")
parser.add_argument("-d", "--darkSample", help="sample is dark", action="store_true")

ops = parser.parse_args()

sample = ops.sampleName
darkSample = ops.darkSample

file  = ROOT.TFile.Open(f"outputs/jets/{sample}_jets.root")
jet_tree = file.Get("Jets")
particle_tree = file.Get("Particles")

#c = ROOT.TCanvas("c","",800,800)
plt = PlotHelper("plots/jets")

def ttree_2_df(tree):
	# Create empty lists to hold the branch data
	data = {branch.GetName(): [] for branch in tree.GetListOfBranches()}

	# Loop over the entries in the TTree and fill the lists
	for entry in range(tree.GetEntries()):
		tree.GetEntry(entry)
		for branch in tree.GetListOfBranches():
			data[branch.GetName()].append(getattr(tree, branch.GetName()))

	# Convert the lists to a Pandas DataFrame
	df = pd.DataFrame(data)


	return df

	


jets = ttree_2_df(jet_tree)
particles = ttree_2_df(particle_tree)

# Print the first few rows of the DataFrame
#print(jets.iloc[0])
#print(len(jets))

#njetsPlot = 10
#njetsPlotted = 0
for i in range(len(jets)):
	particles_in_jet = particles[(particles['jet_index']==i) & (particles['status']==1)]
	n_final = len(particles_in_jet)

	

	# overlap removal, if there is a lepton skip jet
	lepton_jet = False
	for j,p in particles_in_jet.iterrows():
		if abs(p["pid"]) in [11,13,15]:
			lepton_jet = True
	
	if lepton_jet:
		continue
	

	if darkSample:
		
		if jets.iloc[i]['isdark']==1:
			
			particles_displaced = particles_in_jet[(np.sqrt(particles_in_jet['prod_x']**2 + particles_in_jet['prod_y']**2 )<300.) & (particles_in_jet['prod_z']<300.)]
			particles_displaced = particles_displaced[(np.sqrt(particles_displaced['prod_x']**2 + particles_displaced['prod_y']**2 )>1.)]
			particles_prompt = particles_in_jet[(np.sqrt(particles_in_jet['prod_x']**2 + particles_in_jet['prod_y']**2 )<=1.) & (particles_in_jet['prod_z']<=1.)]
				
			n_displaced = len(particles_displaced)
			n_prompt = len(particles_prompt)

			plt.plot1D(f"{sample}_jet_displaced_tracks",";displaced tracks;jets",n_displaced,20,0,20)
			plt.plot1D(f"{sample}_jet_prompt_tracks",";prompt tracks;jets",n_prompt,50,0,50)

		
			plt.plot1D(f"{sample}_jet_displacement",";displacement;darkjets",jets.iloc[i]['displacement'],300,0,1e-8)
			plt.plot1D(f"{sample}_jet_r_inv",";r_inv;darkjets",jets.iloc[i]['r_inv'],30,0,1)


			ndark = len(particles_in_jet[particles_in_jet['isdark']==1])
			plt.plot2D(f"{sample}_jet_displaced_vs_ndark",";ndark;displaced tracks;count",ndark,n_displaced,10,0,10,20,0,20)
			plt.plot2D(f"{sample}_jet_prompt_vs_ndark",";ndark;prompt tracks;count",ndark,n_prompt,10,0,10,50,0,50)
			plt.plot2D(f"{sample}_jet_all_vs_ndark",";ndark;tracks;count",ndark,n_displaced+n_prompt,10,0,10,50,0,50)


			
			if n_final!=0:
				n_final_dark = (particles_in_jet['isdark']==1).sum()
				dark_fraction = float(n_final_dark) / float(n_final)
				plt.plot1D(f"{sample}_jet_n_final_dark",";n_final_dark;darkjets",n_final_dark,8,-0.5,7.5)
				plt.plot1D(f"{sample}_jet_fraction_dark",";fraction_dark;darkjets",dark_fraction,20,0,1)

	else:
		print("running code for non dark sample")
		particles_displaced = particles_in_jet[(np.sqrt(particles_in_jet['prod_x']**2 + particles_in_jet['prod_y']**2 )<300.) & (particles_in_jet['prod_z']<300.)]
		particles_displaced = particles_displaced[(np.sqrt(particles_displaced['prod_x']**2 + particles_displaced['prod_y']**2 )>1.)]
		particles_prompt = particles_in_jet[(np.sqrt(particles_in_jet['prod_x']**2 + particles_in_jet['prod_y']**2 )<=1.) & (particles_in_jet['prod_z']<=1.)]
			
		n_displaced = len(particles_displaced)
		n_prompt = len(particles_prompt)

		plt.plot1D(f"{sample}_jet_displaced_tracks",";displaced tracks;jets",n_displaced,20,0,20)
		plt.plot1D(f"{sample}_jet_prompt_tracks",";prompt tracks;jets",n_prompt,50,0,50)


	#if njetsPlotted<njetsPlot:
	#	for j, p in particles_in_jet.iterrows():
	#		plt.plot2D(f"event_display/{sample}_jet_{i}",";eta;phi;count",p['eta'],p['phi'],20,-2.5,2.5,20,0,np.pi)

	#	njetsPlotted += 1
		
#print(particles.head(100))

fout = ROOT.TFile.Open(f"outputs/jets/hists_{sample}_jets.root","RECREATE")
fout.cd()

plt.drawAll()
