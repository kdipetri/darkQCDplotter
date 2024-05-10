import argparse
from plothelper import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mplhep as hep # ATLAS/CMS plot style
import pyhepmc
import fastjet
import re


# Set CMS plot style
plt.style.use(hep.style.CMS)



parser = argparse.ArgumentParser(usage=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s", "--sampleName", help="sample name", default="higgs_test")
parser.add_argument("-n", "--maxEvents", help="maximum events to process", default=-1)
parser.add_argument("-ja", "--jetAlgorithm", help="Jet algorithm used", default="anti_kt")
parser.add_argument("-jR", "--jetR", help="Jet R parameter", default=0.4)

ops = parser.parse_args()
maxevents = int(ops.maxEvents)
sample = ops.sampleName
jet_algorithm = ops.jetAlgorithm
R = float(ops.jetR)

# Define the pattern to match numbers
pattern = r"m=(\d+)_xio=(\d+)_xil=(\d+)"

# Use re.findall to extract the numbers
matches = re.findall(pattern, sample)

if matches:
	m,xio,xil = matches[0]


# setup fastjet
if jet_algorithm=='anti_kt':
	jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, R)
if jet_algorithm=="kt":
	jetdef = fastjet.JetDefinition(fastjet.kt_algorithm, R)
if jet_algorithm=="cambridge":
	jetdef = fastjet.JetDefinition(fastjet.cambridge_algorithm, R)

sample = ops.sampleName
inFileName  = "inputs/" +sample+".hepmc" # if not ops.inFileName: 

# Create figure and axis
fig, ax = plt.subplots(figsize=(8, 7))


# check particle is in jet
def inJet(particle,jet,R):
    deta = particle.momentum.eta() - jet.eta()
    dphi = particle.momentum.phi() - jet.phi()
    dR = np.sqrt(deta**2 + dphi**2)
    return dR<R


# pyhepmc.open can read most HepMC formats using auto-detection
with pyhepmc.open(inFileName) as f:

	# loop over events
	for ievt, event in enumerate(f):
		print(f"Event {ievt}")
		# for debugging see only 5 events
		#if ievt>1000:break

		# loop over particles
		npart = 0
		nleps = 0
		particles = [] # to cluster in jets
		for particle in event.particles:
			
			# decide which particles to keep
			accept = (particle.status == 1) # final state particle
			accept *= (abs(particle.momentum.eta()) < 2.5) # detector acceptance
			accept *= (abs(particle.momentum.pt()) > 1) # reconstructable pT
			accept *= (abs(particle.pid) not in [12, 14, 16]) # not a neutrino	
			
			
			# plot scatter before filtering dark particles
			if accept:
				if abs(particle.pid) in [4900113,4900111,4900101,4900021]:
					color = 'red'
				else:
					color='gray'
				ax.scatter(particle.momentum.eta(),particle.momentum.phi(),color=color, alpha=0.5,s=5*particle.momentum.pt())
			
			# filter dark particles
			accept *= (abs(particle.pid) not in [4900113,4900111,4900101,4900021]) # not a dark particle
			
			if not accept:
				continue


			
			particles.append( fastjet.PseudoJet(particle.momentum.px,particle.momentum.py,particle.momentum.pz,particle.momentum.e ) )

			npart += 1

			



		cluster = fastjet.ClusterSequence(particles, jetdef)
		inc_jets = cluster.inclusive_jets()

		# count number of jets
		njets = len(inc_jets)
		
		## Check what dark particles are in the jets
		if njets>0:				
			for i,jet in enumerate(inc_jets):
				if jet.pt() < 25 or abs(jet.eta()) > 2.5: continue

				#ndark = 0
				#ndark_stable = 0
				nparticles = 0
				jet_isdark = 0

				for particle in event.particles:
					# decide which particles to keep
					#accept = (particle.status == 1) # final state particle
					accept = (abs(particle.momentum.eta()) < 2.5) # detector acceptance
					accept *= (abs(particle.momentum.pt()) > 1) # reconstructable pT
					accept *= (abs(particle.pid) not in [12, 14, 16]) # not a neutrino
					
					if not accept:
						continue




					if inJet(particle,jet,R):
						
						particle_isdark = 0
						if particle.status ==1:
							nparticles += 1
						if abs(particle.pid) in [4900113,4900111,4900101,4900021]:
							jet_isdark = 1
							particle_isdark = 1
							#ndark +=1

							#if particle.status ==1:
								#ndark_stable +=1
		
				
				if jet_isdark==0:color='b'
				else: color='r'
		
				# Draw the circle
				if nparticles>3:
					circle = plt.Circle((jet.eta(), jet.phi()), R, color=color, fill=False, linewidth=2,linestyle='--')
					ax.add_patch(circle)

		if matches:
			title = f"Event={ievt}, m={m}, xio={xio}, xil={xil}, R={R} \n {jet_algorithm}"
		else:
			title = sample


				
		# Set plot title and labels
		ax.set_title(title)
		ax.set_xlabel(r'$\eta$')
		ax.set_ylabel(r'$\phi$')
		# Set axis limits
		ax.set_xlim(-3, 3)  # Set X-axis limits
		ax.set_ylim(-3.5, 3.5)  # Set Y-axis limits


		# Add legend
		ax.legend()

		plt.savefig(f"./plots/jets/event_display/{jet_algorithm}/{sample}_event_{ievt}")
		# Clear the plot
		ax.clear()


		if maxevents!=-1 and ievt > maxevents : break 









