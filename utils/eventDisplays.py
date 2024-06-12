
'''
Author: Karri DiPetrillo
Date: June 11, 2024
Purpose: make event displays
'''

import pyhepmc
import numpy as np
import argparse
#import fastjet
import pyjet
from particle import Particle
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def decaysToSelf(particle) :

	if particle.end_vertex: 
		for child in particle.end_vertex.particles_out:
			if child.pid==particle.pid : 
				return True
				break
	
	return False

charged_pids = [543, 5334, 5314, 20413, 4122, 4124]
neutral_pids = [35, 5212, 5214, 551]
# Get charge
def isCharged(particle):
	if abs(particle.pid) > 1000000 : return False
	if abs(particle.pid) > 5122 and abs(particle.pid) < 5555 : return False # bottom
	if abs(particle.pid) > 4122 and abs(particle.pid) < 4445 : return False # bottom
	if abs(particle.pid) in charged_pids : return True
	if abs(particle.pid) in neutral_pids : return False
	part = Particle.from_pdgid(particle.pid)
	if part.charge !=0 : return True
	else : return False

def higgsParent(particle):

	if particle.production_vertex: 
		for parent in particle.production_vertex.particles_in:
			print(parent.pid)
			if parent.pid==25: 
				return True
			else : return higgsParent(parent)

	return False

if __name__ == "__main__":
	
	# user options
	parser = argparse.ArgumentParser(usage=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-s", "--sampleName", help="Sample name", default="vector_m-1_ctau-10")
	parser.add_argument("-n", "--maxEvents", help="maximum events to process", default=10)
	ops = parser.parse_args()

	# setup fastjet
	#jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4)

	sample = ops.sampleName
	inFileName  = "inputs/" +sample+".hepmc" # if not ops.inFileName: 
	outFileName = "outputs/displays_"+sample+".root" # if not ops.inFileName: 

	maxevents = ops.maxEvents
		

	# pyhepmc.open can read most HepMC formats using auto-detection
	with pyhepmc.open(inFileName) as f:

		# loop over events
		for ievt, event in enumerate(f):

			print("Event ", ievt)

			if maxevents!=-1 and ievt > maxevents : break 

			# loop over particles
			npart = 0
			ncharged = 0
			ndisplaced = 0
			nleps = 0
			particles = [] # to cluster in jets

			# start with eta phi plots
			from_scalar_e = []
			from_scalar_mu = []
			from_scalar_gamma = []
			from_scalar_chad = []
			from_scalar_nhad = []
			from_scalar_inv = []

			from_other_e = []
			from_other_mu = []
			from_other_gamma = []
			from_other_chad = []
			from_other_nhad = []

			higgs_pt = 0
			higgs_eta = -99
			higgs_phi = -99

			jets = []

			for particle in event.particles:
				
				# check if there's a higgs
				if abs(particle.pid)==25: # just take the last one # and particle.status==62:

					higgs_pt  =  particle.momentum.pt() 
					higgs_eta =  particle.momentum.eta() 
					higgs_phi =  particle.momentum.phi() 


				# check for dark mesons
				#if abs(particle.pid) in [4900111, 4900113] and not decaysToSelf(particle): 
				#
				#    prodvtx = particle.production_vertex
				#    prodvec = prodvtx.position
				#    prodrxy   = (prodvec.x**2 + prodvec.y**2)**0.5 
				#
				#    plt.plot1D("{}_dark_pt".format(sample)   ,";pt;dark"  , particle.momentum.pt(), 100, 0, 500)
				#    plt.plot1D("{}_dark_eta".format(sample)  ,";eta;dark" , particle.momentum.eta(), 100, -10, 10)
				#    plt.plot1D("{}_dark_mass".format(sample) ,";mass;dark", particle.momentum.m()  , 100, 0,60)
				#    plt.plot1D("{}_dark_status".format(sample) ,";status;dark", particle.status  , 100,-0.5,99.5)
				#    plt.plot1D("{}_dark_prodrxy".format(sample) ,";prodrxy;dark", prodrxy  , 100,0,1000)
				#
				#    if particle.end_vertex: 
				#        decayvec = particle.end_vertex.position
				#        decayrxy = (decayvec.x**2 + decayvec.y**2)**0.5 
				#        #for child in particle.end_vertex.particles_out:
				#        #    print("     ", child.id, child.pid, child.status, child.momentum.pt(), child.momentum.eta(), child.momentum.phi())
				#
				#        plt.plot1D("{}_dark_decayrxy".format(sample) ,";decayrxy;dark", decayrxy  , 100,0,1000)



				# keep only stable reconstructable particles
				if particle.status != 1: continue # final state particle
				if abs(particle.momentum.eta()) > 2.5: continue # detector acceptance
				if abs(particle.momentum.pt()) < 1: continue # reconstructable pT

				prodvtx = particle.production_vertex
				prodvec = prodvtx.position
				prodrxy = (prodvec.x**2 + prodvec.y**2)**0.5 
				
				charged = isCharged(particle)
				charged *= (prodrxy < 300)
				displaced = (prodrxy > 0.5 and charged)

				# now do some counting
				pid = abs(particle.pid)
				pt  = particle.momentum.pt()
				eta = particle.momentum.eta()
				phi = particle.momentum.phi()
				fromHiggs = higgsParent(particle)

				accept =  (abs(particle.pid) not in [12, 14, 16]) # not a neutrino
				accept *= (abs(particle.pid) < 1000000)# not a dark sector particle     
				accept *= (abs(prodrxy)<1000) # produced within 1 m

				if not accept:
					continue

				npart += 1
				if charged : ncharged+=1
				if displaced : ndisplaced += 1


				if fromHiggs: 
					
					if      pid==11 : from_scalar_e.append( (pt,eta,phi) )
					elif 	pid==13 : from_scalar_mu.append( (pt,eta,phi) ) 
					elif 	pid==22 : from_scalar_gamma.append( (pt,eta,phi) ) 
					elif 	charged : from_scalar_chad.append( (pt,eta,phi) )
					else 			: from_scalar_nhad.append( (pt,eta,phi) )

				else : 
					if      pid==11 : from_other_e.append( (pt,eta,phi) )
					elif 	pid==13 : from_other_mu.append( (pt,eta,phi) ) 
					elif 	pid==22 : from_other_gamma.append( (pt,eta,phi) ) 
					elif 	charged : from_other_chad.append( (pt,eta,phi) )
					else 			: from_other_nhad.append( (pt,eta,phi) )

				#particles.append( fastjet.PseudoJet(particle.momentum.px,particle.momentum.py,particle.momentum.pz,particle.momentum.e ) )
				# should do proper overlap removal
				if particle.pid not in [11,13,15] or particle.momentum.pt()<20 : # remove high pT leptons
					particles.append( (particle.momentum.e,particle.momentum.px,particle.momentum.py,particle.momentum.pz) )
			

			# format structured arrays for plotting
			dt=np.dtype([('pt', 'f8'), ('eta', 'f8'), ('phi', 'f8')])

			arr_from_scalar_e		 = np.array(from_scalar_e, dt) 
			arr_from_scalar_mu		 = np.array(from_scalar_mu, dt) 
			arr_from_scalar_gamma	 = np.array(from_scalar_gamma, dt) 
			arr_from_scalar_chad	 = np.array(from_scalar_chad, dt) 
			arr_from_scalar_nhad	 = np.array(from_scalar_nhad, dt) 
			#arr_from_scalar_inv		 = np.array(from_scalar_inv, dt) 

			arr_from_other_e		 = np.array(from_other_e, dt) 
			arr_from_other_mu		 = np.array(from_other_mu, dt) 
			arr_from_other_gamma	 = np.array(from_other_gamma, dt) 
			arr_from_other_chad		 = np.array(from_other_nhad, dt) 
			arr_from_other_nhad		 = np.array(from_other_chad, dt) 

			# cluster jets
	
			# format structured array for pyjet input 
			dtjet=np.dtype([('E', 'f8'), ('px', 'f8'), ('py', 'f8'), ('pz', 'f8')])
			arr = np.array(particles, dtjet)

			sequence = pyjet.cluster(arr, R=1.0, p=-1, ep=True)
			inc_jets = sequence.inclusive_jets()  # list of PseudoJets

			njets = 0
			for jet in inc_jets:

				if jet.pt < 25 or abs(jet.eta) > 2.5: continue

				jets.append( ( jet.pt, jet.eta, jet.phi ) ) 

				njets += 1

			arr_jets = np.array(jets, dt)

			#plt.plot1D("{}_nparticles".format(sample)  ,";n_{particles};events"   , npart, 50, 0, 250)
			#plt.plot1D("{}_ncharged".format(sample)  ,";n_{charged};events"   , ncharged, 50, 0, 200)
			#plt.plot1D("{}_ndisplaced".format(sample)  ,";n_{displaced};events"   , ndisplaced, 50, 0, 50)
			#plt.plot1D("{}_nleptons".format(sample)    ,";n_{leptons};events"     , nleps,  5, -0.5,  4.5)
			#plt.plot1D("{}_njets".format(sample)       ,";n_{jets};events"        , njets, 13, -0.5, 12.5)


			# now make event display
			fig = plt.figure(figsize=(8,8))
			ax = plt.gca()

			# Plot parameters
			ax.set_xlim(-3.5,3.5)
			ax.set_ylim(-3.5,3.5)
			ax.set_xlabel(r'$\eta$', fontsize=20)
			ax.set_ylabel(r'$\phi$', fontsize=20)
			ax.tick_params(axis='both', which='major', labelsize=16)
			
			# Function that sets the scaling for markers
			# Two methods for the moment: use tanh or square.
			# Scaling using just the energy is also an option
			def scale(particles, higgs_pt, method=0):
				"""Just to scale to a reasonable dot size"""
				energies = particles["pt"]
				if len(energies) == 0: return []
				if method == 0:
					e_normed = energies*5
				elif method == 1: 
					e_normed = 500000.*np.square(energies/higgs_pt)
				elif method == 2:
					e_normed = 1000.*np.tanh(energies/higgs_pt)
				else:
					e_normed = 2500.*energies/higgs_pt
				return e_normed

			# Add scatters to figure
			ax.scatter(arr_from_scalar_e["eta"], arr_from_scalar_e["phi"],
					   s=scale(arr_from_scalar_e,higgs_pt),
					   c='xkcd:blue', marker='p')
			ax.scatter(arr_from_scalar_mu["eta"], arr_from_scalar_mu["phi"],
					   s=scale(arr_from_scalar_mu,higgs_pt),
					   c='xkcd:blue', marker='^')
			ax.scatter(arr_from_scalar_gamma["eta"], arr_from_scalar_gamma["phi"],
					   s=scale(arr_from_scalar_gamma,higgs_pt),
					   c='xkcd:blue', marker='s')
			ax.scatter(arr_from_scalar_chad["eta"], arr_from_scalar_chad["phi"],
					   s=scale(arr_from_scalar_chad,higgs_pt),
					   c='xkcd:blue', marker='o')
			ax.scatter(arr_from_scalar_nhad["eta"], arr_from_scalar_nhad["phi"],
					   s=scale(arr_from_scalar_nhad,higgs_pt),
					   edgecolors='xkcd:blue', marker='o', facecolors='none')
			ax.scatter(arr_from_other_e["eta"], arr_from_other_e["phi"],
					   s=scale(arr_from_other_e,higgs_pt),
					   c='xkcd:magenta', marker='p')
			ax.scatter(arr_from_other_mu["eta"], arr_from_other_mu["phi"],
					   s=scale(arr_from_other_mu,higgs_pt),
					   c='xkcd:magenta', marker='^')
			ax.scatter(arr_from_other_gamma["eta"], arr_from_other_gamma["phi"],
					   s=scale(arr_from_other_gamma,higgs_pt),
					   c='xkcd:magenta', marker='s')
			ax.scatter(arr_from_other_chad["eta"], arr_from_other_chad["phi"],
					   s=scale(arr_from_other_chad,higgs_pt),
					   c='xkcd:magenta', marker='o')
			ax.scatter(arr_from_other_nhad["eta"], arr_from_other_nhad["phi"],
					   s=scale(arr_from_other_nhad,higgs_pt),
					   edgecolors='xkcd:magenta', marker='o', facecolors='none')

			ax.scatter(arr_jets["eta"], arr_jets["phi"],
					   s=2500, edgecolors='xkcd:gray', marker='o',facecolors='none')

			# Add the higgs mediator to the plot
			ax.scatter(higgs_eta, higgs_phi,
					   s=5*higgs_pt, marker='*', c='xkcd:goldenrod')

			# Legend 1 is particle type
			line0 = ax.scatter([-100], [-100], label='higgs mediator', marker='*', c='xkcd:goldenrod')
			line1 = ax.scatter([-100], [-100], label='$e$',marker='p', c='xkcd:black')
			line2 = ax.scatter([-100], [-100], label='$\mu$', marker='^', c='xkcd:black')
			line3 = ax.scatter([-100], [-100], label='$\gamma$', marker='s', c='xkcd:black')
			line4 = ax.scatter([-100], [-100], label='charged hadron', marker='o', c='xkcd:black')
			line5 = ax.scatter([-100], [-100], label='neutral hadron', marker='o', edgecolors='xkcd:black', facecolors='none')
			line6 = ax.scatter([-100], [-100], label='AK4 jets', marker='o',
							   facecolors='none', edgecolors='xkcd:gray')
			#line8 = ax.scatter([-100], [-100], label='AK15 jets', marker='o',
			#				   facecolors='none', edgecolors='xkcd:bright yellow')
			first_legend = plt.legend(handles=[line0, line6, line1, line2, line3, line4, line5],
									  loc='upper right', fontsize=12)
			ax.add_artist(first_legend)

			# Legend 2 is about particle origin
			blue_patch = mpatches.Patch(color='xkcd:blue', label='from higgs')
			magenta_patch = mpatches.Patch(color='xkcd:magenta', label='not from higgs')
			plt.legend(handles=[blue_patch, magenta_patch],loc='upper left')

			# build a rectangle in axes coords
			left, width = .0, 1.
			bottom, height = .0, 1.
			center = left + width/2.
			right = left + width
			top = bottom + height

			# axes coordinates are 0,0 is bottom left and 1,1 is upper right
			p = mpatches.Rectangle((left, bottom), width, height,
				fill=False, transform=ax.transAxes, clip_on=False)

			ax.add_patch(p)
			# Print event number
			ax.text(left, top, 'Event %d'%ievt, horizontalalignment='left',
					verticalalignment='bottom', transform=ax.transAxes, fontsize=12)
			## Print sample details
			#ax.text(right, top, 'mMed=%d$\,$GeV,mDark=%d$\,$GeV,T=%d$\,$K,'
			#        '%s'%(mMed,mDark,temp,decayMode), horizontalalignment='right',
			#        verticalalignment='bottom', transform=ax.transAxes, fontsize=12)
			## Print details about cuts
			#ax.text(left+0.02, bottom+0.01, 'Final particles have $P_{T}>1\,$GeV',
			#        horizontalalignment='left', verticalalignment='bottom',
			#        transform=ax.transAxes, fontsize=12)
			# Print details of higgs mediator
			if higgs_pt > 0 : 
				ax.text(left+0.02, bottom+0.05, 'Higgs $p_{T}=%d\,$GeV'%(higgs_pt),
						horizontalalignment='left', verticalalignment='bottom',
						transform=ax.transAxes, fontsize=12)

			fig.savefig('plots/displays/%s-evt%i.png'%(sample, ievt))

		# save histos to file
		#fout = ROOT.TFile.Open(outFileName,"RECREATE")
		#fout.cd()
		#plt.drawAll()

