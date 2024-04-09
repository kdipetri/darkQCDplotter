'''
Author: Matias Mantinan
Date: April 8, 2024
Purpose: Analyze jets 
'''

import pyhepmc
import numpy as np
import argparse
import fastjet
import math
import ctypes
from plothelper import *


def inJet(particle,jet,R):
    deta = particle.momentum.eta() - jet.eta()
    dphi = particle.momentum.phi() - jet.phi()
    dR = math.sqrt(deta**2 + dphi**2)
    return dR<R

if __name__ == "__main__":
    
    # user options
    parser = argparse.ArgumentParser(usage=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-s", "--sampleName", help="Sample name", default="wjets")
    parser.add_argument("-i", "--inFileName", help="Input file name", default="wjets.hepmc")
    parser.add_argument("-o", "--outFileName", help="Output file name", default="wjets.root")
    parser.add_argument("-n", "--maxEvents", help="maximum events to process", default=-1)
    ops = parser.parse_args()

    # setup plotter
    plt = PlotHelper()

    # setup fastjet
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4)

    sample = ops.sampleName
    inFileName  = "inputs/" +sample+".hepmc" # if not ops.inFileName: 
    outFileName = "outputs/"+sample+"_jets.root" # if not ops.inFileName: 




    maxevents = ops.maxEvents

    # pyhepmc.open can read most HepMC formats using auto-detection
    with pyhepmc.open(inFileName) as f:




        # Create output file to save the jets
        output_file = ROOT.TFile(outFileName, "RECREATE")

        # Create jets TTrees
        jet_tree = ROOT.TTree("Jets", "Jet properties")
        jet_pt = ctypes.c_double()
        jet_eta = ctypes.c_double()
        jet_phi = ctypes.c_double()
        jet_mass = ctypes.c_double()
        jet_r_inv = ctypes.c_double()
        jet_nparticles = ctypes.c_double()
        event_index = ctypes.c_int()
        jet_tree.Branch("pt", jet_pt, "pt/D")
        jet_tree.Branch("eta", jet_eta, "eta/D")
        jet_tree.Branch("phi", jet_phi, "phi/D")
        jet_tree.Branch("mass", jet_mass, "mass/D")
        jet_tree.Branch("r_inv",jet_r_inv,"r_inv/D")
        jet_tree.Branch("nparticles",jet_nparticles,"nparticles/I")
        jet_tree.Branch("event_index", event_index, "event_index/I")

        # Create a TTree to store constituent particles
        particle_tree = ROOT.TTree("Particles", "Particle properties")
        particle_pt = ctypes.c_double()
        particle_eta = ctypes.c_double()
        particle_phi = ctypes.c_double()
        particle_mass = ctypes.c_double()
        jet_index = ctypes.c_int()
        particle_id = ctypes.c_int()
        particle_tree.Branch("pt", particle_pt, "pt/D")
        particle_tree.Branch("eta", particle_eta, "eta/D")
        particle_tree.Branch("phi", particle_phi, "phi/D")
        particle_tree.Branch("mass", particle_mass, "mass/D")
        particle_tree.Branch("pid", particle_id, "pid/I")
        particle_tree.Branch("jet_index", jet_index, "jet_index/I")
        particle_tree.Branch("event_index", event_index, "event_index/I")




        # loop over events
        for ievt, event in enumerate(f):
            # for debugging see only 5 events
            if ievt>10000: 
                break

            # loop over particles
            npart = 0
            nleps = 0
            particles = [] # to cluster in jets
            for particle in event.particles:
                
                # decide which particles to keep
                accept = (particle.status == 1) # final state particle
                accept *= (abs(particle.momentum.eta()) < 2.5) # detector acceptance
                accept *= (abs(particle.momentum.pt()) > 1) # reconstructable pT
                accept *= (abs(particle.momentum.pt()) > 1) # reconstructable pT
                accept *= (abs(particle.pid) not in [12, 14, 16]) # not a neutrino
                accept *= (abs(particle.pid) not in [4900113,4900111,4900101,4900021]) # not a dark particle
                if not accept:
                    continue



                #print(particle) 
                #print(particle.id, particle.pid, particle.status, particle.momentum.pt(), particle.momentum.eta(), particle.momentum.phi() )

                
                particles.append( fastjet.PseudoJet(particle.momentum.px,particle.momentum.py,particle.momentum.pz,particle.momentum.e ) )

                #print(particle.pid)
                #print(particles[-1])
                npart += 1

                



            # jet level plots 
            cluster = fastjet.ClusterSequence(particles, jetdef)



            inc_jets = cluster.inclusive_jets()
            njets = 0

            for jet in inc_jets:
                #jet_constituents = cluster.constituents(jet)
                #for particle in jet_constituents:
                #    print(particle.pid())
                if jet.pt() < 25 or abs(jet.eta()) > 2.5: continue

                #print("pt:", jet.pt() ,"eta:", jet.eta()) 
                #plt.plot1D("{}_jet_pt".format(sample)  ,";pt;truth jets" , jet.pt(), 50, 0, 300)
                #plt.plot1D("{}_jet_eta".format(sample),";eta;truth jets", jet.eta(), 50, -3, 3)

                njets += 1



            

            ## Check what dark particles are in the jets
            if njets>0:

                    
                    
                for jet_index,jet in enumerate(inc_jets):
                    if jet.pt() < 25 or abs(jet.eta()) > 2.5: continue



                    ndark = 0
                    ndark_stable = 0
                    nparticles = 0
                    for particle in event.particles:
                        # decide which particles to keep
                        #accept = (particle.status == 1) # final state particle
                        accept = (abs(particle.momentum.eta()) < 2.5) # detector acceptance
                        accept *= (abs(particle.momentum.pt()) > 1) # reconstructable pT
                        accept *= (abs(particle.momentum.pt()) > 1) # reconstructable pT
                        accept *= (abs(particle.pid) not in [12, 14, 16]) # not a neutrino
                        
                        if not accept:
                            continue
    
                        
                        if inJet(particle,jet,0.4):
            
                            particle_pt = particle.momentum.pt()
                            particle_eta = particle.momentum.eta()
                            particle_phi = particle.momentum.phi()
                            particle_mass = particle.momentum.m()
                            particle_id = particle.pid
                            event_index = ievt
                            jet_index = jet_index

                            particle_tree.Fill()
            
                            nparticles += 1
                            if abs(particle.pid) in [4900113,4900111,4900101,4900021]:
                                ndark +=1
                                print(f"{particle.pid} is in jet")
                                if particle.status ==1:
                                    ndark_stable +=1

                    #print(f"{nparticles} particles in jets")

                    if ndark == 0 :
                        r_inv = 0
                    else:   
                        r_inv = ndark_stable / (ndark_stable + ndark)
                        #print(f"jet has r_inv = {r_inv}")




                    jet_pt = jet.pt()
                    jet_eta = jet.eta()
                    jet_phi = jet.phi()
                    jet_mass = jet.m()
                    jet_r_inv = r_inv
                    jet_nparticles = nparticles
                                        

                    jet_tree.Fill()




            if maxevents!=-1 and ievt > maxevents : break 




        # Write the trees to the ROOT file
        output_file.Write()

        # Close the ROOT file
        output_file.Close()
    
        # save histos to file
        #fout = ROOT.TFile.Open(outFileName,"RECREATE")
        #fout.cd()

