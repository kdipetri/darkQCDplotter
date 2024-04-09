'''
Author: Karri DiPetrillo
Date: February 26, 2024
Purpose: Loop over hepmc output file
'''

import pyhepmc
import numpy as np
import argparse
import fastjet

from plothelper import *

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
    outFileName = "outputs/"+sample+".root" # if not ops.inFileName: 

    maxevents = ops.maxEvents

    # pyhepmc.open can read most HepMC formats using auto-detection
    with pyhepmc.open(inFileName) as f:

        # loop over events
        for ievt, event in enumerate(f):

            # loop over particles
            npart = 0
            nleps = 0
            particles = [] # to cluster in jets
            for particle in event.particles:
                
                # check if there's a Z
                if abs(particle.pid)==23 and particle.status==62: 

                    plt.plot1D("{}_z_pt".format(sample)   ,";pt;Zs"  , particle.momentum.pt(), 100, 0, 500)
                    plt.plot1D("{}_z_eta".format(sample)  ,";eta;Zs" , particle.momentum.eta(), 100, -10, 10)
                    plt.plot1D("{}_z_mass".format(sample) ,";mass;Zs", particle.momentum.m()  , 100, 60,125)

                # check if there's a W
                if abs(particle.pid)==24 and particle.status==62: 

                    #print(particle.id, particle.pid, particle.status, particle.momentum.pt(), particle.momentum.eta(), particle.momentum.phi())

                    plt.plot1D("{}_w_pt".format(sample)   ,";pt;Ws"  , particle.momentum.pt(), 100, 0, 500)
                    plt.plot1D("{}_w_eta".format(sample)  ,";eta;Ws" , particle.momentum.eta(), 100, -10, 10)
                    plt.plot1D("{}_w_mass".format(sample) ,";mass;Ws", particle.momentum.m()  , 100, 50,110)

                # check if there's a top
                if abs(particle.pid)==6 and particle.status==62: 

                    #print(particle.id, particle.pid, particle.status, particle.momentum.pt(), particle.momentum.eta(), particle.momentum.phi())

                    plt.plot1D("{}_t_pt".format(sample)   ,";pt;tops"  , particle.momentum.pt(), 100, 0, 500)
                    plt.plot1D("{}_t_eta".format(sample)  ,";eta;tops" , particle.momentum.eta(), 100, -10, 10)
                    plt.plot1D("{}_t_mass".format(sample) ,";mass;tops", particle.momentum.m()  , 100, 140,200)

                # decide which particles to keep
                accept = (particle.status == 1) # final state particle
                accept *= (abs(particle.momentum.eta()) < 2.5) # detector acceptance
                accept *= (abs(particle.momentum.pt()) > 1) # reconstructable pT
                accept *= (abs(particle.momentum.pt()) > 1) # reconstructable pT
                accept *= (abs(particle.pid) not in [12, 14, 16]) # not a neutrino
                if not accept:
                    continue

                #print(particle) 
                #print(particle.id, particle.pid, particle.status, particle.momentum.pt(), particle.momentum.eta(), particle.momentum.phi() )

                plt.plot1D("{}_particle_pid".format(sample) ,";pid;stable particles", abs(particle.pid), 250, 0, 250)
                plt.plot1D("{}_particle_pt".format(sample)  ,";pt;stable particles" , particle.momentum.pt(), 100, 0, 100)
                plt.plot1D("{}_particle_eta".format(sample),";eta;stable particles", particle.momentum.eta(), 100, -3, 3)

                particles.append( fastjet.PseudoJet(particle.momentum.px,particle.momentum.py,particle.momentum.pz,particle.momentum.e) )

                npart += 1

                # decide which leptons to keep
                acceptlep = abs(particle.pid) in [11, 13] 
                acceptlep *= (particle.momentum.pt()>20) # pT

                if not acceptlep:
                    continue

                nleps+=1
                plt.plot1D("{}_lepton_pt".format(sample)  ,";pt;truth leptons", particle.momentum.pt(), 50, 0, 200)
                plt.plot1D("{}_lepton_eta".format(sample),";eta;truth leptons", particle.momentum.eta(), 50, -3, 3)
            
            # jet level plots 
            cluster = fastjet.ClusterSequence(particles, jetdef)

            inc_jets = cluster.inclusive_jets()
            njets = 0
            for jet in inc_jets:

                if jet.pt() < 25 or abs(jet.eta()) > 2.5: continue

                #print("pt:", jet.pt() ,"eta:", jet.eta()) 
                plt.plot1D("{}_jet_pt".format(sample)  ,";pt;truth jets" , jet.pt(), 50, 0, 300)
                plt.plot1D("{}_jet_eta".format(sample),";eta;truth jets", jet.eta(), 50, -3, 3)

                njets += 1

            plt.plot1D("{}_nparticles".format(sample)  ,";n_{particles};events"   , npart, 50, 0, 250)
            plt.plot1D("{}_nleptons".format(sample)    ,";n_{leptons};events"     , nleps,  5, -0.5,  4.5)
            plt.plot1D("{}_njets".format(sample)       ,";n_{jets};events"        , njets, 13, -0.5, 12.5)


            if maxevents!=-1 and ievt > maxevents : break 


        # save histos to file
        fout = ROOT.TFile.Open(outFileName,"RECREATE")
        fout.cd()
        plt.drawAll()

