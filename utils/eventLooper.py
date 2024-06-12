'''
Author: Karri DiPetrillo
Date: February 26, 2024
Purpose: Loop over hepmc output file
'''

import pyhepmc
import numpy as np
import argparse
#import fastjet
import pyjet
from particle import Particle
from plothelper import *


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
    #jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4)

    sample = ops.sampleName
    inFileName  = "inputs/" +sample+".hepmc" # if not ops.inFileName: 
    outFileName = "outputs/"+sample+".root" # if not ops.inFileName: 

    maxevents = ops.maxEvents
        


    # pyhepmc.open can read most HepMC formats using auto-detection
    with pyhepmc.open(inFileName) as f:

        # loop over events
        for ievt, event in enumerate(f):

            if ievt % 1000 == 0: print("Event ", ievt)

            #if ievt > 10 : break # for debugging
            #print(ievt)

            # loop over particles
            npart = 0
            ncharged = 0
            ndisplaced = 0
            nleps = 0
            particles = [] # to cluster in jets

            for particle in event.particles:
                
                # check if there's a higgs
                if abs(particle.pid)==25 and particle.status==62:

                    plt.plot1D("{}_h_pt".format(sample)   ,";pt;Hs"  , particle.momentum.pt(), 100, 0, 500)
                    plt.plot1D("{}_h_eta".format(sample)  ,";eta;Hs" , particle.momentum.eta(), 100, -10, 10)
                    plt.plot1D("{}_h_mass".format(sample) ,";mass;Hs", particle.momentum.m()  , 100, 100,150)

                # check for dark mesons
                if abs(particle.pid) in [4900111, 4900113] and not decaysToSelf(particle): 

                    #print(particle.id, particle.pid, particle.status, particle.momentum.pt(), particle.momentum.eta(), particle.momentum.phi())
                    prodvtx = particle.production_vertex
                    prodvec = prodvtx.position
                    prodrxy   = (prodvec.x**2 + prodvec.y**2)**0.5 


                    plt.plot1D("{}_dark_pt".format(sample)   ,";pt;dark"  , particle.momentum.pt(), 100, 0, 500)
                    plt.plot1D("{}_dark_eta".format(sample)  ,";eta;dark" , particle.momentum.eta(), 100, -10, 10)
                    plt.plot1D("{}_dark_mass".format(sample) ,";mass;dark", particle.momentum.m()  , 100, 0,60)
                    plt.plot1D("{}_dark_status".format(sample) ,";status;dark", particle.status  , 100,-0.5,99.5)
                    plt.plot1D("{}_dark_prodrxy".format(sample) ,";prodrxy;dark", prodrxy  , 100,0,1000)

                    if particle.end_vertex: 
                        decayvec = particle.end_vertex.position
                        decayrxy = (decayvec.x**2 + decayvec.y**2)**0.5 
                        #for child in particle.end_vertex.particles_out:
                        #    print("     ", child.id, child.pid, child.status, child.momentum.pt(), child.momentum.eta(), child.momentum.phi())

                        plt.plot1D("{}_dark_decayrxy".format(sample) ,";decayrxy;dark", decayrxy  , 100,0,1000)

                # check if there's a Z
                if abs(particle.pid)==23 and particle.status==62: 

                    plt.plot1D("{}_z_pt".format(sample)   ,";pt;Zs"  , particle.momentum.pt(), 100, 0, 500)
                    plt.plot1D("{}_z_eta".format(sample)  ,";eta;Zs" , particle.momentum.eta(), 100, -10, 10)
                    plt.plot1D("{}_z_mass".format(sample) ,";mass;Zs", particle.momentum.m()  , 100, 60,125)

                # check if there's a W
                if abs(particle.pid)==24 and particle.status==62: 


                    plt.plot1D("{}_w_pt".format(sample)   ,";pt;Ws"  , particle.momentum.pt(), 100, 0, 500)
                    plt.plot1D("{}_w_eta".format(sample)  ,";eta;Ws" , particle.momentum.eta(), 100, -10, 10)
                    plt.plot1D("{}_w_mass".format(sample) ,";mass;Ws", particle.momentum.m()  , 100, 50,110)

                # check if there's a top
                if abs(particle.pid)==6 and particle.status==62: 

                    #print(particle.id, particle.pid, particle.status, particle.momentum.pt(), particle.momentum.eta(), particle.momentum.phi())

                    plt.plot1D("{}_t_pt".format(sample)   ,";pt;tops"  , particle.momentum.pt(), 100, 0, 500)
                    plt.plot1D("{}_t_eta".format(sample)  ,";eta;tops" , particle.momentum.eta(), 100, -10, 10)
                    plt.plot1D("{}_t_mass".format(sample) ,";mass;tops", particle.momentum.m()  , 100, 140,200)

                # keep only stable reconstructable particles
                if particle.status != 1: continue # final state particle
                if abs(particle.momentum.eta()) > 2.5: continue # detector acceptance
                if abs(particle.momentum.pt()) < 1: continue # reconstructable pT

                prodvtx = particle.production_vertex
                prodvec = prodvtx.position
                rxy   = (prodvec.x**2 + prodvec.y**2)**0.5 
                
                #print(charged)
                charged = isCharged(particle)
                charged *= (rxy < 300)
                displaced = (rxy > 0.5 and charged)

                # now do some counting
                accept =  (abs(particle.pid) not in [12, 14, 16]) # not a neutrino
                accept *= (abs(particle.pid) < 1000000)# not a dark sector particle     
                accept *= (abs(rxy)<1000) # produced within 1 m

                if not accept:
                    continue

                npart += 1
                if charged : ncharged+=1
                if displaced : ndisplaced += 1


                #print(particle) 
                #print(particle.id, particle.pid, particle.status, particle.momentum.pt(), particle.momentum.eta(), particle.momentum.phi() )

                plt.plot1D("{}_particle_pid".format(sample) ,";pid;stable particles", abs(particle.pid), 250, 0, 250)
                plt.plot1D("{}_particle_pt".format(sample)  ,";pt;stable particles" , particle.momentum.pt(), 100, 0, 100)
                plt.plot1D("{}_particle_eta".format(sample),";eta;stable particles", particle.momentum.eta(), 100, -3, 3)
                plt.plot1D("{}_particle_prodrxy".format(sample),";rxy;stable particles", rxy, 100, 0,1000)
                # separate into charged, displaced, etc

                if charged : 
                    plt.plot1D("{}_charged_pid".format(sample) ,";pid;charged particles", abs(particle.pid), 250, 0, 250)
                    plt.plot1D("{}_charged_pt".format(sample)  ,";pt;charged particles" , particle.momentum.pt(), 100, 0, 100)
                    plt.plot1D("{}_charged_eta".format(sample),";eta;charged particles", particle.momentum.eta(), 100, -3, 3)
                    plt.plot1D("{}_charged_prodrxy".format(sample),";rxy;charged particles", rxy, 100, 0,1000)

                if displaced:
                    plt.plot1D("{}_displaced_pid".format(sample) ,";pid;displaced particles", abs(particle.pid), 250, 0, 250)
                    plt.plot1D("{}_displaced_pt".format(sample)  ,";pt;displaced particles" , particle.momentum.pt(), 100, 0, 100)
                    plt.plot1D("{}_displaced_eta".format(sample),";eta;displaced particles", particle.momentum.eta(), 100, -3, 3)
                    plt.plot1D("{}_displaced_prodrxy".format(sample),";rxy;displaced particles", rxy, 100, 0,300)

                #particles.append( fastjet.PseudoJet(particle.momentum.px,particle.momentum.py,particle.momentum.pz,particle.momentum.e ) )
                # should do proper overlap removal
                if particle.pid not in [11,13,15] or particle.momentum.pt()<20 : # remove high pT leptons
                    particles.append( (particle.momentum.e,particle.momentum.px,particle.momentum.py,particle.momentum.pz) )



                # decide which leptons to keep
                acceptlep = abs(particle.pid) in [11, 13] 
                acceptlep *= (particle.momentum.pt()>20) # pT

                if not acceptlep:
                    continue

                nleps+=1
                plt.plot1D("{}_lepton_pt".format(sample)  ,";pt;truth leptons", particle.momentum.pt(), 50, 0, 200)
                plt.plot1D("{}_lepton_eta".format(sample),";eta;truth leptons", particle.momentum.eta(), 50, -3, 3)
            
            # jet level plots 
            #cluster = fastjet.ClusterSequence(particles, jetdef)
            #inc_jets = cluster.inclusive_jets()
    
            # format structured array for pyjet input 
            dt=np.dtype([('E', 'f8'), ('px', 'f8'), ('py', 'f8'), ('pz', 'f8')])
            arr = np.array(particles, dt)

            sequence = pyjet.cluster(arr, R=1.0, p=-1, ep=True)
            inc_jets = sequence.inclusive_jets()  # list of PseudoJets

            njets = 0
            for jet in inc_jets:

                if jet.pt < 25 or abs(jet.eta) > 2.5: continue

                #print("pt:", jet.pt() ,"eta:", jet.eta()) 
                plt.plot1D("{}_jet_pt".format(sample)  ,";pt;truth jets" , jet.pt, 50, 0, 300)
                plt.plot1D("{}_jet_eta".format(sample),";eta;truth jets", jet.eta, 50, -3, 3)

                njets += 1

            plt.plot1D("{}_nparticles".format(sample)  ,";n_{particles};events"   , npart, 50, 0, 250)
            plt.plot1D("{}_ncharged".format(sample)  ,";n_{charged};events"   , ncharged, 50, 0, 200)
            plt.plot1D("{}_ndisplaced".format(sample)  ,";n_{displaced};events"   , ndisplaced, 50, 0, 50)
            plt.plot1D("{}_nleptons".format(sample)    ,";n_{leptons};events"     , nleps,  5, -0.5,  4.5)
            plt.plot1D("{}_njets".format(sample)       ,";n_{jets};events"        , njets, 13, -0.5, 12.5)


            if maxevents!=-1 and ievt > maxevents : break 


        # save histos to file
        fout = ROOT.TFile.Open(outFileName,"RECREATE")
        fout.cd()
        plt.drawAll()

