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
from plothelper import *
from array import array


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
    outFileName = "outputs/jets/"+sample +"_jets.root" # if not ops.inFileName: 




    maxevents = ops.maxEvents

    # pyhepmc.open can read most HepMC formats using auto-detection
    with pyhepmc.open(inFileName) as f:




        # Create output file to save the jets
        output_file = ROOT.TFile(outFileName, "RECREATE")

        # Create jets TTrees
        jet_tree = ROOT.TTree("Jets", "Jet properties")
        jet_pt = array('d',[0])
        jet_eta = array('d',[0])
        jet_phi = array('d',[0])
        jet_mass = array('d',[0])
        jet_r_inv = array('d',[0])
        jet_displacement = array('d',[0])
        jet_nparticles = array('i',[0])
        jet_ndark = array('i',[0])
        jet_isdark = array('i',[0])
        event_index = array('i',[0])
        
        jet_tree.Branch("pt", jet_pt, "pt/D")
        jet_tree.Branch("eta", jet_eta, "eta/D")
        jet_tree.Branch("phi", jet_phi, "phi/D")
        jet_tree.Branch("mass", jet_mass, "mass/D")
        jet_tree.Branch("r_inv",jet_r_inv,"r_inv/D")
        jet_tree.Branch("displacement",jet_displacement,"displacement/D")
        jet_tree.Branch("nparticles",jet_nparticles,"nparticles/I")
        jet_tree.Branch("ndark",jet_ndark,"ndark/I")
        jet_tree.Branch("isdark",jet_isdark,"isdark/I")
        jet_tree.Branch("event_index", event_index, "event_index/I")

        # Create a TTree to store constituent particles
        particle_tree = ROOT.TTree("Particles", "Particle properties")
        particle_pt = array('d',[0])
        particle_eta = array('d',[0])
        particle_phi = array('d',[0])
        particle_mass = array('d',[0])
        particle_prod_x = array('d',[0])
        particle_prod_y = array('d',[0])
        particle_prod_z = array('d',[0])
        particle_prod_t = array('d',[0])
        particle_end_x = array('d',[0])
        particle_end_y = array('d',[0])
        particle_end_z = array('d',[0])
        particle_end_t = array('d',[0])
        particle_displacement = array('d',[0])
        particle_isdark = array('i',[0])
        jet_index = array('i',[0])
        particle_id = array('i',[0])
        particle_status = array('i',[0])
        
        

        particle_tree.Branch("pt", particle_pt, "pt/D")
        particle_tree.Branch("eta", particle_eta, "eta/D")
        particle_tree.Branch("phi", particle_phi, "phi/D")
        particle_tree.Branch("mass", particle_mass, "mass/D")
        particle_tree.Branch("prod_x",particle_prod_x,"prod_x/D")
        particle_tree.Branch("prod_y",particle_prod_y,"prod_y/D")
        particle_tree.Branch("prod_z",particle_prod_z,"prod_z/D")
        particle_tree.Branch("prod_t",particle_prod_t,"prod_t/D")
        particle_tree.Branch("displacement",particle_displacement,"displacement/D")
        particle_tree.Branch("end_x",particle_end_x,"end_x/D")
        particle_tree.Branch("end_y",particle_end_y,"end_y/D")
        particle_tree.Branch("end_z",particle_end_z,"end_z/D")
        particle_tree.Branch("end_t",particle_end_t,"end_t/D")
        particle_tree.Branch("isdark",particle_isdark,"isdark/I")
        particle_tree.Branch("pid", particle_id, "pid/I")
        particle_tree.Branch("status", particle_status, "status/I")
        particle_tree.Branch("jet_index", jet_index, "jet_index/I")
        particle_tree.Branch("event_index", event_index, "event_index/I")


        jet_index[0] = 0

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

                    
                    
                for i,jet in enumerate(inc_jets):
                    #if i>100: break # lets start with the first 10 jets
                    if jet.pt() < 25 or abs(jet.eta()) > 2.5: continue

                    #print(f"jet number {jet_index[0]}")


                    ndark = 0
                    ndark_stable = 0
                    nparticles = 0
                    jet_isdark[0] = 0
                    jet_displacement[0] = 0.


                    for particle in event.particles:
                        # decide which particles to keep
                        #accept = (particle.status == 1) # final state particle
                        accept = (abs(particle.momentum.eta()) < 2.5) # detector acceptance
                        accept *= (abs(particle.momentum.pt()) > 1) # reconstructable pT
                        accept *= (abs(particle.pid) not in [12, 14, 16]) # not a neutrino
                        
                        if not accept:
                            continue
    



                        if inJet(particle,jet,0.4):
                            
                            particle_isdark[0] = 0
                            nparticles += 1
                            if abs(particle.pid) in [4900113,4900111,4900101,4900021]:
                                jet_isdark[0] = 1
                                particle_isdark[0] = 1
                                ndark +=1

                                if particle.status ==1:
                                    ndark_stable +=1

                            
                            
                            
                            particle_pt[0] = particle.momentum.pt()
                            particle_eta[0] = particle.momentum.eta()
                            particle_phi[0] = particle.momentum.phi()
                            particle_mass[0] = particle.momentum.m()
                            particle_id[0] = particle.pid
                            particle_status[0] = particle.status
                            event_index[0] = ievt
                            

                            
                            prod_vertex = particle.production_vertex
                            end_vertex = particle.end_vertex



                            if not (prod_vertex is None):
                                particle_prod_x[0] = prod_vertex.position.x
                                particle_prod_y[0] = prod_vertex.position.y
                                particle_prod_z[0] = prod_vertex.position.z
                                particle_prod_t[0] = prod_vertex.position.y
                            else:
                                particle_prod_x[0] = 0
                                particle_prod_y[0] = 0
                                particle_prod_z[0] = 0
                                particle_prod_t[0] = 0
                                

                            if not (end_vertex is None):
                                particle_end_x[0] = end_vertex.position.x
                                particle_end_y[0] = end_vertex.position.y
                                particle_end_z[0] = end_vertex.position.z
                                particle_end_t[0] = end_vertex.position.y
                            else:
                                particle_end_x[0] = 0
                                particle_end_y[0] = 0
                                particle_end_z[0] = 0
                                particle_end_t[0] = 0
                            
                            displacement = np.sqrt(particle_end_x[0]**2 + particle_end_y[0]**2 +particle_end_z[0]**2)
                            particle_displacement[0] = displacement

                            particle_tree.Fill()

                            if particle_id[0] in [4900101,4900021] and particle_status[0]!=1:
                                displacement = np.sqrt(particle_end_x[0]**2 + particle_end_y[0]**2 +particle_end_z[0]**2)
                                if jet_displacement[0] == 0.:
                                    jet_displacement[0] = displacement
                                else:   
                                    jet_displacement[0] = min(jet_displacement[0],displacement)



                            

                    #print(f"{nparticles} particles in jets")

                    if ndark == 0 :
                        r_inv = 0.
                    else:   
                        r_inv = float(ndark_stable) / float(ndark)
                        #print(f"jet has r_inv = {r_inv}")




                    jet_pt[0] = jet.pt()
                    jet_eta[0] = jet.eta()
                    jet_phi[0] = jet.phi()
                    jet_mass[0] = jet.m()
                    jet_r_inv[0] = r_inv
                    jet_nparticles[0] = nparticles
                    jet_ndark[0] = ndark
                    

                    if nparticles > 0:
                        jet_tree.Fill()

                        jet_index[0] +=1

                    


            if maxevents!=-1 and ievt > maxevents : break 




        # Write the trees to the ROOT file
        output_file.Write()

        # Close the ROOT file
        output_file.Close()
    
        # save histos to file
        #fout = ROOT.TFile.Open(outFileName,"RECREATE")
        #fout.cd()

