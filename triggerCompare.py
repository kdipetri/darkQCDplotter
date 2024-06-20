'''
Author: Matias Mantinan
Date: May 29, 2024
Purpose: Compare triggers 
'''


import pyhepmc
import numpy as np
import fastjet
import matplotlib.pyplot as plt
import mplhep as hep # ATLAS/CMS plot style
import glob
import multiprocessing as mp


# Set CMS plot style
plt.style.use(hep.style.CMS)

# setup fastjet
jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4)


def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=50, fill='█'):
    """
    Call in a loop to create a progress bar.
    @param iteration: Required: current iteration (int)
    @param total: Required: total iterations (int)
    @param prefix: Optional: prefix string (str)
    @param suffix: Optional: suffix string (str)
    @param decimals: Optional: positive number of decimals in percent complete (int)
    @param length: Optional: character length of bar (int)
    @param fill: Optional: bar fill character (str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end='\r')
    # Print New Line on Complete
    if iteration == total: 
        print()



def plot_event(event,ievt,trigger,sample):
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(8, 7))

    for particle in event.particles:
        
        # decide which particles to keep
        accept = (abs(particle.momentum.eta()) < 2.5) # detector acceptance
        accept *= (abs(particle.momentum.pt()) > 1) # reconstructable pT
        accept *= (abs(particle.pid) not in [12, 14, 16]) # not a neutrino	

        if accept and particle.pid == 25:
            ax.scatter(particle.momentum.eta(),particle.momentum.phi(),color="blue", alpha=0.5,s=5*particle.momentum.pt())


        accept *= (particle.status == 1) # final state particle
        
        
        # plot scatter before filtering dark particles
        if accept:
            if abs(particle.pid) in [4900113,4900111,4900101,4900021]:
                color = 'red'
            else:
                color='gray'
            ax.scatter(particle.momentum.eta(),particle.momentum.phi(),color=color, alpha=0.5,s=5*particle.momentum.pt())
        

    # Set plot title and labels
    ax.set_title("Event display")
    ax.set_xlabel(r'$\eta$')
    ax.set_ylabel(r'$\phi$')
    # Set axis limits
    ax.set_xlim(-3, 3)  # Set X-axis limits
    ax.set_ylim(-3.5, 3.5)  # Set Y-axis limits


    # Add legend
    ax.legend()

    sample = sample.replace(".0","")

    plt.savefig(f"./plots/trigger/event_display/{trigger}/{sample}_event_{ievt}")
    # Clear the plot
    ax.clear()



def xs(sample):
    if "qcd" in sample: return 1e8
    if "wjets" in sample: return 1e5
    if "zjets" in sample: return 3e4
    if "ttbar" in sample: return 5e2
    if "higgs" in sample: return 1e2
    else: return 1e2

def compare_triggers(sample, maxevents=-1):

    inFileName  = "inputs/" +sample+".hepmc" # if not ops.inFileName: 



    # pyhepmc.open can read most HepMC formats using auto-detection
    with pyhepmc.open(inFileName) as f:

        events_missing_trigger = 0 # events have missing pt > 200 GeV
        events_dilepton_trigger = 0 # events that DONT have two electrons with at least 15 GeV
        events_HT_trigger = 0 # events that have jet with more than 1 TeV
        events_single_lepton = 0 # single lepton trigger


        # loop over events
        for ievt, event in enumerate(f):

            # loop over particles
            npart = 0
            particles = [] # to cluster in jets
            
            # missing pt comming from neutrinos and dark 
            missing_px = 0
            missing_py = 0

            leptons_15 = 0 # number of leptons with pt more than 15 GeV
            leptons_27 = 0 # number of leptons with pt more than 27 GeV

            for particle in event.particles:

                # check invisible particles missing pT
                accept = (particle.status == 1) # final state particle
                accept *= (abs(particle.momentum.eta()) < 2.5) # detector acceptance
                accept *= (abs(particle.momentum.pt()) > 1) # reconstructable pT
                if accept and abs(particle.pid) in [12, 14, 16,4900113,4900111,4900101,4900021]:
                    # add px and py in the corresponding place
                    missing_px += particle.momentum.px
                    missing_py += particle.momentum.py
                    
                        


                # decide which particles to keep
                accept = (particle.status == 1) # final state particle
                accept *= (abs(particle.momentum.eta()) < 2.5) # detector acceptance
                accept *= (abs(particle.momentum.pt()) > 1) # reconstructable pT
                accept *= (abs(particle.momentum.pt()) > 1) # reconstructable pT
                accept *= (abs(particle.pid) not in [12, 14, 16]) # not a neutrino
                if not accept:
                    continue



                particles.append( fastjet.PseudoJet(particle.momentum.px,particle.momentum.py,particle.momentum.pz,particle.momentum.e) )

                npart += 1


                if ((abs(particle.pid) in [11,13]) and (particle.momentum.pt()>27)):
                    leptons_27 +=1

                if ((abs(particle.pid) in [11,13]) and (particle.momentum.pt()>15)):
                    leptons_15 +=1


            # check lepton trigger
            if leptons_15>1:
                events_dilepton_trigger += 1

            
            missing_pt = np.sqrt(missing_px**2 + missing_py**2)

            if missing_pt > 200:
                events_missing_trigger += 1
            #    print("Plotting missing_pt event")
            #    plot_event(event,ievt,"missing_pt",sample)




            if leptons_27>0:
                events_single_lepton +=1


            # jet clustering check HT trigger
            cluster = fastjet.ClusterSequence(particles, jetdef)

            inc_jets = cluster.inclusive_jets()
            for jet in inc_jets:
                if jet.pt()>1000:
                    events_HT_trigger +=1
                    break



            if maxevents!=-1 and ievt >= maxevents-1 : break 


        nevents = ievt +1

        return events_missing_trigger,events_dilepton_trigger,events_HT_trigger, events_single_lepton,nevents



def analyse_sample(sample,num_samples,sample_index):
    print(f"Working on file {sample}, {sample_index}/{num_samples}")
    
    #print_progress_bar(sample_index + 1, num_samples, prefix='Progress:', suffix='Complete', length=50)

    f = open("triggers_output","a")

    events_missing_trigger,events_dilepton_trigger,events_HT_trigger,events_single_lepton,nevents = compare_triggers(sample)

    f.write(f"-------------------------------\n")
    f.write(f"----{sample}---\n")
    f.write(f"-------------------------------\n")
    f.write(f"nevents: {nevents}\n")
    f.write(f"Events missing pt trigger: {events_missing_trigger} -> {events_missing_trigger *xs(sample) / nevents}[nb]\n")
    f.write(f"Events dilepton trigger: {events_dilepton_trigger} -> {events_dilepton_trigger *xs(sample) / nevents}[nb]\n")
    f.write(f"Events HT trigger: {events_HT_trigger} -> {events_HT_trigger *xs(sample) / nevents}[nb]\n")
    f.write(f"Events single lepton trigger: {events_single_lepton} -> {events_single_lepton *xs(sample) / nevents}[nb]\n")
    f.flush() # to see the output in real time

    f.close()

def main():
    #samples = ["qcd","wjets","zjets","ttbar","higgs_portal_m=5_xio=1_xil=1_ctauMin","higgs_portal_m=10_xio=1_xil=1_ctauMin","higgs_portal_m=15_xio=1_xil=1_ctauMin"]
    #samples = ["higgs_portal_m=5_xio=1_xil=1_ctauMin","higgs_portal_m=10_xio=1_xil=1_ctauMin","higgs_portal_m=15_xio=1_xil=1_ctauMin"]
    input_path = "inputs"
    samples = []
    for hepmc_file in glob.glob(input_path+"/*.hepmc"):
        hepmc_file = hepmc_file.split("/")[-1] # get rid of the path
        samples.append(hepmc_file[:-6]) # get rid of the extention


    f = open("triggers_output","a")
    f.write(f"-------------------------------\n")
    f.write(f"---- Starting a run ---\n")
    f.write(f"-------------------------------\n")
    f.flush() # to see the output in real time

    f.close()


    
    num_samples = len(samples)
    args = [(sample,num_samples, i) for i,sample in enumerate(samples)]
    # Determine the number of available CPUs
    num_cpus = mp.cpu_count()
    
    # Create a pool of worker processes
    with mp.Pool(num_cpus) as pool:
        # Run the analysis function on each data file in parallel using starmap
        result_files = pool.starmap(analyse_sample, args)
    
    #for sample_index,sample in enumerate(samples):
    #    analyse_sample(sample,num_samples,sample_index)




if __name__ == "__main__":
    main()
