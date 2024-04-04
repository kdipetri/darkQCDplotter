import argparse
from plothelper import *


parser = argparse.ArgumentParser(usage=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s", "--dataName", help="Data name: background, signal, all", default="background")

ops = parser.parse_args()
data = ops.dataName

fbackgorund = ROOT.TFile.Open("outputs/background.root")
fsignal = ROOT.TFile.Open("outputs/signal.root")


c = ROOT.TCanvas("c","",800,800)

def col(sample):
	if "qcd" in sample: return ROOT.kGray
	if "wjets" in sample: return ROOT.kBlue+1
	if "zjets" in sample: return ROOT.kRed+1
	if "ttbar" in sample: return ROOT.kOrange+1
	if "higgs_portal_m=5" in sample: return ROOT.kMagenta+2
	if "higgs_portal_m=10" in sample: return ROOT.kBlack
	if "higgs_portal_m=15" in sample: return ROOT.kGreen+2

def xs(sample):
	if "qcd" in sample: return 1e8
	if "wjets" in sample: return 1e5
	if "zjets" in sample: return 3e4
	if "ttbar" in sample: return 5e2
	if "higgs" in sample: return 1e2

def label(sample):
	if "qcd" in sample: return "QCD"
	if "wjets" in sample: return "W+jets"
	if "zjets" in sample: return "Z+jets"
	if "ttbar" in sample: return "t #bar{t}"
	if "higgs_portal_m=5" in sample: return "Higgs Portal m=5, 100pb"
	if "higgs_portal_m=10" in sample: return "Higgs Portal m=10, 100pb"
	if "higgs_portal_m=15" in sample: return "Higgs Portal m=15, 100pb"

def legend(xmin=0.55,ymin=0.75,xmax=0.9,ymax=0.92):
	leg = ROOT.TLegend(xmin,ymin,xmax,ymax)
	leg.SetBorderSize(0)
	leg.SetTextFont(20)
	leg.SetTextSize(0.03)
	return leg

def getHist(sample,dist,f):

	hist = f.Get(sample+"_"+dist)

	hist.SetLineColor(col(sample))
	hist.SetLineWidth(3)

	hist.SetFillColorAlpha(col(sample),0.7)
	hist.SetFillStyle(1001)	

	hist.SetMarkerColor(col(sample))
	hist.SetMarkerSize(0)	

	hist.GetXaxis().SetNdivisions(505)
	hist.GetYaxis().SetNdivisions(505)

	hist.Scale(xs(sample)/hist.Integral(0,-1))
	hist.SetDirectory(0)

	
	
	return hist

def stack1D(samples,dist,f):

	c.cd()

	leg = legend()

	hstack = ROOT.THStack(dist,"")

	for i,sample in enumerate(samples): 
		try:
			hist = getHist(sample,dist,f) 
		except:
			print("Error geting hist")
			continue

		if not hist: 
			print(sample+"_"+dist, "not found")
			continue

		hstack.Add(hist)

		leg.AddEntry(hist,label(sample),"f")

		xtitle = hist.GetXaxis().GetTitle()
		ytitle = hist.GetYaxis().GetTitle()


	hstack.Draw("h")#"h" if i==0 else "hsame")
	hstack.GetXaxis().SetTitle(xtitle)
	hstack.GetYaxis().SetTitle(ytitle)
	hstack.GetXaxis().SetNdivisions(505)
	hstack.GetYaxis().SetNdivisions(505)
	#hstack.SetMinimum(.0001)

	
	leg.Draw()
	c.SetLogy()
	c.SaveAs("plots/stack_"+data+"_"+dist+".png")


def compare1D(samples,dist,f):

	c.cd()

	leg = legend()

	hstack = ROOT.THStack(dist,"")	

	for i,sample in enumerate(samples): 
		try:
			hist = getHist(sample,dist,f) 
		except:
			print(f"Error geting hist for {sample}")
			continue


		if not hist: 
			print(sample+"_"+dist, "not found")
			continue

		#hist.Scale(1.0/hist.Integral(0,-1))
		hist.SetFillStyle(0)
		hstack.Add(hist)

		leg.AddEntry(hist,label(sample),"l")

		xtitle = hist.GetXaxis().GetTitle()
		ytitle = hist.GetYaxis().GetTitle()


	hstack.Draw("h nostack")#"h" if i==0 else "hsame")
	hstack.GetXaxis().SetTitle(xtitle)
	hstack.GetYaxis().SetTitle(ytitle)
	hstack.GetXaxis().SetNdivisions(505)
	hstack.GetYaxis().SetNdivisions(505)
	#hstack.SetMaximum(ymin)
	#hstack.SetMinimum(ymax)
	leg.Draw()
	c.SetLogy()
	c.SaveAs("plots/compare_"+data+"_"+dist+".png")

def plot1D(sample,dist,f):

	c.cd()

	leg = legend(0.75,0.85,0.9,0.9)

	hstack = ROOT.THStack(dist,"")

	hist = getHist(sample,dist,f) 
	if not hist: 
		print(sample+"_"+dist, "not found")
		return

	#hist.Scale(1.0/hist.Integral(0,-1))

	leg.AddEntry(hist,label(sample),"l")

	hist.Draw("h")
	#hstack.SetMaximum(10)
	#hstack.SetMinimum(0.001)
	leg.Draw()
	c.SetLogy()
	c.SaveAs("plots/clean_"+dist+".png")



def plotSignalAndBackground(background,signal,dist,ymin=None,ymax=None,xmin=None,xmax=None):

	c.cd()

	leg = legend()

	hstack = ROOT.THStack(dist,"")

	for i,sample in enumerate(background): 
		try:
			hist = getHist(sample,dist,fbackgorund) 
		except:
			print("Error geting hist")
			continue

		if not hist: 
			print(sample+"_"+dist, "not found")
			continue

		hstack.Add(hist)

		leg.AddEntry(hist,label(sample),"f")

		xtitle = hist.GetXaxis().GetTitle()
		ytitle = hist.GetYaxis().GetTitle()



	hstack.Draw("h")#"h" if i==0 else "hsame")


	for sample in signal:
		try:
			hist = getHist(sample,dist,fsignal) 
		except:
			print("Error geting hist")
			continue

		if not hist: 
			print(sample+"_"+dist, "not found")
			continue

		hist.Draw("same")
		leg.AddEntry(hist,label(sample),"f")


	hstack.GetXaxis().SetTitle(xtitle)
	hstack.GetYaxis().SetTitle(ytitle)
	hstack.GetXaxis().SetNdivisions(505)
	hstack.GetYaxis().SetNdivisions(505)

	if ymin:
		hstack.SetMinimum(ymin)
	if ymax:
		hstack.SetMaximum(ymax)

	leg.Draw()
	c.SetLogy()
	c.SaveAs("plots/signal_and_background_"+dist+".png")




def runPlots(samples,f):
	compare1D(samples,"njets",f)
	compare1D(samples,"nleptons",f)
	compare1D(samples,"nparticles",f)
	compare1D(samples,"particle_eta",f)
	compare1D(samples,"particle_pt",f)
	compare1D(samples,"particle_pid",f)
	compare1D(samples,"lepton_eta",f)
	compare1D(samples,"lepton_pt",f)
	compare1D(samples,"jet_pt",f)
	compare1D(samples,"jet_eta",f)


	stack1D(samples,"njets",f)
	stack1D(samples,"nleptons",f)
	stack1D(samples,"nparticles",f)
	stack1D(samples,"particle_eta",f)
	stack1D(samples,"particle_pt",f)
	stack1D(samples,"particle_pid",f)
	stack1D(samples,"lepton_eta",f)
	stack1D(samples,"lepton_pt",f)
	stack1D(samples,"jet_pt",f)
	stack1D(samples,"jet_eta",f)





background = ["ttbar","zjets","wjets","qcd"]
signal = ["higgs_portal_m=5_xio=1_xil=1_ctauMin","higgs_portal_m=10_xio=1_xil=1_ctauMin","higgs_portal_m=15_xio=1_xil=1_ctauMin"]

plot1D("wjets","w_mass",fbackgorund)
plot1D("wjets","w_pt",fbackgorund)
plot1D("wjets","w_eta",fbackgorund)
plot1D("zjets","z_mass",fbackgorund)
plot1D("zjets","z_pt",fbackgorund)
plot1D("zjets","z_eta",fbackgorund)
plot1D("ttbar","t_mass",fbackgorund)
plot1D("ttbar","t_pt",fbackgorund)
plot1D("ttbar","t_eta",fbackgorund)

runPlots(background,fbackgorund)
runPlots(signal,fsignal)

plotSignalAndBackground(background,signal,"njets",ymin=1e-2)
plotSignalAndBackground(background,signal,"nleptons",ymin=1e-2)
plotSignalAndBackground(background,signal,"nparticles",ymin=1e-2)
plotSignalAndBackground(background,signal,"particle_eta",ymin=1e-2)
plotSignalAndBackground(background,signal,"particle_pt",ymin=1e-2)
plotSignalAndBackground(background,signal,"particle_pid",ymin=1e-2)
plotSignalAndBackground(background,signal,"lepton_eta",ymin=1e-2,ymax=1e5)
plotSignalAndBackground(background,signal,"lepton_pt",ymin=1e-2)
plotSignalAndBackground(background,signal,"jet_eta",ymin=1e-2,ymax=1e8)
plotSignalAndBackground(background,signal,"jet_pt",ymin=1e-3,ymax=1e9)


