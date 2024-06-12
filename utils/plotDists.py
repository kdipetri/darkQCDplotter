

from plothelper import *



f = ROOT.TFile.Open("outputs/all.root")
fout = ROOT.TFile.Open("outputs/plots.root","RECREATE")



c = ROOT.TCanvas("c","",800,800)

def col(sample):
	if "qcd" in sample: return ROOT.kGray
	if "wjets" in sample: return ROOT.kViolet+1
	if "zjets" in sample: return ROOT.kAzure+1
	if "ttbar" in sample: return ROOT.kOrange+1
	return ROOT.kRed+1

def xs(sample):
	if "qcd" in sample: return 1e8
	if "wjets" in sample: return 1e5
	if "zjets" in sample: return 3e4
	if "ttbar" in sample: return 5e2
	return 100

def label(sample):
	if "qcd" in sample: return "QCD"
	if "wjets" in sample: return "W+jets"
	if "zjets" in sample: return "Z+jets"
	if "ttbar" in sample: return "t#bar{t}"
	return "signal"

def legend(xmin=0.7,ymin=0.75,xmax=0.9,ymax=0.92):
	leg = ROOT.TLegend(xmin,ymin,xmax,ymax)
	leg.SetBorderSize(0)
	leg.SetTextFont(42)
	leg.SetTextSize(0.05)
	return leg

def getHist(sample,dist,sig=False):

	hist = f.Get(sample+"_"+dist)

	if not hist: 
		print(sample+"_"+dist, "not found")
		return 

	hist.SetLineColor(col(sample))
	hist.SetLineWidth(3)
	if not sig: 
		hist.SetFillColorAlpha(col(sample),0.7)
		hist.SetFillStyle(1001)	

	hist.SetMarkerColor(col(sample))
	hist.SetMarkerSize(0)	

	hist.GetXaxis().SetNdivisions(505)
	hist.GetYaxis().SetNdivisions(505)

	hist.Scale(xs(sample)/hist.Integral(0,-1))
	hist.SetDirectory(0)

	return hist

def stack1D(dist):

	c.cd()

	samples = ["ttbar","zjets","wjets","qcd"]

	leg = legend()

	hstack = ROOT.THStack(dist,"")

	for i,sample in enumerate(samples): 

		hist = getHist(sample,dist) 
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
	hstack.SetMinimum(1)

	signals = ["vector_m-1_ctau-10"]
	for i,signal in enumerate(signals):
		hist = getHist(signal,dist,sig=True) 
		if not hist: 
			print(sample+"_"+dist, "not found")
			continue

		hist.Draw("hsame")

		leg.AddEntry(hist,label(signal),"l")

	leg.Draw()
	c.SetLogy()
	c.SaveAs("plots/compare/stack_"+dist+".png")


def compare1D(dist):

	c.cd()

	leg = legend()

	hstack = ROOT.THStack(dist,"")

	samples = ["qcd","wjets","zjets","ttbar","vector_m-1_ctau-10"]
	

	for i,sample in enumerate(samples): 

		hist = getHist(sample,dist) 
		if not hist: 
			print(sample+"_"+dist, "not found")
			continue

		hist.Scale(1.0/hist.Integral(0,-1))
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
	hstack.SetMaximum(10)
	hstack.SetMinimum(0.001)
	leg.Draw()
	c.SetLogy()
	c.SaveAs("plots/compare/compare_"+dist+".png")

def plot1D(sample,dist):

	c.cd()

	leg = legend(0.75,0.85,0.9,0.9)

	hstack = ROOT.THStack(dist,"")

	hist = getHist(sample,dist) 
	if not hist: 
		print(sample+"_"+dist, "not found")
		return

	hist.Scale(1.0/hist.Integral(0,-1))

	leg.AddEntry(hist,label(sample),"l")

	hist.Draw("h")
	hstack.SetMaximum(10)
	hstack.SetMinimum(0.001)
	leg.Draw()
	c.SetLogy()
	c.SaveAs("plots/compare/clean_"+dist+".png")


plot1D("wjets","w_mass")
plot1D("wjets","w_pt")
plot1D("wjets","w_eta")
plot1D("zjets","z_mass")
plot1D("zjets","z_pt")
plot1D("zjets","z_eta")
plot1D("ttbar","t_mass")
plot1D("ttbar","t_pt")
plot1D("ttbar","t_eta")

plot1D("vector_m-1_ctau-10","h_pt")
plot1D("vector_m-1_ctau-10","h_eta")
plot1D("vector_m-1_ctau-10","h_mass")

plot1D("vector_m-1_ctau-10","dark_pt")
plot1D("vector_m-1_ctau-10","dark_eta")
plot1D("vector_m-1_ctau-10","dark_mass")
plot1D("vector_m-1_ctau-10","dark_decayrxy")
plot1D("vector_m-1_ctau-10","dark_prodrxy")
plot1D("vector_m-1_ctau-10","dark_status")

compare1D("njets")
compare1D("nleptons")
compare1D("nparticles")
compare1D("ncharged")
compare1D("ndisplaced")
compare1D("particle_eta")
compare1D("particle_pt")
compare1D("particle_pid")
compare1D("particle_prodrxy")
compare1D("charged_eta")
compare1D("charged_pt")
compare1D("charged_pid")
compare1D("charged_prodrxy")
compare1D("displaced_eta")
compare1D("displaced_pt")
compare1D("displaced_pid")
compare1D("displaced_prodrxy")
compare1D("lepton_eta")
compare1D("lepton_pt")
compare1D("jet_pt")
compare1D("jet_eta")


stack1D("njets")
stack1D("nleptons")
stack1D("nparticles")
stack1D("ncharged")
stack1D("ndisplaced")
stack1D("particle_eta")
stack1D("particle_pt")
stack1D("particle_pid")
stack1D("particle_prodrxy")
stack1D("charged_eta")
stack1D("charged_pt")
stack1D("charged_pid")
stack1D("charged_prodrxy")
stack1D("displaced_eta")
stack1D("displaced_pt")
stack1D("displaced_pid")
stack1D("displaced_prodrxy")
stack1D("lepton_eta")
stack1D("lepton_pt")
stack1D("jet_pt")
stack1D("jet_eta")
