

from plothelper import *



f = ROOT.TFile.Open("outputs/all.root")
fout = ROOT.TFile.Open("outputs/plots.root","RECREATE")



c = ROOT.TCanvas("c","",800,800)

def col(sample):
	if "qcd" in sample: return ROOT.kGray
	if "wjets" in sample: return ROOT.kBlue+1
	if "zjets" in sample: return ROOT.kRed+1
	if "ttbar" in sample: return ROOT.kOrange+1

def xs(sample):
	if "qcd" in sample: return 1e8
	if "wjets" in sample: return 1e5
	if "zjets" in sample: return 3e4
	if "ttbar" in sample: return 5e2

def label(sample):
	if "qcd" in sample: return "QCD"
	if "wjets" in sample: return "W+jets"
	if "zjets" in sample: return "Z+jets"
	if "ttbar" in sample: return "t#bar{t}"

def legend(xmin=0.7,ymin=0.75,xmax=0.9,ymax=0.92):
	leg = ROOT.TLegend(xmin,ymin,xmax,ymax)
	leg.SetBorderSize(0)
	leg.SetTextFont(42)
	leg.SetTextSize(0.05)
	return leg

def getHist(sample,dist):

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
	hstack.SetMinimum(10)
	leg.Draw()
	c.SetLogy()
	c.SaveAs("plots/stack_"+dist+".png")


def compare1D(dist):

	c.cd()

	leg = legend()

	hstack = ROOT.THStack(dist,"")

	samples = ["qcd","wjets","zjets","ttbar",]
	

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
	c.SaveAs("plots/compare_"+dist+".png")

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
	c.SaveAs("plots/clean_"+dist+".png")


plot1D("wjets","w_mass")
plot1D("wjets","w_pt")
plot1D("wjets","w_eta")
plot1D("zjets","z_mass")
plot1D("zjets","z_pt")
plot1D("zjets","z_eta")
plot1D("ttbar","t_mass")
plot1D("ttbar","t_pt")
plot1D("ttbar","t_eta")

compare1D("njets")
compare1D("nleptons")
compare1D("nparticles")
compare1D("particle_eta")
compare1D("particle_pt")
compare1D("particle_pid")
compare1D("lepton_eta")
compare1D("lepton_pt")
compare1D("jet_pt")
compare1D("jet_eta")


stack1D("njets")
stack1D("nleptons")
stack1D("nparticles")
stack1D("particle_eta")
stack1D("particle_pt")
stack1D("particle_pid")
stack1D("lepton_eta")
stack1D("lepton_pt")
stack1D("jet_pt")
stack1D("jet_eta")
