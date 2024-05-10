import argparse
from plothelper import *



fsignal = ROOT.TFile.Open("outputs/jets/signal_jets.root")
fbackground = ROOT.TFile.Open("outputs/jets/background_jets.root")


c = ROOT.TCanvas("c","",800,800)

def col(sample):
	if "qcd" in sample: return ROOT.kBlue+1
	if "ttbar" in sample: return ROOT.kRed+1
	if "wjets" in sample: return ROOT.kGray+1
	if "zjets" in sample: return ROOT.kOrange
	if "higgs_portal_m=5" in sample: return ROOT.kBlue+1
	if "higgs_portal_m=10" in sample: return ROOT.kRed+1
	if "higgs_portal_m=15" in sample: return ROOT.kGreen+2


def xs(sample):
	if "qcd" in sample: return 1e8
	if "wjets" in sample: return 1e5
	if "zjets" in sample: return 3e4
	if "ttbar" in sample: return 5e2
	if "higgs" in sample: return 1e2

def label(sample):
	if "qcd" in sample: return "qcd"
	if "ttbar" in sample: return "ttbar"
	if "wjets" in sample: return "wjets"
	if "zjets" in sample: return "zjets"
	if "higgs_portal_m=5" in sample: return "Higgs Portal m=5"
	if "higgs_portal_m=10" in sample: return "Higgs Portal m=10"
	if "higgs_portal_m=15" in sample: return "Higgs Portal m=15"

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
	c.SaveAs("plots/jets/stack_jets_"+dist+".png")


def compare1D(samples,dist,f,name):

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
	c.SaveAs("plots/jets/compare_jets_"+name+"_"+dist+".png")




def plot1D(sample,dist,f):

	c.cd()

	leg = legend(0.65,0.85,0.9,0.9)

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
	c.SaveAs("plots/jets/clean_jets_"+sample+"_"+dist+".png")








def plotSignalAndBackground(background,signal,dist,fsignal,fbackground,ymin=None,ymax=None,xmin=None,xmax=None):

	c.cd()

	leg = legend()

	hstack = ROOT.THStack(dist,"")

	for i,sample in enumerate(background): 
		try:
			hist = getHist(sample,dist,fbackground) 
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
	c.SaveAs("plots/jets/signal_and_background_jets_"+dist+".png")







signal = ["higgs_portal_m=5_xio=1_xil=1_ctauMin","higgs_portal_m=10_xio=1_xil=1_ctauMin","higgs_portal_m=15_xio=1_xil=1_ctauMin"]
background = ["ttbar","zjets","wjets","qcd"]

for sample in signal:
	plot1D(sample,"jet_displacement",fsignal)
	plot1D(sample,"jet_r_inv",fsignal)
	plot1D(sample,"jet_fraction_dark",fsignal)
	plot1D(sample,"jet_prompt_tracks",fsignal)
	plot1D(sample,"jet_displaced_tracks",fsignal)

for sample in background:
	plot1D(sample,"jet_prompt_tracks",fbackground)
	plot1D(sample,"jet_displaced_tracks",fbackground)
	
	
compare1D(signal,"jet_displacement",fsignal,"signal")
compare1D(signal,"jet_r_inv",fsignal,"signal")
compare1D(signal,"jet_fraction_dark",fsignal,"signal")
compare1D(signal,"jet_prompt_tracks",fsignal,"signal")
compare1D(signal,"jet_displaced_tracks",fsignal,"signal")

compare1D(background,"jet_prompt_tracks",fbackground,"background")
compare1D(background,"jet_displaced_tracks",fbackground,"background")

plotSignalAndBackground(background,signal,"jet_prompt_tracks",fsignal,fbackground,ymin=1e-2,ymax=1e8)
plotSignalAndBackground(background,signal,"jet_displaced_tracks",fsignal,fbackground,ymin=1e-2,ymax=1e8)


