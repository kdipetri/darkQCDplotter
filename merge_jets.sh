rm outputs/jets/signal_jets.root
hadd outputs/jets/signal_jets.root outputs/jets/hists_higgs_portal_m=5_xio=1_xil=1_ctauMin_jets.root outputs/jets/hists_higgs_portal_m=10_xio=1_xil=1_ctauMin_jets.root outputs/jets/hists_higgs_portal_m=15_xio=1_xil=1_ctauMin_jets.root


rm outputs/jets/background_jets.root
hadd outputs/jets/background_jets.root outputs/jets/hists_qcd_jets.root outputs/jets/hists_ttbar_jets.root outputs/jets/hists_wjets_jets.root outputs/jets/hists_zjets_jets.root