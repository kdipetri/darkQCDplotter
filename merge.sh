rm outputs/background.root
rm outputs/signal.root
hadd outputs/background.root outputs/qcd.root outputs/wjets.root outputs/zjets.root outputs/ttbar.root 
hadd outputs/signal.root outputs/higgs_portal_m=5_xio=1_xil=1_ctauMin.root outputs/higgs_portal_m=10_xio=1_xil=1_ctauMin.root outputs/higgs_portal_m=15_xio=1_xil=1_ctauMin.root 