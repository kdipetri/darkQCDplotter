#python3 jetFinder.py -s higgs_portal_m=5_xio=1_xil=1_ctauMin
#python3 jetFinder.py -s higgs_portal_m=10_xio=1_xil=1_ctauMin
#python3 jetFinder.py -s higgs_portal_m=15_xio=1_xil=1_ctauMin

#python3 jetFinder.py -s qcd
#python3 jetFinder.py -s ttbar
#python3 jetFinder.py -s wjets
#python3 jetFinder.py -s zjets


python3 jetLooper.py -s higgs_portal_m=5_xio=1_xil=1_ctauMin -d
python3 jetLooper.py -s higgs_portal_m=10_xio=1_xil=1_ctauMin -d
python3 jetLooper.py -s higgs_portal_m=15_xio=1_xil=1_ctauMin -d

#python3 jetLooper.py -s qcd
#python3 jetLooper.py -s ttbar
#python3 jetLooper.py -s wjets
#python3 jetLooper.py -s zjets