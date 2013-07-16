{
#include <stdexcept>

  // Give full path to libSusyEvent.so or have the directory in LD_LIBRARY_PATH
  // With the default setting, the library exists in SusyAnalysis/SusyNtuplizer/macro after issuing the make command in the macro directory.
  gSystem->Load("libSusyEvent.so");
  gSystem->Load("libSigmaScan.so");
  // cmsenv needs to be issued before using this tool
  gSystem->AddIncludePath("-I" + TString(gSystem->Getenv("CMSSW_BASE")) + "/src");
  gROOT->LoadMacro("SigmaScan.cc+");

  // Create the main object
  SigmaScan sigmaScan;

  /* COMMON BLOCK */

  // WORKSPACE:
  //// Directory where all the products of step 1 and 2 go into.
  //// This parameter should be unique to the signal and background considered.
  sigmaScan.addParam("WORKSPACE", "~/scratch0/gsq_W_mu_incl");

  // NCPU:
  //// Number of parallel threads to run in the tree prduction and scan. (optional - defaults to 1)
  sigmaScan.addParam("NCPU", "4");

  /* COMMON BLOCK */

  /* STEP1 BLOCK */

  // SIGNAL.LIST:
  //// The list (one file per line) of paths (URLs) to the signal sample files.
  sigmaScan.addParam("SIGNAL.Spectra_gsq_W.LIST", "Spectra_gsq_W.list");
  // SIGNAL.GRIDPARAMPATTERN:
  //// Perl regular expression (see documentations online) specifying which part of the file name of the signal sources should be taken as the unique identifier of the signal point.
  //// This identifier is handled as a whole string throughout the execution of the tool (the string will be stored in the end product tree).
  //// Therefore the format only matters in your macro that will run over the product of the scan.
  //// The pattern should have at least one matching sequence (...), and no wrapped sequences. Each matched sequence is concatenated with an underscore.
  //// Note that files whose names do not match the pattern will not be processed.
  //// Example: pattern="tree_naturalHiggsinoNLSPout_mst_([0-9]+)_M3_5025_mu_([0-9]+)[.]root" file="tree_naturalHiggsinoNLSPout_mst_100_M3_5025_mu_125.root" -> gridParams=100_125
  sigmaScan.addParam("SIGNAL.Spectra_gsq_W.GRIDPARAMPATTERN", "tree_([0-9]+_[0-9]+)_375[.]root");
  // SIGNAL.DECAY:
  //// Gen-level filter to apply on the signal samples. Only events matching the specific decay pattern will be used.
  //// Filter can be an empty string to accept all events.
  //// Particles should be specified by PDG IDs. Charge is ignored in the matching. Only for generic jets and b-jets, symbols j and b are used respectively.
  //// The first particle in each "word" is the final state particle. Users can give as many decay mothers as wanted.
  //// The example below selects events with >=1 photon from a neutralino and >=1 electron from either a W from a chargino or a Z from a neutralino.
  //// Photon kinematics is constrained to Pt>40GeV and |eta|<2.5. The default thresholds are Pt>10GeV and |eta|<3.
  //// Other examples: 22[40.,1.5]<25<1000025 && 22[25.,1.5]<25<1000025 && b<6<1000006
  //// (>=1 barrel photon with Pt>40GeV from a higgs from a higgsino, >=1 barrel photon with Pt>25GeV (guaranteed to not overlap with the leading) from a higgs from a higgsino, and a b-jet from a top from a stop)
  sigmaScan.addParam("SIGNAL.Spectra_gsq_W.DECAY", "22[40.,2.5]<1000022 && (13<24<1000024 || 13<23<1000022)");

  // BACKGROUND.XYZ.LIST:
  //// The list (one file per line) of paths (URLs) to the background sample files.
  //// Unlike the signal samples, background trees will be handled as a single TChain (except in the case of multi-thread processing, which is to be implemented).
  //// Background samples do not have to be declared beforehand; a directory will be automatically created at the first appearance of the sample in the configuration.
  sigmaScan.addParam("BACKGROUND.WGToLNuG.LIST", "WGToLNuG.list");
  // BACKGROUND.XYZ.DECAY:
  //// Gen-level filter to apply on the background samples. See above for the syntax description.
  //// The example here selects events with any non-hadronic photon and an electron from a W.
  sigmaScan.addParam("BACKGROUND.WGToLNuG.DECAY", "22 && 13<24");

  sigmaScan.addParam("BACKGROUND.ZGToLLG.LIST", "ZGToLLG.list");
  sigmaScan.addParam("BACKGROUND.ZGToLLG.DECAY", "22 && 13<23");

  /* STEP1 BLOCK */

  // Run the first step = produce ntuples for all signal and background samples.
  try{
    sigmaScan.produceTrees();
  }
  catch(exception& e){
    cout << e.what() << endl;
    return;
  }

  /* STEP2 BLOCK */

  // JOBNAME:
  //// Name of the scan to be run. This parameter will be used as the file name of the scan product.
  sigmaScan.addParam("SCAN.NAME", "met_ph30_mu30_mt100");
  // METHOD:
  //// Scan method. Only CutBased is implemented as of December 2012. Possible future additions are: BDT, Likelihood
  //// In the CutBased scan, the user will specify one or more variables to cut on to find the optimal set of cuts that maximizes the signal/sqrt(signal+bkg).
  //// Simple rectangular cuts will be used, and yields below and above each cut will be tested.
  sigmaScan.addParam("SCAN.METHOD", "CutBased");
  // INTEGRATEDLUMI:
  //// Integrated luminosity (in /pb) to normalize all samples to.
  sigmaScan.addParam("INTEGRATEDLUMI", "20000.");

  // SIGNAL.NAME:
  //// Name of the signal sample
  sigmaScan.addParam("SIGNAL.NAME", "Spectra_gsq_W");
  // SIGNAL.XSECFILE:
  //// A file that contains the cross section of the each signal point.
  //// Each line has to have the format "$GRIDPARAMS $XSEC $NEVENTS" where GRIDPARAMS is the same parameter explained above. XSEC is in pb.
  //// $GRIDPARAMS is the parameter generated from the file name (see SIGNAL.GRIDPARAMPATTERN above).
  sigmaScan.addParam("SIGNAL.XSECFILE", "Spectra_gsq_W_xsec.txt");
  // SIGNAL.EFFICIENCY:
  //// A weight to give to the signal samples uniformly. The value can come from e.g. object ID efficiencies.
  sigmaScan.addParam("SIGNAL.EFFICIENCY", "1.");

  // BACKGROUND.XYZ.XSEC:
  //// Cross section (in pb) of the sample XYZ.
  sigmaScan.addParam("BACKGROUND.WGToLNuG.XSEC", "461.6");
  // BACKGROUND.XYZ.NEVENTS:
  //// Total number of events in the sample.
  sigmaScan.addParam("BACKGROUND.WGToLNuG.NEVENTS", "5000000");
  // BACKGROUND.XYZ.FAKERATE"
  //// A weight to give to the sample.
  sigmaScan.addParam("BACKGROUND.WGToLNuG.FAKERATE", "1");

  sigmaScan.addParam("BACKGROUND.ZGToLLG.XSEC", "132.6");
  sigmaScan.addParam("BACKGROUND.ZGToLLG.NEVENTS", "6000000");
  sigmaScan.addParam("BACKGROUND.ZGToLLG.FAKERATE", "1");

  // SCAN.STEPS:
  //// CutBased spedcific configuration. Number of steps of the cut-value scan. (Each cut value needs to be given a range)
  sigmaScan.addParam("SCAN.STEPS", "10");

  // VARIABLES
  //// Syntax:
  //// variable[x:] -> Apply a fixed cut variable>x before starting the scan
  //// variable[:x] -> Apply a fixed cut variable<x before starting the scan
  //// variable[x:y] -> (CutBased specific) Scan the variable from x to y to find the optimal cut
  //// List of available variables is in CutVarTreeProducer.h. If you wish to add one (or more), edit the list and give the implementation in the CutVarTreeProducer::produce() function.
  sigmaScan.addVariable("met[50:200]");
  sigmaScan.addVariable("photon1Pt[30:]");
  sigmaScan.addVariable("muon1Pt[30:]");
  sigmaScan.addVariable("mtMuon1[100:]");

  /* STEP2 BLOCK */

  try{
    sigmaScan.scan();
  }
  catch(exception& e){
    cout << e.what() << endl;
    return;
  }
}
