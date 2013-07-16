#ifndef CutVarTreeProducer_h
#define CutVarTreeProducer_h

#include "SusyAnalysis/SusyNtuplizer/src/SusyEvent.h"

#include "TString.h"
#include "TTree.h"

#include <vector>
#include <set>

TString const CutVariables[] = {
  "met",
  "ht",
  "recoMet",
  "photon1Pt",
  "photon2Pt",
  "dPhiPhoton1Photon2",
  "dPhiPhoton1Met",
  "dPhiDiphotonMet",
  "dPtPhoton1Met",
  "dPtPhoton1Photon2",
  "mtPhoton1",
  "mtPhoton2",
  "electron1Pt",
  "electron2Pt",
  "dPhiPhoton1Electron1",
  "dPhiElectron1Electron2",
  "dPtElectron1Met",
  "mtElectron1",
  "mtElectron2",
  "muon1Pt",
  "muon2Pt",
  "dPhiPhoton1Muon1",
  "dPhiMuon1Muon2",
  "dPtMuon1Met",
  "mtMuon1",
  "dPhiElectron1Muon1",
  "nuPt",
  "dMetNuPt",
  "dRecoMetNuPt",
  "nuPz",
  "wP",
  "nJets",
  "nPhotons",
  "nElectrons",
  "nMuons"
};

class GenFilter {
 public:
  GenFilter(TString const&);
  ~GenFilter();
  TString toString() const;
  bool pass(std::vector<susy::Particle> const&, std::set<short>* = 0) const;
  enum Types {
    Atomic,
    LogicOR,
    LogicAND,
    Wrapping,
    nTypes
  };
 private:
  TString expr_;
  float ptThreshold_;
  float maxEta_;
  bool evalTo_;
  std::vector<unsigned> decayChain_;
  std::vector<GenFilter*> subfilters_;
  Types type_;
};

class CutVarTreeProducer {
 public:
  CutVarTreeProducer(TString const&, unsigned = 1);
  void produce(TTree&, TString const&);
  unsigned nRunning();
 private:
  void run() const;

  GenFilter const filter_;
  unsigned const nThreads_;
  std::map<TString, TObjArray> jobStack_;
};

#endif
