#include "CutVarTreeProducer.h"

#include "SusyAnalysis/SusyNtuplizer/src/SusyEvent.h"

#include <iostream>
#include <stdexcept>
#include <map>
#include <limits>

#include "TChain.h"
#include "TPRegexp.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TArrayI.h"
#include "TThread.h"

GenFilter::GenFilter(TString const& _expr) :
  ptThreshold_(10.),
  maxEta_(3.),
  evalTo_(true),
  decayChain_(0),
  subfilters_(0),
  type_(nTypes)
{
  TString expr(_expr.Strip(TString::kBoth));

  if(expr.Length() == 0) return;

  // Example expression: (neutralino -> Z -> e OR chargino -> W -> e) AND neutralino -> gamma
  // (11[20.,3.0]<23<1000022 || 11[20.,3.0]<24<1000024) && 22[30.,2.5]<1000022
  TPRegexp chainPat("(![ ]*|)([0-9jb]+[ ]*(?:\\[[ ]*[0-9.]+[ ]*,[ ]*[0-9.]+[ ]*\\][ ]*|)(?:<[ ]*[0-9]+[ ]*)*)");
  //                 1      12                                                                              2
  TPRegexp openPat("(?:![ ]*|)\\(");
  //
  TPRegexp particlePat("([0-9jb]+)[ ]*(?:\\[[ ]*([0-9.]+)[ ]*,[ ]*([0-9.]+)[ ]*\\]|)");
  //                    1        1              2       2         3       3
  TPRegexp wrappingPat("^(?:![ ]*|)\\(.*\\)$");
  //
  unsigned pos(0);

  if(wrappingPat.MatchB(expr)){
    pos = expr.Index("(", pos) + 1;
    unsigned nOpen(1);
    unsigned nClose(0);
    while(pos < unsigned(expr.Length())){
      unsigned open(expr.Index("(", pos));
      unsigned close(expr.Index(")", pos));
      if(open < close){
        ++nOpen;
        pos = open + 1;
      }
      else if(close < open){
        ++nClose;
        pos = close + 1;
      }
      else pos = -1;
      if(nClose == nOpen) break;
    }

    if(pos == unsigned(-1)) throw std::invalid_argument(("Invalid filter expression " + _expr).Data());
    if(pos == unsigned(expr.Length())){
      type_ = Wrapping;
      unsigned open(expr.Index("("));
      evalTo_ = unsigned(expr.Index("!")) > open;
      subfilters_.push_back(new GenFilter(expr(open + 1, pos - open - 2)));

      return;
    }
  }

  if(chainPat.MatchB(expr)){
    TArrayI positions;
    chainPat.Match(expr, "", 0, 3, &positions);

    if(positions[1] - positions[0] == expr.Length()){
      type_ = Atomic;
      evalTo_ = positions[2] == positions[3];

      TString chain(expr(positions[4], positions[5]));
      TObjArray* chainElems(chain.Tokenize("<"));
      decayChain_.resize(chainElems->GetEntries());
      for(int iP(0); iP < chainElems->GetEntries(); ++iP){
        TObjArray* partProps(particlePat.MatchS(chainElems->At(iP)->GetName(), ""));
        TString pdgStr(partProps->At(1)->GetName());
        int pdgId(0);
        if(pdgStr == "j")
          pdgId = 100;
        else if(pdgStr == "b")
          pdgId = 500;
        else
          pdgId = pdgStr.Atoi();

        decayChain_[iP] = pdgId;

        if(iP == 0 && partProps->GetEntries() == 4){
          ptThreshold_ = TString(partProps->At(2)->GetName()).Atof();
          maxEta_ = TString(partProps->At(3)->GetName()).Atof();
        }

        delete partProps;
      }
      delete chainElems;

      return;
    }
  }

  pos = 0;
  bool noOperand(true);
  bool noSubexpr(false);
  while(pos < unsigned(expr.Length())){
    unsigned amp(expr.Index("&", pos));
    unsigned bar(expr.Index("|", pos));
    unsigned open(expr.Index(openPat, pos));
    unsigned dec(expr.Index(chainPat, pos));

    unsigned min(expr.Length());
    if(amp < min) min = amp;
    if(bar < min) min = bar;
    if(open < min) min = open;
    if(dec < min) min = dec;

    if(min == unsigned(expr.Length()) || min != pos) break;

    if(min == dec){
      if(noSubexpr) throw std::invalid_argument(("Invalid filter expression " + _expr).Data());

      TArrayI positions;
      chainPat.Match(expr, "", min, 3, &positions);
      subfilters_.push_back(new GenFilter(expr(positions[0], positions[1] - positions[0])));
      pos = positions[1];

      noSubexpr = true;
      noOperand = false;
    }
    else if(min == open){
      if(noSubexpr) throw std::invalid_argument(("Invalid filter expression " + _expr).Data());
      
      pos = expr.Index("(", pos) + 1;
      unsigned nOpen(1);
      unsigned nClose(0);
      while(pos < unsigned(expr.Length())){
        open = expr.Index("(", pos);
        unsigned close(expr.Index(")", pos));
        if(open < close){
          ++nOpen;
          pos = open + 1;
        }
        else if(close < open){
          ++nClose;
          pos = close + 1;
        }
        else pos = -1;
        if(nClose == nOpen) break;
      }
      if(pos == unsigned(-1)) throw std::invalid_argument(("Invalid filter expression " + _expr).Data());

      subfilters_.push_back(new GenFilter(expr(min, pos - min)));

      noSubexpr = true;
      noOperand = false;
    }
    else if(min == amp){
      if(noOperand) throw std::invalid_argument(("Invalid filter expression " + _expr).Data());

      if(expr(amp + 1) != '&') throw std::invalid_argument(("Invalid filter expression " + _expr).Data());
      if(type_ != nTypes && type_ != LogicAND)
        throw std::invalid_argument(("Invalid filter expression " + _expr).Data());
      type_ = LogicAND;
      pos = amp + 2;

      noSubexpr = false;
      noOperand = true;
    }
    else if(min == bar){
      if(noOperand) throw std::invalid_argument(("Invalid filter expression " + _expr).Data());

      if(expr(bar + 1) != '|') throw std::invalid_argument(("Invalid filter expression " + _expr).Data());
      if(type_ != nTypes && type_ != LogicOR)
        throw std::invalid_argument(("Invalid filter expression " + _expr).Data());
      type_ = LogicOR;
      pos = bar + 2;

      noSubexpr = false;
      noOperand = true;
    }

    while(expr(pos) == ' ' && pos < unsigned(expr.Length())) ++pos;
  }
}

GenFilter::~GenFilter()
{
  for(unsigned iF(0); iF < subfilters_.size(); ++iF)
    delete subfilters_[iF];
}

TString
GenFilter::toString() const
{
  TString result("");
  if(!evalTo_) result += "!";
  switch(type_){
  case Atomic:
    for(unsigned iD(0); iD < decayChain_.size(); ++iD){
      int pdgId(decayChain_[iD]);
      TString pdgStr;
      if(pdgId == 100)
        pdgStr = "j";
      else if(pdgId == 500)
        pdgStr = "b";
      else
        pdgStr = TString::Format("%d", pdgId);
        
      result += pdgStr;

      if(iD==0){
        result += "[pt";
        result += TString::Format("%.1f", ptThreshold_);
        result += ",eta";
        result += TString::Format("%.1f", maxEta_);
        result += "]";
      }
      if(iD != decayChain_.size() - 1)
        result += " < ";
    }
    break;
  case Wrapping:
    result += subfilters_[0]->toString();
    break;
  default:
    result += "(";
    for(unsigned iF(0); iF < subfilters_.size(); ++iF){
      result += subfilters_[iF]->toString();
      if(type_ == LogicAND && iF < subfilters_.size() - 1)
        result += " && ";
      else if(type_ == LogicOR && iF < subfilters_.size() - 1)
        result += " || ";
    }
    result += ")";
    break;
  }

  return result;
}

bool
GenFilter::pass(std::vector<susy::Particle> const& _genParticles, std::set<short>* _indexList/* = 0*/) const
{
  if(expr_ == "") return true;

  std::set<short>* indexList(0);
  if(_indexList) indexList = _indexList;
  else indexList = new std::set<short>;

  bool result(false);

  switch(type_){
  case Atomic:
    {
      unsigned nP(_genParticles.size());
      std::vector<bool> isHadron(_genParticles.size(), false);

      if(decayChain_[0] != 100 && decayChain_[0] != 500){
        for(unsigned iP(0); iP < nP; ++iP){
          susy::Particle const& particle(_genParticles[iP]);
          if(particle.status != 1) continue;
          if(particle.motherIndex == -1) continue;
          if((std::abs(particle.pdgId) / 100) % 10 != 0)
            isHadron[particle.motherIndex]  = true;
        }
      }

      for(unsigned iP(0); iP < nP; ++iP){
        short index(iP);
        susy::Particle const* particle(&_genParticles[index]);
        if(particle->status != 1) continue;
        if(decayChain_[0] == 100){
          if(std::abs(particle->pdgId) < 100) continue;
          while((std::abs(particle->pdgId) / 100) % 10 != 0 && particle->motherIndex != -1){
            index = particle->motherIndex;
            particle = &_genParticles[index];
          }
        }
        else if(decayChain_[0] == 500){
          while(std::abs(particle->pdgId) != 5 && particle->motherIndex != -1){
            index = particle->motherIndex;
            particle = &(_genParticles[index]);
          }

          if(particle->motherIndex == -1) continue;
        }
        else{
          if(std::abs(particle->pdgId) != decayChain_[0]) continue;
          if(particle->motherIndex != -1 && isHadron[particle->motherIndex]) continue;
        }

        if(indexList->find(index) != indexList->end()) continue;

        if(particle->momentum.Pt() < ptThreshold_) continue;
        if(std::abs(particle->momentum.Eta()) > maxEta_) continue;

        // radiating particles
        while(particle->motherIndex != -1 && _genParticles[particle->motherIndex].pdgId == particle->pdgId)
          particle = &_genParticles[particle->motherIndex];

        unsigned nD(decayChain_.size());

        unsigned iM(1);
        for(;iM < nD; ++iM){
          if(particle->motherIndex == -1) break;
          particle = &_genParticles[particle->motherIndex];

          if(decayChain_[iM] != std::abs(particle->pdgId)) break;
        }

        if(iM == nD){
          indexList->insert(index);
          result = true;
          break;
        }
      }
    }
    break;
  case LogicAND:
    {
      result = true;
      for(unsigned iF(0); iF < subfilters_.size(); ++iF){
        if(!subfilters_[iF]->pass(_genParticles, indexList)){
          result = false;
          break;
        }
      }
    }
    break;
  case LogicOR:
    {
      for(unsigned iF(0); iF < subfilters_.size(); ++iF){
        if(subfilters_[iF]->pass(_genParticles, indexList)){
          result = true;
          break;
        }
      }
    }
    break;
  case Wrapping:
    result = subfilters_[0]->pass(_genParticles, indexList);
    break;
  default:
    break;
  }

  if(!_indexList) delete indexList;

  return result ^ !evalTo_;
}




CutVarTreeProducer::CutVarTreeProducer(TString const& _expr, unsigned _nThreads/* = 1*/) :
  filter_(_expr),
  nThreads_(_nThreads),
  jobStack_()
{
  std::cout << filter_.toString() << std::endl;
}

void
CutVarTreeProducer::produce(TTree& _susyTree, TString const& _outputName)
{
  if(nThreads_ > 1){
    while(nRunning() >= nThreads_) TThread::Sleep(5);

    TThread::Lock();

    TObjArray& arr(jobStack_[_outputName]);
    arr.SetOwner();
    arr.Add(&_susyTree);

    TThread::UnLock();

    TThread* thread(new TThread(_outputName, TThread::VoidFunc_t(&CutVarTreeProducer::run), reinterpret_cast<void*>(this)));
    thread->Run();
  }
  else{
    TObjArray& arr(jobStack_[_outputName]);
    arr.SetOwner();
    arr.Add(&_susyTree);

    run();

    TFile* file(0);
    if(!_susyTree.InheritsFrom(TChain::Class()))
      file = _susyTree.GetCurrentFile();

    jobStack_.erase(_outputName);
    delete file;
  }
}

unsigned
CutVarTreeProducer::nRunning()
{
  if(nThreads_ <= 1) return 0;

  unsigned result(0);

  std::map<TString, TObjArray>::const_iterator jItr(jobStack_.begin());
  while(jItr != jobStack_.end()){
    TString name(jItr->first);
    TThread* thread(TThread::GetThread(name));
    if(!thread || thread->GetState() == TThread::kTerminatedState || thread->GetState() == TThread::kFinishedState || thread->GetState() == TThread::kCanceledState){
      TThread::Lock();

      TThread::Delete(thread);

      TTree* tree(static_cast<TTree*>(jItr->second.At(0)));
      TFile* file(0);
      if(!tree->InheritsFrom(TChain::Class()))
        file = tree->GetCurrentFile();

      jobStack_.erase(name); // erases the ttree / tchain too (TObjArray is SetOwned)
      jItr = jobStack_.upper_bound(name);

      delete file;

      TThread::UnLock();
    }
    else{
      ++jItr;
      ++result;
    }
  }

  return result;
}

void
CutVarTreeProducer::run() const
{
  if(jobStack_.size() == 0) return;

  std::map<TString, int> ints;
  std::map<TString, float> floats;

  susy::Event* event(new susy::Event());

  TString outputName;
  if(nThreads_ > 1) outputName = TThread::Self()->GetName();
  else outputName = jobStack_.begin()->first;

  std::map<TString, TObjArray>::const_iterator jItr(jobStack_.find(outputName));

  if(jItr == jobStack_.end()){
    std::cerr << "No job info for " << outputName << " found in the stack" << std::endl;
    return;
  }

  TObjArray const& arr(jItr->second);
  TTree* susyTree(static_cast<TTree*>(arr.At(0)));

  TThread::Lock();

  susyTree->GetCurrentFile()->cd();

  susyTree->SetBranchAddress("susyEvent", &event);

  unsigned nV(sizeof(CutVariables) / sizeof(TString)); 

  TFile outputFile(outputName, "recreate");
  outputFile.cd();
  TTree* outTree(new TTree("eventTree", "simplified events"));

  for(unsigned iV(0); iV < nV; ++iV){
    TString const& var(CutVariables[iV]);
    if(var(0) == 'n')
      outTree->Branch(var, &(ints[var]), var + "/I");
    else
      outTree->Branch(var, &(floats[var]), var + "/F");
  }

  // not unlocking on purpose

  std::map<TString, int>::iterator intsEnd(ints.end());
  std::map<TString, float>::iterator floatsEnd(floats.end());

  long iEntry(0);

  while(susyTree->GetEntry(iEntry++) != 0){

    // locks at the end of the loop block

    std::vector<susy::Particle>& genParticles(event->genParticles);

    bool passGenFilter(filter_.pass(genParticles));

    if(!passGenFilter) continue;

    TThread::UnLock();

    for(std::map<TString, int>::iterator intsItr(ints.begin()); intsItr != intsEnd; ++intsItr) intsItr->second = 0;
    for(std::map<TString, float>::iterator floatsItr(floats.begin()); floatsItr != floatsEnd; ++floatsItr) floatsItr->second = 0;

    unsigned nP(genParticles.size());
    std::vector<bool> isHadron(genParticles.size(), false);
    for(unsigned iP(0); iP < nP; ++iP){
      susy::Particle& particle(genParticles[iP]);
      if(particle.status != 1) continue;
      if(particle.motherIndex == -1) continue;
      if((std::abs(particle.pdgId) / 100) % 10 != 0)
        isHadron[particle.motherIndex]  = true;
    }

    susy::Particle* electrons[] = {0, 0};
    susy::Particle* muons[] = {0, 0};
    susy::Particle* photons[] = {0, 0};
    TVector2 metV(0., 0.);

    for(unsigned iP(0); iP < nP; ++iP){
      susy::Particle& particle(genParticles[iP]);
      if(particle.status != 1) continue;
      if(particle.motherIndex == -1) continue;
      if(isHadron[particle.motherIndex]) continue;

      int pdgId(std::abs(particle.pdgId));
      float pt(particle.momentum.Pt());

      switch(pdgId){
      case 11:
        if(pt > floats["electron1Pt"]){
          floats["electron2Pt"] = floats["electron1Pt"];
          floats["electron1Pt"] = pt;
          electrons[0] = &particle;
        }
        else if(pt > floats["electron2Pt"]){
          floats["electron2Pt"] = pt;
          electrons[1] = &particle;
        }

        ints["nElectrons"] += 1;
        break;
      case 13:
        if(pt > floats["muon1Pt"]){
          floats["muon2Pt"] = floats["muon1Pt"];
          floats["muon1Pt"] = pt;
          muons[0] = &particle;
        }
        else if(pt > floats["muon2Pt"]){
          floats["muon2Pt"] = pt;
          muons[1] = &particle;
        }

        ints["nMuons"] += 1;
        break;
      case 22:
        if(pt > floats["photon1Pt"]){
          floats["photon2Pt"] = floats["photon1Pt"];
          floats["photon1Pt"] = pt;
          photons[0] = &particle;
        }
        else if(pt > floats["photon2Pt"]){
          floats["photon2Pt"] = pt;
          photons[1] = &particle;
        }

        ints["nPhotons"] += 1;
        break;
      case 12:
      case 14:
      case 16:
      case 1000022:
      case 1000039:
        metV += TVector2(particle.momentum.X(), particle.momentum.Y());
        break;
      default:
        break;
      }
    }

    floats["met"] = metV.Mod();

    TVector2 photon1Pt;
    TVector2 diphotonPt;
    TVector2 electron1Pt;
    TVector2 muon1Pt;

    if(photons[0]) photon1Pt.Set(photons[0]->momentum.X(), photons[0]->momentum.Y());
    if(photons[0] && photons[1]) diphotonPt.Set(photons[0]->momentum.X() + photons[1]->momentum.X(), photons[0]->momentum.Y() + photons[1]->momentum.Y());
    if(electrons[0]) electron1Pt.Set(electrons[0]->momentum.X(), electrons[0]->momentum.Y());
    if(muons[0]) muon1Pt.Set(muons[0]->momentum.X(), muons[0]->momentum.Y());

    if(photons[0]){
      floats["dPhiPhoton1Met"] = photon1Pt.DeltaPhi(metV);
      floats["dPtPhoton1Met"] = photons[0]->momentum.Pt() - metV.Mod();
      if(photons[1]){
        floats["dPhiPhoton1Photon2"] = photons[0]->momentum.DeltaPhi(photons[1]->momentum);
        floats["dPhiDiphotonMet"] = diphotonPt.DeltaPhi(metV);
        floats["dPtPhoton1Photon2"] = photon1Pt.Mod() - photons[1]->momentum.Pt();
      }
    }

    if(photons[0] && electrons[0]) floats["dPhiPhoton1Electron1"] = photons[0]->momentum.DeltaPhi(electrons[0]->momentum);
    if(photons[0] && muons[0]) floats["dPhiPhoton1Muon1"] = photons[0]->momentum.DeltaPhi(muons[0]->momentum);

    if(electrons[0]){
      floats["mtElectron1"] = std::sqrt((floats["met"] + electron1Pt.Mod()) * (floats["met"] + electron1Pt.Mod()) - (electron1Pt + metV).Mod2());
      floats["dPtElectron1Met"] = electron1Pt.Mod() - metV.Mod();
      if(electrons[1]) floats["dPhiElectron1Electron2"] = electrons[0]->momentum.DeltaPhi(electrons[1]->momentum);
    }

    if(muons[0]){
      floats["mtMuon1"] = std::sqrt((floats["met"] + muon1Pt.Mod()) * (floats["met"] + muon1Pt.Mod()) - (muon1Pt + metV).Mod2());
      floats["dPtMuon1Met"] = muon1Pt.Mod() - metV.Mod();
      if(muons[1]) floats["dPhiMuon1Muon2"] = muons[0]->momentum.DeltaPhi(muons[1]->momentum);
    }

    if(electrons[0] && muons[0]) floats["dPhiElectron1Muon1"] = electrons[0]->momentum.DeltaPhi(muons[0]->momentum);

    TVector2* leptonPt(0);
    double leptonPz(0.);
    double leptonP(0.);
    if(electron1Pt.Mod() > muon1Pt.Mod()){
      leptonPt = &electron1Pt;
      leptonPz = electrons[0]->momentum.Z();
      leptonP = electrons[0]->momentum.P();
    }
    else{
      leptonPt = &muon1Pt;
      leptonPz = muons[0]->momentum.Z();
      leptonP = muons[0]->momentum.P();
    }

    TVector2 nuPtV(photon1Pt + *leptonPt);
    nuPtV *= -1.;
    floats["nuPt"] = nuPtV.Mod();
    floats["dMetNuPt"] = (metV - nuPtV).Mod();

    double A(80.4 * 80.4 / 2. + *leptonPt * metV);
    double B(A * A - leptonPt->Mod2() * metV.Mod2());
    if(B > 0.){
      if(A * leptonPz > 0.) floats["nuPz"] = (A * leptonPz - leptonP * std::sqrt(B)) / leptonPt->Mod2();
      else floats["nuPz"] = (A * leptonPz + leptonP * std::sqrt(B)) / leptonPt->Mod2();

      TVector3 wp(leptonPt->X() + metV.X(), leptonPt->Y() + metV.Y(), leptonPz + floats["nuPz"]);
      floats["wP"] = wp.Mag();
    }
    else{
      floats["nuPz"] = std::numeric_limits<float>::max();
      floats["wP"] = std::numeric_limits<float>::max();
    }

    std::vector<TLorentzVector*> jets(0);
    for(unsigned iP(0); iP < nP; ++iP){
      susy::Particle& particle(genParticles[iP]);
      if(particle.status != 1) continue;

      if(electrons[0]){
        if(&particle == electrons[0]) continue;
        if(electrons[0]->momentum.DeltaR(particle.momentum) < 0.4) continue;
        if(electrons[1]){
          if(&particle == electrons[1]) continue;
          if(electrons[1]->momentum.DeltaR(particle.momentum) < 0.4) continue;
        }
      }
      if(muons[0]){
        if(&particle == muons[0]) continue;
        if(muons[0]->momentum.DeltaR(particle.momentum) < 0.4) continue;
        if(muons[1]){
          if(&particle == muons[1]) continue;
          if(muons[1]->momentum.DeltaR(particle.momentum) < 0.4) continue;
        }
      }
      if(photons[0]){
        if(&particle == photons[0]) continue;
        if(photons[0]->momentum.DeltaR(particle.momentum) < 0.4) continue;
        if(photons[1]){
          if(&particle == photons[1]) continue;
          if(photons[1]->momentum.DeltaR(particle.momentum) < 0.4) continue;
        }
      }
    }

    susy::PFJetCollection& pfJets(event->pfJets["ak5"]);
    unsigned nJ(pfJets.size());
    for(unsigned iJ(0); iJ < nJ; ++iJ){
      susy::PFJet& jet(pfJets[iJ]);
      if(jet.momentum.Pt() < 30.) continue;
      if(std::abs(jet.momentum.Eta()) > 3.) continue;
      if(photons[0]){
        if(photons[0]->momentum.DeltaR(jet.momentum) < 0.5) continue;
        if(photons[1] && photons[1]->momentum.DeltaR(jet.momentum) < 0.5) continue;
      }
      if(electrons[0]){
        if(electrons[0]->momentum.DeltaR(jet.momentum) < 0.5) continue;
        if(electrons[1] && electrons[1]->momentum.DeltaR(jet.momentum) < 0.5) continue;
      }
      if(muons[0]){
        if(muons[0]->momentum.DeltaR(jet.momentum) < 0.5) continue;
        if(muons[1] && muons[1]->momentum.DeltaR(jet.momentum) < 0.5) continue;
      }

      floats["ht"] += jet.momentum.Pt();
      ints["nJets"] += 1;
    }

    TVector2& recoMetV(event->metMap["pfMet"].mEt);
    floats["recoMet"] = recoMetV.Mod();
    floats["dRecoMetNuPt"] = (recoMetV - nuPtV).Mod();

    TThread::Lock();

    outTree->Fill();

    // not unlocking on purpose
  }

  // thread still locked

  susyTree->ResetBranchAddresses();
  delete event;

  outputFile.cd();
  outputFile.Write();

  TThread::UnLock();
}

