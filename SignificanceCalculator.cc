#include "SignificanceCalculator.h"

#include <stdexcept>
#include <vector>
#include <set>
#include <cmath>
#include <iostream>

#include "TObjArray.h"
#include "TFile.h"
#include "TPRegexp.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TObjString.h"
#include "TThread.h"
#include "TChain.h"

SignificanceCalculator::SignificanceCalculator(TString const& _outputName, std::vector<TString> const& _variables, unsigned _nThreads) :
  backgrounds_(),
  gridParams_(new TString()),
  significance_(0.),
  signalEfficiency_(0.),
  variables_(_variables),
  fixedCuts_(0),
  distributions_(_variables.size(), 0),
  defaultRanges_(_variables.size(), std::pair<double, double>(0., 0.)),
  output_(0),
  nThreads_(_nThreads),
  jobStack_(),
  initialized_(false),
  verbose_(false)
{
  if(variables_.size() == 0)
    throw std::invalid_argument("No variables given");

  TFile* outputFile(new TFile(_outputName, "recreate"));
  outputFile->cd();
  output_ = new TTree("sigmaTree", "Signal significance");
  output_->Branch("gridParams", "TString", &gridParams_);
  output_->Branch("significance", &significance_, "significance/F");
  output_->Branch("signalEfficiency", &signalEfficiency_, "signalEfficiency/F");

  for(unsigned iV(0); iV < variables_.size(); ++iV){
    TString const& varName(variables_[iV]);
    distributions_[iV] = new TCanvas(varName, varName);
    output_->Branch("dist_" + varName, "TCanvas", &(distributions_[iV]));
  }
}

SignificanceCalculator::~SignificanceCalculator()
{
  TFile* outputFile(output_->GetCurrentFile());
  outputFile->cd();
  output_->Write();
  delete outputFile;

  delete gridParams_;
  for(unsigned iV(0); iV < variables_.size(); ++iV)
    delete distributions_[iV];
}

unsigned
SignificanceCalculator::nRunning()
{
  if(nThreads_ <= 1) return 0;

  unsigned result(0);

  TThread::Lock();

  std::map<TString, std::pair<TTree*, double> >::iterator jItr(jobStack_.begin());
  while(jItr != jobStack_.end()){
    TString gridParams(jItr->first);
    TThread* thread(TThread::GetThread(gridParams));
    if(!thread || thread->GetState() == TThread::kTerminatedState || thread->GetState() == TThread::kFinishedState || thread->GetState() == TThread::kCanceledState){
      TThread::Delete(thread);

      TTree* tree(jItr->second.first);
      if(tree->InheritsFrom(TChain::Class()))
        delete tree;
      else
        delete tree->GetCurrentFile();

      jobStack_.erase(jItr);
      jItr = jobStack_.upper_bound(gridParams);
    }
    else{
      ++jItr;
      ++result;
    }
  }
    
  TThread::UnLock();

  return result;
}



CutBasedCalculator::CutBasedCalculator(TString const& _outputName, std::vector<TString> const& _variables, std::vector<TString> const& _varRanges, unsigned _nThreads/* = 1*/) :
  SignificanceCalculator(_outputName, _variables, _nThreads),
  bounds_(0),
  floatingCuts_(0),
  cutValsI_(0),
  cutValsF_(0),
  upperCut_(0),
  signalYield_(0.),
  bkgYield_(0.),
  bkgYields_()
{
  if(variables_.size() != _varRanges.size())
    throw std::invalid_argument("Number of ranges does not match the number of variables");

  std::set<TString> intVars;

  TPRegexp rangePat("\\[([0-9.-]*)(|[:]([0-9.-]*))\\]");
  for(unsigned iV(0); iV < variables_.size(); ++iV){
    TString const& varName(variables_[iV]);
    bool isInt = varName(0) == 'n';
    if(isInt) intVars.insert(varName);

    TObjArray* matches(rangePat.MatchS(_varRanges[iV]));
    if(matches->GetEntries() == 0){
      delete matches;
      throw std::invalid_argument(("Invalid variable configuration " + _varRanges[iV]).Data());
    }

    TString low(matches->At(1)->GetName());
    TString colon(matches->At(2)->GetName());
    TString high(matches->GetEntries() > 3 ? matches->At(3)->GetName() : "");

    delete matches;

    if(low == "" && high == "")
      throw std::invalid_argument(("Invalid variable configuration " + _varRanges[iV]).Data());

    if(colon == ""){
      if(!isInt)
        throw std::invalid_argument(("Invalid variable configuration " + _varRanges[iV]).Data());
      fixedCuts_.push_back(TCut(varName + "==" + low));
    }
    else if(low == ""){
      fixedCuts_.push_back(TCut(varName + "<" + high));
    }
    else if(high == ""){
      fixedCuts_.push_back(TCut(varName + ">" + low));
    }
    else{
      floatingCuts_.push_back(varName + "?");
      bounds_.push_back(std::pair<float, float>(low.Atof(), high.Atof()));
    }
  }

  upperCut_.resize(floatingCuts_.size());
  cutValsI_.resize(intVars.size());
  cutValsF_.resize(floatingCuts_.size() - intVars.size());

  unsigned iI(0);
  unsigned iF(0);
  for(unsigned iB(0); iB < floatingCuts_.size(); ++iB){
    TString varName(floatingCuts_[iB]);
    varName.Remove(varName.Length() - 1);

    TString leaf;
    void* addr(0);
    if(intVars.find(varName) != intVars.end()){
      leaf = "cut/I";
      addr = reinterpret_cast<void*>(&(cutValsI_[iI++]));
    }
    else{
      leaf = "cut/F";
      addr = reinterpret_cast<void*>(&(cutValsF_[iF++]));
    }

    output_->Branch("cut_" + varName, addr, leaf);
    output_->Branch("isUpper_" + varName, &(upperCut_[iB]), "isUpper/b");
  }

  output_->Branch("signalYield", &signalYield_, "signalYield/F");
  output_->Branch("bkgYield", &bkgYield_, "bkgYield/F");
}

void
CutBasedCalculator::configure(TMap const& _config)
{
  TObjString* stepsStr(dynamic_cast<TObjString*>(_config.GetValue("STEPS")));
  if(!stepsStr)
    throw std::runtime_error("Parameter SCAN.STEPS not set");

  steps_ = stepsStr->GetString().Atoi();
}

void
CutBasedCalculator::initialize()
{
  unsigned nFloat(floatingCuts_.size());

  TCut baseCut;
  for(unsigned iC(0); iC < fixedCuts_.size(); ++iC)
    baseCut *= fixedCuts_[iC];

  std::vector<float> cutVals(nFloat, 0.);
  std::vector<bool> isInt(nFloat, false);
  unsigned nY(1);

  for(unsigned iC(0); iC < nFloat; ++iC){
    cutVals[iC] = bounds_[iC].first;
    nY *= 2;
    if(floatingCuts_[iC](0) == 'n') isInt[iC] = true;
  }

  for(std::map<TString, std::pair<TTree*, double> >::iterator bItr(backgrounds_.begin()); bItr != backgrounds_.end(); ++bItr){
    TTree* sampleTree(bItr->second.first);
    double weight(bItr->second.second);

    for(unsigned iV(0); iV < variables_.size(); ++iV){
      TString const& varName(variables_[iV]);
      double max(sampleTree->GetMaximum(varName));
      double min(sampleTree->GetMinimum(varName));
      if(max > defaultRanges_[iV].second) defaultRanges_[iV].second = max;
      if(min < defaultRanges_[iV].first) defaultRanges_[iV].first = min;
    }

    bool done(false);
    while(!done){
      for(unsigned iY(0); iY < nY; ++iY){
        TCut cut(baseCut);
        for(unsigned iC(0); iC < floatingCuts_.size(); ++iC){
          TString cutStr(floatingCuts_[iC]);
          if(((iY >> iC) & 1) == 1) cutStr.ReplaceAll("?", ">" + TString::Format("%.1f", cutVals[iC]));
          else cutStr.ReplaceAll("?", "<" + TString::Format("%.1f", cutVals[iC]));
          cut *= cutStr;
        }

        bkgYields_[TString(cut)] += sampleTree->GetEntries(cut) * weight;
      }

      done = true;
      for(unsigned iC(0); iC < floatingCuts_.size(); ++iC){
        if(isInt[iC])
          cutVals[iC] += 1.;
        else
          cutVals[iC] += (bounds_[iC].second - bounds_[iC].first) / steps_;

        if(cutVals[iC] >= bounds_[iC].second) cutVals[iC] = bounds_[iC].first;
        else{
          done = false;
          break;
        }
      }
    }
  }

  initialized_ = true;
}

void
CutBasedCalculator::calculate(TString const& _gridParams, TTree* _signal, double _weight)
{
  if(!initialized_) return;

  if(nThreads_ > 1){
    while(nRunning() >= nThreads_) TThread::Sleep(5);

    TThread::Lock();

    std::pair<TTree*, double>& job(jobStack_[_gridParams]);
    job.first = _signal;
    job.second = _weight;

    TThread::UnLock();

    TThread* thread(new TThread(_gridParams, TThread::VoidFunc_t(&CutBasedCalculator::run), reinterpret_cast<void*>(this)));
    thread->Run();
  }
  else{
    std::pair<TTree*, double>& job(jobStack_[_gridParams]);
    job.first = _signal;
    job.second = _weight;

    run();

    if(_signal->InheritsFrom(TChain::Class()))
      delete _signal;
    else
      delete _signal->GetCurrentFile();
  }
}

void
CutBasedCalculator::run()
{
  if(jobStack_.size() == 0) return;

  TString gridParams;
  TTree* signal(0);
  double weight(1.);

  if(nThreads_ > 1){
    gridParams = TThread::Self()->GetName();
    std::map<TString, std::pair<TTree*, double> >::const_iterator jItr(jobStack_.find(gridParams));
    if(jItr == jobStack_.end()){
      std::cerr << "No job info for " << gridParams << " found in the stack" << std::endl;
      return;
    }
    signal = jItr->second.first;
    weight = jItr->second.second;
  }
  else{
    std::map<TString, std::pair<TTree*, double> >::const_iterator jItr(jobStack_.begin());
    gridParams = jItr->first;
    signal = jItr->second.first;
    weight = jItr->second.second;
  }

  double totalSignal(signal->GetEntries());

  if(totalSignal == 0) return;

  TCut baseCut;
  for(unsigned iC(0); iC < fixedCuts_.size(); ++iC)
    baseCut *= fixedCuts_[iC];

  unsigned nFloat(floatingCuts_.size());

  std::vector<float> cutVals(nFloat, 0.);
  std::vector<bool> isInt(nFloat, false);

  unsigned nY(1);
  for(unsigned iC(0); iC < nFloat; ++iC){
    cutVals[iC] = bounds_[iC].first;
    nY *= 2;
    if(floatingCuts_[iC](0) == 'n') isInt[iC] = true;
  }

  std::vector<int> cutValsI(cutValsI_.size(), 0);
  std::vector<float> cutValsF(cutValsF_.size(), 0);

  double significance(0.);
  float signalYield(0.);
  float bkgYield(0.);
  std::vector<unsigned char> upperCut(nFloat, 0);

  bool done(false);
  while(!done){
    for(unsigned iY(0); iY < nY; ++iY){
      TCut cut(baseCut);

      for(unsigned iC(0); iC < nFloat; ++iC){
        TString cutStr(floatingCuts_[iC]);
        if(((iY >> iC) & 1) == 1)
          cutStr.ReplaceAll("?", ">" + TString::Format("%.1f", cutVals[iC]));
        else
          cutStr.ReplaceAll("?", "<" + TString::Format("%.1f", cutVals[iC]));

        cut *= cutStr;
      }

      std::map<TString, float>::const_iterator byItr(bkgYields_.find(TString(cut)));
      if(byItr == bkgYields_.end())
        throw std::runtime_error(("Background yield for " + TString(cut) + " not calculated").Data());
      float localBkg(byItr->second);

      TThread::Lock();

      float localSignal(signal->GetEntries(cut) * weight);

      TThread::UnLock();

      double localSigma(localSignal / std::sqrt(localSignal + localBkg));

      if(localSigma > significance){
        significance = localSigma;
        signalYield = localSignal;
        bkgYield = localBkg;
        unsigned iI(0);
        unsigned iF(0);

        for(unsigned iC(0); iC < nFloat; ++iC){
          upperCut[iC] = iY >> iC & 1;
          if(isInt[iC]) cutValsI[iI++] = cutVals[iC];
          else cutValsF[iF++] = cutVals[iC];
        }
      }
    }

    done = true;

    for(unsigned iC(0); iC < nFloat; ++iC){
      if(isInt[iC])
        cutVals[iC] += 1.;
      else
        cutVals[iC] += (bounds_[iC].second - bounds_[iC].first) / steps_;

      if(cutVals[iC] >= bounds_[iC].second) cutVals[iC] = bounds_[iC].first;
      else{
        done = false;
        break;
      }
    }
  }

  TThread::Lock();

  output_->GetDirectory()->cd();

  std::vector<THStack*> stacks;
  std::vector<TLegend*> legends;
  for(unsigned iV(0); iV < variables_.size(); ++iV){
    TString const& varName(variables_[iV]);
    TCanvas* canvas(distributions_[iV]);
    canvas->cd();
    canvas->Clear();

    std::pair<double, double> range(defaultRanges_[iV]);
    double max(signal->GetMaximum(varName));
    double min(signal->GetMinimum(varName));
    if(max > range.second) range.second = max;
    if(min < range.first) range.first = min;
    TString rangeStr("(100," + TString::Format("%f", range.first) + "," + TString::Format("%f", range.second) + ")");

    TCut nMinus1;
    for(unsigned iC(0); iC < fixedCuts_.size(); ++iC)
      if(!TString(fixedCuts_[iC]).Contains(varName)) nMinus1 *= fixedCuts_[iC];

    THStack* stack = new THStack(varName, varName);
    TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);

    unsigned iS(0);
    TString hName;
    for(std::map<TString, std::pair<TTree*, double> >::const_iterator bItr(backgrounds_.begin()); bItr != backgrounds_.end(); ++bItr){
      TString const& sampleName(bItr->first);
      TTree* tree(bItr->second.first);

      hName = "h_" + varName + "_" + sampleName;
      tree->Draw(varName + ">>" + hName + rangeStr, nMinus1);

      TH1F* histo = static_cast<TH1F*>(gDirectory->Get(hName));
      histo->Scale(bItr->second.second);
      histo->SetLineColor(2 + iS);
      histo->SetFillColor(2 + iS);
      histo->SetFillStyle(1001);
      ++iS;

      stack->Add(histo);
      legend->AddEntry(histo, sampleName);
    }

    hName = "h_" + varName;
    signal->Draw(varName + ">>" + hName + rangeStr, nMinus1);

    TH1F* histo = static_cast<TH1F*>(gDirectory->Get(hName));
    histo->Scale(weight);
    histo->SetLineColor(kBlack);

    legend->AddEntry(histo, gridParams);

    stack->SetTitle(varName);

    stack->Draw();
    histo->Draw("same");
    legend->Draw();

    stacks.push_back(stack);
    legends.push_back(legend);
  }

  *gridParams_ = gridParams;
  significance_ = significance;
  signalEfficiency_ = signalYield / totalSignal;
  cutValsI_ = cutValsI;
  cutValsF_ = cutValsF;
  upperCut_ = upperCut;
  signalYield_ = signalYield;
  bkgYield_ = bkgYield;

  output_->Fill();

  for(unsigned iS(0); iS < stacks.size(); ++iS){
    delete stacks[iS];
    delete legends[iS];
  }

  TThread::UnLock();
}
