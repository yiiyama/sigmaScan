#ifndef SignificanceCalculator_h
#define SignificanceCalculator_h

#include <map>
#include <utility>

#include "TString.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TMap.h"

class SignificanceCalculator {
 public:
  SignificanceCalculator(TString const&, std::vector<TString> const&, unsigned);
  virtual ~SignificanceCalculator();
  virtual void configure(TMap const&) {}
  virtual void setBackground(TString const& _name, std::pair<TTree*, double> const& _treeAndWeight) { if(!initialized_) backgrounds_[_name] = _treeAndWeight; };
  virtual void initialize() = 0;
  virtual void calculate(TString const&, TTree*, double) = 0;
  virtual void setVerbose(bool _v) { verbose_ = _v; }
  virtual unsigned nRunning();
 protected:
  virtual void run() = 0;

  std::map<TString, std::pair<TTree*, double> > backgrounds_;
  TString* gridParams_;
  float significance_;
  float signalEfficiency_;
  std::vector<TString> variables_;
  std::vector<TCut> fixedCuts_;
  std::vector<TCanvas*> distributions_;
  std::vector<std::pair<double, double> > defaultRanges_;
  TTree* output_;
  unsigned nThreads_;
  std::map<TString, std::pair<TTree*, double> > jobStack_;
  bool initialized_;
  bool verbose_;
};

class CutBasedCalculator : public SignificanceCalculator {
 public:
  CutBasedCalculator(TString const&, std::vector<TString> const&, std::vector<TString> const&, unsigned = 1);
  ~CutBasedCalculator() {}
  void configure(TMap const&);
  void initialize();
  void calculate(TString const&, TTree*, double);
 private:
  void run();

  int steps_;
  std::vector<std::pair<float, float> > bounds_;
  std::vector<TString> floatingCuts_;
  std::vector<int> cutValsI_;
  std::vector<float> cutValsF_;
  std::vector<unsigned char> upperCut_;
  float signalYield_;
  float bkgYield_;
  std::map<TString, float> bkgYields_;
};

#endif
