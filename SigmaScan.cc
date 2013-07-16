#include "CutVarTreeProducer.h"
#include "SignificanceCalculator.h"

#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <map>

#include "TMap.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TPRegexp.h"
#include "TFile.h"
#include "TTree.h"
#include "TThread.h"

class SigmaScan {
public:
  SigmaScan() : config_() { config_.SetOwner(); }
  ~SigmaScan() {}

  void addParam(TString const&, TString const&);
  void addVariable(TString const&);
  void produceTrees() const;
  void scan() const;
private:
  TMap config_;
};

void
SigmaScan::addParam(TString const& _name, TString const& _value)
{
  TString name(_name);
  int dot(0);
  TMap* config(&config_);
  while((dot = name.Index(".")) > 0){
    TString subConfigName(name(0, dot));
    TMap* subConfig(dynamic_cast<TMap*>(config->GetValue(subConfigName)));
    if(!subConfig){
      TObjString* key(new TObjString(subConfigName));
      subConfig = new TMap;
      config->Add(key, subConfig);
    }
    config = subConfig;
    name.Remove(0, dot + 1);
  }
  config->Add(new TObjString(name), new TObjString(_value));
}

void
SigmaScan::addVariable(TString const& _nameAndRange)
{
  TObjArray* variables(dynamic_cast<TObjArray*>(config_.GetValue("VARIABLES")));
  if(!variables){
    variables = new TObjArray;
    variables->SetOwner();
    config_.Add(new TObjString("VARIABLES"), variables);
  }

  variables->Add(new TObjString(_nameAndRange));
}

void
SigmaScan::produceTrees() const
{
  TObjString* confStr(0);
  TMap* subConfig(0);

  confStr = dynamic_cast<TObjString*>(config_.GetValue("WORKSPACE"));
  if(!confStr){
    std::cerr << "Parameter WORKSPACE not set" << std::endl;
    return;
  }
  TString const& workSpace(confStr->GetString());

  unsigned nCPU(1);
  confStr = dynamic_cast<TObjString*>(config_.GetValue("NCPU"));
  if(confStr)
    nCPU = confStr->GetString().Atoi();

  void* dirp(gSystem->OpenDirectory(workSpace));
  if(dirp)
    std::cout << "Reusing existing directory " << workSpace << std::endl;
  else{
    if(gSystem->mkdir(workSpace, true) != 0){
      std::cerr << "Cannot create workspace " << workSpace << std::endl;
      return;
    }
  }
  gSystem->FreeDirectory(dirp);

  subConfig = dynamic_cast<TMap*>(config_.GetValue("SIGNAL"));
  if(subConfig){
    TMap const& sigConfSet(*subConfig);

    TMapIter* sigItr(static_cast<TMapIter*>(sigConfSet.MakeIterator()));
    TObject* key(0);
    while((key = sigItr->Next())){
      TString sigName(key->GetName());

      TString sigDir(workSpace + "/" + sigName);
      dirp = gSystem->OpenDirectory(sigDir);
      if(dirp)
        std::cout << "Directory " << sigDir << " already exists and will not be reprocessed." << std::endl;
      else{
        std::cout << "Producing " << sigName << " trees.." << std::endl;
        if(gSystem->mkdir(sigDir) != 0){
          std::cerr << "Cannot create directory " << sigDir << std::endl;
          return;
        }

        subConfig = dynamic_cast<TMap*>(sigConfSet.GetValue(key));
        if(!subConfig){
          std::cerr << "BACKGROUND." << sigName << " is not a ParameterSet" << std::endl;
          return;
        }
        TMap const& sigConf(*subConfig);

        confStr = dynamic_cast<TObjString*>(sigConf.GetValue("LIST"));
        if(!confStr){
          std::cerr << "Parameter SIGNAL.LIST not set" << std::endl;
          return;
        }
        TString const& signalList(confStr->GetString());

        confStr = dynamic_cast<TObjString*>(sigConf.GetValue("DECAY"));
        if(!confStr){
          std::cerr << "Parameter SIGNAL.DECAY not set" << std::endl;
          return;
        }
        TString const& signalDecay(confStr->GetString());

        confStr = dynamic_cast<TObjString*>(sigConf.GetValue("GRIDPARAMPATTERN"));
        if(!confStr){
          std::cerr << "Parameter SIGNAL.GRIDPARAMPATTERN not set" << std::endl;
          return;
        }
        TString const& paramPatStr(confStr->GetString());
        TPRegexp paramPat(paramPatStr);

        CutVarTreeProducer producer(signalDecay, nCPU);

        std::ifstream list(signalList);
        if(!list.is_open()){
          std::cerr << "Cannot open file list" << std::endl;
          return;
        }

        std::string buf;
        TString line;
        while(true){
          std::getline(list, buf);
          if(!list.good()) break;

          line = buf;

          TObjArray* matches(paramPat.MatchS(gSystem->BaseName(line)));
          if(matches->GetEntries() < 2){
            delete matches;
            std::cerr << "Cannot extract gridParams from file name " << line << " using pattern " << paramPatStr << std::endl;
            continue;
          }

          TString gridParams;
          int iP(1);
          while(true){
            gridParams += matches->At(iP)->GetName();
            if(++iP == matches->GetEntries()) break;
            else gridParams += "_";
          }
          TString outputName(sigDir + "/signal_" + gridParams);
          outputName += ".root";

          delete matches;

          TThread::Lock();

          TFile* source(TFile::Open(line));

          TTree* tree(0);
          if(!source || source->IsZombie() || (tree = dynamic_cast<TTree*>(source->Get("susyTree"))) == 0){
            delete source;
            std::cerr << "Cannot process signal file " << line << std::endl;
            TThread::UnLock();
            return;
          }

          TThread::UnLock();

          std::cout << "Signal: " << gridParams << std::endl;

          producer.produce(*tree, outputName);
        }

        while(producer.nRunning() > 0) TThread::Sleep(5);

      }
      gSystem->FreeDirectory(dirp);
    }
    delete sigItr;
  }

  subConfig = dynamic_cast<TMap*>(config_.GetValue("BACKGROUND"));
  if(subConfig){
    TMap const& bkgConfSet(*subConfig);

    TMapIter* bkgItr(static_cast<TMapIter*>(bkgConfSet.MakeIterator()));
    TObject* key(0);
    while((key = bkgItr->Next())){
      TString bkgName(key->GetName());

      TString bkgDir(workSpace + "/" + bkgName);
      dirp = gSystem->OpenDirectory(bkgDir);
      if(dirp)
        std::cout << "Directory " << bkgDir << " already exists and will not be reprocessed." << std::endl;
      else{
	std::cout << "Producing " << bkgName << " trees.." << std::endl;
	if(gSystem->mkdir(bkgDir) != 0){
	  std::cerr << "Cannot create directory " << bkgDir << std::endl;
	  return;
	}

	subConfig = dynamic_cast<TMap*>(bkgConfSet.GetValue(key));
	if(!subConfig){
	  std::cerr << "BACKGROUND." << bkgName << " is not a ParameterSet" << std::endl;
	  return;
	}
	TMap const& bkgConf(*subConfig);

	confStr = dynamic_cast<TObjString*>(bkgConf.GetValue("LIST"));
	if(!confStr){
	  std::cerr << "Parameter BACKGROUND." << bkgName << ".LIST not set" << std::endl;
	  return;
	}
	TString const& bkgList(confStr->GetString());

	confStr = dynamic_cast<TObjString*>(bkgConf.GetValue("DECAY"));
	if(!confStr){
	  std::cerr << "Parameter BACKGROUND." << bkgName << ".DECAY not set" << std::endl;
	  return;
	}
	TString const& bkgDecay(confStr->GetString());

        CutVarTreeProducer producer(bkgDecay, nCPU);

	std::ifstream list(bkgList);
	if(!list.is_open()){
	  std::cerr << "Cannot open list " << bkgList << std::endl;
	  return;
	}

        unsigned iFile(0);
	std::string buf;
        TString line;
        while(true){
          std::getline(list, buf);
          if(!list.good()) break;

          line = buf;

          TThread::Lock();

          TFile* source(TFile::Open(line));
          TTree* tree(0);
          if(!source || source->IsZombie() || (tree = dynamic_cast<TTree*>(source->Get("susyTree"))) == 0){
            delete source;
            std::cerr << "Cannot process signal file " << line << std::endl;
            TThread::UnLock();
            return;
          }

          TThread::UnLock();

          TString prodId(bkgName + "_" + TString::Format("%d", iFile));
          TString outputName(bkgDir + "/bkg_" + prodId + ".root");

          std::cout << "BKG: " << prodId << std::endl;

          producer.produce(*tree, outputName);

          ++iFile;
        }

        while(producer.nRunning() > 0) TThread::Sleep(5);

      }
      gSystem->FreeDirectory(dirp);
    }
    delete bkgItr;
  }
}

void
SigmaScan::scan() const
{
  TObjString* confStr(0);
  TMap* subConfig(0);
  TObjArray* listConfig(0);

  confStr = dynamic_cast<TObjString*>(config_.GetValue("WORKSPACE"));
  if(!confStr){
    std::cerr << "Parameter WORKSPACE not set" << std::endl;
    return;
  }
  TString const& workSpace(confStr->GetString());

  unsigned nCPU(1);
  confStr = dynamic_cast<TObjString*>(config_.GetValue("NCPU"));
  if(confStr)
    nCPU = confStr->GetString().Atoi();

  subConfig = dynamic_cast<TMap*>(config_.GetValue("SCAN"));
  if(!subConfig){
    std::cerr << "ParameterSet SCAN not set" << std::endl;
    return;
  }
  TMap const& scanConfigs(*subConfig);

  confStr = dynamic_cast<TObjString*>(scanConfigs.GetValue("NAME"));
  if(!confStr){
    std::cerr << "Parameter SCAN.NAME not set" << std::endl;
    return;
  }
  TString const& scanName(confStr->GetString());

  confStr = dynamic_cast<TObjString*>(scanConfigs.GetValue("METHOD"));
  if(!confStr){
    std::cerr << "Parameter SCAN.METHOD not set" << std::endl;
    return;
  }
  TString const& calculationMethod(confStr->GetString());

  confStr = dynamic_cast<TObjString*>(config_.GetValue("INTEGRATEDLUMI"));
  if(!confStr){
    std::cerr << "Parameter INTEGRATEDLUMI not set" << std::endl;
    return;
  }
  double intL(confStr->GetString().Atof());

  listConfig = dynamic_cast<TObjArray*>(config_.GetValue("VARIABLES"));
  if(!listConfig){
    std::cerr << "ParameterSet VARIABLES not set" << std::endl;
    return;
  }
  TObjArray const& varConfs(*listConfig);

  std::vector<TString> variables(0);
  std::vector<TString> varRanges(0);
  TPRegexp varPat("([a-zA-Z0-9]+)(\\[[0-9.-]*(?:|[:][0-9.-]*)\\])");
  for(int iV(0); iV < varConfs.GetEntries(); ++iV){
    TString varConf(varConfs[iV]->GetName());
    TObjArray* matches(varPat.MatchS(varConf));
    if(matches->GetEntries() > 0){
      TString varName(matches->At(1)->GetName());
      bool exists(false);
      for(unsigned i(0); i < sizeof(CutVariables)/sizeof(TString); ++i){
        if(CutVariables[i] == varName){
          exists = true;
          break;
        }
      }
      if(!exists){
        std::cerr << "Variable " << varName << " does not exist" << std::endl;
        continue;
      }
      variables.push_back(varName);
      varRanges.push_back(matches->At(2)->GetName());
    }
    delete matches;
  }

  subConfig = dynamic_cast<TMap*>(config_.GetValue("SIGNAL"));
  if(!subConfig){
    std::cerr << "ParameterSet SIGNAL not set" << std::endl;
    return;
  }
  TMap const& signalConf(*subConfig);

  confStr = dynamic_cast<TObjString*>(signalConf.GetValue("NAME"));
  if(!confStr){
    std::cerr << "Parameter SIGNAL.NAME not set" << std::endl;
    return;
  }
  TString const& signalName(confStr->GetString());

  TString signalDir(workSpace + "/" + signalName);
  void* signalDirP(gSystem->OpenDirectory(signalDir));
  if(!signalDirP){
    std::cerr << "Signal files do not exist in " << workSpace << std::endl;
    return;
  }

  confStr = dynamic_cast<TObjString*>(signalConf.GetValue("XSECFILE"));
  if(!confStr){
    std::cerr << "Parameter SIGNAL.XSECFILE not set" << std::endl;
    return;
  }
  TString const& xsecFileName(confStr->GetString());

  confStr = dynamic_cast<TObjString*>(signalConf.GetValue("EFFICIENCY"));
  if(!confStr){
    std::cerr << "Parameter SIGNAL.EFFICIENCY not set" << std::endl;
    return;
  }
  double signalEfficiency(confStr->GetString().Atof());

  std::map<TString, double> signalWeights;
  std::ifstream xsecFile(xsecFileName);
  if(!xsecFile.is_open()){
    std::cerr << "Cannot open cross section file " << xsecFileName << std::endl;
    return;
  }

  std::stringstream bufstream;
  TString gridParams;
  double xsec;
  int N;
  std::string buf;
  while(true){
    std::getline(xsecFile, buf);
    if(!xsecFile.good()) break;
    bufstream.clear();
    bufstream.str("");
    bufstream << buf;
    bufstream >> gridParams >> xsec >> N;
    if(gridParams == "" || xsec == 0. || N == 0){
      std::cerr << "Incorrect format used in the signal cross section file" << std::endl;
      return;
    }
    signalWeights[gridParams] = intL / (N / xsec) * signalEfficiency;
  }
  xsecFile.close();

  subConfig = dynamic_cast<TMap*>(config_.GetValue("BACKGROUND"));
  if(!subConfig){
    std::cerr << "ParameterSet BACKGROUND not set" << std::endl;
    return;
  }
  TMap const& bkgConfSet(*subConfig);

  std::map<TString, std::pair<TTree*, double> > bkgSamples;

  TMapIter* bkgItr(static_cast<TMapIter*>(bkgConfSet.MakeIterator()));
  TObject* key(0);
  while((key = bkgItr->Next())){
    TString bkgName(key->GetName());

    TChain* chain(new TChain("eventTree"));
    if(chain->Add(workSpace + "/" + bkgName + "/bkg_*.root") == 0){
      std::cerr << "No files for " << bkgName << std::endl;
      continue;
    }

    subConfig = dynamic_cast<TMap*>(bkgConfSet.GetValue(key));
    if(!subConfig){
      std::cerr << "BACKGROUND." << bkgName << " is not a ParameterSet" << std::endl;
      return;
    }
    TMap const& bkgConf(*subConfig);

    confStr = dynamic_cast<TObjString*>(bkgConf.GetValue("XSEC"));
    if(!confStr){
      std::cerr << "Parameter BACKGROUND." << bkgName << ".XSEC not set" << std::endl;
      return;
    }
    double bkgXSec(confStr->GetString().Atof());

    confStr = dynamic_cast<TObjString*>(bkgConf.GetValue("NEVENTS"));
    if(!confStr){
      std::cerr << "Parameter BACKGROUND." << bkgName << ".NEVENTS not set" << std::endl;
      return;
    }
    double bkgEvents(confStr->GetString().Atof());

    confStr = dynamic_cast<TObjString*>(bkgConf.GetValue("FAKERATE"));
    if(!confStr){
      std::cerr << "Parameter BACKGROUND." << bkgName << ".FAKERATE not set" << std::endl;
      return;
    }
    double fakeRate(confStr->GetString().Atof());

    double weight(intL / (bkgEvents / bkgXSec) * fakeRate);

    bkgSamples[bkgName] = std::pair<TTree*, double>(chain, weight);
  }

  delete bkgItr;

  SignificanceCalculator* calculator(0);
  try{
    if(calculationMethod == "CutBased")
      calculator = new CutBasedCalculator(workSpace + "/" + scanName + ".root", variables, varRanges, nCPU);
    else{
      std::cerr << "Calculation method " << calculationMethod << " not implemented" << std::endl;
      return;
    }
    calculator->configure(scanConfigs);
  }
  catch(std::exception& e){
    std::cerr << "Cannot configure significance calculator: " << e.what() << std::endl;
    delete calculator;
    return;
  }

  for(std::map<TString, std::pair<TTree*, double> >::iterator bItr(bkgSamples.begin()); bItr != bkgSamples.end(); ++bItr)
    calculator->setBackground(bItr->first, bItr->second);

  std::cout << "Calculating background yields" << std::endl;

  calculator->initialize();

  std::cout << " Done." << std::endl;

  TPRegexp namePat("signal_([0-9a-zA-Z_]+).root");
  TString fileName;
  while((fileName = gSystem->GetDirEntry(signalDirP)) != ""){
    TObjArray* matches(namePat.MatchS(fileName));
    if(matches->GetEntries() == 0){
      delete matches;
      continue;
    }

    gridParams = matches->At(1)->GetName();
    delete matches;

    if(signalWeights.find(gridParams) == signalWeights.end()){
      std::cerr << "No weight inforamtion found for signal point " << gridParams << ", skipping." << std::endl;
      continue;
    }

    TThread::Lock();

    TFile* source(TFile::Open(signalDir + "/" + fileName));
    TTree* signal(0);
    if(!source || source->IsZombie() || (signal = dynamic_cast<TTree*>(source->Get("eventTree"))) == 0){
      delete source;
      std::cerr << "Cannot process event file " << fileName << std::endl;
      continue;
    }

    TThread::UnLock();

    std::cout << "Calculating " << gridParams << std::endl;

    try{
      calculator->calculate(gridParams, signal, signalWeights[gridParams]);
    }
    catch(std::exception& e){
      std::cerr << "Signal point " << gridParams << " skipped due to errors: " << e.what() << std::endl;
    }
  }

  while(calculator->nRunning() > 0) TThread::Sleep(5);

  delete calculator;

  for(std::map<TString, std::pair<TTree*, double> >::iterator bItr(bkgSamples.begin()); bItr != bkgSamples.end(); ++bItr)
    delete bItr->second.first;
}
