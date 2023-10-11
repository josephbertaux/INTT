#ifndef INTT_INTTCOMBINEDRAWDATACONVERTEREB_H
#define INTT_INTTCOMBINEDRAWDATACONVERTEREB_H

#include <fun4all/SubsysReco.h>

#include <Rtypes.h> // for ROOT data types

#include <map>
#include <string>
#include <tuple>
#include <vector>

class PHCompositeNode;
class TTree;
class TFile;

class InttEvent;

class InttCombinedRawDataConverterEB : public SubsysReco
{
 public:
  InttCombinedRawDataConverterEB(std::string const& name = "InttRawDataConverterEB");

  int SetOutputFile(std::string const&);
  int WriteOutputFile();

  int Init(PHCompositeNode*) override;
  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

 private:
  std::string m_InttRawNodeName = "INTTRAWHIT";

  TFile* file = nullptr;
  TTree* tree = nullptr;

  InttEvent* inttEvt;

  Int_t n_evt = 0;
  Int_t num_hits = 0;

  typedef std::map<std::string, std::vector<Int_t>*> Branches_i_t;
  typedef std::map<std::string, std::vector<Long64_t>*> Branches_l_t;
  typedef std::map<std::string, std::vector<Double_t>*> Branches_d_t;
  Branches_i_t branches_i;
  Branches_l_t branches_l;
  Branches_d_t branches_d;
};

#endif  // INTT_INTTCOMBINEDRAWDATACONVERTEREB_H
