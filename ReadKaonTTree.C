/**
  * From Jason Brylawskyj
  *
  * \Update to analyze 3p85_FXT_2018
  * \author Ding Chen
  * \date Feb 24, 2020
  *
  * \Upadate to analyze 26p5_FXT_2018 data, new centrality def, new kaon TTree
  * \author Ding Chen
  * \date Jun 29, 2020
  *
  * \Update to do the systematic analysis 26p5_FXT_2018
  * \author Ding Chen
  * \date Sep 28, 2020
  */

// C++ classes
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include <map>

// ROOT classes
#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStreamerElement.h>
#include <TStyle.h>
#include "Riostream.h"
#include "TSystemDirectory.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TChain.h"
#include "TLegend.h"
#include "TFractionFitter.h"
#include "TVirtualFitter.h"
#include "TCut.h"
#include "TObject.h"
#include "TLine.h"
#include "TLeaf.h"
#include "TGaxis.h"
#include "TSpectrum.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TSpline.h"
#include "TPaveText.h"
#include "Fit/FitResult.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom3.h"
#include "TDatime.h"
#include "Math/MinimizerOptions.h"

using namespace std;
const Int_t _Ncentralities = 9; // 9 centrality bins
const double _d_K_m        = 0.493677;
const double _m_phi        = 1.019461;
// ======================== (1) prerequsites for read kaon TTree ===============
struct st_track // Useful Physics information of a track (from kaon TTree)
{
  bool  b_pos_charge;
  bool  b_bad_TOF;

  int   runNumber;
  int   eventNumber;
  int   nGoodTracks;
  int   nFXTMult;

  double  px;
  double  py;
  double  pz;

  double  x_vtx;
  double  y_vtx;
  double  z_vtx;
  double  r_vtx;
  double  DCA_r;

  double  tofBeta;
  double  mass2;

  double  nHitsDedx;
  double  nHitsFit;
  double  nHitsRatio;

  double  d_TPCnSigmaKaon;
  double  d_TOFnSigmaKaon;
};
struct st_event // Positive and negative tracks of an event
{
  vector<st_track> v_trk_pl;
  vector<st_track> v_trk_mi;
};
TChain * t_K   = new TChain("t_K");

// =============================== (2) Analysis Start ==========================
void ReadKaonTTree( string FileName,
                   TString outFile = "test",
                    // double inputParameter1 = 0.
                    Int_t   inputp2 = 0, // sysErr cut Indexes 0-15
                    Int_t   inputp3 = 0, // sysErr cut variations, each systematic check has 2 or 3 vertions
                    Int_t   inputp4 = 0 // Iteration of the analysis is. In this analysis, 2 iterations is enough
){
  // ----------------------------- read TTree ----------------------------------
  if( FileName.find(".list") != string::npos ||
      FileName.find(".lis") != string::npos ) {
    std::ifstream inputStream( FileName.c_str() );
    if(!inputStream) {
      cout << "ERROR: Cannot open list file " << FileName << endl;
    }
    Int_t nFile = 0;
    string file;
    while(getline(inputStream, file)) {
      if(file.find(".root") != string::npos) {
        TFile* ftmp = TFile::Open(file.c_str());
        if(ftmp && !ftmp->IsZombie() && ftmp->GetNkeys()) {
          cout << " Read in root file " << file << endl;
          t_K->Add(file.c_str());
          ++nFile;
        } //if(ftmp && !ftmp->IsZombie() && ftmp->GetNkeys())
        if (ftmp) {
          ftmp->Close();
        } //if (ftmp)
      } //if(file.find(".picoDst.root") != string::npos)
    } //while (getline(inputStream, file))
    cout << " Total " << nFile << " files have been read in. " << endl;
  } //if(FileName.find(".list") != string::npos || FileName.find(".lis" != string::npos))
  else if(FileName.find(".root") != string::npos) {
    t_K->Add(FileName.c_str());
  }
  else {
    cout << " No good input file to read ... " << endl;
  }
  if( t_K->GetEntries() == 0 ){cout << "Error:Could not find 't_K' in file: " << FileName << endl; return;}

  Int_t sys_cutN = inputp2; // sysErr cut Indexes 0-15
  Int_t sys_varN = inputp3; // sysErr cut variations, each systematic check has 2 or 3 vertions
  Int_t sys_iterN = inputp4; // Iteration of the analysis is. In this analysis, 2 iterations is enough
  string sys_object[17]  = {"primary", "etaGap", "etaRange",
                            "vz", "vr", "dedx", "dca",
                            "nHitsFit", "ratio", "nSigK", "mass2",
                            "pT", "dipAngle", "vtxDiff", "mthdDiff",
                            "binning",
                            "TPCpid"};
  std::cout << "sys_cutN == "<< sys_cutN <<": "<< sys_object[sys_cutN] << std::endl;
  // systematic cut Levels
  double d_zvtxCutLvlLow,d_zvtxCutLvlHigh,d_rvtxCutLvl,d_nHitsDedxCutLvl,d_DCACutLvl,d_nHitsFitCutLvl,d_ratioCutLvl;
  double d_nSigmaKaonCutLvl,d_KaonM2CutLvllow,d_KaonM2CutLvlHigh,d_KaonpTCutLvlLow;
  double d_dipAngleCutLvl,d_vtxDiffCutLvl;
  bool b_TPCpid = false;
  // default cut level - primary
  d_zvtxCutLvlLow = 198.0, d_zvtxCutLvlHigh = 202.0;
  d_rvtxCutLvl = 2.0;
  d_nHitsDedxCutLvl = 0;
  d_DCACutLvl = 3.0;
  d_nHitsFitCutLvl = 15;
  d_ratioCutLvl = 0.51;
  d_nSigmaKaonCutLvl=2.0, d_KaonM2CutLvllow=0.16,d_KaonM2CutLvlHigh=0.32, d_KaonpTCutLvlLow=0.2 ;
  d_dipAngleCutLvl=0.04,d_vtxDiffCutLvl=0.5;
  // # Systematic Analysis
  if(sys_cutN == 3){
    // sys_cutN == 3; // vz
    if(sys_varN == 1){
      d_zvtxCutLvlLow = 198.4, d_zvtxCutLvlHigh = 201.6;
    } else if(sys_varN == 2){
      d_zvtxCutLvlLow = 197.6, d_zvtxCutLvlHigh = 202.4;
    }
  } else if(sys_cutN == 4){
    // sys_cutN == 4; // vr
    if(sys_varN == 1){
      d_rvtxCutLvl = 1.6;
    } else if(sys_varN == 2){
      d_rvtxCutLvl = 2.4;
    }
  } else if(sys_cutN == 5){
    // sys_cutN == 5; // dedx
    if(sys_varN == 1){
      d_nHitsDedxCutLvl = 10;
    } else if(sys_varN == 2){
      d_nHitsDedxCutLvl = 20;
    }
  } else if(sys_cutN == 6){
    // sys_cutN == 6; // dca
    if(sys_varN == 1){
      d_DCACutLvl = 1.0;
    } else if(sys_varN == 2){
      d_DCACutLvl = -999.0;
       // no DCA cut
    }
  } else if(sys_cutN == 7){
    // sys_cutN == 7; // nHitsFit
    if(sys_varN == 1){
      d_nHitsFitCutLvl = 10;
    } else if(sys_varN == 2){
      d_nHitsFitCutLvl = 20;
    }
  } else if(sys_cutN == 8){
    // sys_cutN == 8; // ratio
    if(sys_varN == 1){
      d_ratioCutLvl = 0.45;
    } else if(sys_varN == 2){
      d_ratioCutLvl = 0.55;
    }
  } else if(sys_cutN == 9){
    // sys_cutN == 9; // nSigmaKaon
    if(sys_varN == 1){
      d_nSigmaKaonCutLvl = 1.8;
    } else if(sys_varN == 2){
      d_nSigmaKaonCutLvl = 2.2;
    }
  } else if(sys_cutN == 10){
    // sys_cutN == 10; // Mass2
    if(sys_varN == 1){
      d_KaonM2CutLvllow     = 0.17;
      d_KaonM2CutLvlHigh    = 0.31;
    } else if(sys_varN == 2){
      d_KaonM2CutLvllow     = 0.15;
      d_KaonM2CutLvlHigh    = 0.33;
    }
  } else if(sys_cutN == 11){
    // sys_cutN == 11; // pTlow
    if(sys_varN == 1){
      d_KaonpTCutLvlLow     = 0.0;
    } else if(sys_varN == 2){
      d_KaonpTCutLvlLow     = 0.4;
    }
  } else if(sys_cutN == 12){
    // sys_cutN == 12; // dip angle
    d_dipAngleCutLvl = 0.0; // no dip angle cut
  } else if(sys_cutN == 13){
    // sys_cutN == 13; // vtxDiff
    if(sys_varN == 1){
      d_vtxDiffCutLvl     = 0.2;
    } else if(sys_varN == 2){
      d_vtxDiffCutLvl     = 1.0;
    }
  } else if(sys_cutN == 16){
    // sys_cutN == 16; // TPCpid
    b_TPCpid = true;
    if(sys_varN == 0){
      d_nSigmaKaonCut = 2.0;
    } else if(sys_varN == 1){
      d_nSigmaKaonCut = 3.0;
    } else if(sys_varN == 2){
      d_nSigmaKaonCut = 4.0;
    }
  }

  // multi job for ME using script
  // Double_t d_entryRange = inputParameter1;
  int N_entries = t_K -> GetEntries();
  // Get the number of entries in the TTree

  // int N_entries = (((d_entryRange+1)*5550) < (t_K -> GetEntries())) ? ((d_entryRange+1)*5550) :  (t_K -> GetEntries());

  // ---------------------------- Output file and Hists ------------------------
  outFile.Prepend(Form("_var%d_iter%d_", sys_varN, sys_iterN));
  outFile.Prepend(sys_object[sys_cutN]);
  outFile.Prepend("sys_");
  outFile.Append(".readKTree.result.root");
  TFile * tf_out = new TFile(outFile,"RECREATE");
  TH1D * hist_dip_angle = new TH1D("hist_dip_angle","hist_dip_angle",1000,-1,1.0);
  TH1D * hist_SE_mass_Phi     = new TH1D("hist_SE_mass_Phi","Same event invariant mass",200,0.9,1.1);
  TH1D * hist_ME_mass_Phi     = new TH1D("hist_ME_mass_Phi","Mixed event invariant mass",200,0.9,1.1);
  // ------------------------------- QA plots ----------------------------------
  TH1D *hist_SE_PhiMeson_pT  = new TH1D("hist_SE_PhiMeson_pT","pT distribution of #phi",200,0.0,10);
  TH1D *hist_SE_PhiMeson_mT  = new TH1D("hist_SE_PhiMeson_mT","mT distribution of #phi",200,0.0,10);
  TH1D *hist_SE_PhiMeson_rap  = new TH1D("hist_SE_PhiMeson_rap","y distribution of #phi",200,-10.,10);
  TH1D *hist_PhiMeson_tight_rap  = new TH1D("hist_PhiMeson_tight_rap","y distribution with tight cuts of #phi",200,-10.0,10.0);
  TH2D *hist_SE_pt_y_PhiMeson = new TH2D("hist_SE_pt_y_PhiMeson","p_{T} [GeV/c] vs. y of #phi",500,-3.0,0.5,500,0.0,3.5);
  TH2D * h2_TOF_beta_pq                  = new TH2D("h2_TOF_beta_pq","1/#beta vs. pq",500,-3,3,500,0,3);
  TH2D * h2_dEdx_pq                      = new TH2D("h2_dEdx_pq","dEdx vs. pq",500,-3,3,500,2e-6,10e-6);
  TH3D * h3_dEdx_pq_vs_nsig              = new TH3D("h3_dEdx_pq_vs_nsig","nSigmaKaon vs. dEdx vs. pq",500,-3,3,500,2e-6,10e-6,100,0.0,10.0);
  // waiting for dedx value in the K TTree for plots that need dedx info
  TH1D *hist_KaonPlus_pT  = new TH1D("hist_KaonPlus_pT","pT distribution of K^{+}",500,0.0,3.5);
  TH1D *hist_KaonMinus_pT = new TH1D("hist_KaonMinus_pT","pT distribution of K^{-}",500,0.0,3.5);
  TH1D *hist_KaonPlus_rap  = new TH1D("hist_KaonPlus_rap","y distribution of K^{+}",500,-3.0,0.5);
  TH1D *hist_KaonMinus_rap = new TH1D("hist_KaonMinus_rap","y distribution of K^{-}",500,-3.0,0.5);
  TH1D *hist_KaonPlus_eta  = new TH1D("hist_KaonPlus_eta","#eta distribution of K^{+}",500,-3.0,0.5);
  TH1D *hist_KaonMinus_eta = new TH1D("hist_KaonMinus_eta","#eta distribution of K^{-}",500,-3.0,0.5);
  TH2D *hist_pt_y_kaonPlus = new TH2D("hist_pt_y_kaonPlus","p_{T} [GeV/c] vs. y of K^{+}",500,-3.0,0.5,500,0.0,3.5);
  TH2D *hist_pt_y_kaonMinus = new TH2D("hist_pt_y_kaonMinus","p_{T} [GeV/c] vs. y of K^{-}",500,-3.0,0.5,500,0.0,3.5);
  TH2D *hist_pt_eta_kaonPlus = new TH2D("hist_pt_eta_kaonPlus","p_{T} [GeV/c] vs. #eta of K^{+}",500,-3.0,0.5,500,0.0,3.5);
  TH2D *hist_pt_eta_kaonMinus = new TH2D("hist_pt_eta_kaonMinus","p_{T} [GeV/c] vs. #eta of K^{-}",500,-3.0,0.5,500,0.0,3.5);
  // ----- InvMass plots in different centrality and pT or y bins --------------
  double ptSetA[3]  = {0.4, 1.2, 2.0};
  double ptSetB[5]  = {0.4, 0.7, 1.0, 1.4, 2.0};
  double ptSetC[11] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6, 2.0, 2.5, 3.0, 4.0};

  double rapSetA[5]  = {-2.0, -1.5, -1.0, -0.5, 0};

  double centSetA[5]  = {0, 10, 40, 60, 80}; // %
  double centSetB[10]  = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80}; // %
  // Int_t cenSection[_Ncentralities]={11,22,37,57,82,113,151,174,245};//10,17,28,41,57,77,100,127,160,245 version 0 cent
  Int_t cenSection[9]={6,12,22,39,64,100,154,191,241}; // From UC Davis, cut on nFXTMult
  // pt SetA, cent SetA
  TH1D *mHist_SE_InvM_ptSetA_centSetA[2][6];
  TH1D *mHist_ME_InvM_ptSetA_centSetA[2][6];
  for(int pt=0; pt<2; pt++)
  {
    for(int cent=0; cent<6;cent++){
      mHist_SE_InvM_ptSetA_centSetA[pt][cent] = new TH1D(Form("Hist_SE_InvM_ptSetA%d_centSetA%d",pt,cent),
      Form("Hist_SE_InvM_ptSetA%d_centSetA%d",pt,cent),
      200,0.9,1.1);
      mHist_SE_InvM_ptSetA_centSetA[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_ME_InvM_ptSetA_centSetA[pt][cent] = new TH1D(Form("Hist_ME_InvM_ptSetA%d_centSetA%d",pt,cent),
      Form("Hist_ME_InvM_ptSetA%d_centSetA%d",pt,cent),
      200,0.9,1.1);
      mHist_ME_InvM_ptSetA_centSetA[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
    }
  }
  // pt SetA, cent SetB
  TH1D *mHist_SE_InvM_ptSetA_centSetB[2][9];
  TH1D *mHist_ME_InvM_ptSetA_centSetB[2][9];
  for(int pt=0; pt<2; pt++)
  {
    for(int cent=0; cent<9;cent++){
      mHist_SE_InvM_ptSetA_centSetB[pt][cent] = new TH1D(Form("Hist_SE_InvM_ptSetA%d_centSetB%d",pt,cent),
      Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetB[cent],centSetB[cent+1]),
      200,0.9,1.1);
      mHist_SE_InvM_ptSetA_centSetB[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_ME_InvM_ptSetA_centSetB[pt][cent] = new TH1D(Form("Hist_ME_InvM_ptSetA%d_centSetB%d",pt,cent),
      Form("ME, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetB[cent],centSetB[cent+1]),
      200,0.9,1.1);
      mHist_ME_InvM_ptSetA_centSetB[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
    }
  }
  // pt SetB, cent SetA
  TH1D *mHist_SE_InvM_ptSetB_centSetA[4][6];
  TH1D *mHist_ME_InvM_ptSetB_centSetA[4][6];
  for(int pt=0; pt<4; pt++)
  {
    for(int cent=0; cent<6;cent++){
      mHist_SE_InvM_ptSetB_centSetA[pt][cent] = new TH1D(Form("Hist_SE_InvM_ptSetB%d_centSetA%d",pt,cent),
      Form("Hist_SE_InvM_ptSetA%d_centSetA%d",pt,cent),
      200,0.9,1.1);
      mHist_SE_InvM_ptSetB_centSetA[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_ME_InvM_ptSetB_centSetA[pt][cent] = new TH1D(Form("Hist_ME_InvM_ptSetB%d_centSetA%d",pt,cent),
      Form("Hist_ME_InvM_ptSetA%d_centSetA%d",pt,cent),
      200,0.9,1.1);
      mHist_ME_InvM_ptSetB_centSetA[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
    }
  }
  // pt SetB, cent SetB
  TH1D *mHist_SE_InvM_ptSetB_centSetB[4][9];
  TH1D *mHist_ME_InvM_ptSetB_centSetB[4][9];
  for(int pt=0; pt<4; pt++)
  {
    for(int cent=0; cent<9;cent++){
      mHist_SE_InvM_ptSetB_centSetB[pt][cent] = new TH1D(Form("Hist_SE_InvM_ptSetB%d_centSetB%d",pt,cent),
      Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[pt],ptSetB[pt+1],centSetB[cent],centSetB[cent+1]),
      200,0.9,1.1);
      mHist_SE_InvM_ptSetB_centSetB[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_ME_InvM_ptSetB_centSetB[pt][cent] = new TH1D(Form("Hist_ME_InvM_ptSetB%d_centSetB%d",pt,cent),
      Form("ME, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[pt],ptSetB[pt+1],centSetB[cent],centSetB[cent+1]),
      200,0.9,1.1);
      mHist_ME_InvM_ptSetB_centSetB[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
    }
  }
  // pt SetC, cent 0-60%, 0-80%
  TH1D *mHist_SE_InvM_ptSetC_centAll[10][2];
  TH1D *mHist_ME_InvM_ptSetC_centAll[10][2];
  for(int pt=0; pt<10; pt++)
  {
    for(int cent=0; cent<2;cent++){
      mHist_SE_InvM_ptSetC_centAll[pt][cent] = new TH1D(Form("Hist_SE_InvM_ptSetC%d_centAll%d",pt,cent),
      Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetC[pt],ptSetC[pt+1],centSetA[0],centSetA[cent+3]),
      200,0.9,1.1);
      mHist_SE_InvM_ptSetC_centAll[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_ME_InvM_ptSetC_centAll[pt][cent] = new TH1D(Form("Hist_ME_InvM_ptSetC%d_centAll%d",pt,cent),
      Form("ME, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetC[pt],ptSetC[pt+1],centSetA[0],centSetA[cent+3]),
      200,0.9,1.1);
      mHist_ME_InvM_ptSetC_centAll[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
    }
  }
  // rap SetA, cent SetA
  TH1D *mHist_SE_InvM_rapSetA_centSetA[4][6];
  TH1D *mHist_ME_InvM_rapSetA_centSetA[4][6];
  for(int rap=0; rap<4; rap++)
  {
    for(int cent=0; cent<6;cent++){
      mHist_SE_InvM_rapSetA_centSetA[rap][cent] = new TH1D(Form("Hist_SE_InvM_rapSetA%d_centSetA%d",rap,cent),
      Form("Hist_SE_InvM_rapSetA%d_centSetA%d",rap,cent),
      200,0.9,1.1);
      mHist_SE_InvM_rapSetA_centSetA[rap][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_ME_InvM_rapSetA_centSetA[rap][cent] = new TH1D(Form("Hist_ME_InvM_rapSetA%d_centSetA%d",rap,cent),
      Form("Hist_ME_InvM_rapSetA%d_centSetA%d",rap,cent),
      200,0.9,1.1);
      mHist_ME_InvM_rapSetA_centSetA[rap][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
    }
  }
  // rap SetA, cent SetB
  TH1D *mHist_SE_InvM_rapSetA_centSetB[4][9];
  TH1D *mHist_ME_InvM_rapSetA_centSetB[4][9];
  for(int rap=0; rap<4; rap++)
  {
    for(int cent=0; cent<9;cent++){
      mHist_SE_InvM_rapSetA_centSetB[rap][cent] = new TH1D(Form("Hist_SE_InvM_rapSetA%d_centSetB%d",rap,cent),
      Form("SE, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[rap],rapSetA[rap+1],centSetB[cent],centSetB[cent+1]),
      200,0.9,1.1);
      mHist_SE_InvM_rapSetA_centSetB[rap][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_ME_InvM_rapSetA_centSetB[rap][cent] = new TH1D(Form("Hist_ME_InvM_rapSetA%d_centSetB%d",rap,cent),
      Form("ME, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[rap],rapSetA[rap+1],centSetB[cent],centSetB[cent+1]),
      200,0.9,1.1);
      mHist_ME_InvM_rapSetA_centSetB[rap][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
    }
  }
  // ============== (3) First TTree loop to get vector of events ===============
  int pre_runNumber   = -999;
  int pre_eventNumber = -999;
  // Identify same/differnt event(s)
  // vector of tracks and events
  vector<st_track> v_trk_pl;
  vector<st_track> v_trk_mi;
  vector<st_event> v_evt;
  int N_events = 0; //count # of events
  for( int i_entries = 0 /*(d_entryRange*5550)*/; i_entries< N_entries; i_entries++){
    // get kaon TTree info
    t_K->GetEntry(i_entries);
    bool    b_pos_charge     = t_K-> GetLeaf("b_pos_charge")->GetValue(0);
    bool    b_bad_TOF        = t_K-> GetLeaf("b_bad_TOF")->GetValue(0);

    int     runNumber        = t_K-> GetLeaf("runNumber")->GetValue(0);
    int     eventNumber      = t_K-> GetLeaf("eventNumber")->GetValue(0);
    int     nGoodTracks      = t_K-> GetLeaf("nGoodTracks")->GetValue(0);
    int     nFXTMult         = t_K-> GetLeaf("nFXTMult")->GetValue(0);

    double  px               = t_K-> GetLeaf("px")->GetValue(0);
    double  py               = t_K-> GetLeaf("py")->GetValue(0);
    double  pz               = t_K-> GetLeaf("pz")->GetValue(0);

    double  x_vtx            = t_K-> GetLeaf("x_vtx")->GetValue(0);
    double  y_vtx            = t_K-> GetLeaf("y_vtx")->GetValue(0);
    double  z_vtx            = t_K-> GetLeaf("z_vtx")->GetValue(0);
    double  r_vtx            = t_K-> GetLeaf("r_vtx")->GetValue(0);
    double  DCA_r            = t_K-> GetLeaf("DCA_r")->GetValue(0);

    double  tofBeta          = t_K-> GetLeaf("tofBeta")->GetValue(0);
    double  mass2            = t_K-> GetLeaf("mass2")->GetValue(0);

    double  nHitsDedx        = t_K-> GetLeaf("nHitsDedx")->GetValue(0);
    double  nHitsFit         = t_K-> GetLeaf("nHitsFit")->GetValue(0);
    double  nHitsRatio       = t_K-> GetLeaf("nHitsRatio")->GetValue(0);

    double  d_TPCnSigmaKaon  = t_K-> GetLeaf("d_TPCnSigmaKaon")->GetValue(0);
    double  d_TOFnSigmaKaon  = t_K-> GetLeaf("d_TOFnSigmaKaon")->GetValue(0);

    // QA cuts
    if((z_vtx < d_zvtxCutLvlLow) || (z_vtx > d_zvtxCutLvlHigh)) continue; // sys_cutN == 3; // vz
    if(r_vtx > d_rvtxCutLvl) continue; // sys_cutN == 4; // vr
    if(nHitsDedx <= d_nHitsDedxCutLvl) continue; // sys_cutN == 5; // dedx
    if(DCA_r >= d_DCACutLvl) continue; // sys_cutN == 6; // dca
    if(nHitsFit < d_nHitsFitCutLvl) continue; // sys_cutN == 7; // nHitsFit
    if(nHitsRatio < d_ratioCutLvl) continue; // sys_cutN == 8; // ratio
    if(d_TPCnSigmaKaon >= d_nSigmaKaonCutLvl) continue; // sys_cutN == 9; // nSigmaKaon // sys_cutN == 16; // TPCpid
    if(b_TPCpid == false && (mass2 <= d_KaonM2CutLvllow) || (mass2 >= d_KaonM2CutLvlHigh)) continue; // sys_cutN == 10; // Mass2
    if(sqrt(px*px+py*py) <= d_KaonpTCutLvlLow) continue; // sys_cutN == 11; // pTlow
    // if(d_dip_angle <= d_dipAngleCutLvl) continue; // sys_cutN == 12; // dip angle
    // if(d_vtxDiff >= d_vtxDiffCutLvl) continue; // sys_cutN == 13; // vtxDiff
    double d_inv_tofBeta     = -999.0;
    double pT  = sqrt(px*px+py*py);
    double mom = sqrt(px*px+py*py+pz*pz);
    double E   = sqrt((px*px+py*py+pz*pz)+_d_K_m*_d_K_m);
    double y   = ((E-pz) != 0.0) ? 0.5*TMath::Log( (E + pz) / (E - pz) ) : -999.0;
    double eta = ((mom - pz) != 0.0) ? 0.5*TMath::Log( (mom + pz) / (mom - pz) ) : -999.0;
    // double mass2      = mom*mom*((1.0/(tofBeta*tofBeta))-1.0);
    double d_q            = 1.;
    if(!b_pos_charge) d_q = -1.;
    double d_pq           = fabs(mom) * d_q;
    if(tofBeta != 0.0){
      d_inv_tofBeta = 1.0 / tofBeta;
      h2_TOF_beta_pq  -> Fill(d_pq,d_inv_tofBeta);
    }
    // Get physics values
    // Fill QA plots
    // h2_dEdx_pq      -> Fill(d_pq,d_dEdx);
    // h3_dEdx_pq_vs_nsig      -> Fill(d_pq,d_dEdx,d_TPCnSigmaKaon);

    // push back last event info to event vector when a new event show up
    bool b_new_event = false;
    if((pre_runNumber != runNumber)||(pre_eventNumber != eventNumber)){
      N_events++;
      b_new_event     = true;
      pre_runNumber   = runNumber;
      pre_eventNumber = eventNumber;
      st_event evt; // an event contains a vector of K+ tracks and K- tracks info
      evt.v_trk_pl = v_trk_pl;
      evt.v_trk_mi = v_trk_mi;
      v_evt.push_back(evt);
      v_trk_pl.clear(); // reset for tracks of next event
      v_trk_mi.clear(); // reset for tracks of next event
    }
    // push back K+/- tracks info
    st_track trk;
    trk.b_pos_charge    = b_pos_charge;
    trk.b_bad_TOF       = b_bad_TOF;

    trk.runNumber       = runNumber;
    trk.eventNumber     = eventNumber;
    trk.nGoodTracks     = nGoodTracks;
    trk.nFXTMult        = nFXTMult;

    trk.px              = px;
    trk.py              = py;
    trk.pz              = pz;

    trk.x_vtx           = x_vtx;
    trk.y_vtx           = y_vtx;
    trk.z_vtx           = z_vtx;
    trk.r_vtx           = r_vtx;
    trk.DCA_r           = DCA_r;

    trk.tofBeta         = tofBeta;
    trk.mass2           = mass2;

    trk.nHitsDedx       = nHitsDedx;
    trk.nHitsFit        = nHitsFit;
    trk.nHitsRatio      = nHitsRatio;

    trk.d_TPCnSigmaKaon = d_TPCnSigmaKaon;
    trk.d_TOFnSigmaKaon = d_TOFnSigmaKaon;

    if(b_pos_charge){
      v_trk_pl.push_back(trk);
      hist_KaonPlus_pT->Fill(pT);
      if(y!=-999.0){
        hist_KaonPlus_rap->Fill(y);
        hist_pt_y_kaonPlus->Fill(y,pT);
      }
      if(eta!=-999.0){
        hist_KaonPlus_eta->Fill(eta);
        hist_pt_eta_kaonPlus->Fill(eta,pT);
      }
    }
    else{
      v_trk_mi.push_back(trk);
      hist_KaonMinus_pT->Fill(pT);
      if(y!=-999.0){
        hist_KaonMinus_rap->Fill(y);
        hist_pt_y_kaonMinus->Fill(y,pT);
      }
      if(eta!=-999.0){
        hist_KaonMinus_eta->Fill(eta);
        hist_pt_eta_kaonMinus->Fill(eta,pT);
      }
    } // push back tracks and fill QA hists
  }
  cout<<endl<<" LOOPING OVER EVENTS "<< v_evt.size()<<endl<<endl;
  // ======================= (4) Same Event K+/- pairs loop ====================
  for(unsigned int i = 0; i < v_evt.size(); i++){
    st_event         evt      = v_evt[i];
    vector<st_track> v_trk_pl = v_evt[i].v_trk_pl;
    vector<st_track> v_trk_mi = v_evt[i].v_trk_mi;
    // --------------------- (4.1) Positive track loop -------------------------
    for(unsigned int j = 0; j < v_trk_pl.size(); j++){
      st_track trk0 = v_trk_pl[j]; // j-th of positive track of i-th event
      bool    b_pos_charge0    = trk0.b_pos_charge;
      bool    b_bad_TOF0       = trk0.b_bad_TOF;

      int     runNumber0       = trk0.runNumber;
      int     eventNumber0     = trk0.eventNumber;
      int     nGoodTracks0     = trk0.nGoodTracks;
      int     nFXTMult0        = trk0.nFXTMult;

      double  px0              = trk0.px;
      double  py0              = trk0.py;
      double  pz0              = trk0.pz;

      double  x_vtx0           = trk0.x_vtx;
      double  y_vtx0           = trk0.y_vtx;
      double  z_vtx0           = trk0.z_vtx;
      double  r_vtx0           = trk0.r_vtx;
      double  DCA_r0           = trk0.DCA_r;

      double  d_tofBeta0       = trk0.tofBeta;
      double  mass2_0          = trk0.mass2;

      double  d_TPCnSigmaKaon0 = trk0.d_TPCnSigmaKaon;
      double  d_TOFnSigmaKaon0 = trk0.d_TOFnSigmaKaon;
      // --------------------- (4.2) Negative track loop -----------------------
      for(unsigned int k = 0; k < v_trk_mi.size(); k++){
        st_track trk1 = v_trk_mi[k]; // k-th of positive track of i-th event
        bool    b_pos_charge1    = trk1.b_pos_charge;
        bool    b_bad_TOF1       = trk1.b_bad_TOF;

        int     runNumber1       = trk1.runNumber;
        int     eventNumber1     = trk1.eventNumber;
        int     nGoodTracks1     = trk1.nGoodTracks;
        int     nFXTMult1        = trk1.nFXTMult;

        double  px1              = trk1.px;
        double  py1              = trk1.py;
        double  pz1              = trk1.pz;

        double  x_vtx1           = trk1.x_vtx;
        double  y_vtx1           = trk1.y_vtx;
        double  z_vtx1           = trk1.z_vtx;
        double  r_vtx1           = trk1.r_vtx;
        double  DCA_r1           = trk1.DCA_r;

        double  d_tofBeta1       = trk1.tofBeta;
        double  mass2_1          = trk1.mass2;

        double  d_TPCnSigmaKaon1 = trk1.d_TPCnSigmaKaon;
        double  d_TOFnSigmaKaon1 = trk1.d_TOFnSigmaKaon;
        // K+ Variables
        double d_M0 = _d_K_m;
        double d_E0   = sqrt((px0*px0+py0*py0+pz0*pz0)+_d_K_m*_d_K_m);
        double d_mom0 = sqrt(px0*px0+py0*py0+pz0*pz0);
        double d_y0   = ((d_E0-pz0) != 0.0) ? 0.5*TMath::Log( (d_E0 + pz0) / (d_E0 - pz0) ) : -999.0;
        double eta0   = ((d_mom0 - pz0) != 0.0) ? 0.5*TMath::Log( (d_mom0 + pz0) / (d_mom0 - pz0) ) : -999.0;
        double d_pT0  = sqrt(px0*px0+py0*py0);
        double d_mT0  = sqrt(d_pT0*d_pT0 + d_M0*d_M0);
        double mass2_0 = d_mom0*d_mom0*((1.0/(d_tofBeta0*d_tofBeta0))-1.0);
        // K- Variables
        double d_M1 = _d_K_m;
        double d_E1   = sqrt((px1*px1+py1*py1+pz1*pz1)+d_M1*d_M1);
        double d_mom1 = sqrt(px1*px1+py1*py1+pz1*pz1);
        double d_y1   = ((d_E1-pz1) != 0.0) ? 0.5*TMath::Log( (d_E1 + pz1) / (d_E1 - pz1) ) : -999.0;
        double eta1   = ((d_mom1 - pz1) != 0.0) ? 0.5*TMath::Log( (d_mom1 + pz1) / (d_mom1 - pz1) ) : 1.0;
        double d_pT1  = sqrt(px1*px1+py1*py1);
        double d_mT1  = sqrt(d_pT1*d_pT1 + d_M1*d_M1);
        double mass2_1 = d_mom1*d_mom1*((1.0/(d_tofBeta1*d_tofBeta1))-1.0);
        // SE cuts
        if(runNumber0 != runNumber1) continue;
        if(eventNumber0 != eventNumber1) continue;
        if(nFXTMult0 != nFXTMult1) continue;
        // phi Variables
        double d_dip_angle = TMath::ACos((d_pT0*d_pT1+pz0*pz1) / (d_mom0*d_mom1) );
        double d_Phi_pT = sqrt(px0*px0 + py0*py0 +px1*px1 +py1+py1 + 2.*px0*px1 + 2.*py0*py1);
        double d_mT_phi = sqrt(d_Phi_pT*d_Phi_pT + _m_phi*_m_phi );
        double d_phi_pz = pz0+pz1;
        double d_phi_E  = d_E0+d_E1;
        double d_phi_y  = ((d_phi_E - d_phi_pz) != 0.0) ?  0.5*TMath::Log( (d_phi_E + d_phi_pz) / (d_phi_E - d_phi_pz) ) : -9999;
        double d_inv_m  = sqrt(  d_M0*d_M0
                              + d_M1*d_M1
                              + 2.0 *d_E0*d_E1
                              - 2.0 *(px0*px1+py0*py1+pz0*pz1) );
        hist_dip_angle         ->Fill(d_dip_angle);
        hist_SE_mass_Phi    ->Fill(d_inv_m);
        hist_SE_PhiMeson_pT ->Fill(d_Phi_pT);
        hist_SE_PhiMeson_mT ->Fill(d_mT_phi);
        hist_SE_PhiMeson_rap ->Fill(d_phi_y);
        hist_SE_pt_y_PhiMeson ->Fill(d_phi_y,d_Phi_pT);
        // invM cut
        // QA cuts
        // if((d_inv_m <= 0.9) || (d_inv_m >= 1.1)) continue;
        // if((z_vtx < d_zvtxCutLvlLow) || (z_vtx > d_zvtxCutLvlHigh)) continue; // sys_cutN == 3; // vz
        // if(r_vtx > d_rvtxCutLvl) continue; // sys_cutN == 4; // vr
        // if(nHitsDedx <= d_nHitsDedxCutLvl) continue; // sys_cutN == 5; // dedx
        // if(DCA_r >= d_DCACutLvl) continue; // sys_cutN == 6; // dca
        // if(nHitsFit < d_nHitsFitCutLvl) continue; // sys_cutN == 7; // nHitsFit
        // if(nHitsRatio < d_ratioCutLvl) continue; // sys_cutN == 8; // ratio
        // if(d_TPCnSigmaKaon >= d_nSigmaKaonCutLvl) continue; // sys_cutN == 9; // nSigmaKaon // sys_cutN == 16; // TPCpid
        // if(b_TPCpid == false && (mass2 <= d_KaonM2CutLvllow) || (mass2 >= d_KaonM2CutLvlHigh)) continue; // sys_cutN == 10; // Mass2
        // if(sqrt(px*px+py*py) <= d_KaonpTCutLvlLow) continue; // sys_cutN == 11; // pTlow
        if(d_dip_angle < d_dipAngleCutLvl) continue; // sys_cutN == 12; // dip angle
        // if(d_vtxDiff >= d_vtxDiffCutLvl) continue; // sys_cutN == 13; // vtxDiff

        // if(d_dip_angle<0.04) continue; // dip-angle cut
        // -------------------- (4.2.1) centrality def -------------------------
        // Same Event nGoodTracks1 == nGoodTracks0
        Int_t centrality = 0;
        bool a_b_cent[_Ncentralities]={false};
        bool b_pileup   = (nFXTMult1 >= 241);
        bool b_low_mult = (nFXTMult1 < 2);
        a_b_cent[0]     = (nFXTMult1 >= cenSection[7] && nFXTMult1 < cenSection[8]); // 0 - 5%, cent 1
        a_b_cent[1]     = (nFXTMult1 >= cenSection[6] && nFXTMult1 < cenSection[7]); // 5 - 10%, cent 2
        a_b_cent[2]     = (nFXTMult1 >= cenSection[5] && nFXTMult1 < cenSection[6]); // 10 - 20%, cent 3
        a_b_cent[3]     = (nFXTMult1 >= cenSection[4]  && nFXTMult1 < cenSection[5]); // 20 - 30%, cent 4
        a_b_cent[4]     = (nFXTMult1 >= cenSection[3]  && nFXTMult1 < cenSection[4]); // 30 - 40%, cent 5
        a_b_cent[5]     = (nFXTMult1 >= cenSection[2]  && nFXTMult1 < cenSection[3]); // 40 - 50%, cent 6
        a_b_cent[6]     = (nFXTMult1 >= cenSection[1]  && nFXTMult1 < cenSection[2]); // 50 - 60%, cent 7
        a_b_cent[7]     = (nFXTMult1 >= cenSection[0]  && nFXTMult1 < cenSection[1]); // 60 - 70%, cent 8
        a_b_cent[8]     = (nFXTMult1 >= 2  && nFXTMult1 < cenSection[0]); // 70 - 80%, cent 9
        for(int i=0;i<_Ncentralities;i++){
          if(a_b_cent[i]) centrality = i+1;
        }
        // -------------------- (4.2.2) Fill SE InvM plots -------------------------
        for(int pt=0; pt<2; pt++)
        {// pt SetA, cent SetA
          if(d_Phi_pT >= ptSetA[pt] && d_Phi_pT <= ptSetA[pt+1]){
            if(centrality >= 1 && centrality <= 2) mHist_SE_InvM_ptSetA_centSetA[pt][0]->Fill(d_inv_m); // 0-10%
            if(centrality >= 3 && centrality <= 5) mHist_SE_InvM_ptSetA_centSetA[pt][1]->Fill(d_inv_m); // 10-40%
            if(centrality >= 6 && centrality <= 7) mHist_SE_InvM_ptSetA_centSetA[pt][2]->Fill(d_inv_m); // 40-60%
            if(centrality >= 6 && centrality <= 9) mHist_SE_InvM_ptSetA_centSetA[pt][3]->Fill(d_inv_m); // 40-80%
            if(centrality >= 1 && centrality <= 7) mHist_SE_InvM_ptSetA_centSetA[pt][4]->Fill(d_inv_m); // 0-60%
            if(centrality >= 1 && centrality <= 9) mHist_SE_InvM_ptSetA_centSetA[pt][5]->Fill(d_inv_m); // 0-80%
          }
        }
        for(int pt=0; pt<2; pt++)
        {// pt SetA, cent SetB
          for(int cent=0; cent<9;cent++){
            if(d_Phi_pT >= ptSetA[pt] && d_Phi_pT <= ptSetA[pt+1]){
              if(centrality == cent+1 ) mHist_SE_InvM_ptSetA_centSetB[pt][cent]->Fill(d_inv_m);
            }
          }
        }
        for(int i=0; i<4; i++)
        {// pt SetB, cent SetA
          if(d_Phi_pT >= ptSetB[i] && d_Phi_pT <= ptSetB[i+1]){
            if(centrality >= 1 && centrality <= 2) mHist_SE_InvM_ptSetB_centSetA[i][0]->Fill(d_inv_m); // 0-10%
            if(centrality >= 3 && centrality <= 5) mHist_SE_InvM_ptSetB_centSetA[i][1]->Fill(d_inv_m); // 10-40%
            if(centrality >= 6 && centrality <= 7) mHist_SE_InvM_ptSetB_centSetA[i][2]->Fill(d_inv_m); // 40-60%
            if(centrality >= 6 && centrality <= 9) mHist_SE_InvM_ptSetB_centSetA[i][3]->Fill(d_inv_m); // 40-80%
            if(centrality >= 1 && centrality <= 7) mHist_SE_InvM_ptSetB_centSetA[i][4]->Fill(d_inv_m); // 0-60%
            if(centrality >= 1 && centrality <= 9) mHist_SE_InvM_ptSetB_centSetA[i][5]->Fill(d_inv_m); // 0-80%
          }

          // rap SetA, cent SetA
          if(d_phi_y >= rapSetA[i] && d_phi_y <= rapSetA[i+1]){
            if(centrality >= 1 && centrality <= 2) mHist_SE_InvM_rapSetA_centSetA[i][0]->Fill(d_inv_m); // 0-10%
            if(centrality >= 3 && centrality <= 5) mHist_SE_InvM_rapSetA_centSetA[i][1]->Fill(d_inv_m); // 10-40%
            if(centrality >= 6 && centrality <= 7) mHist_SE_InvM_rapSetA_centSetA[i][2]->Fill(d_inv_m); // 40-60%
            if(centrality >= 6 && centrality <= 9) mHist_SE_InvM_rapSetA_centSetA[i][3]->Fill(d_inv_m); // 40-80%
            if(centrality >= 1 && centrality <= 7) mHist_SE_InvM_rapSetA_centSetA[i][4]->Fill(d_inv_m); // 0-60%
            if(centrality >= 1 && centrality <= 9) mHist_SE_InvM_rapSetA_centSetA[i][5]->Fill(d_inv_m); // 0-80%
          }
          for(int cent=0; cent<9;cent++){
            if(d_phi_y >= rapSetA[i] && d_phi_y <= rapSetA[i+1]){
              if(centrality == cent+1 ) mHist_SE_InvM_rapSetA_centSetB[i][cent]->Fill(d_inv_m);
            }
          }
        }
        for(int pt=0; pt<4; pt++)
        {// pt SetB, cent SetB
          for(int cent=0; cent<9;cent++){
            if(d_Phi_pT >= ptSetB[pt] && d_Phi_pT <= ptSetB[pt+1]){
              if(centrality == cent+1 ) mHist_SE_InvM_ptSetB_centSetB[pt][cent]->Fill(d_inv_m);
            }
          }
        }
        for(int pt=0; pt<10; pt++)
        {// pt SetC, cent 0-60%, 0-80%
          for(int cent=0; cent<2;cent++){
            if(d_Phi_pT >= ptSetC[pt] && d_Phi_pT <= ptSetC[pt+1]){
              if(centrality >= 1 && centrality <= cent*2 + 7 ) mHist_SE_InvM_ptSetC_centAll[pt][cent]->Fill(d_inv_m);
            }
          }
        }

      }
    }
  }
  // ====================== (5) Mixed Event K+/- pairs loop ====================
  int j_start = 0;
  unsigned long long pre_event_pl_runnumber = 9999;
  unsigned long long this_mixed_event_pl_runnumber = 9999;
  bool b_found_mixed_evt = false;
  unsigned int i_found_mixed_evt = 0;
  bool b_half  = false;
  bool b_quart = false;
  for(int i=0;i<N_entries;i++){
    if(( i > ((double) N_entries/4.0))&&(!b_quart)&&(!b_half)) {cout<<" 0.25 done"<<endl; b_quart = true;}
    if(( i > ((double) (3.0*N_entries)/4.0))&&(!b_quart)&&(b_half)) {cout<<" 0.25 done"<<endl; b_quart = true;}
    if(( i > ((double) N_entries/2.0))&&(!b_half)) {cout<<" halfway"<<endl; b_half = true; b_quart = false;}
    t_K -> GetEntry(i);
    bool b_pos_charge0 = t_K -> GetLeaf("b_pos_charge")->GetValue(0);
    bool b_bad_TOF0    = t_K -> GetLeaf("b_pos_charge")->GetValue(0);

    unsigned int i_runNumber0   = t_K-> GetLeaf("runNumber")->GetValue(0);
    unsigned int i_eventNumber0 = t_K-> GetLeaf("eventNumber")->GetValue(0);
    unsigned int nGoodTracks0   = t_K-> GetLeaf("nGoodTracks")->GetValue(0);
    unsigned int nFXTMult0      = t_K-> GetLeaf("nFXTMult")->GetValue(0);

    double px0                  = t_K-> GetLeaf("px")->GetValue(0);
    double py0                  = t_K-> GetLeaf("py")->GetValue(0);
    double pz0                  = t_K-> GetLeaf("pz")->GetValue(0);

    double d_xvtx0              = t_K-> GetLeaf("x_vtx")->GetValue(0);
    double d_yvtx0              = t_K-> GetLeaf("y_vtx")->GetValue(0);
    double d_zvtx0              = t_K-> GetLeaf("z_vtx")->GetValue(0);
    double d_rvtx0              = t_K-> GetLeaf("r_vtx")->GetValue(0);
    double DCA_r0               = t_K-> GetLeaf("DCA_r")->GetValue(0);

    double  d_tofBeta0          = t_K-> GetLeaf("tofBeta")->GetValue(0);
    double  mass2_0             = t_K-> GetLeaf("mass2")->GetValue(0);

    double nHitsDedx0           = t_K-> GetLeaf("nHitsDedx")->GetValue(0);
    double nHitsFit0            = t_K-> GetLeaf("nHitsFit")->GetValue(0);
    double nHitsRatio0          = t_K-> GetLeaf("nHitsRatio")->GetValue(0);

    double  d_TPCnSigmaKaon0    = t_K-> GetLeaf("d_TPCnSigmaKaon")->GetValue(0);
    double  d_TOFnSigmaKaon0    = t_K-> GetLeaf("d_TOFnSigmaKaon")->GetValue(0);
    // QA cuts
    if((d_zvtx0 < d_zvtxCutLvlLow) || (d_zvtx0 > d_zvtxCutLvlHigh)) continue; // sys_cutN == 3; // vz
    if(d_rvtx0 > d_rvtxCutLvl) continue; // sys_cutN == 4; // vr
    if(nHitsDedx0 <= d_nHitsDedxCutLvl) continue; // sys_cutN == 5; // dedx
    if(DCA_r0 >= d_DCACutLvl) continue; // sys_cutN == 6; // dca
    if(nHitsFit0 < d_nHitsFitCutLvl) continue; // sys_cutN == 7; // nHitsFit
    if(nHitsRatio0 < d_ratioCutLvl) continue; // sys_cutN == 8; // ratio
    if(d_TPCnSigmaKaon0 >= d_nSigmaKaonCutLvl) continue; // sys_cutN == 9; // nSigmaKaon // sys_cutN == 16; // TPCpid
    if(b_TPCpid == false && (mass2_0 <= d_KaonM2CutLvllow) || (mass2_0 >= d_KaonM2CutLvlHigh)) continue; // sys_cutN == 10; // Mass2
    if(sqrt(px0*px0+py0*py0) <= d_KaonpTCutLvlLow) continue; // sys_cutN == 11; // pTlow
    // if(d_dip_angle <= d_dipAngleCutLvl) continue; // sys_cutN == 12; // dip angle
    // if(d_vtxDiff >= d_vtxDiffCutLvl) continue; // sys_cutN == 13; // vtxDiff
    double d_pT0    = sqrt(px0*px0+py0*py0);
    double d_mom0   = sqrt(px0*px0+py0*py0+pz0*pz0);
    double eta0     = ((d_mom0 - pz0) != 0.0) ? 0.5*TMath::Log( (d_mom0 + pz0) / (d_mom0 - pz0) ) : 1.0;
    double mass2_0  = d_mom0*d_mom0*((1.0/(d_tofBeta0*d_tofBeta0))-1.0);
    unsigned long long event_pl_runnumber0 = (i_runNumber0-19150000)*10000000+i_eventNumber0;
    if(event_pl_runnumber0 != pre_event_pl_runnumber) {
      j_start = i+1; // ME start from the next entry in kaon TTree
      pre_event_pl_runnumber = event_pl_runnumber0;
      b_found_mixed_evt = false;
    }
    // K0 track event centrality
    Int_t centrality0 = 0;
    bool a_b_cent0[_Ncentralities]={false};
    bool b_pileup0   = (nFXTMult0 > 241);
    bool b_low_mult0 = (nFXTMult0 < 2);
    a_b_cent0[0]     = (nFXTMult0 >= cenSection[7] && nFXTMult0 < cenSection[8]); // 0 - 5%, cent 1
    a_b_cent0[1]     = (nFXTMult0 >= cenSection[6] && nFXTMult0 < cenSection[7]); // 5 - 10%, cent 2
    a_b_cent0[2]     = (nFXTMult0 >= cenSection[5] && nFXTMult0 < cenSection[6]); // 10 - 20%, cent 3
    a_b_cent0[3]     = (nFXTMult0 >= cenSection[4]  && nFXTMult0 < cenSection[5]); // 20 - 30%, cent 4
    a_b_cent0[4]     = (nFXTMult0 >= cenSection[3]  && nFXTMult0 < cenSection[4]); // 30 - 40%, cent 5
    a_b_cent0[5]     = (nFXTMult0 >= cenSection[2]  && nFXTMult0 < cenSection[3]); // 40 - 50%, cent 6
    a_b_cent0[6]     = (nFXTMult0 >= cenSection[1]  && nFXTMult0 < cenSection[2]); // 50 - 60%, cent 7
    a_b_cent0[7]     = (nFXTMult0 >= cenSection[0]  && nFXTMult0 < cenSection[1]); // 60 - 70%, cent 8
    a_b_cent0[8]     = (nFXTMult0 >= 2  && nFXTMult0 < cenSection[0]); // 70 - 80%, cent 9
    for(int i=0;i<_Ncentralities;i++){
      if(a_b_cent0[i]) centrality0 = i+1;
    }
    // ====================== (5.1) Mixed Event K- loop ========================
    for(int j=j_start;j<N_entries;j++){
      t_K -> GetEntry(j);
      bool b_pos_charge1 = t_K -> GetLeaf("b_pos_charge")->GetValue(0);
      bool b_bad_TOF1    = t_K -> GetLeaf("b_pos_charge")->GetValue(0);

      unsigned int i_runNumber1   = t_K-> GetLeaf("runNumber")->GetValue(0);
      unsigned int i_eventNumber1 = t_K-> GetLeaf("eventNumber")->GetValue(0);
      unsigned int nGoodTracks1   = t_K-> GetLeaf("nGoodTracks")->GetValue(0);
      unsigned int nFXTMult1      = t_K-> GetLeaf("nFXTMult")->GetValue(0);

      double px1                  = t_K-> GetLeaf("px")->GetValue(0);
      double py1                  = t_K-> GetLeaf("py")->GetValue(0);
      double pz1                  = t_K-> GetLeaf("pz")->GetValue(0);

      double d_xvtx1              = t_K-> GetLeaf("x_vtx")->GetValue(0);
      double d_yvtx1              = t_K-> GetLeaf("y_vtx")->GetValue(0);
      double d_zvtx1              = t_K-> GetLeaf("z_vtx")->GetValue(0);
      double d_rvtx1              = t_K-> GetLeaf("r_vtx")->GetValue(0);
      double DCA_r1               = t_K-> GetLeaf("DCA_r")->GetValue(0);

      double  d_tofBeta1          = t_K-> GetLeaf("tofBeta")->GetValue(0);
      double  mass2_1             = t_K-> GetLeaf("mass2")->GetValue(0);

      double nHitsDedx1           = t_K-> GetLeaf("nHitsDedx")->GetValue(0);
      double nHitsFit1            = t_K-> GetLeaf("nHitsFit")->GetValue(0);
      double nHitsRatio1          = t_K-> GetLeaf("nHitsRatio")->GetValue(0);

      double  d_TPCnSigmaKaon1    = t_K-> GetLeaf("d_TPCnSigmaKaon")->GetValue(0);
      double  d_TOFnSigmaKaon1    = t_K-> GetLeaf("d_TOFnSigmaKaon")->GetValue(0);

      // QA cuts
      if((d_zvtx1 < d_zvtxCutLvlLow) || (d_zvtx1 > d_zvtxCutLvlHigh)) continue; // sys_cutN == 3; // vz
      if(d_rvtx1 > d_rvtxCutLvl) continue; // sys_cutN == 4; // vr
      if(nHitsDedx1 <= d_nHitsDedxCutLvl) continue; // sys_cutN == 5; // dedx
      if(DCA_r1 >= d_DCACutLvl) continue; // sys_cutN == 6; // dca
      if(nHitsFit1 < d_nHitsFitCutLvl) continue; // sys_cutN == 7; // nHitsFit
      if(nHitsRatio1 < d_ratioCutLvl) continue; // sys_cutN == 8; // ratio
      if(d_TPCnSigmaKaon1 >= d_nSigmaKaonCutLvl) continue; // sys_cutN == 9; // nSigmaKaon // sys_cutN == 16; // TPCpid
      if(b_TPCpid == false && (mass2_1 <= d_KaonM2CutLvllow) || (mass2_1 >= d_KaonM2CutLvlHigh)) continue; // sys_cutN == 10; // Mass2
      if(sqrt(px1*px1+py1*py1) <= d_KaonpTCutLvlLow) continue; // sys_cutN == 11; // pTlow
      // if(d_dip_angle <= d_dipAngleCutLvl) continue; // sys_cutN == 12; // dip angle
      // if(d_vtxDiff >= d_vtxDiffCutLvl) continue; // sys_cutN == 13; // vtxDiff

      unsigned long long event_pl_runnumber1 = (i_runNumber1-19150000)*10000000+i_eventNumber1;
      if(event_pl_runnumber0 == event_pl_runnumber1) continue; // skip SE
      if(b_found_mixed_evt && (event_pl_runnumber1 != this_mixed_event_pl_runnumber) && i_found_mixed_evt >= 5 ){
        break; // found 5 ME
      }
      if(fabs(d_zvtx1 - d_zvtx0)>d_vtxDiffCutLvl) continue; //bad zvtx
      // if(d_vtxDiff >= d_vtxDiffCutLvl) continue; // sys_cutN == 13; // vtxDiff
      if(b_pos_charge0 == b_pos_charge1) continue; // require opposite charge +-/-+ mixed events
      Int_t centrality1 = 0;
      bool a_b_cent1[_Ncentralities]={false};
      bool b_pileup1   = (nFXTMult1 > 241);
      bool b_low_mult1 = (nFXTMult1 < 2);
      a_b_cent1[0]     = (nFXTMult1 >= cenSection[7] && nFXTMult1 < cenSection[8]); // 0 - 5%, cent 1
      a_b_cent1[1]     = (nFXTMult1 >= cenSection[6] && nFXTMult1 < cenSection[7]); // 5 - 10%, cent 2
      a_b_cent1[2]     = (nFXTMult1 >= cenSection[5] && nFXTMult1 < cenSection[6]); // 10 - 20%, cent 3
      a_b_cent1[3]     = (nFXTMult1 >= cenSection[4]  && nFXTMult1 < cenSection[5]); // 20 - 30%, cent 4
      a_b_cent1[4]     = (nFXTMult1 >= cenSection[3]  && nFXTMult1 < cenSection[4]); // 30 - 40%, cent 5
      a_b_cent1[5]     = (nFXTMult1 >= cenSection[2]  && nFXTMult1 < cenSection[3]); // 40 - 50%, cent 6
      a_b_cent1[6]     = (nFXTMult1 >= cenSection[1]  && nFXTMult1 < cenSection[2]); // 50 - 60%, cent 7
      a_b_cent1[7]     = (nFXTMult1 >= cenSection[0]  && nFXTMult1 < cenSection[1]); // 60 - 70%, cent 8
      a_b_cent1[8]     = (nFXTMult1 >= 2  && nFXTMult1 < cenSection[0]); // 70 - 80%, cent 9
      for(int i=0;i<_Ncentralities;i++){
        if(a_b_cent1[i]) centrality1 = i+1;
      }
      if(centrality1!=centrality0) continue; // require same centrality of the two mixed events
      // Found ME
      b_found_mixed_evt = true;
      i_found_mixed_evt++;
      this_mixed_event_pl_runnumber = event_pl_runnumber1;
      double d_pT1   = sqrt(px1*px1+py1*py1);
      double d_mom1  = sqrt(px1*px1+py1*py1+pz1*pz1);
      double eta1    = ((d_mom1 - pz1) != 0.0) ? 0.5*TMath::Log( (d_mom1 + pz1) / (d_mom1 - pz1) ) : 1.0;
      double mass2_1 = d_mom1*d_mom1*((1.0/(d_tofBeta1*d_tofBeta1))-1.0);

      double d_M1  = _d_K_m;
      double d_M0  = _d_K_m;
      double d_E0  = sqrt((px0*px0+py0*py0+pz0*pz0)+d_M0*d_M0);
      double d_E1  = sqrt((px1*px1+py1*py1+pz1*pz1)+d_M1*d_M1);
      double d_y0  = 0.5*TMath::Log((d_E0 + pz0)/(d_E0 - pz0));
      double d_y1  = 0.5*TMath::Log((d_E1 + pz1)/(d_E1 - pz1));
      double d_mT0 = sqrt(d_pT0*d_pT0 + d_M0*d_M0);
      double d_mT1 = sqrt(d_pT1*d_pT1 + d_M1*d_M1);
      // phi physics variables
      double d_dip_angle = TMath::ACos((d_pT0*d_pT1+pz0*pz1) / (d_mom0*d_mom1) );
      double d_Phi_pT = sqrt(px0*px0 + py0*py0 +px1*px1 +py1+py1 + 2.*px0*px1 + 2.*py0*py1);
      double d_mT_phi = sqrt(d_Phi_pT*d_Phi_pT + _m_phi*_m_phi );
      double d_phi_pz = pz0+pz1;
      double d_phi_E  = d_E0+d_E1;
      double d_phi_y  = ((d_phi_E - d_phi_pz) != 0.0) ?  0.5*TMath::Log( (d_phi_E + d_phi_pz) / (d_phi_E - d_phi_pz) ) : -9999;
      double d_inv_m  = sqrt( d_M0*d_M0
                             +d_M1*d_M1
                             +2.0*d_E0*d_E1
                             -2.0*(px0*px1+py0*py1+pz0*pz1) );
      hist_ME_mass_Phi->Fill(d_inv_m);
      if((d_inv_m <= 0.9) || (d_inv_m >= 1.1)) continue;
      if(d_dip_angle < d_dipAngleCutLvl) continue; // dip-angle cut
      // if(d_dip_angle < d_dipAngleCutLvl) continue; // sys_cutN == 12; // dip angle
      // -------------------- (5.1.2) Fill ME InvM plots -------------------------
      for(int pt=0; pt<2; pt++)
      {// pt SetA, cent SetA
        if(d_Phi_pT >= ptSetA[pt] && d_Phi_pT <= ptSetA[pt+1]){
          if(centrality1 >= 1 && centrality1 <= 2) mHist_ME_InvM_ptSetA_centSetA[pt][0]->Fill(d_inv_m); // 0-10%
          if(centrality1 >= 3 && centrality1 <= 5) mHist_ME_InvM_ptSetA_centSetA[pt][1]->Fill(d_inv_m); // 10-40%
          if(centrality1 >= 6 && centrality1 <= 7) mHist_ME_InvM_ptSetA_centSetA[pt][2]->Fill(d_inv_m); // 40-60%
          if(centrality1 >= 6 && centrality1 <= 9) mHist_ME_InvM_ptSetA_centSetA[pt][3]->Fill(d_inv_m); // 40-80%
          if(centrality1 >= 1 && centrality1 <= 7) mHist_ME_InvM_ptSetA_centSetA[pt][4]->Fill(d_inv_m); // 0-60%
          if(centrality1 >= 1 && centrality1 <= 9) mHist_ME_InvM_ptSetA_centSetA[pt][5]->Fill(d_inv_m); // 0-80%
        }
      }
      for(int pt=0; pt<2; pt++)
      {// pt SetA, cent SetB
        for(int cent=0; cent<9;cent++){
          if(d_Phi_pT >= ptSetA[pt] && d_Phi_pT <= ptSetA[pt+1]){
            if(centrality1 == cent+1 ) mHist_ME_InvM_ptSetA_centSetB[pt][cent]->Fill(d_inv_m);
          }
        }
      }
      for(int i=0; i<4; i++)
      {// pt SetB, cent SetA
        if(d_Phi_pT >= ptSetB[i] && d_Phi_pT <= ptSetB[i+1]){
          if(centrality1 >= 1 && centrality1 <= 2) mHist_ME_InvM_ptSetB_centSetA[i][0]->Fill(d_inv_m); // 0-10%
          if(centrality1 >= 3 && centrality1 <= 5) mHist_ME_InvM_ptSetB_centSetA[i][1]->Fill(d_inv_m); // 10-40%
          if(centrality1 >= 6 && centrality1 <= 7) mHist_ME_InvM_ptSetB_centSetA[i][2]->Fill(d_inv_m); // 40-60%
          if(centrality1 >= 6 && centrality1 <= 9) mHist_ME_InvM_ptSetB_centSetA[i][3]->Fill(d_inv_m); // 40-80%
          if(centrality1 >= 1 && centrality1 <= 7) mHist_ME_InvM_ptSetB_centSetA[i][4]->Fill(d_inv_m); // 0-60%
          if(centrality1 >= 1 && centrality1 <= 9) mHist_ME_InvM_ptSetB_centSetA[i][5]->Fill(d_inv_m); // 0-80%
        }

        // rap SetA, cent SetA
        if(d_phi_y >= rapSetA[i] && d_phi_y <= rapSetA[i+1]){
          if(centrality1 >= 1 && centrality1 <= 2) mHist_ME_InvM_rapSetA_centSetA[i][0]->Fill(d_inv_m); // 0-10%
          if(centrality1 >= 3 && centrality1 <= 5) mHist_ME_InvM_rapSetA_centSetA[i][1]->Fill(d_inv_m); // 10-40%
          if(centrality1 >= 6 && centrality1 <= 7) mHist_ME_InvM_rapSetA_centSetA[i][2]->Fill(d_inv_m); // 40-60%
          if(centrality1 >= 6 && centrality1 <= 9) mHist_ME_InvM_rapSetA_centSetA[i][3]->Fill(d_inv_m); // 40-80%
          if(centrality1 >= 1 && centrality1 <= 7) mHist_ME_InvM_rapSetA_centSetA[i][4]->Fill(d_inv_m); // 0-60%
          if(centrality1 >= 1 && centrality1 <= 9) mHist_ME_InvM_rapSetA_centSetA[i][5]->Fill(d_inv_m); // 0-80%
        }
        for(int cent=0; cent<9;cent++){
          if(d_phi_y >= rapSetA[i] && d_phi_y <= rapSetA[i+1]){
            if(centrality1 == cent+1 ) mHist_ME_InvM_rapSetA_centSetB[i][cent]->Fill(d_inv_m);
          }
        }
      }
      for(int pt=0; pt<4; pt++)
      {// pt SetB, cent SetB
        for(int cent=0; cent<9;cent++){
          if(d_Phi_pT >= ptSetB[pt] && d_Phi_pT <= ptSetB[pt+1]){
            if(centrality1 == cent+1 ) mHist_ME_InvM_ptSetB_centSetB[pt][cent]->Fill(d_inv_m);
          }
        }
      }
      for(int pt=0; pt<10; pt++)
      {// pt SetC, cent 0-60%, 0-80%
        for(int cent=0; cent<2;cent++){
          if(d_Phi_pT >= ptSetC[pt] && d_Phi_pT <= ptSetC[pt+1]){
            if(centrality1 >= 1 && centrality1 <= cent*2 + 7 ) mHist_ME_InvM_ptSetC_centAll[pt][cent]->Fill(d_inv_m);
          }
        }
      }

    }
  }
  // -------------------------------- Set titles -------------------------------
  for(int i=0; i<4; i++)
  {// pt SetB, cent SetA
    mHist_SE_InvM_ptSetB_centSetA[i][0]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[1]));
    mHist_SE_InvM_ptSetB_centSetA[i][1]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[1],centSetA[2]));
    mHist_SE_InvM_ptSetB_centSetA[i][2]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[2],centSetA[3]));
    mHist_SE_InvM_ptSetB_centSetA[i][3]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[2],centSetA[4]));
    mHist_SE_InvM_ptSetB_centSetA[i][4]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[3]));
    mHist_SE_InvM_ptSetB_centSetA[i][5]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[4]));

    mHist_ME_InvM_ptSetB_centSetA[i][0]->SetTitle(Form("ME, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[1]));
    mHist_ME_InvM_ptSetB_centSetA[i][1]->SetTitle(Form("ME, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[1],centSetA[2]));
    mHist_ME_InvM_ptSetB_centSetA[i][2]->SetTitle(Form("ME, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[2],centSetA[3]));
    mHist_ME_InvM_ptSetB_centSetA[i][3]->SetTitle(Form("ME, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[2],centSetA[4]));
    mHist_ME_InvM_ptSetB_centSetA[i][4]->SetTitle(Form("ME, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[3]));
    mHist_ME_InvM_ptSetB_centSetA[i][5]->SetTitle(Form("ME, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[4]));
    // rap SetA, cent SetA
    mHist_SE_InvM_rapSetA_centSetA[i][0]->SetTitle(Form("SE, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[1]));
    mHist_SE_InvM_rapSetA_centSetA[i][1]->SetTitle(Form("SE, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[1],centSetA[2]));
    mHist_SE_InvM_rapSetA_centSetA[i][2]->SetTitle(Form("SE, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[2],centSetA[3]));
    mHist_SE_InvM_rapSetA_centSetA[i][3]->SetTitle(Form("SE, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[2],centSetA[4]));
    mHist_SE_InvM_rapSetA_centSetA[i][4]->SetTitle(Form("SE, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[3]));
    mHist_SE_InvM_rapSetA_centSetA[i][5]->SetTitle(Form("SE, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[4]));

    mHist_ME_InvM_rapSetA_centSetA[i][0]->SetTitle(Form("ME, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[1]));
    mHist_ME_InvM_rapSetA_centSetA[i][1]->SetTitle(Form("ME, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[1],centSetA[2]));
    mHist_ME_InvM_rapSetA_centSetA[i][2]->SetTitle(Form("ME, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[2],centSetA[3]));
    mHist_ME_InvM_rapSetA_centSetA[i][3]->SetTitle(Form("ME, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[2],centSetA[4]));
    mHist_ME_InvM_rapSetA_centSetA[i][4]->SetTitle(Form("ME, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[3]));
    mHist_ME_InvM_rapSetA_centSetA[i][5]->SetTitle(Form("ME, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[4]));
  }
  for(int pt=0; pt<2; pt++)
  {// pt SetA, cent SetA
    mHist_SE_InvM_ptSetA_centSetA[pt][0]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[1]));
    mHist_SE_InvM_ptSetA_centSetA[pt][1]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[1],centSetA[2]));
    mHist_SE_InvM_ptSetA_centSetA[pt][2]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[2],centSetA[3]));
    mHist_SE_InvM_ptSetA_centSetA[pt][3]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[2],centSetA[4]));
    mHist_SE_InvM_ptSetA_centSetA[pt][4]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[3]));
    mHist_SE_InvM_ptSetA_centSetA[pt][5]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[4]));

    mHist_ME_InvM_ptSetA_centSetA[pt][0]->SetTitle(Form("ME, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[1]));
    mHist_ME_InvM_ptSetA_centSetA[pt][1]->SetTitle(Form("ME, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[1],centSetA[2]));
    mHist_ME_InvM_ptSetA_centSetA[pt][2]->SetTitle(Form("ME, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[2],centSetA[3]));
    mHist_ME_InvM_ptSetA_centSetA[pt][3]->SetTitle(Form("ME, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[2],centSetA[4]));
    mHist_ME_InvM_ptSetA_centSetA[pt][4]->SetTitle(Form("ME, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[3]));
    mHist_ME_InvM_ptSetA_centSetA[pt][5]->SetTitle(Form("ME, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[4]));
  }

  hist_SE_mass_Phi->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
  hist_ME_mass_Phi->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
  hist_SE_PhiMeson_pT->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hist_SE_PhiMeson_mT->GetXaxis()->SetTitle("m_{T} [GeV/c^{2}]");
  hist_SE_PhiMeson_rap->GetXaxis()->SetTitle("y");
  hist_PhiMeson_tight_rap->GetXaxis()->SetTitle("y");
  hist_SE_pt_y_PhiMeson->GetXaxis()->SetTitle("y");
  hist_SE_pt_y_PhiMeson->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  hist_dip_angle->GetXaxis()->SetTitle("radian");
  hist_KaonPlus_pT->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hist_KaonMinus_pT->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hist_KaonPlus_rap->GetXaxis()->SetTitle("y");
  hist_KaonMinus_rap->GetXaxis()->SetTitle("y");
  hist_KaonPlus_eta->GetXaxis()->SetTitle("#eta");
  hist_KaonMinus_eta->GetXaxis()->SetTitle("#eta");
  hist_pt_y_kaonPlus->GetXaxis()->SetTitle("y");
  hist_pt_y_kaonPlus->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  hist_pt_y_kaonMinus->GetXaxis()->SetTitle("y");
  hist_pt_y_kaonMinus->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  h2_TOF_beta_pq->GetXaxis()->SetTitle("pq [(GeV/c)C]");
  h2_TOF_beta_pq->GetYaxis()->SetTitle("1/#beta");
  h2_dEdx_pq->GetXaxis()->SetTitle("pq [(GeV/c)C]");
  h2_dEdx_pq->GetYaxis()->SetTitle("dEdx [GeV/cm]");
  h3_dEdx_pq_vs_nsig->GetXaxis()->SetTitle("pq [(GeV/c)C]");
  h3_dEdx_pq_vs_nsig->GetYaxis()->SetTitle("dEdx [GeV/cm]");
  h3_dEdx_pq_vs_nsig->GetZaxis()->SetTitle("|nSigmaKaon|");
  tf_out->Write();
}
