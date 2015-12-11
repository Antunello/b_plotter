#include "IPxDStandaloneTool.h"

#include <TH2.h>
#include <TFile.h>
#include <TROOT.h>

#include <iostream>

#include "Utils.h"

IPxDStandaloneTool::IPxDStandaloneTool() {

  m_runMode = 0;
}

IPxDStandaloneTool::~IPxDStandaloneTool() {

}

int IPxDStandaloneTool::initTrainingMode(std::string fname, std::string jetCollection, std::vector<std::string> grades) { 

  if(m_runMode != 0) return 1;
  m_runMode = 1;

  TDirectory *gTemp = gDirectory;

  m_file = TFile::Open(fname.c_str(), "RECREATE");
  if(!m_file || m_file->IsZombie()) {
    std::cout << "Error in IPxDStandaloneTool::initTrainingMode(): could not open file " << fname << std::endl;
    return 1;
  }

  int nGrades = grades.size();
  m_trainings.resize(nGrades);
  for(int iGrade=0; iGrade<nGrades; iGrade++) {
    Training *t = new Training();
    t->name = grades[iGrade];
    for(int iHypo=0; iHypo<3; iHypo++) {
      Hypothesis *h = new Hypothesis;
      int hypo;
      switch(iHypo) {
      case 0: hypo = 0; h->name = "U"; break;
      //case 0: hypo = 0; h->name = "L"; break;
      case 1: hypo = 4; h->name = "C"; break;
      case 2: hypo = 5; h->name = "B"; break;
      }
      m_file->mkdir(("IP2D/"+jetCollection+"/"+h->name+"/"+t->name).c_str());
      m_file->cd(("IP2D/"+jetCollection+"/"+h->name+"/"+t->name).c_str());
      h->h_2D = new TH1D("SipA0", "", 400, -40., 60.);
      h->h_2D->Sumw2();
      h->h_2D->SetDirectory(gDirectory);
      const int nbx  = 35;
      const int nby  = 20;
      double xbi[nbx+1] = {-40.0,-20.0,-15.0,-10.0,-8.0,-6.0,-5.5,-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,
			 0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,8.0,10.0,15.0,20.0,40.0,60.0};
      double ybi[nby+1] = {-40.,-20.,-12.,-8.,-4.,-3.,-2.,-1.,-0.5,0.,0.5,1.,2.,3.,4.,6.,8.,12.,20.,40.,60.};
      m_file->mkdir(("IP3D/"+jetCollection+"/"+h->name+"/"+t->name).c_str());
      m_file->cd(("IP3D/"+jetCollection+"/"+h->name+"/"+t->name).c_str());
      h->h_3D = new TH2D("Sip3D", "", nbx, xbi, nby, ybi);
      h->h_3D->Sumw2();
      h->h_3D->SetDirectory(gDirectory);
      t->h[hypo] = h;
    }
    m_trainings[iGrade] = t;
  }

  gTemp->cd();

  return 0;
}

int IPxDStandaloneTool::initEvaluationMode(std::string fname, std::string jetCollection, std::vector<std::string> grades) {

  if(m_runMode != 0) return 1;
  m_runMode = 2;

  m_file = TFile::Open(fname.c_str(), "READ");
  if(!m_file || m_file->IsZombie()) {
    std::cout << "Error in IPxDStandaloneTool::initEvaluationMode(): could not open file " << fname << " for reading" << std::endl;
    return 1;
  }

  int nGrades = grades.size();
  m_trainings.resize(nGrades);
  for(int iGrade=0; iGrade<nGrades; iGrade++) {
    Training *t = new Training();
    t->name = grades[iGrade];
    for(int iHypo=0; iHypo<3; iHypo++) {
      Hypothesis *h = new Hypothesis;
      int hypo;
      switch(iHypo) {
      case 0: hypo = 0; h->name = "U"; break;
      //case 0: hypo = 0; h->name = "L"; break;
      case 1: hypo = 4; h->name = "C"; break;
      case 2: hypo = 5; h->name = "B"; break;
      }
      h->h_2D = Utils::readObject<TH1>(m_file, "IP2D/"+jetCollection+"/"+h->name+"/"+t->name+"/SipA0");
      h->h_3D = Utils::readObject<TH2>(m_file, "IP3D/"+jetCollection+"/"+h->name+"/"+t->name+"/Sip3D");
      h->I_2D = h->h_2D->Integral(0, h->h_2D->GetNbinsX()+1);
      h->I_3D = h->h_3D->Integral(0, h->h_3D->GetNbinsX()+1, 0, h->h_3D->GetNbinsY()+1);
      t->h[hypo] = h;
    }
    m_trainings[iGrade] = t;
  }

  m_file->Clear();
  m_file->Close();
  delete m_file;

  return 0;
}


void IPxDStandaloneTool::fillTrack(int grade, float d0sig, float z0sig, int hypothesis) {

  fillTrack(2, grade, d0sig, z0sig, hypothesis);
  fillTrack(3, grade, d0sig, z0sig, hypothesis);
}

void IPxDStandaloneTool::fillTrack(int nD, int grade, float d0sig, float z0sig, int hypothesis) {

  if(m_runMode != 1) return;
  if(grade < 0 || grade >= (int)m_trainings.size()) return;
  if(hypothesis != 0 && hypothesis != 4 && hypothesis != 5) return;
  if(nD < 2 || nD > 3) return;

  Training *t = m_trainings[grade];
  Hypothesis *h = t->h[hypothesis];

  if(nD == 2) h->h_2D->Fill(d0sig);
  if(nD == 3) h->h_3D->Fill(d0sig, z0sig);
}

void IPxDStandaloneTool::write() {

  if(m_runMode != 1) return;

  m_file->Write();
  m_file->Close();
  delete m_file;
}


float IPxDStandaloneTool::getTrackProb(int nD, int grade, float d0sig, float z0sig, int hypothesis) {

  if(m_runMode != 2) return -999;
  if(grade < 0 || grade >= (int)m_trainings.size()) return -999;
  if(hypothesis != 0 && hypothesis != 4 && hypothesis != 5) return -999;
  if(nD < 2 || nD > 3) return -999;

  Training *t = m_trainings[grade];
  Hypothesis *h = t->h[hypothesis];

  float prob = -999;
  if(nD == 2) {
    prob = h->h_2D->GetBinContent(h->h_2D->FindBin(d0sig)) / h->I_2D;
  }
  if(nD == 3) {
    prob = h->h_3D->GetBinContent(h->h_3D->FindBin(d0sig, z0sig)) / h->I_3D;
  }

  return prob;
}

float IPxDStandaloneTool::getTrackLLR(int nD, int grade, float d0sig, float z0sig, int hypoS, int hypoB) {

  if(m_runMode != 2) return -999;

  float pS = getTrackProb(nD, grade, d0sig, z0sig, hypoS);
  float pB = getTrackProb(nD, grade, d0sig, z0sig, hypoB);

  float llR = -999;
  if(pB > 0 && pS > 0) {
    llR = log(pS/pB);
  }

  return llR;
}

float IPxDStandaloneTool::getJetLLR(int nD, std::vector<int> grade, std::vector<float> d0sig, std::vector<float> z0sig, int hypoS, int hypoB) {

  float llR = 0;

  for(unsigned int i=0; i<grade.size(); i++) {
    float illR = getTrackLLR(nD, grade[i], d0sig[i], z0sig[i], hypoS, hypoB);
    if(illR != -999) {
      llR += illR;
    }
  }

  return llR;
}

float IPxDStandaloneTool::getJetProb(int nD, std::vector<int> grade, std::vector<float> d0sig, std::vector<float> z0sig, int hypothesis) {

  float prob = 1;

  for(unsigned int i=0; i<grade.size(); i++) {
    float iProb = getTrackProb(nD, grade[i], d0sig[i], z0sig[i], hypothesis);
    if(prob != -999) {
      prob *= iProb;
    }
  }

  return prob;
}


