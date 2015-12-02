#include "IPxDTrainingTool.h"

#include <TFile.h>
#include <iostream>

IPxDTrainingTool::IPxDTrainingTool() {

  m_runMode = 0;
}

IPxDTrainingTool::~IPxDTrainingTool() {

}

void IPxDTrainingTool::initTrainingMode(int nGrades) {

  if(m_runMode != 0) return;
  m_runMode = 1;

  m_trainings.resize(nGrades);
  for(int iGrade=0; iGrade<nGrades; iGrade++) {
    Training *t = new Training();
    t->name = TString::Format("Grade%d", iGrade);
    for(int iHypo=0; iHypo<3; iHypo++) {
      Hypothesis *h = new Hypothesis;
      int hypo;
      switch(iHypo) {
      case 0: hypo = 0; h->name = "L"; break;
      case 1: hypo = 4; h->name = "C"; break;
      case 2: hypo = 5; h->name = "B"; break;
      }
      TString suffix = t->name+"_"+h->name;
      h->h_2D = new TH1D("h_2D_"+suffix, "", 400,-40.,60.);
      h->h_2D->Sumw2();
      h->h_2D->SetDirectory(0);
      const int nbx  = 35;
      const int nby  = 20;
      double xbi[nbx+1] = {-40.0,-20.0,-15.0,-10.0,-8.0,-6.0,-5.5,-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,
			 0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,8.0,10.0,15.0,20.0,40.0,60.0};
      double ybi[nby+1] = {-40.,-20.,-12.,-8.,-4.,-3.,-2.,-1.,-0.5,0.,0.5,1.,2.,3.,4.,6.,8.,12.,20.,40.,60.};
      h->h_3D = new TH2D("h_3D_"+suffix, "", nbx, xbi, nby, ybi);
      h->h_3D->Sumw2();
      h->h_3D->SetDirectory(0);
      t->h[hypo] = h;
    }
    m_trainings[iGrade] = t;
  }
}

void IPxDTrainingTool::initEvaluationModeDB(TString fname, TString jetCollection, std::vector<TString> grades) {

//  if(m_runMode != 0) return;
  m_runMode = 2;

  TFile *F = new TFile(fname, "read");
  if(!F || F->IsZombie()) {
    m_runMode = 0;
    std::cout << "Error: could not open file " << fname << " for reading" << std::endl;
    return;
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
      case 1: hypo = 4; h->name = "C"; break;
      case 2: hypo = 5; h->name = "B"; break;
      }
      h->h_2D = (TH1*)F->Get("IP2D/"+jetCollection+"/"+h->name+"/"+t->name+"/SipA0")->Clone();
      h->h_3D = (TH2*)F->Get("IP3D/"+jetCollection+"/"+h->name+"/"+t->name+"/Sip3D")->Clone();
      t->h[hypo] = h;
    }
    m_trainings[iGrade] = t;
  }

  F->Clear();
  F->Close();
  delete F;
}

void IPxDTrainingTool::initEvaluationMode(int nGrades, TString fname) {

  if(m_runMode != 0) return;
  m_runMode = 2;

  TFile *F = new TFile(fname, "read");
  if(!F || F->IsZombie()) {
    m_runMode = 0;
    std::cout << "Error: could not open file " << fname << " for reading" << std::endl;
    return;
  }

  m_trainings.resize(nGrades);
  for(int iGrade=0; iGrade<nGrades; iGrade++) {
    Training *t = new Training();
    t->name = TString::Format("Grade%d", iGrade);
    for(int iHypo=0; iHypo<3; iHypo++) {
      Hypothesis *h = new Hypothesis;
      int hypo;
      switch(iHypo) {
      case 0: hypo = 0; h->name = "L"; break;
      case 1: hypo = 4; h->name = "C"; break;
      case 2: hypo = 5; h->name = "B"; break;
      }
      TString suffix = t->name+"_"+h->name;
      h->h_2D = (TH1*)F->Get("h_2D_"+suffix)->Clone();
      h->h_3D = (TH2*)F->Get("h_3D_"+suffix)->Clone();
      t->h[hypo] = h;
    }
    m_trainings[iGrade] = t;
  }

  F->Clear();
  F->Close();
  delete F;
}

void IPxDTrainingTool::fillTrack(int grade, float d0sig, float z0sig, int hypothesis) {

  if(m_runMode != 1) return;
  if(grade < 0 || grade >= (int)m_trainings.size()) return;
  if(hypothesis != 0 && hypothesis != 4 && hypothesis != 5) return;

  Training *t = m_trainings[grade];
  Hypothesis *h = t->h[hypothesis];

  h->h_2D->Fill(d0sig);
  h->h_3D->Fill(d0sig, z0sig);
}

void IPxDTrainingTool::write(TString fname) {

  if(m_runMode != 1) return;

  TFile *F = new TFile(fname, "recreate");
  if(!F || F->IsZombie()) {
    m_runMode = 0;
    std::cout << "Error: could not open file " << fname << " for writing" << std::endl;
    return;
  }

  for(unsigned int iGrade=0; iGrade<m_trainings.size(); iGrade++) {
    Training *t = m_trainings[iGrade];
    for(int iHypo=0; iHypo<3; iHypo++) {
      int hypo;
      switch(iHypo) {
      case 0: hypo = 0; break;
      case 1: hypo = 4; break;
      case 2: hypo = 5; break;
      }
      Hypothesis *h = t->h[hypo];
      h->h_2D->Write();
      h->h_3D->Write();
    }
  }

  F->Close();
  delete F;
}


float IPxDTrainingTool::getTrackProb(int nD, int grade, float d0sig, float z0sig, int hypothesis) {

  if(m_runMode != 2) return -999;
  if(grade < 0 || grade >= (int)m_trainings.size()) return -999;
  if(hypothesis != 0 && hypothesis != 4 && hypothesis != 5) return -999;
 	
	if(nD < 2 || nD > 3) return -999;

  Training *t = m_trainings[grade];
  Hypothesis *h = t->h[hypothesis];

  float prob = -999;
  if(nD == 2) {
    prob = h->h_2D->GetBinContent(h->h_2D->FindBin(d0sig)) /
      h->h_2D->Integral(0, h->h_2D->GetNbinsX()+1);
  }
  if(nD == 3) {
    prob = h->h_3D->GetBinContent(h->h_3D->FindBin(d0sig, z0sig)) /
      h->h_3D->Integral(0, h->h_3D->GetNbinsX()+1, 0, h->h_3D->GetNbinsY()+1);
  }
  return prob;
}

std::vector<float> IPxDTrainingTool::getTrackProb(int nD, int grade, float d0sig, float z0sig) {

  std::vector<float> probs;
  if(m_runMode != 2) return probs;

  for(int hypo=0; hypo<3; hypo++) {
    probs.push_back(getTrackProb(nD, grade, d0sig, z0sig, hypo));
  }

  return probs;
}

float IPxDTrainingTool::getTrackLLR(int nD, int grade, float d0sig, float z0sig, int hypoS, int hypoB) {

  if(m_runMode != 2) return -999;
  float pS = getTrackProb(nD, grade, d0sig, z0sig, hypoS);
  float pB = getTrackProb(nD, grade, d0sig, z0sig, hypoB);

  float llR = -999;
  if(pB > 0 && pS > 0) {
    llR = log(pS/pB);
  }

  return llR;
}

float IPxDTrainingTool::getJetLLR(int nD, std::vector<int> grade, std::vector<float> d0sig, std::vector<float> z0sig, int hypoS, int hypoB) {
	m_runMode = 2;
  float llR = 0;
  for(unsigned int i=0; i<grade.size(); i++) {
    float illR = getTrackLLR(nD, grade[i], d0sig[i], z0sig[i], hypoS, hypoB);
    if(illR != -999) {
     	llR += illR;
    }
  }
  return llR;
}


