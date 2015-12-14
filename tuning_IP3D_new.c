#include <string>
#include <iostream>
#include <sstream>
#include <dirent.h>
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TApplication.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TSystem.h"
#include <algorithm>
#include "TStyle.h"
#include "TColor.h"
#include "TProfile.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "IPxDStandaloneTool.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include <math.h>

std::vector<std::string> get_track_grades(){

	std::vector<std::string> grades;
	
	grades.push_back("0HitIn0HitNInExp2");
	grades.push_back("0HitIn0HitNInExpIn");
	grades.push_back("0HitIn0HitNInExpNIn");
	grades.push_back("0HitIn0HitNIn");
	grades.push_back("0HitInExp");
	grades.push_back("0HitIn");
	grades.push_back("0HitNInExp");
	grades.push_back("0HitNIn");
	grades.push_back("InANDNInShared");
	grades.push_back("PixShared");
	grades.push_back("SctShared");
	grades.push_back("InANDNInSplit");
	grades.push_back("PixSplit");
	grades.push_back("Good");

	return grades;
}



void initBranches(TChain* myChain){

	myChain->SetBranchStatus("*",0);
	myChain->SetBranchStatus("bH_x",1);
	myChain->SetBranchStatus("bH_y",1);
	myChain->SetBranchStatus("bH_z",1);
	myChain->SetBranchStatus("bH_pt",1);
	myChain->SetBranchStatus("bH_phi",1);
	myChain->SetBranchStatus("bH_eta",1);
	myChain->SetBranchStatus("njets",1);
	myChain->SetBranchStatus("jet_pt",1);
	myChain->SetBranchStatus("jet_E",1);
	myChain->SetBranchStatus("jet_truthMatch",1);
	myChain->SetBranchStatus("jet_eta",1);
	myChain->SetBranchStatus("jet_phi",1);
	myChain->SetBranchStatus("jet_truthflav",1);
	myChain->SetBranchStatus("jet_LabDr_HadF",1);
	myChain->SetBranchStatus("jet_trk_ip3d_llr",1);
	myChain->SetBranchStatus("jet_ip3d_llr",1);
	myChain->SetBranchStatus("jet_trk_vtx_X",1);
	myChain->SetBranchStatus("jet_trk_vtx_Y",1);
	myChain->SetBranchStatus("jet_trk_d0",1);
	myChain->SetBranchStatus("jet_trk_z0",1);
	myChain->SetBranchStatus("jet_trk_ip3d_d0",1);
	myChain->SetBranchStatus("jet_trk_ip3d_z0",1);
	myChain->SetBranchStatus("jet_trk_d0_truth",1);
	myChain->SetBranchStatus("jet_trk_z0_truth",1);
	myChain->SetBranchStatus("jet_trk_ip3d_d0sig",1);
	myChain->SetBranchStatus("jet_trk_ip3d_z0sig",1);
	myChain->SetBranchStatus("jet_trk_pt",1);
	myChain->SetBranchStatus("jet_trk_pt",1);
	myChain->SetBranchStatus("jet_trk_eta",1);
	myChain->SetBranchStatus("jet_trk_theta",1);
	myChain->SetBranchStatus("jet_trk_phi",1);
	myChain->SetBranchStatus("jet_trk_algo",1);
	myChain->SetBranchStatus("jet_jf_nvtx",1);
 	myChain->SetBranchStatus("bH_Lxy",1);
	myChain->SetBranchStatus("jet_jf_ntrkAtVx",1);
	myChain->SetBranchStatus("jet_sv1_ntrkv",1);
	myChain->SetBranchStatus("jet_sv1_Nvtx",1);
	myChain->SetBranchStatus("jet_ip3d_ntrk",1);
	myChain->SetBranchStatus("jet_trk_ip3d_grade",1);
	myChain->SetBranchStatus("jet_trk_nInnHits",1);
	myChain->SetBranchStatus("jet_trk_orig",1);
	myChain->SetBranchStatus("jet_trk_nNextToInnHits",1);
	myChain->SetBranchStatus("bH_nBtracks",1);
	myChain->SetBranchStatus("bH_nCtracks",1);
	myChain->SetBranchStatus("bH_dRjet",1);
	myChain->SetBranchStatus("jet_aliveAfterOR",1);
}


void tuning_IP3D(std::string inputFolder, double n_cut){

	std::string chain_name = "bTag_AntiKt4EMTopoJets";
	TChain* myChain = new TChain(chain_name.c_str());
	float eta_cut = 2.5;
	DIR* dir;
	dirent* pdir;
	dir = opendir(inputFolder.c_str());
	while (pdir = readdir(dir)){
		std::string foldName = pdir->d_name;
		if(foldName.find("mc")==std::string::npos) continue;
		//cout << pdir->d_name << endl;
		DIR* dir2;
		dirent* pdir2;
		dir2 = opendir((inputFolder+"/"+foldName).c_str());
		while (pdir2 = readdir(dir2)){
			std::string fName=pdir2->d_name;
			if(fName.find("root")==std::string::npos) continue;
			myChain->Add( (inputFolder+"/"+foldName+"/"+fName).c_str() );
		}	

	}	

	std::vector<std::string> grades = get_track_grades();

	IPxDStandaloneTool *tool = new IPxDStandaloneTool();
//	tool->initTrainingMode(14);
	tool->initTrainingMode("ip3d_tuning_new.root","AntiKt4EMTopo" ,grades);

	std::cout<<"Chain Entries:"<<myChain->GetEntries()<<std::endl;
	initBranches(myChain);
	
	std::vector<float> *jet_pt = 0;
	std::vector<int> *jet_truthMatch = 0;
	std::vector<int> *jet_truthflav = 0;
	std::vector<int> *jet_LabDr_HadF = 0;
	Int_t njets = 0;
	Int_t eventnb = 0;
	std::vector<std::vector<float> > *jet_trk_ip3d_d0sig = 0;
	std::vector<std::vector<float> > *jet_trk_pt = 0;
	std::vector<std::vector<int> > *jet_trk_ip3d_grade = 0;
	std::vector<std::vector<float> > *jet_trk_ip3d_z0sig = 0;
	std::vector<int> *jet_btag_ntrk = 0;
	std::vector<int> *jet_aliveAfterOR =0;   

	TBranch *b_njets, *b_eventnb, *b_jet_pt,*b_jet_trk_pt, *b_jet_truthMatch,*b_jet_truthflav,*b_jet_LabDr_HadF, *b_jet_trk_ip3d_d0sig, *b_jet_trk_ip3d_z0sig, *b_jet_trk_ip3d_grade, *b_jet_btag_ntrk, *b_jet_aliveAfterOR;

	myChain->SetBranchAddress("njets", &njets, &b_njets);
	myChain->SetBranchAddress("eventnb", &eventnb, &b_eventnb);
  myChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
  myChain->SetBranchAddress("jet_trk_pt", &jet_trk_pt, &b_jet_trk_pt);
	myChain->SetBranchAddress("jet_truthMatch", &jet_truthMatch, &b_jet_truthMatch);
	myChain->SetBranchAddress("jet_truthflav", &jet_truthflav, &b_jet_truthflav);
	myChain->SetBranchAddress("jet_LabDr_HadF", &jet_LabDr_HadF, &b_jet_LabDr_HadF);
	myChain->SetBranchAddress("jet_trk_ip3d_z0sig", &jet_trk_ip3d_z0sig, &b_jet_trk_ip3d_z0sig);
	myChain->SetBranchAddress("jet_trk_ip3d_d0sig", &jet_trk_ip3d_d0sig, &b_jet_trk_ip3d_d0sig);
	myChain->SetBranchAddress("jet_trk_ip3d_grade", &jet_trk_ip3d_grade, &b_jet_trk_ip3d_grade);
  myChain->SetBranchAddress("jet_btag_ntrk", &jet_btag_ntrk, &b_jet_btag_ntrk);
	myChain->SetBranchAddress("jet_aliveAfterOR", &jet_aliveAfterOR, &b_jet_aliveAfterOR);

	Long64_t nentries = myChain->GetEntriesFast();
	for(Long64_t jentry = 0; jentry<nentries; jentry++){

		Long64_t ientry = myChain->LoadTree(jentry);
		if(ientry<0) break;
		if(ientry%1000 == 0)std::cout<<ientry<<std::endl;
		b_njets->GetEntry(ientry);
		b_eventnb->GetEntry(ientry);
		for(int i = 0; i<njets; i++){
			b_jet_pt->GetEntry(ientry);
			b_jet_trk_pt->GetEntry(ientry);
			b_jet_truthMatch->GetEntry(ientry);
			b_jet_aliveAfterOR->GetEntry(ientry);
  		if (jet_truthMatch  ->at(i)!=1) continue;
  		if (jet_aliveAfterOR->at(i)!=1) continue;
  		if (jet_pt->at(i)<25e3)         continue;
			//b_jet_truthflav->GetEntry(ientry);
			b_jet_LabDr_HadF->GetEntry(ientry);
			b_jet_trk_ip3d_d0sig->GetEntry(ientry);
			b_jet_trk_ip3d_z0sig->GetEntry(ientry);
			b_jet_trk_ip3d_grade->GetEntry(ientry);
			b_jet_btag_ntrk->GetEntry(ientry);

			Int_t idx_trkdz0sig[jet_trk_ip3d_d0sig->at(i).size()];
			float trk_dz0sig_vec[jet_trk_ip3d_d0sig->at(i).size()];
			for(int l=0; l<jet_trk_ip3d_d0sig->at(i).size();l++){
				if(jet_trk_ip3d_d0sig->at(i)[l] == -999)trk_dz0sig_vec[l] = -999;
				else trk_dz0sig_vec[l] =/*(jet_trk_ip3d_z0sig->at(i)[l]/abs(jet_trk_ip3d_z0sig->at(i)[l]))*/jet_trk_pt->at(i)[l]* sqrt((pow(jet_trk_ip3d_d0sig->at(i)[l],2) + pow(jet_trk_ip3d_z0sig->at(i)[l],2)));
			}
			TMath::Sort(int (jet_trk_ip3d_d0sig->at(i).size()), trk_dz0sig_vec, idx_trkdz0sig);
	
			Int_t idx_trkpt[jet_trk_ip3d_d0sig->at(i).size()];
			TMath::Sort(int (jet_trk_ip3d_d0sig->at(i).size()), &(jet_trk_pt->at(i)[0]), idx_trkpt);
//			double r_cut = ( 0.028 * jet_pt->at(i) - 5.46001e-8* pow(jet_pt->at(i),2) + 2.85982e-14 * pow(jet_pt->at(i),3) + 5.69244e-21 * pow(jet_pt->at(i),4) + 1.75064e-28*pow(jet_pt->at(i),5) );
			//double r_cut = 0.02*jet_pt->at(i); // ( 0.0553861 * jet_pt->at(i) - 5.46001e-8* pow(jet_pt->at(i),2) + 2.85982e-14 * pow(jet_pt->at(i),3) + 5.69244e-21 * pow(jet_pt->at(i),4) + 1.75064e-28*pow(jet_pt->at(i),5) );
	

//			double r_cut =		0.0710166e6 - 7.4378e-1 * pow(jet_pt->at(i),1) + 3.3711e-6 * pow(jet_pt->at(i),2) - 7.68983e-12 * pow(jet_pt->at(i),3) + 9.94882e-18 * pow(jet_pt->at(i),4) - 7.80337e-24 * pow(jet_pt->at(i),5) +	3.78802e-30 * pow(jet_pt->at(i),6) - 1.11231e-36 * pow(jet_pt->at(i),7) +	1.8105e-43* pow(jet_pt->at(i),8) -1 .25372e-50 * pow(jet_pt->at(i),9);
			//double n_cut=1.0;
			double r_cut =	0.5*(	0.0710166 +
										-7.4378e-1 * pow(jet_pt->at(i)*1e-6,1) +
										3.3711 * pow(jet_pt->at(i)*1e-6,2) +
										-7.68983 * pow(jet_pt->at(i)*1e-6,3) +
										9.94882 * pow(jet_pt->at(i)*1e-6,4) +
										-7.80337 * pow(jet_pt->at(i)*1e-6,5) +
										3.78802 * pow(jet_pt->at(i)*1e-6,6) +
										-1.11231 * pow(jet_pt->at(i)*1e-6,7) +
										1.8105e-1 * pow(jet_pt->at(i)*1e-6,8) +
										-1.25372e-2 * pow(jet_pt->at(i)*1e-6,9));
			for(int ok=0; ok<jet_btag_ntrk->at(i); ok++){
				int k = idx_trkdz0sig[ok];
				//int k = idx_trkpt[ok];
				/*if(jet_trk_pt->at(i)[k] > r_cut *jet_pt->at(i))*/
				if(ok < n_cut )tool->fillTrack(jet_trk_ip3d_grade->at(i)[k],jet_trk_ip3d_d0sig->at(i)[k],jet_trk_ip3d_z0sig->at(i)[k], jet_LabDr_HadF->at(i));

			}

		}

	}
	tool->write();

}
