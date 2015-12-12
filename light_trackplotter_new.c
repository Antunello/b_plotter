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
#include "THStack.h"
#include "TLorentzVector.h"
#include "IPxDStandaloneTool.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TrackPlotting/TrackPlotting.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

//#include "AtlasStyle.C"
//#include "AtlasUtils.C"

#define BJETS 5
#define CJETS 4
#define LJETS 0


#define TRKORIG_FAKE = -1
#define TRKORIG_B = 0
#define TRKORIG_C = 1
#define TRKORIG_FRAG = 2
#define TRKORIG_GEANT = 3

static const float MV1_CUT = 0.945487725; 
static const float MV1c_CUT = 0.779833333333;
static const float MV2c00_CUT = 0.0308333333333;
static const float MV2c10_CUT = -0.00416666666667;
static const float MV2c20_CUT = -0.0215;
static const float IP3D_CUT = 2.007;
static const float IP3DSV1_CUT = 4.3625;
static const float MVb_CUT = -0.120991666667;
static const float SV1_CUT = 0.354833333333;
static const float JF_CUT = -1.6125;


struct track_data{
	int trk_orig;
	double trk_pt;
	int trk_grade;
};


struct by_pt { 
    bool operator()(track_data const &a, track_data const &b) { 
        return a.trk_pt > b.trk_pt;
    }
};

void FicoPlot();

TH1D* ratioplot(TH2D* map_2d, std::string hname);
void FicoPlot(){

	const int NRGBs=5;
	const int NCont =255;
	double stops[NRGBs] = {0.00,0.34,0.61,1.00,1.00};
	double red[NRGBs] = {0.00,0.00,0.87,1.00,1.00};
	double green[NRGBs] = {0.00,1.00,0.80,0.20,0.00};
	double blue[NRGBs] ={1.00,1.00,0.20,0.00,0.00};
 
	TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,NCont);
	gStyle->SetNumberContours(NCont);
	gPad->SetLogz(0);
	gStyle->SetOptStat(0);
//	gPad->SetRightMargin(0.15);
//	  gStyle->SetPaintTextFormat("4.0f");
//
//	  }
//
}


TH1D* ratioplot(TH2D* map_2d, std::string hname){


	TH1D* ratio = new TH1D(hname.c_str(),hname.c_str(), map_2d->GetYaxis()->GetNbins(), 0, map_2d->GetYaxis()->GetBinCenter(map_2d->GetYaxis()->GetNbins() -1) );
	for(int y = 0; y<map_2d->GetYaxis()->GetNbins(); y++){
		int sd0_c=0, pd0_c = 0;
		for(int x = 0; x<map_2d->GetXaxis()->GetNbins(); x++){
			if(map_2d->GetXaxis()->GetBinCenter(x+1)>0) pd0_c += map_2d->GetBinContent(x+1,y+1);
			else sd0_c += map_2d->GetBinContent(x+1,y+1);
		}
		if(sd0_c != 0 && pd0_c != 0)ratio->SetBinContent(y+1,(double)pd0_c/(double)(sd0_c+pd0_c));
	}
	ratio->Rebin(4);
	ratio->Scale(0.25);
	return ratio;
}

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

struct ip3d_llr_ptbin{
	std::string ID;
	int nbin;
	double pt_limit;
	double ip3d_cut;
	double seventy_ip3d_cut;
	std::vector<double> v_flat_cut;

	std::vector<double> pt_bins, pt_bin_width;
	std::vector<TH1D*> b_ip3d_llr_ptbin;
	std::vector<TH1D*> l_ip3d_llr_ptbin;
	TH1D* b_eff;
	TH1D* b_eff_flat;
	TH1D* l_rej;
	TH1D* l_rej_flat;
	

	void Init(int nb, double pt_l, std::string id);
	void InitHistos();
	void set_flat_cut();

};

void ip3d_llr_ptbin::Init(int nb, double pt_l, std::string id){
	nbin = nb;
	pt_limit = pt_l;
	ID = id;

	for(int i=0; i<nb+1; i++){
		double ptbin =( (double)(i) )*(pt_limit / ( (double)(nb) ) );
		pt_bins.push_back( ptbin );
		pt_bin_width.push_back( pt_limit / ( (double)(nb)) );
	}
}

void ip3d_llr_ptbin::InitHistos(){
	//int nbin = 4;	

	//double pt_limit = 3e6;
	double pt_step = pt_limit/((double)nbin);

	b_eff = new TH1D((ID+"_b_eff").c_str(),(ID+"_b_eff;jet_pt (MeV);b_eff").c_str(), nbin,0,pt_limit);
	b_eff_flat = new TH1D((ID+"_b_eff_flat").c_str(),(ID+"_b_eff_flat;jet_pt (MeV);b_eff_flat").c_str(), nbin,0,pt_limit);
	l_rej = new TH1D((ID+"_l_rej").c_str(),(ID+"_l_rej;jet_pt (MeV);l_rej").c_str(), nbin,0,pt_limit);
	l_rej_flat = new TH1D((ID+"_l_rej_flat").c_str(),(ID+"_l_rej_flat;jet_pt (MeV);l_rej_flat").c_str(), nbin,0,pt_limit);


	for(int i=0; i<nbin; i++){

		//float flat_cut = 30;	
		std::stringstream ssb;
		ssb.str(std::string());
		ssb<<ID<<"_b_IP3D_LLR_ptbin"<<i<<"_minpt_"<<i*pt_step<<"_maxpt_"<<(i+1)*pt_step;

		std::stringstream ssl;
		ssl.str(std::string());
		ssl<<ID<<"_l_IP3D_LLR_ptbin"<<i<<"_minpt_"<<i*pt_step<<"_maxpt_"<<(i+1)*pt_step;

		TH1D *b_hist = new TH1D(ssb.str().c_str(),(ssb.str()+";ip3d_llr;entries").c_str(),3000, -12, 30);
		TH1D *l_hist = new TH1D(ssl.str().c_str(),(ssl.str()+";ip3d_llr;entries").c_str(),3000, -12, 30);

		b_ip3d_llr_ptbin.push_back(b_hist);
		l_ip3d_llr_ptbin.push_back(l_hist);


	}


}

void ip3d_llr_ptbin::set_flat_cut(){

	for(int i=0; i<b_ip3d_llr_ptbin.size(); i++ ){

		//std::cout<<"LLR Mean"<<b_ip3d_llr_ptbin[i]->GetMean() <<std::endl;
		for(int j=1; j<=b_ip3d_llr_ptbin[i]->GetNbinsX(); j++){
			double flat_cut = 0.;
			double eff_b = b_ip3d_llr_ptbin[i]->Integral(j,b_ip3d_llr_ptbin[i]->GetNbinsX()) / b_ip3d_llr_ptbin[i]->Integral();
			//std::cout<<"eb: "<<eb<<std::endl
			if(eff_b>0.69 && eff_b<0.71){
				//std::cout<<"eb: "<<eff_b<<std::endl;
				flat_cut = b_ip3d_llr_ptbin[i]->GetBinCenter(j); 
				//std::cout<<"FLAT CUT: "<<flat_cut<<"\tEB"<<eff_b<<std::endl;
				v_flat_cut.push_back(flat_cut);
				break;	
			}

		}
	}

}


TH1D* GetRatio(TGraphAsymmErrors* new_roc, TGraphAsymmErrors* old_roc){
	const int n_bin = 110;
	const double max = 1.1; 
	const int graph_size = new_roc->GetN();

	double* nr_x = new_roc->GetX(); 
	double* nr_y = new_roc->GetY(); 
	double* or_x = old_roc->GetX(); 
	double* or_y = old_roc->GetY(); 

	TH1D* ratio = new TH1D("ratio","ratio; eff_b; rejection_ratio", n_bin, 0, max );

	for(int i = 0; i <= n_bin; i++){
		std::vector<double> ny;
		std::vector<double> oy;
	
		double step_size =(double)(max/(double)n_bin);
		for(int j = 0; j < graph_size; j++){
			if(nr_x[j]>= (i-1)*step_size && nr_x[j]< i*step_size ) ny.push_back(nr_y[j]);
			if(or_x[j]>= (i-1)*step_size && or_x[j]< i*step_size ) oy.push_back(or_y[j]);
		}
		double ny_bin = 0;
		double oy_bin = 0;	
		for(int k  = 0; k<ny.size(); k++) ny_bin += ((double)ny[k]/(double)ny.size());
		for(int k  = 0; k<oy.size(); k++) oy_bin += ((double)oy[k]/(double)oy.size());
		
		//std::cout<<ny_bin<<"\t"<<oy_bin<<std::endl;
		if(oy_bin!=0 && ny_bin!=0) ratio->SetBinContent(i, ny_bin/oy_bin);
	}
	ratio->SetLineStyle(2);
	ratio->SetLineColor(kRed);
	ratio->SetLineWidth(2);
	ratio->Draw();
	double axis_size = 0.06;
	ratio->GetXaxis()->SetTitleSize(axis_size);
	ratio->GetYaxis()->SetTitleSize(axis_size);
	ratio->GetXaxis()->SetLabelSize(axis_size);
	ratio->GetYaxis()->SetLabelSize(axis_size);
	ratio->GetYaxis()->SetTitleOffset(0.7);
	ratio->GetXaxis()->SetRangeUser(0.4,1.0);
	//ratio->GetYaxis()->SetRangeUser(0.,7.);
	return ratio;
}
struct grades_histo{

	std::vector<TH1D*> grades_jpt;
	std::vector<TH1D*> grades_bHLxy;

	std::vector<TH1D*> ratio_grades_jpt;
	std::vector<TH1D*> ratio_grades_bHLxy;

	std::vector<TH1D*> grades_d0sig;
	//std::vector<TString> grade_string = get_track_grades();

	void initHistos(std::string ID);
	void getRatioHistos();
};

void grades_histo::initHistos(std::string ID){
	std::vector<std::string> grade_string = get_track_grades();

	int grade_color[14] = {kBlue+2, kRed, kGreen+2, kAzure-4, kOrange-3, kPink+7, kMagenta, kViolet-6, kSpring, kTeal, kYellow, kGray, kBlack, kCyan };
	for(int i = 0; i<grade_string.size(); i++){
		std::stringstream ss;
		ss.str(std::string());
		ss<<ID<<"_"<<grade_string[i]<<"_jet_pt";
		TH1D* grade_jpt = new TH1D(ss.str().c_str(),(ss.str() + "; jet_pt (MeV); ").c_str(), 300, 0, 3e6);
		TH1D* ratio_grade_jpt = new TH1D(("ratio_"+ss.str()).c_str(),(("ratio_"+ss.str()) + "; jet_pt (MeV); ").c_str(), 300, 0, 3e6);

		grade_jpt->SetFillColor(grade_color[i]);
		ratio_grade_jpt->SetFillColor(grade_color[i]);
		grades_jpt.push_back(grade_jpt);
		ratio_grades_jpt.push_back(ratio_grade_jpt);
		ss.str(std::string());
		ss<<ID<<"_"<<grade_string[i]<<"_bH_Lxy";
		TH1D* grade_bHLxy = new TH1D(ss.str().c_str(),(ss.str() + "; bH_Lxy (mm); ").c_str(), 300, 0, 300);
		TH1D* ratio_grade_bHLxy = new TH1D(("ratio_"+ss.str()).c_str(),(("ratio_"+ss.str()) + "; bH_Lxy (mm); ").c_str(), 300, 0, 300);
		grade_bHLxy->SetFillColor(grade_color[i]);
		ratio_grade_bHLxy->SetFillColor(grade_color[i]);
		grades_bHLxy.push_back(grade_bHLxy);
		ratio_grades_bHLxy.push_back(ratio_grade_bHLxy);

		ss.str(std::string());
		ss<<ID<<"_"<<grade_string[i]<<"_d0sig";
		TH1D* grade_d0sig = new TH1D(ss.str().c_str(),(ss.str()+";d0sig;").c_str(),300,-40,60);
		grades_d0sig.push_back(grade_d0sig);
	}

}

void grades_histo::getRatioHistos(){

	std::vector<std::string> grade_string = get_track_grades();
	for(int j = 1; j <grades_jpt[0]->GetNbinsX(); j++ ){
		double binsum = 0;
		for(int i = 0; i< grade_string.size(); i++) binsum += grades_jpt[i]->GetBinContent(j);
		if(binsum != 0)for(int i = 0; i< grade_string.size(); i++)	ratio_grades_jpt[i]->SetBinContent(j,grades_jpt[i]->GetBinContent(j)/binsum);
	}

	for(int j = 1; j <grades_bHLxy[0]->GetNbinsX(); j++ ){
		double binsum = 0;
		for(int i = 0; i< grade_string.size(); i++) binsum += grades_bHLxy[i]->GetBinContent(j);
		if(binsum!=0)for(int i = 0; i< grade_string.size(); i++)	ratio_grades_bHLxy[i]->SetBinContent(j,grades_bHLxy[i]->GetBinContent(j)/binsum);
	}
}


void light_trackplotter(std::string inputFolder, int jet_flav, double ts, float min_pt = 0., float max_pt = 3e+6){
	TH1D* h_b_jet_pt = new TH1D("b_jet_pt", "b_jet_pt", 1000, 0, 2e6);
	TH1D* h_l_jet_pt = new TH1D("l_jet_pt", "l_jet_pt", 1000, 0, 2e6);

	TH1D* BC_grades_histo = new TH1D("BC_grades","BC_grades;grades;",14,-0.5,13.4);


	TH2D* B_IP3DtrkVsJetPt = new TH2D("B_IP3DtrkVsJetPt","B_IP3DtrkVsJetPt;Jet Pt (MeV); IP3D trk", 100, 0, 3e+6, 31, -0.5, 30.5); B_IP3DtrkVsJetPt->Sumw2();
	TH2D* BC_IP3DtrkVsJetPt = new TH2D("BC_IP3DtrkVsJetPt","BC_IP3DtrkVsJetPt;Jet Pt (MeV); IP3D trk", 100, 0, 3e+6, 31, -0.5, 30.5); BC_IP3DtrkVsJetPt->Sumw2();
	TH2D* BC_trkVsJetPt = new TH2D("BC_trkVsJetPt","BC_trkVsJetPt;Jet Pt (MeV); trk", 100, 0, 3e+6, 31, -0.5, 30.5); BC_trkVsJetPt->Sumw2();
	TH2D* BC_IP3DtrkVsbHLxy = new TH2D("BC_IP3DtrkVsbHLxy","BC_IP3DtrkVsbHLxy; bH_Lxy; IP3D trk", 300, -0.5, 300.5 , 31, -0.5, 30.5); BC_IP3DtrkVsbHLxy->Sumw2();
	TH2D* BC_trkVsbHLxy = new TH2D("BC_trkVsbHLxy","BC_trkVsbHLxy; bH_Lxy; trk", 300, -0.5, 300.5 , 31, -0.5, 30.5); BC_trkVsbHLxy->Sumw2();
	TH1D* B_IP3Dtrk = new TH1D("B_IP3Dtrk","B_IP3Dtrk; IP3D trk", 31, -0.5, 30.5);// B_IP3DtrkVsJetPt->Sumw2();

	grades_histo Bjet_grades, BC_grades, frag_grades, light_grades;
	BC_grades.initHistos("BC_");
	frag_grades.initHistos("frag_");
	light_grades.initHistos("light_");
	Bjet_grades.initHistos("Bjet_");
//	std::vector<TH1D*> BC_d0sig_grades, frag_d0sig_grades, light_d0sig_grades;
//	init_d0sig_grades(BC_d0sig_grades, "BC");
//	init_d0sig_grades(frag_d0sig_grades, "frag");
//	init_d0sig_grades(light_d0sig_grades, "light");

	TH2D* bHpt_vs_jet_pt = new TH2D("bHpt_vs_jet_pt","bHpt_vs_jet_pt;jet_pt [MeV];bH_pt [MeV]", 300, 0, 3e6, 300, 0, 3e6);
	TH2D* bHptratio_vs_jet_pt = new TH2D("bHptratio_vs_jet_pt","bHptratio_vs_jet_pt;jet_pt [MeV];bH_pt [MeV]", 300, 0, 3e6, 300, 0, 10);

	TH1D* BC_d0sig = new TH1D("BC_d0sig","BC_d0sig;d0sig",300, -150,150);
	TH2D* BC_sumtrkpt_jetpt = new TH2D("BC_sumtrkpt_jetpt","BC_sumtrkpt_jetpt;jet_pt (MeV); sumtrk_pt (MeV)",300, 0, 3e6, 300, 0 ,3e6); BC_sumtrkpt_jetpt->Sumw2();
	TH2D* BC_sumtrkpt_bHpt = new TH2D("BC_sumtrkpt_bHpt","BC_sumtrkpt_bHpt;jet_pt (MeV); sumtrk_pt (MeV)",300, 0, 3e6, 300, 0 ,3e6); BC_sumtrkpt_bHpt->Sumw2();
	TH1D* BC_z0sig = new TH1D("BC_z0sig","BC_z0sig;z0sig",300, -150,150);
	TH1D* BC_dz0sig = new TH1D("BC_dz0sig","BC_dz0sig;dz0sig",300, 0,300);
	TH1D* BC_trkpt = new TH1D("BC_trkpt","BC_trkpt;trkpt",300, 0,3e6);
	TH2D* BC_d0sig_vs_z0sig = new TH2D("BC_d0sig_vs_z0sig","BC_d0sig_vs_z0sig;d0sig;z0sig",300,-150,150,300,-150,150);
	TH1D* frag_d0sig = new TH1D("frag_d0sig","frag_d0sig;d0sig",300, -150,150);
	TH2D* frag_sumtrkpt_jetpt = new TH2D("frag_sumtrkpt_jetpt","frag_sumtrkpt_jetpt;jet_pt (MeV); sumtrk_pt (MeV)",300, 0, 3e6, 300, 0 ,3e6); frag_sumtrkpt_jetpt->Sumw2();
	TH2D* frag_sumtrkpt_bHpt = new TH2D("frag_sumtrkpt_bHpt","frag_sumtrkpt_bHpt;jet_pt (MeV); sumtrk_pt (MeV)",300, 0, 3e6, 300, 0 ,3e6); frag_sumtrkpt_bHpt->Sumw2();
	TH1D* frag_z0sig = new TH1D("frag_z0sig","frag_z0sig;z0sig",300, -150,150);
	TH1D* frag_dz0sig = new TH1D("frag_dz0sig","frag_dz0sig;dz0sig",300, 0,300);
	TH1D* frag_trkpt = new TH1D("frag_trkpt","frag_trkpt;trkpt",300, 0,3e6);
	TH2D* frag_d0sig_vs_z0sig = new TH2D("frag_d0sig_vs_z0sig","frag_d0sig_vs_z0sig;d0sig;z0sig",300,-150,150,300,-150,150);

	TH1D* BC_d0sig_custom = new TH1D("BC_d0sig_custom","BC_d0sig_custom;d0sig_custom",300, -150,150);
	TH1D* BC_z0sig_custom = new TH1D("BC_z0sig_custom","BC_z0sig_custom;z0sig_custom",300, -150,150);
	TH2D* BC_d0sig_vs_z0sig_custom = new TH2D("BC_d0sig_vs_z0sig_custom","BC_d0sig_vs_z0sig_custom;d0sig_custom;z0sig_custom",300,-150,150,300,-150,150);
	TH1D* frag_d0sig_custom = new TH1D("frag_d0sig_custom","frag_d0sig_custom;d0sig_custom",300, -150,150);
	TH1D* frag_z0sig_custom = new TH1D("frag_z0sig_custom","frag_z0sig_custom;z0sig_custom",300, -150,150);
	TH2D* frag_d0sig_vs_z0sig_custom = new TH2D("frag_d0sig_vs_z0sig_custom","frag_d0sig_vs_z0sig_custom;d0sig_custom;z0sig_custom",300,-150,150,300,-150,150);


	TH2D* BC_leadOrder_vs_jetPt = new TH2D("BC_leadOrder_vs_jetPt","BC_leadOrder_vs_jetPt;Jet Pt (MeV); Lead Order", 300, 0, 3e+6, 151, -.5,150.5);
	TH2D* BC_leadOrder_vs_bHLxy = new TH2D("BC_leadOrder_vs_bHLxy","BC_leadOrder_vs_bHLxy;bH_Lxy (mm); Lead Order", 100, 0, 300, 151, -.5,150.5);
	TH2D* BC_leadOrder_vs_tvSize = new TH2D("BC_leadOrder_vs_tvSize","BC_leadOrder_vs_tvSize;tvSize; Lead Order", 151, -0.5, 151, 151, -.5,150.5);
	TH2D* BC_tvSize_vs_jetPt = new TH2D("BC_tvSize_vs_jetPt","BC_tvSize_vs_jetPt;Jet Pt (MeV); tvSize", 300, 0, 3e+6, 151, -.5,150.5);
	TH2D* BC_tvSize_vs_bHLxy = new TH2D("BC_tvSize_vs_bHLxy","BC_tvSize_vs_bHLxy;bH_Lxy (mm); tvSize", 100, 0, 300, 151, -.5,150.5);
	TH2D* B_bHLxyVsJetPt = new TH2D("B_bHLxyVsJetPt","B_bHLxyVsJetPt;Jet Pt (MeV); bH_Lxy", 100, 0, 3e6, 301, -0.5,300.5);
	TH2D* B_bHLxyVsJetPt_01 = new TH2D("B_bHLxyVsJetPt_01","B_bHLxyVsJetPt_01;Jet Pt (MeV); bH_Lxy", 100, 0, 3e6, 301, -0.5,300.5);
	TH2D* B_bHLxyVsJetPt_gt1 = new TH2D("B_bHLxyVsJetPt_gt1","B_bHLxyVsJetPt_gt1;Jet Pt (MeV); bH_Lxy", 100, 0, 3e6, 301, -0.5,300.5);

	TH2D* B_nPixHitsvsnSiHits = new TH2D("B_nPixHitsvsnSiHits","B_nPixHitsvsnSiHits; nPixHits; nSiHits", 11, -0.5, 10.5, 21, -0.5,20.2);
	TH2D* B_IP3D_nPixHitsvsnSiHits = new TH2D("B_IP3D_nPixHitsvsnSiHits","B_IP3D_nPixHitsvsnSiHits; nPixHits; nSiHits", 11, -0.5, 10.5, 21, -0.5,20.2);


	TH2D* frag_IP3DtrkVsJetPt = new TH2D("frag_IP3DtrkVsJetPt","frag_IP3DtrkVsJetPt;Jet Pt (MeV); IP3D trk", 100, 0, 3e+6, 31, -0.5, 30.5); frag_IP3DtrkVsJetPt->Sumw2();

	TH2D* Btruth_tracks_vs_jetpt = new TH2D("Btruth_tracks_vs_jetpt","Btruth_tracks_vs_jetpt;Jet Pt (MeV); IP3D trk", 100, 0, 3e+6, 31, -0.5, 30.5); Btruth_tracks_vs_jetpt->Sumw2();
	TH2D* Btruth_tracks_vs_bHLxy = new TH2D("Btruth_tracks_vs_bHLxy","Btruth_tracks_vs_bHLxy;bH_Lxy; IP3D trk", 300, -0.5, 300.5, 31, -0.5, 30.5); Btruth_tracks_vs_bHLxy->Sumw2();

	TH2D* BC_trkpt_vs_jetpt = new TH2D("BC_trkpt_vs_jetpt","BC_trkpt_vs_jetpt;Jet Pt (MeV); trk_pt/jet_pt  ", 100, 0, 3e+6, 300, 0, 1); BC_trkpt_vs_jetpt->Sumw2();
	TH2D* BC_trkpt_vs_jetpt_log = new TH2D("BC_trkpt_vs_jetpt_log","BC_trkpt_vs_jetpt_log;Jet Pt (MeV); log10(trk_pt/jet_pt)  ", 100, 0, 3e+6, 100, -3, 1); BC_trkpt_vs_jetpt_log->Sumw2();
	TH2F* B_npixVsbH_Lxy = new TH2F("B_npixVsbH_Lxy","B_npixVsbH_Lxy; bH_Lxy; pixel hit", 201, -0.5, 200.5,11, -0.5, 10.5);
	TH2F* B_npixVsjetpt = new TH2F("B_npixVsjetpt","B_npixVsjetpt; jet pt; pixel hit", 100, 0, 3e6, 11, -0.5, 10.5);
	TH2D* frag_trkpt_vs_jetpt = new TH2D("frag_trkpt_vs_jetpt","frag_trkpt_vs_jetpt;Jet Pt (MeV); trk_pt/jet_pt ", 100, 0, 3e+6, 300, 0, 1); frag_trkpt_vs_jetpt->Sumw2();
	TH2D* frag_trkpt_vs_jetpt_log = new TH2D("frag_trkpt_vs_jetpt_log","frag_trkpt_vs_jetpt_log;Jet Pt (MeV); log10(trk_pt/jet_pt) ", 100, 0, 3e+6, 100, -3, 1);frag_trkpt_vs_jetpt_log->Sumw2();

	TH2F* B_IP3D_npixVsbH_Lxy = new TH2F("B_IP3D_npixVsbH_Lxy","B_IP3D_npixVsbH_Lxy; bH_Lxy; pixel hit", 201, -0.5, 200.5,11, -0.5, 10.5);
	TH2F* B_IP3D_npixVsjetpt = new TH2F("B_IP3D_npixVsjetpt","B_IP3D_npixVsjetpt; jet pt; pixel hit", 100, 0, 3e6, 11, -0.5, 10.5);

	TH2D* B_dRjetVsjetpt = new TH2D("B_dRjetVsjetpt","B_dRjetVsjetpt;jet pt (MeV); dR", 100, 0, 3e6, 100, 0,0.6);
	TH2D* B_dRjetVsjetpt_sel = new TH2D("B_dRjetVsjetpt_sel","B_dRjetVsjetpt_sel;jet pt (MeV); dR", 100, 0, 3e6, 100, 0,0.6);
	TH2D* B_dRjetVsIP3Dtrk = new TH2D("B_dRjetVsIP3Dtrk","B_dRjetVsIP3Dtrk;IP3Dtrk; dR", 31, -0.5, 31.5, 100, 0,0.6);

	TH2D* BC_dRtrkjetVsjetpt = new TH2D("BC_dRtrkjetVsjetpt","BC_dRtrkjetVsjetpt; jet pt [MeV]; dR (trk -jet)", 100, 0, 3e6, 100, 0, 0.6);
	TH2D* frag_dRtrkjetVsjetpt = new TH2D("frag_dRtrkjetVsjetpt","frag_dRtrkjetVsjetpt; jet pt [MeV]; dR (trk -jet)", 100, 0, 3e6, 100, 0, 0.6);
	TH2D* light_dRtrkjetVsjetpt = new TH2D("light_dRtrkjetVsjetpt","light_dRtrkjetVsjetpt; jet pt [MeV]; dR (trk -jet)", 100, 0, 3e6, 100, 0, 0.6);

	TH2D* BC_d0vsjetpt = new TH2D("BC_d0vsjetpt","BC_d0vsjetpt; jetpt; d0", 100, 0, 3e6, 100, -3, 3);
	TH2D* BC_d0sigvsjetpt = new TH2D("BC_d0sigvsjetpt","BC_d0sigvsjetpt; jetpt; d0", 100, 0, 3e6, 100, -30, 30);
	TH2D* B_z0sinthetavsjetpt = new TH2D("B_z0sinthetavsjetpt","B_z0sinthetavsjetpt; jetpt; z0sintheta", 100, 0, 3e6, 100, -3, 3);

	TH2D* BC_d0vsjetpt_ip3d_sel = new TH2D("BC_d0vsjetpt_ip3d_sel","BC_d0vsjetpt_ip3d_sel; jetpt; d0", 100, 0, 3e6, 100, -3, 3);
	TH2D* B_z0sinthetavsjetpt_ip3d_sel = new TH2D("B_z0sinthetavsjetpt_ip3d_sel","B_z0sinthetavsjetpt_ip3d_sel; jetpt; z0sintheta", 100, 0, 3e6, 100, -3, 3);

	TH1F* ip3d_llr_distro_b = new TH1F("ip3d_llr_distro_b","ip3d_llr_distro_b;ip3d_llr",3000,-12,30);
	TH1F* ip3d_llr_distro_l = new TH1F("ip3d_llr_distro_l","ip3d_llr_distro_l;ip3d_llr",3000,-12,30);

	TH2F* B_UT_trkptVsjetpt = new TH2F(" B_UT_trkptVsjetpt"," B_UT_trkptVsjetpt; jet pt; trk pt", 100, 0, 3e6, 100, 0, 1e6);
	TH2F* B_UT_d0Vsjetpt = new TH2F(" B_UT_d0Vsjetpt"," B_UT_d0Vsjetpt; jet pt; d0", 100, 0, 3e6, 100,-10, 10);
	TH2F* B_UT_z0sinthetaVsjetpt = new TH2F(" B_UT_z0sinthetaVsjetpt"," B_UT_z0sinthetaVsjetpt; jet pt; z0sintheta", 100, 0, 3e6, 100,-10, 10);
	TH2F* B_UT_trkVtxVsjetpt = new TH2F("B_UT_trkVtxVsjetpt","B_UT_trkVtxVsjetpt; jet pt; trk Lxy", 100, 0, 3e6, 100, 0, 500);
	TH2F* B_trkVtxVsjetpt = new TH2F("B_trkVtxVsjetpt","B_trkVtxVsjetpt; jet pt; trk Lxy", 100, 0, 3e6, 100, 0, 500);
	TH2F* B_UT_ip3dgradeVsjetpt = new TH2F("B_UT_ip3dgradeVsjetpt","B_UT_ip3dgradeVsjetpt; jet pt; ip3d grade", 100, 0, 3e6, 31, -10.5, 20.5);
	TH2F* B_UT_npixVsjetpt = new TH2F("B_UT_npixVsjetpt","B_UT_npixVsjetpt; jet pt; pixel hit", 100, 0, 3e6, 11, -0.5, 10.5);
	TH2F* B_UT_npixVsbH_Lxy = new TH2F("B_UT_npixVsbH_Lxy","B_UT_npixVsbH_Lxy; bH_Lxy; pixel hit", 201, -0.5, 200.5,11, -0.5, 10.5);
	TH2F* B_UTnoC_npixVsjetpt = new TH2F("B_UTnoC_npixVsjetpt","B_UTnoC_npixVsjetpt; jet pt; pixel hit", 100, 0, 3e6, 11, -0.5, 10.5);
	TH2F* B_UTnoC_npixVsbH_Lxy = new TH2F("B_UTnoC_npixVsbH_Lxy","B_UTnoC_npixVsbH_Lxy; bH_Lxy; pixel hit", 201, -0.5, 200.5,11, -0.5, 10.5);
	TH2F* frag_UT_npixVsjetpt = new TH2F("frag_UT_npixVsjetpt","frag_UT_npixVsjetpt; jet pt; pixel hit", 100, 0, 3e6, 11, -0.5, 10.5);
	TH2F* B_1pixHitIP3D_etaphimap = new TH2F("B_1pixHitIP3D_etaphimap","B_1pixHitIP3D_etaphimap",100, -3,3,100, 0,3.2);
	TH2F* B_0pixHit_etaphimap = new TH2F("B_0pixHit_etaphimap","B_0pixHit_etaphimap",100, -3,3,100, 0,3.2);
	TH2F* B_1pixHit_etaphimap = new TH2F("B_1pixHit_etaphimap","B_1pixHit_etaphimap",100, -3,3,100, 0,3.2);
	TH2F* frag_0pixHit_etaphimap = new TH2F("frag_0pixHit_etaphimap","frag_0pixHit_etaphimap",100, -3,3,100, 0,3.2);
	TH1F* custom_ip3d_llr_distro_b = new TH1F("custom_ip3d_llr_distro_b","custom_ip3d_llr_distro_b;ip3d_llr",3000,-12,30);
	TH1F* custom_ip3d_llr_distro_l = new TH1F("custom_ip3d_llr_distro_l","custom_ip3d_llr_distro_l;ip3d_llr",3000,-12,30);


	std::cout<<"Initialiasing ip3d_llr_ptbin parameters ...."<<std::endl;
	ip3d_llr_ptbin cust_ip3d_llr, std_ip3d_llr;
	cust_ip3d_llr.Init(100,3e6,"custom" );
	std_ip3d_llr.Init(100,3e6,"std" );

	cust_ip3d_llr.ip3d_cut = 7.507;
	std_ip3d_llr.ip3d_cut = 2.003;

	std::cout<<"Initialiasing ip3d_llr_ptbin histos ...."<<std::endl;

	cust_ip3d_llr.InitHistos();
	std_ip3d_llr.InitHistos();
	std::cout<<"DONE"<<std::endl;
		



	TH2D* vtx_B = new TH2D("vtx_B","vtx_B",200, -50, 50, 200, -50, 50);
	TH2D* vtx_frag = new TH2D("vtx_frag","vtx_frag",200, -50, 50, 200, -50, 50);

	int light_flav = 0;
	//float jpt_frac = 0.001;


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
	//tool->initTrainingMode(14);
	//tool->initEvaluationModeDB( "/afs/cern.ch/atlas/groups/perf-flavtag/ReferenceHistograms/BTagCalibRUN2-08-05.root","AntiKt4EMTopo", grades);
	//tool->initEvaluationMode("/afs/cern.ch/user/a/amiucci/B_tagging/xAODAthena/btagIBLAnalysis/macros/ip3d_tuning_new.root","AntiKt4EMTopo",grades);
	tool->initEvaluationMode("/afs/cern.ch/user/a/amiucci/B_tagging/xAODAthena/btagIBLAnalysis/macros/BTagCalibRUN2-08-15_cand.root","AntiKt4EMTopo",grades);
	

	std::cout<<"Chain Entries:"<<myChain->GetEntries()<<std::endl;
	initBranches(myChain);

	
	std::vector<float> *bH_x = 0;
	std::vector<float> *bH_y = 0;
	std::vector<float> *bH_z = 0;
	Double_t PVx = 0;
	Double_t PVy = 0;
	Double_t PVz = 0;
	std::vector<float> *bH_phi = 0;
	std::vector<float> *bH_eta = 0;
	std::vector<int> *bH_nBtracks = 0;
	std::vector<int> *bH_nCtracks = 0;
	std::vector<float> *bH_pt = 0;
	std::vector<float> *jet_pt = 0;
	std::vector<float> *jet_dRiso = 0;
	std::vector<float> *jet_E = 0;
	std::vector<float> *jet_eta = 0;
	std::vector<float> *jet_phi = 0;
	std::vector<int> *jet_truthMatch = 0;
	std::vector<int> *jet_truthflav = 0;
	std::vector<int> *jet_LabDr_HadF = 0;
	Int_t njets = 0;
	Int_t eventnb = 0;
	std::vector<float> *jet_ip3d_llr = 0;	
	std::vector< std::vector<float> > *jet_trk_ip3d_llr = 0;	
	std::vector<float> *bH_Lxy = 0;	
	std::vector<float> *jet_jf_ntrkAtVx = 0;
	std::vector<int> *jet_sv1_ntrkv = 0;
	std::vector<int> *jet_ip3d_ntrk = 0;
	std::vector<int> *jet_jf_nvtx = 0;
	std::vector<int> *jet_sv1_Nvtx = 0;
	std::vector<std::vector<float> > *jet_trk_d0 = 0;
	std::vector<std::vector<float> > *jet_trk_d0_truth = 0;
	std::vector<std::vector<float> > *jet_trk_ip3d_d0 = 0;
	std::vector<std::vector<float> > *jet_trk_ip3d_d0sig = 0;
	std::vector<std::vector<int> > *jet_trk_ip3d_grade = 0;
	std::vector<std::vector<int> > *jet_trk_orig = 0;
	std::vector<std::vector<int> > *jet_trk_algo = 0;
	std::vector<std::vector<float> > *jet_trk_z0 = 0;
	std::vector<std::vector<float> > *jet_trk_z0_truth = 0;
	std::vector<std::vector<float> > *jet_trk_ip3d_z0 = 0;
	std::vector<std::vector<float> > *jet_trk_ip3d_z0sig = 0;
	std::vector<std::vector<float> > *jet_trk_vtx_X = 0;
	std::vector<std::vector<float> > *jet_trk_vtx_Y = 0;
	std::vector<std::vector<float> > *jet_trk_eta = 0;
	std::vector<std::vector<float> > *jet_trk_theta = 0;
	std::vector<std::vector<float> > *jet_trk_phi = 0;
	std::vector<std::vector<float> > *jet_trk_pt = 0;
	std::vector<std::vector<int> > *jet_trk_nInnHits = 0;
	std::vector<std::vector<int> > *jet_trk_nPixHits = 0;
	std::vector<std::vector<int> > *jet_trk_nSCTHits = 0;
	std::vector<std::vector<int> > *jet_trk_nNextToInnHits = 0;
	std::vector<int> *jet_btag_ntrk = 0;
	std::vector<int> *jet_aliveAfterOR =0;   
	std::vector<float> *bH_dRjet = 0;

	TBranch *b_PVx, *b_PVy, *b_PVz, *b_bH_dRjet, *b_jet_btag_ntrk, *b_eventnb,*b_jet_aliveAfterOR, *b_njets, *b_jet_pt,*b_jet_E, *b_jet_eta, *b_jet_phi, *b_jet_truthMatch,*b_jet_LabDr_HadF, *b_jet_truthflav, *b_jet_ip3d_llr, *b_jet_trk_ip3d_llr;
	TBranch *b_jet_trk_d0, *b_jet_trk_ip3d_d0, *b_jet_trk_d0_truth,*b_jet_trk_nPixHits,*b_jet_trk_nSCTHits, *b_jet_trk_ip3d_d0sig, *b_jet_trk_vtx_X, *b_jet_trk_vtx_Y, *b_jet_trk_pt;
	TBranch *b_jet_trk_z0, *b_jet_trk_ip3d_z0, *b_jet_trk_z0_truth, *b_jet_trk_ip3d_z0sig, *b_jet_dRiso, *b_jet_trk_phi, *b_jet_trk_eta, *b_jet_trk_theta;
	TBranch *b_bH_Lxy, *b_bH_x, *b_bH_y, *b_bH_z, *b_bH_pt, *b_bH_eta, *b_bH_nBtracks, *b_bH_nCtracks, *b_bH_phi, *b_jet_jf_ntrkAtVx, *b_jet_sv1_Nvtx, *b_jet_jf_nvtx,*b_jet_sv1_ntrkv,*b_jet_ip3d_ntrk;
	TBranch *b_jet_trk_nInnHits, *b_jet_trk_nNextToInnHits, *b_jet_trk_ip3d_grade,*b_jet_trk_algo, *b_jet_trk_orig;
	myChain->SetBranchAddress("njets", &njets, &b_njets);
	myChain->SetBranchAddress("eventnb", &eventnb, &b_eventnb);
  myChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
  myChain->SetBranchAddress("jet_E", &jet_E, &b_jet_E);
  myChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
  myChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
	myChain->SetBranchAddress("jet_truthMatch", &jet_truthMatch, &b_jet_truthMatch);
	myChain->SetBranchAddress("jet_truthflav", &jet_truthflav, &b_jet_truthflav);
	myChain->SetBranchAddress("jet_LabDr_HadF", &jet_LabDr_HadF, &b_jet_LabDr_HadF);
	myChain->SetBranchAddress("jet_ip3d_llr", &jet_ip3d_llr, &b_jet_ip3d_llr);
	myChain->SetBranchAddress("jet_trk_ip3d_llr", &jet_trk_ip3d_llr, &b_jet_trk_ip3d_llr);
	myChain->SetBranchAddress("jet_trk_vtx_X", &jet_trk_vtx_X, &b_jet_trk_vtx_X);
	myChain->SetBranchAddress("jet_trk_vtx_Y", &jet_trk_vtx_Y, &b_jet_trk_vtx_Y);
	myChain->SetBranchAddress("jet_trk_d0", &jet_trk_d0, &b_jet_trk_d0);
	myChain->SetBranchAddress("jet_trk_d0_truth", &jet_trk_d0_truth, &b_jet_trk_d0_truth);
	myChain->SetBranchAddress("jet_trk_ip3d_d0", &jet_trk_ip3d_d0, &b_jet_trk_ip3d_d0);
	myChain->SetBranchAddress("jet_trk_ip3d_d0sig", &jet_trk_ip3d_d0sig, &b_jet_trk_ip3d_d0sig);
	myChain->SetBranchAddress("jet_ip3d_ntrk", &jet_ip3d_ntrk, &b_jet_ip3d_ntrk);
	myChain->SetBranchAddress("jet_trk_theta", &jet_trk_theta, &b_jet_trk_theta);
	myChain->SetBranchAddress("jet_trk_z0", &jet_trk_z0, &b_jet_trk_z0);
	myChain->SetBranchAddress("jet_trk_z0_truth", &jet_trk_z0_truth, &b_jet_trk_z0_truth);
	myChain->SetBranchAddress("jet_trk_ip3d_z0", &jet_trk_ip3d_z0, &b_jet_trk_ip3d_z0);
	myChain->SetBranchAddress("jet_trk_ip3d_z0sig", &jet_trk_ip3d_z0sig, &b_jet_trk_ip3d_z0sig);
	myChain->SetBranchAddress("jet_trk_ip3d_grade", &jet_trk_ip3d_grade, &b_jet_trk_ip3d_grade);
	myChain->SetBranchAddress("jet_trk_orig", &jet_trk_orig, &b_jet_trk_orig);
	myChain->SetBranchAddress("jet_trk_algo", &jet_trk_algo, &b_jet_trk_algo);
	myChain->SetBranchAddress("jet_trk_pt", &jet_trk_pt, &b_jet_trk_pt);
	myChain->SetBranchAddress("jet_trk_eta", &jet_trk_eta, &b_jet_trk_eta);
	myChain->SetBranchAddress("jet_dRiso", &jet_dRiso, &b_jet_dRiso);
	myChain->SetBranchAddress("jet_trk_phi", &jet_trk_phi, &b_jet_trk_phi);
	myChain->SetBranchAddress("PVx", &PVx, &b_PVx);
	myChain->SetBranchAddress("PVy", &PVy, &b_PVy);
	myChain->SetBranchAddress("PVz", &PVz, &b_PVz);
	myChain->SetBranchAddress("bH_Lxy", &bH_Lxy, &b_bH_Lxy);
	myChain->SetBranchAddress("bH_y", &bH_y, &b_bH_y);
	myChain->SetBranchAddress("bH_nBtracks", &bH_nBtracks, &b_bH_nBtracks);
	myChain->SetBranchAddress("bH_nCtracks", &bH_nCtracks, &b_bH_nCtracks);
	myChain->SetBranchAddress("bH_dRjet", &bH_dRjet, &b_bH_dRjet);
	myChain->SetBranchAddress("bH_x", &bH_x, &b_bH_x);
	myChain->SetBranchAddress("bH_z", &bH_z, &b_bH_z);
	myChain->SetBranchAddress("bH_eta", &bH_eta, &b_bH_eta);
	myChain->SetBranchAddress("bH_phi", &bH_phi, &b_bH_phi);
	myChain->SetBranchAddress("bH_pt", &bH_pt, &b_bH_pt);
	myChain->SetBranchAddress("jet_jf_ntrkAtVx", &jet_jf_ntrkAtVx, &b_jet_jf_ntrkAtVx);
	myChain->SetBranchAddress("jet_jf_nvtx", &jet_jf_nvtx, &b_jet_jf_nvtx);
	myChain->SetBranchAddress("jet_sv1_ntrkv", &jet_sv1_ntrkv, &b_jet_sv1_ntrkv);
	myChain->SetBranchAddress("jet_sv1_Nvtx", &jet_sv1_Nvtx, &b_jet_sv1_Nvtx);	
	myChain->SetBranchAddress("jet_trk_nNextToInnHits", &jet_trk_nNextToInnHits, &b_jet_trk_nNextToInnHits);
	myChain->SetBranchAddress("jet_trk_nInnHits", &jet_trk_nInnHits, &b_jet_trk_nInnHits);
	myChain->SetBranchAddress("jet_trk_nPixHits", &jet_trk_nPixHits, &b_jet_trk_nPixHits);
	myChain->SetBranchAddress("jet_trk_nSCTHits", &jet_trk_nSCTHits, &b_jet_trk_nSCTHits);
  myChain->SetBranchAddress("jet_btag_ntrk", &jet_btag_ntrk, &b_jet_btag_ntrk);
	myChain->SetBranchAddress("jet_aliveAfterOR", &jet_aliveAfterOR, &b_jet_aliveAfterOR);

	Long64_t nentries = myChain->GetEntriesFast();

	int n_b_tracks = 0;
	int n_b_sel_tracks =0;

	for(Long64_t jentry =0; jentry<nentries; jentry++){
		Long64_t ientry = myChain->LoadTree(jentry);
		if(ientry<0) break;
		if(ientry%1000 == 0)std::cout<<ientry<<"\t"<<jentry<<std::endl;
		b_njets->GetEntry(ientry);
		b_eventnb->GetEntry(ientry);
		for(int i = 0; i < njets; i++){
			b_jet_pt->GetEntry(ientry);
			b_bH_pt->GetEntry(ientry);
			b_jet_truthMatch->GetEntry(ientry);
			b_jet_aliveAfterOR->GetEntry(ientry);
  		if (jet_truthMatch  ->at(i)!=1) continue;
  		if (jet_aliveAfterOR->at(i)!=1) continue;
  		if (jet_pt->at(i)<25e3)         continue;
			b_jet_truthflav->GetEntry(ientry);
			b_jet_LabDr_HadF->GetEntry(ientry);
			b_jet_ip3d_llr->GetEntry(ientry);
			
			//jpt_frac = ( 0.0553861 - 5.46001e-8* pow(jet_pt->at(i),1) + 2.85982e-14 * pow(jet_pt->at(i),2) + 5.69244e-21 * pow(jet_pt->at(i),3) + 1.75064e-28*pow(jet_pt->at(i),4) );
			//jpt_frac = ( 0.028 - 5.46001e-8* pow(jet_pt->at(i),1) + 2.85982e-14 * pow(jet_pt->at(i),2) + 5.69244e-21 * pow(jet_pt->at(i),3) + 1.75064e-28*pow(jet_pt->at(i),4) );
			//double jpt_frac = jf;
			long double jpt_frac =	0.5*	(0.0980166 +
										-7.4378e-1 * pow(jet_pt->at(i)*1e-6,1) +
										3.3711 * pow(jet_pt->at(i)*1e-6,2) +
										-7.68983 * pow(jet_pt->at(i)*1e-6,3) +
										9.94882 * pow(jet_pt->at(i)*1e-6,4) +
										-7.80337 * pow(jet_pt->at(i)*1e-6,5) +
										3.78802 * pow(jet_pt->at(i)*1e-6,6) +
										-1.11231 * pow(jet_pt->at(i)*1e-6,7) +
										1.8105e-1 * pow(jet_pt->at(i)*1e-6,8) +
										-1.25372e-2 * pow(jet_pt->at(i)*1e-6,9));
			b_jet_trk_pt->GetEntry(ientry);
			b_jet_trk_ip3d_z0sig->GetEntry(ientry);
			b_jet_trk_ip3d_d0sig->GetEntry(ientry);
			b_bH_eta->GetEntry(ientry);
			b_bH_phi->GetEntry(ientry);
			b_jet_trk_eta->GetEntry(ientry);
			b_jet_trk_phi->GetEntry(ientry);	
			b_jet_phi->GetEntry(ientry);
			b_jet_eta->GetEntry(ientry);
			b_jet_dRiso->GetEntry(ientry);			


			Int_t idx_trkpt[jet_trk_pt->at(i).size()];
			float *trk_pt_vec = &(jet_trk_pt->at(i)[0]);
			TMath::Sort(int (jet_trk_pt->at(i).size()), trk_pt_vec, idx_trkpt);

			Int_t idx_trkdz0sig[jet_trk_ip3d_d0sig->at(i).size()];
			float trk_dz0sig_vec[jet_trk_ip3d_d0sig->at(i).size()];
			for(int l=0; l<jet_trk_ip3d_d0sig->at(i).size();l++){
				if(jet_trk_ip3d_d0sig->at(i)[l] == -999)trk_dz0sig_vec[l] = -999;
				else trk_dz0sig_vec[l] =/*(jet_trk_ip3d_z0sig->at(i)[l]/abs(jet_trk_ip3d_z0sig->at(i)[l]))*/ jet_trk_pt->at(i)[l] * sqrt((pow(jet_trk_ip3d_d0sig->at(i)[l],2) + pow(jet_trk_ip3d_z0sig->at(i)[l],2)));
			}
			TMath::Sort(int (jet_trk_ip3d_d0sig->at(i).size()), trk_dz0sig_vec, idx_trkdz0sig);
			//std::cout<<(long double)jpt_frac<<"\t"<< jet_pt->at(i)<<std::endl;
			//if(jet_truthflav->at(i) == jet_flav){
			if(jet_LabDr_HadF->at(i) == jet_flav){
				b_jet_eta->GetEntry(ientry);
				h_b_jet_pt->Fill(jet_pt->at(i));
				if(jet_truthMatch->at(i) ==1 && jet_aliveAfterOR->at(i) ==1 && jet_pt->at(i) > 25e3 &&  abs(jet_eta->at(i)) < eta_cut && jet_pt->at(i)>min_pt && jet_pt->at(i)<max_pt){

					bHpt_vs_jet_pt->Fill(jet_pt->at(i),bH_pt->at(i));					
					bHptratio_vs_jet_pt->Fill(jet_pt->at(i),bH_pt->at(i)/jet_pt->at(i));					

					b_jet_trk_d0->GetEntry(ientry);
					b_jet_trk_z0->GetEntry(ientry);
					b_jet_trk_ip3d_d0->GetEntry(ientry);
					b_jet_trk_orig->GetEntry(ientry);
					b_jet_trk_algo->GetEntry(ientry);
					b_jet_trk_vtx_X->GetEntry(ientry);	
					b_jet_trk_vtx_Y->GetEntry(ientry);	
					b_bH_nBtracks->GetEntry(ientry);
					b_bH_nCtracks->GetEntry(ientry);
					b_bH_dRjet->GetEntry(ientry);
					b_bH_Lxy->GetEntry(ientry);
					b_jet_trk_ip3d_grade->GetEntry(ientry);
					b_jet_trk_ip3d_z0->GetEntry(ientry);
					b_jet_trk_theta->GetEntry(ientry);
					b_jet_ip3d_llr->GetEntry(ientry);
					b_jet_btag_ntrk->GetEntry(ientry);
					b_PVx->GetEntry(ientry);
					b_PVy->GetEntry(ientry);
					b_PVz->GetEntry(ientry);
					b_jet_trk_nPixHits->GetEntry(ientry);
					b_jet_trk_nSCTHits->GetEntry(ientry);

	
					std::vector<int> ip3d_trk_grade_sel;
					std::vector<float> ip3d_d0sig_sel;	
					std::vector<float> ip3d_z0sig_sel;	
					std::vector<float> d0 = jet_trk_d0->at(i);
					std::vector<float> z0 = jet_trk_z0->at(i);
					std::vector<float> ip3d_d0 = jet_trk_ip3d_d0->at(i);
					std::vector<float> ip3d_z0 = jet_trk_ip3d_z0->at(i);
					std::vector<float> theta = jet_trk_theta->at(i);
					std::vector<float> trk_pt = jet_trk_pt->at(i);
					std::vector<float> vtx_X = jet_trk_vtx_X->at(i);
					std::vector<float> vtx_Y = jet_trk_vtx_Y->at(i);
					std::vector<int> trk_orig = jet_trk_orig->at(i);
					std::vector<int> trk_algo = jet_trk_algo->at(i);
					
					ip3d_llr_distro_b->Fill(jet_ip3d_llr->at(i));
					for(int m = 0; m<std_ip3d_llr.nbin; m++){
						if( jet_pt->at(i) >= std_ip3d_llr.pt_bins[m] && jet_pt->at(i) < std_ip3d_llr.pt_bins[m+1])std_ip3d_llr.b_ip3d_llr_ptbin[m]->Fill(jet_ip3d_llr->at(i));
					}
					
					int n_b_ip3d_trk = 0;
					int n_bc_ip3d_trk = 0;
					int n_bc_trk = 0;
					int n_frag_ip3d_trk = 0;
				
	
					/*std::vector<track_data> tVector;
					for(int k = 0; k<trk_pt.size();k++){
						track_data tD;
						tD.trk_pt = trk_pt[k];
						tD.trk_orig = trk_orig[k];
						tD.trk_grade = jet_trk_ip3d_grade->at(i)[k]; 
						tVector.push_back(tD);
					}
					std::sort(tVector.begin(),tVector.end(),by_pt());
					for(unsigned int k=0; k<tVector.size(); k++){
						//std::cout<<tVector[k].trk_pt<<std::endl;
					}*/


	
					//for(unsigned int k =0; k<trk_pt.size(); k++){
					int trk_size = trk_pt.size();
					//std::cout << std::endl;
					double BC_sumtrk_pt = 0.;
					double frag_sumtrk_pt = 0.;
					for(unsigned int ok =0; ok<trk_size; ok++){
						//int k = idx_trkpt[ok];
						int k = idx_trkdz0sig[ok];
						double d0_sin = TMath::Sin( (bH_phi->at(i) - jet_trk_phi->at(i)[k]) * (TMath::Pi()/180.));
						double z0_sin = TMath::Sin( (bH_eta->at(i) - jet_trk_eta->at(i)[k]) * (TMath::Pi()/180.));
						//std::cout<<trk_dz0sig_vec[k]<<"\t"<<k<<"\t"<<ok<<std::endl;
						//int k = ok;
						if(jet_trk_ip3d_grade->at(i)[k] >-1 &&( jet_trk_orig->at(i)[k] == 0 || jet_trk_orig->at(i)[k] == 1)) {
							BC_leadOrder_vs_jetPt->Fill(jet_pt->at(i),ok);
							BC_leadOrder_vs_bHLxy->Fill(bH_Lxy->at(i),ok);
							BC_tvSize_vs_jetPt->Fill(jet_pt->at(i), trk_size);
							BC_tvSize_vs_bHLxy->Fill(bH_Lxy->at(i), trk_size);
							BC_leadOrder_vs_tvSize->Fill(trk_size, ok);
							//if(jet_trk_orig->at(i)[k] == 0 ){
							BC_d0sig->Fill(jet_trk_ip3d_d0sig->at(i)[k]);
							BC_z0sig->Fill(jet_trk_ip3d_z0sig->at(i)[k]);
							BC_dz0sig->Fill(trk_dz0sig_vec[k]);
							BC_trkpt->Fill(jet_trk_pt->at(i)[k]);
							BC_d0sig_vs_z0sig->Fill(jet_trk_ip3d_d0sig->at(i)[k],jet_trk_ip3d_z0sig->at(i)[k]);
							double d0_sig_custom = (d0_sin/TMath::Abs(d0_sin))*abs(jet_trk_ip3d_d0sig->at(i)[k])*(jet_trk_d0->at(i)[k]/abs(jet_trk_d0->at(i)[k]));
							double z0_sig_custom = (z0_sin/TMath::Abs(z0_sin))*abs(jet_trk_ip3d_z0sig->at(i)[k])*((jet_trk_z0->at(i)[k] -PVz) /abs(jet_trk_z0->at(i)[k] -PVz));
							BC_d0sig_custom->Fill(d0_sig_custom);
							BC_z0sig_custom->Fill(z0_sig_custom);
							BC_d0sig_vs_z0sig_custom->Fill((double)d0_sig_custom, (double)z0_sig_custom);
							BC_grades_histo->Fill(jet_trk_ip3d_grade->at(i)[k]);
						}
						int ta = trk_algo[k];
						bool unsel_trk = false;					
						B_nPixHitsvsnSiHits->Fill(jet_trk_nPixHits->at(i)[k], jet_trk_nPixHits->at(i)[k]+jet_trk_nSCTHits->at(i)[k]);
						if(jet_trk_ip3d_grade->at(i)[k] >-1 && ( trk_orig[k] == 0 || trk_orig[k] == 1)) B_IP3D_nPixHitsvsnSiHits->Fill(jet_trk_nPixHits->at(i)[k],jet_trk_nPixHits->at(i)[k]+jet_trk_nSCTHits->at(i)[k]);
						if(trk_orig[k] == 0 || trk_orig[k] == 1){
								unsel_trk = true;
								BC_d0vsjetpt->Fill(jet_pt->at(i), d0[k] );
								BC_d0sigvsjetpt->Fill(jet_pt->at(i), jet_trk_ip3d_d0sig->at(i)[k] );
								B_z0sinthetavsjetpt->Fill(jet_pt->at(i), (z0[k]-(PVz))*sin(theta[k]));

						}	

						if(trk_orig[k] == 0 || trk_orig[k] == 1){
							B_npixVsjetpt->Fill(jet_pt->at(i), (jet_trk_nPixHits->at(i))[k]);
							B_npixVsbH_Lxy->Fill(bH_Lxy->at(i), (jet_trk_nPixHits->at(i))[k]);
							n_bc_trk++;
							BC_sumtrk_pt+=trk_pt[k];
						}

						bool lead_cut = (ok < ts);
						//if( jet_trk_ip3d_grade->at(i)[k] >-1 && abs(ip3d_d0[k]) < 1  && trk_pt[k] >= jpt_frac*jet_pt->at(i)){
						bool nPixCut = true;// ((jet_trk_nPixHits->at(i))[k]>=2 && jet_pt->at(i) < 4e5) || ((jet_trk_nPixHits->at(i))[k]>=0 && jet_pt->at(i) >= 4e5);

						if(nPixCut && lead_cut &&jet_trk_ip3d_grade->at(i)[k] >-1){
                 Bjet_grades.grades_jpt[jet_trk_ip3d_grade->at(i)[k]]->Fill(jet_pt->at(i));
                 Bjet_grades.grades_bHLxy[jet_trk_ip3d_grade->at(i)[k]]->Fill(bH_Lxy->at(i));
                 Bjet_grades.grades_d0sig[jet_trk_ip3d_grade->at(i)[k]]->Fill(jet_trk_ip3d_d0sig->at(i)[k]);
						}
	
						if(jet_trk_ip3d_grade->at(i)[k] >-1 &&( trk_orig[k] == 0 || trk_orig[k] == 1)) {
							if(trk_orig[k]==0) n_b_ip3d_trk++;
							n_bc_ip3d_trk++;
							/*if(jet_dRiso->at(i)>0.8)*/BC_dRtrkjetVsjetpt->Fill(jet_pt->at(i), sqrt( pow(jet_trk_phi->at(i)[k]-jet_phi->at(i),2) + pow(jet_trk_eta->at(i)[k]-jet_eta->at(i),2) ));
							B_IP3D_npixVsjetpt->Fill(jet_pt->at(i), (jet_trk_nPixHits->at(i))[k]);
							B_IP3D_npixVsbH_Lxy->Fill(bH_Lxy->at(i), (jet_trk_nPixHits->at(i))[k]);
							BC_trkpt_vs_jetpt->Fill(jet_pt->at(i),(trk_pt[k]/jet_pt->at(i)));
							BC_trkpt_vs_jetpt_log->Fill(jet_pt->at(i),log10(trk_pt[k]/jet_pt->at(i)));
							vtx_B->Fill(vtx_X[k],vtx_Y[k]);
							BC_d0vsjetpt_ip3d_sel->Fill(jet_pt->at(i), d0[k] );
							B_z0sinthetavsjetpt_ip3d_sel->Fill(jet_pt->at(i), (z0[k]-(PVz))*sin(theta[k]));
							B_trkVtxVsjetpt->Fill(jet_pt->at(i), sqrt(pow(vtx_X[k],2)+pow(vtx_Y[k],2)));
							if(nPixCut && lead_cut){	
								BC_grades.grades_jpt[jet_trk_ip3d_grade->at(i)[k]]->Fill(jet_pt->at(i));	
								BC_grades.grades_bHLxy[jet_trk_ip3d_grade->at(i)[k]]->Fill(bH_Lxy->at(i));	
								BC_grades.grades_d0sig[jet_trk_ip3d_grade->at(i)[k]]->Fill(jet_trk_ip3d_d0sig->at(i)[k]);
							}
							if((jet_trk_nPixHits->at(i))[k] == 1) B_1pixHitIP3D_etaphimap->Fill((jet_trk_eta->at(i))[k],(jet_trk_phi->at(i))[k]);
							unsel_trk = false;
						}else if(( trk_orig[k] == 0 || trk_orig[k] == 1)){
							B_UT_trkptVsjetpt->Fill(jet_pt->at(i),trk_pt[k]);
							B_UT_d0Vsjetpt->Fill(jet_pt->at(i), d0[k]);	
							B_UT_z0sinthetaVsjetpt->Fill(jet_pt->at(i),( z0[k]-PVz) * sin(theta[k]));
							B_UT_trkVtxVsjetpt->Fill(jet_pt->at(i), sqrt(pow(vtx_X[k],2)+pow(vtx_Y[k],2)));	
							B_UT_ip3dgradeVsjetpt->Fill(jet_pt->at(i), (jet_trk_ip3d_grade->at(i))[k]);
							B_UT_npixVsjetpt->Fill(jet_pt->at(i), (jet_trk_nPixHits->at(i))[k]);
							B_UT_npixVsbH_Lxy->Fill(bH_Lxy->at(i), (jet_trk_nPixHits->at(i))[k]);
							if(trk_orig[k] == 0) B_UTnoC_npixVsjetpt->Fill(jet_pt->at(i), (jet_trk_nPixHits->at(i))[k]);
							if(trk_orig[k] == 0) B_UTnoC_npixVsbH_Lxy->Fill(bH_Lxy->at(i), (jet_trk_nPixHits->at(i))[k]);
							if((jet_trk_nPixHits->at(i))[k]) B_0pixHit_etaphimap->Fill((jet_trk_eta->at(i))[k],(jet_trk_phi->at(i))[k]);
							if((jet_trk_nPixHits->at(i))[k]) B_1pixHit_etaphimap->Fill((jet_trk_eta->at(i))[k],(jet_trk_phi->at(i))[k]);
						}
					
/*						if(unsel_trk && trk_pt[k] > 1e3 && abs( (z0[k]-(PVz))*sin(theta[k]) ) < 1.5 && abs(d0[k]) < 1){
							B_UT_trkptVsjetpt->Fill(jet_pt->at(i),trk_pt[k]);
							B_UT_d0Vsjetpt->Fill(jet_pt->at(i), d0[k]);	
							B_UT_z0sinthetaVsjetpt->Fill(jet_pt->at(i),( z0[k]-PVz) * sin(theta[k]));
							B_UT_trkVtxVsjetpt->Fill(jet_pt->at(i), sqrt(pow(vtx_X[k],2)+pow(vtx_Y[k],2)));	
							B_UT_ip3dgradeVsjetpt->Fill(jet_pt->at(i), (jet_trk_ip3d_grade->at(i))[k]);
							B_UT_npixVsjetpt->Fill(jet_pt->at(i), (jet_trk_nPixHits->at(i))[k]);
						}
	*/
						if(jet_trk_ip3d_grade->at(i)[k] >-1 && trk_orig[k] == 2 ){
							n_frag_ip3d_trk++;
							/*if(jet_dRiso->at(i)>0.8)*/frag_dRtrkjetVsjetpt->Fill(jet_pt->at(i), sqrt( pow(jet_trk_phi->at(i)[k]-jet_phi->at(i),2) + pow(jet_trk_eta->at(i)[k]-jet_eta->at(i),2) ));
							frag_trkpt_vs_jetpt->Fill(jet_pt->at(i), (trk_pt[k]/jet_pt->at(i)));
							frag_trkpt_vs_jetpt_log->Fill(jet_pt->at(i), log10(trk_pt[k]/jet_pt->at(i)));
							if(nPixCut && lead_cut){
								frag_grades.grades_jpt[jet_trk_ip3d_grade->at(i)[k]]->Fill(jet_pt->at(i));	
								frag_grades.grades_bHLxy[jet_trk_ip3d_grade->at(i)[k]]->Fill(bH_Lxy->at(i));	
								frag_grades.grades_d0sig[jet_trk_ip3d_grade->at(i)[k]]->Fill(jet_trk_ip3d_d0sig->at(i)[k]);
							}
							vtx_frag->Fill(vtx_X[k],vtx_Y[k]);
							frag_d0sig->Fill(jet_trk_ip3d_d0sig->at(i)[k]);
							frag_z0sig->Fill(jet_trk_ip3d_z0sig->at(i)[k]);
							frag_dz0sig->Fill(trk_dz0sig_vec[k]);
							frag_trkpt->Fill(jet_trk_pt->at(i)[k]);
							frag_d0sig_vs_z0sig->Fill(jet_trk_ip3d_d0sig->at(i)[k],jet_trk_ip3d_z0sig->at(i)[k]);
							//frag_d0sig_custom->Fill((d0_sin/TMath::Abs(d0_sin))*abs(jet_trk_ip3d_d0sig->at(i)[k]));
							//frag_z0sig_custom->Fill((z0_sin/TMath::Abs(z0_sin))*abs(jet_trk_ip3d_z0sig->at(i)[k]));
							//frag_d0sig_vs_z0sig_custom->Fill(abs(jet_trk_ip3d_d0sig->at(i)[k]),abs(jet_trk_ip3d_z0sig->at(i)[k]));
							double d0_sig_custom = (d0_sin/TMath::Abs(d0_sin))*abs(jet_trk_ip3d_d0sig->at(i)[k])*(jet_trk_d0->at(i)[k]/abs(jet_trk_d0->at(i)[k]));
							double z0_sig_custom = (z0_sin/TMath::Abs(z0_sin))*abs(jet_trk_ip3d_z0sig->at(i)[k])*((jet_trk_z0->at(i)[k] -PVz) /abs(jet_trk_z0->at(i)[k] -PVz));
							frag_d0sig_custom->Fill(d0_sig_custom);
							frag_z0sig_custom->Fill(z0_sig_custom);
							frag_d0sig_vs_z0sig_custom->Fill((double)d0_sig_custom, (double)z0_sig_custom);
							frag_sumtrk_pt+=trk_pt[k];
						}else if( trk_orig[k] == 2 ){
							frag_UT_npixVsjetpt->Fill(jet_pt->at(i), (jet_trk_nPixHits->at(i))[k]);
							frag_sumtrk_pt+=trk_pt[k];
							if((jet_trk_nPixHits->at(i))[k] == 0) frag_0pixHit_etaphimap->Fill((jet_trk_eta->at(i))[k],(jet_trk_phi->at(i))[k]  );
						}
						//bool lead_cut = (trk_pt.size() > 10 && ok < ts * trk_pt.size());
						if(  nPixCut && jet_trk_ip3d_grade->at(i)[k] >-1  /* &&  trk_pt[k] >= jpt_frac*jet_pt->at(i) */ ){
							n_b_tracks++;
							if(lead_cut){
								n_b_sel_tracks++;
								ip3d_d0sig_sel.push_back( (jet_trk_ip3d_d0sig->at(i))[k] );
								ip3d_z0sig_sel.push_back( (jet_trk_ip3d_z0sig->at(i))[k] );
								ip3d_trk_grade_sel.push_back( (jet_trk_ip3d_grade->at(i))[k] );
								//if(trk_orig[k] == 2)frag_d0sig_grades[jet_trk_ip3d_grade->at(i)[k]]->Fill(jet_trk_ip3d_d0sig->at(i)[k]);
							}
						}
					}
					float custom_ip3d_llr = tool->getJetLLR(3,ip3d_trk_grade_sel ,ip3d_d0sig_sel, ip3d_z0sig_sel, jet_flav, light_flav);

					if(ip3d_d0sig_sel.size()  != 0){
						custom_ip3d_llr_distro_b->Fill(custom_ip3d_llr);
						for(int m = 0; m<cust_ip3d_llr.nbin; m++){
							if( jet_pt->at(i) >= cust_ip3d_llr.pt_bins[m] && jet_pt->at(i) < cust_ip3d_llr.pt_bins[m+1])cust_ip3d_llr.b_ip3d_llr_ptbin[m]->Fill(custom_ip3d_llr);
						}
					}
					//else std::cout<<"Numb Trk in Jet higher than cut: "<< ip3d_d0sig_sel.size()  <<std::endl;
					if(BC_sumtrk_pt){
						BC_sumtrkpt_jetpt->Fill(jet_pt->at(i), BC_sumtrk_pt);
						BC_sumtrkpt_bHpt->Fill(bH_pt->at(i), BC_sumtrk_pt);

					}
					if(frag_sumtrk_pt){
						frag_sumtrkpt_jetpt->Fill(jet_pt->at(i), frag_sumtrk_pt);
						frag_sumtrkpt_bHpt->Fill(bH_pt->at(i), frag_sumtrk_pt);

					}
					Btruth_tracks_vs_jetpt->Fill(jet_pt->at(i), bH_nBtracks->at(i) + bH_nCtracks->at(i));
					Btruth_tracks_vs_bHLxy->Fill( bH_Lxy->at(i) , bH_nBtracks->at(i) + bH_nCtracks->at(i));

					B_dRjetVsjetpt->Fill(jet_pt->at(i), bH_dRjet->at(i));
					if(n_b_ip3d_trk == 0)B_dRjetVsjetpt_sel->Fill(jet_pt->at(i), bH_dRjet->at(i));
					B_dRjetVsIP3Dtrk->Fill(n_b_ip3d_trk, bH_dRjet->at(i));
					B_IP3Dtrk->Fill(n_b_ip3d_trk);
					B_bHLxyVsJetPt->Fill(jet_pt->at(i),bH_Lxy->at(i));
					if(n_bc_ip3d_trk <2) B_bHLxyVsJetPt_01->Fill(jet_pt->at(i),bH_Lxy->at(i));
					if(n_bc_ip3d_trk >=2) B_bHLxyVsJetPt_gt1->Fill(jet_pt->at(i),bH_Lxy->at(i));
					B_IP3DtrkVsJetPt->Fill(jet_pt->at(i),n_b_ip3d_trk);	
					BC_IP3DtrkVsJetPt->Fill(jet_pt->at(i),n_bc_ip3d_trk);	
					BC_trkVsJetPt->Fill(jet_pt->at(i),n_bc_trk);	
					BC_IP3DtrkVsbHLxy->Fill(bH_Lxy->at(i),n_bc_ip3d_trk);	
					BC_trkVsbHLxy->Fill(bH_Lxy->at(i),n_bc_trk);	
					frag_IP3DtrkVsJetPt->Fill(jet_pt->at(i),n_frag_ip3d_trk);	
				}
			//}else if(jet_truthflav->at(i) == light_flav){
			}else if(jet_LabDr_HadF->at(i) == light_flav){
				h_l_jet_pt->Fill(jet_pt->at(i));
				b_jet_truthMatch->GetEntry(ientry);
				//b_jet_trk_phi->GetEntry(ientry);
				//b_jet_trk_eta->GetEntry(ientry);
				//b_jet_phi->GetEntry(ientry);
				//b_jet_eta->GetEntry(ientry);
				b_jet_pt->GetEntry(ientry);
				b_bH_Lxy->GetEntry(ientry);
				b_jet_aliveAfterOR->GetEntry(ientry);
				if(jet_truthMatch->at(i) ==1 && jet_aliveAfterOR->at(i) ==1 && abs(jet_eta->at(i)) < eta_cut && jet_pt->at(i)>min_pt && jet_pt->at(i)<max_pt){
				
					b_jet_trk_d0->GetEntry(ientry);
					b_jet_trk_ip3d_d0->GetEntry(ientry);
					//b_jet_trk_pt->GetEntry(ientry);
					b_jet_trk_ip3d_d0sig->GetEntry(ientry);
					b_jet_trk_ip3d_z0sig->GetEntry(ientry);
					b_jet_trk_ip3d_d0->GetEntry(ientry);
					b_jet_trk_ip3d_z0->GetEntry(ientry);
					b_jet_trk_algo->GetEntry(ientry);
					b_jet_trk_ip3d_grade->GetEntry(ientry);
					b_jet_trk_theta->GetEntry(ientry);
					b_jet_btag_ntrk->GetEntry(ientry);
					b_jet_trk_nPixHits->GetEntry(ientry);
					std::vector<int> ip3d_trk_grade_sel;
					std::vector<float> ip3d_d0sig_sel;	
					std::vector<float> ip3d_z0sig_sel;	
					std::vector<float> d0 = jet_trk_d0->at(i);
					std::vector<float> ip3d_z0 = jet_trk_ip3d_z0->at(i);
					std::vector<float> ip3d_d0 = jet_trk_ip3d_d0->at(i);
					std::vector<float> trk_pt = jet_trk_pt->at(i);
					std::vector<float> theta = jet_trk_theta->at(i);
					std::vector<int> trk_algo = jet_trk_algo->at(i);

					int n_b_ip3d_trk = 0;
					int n_frag_ip3d_trk = 0;
					
					ip3d_llr_distro_l->Fill(jet_ip3d_llr->at(i));
					for(int m = 0; m<std_ip3d_llr.nbin; m++){
						if( jet_pt->at(i) >= std_ip3d_llr.pt_bins[m] && jet_pt->at(i) < std_ip3d_llr.pt_bins[m+1])std_ip3d_llr.l_ip3d_llr_ptbin[m]->Fill(jet_ip3d_llr->at(i));
					}
					//for(unsigned int k =0; k<trk_pt.size(); k++){
					for(unsigned int ok =0; ok<trk_pt.size(); ok++){
						//int k = idx_trkpt[ok];
						int k = idx_trkdz0sig[ok];
						int ta = trk_algo[k];
						bool nPixCut = true;// ((jet_trk_nPixHits->at(i))[k]>=2 && jet_pt->at(i) < 4e5) || ((jet_trk_nPixHits->at(i))[k]>=0 && jet_pt->at(i) >= 4e5);
						//bool lead_cut = (trk_pt.size() > 10 &&ok < ts * trk_pt.size());
						bool lead_cut = (ok < ts);
						/*if(jet_dRiso->at(i)>0.8)*/light_dRtrkjetVsjetpt->Fill(jet_pt->at(i), sqrt( pow(jet_trk_phi->at(i)[k]-jet_phi->at(i),2) + pow(jet_trk_eta->at(i)[k]-jet_eta->at(i),2) ));

						if(nPixCut && jet_trk_ip3d_grade->at(i)[k] >-1  && lead_cut/* && trk_pt[k] >= jpt_frac*jet_pt->at(i)*/  ){
							
							light_grades.grades_jpt[jet_trk_ip3d_grade->at(i)[k]]->Fill(jet_pt->at(i));	
							light_grades.grades_bHLxy[jet_trk_ip3d_grade->at(i)[k]]->Fill(bH_Lxy->at(i));	
							light_grades.grades_d0sig[jet_trk_ip3d_grade->at(i)[k]]->Fill(jet_trk_ip3d_d0sig->at(i)[k]);
	//						light_d0sig_grades[jet_trk_ip3d_grade->at(i)[k]]->Fill(jet_trk_ip3d_d0sig->at(i)[k]);
							ip3d_d0sig_sel.push_back( (jet_trk_ip3d_d0sig->at(i))[k] );
							ip3d_z0sig_sel.push_back( (jet_trk_ip3d_z0sig->at(i))[k] );
							ip3d_trk_grade_sel.push_back( (jet_trk_ip3d_grade->at(i))[k] );
						}//else std::cout<<"CAZZO"<<std::endl;
					}	
					float custom_ip3d_llr = tool->getJetLLR(3,ip3d_trk_grade_sel ,ip3d_d0sig_sel, ip3d_z0sig_sel, jet_flav, light_flav);

					if(ip3d_d0sig_sel.size() != 0){
						custom_ip3d_llr_distro_l->Fill(custom_ip3d_llr);
						for(int m = 0; m<cust_ip3d_llr.nbin; m++){
							if( jet_pt->at(i) >= cust_ip3d_llr.pt_bins[m] && jet_pt->at(i) < cust_ip3d_llr.pt_bins[m+1])cust_ip3d_llr.l_ip3d_llr_ptbin[m]->Fill(custom_ip3d_llr);
						}
					}
				}
			}
		}
	}
/*	TCanvas *c1 = new TCanvas("c1","c1",1200,800);
	FicoPlot();
	c1->Divide(2,2);
	c1->cd(1);
	B_IP3DtrkVsJetPt->Draw("colz");
	//BC_trkpt_vs_jetpt->Draw("colz");
	c1->cd(2);
	//TProfile* prof = B_IP3DtrkVsJetPt->ProfileX("_pfx",1,-1,"s");
	//prof->SetLineColor(kGreen+3);
	//Btruth_tracks_vs_jetpt->ProfileX("_pfx",1,-1,"s")->Draw();
	//Btruth_tracks_vs_jetpt->Draw("colz");
	//B_IP3Dtrk->Draw();
	B_dRjetVsjetpt->Draw("colz");
	//frag_trkpt_vs_jetpt->Draw("colz");
	//prof->Draw("same");
	c1->cd(3);
	//frag_IP3DtrkVsJetPt->Draw("colz");
//	TProfile* prof = Btruth_tracks_vs_jetpt->ProfileX("_pfx",1,-1,"s");
//	prof->SetLineColor(kGreen+3);
//	prof->Draw();
//	B_IP3DtrkVsJetPt->ProfileX("_pfx",1,-1,"s")->Draw("same");
	B_dRjetVsIP3Dtrk->Draw("colz");
	c1->cd(4);
	//TProfile* prof_frag = frag_IP3DtrkVsJetPt->ProfileX("_pfx",1,-1,"s");
	//prof_frag->SetLineColor(kRed);
	//prof_frag->Draw();
	B_dRjetVsjetpt_sel->Draw("colz");
*/
	TCanvas *c2 = new TCanvas("c2","c2",800,600);
	c2->Divide(3,1);
	c2->cd(1);

	gPad->SetPad(0.00,0.35,0.4,1.0);
	gPad->SetBottomMargin(0);
	std::vector<double> eff_b;
	std::vector<double> rej_l;

	double web[2] ={0};
	double wrl[2] ={0};

	for(int i=1; i<=ip3d_llr_distro_b->GetNbinsX(); i++){
		double eb = ip3d_llr_distro_b->Integral(i,ip3d_llr_distro_b->GetNbinsX()) / ip3d_llr_distro_b->Integral();
		double el = ip3d_llr_distro_l->Integral(i,ip3d_llr_distro_l->GetNbinsX()) / ip3d_llr_distro_l->Integral();

		if(eb>0.69 && eb <0.71) std_ip3d_llr.seventy_ip3d_cut = ip3d_llr_distro_b->GetBinCenter(i);
	}
	for(int i=1; i<=ip3d_llr_distro_b->GetNbinsX(); i++){
		double eb = ip3d_llr_distro_b->Integral(i,ip3d_llr_distro_b->GetNbinsX()) / ip3d_llr_distro_b->Integral();
		double el = ip3d_llr_distro_l->Integral(i,ip3d_llr_distro_l->GetNbinsX()) / ip3d_llr_distro_l->Integral();

		if(abs(ip3d_llr_distro_b->GetBinCenter(i)-std_ip3d_llr.seventy_ip3d_cut) < 0.01){
			web[1]=eb;
			wrl[1]=1/el;
			std::cout<<"IP3D: "<<ip3d_llr_distro_b->GetBinCenter(i)<<"\tstd_ip3d_beff: "<<web[1]<<"\tstd_ip3d_lrej: "<<wrl[1]<<std::endl;
		}
		eff_b.push_back(eb);
		if(el>0)rej_l.push_back(1/el);
		else rej_l.push_back(-1);

	}

	std::vector<double> c_eff_b;
	std::vector<double> c_rej_l;
	for(int i=1; i<=custom_ip3d_llr_distro_b->GetNbinsX(); i++){
		double eb = custom_ip3d_llr_distro_b->Integral(i,custom_ip3d_llr_distro_b->GetNbinsX()) / custom_ip3d_llr_distro_b->Integral();
		double el = custom_ip3d_llr_distro_l->Integral(i,custom_ip3d_llr_distro_l->GetNbinsX()) / custom_ip3d_llr_distro_l->Integral();
		if(eb>0.69 && eb <0.71) cust_ip3d_llr.seventy_ip3d_cut = ip3d_llr_distro_b->GetBinCenter(i);
	}

	for(int i=1; i<=custom_ip3d_llr_distro_b->GetNbinsX(); i++){
		double eb = custom_ip3d_llr_distro_b->Integral(i,custom_ip3d_llr_distro_b->GetNbinsX()) / custom_ip3d_llr_distro_b->Integral();
		double el = custom_ip3d_llr_distro_l->Integral(i,custom_ip3d_llr_distro_l->GetNbinsX()) / custom_ip3d_llr_distro_l->Integral();
		if(abs(custom_ip3d_llr_distro_b->GetBinCenter(i)-cust_ip3d_llr.seventy_ip3d_cut) < 0.01){
			web[0]=eb;
			wrl[0]=1/el;
			std::cout<<"IP3D: "<<custom_ip3d_llr_distro_b->GetBinCenter(i)<<"\tcustom_ip3d_beff: "<<web[0]<<"\tcustom_ip3d_lrej: "<<wrl[0]<<std::endl;
		}

		c_eff_b.push_back(eb);
		if(el>0)c_rej_l.push_back(1/el);
		else c_rej_l.push_back(-1);
	}
	gPad->SetLogy();
	TGraphAsymmErrors* ROC_graph = new TGraphAsymmErrors(eff_b.size(),&(eff_b[0]),&(rej_l[0]));
	TGraphAsymmErrors* c_ROC_graph = new TGraphAsymmErrors(c_eff_b.size(),&(c_eff_b[0]),&(c_rej_l[0]));
	TGraph* WP = new TGraph(2,web,wrl);
	std::stringstream ss;
	ss.str(std::string());
	ss<<"custom_ip3d_llr_customnPix_trksize_"<<ts<<"_minpt"<<min_pt<<"_maxpt"<<max_pt<<";eff_b;rej_l";	
	std::stringstream s;
	s.str(std::string());
	s<<"custom_ip3d_llr_customnPix_trksize_"<<ts<<"_minpt"<<min_pt<<"_maxpt"<<max_pt;	
	std::stringstream ls;
	ls<<"ip3d_llr_minpt"<<min_pt<<"_maxpt"<<max_pt;	
	ROC_graph->SetName(ls.str().c_str());
	ROC_graph->SetTitle("ip3d_llr; b-tagging efficiency;Light-jet rejection");
	c_ROC_graph->SetTitle(ss.str().c_str());
	c_ROC_graph->SetName(s.str().c_str());
	c_ROC_graph->SetMarkerColor(kRed);
	c_ROC_graph->SetMarkerStyle(21);
	c_ROC_graph->SetMarkerSize(0.2);
	ROC_graph->SetMarkerColor(kBlack);
	ROC_graph->SetMarkerStyle(21);
	ROC_graph->SetMarkerSize(0.2);
	ROC_graph->Draw("AP");
	ROC_graph->GetYaxis()->SetRangeUser(1,1e+6);	
	ROC_graph->GetXaxis()->SetRangeUser(0.4,1.0);	
	ROC_graph->Draw("AP");
	c_ROC_graph->Draw("Psame");
	
	WP->SetMarkerStyle(30);
	WP->SetMarkerColor(kOrange);
	WP->SetMarkerSize(2);
	WP->Draw("Psame");
	gPad->BuildLegend();
	c2->cd(2);
	gPad->SetPad(0.00,0.00,0.4,0.35);
	gPad->SetTopMargin(0);
	gPad->SetBottomMargin(0.28);
	gPad->SetGridx();
	gPad->SetGridy();
	gStyle->SetOptStat(0);
	TH1D* ratio = GetRatio(c_ROC_graph, ROC_graph);
	c2->cd(3);
	gPad->SetPad(0.4,0.0,1.0,1.0);
	ip3d_llr_distro_b->SetLineColor(kRed);
	custom_ip3d_llr_distro_b->SetLineColor(kRed+2);
	custom_ip3d_llr_distro_b->SetFillColor(kRed+2);
	custom_ip3d_llr_distro_b->SetFillStyle(3004);

	ip3d_llr_distro_l->SetLineColor(kBlue);
	custom_ip3d_llr_distro_l->SetLineColor(kBlue+2);
	custom_ip3d_llr_distro_l->SetFillColor(kBlue+2);
	custom_ip3d_llr_distro_l->SetFillStyle(3003);
//	ip3d_llr_distro_l->Draw("same");
//
	THStack *llr_stack = new THStack();


	llr_stack->Add(ip3d_llr_distro_l);
	llr_stack->Add(ip3d_llr_distro_b);
	llr_stack->Add(custom_ip3d_llr_distro_b);
	llr_stack->Add(custom_ip3d_llr_distro_l);
	llr_stack->Draw("nostack");
	//ip3d_llr_distro_l->Draw();
	//ip3d_llr_distro_b->Draw("same");
	gPad->SetLogy();
	//custom_ip3d_llr_distro_b->Draw("same");
	//custom_ip3d_llr_distro_l->Draw("same");

	TLegend* leg = new TLegend(0.4,0.7,0.9,0.9);
	leg->AddEntry(custom_ip3d_llr_distro_b, "Custom Ip3d B-jets");
	leg->AddEntry(custom_ip3d_llr_distro_l, "Custom Ip3d L-jets");
	leg->AddEntry(ip3d_llr_distro_b, "Ip3d B-jets");
	leg->AddEntry(ip3d_llr_distro_l, "Ip3d L-jets");
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->Draw();


	TCanvas *c3 = new TCanvas("eff_plot","eff_plot", 1200, 400); 
	c3->Divide(3,2);
	c3->cd(1);
	gPad->SetPad(0.0,0.00,0.33,1.0);

	cust_ip3d_llr.set_flat_cut();
	std_ip3d_llr.set_flat_cut();

	for(int m=0; m<cust_ip3d_llr.nbin; m++){
		double num = cust_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral(cust_ip3d_llr.b_ip3d_llr_ptbin[m]->FindBin(cust_ip3d_llr.seventy_ip3d_cut),cust_ip3d_llr.b_ip3d_llr_ptbin[m]->GetNbinsX());
		double den = cust_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral();
		double eb = num/den;

//		double eb = cust_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral(cust_ip3d_llr.b_ip3d_llr_ptbin[m]->FindBin(cust_ip3d_llr.seventy_ip3d_cut),cust_ip3d_llr.b_ip3d_llr_ptbin[m]->GetNbinsX()) /cust_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral();
		double err_eb = sqrt( pow(sqrt(num)/den,2) + pow( (num/pow(den,2)) *sqrt(den) ,2) );//cust_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral(cust_ip3d_llr.b_ip3d_llr_ptbin[m]->FindBin(cust_ip3d_llr.seventy_ip3d_cut),cust_ip3d_llr.b_ip3d_llr_ptbin[m]->GetNbinsX()) /cust_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral();
		if(eb>-0.1){
			cust_ip3d_llr.b_eff->SetBinContent(m+1,eb);	
			cust_ip3d_llr.b_eff->SetBinError(m+1,err_eb);
		}
		//std::cout<<"eb: "<<eb<<std::endl;
	}
	for(int m=0; m<std_ip3d_llr.nbin; m++){
		double num = std_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral(std_ip3d_llr.b_ip3d_llr_ptbin[m]->FindBin(std_ip3d_llr.seventy_ip3d_cut),std_ip3d_llr.b_ip3d_llr_ptbin[m]->GetNbinsX());
		double den = std_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral();
		double eb = num/den;//std_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral(std_ip3d_llr.b_ip3d_llr_ptbin[m]->FindBin(std_ip3d_llr.seventy_ip3d_cut),std_ip3d_llr.b_ip3d_llr_ptbin[m]->GetNbinsX()) /std_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral();
		double err_eb = sqrt( pow(sqrt(num)/den,2) + pow( (num/pow(den,2)) *sqrt(den) ,2) );//cust_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral(cust_ip3d_llr.b_ip3d_llr_ptbin[m]->FindBin(cust_ip3d_llr.seventy_ip3d_cut),cust_ip3d_llr.b_ip3d_llr_ptbin[m]->GetNbinsX()) /cust_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral();
		if(eb>-0.1){
			std_ip3d_llr.b_eff->SetBinContent(m+1,eb);	
			std_ip3d_llr.b_eff->SetBinError(m+1,err_eb);
		}
		//std::cout<<"eb: "<<eb<<std::endl;
	}

	for(int m=0; m<cust_ip3d_llr.nbin; m++){
		double num = cust_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral(cust_ip3d_llr.b_ip3d_llr_ptbin[m]->FindBin(cust_ip3d_llr.v_flat_cut[m]),cust_ip3d_llr.b_ip3d_llr_ptbin[m]->GetNbinsX());
		double den = cust_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral();

		double eb = num/den;//cust_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral(cust_ip3d_llr.b_ip3d_llr_ptbin[m]->FindBin(cust_ip3d_llr.v_flat_cut[m]),cust_ip3d_llr.b_ip3d_llr_ptbin[m]->GetNbinsX()) /cust_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral();
		double err_eb = sqrt( pow(sqrt(num)/den,2) + pow( (num/pow(den,2)) *sqrt(den) ,2) );//cust_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral(cust_ip3d_llr.b_ip3d_llr_ptbin[m]->FindBin(cust_ip3d_llr.seventy_ip3d_cut),cust_ip3d_llr.b_ip3d_llr_ptbin[m]->GetNbinsX()) /cust_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral();
		if(eb>-0.1){
			cust_ip3d_llr.b_eff_flat->SetBinContent(m+1,eb);	
			cust_ip3d_llr.b_eff_flat->SetBinError(m+1,err_eb);
		}
		//std::cout<<"eb: "<<eb<<std::endl;
	}
	for(int m=0; m<std_ip3d_llr.nbin; m++){
		
		double num = std_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral(std_ip3d_llr.b_ip3d_llr_ptbin[m]->FindBin(std_ip3d_llr.v_flat_cut[m]),std_ip3d_llr.b_ip3d_llr_ptbin[m]->GetNbinsX());
		double den = std_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral();


		double eb = num/den;//std_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral(std_ip3d_llr.b_ip3d_llr_ptbin[m]->FindBin(std_ip3d_llr.v_flat_cut[m]),std_ip3d_llr.b_ip3d_llr_ptbin[m]->GetNbinsX()) /std_ip3d_llr.b_ip3d_llr_ptbin[m]->Integral();
		double err_eb = sqrt( pow(sqrt(num)/den,2) + pow( (num/pow(den,2)) *sqrt(den) ,2) );
		if(eb>-0.1){
			std_ip3d_llr.b_eff_flat->SetBinContent(m+1,eb);	
			std_ip3d_llr.b_eff_flat->SetBinError(m+1,err_eb);
		}
		//std::cout<<"eb: "<<eb<<std::endl;
	}

	cust_ip3d_llr.b_eff->Draw();
	cust_ip3d_llr.b_eff->GetYaxis()->SetRangeUser(0,1);
	cust_ip3d_llr.b_eff->SetLineColor(kRed);
	cust_ip3d_llr.b_eff->SetMarkerColor(kRed);
	cust_ip3d_llr.b_eff->SetMarkerStyle(22);
	cust_ip3d_llr.b_eff->SetMarkerSize(0.5);
	cust_ip3d_llr.b_eff->Rebin(4);
	cust_ip3d_llr.b_eff->Scale(0.25);
	cust_ip3d_llr.b_eff->GetXaxis()->SetRangeUser(0,2e6);//Scale(0.25);
	cust_ip3d_llr.b_eff->Draw("PLE");	
	std_ip3d_llr.b_eff->SetLineColor(kBlack);
	//std_ip3d_llr.b_eff->SetLineStyle(2);
	std_ip3d_llr.b_eff->SetMarkerColor(kBlack);
	std_ip3d_llr.b_eff->SetMarkerStyle(21);
	std_ip3d_llr.b_eff->SetMarkerSize(0.5);
	std_ip3d_llr.b_eff->Rebin(4);
	std_ip3d_llr.b_eff->Scale(0.25);
	std_ip3d_llr.b_eff->Draw("PLEsame");
	cust_ip3d_llr.b_eff_flat->SetLineColor(kRed);
	cust_ip3d_llr.b_eff_flat->SetMarkerColor(kRed);
	cust_ip3d_llr.b_eff_flat->SetMarkerStyle(26);
	cust_ip3d_llr.b_eff_flat->SetLineStyle(2);
	cust_ip3d_llr.b_eff_flat->SetMarkerSize(0.5);
	cust_ip3d_llr.b_eff_flat->Rebin(4);
	cust_ip3d_llr.b_eff_flat->Scale(0.25);
	cust_ip3d_llr.b_eff_flat->Draw("PLEsame");	
	std_ip3d_llr.b_eff_flat->SetLineColor(kBlack);
	std_ip3d_llr.b_eff_flat->SetLineStyle(2);
	std_ip3d_llr.b_eff_flat->SetMarkerColor(kBlack);
	std_ip3d_llr.b_eff_flat->SetMarkerStyle(25);
	std_ip3d_llr.b_eff_flat->SetMarkerSize(0.5);
	std_ip3d_llr.b_eff_flat->Rebin(4);
	std_ip3d_llr.b_eff_flat->Scale(0.25);
	std_ip3d_llr.b_eff_flat->Draw("PLEsame");

	c3->cd(2);
	gPad->SetPad(0,0,0,0);

	c3->cd(3);
	gPad->SetPad(0.33,1.00,0.66,0.33);
	gPad->SetLogy();
	gPad->SetBottomMargin(0);
	gPad->SetLeftMargin(0.2);



	for(int m=0;m<cust_ip3d_llr.nbin; m++){
		double den = cust_ip3d_llr.l_ip3d_llr_ptbin[m]->Integral(cust_ip3d_llr.l_ip3d_llr_ptbin[m]->FindBin(cust_ip3d_llr.seventy_ip3d_cut),cust_ip3d_llr.l_ip3d_llr_ptbin[m]->GetNbinsX());
		double num = cust_ip3d_llr.l_ip3d_llr_ptbin[m]->Integral();
		double lr = num/den;
		double err_lr = sqrt( pow(sqrt(num)/den,2) + pow( (num/pow(den,2)) *sqrt(den) ,2) );
		double el = den/num;//cust_ip3d_llr.l_ip3d_llr_ptbin[m]->Integral(cust_ip3d_llr.l_ip3d_llr_ptbin[m]->FindBin(cust_ip3d_llr.seventy_ip3d_cut),cust_ip3d_llr.l_ip3d_llr_ptbin[m]->GetNbinsX()) /cust_ip3d_llr.l_ip3d_llr_ptbin[m]->Integral();
		if(el > 0){
			cust_ip3d_llr.l_rej->SetBinContent(m+1,lr);	
			cust_ip3d_llr.l_rej->SetBinError(m+1,err_lr);
		}
	}
	for(int m=0; m<std_ip3d_llr.nbin; m++){
		double den = std_ip3d_llr.l_ip3d_llr_ptbin[m]->Integral(std_ip3d_llr.l_ip3d_llr_ptbin[m]->FindBin(std_ip3d_llr.seventy_ip3d_cut),std_ip3d_llr.l_ip3d_llr_ptbin[m]->GetNbinsX());
		double num = std_ip3d_llr.l_ip3d_llr_ptbin[m]->Integral();
		double lr = num/den;
		double err_lr = sqrt( pow(sqrt(num)/den,2) + pow( (num/pow(den,2)) *sqrt(den) ,2) );
		double el = den/num;
		//double el = std_ip3d_llr.l_ip3d_llr_ptbin[m]->Integral(std_ip3d_llr.l_ip3d_llr_ptbin[m]->FindBin(std_ip3d_llr.seventy_ip3d_cut),std_ip3d_llr.l_ip3d_llr_ptbin[m]->GetNbinsX()) /std_ip3d_llr.l_ip3d_llr_ptbin[m]->Integral();
		if(el > 0){
			std_ip3d_llr.l_rej->SetBinContent(m+1,lr);	
			std_ip3d_llr.l_rej->SetBinError(m+1,err_lr);
		}
	}
	

	for(int m=0;m<cust_ip3d_llr.nbin; m++){
		double den = cust_ip3d_llr.l_ip3d_llr_ptbin[m]->Integral(cust_ip3d_llr.l_ip3d_llr_ptbin[m]->FindBin(cust_ip3d_llr.v_flat_cut[m]),cust_ip3d_llr.l_ip3d_llr_ptbin[m]->GetNbinsX());
		double num = cust_ip3d_llr.l_ip3d_llr_ptbin[m]->Integral();
		double lr = num/den;
		double err_lr = sqrt( pow(sqrt(num)/den,2) + pow( (num/pow(den,2)) *sqrt(den) ,2) );
		double el = den/num;//cust_ip3d_llr.l_ip3d_llr_ptbin[m]->Integral(cust_ip3d_llr.l_ip3d_llr_ptbin[m]->FindBin(cust_ip3d_llr.v_flat_cut[m]),cust_ip3d_llr.l_ip3d_llr_ptbin[m]->GetNbinsX()) /cust_ip3d_llr.l_ip3d_llr_ptbin[m]->Integral();
		//double el = cust_ip3d_llr.l_ip3d_llr_ptbin[m]->Integral(cust_ip3d_llr.l_ip3d_llr_ptbin[m]->FindBin(cust_ip3d_llr.v_flat_cut[m]),cust_ip3d_llr.l_ip3d_llr_ptbin[m]->GetNbinsX()) /cust_ip3d_llr.l_ip3d_llr_ptbin[m]->Integral();
		if(el > 0){
			cust_ip3d_llr.l_rej_flat->SetBinContent(m+1,lr);	
			cust_ip3d_llr.l_rej_flat->SetBinError(m+1,err_lr);	
		}
	}
	for(int m=0; m<std_ip3d_llr.nbin; m++){
		double den = std_ip3d_llr.l_ip3d_llr_ptbin[m]->Integral(std_ip3d_llr.l_ip3d_llr_ptbin[m]->FindBin(std_ip3d_llr.v_flat_cut[m]),std_ip3d_llr.l_ip3d_llr_ptbin[m]->GetNbinsX());
		double num = std_ip3d_llr.l_ip3d_llr_ptbin[m]->Integral();
		double lr = num/den;
		double err_lr = sqrt( pow(sqrt(num)/den,2) + pow( (num/pow(den,2)) *sqrt(den) ,2) );
		double el = den/num;//std_ip3d_llr.l_ip3d_llr_ptbin[m]->Integral(std_ip3d_llr.l_ip3d_llr_ptbin[m]->FindBin(std_ip3d_llr.v_flat_cut[m]),std_ip3d_llr.l_ip3d_llr_ptbin[m]->GetNbinsX()) /std_ip3d_llr.l_ip3d_llr_ptbin[m]->Integral();
		//double el = std_ip3d_llr.l_ip3d_llr_ptbin[m]->Integral(std_ip3d_llr.l_ip3d_llr_ptbin[m]->FindBin(std_ip3d_llr.v_flat_cut[m]),std_ip3d_llr.l_ip3d_llr_ptbin[m]->GetNbinsX()) /std_ip3d_llr.l_ip3d_llr_ptbin[m]->Integral();
		if(el > 0){
			std_ip3d_llr.l_rej_flat->SetBinContent(m+1,lr);	
			std_ip3d_llr.l_rej_flat->SetBinError(m+1,err_lr);	

		}
	}


	cust_ip3d_llr.l_rej->Draw();
	cust_ip3d_llr.l_rej->GetYaxis()->SetRangeUser(1,1e4);
	cust_ip3d_llr.l_rej->SetLineColor(kRed);
	cust_ip3d_llr.l_rej->SetMarkerColor(kRed);
	cust_ip3d_llr.l_rej->SetMarkerStyle(22);
	cust_ip3d_llr.l_rej->SetMarkerSize(0.5);
	cust_ip3d_llr.l_rej->Rebin(4);
	cust_ip3d_llr.l_rej->Scale(0.25);
	cust_ip3d_llr.l_rej->GetYaxis()->SetTitleSize(0.05);
	cust_ip3d_llr.l_rej->GetYaxis()->SetLabelSize(0.05);
	cust_ip3d_llr.l_rej->GetXaxis()->SetLabelSize(0.0);
	cust_ip3d_llr.l_rej->GetXaxis()->SetRangeUser(0,2e6);//GetYaxis()->SetLabelSize(0.05);
	cust_ip3d_llr.l_rej->GetYaxis()->SetTitle("Light-jet rejection (fix cut)");
	cust_ip3d_llr.l_rej->SetTitle("");
	cust_ip3d_llr.l_rej->Draw("PLE");	
	
	std_ip3d_llr.l_rej->SetLineColor(kBlack);
	//std_ip3d_llr.l_rej->SetLineStyle(2);
	std_ip3d_llr.l_rej->SetMarkerColor(kBlack);
	std_ip3d_llr.l_rej->SetMarkerStyle(21);
	std_ip3d_llr.l_rej->SetMarkerSize(0.5);
	std_ip3d_llr.l_rej->Rebin(4);
	std_ip3d_llr.l_rej->Scale(0.25);
	std_ip3d_llr.l_rej->Draw("PLEsame");

	TLegend* leg_fix = new TLegend(0.6,0.9,0.9,0.7);
	leg_fix->SetBorderSize(0);
	leg_fix->SetFillStyle(0);
	leg_fix->AddEntry(cust_ip3d_llr.l_rej,"new tuning","PL");
	leg_fix->AddEntry(std_ip3d_llr.l_rej,"standard tuning","PL");
	leg_fix->Draw();

	c3->cd(4);
	gPad->SetTopMargin(0);
	gPad->SetBottomMargin(0.2);
	gPad->SetLeftMargin(0.2);
	gPad->SetPad(0.33,0.33,0.66,0.0);
	TH1D* l_rej_ratio_fix = new TH1D(*cust_ip3d_llr.l_rej);
	l_rej_ratio_fix->Divide(std_ip3d_llr.l_rej);
	l_rej_ratio_fix->SetLineColor(kRed);
	l_rej_ratio_fix->SetMarkerColor(kRed);
	l_rej_ratio_fix->SetMarkerStyle(22);
	l_rej_ratio_fix->SetMarkerSize(0.5);
	l_rej_ratio_fix->GetXaxis()->SetTitleSize(0.1);
	l_rej_ratio_fix->GetXaxis()->SetLabelSize(0.1);
	l_rej_ratio_fix->GetXaxis()->SetRangeUser(0,2e6);//GetXaxis()->SetLabelSize(0.1);
	l_rej_ratio_fix->GetYaxis()->SetTitleSize(0.1);
	l_rej_ratio_fix->GetYaxis()->SetTitleOffset(0.5);
	l_rej_ratio_fix->GetYaxis()->SetLabelSize(0.1);
	l_rej_ratio_fix->GetYaxis()->SetTitle("ratio");
	l_rej_ratio_fix->SetTitle("");
	l_rej_ratio_fix->Draw("PLE");
	gPad->SetGridy();

	c3->cd(5);
	gPad->SetPad(0.66,1.00,1.0,0.33);
	gPad->SetLogy();
	gPad->SetBottomMargin(0);
	gPad->SetLeftMargin(0.2);
	//cust_ip3d_llr.l_rej_flat->SetLineStyle(2);
	cust_ip3d_llr.l_rej_flat->SetLineColor(kRed);
	cust_ip3d_llr.l_rej_flat->SetMarkerColor(kRed);
	cust_ip3d_llr.l_rej_flat->SetMarkerStyle(22);
	cust_ip3d_llr.l_rej_flat->SetMarkerSize(0.5);
	cust_ip3d_llr.l_rej_flat->Rebin(4);
	cust_ip3d_llr.l_rej_flat->Scale(0.25);
	cust_ip3d_llr.l_rej_flat->GetYaxis()->SetTitleSize(0.05);
	cust_ip3d_llr.l_rej_flat->GetYaxis()->SetLabelSize(0.05);
	cust_ip3d_llr.l_rej_flat->GetXaxis()->SetLabelSize(0.0);
	cust_ip3d_llr.l_rej_flat->GetXaxis()->SetRangeUser(0,2e6);//->SetLabelSize(0.05);
	cust_ip3d_llr.l_rej_flat->GetYaxis()->SetTitle("Light-jet rejection (flat cut)");
	cust_ip3d_llr.l_rej_flat->SetTitle("");
	cust_ip3d_llr.l_rej_flat->Draw("PLE");	
	std_ip3d_llr.l_rej_flat->SetLineColor(kBlack);
	//std_ip3d_llr.l_rej_flat->SetLineStyle(2);
	std_ip3d_llr.l_rej_flat->SetMarkerColor(kBlack);
	std_ip3d_llr.l_rej_flat->SetMarkerStyle(21);
	std_ip3d_llr.l_rej_flat->SetMarkerSize(0.5);
	std_ip3d_llr.l_rej_flat->Rebin(4);
	std_ip3d_llr.l_rej_flat->Scale(0.25);
	std_ip3d_llr.l_rej_flat->Draw("PLEsame");
	TLegend* leg_flat = new TLegend(0.6,0.9,0.9,0.7);
	leg_flat->SetBorderSize(0);
	leg_flat->SetFillStyle(0);
	leg_flat->AddEntry(cust_ip3d_llr.l_rej_flat,"new tuning","PL");
	leg_flat->AddEntry(std_ip3d_llr.l_rej_flat,"standard tuning","PL");
	leg_flat->Draw();
	c3->cd(6);
	gPad->SetTopMargin(0);
	gPad->SetBottomMargin(0.2);
	gPad->SetLeftMargin(0.2);
	gPad->SetPad(0.66,0.33,1.0,0.0);
	TH1D* l_rej_ratio_flat = new TH1D(*cust_ip3d_llr.l_rej_flat);
	l_rej_ratio_flat->Divide(std_ip3d_llr.l_rej_flat);
	l_rej_ratio_flat->SetLineColor(kRed);
	l_rej_ratio_flat->SetMarkerColor(kRed);
	l_rej_ratio_flat->SetMarkerStyle(22);
	l_rej_ratio_flat->SetMarkerSize(0.5);
	l_rej_ratio_flat->GetXaxis()->SetTitleSize(0.1);
	l_rej_ratio_flat->GetXaxis()->SetLabelSize(0.1);
	l_rej_ratio_flat->GetXaxis()->SetRangeUser(0,2e6);//)LabelSize(0.1);
	l_rej_ratio_flat->GetYaxis()->SetTitleSize(0.1);
	l_rej_ratio_flat->GetYaxis()->SetTitleOffset(0.5);
	l_rej_ratio_flat->GetYaxis()->SetLabelSize(0.1);
	l_rej_ratio_flat->GetYaxis()->SetTitle("ratio");
	l_rej_ratio_flat->SetTitle("");
	l_rej_ratio_flat->Draw("PLE");
	gPad->SetGridy();

	TCanvas *c55 = new TCanvas("IP3D_LLR_pt","IP3D_LLr_pt",1200,400);

	std::vector<double> custLlr_b_mean, custLlr_b_rms, custLlr_l_mean, custLlr_l_rms;
	std::vector<double> stdLlr_b_mean, stdLlr_b_rms, stdLlr_l_mean, stdLlr_l_rms;

	for(int i=0; i<cust_ip3d_llr.nbin; i++){
		std::cout<<"MEAN: "<<std_ip3d_llr.b_ip3d_llr_ptbin[i]->GetMean()<<std::endl;
		std::cout<<"PT: "<<std_ip3d_llr.pt_bins[i]<<std::endl;
		custLlr_b_mean.push_back(cust_ip3d_llr.b_ip3d_llr_ptbin[i]->GetMean());
		custLlr_b_rms.push_back(cust_ip3d_llr.b_ip3d_llr_ptbin[i]->GetRMS());
		custLlr_l_mean.push_back(cust_ip3d_llr.l_ip3d_llr_ptbin[i]->GetMean());
		custLlr_l_rms.push_back(cust_ip3d_llr.l_ip3d_llr_ptbin[i]->GetRMS());
	}

	TGraphErrors *b_cust_ip3d_llr_pt_graph = new TGraphErrors((int)cust_ip3d_llr.nbin,&(cust_ip3d_llr.pt_bins[0]), &(custLlr_b_mean[0]), &(cust_ip3d_llr.pt_bin_width[0]), &(custLlr_b_rms[0]));
	TGraphErrors *l_cust_ip3d_llr_pt_graph = new TGraphErrors((int)cust_ip3d_llr.nbin,&(cust_ip3d_llr.pt_bins[0]), &(custLlr_l_mean[0]), &(cust_ip3d_llr.pt_bin_width[0]), &(custLlr_l_rms[0]));
	TGraph *cust_flat_cut = new TGraph(cust_ip3d_llr.nbin,&(cust_ip3d_llr.pt_bins[0]), &(cust_ip3d_llr.v_flat_cut[0]));

	for(int i=0; i<std_ip3d_llr.nbin; i++){
		stdLlr_b_mean.push_back(std_ip3d_llr.b_ip3d_llr_ptbin[i]->GetMean());
		stdLlr_b_rms.push_back(std_ip3d_llr.b_ip3d_llr_ptbin[i]->GetRMS());
		stdLlr_l_mean.push_back(std_ip3d_llr.l_ip3d_llr_ptbin[i]->GetMean());
		stdLlr_l_rms.push_back(std_ip3d_llr.l_ip3d_llr_ptbin[i]->GetRMS());
	}

	TGraphErrors *b_std_ip3d_llr_pt_graph = new TGraphErrors(std_ip3d_llr.nbin,&(std_ip3d_llr.pt_bins[0]), &(stdLlr_b_mean[0]), &(std_ip3d_llr.pt_bin_width[0]), &(stdLlr_b_rms[0]));
	TGraphErrors *l_std_ip3d_llr_pt_graph = new TGraphErrors(std_ip3d_llr.nbin,&(std_ip3d_llr.pt_bins[0]), &(stdLlr_l_mean[0]), &(std_ip3d_llr.pt_bin_width[0]), &(stdLlr_l_rms[0]));
	TGraph *std_flat_cut = new TGraph(std_ip3d_llr.nbin,&(std_ip3d_llr.pt_bins[0]), &(std_ip3d_llr.v_flat_cut[0]));
	c55->Divide(2,1);
	c55->cd(1);
	//TMultiGraph* cust_multigraph = new TMultiGraph();
	b_cust_ip3d_llr_pt_graph->SetFillColor(kRed);
	l_cust_ip3d_llr_pt_graph->SetFillColor(kBlue);
	b_cust_ip3d_llr_pt_graph->SetFillStyle(3004);
	l_cust_ip3d_llr_pt_graph->SetFillStyle(3003);
	b_cust_ip3d_llr_pt_graph->Draw("api2");
	b_cust_ip3d_llr_pt_graph->GetYaxis()->SetRangeUser(-12,30);
	b_cust_ip3d_llr_pt_graph->GetYaxis()->SetTitle("IP3D_LLR");
	b_cust_ip3d_llr_pt_graph->GetXaxis()->SetTitle("jet_pt [MeV]");
	
	l_cust_ip3d_llr_pt_graph->Draw("pi2");
	cust_flat_cut->Draw("l");
	TLegend *leg_cust = new TLegend(0.5,0.7,0.9,0.9);
	leg_cust->SetBorderSize(0);
	leg_cust->SetFillStyle(0);
	leg_cust->AddEntry(b_cust_ip3d_llr_pt_graph,"b-jet IP3D custom","F");
	leg_cust->AddEntry(l_cust_ip3d_llr_pt_graph,"light-jet IP3D custom","F");
	leg_cust->AddEntry(cust_flat_cut,"Flat Cut","l");
	leg_cust->Draw();
	//cust_multigraph->Draw("a");
	c55->cd(2);
	//TMultiGraph* cust_multigraph = new TMultiGraph();
	b_std_ip3d_llr_pt_graph->SetFillColor(kRed);
	l_std_ip3d_llr_pt_graph->SetFillColor(kBlue);
	b_std_ip3d_llr_pt_graph->SetFillStyle(3004);
	l_std_ip3d_llr_pt_graph->SetFillStyle(3003);
	b_std_ip3d_llr_pt_graph->Draw("api2");
	b_std_ip3d_llr_pt_graph->GetYaxis()->SetRangeUser(-12,30);
	b_std_ip3d_llr_pt_graph->GetYaxis()->SetTitle("IP3D_LLR");
	b_std_ip3d_llr_pt_graph->GetXaxis()->SetTitle("jet_pt [MeV]");
	l_std_ip3d_llr_pt_graph->Draw("pi2");
	std_flat_cut->Draw("l");
	TLegend *leg_std = new TLegend(0.5,0.7,0.9,0.9);
	leg_std->SetBorderSize(0);
	leg_std->SetFillStyle(0);
	leg_std->AddEntry(b_std_ip3d_llr_pt_graph,"b-jet IP3D std","F");
	leg_std->AddEntry(l_std_ip3d_llr_pt_graph,"light-jet IP3D std","F");
	leg_std->AddEntry(std_flat_cut,"Flat Cut","l");
	leg_std->Draw();
	//cust_multigraph->Draw("a");
	TFile *output_file = new TFile("IP3D_ROC_new.root","update");
	ROC_graph->Write();
	c_ROC_graph->Write();
	output_file->Close();
	TFile *output_sum = new TFile("IP3D_summary_new.root","recreate");
	ROC_graph->Write();
	c_ROC_graph->Write();

	B_nPixHitsvsnSiHits->Write();
	B_IP3D_nPixHitsvsnSiHits->Write();

	bHpt_vs_jet_pt->Write();
	bHpt_vs_jet_pt->ProfileX("_pfx",1,-1,"s")->Write();
	TH1D* bHpt_ratio_jetpt = new TH1D("bHpt_ratio_jetpt","bHpt_ratio_jetpt;jet_pt (MeV);ratio", 300, 0, 3e6);
	for(int k=1;k<=bHpt_vs_jet_pt->ProfileX("_pfx",1,-1,"s")->GetNbinsX();k++) bHpt_ratio_jetpt->SetBinContent(k,bHpt_vs_jet_pt->ProfileX("_pfx",1,-1,"s")->GetBinContent(k)/(bHpt_vs_jet_pt->ProfileX("_pfx",1,-1,"s")->GetBinCenter(k)+1));
	bHpt_ratio_jetpt->Write();
	bHptratio_vs_jet_pt->Write();
	bHptratio_vs_jet_pt->ProfileX("_pfx",1,-1,"s")->Write();

	std::cout<<"Number of b+c tracks:\t"<<n_b_tracks<<std::endl;
	std::cout<<"Number of b+c selected tracks:\t"<<n_b_sel_tracks<<std::endl;
	std::cout<<"Lead Cut efficiency:\t"<<((double)(n_b_sel_tracks)/(double)(n_b_tracks))<<std::endl;



	THStack *frag_stack_grades_jpt = new THStack("frag_stack_grades_jpt","frag_stack_grades_jpt; jet_pt (MeV)");
	THStack *frag_stack_grades_bHLxy = new THStack("frag_stack_grades_bH_Lxy","frag_stack_grades_bH_Lxy; bH_Lxy (mm)");
	THStack *BC_stack_grades_jpt = new THStack("BC_stack_grades_jpt","BC_stack_grades_jpt; jet_pt (MeV)");
	THStack *BC_stack_grades_bHLxy = new THStack("BC_stack_grades_bH_Lxy","BC_stack_grades_bH_Lxy; bH_Lxy (mm)");

	THStack *ratio_frag_stack_grades_jpt = new THStack("ratio_frag_stack_grades_jpt","ratio_frag_stack_grades_jpt; jet_pt (MeV)");
	THStack *ratio_frag_stack_grades_bHLxy = new THStack("ratio_frag_stack_grades_bH_Lxy","ratio_frag_stack_grades_bH_Lxy; bH_Lxy (mm)");
	THStack *ratio_BC_stack_grades_jpt = new THStack("ratio_BC_stack_grades_jpt","ratio_BC_stack_grades_jpt; jet_pt (MeV)");
	THStack *ratio_BC_stack_grades_bHLxy = new THStack("ratio_BC_stack_grades_bH_Lxy","ratio_BC_stack_grades_bH_Lxy; bH_Lxy (mm)");

	frag_grades.getRatioHistos();
	BC_grades.getRatioHistos();

	for(int grade=0; grade<14; grade++){

		frag_stack_grades_jpt->Add(frag_grades.grades_jpt[grade]);
		frag_stack_grades_bHLxy->Add(frag_grades.grades_bHLxy[grade]);
		BC_stack_grades_jpt->Add(BC_grades.grades_jpt[grade]);
		BC_stack_grades_bHLxy->Add(BC_grades.grades_bHLxy[grade]);

		ratio_frag_stack_grades_jpt->Add(frag_grades.ratio_grades_jpt[grade]);
		ratio_frag_stack_grades_bHLxy->Add(frag_grades.ratio_grades_bHLxy[grade]);
		ratio_BC_stack_grades_jpt->Add(BC_grades.ratio_grades_jpt[grade]);
		ratio_BC_stack_grades_bHLxy->Add(BC_grades.ratio_grades_bHLxy[grade]);
	}


	frag_stack_grades_jpt->Write();
	frag_stack_grades_bHLxy->Write();

	BC_d0sig->Write();
	BC_z0sig->Write();
	BC_dz0sig->Write();
	BC_trkpt->Write();
	BC_d0sig_vs_z0sig->Write();
	BC_d0sig_custom->Write();
	BC_z0sig_custom->Write();
	BC_d0sig_vs_z0sig_custom->Write();
	frag_d0sig->Write();
	frag_z0sig->Write();
	frag_dz0sig->Write();
	frag_trkpt->Write();
	frag_d0sig_vs_z0sig->Write();
	frag_d0sig_custom->Write();
	frag_z0sig_custom->Write();
	frag_d0sig_vs_z0sig_custom->Write();
	BC_sumtrkpt_jetpt->Write();
	BC_sumtrkpt_jetpt->ProfileX("_pfx",1,-1,"s")->Write();
	BC_sumtrkpt_bHpt->Write();
	BC_sumtrkpt_bHpt->ProfileX("_pfx",1,-1,"s")->Write();
	frag_sumtrkpt_jetpt->Write();
	frag_sumtrkpt_jetpt->ProfileX("_pfx",1,-1,"s")->Write();
	frag_sumtrkpt_bHpt->Write();
	frag_sumtrkpt_bHpt->ProfileX("_pfx",1,-1,"s")->Write();
	TH1D* BC_sumtrkpt_ratio_jetpt = new TH1D("BC_sumtrkpt_ratio_jetpt","BC_sumtrkpt_ratio_jetpt;jet_pt (MeV);ratio", 300, 0, 3e6);
	for(int k=1;k<=BC_sumtrkpt_jetpt->ProfileX("_pfx",1,-1,"s")->GetNbinsX();k++) BC_sumtrkpt_ratio_jetpt->SetBinContent(k,BC_sumtrkpt_jetpt->ProfileX("_pfx",1,-1,"s")->GetBinContent(k)/(BC_sumtrkpt_jetpt->ProfileX("_pfx",1,-1,"s")->GetBinCenter(k)+1));
	BC_sumtrkpt_ratio_jetpt->Write();
	TH1D* BC_sumtrkpt_ratio_bHpt = new TH1D("BC_sumtrkpt_ratio_bHpt","BC_sumtrkpt_ratio_bHpt;bH_pt (MeV);ratio", 300, 0, 3e6);
	for(int k=1;k<=BC_sumtrkpt_bHpt->ProfileX("_pfx",1,-1,"s")->GetNbinsX();k++) BC_sumtrkpt_ratio_bHpt->SetBinContent(k,BC_sumtrkpt_bHpt->ProfileX("_pfx",1,-1,"s")->GetBinContent(k)/(BC_sumtrkpt_bHpt->ProfileX("_pfx",1,-1,"s")->GetBinCenter(k)+1));
	BC_sumtrkpt_ratio_bHpt->Write();
	BC_stack_grades_jpt->Write();
	TH1D* frag_sumtrkpt_ratio_jetpt = new TH1D("frag_sumtrkpt_ratio_jetpt","frag_sumtrkpt_ratio_jetpt;jet_pt (MeV);ratio", 300, 0, 3e6);
	for(int k=1;k<=frag_sumtrkpt_jetpt->ProfileX("_pfx",1,-1,"s")->GetNbinsX();k++) frag_sumtrkpt_ratio_jetpt->SetBinContent(k,frag_sumtrkpt_jetpt->ProfileX("_pfx",1,-1,"s")->GetBinContent(k)/(frag_sumtrkpt_jetpt->ProfileX("_pfx",1,-1,"s")->GetBinCenter(k)+1));
	frag_sumtrkpt_ratio_jetpt->Write();
	TH1D* frag_sumtrkpt_ratio_bHpt = new TH1D("frag_sumtrkpt_ratio_bHpt","frag_sumtrkpt_ratio_bHpt;bH_pt (MeV);ratio", 300, 0, 3e6);
	for(int k=1;k<=frag_sumtrkpt_bHpt->ProfileX("_pfx",1,-1,"s")->GetNbinsX();k++) frag_sumtrkpt_ratio_bHpt->SetBinContent(k,frag_sumtrkpt_bHpt->ProfileX("_pfx",1,-1,"s")->GetBinContent(k)/(frag_sumtrkpt_bHpt->ProfileX("_pfx",1,-1,"s")->GetBinCenter(k)+1));
	frag_sumtrkpt_ratio_bHpt->Write();
	BC_stack_grades_bHLxy->Write();
	BC_leadOrder_vs_tvSize->Write();
	BC_leadOrder_vs_tvSize->ProfileX("_pfx",1,-1,"s")->Write();
	BC_leadOrder_vs_jetPt->Write();
	BC_leadOrder_vs_jetPt->ProfileX("_pfx",1,-1,"s")->Write();
	BC_leadOrder_vs_bHLxy->Write();
	BC_leadOrder_vs_bHLxy->ProfileX("_pfx",1,-1,"s")->Write();
	BC_tvSize_vs_jetPt->Write();
	BC_tvSize_vs_jetPt->ProfileX("_pfx",1,-1,"s")->Write();
	BC_tvSize_vs_bHLxy->Write();
	BC_tvSize_vs_bHLxy->ProfileX("_pfx",1,-1,"s")->Write();
	ratio_frag_stack_grades_jpt->Write();
	ratio_frag_stack_grades_bHLxy->Write();
	ratio_BC_stack_grades_jpt->Write();
	ratio_BC_stack_grades_bHLxy->Write();
	B_IP3DtrkVsJetPt->Write();
	BC_IP3DtrkVsJetPt->Write();
	BC_IP3DtrkVsJetPt->ProfileX("_pfx",1,-1,"s")->Write();
	BC_trkVsJetPt->Write();
	BC_trkVsJetPt->ProfileX("_pfx",1,-1,"s")->Write();
	BC_IP3DtrkVsbHLxy->Write();
	BC_IP3DtrkVsbHLxy->ProfileX("_pfx",1,-1,"s")->Write();
	BC_trkVsbHLxy->Write();
	BC_trkVsbHLxy->ProfileX("_pfx",1,-1,"s")->Write();
	BC_trkpt_vs_jetpt->Write();
	BC_trkpt_vs_jetpt->ProfileX("_pfx",1,-1,"s")->Write();
	BC_trkpt_vs_jetpt_log->Write();
	BC_trkpt_vs_jetpt_log->ProfileX("_pfx",1,-1,"s")->Write();
	TProfile* prof = B_IP3DtrkVsJetPt->ProfileX("_pfx",1,-1,"s");
	Btruth_tracks_vs_jetpt->ProfileX("_pfx",1,-1,"s")->Write();
	Btruth_tracks_vs_jetpt->Write();
	Btruth_tracks_vs_bHLxy->ProfileX("_pfx",1,-1,"s")->Write();
	Btruth_tracks_vs_bHLxy->Write();
	B_IP3Dtrk->Write();
	B_dRjetVsjetpt->Write();
	frag_trkpt_vs_jetpt->Write();
	frag_trkpt_vs_jetpt->ProfileX("_pfx",1,-1,"s")->Write();
	frag_trkpt_vs_jetpt_log->Write();
	frag_trkpt_vs_jetpt_log->ProfileX("_pfx",1,-1,"s")->Write();
	prof->Write();
	frag_IP3DtrkVsJetPt->Write();
	B_dRjetVsIP3Dtrk->Write();
	TProfile* prof_frag = frag_IP3DtrkVsJetPt->ProfileX("_pfx",1,-1,"s");
	prof_frag->Write();
	B_dRjetVsjetpt_sel->Write();
	BC_d0vsjetpt->Write();
	BC_d0sigvsjetpt->Write();
	B_z0sinthetavsjetpt->Write();
	BC_d0vsjetpt_ip3d_sel->Write();
	B_z0sinthetavsjetpt_ip3d_sel->Write();

	BC_dRtrkjetVsjetpt->Write();
	frag_dRtrkjetVsjetpt->Write();
	light_dRtrkjetVsjetpt->Write();
	BC_dRtrkjetVsjetpt->ProfileX("_pfx",1,-1,"s")->Write();
	frag_dRtrkjetVsjetpt->ProfileX("_pfx",1,-1,"s")->Write();
	light_dRtrkjetVsjetpt->ProfileX("_pfx",1,-1,"s")->Write();
	

	B_1pixHitIP3D_etaphimap->Write();
	B_0pixHit_etaphimap->Write();	
	B_1pixHit_etaphimap->Write();	
	frag_0pixHit_etaphimap->Write();
	frag_UT_npixVsjetpt->Write();
	B_npixVsjetpt->Write();
	B_npixVsjetpt->ProfileX("_pfx",1,-1,"s")->Write();
	B_npixVsbH_Lxy->Write();
	B_npixVsbH_Lxy->ProfileX("_pfx",1,-1,"s")->Write();
	B_IP3D_npixVsjetpt->Write();
	B_IP3D_npixVsjetpt->ProfileX("_pfx",1,-1,"s")->Write();
	B_IP3D_npixVsbH_Lxy->Write();
	B_IP3D_npixVsbH_Lxy->ProfileX("_pfx",1,-1,"s")->Write();
	B_trkVtxVsjetpt->Write();
	B_trkVtxVsjetpt->ProfileX("_pfx",1,-1,"s")->Write();
	B_UT_npixVsjetpt->Write();
	B_UT_npixVsjetpt->ProfileX("_pfx",1,-1,"s")->Write();
	B_UT_npixVsbH_Lxy->Write();
	B_UT_npixVsbH_Lxy->ProfileX("_pfx",1,-1,"s")->Write();
	B_UTnoC_npixVsjetpt->Write();
	B_UTnoC_npixVsjetpt->ProfileX("_pfx",1,-1,"s")->Write();
	B_UTnoC_npixVsbH_Lxy->Write();
	B_UTnoC_npixVsbH_Lxy->ProfileX("_pfx",1,-1,"s")->Write();
	B_UT_trkptVsjetpt->Write();//Fill(jet_pt->at(i),trk_pt[k]);
	B_UT_d0Vsjetpt->Write();//Fill(jet_pt->at(i), d0[k]);	
	B_UT_z0sinthetaVsjetpt->Write();//Fill(jet_pt->at(i), z0[k] * sin(theta[k]));
	B_UT_trkVtxVsjetpt->Write();//Fill(jet_pt->at(i), sqrt(pow(vtx_X[k],2)+pow(vtx_Y[k],2)));	
	B_UT_trkVtxVsjetpt->ProfileX("_pfx",1,-1,"s")->Write();//Fill(jet_pt->at(i), sqrt(pow(vtx_X[k],2)+pow(vtx_Y[k],2)));	
	B_UT_ip3dgradeVsjetpt->Write();//Fill(jet_pt->at(i), (jet_trk_ip3d_grade->at(i))[k]);
	B_bHLxyVsJetPt->Write();
	B_bHLxyVsJetPt->ProfileX("_pfx",1,-1,"s")->Write();
	B_bHLxyVsJetPt_01->Write();
	B_bHLxyVsJetPt_01->ProfileX("_pfx",1,-1,"s")->Write();
	B_bHLxyVsJetPt_gt1->Write();
	B_bHLxyVsJetPt_gt1->ProfileX("_pfx",1,-1,"s")->Write();

	BC_grades_histo->Scale(1/BC_grades_histo->Integral());
	//BC_grades_histo->Draw("goff");
	

	std::vector<std::string> tg = get_track_grades();
	for(int i=0; i<=14; i++)
		BC_grades_histo->GetXaxis()->SetBinLabel(i,tg[i-1].c_str());
	for(int i=0; i<14; i++){
		BC_grades.grades_d0sig[i]->Scale(1/(BC_grades.grades_d0sig[i]->Integral()));
		BC_grades.grades_d0sig[i]->SetLineColor(kRed);
		Bjet_grades.grades_d0sig[i]->Scale(1/(Bjet_grades.grades_d0sig[i]->Integral()));
		Bjet_grades.grades_d0sig[i]->SetLineColor(kBlack);
		frag_grades.grades_d0sig[i]->Scale(1/(frag_grades.grades_d0sig[i]->Integral()));
		frag_grades.grades_d0sig[i]->SetLineColor(kGreen+2);
		light_grades.grades_d0sig[i]->Scale(1/(light_grades.grades_d0sig[i]->Integral()));
		light_grades.grades_d0sig[i]->SetLineColor(kBlue+2);
		THStack* grade_d0sig = new THStack(("d0sig_"+(std::string)tg[i]).c_str(), ("d0sig_"+(std::string)tg[i]+";d0sig;").c_str());
		grade_d0sig->Add(BC_grades.grades_d0sig[i]);
		grade_d0sig->Add(Bjet_grades.grades_d0sig[i]);
		grade_d0sig->Add(frag_grades.grades_d0sig[i]);
		grade_d0sig->Add(light_grades.grades_d0sig[i]);
		grade_d0sig->Write();
	}
	BC_grades_histo->Write();
	h_b_jet_pt->Scale(1./(double)h_b_jet_pt->Integral());	
	h_l_jet_pt->Scale(1./(double)h_l_jet_pt->Integral());	
	h_b_jet_pt->Write();
	h_l_jet_pt->Write();
	c2->Write();
	c3->Write();
	c55->Write();
	output_sum->Close();
}


//TH1D* flat_efficiency_plot(ip3d_llr_ptbin, int nbin = 20, double eff = 0.7, double max_pt=3e+6){

//	std::vector<TH1D*> v_light_llr, v_b_llr;
//	for(int i = 0; i<nbin; i++){
//		for
//	}
//
//}

	

