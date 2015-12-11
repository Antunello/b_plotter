#ifndef IPxDSTANDALONETOOL_H
#define IPxDSTANDALONETOOL_H

#include <vector>
#include <map>
#include <string>

class TFile;
class TH2;
class TH1;

class IPxDStandaloneTool {

public: 
  IPxDStandaloneTool();
  ~IPxDStandaloneTool();

  int initTrainingMode(std::string fname, std::string jetColl, std::vector<std::string> grades);
  void fillTrack(int grade, float d0sig, float z0sig, int hypothesis);
  void fillTrack(int nD, int grade, float d0sig, float z0sig, int hypothesis);
  void write();

  int initEvaluationMode(std::string fname, std::string jetColl, std::vector<std::string> grades);
  float getTrackLLR(int nD, int grade, float d0sig, float z0sig, int hypoS, int hypoB);
  float getTrackProb(int nD, int grade, float d0sig, float z0sig, int hypothesis);
  float getJetLLR(int nD, std::vector<int> grade, std::vector<float> d0sig, std::vector<float> z0sig, int hypoS, int hypoB);
  float getJetProb(int nD, std::vector<int> grade, std::vector<float> d0sig, std::vector<float> z0sig, int hypothesis);

private:

  struct Hypothesis {
    std::string name;
    TH1 *h_2D;
    TH2 *h_3D;
    double I_2D;
    double I_3D;
  };

  struct Training {
    std::string name;
    std::map<int, Hypothesis *>h; // 0: L - 4: C - 5: B
  };

  std::vector<Training*> m_trainings;

  TFile *m_file;

  int m_runMode;
};


#endif
