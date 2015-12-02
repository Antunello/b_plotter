#ifndef IPxDTRAININGTOOL_H
#define IPxDTRAININGTOOL_H

#include <TH2.h>
#include <TString.h>

#include <vector>
#include <map>

class IPxDTrainingTool {

public: 
  IPxDTrainingTool();
  ~IPxDTrainingTool();

  void initTrainingMode(int nGrades);
  void fillTrack(int grade, float d0sig, float z0sig, int hypothesis);
  void write(TString fname);

  void initEvaluationMode(int nGrades, TString fname);
  void initEvaluationModeDB(TString fname, TString jetColl, std::vector<TString> grades);
  float getTrackProb(int nD, int grade, float d0sig, float z0sig, int hypothesis);
  std::vector<float> getTrackProb(int nD, int grade, float d0sig, float z0sig);
  float getTrackLLR(int nD, int grade, float d0sig, float z0sig, int hypoS, int hypoB);
  float getJetLLR(int nD, std::vector<int> grade, std::vector<float> d0sig, std::vector<float> z0sig, int hypoS, int hypoB);

private:

  struct Hypothesis {
    TString name;
    TH1 *h_2D;
    TH2 *h_3D;
  };

  struct Training {
    TString name;
    std::map<int, Hypothesis *>h; // 0: L - 4: C - 5: B
  };

  std::vector<Training*> m_trainings;

  int m_runMode;
};


#endif
