#ifndef JetLearning_JetLearningAnalysis_H
#define JetLearning_JetLearningAnalysis_H

#include <EventLoop/Algorithm.h>
// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"

#ifndef __CINT__
#include "xAODCore/ShallowCopy.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODCaloEvent/CaloClusterContainer.h"
#include "xAODJet/JetConstituentVector.h"
#endif

//JetReclustering
#include <JetReclustering/JetReclusteringTool.h>

#include "InDetTrackSelectionTool/InDetTrackSelectionTool.h"

#include <TTree.h>

class JetLearningAnalysis : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;
  bool m_debug = false;
  bool m_doLC = false;
  bool m_doClusterInfo = false;
  bool m_doJetReclustering = false;
  bool m_doTracks = false;

  InDet::InDetTrackSelectionTool* m_trackSelectionTool; //!
  JetReclusteringTool* m_jetReclusteringTool; //!

  std::string m_jets = "AntiKt4NoAreaJets",
              m_truth_jets = "AntiKt4TruthJets",
              m_clusts = "CaloCalTopoClusters",
              m_vertices       = "PrimaryVertices";

  int m_mu;
  int m_NPV;
  float m_rho;
  float m_sigmarho;
  int m_eventNumber;
  float m_eventWeight;

  std::vector<float> m_j0pt; //!
  std::vector<float> m_jnoarea0pt; //!
  std::vector<float> m_j0eta; //!
  std::vector<float> m_j0phi; //!
  std::vector<float> m_j0m; //!
  std::vector<float> m_j0area; //!
  std::vector<bool> m_j0isPU; //!
  std::vector<float> m_j0Rpt; //!
  std::vector<float> m_j0deltaRpt; //!
  std::vector<float> m_j0JVF; //!
  std::vector<float> m_j0corrJVF; //!
  std::vector<float> m_j0JVT; //!

  std::vector<std::vector<float> > m_j0_trkpt; //!
  std::vector<std::vector<float> > m_j0_trketa; //!
  std::vector<std::vector<float> > m_j0_trkphi; //!
  std::vector<std::vector<float> > m_j0_trkm; //!
  std::vector<std::vector<bool> > m_j0_trkisHS; //!
  std::vector<std::vector<bool> > m_j0_trkisPU; //!

  std::vector<std::vector<float> > m_j0_clpt; //!
  std::vector<std::vector<float> > m_j0_cleta; //!
  std::vector<std::vector<float> > m_j0_clphi; //!
  std::vector<std::vector<float> > m_j0_clm; //!
  std::vector<std::vector<float> > m_j0_cllatwidth; //!
  std::vector<std::vector<float> > m_j0_cllongwidth; //!

  std::vector<float> m_j0sumpt5; //!
  std::vector<float> m_j0sumpt10; //!
  std::vector<float> m_j0sumpt15; //!
  std::vector<float> m_j0sumpt20; //!
  std::vector<float> m_j0sumpt25; //!
  std::vector<float> m_j0sumpt30; //!
  std::vector<float> m_j0sumpt35; //!
  std::vector<float> m_j0sumpt40; //!

  std::vector<float> m_clpt; //!
  std::vector<float> m_cleta; //!
  std::vector<float> m_clphi; //!
  std::vector<float> m_clm; //!
  std::vector<float> m_cllatwidth; //!
  std::vector<float> m_cllongwidth; //!

  std::vector<float> m_tj0pt; //!
  std::vector<float> m_tj0eta; //!
  std::vector<float> m_tj0phi; //!
  std::vector<float> m_tj0m; //!
  std::vector<float> m_tj0mindr; //!

private:
  xAOD::TEvent *m_event; //!
  xAOD::TStore *m_store;  //!
  TTree *m_tree; //!

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!


  EL::StatusCode SetRho(DataVector<xAOD::Jet_v1> jets,float event_rho);
  EL::StatusCode FindRho(const xAOD::CaloClusterContainer* in_clusters,float& rho, float& sigma);
  EL::StatusCode SetMinDR(DataVector<xAOD::Jet_v1> jets,float threshold=5000);
  EL::StatusCode FindTruthMatch(DataVector<xAOD::Jet_v1> jets, DataVector<xAOD::Jet_v1> tjets);
  EL::StatusCode CalcAnnulus(const xAOD::Jet_v1*);

  // this is a standard constructor
  JetLearningAnalysis ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  // this is needed to distribute the algorithm to the workers
  ClassDef(JetLearningAnalysis, 1);
};

#endif
