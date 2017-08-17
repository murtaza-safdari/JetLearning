#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <JetLearning/JetLearningAnalysis.h>
#include "xAODEventInfo/EventInfo.h"
#include "EventLoop/OutputStream.h"

//Tracking
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "InDetTrackSelectionTool/InDetTrackSelectionTool.h"

//JetReclustering
#include <JetReclustering/JetReclusteringTool.h>

//FastJet
#include <xAODJet/FastJetLinkBase.h>
#include "JetInterface/IPseudoJetGetter.h"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include <fastjet/tools/Subtractor.hh>
#include "fastjet/PseudoJet.hh"

//CaloClusters
#include "xAODCaloEvent/CaloClusterContainer.h"
#include "xAODCaloEvent/CaloClusterChangeSignalState.h"

//Jets
#include "xAODJet/JetAttributes.h"
#include "xAODBase/IParticleHelpers.h"
#include <xAODBase/IParticleContainer.h>
#include <xAODBase/IParticle.h>

//Truth
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthEvent.h"
/*#include "xAODTruth/TruthVertexContainer.h"
#include "xAODTruth/TruthPileupEventContainer.h"*/

//xAOD stuff
#include "xAODCore/ShallowAuxContainer.h"
#include "xAODJet/JetConstituentVector.h"
#include "FourMomUtils/xAODP4Helpers.h"

// xAH includes
#include "xAODAnaHelpers/HelperFunctions.h"
//#include "xAODAnaHelpers/tools/ReturnCheck.h"
namespace HF = HelperFunctions;

#define ARRAY_INIT {}
// this is needed to distribute the algorithm to the workers
ClassImp(JetLearningAnalysis)

JetLearningAnalysis :: JetLearningAnalysis () :
  m_tree(new TTree("oTree", "output tree")),
  m_trackSelectionTool{nullptr},
  m_jetReclusteringTool{nullptr},
  m_mu(-999),
  m_rho(-99),
  m_sigmarho(-99),
  m_NPV(-99),
  m_eventNumber(-999),
  m_eventWeight(-999),
  m_j0pt{ARRAY_INIT},
  m_jnoarea0pt{ARRAY_INIT},
  m_j0eta{ARRAY_INIT},
  m_j0phi{ARRAY_INIT},
  m_j0m{ARRAY_INIT},
  m_j0area{ARRAY_INIT},
  m_j0isPU{ARRAY_INIT},
  m_j0Rpt{ARRAY_INIT},
  m_j0deltaRpt{ARRAY_INIT},
  m_j0JVF{ARRAY_INIT},
  m_j0corrJVF{ARRAY_INIT},
  m_j0JVT{ARRAY_INIT},
  m_j0_clpt{ARRAY_INIT},
  m_j0_cleta{ARRAY_INIT},
  m_j0_clphi{ARRAY_INIT},
  m_j0_clm{ARRAY_INIT},
  m_j0_cllatwidth{ARRAY_INIT},
  m_j0_cllongwidth{ARRAY_INIT},
  m_j0_trkpt{ARRAY_INIT},
  m_j0_trketa{ARRAY_INIT},
  m_j0_trkphi{ARRAY_INIT},
  m_j0_trkm{ARRAY_INIT},
  m_j0_trkisHS{ARRAY_INIT},
  m_j0_trkisPU{ARRAY_INIT},
  m_j0sumpt5{ARRAY_INIT},
  m_j0sumpt10{ARRAY_INIT},
  m_j0sumpt15{ARRAY_INIT},
  m_j0sumpt20{ARRAY_INIT},
  m_j0sumpt25{ARRAY_INIT},
  m_j0sumpt30{ARRAY_INIT},
  m_j0sumpt35{ARRAY_INIT},
  m_j0sumpt40{ARRAY_INIT},
  m_clpt{ARRAY_INIT},
  m_cleta{ARRAY_INIT},
  m_clphi{ARRAY_INIT},
  m_clm{ARRAY_INIT},
  m_cllatwidth{ARRAY_INIT},
  m_cllongwidth{ARRAY_INIT},
  m_tj0pt{ARRAY_INIT},
  m_tj0eta{ARRAY_INIT},
  m_tj0phi{ARRAY_INIT},
  m_tj0m{ARRAY_INIT},
  m_tj0mindr{ARRAY_INIT}
{}


EL::StatusCode JetLearningAnalysis :: setupJob (EL::Job& job)
{
  job.useXAOD();
  //xAOD::Init(); // call before opening first file
  xAOD::Init("JetLearningAnalysis").ignore(); // call before opening first file

  EL::OutputStream out ("outputTree");
  job.outputAdd (out);
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode JetLearningAnalysis :: histInitialize ()
{
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode JetLearningAnalysis :: fileExecute ()
{
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode JetLearningAnalysis :: changeInput (bool firstFile)
{
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode JetLearningAnalysis :: initialize ()
{
  m_trackSelectionTool = new InDet::InDetTrackSelectionTool("TrackSelection");
  ANA_CHECK(m_trackSelectionTool->setProperty("CutLevel", "Loose" )); // set tool to apply the pre-defined "Loose" cuts        
  ANA_CHECK(m_trackSelectionTool->setProperty("maxAbsEta", 4.5 )); // set |eta| cut
  ANA_CHECK(m_trackSelectionTool->initialize());

  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();

  if(m_doJetReclustering){
    m_jetReclusteringTool = new JetReclusteringTool(m_jets+std::to_string(std::rand()));
    ANA_CHECK(m_jetReclusteringTool->setProperty("InputJetContainer",  m_clusts));
    ANA_CHECK(m_jetReclusteringTool->setProperty("OutputJetContainer", m_jets));
    ANA_CHECK(m_jetReclusteringTool->setProperty("ReclusterRadius",    0.4));
    ANA_CHECK(m_jetReclusteringTool->setProperty("ReclusterAlgorithm", "AntiKt"));
    ANA_CHECK(m_jetReclusteringTool->setProperty("InputJetPtMin",      0));
    ANA_CHECK(m_jetReclusteringTool->setProperty("RCJetPtMin",         0.1));
    ANA_CHECK(m_jetReclusteringTool->setProperty("RCJetPtFrac",        0));
    ANA_CHECK(m_jetReclusteringTool->setProperty("DoArea",        true));
    ANA_CHECK(m_jetReclusteringTool->initialize());
  }

  TFile *file = wk()->getOutputFile ("outputTree");
  m_tree->SetDirectory (file);

  m_tree->Branch ("mu",              &m_mu, "mu/I");
  m_tree->Branch ("NPV",              &m_NPV, "NPV/I");
  m_tree->Branch ("rho",              &m_rho, "rho/F");
  m_tree->Branch ("sigmarho",              &m_sigmarho, "sigmarho/F");
  m_tree->Branch ("event_number",              &m_eventNumber, "mu/I");
  m_tree->Branch ("event_weight",              &m_eventWeight, "mu/F");

  m_tree->Branch("j0pt","std::vector<float>", &m_j0pt);
  m_tree->Branch("jnoarea0pt","std::vector<float>", &m_jnoarea0pt);
  m_tree->Branch("j0eta","std::vector<float>", &m_j0eta);
  m_tree->Branch("j0phi","std::vector<float>", &m_j0phi);
  m_tree->Branch("j0m","std::vector<float>", &m_j0m);
  m_tree->Branch("j0area","std::vector<float>", &m_j0area);
  m_tree->Branch("j0isPU","std::vector<bool>", &m_j0isPU);
  m_tree->Branch("j0Rpt","std::vector<float>", &m_j0Rpt);
  m_tree->Branch("j0deltaRpt","std::vector<float>", &m_j0deltaRpt);
  m_tree->Branch("j0JVF","std::vector<float>", &m_j0JVF);
  m_tree->Branch("j0corrJVF","std::vector<float>", &m_j0corrJVF);
  m_tree->Branch("j0JVT","std::vector<float>", &m_j0JVT);

  m_tree->Branch("j0_clpt","std::vector<std::vector<float> >", &m_j0_clpt);
  m_tree->Branch("j0_cleta","std::vector<std::vector<float> >", &m_j0_cleta);
  m_tree->Branch("j0_clphi","std::vector<std::vector<float> >", &m_j0_clphi);
  m_tree->Branch("j0_clm","std::vector<std::vector<float> >", &m_j0_clm);
  m_tree->Branch("j0_cllatwidth","std::vector<std::vector<float> >", &m_j0_cllatwidth);
  m_tree->Branch("j0_cllongwidth","std::vector<std::vector<float> >", &m_j0_cllongwidth);

  m_tree->Branch("j0_trkpt","std::vector<std::vector<float> >", &m_j0_trkpt);
  m_tree->Branch("j0_trketa","std::vector<std::vector<float> >", &m_j0_trketa);
  m_tree->Branch("j0_trkphi","std::vector<std::vector<float> >", &m_j0_trkphi);
  m_tree->Branch("j0_trkm","std::vector<std::vector<float> >", &m_j0_trkm);
  m_tree->Branch("j0_trkisHS","std::vector<std::vector<bool> >", &m_j0_trkisHS);
  m_tree->Branch("j0_trkisPU","std::vector<std::vector<bool> >", &m_j0_trkisPU);

  m_tree->Branch("j0sumpt5","std::vector<float>", &m_j0sumpt5);
  m_tree->Branch("j0sumpt10","std::vector<float>", &m_j0sumpt10);
  m_tree->Branch("j0sumpt15","std::vector<float>", &m_j0sumpt15);
  m_tree->Branch("j0sumpt20","std::vector<float>", &m_j0sumpt20);
  m_tree->Branch("j0sumpt25","std::vector<float>", &m_j0sumpt25);
  m_tree->Branch("j0sumpt30","std::vector<float>", &m_j0sumpt30);
  m_tree->Branch("j0sumpt35","std::vector<float>", &m_j0sumpt35);
  m_tree->Branch("j0sumpt40","std::vector<float>", &m_j0sumpt40);

  m_tree->Branch("clpt","std::vector<float>", &m_clpt);
  m_tree->Branch("cleta","std::vector<float>", &m_cleta);
  m_tree->Branch("clphi","std::vector<float>", &m_clphi);
  m_tree->Branch("clm","std::vector<float>", &m_clm);
  m_tree->Branch("cllatwidth","std::vector<float>", &m_cllatwidth);
  m_tree->Branch("cllongwidth","std::vector<float>", &m_cllongwidth);

  m_tree->Branch("tj0pt","std::vector<float>", &m_tj0pt);
  m_tree->Branch("tj0eta","std::vector<float>", &m_tj0eta);
  m_tree->Branch("tj0phi","std::vector<float>", &m_tj0phi);
  m_tree->Branch("tj0m","std::vector<float>", &m_tj0m);
  m_tree->Branch("tj0mindr","std::vector<float>", &m_tj0mindr);


  // as a check, let's see the number of events in our xAOD
  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode JetLearningAnalysis :: execute ()
{
  xAOD::TEvent* event = wk()->xaodEvent();
  //Event info
  std::string m_eventInfo = "EventInfo";
  const xAOD::EventInfo*                        eventInfo     (nullptr);
  if(!m_eventInfo.empty()) ANA_CHECK(HF::retrieve(eventInfo,    m_eventInfo,        m_event, m_store, m_debug));
  //TopoClusters info
  const xAOD::CaloClusterContainer*                        in_clusters     (nullptr);
  if(!m_clusts.empty())
    ANA_CHECK(HF::retrieve(in_clusters,     m_clusts,       m_event, m_store, m_debug));
  //Vertices
  const xAOD::VertexContainer*                  vertices      (nullptr);
  if(!m_vertices.empty()) ANA_CHECK(HF::retrieve(vertices,    m_vertices,       m_event, m_store, m_debug));
  //Tracks
  std::string m_tracks = "InDetTrackParticles";
  const xAOD::TrackParticleContainer*                  tracks      (nullptr);
  if(!m_tracks.empty()) ANA_CHECK(HF::retrieve(tracks,    m_tracks,       m_event, m_store, m_debug));
  //Truth jets
  const xAOD::JetContainer*                     truth_jets       (nullptr);
  if(!m_truth_jets.empty()) ANA_CHECK(HF::retrieve(truth_jets,    m_truth_jets,       m_event, m_store, m_debug));
  
  //LC Weights
  CaloClusterChangeSignalStateList stateHelperList;
  for(const auto clust: *in_clusters){
    if(m_doLC) stateHelperList.add(clust,xAOD::CaloCluster::State(1)); //default is calibrated but we can make it explicit anyway
    else stateHelperList.add(clust,xAOD::CaloCluster::State(0));
  }

  //Event Info
  m_mu = eventInfo->averageInteractionsPerCrossing();
  m_NPV = 0;
  for (auto *ivert : *vertices){
    if ((ivert)->nTrackParticles() >= 2) ++m_NPV;
  }
  m_eventNumber = eventInfo->eventNumber();
  m_eventWeight = eventInfo->mcEventWeight();

  float event_rho,event_sigma;
  FindRho(in_clusters,event_rho,event_sigma);
  m_rho = event_rho;
  m_sigmarho = event_sigma;

  //Jets
  if(SetMinDR(HF::sort_container_pt(truth_jets)) != EL::StatusCode::SUCCESS)
    Error("JetLearningAnalysis::execute()","Error in SetMinDR (only one jet in event?)");
  if(m_doJetReclustering) m_jetReclusteringTool->execute();

  //TopoClusters info
  const xAOD::JetContainer*                        jets     (nullptr);
  if(!m_jets.empty())
    ANA_CHECK(HF::retrieve(jets,     m_jets,       m_event, m_store, m_debug));
  if(FindTruthMatch(HF::sort_container_pt(jets), HF::sort_container_pt(truth_jets)) != EL::StatusCode::SUCCESS)
    Error("JetLearningAnalysis::execute()","Error in FindTruthMatch");

  m_j0pt.clear();
  m_jnoarea0pt.clear();
  m_j0eta.clear();
  m_j0phi.clear();
  m_j0m.clear();
  m_j0area.clear();
  m_j0isPU.clear();

  m_j0Rpt.clear();
  m_j0deltaRpt.clear();
  m_j0JVF.clear();
  m_j0corrJVF.clear();
  m_j0JVT.clear();

  m_j0_clpt.clear();
  m_j0_cleta.clear();
  m_j0_clphi.clear();
  m_j0_clm.clear();
  m_j0_cllatwidth.clear();
  m_j0_cllongwidth.clear();

  m_j0_trkpt.clear();
  m_j0_trketa.clear();
  m_j0_trkphi.clear();
  m_j0_trkm.clear();
  m_j0_trkisHS.clear();
  m_j0_trkisPU.clear();

  m_j0sumpt5.clear();
  m_j0sumpt10.clear();
  m_j0sumpt15.clear();
  m_j0sumpt20.clear();
  m_j0sumpt25.clear();
  m_j0sumpt30.clear();
  m_j0sumpt35.clear();
  m_j0sumpt40.clear();

  m_tj0pt.clear();
  m_tj0eta.clear();
  m_tj0phi.clear();
  m_tj0m.clear();
  m_tj0mindr.clear();

  // All clusters info
  m_clpt.clear();
  m_cleta.clear();
  m_clphi.clear();
  m_clm.clear();
  m_cllatwidth.clear();
  m_cllongwidth.clear();
  /*
  for(const auto clust: *in_clusters){
    m_clpt.push_back(clust->pt()*0.001);
    m_cleta.push_back(clust->eta());
    m_clphi.push_back(clust->phi());
    m_clm.push_back(clust->m());
  }
  */

  //PV
  const xAOD::Vertex *pv = 0;
  for(const auto& vx : *vertices) if(vx->vertexType()==1) {pv = vx;break;}
  int n_putracks = 0;
  if(m_doTracks){
    for (const auto& trk : *tracks) {
      if(m_trackSelectionTool->accept(*trk,pv) && trk->vertex() && trk->vertex()!=pv && trk->pt()<30e3) n_putracks++;
    }
  }
  if (!n_putracks) n_putracks++;

  DataVector<xAOD::Jet_v1> sorted_truth_jets = HF::sort_container_pt(truth_jets);
  static SG::AuxElement::ConstAccessor< int > truth_match_i("truth_match_i");
  static SG::AuxElement::ConstAccessor< float > minDR("minDR");
  static SG::AuxElement::ConstAccessor< float > sumpT5("sumpT5");
  static SG::AuxElement::ConstAccessor< float > sumpT10("sumpT10");
  static SG::AuxElement::ConstAccessor< float > sumpT15("sumpT15");
  static SG::AuxElement::ConstAccessor< float > sumpT20("sumpT20");
  static SG::AuxElement::ConstAccessor< float > sumpT25("sumpT25");
  static SG::AuxElement::ConstAccessor< float > sumpT30("sumpT30");
  static SG::AuxElement::ConstAccessor< float > sumpT35("sumpT35");
  static SG::AuxElement::ConstAccessor< float > sumpT40("sumpT40");

  for(const auto& jet : *jets){
    float area(-99.0),pt_unsub(-99.0),pt_sub(-99.0);
    jet->getAttribute("ActiveArea", area);
    if(m_doJetReclustering){
      pt_unsub = jet->pt();
      pt_sub = pt_unsub - event_rho*area;
    }
    else{
      pt_sub = jet->pt();
      pt_unsub = pt_unsub + event_rho*area;
    }

    if(pt_sub<5) continue; //only store jets with pT>5 !!!MeV!!!
    m_jnoarea0pt.push_back(pt_unsub*0.001);
    m_j0eta.push_back(jet->eta());
    m_j0phi.push_back(jet->phi());
    m_j0m.push_back(jet->m()*0.001);
    m_j0area.push_back(area);
    m_j0pt.push_back(pt_sub*0.001);
    //std::cout << "Jet" << std::endl;
    ///std::cout << jet->pt()*0.001 << ";" << jet->eta() << ";" << jet->phi() << std::endl;

    bool isPU;
    jet->getAttribute("isPU", isPU);
    m_j0isPU.push_back(isPU);
    int truthmatch = truth_match_i(*jet);
    if(truthmatch < 0){
      m_tj0pt.push_back(-99);
      m_tj0eta.push_back(-99);
      m_tj0phi.push_back(-99);
      m_tj0m.push_back(-99);
      m_tj0mindr.push_back(-99);
    }
    else{
      auto tjet = sorted_truth_jets.at(truthmatch);
      float mindr = minDR(*tjet);
      //if(mindr<0.6) continue; //only store jets matched to isolated truth jets

      m_tj0pt.push_back(tjet->pt()/1000.);
      m_tj0eta.push_back(tjet->eta());
      m_tj0phi.push_back(tjet->phi());
      m_tj0m.push_back(tjet->m()/1000.);
      m_tj0mindr.push_back(mindr);
    }

    if(m_doTracks){
      std::vector<const xAOD::IParticle*> jettracks;
      jet->getAssociatedObjects<xAOD::IParticle>(xAOD::JetAttribute::GhostTrack,jettracks);
      double ptsum_all = 0;
      double ptsum_pv = 0;
      double ptsum_pileup = 0;
      std::vector<float> trkpt;
      std::vector<float> trketa;
      std::vector<float> trkphi;
      std::vector<float> trkm;
      std::vector<bool> trkisHS;
      std::vector<bool> trkisPU;
      for (size_t i = 0; i < jettracks.size(); i++) {
        const xAOD::TrackParticle* trk = static_cast<const xAOD::TrackParticle*>(jettracks[i]);
        bool accept = (trk->pt()>500 && m_trackSelectionTool->accept(*trk,pv));
        if (accept) ptsum_all += trk->pt();
        if (accept && ((!trk->vertex() && fabs((trk->z0()+trk->vz()-pv->z())*sin(trk->theta()))<3.) || trk->vertex()==pv)) ptsum_pv += trk->pt();
        if (accept && trk->vertex() && trk->vertex()!=pv) ptsum_pileup += trk->pt();
        if (accept){
          trkpt.push_back(trk->pt());
          trketa.push_back(trk->eta());
          trkphi.push_back(trk->phi());
          trkm.push_back(trk->m());
          trkisHS.push_back(accept && ((!trk->vertex() && fabs((trk->z0()+trk->vz()-pv->z())*sin(trk->theta()))<3.) || trk->vertex()==pv));
          trkisPU.push_back(accept && trk->vertex() && trk->vertex()!=pv);
        }
      }
      double JVF = ptsum_all>0 ? ptsum_pv/ptsum_all : -1;
      double Rpt = ptsum_pv/jet->pt();
      double corrJVF = ptsum_pv+ptsum_pileup>0 ? ptsum_pv/(ptsum_pv+100*ptsum_pileup/n_putracks) : -1;
      //double JVT = corrJVF>=0 ? likelihood->Interpolate(corrJVF,std::min(Rpt,(float)1.)) : -0.1;
      double JVT = -2;
      //std::cout << Rpt << ";" << JVF << ";" << corrJVF << std::endl;

      //std::cout << "pv sum: " << Rpt*jet->pt() << std::endl;
      std::vector<float> sumPtTrkvec;
      jet->getAttribute( xAOD::JetAttribute::SumPtTrkPt500, sumPtTrkvec );
      std::sort(sumPtTrkvec.begin(),sumPtTrkvec.end()); //sorted, low to high
      //for(int si=0; si<sumPtTrkvec.size(); ++si) std::cout << "Vertex " << si << "; sum " << sumPtTrkvec[si] << std::endl;
      float maxRpt = sumPtTrkvec[sumPtTrkvec.size()-1];
      float medianRpt = sumPtTrkvec[sumPtTrkvec.size()/2];
      float deltaRpt = maxRpt-medianRpt;
      deltaRpt = deltaRpt/jet->pt();
      m_j0deltaRpt.push_back(deltaRpt);
      m_j0Rpt.push_back(Rpt);
      m_j0JVF.push_back(JVF);
      m_j0corrJVF.push_back(corrJVF);
      m_j0JVT.push_back(JVT);
      m_j0_trkpt.push_back(trkpt);
      m_j0_trketa.push_back(trketa);
      m_j0_trkphi.push_back(trkphi);
      m_j0_trkm.push_back(trkm);
      m_j0_trkisHS.push_back(trkisHS);
      m_j0_trkisPU.push_back(trkisPU);
    }

    if(m_doClusterInfo){
      std::vector<float> clpt;
      std::vector<float> cleta;
      std::vector<float> clphi;
      std::vector<float> clm;
      std::vector<float> cllatwidth;
      std::vector<float> cllongwidth;
      for(auto constit: jet->getConstituents()){
        clpt.push_back(constit->pt()*0.001);
        cleta.push_back(constit->eta());
        clphi.push_back(constit->phi());
        clm.push_back(constit->m());
        cllatwidth.push_back(constit->auxdata<float>("LATERAL"));
        cllongwidth.push_back(constit->auxdata<float>("LONGITUDINAL"));
        //std::cout << constit->pt()*0.001 << ";" << constit->eta() << ";" << constit->phi() << std::endl;
      }
      m_j0_clpt.push_back(clpt);
      m_j0_cleta.push_back(cleta);
      m_j0_clphi.push_back(clphi);
      m_j0_clm.push_back(clm);
      m_j0_cllatwidth.push_back(cllatwidth);
      m_j0_cllongwidth.push_back(cllongwidth);
    }

    CalcAnnulus(jet);
    m_j0sumpt5.push_back(sumpT5(*jet)/1000.);
    m_j0sumpt10.push_back(sumpT10(*jet)/1000.);
    m_j0sumpt15.push_back(sumpT15(*jet)/1000.);
    m_j0sumpt20.push_back(sumpT20(*jet)/1000.);
    m_j0sumpt25.push_back(sumpT25(*jet)/1000.);
    m_j0sumpt30.push_back(sumpT30(*jet)/1000.);
    m_j0sumpt35.push_back(sumpT35(*jet)/1000.);
    m_j0sumpt40.push_back(sumpT40(*jet)/1000.);

  }

  m_tree->Fill();
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode JetLearningAnalysis :: postExecute ()
{
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode JetLearningAnalysis :: finalize ()
{
  xAOD::TEvent* event = wk()->xaodEvent();
  if(m_trackSelectionTool) delete m_trackSelectionTool;
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode JetLearningAnalysis :: histFinalize ()
{
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode JetLearningAnalysis :: SetRho(DataVector<xAOD::Jet_v1> jets,float event_rho){
  static SG::AuxElement::Decorator< float > rho("rho");
  for(int iJ=0; iJ<jets.size(); iJ++) {
    DataVector<xAOD::Jet_v1>::ElementProxy jet = jets.at(iJ);
    rho(*jet) = event_rho;
  }
}

EL::StatusCode JetLearningAnalysis :: FindRho(const xAOD::CaloClusterContainer* in_clusters,float& rho, float& sigma){
  //CaloClusterChangeSignalStateList stateHelperList;
  std::vector<fastjet::PseudoJet> inputConst;
  for(const auto clust: *in_clusters){
    //if(m_doLC) stateHelperList.add(clust,xAOD::CaloCluster::State(1)); //default is calibrated but we can make it explicit anyway
    //else stateHelperList.add(clust,xAOD::CaloCluster::State(0));
    fastjet::PseudoJet test;
    test = fastjet::PseudoJet(clust->p4());
    if(clust->e() >= 0) inputConst.push_back(test);
  }

  fastjet::Selector jselector = fastjet::SelectorAbsRapRange(0.0,2.1);
  fastjet::JetAlgorithm algo = fastjet::kt_algorithm;
  float jetR = 0.4;
  fastjet::JetDefinition jetDef(algo, jetR,fastjet::E_scheme, fastjet::Best);
  fastjet::AreaDefinition area_def(fastjet::voronoi_area, fastjet::VoronoiAreaSpec(0.9));

  fastjet::JetMedianBackgroundEstimator bge(jselector,jetDef,area_def);
  bge.set_particles(inputConst);

  rho = bge.rho();
  sigma = bge.sigma();
}

float deltaR(DataVector<xAOD::Jet_v1>::ElementProxy jet1, DataVector<xAOD::Jet_v1>::ElementProxy jet2){
  fastjet::PseudoJet pjet1 = fastjet::PseudoJet(jet1->p4());
  fastjet::PseudoJet pjet2 = fastjet::PseudoJet(jet2->p4());
  float deltaphi = pjet1.delta_phi_to(pjet2);
  float deltaeta = jet1->eta() - jet2->eta();
  float dR =  sqrt(pow(deltaphi,2) + pow(deltaeta,2));
  return dR;
}

EL::StatusCode JetLearningAnalysis :: FindTruthMatch(DataVector<xAOD::Jet_v1> jets, DataVector<xAOD::Jet_v1> tjets){
  std::vector<int> matched;
  float MAXJETTRUTHMATCHDR = 0.3;
  float MINPUJETTRUTHMATCHDR = 0.6;
  static SG::AuxElement::Decorator< int > truth_match_i("truth_match_i");
  static SG::AuxElement::Decorator< bool > isPU("isPU");
  for(int iJ=0; iJ<jets.size(); iJ++) {
    DataVector<xAOD::Jet_v1>::ElementProxy jet = jets.at(iJ);
    float mindR= 999.99;
    float maxPt=-999.99;
    //int minDRindex =-1;
    int maxPtIndex =-1;
    for(int iTrueJ=0; iTrueJ<tjets.size(); iTrueJ++){
      DataVector<xAOD::Jet_v1>::ElementProxy tjet = tjets.at(iTrueJ);

      float dR = deltaR(jet,tjet);
      if(dR<mindR){ mindR = dR;}
      if(std::find(matched.begin(), matched.end(), iTrueJ) != matched.end())
        continue; //if true jet has already been matched, skip it -> truth jet is matched to at most one reco jet 
      if(dR < MAXJETTRUTHMATCHDR && maxPt < tjet->pt())
        { maxPt = tjet->pt(); maxPtIndex = iTrueJ;} //match to highest pT truth jet within MAXDR, not closest
    }//true jets

    if(maxPtIndex != -1){
      matched.push_back(maxPtIndex);
      /*if (!(jet(maxPtIndex, TruthJetType).Exists(JetType+"_match")))
        jet(maxPtIndex, TruthJetType).Add(JetType+"_match", thejet, true);*/ //add link from truth jet to jet as well eventually
    }
    truth_match_i(*jet) = maxPtIndex; //if no match truth_match_i == -1
    isPU(*jet) = mindR>MINPUJETTRUTHMATCHDR; //mindR to any truth jet
  }//jet loop
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode JetLearningAnalysis :: SetMinDR(DataVector<xAOD::Jet_v1> jets,float threshold){
  static SG::AuxElement::Decorator< float > minDR("minDR");
  for(int iJ=0; iJ<jets.size(); iJ++) {
    DataVector<xAOD::Jet_v1>::ElementProxy jet1 = jets.at(iJ);
    float mindR= 999.99;
    for(int jJ=0; jJ<jets.size(); jJ++){
      if(iJ==jJ) continue;
      DataVector<xAOD::Jet_v1>::ElementProxy jet2 = jets.at(jJ);
      if(jet2->pt() < threshold) continue;

      float dR = deltaR(jet1,jet2);

      if(dR<mindR){ mindR = dR;}
    } //inner jet loop

    //if(mindR < 999.99){
    minDR(*jet1) = mindR; 
    //}
    //else{
    //  return EL::StatusCode::FAILURE;
    //}
  }//jet loop
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode JetLearningAnalysis :: CalcAnnulus(const xAOD::Jet_v1* jet){
  static SG::AuxElement::Decorator< float > sumpT5("sumpT5");
  static SG::AuxElement::Decorator< float > sumpT10("sumpT10");
  static SG::AuxElement::Decorator< float > sumpT15("sumpT15");
  static SG::AuxElement::Decorator< float > sumpT20("sumpT20");
  static SG::AuxElement::Decorator< float > sumpT25("sumpT25");
  static SG::AuxElement::Decorator< float > sumpT30("sumpT30");
  static SG::AuxElement::Decorator< float > sumpT35("sumpT35");
  static SG::AuxElement::Decorator< float > sumpT40("sumpT40");

  float sumpTs[8] = {0,0,0,0,0,0,0,0};
  for(auto constit: jet->getConstituents()){
    float dR = xAOD::P4Helpers::deltaR(constit->rawConstituent(),jet);
    for(int i=0; i<8; i++){
      if(dR<float(i+1)*0.05){
        sumpTs[i]+=constit->pt();
      }
    }
  }
  sumpT5(*jet) = sumpTs[0];
  sumpT10(*jet) = sumpTs[1];
  sumpT15(*jet) = sumpTs[2];
  sumpT20(*jet) = sumpTs[3];
  sumpT25(*jet) = sumpTs[4];
  sumpT30(*jet) = sumpTs[5];
  sumpT35(*jet) = sumpTs[6];
  sumpT40(*jet) = sumpTs[7];
}
