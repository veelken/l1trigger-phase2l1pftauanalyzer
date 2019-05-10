#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/L1TrackAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"                              // JetMCTagUtils::genTauDecayMode()

#include "L1Trigger/TallinnL1PFTauAnalyzer/interface/GenChargedHadronToTrackMatch.h" // GenChargedHadronToOfflineTrackMatch, GenChargedHadronToL1TrackMatch

#include <algorithm> // std::sort()

enum { kQualityCuts_disabled, kQualityCuts_enabled };

L1TrackAnalyzer::L1TrackAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_genTaus_ = cfg.getParameter<edm::InputTag>("srcGenTaus");
  token_genTaus_ = consumes<reco::GenJetCollection>(src_genTaus_);

  src_offlineVertices_ = cfg.getParameter<edm::InputTag>("srcOfflineVertices");
  token_offlineVertices_ = consumes<reco::VertexCollection>(src_offlineVertices_);

  src_offlineTracks_ = cfg.getParameter<edm::InputTag>("srcOfflineTracks");
  token_offlineTracks_ = consumes<reco::TrackCollection>(src_offlineTracks_);

  src_offlinePFCands_ = cfg.getParameter<edm::InputTag>("srcOfflinePFCands");
  token_offlinePFCands_ = consumes<reco::PFCandidateCollection>(src_offlinePFCands_);

  src_l1Vertices_ = cfg.getParameter<edm::InputTag>("srcL1Vertices");
  token_l1Vertices_ = consumes<l1t::VertexCollection>(src_l1Vertices_);

  src_l1Tracks_ = cfg.getParameter<edm::InputTag>("srcL1Tracks");
  token_l1Tracks_ = consumes<l1t::TrackCollection>(src_l1Tracks_);

  src_l1PFVertex_z_ = cfg.getParameter<edm::InputTag>("srcL1PFVertex_z");
  token_l1PFVertex_z_ = consumes<float>(src_l1PFVertex_z_);

  src_l1PFCands_ = cfg.getParameter<edm::InputTag>("srcL1PFCands");
  token_l1PFCands_ = consumes<l1t::PFCandidateCollection>(src_l1PFCands_);

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

L1TrackAnalyzer::~L1TrackAnalyzer()
{
  for ( auto efficiencyPlot : efficiencyPlots_offlineTracks_woQualityCuts_ ) 
  {
    delete efficiencyPlot;
  }

  for ( auto efficiencyPlot : efficiencyPlots_offlineTracks_wQualityCuts_ ) 
  {
    delete efficiencyPlot;
  }

  for ( auto efficiencyPlot : efficiencyPlots_offlinePFCandTracks_woQualityCuts_ ) 
  {
    delete efficiencyPlot;
  }

  for ( auto efficiencyPlot : efficiencyPlots_offlinePFCandTracks_wQualityCuts_ ) 
  {
    delete efficiencyPlot;
  }
  
  for ( auto efficiencyPlot : efficiencyPlots_l1Tracks_woQualityCuts_ ) 
  {
    delete efficiencyPlot;
  }

  for ( auto efficiencyPlot : efficiencyPlots_l1Tracks_wQualityCuts_ ) 
  {
    delete efficiencyPlot;
  }

  for ( auto efficiencyPlot : efficiencyPlots_l1PFCandTracks_woQualityCuts_ ) 
  {
    delete efficiencyPlot;
  }

  for ( auto efficiencyPlot : efficiencyPlots_l1PFCandTracks_wQualityCuts_ ) 
  {
    delete efficiencyPlot;
  }
}

void L1TrackAnalyzer::beginJob()
{
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    throw cms::Exception("L1TrackAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());

  std::vector<std::string> genTau_decayModes = { "oneProng0Pi0", "oneProng1Pi0", "oneProng2Pi0", "threeProng0Pi0", "threeProng1Pi0", "all" };
  std::vector<double> genChargedHadron_absEtaRanges = { 1.0, 1.4 };
  for ( auto genTau_decayMode : genTau_decayModes )
  {
    for ( auto genChargedHadron_absEtaRange : genChargedHadron_absEtaRanges )
    {
      efficiencyPlotEntryType* efficiencyPlot_offlineTracks_woQualityCuts = new efficiencyPlotEntryType(
        "offlineTrack_woQualityCuts", 1., genChargedHadron_absEtaRange, genTau_decayMode); 
      efficiencyPlot_offlineTracks_woQualityCuts->bookHistograms(dqmStore);
      efficiencyPlots_offlineTracks_woQualityCuts_.push_back(efficiencyPlot_offlineTracks_woQualityCuts);
            
      efficiencyPlotEntryType* efficiencyPlot_offlineTracks_wQualityCuts = new efficiencyPlotEntryType(
        "offlineTrack_wQualityCuts", 1., genChargedHadron_absEtaRange, genTau_decayMode); 
      efficiencyPlot_offlineTracks_wQualityCuts->bookHistograms(dqmStore);
      efficiencyPlots_offlineTracks_wQualityCuts_.push_back(efficiencyPlot_offlineTracks_wQualityCuts);
      
      efficiencyPlotEntryType* efficiencyPlot_offlinePFCandTracks_woQualityCuts = new efficiencyPlotEntryType(
        "offlinePFCandTrack_woQualityCuts", 1., genChargedHadron_absEtaRange, genTau_decayMode); 
      efficiencyPlot_offlinePFCandTracks_woQualityCuts->bookHistograms(dqmStore);
      efficiencyPlots_offlinePFCandTracks_woQualityCuts_.push_back(efficiencyPlot_offlinePFCandTracks_woQualityCuts);
            
      efficiencyPlotEntryType* efficiencyPlot_offlinePFCandTracks_wQualityCuts = new efficiencyPlotEntryType(
        "offlinePFCandTrack_wQualityCuts", 1., genChargedHadron_absEtaRange, genTau_decayMode); 
      efficiencyPlot_offlinePFCandTracks_wQualityCuts->bookHistograms(dqmStore);
      efficiencyPlots_offlinePFCandTracks_wQualityCuts_.push_back(efficiencyPlot_offlinePFCandTracks_wQualityCuts);
      
      efficiencyPlotEntryType* efficiencyPlot_l1Tracks_woQualityCuts = new efficiencyPlotEntryType(
        "l1Track_woQualityCuts", 1., genChargedHadron_absEtaRange, genTau_decayMode); 
      efficiencyPlot_l1Tracks_woQualityCuts->bookHistograms(dqmStore);
      efficiencyPlots_l1Tracks_woQualityCuts_.push_back(efficiencyPlot_l1Tracks_woQualityCuts);
            
      efficiencyPlotEntryType* efficiencyPlot_l1Tracks_wQualityCuts = new efficiencyPlotEntryType(
        "l1Track_wQualityCuts", 1., genChargedHadron_absEtaRange, genTau_decayMode); 
      efficiencyPlot_l1Tracks_wQualityCuts->bookHistograms(dqmStore);
      efficiencyPlots_l1Tracks_wQualityCuts_.push_back(efficiencyPlot_l1Tracks_wQualityCuts);
      
      efficiencyPlotEntryType* efficiencyPlot_l1PFCandTracks_woQualityCuts = new efficiencyPlotEntryType(
        "l1PFCandTrack_woQualityCuts", 1., genChargedHadron_absEtaRange, genTau_decayMode); 
      efficiencyPlot_l1PFCandTracks_woQualityCuts->bookHistograms(dqmStore);
      efficiencyPlots_l1PFCandTracks_woQualityCuts_.push_back(efficiencyPlot_l1PFCandTracks_woQualityCuts);
            
      efficiencyPlotEntryType* efficiencyPlot_l1PFCandTracks_wQualityCuts = new efficiencyPlotEntryType(
        "l1PFCandTrack_wQualityCuts", 1., genChargedHadron_absEtaRange, genTau_decayMode); 
      efficiencyPlot_l1PFCandTracks_wQualityCuts->bookHistograms(dqmStore);
      efficiencyPlots_l1PFCandTracks_wQualityCuts_.push_back(efficiencyPlot_l1PFCandTracks_wQualityCuts);
    }
  }
}

namespace
{
  // auxiliary function to select offline reconstructed tracks of "good quality" 
  bool passesOfflineTrackQualityCuts(const reco::Track& track, const reco::Vertex* primaryVertex)
  {    
    bool passesQualityCuts = ( primaryVertex                                             && 
			       TMath::Abs(track.dz(primaryVertex->position()))  <   0.4  &&
			       track.hitPattern().numberOfValidPixelHits()      >=  0    &&
			       track.hitPattern().numberOfValidHits()           >=  3    &&
			       TMath::Abs(track.dxy(primaryVertex->position())) <   0.03 &&
			       track.normalizedChi2()                           < 100.   ) ? true : false;
    return passesQualityCuts;
  }

  // auxiliary function for matches between generator-level charged hadrons produced in tau decays and offline reconstructed tracks
  std::vector<GenChargedHadronToOfflineTrackMatch_and_genTau_decayMode> 
  getGenChargedHadronToOfflineTrackMatches(const std::vector<GenChargedHadron_and_genTau_decayMode>& genTauChargedHadrons, std::vector<const reco::Track*> recTracks)
  {
    std::vector<GenChargedHadronToOfflineTrackMatch_and_genTau_decayMode> genChargedHadronToTrackMatches;
    const double dRmatch = 0.05;
    for ( auto genTauChargedHadron : genTauChargedHadrons )
    {
      for ( auto recTrack : recTracks )
      {
	double dR = deltaR(genTauChargedHadron.genChargedHadron_eta(), genTauChargedHadron.genChargedHadron_phi(), recTrack->eta(), recTrack->phi());
	if ( dR < dRmatch ) 
        {
	  genChargedHadronToTrackMatches.push_back(GenChargedHadronToOfflineTrackMatch_and_genTau_decayMode(
	    GenChargedHadronToOfflineTrackMatch(genTauChargedHadron.genChargedHadron(), recTrack),
	    genTauChargedHadron.genTau_decayMode(),
	    dR));
	}
      }
      // CV: add "empty" match to allow for genParticles without recTrack match (tracking inefficiency)
      genChargedHadronToTrackMatches.push_back(GenChargedHadronToOfflineTrackMatch_and_genTau_decayMode(
	GenChargedHadronToOfflineTrackMatch(genTauChargedHadron.genChargedHadron(), nullptr), 
	genTauChargedHadron.genTau_decayMode(),
	1.e+3));
    }
    return genChargedHadronToTrackMatches;
  }

  // auxiliary function to select "good quality" tracks reconstructed on L1 trigger level
  bool passesL1TrackQualityCuts(const l1t::Track& track, double primaryVertex_z)
  {    
    //bool passesQualityCuts = ( TMath::Abs(track.vertex().z() - primaryVertex_z) < 0.2 ) ? true : false;
    bool passesQualityCuts = true; // CV: dz cut not implemented for l1t::Track yet
    return passesQualityCuts;
  }

  bool passesL1TrackQualityCuts(const l1t::PFTrack& track, double primaryVertex_z)
  {    
    bool passesQualityCuts = ( TMath::Abs(track.vertex().z() - primaryVertex_z) < 0.2 ) ? true : false;
    return passesQualityCuts;
  }

  // auxiliary function for matches between generator-level charged hadrons produced in tau decays and tracks reconstructed on L1 trigger level
  std::vector<GenChargedHadronToL1TrackMatch_and_genTau_decayMode> 
  getGenChargedHadronToL1TrackMatches(const std::vector<GenChargedHadron_and_genTau_decayMode>& genTauChargedHadrons, std::vector<const l1t::Track*> recTracks)
  {
    std::vector<GenChargedHadronToL1TrackMatch_and_genTau_decayMode> genChargedHadronToTrackMatches;
    const double dRmatch = 0.05;
    for ( auto genTauChargedHadron : genTauChargedHadrons )
    {
      for ( auto recTrack : recTracks )
      {
	const unsigned nParam = 4;
	double recTrack_eta = recTrack->getMomentum(nParam).eta();
	double recTrack_phi = recTrack->getMomentum(nParam).phi();
	double dR = deltaR(genTauChargedHadron.genChargedHadron_eta(), genTauChargedHadron.genChargedHadron_phi(), recTrack_eta, recTrack_phi);
	if ( dR < dRmatch ) 
        {
	  genChargedHadronToTrackMatches.push_back(GenChargedHadronToL1TrackMatch_and_genTau_decayMode(
	    GenChargedHadronToL1TrackMatch(genTauChargedHadron.genChargedHadron(), recTrack),
	    genTauChargedHadron.genTau_decayMode(),
	    dR));
	}
      }
      // CV: add "empty" match to allow for genParticles without recTrack match (tracking inefficiency)
      genChargedHadronToTrackMatches.push_back(GenChargedHadronToL1TrackMatch_and_genTau_decayMode(
	GenChargedHadronToL1TrackMatch(genTauChargedHadron.genChargedHadron(), nullptr),
	genTauChargedHadron.genTau_decayMode(),
	1.e+3));
    }
    return genChargedHadronToTrackMatches;
  }

  // auxiliary function to sort (ranke) matches (needed for cleaning)
  template<class T>
  bool isLowerDeltaR(const T* genChargedHadronToTrackMatch1, const T* genChargedHadronToTrackMatch2)
  {
    return (genChargedHadronToTrackMatch1->dR() < genChargedHadronToTrackMatch2->dR());
  }

  // auxiliary function to clean matches (in order to guarantee that each generator-level charged hadron and each reconstructed track is used in the matching only once)
  template <class T>
  std::vector<const T*> 
  cleanGenChargedHadronToTrackMatches(const std::vector<T>& genChargedHadronToTrackMatches)
  {
    std::list<const T*> matches_remaining;
    for ( auto match : genChargedHadronToTrackMatches )
    {
      matches_remaining.push_back(&match);
    }

    std::vector<const T*> matches_selected;
    while ( matches_remaining.size() > 0 ) 
    {
      matches_remaining.sort(isLowerDeltaR<T>);
      const T* bestMatch = matches_remaining.front();
      assert(bestMatch);
      matches_selected.push_back(bestMatch);
      for ( auto match : matches_remaining ) 
      {
        if ( isOverlap(match->genChargedHadronToTrackMatch(), bestMatch->genChargedHadronToTrackMatch()) )
	{
          matches_remaining.remove(match);
        }
      }      
    }
    return matches_selected;
  }
}
    
void L1TrackAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<reco::GenJetCollection> genTaus;
  evt.getByToken(token_genTaus_, genTaus);

  std::vector<GenChargedHadron_and_genTau_decayMode> genTauChargedHadrons;
  for ( reco::GenJetCollection::const_iterator genTau = genTaus->begin(); genTau != genTaus->end(); ++genTau ) 
  {
    std::string genTau_decayMode = JetMCTagUtils::genTauDecayMode(*genTau);
    std::vector<const reco::GenParticle*> genTau_daughters = genTau->getGenConstituents();
    for ( auto genTau_daughter : genTau_daughters ) 
    {
      if ( genTau_daughter->charge() != 0 ) 
      {
	genTauChargedHadrons.push_back(GenChargedHadron_and_genTau_decayMode(genTau_daughter, genTau_decayMode));
      }
    }
  }

  const double evtWeight = 1.;

  //-----------------------------------------------------------------------------
  // process offline reconstructed vertices
  const reco::Vertex* offlinePrimaryVertex = nullptr;
  if ( src_offlineVertices_.label() != "" ) 
  {
    edm::Handle<reco::VertexCollection> offlineVertices;
    evt.getByToken(token_offlineVertices_, offlineVertices);    
    if ( offlineVertices->size() > 0 ) offlinePrimaryVertex = &offlineVertices->at(0);
  }
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // process offline reconstructed tracks
  if ( src_offlineTracks_.label() != "" )
  {
    edm::Handle<reco::TrackCollection> offlineTracks;
    evt.getByToken(token_offlineTracks_, offlineTracks);
  
    for ( int idxQualityCuts = kQualityCuts_disabled; idxQualityCuts <= kQualityCuts_enabled; ++idxQualityCuts )
    {
      std::vector<const reco::Track*> selectedOfflineTracks;
      for ( reco::TrackCollection::const_iterator recTrack = offlineTracks->begin(); recTrack != offlineTracks->end(); ++recTrack )
      {
        if (  idxQualityCuts == kQualityCuts_disabled                                                                    ||
	     (idxQualityCuts == kQualityCuts_enabled  && passesOfflineTrackQualityCuts(*recTrack, offlinePrimaryVertex)) )
	{
	  selectedOfflineTracks.push_back(&(*recTrack));
	}
      }
      
      std::vector<GenChargedHadronToOfflineTrackMatch_and_genTau_decayMode> genChargedHadronToTrackMatches = getGenChargedHadronToOfflineTrackMatches(
	genTauChargedHadrons,
	selectedOfflineTracks);     
      std::vector<const GenChargedHadronToOfflineTrackMatch_and_genTau_decayMode*> cleanedGenChargedHadronToTrackMatches = cleanGenChargedHadronToTrackMatches(
        genChargedHadronToTrackMatches);     

      const std::vector<efficiencyPlotEntryType*>* efficiencyPlots = nullptr;
      if      ( idxQualityCuts == kQualityCuts_disabled )  efficiencyPlots = &efficiencyPlots_offlineTracks_woQualityCuts_;
      else if ( idxQualityCuts == kQualityCuts_enabled  )  efficiencyPlots = &efficiencyPlots_offlineTracks_wQualityCuts_;
      else assert(0);
      for ( auto efficiencyPlot : *efficiencyPlots )
      {
	efficiencyPlot->fillHistograms(cleanedGenChargedHadronToTrackMatches, evtWeight);
      }
    }
  }
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // process tracks associated to offline reconstructed PF candidates
  if ( src_offlinePFCands_.label() != "" )
  {
    edm::Handle<reco::PFCandidateCollection> offlinePFCands;
    evt.getByToken(token_offlinePFCands_, offlinePFCands);
  
    for ( int idxQualityCuts = kQualityCuts_disabled; idxQualityCuts <= kQualityCuts_enabled; ++idxQualityCuts )
    {
      std::vector<const reco::Track*> selectedOfflineTracks;
      for ( reco::PFCandidateCollection::const_iterator recPFCand = offlinePFCands->begin(); recPFCand != offlinePFCands->end(); ++recPFCand )
      {
	if ( recPFCand->charge() != 0 )
	{
	  const reco::Track* recTrack = nullptr;
	  if      ( recPFCand->trackRef().isNonnull()    ) recTrack = recPFCand->trackRef().get();
	  else if ( recPFCand->gsfTrackRef().isNonnull() ) recTrack = recPFCand->gsfTrackRef().get();
	  if ( !recTrack ) continue;

          if (  idxQualityCuts == kQualityCuts_disabled                                                                     ||
	       (idxQualityCuts == kQualityCuts_enabled  && passesOfflineTrackQualityCuts(*recTrack, offlinePrimaryVertex)) )
  	  {
	    selectedOfflineTracks.push_back(&(*recTrack));
	  }
	}
      }
      
      std::vector<GenChargedHadronToOfflineTrackMatch_and_genTau_decayMode> genChargedHadronToTrackMatches = getGenChargedHadronToOfflineTrackMatches(
	genTauChargedHadrons,
	selectedOfflineTracks);     
      std::vector<const GenChargedHadronToOfflineTrackMatch_and_genTau_decayMode*> cleanedGenChargedHadronToTrackMatches = cleanGenChargedHadronToTrackMatches(
        genChargedHadronToTrackMatches);     

      const std::vector<efficiencyPlotEntryType*>* efficiencyPlots = nullptr;
      if      ( idxQualityCuts == kQualityCuts_disabled )  efficiencyPlots = &efficiencyPlots_offlinePFCandTracks_woQualityCuts_;
      else if ( idxQualityCuts == kQualityCuts_enabled  )  efficiencyPlots = &efficiencyPlots_offlinePFCandTracks_wQualityCuts_;
      else assert(0);
      for ( auto efficiencyPlot : *efficiencyPlots )
      {
	efficiencyPlot->fillHistograms(cleanedGenChargedHadronToTrackMatches, evtWeight);
      }
    }
  }
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // process tracks reconstructed on L1 trigger level
  if ( src_l1Vertices_.label() != "" && src_l1Tracks_.label() != "" )
  {
    edm::Handle<l1t::VertexCollection> l1Vertices;
    evt.getByToken(token_l1Vertices_, l1Vertices);
    if ( l1Vertices->size() > 0 ) 
    {
      double l1PrimaryVertex_z = l1Vertices->at(0).z0();

      edm::Handle<l1t::TrackCollection> l1Tracks;
      evt.getByToken(token_l1Tracks_, l1Tracks);
  
      for ( int idxQualityCuts = kQualityCuts_disabled; idxQualityCuts <= kQualityCuts_enabled; ++idxQualityCuts )
      {
        std::vector<const l1t::Track*> selectedL1Tracks;
        for ( l1t::TrackCollection::const_iterator recTrack = l1Tracks->begin(); recTrack != l1Tracks->end(); ++recTrack )
        {
          if (  idxQualityCuts == kQualityCuts_disabled                                                            ||
	       (idxQualityCuts == kQualityCuts_enabled  && passesL1TrackQualityCuts(*recTrack, l1PrimaryVertex_z)) )
  	  {
	    selectedL1Tracks.push_back(&(*recTrack));
	  }
	}
      
        std::vector<GenChargedHadronToL1TrackMatch_and_genTau_decayMode> genChargedHadronToTrackMatches = getGenChargedHadronToL1TrackMatches(
  	  genTauChargedHadrons,
          selectedL1Tracks);     
        std::vector<const GenChargedHadronToL1TrackMatch_and_genTau_decayMode*> cleanedGenChargedHadronToTrackMatches = cleanGenChargedHadronToTrackMatches(
          genChargedHadronToTrackMatches);     

        const std::vector<efficiencyPlotEntryType*>* efficiencyPlots = nullptr;
        if      ( idxQualityCuts == kQualityCuts_disabled )  efficiencyPlots = &efficiencyPlots_l1Tracks_woQualityCuts_;
        else if ( idxQualityCuts == kQualityCuts_enabled  )  efficiencyPlots = &efficiencyPlots_l1Tracks_wQualityCuts_;
        else assert(0);
        for ( auto efficiencyPlot : *efficiencyPlots )
        {
	  efficiencyPlot->fillHistograms(cleanedGenChargedHadronToTrackMatches, evtWeight);
        }
      }
    }
  }
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // process tracks associated to PF candidates reconstructed on L1 trigger level
  if ( src_l1PFVertex_z_.label() != "" && src_l1PFCands_.label() != "" )
  {
    edm::Handle<float> l1PrimaryVertex_z;
    evt.getByToken(token_l1PFVertex_z_, l1PrimaryVertex_z);

    edm::Handle<l1t::PFCandidateCollection> l1PFCands;
    evt.getByToken(token_l1PFCands_, l1PFCands);
  
    for ( int idxQualityCuts = kQualityCuts_disabled; idxQualityCuts <= kQualityCuts_enabled; ++idxQualityCuts )
    {
      std::vector<const l1t::Track*> selectedL1Tracks;
      for ( l1t::PFCandidateCollection::const_iterator recPFCand = l1PFCands->begin(); recPFCand != l1PFCands->end(); ++recPFCand )
      {
	if ( recPFCand->charge() != 0 )
	{
	  const l1t::PFTrack* recPFTrack = nullptr;
	  if ( recPFCand->pfTrack().isNonnull() ) recPFTrack = recPFCand->pfTrack().get();
	  if ( !recPFTrack ) continue;

	  const l1t::Track* recTrack = nullptr;
	  if ( recPFTrack->track().isNonnull() ) recTrack = recPFTrack->track().get();
	  if ( !recTrack ) continue;

          if (  idxQualityCuts == kQualityCuts_disabled                                                               ||
	       (idxQualityCuts == kQualityCuts_enabled  && passesL1TrackQualityCuts(*recPFTrack, *l1PrimaryVertex_z)) ) // CV: apply dz cut using l1t::PFTrack
	  {
	    selectedL1Tracks.push_back(&(*recTrack));
	  }
	}
      }
      
      std::vector<GenChargedHadronToL1TrackMatch_and_genTau_decayMode> genChargedHadronToTrackMatches = getGenChargedHadronToL1TrackMatches(
	genTauChargedHadrons,
	selectedL1Tracks);     
      std::vector<const GenChargedHadronToL1TrackMatch_and_genTau_decayMode*> cleanedGenChargedHadronToTrackMatches = cleanGenChargedHadronToTrackMatches(
        genChargedHadronToTrackMatches);     

      const std::vector<efficiencyPlotEntryType*>* efficiencyPlots = nullptr;
      if      ( idxQualityCuts == kQualityCuts_disabled )  efficiencyPlots = &efficiencyPlots_l1PFCandTracks_woQualityCuts_;
      else if ( idxQualityCuts == kQualityCuts_enabled  )  efficiencyPlots = &efficiencyPlots_l1PFCandTracks_wQualityCuts_;
      else assert(0);
      for ( auto efficiencyPlot : *efficiencyPlots )
      {
	efficiencyPlot->fillHistograms(cleanedGenChargedHadronToTrackMatches, evtWeight);
      }
    }
  }
  //-----------------------------------------------------------------------------
}

void L1TrackAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(L1TrackAnalyzer);


