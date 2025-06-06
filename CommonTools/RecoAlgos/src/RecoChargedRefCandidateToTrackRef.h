#ifndef CommonTools_RecoAlgos_TrackToRefCandidate_h
#define CommonTools_RecoAlgos_TrackToRefCandidate_h
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ParameterSet/interface/ParameterSetfwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedRefCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedRefCandidateFwd.h"

namespace edm {
  class EventSetup;
  class ParameterSet;
  class ConsumesCollector;
}  // namespace edm

namespace converter {
  struct RecoChargedRefCandidateToTrackRef {
    typedef reco::RecoChargedRefCandidate value_type;
    typedef reco::RecoChargedRefCandidateCollection Components;
    typedef reco::TrackRef Candidate;
    RecoChargedRefCandidateToTrackRef(const edm::ParameterSet& cfg, const edm::ConsumesCollector&) {}
    void beginFirstRun(const edm::EventSetup&) {}
    void convert(const reco::RecoChargedRefCandidateRef& c, reco::TrackRef& trkRef) const { trkRef = c->track(); }
    static void fillPSetDescription(edm::ParameterSetDescription& desc) {}
  };
}  // namespace converter

#endif
