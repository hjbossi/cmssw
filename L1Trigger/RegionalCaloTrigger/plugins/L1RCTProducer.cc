#include "L1Trigger/RegionalCaloTrigger/interface/L1RCTProducer.h"

// RunInfo stuff
#include "CondFormats/RunInfo/interface/RunInfo.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <vector>
using std::vector;
#include <iostream>

using std::cout;
using std::endl;
namespace {
  constexpr int crateFED[18][6] = {{613, 614, 603, 702, 718, 1118},
                                   {611, 612, 602, 700, 718, 1118},
                                   {627, 610, 601, 716, 722, 1122},
                                   {625, 626, 609, 714, 722, 1122},
                                   {623, 624, 608, 712, 722, 1122},
                                   {621, 622, 607, 710, 720, 1120},
                                   {619, 620, 606, 708, 720, 1120},
                                   {617, 618, 605, 706, 720, 1120},
                                   {615, 616, 604, 704, 718, 1118},
                                   {631, 632, 648, 703, 719, 1118},
                                   {629, 630, 647, 701, 719, 1118},
                                   {645, 628, 646, 717, 723, 1122},
                                   {643, 644, 654, 715, 723, 1122},
                                   {641, 642, 653, 713, 723, 1122},
                                   {639, 640, 652, 711, 721, 1120},
                                   {637, 638, 651, 709, 721, 1120},
                                   {635, 636, 650, 707, 721, 1120},
                                   {633, 634, 649, 705, 719, 1118}};
}
L1RCTProducer::L1RCTProducer(const edm::ParameterSet &conf)
    : rctLookupTables(new L1RCTLookupTables),
      rct(new L1RCT(rctLookupTables.get())),
      useEcal(conf.getParameter<bool>("useEcal")),
      useHcal(conf.getParameter<bool>("useHcal")),
      ecalDigis(conf.getParameter<std::vector<edm::InputTag>>("ecalDigis")),
      hcalDigis(conf.getParameter<std::vector<edm::InputTag>>("hcalDigis")),
      bunchCrossings(conf.getParameter<std::vector<int>>("BunchCrossings")),
      getFedsFromOmds(conf.getParameter<bool>("getFedsFromOmds")),
      queryDelayInLS(conf.getParameter<unsigned int>("queryDelayInLS")),
      queryIntervalInLS(conf.getParameter<unsigned int>("queryIntervalInLS")),
      conditionsLabel(conf.getParameter<std::string>("conditionsLabel")),
      fedUpdatedMask(nullptr),

      rctParamsToken_(esConsumes<edm::Transition::BeginRun>(edm::ESInputTag("", conditionsLabel))),
      emScaleToken_(esConsumes<edm::Transition::BeginRun>(edm::ESInputTag("", conditionsLabel))),
      ecalScaleToken_(esConsumes<edm::Transition::BeginRun>(edm::ESInputTag("", conditionsLabel))),
      hcalScaleToken_(esConsumes<edm::Transition::BeginRun>(edm::ESInputTag("", conditionsLabel))),

      beginRunRunInfoToken_(esConsumes<edm::Transition::BeginRun>()),

      beginRunChannelMaskToken_(esConsumes<edm::Transition::BeginRun>()),
      beginRunHotChannelMaskToken_(esConsumes<edm::Transition::BeginRun>()) {
  produces<L1CaloEmCollection>();
  produces<L1CaloRegionCollection>();

  if (getFedsFromOmds) {
    beginLumiChannelMaskToken_ = esConsumes<edm::Transition::BeginLuminosityBlock>();
    beginLumiHotChannelMaskToken_ = esConsumes<edm::Transition::BeginLuminosityBlock>();
    beginLumiRunInfoToken_ = esConsumes<edm::Transition::BeginLuminosityBlock>();
    omdsRunInfoToken_ = esConsumes<edm::Transition::BeginLuminosityBlock>(edm::ESInputTag("", "OmdsFedVector"));
  }

  for (unsigned int ihc = 0; ihc < hcalDigis.size(); ihc++) {
    consumes<edm::SortedCollection<HcalTriggerPrimitiveDigi, edm::StrictWeakOrdering<HcalTriggerPrimitiveDigi>>>(
        hcalDigis[ihc]);
  }

  for (unsigned int iec = 0; iec < ecalDigis.size(); iec++) {
    consumes<edm::SortedCollection<EcalTriggerPrimitiveDigi, edm::StrictWeakOrdering<EcalTriggerPrimitiveDigi>>>(
        ecalDigis[iec]);
  }
}

void L1RCTProducer::beginRun(edm::Run const &run, const edm::EventSetup &eventSetup) {
  //  std::cout << "getFedsFromOmds is " << getFedsFromOmds << std::endl;

  updateConfiguration(eventSetup);

  // list of RCT channels to mask
  L1RCTChannelMask const &channelMask = eventSetup.getData(beginRunChannelMaskToken_);

  // list of Noisy RCT channels to mask
  L1RCTNoisyChannelMask const &hotChannelMask = eventSetup.getData(beginRunHotChannelMaskToken_);

  updateFedVector(channelMask, hotChannelMask, getFedVectorFromRunInfo(beginRunRunInfoToken_, eventSetup));
}

void L1RCTProducer::beginLuminosityBlock(edm::LuminosityBlock const &lumiSeg, const edm::EventSetup &context) {
  // check LS number every LS, if the checkOMDS flag is set AND it's the right
  // LS, update the FED vector from OMDS can pass the flag as the bool??  but
  // only check LS number if flag is true anyhow
  if (getFedsFromOmds) {
    throw cms::Exception("L1RCTProducer Configuration")
        << "L1RCTProducer is being run with the configuration parameter getFedsFromOmds set true. "
        << "Underlying Framework changes have broken the implementation of that option. "
        << "It was not fixed because we believe this option is no longer used or needed. "
        << "If you actually need this option, please report this failure to the Framework. "
        << "For more details see GitHub Issue 43697.";
    /*
    unsigned int nLumi = lumiSeg.luminosityBlock();  // doesn't even need the (unsigned int) cast
                                                     // because LuminosityBlockNumber_t is already
                                                     // an unsigned int
    // LS count starts at 1, want to be able to delay 0 LS's intuitively
    if (((nLumi - 1) == queryDelayInLS) ||
        (queryIntervalInLS > 0 &&
         nLumi % queryIntervalInLS == 0))  // to guard against problems if online DQM crashes; every 100
                                           // LS is ~20-30 minutes, not too big a load, hopefully not too
                                           // long between
    {
      //	  std::cout << "Lumi section for this FED vector update is " <<
      // nLumi << std::endl;

      // list of RCT channels to mask
      L1RCTChannelMask const &channelMask = context.getData(beginLumiChannelMaskToken_);

      // list of Noisy RCT channels to mask
      L1RCTNoisyChannelMask const &hotChannelMask = context.getData(beginLumiHotChannelMaskToken_);

      updateFedVector(channelMask, hotChannelMask, getFedVectorFromOmds(context));
    } else if (queryIntervalInLS <= 0) {
      // don't do interval checking... cout message??
    }
    */
  }
}

void L1RCTProducer::updateConfiguration(const edm::EventSetup &eventSetup) {
  // Refresh configuration information every event
  // Hopefully, this does not take too much time
  // There should be a call back function in future to
  // handle changes in configuration
  // parameters to configure RCT (thresholds, etc)
  const L1RCTParameters *r = &eventSetup.getData(rctParamsToken_);

  // SCALES

  // energy scale to convert eGamma output
  const L1CaloEtScale *s = &eventSetup.getData(emScaleToken_);

  // get energy scale to convert input from ECAL
  const L1CaloEcalScale *e = &eventSetup.getData(ecalScaleToken_);

  // get energy scale to convert input from HCAL
  const L1CaloHcalScale *h = &eventSetup.getData(hcalScaleToken_);

  // set scales
  rctLookupTables->setEcalScale(e);
  rctLookupTables->setHcalScale(h);

  rctLookupTables->setRCTParameters(r);
  rctLookupTables->setL1CaloEtScale(s);
}

void L1RCTProducer::updateFedVector(const L1RCTChannelMask &channelMask,
                                    const L1RCTNoisyChannelMask &hotChannelMask,
                                    const std::vector<int> &Feds)
// http://cmslxr.fnal.gov/lxr/source/FWCore/Framework/interface/EventSetup.h
{
  rctLookupTables->setNoisyChannelMask(&hotChannelMask);

  // Update the channel mask according to the FED VECTOR
  // This is the beginning of run. We delete the old
  // create the new and set it in the LUTs

  fedUpdatedMask = std::make_unique<L1RCTChannelMask>();
  // copy a constant object
  for (int i = 0; i < 18; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 28; k++) {
        fedUpdatedMask->ecalMask[i][j][k] = channelMask.ecalMask[i][j][k];
        fedUpdatedMask->hcalMask[i][j][k] = channelMask.hcalMask[i][j][k];
      }
      for (int k = 0; k < 4; k++) {
        fedUpdatedMask->hfMask[i][j][k] = channelMask.hfMask[i][j][k];
      }
    }
  }

  // so can create/initialize/assign const quantity in one line accounting for
  // if statement wikipedia says this is exactly what it's for:
  // http://en.wikipedia.org/wiki/%3F:#C.2B.2B

  //   std::cout << "Contents of ";
  //   std::cout << (getFromOmds ? "OMDS RunInfo" : "standard RunInfo");
  //   std::cout << " FED vector" << std::endl;
  //   printFedVector(Feds);

  bool useUpgradedHF = false;

  std::vector<int> caloFeds;  // pare down the feds to the interesting ones
  // is this unneccesary?
  // Mike B : This will decrease the find speed so better do it
  for (std::vector<int>::const_iterator cf = Feds.begin(); cf != Feds.end(); ++cf) {
    int fedNum = *cf;
    if ((fedNum > 600 && fedNum < 724) || fedNum == 1118 || fedNum == 1120 || fedNum == 1122)
      caloFeds.push_back(fedNum);

    if (fedNum == 1118 || fedNum == 1120 || fedNum == 1122)
      useUpgradedHF = true;
  }

  for (int cr = 0; cr < 18; ++cr) {
    for (crateSection cs = c_min; cs <= c_max; cs = crateSection(cs + 1)) {
      bool fedFound = false;

      // Try to find the FED
      std::vector<int>::iterator fv = std::find(caloFeds.begin(), caloFeds.end(), crateFED[cr][cs]);
      if (fv != caloFeds.end())
        fedFound = true;

      if (!fedFound) {
        int eta_min = 0;
        int eta_max = 0;
        bool phi_even[2] = {false};  //, phi_odd = false;
        bool ecal = false;

        switch (cs) {
          case ebEvenFed:
            eta_min = minBarrel;
            eta_max = maxBarrel;
            phi_even[0] = true;
            ecal = true;
            break;

          case ebOddFed:
            eta_min = minBarrel;
            eta_max = maxBarrel;
            phi_even[1] = true;
            ecal = true;
            break;

          case eeFed:
            eta_min = minEndcap;
            eta_max = maxEndcap;
            phi_even[0] = true;
            phi_even[1] = true;
            ecal = true;
            break;

          case hbheFed:
            eta_min = minBarrel;
            eta_max = maxEndcap;
            phi_even[0] = true;
            phi_even[1] = true;
            ecal = false;
            break;

          case hfFed:
            if (useUpgradedHF)
              break;

            eta_min = minHF;
            eta_max = maxHF;

            phi_even[0] = true;
            phi_even[1] = true;
            ecal = false;
            break;

          case hfFedUp:
            if (!useUpgradedHF)
              break;

            eta_min = minHF;
            eta_max = maxHF;

            phi_even[0] = true;
            phi_even[1] = true;
            ecal = false;
            break;

          default:
            break;
        }
        for (int ieta = eta_min; ieta <= eta_max; ++ieta) {
          if (ieta <= 28)  // barrel and endcap
            for (int even = 0; even <= 1; even++) {
              if (phi_even[even]) {
                if (ecal)
                  fedUpdatedMask->ecalMask[cr][even][ieta - 1] = true;
                else
                  fedUpdatedMask->hcalMask[cr][even][ieta - 1] = true;
              }
            }
          else
            for (int even = 0; even <= 1; even++)
              if (phi_even[even])
                fedUpdatedMask->hfMask[cr][even][ieta - 29] = true;
        }
      }
    }
  }

  rctLookupTables->setChannelMask(fedUpdatedMask.get());
}

const std::vector<int> L1RCTProducer::getFedVectorFromRunInfo(const edm::ESGetToken<RunInfo, RunInfoRcd> &token,
                                                              const edm::EventSetup &eventSetup) const {
  //  std::cout << "Getting FED vector from standard RunInfo object" <<
  //  std::endl;
  // get FULL FED vector from RUNINFO
  return eventSetup.getData(token).m_fed_in;
}

const std::vector<int> L1RCTProducer::getFedVectorFromOmds(const edm::EventSetup &eventSetup) const {
  //  std::cout << "Getting FED vector from my specific ES RunInfo object" <<
  //  std::endl;

  // get FULL FED vector from RunInfo object specifically created to have OMDS
  // fed vector
  edm::ESHandle<RunInfo> sum = eventSetup.getHandle(omdsRunInfoToken_);
  if (sum.isValid()) {
    return sum->m_fed_in;
  } else {
    return getFedVectorFromRunInfo(beginLumiRunInfoToken_, eventSetup);
  }
}

void L1RCTProducer::produce(edm::Event &event, const edm::EventSetup &eventSetup) {
  std::unique_ptr<L1CaloEmCollection> rctEmCands(new L1CaloEmCollection);
  std::unique_ptr<L1CaloRegionCollection> rctRegions(new L1CaloRegionCollection);

  if (!(ecalDigis.size() == hcalDigis.size() && hcalDigis.size() == bunchCrossings.size()))
    throw cms::Exception("BadInput") << "From what I see the number of your your ECAL input digi "
                                        "collections.\n"
                                     << "is different from the size of your HCAL digi input collections\n"
                                     << "or the size of your BX factor collection"
                                     << "They must be the same to correspond to the same Bxs\n"
                                     << "It does not matter if one of them is empty\n";

  // loop through and process each bx
  for (unsigned short sample = 0; sample < bunchCrossings.size(); sample++) {
    edm::Handle<EcalTrigPrimDigiCollection> ecal;
    edm::Handle<HcalTrigPrimDigiCollection> hcal;

    EcalTrigPrimDigiCollection ecalIn;
    HcalTrigPrimDigiCollection hcalIn;

    if (useHcal && event.getByLabel(hcalDigis[sample], hcal))
      hcalIn = *hcal;

    if (useEcal && event.getByLabel(ecalDigis[sample], ecal))
      ecalIn = *ecal;

    rct->digiInput(ecalIn, hcalIn);
    rct->processEvent();

    // Stuff to create
    for (int j = 0; j < 18; j++) {
      L1CaloEmCollection isolatedEGObjects = rct->getIsolatedEGObjects(j);
      L1CaloEmCollection nonisolatedEGObjects = rct->getNonisolatedEGObjects(j);
      for (int i = 0; i < 4; i++) {
        isolatedEGObjects.at(i).setBx(bunchCrossings[sample]);
        nonisolatedEGObjects.at(i).setBx(bunchCrossings[sample]);
        rctEmCands->push_back(isolatedEGObjects.at(i));
        rctEmCands->push_back(nonisolatedEGObjects.at(i));
      }
    }

    for (int i = 0; i < 18; i++) {
      std::vector<L1CaloRegion> regions = rct->getRegions(i);
      for (int j = 0; j < 22; j++) {
        regions.at(j).setBx(bunchCrossings[sample]);
        rctRegions->push_back(regions.at(j));
      }
    }
  }

  // putting stuff back into event
  event.put(std::move(rctEmCands));
  event.put(std::move(rctRegions));
}

// print contents of (FULL) FED vector
void L1RCTProducer::printFedVector(const std::vector<int> &fedVector) {
  std::cout << "Contents of given fedVector: ";
  std::copy(fedVector.begin(), fedVector.end(), std::ostream_iterator<int>(std::cout, ", "));
  std::cout << std::endl;
}

// print contents of RCT channel mask fedUpdatedMask
void L1RCTProducer::printUpdatedFedMask() {
  if (fedUpdatedMask != nullptr) {
    fedUpdatedMask->print(std::cout);
  } else {
    std::cout << "Trying to print contents of fedUpdatedMask, but it doesn't exist!" << std::endl;
  }
}

// print contents of RCT channel mask fedUpdatedMask
void L1RCTProducer::printUpdatedFedMaskVerbose() {
  if (fedUpdatedMask != nullptr) {
    // print contents of fedvector
    std::cout << "Contents of fedUpdatedMask: ";
    //       std::copy(fedUpdatedMask.begin(), fedUpdatedMask.end(),
    //       std::ostream_iterator<int>(std::cout, ", "));
    std::cout << "--> ECAL mask: " << std::endl;
    for (int i = 0; i < 18; i++) {
      for (int j = 0; j < 2; j++) {
        for (int k = 0; k < 28; k++) {
          std::cout << fedUpdatedMask->ecalMask[i][j][k] << ", ";
        }
      }
    }
    std::cout << "--> HCAL mask: " << std::endl;
    for (int i = 0; i < 18; i++) {
      for (int j = 0; j < 2; j++) {
        for (int k = 0; k < 28; k++) {
          std::cout << fedUpdatedMask->hcalMask[i][j][k] << ", ";
        }
      }
    }
    std::cout << "--> HF mask: " << std::endl;
    for (int i = 0; i < 18; i++) {
      for (int j = 0; j < 2; j++) {
        for (int k = 0; k < 4; k++) {
          std::cout << fedUpdatedMask->hfMask[i][j][k] << ", ";
        }
      }
    }

    std::cout << std::endl;
  } else {
    // print error message
    std::cout << "Trying to print contents of fedUpdatedMask, but it doesn't exist!" << std::endl;
  }
}
