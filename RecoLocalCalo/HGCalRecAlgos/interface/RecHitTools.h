#ifndef __RecoLocalCalo_HGCalRecAlgos_RecHitTools_h__
#define __RecoLocalCalo_HGCalRecAlgos_RecHitTools_h__

#include <array>
#include <cmath>
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/ForwardDetId/interface/HFNoseDetId.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

class CaloGeometry;
class CaloSubdetectorGeometry;
class DetId;

namespace edm {
  class Event;
  class EventSetup;
}  // namespace edm

namespace hgcal {
  class RecHitTools {
  public:
    struct siliconWaferInfo {
      int32_t type, partialType, orientation, placementIndex, cassette;
      siliconWaferInfo(int32_t t = 0, int32_t p = 0, int32_t o = 0, int32_t i = 0, int32_t c = 0)
          : type(t), partialType(p), orientation(o), placementIndex(i), cassette(c) {}
    };
    struct scintillatorTileInfo {
      int32_t type, sipm, cassette;
      scintillatorTileInfo(int32_t t = 0, int32_t s = 0, int32_t c = 0) : type(t), sipm(s), cassette(c) {}
    };
    RecHitTools()
        : geom_(nullptr),
          eeOffset_(0),
          fhOffset_(0),
          bhFirstLayer_(0),
          bhOffset_(0),
          fhLastLayer_(0),
          noseLastLayer_(0),
          hcalBarrelFirstLayer_(1),
          hcalBarrelLastLayer_(4),
          ecalBarrelFirstLayer_(0),
          ecalBarrelLastLayer_(0),
          geometryType_(0) {}
    ~RecHitTools() {}

    void setGeometry(CaloGeometry const&);
    const CaloSubdetectorGeometry* getSubdetectorGeometry(const DetId& id) const;

    GlobalPoint getPosition(const DetId& id) const;
    GlobalPoint getPositionLayer(int layer, bool nose = false, bool barrel = false) const;
    // zside returns +/- 1
    int zside(const DetId& id) const;

    std::float_t getSiThickness(const DetId&) const;
    std::float_t getRadiusToSide(const DetId&) const;
    int getSiThickIndex(const DetId&) const;

    std::pair<float, float> getScintDEtaDPhi(const DetId&) const;

    unsigned int getLayer(DetId::Detector type, bool nose = false) const;
    unsigned int getLayer(ForwardSubdetector type) const;
    unsigned int getLayer(const DetId&) const;
    unsigned int getLayerWithOffset(const DetId&) const;
    int getCellType(const DetId& id) const;
    std::pair<int, int> getWafer(const DetId&) const;
    std::pair<int, int> getCell(const DetId&) const;

    bool isHalfCell(const DetId&) const;

    bool isSilicon(const DetId&) const;
    bool isScintillator(const DetId&) const;
    bool isScintillatorFine(const DetId& id) const;
    bool isBarrel(const DetId&) const;

    bool isOnlySilicon(const unsigned int layer) const;

    // 4-vector helper functions using GlobalPoint
    float getEta(const GlobalPoint& position, const float& vertex_z = 0.) const;
    float getPhi(const GlobalPoint& position) const;
    float getPt(const GlobalPoint& position, const float& hitEnergy, const float& vertex_z = 0.) const;

    // 4-vector helper functions using DetId
    float getEta(const DetId& id, const float& vertex_z = 0.) const;
    float getPhi(const DetId& id) const;
    float getPt(const DetId& id, const float& hitEnergy, const float& vertex_z = 0.) const;
    int getScintMaxIphi(const DetId& id) const;

    inline const CaloGeometry* getGeometry() const { return geom_; };
    unsigned int lastLayerEE(bool nose = false) const { return (nose ? HFNoseDetId::HFNoseLayerEEmax : fhOffset_); }
    unsigned int lastLayerFH() const { return fhLastLayer_; }
    unsigned int firstLayerBH() const { return bhFirstLayer_; }
    unsigned int lastLayerBH() const { return bhLastLayer_; }
    unsigned int lastLayer(bool nose = false) const { return (nose ? noseLastLayer_ : bhLastLayer_); }
    unsigned int lastLayerECAL() const { return ecalBarrelLastLayer_; }
    unsigned int lastLayerBarrel() const { return hcalBarrelLastLayer_; }
    std::pair<uint32_t, uint32_t> firstAndLastLayer(DetId::Detector det, int subdet) const;
    unsigned int maxNumberOfWafersPerLayer(bool nose = false) const {
      return (nose ? maxNumberOfWafersNose_ : maxNumberOfWafersPerLayer_);
    }
    inline int getScintMaxIphi() const { return bhMaxIphi_; }
    inline int getGeometryType() const { return geometryType_; }
    bool maskCell(const DetId& id, int corners = 3) const;

    // Informaion of the wafer/tile
    siliconWaferInfo getWaferInfo(const DetId& id) const;
    scintillatorTileInfo getTileInfo(const DetId& id) const;
    int getWaferTypes(DetId::Detector det, int subdet = ForwardSubdetector::ForwardEmpty) const;
    std::vector<double> getSiThickness(DetId::Detector det, int subdet = ForwardSubdetector::ForwardEmpty) const;

  private:
    const CaloGeometry* geom_;
    unsigned int eeOffset_, fhOffset_, bhFirstLayer_, bhLastLayer_, bhOffset_, fhLastLayer_, noseLastLayer_;
    unsigned int hcalBarrelFirstLayer_, hcalBarrelLastLayer_, ecalBarrelFirstLayer_, ecalBarrelLastLayer_;
    unsigned int maxNumberOfWafersPerLayer_, maxNumberOfWafersNose_;
    int geometryType_;
    int bhMaxIphi_;
  };
}  // namespace hgcal

#endif
