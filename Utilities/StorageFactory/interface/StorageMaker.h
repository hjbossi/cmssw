#ifndef STORAGE_FACTORY_STORAGE_MAKER_H
#define STORAGE_FACTORY_STORAGE_MAKER_H

#include "Utilities/StorageFactory/interface/IOTypes.h"
#include <string>
#include <memory>

namespace edm::storage {
  class Storage;
  class StorageMaker {
  public:
    struct AuxSettings {
      unsigned int timeout = 0;
      unsigned int debugLevel = 0;

      AuxSettings &setDebugLevel(unsigned int iLevel) {
        debugLevel = iLevel;
        return *this;
      }

      AuxSettings &setTimeout(unsigned int iTime) {
        timeout = iTime;
        return *this;
      }
    };

    enum class UseLocalFile { kYes, kCheckFromPath, kNo };

    StorageMaker() = default;
    virtual ~StorageMaker() = default;
    // implicit copy constructor
    // implicit assignment operator

    virtual void stagein(const std::string &proto, const std::string &path, const AuxSettings &aux) const;
    virtual std::unique_ptr<Storage> open(const std::string &proto,
                                          const std::string &path,
                                          int mode,
                                          const AuxSettings &aux) const = 0;
    virtual bool check(const std::string &proto,
                       const std::string &path,
                       const AuxSettings &aux,
                       IOOffset *size = nullptr) const;

    // If the file being read in the end is a local file or not
    virtual UseLocalFile usesLocalFile() const = 0;
  };
}  // namespace edm::storage
#endif  // STORAGE_FACTORY_STORAGE_MAKER_H
