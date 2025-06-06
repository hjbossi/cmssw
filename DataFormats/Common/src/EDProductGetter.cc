// -*- C++ -*-
//
// Package:     EDProduct
// Class  :     EDProductGetter
//
// Implementation:
//     <Notes on implementation>
//
// Original Author:  Chris Jones
//         Created:  Tue Nov  1 15:06:41 EST 2005
//

// system include files
#include <atomic>

// user include files
#include "DataFormats/Common/interface/EDProductGetter.h"
#include "DataFormats/Provenance/interface/ProductID.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/Likely.h"

namespace edm {
  //
  // constants, enums and typedefs
  //

  //
  // static data member definitions
  //

  //
  // constructors and destructor
  //
  EDProductGetter::EDProductGetter() {}

  // EDProductGetter::EDProductGetter(EDProductGetter const& rhs)
  // {
  //    // do actual copying here;
  // }

  EDProductGetter::~EDProductGetter() {}

  //
  // assignment operators
  //
  // EDProductGetter const& EDProductGetter::operator=(EDProductGetter const& rhs)
  // {
  //   //An exception safe implementation is
  //   EDProductGetter temp(rhs);
  //   swap(rhs);
  //
  //   return *this;
  // }

  //
  // member functions
  //

  //
  // const member functions
  //

  //
  // static member functions
  //

  EDProductGetter const* mustBeNonZero(EDProductGetter const* prodGetter,
                                       std::string refType,
                                       ProductID const& productID) {
    if (prodGetter != nullptr)
      return prodGetter;
    throw Exception(errors::InvalidReference, refType)
        << "Attempt to construct a " << refType << " with ProductID " << productID << "\n"
        << "but with a null pointer to a product getter.\n"
        << "The product getter pointer passed to the constructor must refer\n"
        << "to a real getter, such as an EventPrincipal.\n";
  }

  thread_local EDProductGetter const* s_productGetter = nullptr;
  static std::atomic<EDProductGetter const*> s_multiThreadProductGetter{nullptr};
  EDProductGetter const* EDProductGetter::switchProductGetter(EDProductGetter const* iNew) {
    //std::cout <<"switch from "<<s_productGetter<<" to "<<iNew<<std::endl;
    EDProductGetter const* old = s_productGetter;
    s_productGetter = iNew;
    return old;
  }

  void EDProductGetter::setMultiThreadProductGetter(EDProductGetter const* prodGetter) {
    EDProductGetter const* expected = nullptr;
    while (not s_multiThreadProductGetter.compare_exchange_strong(expected, prodGetter, std::memory_order_acq_rel)) {
      expected = nullptr;
    };
  }

  void EDProductGetter::unsetMultiThreadProductGetter() {
    s_multiThreadProductGetter.store(nullptr, std::memory_order_release);
  }

  void EDProductGetter::assignEDProductGetter(EDProductGetter const*& iGetter) {
    //std::cout <<"assign "<<s_productGetter<<std::endl;
    if LIKELY (s_productGetter != nullptr) {
      iGetter = s_productGetter;
      return;
    }
    iGetter = s_multiThreadProductGetter.load(std::memory_order_acquire);
  }

}  // namespace edm
