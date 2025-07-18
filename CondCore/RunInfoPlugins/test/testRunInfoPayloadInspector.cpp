#include <iostream>
#include <sstream>
#include "CondCore/Utilities/interface/PayloadInspector.h"
#include "CondCore/RunInfoPlugins/plugins/RunInfo_PayloadInspector.cc"
#include "CondCore/RunInfoPlugins/plugins/LHCInfoPerFill_PayloadInspector.cc"
#include "CondCore/RunInfoPlugins/plugins/LHCInfoPerLS_PayloadInspector.cc"
#include "FWCore/PluginManager/interface/PluginManager.h"
#include "FWCore/PluginManager/interface/standard.h"
#include "FWCore/ServiceRegistry/interface/ServiceRegistry.h"

int main(int argc, char** argv) {
  Py_Initialize();
  edmplugin::PluginManager::Config config;
  edmplugin::PluginManager::configure(edmplugin::standard::config());

  std::vector<edm::ParameterSet> psets;
  edm::ParameterSet pSet;
  pSet.addParameter("@service_type", std::string("SiteLocalConfigService"));
  psets.push_back(pSet);
  edm::ServiceToken servToken(edm::ServiceRegistry::createSet(psets));
  edm::ServiceRegistry::Operate operate(servToken);

  std::string connectionString("frontier://FrontierProd/CMS_CONDITIONS");

  std::string tag = "runinfo_31X_hlt";
  cond::Time_t start = static_cast<unsigned long long>(311950);
  cond::Time_t end = static_cast<unsigned long long>(312237);

  std::cout << "## RunInfo testing" << std::endl;

  RunInfoTest histo0;
  histo0.process(connectionString, PI::mk_input(tag, end, end));
  std::cout << histo0.data() << std::endl;

  RunInfoParameters histo1;
  histo1.process(connectionString, PI::mk_input(tag, end, end));
  std::cout << histo1.data() << std::endl;

  RunInfoBFieldHistory histo2;
  histo2.process(connectionString, PI::mk_input(tag, start, end));
  std::cout << histo2.data() << std::endl;

  std::cout << "## LHCInfoPerFill testing" << std::endl;

  tag = "LHCInfoPerFill_duringFill_hlt_v1";
  start = static_cast<unsigned long long>(1686852700471354);
  end = static_cast<unsigned long long>(1686852700471354);

  LHCInfoPerFill_Display histo3;
  histo3.process(connectionString, PI::mk_input(tag, end, end));
  std::cout << histo3.data() << std::endl;

  std::cout << "## LHCInfoPerLS testing" << std::endl;

  tag = "LHCInfoPerLS_duringFill_hlt_v1";
  start = static_cast<unsigned long long>(1690138350452868);
  end = static_cast<unsigned long long>(1690138350452868);

  LHCInfoPerLS_Display histo4;
  histo4.process(connectionString, PI::mk_input(tag, end, end));
  std::cout << histo4.data() << std::endl;

  Py_Finalize();
}
