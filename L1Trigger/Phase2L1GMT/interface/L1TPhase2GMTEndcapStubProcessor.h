
#ifndef L1TPHASE2GMTENDCAPSTUBPROCESSOR
#define L1TPHASE2GMTENDCAPSTUBPROCESSOR

#include "DataFormats/L1TMuonPhase2/interface/MuonStub.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "DataFormats/MuonData/interface/MuonDigiCollection.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "L1Trigger/L1TMuon/interface/GeometryTranslator.h"
#include "DataFormats/RPCDigi/interface/RPCDigi.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "L1Trigger/L1TMuon/interface/MuonTriggerPrimitive.h"
#include "L1Trigger/L1TTwinMux/interface/RPCHitCleaner.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/GEMDigi/interface/GEMPadDigiCluster.h"
#include "DataFormats/GEMDigi/interface/ME0TriggerDigi.h"

class L1TPhase2GMTEndcapStubProcessor {
public:
  L1TPhase2GMTEndcapStubProcessor();
  L1TPhase2GMTEndcapStubProcessor(const edm::ParameterSet&);
  ~L1TPhase2GMTEndcapStubProcessor();

  l1t::MuonStubCollection makeStubs(const MuonDigiCollection<CSCDetId, CSCCorrelatedLCTDigi>& csc,
                                    const MuonDigiCollection<RPCDetId, RPCDigi>& rpc,
				    const MuonDigiCollection<GEMDetId, ME0TriggerDigi>& me0,
				    const MuonDigiCollection<GEMDetId, GEMPadDigiCluster>& gem,
                                    const L1TMuon::GeometryTranslator* t,
                                    const edm::EventSetup& iSetup);

private:
  l1t::MuonStub buildCSCOnlyStub(const CSCDetId&,
                                 const CSCCorrelatedLCTDigi&,
                                 const L1TMuon::GeometryTranslator*,
                                 unsigned int);
  l1t::MuonStub buildRPCOnlyStub(const RPCDetId&, const RPCDigi&, const L1TMuon::GeometryTranslator*);
  l1t::MuonStub buildME0OnlyStub(const GEMDetId&, const ME0TriggerDigi&, const L1TMuon::GeometryTranslator*, unsigned int);
  l1t::MuonStub buildGEMOnlyStub(const GEMDetId&, const GEMPadDigiCluster&, const L1TMuon::GeometryTranslator*, unsigned int);
  l1t::MuonStubCollection clusterRPCStubs(const l1t::MuonStubCollection&);

  int minBX_;
  int maxBX_;
  double coord1LSB_;
  double coord2LSB_;
  double eta1LSB_;
  double eta2LSB_;
  double etaMatch_;
  double phiMatch_;
  bool verbose_;
};

#endif
