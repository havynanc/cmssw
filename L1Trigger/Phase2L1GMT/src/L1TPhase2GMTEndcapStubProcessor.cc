#include "L1Trigger/Phase2L1GMT/interface/L1TPhase2GMTEndcapStubProcessor.h"
#include "L1Trigger/L1TMuon/interface/MuonTriggerPrimitive.h"
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>

L1TPhase2GMTEndcapStubProcessor::L1TPhase2GMTEndcapStubProcessor() : minBX_(-3), maxBX_(3) {}

L1TPhase2GMTEndcapStubProcessor::L1TPhase2GMTEndcapStubProcessor(const edm::ParameterSet& iConfig)
    : minBX_(iConfig.getParameter<int>("minBX")),
      maxBX_(iConfig.getParameter<int>("maxBX")),
      coord1LSB_(iConfig.getParameter<double>("coord1LSB")),
      coord2LSB_(iConfig.getParameter<double>("coord2LSB")),
      eta1LSB_(iConfig.getParameter<double>("eta1LSB")),
      eta2LSB_(iConfig.getParameter<double>("eta2LSB")),
      etaMatch_(iConfig.getParameter<double>("etaMatch")),
      phiMatch_(iConfig.getParameter<double>("phiMatch")),
      verbose_(iConfig.getParameter<unsigned int>("verbose")) {}

L1TPhase2GMTEndcapStubProcessor::~L1TPhase2GMTEndcapStubProcessor() {}

l1t::MuonStub L1TPhase2GMTEndcapStubProcessor::buildCSCOnlyStub(const CSCDetId& detid,
                                                                const CSCCorrelatedLCTDigi& digi,
                                                                const L1TMuon::GeometryTranslator* translator,
                                                                unsigned int tag) {
  int endcap = detid.endcap();
  int station = detid.station();
  int chamber = detid.chamber();
  int ring = detid.ring();

  L1TMuon::TriggerPrimitive primitive(detid, digi);

  const GlobalPoint& gp = translator->getGlobalPoint(primitive);

  int phi = int(gp.phi().value() / coord1LSB_);
  while (phi < -std::pow(2, Phase2L1GMT::BITSSTUBCOORD-1)) {
    phi = phi + std::pow(2, Phase2L1GMT::BITSSTUBCOORD);
  }
  while (phi >= std::pow(2, Phase2L1GMT::BITSSTUBCOORD-1)) {
    phi = phi - std::pow(2, Phase2L1GMT::BITSSTUBCOORD);
  }
  int eta1 = int(gp.eta() / eta1LSB_);

  int etaRegion = 0;
  // endcap: 1=forward (+Z), 2=backward(-Z) from CSCDetId
  int sign = endcap == 2 ? -1 : 1;

  if (((station == 1) && (ring == 3)) || ((station == 2) && (ring == 2)) || ((station == 3) && (ring == 2)))
    etaRegion = sign * 3;
  else if (((station == 1) && (ring == 2)) || ((station == 4) && (ring == 2)))
    etaRegion = sign * 4;
  else if (ring == 1)
    etaRegion = sign * 5;

  int sector = fabs(chamber);

  int bx = digi.getBX() - 8;
  int quality = 1;

  uint tfLayer = 0;
  if ((station == 1) && (ring == 1))  //ME1/1
    tfLayer = 1;
  else if (station == 4)  //ME4/x
    tfLayer = 3;
  else if (station == 3)  //ME3/x
    tfLayer = 5;
  else if (station == 2)  //ME2/x
    tfLayer = 6;
  else if ((station == 1) && ((ring == 2) || (ring == 3))) //ME1/2 and ME1/3
    tfLayer = 10;

  l1t::MuonStub stub(etaRegion, sector, tfLayer, tfLayer, phi, 0, tag, bx, quality, eta1, 0, 1, 0);

  stub.setOfflineQuantities(gp.phi().value(), 0.0, gp.eta(), 0.0);
  return stub;
}

l1t::MuonStub L1TPhase2GMTEndcapStubProcessor::buildRPCOnlyStub(const RPCDetId& detid,
                                                                const RPCDigi& digi,
                                                                const L1TMuon::GeometryTranslator* translator) {
  L1TMuon::TriggerPrimitive primitive(detid, digi);
  const GlobalPoint& gp = translator->getGlobalPoint(primitive);

  int phi2 = int(gp.phi().value() / coord2LSB_);
  while (phi2 < -std::pow(2, Phase2L1GMT::BITSSTUBCOORD-1)) {
    phi2 = phi2 + std::pow(2, Phase2L1GMT::BITSSTUBCOORD);
  }
  while (phi2 >= std::pow(2, Phase2L1GMT::BITSSTUBCOORD-1)) {
    phi2 = phi2 - std::pow(2, Phase2L1GMT::BITSSTUBCOORD);
  }
  int eta2 = int(gp.eta() / eta2LSB_);

  int etaRegion = 0;
  int sector = (detid.sector() - 1) * 6 + detid.subsector();
  int station = detid.station();
  bool tag = detid.trIndex();
  int bx = digi.bx();
  int quality = 2;

  int ring = detid.ring();
  int sign = eta2 < 0 ? -1 : 1;

  uint tfLayer = 0;

  if (((station == 1) && (ring == 3)) || (station == 2) || ((station == 3) && (ring == 2)) || ((station == 3) && (ring == 3))) //RE1/3, RE2/x, RE3/2, RE3/3
    etaRegion = sign * 3;
  else if (((station == 1) && (ring == 2)) || ((station == 4) && (ring == 2)) || ((station == 4) && (ring == 3))) //RE1/2, RE4/2, RE4/3
    etaRegion = sign * 4;
  else if (ring == 1) //REx/1
    etaRegion = sign * 5;

  if (station == 4) //RE4/x
    tfLayer = 2;
  else if (station == 3) //RE3/x
    tfLayer = 4;
  else if (station == 2) //RE2/x
    tfLayer = 7;
  else if ((station == 1) && (ring == 3)) //RE1/3
    tfLayer = 8;
  else if ((station == 1) && (ring == 2)) { //RE1/2
    // In the barrel, layer() differentiates between inner (1) and outer (2) RPCs.
    // Currently layer()=1 always in the encap RPCs but this is likely how the two RE1/2's are supposed to be differentiated.
    // For now, the outer RE1/2 is absorbed into layer 9 since layer()==2 is never satisfied.
    if (detid.layer() == 1)
      tfLayer = 9;
    else if (detid.layer() == 2)
      tfLayer = 8;
  }

  l1t::MuonStub stub(etaRegion, sector, tfLayer, tfLayer, phi2, 0, tag, bx, quality, eta2, 0, 2, 0);
  stub.setOfflineQuantities(0.0, gp.phi().value(), 0.0, gp.eta()); //offline quantities will be an independent reminder of csc or rpc
  //std::cout<<"RPC detID: " << "station:" << detid.station() << " sector:" << detid.sector() << " subsector:" << detid.subsector() << " etaRegion:"<< etaRegion << " tfLayer" << tfLayer << " layer:" << detid.layer() << " roll:" << detid.roll() << " trIndex:" << detid.trIndex() << " chamberId:" << detid.chamberId() << std::endl;
  
  return stub;
}

l1t::MuonStub L1TPhase2GMTEndcapStubProcessor::buildME0OnlyStub(const GEMDetId& detid,
                                                                const ME0TriggerDigi& digi,
                                                                const L1TMuon::GeometryTranslator* translator,
                                                                unsigned int tag) {

  L1TMuon::TriggerPrimitive primitive(detid, digi);
  const GlobalPoint& gp = translator->getGlobalPoint(primitive);

  int phi = int(gp.phi().value() / coord1LSB_);
  while (phi < -std::pow(2, Phase2L1GMT::BITSSTUBCOORD-1)) {
    phi = phi + std::pow(2, Phase2L1GMT::BITSSTUBCOORD);
  }
  while (phi >= std::pow(2, Phase2L1GMT::BITSSTUBCOORD-1)) {
    phi = phi - std::pow(2, Phase2L1GMT::BITSSTUBCOORD);
  }
  int eta1 = int(gp.eta() / eta1LSB_);

  int sign = detid.region();
  int chamber = detid.chamber();
  int etaRegion = sign * 5;
  int bx = digi.getBX() - 8;  // ME0TriggerDigi bx is centered at BX8, similar to CSC
  int quality = 1;
  uint tfLayer = 10;

  l1t::MuonStub stub(etaRegion, chamber, tfLayer, tfLayer, phi, 0, tag, bx, quality, eta1, 0, 1, 0);
  stub.setOfflineQuantities(gp.phi().value(), 0.0, gp.eta(), 0.0);

  return stub;
}

l1t::MuonStub L1TPhase2GMTEndcapStubProcessor::buildGEMOnlyStub(const GEMDetId& detid,
                                                                const GEMPadDigiCluster& digi,
                                                                const L1TMuon::GeometryTranslator* translator,
                                                                unsigned int tag) {
  L1TMuon::TriggerPrimitive primitive(detid, digi);
  const GlobalPoint& gp = translator->getGlobalPoint(primitive);

  int phi = int(gp.phi().value() / coord1LSB_);
  while (phi < -std::pow(2, Phase2L1GMT::BITSSTUBCOORD-1)) {
    phi = phi + std::pow(2, Phase2L1GMT::BITSSTUBCOORD);
  }
  while (phi >= std::pow(2, Phase2L1GMT::BITSSTUBCOORD-1)) {
    phi = phi - std::pow(2, Phase2L1GMT::BITSSTUBCOORD);
  }
  int eta1 = int(gp.eta() / eta1LSB_);

  int sign = detid.region();
  int chamber = detid.chamber();
  int station = detid.station();
  int etaRegion = sign * 5;
  int bx = digi.bx();
  int quality = 1;

  uint tfLayer = 0;
  if (station == 1)
    tfLayer = 0;
  else if (station == 2)
    tfLayer = 7;

  l1t::MuonStub stub(etaRegion, chamber, tfLayer, tfLayer, phi, 0, tag, bx, quality, eta1, 0, 1, 0);
  stub.setOfflineQuantities(gp.phi().value(), 0.0, gp.eta(), 0.0);

  return stub;
}


l1t::MuonStubCollection L1TPhase2GMTEndcapStubProcessor::clusterRPCStubs(const l1t::MuonStubCollection& rpcStubs) {

  l1t::MuonStubCollection out;
  l1t::MuonStubCollection inRPC = rpcStubs;

  while (!inRPC.empty()) {
    l1t::MuonStubCollection freeRPC;

    int nRPC = 1;
    float phiF = inRPC[0].offline_coord2();
    float etaF = inRPC[0].offline_eta2();
    int phi = inRPC[0].coord1();
    int eta = inRPC[0].eta1();

    for (unsigned i = 1; i < inRPC.size(); ++i) {
      if (fabs(deltaPhi(inRPC[0].offline_coord2(), inRPC[i].offline_coord2())) < phiMatch_ &&
          inRPC[0].depthRegion() == inRPC[i].depthRegion() &&
          fabs(inRPC[0].offline_eta2() - inRPC[i].offline_eta2()) < etaMatch_ &&
          inRPC[0].bxNum() == inRPC[i].bxNum()) {
        phiF += inRPC[i].offline_coord2();
        etaF += inRPC[i].offline_eta2();
        phi += inRPC[i].coord1();
        eta += inRPC[i].eta1();
        nRPC++;
      } else {
        freeRPC.push_back(inRPC[i]);
      }
    }
    l1t::MuonStub stub(inRPC[0].etaRegion(),
                       inRPC[0].phiRegion(),
                       inRPC[0].depthRegion(),
                       inRPC[0].tfLayer(),
                       phi / nRPC,
		       0,
                       0,
                       inRPC[0].bxNum(),
                       2,
                       eta / nRPC,
		       0,
                       2,
                       0);
    stub.setOfflineQuantities(0.0, phiF / nRPC, 0.0, etaF / nRPC);
    out.push_back(stub);
    inRPC = freeRPC;
  };
  return out;

}



l1t::MuonStubCollection L1TPhase2GMTEndcapStubProcessor::makeStubs(
    const MuonDigiCollection<CSCDetId, CSCCorrelatedLCTDigi>& csc,
    const MuonDigiCollection<RPCDetId, RPCDigi>& cleaned,
    const MuonDigiCollection<GEMDetId, ME0TriggerDigi>& me0,
    const MuonDigiCollection<GEMDetId, GEMPadDigiCluster>& gem,
    const L1TMuon::GeometryTranslator* t,
    const edm::EventSetup& iSetup) {

  l1t::MuonStubCollection cscStubs;
  auto chamber = csc.begin();
  auto chend = csc.end();
  for (; chamber != chend; ++chamber) {
    auto digi = (*chamber).second.first;
    auto dend = (*chamber).second.second;
    unsigned int tag = 0;
    for (; digi != dend; ++digi) {
      l1t::MuonStub stub = buildCSCOnlyStub((*chamber).first, *digi, t, tag);
      tag = tag + 1;
      if (stub.bxNum() >= minBX_ && stub.bxNum() <= maxBX_)
        cscStubs.push_back(stub);
    }
  }

  l1t::MuonStubCollection rpcStubs;
  auto rpcchamber = cleaned.begin();
  auto rpcchend = cleaned.end();
  for (; rpcchamber != rpcchend; ++rpcchamber) {
    if ((*rpcchamber).first.region() == 0)
      continue;
    auto digi = (*rpcchamber).second.first;
    auto dend = (*rpcchamber).second.second;
    for (; digi != dend; ++digi) {
      l1t::MuonStub stub = buildRPCOnlyStub((*rpcchamber).first, *digi, t);
      if (stub.bxNum() >= minBX_ && stub.bxNum() <= maxBX_)
        rpcStubs.push_back(stub);
    }
  }

  l1t::MuonStubCollection me0Stubs;
  auto me0chamber = me0.begin();
  auto me0chend = me0.end();
  for (; me0chamber != me0chend; ++me0chamber) {
    auto digi = (*me0chamber).second.first;
    auto dend = (*me0chamber).second.second;
    unsigned int tag = 0;
    for (; digi != dend; ++digi) {
      l1t::MuonStub stub = buildME0OnlyStub((*me0chamber).first, *digi, t, tag);
      tag = tag + 1; //need to understand this for sure
      if (stub.bxNum() >= minBX_ && stub.bxNum() <= maxBX_)
        me0Stubs.push_back(stub);
    }
  }

  l1t::MuonStubCollection gemStubs;
  auto gemchamber = gem.begin();
  auto gemchend = gem.end();
  for (; gemchamber != gemchend; ++gemchamber) {
    auto digi = (*gemchamber).second.first;
    auto dend = (*gemchamber).second.second;
    unsigned int tag = 0;
    for (; digi != dend; ++digi) {
      l1t::MuonStub stub = buildGEMOnlyStub((*gemchamber).first, *digi, t, tag);
      tag = tag + 1; //need to understand this for sure
      if (stub.bxNum() >= minBX_ && stub.bxNum() <= maxBX_)
        gemStubs.push_back(stub);
    }
  }


  l1t::MuonStubCollection clusteredRPCStubs = clusterRPCStubs(rpcStubs);
  l1t::MuonStubCollection combinedStubs;
  for (const auto& stub : cscStubs) {
    combinedStubs.push_back(stub);
  }
  for (const auto& stub : clusteredRPCStubs) {
    combinedStubs.push_back(stub);
  }
  for (const auto& stub : me0Stubs) {
    combinedStubs.push_back(stub);
  }
  for (const auto& stub : gemStubs) {
    combinedStubs.push_back(stub);
  }

  if (verbose_) {
    edm::LogInfo("EndcapStub") << "CSC Stubs";
    for (const auto& stub : cscStubs)
      edm::LogInfo("EndcapStub") << "CSC Stub bx=" << stub.bxNum() << " TF=" << stub.tfLayer()
                                 << " etaRegion=" << stub.etaRegion() << " phiRegion=" << stub.phiRegion()
                                 << " depthRegion=" << stub.depthRegion() << "  coord1=" << stub.offline_coord1() << ","
                                 << stub.coord1() << " coord2=" << stub.offline_coord2() << "," << stub.coord2()
                                 << " eta1=" << stub.offline_eta1() << "," << stub.eta1()
                                 << " eta2=" << stub.offline_eta2() << "," << stub.eta2()
                                 << " quality=" << stub.quality() << " etaQuality=" << stub.etaQuality();

    edm::LogInfo("EndcapStub") << "RPC Stubs";
    for (const auto& stub : rpcStubs)
      edm::LogInfo("EndcapStub") << "RPC Stub bx=" << stub.bxNum() << " TF=" << stub.tfLayer()
                                 << " etaRegion=" << stub.etaRegion() << " phiRegion=" << stub.phiRegion()
                                 << " depthRegion=" << stub.depthRegion() << "  coord1=" << stub.offline_coord1() << ","
                                 << stub.coord1() << " coord2=" << stub.offline_coord2() << "," << stub.coord2()
                                 << " eta1=" << stub.offline_eta1() << "," << stub.eta1()
                                 << " eta2=" << stub.offline_eta2() << "," << stub.eta2()
                                 << " quality=" << stub.quality() << " etaQuality=" << stub.etaQuality();

    for (const auto& stub : combinedStubs)
      edm::LogInfo("EndcapStub") << "Combined Stub bx=" << stub.bxNum() << " TF=" << stub.tfLayer()
                                 << " etaRegion=" << stub.etaRegion() << " phiRegion=" << stub.phiRegion()
                                 << " depthRegion=" << stub.depthRegion() << "  coord1=" << stub.offline_coord1() << ","
                                 << stub.coord1() << " coord2=" << stub.offline_coord2() << "," << stub.coord2()
                                 << " eta1=" << stub.offline_eta1() << "," << stub.eta1()
                                 << " eta2=" << stub.offline_eta2() << "," << stub.eta2()
                                 << " quality=" << stub.quality() << " etaQuality=" << stub.etaQuality();
  }

  return combinedStubs;
}
