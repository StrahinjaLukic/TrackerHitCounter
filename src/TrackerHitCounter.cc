#include "DD4hep/Detector.h"
#include <DD4hep/DetType.h>
#include <TrackerHitCounter.h>
#include "DDRec/DetectorData.h"
#include <IMPL/LCCollectionVec.h>
#include <EVENT/SimTrackerHit.h>
#include <UTIL/CellIDDecoder.h>
#include <vector>

using namespace marlin ;


TrackerHitCounter aTrackerHitCounter ;


TrackerHitCounter::TrackerHitCounter() : Processor("TrackerHitCounter"),
    hitCounters(), m_trkHitCollNames()
{

  // modify processor description
  _description = "TrackerHitCounter counts SimTrackerHits in tracker detector elements "
          "and reports the number of hits per unit area." ;

  StringVec defaultTrkHitCollections;
  defaultTrkHitCollections.push_back(std::string("VXDCollection"));
  defaultTrkHitCollections.push_back(std::string("SITCollection"));
  defaultTrkHitCollections.push_back(std::string("FTDCollection"));
  defaultTrkHitCollections.push_back(std::string("TPCCollection"));
  defaultTrkHitCollections.push_back(std::string("SETCollection"));

  registerProcessorParameter("TrkHitCollections" ,
                             "Tracker hit collections that will be analysed",
                             m_trkHitCollNames ,
                             defaultTrkHitCollections ) ;
}

void TrackerHitCounter::init() { 

  // usually a good idea to
  printParameters() ;

  // v01-19-05 namespace
  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();

  const std::vector< dd4hep::DetElement > &detElements = theDetector.detectors("tracker", true);

  for(dd4hep::DetElement element : detElements) {

    streamlog_out(MESSAGE) << "*********************************************************************\n";
          streamlog_out(MESSAGE) << "Detector element \'" << element.name() << "\' of type \'"
        << element.type() << "\':\n";

    hitCounters[element.id()] = new SystemHitCounter;

    try {
      streamlog_out(DEBUG) << "Trying ZPlanarData.\n";
      dd4hep::rec::ZPlanarData* layering = element.extension<dd4hep::rec::ZPlanarData>() ;

      // We assume that layer numbering always starts from 0.
      int ilay = 0;
      for ( auto layer : layering->layers ) {
          streamlog_out(MESSAGE) << "  Layer " << ilay+1 << ":\n";
          double ls = layer.zHalfSensitive*2;
          double ws = layer.widthSensitive;
          int nModules = layer.ladderNumber;
          double area = ls*ws*nModules;
          (*hitCounters.at(element.id()))[ilay] = new LayerHitCounter(area);
          streamlog_out(MESSAGE) << "    Length of sensitive area: " << ls/dd4hep::mm << " mm\n";
          streamlog_out(MESSAGE) << "    Width of sensitive area: " << ws/dd4hep::mm << " mm\n";
          streamlog_out(MESSAGE) << "    Number of ladders: " << nModules << "\n";
          streamlog_out(MESSAGE) << "    Total sensitive area: " << area/dd4hep::cm2 << " cm^2\n";
          //streamlog_out(MESSAGE) << "    Sensors per ladder: " << layer.sensorsPerLadder << "\n";
          ilay++;
      }
    }
    catch ( std::exception &e) {
      streamlog_out(DEBUG) << "Caught exception " << e.what() << std::endl;
      try {
            streamlog_out(DEBUG) << "Trying ZDiskPetalsData.\n";
            dd4hep::rec::ZDiskPetalsData *layering = element.extension<dd4hep::rec::ZDiskPetalsData>() ;
            // We assume that layer numbering always starts from 0.
            int ilay = 0;
            for ( auto layer : layering->layers ) {
                streamlog_out(MESSAGE) << "  Layer " << ilay+1 << ":\n";
                double ls = layer.lengthSensitive;
                double wsi = layer.widthInnerSensitive;
                double wso = layer.widthOuterSensitive;
                int nModules = layer.petalNumber;
                double area = ls*(wsi+wso)*nModules/2;
                (*hitCounters.at(element.id()))[ilay] = new LayerHitCounter(area);
                streamlog_out(MESSAGE) << "    Length of sensitive area: " << ls/dd4hep::mm << " mm\n";
                streamlog_out(MESSAGE) << "    Inner width of sensitive area: " << wsi/dd4hep::mm << " mm\n";
                streamlog_out(MESSAGE) << "    Outer width of sensitive area: " << wso/dd4hep::mm << " mm\n";
                streamlog_out(MESSAGE) << "    Number of petals: " << nModules << "\n";
                streamlog_out(MESSAGE) << "    Total sensitive area: " << area/dd4hep::cm2 << " cm^2\n";
                //streamlog_out(MESSAGE) << "    Sensors per petal: " << layer.sensorsPerPetal << "\n";
                ilay++;
            }
          }
      catch ( std::exception &e1) {
        streamlog_out(DEBUG) << "Caught exception " << e1.what() << std::endl;
        streamlog_out(MESSAGE) << "  No layering extension in the "
                "detector element \'" << element.name() << "\'.\nTotal hits will be counted.\n";
        (*hitCounters.at(element.id()))[0] = new LayerHitCounter(-1.);
      }
    }


    streamlog_out(MESSAGE) << "\n    Added " << hitCounters.at(element.id())->size()
            << " hit counters for subsystem \'" << element.name() << "\'.\n\n";
  }
}


// Performed at the end of each run (file)
void TrackerHitCounter::processRunHeader( LCRunHeader*) {
  _nRun++;
  streamlog_out(MESSAGE) << "Processing run " << _nRun << "\n";

  for (auto sysCtr : hitCounters) {
      for (auto layerCtr : *sysCtr.second) {
          layerCtr.second->MarkRun();
      }
  }
}


void TrackerHitCounter::processEvent( LCEvent * evt) {

  _nEvt++;

  for (auto collname : m_trkHitCollNames) {

    streamlog_out(DEBUG) << "Looking into collection " << collname << "\n";
    LCCollectionVec *col = dynamic_cast<LCCollectionVec*>(evt->getCollection(collname));
    CellIDDecoder<SimTrackerHit> decoder(col);

    for (auto lcobj : (*col)) {

      SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(lcobj);

      if(!hit) {
        streamlog_out(WARNING) << "Collection " << collname << " contains objects "
            "of type other than SimTrackerHit!\nSkipping collection.\n";
        break;
      }

      int nsys = decoder(hit)["system"];
      streamlog_out(DEBUG) << "Found hit belonging to system #" << nsys << "\n";
      HitCtrMapIter hcmit = hitCounters.find(nsys);

      if (hcmit == hitCounters.end()) {
        streamlog_out(WARNING) << "Hit belongs to a system that is not analysed.\n";
      }
      else {
        SystemHitCounter *shc = hcmit->second;
        if (shc->size() == 1) {
          shc->at(0)->addHit();
        }
        else {
          int nLayer = decoder(hit)["layer"];
          HitCtrLayerIter hclit = shc->find(nLayer);
          if (hclit == shc->end()) {
            streamlog_out(ERROR) << "Hit in system ID=" << hcmit->first
                 << " belongs to layer number " << nLayer << ". Out of range!\n"
                    "  This should only happen if the xml detector description\n"
                    "  is different than the one used in the simulation.\n";
          }
          else {
            hclit->second->addHit();
          }
        }
      }

    } // Loop over hits in collection

  } // Loop over collections

}


void TrackerHitCounter::check( LCEvent *  ) { }


void TrackerHitCounter::end() {

  streamlog_out(MESSAGE) << "******************************************************\n";
  streamlog_out(MESSAGE) << "REPORT: \n";

  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();
  const std::vector< dd4hep::DetElement > &detElements = theDetector.detectors("tracker", true);

  for(dd4hep::DetElement element : detElements) {

    streamlog_out(MESSAGE) << "-----------------------------------------\n"
        "Subsystem : " << element.name() << "\n";
    auto hitctr = hitCounters.at(element.id());


    if (hitctr->size() == 1) {

      auto layerctr = hitctr->at(0);
      streamlog_out(MESSAGE) << "  Total: " << layerctr->getNHits() << " hits.\n";
      streamlog_out(MESSAGE) << "    (" << layerctr->getHitsPerRun() << " +- "
              << layerctr->getStDevHitsPerRun() << ") hits/run.\n";

      if (layerctr->isAreaAvailable()) {
        streamlog_out(MESSAGE) << "    " << layerctr->getHitsPerCm2() << " hits/cm^2.\n";
        streamlog_out(MESSAGE) << "    " << layerctr->getHitsPerCm2()/_nEvt << " hits/cm^2/event.\n";
        streamlog_out(MESSAGE) << "    (" << layerctr->getHitsPerCm2PerRun() << " +- "
                << layerctr->getStDevHitsPerCm2PerRun() << ") hits/cm^2/run.\n";
      }
      else {
        streamlog_out(MESSAGE) << "    " << layerctr->getNHits()/_nEvt << " hits/event.\n";
      }
    }
    else {

      for (auto layerpair : *hitctr ) {

        auto layerctr = layerpair.second;
        streamlog_out(MESSAGE) << "  Layer " << layerpair.first+1 << ": "
            << layerctr->getNHits() << " hits.\n";
        streamlog_out(MESSAGE) << "    (" << layerctr->getHitsPerRun() << " +- "
                  << layerctr->getStDevHitsPerRun() << ") hits/run.\n";

        if (layerctr->isAreaAvailable()) {
          streamlog_out(MESSAGE) << "    " << layerctr->getHitsPerCm2() << " hits/cm^2.\n";
          streamlog_out(MESSAGE) << "    " << layerctr->getHitsPerCm2()/_nEvt << " hits/cm^2/event.\n";
          streamlog_out(MESSAGE) << "    (" << layerctr->getHitsPerCm2PerRun() << " +- "
                  << layerctr->getStDevHitsPerCm2PerRun() << ") hits/cm^2/run.\n";
        }
        else {
          streamlog_out(MESSAGE) << "    " << layerctr->getNHits()/_nEvt << " hits/event.\n";
        }
      }
    }

    streamlog_out(MESSAGE) << "\n";
  }
  streamlog_out(MESSAGE) << "Analysed a total of " << _nEvt << " events in " << _nRun << " runs.\n";
  streamlog_out(MESSAGE) << "******************************************************\n";


  // Clean up
  for (auto systempair : hitCounters) {
        for (auto layerpair : *(systempair.second) ) {
            delete layerpair.second;
        }
        delete systempair.second;
    }
    hitCounters.clear();

}

