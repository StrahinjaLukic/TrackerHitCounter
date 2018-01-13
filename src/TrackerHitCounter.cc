#include "DD4hep/Detector.h"
#include <DD4hep/DetType.h>
#include <TrackerHitCounter.h>
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"
#include <IMPL/LCCollectionVec.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <vector>

using namespace marlin ;

// Todo: update to v00-05:
// DD4hep -> dd4hep
// DD4hep::DDRec -> dd4hep::rec
// DD4hep::Geometry::LCDD -> dd4hep::Detector

TrackerHitCounter aTrackerHitCounter ;


TrackerHitCounter::TrackerHitCounter() : Processor("TrackerHitCounter"),
    hitCounters(), m_trkHitCollNames()
{

  // modify processor description
  _description = "TrackerHitCounter counts hits in tracker detector elements "
          "and reports number of hits per unit area." ;

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
          streamlog_out(MESSAGE) << "  Layer " << ilay << ":\n";
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
                streamlog_out(MESSAGE) << "  Layer " << ilay << ":\n";
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
                "detector element \'" << element.name() << "\'. Total hits will be counted.\n";
        (*hitCounters.at(element.id()))[1] = new LayerHitCounter();
      }
    }


    streamlog_out(MESSAGE) << "\n    Added " << hitCounters.at(element.id())->size()
            << " hit counters for ID=" << element.id() << ".\n\n";
  }
}


void TrackerHitCounter::processRunHeader( LCRunHeader*) { }


void TrackerHitCounter::processEvent( LCEvent * evt) {

  for (auto collname : m_trkHitCollNames) {

    streamlog_out(MESSAGE) << "Looking into collection " << collname << "\n";
    LCCollectionVec *col = dynamic_cast<LCCollectionVec*>(evt->getCollection(collname));
    CellIDDecoder<SimTrackerHit> decoder(col);

    for (auto lcobj : (*col)) {

      SimTrackerHitImpl* hit = dynamic_cast<SimTrackerHitImpl*>(lcobj);

      if(!hit) {
        streamlog_out(WARNING) << "Collection " << collname << " does not contain "
            "SimTrackerHits! Skipping collection.\n";
        break;
      }

      int nsys = decoder(hit)["system"];
      streamlog_out(MESSAGE) << "Found hit decoded to system #" << nsys << "\n";
      HitCtrMapIter hcmit = hitCounters.find(nsys);

      if (hcmit == hitCounters.end()) {
        streamlog_out(WARNING) << "Hit belongs to system that is not analysed.\n";
      }
      else {
        int nLayer = decoder(hit)["layer"];
        HitCtrLayerIter hclit = hcmit->second->find(nLayer);
        if (hclit == hcmit->second->end()) {
          streamlog_out(ERROR) << "Hit in system ID=" << hcmit->first
               << " belongs to layer number " << nLayer << ". Out of range!\n"
                  "  This should only happen if the xml detector description\n"
                  "  is different than the one used in the simulation.\n";
        }
        else {
          hclit->second->addHit();
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

    for (auto layerpair : *hitctr ) {

      int nhits = layerpair.second->getNHits();
      double hitspercm2 = layerpair.second->getHitsPerCm2();
      streamlog_out(MESSAGE) << "  Layer " << layerpair.first << ": " << nhits << " hits. ("
          << hitspercm2 << " hits/cm^2).\n";
    }

    streamlog_out(MESSAGE) << "\n";
  }
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

