#include "DD4hep/Detector.h"
#include <DD4hep/DetType.h>
#include <TrackerHitCounter.h>
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"
#include <vector>

using namespace marlin ;

// Todo: update to v00-05:
// DD4hep -> dd4hep
// DD4hep::DDRec -> dd4hep::rec
// DD4hep::Geometry::LCDD -> dd4hep::Detector

TrackerHitCounter aTrackerHitCounter ;


TrackerHitCounter::TrackerHitCounter() : Processor("TrackerHitCounter"),
    hitCounters()
{

  // modify processor description
  _description = "TrackerHitCounter counts hits in tracker detector elements "
          "and reports number of hits per unit area." ;

}

void TrackerHitCounter::init() { 

  // usually a good idea to
  // printParameters() ;

  // v01-19-05 namespace
  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();

  const std::vector< dd4hep::DetElement > &detElements = theDetector.detectors("tracker", true);

  for(dd4hep::DetElement element : detElements) {

    streamlog_out(MESSAGE) << "*********************************************************************\n";
          streamlog_out(MESSAGE) << "Detector element \'" << element.name() << "\' of type \'"
        << element.type() << "\':\n";

    hitCounters[element.volumeID()] = new SystemHitCounter;

    try {
      streamlog_out(DEBUG) << "Trying ZPlanarData.\n";
      dd4hep::rec::ZPlanarData* layering = element.extension<dd4hep::rec::ZPlanarData>() ;

      int ilay = 1;
      for ( auto layer : layering->layers ) {
          streamlog_out(MESSAGE) << "  Layer " << ilay << ":\n";
          double ls = layer.zHalfSensitive*2;
          double ws = layer.widthSensitive;
          int nModules = layer.ladderNumber;
          double area = ls*ws*nModules;
          (*hitCounters.at(element.volumeID()))[ilay] = new LayerHitCounter(area);
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
            int ilay = 1;
            for ( auto layer : layering->layers ) {
                streamlog_out(MESSAGE) << "  Layer " << ilay << ":\n";
                double ls = layer.lengthSensitive;
                double wsi = layer.widthInnerSensitive;
                double wso = layer.widthOuterSensitive;
                int nModules = layer.petalNumber;
                double area = ls*(wsi+wso)*nModules/2;
                (*hitCounters.at(element.volumeID()))[ilay] = new LayerHitCounter(area);
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
        (*hitCounters.at(element.volumeID()))[1] = new LayerHitCounter();
      }
    }


    streamlog_out(MESSAGE) << "\n    Added " << hitCounters.at(element.volumeID())->size()
            << " hit counters for ID=" << element.volumeID() << ".\n\n";
  }
}


void TrackerHitCounter::processRunHeader( LCRunHeader*) { } 

void TrackerHitCounter::processEvent( LCEvent * ) {



}

void TrackerHitCounter::check( LCEvent *  ) { }

void TrackerHitCounter::end() {

    for (auto systempair : hitCounters) {
        for (auto modulepair : *(systempair.second) ) {
            delete modulepair.second;
        }
        delete systempair.second;
    }
    hitCounters.clear();

}

