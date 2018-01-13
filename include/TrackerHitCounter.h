#ifndef TrackerHitCounter_h
#define TrackerHitCounter_h 1

#include "marlin/Processor.h"
#include "DD4hep/DD4hepUnits.h"
#include <string>

/** TrackerHitCounter: Marlin processor that counts SimTrackerHits in tracker detector
 *     elements and reports the number of hits per unit area.
 *
 *  S. Lukić, Jan 2018
 */

class LayerHitCounter {
public:
    LayerHitCounter() :
        _area(1.), _nHits(0)
    {}

    LayerHitCounter(double area) :
        _area(area), _nHits(0)
    {}

    virtual ~LayerHitCounter() {}

    int getNHits() const { return _nHits; }
    double getHitsPerCm2() const { return double(_nHits)/(_area/dd4hep::cm2); }
    bool isAreaAvailable() const { return (_area > 0.); }
    void addHit() { _nHits++; }
    void addHits(int newHits) { _nHits += newHits; }


private:
    const double _area;
    int _nHits;
};

typedef std::map<int, LayerHitCounter*> SystemHitCounter;
typedef SystemHitCounter::iterator HitCtrLayerIter;
typedef std::map<int, SystemHitCounter*> HitCtrMap;
typedef HitCtrMap::iterator HitCtrMapIter;

class TrackerHitCounter : public marlin::Processor {

public:

  virtual Processor*  newProcessor() { return new TrackerHitCounter ; }


  TrackerHitCounter() ;

  /** Read and report the areas of tracker elements.
   *  Construct the corresponding hit counters.
   */
  virtual void init() ;

  /// Count runs
  virtual void processRunHeader( LCRunHeader* run ) ;

  /// Count hits in subsystems
  virtual void processEvent( LCEvent * evt ) ;

  /// do nothing
  virtual void check( LCEvent * evt ) ;

  /// report hits and clean up.
  virtual void end() ;


protected:

  HitCtrMap hitCounters{} ;

  StringVec m_trkHitCollNames{} ;

  int _nRun{} ;
  int _nEvt{} ;
} ;


#endif // TrackerHitCounter_h



