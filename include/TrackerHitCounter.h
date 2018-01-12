#ifndef AreaReader_h
#define AreaReader_h 1

#include "marlin/Processor.h"
#include <string>

/** AreaReader: Reports basic geometrical information on tracker detector
 *              elements.
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

    double getHitsPerPerCm2() const { return double(_nHits)/_area; }
    void addHits(int newHits) { _nHits += newHits; }

private:
    const double _area;
    int _nHits;
};

typedef std::map<int, LayerHitCounter*> SystemHitCounter;


class TrackerHitCounter : public marlin::Processor {

public:

  virtual Processor*  newProcessor() { return new TrackerHitCounter ; }


  TrackerHitCounter() ;

  /** Read and report the areas of tracker elements.
   */
  virtual void init() ;

  /// do nothing
  virtual void processRunHeader( LCRunHeader* run ) ;

  /// Count hits
  virtual void processEvent( LCEvent * evt ) ; 

  /// do nothing
  virtual void check( LCEvent * evt ) ; 

  /// report hits and clean up.
  virtual void end() ;


protected:

  std::map<int, SystemHitCounter*> hitCounters;

} ;

#endif



