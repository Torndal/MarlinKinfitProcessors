#ifndef ZHHKinFit_h
#define ZHHKinFit_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <TFile.h>
#include <TTree.h>
#include "TopEventILC.h"

using namespace lcio ;
using namespace marlin ;


/**  An example processor for a kinematic fit
 *   
 *   ... testing a ZHH -> 2leptons + 4jets hypothesis
 *   with energy and momentum conservation
 *   and an equal mass constraint
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Needs 4 reconstructed jets
 *
 *  <h4>Output</h4> 
 *  A histogram.
 */

class ZHHKinFit : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new ZHHKinFit ; }

  ZHHKinFit(const ZHHKinFit&) = delete ;
  ZHHKinFit& operator=(const ZHHKinFit&) = delete ;  
  
  ZHHKinFit() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  double JetEnergyResolution(double E);

 protected:

  /** Input collection name.
   */
  std::string _jetcolName{}, _colLeptons, _name{},_OutputCollection{} ;
  // Output collection name.
  std::string _OutCol;
  /** Input parameter: center of mass energy.
   */
  float _ecm{}, _isrpzmax{};
  int _fitISR{}, _ifitter{}, _ievttrace{};
  bool _traceall{};

  double b{}, ISRPzMaxB{};

 
  float prob{}, chi2{}, bestchi2{}, bestprob{}, bestnit{}, bestmass1{}, bestmass2{}, beststartmass1{}, beststartmass2{}, bestphotonenergy{}, startmass1{}, startmass2{}, variable{};
  float momentum[3]{}, energy{};
  TLorentsVector bestfitlep1{}, bestfitlep2{}, bestfitjet1{}, bestfit2jet{}, bestfit3jet{}, bestfit4jet{};

  int _nRun{}, _nEvt{}, nit{};
 
  int bestperm{}, errorflag{};
  
  TopEventILC* topevent{};
 

  //output
  // TTree *outTree;
         
} ;

#endif
