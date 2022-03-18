#include "ZHHKinFit.h"
#include <iostream>
#include <vector>
#include <string>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/ITupleFactory.h>
#include <AIDA/ITuple.h>
#endif

#include "UTIL/LCRelationNavigator.h"
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/ReconstructedParticle.h>
//#include <root/TLorentzVector.h>
#include <CLHEP/Vector/LorentzVector.h>
#include "JetFitObject.h"
#include "LeptonFitObject.h"
#include "ISRPhotonFitObject.h"
#include "MomentumConstraint.h"
#include "OPALFitterGSL.h"
#include "NewFitterGSL.h"
#include "TextTracer.h"
#include "NewtonFitterGSL.h"
#include "FourJetPairing.h"
#include "MassConstraint.h"

using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace CLHEP ;


ZHHKinFit aZHHKinFit ;

// function to define the jet energy resolution (in GeV)
double ZHHKinFit::JetEnergyResolution(double E)
{
  // examples here derived by Benjamin Hermberg from e+e- -> udsc:
  
  // 1) comparing jet-level to quark-level energies 
  //    (using MarlinReco/Analysis/RecoMCTruthLink/QuarkJetPairing.cc)
  // double result = std::sqrt(pow(0.6908,2)*(E)+(pow(0.02596,2)*pow(E,2))); 
  
  // 2) 120%/sqrt(E), gives best convergence of 5C fit on e+e- -> udsc
  double result = 1.2*std::sqrt(E);  
  return result;      
}

ZHHKinFit::ZHHKinFit() : Processor("ZHHKinFit") {
  
  // modify processor description
  _description = "ZHHKinFit does a kinematic fit on 2 lepton + 4 jet events" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "JetCollectionName" , 
			   "Name of the Jet collection"  ,
			   _jetcolName ,
			   std::string("Durham4Jets") ) ;
                           
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "InputIsoLeptonsCollection",
			   "Name of collection with the selected isolated lepton",
			   _colLeptons,
			   std::string("ISOLeptons") );

  registerProcessorParameter( "ECM" ,
                              "Center-of-Mass Energy in GeV",
                              _ecm,
                              (float)500.);
                              
  registerProcessorParameter( "FitISR" ,
                              "0: Fit hypothesis without ISR   1: Fit hypothesis including ISR",
                              _fitISR,
                              (int) 1);
                              
  registerProcessorParameter( "ISRPzMax" ,
                              "Maximum possible energy for a single ISR photon",
                              _isrpzmax,
                              (float)125.6);

  registerProcessorParameter( "SigmaEnergyScaleFactor" ,
			      "Factor for scaling up energy error",
			      m_SigmaEnergyScaleFactor,
			      float(1.0f));

  registerProcessorParameter( "SigmaAnglesScaleFactor" ,
                              "Factor for scaling up angular errors",
                              m_SigmaAnglesScaleFactor,
                              float(1.0f));

  registerProcessorParameter( "fitter" ,
                              "0 = OPALFitter, 1 = NewFitter, 2 = NewtonFitter",
                              _ifitter,
                              (int)0);
                              
  registerProcessorParameter( "traceall" ,
                              "set true if every event should be traced",
                              _traceall,
                              (bool)false);
                              
  registerProcessorParameter( "ievttrace" ,
                              "number of individual event to be traced",
                              _ievttrace,
                              (int)0);
   
  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "FitOutputCollection",
			    "Output Fit Collection" ,
			    _OutCol,
			    std::string("FitReco")) ;

  //registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
  //			    "PairedOutputCollection",
  //			    "Output Paired Reco Collection" ,
  //			    _OutPairedCol,
  //			    std::string("PairedReco")) ;
}


void ZHHKinFit::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  //DDMarlinCED::init(this);

  //m_Bfield = MarlinUtil::getBzAtOrigin();
  m_Bfield = 3.5 ;
  streamlog_out(DEBUG0) << " BField =  "<< m_Bfield << " Tesla" << std::endl ;
  c = 2.99792458e8;
  mm2m = 1e-3;
  eV2GeV = 1e-9;
  eB = m_Bfield * c * mm2m * eV2GeV;

  b = (double) 0.00464564*( std::log(_ecm*_ecm*3814714.)-1. );
  //= 2*alpha/pi*( ln(s/m_e^2)-1 )
  ISRPzMaxB = std::pow((double)_isrpzmax,b);
  
}

void ZHHKinFit::processRunHeader( LCRunHeader* ) { 
  _nRun++ ;
} 

void ZHHKinFit::processEvent( LCEvent * evt ) { 

    
  streamlog_out(DEBUG) 
    << " processing event " << evt->getEventNumber() 
    << "  in run "          << evt->getRunNumber() 
    << std::endl ;
  // this gets called for every event 
  // usually the working horse ...

#ifdef MARLIN_USE_AIDA
    
  // define a histogram pointer
  static AIDA::IHistogram1D* hRecHMassBest ;    
  static AIDA::IHistogram1D* hRecHMassAll ;    
  static AIDA::IHistogram1D* hRecHMassNoFitBest ;    
  static AIDA::IHistogram1D* hRecHMassNoFitAll ;    
  static AIDA::IHistogram1D* hTestHMassNoFitAll ;    
  static AIDA::IHistogram1D* hFitProbBest ;    
  static AIDA::IHistogram1D* hFitProbAll ;    
  static AIDA::IHistogram1D* hNItBest ;    
  static AIDA::IHistogram1D* hNItAll ;    
  static AIDA::IHistogram1D* hPhotonEnergy ;    
  static AIDA::IHistogram1D* hJetMass ;    
  static AIDA::IHistogram1D* hFitError;    
             
  streamlog_out(DEBUG) 
    << " processing event " << evt->getEventNumber() 
    << "  in run "          << evt->getRunNumber() 
    << std::endl ;
  
  if( isFirstEvent() ) { 
    
    hRecHMassBest = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecHMassBest", "M_W", 200, 0., 200. ) ; 
    hRecHMassAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecHMassAll", "M_W", 200, 0., 200. ) ; 
    hRecHMassNoFitBest = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecHMassNoFitBest", "M_W", 200, 0., 200. ) ; 
    hRecHMassNoFitAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecHMassNoFitAll", "M_W", 200, 0., 200. ) ; 
    hTestHMassNoFitAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hTestHMassNoFitAll", "M_W", 200, 0., 200. ) ; 
    hFitProbBest = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hFitProb", "fit probability", 100, 0., 1. ) ; 
    hFitProbAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hFitProbAll", "fit probability", 100, 0., 1. ) ; 
    hNItBest = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hNItBest", "number of iterations", 200, 0., 200. ) ; 
    hNItAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hNItAll", "number of iterations", 200, 0., 200. ) ; 
    hPhotonEnergy = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPhotonEnergy", "ISR photon energy", 200, 0., 400. ) ; 
    hJetMass = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hJetMass", "Jet Mass", 200, 0., 100. ) ; 
    if (_ifitter == 1) {
      hFitError = 
        AIDAProcessor::histogramFactory(this)->
        createHistogram1D( "hFitError", "Error flag", 100, -0.5, 99.5 ) ; 
    }
    else {    
      hFitError = 
        AIDAProcessor::histogramFactory(this)->
        createHistogram1D( "hFitError", "Error flag", 10, -0.5, 9.5 ) ; 
    }    

  }

#endif   
   
  streamlog_out(DEBUG) 
    << " processing event " << evt->getEventNumber() 
    << "  in run "          << evt->getRunNumber() 
    << std::endl ;
  
  
  HepLorentzVector lvec;
  
    
  // fill histogram from LCIO data :

  //////////////////   JETS ///////////////////////////
   
  LCCollection* jetcol = evt->getCollection( _jetcolName ) ;
  LCCollection* lepcol = evt->getCollection( _colLeptons ) ;
  if (jetcol != 0 && lepcol != 0) {
  
    int nJETS = jetcol->getNumberOfElements()  ;
    int nLEPS = lepcol->getNumberOfElements()  ;
    streamlog_out(DEBUG) 
      << " found " << nJETS
      << " jets in event " << evt->getEventNumber() 
      << "  in run "          << evt->getRunNumber() 
      << std::endl ;
                      
    if (nJETS != 4) return;               
    if (nLEPS != 2) return;               
                   
    float yminus = jetcol ->parameters().getFloatVal( "YMinus");              
    streamlog_out(DEBUG)  << " yminus = " << yminus << std::endl ;
                     
    float yplus = jetcol ->parameters().getFloatVal( "YPlus");              
    streamlog_out(DEBUG)  << " yplus = " << yplus << std::endl ;
                            
    // original fit objects - save for next permutation
    JetFitObject* j1 = 0;
    JetFitObject* j2 = 0;
    JetFitObject* j3 = 0;
    JetFitObject* j4 = 0;
    //orginal lepton fit objects
    LeptonFitObject* l1 = 0; 
    LeptonFitObject* l2 = 0;
        
    // angular resolutions - optimised for best convergence of 5C fit on e+e- -> udsc
    //double errtheta = 0.1;   //   100mrad
    //double errphi = 0.1;     //   100mrad
       
    for(int i=0; i< nLEPS ; i++){
      float parameters[ 3 ]{ 0.0 } , errors[ 3 ]{ 0.0 };
      ReconstructedParticle* l = dynamic_cast<ReconstructedParticle*>( lepcol->getElementAt( i ) ) ;
      if (l) {
	streamlog_out(DEBUG)
          << " found lepton in event " << evt->getEventNumber()
          << "  in run "          << evt->getRunNumber()
          << std::endl ;
	//lvec = HepLorentzVector ((l->getMomentum())[0],(l->getMomentum())[1],(l->getMomentum())[2],l->getEnergy());
	//double invPt = 1./sqrt(lvec.px()*lvec.px()+lvec.py()*lvec.py());

	getLeptonParameters( l , parameters , errors );
        if (i == 0) {
          //l1 = new LeptonFitObject (invPt, lvec.theta(), lvec.phi(),
          //                       JetEnergyResolution(lvec.e()), errtheta, errphi, lvec.m());
	  l1 = new LeptonFitObject (parameters[ 0 ] , parameters[ 1 ] , parameters[ 2 ],
				    errors[ 0 ] , errors[ 1 ] , errors[ 2 ] , l->getMass() );
          l1->setName("Lepton1");
          streamlog_out(DEBUG)  << " start four-vector of first  lepton: " << *l1  << std::endl ;
        }
        else if (i == 1) {
	  //l2 = new LeptonFitObject (invPt, lvec.theta(), lvec.phi(),
          //                       JetEnergyResolution(lvec.e()), errtheta, errphi, lvec.m());
          l2 = new LeptonFitObject (parameters[ 0 ] , parameters[ 1 ] , parameters[ 2 ],
				    errors[ 0 ] , errors[ 1 ] , errors[ 2 ] , l->getMass() );
          l2->setName("Lepton2");
          streamlog_out(DEBUG) << " start four-vector of second  lepton: " << *l2  << std::endl ;
        }
      }
    }
   
    ReconstructedParticle* jrps[4];

    for(int i=0; i< nJETS ; i++){
      ReconstructedParticle* j = dynamic_cast<ReconstructedParticle*>( jetcol->getElementAt( i ) ) ;      
      if (j) {
	jrps[i] = j;
	streamlog_out(DEBUG) 
	  << " found jet in event " << evt->getEventNumber() 
	  << "  in run "          << evt->getRunNumber() 
	  << std::endl ;
	//lvec = HepLorentzVector ((j->getMomentum())[0],(j->getMomentum())[1],(j->getMomentum())[2],j->getEnergy()); 
	double sigmaE , sigmaTheta , sigmaPhi;
	TLorentzVector jetFourMomentum( j->getMomentum(), j->getEnergy() );
	std::vector< float > jetCovMat( 10 , 0.0 );
	for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element ) {
	  jetCovMat[ i_Element ] = j->getCovMatrix()[ i_Element ];
	}
	getJetResolutions( jetFourMomentum , jetCovMat , sigmaE , sigmaTheta , sigmaPhi );
	if (i == 0) { 
	  j1 = new JetFitObject( jetFourMomentum.E() , jetFourMomentum.Theta() , jetFourMomentum.Phi(),
				 m_SigmaEnergyScaleFactor * sigmaE , m_SigmaAnglesScaleFactor * sigmaTheta , m_SigmaAnglesScaleFactor * sigmaPhi , jetFourMomentum.M() );
	  j1->setName("Jet1");
	  streamlog_out(DEBUG)  << " start four-vector of first  jet: " << *j1  << std::endl ;
	}
	else if (i == 1) { 
	  j2 = new JetFitObject( jetFourMomentum.E() , jetFourMomentum.Theta() , jetFourMomentum.Phi(),
                                 m_SigmaEnergyScaleFactor * sigmaE , m_SigmaAnglesScaleFactor * sigmaTheta , m_SigmaAnglesScaleFactor * sigmaPhi , jetFourMomentum.M() );
	  j2->setName("Jet2");
	  streamlog_out(DEBUG) << " start four-vector of second  jet: " << *j2  << std::endl ;
	}
	else if (i == 2) {
	  j3 = new JetFitObject( jetFourMomentum.E() , jetFourMomentum.Theta() , jetFourMomentum.Phi(),
                                 m_SigmaEnergyScaleFactor * sigmaE , m_SigmaAnglesScaleFactor * sigmaTheta , m_SigmaAnglesScaleFactor * sigmaPhi , jetFourMomentum.M() );
	  j3->setName("Jet3");
	  streamlog_out(DEBUG) << " start four-vector of third  jet: " << *j3  << std::endl ;
	}
	else if (i == 3) { 
	  j4 = new JetFitObject( jetFourMomentum.E() , jetFourMomentum.Theta() , jetFourMomentum.Phi(),
                                 m_SigmaEnergyScaleFactor * sigmaE , m_SigmaAnglesScaleFactor * sigmaTheta , m_SigmaAnglesScaleFactor * sigmaPhi , jetFourMomentum.M() );
	  j4->setName("Jet4");
	  streamlog_out(DEBUG) << " start four-vector of fourth  jet: " << *j4  << std::endl ;
	}
#ifdef MARLIN_USE_AIDA
	hJetMass->fill(j->getMass());
#endif           
      }
    }
#ifdef MARLIN_USE_AIDA
    double en, px, py, pz, mass; 
    if (jrps[0] && jrps[1] && jrps[2] && jrps[3]) {
      en = jrps[0]->getEnergy()     +jrps[1]->getEnergy();
      px = jrps[0]->getMomentum()[0]+jrps[1]->getMomentum()[0];
      py = jrps[0]->getMomentum()[1]+jrps[1]->getMomentum()[1];
      pz = jrps[0]->getMomentum()[2]+jrps[1]->getMomentum()[2];
      mass = en*en-px*px-py*py-pz*pz;
      if (mass >= 0) hTestHMassNoFitAll->fill( std::sqrt(mass) ) ;
      en = jrps[2]->getEnergy()     +jrps[3]->getEnergy();
      px = jrps[2]->getMomentum()[0]+jrps[3]->getMomentum()[0];
      py = jrps[2]->getMomentum()[1]+jrps[3]->getMomentum()[1];
      pz = jrps[2]->getMomentum()[2]+jrps[3]->getMomentum()[2];
      mass = en*en-px*px-py*py-pz*pz;
      if (mass >= 0) hTestHMassNoFitAll->fill( std::sqrt(mass) ) ;
      en = jrps[0]->getEnergy()     +jrps[2]->getEnergy();
      px = jrps[0]->getMomentum()[0]+jrps[2]->getMomentum()[0];
      py = jrps[0]->getMomentum()[1]+jrps[2]->getMomentum()[1];
      pz = jrps[0]->getMomentum()[2]+jrps[2]->getMomentum()[2];
      mass = en*en-px*px-py*py-pz*pz;
      if (mass >= 0) hTestHMassNoFitAll->fill( std::sqrt(mass) ) ;
      en = jrps[1]->getEnergy()     +jrps[3]->getEnergy();
      px = jrps[1]->getMomentum()[0]+jrps[3]->getMomentum()[0];
      py = jrps[1]->getMomentum()[1]+jrps[3]->getMomentum()[1];
      pz = jrps[1]->getMomentum()[2]+jrps[3]->getMomentum()[2];
      mass = en*en-px*px-py*py-pz*pz;
      if (mass >= 0) hTestHMassNoFitAll->fill( std::sqrt(mass) ) ;
      en = jrps[0]->getEnergy()     +jrps[3]->getEnergy();
      px = jrps[0]->getMomentum()[0]+jrps[3]->getMomentum()[0];
      py = jrps[0]->getMomentum()[1]+jrps[3]->getMomentum()[1];
      pz = jrps[0]->getMomentum()[2]+jrps[3]->getMomentum()[2];
      mass = en*en-px*px-py*py-pz*pz;
      if (mass >= 0) hTestHMassNoFitAll->fill( std::sqrt(mass) ) ;
      en = jrps[1]->getEnergy()     +jrps[2]->getEnergy();
      px = jrps[1]->getMomentum()[0]+jrps[2]->getMomentum()[0];
      py = jrps[1]->getMomentum()[1]+jrps[2]->getMomentum()[1];
      pz = jrps[1]->getMomentum()[2]+jrps[2]->getMomentum()[2];
      mass = en*en-px*px-py*py-pz*pz;
      if (mass >= 0) hTestHMassNoFitAll->fill( std::sqrt(mass) ) ;
    }  
#endif           
       
    const int NJETS = 4;
    const int NLEPS = 2;
    streamlog_out(DEBUG)  << "*j1" << *j1  << "*j2" << *j2  << "*j3" << *j3  << "*j4" << *j4  << std::endl ;
    streamlog_out(DEBUG)  << "*l1" << *l1  << "*l2" << *l2  << std::endl ;
    
       // these get changed by the fit -> reset after each permutation!
       JetFitObject fitjets[NJETS] = {*j1, *j2, *j3, *j4};
       LeptonFitObject fitleps[NLEPS] = {*l1, *l2};
       for (int i = 0; i < NJETS; ++i)
         streamlog_out(DEBUG)  << "fitjets[ " << i << "]: " << fitjets[i]  << std::endl ;
       for (int i = 0; i < NLEPS; ++i)
         streamlog_out(DEBUG)  << "fitleps[ " << i << "]: " << fitleps[i]  << std::endl ;
 
       // these point allways to the fitjets array, which gets reset.
       JetFitObject *jets[NJETS];
       for (int i = 0; i < NJETS; ++i) jets[i] = &fitjets[i];
       for (int i = 0; i < NJETS; ++i)
         streamlog_out(DEBUG)  << "start four-vector of jets[ " << i << "]: " << *(jets[i])  << std::endl ;
       LeptonFitObject *leps[NLEPS];
       for (int i = 0; i < NLEPS; ++i) leps[i] = &fitleps[i];
       for (int i = 0; i < NLEPS; ++i)
         streamlog_out(DEBUG)  << "start four-vector of leps[ " << i << "]: " << *(leps[i])  << std::endl ;

       FourJetPairing pairing (jets);
       JetFitObject *permutedjets[NJETS];

       bestprob = 0.;
       bestnit = 0;
       bestmass1 = 0., bestmass2 = 0.;
       beststartmass1 = 0., beststartmass2 = 0.;
       startmass1 = 0., startmass2 = 0.;
       bestphotonenergy = 0.;

       for (int iperm = 0; iperm < pairing.getNPerm(); iperm++) {

         streamlog_out(DEBUG) 
                       << " ================================================= "  
                       << std::endl ;
         streamlog_out(DEBUG) 
                       << " iperm = " << iperm 
                       << std::endl ;

         // important: (re-)set fitjets array!
         fitjets[0] = *j1;
         fitjets[1] = *j2;
         fitjets[2] = *j3;
         fitjets[3] = *j4;

	 fitleps[0] = *l1;
	 fitleps[1] = *l2;


         pairing.nextPermutation (permutedjets);
         for (int i = 0; i < NJETS; ++i)
            streamlog_out(DEBUG)  << "start four-vector of jet " << i << ": " << *(permutedjets[i])  << std::endl ;
                              
         //MomentumConstraint pxc (1, 0);
         // crossing angle 14 mrad = 7/500
         MomentumConstraint pxc (0, 1, 0, 0, 7.0);
         pxc.setName("sum(p_x)");
         for (int i = 0; i < NJETS; ++i)
	   pxc.addToFOList (*(permutedjets[i]));
	 for (int i = 0; i < NLEPS; ++i)
	   pxc.addToFOList (*(leps[i]));

         MomentumConstraint pyc (0, 0, 1);
         pyc.setName("sum(p_y)");
         for (int i = 0; i < NJETS; ++i)
            pyc.addToFOList (*(permutedjets[i]));
	 for (int i = 0; i < NLEPS; ++i)
	   pyc.addToFOList (*(leps[i]));
        
         MomentumConstraint pzc (0, 0, 0, 1);
         pzc.setName("sum(p_z)");
         for (int i = 0; i < NJETS; ++i)
            pzc.addToFOList (*(permutedjets[i]));
	 for (int i = 0; i < NLEPS; ++i)
	   pzc.addToFOList (*(leps[i]));

         streamlog_out(DEBUG) << "ECM = " << _ecm  << std::endl ; 
	 MomentumConstraint ec(1, 0, 0, 0, _ecm);
         ec.setName("sum(E)");
         for (int i = 0; i < NJETS; ++i)
	   ec.addToFOList (*(permutedjets[i]));
	 for (int i = 0; i < NLEPS; ++i)
	   ec.addToFOList (*(leps[i]));

	 streamlog_out(DEBUG)  << "Value of pxc before fit: " << pxc.getValue() << std::endl ;
	 streamlog_out(DEBUG)  << "Value of pyc before fit: " << pyc.getValue() << std::endl ;
	 streamlog_out(DEBUG)  << "Value of pzc before fit: " << pzc.getValue() << std::endl ;
	 streamlog_out(DEBUG)  << "Value of ec before fit: " << ec.getValue() << std::endl ;
                      

         // ISR Photon initialized with missing p_z
         ISRPhotonFitObject *photon = new ISRPhotonFitObject (0., 0., -pzc.getValue(), b, ISRPzMaxB);
        
	 if(_fitISR){
	   streamlog_out(DEBUG)  << "start four-vector of ISR photon: " << *(photon) << std::endl ;
                      
	   pxc.addToFOList (*(photon));
	   pyc.addToFOList (*(photon));
	   pzc.addToFOList (*(photon));
	   ec.addToFOList  (*(photon));
         }

         MassConstraint w(0.);
         w.addToFOList (*(permutedjets[0]), 1);
         w.addToFOList (*(permutedjets[1]), 1);
         w.addToFOList (*(permutedjets[2]), 2);
         w.addToFOList (*(permutedjets[3]), 2);

         startmass1 = w.getMass(1);
         startmass2 = w.getMass(2);
         
	 streamlog_out(DEBUG) << "start mass of H 1: " << startmass1 << std::endl ;
	 streamlog_out(DEBUG) << "start mass of H 2: " << startmass2 << std::endl ;
                      
#ifdef MARLIN_USE_AIDA
         hRecHMassNoFitAll->fill( 0.5*(startmass1 + startmass2) ) ;
         //hRecHMassNoFitAll->fill( startmass2 ) ;
#endif
         BaseFitter *pfitter;
         if (_ifitter == 1) {
           pfitter = new NewFitterGSL();
           if (evt->getEventNumber()== _ievttrace || _traceall) (dynamic_cast<NewFitterGSL*>(pfitter))->setDebug (1);
         }
         else if (_ifitter == 2) {
           pfitter = new NewtonFitterGSL();
           if (evt->getEventNumber()== _ievttrace || _traceall) (dynamic_cast<NewtonFitterGSL*>(pfitter))->setDebug (1);
         }
         else {
           // OPALFitter has no method setDebug !
           pfitter = new OPALFitterGSL();
         }
         BaseFitter &fitter = *pfitter;
  
         TextTracer tracer (std::cout);
         if (evt->getEventNumber()== _ievttrace || _traceall) fitter.setTracer (tracer);
         
         for (int i = 0; i < NLEPS; ++i)
	   fitter.addFitObject (*(leps[i]));
         for (int i = 0; i < NJETS; ++i)
	   fitter.addFitObject (*(permutedjets[i]));
         if(_fitISR){
	   fitter.addFitObject (*(photon));
         }
         fitter.addConstraint (pxc);
         fitter.addConstraint (pyc);
         fitter.addConstraint (pzc);
         fitter.addConstraint (ec);
         fitter.addConstraint (w);

         prob = fitter.fit();
         chi2 = fitter.getChi2();
         nit = fitter.getIterations();

         streamlog_out(DEBUG) << "fit probability = " << prob << std::endl ;  
         streamlog_out(DEBUG) << "fit chi2 = " << chi2  << std::endl ; 
         streamlog_out(DEBUG) << "error code: " << fitter.getError() << std::endl ;
                      
         for (int i = 0; i < NJETS; ++i) {
	   streamlog_out(DEBUG)  << "final four-vector of jet " << i << ": " << *(permutedjets[i]) << std::endl ;
	 }
         if(_fitISR){
	   streamlog_out(DEBUG)  << "final four-vector of ISR photon: " << *(photon) << std::endl ;
	 }

         streamlog_out(DEBUG)  << "final mass of H 1: " << w.getMass(1) << std::endl ;
	 streamlog_out(DEBUG)  << "final mass of H 2: " << w.getMass(2) << std::endl ;
                      
         hFitError->fill( fitter.getError() ) ;
         if (fitter.getError() == 0) {
#ifdef MARLIN_USE_AIDA
           hFitProbAll->fill( prob ) ;
           hNItAll->fill( nit ) ;
           hRecHMassAll->fill( 0.5*(w.getMass(1)+w.getMass(2)) ) ;
           //hRecHMassAll->fill( w.getMass(2)) ;
#endif
	   //TLorentzVector recoJetFourMomentum( recoJet->getMomentum()[ 0 ] , recoJet->getMomentum()[ 1 ] , recoJet->getMomentum()[ 2 ] , recoJet->getEnergy() );
	   //           if (prob > bestprob && w.getMass(1) > 50 && w.getMass(1) < 110) {
           if (prob > bestprob) {
	     bestchi2 = chi2;
             bestprob = prob;
             bestnit  = nit;
             bestmass1 = w.getMass(1);
             bestmass2 = w.getMass(2);
             beststartmass1 = startmass1;
             beststartmass2 = startmass2;
             bestphotonenergy = photon->getE();
	     bestfitlep1.SetPxPyPzE(leps[0]->getPx(),leps[0]->getPy(),leps[0]->getPz(),leps[0]->getE());
	     bestfitlep2.SetPxPyPzE(leps[1]->getPx(),leps[1]->getPy(),leps[1]->getPz(),leps[1]->getE());
	     bestfitjet1.SetPxPyPzE(permutedjets[0]->getPx(),permutedjets[0]->getPy(),permutedjets[0]->getPz(),permutedjets[0]->getE());
	     bestfitjet2.SetPxPyPzE(permutedjets[1]->getPx(),permutedjets[1]->getPy(),permutedjets[1]->getPz(),permutedjets[1]->getE());
	     bestfitjet3.SetPxPyPzE(permutedjets[2]->getPx(),permutedjets[2]->getPy(),permutedjets[2]->getPz(),permutedjets[2]->getE());
	     bestfitjet4.SetPxPyPzE(permutedjets[3]->getPx(),permutedjets[3]->getPy(),permutedjets[3]->getPz(),permutedjets[3]->getE());

           }
         }
         else {
	   streamlog_out(DEBUG) << "FIT ERROR = " << fitter.getError() << ", not filling histograms!"  << std::endl ;
	 }
         delete photon;
         streamlog_out(DEBUG) << "end permutation " << std::endl ;
       }

       //4-momentum of KinFit objects after fit for best combination found in the fit
       LCCollectionVec *OutputCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
       /*ReconstructedParticleImpl* FitLep1 = new ReconstructedParticleImpl;
       ReconstructedParticleImpl* FitLep2 = new ReconstructedParticleImpl;
       ReconstructedParticleImpl* FitJet1 = new ReconstructedParticleImpl;
       ReconstructedParticleImpl* FitJet2 = new ReconstructedParticleImpl;
       ReconstructedParticleImpl* FitJet3 = new ReconstructedParticleImpl;
       ReconstructedParticleImpl* FitJet4 = new ReconstructedParticleImpl;*/
       //4-momentum of KinFit objects before fit for best combination found in the fit
       //LCCollectionVec *OutputCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
       /*ReconstructedParticleImpl* RecoLep1 = new ReconstructedParticleImpl;
       ReconstructedParticleImpl* RecoLep2 = new ReconstructedParticleImpl;
       ReconstructedParticleImpl* RecoJet1 = new ReconstructedParticleImpl;
       ReconstructedParticleImpl* RecoJet2 = new ReconstructedParticleImpl;
       ReconstructedParticleImpl* RecoJet3 = new ReconstructedParticleImpl;
       ReconstructedParticleImpl* RecoJet4 = new ReconstructedParticleImpl;*/


       OutputCol->parameters().setValue("bestchisq", (float)bestchi2);
       streamlog_out(DEBUG) << " chi2:   " << bestchi2 << std::endl ;
       OutputCol->parameters().setValue("bestprob", (float)bestprob);
       streamlog_out(DEBUG) << " prob:   " << bestprob << std::endl ;
       OutputCol->parameters().setValue("bestmass1", (float)bestmass1);
       OutputCol->parameters().setValue("bestmass2", (float)bestmass2);

       streamlog_out(DEBUG) << "==============  end of permutations for event " << evt->getEventNumber() <<  " ==============" << std::endl ;
       streamlog_out(DEBUG)  << "best mass of H 1: " << bestmass1 << std::endl ;
       streamlog_out(DEBUG)  << "best mass of H 2: " << bestmass2 << std::endl ;


#ifdef MARLIN_USE_AIDA
       if (bestprob > 0) {
         hFitProbBest->fill( bestprob ) ;
         hNItBest->fill( bestnit ) ;
         hRecHMassBest->fill( 0.5*(bestmass1+bestmass2) ) ;
         //hRecHMassBest->fill( bestmass2 ) ;
         hRecHMassNoFitBest->fill( 0.5*(beststartmass1+beststartmass2) ) ;
         //hRecHMassNoFitBest->fill( beststartmass2 ) ;
         hPhotonEnergy->fill( _fitISR ? bestphotonenergy : 0. );
       } 
#endif

       delete j1;
       delete j2;
       delete j3;
       delete j4;
       delete l1;
       delete l2;

       evt->addCollection(OutputCol, _OutCol.c_str() );
  }



  _nEvt ++ ;
}

void ZHHKinFit::getJetResolutions(TLorentzVector jetFourMomentum , std::vector<float> jetCovMat , double &sigmaE , double &sigmaTheta , double &sigmaPhi )
{
  float Px , Py , Pz , P2 , Pt , Pt2;
  float dTheta_dPx , dTheta_dPy , dTheta_dPz , dPhi_dPx , dPhi_dPy;
  float sigmaPx2 , sigmaPy2 , sigmaPz2 , sigmaPxPy , sigmaPxPz , sigmaPyPz;

  Px= jetFourMomentum.Px();
  Py= jetFourMomentum.Py();
  Pz= jetFourMomentum.Pz();
  P2= ( jetFourMomentum.Vect() ).Mag2();
  Pt2= pow( Px , 2 ) + pow( Py , 2 );
  Pt= sqrt( Pt2 );
  sigmaPx2= jetCovMat[ 0 ];
  sigmaPxPy= jetCovMat[ 1 ];
  sigmaPy2= jetCovMat[ 2 ];
  sigmaPxPz= jetCovMat[ 3 ];
  sigmaPyPz= jetCovMat[ 4 ];
  sigmaPz2= jetCovMat[ 5 ];

  dTheta_dPx= Px * Pz / ( P2 * Pt );
  dTheta_dPy= Py * Pz / ( P2 * Pt );
  dTheta_dPz= -Pt / P2;
  dPhi_dPx= -Py / Pt2;
  dPhi_dPy= Px / Pt2;

  sigmaE= std::sqrt( jetCovMat[ 9 ] );
  sigmaTheta= std::sqrt( std::fabs( std::pow( dTheta_dPx , 2 ) * sigmaPx2 + std::pow( dTheta_dPy , 2 ) * sigmaPy2 + std::pow( dTheta_dPz , 2 ) * sigmaPz2 +
				    2.0 * dTheta_dPx * dTheta_dPy * sigmaPxPy + 2.0 * dTheta_dPx * dTheta_dPz * sigmaPxPz + 2.0 * dTheta_dPy * dTheta_dPz * sigmaPyPz ) );
  sigmaPhi= std::sqrt( std::fabs( std::pow( dPhi_dPx , 2 ) * sigmaPx2 + std::pow( dPhi_dPy , 2 ) * sigmaPy2 + 2.0 * dPhi_dPx * dPhi_dPy * sigmaPxPy ) );
  streamlog_out(DEBUG6) << "E= " << jetFourMomentum.E() << std::endl ;
  streamlog_out(DEBUG6) << "Theta= " << jetFourMomentum.Theta() << std::endl ;
  streamlog_out(DEBUG6) << "Phi= " << jetFourMomentum.Phi() << std::endl ;
  streamlog_out(DEBUG6) << "sigmaE= " << sigmaE << std::endl ;
  streamlog_out(DEBUG6) << "sigmaTheta= " << sigmaTheta << std::endl ;
  streamlog_out(DEBUG6) << "sigmaPhi= " << sigmaPhi << std::endl ;

}

void ZHHKinFit::getLeptonParameters( ReconstructedParticle* lepton , float (&parameters)[ 3 ] , float (&errors)[ 3 ] )
{
  TrackVec trackVec = lepton->getTracks();
  if ( trackVec.size() != 1 )
    {
      streamlog_out(DEBUG4)  << "Number of tracks for lepton is not exactly ONE!!! (nTracks = " << trackVec.size() << " ) " << std::endl ;
      //streamlog_out(DEBUG4) << "Input lepton" << *(lepton) << std::endl;
      TLorentzVector leptonFourMomentum( lepton->getMomentum() , lepton->getEnergy() );
      float Px= leptonFourMomentum.Px();
      float Py= leptonFourMomentum.Py();
      float Pz= leptonFourMomentum.Pz();
      float P2= ( leptonFourMomentum.Vect() ).Mag2();
      float Pt2= std::pow( Px , 2 ) + std::pow( Py , 2 );
      float Pt= std::sqrt( Pt2 );

      float sigmaPx2= lepton->getCovMatrix()[ 0 ];
      float sigmaPxPy= lepton->getCovMatrix()[ 1 ];
      float sigmaPy2= lepton->getCovMatrix()[ 2 ];
      float sigmaPxPz= lepton->getCovMatrix()[ 3 ];
      float sigmaPyPz= lepton->getCovMatrix()[ 4 ];
      float sigmaPz2= lepton->getCovMatrix()[ 5 ];

      float dInvPt_dPx= - Px / ( Pt * Pt2 );
      float dInvPt_dPy= - Py / ( Pt * Pt2 );
      float dTheta_dPx= Px * Pz / ( P2 * Pt );
      float dTheta_dPy= Py * Pz / ( P2 * Pt );
      float dTheta_dPz= -Pt / P2;
      float dPhi_dPx= -Py / Pt2;
      float dPhi_dPy= Px / Pt2;

      parameters[ 0 ] = 1.0 / std::sqrt( std::pow( leptonFourMomentum.Px() , 2 ) + std::pow( leptonFourMomentum.Py() , 2 ) );
      parameters[ 1 ] = leptonFourMomentum.Theta();
      parameters[ 2 ] = leptonFourMomentum.Phi();
      errors[ 0 ]= std::sqrt( std::pow( dInvPt_dPx , 2 ) * sigmaPx2 + std::pow( dInvPt_dPy , 2 ) * sigmaPy2 + 2.0 * dInvPt_dPx * dInvPt_dPy * sigmaPxPy );
      errors[ 1 ]= std::sqrt( std::fabs( std::pow( dTheta_dPx , 2 ) * sigmaPx2 + std::pow( dTheta_dPy , 2 ) * sigmaPy2 + std::pow( dTheta_dPz , 2 ) * sigmaPz2 +
					 2.0 * dTheta_dPx * dTheta_dPy * sigmaPxPy + 2.0 * dTheta_dPx * dTheta_dPz * sigmaPxPz + 2.0 * dTheta_dPy * dTheta_dPz * sigmaPyPz ) );
      errors[ 2 ]= std::sqrt( std::fabs( std::pow( dPhi_dPx , 2 ) * sigmaPx2 + std::pow( dPhi_dPy , 2 ) * sigmaPy2 + 2.0 * dPhi_dPx * dPhi_dPy * sigmaPxPy ) );
    }
  else
    {
      streamlog_out(DEBUG4)  << "Lepton has exactly ONE track:" << std::endl ;
      //streamlog_out(DEBUG4) << *lepton << std::endl;
      //streamlog_out(DEBUG4) << *trackVec[ 0 ] << std::endl;
      float Omega= trackVec[ 0 ]->getOmega();
      float tanLambda= trackVec[ 0 ]->getTanLambda();
      float Theta= 2.0 * atan( 1.0 ) - atan( tanLambda );//atan( 1.0 / tanLambda );
      float Phi= trackVec[ 0 ]->getPhi();

      float sigmaOmega= std::sqrt( trackVec[ 0 ]->getCovMatrix()[ 5 ] );
      float sigmaTanLambda= std::sqrt( trackVec[ 0 ]->getCovMatrix()[ 14 ] );
      float sigmaPhi= std::sqrt( trackVec[ 0 ]->getCovMatrix()[ 2 ] );

      float dTheta_dTanLambda= -1.0 / ( 1.0 + std::pow( tanLambda , 2 ) );

      parameters[ 0 ]= Omega / eB;
      parameters[ 1 ]= Theta;
      parameters[ 2 ]= Phi;
      errors[ 0 ]= sigmaOmega / eB;
      errors[ 1 ]= std::fabs( dTheta_dTanLambda ) * sigmaTanLambda;
      errors[ 2 ]= sigmaPhi;
    }
  streamlog_out(DEBUG6) << "Inverse pT= " << parameters[ 0 ] << std::endl ;
  streamlog_out(DEBUG6) << "Theta= " << parameters[ 1 ] << std::endl ;
  streamlog_out(DEBUG6) << "Phi= " << parameters[ 2 ] << std::endl ;
  streamlog_out(DEBUG6) << "SigmaInverse pT= " << errors[ 0 ] << std::endl ;
  streamlog_out(DEBUG6) << "SigmaTheta= " << errors[ 1 ] << std::endl ;
  streamlog_out(DEBUG6) << "SigmaPhi= " << errors[ 2 ] << std::endl ;
}



void ZHHKinFit::check( LCEvent* ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ZHHKinFit::end(){ 

}
