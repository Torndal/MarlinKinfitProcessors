#include "ZHllqq5CFit.h"

using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace CLHEP ;

ZHllqq5CFit aZHllqq5CFit ;

ZHllqq5CFit::ZHllqq5CFit() :
Processor("ZHllqq5CFit"),
m_nRun(0),
m_nEvt(0),
m_nRunSum(0),
m_nEvtSum(0),
m_Bfield(0.0),
c(0.0),
mm2m(0.0),
eV2GeV(0.0),
eB(0.0)
{

//	modify processor description
	_description = "ZHllqq5CFit does a fit on 2 jet events (Px, Py, Pz, E, M12 = MZ)" ;

//	register steering parameters: name, description, class-variable, default value

	registerInputCollection( 	LCIO::RECONSTRUCTEDPARTICLE,
					"isolatedlaptonCollection" ,
					"Name of the Isolated Lepton collection"  ,
					m_inputIsolatedlaptonCollection ,
					std::string("ISOLeptons")
				);

	registerInputCollection( 	LCIO::RECONSTRUCTEDPARTICLE,
					"JetCollectionName" ,
					"Name of the Jet collection"  ,
					m_inputJetCollection ,
					std::string("Durham_2Jets")
				);

	registerInputCollection( 	LCIO::VERTEX,
					"SLDVertexCollection" ,
					"Name of Semi-Leptonic Decay Vertices Collection"  ,
					m_inputSLDVertexCollection ,
					std::string("SemiLeptonicDecayVertex")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"JetSLDRelationCollection",
					"Name of the Jet-SemiLeptonicDecay Relation collection",
					m_inputJetSLDLink,
					std::string("JetSLDLinkName")
				);

	registerOutputCollection( 	LCIO::RECONSTRUCTEDPARTICLE,
					"FitOutputColection" ,
					"Name of Output Fit collection"  ,
					m_outputFitcollection ,
					std::string("FitReco")
				);

	registerOutputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"outputJetCollection",
					"Name of output jet collection",
					m_outputJetCollection,
					std::string("Durham_2JetsKinFit")
				);

	registerInputCollection( 	LCIO::RECONSTRUCTEDPARTICLE,
					"TrueJets" ,
					"Name of the TrueJetCollection input collection",
					_trueJetCollectionName ,
					std::string("TrueJets")
				);

	registerInputCollection( 	LCIO::RECONSTRUCTEDPARTICLE,
					"FinalColourNeutrals" ,
					"Name of the FinalColourNeutralCollection input collection"  ,
					_finalColourNeutralCollectionName ,
					std::string("FinalColourNeutrals")
				);

	registerInputCollection( 	LCIO::RECONSTRUCTEDPARTICLE,
					"InitialColourNeutrals" ,
					"Name of the InitialColourNeutralCollection input collection"  ,
					_initialColourNeutralCollectionName ,
					std::string("InitialColourNeutrals")
				);

	registerInputCollection( 	LCIO::LCRELATION,
					"TrueJetPFOLink" ,
					"Name of the TrueJetPFOLink input collection"  ,
					_trueJetPFOLink,
					std::string("TrueJetPFOLink")
				);

	registerInputCollection( 	LCIO::LCRELATION,
					"TrueJetMCParticleLink" ,
					"Name of the TrueJetMCParticleLink input collection"  ,
					_trueJetMCParticleLink,
					std::string("TrueJetMCParticleLink")
				);

	registerInputCollection( 	LCIO::LCRELATION,
					"FinalElementonLink rueJetMCParticleLink" ,
					"Name of the  FinalElementonLink input collection"	,
					_finalElementonLink,
					std::string("FinalElementonLink")
				);

	registerInputCollection( 	LCIO::LCRELATION,
					"InitialElementonLink" ,
					"Name of the  InitialElementonLink input collection"  ,
					_initialElementonLink,
					std::string("InitialElementonLink")
				);

	registerInputCollection( 	LCIO::LCRELATION,
					"FinalColourNeutralLink" ,
					"Name of the  FinalColourNeutralLink input collection"  ,
					_finalColourNeutralLink,
					std::string("FinalColourNeutralLink")
				);

	registerInputCollection( 	LCIO::LCRELATION,
					"InitialColourNeutralLink" ,
					"Name of the  InitialColourNeutralLink input collection"  ,
					_initialColourNeutralLink,
					std::string("InitialColourNeutralLink")
				);

	registerProcessorParameter(	"nJets",
					"Number of jet should be in the event",
					m_nAskedJets,
					int(2)
				);

	registerProcessorParameter(	"nIsoLeps",
					"Number of Isolated Leptons should be in the event",
					m_nAskedIsoLeps,
					int(2)
				);

	registerProcessorParameter(	"ECM" ,
					"Center-of-Mass Energy in GeV",
					m_ECM,
					float(250.f)
				);

	registerProcessorParameter(	"ISRPzMax" ,
					"Maximum possible energy for a single ISR photon",
					m_isrpzmax,
					float(35.f)
				);

	registerProcessorParameter(	"includeISR",
					"Include ISR in fit hypothesis; false: without ISR , true: with ISR",
					m_fitISR,
					bool(true)
				);

	registerProcessorParameter(	"fitter" ,
					"0 = OPALFitter, 1 = NewFitter, 2 = NewtonFitter",
					m_fitter,
					int(0)
				);

	registerProcessorParameter(	"outputFilename",
					"name of output root file",
					m_outputFile,
					std::string("")
				);

	registerProcessorParameter(	"traceall" ,
					"set true if every event should be traced",
					m_traceall,
					(bool)false
				);

	registerProcessorParameter(	"ievttrace" ,
					"number of individual event to be traced",
					m_ievttrace,
					(int)0
				);

}

void ZHllqq5CFit::init()
{
//	usually a good idea to
	streamlog_out(DEBUG) << "   init called  " << std::endl;
	this->Clear();
	m_nRun = 0;
	m_nEvt = 0;
	m_nRunSum = 0;
	m_nEvtSum = 0;
	DDMarlinCED::init(this);

	m_Bfield = MarlinUtil::getBzAtOrigin();
	streamlog_out(DEBUG0) << " BField =  "<< m_Bfield << " Tesla" << std::endl ;
	c = 2.99792458e8;
	mm2m = 1e-3;
	eV2GeV = 1e-9;
	eB = m_Bfield * c * mm2m * eV2GeV;


	b = (double) 0.00464564 * ( std::log( m_ECM * m_ECM * 3814714. ) - 1. );
//	  = 2*alpha/pi*( ln(s/m_e^2)-1 )
	ISRPzMaxB = std::pow((double)m_isrpzmax,b);

	printParameters();

	m_pTFile = new TFile(m_outputFile.c_str(),"recreate");

	m_pTTree = new TTree("eventTree","eventTree");
	m_pTTree->SetDirectory(m_pTFile);
	m_pTTree->Branch("run", &m_nRun, "run/I");
	m_pTTree->Branch("event", &m_nEvt, "event/I");
	m_pTTree->Branch("nJets",&m_nJets,"nJets/I") ;
	m_pTTree->Branch("nIsoLeptons",&m_nIsoLeps,"nIsoLeptons/I") ;
	m_pTTree->Branch("nSLDecayBHadron",&m_nSLDecayBHadron,"nSLDecayBHadron/I") ;
	m_pTTree->Branch("nSLDecayCHadron",&m_nSLDecayCHadron,"nSLDecayCHadron/I") ;
	m_pTTree->Branch("nSLDecayTauLepton",&m_nSLDecayTauLepton,"nSLDecayTauLepton/I") ;
	m_pTTree->Branch("nSLDecayTotal",&m_nSLDecayTotal,"nSLDecayTotal/I") ;
	m_pTTree->Branch( "FitErrorCode_woNu" , &m_FitErrorCode_woNu , "FitErrorCode_woNu/I" );
	m_pTTree->Branch( "ZMassBeforeFit_woNu" , &m_ZMassBeforeFit_woNu , "ZMassBeforeFit_woNu/F" );
	m_pTTree->Branch( "HMassBeforeFit_woNu" , &m_HMassBeforeFit_woNu , "HMassBeforeFit_woNu/F" );
	m_pTTree->Branch( "ZMassAfterFit_woNu" , &m_ZMassAfterFit_woNu , "ZMassAfterFit_woNu/F" );
	m_pTTree->Branch( "HMassAfterFit_woNu" , &m_HMassAfterFit_woNu , "HMassAfterFit_woNu/F" );
	m_pTTree->Branch( "FitProbability_woNu" , &m_FitProbability_woNu , "FitProbability_woNu/F" );
	m_pTTree->Branch( "pullJetEnergy_woNu" , &m_pullJetEnergy_woNu );
	m_pTTree->Branch( "pullJetTheta_woNu" , &m_pullJetTheta_woNu );
	m_pTTree->Branch( "pullJetPhi_woNu" , &m_pullJetPhi_woNu );
	m_pTTree->Branch( "pullLeptonInvPt_woNu" , &m_pullLeptonInvPt_woNu );
	m_pTTree->Branch( "pullLeptonTheta_woNu" , &m_pullLeptonTheta_woNu );
	m_pTTree->Branch( "pullLeptonPhi_woNu" , &m_pullLeptonPhi_woNu );
	m_pTTree->Branch( "normalizedResidualJetEnergy_woNu" , &m_normalizedResidualJetEnergy_woNu );
	m_pTTree->Branch( "normalizedResidualJetTheta_woNu" , &m_normalizedResidualJetTheta_woNu );
	m_pTTree->Branch( "normalizedResidualJetPhi_woNu" , &m_normalizedResidualJetPhi_woNu );
	m_pTTree->Branch( "normalizedResidualLeptonInvPt_woNu" , &m_normalizedResidualLeptonInvPt_woNu );
	m_pTTree->Branch( "normalizedResidualLeptonTheta_woNu" , &m_normalizedResidualLeptonTheta_woNu );
	m_pTTree->Branch( "normalizedResidualLeptonPhi_woNu" , &m_normalizedResidualLeptonPhi_woNu );
	m_pTTree->Branch( "FitErrorCode_wNu" , &m_FitErrorCode_wNu , "FitErrorCode_wNu/I" );
	m_pTTree->Branch( "ZMassBeforeFit_wNu" , &m_ZMassBeforeFit_wNu , "ZMassBeforeFit_wNu/F" );
	m_pTTree->Branch( "HMassBeforeFit_wNu" , &m_HMassBeforeFit_wNu , "HMassBeforeFit_wNu/F" );
	m_pTTree->Branch( "ZMassAfterFit_wNu" , &m_ZMassAfterFit_wNu , "ZMassAfterFit_wNu/F" );
	m_pTTree->Branch( "HMassAfterFit_wNu" , &m_HMassAfterFit_wNu , "HMassAfterFit_wNu/F" );
	m_pTTree->Branch( "FitProbability_wNu" , &m_FitProbability_wNu , "FitProbability_wNu/F" );
	m_pTTree->Branch( "pullJetEnergy_wNu" , &m_pullJetEnergy_wNu );
	m_pTTree->Branch( "pullJetTheta_wNu" , &m_pullJetTheta_wNu );
	m_pTTree->Branch( "pullJetPhi_wNu" , &m_pullJetPhi_wNu );
	m_pTTree->Branch( "pullLeptonInvPt_wNu" , &m_pullLeptonInvPt_wNu );
	m_pTTree->Branch( "pullLeptonTheta_wNu" , &m_pullLeptonTheta_wNu );
	m_pTTree->Branch( "pullLeptonPhi_wNu" , &m_pullLeptonPhi_wNu );
	m_pTTree->Branch( "normalizedResidualJetEnergy_wNu" , &m_normalizedResidualJetEnergy_wNu );
	m_pTTree->Branch( "normalizedResidualJetTheta_wNu" , &m_normalizedResidualJetTheta_wNu );
	m_pTTree->Branch( "normalizedResidualJetPhi_wNu" , &m_normalizedResidualJetPhi_wNu );
	m_pTTree->Branch( "normalizedResidualLeptonInvPt_wNu" , &m_normalizedResidualLeptonInvPt_wNu );
	m_pTTree->Branch( "normalizedResidualLeptonTheta_wNu" , &m_normalizedResidualLeptonTheta_wNu );
	m_pTTree->Branch( "normalizedResidualLeptonPhi_wNu" , &m_normalizedResidualLeptonPhi_wNu );
	m_pTTree->Branch( "FitErrorCode" , &m_FitErrorCode , "FitErrorCode/I" );
	m_pTTree->Branch( "ZMassBeforeFit" , &m_ZMassBeforeFit , "ZMassBeforeFit/F" );
	m_pTTree->Branch( "HMassBeforeFit" , &m_HMassBeforeFit , "HMassBeforeFit/F" );
	m_pTTree->Branch( "ZMassAfterFit" , &m_ZMassAfterFit , "ZMassAfterFit/F" );
	m_pTTree->Branch( "HMassAfterFit" , &m_HMassAfterFit , "HMassAfterFit/F" );
	m_pTTree->Branch( "FitProbability" , &m_FitProbability , "FitProbability/F" );
	m_pTTree->Branch( "pullJetEnergy" , &m_pullJetEnergy );
	m_pTTree->Branch( "pullJetTheta" , &m_pullJetTheta );
	m_pTTree->Branch( "pullJetPhi" , &m_pullJetPhi );
	m_pTTree->Branch( "pullLeptonInvPt" , &m_pullLeptonInvPt );
	m_pTTree->Branch( "pullLeptonTheta" , &m_pullLeptonTheta );
	m_pTTree->Branch( "pullLeptonPhi" , &m_pullLeptonPhi );
	m_pTTree->Branch( "normalizedResidualJetEnergy" , &m_normalizedResidualJetEnergy );
	m_pTTree->Branch( "normalizedResidualJetTheta" , &m_normalizedResidualJetTheta );
	m_pTTree->Branch( "normalizedResidualJetPhi" , &m_normalizedResidualJetPhi );
	m_pTTree->Branch( "normalizedResidualLeptonInvPt" , &m_normalizedResidualLeptonInvPt );
	m_pTTree->Branch( "normalizedResidualLeptonTheta" , &m_normalizedResidualLeptonTheta );
	m_pTTree->Branch( "normalizedResidualLeptonPhi" , &m_normalizedResidualLeptonPhi );

	streamlog_out(DEBUG) << "   init finished  " << std::endl;

}

void ZHllqq5CFit::Clear()
{
	streamlog_out(DEBUG) << "   Clear called  " << std::endl;

	m_nJets = 0;
	m_nIsoLeps = 0;
	m_nSLDecayBHadron = 0;
	m_nSLDecayCHadron = 0;
	m_nSLDecayTauLepton = 0;
	m_nSLDecayTotal = 0;
	m_FitErrorCode_woNu = 0;
	m_ZMassBeforeFit_woNu = 0.0;
	m_HMassBeforeFit_woNu = 0.0;
	m_ZMassAfterFit_woNu = 0.0;
	m_HMassAfterFit_woNu = 0.0;
	m_FitProbability_woNu = 0.0;
	m_pullJetEnergy_woNu.clear();
	m_pullJetTheta_woNu.clear();
	m_pullJetPhi_woNu.clear();
	m_pullLeptonInvPt_woNu.clear();
	m_pullLeptonTheta_woNu.clear();
	m_pullLeptonPhi_woNu.clear();
	m_normalizedResidualJetEnergy_woNu.clear();
	m_normalizedResidualJetTheta_woNu.clear();
	m_normalizedResidualJetPhi_woNu.clear();
	m_normalizedResidualLeptonInvPt_woNu.clear();
	m_normalizedResidualLeptonTheta_woNu.clear();
	m_normalizedResidualLeptonPhi_woNu.clear();
	m_FitErrorCode_wNu = 0;
	m_ZMassBeforeFit_wNu = 0.0;
	m_HMassBeforeFit_wNu = 0.0;
	m_ZMassAfterFit_wNu = 0.0;
	m_HMassAfterFit_wNu = 0.0;
	m_FitProbability_wNu = 0.0;
	m_pullJetEnergy_wNu.clear();
	m_pullJetTheta_wNu.clear();
	m_pullJetPhi_wNu.clear();
	m_pullLeptonInvPt_wNu.clear();
	m_pullLeptonTheta_wNu.clear();
	m_pullLeptonPhi_wNu.clear();
	m_normalizedResidualJetEnergy_wNu.clear();
	m_normalizedResidualJetTheta_wNu.clear();
	m_normalizedResidualJetPhi_wNu.clear();
	m_normalizedResidualLeptonInvPt_wNu.clear();
	m_normalizedResidualLeptonTheta_wNu.clear();
	m_normalizedResidualLeptonPhi_wNu.clear();
	m_FitErrorCode = 0;
	m_ZMassBeforeFit = 0.0;
	m_HMassBeforeFit = 0.0;
	m_ZMassAfterFit = 0.0;
	m_HMassAfterFit = 0.0;
	m_FitProbability = 0.0;
	m_pullJetEnergy.clear();
	m_pullJetTheta.clear();
	m_pullJetPhi.clear();
	m_pullLeptonInvPt.clear();
	m_pullLeptonTheta.clear();
	m_pullLeptonPhi.clear();
	m_normalizedResidualJetEnergy.clear();
	m_normalizedResidualJetTheta.clear();
	m_normalizedResidualJetPhi.clear();
	m_normalizedResidualLeptonInvPt.clear();
	m_normalizedResidualLeptonTheta.clear();
	m_normalizedResidualLeptonPhi.clear();

}

void ZHllqq5CFit::processRunHeader()
{
	m_nRun++ ;
}

void ZHllqq5CFit::processEvent( EVENT::LCEvent *pLCEvent )
{
	this->Clear();
	m_nRun = pLCEvent->getRunNumber();
	m_nEvt = pLCEvent->getEventNumber();
	streamlog_out(MESSAGE) << "	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	///////////////////////////////////////////////// processing event " << m_nEvt << " in run " << m_nRun << " /////////////////////////////////////////////////////////////////////" << std::endl ;
	streamlog_out(MESSAGE) << "	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;

	LCCollection *inputJetCollection{};
	LCCollection *inputLeptonCollection{};
	LCCollection *inputSLDecayCollection{};
	IMPL::LCCollectionVec* outputJetCollection(NULL);
	outputJetCollection = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	outputJetCollection->setSubset( true );
	try
	{
		streamlog_out(DEBUG0) << "	getting jet collection: " << m_inputJetCollection << std::endl ;
		inputJetCollection = pLCEvent->getCollection( m_inputJetCollection );
		streamlog_out(DEBUG0) << "	getting isolated lepton collection: " << m_inputIsolatedlaptonCollection << std::endl ;
		inputLeptonCollection = pLCEvent->getCollection( m_inputIsolatedlaptonCollection );
		streamlog_out(DEBUG0) << "	getting semi-leptonic vertex collection: " << m_inputSLDVertexCollection << std::endl ;
		inputSLDecayCollection = pLCEvent->getCollection( m_inputSLDVertexCollection );
		m_nCorrectedSLD = inputSLDecayCollection->getNumberOfElements();
		m_nSLDecayBHadron = inputSLDecayCollection->getParameters().getIntVal( "nBHadronSLD_found" );
		m_nSLDecayCHadron = inputSLDecayCollection->getParameters().getIntVal( "nCHadronSLD_found" );
		m_nSLDecayTauLepton = inputSLDecayCollection->getParameters().getIntVal( "nTauLeptonSLD_found" );
		m_nSLDecayTotal = inputSLDecayCollection->getParameters().getIntVal( "nTotalSLD_found" );
		m_nJets = inputJetCollection->getNumberOfElements();
		m_nIsoLeps = inputLeptonCollection->getNumberOfElements();
		streamlog_out(DEBUG8) << "	Number of jets: " << m_nJets << std::endl ;
		streamlog_out(DEBUG8) << "	Number of isolatedLeptons: " << m_nIsoLeps << std::endl ;
		streamlog_out(DEBUG8) << "	Number of found semi-leptonic decays: " << m_nSLDecayTotal << std::endl ;
		streamlog_out(DEBUG8) << "	Number of corrected semi-leptonic decays: " << m_nCorrectedSLD << std::endl ;
		if ( m_nJets != m_nAskedJets || m_nIsoLeps != m_nAskedIsoLeps || m_nCorrectedSLD != m_nSLDecayTotal ) return;
		bool traceEvent = false;
		if ( pLCEvent->getEventNumber() == m_ievttrace || m_traceall ) traceEvent = true;

		std::vector< ReconstructedParticle* > Leptons{};
		for ( int i_lep = 0 ; i_lep < m_nIsoLeps ; ++i_lep )
		{
			ReconstructedParticle* lepton = dynamic_cast<ReconstructedParticle*>( inputLeptonCollection->getElementAt( i_lep ) );
			Leptons.push_back( lepton );
		}
		ReconstructedParticle* jet1 = dynamic_cast<ReconstructedParticle*>( inputJetCollection->getElementAt( 0 ) );
		TLorentzVector jet1tlv( jet1->getMomentum() , jet1->getEnergy() );
		std::vector< float > jet1initialCovMat( 10 , 0.0 );
		for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element )
		{
			jet1initialCovMat[ i_Element ] = jet1->getCovMatrix()[ i_Element ];
		}
		ReconstructedParticle* jet2 = dynamic_cast<ReconstructedParticle*>( inputJetCollection->getElementAt( 1 ) );
		TLorentzVector jet2tlv( jet2->getMomentum() , jet2->getEnergy() );
		std::vector< float > jet2initialCovMat( 10 , 0.0 );
		for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element )
		{
			jet2initialCovMat[ i_Element ] = jet2->getCovMatrix()[ i_Element ];
		}

		float fitProbability = 0.0;
		float fitProbabilityBestFit_wNu = 0.0;
		float fitProbability_woNu = 0.0;
		int fitErrorCode = 0;
		float fitOutputs_temp[ 18 ]{ 0.0 };
		std::vector< TLorentzVector > fittedObjects_temp{};
		float pull_temp[ 12 ]{ 0.0 };
		float fitOutputs_woNu[ 18 ]{ 0.0 };
		std::vector< TLorentzVector > fittedObjects_woNu{};
		float pull_woNu[ 12 ]{ 0.0 };
		float fitOutputs_wNu[ 18 ]{ 0.0 };
		std::vector< TLorentzVector > fittedObjects_wNu{};
		float pull_wNu[ 12 ]{ 0.0 };
		float fitOutputs[ 18 ]{ 0.0 };
		std::vector< TLorentzVector > fittedObjects{};
		float pull[ 12 ]{ 0.0 };

		TLorentzVector Nu1tlv( 0.0 , 0.0 , 0.0 , 0.0 );
		TLorentzVector Nu2tlv( 0.0 , 0.0 , 0.0 , 0.0 );
		std::vector< float > nu1CovMat( 10 , 0.0 );
		std::vector< float > nu2CovMat( 10 , 0.0 );

		m_FitErrorCode_woNu = performFIT( jet1tlv , jet1initialCovMat , jet2tlv , jet2initialCovMat , Leptons , fitProbability_woNu , fitOutputs_temp , fittedObjects_temp , pull_temp , traceEvent );
		if ( m_FitErrorCode_woNu == 0 )
		{
			for ( unsigned int i = 0 ; i < sizeof( fitOutputs_woNu ) / sizeof( fitOutputs_woNu[ 0 ] ) ; ++i ) fitOutputs_woNu[ i ] = fitOutputs_temp[ i ];
			for ( unsigned int i = 0 ; i < sizeof( pull_woNu ) / sizeof( pull_woNu[ 0 ] ) ; ++i ) pull_woNu[ i ] = pull_temp[ i ];
			for ( unsigned int i = 0 ; i < fittedObjects_woNu.size()  ; ++i ) fittedObjects_woNu.push_back( fittedObjects_temp[ i ] );
			m_ZMassBeforeFit_woNu = fitOutputs_woNu[ 2 ];
			m_HMassBeforeFit_woNu = fitOutputs_woNu[ 3 ];
			m_ZMassAfterFit_woNu = fitOutputs_woNu[ 4 ];
			m_HMassAfterFit_woNu = fitOutputs_woNu[ 5 ];
			m_FitProbability_woNu = fitProbability_woNu;
		}

		LCRelationNavigator JetSLDNav( pLCEvent->getCollection( m_inputJetSLDLink ) );

		TLorentzVector jet1FourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
		std::vector< float > jet1CovMat( 10 , 0.0 );
		pfoVector jet1ZeroNeutrinos;
		getSLDsInJet( JetSLDNav , jet1 , jet1ZeroNeutrinos );
		std::vector< int > jet1nSLDSolutions{};
		int njet1SLDSolutions = 1;
		streamlog_out(DEBUG6) << "	Number of semi-leptonic decays in jet1: " << jet1ZeroNeutrinos.size() << std::endl ;
		for ( unsigned int i_sld = 0 ; i_sld < jet1ZeroNeutrinos.size() ; ++i_sld )
		{
			int nNeutrinos = jet1ZeroNeutrinos[ i_sld ]->getParticles().size() + 1;
			streamlog_out(DEBUG6) << "	Number of Neutrinos for semi-leptonic decay[ " << i_sld << " ] in jet1: " << nNeutrinos << std::endl ;
			jet1nSLDSolutions.push_back( nNeutrinos );
			njet1SLDSolutions *= nNeutrinos;
		}

		TLorentzVector jet2FourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
		std::vector< float > jet2CovMat( 10 , 0.0 );
		pfoVector jet2ZeroNeutrinos;
		getSLDsInJet( JetSLDNav , jet2 , jet2ZeroNeutrinos );
		std::vector< int > jet2nSLDSolutions{};
		int njet2SLDSolutions = 1;
		streamlog_out(DEBUG6) << "	Number of semi-leptonic decays in jet2: " << jet2ZeroNeutrinos.size() << std::endl ;
		for ( unsigned int i_sld = 0 ; i_sld < jet2ZeroNeutrinos.size() ; ++i_sld )
		{
			int nNeutrinos = jet2ZeroNeutrinos[ i_sld ]->getParticles().size() + 1;
			streamlog_out(DEBUG6) << "	Number of Neutrinos for semi-leptonic decay[ " << i_sld << " ] in jet2: " << nNeutrinos << std::endl ;
			jet2nSLDSolutions.push_back( nNeutrinos );
			njet2SLDSolutions *= nNeutrinos;
		}

		std::vector<int> jet1SLDCombination( jet1ZeroNeutrinos.size() , 0 );
		std::vector<int> jet2SLDCombination( jet2ZeroNeutrinos.size() , 0 );
		std::vector<int> jet1FinalNuSolutions( jet1ZeroNeutrinos.size() , 0 );
		std::vector<int> jet2FinalNuSolutions( jet2ZeroNeutrinos.size() , 0 );
		for ( int i_jet1 = 0 ; i_jet1 < njet1SLDSolutions ; ++i_jet1 )
		{
			streamlog_out(DEBUG6) << "" << std::endl ;
			jet1FourMomentum = jet1tlv;
			for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element ) jet1CovMat[ i_Element ] = jet1initialCovMat[ i_Element ];
			getSLDCombination( jet1nSLDSolutions , i_jet1 , jet1SLDCombination );
			streamlog_out(DEBUG6) << "	Preparing Jet1 for kinematic fit with " << jet1ZeroNeutrinos.size() << " semi-leptonic decays:" << std::endl ;
			TLorentzVector Nu1FourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
			std::vector< float > Nu1CovMat( 10 , 0.0 );
			for ( unsigned int i_sld1 = 0 ; i_sld1 < jet1ZeroNeutrinos.size() ; ++i_sld1 )
			{
				streamlog_out(DEBUG6) << "		Adding solution [" << jet1SLDCombination[ i_sld1 ] << "] From semi-leptonic decay[ " << i_sld1 << " ] to Jet1" << std::endl ;
				ReconstructedParticle* Neutrino1 = NULL;
				if ( jet1SLDCombination[ i_sld1 ] == 0 )
				{
					Neutrino1 = jet1ZeroNeutrinos[ i_sld1 ];
				}
				else
				{
					Neutrino1 = jet1ZeroNeutrinos[ i_sld1 ]->getParticles()[ jet1SLDCombination[ i_sld1 ] - 1 ];
				}
				streamlog_out(DEBUG1) << *Neutrino1 << std::endl;
				Nu1FourMomentum += TLorentzVector( Neutrino1->getMomentum() , Neutrino1->getEnergy() );
				for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element ) Nu1CovMat[ i_Element ] += Neutrino1->getCovMatrix()[ i_Element ];
			}
			jet1FourMomentum += Nu1FourMomentum;
			for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element ) jet1CovMat[ i_Element ] += Nu1CovMat[ i_Element ];
			for ( int i_jet2 = 0 ; i_jet2 < njet2SLDSolutions ; ++i_jet2 )
	 		{
				jet2FourMomentum = jet2tlv;
				for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element ) jet2CovMat[ i_Element ] = jet2initialCovMat[ i_Element ];
	 			getSLDCombination( jet2nSLDSolutions , i_jet2 , jet2SLDCombination );
				streamlog_out(DEBUG6) << "	Preparing Jet2 for kinematic fit with " << jet2ZeroNeutrinos.size() << " semi-leptonic decays:" << std::endl ;
				TLorentzVector Nu2FourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
				std::vector< float > Nu2CovMat( 10 , 0.0 );
				for ( unsigned int i_sld2 = 0 ; i_sld2 < jet2ZeroNeutrinos.size() ; ++i_sld2 )
				{
					streamlog_out(DEBUG6) << "		Adding solution [" << jet2SLDCombination[ i_sld2 ] << "] From semi-leptonic decay[ " << i_sld2 << " ] to Jet1" << std::endl ;
					ReconstructedParticle* Neutrino2 = NULL;
					if ( jet2SLDCombination[ i_sld2 ] == 0 )
					{
						Neutrino2 = jet2ZeroNeutrinos[ i_sld2 ];
					}
					else
					{
						Neutrino2 = jet2ZeroNeutrinos[ i_sld2 ]->getParticles()[ jet2SLDCombination[ i_sld2 ] - 1 ];
					}
					streamlog_out(DEBUG1) << *Neutrino2 << std::endl;
					Nu2FourMomentum += TLorentzVector( Neutrino2->getMomentum() , Neutrino2->getEnergy() );
					for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element ) Nu2CovMat[ i_Element ] += Neutrino2->getCovMatrix()[ i_Element ];
				}
				jet2FourMomentum += Nu2FourMomentum;
				for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element ) jet2CovMat[ i_Element ] += Nu2CovMat[ i_Element ];
				fitErrorCode = performFIT( jet1FourMomentum , jet1CovMat , jet2FourMomentum , jet2CovMat , Leptons , fitProbability , fitOutputs_temp , fittedObjects_temp , pull_temp , traceEvent );
				if ( fitProbability >= fitProbabilityBestFit_wNu ) m_FitErrorCode_wNu = fitErrorCode;
				if ( fitErrorCode == 0 )
				{
					if ( fitProbability >= fitProbabilityBestFit_wNu )
					{
						for ( unsigned int i = 0 ; i < sizeof( fitOutputs_wNu ) / sizeof( fitOutputs_wNu[ 0 ] ) ; ++i ) fitOutputs_woNu[ i ] = fitOutputs_temp[ i ];
						for ( unsigned int i = 0 ; i < sizeof( pull_wNu ) / sizeof( pull_wNu[ 0 ] ) ; ++i ) pull_woNu[ i ] = pull_temp[ i ];
						for ( unsigned int i = 0 ; i < fittedObjects_wNu.size()  ; ++i ) fittedObjects_wNu.push_back( fittedObjects_temp[ i ] );
						for ( unsigned int i_sld1 = 0 ; i_sld1 < jet1SLDCombination.size() ; ++i_sld1 ) jet1FinalNuSolutions[ i_sld1 ] = jet1SLDCombination[ i_sld1 ];
						for ( unsigned int i_sld2 = 0 ; i_sld2 < jet2SLDCombination.size() ; ++i_sld2 ) jet2FinalNuSolutions[ i_sld2 ] = jet2SLDCombination[ i_sld2 ];
						fitProbabilityBestFit_wNu = fitProbability;
						m_ZMassBeforeFit_wNu = fitOutputs_wNu[ 2 ];
						m_HMassBeforeFit_wNu = fitOutputs_wNu[ 3 ];
						m_ZMassAfterFit_wNu = fitOutputs_wNu[ 4 ];
						m_HMassAfterFit_wNu = fitOutputs_wNu[ 5 ];
						m_FitProbability_wNu = fitProbabilityBestFit_wNu;
					}
				}
	 		}
		}
		std::vector<EVENT::ReconstructedParticle*> jet1Neutrinos{};
		std::vector<EVENT::ReconstructedParticle*> jet2Neutrinos{};
		if ( m_FitProbability_wNu > m_FitProbability_woNu )
		{
			m_FitErrorCode = m_FitErrorCode_wNu;
			m_ZMassBeforeFit = m_ZMassBeforeFit_wNu;
			m_HMassBeforeFit = m_HMassBeforeFit_wNu;
			m_ZMassAfterFit = m_ZMassAfterFit_wNu;
			m_HMassAfterFit = m_HMassAfterFit_wNu;
			m_FitProbability = m_FitProbability_wNu;
			for ( unsigned int i = 0 ; i < sizeof( fitOutputs ) / sizeof( fitOutputs[ 0 ] ) ; ++i ) fitOutputs[ i ] = fitOutputs_wNu[ i ];
			for ( unsigned int i = 0 ; i < sizeof( pull ) / sizeof( pull[ 0 ] ) ; ++i ) pull[ i ] = pull_wNu[ i ];
			for ( unsigned int i = 0 ; i < fittedObjects_wNu.size()  ; ++i ) fittedObjects.push_back( fittedObjects_wNu[ i ] );
			for ( unsigned int i = 0 ; i < jet1FinalNuSolutions.size() ; ++i )
			{
				if ( jet1FinalNuSolutions[ i ] == 0 )
				{
					jet1Neutrinos.push_back( ( EVENT::ReconstructedParticle* )jet1ZeroNeutrinos[ i ] );
				}
				else
				{
					jet1Neutrinos.push_back( ( EVENT::ReconstructedParticle* )jet1ZeroNeutrinos[ i ]->getParticles()[ jet1FinalNuSolutions[ i ] - 1 ] );
				}
			}
			for ( unsigned int i = 0 ; i < jet2FinalNuSolutions.size() ; ++i )
			{
				if ( jet2FinalNuSolutions[ i ] == 0 )
				{
					jet2Neutrinos.push_back( ( EVENT::ReconstructedParticle* )jet2ZeroNeutrinos[ i ] );
				}
				else
				{
					jet2Neutrinos.push_back( ( EVENT::ReconstructedParticle* )jet2ZeroNeutrinos[ i ]->getParticles()[ jet2FinalNuSolutions[ i ] - 1 ] );
				}
			}
		}
		else
		{
			m_FitErrorCode = m_FitErrorCode_woNu;
			m_ZMassBeforeFit = m_ZMassBeforeFit_woNu;
			m_HMassBeforeFit = m_HMassBeforeFit_woNu;
			m_ZMassAfterFit = m_ZMassAfterFit_woNu;
			m_HMassAfterFit = m_HMassAfterFit_woNu;
			m_FitProbability = m_FitProbability_woNu;
			for ( unsigned int i = 0 ; i < sizeof( fitOutputs ) / sizeof( fitOutputs[ 0 ] ) ; ++i ) fitOutputs[ i ] = fitOutputs_woNu[ i ];
			for ( unsigned int i = 0 ; i < sizeof( pull ) / sizeof( pull[ 0 ] ) ; ++i ) pull[ i ] = pull_woNu[ i ];
			for ( unsigned int i = 0 ; i < fittedObjects_woNu.size()  ; ++i ) fittedObjects.push_back( fittedObjects_woNu[ i ] );
		}
		m_pullJetEnergy_woNu.push_back( pull_woNu[ 0 ] );	m_pullJetEnergy_woNu.push_back( pull_woNu[ 1 ] );
		m_pullJetTheta_woNu.push_back( pull_woNu[ 2 ] ); 	m_pullJetTheta_woNu.push_back( pull_woNu[ 3 ] );
		m_pullJetPhi_woNu.push_back( pull_woNu[ 4 ] );		m_pullJetPhi_woNu.push_back( pull_woNu[ 5 ] );
		m_pullLeptonInvPt_woNu.push_back( pull_woNu[ 6 ] );	m_pullLeptonInvPt_woNu.push_back( pull_woNu[ 7 ] );
		m_pullLeptonTheta_woNu.push_back( pull_woNu[ 8 ] );	m_pullLeptonTheta_woNu.push_back( pull_woNu[ 9 ] );
		m_pullLeptonPhi_woNu.push_back( pull_woNu[ 10 ] );	m_pullLeptonPhi_woNu.push_back( pull_woNu[ 11 ] );
		m_pullJetEnergy_wNu.push_back( pull_wNu[ 0 ] );		m_pullJetEnergy_wNu.push_back( pull_wNu[ 1 ] );
		m_pullJetTheta_wNu.push_back( pull_wNu[ 2 ] ); 		m_pullJetTheta_wNu.push_back( pull_wNu[ 3 ] );
		m_pullJetPhi_wNu.push_back( pull_wNu[ 4 ] );		m_pullJetPhi_wNu.push_back( pull_wNu[ 5 ] );
		m_pullLeptonInvPt_wNu.push_back( pull_wNu[ 6 ] );	m_pullLeptonInvPt_wNu.push_back( pull_wNu[ 7 ] );
		m_pullLeptonTheta_wNu.push_back( pull_wNu[ 8 ] );	m_pullLeptonTheta_wNu.push_back( pull_wNu[ 9 ] );
		m_pullLeptonPhi_wNu.push_back( pull_wNu[ 10 ] );	m_pullLeptonPhi_wNu.push_back( pull_wNu[ 11 ] );
		m_pullJetEnergy.push_back( pull[ 0 ] );			m_pullJetEnergy.push_back( pull[ 1 ] );
		m_pullJetTheta.push_back( pull[ 2 ] ); 			m_pullJetTheta.push_back( pull[ 3 ] );
		m_pullJetPhi.push_back( pull[ 4 ] );			m_pullJetPhi.push_back( pull[ 5 ] );
		m_pullLeptonInvPt.push_back( pull[ 6 ] );		m_pullLeptonInvPt.push_back( pull[ 7 ] );
		m_pullLeptonTheta.push_back( pull[ 8 ] );		m_pullLeptonTheta.push_back( pull[ 9 ] );
		m_pullLeptonPhi.push_back( pull[ 10 ] );		m_pullLeptonPhi.push_back( pull[ 11 ] );

		ReconstructedParticleImpl* outJet1 = dynamic_cast<ReconstructedParticleImpl*>( inputJetCollection->getElementAt( 0 ) );
		ReconstructedParticleImpl* outJet2 = dynamic_cast<ReconstructedParticleImpl*>( inputJetCollection->getElementAt( 1 ) );

		outputJetCollection->addElement( outJet1 );
		outputJetCollection->addElement( outJet2 );
		pLCEvent->addCollection( outputJetCollection , m_outputJetCollection.c_str() );
		getNormalizedResiduals( pLCEvent , outJet1 , outJet2 );
	}
	catch(DataNotAvailableException &e)
	{
		streamlog_out(MESSAGE) << "processEvent : Input collections not found in event " << m_nEvt << std::endl;
	}
	m_pTTree->Fill();

}

void ZHllqq5CFit::getNormalizedResiduals( EVENT::LCEvent *pLCEvent , ReconstructedParticleImpl* outJet1 , ReconstructedParticleImpl* outJet2 )
{
	TrueJet_Parser* trueJet	= this;
	trueJet->getall( pLCEvent );

}

void ZHllqq5CFit::getSLDsInJet( LCRelationNavigator JetSLDNav , EVENT::ReconstructedParticle* jet , pfoVector &Neutrinos )
{
	const EVENT::LCObjectVec& SLDVertices = JetSLDNav.getRelatedToObjects( jet );
	Neutrinos.clear();
	for ( unsigned int i_sld = 0 ; i_sld < SLDVertices.size() ; ++i_sld )
	{
		Vertex *sldVertex = (Vertex*) SLDVertices.at( i_sld );
		ReconstructedParticle *sldRP = (ReconstructedParticle *) ( sldVertex->getAssociatedParticle() );
		Neutrinos.push_back( sldRP->getParticles()[ 0 ] );
	}
}

void ZHllqq5CFit::getSLDCombination( std::vector<int> jetnSLDSolutions , int iteration , std::vector<int> &jetSLDCombination )
{
	int k,l,m;
	m = iteration;
	for ( unsigned int i = 0 ; i < jetnSLDSolutions.size() ; ++i )
	{
		k = m / jetnSLDSolutions[ jetnSLDSolutions.size() - i - 1 ];
		l = m - k * jetnSLDSolutions[ jetnSLDSolutions.size() - i - 1 ];
		m = k;
		jetSLDCombination[ jetnSLDSolutions.size() - i - 1 ] = l;
	}
}

int ZHllqq5CFit::performFIT( 	TLorentzVector jet1FourMomentum , std::vector<float> jet1CovMat , TLorentzVector jet2FourMomentum , std::vector<float> jet2CovMat , pfoVector leptons ,
				float &fitProbability , float (&fitOutputs)[ 18 ] , std::vector< TLorentzVector > &fittedObjects , float (&pull)[ 12 ] , bool traceEvent )
{
	JetFitObject *jet[ 2 ];
	float sigmaE , sigmaTheta , sigmaPhi;
	LeptonFitObject *lepton[ 2 ];
	int ErrorCode;
	fittedObjects.clear();
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////													//////
//////					Set JetFitObjects						//////
//////													//////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	float jet1FittedEnergy(0.0),jet1FittedTheta(0.0),jet1FittedPhi(0.0);
	float jet2FittedEnergy(0.0),jet2FittedTheta(0.0),jet2FittedPhi(0.0);
	streamlog_out(DEBUG6) << "		get jet1 resolutions"  << std::endl ;
	getJetResolutions( jet1FourMomentum , jet1CovMat , sigmaE , sigmaTheta , sigmaPhi );
	jet[ 0 ] = new JetFitObject( jet1FourMomentum.E() , jet1FourMomentum.Theta() , jet1FourMomentum.Phi() , sigmaE , sigmaTheta , sigmaPhi , jet1FourMomentum.M() );
	jet[ 0 ]->setName( "jet1" );
	streamlog_out(DEBUG8)  << " start four-vector of jet1: " << *jet[ 0 ]  << std::endl ;
	streamlog_out(DEBUG6) << "		get jet2 resolutions"  << std::endl ;
	getJetResolutions( jet2FourMomentum , jet2CovMat , sigmaE , sigmaTheta , sigmaPhi );
	jet[ 1 ] = new JetFitObject( jet2FourMomentum.E() , jet2FourMomentum.Theta() , jet2FourMomentum.Phi() , sigmaE , sigmaTheta , sigmaPhi , jet2FourMomentum.M() );
	jet[ 1 ]->setName( "jet2" );
	streamlog_out(DEBUG8)  << " start four-vector of jet2: " << *jet[ 1 ]  << std::endl ;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////													//////
//////					Set LeptonFitObjects						//////
//////													//////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	float lepton1FittedEnergy(0.0),lepton1FittedTheta(0.0),lepton1FittedPhi(0.0);
	float lepton2FittedEnergy(0.0),lepton2FittedTheta(0.0),lepton2FittedPhi(0.0);
	float parameters[ 3 ]{ 0.0 } , errors[ 3 ]{ 0.0 };
	streamlog_out(DEBUG6) << "		get lepton1 parameters"  << std::endl ;
	getLeptonParameters( leptons[ 0 ] , parameters , errors );
	lepton[ 0 ] = new LeptonFitObject ( parameters[ 0 ] , parameters[ 1 ] , parameters[ 2 ] , errors[ 0 ] , errors[ 1 ] , errors[ 2 ] , leptons[ 0 ]->getMass() );
	lepton[ 0 ]->setName( "lepton1" );
	streamlog_out(DEBUG8)  << " start four-vector of lepton1: " << *lepton[ 0 ]  << std::endl ;
	streamlog_out(DEBUG6) << "		get lepton2 parameters"  << std::endl ;
	getLeptonParameters( leptons[ 1 ] , parameters , errors );
	lepton[ 1 ] = new LeptonFitObject ( parameters[ 0 ] , parameters[ 1 ] , parameters[ 2 ] , errors[ 0 ] , errors[ 1 ] , errors[ 2 ] , leptons[ 1 ]->getMass() );
	lepton[ 1 ]->setName( "lepton2" );
	streamlog_out(DEBUG8)  << " start four-vector of lepton2: " << *lepton[ 1 ]  << std::endl ;

////	these don't get changed by the fit -> to obtain start values later!
	const int NJETS = 2;
	const int NLEPTONS = 2;

	streamlog_out(DEBUG8)  << "	<<<<<<<<<<<<<<<<<<<<<<<<<< START VALUES >>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	JetFitObject startjets[NJETS] = {*jet[0], *jet[1]};
	for (int i = 0; i < NJETS; ++i)	streamlog_out(DEBUG8)  << "	startjets[" << i << "]: " << startjets[i]  << std::endl;

	LeptonFitObject startleptons[NLEPTONS] = {*lepton[0], *lepton[1]};
	for (int i = 0; i < NLEPTONS; ++i) streamlog_out(DEBUG8)  << "	startleptons[" << i << "]: " << startleptons[i]  << std::endl;

	JetFitObject fitjets[NJETS] = {*jet[0], *jet[1]};
	for (int i = 0; i < NJETS; ++i) streamlog_out(DEBUG8)  << "	fitjets[" << i << "]: " << fitjets[i]  << std::endl ;

	LeptonFitObject fitleptons[NLEPTONS] = {*lepton[0], *lepton[1]};
	for (int i = 0; i < NLEPTONS; ++i) streamlog_out(DEBUG8)  << "	fitleptons[" << i << "]: " << fitleptons[i]  << std::endl ;

////	these get changed by the fit
	streamlog_out(DEBUG8)  << "	<<<<<<<<<<<<<<<<<<<<<<<<<< VALUES TO BE FITTED >>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	JetFitObject *JETs[NJETS];
	for (int i = 0; i < NJETS; ++i)
	{
		JETs[i] = &fitjets[i];
		streamlog_out(DEBUG8)  << "	start four-vector of jet " << i << ": " << *(JETs[i])  << std::endl ;
	}

	LeptonFitObject *LEPTONs[NLEPTONS];
	for (int i = 0; i < NLEPTONS; ++i)
	{
		LEPTONs[i] = &fitleptons[i];
		streamlog_out(DEBUG8)  << "	start four-vector of leptons " << i << ": " << *(LEPTONs[i])  << std::endl ;
	}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////													//////
//////					Set Constraints Before Fit					//////
//////													//////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	float target_p_due_crossing_angle = m_ECM * 0.007; // crossing angle = 14 mrad
	MomentumConstraint pxc ( 0 , 1 , 0 , 0 , target_p_due_crossing_angle );//Factor for: (energy sum, px sum, py sum,pz sum,target value of sum)

	pxc.setName("sum(p_x)");
	for (int i = 0; i < NJETS; ++i) pxc.addToFOList(*(JETs[i]));
	for (int i = 0; i < NLEPTONS; ++i) pxc.addToFOList(*(LEPTONs[i]));

	MomentumConstraint pyc (0, 0, 1, 0, 0);
	pyc.setName("sum(p_y)");
	for (int i = 0; i < NJETS; ++i) pyc.addToFOList(*(JETs[i]));
	for (int i = 0; i < NLEPTONS; ++i) pyc.addToFOList(*(LEPTONs[i]));

	MomentumConstraint pzc (0, 0, 0, 1, 0);
	pzc.setName("sum(p_z)");
	for (int i = 0; i < NJETS; ++i) pzc.addToFOList(*(JETs[i]));
	for (int i = 0; i < NLEPTONS; ++i) pzc.addToFOList(*(LEPTONs[i]));

	double E_lab = 2 * sqrt( std::pow( 0.548579909e-3 , 2 ) + std::pow( m_ECM / 2 , 2 ) + std::pow( target_p_due_crossing_angle , 2 ) + 0. + 0.);
	MomentumConstraint ec(1, 0, 0, 0, E_lab);
	ec.setName("sum(E)");
	for (int i = 0; i < NJETS; ++i) ec.addToFOList(*(JETs[i]));
	for (int i = 0; i < NLEPTONS; ++i) ec.addToFOList(*(LEPTONs[i]));

	streamlog_out(DEBUG8)  << "	Value of E_lab before adding ISR: " << E_lab << std::endl ;
	streamlog_out(DEBUG8)  << "	Value of target_p_due_crossing_angle before adding ISR: " << target_p_due_crossing_angle << std::endl ;
	streamlog_out(DEBUG8)  << "	Value of pxc before adding ISR: " << pxc.getValue() << std::endl ;
	streamlog_out(DEBUG8)  << "	Value of pyc before adding ISR: " << pyc.getValue() << std::endl ;
	streamlog_out(DEBUG8)  << "	Value of pzc before adding ISR: " << pzc.getValue() << std::endl ;
	streamlog_out(DEBUG8)  << "	Value of ec before adding ISR: " << ec.getValue() << std::endl ;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////													//////
//////					Set ISR PhotonFitObjects					//////
//////													//////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	ISRPhotonFitObject *photon = new ISRPhotonFitObject (0., 0., -pzc.getValue(), b, ISRPzMaxB);
	if( m_fitISR )
	{
		streamlog_out(DEBUG8)  << "start four-vector of ISR photon: " << *(photon) << std::endl ;

		pxc.addToFOList(*(photon));
		pyc.addToFOList(*(photon));
		pzc.addToFOList(*(photon));
		ec.addToFOList(*(photon));
	}

	streamlog_out(DEBUG8)  << "	Value of E_lab before fit: " << E_lab << std::endl ;
	streamlog_out(DEBUG8)  << "	Value of target_p_due_crossing_angle before fit: " << target_p_due_crossing_angle << std::endl ;
	streamlog_out(DEBUG8)  << "	Value of pxc after adding ISR before fit: " << pxc.getValue() << std::endl ;
	streamlog_out(DEBUG8)  << "	Value of pyc after adding ISR before fit: " << pyc.getValue() << std::endl ;
	streamlog_out(DEBUG8)  << "	Value of pzc after adding ISR before fit: " << pzc.getValue() << std::endl ;
	streamlog_out(DEBUG8)  << "	Value of ec after adding ISR before fit: " << ec.getValue() << std::endl ;

	SoftGaussMassConstraint z(91.2,2.4952/2);
	z.addToFOList(*(LEPTONs[0]), 1);
	z.addToFOList(*(LEPTONs[1]), 1);

	MassConstraint h(125.);
	h.addToFOList(*(JETs[0]), 1);
	h.addToFOList(*(JETs[1]), 1);

	double ZstartMass = z.getMass(1);
	double HstartMass = h.getMass(1);

	streamlog_out(DEBUG8) << "start mass of Z: " << ZstartMass << std::endl ;
	streamlog_out(DEBUG8) << "start mass of H: " << HstartMass << std::endl ;

	BaseFitter *pfitter;

	if ( m_fitter == 1 )
	{
		pfitter = new NewFitterGSL();
		if ( traceEvent ) (dynamic_cast<NewFitterGSL*>(pfitter))->setDebug( 1 );

		streamlog_out(DEBUG6) << "fit using GSL Fitter"  << std::endl ;
	}
	else if ( m_fitter == 2 )
	{
		pfitter = new NewtonFitterGSL();
		if ( traceEvent ) (dynamic_cast<NewtonFitterGSL*>(pfitter))->setDebug( 1 );

		streamlog_out(DEBUG6) << "fit using Newton GSL Fitter"  << std::endl ;
	}
	else
	{
////		OPALFitter has no method setDebug !
		pfitter = new OPALFitterGSL();

		streamlog_out(DEBUG6) << "fit using OPAL GSL Fitter"  << std::endl ;
		if ( traceEvent ) (dynamic_cast<OPALFitterGSL*>(pfitter))->setDebug( 1 );
	}
	BaseFitter &fitter = *pfitter;

	for (int i = 0; i < NJETS; ++i) fitter.addFitObject( *(JETs[i]) );
	for (int i = 0; i < NLEPTONS; ++i) fitter.addFitObject( *(LEPTONs[i]) );

	if( m_fitISR )
	{
		fitter.addFitObject( *(photon) );
		streamlog_out(DEBUG8) << "ISR added to fit"  << std::endl ;
	}

	fitter.addConstraint( pxc );
	fitter.addConstraint( pyc );
	fitter.addConstraint( pzc );
	fitter.addConstraint( ec );
	fitter.addConstraint( z );

	streamlog_out(DEBUG8) << "constraints added"  << std::endl ;

	fitProbability = fitter.fit();
	double chi2 = fitter.getChi2();
	int nit = fitter.getIterations();
	ErrorCode = fitter.getError();

	streamlog_out(DEBUG8)  << "	<<<<<<<<<<<<<<<<<<<<<<<<<< FIT OUTPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	streamlog_out(DEBUG8) << "	fit probability = " << fitProbability << std::endl ;
	streamlog_out(DEBUG8) << "	fit chi2 = " << chi2  << std::endl ;
	streamlog_out(DEBUG8) << "	error code: " << ErrorCode << std::endl ;

	for (int i = 0; i < NJETS; ++i)
	{
		streamlog_out(DEBUG8)  << "final four-vector of jet " << i << ": " << *(JETs[i]) << std::endl ;
		streamlog_out(DEBUG8)  << "final px of jet " << i << ": " << (JETs[i]) << std::endl ;
	}
	for (int i = 0; i < NLEPTONS; ++i)
	{
		streamlog_out(DEBUG8)  << "final four-vector of lepton " << i << ": " << *(LEPTONs[i]) << std::endl ;
		streamlog_out(DEBUG8)  << "final px of lepton " << i << ": " << (LEPTONs[i]) << std::endl ;
	}
	if( m_fitISR ) streamlog_out(DEBUG9)  << "final four-vector of ISR photon: " << *(photon) << std::endl;

	streamlog_out(DEBUG9)  << "fitter error: " << ErrorCode << std::endl;

	streamlog_out(DEBUG8)  << "Value of pxc after fit: " << pxc.getValue() << std::endl ;
	streamlog_out(DEBUG8)  << "Value of pyc after fit: " << pyc.getValue() << std::endl ;
	streamlog_out(DEBUG8)  << "Value of pzc after fit: " << pzc.getValue() << std::endl ;
	streamlog_out(DEBUG8)  << "Value of ec after fit: " << ec.getValue() << std::endl ;

	double pullJet[3][2];
	double pullLepton[3][2];

	float Zmass_after_fit = 0.0;
	float Hmass_after_fit = 0.0;
	ISREnergy = 0.0;
	ISRmomentum[ 0 ] = 0.0;
	ISRmomentum[ 1 ] = 0.0;
	ISRmomentum[ 2 ] = 0.0;
	ZEnergy = 0.0;
	Zmomentum[ 0 ] = 0.0;
	Zmomentum[ 1 ] = 0.0;
	Zmomentum[ 2 ] = 0.0;
	HEnergy = 0.0;
	Hmomentum[ 0 ] = 0.0;
	Hmomentum[ 1 ] = 0.0;
	Hmomentum[ 2 ] = 0.0;

	if ( ErrorCode <= 0)
	{
////		require successfull error calculation for pulls!
		if ( ErrorCode == 0)
		{
			Zmass_after_fit = z.getMass( 1 );
			Hmass_after_fit = h.getMass( 1 );
			jet1FittedEnergy = JETs[ 0 ]->getParam( 0 );
			jet1FittedTheta = JETs[ 0 ]->getParam( 1 );
			jet1FittedPhi = JETs[ 0 ]->getParam( 2 );
			jet2FittedEnergy = JETs[ 1 ]->getParam( 0 );
			jet2FittedTheta = JETs[ 1 ]->getParam( 1 );
			jet2FittedPhi = JETs[ 1 ]->getParam( 2 );
			lepton1FittedEnergy = LEPTONs[ 0 ]->getParam( 0 );
			lepton1FittedTheta = LEPTONs[ 0 ]->getParam( 1 );
			lepton1FittedPhi = LEPTONs[ 0 ]->getParam( 2 );
			lepton2FittedEnergy = LEPTONs[ 1 ]->getParam( 0 );
			lepton2FittedTheta = LEPTONs[ 1 ]->getParam( 1 );
			lepton2FittedPhi = LEPTONs[ 1 ]->getParam( 2 );
			for (int ifo = 0; ifo < 2; ifo++)
			{
				double start, fitted;
				double errfit, errmea, sigma;
				for (int ipar = 0; ipar < 3; ipar++)
				{
					fitted = JETs[ifo]->getParam(ipar);
					start = startjets[ifo].getParam(ipar);
					errfit = JETs[ifo]->getError(ipar);
					errmea = startjets[ifo].getError(ipar);
					sigma = sqrt(fabs(errmea*errmea-errfit*errfit));
//					sigma = errmea*errmea-errfit*errfit;
					if (sigma > 0)
					{
						sigma = sqrt(sigma);
						pullJet[ipar][ifo] = (fitted - start)/sigma;
					}
					else
					{
						pullJet[ipar][ifo] = -4.5;
					}
				}
			}
			for (int ifo = 0; ifo < 2; ifo++)
			{
				double start, fitted;
				double errfit, errmea, sigma;
				for (int ipar = 0; ipar < 3; ipar++)
				{
					fitted = LEPTONs[ifo]->getParam(ipar);
					start = startleptons[ifo].getParam(ipar);
					errfit = LEPTONs[ifo]->getError(ipar);
					errmea = startleptons[ifo].getError(ipar);
					sigma = sqrt(fabs(errmea*errmea-errfit*errfit));
//					sigma = errmea*errmea-errfit*errfit;
					if (sigma > 0)
					{
						sigma = sqrt(sigma);
						pullLepton[ipar][ifo] = (fitted - start)/sigma;
					}
					else
					{
						pullLepton[ipar][ifo] = -4.5;
					}
				}
			}
		}
		else //if  fitter.getError() is not = 0
		{
			streamlog_out(DEBUG6) << " ERROR CALCULATION FAILED" << std::endl ;
		}
		ISREnergy = photon->getE();
		ISRmomentum[0] = photon->getPx();
		ISRmomentum[1] = photon->getPy();
		ISRmomentum[2] = photon->getPz();
		Zmomentum[0] = LEPTONs[0]->getPx() + LEPTONs[1]->getPx();
		Zmomentum[1] = LEPTONs[0]->getPy() + LEPTONs[1]->getPy();
		Zmomentum[2] = LEPTONs[0]->getPz() + LEPTONs[1]->getPz();
		Hmomentum[0] = JETs[0]->getPx() + JETs[1]->getPx();
		Hmomentum[1] = JETs[0]->getPy() + JETs[1]->getPy();
		Hmomentum[2] = JETs[0]->getPz() + JETs[1]->getPz();
		ZEnergy = LEPTONs[0]->getE() + LEPTONs[1]->getE();
		HEnergy = JETs[0]->getE() + JETs[1]->getE();
	}//end-if fitter.getError() <=0
	else
	{
		streamlog_out(DEBUG8) << "FIT ERROR = " << fitter.getError()
					<< ", not filling histograms!"  << std::endl ;
		streamlog_out(DEBUG8)  << "start mass of Z: " << ZstartMass << std::endl ;
		streamlog_out(DEBUG8)  << "start mass of H: " << HstartMass << std::endl ;
		streamlog_out(DEBUG8)  << "final mass of Z: " << z.getMass(1) << std::endl ;
		streamlog_out(DEBUG8)  << "final mass of H: " << h.getMass(1) << std::endl ;
	}


	streamlog_out(DEBUG9)  << "start mass of Z: " << ZstartMass << std::endl ;
	streamlog_out(DEBUG9)  << "start mass of H: " << HstartMass << std::endl ;
	streamlog_out(DEBUG9)  << "mass of fitted Z: " << Zmass_after_fit << std::endl ;
	streamlog_out(DEBUG9)  << "mass of fitted H: " << Hmass_after_fit << std::endl ;
	streamlog_out(DEBUG9)  << "Error Code: " << ErrorCode << std::endl ;

	fitOutputs[ 0 ] = nit; // number of itterations
	fitOutputs[ 1 ] = chi2;
	fitOutputs[ 2 ] = ZstartMass; // Mass of Z-boson (= Invariant mass of di-lepton) before fit
	fitOutputs[ 3 ] = HstartMass; // Mass of Higgs-boson (= Invariant mass of di-jet) before fit
	fitOutputs[ 4 ] = Zmass_after_fit; // Mass of Z-boson (= Invariant mass of di-lepton) after fit
	fitOutputs[ 5 ] = Hmass_after_fit; // Mass of Higgs-boson (= Invariant mass of di-jet) after fit
	fitOutputs[ 6 ] = jet1FittedEnergy;
	fitOutputs[ 7 ] = jet1FittedTheta;
	fitOutputs[ 8 ] = jet1FittedPhi;
	fitOutputs[ 9 ] = jet2FittedEnergy;
	fitOutputs[ 10 ] = jet2FittedTheta;
	fitOutputs[ 11 ] = jet2FittedPhi;
	fitOutputs[ 12 ] = lepton1FittedEnergy;
	fitOutputs[ 13 ] = lepton1FittedTheta;
	fitOutputs[ 14 ] = lepton1FittedPhi;
	fitOutputs[ 15 ] = lepton2FittedEnergy;
	fitOutputs[ 16 ] = lepton2FittedTheta;
	fitOutputs[ 17 ] = lepton2FittedPhi;
	streamlog_out(DEBUG4)  << "Fit Output is formed" <<  std::endl;
	for ( int i = 0 ; i < 17 ; ++i )
	{
		streamlog_out(DEBUG4)  << "Fit Output[ " << i << " ] = " << fitOutputs[ i ] <<  std::endl;
	}

	pull[ 0 ] = pullJet[ 0 ][ 0 ]; // pull Energy jet 1
	pull[ 1 ] = pullJet[ 0 ][ 1 ]; // pull Energy jet 2
	pull[ 2 ] = pullJet[ 1 ][ 0 ]; // pull Theta jet 1
	pull[ 3 ] = pullJet[ 1 ][ 1 ]; // pull Theta jet 2
	pull[ 4 ] = pullJet[ 2 ][ 0 ]; // pull Phi jet 1
	pull[ 5 ] = pullJet[ 2 ][ 1 ]; // pull Phi jet 2
	pull[ 6 ] = pullLepton[ 0 ][ 0 ]; // pull InvPt lepton 1
	pull[ 7 ] = pullLepton[ 0 ][ 1 ]; // pull InvPt lepton 2
	pull[ 8 ] = pullLepton[ 1 ][ 0 ]; // pull Theta lepton 1
	pull[ 9 ] = pullLepton[ 1 ][ 1 ]; // pull Theta lepton 2
	pull[ 10 ] = pullLepton[ 2 ][ 0 ]; // pull Phi lepton 1
	pull[ 11 ] = pullLepton[ 2 ][ 1 ]; // pull Phi lepton 2
	streamlog_out(DEBUG4)  << "Fit pulls are obtained" <<  std::endl;
	for ( int i = 0 ; i < 12 ; ++i )
	{
		streamlog_out(DEBUG4)  << "pull[ " << i << " ] = " << pull[ i ] <<  std::endl;
	}

	fittedObjects.push_back( TLorentzVector( Zmomentum[ 0 ] , Zmomentum[ 1 ] , Zmomentum[ 2 ] , ZEnergy ) );
	fittedObjects.push_back( TLorentzVector( Hmomentum[ 0 ] , Hmomentum[ 1 ] , Hmomentum[ 2 ] , HEnergy ) );
	fittedObjects.push_back( TLorentzVector( ISRmomentum[ 0 ] , ISRmomentum[ 1 ] , ISRmomentum[ 2 ] , ISREnergy ) );

	streamlog_out(DEBUG4)  << "Fit Objects are created" <<  std::endl;
	for ( int i = 0 ; i < 3 ; ++i )
	{
		streamlog_out(DEBUG4)  << "Fit Object[ " << i << " ] = ( " << fittedObjects[ i ].Px() << "	, " << fittedObjects[ i ].Py() << "	, " << fittedObjects[ i ].Pz() << "	, " << fittedObjects[ i ].E() << "	)" <<  std::endl;
	}

	return ErrorCode;

}


void ZHllqq5CFit::getJetResolutions(	TLorentzVector jetFourMomentum , std::vector<float> jetCovMat , float &sigmaE , float &sigmaTheta , float &sigmaPhi )
{
	float Px , Py , Pz , P2 , Pt , Pt2;
	float dTheta_dPx , dTheta_dPy , dTheta_dPz , dPhi_dPx , dPhi_dPy;
	float sigmaPx2 , sigmaPy2 , sigmaPz2 , sigmaPxPy , sigmaPxPz , sigmaPyPz;

	Px		= jetFourMomentum.Px();
	Py		= jetFourMomentum.Py();
	Pz		= jetFourMomentum.Pz();
	P2		= ( jetFourMomentum.Vect() ).Mag2();
	Pt2		= pow( Px , 2 ) + pow( Py , 2 );
	Pt		= sqrt( Pt2 );
	sigmaPx2	= jetCovMat[ 0 ];
	sigmaPxPy	= jetCovMat[ 1 ];
	sigmaPy2	= jetCovMat[ 2 ];
	sigmaPxPz	= jetCovMat[ 3 ];
	sigmaPyPz	= jetCovMat[ 4 ];
	sigmaPz2	= jetCovMat[ 5 ];

	dTheta_dPx	= Px * Pz / ( P2 * Pt );
	dTheta_dPy	= Py * Pz / ( P2 * Pt );
	dTheta_dPz	= -Pt / P2;
	dPhi_dPx	= -Py / Pt2;
	dPhi_dPy	= Px / Pt2;

	sigmaE		= sqrt( jetCovMat[ 9 ] );
	sigmaTheta	= std::sqrt( std::fabs( std::pow( dTheta_dPx , 2 ) * sigmaPx2 + std::pow( dTheta_dPy , 2 ) * sigmaPy2 + std::pow( dTheta_dPz , 2 ) * sigmaPz2 +
	 					2.0 * dTheta_dPx * dTheta_dPy * sigmaPxPy + 2.0 * dTheta_dPx * dTheta_dPz * sigmaPxPz + 2.0 * dTheta_dPy * dTheta_dPz * sigmaPyPz ) );
	sigmaPhi	= std::sqrt( std::fabs( std::pow( dPhi_dPx , 2 ) * sigmaPx2 + std::pow( dPhi_dPy , 2 ) * sigmaPy2 + 2.0 * dPhi_dPx * dPhi_dPy * sigmaPxPy ) );
	streamlog_out(DEBUG6) << "			E		= " << jetFourMomentum.E() << std::endl ;
	streamlog_out(DEBUG6) << "			Theta		= " << jetFourMomentum.Theta() << std::endl ;
	streamlog_out(DEBUG6) << "			Phi		= " << jetFourMomentum.Phi() << std::endl ;
	streamlog_out(DEBUG6) << "			sigmaE		= " << sigmaE << std::endl ;
	streamlog_out(DEBUG6) << "			sigmaTheta	= " << sigmaTheta << std::endl ;
	streamlog_out(DEBUG6) << "			sigmaPhi	= " << sigmaPhi << std::endl ;

}

void ZHllqq5CFit::getLeptonParameters( ReconstructedParticle* lepton , float (&parameters)[ 3 ] , float (&errors)[ 3 ] )
{
	TrackVec trackVec = lepton->getTracks();
	if ( trackVec.size() != 1 )
	{
		streamlog_out(DEBUG4)  << "Number of tracks for lepton is not exactly ONE!!! (nTracks = " << trackVec.size() << " ) " << std::endl ;
		streamlog_out(DEBUG4) << *lepton << std::endl;
		TLorentzVector leptonFourMomentum( lepton->getMomentum() , lepton->getEnergy() );
		float Px		= leptonFourMomentum.Px();
		float Py		= leptonFourMomentum.Py();
		float Pz		= leptonFourMomentum.Pz();
		float P2		= ( leptonFourMomentum.Vect() ).Mag2();
		float Pt2		= std::pow( Px , 2 ) + std::pow( Py , 2 );
		float Pt		= std::sqrt( Pt2 );

		float sigmaPx2		= lepton->getCovMatrix()[ 0 ];
		float sigmaPxPy		= lepton->getCovMatrix()[ 1 ];
		float sigmaPy2		= lepton->getCovMatrix()[ 2 ];
		float sigmaPxPz		= lepton->getCovMatrix()[ 3 ];
		float sigmaPyPz		= lepton->getCovMatrix()[ 4 ];
		float sigmaPz2		= lepton->getCovMatrix()[ 5 ];

		float dInvPt_dPx	= - Px / ( Pt * Pt2 );
		float dInvPt_dPy	= - Py / ( Pt * Pt2 );
		float dTheta_dPx	= Px * Pz / ( P2 * Pt );
		float dTheta_dPy	= Py * Pz / ( P2 * Pt );
		float dTheta_dPz	= -Pt / P2;
		float dPhi_dPx		= -Py / Pt2;
		float dPhi_dPy		= Px / Pt2;

		parameters[ 0 ] = 1.0 / std::sqrt( std::pow( leptonFourMomentum.Px() , 2 ) + std::pow( leptonFourMomentum.Py() , 2 ) );
		parameters[ 1 ] = leptonFourMomentum.Theta();
		parameters[ 2 ] = leptonFourMomentum.Phi();
		errors[ 0 ]	= std::sqrt( std::pow( dInvPt_dPx , 2 ) * sigmaPx2 + std::pow( dInvPt_dPy , 2 ) * sigmaPy2 + 2.0 * dInvPt_dPx * dInvPt_dPy * sigmaPxPy );
		errors[ 1 ]	= std::sqrt( std::fabs( std::pow( dTheta_dPx , 2 ) * sigmaPx2 + std::pow( dTheta_dPy , 2 ) * sigmaPy2 + std::pow( dTheta_dPz , 2 ) * sigmaPz2 +
		 					2.0 * dTheta_dPx * dTheta_dPy * sigmaPxPy + 2.0 * dTheta_dPx * dTheta_dPz * sigmaPxPz + 2.0 * dTheta_dPy * dTheta_dPz * sigmaPyPz ) );
		errors[ 2 ]	= std::sqrt( std::fabs( std::pow( dPhi_dPx , 2 ) * sigmaPx2 + std::pow( dPhi_dPy , 2 ) * sigmaPy2 + 2.0 * dPhi_dPx * dPhi_dPy * sigmaPxPy ) );
	}
	else
	{
		streamlog_out(DEBUG4)  << "	Lepton has exactly ONE track:" << std::endl ;
		streamlog_out(DEBUG4) << *lepton << std::endl;
		streamlog_out(DEBUG4) << *trackVec[ 0 ] << std::endl;
		float Omega		= trackVec[ 0 ]->getOmega();
		float tanLambda		= trackVec[ 0 ]->getTanLambda();
		float Theta		= 2.0 * atan( 1.0 ) - atan( tanLambda );//atan( 1.0 / tanLambda );
		float Phi		= trackVec[ 0 ]->getPhi();

		float sigmaOmega	= std::sqrt( trackVec[ 0 ]->getCovMatrix()[ 5 ] );
		float sigmaTanLambda	= std::sqrt( trackVec[ 0 ]->getCovMatrix()[ 14 ] );
		float sigmaPhi		= std::sqrt( trackVec[ 0 ]->getCovMatrix()[ 2 ] );

		float dTheta_dTanLambda	= -1.0 / ( 1.0 + std::pow( tanLambda , 2 ) );

		parameters[ 0 ]	= Omega / eB;
		parameters[ 1 ]	= Theta;
		parameters[ 2 ]	= Phi;
		errors[ 0 ]	= sigmaOmega / eB;
		errors[ 1 ]	= std::fabs( dTheta_dTanLambda ) * sigmaTanLambda;
		errors[ 2 ]	= sigmaPhi;
	}
	streamlog_out(DEBUG6) << "			Inverse pT	= " << parameters[ 0 ] << std::endl ;
	streamlog_out(DEBUG6) << "			Theta		= " << parameters[ 1 ] << std::endl ;
	streamlog_out(DEBUG6) << "			Phi		= " << parameters[ 2 ] << std::endl ;
	streamlog_out(DEBUG6) << "			SigmaInverse pT	= " << errors[ 0 ] << std::endl ;
	streamlog_out(DEBUG6) << "			SigmaTheta	= " << errors[ 1 ] << std::endl ;
	streamlog_out(DEBUG6) << "			SigmaPhi	= " << errors[ 2 ] << std::endl ;
}

/*
std::vector<std::vector<float>> ZHllqq5CFit::performOldFIT(EVENT::LCEvent *pLCEvent, TLorentzVector Jet0_Nutlv, TLorentzVector Jet1_Nutlv , std::vector<float> nu1CovMat , std::vector<float> nu2CovMat)
{
//	Setvalues();
	std::vector<std::vector<float>> FitResult{};
	std::vector<float> fitStartValues;
	std::vector<float> fitOutputs;
	std::vector<float> fittedParticles;
	std::vector<float> pulls;
	std::vector<float> constraints;
	std::vector<float> uncertainties;
	std::vector<float> diJetSystem;
	LCCollection *inputErrorFlowCollection = pLCEvent->getCollection( m_inputJetCollection );
	LCCollection *inputLeptonCollection = pLCEvent->getCollection( m_inputIsolatedlaptonCollection );

	const int nJets = 2;
	const int nLeptons = 2;
	double Omega = 0.;
	double Omega_uncert = 0.;
	double TanLambda = 0.;
	double TanLambda_err = 0.;
	double theta = 0.;
	double theta_err = 0.;
	double phi = 0.;
	double phi_err = 0.;
	double invers_pT = 0.;
	double invers_pT_err = 0.;

	float pxc_before_ISR = 0.;
	float pyc_before_ISR = 0.;
	float pzc_before_ISR = 0.;
	float ec_before_ISR = 0.;
	float zc_before_ISR = 0.;

	float pxc_before_fit = 0.;
	float pyc_before_fit = 0.;
	float pzc_before_fit = 0.;
	float ec_before_fit = 0.;
	float zc_before_fit = 0.;

	float pxc_after_fit = 0.;
	float pyc_after_fit = 0.;
	float pzc_after_fit = 0.;
	float ec_after_fit = 0.;
	float zc_after_fit = 0.;

	float Px,Py,Pz,Px2,Py2,Pz2,Pt,Pt2,P,P2;
	float SigPx2,SigPxSigPy,SigPy2,SigPxSigPz,SigPySigPz,SigPz2,SigE2;
	float dth_dpx,dth_dpy,dth_dpz,dphi_dpx,dphi_dpy;


	float jet0_Px,jet0_Py,jet0_Pz,jet0_E,jet0_Theta,jet0_Phi;
	float jet1_Px,jet1_Py,jet1_Pz,jet1_E,jet1_Theta,jet1_Phi;
	float fittedjet0_E,fittedjet0_Theta,fittedjet0_Phi;
	float fittedjet1_E,fittedjet1_Theta,fittedjet1_Phi;
	float lepton0_Px,lepton0_Py,lepton0_Pz,lepton0_E;
	float lepton1_Px,lepton1_Py,lepton1_Pz,lepton1_E;
	float jet0_SigmaTheta,jet0_SigmaPhi,jet0_SigmaE;
	float jet1_SigmaTheta,jet1_SigmaPhi,jet1_SigmaE;
	float lepton0_SigmaTheta,lepton0_SigmaPhi,lepton0_SigmaInvpT;
	float lepton1_SigmaTheta,lepton1_SigmaPhi,lepton1_SigmaInvpT;

	float JetResE,JetResTheta,JetResPhi;
	float sigmaScaleFactor = 1.0;
	float m_JetResThetaCoeff = 1.0;
	float m_JetResPhiCoeff = 1.0;


	JetFitObject *jet[nJets];
	LeptonFitObject *lepton[nLeptons];
	int nSingularCovMatrix = 0;
	HepLorentzVector jetvec;
	HepLorentzVector Nuvec;
	HepLorentzVector leptonvec;
	Hep3Vector jet0_unit;
	Hep3Vector jet1_unit;

	float bestzvalue = 10000.;
	int besterr = 999;
	float bestprob = 0.;
	int bestnit = 0;
	float beststartmassZ = 0., beststartmassH = 0.;
	float startmassZ = 0., startmassH = 0.;
	float bestphotonenergy = 0.;
	float chi2startmassZ = 0.;
	float chi2startmassH = 0.;
	memset(Zmomentum, 0, sizeof(Zmomentum));
	memset(Hmomentum, 0, sizeof(Hmomentum));
	memset(ISRmomentum, 0, sizeof(ISRmomentum));
	float bestmassZ = 0.0, bestmassH = 0.0;
	float Z_Energy=0.;
	float H_Energy=0.;
	float chi2best=0.;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////					Set JetFitObjects
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for (int i_jet = 0; i_jet < nJets; i_jet++)
	{
		ReconstructedParticle *j;
		j = dynamic_cast<ReconstructedParticle*>( inputErrorFlowCollection->getElementAt( i_jet ) );

		if ( i_jet == 0 ) Nuvec = HepLorentzVector( Jet0_Nutlv.Px() , Jet0_Nutlv.Py() , Jet0_Nutlv.Pz() , Jet0_Nutlv.E() );
		if ( i_jet == 1 ) Nuvec = HepLorentzVector( Jet1_Nutlv.Px() , Jet1_Nutlv.Py() , Jet1_Nutlv.Pz() , Jet1_Nutlv.E() );

		Px =	j->getMomentum()[0] + Nuvec.px();
		Px2 =	std::pow( Px , 2 );
		Py =	j->getMomentum()[1] + Nuvec.py();
		Py2 =	std::pow( Py , 2 );
		Pz =	j->getMomentum()[2] + Nuvec.pz();
		Pz2 =	std::pow( Pz , 2 );
		Pt2 =	Px2 + Py2;
		Pt =	std::sqrt( Pt2 );
		P =	std::sqrt( Px2 + Py2 + Pz2 );
		P2 =	std::pow( P , 2 );

		TVector3 startJetMomentum( Px , Py , Pz );
		if ( i_jet == 0 )
		{
			jet0_Px = Px;
			jet0_Py = Py;
			jet0_Pz = Pz;
			jet0_E = j->getEnergy() + Nuvec.e();
			jet0_Theta = startJetMomentum.Theta();
			jet0_Phi = startJetMomentum.Phi();
			SigPx2 =	j->getCovMatrix()[0] + nu1CovMat[ 0 ];
			SigPxSigPy =	j->getCovMatrix()[1] + nu1CovMat[ 1 ];
			SigPy2 =	j->getCovMatrix()[2] + nu1CovMat[ 2 ];
			SigPxSigPz =	j->getCovMatrix()[3] + nu1CovMat[ 3 ];
			SigPySigPz =	j->getCovMatrix()[4] + nu1CovMat[ 4 ];
			SigPz2 =	j->getCovMatrix()[5] + nu1CovMat[ 5 ];
			SigE2 =		j->getCovMatrix()[9] + nu1CovMat[ 9 ];
		}
		else if ( i_jet == 1 )
		{
			jet1_Px = Px;
			jet1_Py = Py;
			jet1_Pz = Pz;
			jet1_E = j->getEnergy();
			jet1_Theta = startJetMomentum.Theta();
			jet1_Phi = startJetMomentum.Phi();
			SigPx2 =	j->getCovMatrix()[0] + nu2CovMat[ 0 ];
			SigPxSigPy =	j->getCovMatrix()[1] + nu2CovMat[ 1 ];
			SigPy2 =	j->getCovMatrix()[2] + nu2CovMat[ 2 ];
			SigPxSigPz =	j->getCovMatrix()[3] + nu2CovMat[ 3 ];
			SigPySigPz =	j->getCovMatrix()[4] + nu2CovMat[ 4 ];
			SigPz2 =	j->getCovMatrix()[5] + nu2CovMat[ 5 ];
			SigE2 =		j->getCovMatrix()[9] + nu2CovMat[ 9 ];
		}


		dth_dpx =	Px * Pz / ( P2 * Pt );
		dth_dpy =	Py * Pz / ( P2 * Pt );
		dth_dpz =	-Pt / P2;

		dphi_dpx =	-Py / Pt2;
		dphi_dpy =	Px / Pt2;

		JetResE =	std::sqrt( SigE2 ) * sigmaScaleFactor;
		JetResTheta =	m_JetResThetaCoeff * std::sqrt( std::fabs( SigPx2 * std::pow( dth_dpx , 2 ) + SigPy2 * std::pow( dth_dpy , 2 ) + SigPz2 * std::pow( dth_dpz , 2 ) + 2 * ( SigPxSigPy * dth_dpx * dth_dpy ) + 2 * ( SigPySigPz * dth_dpy * dth_dpz ) + 2 * ( SigPxSigPz * dth_dpx * dth_dpz ) ) );
		JetResPhi =	m_JetResPhiCoeff * std::sqrt( std::fabs( SigPx2 * std::pow( dphi_dpx , 2 ) + SigPy2 * std::pow( dphi_dpy , 2 ) + 2 * ( SigPxSigPy * dphi_dpx * dphi_dpy ) ) );
		streamlog_out(MESSAGE) << "JET[" << i_jet << "]: JetResEnergy = " << JetResE << " , JetResTheta = " << JetResTheta << " , JetResPhi = " << JetResPhi << std::endl;

		if ( SigPx2 == 0. || SigPxSigPy == 0. || SigPxSigPz == 0. || SigPy2 == 0. || SigPySigPz == 0. || SigPz2 == 0. || SigE2 == 0. ) // Check CovMatrix singularity
		{
			streamlog_out(WARNING) << "Covariance Matrix is singular"<<std::endl;
			streamlog_out(WARNING) << "Setting theta and phi Resolution back to default values "<<std::endl;
			JetResTheta = 0.1;//m_jetThetaError;
			JetResPhi = 0.1;//m_jetPhiError;
			nSingularCovMatrix++;
		}
		jetvec = HepLorentzVector ( ( j->getMomentum() )[0] , ( j->getMomentum() )[1] , ( j->getMomentum() )[2] , j->getEnergy() );
		jetvec += Nuvec;
		if ( i_jet == 0 )
		{
			jet0_SigmaTheta = JetResTheta;
			jet0_SigmaPhi = JetResPhi;
			jet0_SigmaE = JetResE;
			jet0_unit = (Hep3Vector(jetvec)).unit();
		}
		else if ( i_jet == 1 )
		{
			jet1_SigmaTheta = JetResTheta;
			jet1_SigmaPhi = JetResPhi;
			jet1_SigmaE = JetResE;
			jet1_unit = (Hep3Vector(jetvec)).unit();
		}

		jet[i_jet] = new JetFitObject ( jetvec.e(), jetvec.theta() , jetvec.phi() , JetResE , JetResTheta , JetResPhi , jetvec.m() );
		streamlog_out(DEBUG)  << " start four-vector of jet[" << i_jet << "]: " << *jet[i_jet]  << std::endl ;
		if ( i_jet == 0 )
		{
			jet[i_jet]->setName("Jet0");
		}
		else if (i_jet == 1)
		{
			jet[i_jet]->setName("Jet1");
		}
		diJetSystem.push_back(jetvec.e());
		diJetSystem.push_back(jetvec.theta());
		diJetSystem.push_back(jetvec.phi());
		diJetSystem.push_back(JetResE);
		diJetSystem.push_back(JetResTheta);
		diJetSystem.push_back(JetResPhi);
	}
	diJetSystem.push_back( acos( jet0_unit.dot(jet1_unit) ) * 45 / atan(1) );

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////					Set LeptonFitObjects
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for(int i_lep = 0; i_lep < nLeptons ; i_lep++)
	{
		ReconstructedParticle* l = dynamic_cast<ReconstructedParticle*>( inputLeptonCollection->getElementAt( i_lep ) ) ;
		leptonvec = HepLorentzVector ((l->getMomentum())[0],(l->getMomentum())[1],(l->getMomentum())[2],l->getEnergy());
		TrackVec tckvec = l->getTracks();
		if ( tckvec.size() != 1 )
		{
			streamlog_out(DEBUG)  << "Number of tracks for lepton[" << i_lep <<"] is not exactly ONE!!! (nTracks = " << tckvec.size() << " ) " << std::endl ;
			invers_pT = 1/(std::sqrt(pow(leptonvec.px(),2)+pow(leptonvec.py(),2)));
			invers_pT_err = 2*std::sqrt(pow(leptonvec.px(),2)+pow(leptonvec.py(),2))*0.00001;
			theta = leptonvec.theta();
			theta_err = 0.1;//m_jetThetaError;
			phi = leptonvec.phi();
			phi_err = 0.1;//m_jetPhiError;
		}
		else
		{
			streamlog_out(DEBUG)  << "Number of tracks for lepton[" << i_lep <<"] is exactly ONE!!!" << std::endl;
			Omega = tckvec[0]->getOmega();
			Omega_uncert = std::sqrt( std::abs(tckvec[0]->getCovMatrix()[5]) );
			streamlog_out(DEBUG)  << "Omega_uncert = " << Omega_uncert << std::endl;
			TanLambda = tckvec[0]->getTanLambda();
			TanLambda_err = std::sqrt( std::abs(tckvec[0]->getCovMatrix()[14]) );
			streamlog_out(DEBUG)  << "TanLambda_err = " << TanLambda_err << std::endl;
			theta = 2 * atan( 1. ) - atan( TanLambda );
			theta_err = TanLambda_err / ( 1 + pow( TanLambda , 2 ) );
			streamlog_out(DEBUG)  << "theta_err = " << theta_err << std::endl;
			phi = tckvec[0]->getPhi();
			phi_err = std::sqrt( std::abs(tckvec[0]->getCovMatrix()[2]) );
			streamlog_out(DEBUG)  << "phi_err = " << phi_err << std::endl;
			invers_pT = Omega / eB;
			invers_pT_err = std::fabs( 1. / eB ) * Omega_uncert;
			streamlog_out(DEBUG)  << "invers_pT_err = " << invers_pT_err << std::endl;
			if ( i_lep == 0 )
			{
				lepton0_Px = cos( phi ) / invers_pT;
				lepton0_Py = sin( phi ) / invers_pT;
				lepton0_Pz = TanLambda / invers_pT;
				lepton0_E = l->getEnergy();
			}
			else if ( i_lep == 1 )
			{
				lepton1_Px = cos( phi ) / invers_pT;
				lepton1_Py = sin( phi ) / invers_pT;
				lepton1_Pz = TanLambda / invers_pT;
				lepton1_E = l->getEnergy();
			}

		}
		if ( i_lep == 0 )
		{
			lepton0_SigmaTheta = theta_err;
			lepton0_SigmaPhi = phi_err;
			lepton0_SigmaInvpT = invers_pT_err;
		}
		else if ( i_lep == 1 )
		{
			lepton1_SigmaTheta = theta_err;
			lepton1_SigmaPhi = phi_err;
			lepton1_SigmaInvpT = invers_pT_err;
		}
		streamlog_out(DEBUG)  << "Lepton fit object from leptonvec: "
			<< 1/(std::sqrt(pow(leptonvec.px(),2)+pow(leptonvec.py(),2))) <<" +- " << 2*std::sqrt(pow(leptonvec.px(),2)+pow(leptonvec.py(),2))*0.00001 << " , "
			<< leptonvec.theta() <<" +- " << 0.1 m_jetThetaError << " , "
			<< leptonvec.phi() <<" +- " << 0.1 m_jetPhiError << std::endl ;

		streamlog_out(DEBUG)  << "Lepton fit object from track:     "
			<< std::fabs( tckvec[0]->getOmega() / eB ) <<" +- " << std::fabs( 1. / eB ) * std::sqrt( tckvec[0]->getCovMatrix()[5] ) << " , "
			<< 2 * atan( 1. ) - atan( tckvec[0]->getTanLambda() ) <<" +- " << std::abs( std::sqrt( tckvec[0]->getCovMatrix()[14]) ) / ( 1 + pow( tckvec[0]->getTanLambda() , 2 ) ) << " , "
			<< tckvec[0]->getPhi() <<" +- " << std::abs( std::sqrt( tckvec[0]->getCovMatrix()[2] ) ) << std::endl ;

		lepton[i_lep] = new LeptonFitObject (invers_pT , theta , phi , invers_pT_err , theta_err , phi_err, leptonvec.m());
		if (i_lep == 0 )
		{
			lepton[i_lep]->setName("Lepton0");
		}
		else if (i_lep == 1)
		{
			lepton[i_lep]->setName("Lepton1");
		}
		streamlog_out(DEBUG)  << " start four-vector of lepton[" << i_lep <<"]: " << *lepton[i_lep]  << std::endl ;
	}

	const int NJETS = 2;
	streamlog_out(DEBUG) << "*jet[0]: " << *jet[0] << std::endl ;
	streamlog_out(DEBUG) << "*jet[1]: " << *jet[1] << std::endl ;

	const int NLEPTONS = 2;
	streamlog_out(MESSAGE) << "*lepton[0]: " << *lepton[0] << std::endl ;
	streamlog_out(MESSAGE) << "*lepton[1]: " << *lepton[1] << std::endl ;

////	these don't get changed by the fit -> to obtain start values later!
	JetFitObject startjets[NJETS] = {*jet[0], *jet[1]};
	for (int i = 0; i < NJETS; ++i)	streamlog_out(DEBUG)  << "startjets[" << i << "]: " << startjets[i]  << std::endl;

	LeptonFitObject startleptons[NLEPTONS] = {*lepton[0], *lepton[1]};
	for (int i = 0; i < NLEPTONS; ++i) streamlog_out(DEBUG)  << "startleptons[" << i << "]: " << startleptons[i]  << std::endl;

	JetFitObject fitjets[NJETS] = {*jet[0], *jet[1]};
	for (int i = 0; i < NJETS; ++i)
		streamlog_out(DEBUG)  << "fitjets[" << i << "]: " << fitjets[i]  << std::endl ;

	LeptonFitObject fitleptons[NLEPTONS] = {*lepton[0], *lepton[1]};
	for (int i = 0; i < NLEPTONS; ++i)
		streamlog_out(DEBUG)  << "fitleptons[" << i << "]: " << fitleptons[i]  << std::endl ;

////	these get changed by the fit
	JetFitObject *JETs[NJETS];
	for (int i = 0; i < NJETS; ++i)
	{
		JETs[i] = &fitjets[i];
		streamlog_out(MESSAGE)  << "start four-vector of jet " << i << ": " << *(JETs[i])  << std::endl ;
	}

	LeptonFitObject *leptons[NLEPTONS];
	for (int i = 0; i < NLEPTONS; ++i)
	{
		leptons[i] = &fitleptons[i];
		streamlog_out(MESSAGE)  << "start four-vector of leptons " << i << ": " << *(leptons[i])  << std::endl ;
	}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////				set constraints befor fit
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	float target_p_due_crossing_angle = m_ECM * 0.007; // crossing angle = 14 mrad
	MomentumConstraint pxc ( 0 , 1 , 0 , 0 , target_p_due_crossing_angle );//Factor for: (energy sum, px sum, py sum,pz sum,target value of sum)

	pxc.setName("sum(p_x)");
	for (int i = 0; i < NJETS; ++i) pxc.addToFOList(*(JETs[i]));
	for (int i = 0; i < NLEPTONS; ++i) pxc.addToFOList(*(leptons[i]));

	MomentumConstraint pyc (0, 0, 1, 0, 0);
	pyc.setName("sum(p_y)");
	for (int i = 0; i < NJETS; ++i) pyc.addToFOList(*(JETs[i]));
	for (int i = 0; i < NLEPTONS; ++i) pyc.addToFOList(*(leptons[i]));

	MomentumConstraint pzc (0, 0, 0, 1, 0);
	pzc.setName("sum(p_z)");
	for (int i = 0; i < NJETS; ++i) pzc.addToFOList(*(JETs[i]));
	for (int i = 0; i < NLEPTONS; ++i) pzc.addToFOList(*(leptons[i]));

	float E_lab= 2 * sqrt( std::pow( 0.548579909e-3 , 2 ) + std::pow( m_ECM / 2 , 2 ) + std::pow( target_p_due_crossing_angle , 2 ) + 0. + 0.);
	MomentumConstraint ec(1, 0, 0, 0, E_lab);
	ec.setName("sum(E)");
	for (int i = 0; i < NJETS; ++i) ec.addToFOList(*(JETs[i]));
	for (int i = 0; i < NLEPTONS; ++i) ec.addToFOList(*(leptons[i]));

	streamlog_out(DEBUG)  << "Value of E_lab before adding ISR: " << E_lab << std::endl ;
	streamlog_out(DEBUG)  << "Value of target_p_due_crossing_angle before adding ISR: " << target_p_due_crossing_angle << std::endl ;
	streamlog_out(DEBUG)  << "Value of pxc before adding ISR: " << pxc.getValue() << std::endl ;
	streamlog_out(DEBUG)  << "Value of pyc before adding ISR: " << pyc.getValue() << std::endl ;
	streamlog_out(DEBUG)  << "Value of pzc before adding ISR: " << pzc.getValue() << std::endl ;
	streamlog_out(DEBUG)  << "Value of ec before adding ISR: " << ec.getValue() << std::endl ;
	pxc_before_ISR = pxc.getValue();
	pyc_before_ISR = pyc.getValue();
	pzc_before_ISR = pzc.getValue();
	ec_before_ISR = ec.getValue();

////	ISR Photon initialized with missing p_z
	ISRPhotonFitObject *photon = new ISRPhotonFitObject (0., 0., -pzc.getValue(), b, ISRPzMaxB);
	float ISRstartPx = photon->getPx();
	float ISRstartPy = photon->getPy();
	float ISRstartPz = photon->getPz();
////	ISRPhotonFitObject(double px, double py, double pz,
////	double b_, double PzMaxB_, double PzMinB_ = 0.);
	if(m_fitISR)
	{
		streamlog_out(MESSAGE)  << "start four-vector of ISR photon: " << *(photon) << std::endl ;

		pxc.addToFOList(*(photon));
		pyc.addToFOList(*(photon));
		pzc.addToFOList(*(photon));
		ec.addToFOList(*(photon));
	}

	streamlog_out(DEBUG)  << "Value of E_lab before fit: " << E_lab << std::endl ;
	streamlog_out(DEBUG)  << "Value of target_p_due_crossing_angle before fit: " << target_p_due_crossing_angle << std::endl ;
	streamlog_out(DEBUG)  << "Value of pxc after adding ISR before fit: " << pxc.getValue() << std::endl ;
	streamlog_out(DEBUG)  << "Value of pyc after adding ISR before fit: " << pyc.getValue() << std::endl ;
	streamlog_out(DEBUG)  << "Value of pzc after adding ISR before fit: " << pzc.getValue() << std::endl ;
	streamlog_out(DEBUG)  << "Value of ec after adding ISR before fit: " << ec.getValue() << std::endl ;
	pxc_before_fit = pxc.getValue();
	pyc_before_fit = pyc.getValue();
	pzc_before_fit = pzc.getValue();
	ec_before_fit = ec.getValue();

	SoftGaussMassConstraint z(91.2,2.4952/2);
	z.addToFOList(*(leptons[0]), 1);
	z.addToFOList(*(leptons[1]), 1);
	zc_before_ISR = z.getValue();
	zc_before_fit = z.getValue();

	MassConstraint h(125.);
	h.addToFOList(*(JETs[0]), 1);
	h.addToFOList(*(JETs[1]), 1);

	startmassZ = z.getMass(1);
	startmassH = h.getMass(1);

	streamlog_out(MESSAGE) << "start mass of Z: " << startmassZ << std::endl ;
	streamlog_out(MESSAGE) << "start mass of H: " << startmassH << std::endl ;

	Hmass_NoFit = startmassH;

	BaseFitter *pfitter;

	int debug = 0;
	if ( pLCEvent->getEventNumber() == m_ievttrace || m_traceall) debug = 10;

	if (m_fitter == 1)
	{
		pfitter = new NewFitterGSL();
		if (pLCEvent->getEventNumber()== m_ievttrace || m_traceall) (dynamic_cast<NewFitterGSL*>(pfitter))->setDebug(debug);

		streamlog_out(DEBUG) << "ifitter is 1"  << std::endl ;
	}
	else if (m_fitter == 2)
	{
		pfitter = new NewtonFitterGSL();
		if (pLCEvent->getEventNumber()== m_ievttrace || m_traceall) (dynamic_cast<NewtonFitterGSL*>(pfitter))->setDebug(debug);

		streamlog_out(DEBUG) << "ifitter is 2"  << std::endl ;
	}
	else
	{
////		OPALFitter has no method setDebug !
		pfitter = new OPALFitterGSL();

		streamlog_out(DEBUG) << "ifitter is not 1 or 2"  << std::endl ;
		if (pLCEvent->getEventNumber()== m_ievttrace || m_traceall) (dynamic_cast<OPALFitterGSL*>(pfitter))->setDebug(debug);
	}
	BaseFitter &fitter = *pfitter;

	TextTracer tracer (std::cout);
	if (pLCEvent->getEventNumber() == m_ievttrace || m_traceall) fitter.setTracer(tracer);
	for (int i = 0; i < NJETS; ++i) fitter.addFitObject(*(JETs[i]));
	for (int i = 0; i < NLEPTONS; ++i) fitter.addFitObject(*(leptons[i]));

	if(m_fitISR)
	{
		fitter.addFitObject (*(photon));
		streamlog_out(DEBUG) << "ISR added to fit"  << std::endl ;
	}

	fitter.addConstraint(pxc);
	fitter.addConstraint(pyc);
	fitter.addConstraint(pzc);
	fitter.addConstraint(ec);
	fitter.addConstraint(z);

	streamlog_out(DEBUG) << "constraints added"  << std::endl ;

	if (fabs(startmassZ-91.2) + fabs(startmassH-125.) < bestzvalue)
	{
		chi2startmassZ = startmassZ;
		chi2startmassH = startmassH;
		bestzvalue = fabs(startmassZ-91.2) + fabs(startmassH-125.);

		streamlog_out(DEBUG) << "best z value is this..." <<  bestzvalue << std::endl ;
	}

	float prob = fitter.fit();
	double chi2 = fitter.getChi2();
	int nit = fitter.getIterations();

	streamlog_out(DEBUG) << "fit probability = " << prob << std::endl ;
	streamlog_out(DEBUG) << "fit chi2 = " << chi2  << std::endl ;
	streamlog_out(DEBUG) << "error code: " << fitter.getError() << std::endl ;

	for (int i = 0; i < NJETS; ++i)
	{
		streamlog_out(MESSAGE)  << "final four-vector of jet " << i << ": " << *(JETs[i]) << std::endl ;
		streamlog_out(DEBUG)  << "final px of jet " << i << ": " << (JETs[i]) << std::endl ;
	}
	for (int i = 0; i < NLEPTONS; ++i)
	{
		streamlog_out(MESSAGE)  << "final four-vector of lepton " << i << ": " << *(leptons[i]) << std::endl ;
		streamlog_out(DEBUG)  << "final px of lepton " << i << ": " << (leptons[i]) << std::endl ;
	}
	if(m_fitISR) streamlog_out(MESSAGE)  << "final four-vector of ISR photon: " << *(photon) << std::endl;
	int ierr = fitter.getError();
	streamlog_out(MESSAGE)  << "fitter error: " << ierr << std::endl;
	if ((besterr > 0 && ierr < besterr) || ( besterr < 0 && ierr == 0)) besterr = ierr;

	streamlog_out(DEBUG)  << "Value of pxc after fit: " << pxc.getValue() << std::endl ;
	streamlog_out(DEBUG)  << "Value of pyc after fit: " << pyc.getValue() << std::endl ;
	streamlog_out(DEBUG)  << "Value of pzc after fit: " << pzc.getValue() << std::endl ;
	streamlog_out(DEBUG)  << "Value of ec after fit: " << ec.getValue() << std::endl ;
	pxc_after_fit = pxc.getValue();
	pyc_after_fit = pyc.getValue();
	pzc_after_fit = pzc.getValue();
	ec_after_fit = ec.getValue();
	zc_after_fit = z.getValue();


	if (ierr <= 0)
	{
		double pullJet[3][2];
		double pullLepton[3][2];
////		require successfull error calculation for pulls!
		if (ierr == 0)
		{
			for (int ifo = 0; ifo < 2; ifo++)
			{
				double start, fitted;
				double errfit, errmea, sigma;
				for (int ipar = 0; ipar < 3; ipar++)
				{
					fitted = JETs[ifo]->getParam(ipar);
					start = startjets[ifo].getParam(ipar);
					errfit = JETs[ifo]->getError(ipar);
					errmea = startjets[ifo].getError(ipar);
					sigma = errmea*errmea-errfit*errfit;
					if (sigma > 0)
					{
						sigma = sqrt(sigma);
						pullJet[ipar][ifo] = (fitted - start)/sigma;
					}
					else
					{
						pullJet[ipar][ifo] = -4.5;
					}
				}
				if ( ifo == 0 )
				{
					fittedjet0_E = JETs[ifo]->getParam(0);
					fittedjet0_Theta = JETs[ifo]->getParam(1);
					fittedjet0_Phi = JETs[ifo]->getParam(2);
				}
				else
				{
					fittedjet1_E = JETs[ifo]->getParam(0);
					fittedjet1_Theta = JETs[ifo]->getParam(1);
					fittedjet1_Phi = JETs[ifo]->getParam(2);
				}

			}
			for (int ifo = 0; ifo < 2; ifo++)
			{
				double start, fitted;
				double errfit, errmea, sigma;
				for (int ipar = 0; ipar < 3; ipar++)
				{
					fitted = leptons[ifo]->getParam(ipar);
					start = startleptons[ifo].getParam(ipar);
					errfit = leptons[ifo]->getError(ipar);
					errmea = startleptons[ifo].getError(ipar);
					sigma = errmea*errmea-errfit*errfit;
					if (sigma > 0)
					{
						sigma = sqrt(sigma);
						pullLepton[ipar][ifo] = (fitted - start)/sigma;
					}
					else
					{
						pullLepton[ipar][ifo] = -4.5;
					}
				}
			}
		}
		if (prob >= bestprob)
		{
			bestprob = prob;
			streamlog_out(DEBUG)  << "BESTPROB: " << bestprob << std::endl ;
			bestnit  = nit;
			bestmassZ = z.getMass(1);
			bestmassH = h.getMass(1);
			beststartmassZ = startmassZ;
			beststartmassH = startmassH;
			bestphotonenergy = photon->getE();
			ISRmomentum[0] = photon->getPx();
			ISRmomentum[1] = photon->getPy();
			ISRmomentum[2] = photon->getPz();
			Zmomentum[0] = leptons[0]->getPx() + leptons[1]->getPx();
			Zmomentum[1] = leptons[0]->getPy() + leptons[1]->getPy();
			Zmomentum[2] = leptons[0]->getPz() + leptons[1]->getPz();
			Hmomentum[0] = JETs[0]->getPx() + JETs[1]->getPx();
			Hmomentum[1] = JETs[0]->getPy() + JETs[1]->getPy();
			Hmomentum[2] = JETs[0]->getPz() + JETs[1]->getPz();
			Z_Energy = leptons[0]->getE() + leptons[1]->getE();
			H_Energy = JETs[0]->getE() + JETs[1]->getE();
			chi2best = fitter.getChi2();
			errorcode = fitter.getError();
			if (ierr == 0) //if  fitter.getError() is = 0
			{
				hpull_jet_E=pullJet[0][0];
				hpull_jet2_E=pullJet[0][1];
				hpull_jet_th=pullJet[1][0];
				hpull_jet2_th=pullJet[1][1];
				hpull_jet_phi=pullJet[2][0];
				hpull_jet2_phi=pullJet[2][1];
				hpull_lepton_InvpT=pullLepton[0][0];
				hpull_lepton2_InvpT=pullLepton[0][1];
				hpull_lepton_th=pullLepton[1][0];
				hpull_lepton2_th=pullLepton[1][1];
				hpull_lepton_phi=pullLepton[2][0];
				hpull_lepton2_phi=pullLepton[2][1];
			}//end if  fitter.getError() is = 0
			else //if  fitter.getError() is not = 0
			{
				streamlog_out(DEBUG) << " ERROR CALCULATION FAILED for best permutation in event " << pLCEvent->getEventNumber() << std::endl ;
			}
		}
	}//end-if fitter.getError() <=0
	else
	{
		streamlog_out(DEBUG) << "FIT ERROR = " << fitter.getError()
					<< " in event " << pLCEvent->getEventNumber()
					<< ", not filling histograms!"  << std::endl ;
		streamlog_out(DEBUG)  << "start mass of Z: " << startmassZ << std::endl ;
		streamlog_out(DEBUG)  << "start mass of H: " << startmassH << std::endl ;
		streamlog_out(DEBUG)  << "final mass of Z: " << z.getMass(1) << std::endl ;
		streamlog_out(DEBUG)  << "final mass of H: " << h.getMass(1) << std::endl ;
	}

	float Zmass_before_fit=beststartmassZ;
	float Zmass_after_fit=bestmassZ;
	float Hmass_before_fit=beststartmassH;
	float Hmass_after_fit=bestmassH;
	Error_code=errorcode;

	streamlog_out(DEBUG)  << "min chi2 start mass of Z: " << chi2startmassZ << std::endl ;
	streamlog_out(DEBUG)  << "min chi2 start mass of H: " << chi2startmassH << std::endl ;
	streamlog_out(DEBUG)  << "best start mass of Z: " << beststartmassZ << std::endl ;
	streamlog_out(DEBUG)  << "best start mass of H: " << beststartmassH << std::endl ;
	streamlog_out(DEBUG)  << "best mass of Z: " << bestmassZ << std::endl ;
	streamlog_out(DEBUG)  << "best mass of H: " << bestmassH << std::endl ;
	streamlog_out(DEBUG)  << "Error Code: " << errorcode << std::endl ;

	fitStartValues.push_back(jet0_Px);
	fitStartValues.push_back(jet1_Px);
	fitStartValues.push_back(jet0_Py);
	fitStartValues.push_back(jet1_Py);
	fitStartValues.push_back(jet0_Pz);
	fitStartValues.push_back(jet1_Pz);
	fitStartValues.push_back(jet0_E);
	fitStartValues.push_back(jet1_E);
	fitStartValues.push_back(lepton0_Px);
	fitStartValues.push_back(lepton1_Px);
	fitStartValues.push_back(lepton0_Py);
	fitStartValues.push_back(lepton1_Py);
	fitStartValues.push_back(lepton0_Pz);
	fitStartValues.push_back(lepton1_Pz);
	fitStartValues.push_back(lepton0_E);
	fitStartValues.push_back(lepton1_E);
	fitStartValues.push_back(ISRstartPx);
	fitStartValues.push_back(ISRstartPy);
	fitStartValues.push_back(ISRstartPz);
	fitStartValues.push_back(jet0_Theta);
	fitStartValues.push_back(jet1_Theta);
	fitStartValues.push_back(jet0_Phi);
	fitStartValues.push_back(jet1_Phi);

	fitOutputs.push_back(ierr);
	fitOutputs.push_back(prob);
	fitOutputs.push_back(nit);
	fitOutputs.push_back(startmassZ);
	fitOutputs.push_back(startmassH);
	fitOutputs.push_back(beststartmassZ);
	fitOutputs.push_back(beststartmassH);
	fitOutputs.push_back(Zmass_after_fit);
	fitOutputs.push_back(Hmass_after_fit);
	fitOutputs.push_back(chi2startmassZ);
	fitOutputs.push_back(chi2startmassH);
	fitOutputs.push_back(chi2best);
	fitOutputs.push_back(bestphotonenergy);
	fitOutputs.push_back(fittedjet0_E);
	fitOutputs.push_back(fittedjet1_E);
	fitOutputs.push_back(fittedjet0_Theta);
	fitOutputs.push_back(fittedjet1_Theta);
	fitOutputs.push_back(fittedjet0_Phi);
	fitOutputs.push_back(fittedjet1_Phi);

	fittedParticles.push_back(ISRmomentum[0]);
	fittedParticles.push_back(ISRmomentum[1]);
	fittedParticles.push_back(ISRmomentum[2]);
	fittedParticles.push_back(Zmomentum[0]);
	fittedParticles.push_back(Zmomentum[1]);
	fittedParticles.push_back(Zmomentum[2]);
	fittedParticles.push_back(Z_Energy);
	fittedParticles.push_back(Hmomentum[0]);
	fittedParticles.push_back(Hmomentum[1]);
	fittedParticles.push_back(Hmomentum[2]);
	fittedParticles.push_back(H_Energy);

	pulls.push_back(hpull_jet_E);
	pulls.push_back(hpull_jet2_E);
	pulls.push_back(hpull_jet_th);
	pulls.push_back(hpull_jet2_th);
	pulls.push_back(hpull_jet_phi);
	pulls.push_back(hpull_jet2_phi);
	pulls.push_back(hpull_lepton_InvpT);
	pulls.push_back(hpull_lepton2_InvpT);
	pulls.push_back(hpull_lepton_th);
	pulls.push_back(hpull_lepton2_th);
	pulls.push_back(hpull_lepton_phi);
	pulls.push_back(hpull_lepton2_phi);

	constraints.push_back(pxc_before_ISR);
	constraints.push_back(pyc_before_ISR);
	constraints.push_back(pzc_before_ISR);
	constraints.push_back(ec_before_ISR);
	constraints.push_back(pxc_before_fit);
	constraints.push_back(pyc_before_fit);
	constraints.push_back(pzc_before_fit);
	constraints.push_back(ec_before_fit);
	constraints.push_back(pxc_after_fit);
	constraints.push_back(pyc_after_fit);
	constraints.push_back(pzc_after_fit);
	constraints.push_back(ec_after_fit);
	constraints.push_back(zc_before_ISR);
	constraints.push_back(zc_before_fit);
	constraints.push_back(zc_after_fit);

	uncertainties.push_back(jet0_SigmaTheta);
	uncertainties.push_back(jet1_SigmaTheta);
	uncertainties.push_back(jet0_SigmaPhi);
	uncertainties.push_back(jet1_SigmaPhi);
	uncertainties.push_back(jet0_SigmaE);
	uncertainties.push_back(jet1_SigmaE);
	uncertainties.push_back(lepton0_SigmaTheta);
	uncertainties.push_back(lepton1_SigmaTheta);
	uncertainties.push_back(lepton0_SigmaPhi);
	uncertainties.push_back(lepton1_SigmaPhi);
	uncertainties.push_back(lepton0_SigmaInvpT);
	uncertainties.push_back(lepton1_SigmaInvpT);

	FitResult.push_back(fitStartValues);
	FitResult.push_back(fitOutputs);
	FitResult.push_back(fittedParticles);
	FitResult.push_back(pulls);
	FitResult.push_back(constraints);
	FitResult.push_back(uncertainties);
	FitResult.push_back(diJetSystem);
	streamlog_out(DEBUG)  << "FitResult returned to event processor successfully: " << std::endl ;


	delete photon;
	return FitResult;
}
*/
void ZHllqq5CFit::check( LCEvent* )
{
//	nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ZHllqq5CFit::end()
{
//	streamlog_out(MESSAGE) << "# of events: " << m_nEvt << std::endl;
//	streamlog_out(ERROR) << "# of nucorrection: " << correction<< std::endl;
//	streamlog_out(ERROR) << "# of Covariance failed: " << nCo<< std::endl;

	m_pTFile->cd();
	m_pTTree->Write();
	m_pTFile->Close();
	delete m_pTFile;

}
