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

	registerInputCollection(	LCIO::LCRELATION,
					"SLDNeutrinoRelationCollection",
					"Name of the JetSemiLeptonicDecayVertex-Neutrinos Relation collection",
					m_inputSLDNuLink,
					std::string("SLDNuLinkName")
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

	registerInputCollection( 	LCIO::MCPARTICLE,
					"MCParticleCollection" ,
					"Name of the MCParticle collection"  ,
					_MCParticleColllectionName ,
					std::string("MCParticlesSkimmed")
					);

	registerInputCollection( 	LCIO::RECONSTRUCTEDPARTICLE,
					"RecoParticleCollection" ,
					"Name of the ReconstructedParticles input collection"  ,
					_recoParticleCollectionName ,
					std::string("PandoraPFOs")
				);

	registerInputCollection( 	LCIO::LCRELATION,
					"RecoMCTruthLink",
					"Name of the RecoMCTruthLink input collection"  ,
					_recoMCTruthLink,
					std::string("RecoMCTruthLink")
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

	registerProcessorParameter(	"SigmaEnergyScaleFactor" ,
					"Factor for scaling up energy error",
					m_SigmaEnergyScaleFactor,
					float(1.0f)
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

	registerProcessorParameter(	"matchTrueJetWithAngle" ,
					"Matching true jet with reco jet: TRUE = jets with closest angle are matched , FALSE = jets containing same leading particle are matched",
					m_matchTrueJetWithAngle,
					(bool)false
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
	m_pTTree->Branch("nCorrectedSLD",&m_nCorrectedSLD,"nCorrectedSLD/I") ;
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
	m_pTTree->Branch("Sigma_Px2",&m_Sigma_Px2);
	m_pTTree->Branch("Sigma_PxPy",&m_Sigma_PxPy);
	m_pTTree->Branch("Sigma_Py2",&m_Sigma_Py2);
	m_pTTree->Branch("Sigma_PxPz",&m_Sigma_PxPz);
	m_pTTree->Branch("Sigma_PyPz",&m_Sigma_PyPz);
	m_pTTree->Branch("Sigma_Pz2",&m_Sigma_Pz2);
	m_pTTree->Branch("Sigma_PxE",&m_Sigma_PxE);
	m_pTTree->Branch("Sigma_PyE",&m_Sigma_PyE);
	m_pTTree->Branch("Sigma_PzE",&m_Sigma_PzE);
	m_pTTree->Branch("Sigma_E2",&m_Sigma_E2);

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
	m_nCorrectedSLD = 0;
	m_FitErrorCode_woNu = -100;
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
	m_FitErrorCode_wNu = -100;
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
	m_FitErrorCode = -100;
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
	m_Sigma_Px2.clear();
	m_Sigma_PxPy.clear();
	m_Sigma_Py2.clear();
	m_Sigma_PxPz.clear();
	m_Sigma_PyPz.clear();
	m_Sigma_Pz2.clear();
	m_Sigma_PxE.clear();
	m_Sigma_PyE.clear();
	m_Sigma_PzE.clear();
	m_Sigma_E2.clear();

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
		if ( m_nJets == m_nAskedJets && m_nIsoLeps == m_nAskedIsoLeps )//&& m_nCorrectedSLD == m_nSLDecayTotal )
		{
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
			m_Sigma_Px2.push_back( jet1initialCovMat[ 0 ] );
			m_Sigma_PxPy.push_back( jet1initialCovMat[ 1 ] );
			m_Sigma_Py2.push_back( jet1initialCovMat[ 2 ] );
			m_Sigma_PxPz.push_back( jet1initialCovMat[ 3 ] );
			m_Sigma_PyPz.push_back( jet1initialCovMat[ 4 ] );
			m_Sigma_Pz2.push_back( jet1initialCovMat[ 5 ] );
			m_Sigma_PxE.push_back( jet1initialCovMat[ 6 ] );
			m_Sigma_PyE.push_back( jet1initialCovMat[ 7 ] );
			m_Sigma_PzE.push_back( jet1initialCovMat[ 8 ] );
			m_Sigma_E2.push_back( jet1initialCovMat[ 9 ] );
			ReconstructedParticle* jet2 = dynamic_cast<ReconstructedParticle*>( inputJetCollection->getElementAt( 1 ) );
			TLorentzVector jet2tlv( jet2->getMomentum() , jet2->getEnergy() );
			std::vector< float > jet2initialCovMat( 10 , 0.0 );
			for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element )
			{
				jet2initialCovMat[ i_Element ] = jet2->getCovMatrix()[ i_Element ];
			}
			m_Sigma_Px2.push_back( jet2initialCovMat[ 0 ] );
			m_Sigma_PxPy.push_back( jet2initialCovMat[ 1 ] );
			m_Sigma_Py2.push_back( jet2initialCovMat[ 2 ] );
			m_Sigma_PxPz.push_back( jet2initialCovMat[ 3 ] );
			m_Sigma_PyPz.push_back( jet2initialCovMat[ 4 ] );
			m_Sigma_Pz2.push_back( jet2initialCovMat[ 5 ] );
			m_Sigma_PxE.push_back( jet2initialCovMat[ 6 ] );
			m_Sigma_PyE.push_back( jet2initialCovMat[ 7 ] );
			m_Sigma_PzE.push_back( jet2initialCovMat[ 8 ] );
			m_Sigma_E2.push_back( jet2initialCovMat[ 9 ] );

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
			double jetNormalizedResiduals[ 6 ]{ 0.0 };

			TLorentzVector Nu1tlv( 0.0 , 0.0 , 0.0 , 0.0 );
			TLorentzVector Nu2tlv( 0.0 , 0.0 , 0.0 , 0.0 );
			std::vector< float > nu1CovMat( 10 , 0.0 );
			std::vector< float > nu2CovMat( 10 , 0.0 );
			streamlog_out(MESSAGE) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl ;
			streamlog_out(MESSAGE) << "	||||||||||||||||||||||||||||  KINFIT WITHOUT NEUTRINO COORECTION  ||||||||||||||||||||||||||||" << std::endl ;
			streamlog_out(MESSAGE) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl ;

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

			streamlog_out(MESSAGE) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl ;
			streamlog_out(MESSAGE) << "	||||||||||||||||||||||||||||  FEED NEUTRINO COORECTION TO KINFIT  ||||||||||||||||||||||||||||" << std::endl ;
			streamlog_out(MESSAGE) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl ;
			LCRelationNavigator JetSLDNav( pLCEvent->getCollection( m_inputJetSLDLink ) );
			LCRelationNavigator SLDNuNav( pLCEvent->getCollection( m_inputSLDNuLink ) );
			
			pfoVectorVector jet1NuSolutions;
			std::vector<std::vector<int>> jet1SLDCombinations;
			pfoVectorVector jet2NuSolutions;
			std::vector<std::vector<int>> jet2SLDCombinations;
			getNeutrinosInJet( JetSLDNav , SLDNuNav , jet1 , jet1NuSolutions , jet1SLDCombinations );
			getNeutrinosInJet( JetSLDNav , SLDNuNav , jet2 , jet2NuSolutions , jet2SLDCombinations );
			int njet1SLDSolutions = jet1SLDCombinations.size();
			int njet2SLDSolutions = jet2SLDCombinations.size();
			TLorentzVector jet1FourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
			std::vector< float > jet1CovMat( 10 , 0.0 );
			TLorentzVector jet2FourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
			std::vector< float > jet2CovMat( 10 , 0.0 );
			std::vector<int> jet1FinalNuSolutions( jet1NuSolutions.size() , 0 );
			std::vector<int> jet2FinalNuSolutions( jet2NuSolutions.size() , 0 );

			for ( int i_jet1 = 0 ; i_jet1 < njet1SLDSolutions ; ++i_jet1 )
			{
				streamlog_out(DEBUG6) << "" << std::endl ;
				jet1FourMomentum = jet1tlv;
				for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element ) jet1CovMat[ i_Element ] = jet1initialCovMat[ i_Element ];
				TLorentzVector Nu1FourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
				std::vector< float > Nu1CovMat( 10 , 0.0 );
				for ( unsigned int i_sld1 = 0 ; i_sld1 < jet1SLDCombinations[ i_jet1 ].size() ; ++i_sld1 )
				{
					ReconstructedParticle* Neutrino1 = jet1NuSolutions[ i_sld1 ][ jet1SLDCombinations[ i_jet1 ][ i_sld1 ] ];
					streamlog_out(DEBUG1) << *Neutrino1 << std::endl;
					Nu1FourMomentum += TLorentzVector( Neutrino1->getMomentum() , Neutrino1->getEnergy() );
					for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element ) Nu1CovMat[ i_Element ] += Neutrino1->getCovMatrix()[ i_Element ];
				}
				jet1FourMomentum += Nu1FourMomentum;
				for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element ) jet1CovMat[ i_Element ] += Nu1CovMat[ i_Element ];
				for ( int i_jet2 = 0 ; i_jet2 < njet2SLDSolutions ; ++i_jet2 )
				{
					streamlog_out(DEBUG6) << "" << std::endl ;
					jet2FourMomentum = jet2tlv;
					for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element ) jet2CovMat[ i_Element ] = jet2initialCovMat[ i_Element ];
					TLorentzVector Nu2FourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
					std::vector< float > Nu2CovMat( 10 , 0.0 );
					for ( unsigned int i_sld2 = 0 ; i_sld2 < jet2SLDCombinations[ i_jet2 ].size() ; ++i_sld2 )
					{
						ReconstructedParticle* Neutrino2 = jet2NuSolutions[ i_sld2 ][ jet2SLDCombinations[ i_jet2 ][ i_sld2 ] ];
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
							for ( unsigned int i = 0 ; i < sizeof( fitOutputs_wNu ) / sizeof( fitOutputs_wNu[ 0 ] ) ; ++i ) fitOutputs_wNu[ i ] = fitOutputs_temp[ i ];
							for ( unsigned int i = 0 ; i < sizeof( pull_wNu ) / sizeof( pull_wNu[ 0 ] ) ; ++i ) pull_wNu[ i ] = pull_temp[ i ];
							for ( unsigned int i = 0 ; i < fittedObjects_wNu.size()  ; ++i ) fittedObjects_wNu.push_back( fittedObjects_temp[ i ] );
							for ( unsigned int i_sld1 = 0 ; i_sld1 < jet1SLDCombinations[ i_jet1 ].size() ; ++i_sld1 ) jet1FinalNuSolutions[ i_sld1 ] = jet1SLDCombinations[ i_jet1 ][ i_sld1 ];
							for ( unsigned int i_sld2 = 0 ; i_sld2 < jet2SLDCombinations[ i_jet2 ].size() ; ++i_sld2 ) jet2FinalNuSolutions[ i_sld2 ] = jet2SLDCombinations[ i_jet2 ][ i_sld2 ];
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
			streamlog_out(DEBUG0) << " Checking Fit Probability of KinFit with and without Neutrino Correction " << std::endl;
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
					jet1Neutrinos.push_back( ( EVENT::ReconstructedParticle* )jet1NuSolutions[ i ][ jet1FinalNuSolutions[ i ] ] );
				}
				for ( unsigned int i = 0 ; i < jet2FinalNuSolutions.size() ; ++i )
				{
					jet2Neutrinos.push_back( ( EVENT::ReconstructedParticle* )jet2NuSolutions[ i ][ jet2FinalNuSolutions[ i ] ] );
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
			streamlog_out(DEBUG0) << " Fit Probability of KinFit with and without Neutrino Correction was checked " << std::endl;
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
			streamlog_out(DEBUG0) << " Pulls of KinFit was filled " << std::endl;

			bool foundTrueJet_woNu = true;
			bool foundTrueJet_wNu = true;
			ReconstructedParticleImpl* outJet1 = dynamic_cast<ReconstructedParticleImpl*>( inputJetCollection->getElementAt( 0 ) );
			ReconstructedParticleImpl* outJet2 = dynamic_cast<ReconstructedParticleImpl*>( inputJetCollection->getElementAt( 1 ) );
			std::vector< float > outJet1CovMat( 10 , 0.0 ); for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element ) outJet1CovMat[ i_Element ] = outJet1->getCovMatrix()[ i_Element ];
			std::vector< float > outJet2CovMat( 10 , 0.0 ); for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element ) outJet2CovMat[ i_Element ] = outJet2->getCovMatrix()[ i_Element ];
			streamlog_out(DEBUG0) << " Output jets are obtained from input jet collection " << std::endl;
			getNormalizedResiduals( pLCEvent , outJet1 , outJet2 , jetNormalizedResiduals , foundTrueJet_woNu );
			streamlog_out(DEBUG0) << " Normalized Residuals of Jets without Neutrino Correction calculated " << std::endl;
			if ( foundTrueJet_woNu )
			{
				m_normalizedResidualJetEnergy_woNu.push_back( jetNormalizedResiduals[ 0 ] );	m_normalizedResidualJetEnergy_woNu.push_back( jetNormalizedResiduals[ 1 ] );
				m_normalizedResidualJetTheta_woNu.push_back( jetNormalizedResiduals[ 2 ] );	m_normalizedResidualJetTheta_woNu.push_back( jetNormalizedResiduals[ 3 ] );
				m_normalizedResidualJetPhi_woNu.push_back( jetNormalizedResiduals[ 4 ] );	m_normalizedResidualJetPhi_woNu.push_back( jetNormalizedResiduals[ 5 ] );
				streamlog_out(DEBUG0) << " Normalized Residuals of Jets without Neutrino Correction was filled " << std::endl;
				if ( m_FitProbability_wNu > m_FitProbability_woNu )
				{
					streamlog_out(DEBUG0) << " Normalized Residuals of Jets without Neutrino Correction are BEST " << std::endl;
					m_normalizedResidualJetEnergy.push_back( jetNormalizedResiduals[ 0 ] );	m_normalizedResidualJetEnergy.push_back( jetNormalizedResiduals[ 1 ] );
					m_normalizedResidualJetTheta.push_back( jetNormalizedResiduals[ 2 ] );	m_normalizedResidualJetTheta.push_back( jetNormalizedResiduals[ 3 ] );
					m_normalizedResidualJetPhi.push_back( jetNormalizedResiduals[ 4 ] );	m_normalizedResidualJetPhi.push_back( jetNormalizedResiduals[ 5 ] );
				}
			}

			streamlog_out(DEBUG0) << " Adding Neutrino Solutions to Jets " << std::endl;
			double jet1Momentum[ 3 ]{ outJet1->getMomentum()[ 0 ] , outJet1->getMomentum()[ 1 ] , outJet1->getMomentum()[ 2 ] };
			double jet1Energy = outJet1->getEnergy();
			for ( unsigned int i_nu1 = 0 ; i_nu1 < jet1Neutrinos.size() ; ++i_nu1 )
			{
				streamlog_out(DEBUG0) << " Adding Neutrino Solution " << i_nu1 << " to Jet1 " << std::endl;
				outJet1->addParticle( jet1Neutrinos[ i_nu1 ] );
				jet1Momentum[ 0 ] += jet1Neutrinos[ i_nu1 ]->getMomentum()[ 0 ];
				jet1Momentum[ 1 ] += jet1Neutrinos[ i_nu1 ]->getMomentum()[ 1 ];
				jet1Momentum[ 2 ] += jet1Neutrinos[ i_nu1 ]->getMomentum()[ 2 ];
				jet1Energy += jet1Neutrinos[ i_nu1 ]->getEnergy();
				for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element ) outJet1CovMat[ i_Element ] += jet1Neutrinos[ i_nu1 ]->getCovMatrix()[ i_Element ];
			}
			outJet1->setMomentum( jet1Momentum );
			outJet1->setEnergy( jet1Energy );
			outJet1->setCovMatrix( outJet1CovMat );
			streamlog_out(DEBUG0) << " FourMomentum of Jet1 updated" << std::endl;

			double jet2Momentum[ 3 ]{ outJet2->getMomentum()[ 0 ] , outJet2->getMomentum()[ 1 ] , outJet2->getMomentum()[ 2 ] };
			double jet2Energy = outJet2->getEnergy();
			for ( unsigned int i_nu2 = 0 ; i_nu2 < jet2Neutrinos.size() ; ++i_nu2 )
			{
				streamlog_out(DEBUG0) << " Adding Neutrino Solution " << i_nu2 << " to Jet2 " << std::endl;
				outJet2->addParticle( jet2Neutrinos[ i_nu2 ] );
				jet2Momentum[ 0 ] += jet2Neutrinos[ i_nu2 ]->getMomentum()[ 0 ];
				jet2Momentum[ 1 ] += jet2Neutrinos[ i_nu2 ]->getMomentum()[ 1 ];
				jet2Momentum[ 2 ] += jet2Neutrinos[ i_nu2 ]->getMomentum()[ 2 ];
				jet2Energy += jet2Neutrinos[ i_nu2 ]->getEnergy();
				for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element ) outJet2CovMat[ i_Element ] += jet2Neutrinos[ i_nu2 ]->getCovMatrix()[ i_Element ];
			}
			outJet2->setMomentum( jet2Momentum );
			outJet2->setEnergy( jet2Energy );
			outJet2->setCovMatrix( outJet2CovMat );
			streamlog_out(DEBUG0) << " FourMomentum of Jet2 updated" << std::endl;
			getNormalizedResiduals( pLCEvent , outJet1 , outJet2 , jetNormalizedResiduals , foundTrueJet_wNu );
			if ( foundTrueJet_wNu )
			{
				m_normalizedResidualJetEnergy_wNu.push_back( jetNormalizedResiduals[ 0 ] );	m_normalizedResidualJetEnergy_wNu.push_back( jetNormalizedResiduals[ 1 ] );
				m_normalizedResidualJetTheta_wNu.push_back( jetNormalizedResiduals[ 2 ] );	m_normalizedResidualJetTheta_wNu.push_back( jetNormalizedResiduals[ 3 ] );
				m_normalizedResidualJetPhi_wNu.push_back( jetNormalizedResiduals[ 4 ] );	m_normalizedResidualJetPhi_wNu.push_back( jetNormalizedResiduals[ 5 ] );
				streamlog_out(DEBUG0) << " Normalized Residuals of Jets with Neutrino Correction was filled " << std::endl;
				if ( m_FitProbability_wNu <= m_FitProbability_woNu )
				{
					streamlog_out(DEBUG0) << " Normalized Residuals of Jets with Neutrino Correction are BEST " << std::endl;
					m_normalizedResidualJetEnergy.push_back( jetNormalizedResiduals[ 0 ] );	m_normalizedResidualJetEnergy.push_back( jetNormalizedResiduals[ 1 ] );
					m_normalizedResidualJetTheta.push_back( jetNormalizedResiduals[ 2 ] );	m_normalizedResidualJetTheta.push_back( jetNormalizedResiduals[ 3 ] );
					m_normalizedResidualJetPhi.push_back( jetNormalizedResiduals[ 4 ] );	m_normalizedResidualJetPhi.push_back( jetNormalizedResiduals[ 5 ] );
				}
			}

			outputJetCollection->addElement( outJet1 );
			streamlog_out(DEBUG0) << " Jet1 Added to Output Jet collection" << std::endl;
			outputJetCollection->addElement( outJet2 );
			streamlog_out(DEBUG0) << " Jet2 Added to Output Jet collection" << std::endl;
			pLCEvent->addCollection( outputJetCollection , m_outputJetCollection.c_str() );
			streamlog_out(DEBUG0) << " Output Jet collection added to event" << std::endl;
		}
	}
	catch(DataNotAvailableException &e)
	{
		streamlog_out(MESSAGE) << "processEvent : Input collections not found in event " << m_nEvt << std::endl;
	}
	m_pTTree->Fill();

}

void ZHllqq5CFit::getNormalizedResiduals( EVENT::LCEvent *pLCEvent , ReconstructedParticleImpl* recoJet1 , ReconstructedParticleImpl* recoJet2 , double (&jetNormalizedResiduals)[ 6 ] , bool &foundTrueJets )
{
	TrueJet_Parser* trueJet	= this;
	trueJet->getall( pLCEvent );

	ReconstructedParticle* leadingParticleJet1 = NULL;
	ReconstructedParticle* leadingParticleJet2 = NULL;
	float leadingEnergyJet1 = 0.0;
	float leadingEnergyJet2 = 0.0;
	streamlog_out(DEBUG0) << " looking for leading particle in jet 1 with " << ( recoJet1->getParticles() ).size() << " particles" << std::endl;
	for ( unsigned int i_par = 0 ; i_par < ( recoJet1->getParticles() ).size() ; ++i_par )
	{
		ReconstructedParticle* pfo = ( ReconstructedParticle* )recoJet1->getParticles()[ i_par ];
		streamlog_out(DEBUG0) << " checking particle " << i_par << " with Energy = " << pfo->getEnergy() << std::endl;
		streamlog_out(DEBUG3) << *pfo << std::endl;
		if ( abs( pfo->getType() ) == 12 || abs( pfo->getType() ) == 14 || abs( pfo->getType() ) == 16 ) continue;
		if ( pfo->getEnergy() > leadingEnergyJet1 )
		{
			leadingParticleJet1 = pfo;
			leadingEnergyJet1 = pfo->getEnergy();
		}
	}
	streamlog_out(DEBUG0) << " looking for leading particle in jet 2 with " << ( recoJet2->getParticles() ).size() << " particles" << std::endl;
	for ( unsigned int i_par = 0 ; i_par < ( recoJet2->getParticles() ).size() ; ++i_par )
	{
		ReconstructedParticle* pfo = ( ReconstructedParticle* )recoJet2->getParticles()[ i_par ];
		streamlog_out(DEBUG0) << " checking particle " << i_par << " with Energy = " << pfo->getEnergy() << std::endl;
		streamlog_out(DEBUG3) << *pfo << std::endl;
		if ( abs( pfo->getType() ) == 12 || abs( pfo->getType() ) == 14 || abs( pfo->getType() ) == 16 ) continue;
		if ( pfo->getEnergy() > leadingEnergyJet2 )
		{
			leadingParticleJet2 = pfo;
			leadingEnergyJet2 = pfo->getEnergy();
			streamlog_out(DEBUG0) << " So far, the energy of leading particle is: " << leadingEnergyJet2 << " GeV" << std::endl;
		}
	}
	if ( leadingParticleJet1 == NULL || leadingParticleJet2 == NULL )
	{
		foundTrueJets = false;
		return;
	}
	streamlog_out(DEBUG0) << "**************************** Leading particle of jet1 ****************************" << std::endl;
	streamlog_out(DEBUG0) << *leadingParticleJet1 << std::endl;
	streamlog_out(DEBUG0) << "**************************** Leading particle of jet2 ****************************" << std::endl;
	streamlog_out(DEBUG0) << *leadingParticleJet2 << std::endl;
	LCObjectVec jet1vec = reltjreco->getRelatedFromObjects( leadingParticleJet1 );
	streamlog_out(DEBUG0) << jet1vec.size() << " true Jet found for leading particle of jet1" << std::endl;
	LCObjectVec jet2vec = reltjreco->getRelatedFromObjects( leadingParticleJet2 );
	streamlog_out(DEBUG0) << jet2vec.size() << " true Jet found for leading particle of jet2" << std::endl;
	if ( jet1vec.size() == 0 || jet2vec.size() == 0 )
	{
		foundTrueJets = false;
		return;
	}

	int trueJet1_index;
	int trueJet2_index;
/*
	int nTrueHadronicJets = 0;
	std::vector<int> trueHadronicJetIndices; trueHadronicJetIndices.clear();
	for (int i_jet = 0 ; i_jet < trueJet->njets() ; i_jet++ )
	{
		if ( type_jet( i_jet ) == 1 )
		{
			++nTrueHadronicJets;
			trueHadronicJetIndices.push_back( i_jet );
		}
	}
	if ( m_matchTrueJetWithAngle && nTrueHadronicJets != 2 )
	{
		foundTrueJets = false;
		return;
	}
	TVector3 jet1SeenMomentum( pseen( trueHadronicJetIndices[ 0 ] )[ 0 ] , pseen( trueHadronicJetIndices[ 0 ] )[ 1 ] , pseen( trueHadronicJetIndices[ 0 ] )[ 2 ] ); jet1SeenMomentum.SetMag( 1.0 );
	TVector3 jet2SeenMomentum( pseen( trueHadronicJetIndices[ 1 ] )[ 0 ] , pseen( trueHadronicJetIndices[ 1 ] )[ 1 ] , pseen( trueHadronicJetIndices[ 1 ] )[ 2 ] ); jet2SeenMomentum.SetMag( 1.0 );
	TVector3 jet1RecMomentum( recoJet1->getMomentum() ); jet1RecMomentum.SetMag( 1.0 );
	TVector3 jet2RecMomentum( recoJet2->getMomentum() ); jet2RecMomentum.SetMag( 1.0 );
	if ( m_matchTrueJetWithAngle )
	{
		if ( jet1RecMomentum.Dot( jet1SeenMomentum ) + jet2RecMomentum.Dot( jet2SeenMomentum ) > jet1RecMomentum.Dot( jet2SeenMomentum ) + jet2RecMomentum.Dot( jet1SeenMomentum ) )
		{
			trueJet1_index = trueHadronicJetIndices[ 0 ];
			trueJet2_index = trueHadronicJetIndices[ 1 ];
		}
		else
		{
			trueJet1_index = trueHadronicJetIndices[ 1 ];
			trueJet2_index = trueHadronicJetIndices[ 0 ];
		}
	}
	else
	{
	}
*/
	trueJet1_index = jetindex( dynamic_cast<ReconstructedParticle*>( jet1vec[ 0 ] ) );
	trueJet2_index = jetindex( dynamic_cast<ReconstructedParticle*>( jet2vec[ 0 ] ) );

	streamlog_out(DEBUG0) << " true Jet[ " << trueJet1_index << " ] has the leading particle of jet1" << std::endl;
	streamlog_out(DEBUG0) << " true Jet[ " << trueJet2_index << " ] has the leading particle of jet2" << std::endl;
	streamlog_out(DEBUG4) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG4) << "	trueJet1 TYPE:	" << type_jet( trueJet1_index ) << std::endl;
	streamlog_out(DEBUG4) << "	trueJet1 (Px,Py,Pz,E):	" << ptrue( trueJet1_index )[ 0 ] << " , " << ptrue( trueJet1_index )[ 1 ] << " , " << ptrue( trueJet1_index )[ 2 ] << " , " << Etrue( trueJet1_index ) << std::endl;
	streamlog_out(DEBUG4) << "recoJet1:" << std::endl;
	streamlog_out(DEBUG4) << *recoJet1 << std::endl;
	streamlog_out(DEBUG4) << "leadingParticleJet1:" << std::endl;
	streamlog_out(DEBUG4) << *leadingParticleJet1 << std::endl;
	streamlog_out(DEBUG4) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG4) << "	trueJet2 TYPE:	" << type_jet( trueJet2_index ) << std::endl;
	streamlog_out(DEBUG4) << "	trueJet2 (Px,Py,Pz,E):	" << ptrue( trueJet2_index )[ 0 ] << " , " << ptrue( trueJet2_index )[ 1 ] << " , " << ptrue( trueJet2_index )[ 2 ] << " , " << Etrue( trueJet2_index ) << std::endl;
	streamlog_out(DEBUG4) << "recoJet2:" << std::endl;
	streamlog_out(DEBUG4) << *recoJet2 << std::endl;
	streamlog_out(DEBUG4) << "leadingParticleJet2:" << std::endl;
	streamlog_out(DEBUG4) << *leadingParticleJet2 << std::endl;

	TVector3 jet1TrueMomentum( ptrue( trueJet1_index )[ 0 ] , ptrue( trueJet1_index )[ 1 ] , ptrue( trueJet1_index )[ 2 ] );
	double jet1TrueEnergy = Etrue( trueJet1_index );
	TVector3 jet1RecoMomentum( recoJet1->getMomentum() );
	double jet1RecoEnergy = recoJet1->getEnergy();

	double jet1EnergyResidual , jet1ThetaResidual , jet1PhiResidual;
	getJetResiduals( jet1TrueMomentum , jet1TrueEnergy , jet1RecoMomentum , jet1RecoEnergy , jet1EnergyResidual , jet1ThetaResidual , jet1PhiResidual );
	double jet1SigmaE , jet1SigmaTheta , jet1SigmaPhi;
	getJetResolutions( TLorentzVector( jet1RecoMomentum , jet1RecoEnergy ) , recoJet1->getCovMatrix() , jet1SigmaE , jet1SigmaTheta , jet1SigmaPhi );

	TVector3 jet2TrueMomentum( ptrue( trueJet2_index )[ 0 ] , ptrue( trueJet2_index )[ 1 ] , ptrue( trueJet2_index )[ 2 ] );
	double jet2TrueEnergy = Etrue( trueJet2_index );
	TVector3 jet2RecoMomentum( recoJet2->getMomentum() );
	double jet2RecoEnergy = recoJet2->getEnergy();

	double jet2EnergyResidual , jet2ThetaResidual , jet2PhiResidual;
	getJetResiduals( jet2TrueMomentum , jet2TrueEnergy , jet2RecoMomentum , jet2RecoEnergy , jet2EnergyResidual , jet2ThetaResidual , jet2PhiResidual );
	double jet2SigmaE , jet2SigmaTheta , jet2SigmaPhi;
	getJetResolutions( TLorentzVector( jet2RecoMomentum , jet2RecoEnergy ) , recoJet2->getCovMatrix() , jet2SigmaE , jet2SigmaTheta , jet2SigmaPhi );

	jetNormalizedResiduals[ 0 ] = jet1EnergyResidual / jet1SigmaE;
	jetNormalizedResiduals[ 1 ] = jet2EnergyResidual / jet2SigmaE;
	jetNormalizedResiduals[ 2 ] = jet1ThetaResidual / jet1SigmaTheta;
	jetNormalizedResiduals[ 3 ] = jet2ThetaResidual / jet2SigmaTheta;
	jetNormalizedResiduals[ 4 ] = jet1PhiResidual / jet1SigmaPhi;
	jetNormalizedResiduals[ 5 ] = jet2PhiResidual / jet2SigmaPhi;

}

void ZHllqq5CFit::getNeutrinosInJet( LCRelationNavigator JetSLDNav , LCRelationNavigator SLDNuNav , EVENT::ReconstructedParticle* jet , pfoVectorVector &Neutrinos , std::vector<std::vector<int>> &NeutrinoInJetCombinations )
{
	const EVENT::LCObjectVec& SLDVertices = JetSLDNav.getRelatedToObjects( jet );
	Neutrinos.clear();
	NeutrinoInJetCombinations.clear();
	pfoVector NeutrinosOfThisSLD;
	std::vector<std::vector<int>>  neutrinoSolutions;
	std::vector<int> i_NeutrinoSolution;
	std::vector<int> singleCombination;
	for ( unsigned int i_sld = 0 ; i_sld < SLDVertices.size() ; ++i_sld )
	{
		Vertex *sldVertex = (Vertex*) SLDVertices.at( i_sld );
		NeutrinosOfThisSLD.clear();
		i_NeutrinoSolution.clear();
		singleCombination.push_back( 0 );
		const EVENT::LCObjectVec& neutrinos = SLDNuNav.getRelatedToObjects( sldVertex );
		for ( unsigned int i_nu = 0 ; i_nu < neutrinos.size() ; ++i_nu )
		{
			EVENT::ReconstructedParticle* neutrino = (ReconstructedParticle*) neutrinos.at( i_nu );
			NeutrinosOfThisSLD.push_back( neutrino );
			i_NeutrinoSolution.push_back( i_nu );
		}
		Neutrinos.push_back( NeutrinosOfThisSLD );
		neutrinoSolutions.push_back( i_NeutrinoSolution );
	}

	int nCombinations = 1;
	for ( unsigned int i_sld = 0 ; i_sld < neutrinoSolutions.size() ; ++i_sld ) nCombinations *= neutrinoSolutions[ i_sld ].size();
	for ( int i = 0 ; i < nCombinations ; ++i ) NeutrinoInJetCombinations.push_back( singleCombination );
	for ( unsigned int i_sld = 0 ; i_sld < neutrinoSolutions.size() ; ++i_sld )
	{
		int nUpSLD = 1;
		int nDownSLD = 1;
		for ( unsigned int j = 0 ; j < i_sld ; ++j )
		{
			nUpSLD *= neutrinoSolutions[ j ].size();
		}
		for ( unsigned int i = i_sld + 1 ; i < neutrinoSolutions.size() ; ++i )
		{
			nDownSLD *= neutrinoSolutions[ i ].size();
		}
		streamlog_out(MESSAGE) << "Number of Big loops in sld[ " << i_sld << " ] = 	" << nUpSLD << std::endl;
		streamlog_out(MESSAGE) << "Number of Small loops in sld[ " << i_sld << " ] = 	" << nDownSLD << std::endl;
		for ( int i = 0 ; i < nUpSLD ; ++i )
		{
			for ( unsigned int i_nu = 0 ; i_nu < neutrinoSolutions[ i_sld ].size() ; ++i_nu )
			{
				for ( int j = 0 ; j < nDownSLD ; ++j )
				{
					streamlog_out(MESSAGE) << "i = " << i << " , i_nu = " << i_nu << " , j = " << j ;
					streamlog_out(MESSAGE) << " ,	Filling sldcombination[ " << j + i * neutrinoSolutions[ i_sld ].size() * nDownSLD + i_nu * nDownSLD << " ][ " << i_sld << " ] = 	" << neutrinoSolutions[ i_sld ][ i_nu ] << std::endl;
					NeutrinoInJetCombinations[ j + i * neutrinoSolutions[ i_sld ].size() * nDownSLD + i_nu * nDownSLD ][ i_sld ] = neutrinoSolutions[ i_sld ][ i_nu ];
				}
			}
		}
	}
	streamlog_out(MESSAGE) << "Total Number of SLD combinations in jet: " << NeutrinoInJetCombinations.size() << " combinations for " << NeutrinoInJetCombinations[ 0 ].size() << " semi-leptonic decays!" << std::endl;
	streamlog_out(MESSAGE) << "Possible combinations:" << std::endl;
	for ( unsigned int i_comb = 0 ; i_comb < NeutrinoInJetCombinations.size() ; ++i_comb )
	{
		streamlog_out(MESSAGE) << "		combination " << i_comb << ":	";
		for ( unsigned int i_sld = 0 ; i_sld < NeutrinoInJetCombinations[ i_comb ].size() ; ++i_sld )
		{
			streamlog_out(MESSAGE) << " " << NeutrinoInJetCombinations[ i_comb ][ i_sld ] << " ";
		}
		streamlog_out(MESSAGE) << " " << std::endl;
	}
}

int ZHllqq5CFit::performFIT( 	TLorentzVector jet1FourMomentum , std::vector<float> jet1CovMat , TLorentzVector jet2FourMomentum , std::vector<float> jet2CovMat , pfoVector leptons ,
				float &fitProbability , float (&fitOutputs)[ 18 ] , std::vector< TLorentzVector > &fittedObjects , float (&pull)[ 12 ] , bool traceEvent )
{
	JetFitObject *jet[ 2 ];
	double sigmaE , sigmaTheta , sigmaPhi;
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
	jet[ 0 ] = new JetFitObject( jet1FourMomentum.E() , jet1FourMomentum.Theta() , jet1FourMomentum.Phi() , m_SigmaEnergyScaleFactor * sigmaE , sigmaTheta , sigmaPhi , jet1FourMomentum.M() );
	jet[ 0 ]->setName( "jet1" );
	streamlog_out(MESSAGE)  << " start four-vector of jet1: " << *jet[ 0 ]  << std::endl ;
	streamlog_out(DEBUG6) << "		get jet2 resolutions"  << std::endl ;
	getJetResolutions( jet2FourMomentum , jet2CovMat , sigmaE , sigmaTheta , sigmaPhi );
	jet[ 1 ] = new JetFitObject( jet2FourMomentum.E() , jet2FourMomentum.Theta() , jet2FourMomentum.Phi() , m_SigmaEnergyScaleFactor * sigmaE , sigmaTheta , sigmaPhi , jet2FourMomentum.M() );
	jet[ 1 ]->setName( "jet2" );
	streamlog_out(MESSAGE)  << " start four-vector of jet2: " << *jet[ 1 ]  << std::endl ;

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
	streamlog_out(MESSAGE)  << " start four-vector of lepton1: " << *lepton[ 0 ]  << std::endl ;
	streamlog_out(DEBUG6) << "		get lepton2 parameters"  << std::endl ;
	getLeptonParameters( leptons[ 1 ] , parameters , errors );
	lepton[ 1 ] = new LeptonFitObject ( parameters[ 0 ] , parameters[ 1 ] , parameters[ 2 ] , errors[ 0 ] , errors[ 1 ] , errors[ 2 ] , leptons[ 1 ]->getMass() );
	lepton[ 1 ]->setName( "lepton2" );
	streamlog_out(MESSAGE)  << " start four-vector of lepton2: " << *lepton[ 1 ]  << std::endl ;

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
		streamlog_out(DEBUG8)  << "	start four-vector of lepton " << i << ": " << *(LEPTONs[i])  << std::endl ;
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
		streamlog_out(MESSAGE)  << "start four-vector of ISR photon: " << *(photon) << std::endl ;

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

	streamlog_out(MESSAGE)  << "	<<<<<<<<<<<<<<<<<<<<<<<<<< FIT OUTPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

	for (int i = 0; i < NJETS; ++i)
	{
		streamlog_out(MESSAGE)  << "final four-vector of jet " << i << ": " << *(JETs[i]) << std::endl ;
		streamlog_out(DEBUG8)  << "final px of jet " << i << ": " << (JETs[i]) << std::endl ;
	}
	for (int i = 0; i < NLEPTONS; ++i)
	{
		streamlog_out(MESSAGE)  << "final four-vector of lepton " << i << ": " << *(LEPTONs[i]) << std::endl ;
		streamlog_out(DEBUG8)  << "final px of lepton " << i << ": " << (LEPTONs[i]) << std::endl ;
	}
	streamlog_out(MESSAGE) << "	fit probability = " << fitProbability << std::endl ;
	streamlog_out(MESSAGE) << "	fit chi2 = " << chi2  << std::endl ;
	streamlog_out(MESSAGE) << "	error code: " << ErrorCode << std::endl ;
	if( m_fitISR ) streamlog_out(MESSAGE)  << "final four-vector of ISR photon: " << *(photon) << std::endl;

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
//					sigma = sqrt(fabs(errmea*errmea-errfit*errfit));
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
//					sigma = sqrt(fabs(errmea*errmea-errfit*errfit));
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


	streamlog_out(MESSAGE)  << "start mass of Z: " << ZstartMass << std::endl ;
	streamlog_out(MESSAGE)  << "start mass of H: " << HstartMass << std::endl ;
	streamlog_out(MESSAGE)  << "mass of fitted Z: " << Zmass_after_fit << std::endl ;
	streamlog_out(MESSAGE)  << "mass of fitted H: " << Hmass_after_fit << std::endl ;
	streamlog_out(DEBUG9)  << "Error Code: " << ErrorCode << std::endl ;
	streamlog_out(DEBUG9)  << "Fit probability: " << fitProbability << std::endl ;

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


void ZHllqq5CFit::getJetResolutions(	TLorentzVector jetFourMomentum , std::vector<float> jetCovMat , double &sigmaE , double &sigmaTheta , double &sigmaPhi )
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

	sigmaE		= std::sqrt( jetCovMat[ 9 ] );
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

void ZHllqq5CFit::getJetResiduals( TVector3 jetTrueMomentum , double jetTrueEnergy , TVector3 jetRecoMomentum , double jetRecoEnergy , double &jetEnergyResidual , double &jetThetaResidual , double &jetPhiResidual )
{
	jetEnergyResidual = jetRecoEnergy - jetTrueEnergy;
	TVector3 jetTrueMomentumUnit = jetTrueMomentum; jetTrueMomentumUnit.SetMag( 1.0 );
	TVector3 jetTruePtUnit( jetTrueMomentum.Px() , jetTrueMomentum.Py() , 0.0 ); jetTruePtUnit.SetMag( 1.0 );
	TVector3 jetRecoMomentumUnitRotated = jetRecoMomentum; jetRecoMomentumUnitRotated.SetPhi( jetTrueMomentum.Phi() ); jetRecoMomentumUnitRotated.SetMag( 1.0 );
	TVector3 jetRecoPtUnit( jetRecoMomentum.Px() , jetRecoMomentum.Py() , 0.0 ); jetRecoPtUnit.SetMag( 1.0 );
	jetThetaResidual = ( jetRecoMomentum.Theta() >= jetTrueMomentum.Theta() ? acos( jetTrueMomentumUnit.Dot( jetRecoMomentumUnitRotated ) ) : -1.0 * acos( jetTrueMomentumUnit.Dot( jetRecoMomentumUnitRotated ) ) );
	jetPhiResidual = ( jetRecoMomentum.Phi() >= jetTrueMomentum.Phi() ? acos( jetTruePtUnit.Dot( jetRecoPtUnit ) ) : -1.0 * acos( jetTruePtUnit.Dot( jetRecoPtUnit ) ) );
}

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
