#include "ZHllqq5CFit.h"
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
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>

#include <GeometryUtil.h>
#include "TLorentzVector.h"
#include <CLHEP/Vector/LorentzVector.h>
#include "JetFitObject.h"
#include "LeptonFitObject.h"
#include "ISRPhotonFitObject.h"
#include "MomentumConstraint.h"
#include "OPALFitterGSL.h"
#include "NewFitterGSL.h"
#include "TextTracer.h"
#include "NewtonFitterGSL.h"
#include "FourJetZHPairing.h"
#include "MassConstraint.h"
#include "SoftGaussParticleConstraint.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/LCCollectionVec.h"
#include <EVENT/LCCollection.h>
#include "TFile.h"
#include "TH1F.h"
#include "TH2I.h"
#include "TH2F.h"
#include "TTree.h"
#define PI 3.14159265

using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace CLHEP ;


ZHllqq5CFit aZHllqq5CFit ;

// function to define the jet energy resolution (in GeV)
double ZHllqq5CFit::JetEnergyResolution(double E)
{
	// examples here derived by Benjamin Hermberg from e+e- -> udsc:
	// 1) default is 120%/sqrt(E), gives best convergence of 5C fit on e+e- -> udsc
	double result = m_jetEnergyError*std::sqrt(E);

	// 2) comparing jet-level to quark-level energies
	//    (using MarlinReco/Analysis/RecoMCTruthLink/QuarkJetPairing.cc)
	if (m_jetEnergyError == 0 ) result = std::sqrt(pow(0.6908,2)*(E)+(pow(0.02596,2)*pow(E,2)));

	return result;
}

ZHllqq5CFit::ZHllqq5CFit() :
Processor("ZHllqq5CFit"),
m_includeHbb(false),
m_includeHcc(false),
m_includeHgg(false),
m_includeHother(false),
m_nRun(0),
m_nEvt(0),
m_nRunSum(0),
m_nEvtSum(0),
m_Bfield(0.f),
c(0.),
mm2m(0.),
eV2GeV(0.),
eB(0.),
m_ECM(250.f),
m_nSLDecayBHadron(0),
m_nSLDecayCHadron(0),
m_nSLDecayTotal(0),
m_nJets(0),
m_nLeptons(0),
m_HDecayMode(-1),
m_dijet_angle_wNu_bestfit(0.),
m_dijet_angle_woNu(0.),
m_dijet_angle_best(0.),
m_iError_wNu_bestfit(-5),
m_probability_wNu_bestfit(0.),
m_chi2_wNu_bestfit(0.),
m_n_itter_wNu_bestfit(0),
m_startmassZ_wNu_bestfit(0.),
m_startmassH_wNu_bestfit(0.),
m_beststartmassZ_wNu_bestfit(0.),
m_beststartmassH_wNu_bestfit(0.),
m_Zmass_after_fit_wNu_bestfit(0.),
m_Hmass_after_fit_wNu_bestfit(0.),
m_bestphotonenergy_wNu_bestfit(0.),
m_chi2startmassZ_wNu_bestfit(0.),
m_chi2startmassH_wNu_bestfit(0.),
m_iError_woNu(-5),
m_probability_woNu(0.),
m_chi2best_woNu(0.),
m_n_itter_woNu(0),
m_startmassZ_woNu(0.),
m_startmassH_woNu(0.),
m_beststartmassZ_woNu(0.),
m_beststartmassH_woNu(0.),
m_Zmass_after_fit_woNu(0.),
m_Hmass_after_fit_woNu(0.),
m_bestphotonenergy_woNu(0.),
m_chi2startmassZ_woNu(0.),
m_chi2startmassH_woNu(0.),
m_iError_best(-5),
m_probability_best(0.),
m_chi2_best(0.),
m_n_itter_best(0),
m_startmassZ_best(0.),
m_startmassH_best(0.),
m_beststartmassZ_best(0.),
m_beststartmassH_best(0.),
m_Zmass_after_fit_best(0.),
m_Hmass_after_fit_best(0.),
m_bestphotonenergy_best(0.),
m_chi2startmassZ_best(0.),
m_chi2startmassH_best(0.),
chi2best(0.),
errorcode(0),
E_lab(0.),
m_pxc_before_ISR_wNu(0.),
m_pyc_before_ISR_wNu(0.),
m_pzc_before_ISR_wNu(0.),
m_ec_before_ISR_wNu(0.),
m_pxc_before_fit_wNu(0.),
m_pyc_before_fit_wNu(0.),
m_pzc_before_fit_wNu(0.),
m_ec_before_fit_wNu(0.),
m_pxc_after_fit_wNu(0.),
m_pyc_after_fit_wNu(0.),
m_pzc_after_fit_wNu(0.),
m_ec_after_fit_wNu(0.),
m_zc_before_ISR_wNu(0.),
m_zc_before_fit_wNu(0.),
m_zc_after_fit_wNu(0.),
m_pxc_before_ISR_woNu(0.),
m_pyc_before_ISR_woNu(0.),
m_pzc_before_ISR_woNu(0.),
m_ec_before_ISR_woNu(0.),
m_pxc_before_fit_woNu(0.),
m_pyc_before_fit_woNu(0.),
m_pzc_before_fit_woNu(0.),
m_ec_before_fit_woNu(0.),
m_pxc_after_fit_woNu(0.),
m_pyc_after_fit_woNu(0.),
m_pzc_after_fit_woNu(0.),
m_ec_after_fit_woNu(0.),
m_zc_before_ISR_woNu(0.),
m_zc_before_fit_woNu(0.),
m_zc_after_fit_woNu(0.),
m_pxc_before_ISR_best(0.),
m_pyc_before_ISR_best(0.),
m_pzc_before_ISR_best(0.),
m_ec_before_ISR_best(0.),
m_pxc_before_fit_best(0.),
m_pyc_before_fit_best(0.),
m_pzc_before_fit_best(0.),
m_ec_before_fit_best(0.),
m_pxc_after_fit_best(0.),
m_pyc_after_fit_best(0.),
m_pzc_after_fit_best(0.),
m_ec_after_fit_best(0.),
m_zc_before_ISR_best(0.),
m_zc_before_fit_best(0.),
m_zc_after_fit_best(0.),
m_mcHmumuPx(0.),
m_mcHmumuPy(0.),
m_mcHmumuPz(0.),
m_mcHmumuPt(0.),
m_mcHmumuP(0.),
m_mcHmumuE(0.),
m_mcISR1Px(0.),
m_mcISR1Py(0.),
m_mcISR1Pz(0.),
m_mcISR1E(0.),
m_mcISR2Px(0.),
m_mcISR2Py(0.),
m_mcISR2Pz(0.),
m_mcISR2E(0.),
m_mcISRPx(0.),
m_mcISRPy(0.),
m_mcISRPz(0.),
m_mcISRE(0.),
m_ISR_startPx_wNu(0.),
m_ISR_startPy_wNu(0.),
m_ISR_startPz_wNu(0.),
m_ISR_startPx_woNu(0.),
m_ISR_startPy_woNu(0.),
m_ISR_startPz_woNu(0.),
m_TotalSigmaPx2(0.),
m_TotalSigmaPy2(0.),
m_TotalSigmaPz2(0.),
m_TotalSigmaE2(0.),
m_TotalSigmaPxPy(0.),
m_TotalSigmaPxPz(0.),
m_TotalSigmaPxE(0.),
m_TotalSigmaPyPz(0.),
m_TotalSigmaPyE(0.),
m_TotalSigmaPzE(0.),
m_ISRPx_total(0.),
m_ISRPy_total(0.),
m_ISRPz_total(0.),
m_ISRPt_total(0.),
m_ISRE_total(0.),
m_FSRPx_total(0.),
m_FSRPy_total(0.),
m_FSRPz_total(0.),
m_FSRPt_total(0.),
m_FSRE_total(0.),
m_pTFile(NULL),
m_pTTree(NULL),
m_pTTree_0(NULL),
m_pTTree_1(NULL),
m_pTTree_2(NULL),
m_pTTree_3(NULL),
m_pTTree_4(NULL),
m_pTTree_5(NULL),
m_pTTree_6(NULL),
h_HDecayMode(NULL),
h_Zmass_beforefit_woNu(NULL),
h_Hmass_beforefit_woNu(NULL),
h_Zmass_beforefit_wNu(NULL),
h_Hmass_beforefit_wNu(NULL),
h_Zmass_afterfit_woNu(NULL),
h_Hmass_afterfit_woNu(NULL),
h_Zmass_afterfit_wNu(NULL),
h_Hmass_afterfit_wNu(NULL),
h_fitError_wNu(NULL),
h_fitError_woNu(NULL),
h_ErrorCode_wNu_woNu(NULL),
h_fitProbability_wNu(NULL),
h_fitProbability_woNu(NULL),
h_fitProbability_best(NULL),
h_nJets(NULL),
h_nLeptons(NULL),
h_nLeptons_nJets(NULL),
h_ISRE_1mcp_fit(NULL),
h_ISRE_2mcp_fit(NULL),
h_ISRpz_1mcp_fit(NULL),
h_ISRpz_2mcp_fit(NULL),
h_ISR_pzc_wNu(NULL),
h_ISR_pzc_woNu(NULL),
h_ISR_pzc_best(NULL),
h_mcpISRPx(NULL),
h_mcpISRPy(NULL),
h_mcpISRPz(NULL),
h_mcpISRE(NULL),
h_pull_jet_E(NULL),
h_pull_jet_theta(NULL),
h_pull_jet_phi(NULL),
h_pull_lepton_InvPt(NULL),
h_pull_lepton_theta(NULL),
h_pull_lepton_phi(NULL),
h_error_jet_E(NULL),
h_error_jet_Theta(NULL),
h_error_jet_Phi(NULL),
h_error_lepton_InvpT(NULL),
h_error_lepton_Theta(NULL),
h_error_lepton_Phi(NULL),
h_SigmaPx2(NULL),
h_SigmaPy2(NULL),
h_SigmaPz2(NULL),
h_SigmaE2(NULL),
h_SigmaPxPy(NULL),
h_SigmaPxPz(NULL),
h_SigmaPxE(NULL),
h_SigmaPyPz(NULL),
h_SigmaPyE(NULL),
h_SigmaPzE(NULL),
h_fitProbability_diJetAngle(NULL),
h_fitProbability_Ejet(NULL),
h_fitProbability_Thetajet(NULL),
h_fitProbability_Phijet(NULL),
h_fitProbability_SigmaEjet(NULL),
h_fitProbability_SigmaThetajet(NULL),
h_fitProbability_SigmaPhijet(NULL),
h_fitProbability_pullEjet(NULL),
h_fitProbability_pullThetajet(NULL),
h_fitProbability_pullPhijet(NULL),
h_fitProbability_constraintPx(NULL),
h_fitProbability_constraintPy(NULL),
h_fitProbability_constraintPz(NULL),
h_fitProbability_constraintE(NULL),
h_fitProbability_sigmaPx(NULL),
h_fitProbability_sigmaPy(NULL),
h_fitProbability_sigmaPz(NULL),
h_fitProbability_sigmaE(NULL),
h_constraintPx_lowFitProb(NULL),
h_constraintPx_midFitProb(NULL),
h_constraintPx_highFitProb(NULL),
h_constraintPy_lowFitProb(NULL),
h_constraintPy_midFitProb(NULL),
h_constraintPy_highFitProb(NULL),
h_constraintPz_lowFitProb(NULL),
h_constraintPz_midFitProb(NULL),
h_constraintPz_highFitProb(NULL),
h_constraintE_lowFitProb(NULL),
h_constraintE_midFitProb(NULL),
h_constraintE_highFitProb(NULL),
h_constraintPx_uncertaintyPx(NULL),
h_constraintPy_uncertaintyPy(NULL),
h_constraintPz_uncertaintyPz(NULL),
h_constraintE_uncertaintyE(NULL),
h_constraintPx_uncertaintyPx_lowFitProb(NULL),
h_constraintPx_uncertaintyPx_midFitProb(NULL),
h_constraintPx_uncertaintyPx_highFitProb(NULL),
h_constraintPy_uncertaintyPy_lowFitProb(NULL),
h_constraintPy_uncertaintyPy_midFitProb(NULL),
h_constraintPy_uncertaintyPy_highFitProb(NULL),
h_constraintPz_uncertaintyPz_lowFitProb(NULL),
h_constraintPz_uncertaintyPz_midFitProb(NULL),
h_constraintPz_uncertaintyPz_highFitProb(NULL),
h_constraintE_uncertaintyE_lowFitProb(NULL),
h_constraintE_uncertaintyE_midFitProb(NULL),
h_constraintE_uncertaintyE_highFitProb(NULL),
h_pull_jet_theta_jet_phi_lowFitProb(NULL),
h_pull_jet_theta_jet_phi_midFitProb(NULL),
h_pull_jet_theta_jet_phi_highFitProb(NULL),
h_NormalizedResidualTheta_mcQmcJ_ThetaMcJet(NULL),
h_NormalizedResidualPhi_mcQmcJ_PhiMcJet(NULL),
h_NormalizedResidualEnergy_mcQmcJ_EnergyMcJet(NULL),
h_NormalizedResidualTheta_mcQmcJ_EnergyMcJet(NULL),
h_NormalizedResidualPhi_mcQmcJ_EnergyMcJet(NULL),
h_NormalizedResidualTheta_mcJrecJ_ThetaRecJet(NULL),
h_NormalizedResidualPhi_mcJrecJ_PhiRecJet(NULL),
h_NormalizedResidualEnergy_mcJrecJ_EnergyRecJet(NULL),
h_NormalizedResidualTheta_mcJrecJ_EnergyRecJet(NULL),
h_NormalizedResidualPhi_mcJrecJ_EnergyRecJet(NULL),
h_NormalizedResidualTheta_mcQrecJ_ThetaRecJet(NULL),
h_NormalizedResidualPhi_mcQrecJ_PhiRecJet(NULL),
h_NormalizedResidualEnergy_mcQrecJ_EnergyRecJet(NULL),
h_NormalizedResidualTheta_mcQrecJ_EnergyRecJet(NULL),
h_NormalizedResidualPhi_mcQrecJ_EnergyRecJet(NULL),
h_recE_fitE(NULL),
h_recPhifitPhi_recThetafitTheta(NULL),
h_recPhifitPhi_recThetafitTheta_lowfit(NULL),
h_recPhifitPhi_recThetafitTheta_midfit(NULL),
h_recPhifitPhi_recThetafitTheta_highfit(NULL),
h_recPhifitPhi_recThetafitTheta_SigmaThetaPhi(NULL),
h_recPhifitPhi_recThetafitTheta_lowfit_SigmaThetaPhi(NULL),
h_recPhifitPhi_recThetafitTheta_midfit_SigmaThetaPhi(NULL),
h_recPhifitPhi_recThetafitTheta_highfit_SigmaThetaPhi(NULL)
{

//	modify processor description
	_description = "ZHllqq5CFit does a fit on 2 jet events (Px, Py, Pz, E, M12 = MZ)" ;


//	register steering parameters: name, description, class-variable, default value

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"JetCollectionName" ,
				"Name of the Jet collection"  ,
				jetcollection ,
				std::string("Durham_2Jets")
				);

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"ErrorFlowCollection" ,
				"Name of the ErrorFlow collection"  ,
				errorflowcollection ,
				std::string("OutputErrorFlowJets")
				);

	registerInputCollection( LCIO::MCPARTICLE,
				"HiggsDecayMode" ,
				"Higgs decay mode"  ,
				hDecayMode ,
				std::string("HdecayMode")
				);

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"ISOLeptonCollectionName" ,
				"Name of the Isolated Lepton collection"  ,
				leptoncollection ,
				std::string("ISOLeptons")
				);

	registerInputCollection( LCIO::MCPARTICLE,
				"SemiLeptonicDecays",
				"Semi-Leptonic Decays Collection",
				SLDecayCollection,
				std::string("")
				);

	registerInputCollection( LCIO::MCPARTICLE,
				"NeutrinoCorrectionB",
				"Collection of Corrected/Estimated neutrino energies from SemiLeptonic decays of B-Hadron",
				NuEnergyCollectionB,
				std::string("NuCorrectB")
				);

	registerInputCollection( LCIO::MCPARTICLE,
				"NeutrinoCorrectionC",
				"Collection of Corrected/Estimated neutrino energies from SemiLeptonic decays of B-Hadron",
				NuEnergyCollectionC,
				std::string("NuCorrectC")
				);

	registerInputCollection( LCIO::MCPARTICLE,
				"MCParticleCollection" ,
				"Name of the MCParticle collection"  ,
				MCPCollection ,
				std::string("MCParticle")
				);

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"MCJetCollectionwNuName" ,
				"Name of the MCJet collection with neutrinos"  ,
				mcjetcollection_wNu ,
				std::string("MCTruthJets_wNu")
				);

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"MCJetCollectionwoNuName" ,
				"Name of the MCJet collection without neutrinos"  ,
				mcjetcollection_woNu ,
				std::string("MCTruthJets_woNu")
				);

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"MCISRCollectionwNuName" ,
				"Name of the MCISR collection"  ,
				mcISRcol ,
				std::string("mcISRcol")
				);

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"MCFSRCollectionwNuName" ,
				"Name of the MCFSR collection"  ,
				mcFSRcol ,
				std::string("mcFSRcol")
				);

	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"FitOutputColection",
				" Output Fit Colection" ,
				outputFitcollection,
				std::string("FitReco")
				);

	registerProcessorParameter("includeHbb",
				"Include H->bb decays",
				m_includeHbb,
				bool(false)
				);

	registerProcessorParameter("includeHcc",
				"Include H->cc decays",
				m_includeHcc,
				bool(false)
				);

	registerProcessorParameter("includeHgg",
				"Include H->gg decays",
				m_includeHgg,
				bool(false)
				);

	registerProcessorParameter("includeHother",
				"Include other Higgs decay modes",
				m_includeHother,
				bool(false)
				);

	registerProcessorParameter("includeISR",
				"Include ISR in fit hypothesis; false: without ISR , true: with ISR",
				m_fitISR,
				bool(true)
				);

	registerProcessorParameter("includeNuCorrection",
				"Include neutrino correction in fit hypothesis; false: without Nu Correction , true: with Nu Correction",
				m_fitNuE,
				bool(true)
				);

	registerProcessorParameter( "useErrorFlow",
				"If true, use covariance matrix for energy uncertainty. Otherwise 1.2/sqrt(E)",
				m_useErrorFlow,
				bool(true)
				);

	registerProcessorParameter( "ECM" ,
				"Center-of-Mass Energy in GeV",
				m_ECM,
				float(250.f)
				);

	registerProcessorParameter( "ISRPzMax" ,
				"Maximum possible energy for a single ISR photon",
				m_isrpzmax,
				(float)35.
				);

	registerProcessorParameter( "SigmaEnergyScaleFactor",
				"Scale Factor t be applied to jet energy uncertainty",
				sigmaScaleFactor,
				(float) 1.0
				);

	registerProcessorParameter( "errene" ,
				"assumed energy resolution for jets as x/sqrt(E) - if 0, then parametrisation is used",
				m_jetEnergyError,
				(double)1.2
				);

	registerProcessorParameter( "errtheta" ,
				"assumed theta resolution for jet axis",
				m_jetThetaError,
				(double)0.1
				);

	registerProcessorParameter( "errphi" ,
				"assumed phi resolution for jet axis",
				m_jetPhiError,
				(double)0.1
				);

	registerProcessorParameter( "fitter" ,
				"0 = OPALFitter, 1 = NewFitter, 2 = NewtonFitter",
				m_fitter,
				(int)0
				);

	registerProcessorParameter( "traceall" ,
				"set true if every event should be traced",
				m_traceall,
				(bool)false
				);

	registerProcessorParameter( "ievttrace" ,
				"number of individual event to be traced",
				m_ievttrace,
				(int)0
				);

	registerProcessorParameter( "JetResThetaCoeff" ,
				"Coefficient to Theta resolution of the jet",
				m_JetResThetaCoeff,
				(float)1.0
				);

	registerProcessorParameter( "JetResPhiCoeff" ,
				"Coefficient to Phi resolution of the jet",
				m_JetResPhiCoeff,
				(float)1.0
				);

	registerProcessorParameter("outputFilename",
				"name of output file",
				m_outputFile,
				std::string("")
				);

}

void ZHllqq5CFit::init()
{
//	usually a good idea to
	streamlog_out(DEBUG) << "   init called  " << std::endl;
	this->Clear();
	m_Bfield = MarlinUtil::getBzAtOrigin();
	printParameters();
	streamlog_out(DEBUG) << " BField =  "<< m_Bfield << " Tesla" << std::endl ;
	c = 2.99792458e8;
	mm2m = 1e-3;
	eV2GeV = 1e-9;
	eB = m_Bfield * c * mm2m * eV2GeV;

	m_nRun = 0;
	m_nEvt = 0;
	m_nRunSum = 0;
	m_nEvtSum = 0;

	b = (double) 0.00464564 * ( std::log( m_ECM * m_ECM * 3814714. ) - 1. );
//	  = 2*alpha/pi*( ln(s/m_e^2)-1 )
	ISRPzMaxB = std::pow((double)m_isrpzmax,b);

	m_pTFile = new TFile(m_outputFile.c_str(),"recreate");

	m_pTTree = new TTree("eventTree","eventTree");
	m_pTTree->SetDirectory(m_pTFile);
	m_pTTree->Branch("run", &m_nRun, "run/I");
	m_pTTree->Branch("event", &m_nEvt, "event/I");
	m_pTTree->Branch("nJets",&m_nJets,"nJets/I") ;
	m_pTTree->Branch("nLeptons",&m_nLeptons,"nLeptons/I") ;
	m_pTTree->Branch("HDecayMode",&m_HDecayMode,"HDecayMode/I") ;

	m_pTTree_0 = new TTree("ZHlljjTree","ZHlljjTree");
	m_pTTree_0->SetDirectory(m_pTFile);
	m_pTTree_0->Branch("nSLDecayBHadron",&m_nSLDecayBHadron,"nSLDecayBHadron/I") ;
	m_pTTree_0->Branch("nSLDecayCHadron",&m_nSLDecayCHadron,"nSLDecayCHadron/I") ;
	m_pTTree_0->Branch("nSLDecayTotal",&m_nSLDecayTotal,"nSLDecayTotal/I") ;
	m_pTTree_0->Branch("bestNuCombination",&m_bestNuCombination) ;
	m_pTTree_0->Branch("iError_wNu",&m_iError_wNu) ;
	m_pTTree_0->Branch("probability_wNu",&m_probability_wNu) ;
	m_pTTree_0->Branch("chi2_wNu",&m_chi2_wNu) ;
	m_pTTree_0->Branch("n_itter_wNu",&m_n_itter_wNu) ;
	m_pTTree_0->Branch("startmassZ_wNu",&m_startmassZ_wNu) ;
	m_pTTree_0->Branch("startmassH_wNu",&m_startmassH_wNu) ;
	m_pTTree_0->Branch("beststartmassZ_wNu",&m_beststartmassZ_wNu) ;
	m_pTTree_0->Branch("beststartmassH_wNu",&m_beststartmassH_wNu) ;
	m_pTTree_0->Branch("Zmass_after_fit_wNu",&m_Zmass_after_fit_wNu) ;
	m_pTTree_0->Branch("Hmass_after_fit_wNu",&m_Hmass_after_fit_wNu) ;
	m_pTTree_0->Branch("chi2startmassZ_wNu",&m_chi2startmassZ_wNu) ;
	m_pTTree_0->Branch("chi2startmassH_wNu",&m_chi2startmassH_wNu) ;
	m_pTTree_0->Branch("iError_wNu_bestfit",&m_iError_wNu_bestfit,"iError_wNu_bestfit/I") ;
	m_pTTree_0->Branch("probability_wNu_bestfit",&m_probability_wNu_bestfit,"probability_wNu_bestfit/F") ;
	m_pTTree_0->Branch("chi2_wNu_bestfit",&m_chi2_wNu_bestfit,"chi2_wNu_bestfit/F") ;
	m_pTTree_0->Branch("n_itter_wNu_bestfit",&m_n_itter_wNu_bestfit,"n_itter_wNu_bestfit/I") ;
	m_pTTree_0->Branch("startmassZ_wNu_bestfit",&m_startmassZ_wNu_bestfit,"startmassZ_wNu_bestfit/F") ;
	m_pTTree_0->Branch("startmassH_wNu_bestfit",&m_startmassH_wNu_bestfit,"startmassH_wNu_bestfit/F") ;
	m_pTTree_0->Branch("beststartmassZ_wNu_bestfit",&m_beststartmassZ_wNu_bestfit,"beststartmassZ_wNu_bestfit/F") ;
	m_pTTree_0->Branch("beststartmassH_wNu_bestfit",&m_beststartmassH_wNu_bestfit,"beststartmassH_wNu_bestfit/F") ;
	m_pTTree_0->Branch("Zmass_after_fit_wNu_bestfit",&m_Zmass_after_fit_wNu_bestfit,"Zmass_after_fit_wNu_bestfit/F") ;
	m_pTTree_0->Branch("Hmass_after_fit_wNu_bestfit",&m_Hmass_after_fit_wNu_bestfit,"Hmass_after_fit_wNu_bestfit/F") ;
	m_pTTree_0->Branch("chi2startmassZ_wNu_bestfit",&m_chi2startmassZ_wNu_bestfit,"chi2startmassZ_wNu_bestfit/F") ;
	m_pTTree_0->Branch("chi2startmassH_wNu_bestfit",&m_chi2startmassH_wNu_bestfit,"chi2startmassH_wNu_bestfit/F") ;
	m_pTTree_0->Branch("bestphotonenergy_wNu_bestfit",&m_bestphotonenergy_wNu_bestfit,"bestphotonenergy_wNu_bestfit/F") ;
	m_pTTree_0->Branch("iError_woNu",&m_iError_woNu,"iError_woNu/I") ;
	m_pTTree_0->Branch("probability_woNu",&m_probability_woNu,"probability_woNu/F") ;
	m_pTTree_0->Branch("chi2best_woNu",&m_chi2best_woNu,"chi2best_woNu/F") ;
	m_pTTree_0->Branch("n_itter_woNu",&m_n_itter_woNu,"n_itter_woNu/I") ;
	m_pTTree_0->Branch("startmassZ_woNu",&m_startmassZ_woNu,"startmassZ_woNu/F") ;
	m_pTTree_0->Branch("startmassH_woNu",&m_startmassH_woNu,"startmassH_woNu/F") ;
	m_pTTree_0->Branch("beststartmassZ_woNu",&m_beststartmassZ_woNu,"beststartmassZ_woNu/F") ;
	m_pTTree_0->Branch("beststartmassH_woNu",&m_beststartmassH_woNu,"beststartmassH_woNu/F") ;
	m_pTTree_0->Branch("Zmass_after_fit_woNu",&m_Zmass_after_fit_woNu,"Zmass_after_fit_woNu/F") ;
	m_pTTree_0->Branch("Hmass_after_fit_woNu",&m_Hmass_after_fit_woNu,"Hmass_after_fit_woNu/F") ;
	m_pTTree_0->Branch("chi2startmassZ_woNu",&m_chi2startmassZ_woNu,"chi2startmassZ_woNu/F") ;
	m_pTTree_0->Branch("chi2startmassH_woNu",&m_chi2startmassH_woNu,"chi2startmassH_woNu/F") ;
	m_pTTree_0->Branch("bestphotonenergy_woNu",&m_bestphotonenergy_woNu,"bestphotonenergy_woNu/F") ;
	m_pTTree_0->Branch("iError_best",&m_iError_best,"iError_best/I") ;
	m_pTTree_0->Branch("probability_best",&m_probability_best,"probability_best/F") ;
	m_pTTree_0->Branch("chi2_best",&m_chi2_best,"chi2_best/F") ;
	m_pTTree_0->Branch("n_itter_best",&m_n_itter_best,"n_itter_best/I") ;
	m_pTTree_0->Branch("startmassZ_best",&m_startmassZ_best,"startmassZ_best/F") ;
	m_pTTree_0->Branch("startmassH_best",&m_startmassH_best,"startmassH_best/F") ;
	m_pTTree_0->Branch("beststartmassZ_best",&m_beststartmassZ_best,"beststartmassZ_best/F") ;
	m_pTTree_0->Branch("beststartmassH_best",&m_beststartmassH_best,"beststartmassH_best/F") ;
	m_pTTree_0->Branch("Zmass_after_fit_best",&m_Zmass_after_fit_best,"Zmass_after_fit_best/F") ;
	m_pTTree_0->Branch("Hmass_after_fit_best",&m_Hmass_after_fit_best,"Hmass_after_fit_best/F") ;
	m_pTTree_0->Branch("chi2startmassZ_best",&m_chi2startmassZ_best,"chi2startmassZ_best/F") ;
	m_pTTree_0->Branch("chi2startmassH_best",&m_chi2startmassH_best,"chi2startmassH_best/F") ;
	m_pTTree_0->Branch("bestphotonenergy_best",&m_bestphotonenergy_best,"bestphotonenergy_best/F") ;
	m_pTTree_1 = new TTree("pulls","pulls");
	m_pTTree_1->SetDirectory(m_pTFile);
	m_pTTree_1->Branch("probability_best",&m_probability_best,"probability_best/F") ;
	m_pTTree_1->Branch("nSLDecayTotal",&m_nSLDecayTotal,"nSLDecayTotal/I") ;
	m_pTTree_1->Branch("ISRE_total",&m_ISRE_total,"ISRE_total/F") ;
	m_pTTree_1->Branch("FSRE_total",&m_FSRE_total,"FSRE_total/F") ;
	m_pTTree_1->Branch("pull_jet_E_wNu",&m_pull_jet_E_wNu) ;
	m_pTTree_1->Branch("pull_jet_th_wNu",&m_pull_jet_th_wNu) ;
	m_pTTree_1->Branch("pull_jet_phi_wNu",&m_pull_jet_phi_wNu) ;
	m_pTTree_1->Branch("pull_lepton_InvpT_wNu",&m_pull_lepton_InvpT_wNu) ;
	m_pTTree_1->Branch("pull_lepton_th_wNu",&m_pull_lepton_th_wNu) ;
	m_pTTree_1->Branch("pull_lepton_phi_wNu",&m_pull_lepton_phi_wNu) ;
	m_pTTree_1->Branch("pull_jet_E_wNu_bestfit",&m_pull_jet_E_wNu_bestfit);
	m_pTTree_1->Branch("pull_jet_th_wNu_bestfit",&m_pull_jet_th_wNu_bestfit);
	m_pTTree_1->Branch("pull_jet_phi_wNu_bestfit",&m_pull_jet_phi_wNu_bestfit);
	m_pTTree_1->Branch("pull_lepton_InvpT_wNu_bestfit",&m_pull_lepton_InvpT_wNu_bestfit);
	m_pTTree_1->Branch("pull_lepton_th_wNu_bestfit",&m_pull_lepton_th_wNu_bestfit);
	m_pTTree_1->Branch("pull_lepton_phi_wNu_bestfit",&m_pull_lepton_phi_wNu_bestfit);
	m_pTTree_1->Branch("pull_jet_E_woNu",&m_pull_jet_E_woNu) ;
	m_pTTree_1->Branch("pull_jet_th_woNu",&m_pull_jet_th_woNu) ;
	m_pTTree_1->Branch("pull_jet_phi_woNu",&m_pull_jet_phi_woNu) ;
	m_pTTree_1->Branch("pull_lepton_InvpT_woNu",&m_pull_lepton_InvpT_woNu) ;
	m_pTTree_1->Branch("pull_lepton_th_woNu",&m_pull_lepton_th_woNu) ;
	m_pTTree_1->Branch("pull_lepton_phi_woNu",&m_pull_lepton_phi_woNu) ;
	m_pTTree_1->Branch("pull_jet_E_best",&m_pull_jet_E_best);
	m_pTTree_1->Branch("pull_jet_th_best",&m_pull_jet_th_best);
	m_pTTree_1->Branch("pull_jet_phi_best",&m_pull_jet_phi_best);
	m_pTTree_1->Branch("pull_lepton_InvpT_best",&m_pull_lepton_InvpT_best);
	m_pTTree_1->Branch("pull_lepton_th_best",&m_pull_lepton_th_best);
	m_pTTree_1->Branch("pull_lepton_phi_best",&m_pull_lepton_phi_best);
	m_pTTree_2 = new TTree("constraints","constraints");
	m_pTTree_2->SetDirectory(m_pTFile);
	m_pTTree_2->Branch("nSLDecayTotal",&m_nSLDecayTotal,"nSLDecayTotal/I") ;
	m_pTTree_2->Branch("probability_best",&m_probability_best,"probability_best/F") ;
	m_pTTree_2->Branch("ISRE_total",&m_ISRE_total,"ISRE_total/F") ;
	m_pTTree_2->Branch("FSRE_total",&m_FSRE_total,"FSRE_total/F") ;
	m_pTTree_2->Branch("pxc_before_ISR_wNu",&m_pxc_before_ISR_wNu,"pxc_before_ISR_wNu/F") ;
	m_pTTree_2->Branch("pyc_before_ISR_wNu",&m_pyc_before_ISR_wNu,"pyc_before_ISR_wNu/F") ;
	m_pTTree_2->Branch("pzc_before_ISR_wNu",&m_pzc_before_ISR_wNu,"pzc_before_ISR_wNu/F") ;
	m_pTTree_2->Branch("ec_before_ISR_wNu",&m_ec_before_ISR_wNu,"ec_before_ISR_wNu/F") ;
	m_pTTree_2->Branch("zc_before_ISR_wNu",&m_zc_before_ISR_wNu,"zc_before_ISR_wNu/F") ;
	m_pTTree_2->Branch("pxc_before_fit_wNu",&m_pxc_before_fit_wNu,"pxc_before_fit_wNu/F") ;
	m_pTTree_2->Branch("pyc_before_fit_wNu",&m_pyc_before_fit_wNu,"pyc_before_fit_wNu/F") ;
	m_pTTree_2->Branch("pzc_before_fit_wNu",&m_pzc_before_fit_wNu,"pzc_before_fit_wNu/F") ;
	m_pTTree_2->Branch("ec_before_fit_wNu",&m_ec_before_fit_wNu,"ec_before_fit_wNu/F") ;
	m_pTTree_2->Branch("zc_before_fit_wNu",&m_zc_before_fit_wNu,"zc_before_fit_wNu/F") ;
	m_pTTree_2->Branch("pxc_after_fit_wNu",&m_pxc_after_fit_wNu,"pxc_after_fit_wNu/F") ;
	m_pTTree_2->Branch("pyc_after_fit_wNu",&m_pyc_after_fit_wNu,"pyc_after_fit_wNu/F") ;
	m_pTTree_2->Branch("pzc_after_fit_wNu",&m_pzc_after_fit_wNu,"pzc_after_fit_wNu/F") ;
	m_pTTree_2->Branch("ec_after_fit_wNu",&m_ec_after_fit_wNu,"ec_after_fit_wNu/F") ;
	m_pTTree_2->Branch("zc_after_fit_wNu",&m_zc_after_fit_wNu,"zc_after_fit_wNu/F") ;
	m_pTTree_2->Branch("pxc_before_ISR_woNu",&m_pxc_before_ISR_woNu,"pxc_before_ISR_woNu/F") ;
	m_pTTree_2->Branch("pyc_before_ISR_woNu",&m_pyc_before_ISR_woNu,"pyc_before_ISR_woNu/F") ;
	m_pTTree_2->Branch("pzc_before_ISR_woNu",&m_pzc_before_ISR_woNu,"pzc_before_ISR_woNu/F") ;
	m_pTTree_2->Branch("ec_before_ISR_woNu",&m_ec_before_ISR_woNu,"ec_before_ISR_woNu/F") ;
	m_pTTree_2->Branch("zc_before_ISR_woNu",&m_zc_before_ISR_woNu,"zc_before_ISR_woNu/F") ;
	m_pTTree_2->Branch("pxc_before_fit_woNu",&m_pxc_before_fit_woNu,"pxc_before_fit_woNu/F") ;
	m_pTTree_2->Branch("pyc_before_fit_woNu",&m_pyc_before_fit_woNu,"pyc_before_fit_woNu/F") ;
	m_pTTree_2->Branch("pzc_before_fit_woNu",&m_pzc_before_fit_woNu,"pzc_before_fit_woNu/F") ;
	m_pTTree_2->Branch("ec_before_fit_woNu",&m_ec_before_fit_woNu,"ec_before_fit_woNu/F") ;
	m_pTTree_2->Branch("zc_before_fit_woNu",&m_zc_before_fit_woNu,"zc_before_fit_woNu/F") ;
	m_pTTree_2->Branch("pxc_after_fit_woNu",&m_pxc_after_fit_woNu,"pxc_after_fit_woNu/F") ;
	m_pTTree_2->Branch("pyc_after_fit_woNu",&m_pyc_after_fit_woNu,"pyc_after_fit_woNu/F") ;
	m_pTTree_2->Branch("pzc_after_fit_woNu",&m_pzc_after_fit_woNu,"pzc_after_fit_woNu/F") ;
	m_pTTree_2->Branch("ec_after_fit_woNu",&m_ec_after_fit_woNu,"ec_after_fit_woNu/F") ;
	m_pTTree_2->Branch("zc_after_fit_woNu",&m_zc_after_fit_woNu,"zc_after_fit_woNu/F") ;
	m_pTTree_2->Branch("pxc_before_ISR_best",&m_pxc_before_ISR_best,"pxc_before_ISR_best/F") ;
	m_pTTree_2->Branch("pyc_before_ISR_best",&m_pyc_before_ISR_best,"pyc_before_ISR_best/F") ;
	m_pTTree_2->Branch("pzc_before_ISR_best",&m_pzc_before_ISR_best,"pzc_before_ISR_best/F") ;
	m_pTTree_2->Branch("ec_before_ISR_best",&m_ec_before_ISR_best,"ec_before_ISR_best/F") ;
	m_pTTree_2->Branch("zc_before_ISR_best",&m_zc_before_ISR_best,"zc_before_ISR_best/F") ;
	m_pTTree_2->Branch("pxc_before_fit_best",&m_pxc_before_fit_best,"pxc_before_fit_best/F") ;
	m_pTTree_2->Branch("pyc_before_fit_best",&m_pyc_before_fit_best,"pyc_before_fit_best/F") ;
	m_pTTree_2->Branch("pzc_before_fit_best",&m_pzc_before_fit_best,"pzc_before_fit_best/F") ;
	m_pTTree_2->Branch("ec_before_fit_best",&m_ec_before_fit_best,"ec_before_fit_best/F") ;
	m_pTTree_2->Branch("zc_before_fit_best",&m_zc_before_fit_best,"zc_before_fit_best/F") ;
	m_pTTree_2->Branch("pxc_after_fit_best",&m_pxc_after_fit_best,"pxc_after_fit_best/F") ;
	m_pTTree_2->Branch("pyc_after_fit_best",&m_pyc_after_fit_best,"pyc_after_fit_best/F") ;
	m_pTTree_2->Branch("pzc_after_fit_best",&m_pzc_after_fit_best,"pzc_after_fit_best/F") ;
	m_pTTree_2->Branch("ec_after_fit_best",&m_ec_after_fit_best,"ec_after_fit_best/F") ;
	m_pTTree_2->Branch("zc_after_fit_best",&m_zc_after_fit_best,"zc_after_fit_best/F") ;
	m_pTTree_2->Branch("mcHmumuPx",&m_mcHmumuPx,"mcHmumuPx/F") ;
	m_pTTree_2->Branch("mcHmumuPy",&m_mcHmumuPy,"mcHmumuPy/F") ;
	m_pTTree_2->Branch("mcHmumuPz",&m_mcHmumuPz,"mcHmumuPz/F") ;
	m_pTTree_2->Branch("mcHmumuPt",&m_mcHmumuPt,"mcHmumuPt/F") ;
	m_pTTree_2->Branch("mcHmumuP",&m_mcHmumuP,"mcHmumuP/F") ;
	m_pTTree_2->Branch("mcHmumuE",&m_mcHmumuE,"mcHmumuE/F") ;
	m_pTTree_2->Branch("mcISR1Px",&m_mcISR1Px,"mcISR1Px/F") ;
	m_pTTree_2->Branch("mcISR1Py",&m_mcISR1Py,"mcISR1Py/F") ;
	m_pTTree_2->Branch("mcISR1Pz",&m_mcISR1Pz,"mcISR1Pz/F") ;
	m_pTTree_2->Branch("mcISR1E",&m_mcISR1E,"mcISR1E/F") ;
	m_pTTree_2->Branch("mcISR2Px",&m_mcISR2Px,"mcISR2Px/F") ;
	m_pTTree_2->Branch("mcISR2Py",&m_mcISR2Py,"mcISR2Py/F") ;
	m_pTTree_2->Branch("mcISR2Pz",&m_mcISR2Pz,"mcISR2Pz/F") ;
	m_pTTree_2->Branch("mcISR2E",&m_mcISR2E,"mcISR2E/F") ;
	m_pTTree_2->Branch("mcISRPx",&m_mcISRPx,"mcISRPx/F") ;
	m_pTTree_2->Branch("mcISRPy",&m_mcISRPy,"mcISRPy/F") ;
	m_pTTree_2->Branch("mcISRPz",&m_mcISRPz,"mcISRPz/F") ;
	m_pTTree_2->Branch("mcISRE",&m_mcISRE,"mcISRE/F") ;
	m_pTTree_3 = new TTree("StartFitObjects","StartFitObjects");
	m_pTTree_3->SetDirectory(m_pTFile);
	m_pTTree_3->Branch("fitError",&m_iError_woNu) ;
	m_pTTree_3->Branch("nSLDecayTotal",&m_nSLDecayTotal,"nSLDecayTotal/I") ;
	m_pTTree_3->Branch("probability_best",&m_probability_best,"probability_best/F") ;
	m_pTTree_3->Branch("ISRE_total",&m_ISRE_total,"ISRE_total/F") ;
	m_pTTree_3->Branch("FSRE_total",&m_FSRE_total,"FSRE_total/F") ;
	m_pTTree_3->Branch("mcHmumuPx",&m_mcHmumuPx,"mcHmumuPx/F") ;
	m_pTTree_3->Branch("mcHmumuPy",&m_mcHmumuPy,"mcHmumuPy/F") ;
	m_pTTree_3->Branch("mcHmumuPz",&m_mcHmumuPz,"mcHmumuPz/F") ;
	m_pTTree_3->Branch("mcHmumuPt",&m_mcHmumuPt,"mcHmumuPt/F") ;
	m_pTTree_3->Branch("mcHmumuP",&m_mcHmumuP,"mcHmumuP/F") ;
	m_pTTree_3->Branch("mcHmumuE",&m_mcHmumuE,"mcHmumuE/F") ;
	m_pTTree_3->Branch("jet_startPx_wNu",&m_jet_startPx_wNu);
	m_pTTree_3->Branch("jet_startPy_wNu",&m_jet_startPy_wNu);
	m_pTTree_3->Branch("jet_startPz_wNu",&m_jet_startPz_wNu);
	m_pTTree_3->Branch("jet_startE_wNu",&m_jet_startE_wNu);
	m_pTTree_3->Branch("jet_startTheta_wNu",&m_jet_startTheta_wNu);
	m_pTTree_3->Branch("jet_startPhi_wNu",&m_jet_startPhi_wNu);
	m_pTTree_3->Branch("dijet_angle_wNu",&m_dijet_angle_wNu);
	m_pTTree_3->Branch("jet_startPx_wNu_bestfit",&m_jet_startPx_wNu_bestfit);
	m_pTTree_3->Branch("jet_startPy_wNu_bestfit",&m_jet_startPy_wNu_bestfit);
	m_pTTree_3->Branch("jet_startPz_wNu_bestfit",&m_jet_startPz_wNu_bestfit);
	m_pTTree_3->Branch("jet_startE_wNu_bestfit",&m_jet_startE_wNu_bestfit);
	m_pTTree_3->Branch("jet_startTheta_wNu_bestfit",&m_jet_startTheta_wNu_bestfit);
	m_pTTree_3->Branch("jet_startPhi_wNu_bestfit",&m_jet_startPhi_wNu_bestfit);
	m_pTTree_3->Branch("jet_fittedE_wNu_bestfit",&m_jet_fittedE_wNu_bestfit);
	m_pTTree_3->Branch("jet_fittedTheta_wNu_bestfit",&m_jet_fittedTheta_wNu_bestfit);
	m_pTTree_3->Branch("jet_fittedPhi_wNu_bestfit",&m_jet_fittedPhi_wNu_bestfit);
	m_pTTree_3->Branch("dijet_angle_wNu_bestfit",&m_dijet_angle_wNu_bestfit,"dijet_angle_wNu_bestfit/F");
	m_pTTree_3->Branch("jet_startPx_woNu",&m_jet_startPx_woNu);
	m_pTTree_3->Branch("jet_startPy_woNu",&m_jet_startPy_woNu);
	m_pTTree_3->Branch("jet_startPz_woNu",&m_jet_startPz_woNu);
	m_pTTree_3->Branch("jet_startE_woNu",&m_jet_startE_woNu);
	m_pTTree_3->Branch("jet_startTheta_woNu",&m_jet_startTheta_woNu);
	m_pTTree_3->Branch("jet_startPhi_woNu",&m_jet_startPhi_woNu);
	m_pTTree_3->Branch("jet_fittedE_woNu",&m_jet_fittedE_woNu);
	m_pTTree_3->Branch("jet_fittedTheta_woNu",&m_jet_fittedTheta_woNu);
	m_pTTree_3->Branch("jet_fittedPhi_woNu",&m_jet_fittedPhi_woNu);
	m_pTTree_3->Branch("dijet_angle_woNu",&m_dijet_angle_woNu,"dijet_angle_woNu/F");
	m_pTTree_3->Branch("jet_startPx_best",&m_jet_startPx_best);
	m_pTTree_3->Branch("jet_startPy_best",&m_jet_startPy_best);
	m_pTTree_3->Branch("jet_startPz_best",&m_jet_startPz_best);
	m_pTTree_3->Branch("jet_startE_best",&m_jet_startE_best);
	m_pTTree_3->Branch("jet_startTheta_best",&m_jet_startTheta_best);
	m_pTTree_3->Branch("jet_startPhi_best",&m_jet_startPhi_best);
	m_pTTree_3->Branch("jet_fittedE_best",&m_jet_fittedE_best);
	m_pTTree_3->Branch("jet_fittedTheta_best",&m_jet_fittedTheta_best);
	m_pTTree_3->Branch("jet_fittedPhi_best",&m_jet_fittedPhi_best);
	m_pTTree_3->Branch("dijet_angle_best",&m_dijet_angle_best,"dijet_angle_best/F");
	m_pTTree_3->Branch("lepton_startPx_wNu",&m_lepton_startPx_wNu);
	m_pTTree_3->Branch("lepton_startPy_wNu",&m_lepton_startPy_wNu);
	m_pTTree_3->Branch("lepton_startPz_wNu",&m_lepton_startPz_wNu);
	m_pTTree_3->Branch("lepton_startE_wNu",&m_lepton_startE_wNu);
	m_pTTree_3->Branch("lepton_startPx_wNu_bestfit",&m_lepton_startPx_wNu_bestfit);
	m_pTTree_3->Branch("lepton_startPy_wNu_bestfit",&m_lepton_startPy_wNu_bestfit);
	m_pTTree_3->Branch("lepton_startPz_wNu_bestfit",&m_lepton_startPz_wNu_bestfit);
	m_pTTree_3->Branch("lepton_startE_wNu_bestfit",&m_lepton_startE_wNu_bestfit);
	m_pTTree_3->Branch("lepton_startPx_woNu",&m_lepton_startPx_woNu);
	m_pTTree_3->Branch("lepton_startPy_woNu",&m_lepton_startPy_woNu);
	m_pTTree_3->Branch("lepton_startPz_woNu",&m_lepton_startPz_woNu);
	m_pTTree_3->Branch("lepton_startE_woNu",&m_lepton_startE_woNu);
	m_pTTree_3->Branch("lepton_startPx_best",&m_lepton_startPx_best);
	m_pTTree_3->Branch("lepton_startPy_best",&m_lepton_startPy_best);
	m_pTTree_3->Branch("lepton_startPz_best",&m_lepton_startPz_best);
	m_pTTree_3->Branch("lepton_startE_best",&m_lepton_startE_best);
	m_pTTree_3->Branch("ISR_startPx_wNu",&m_ISR_startPx_wNu,"ISR_startPx_wNu/F") ;
	m_pTTree_3->Branch("ISR_startPy_wNu",&m_ISR_startPy_wNu,"ISR_startPy_wNu/F") ;
	m_pTTree_3->Branch("ISR_startPz_wNu",&m_ISR_startPz_wNu,"ISR_startPz_wNu/F") ;
	m_pTTree_3->Branch("ISR_startPx_woNu",&m_ISR_startPx_woNu,"ISR_startPx_woNu/F") ;
	m_pTTree_3->Branch("ISR_startPy_woNu",&m_ISR_startPy_woNu,"ISR_startPy_woNu/F") ;
	m_pTTree_3->Branch("ISR_startPz_woNu",&m_ISR_startPz_woNu,"ISR_startPz_woNu/F") ;
	m_pTTree_3->Branch("jet_SigmaTheta",&m_jet_SigmaTheta);
	m_pTTree_3->Branch("jet_SigmaPhi",&m_jet_SigmaPhi);
	m_pTTree_3->Branch("jet_SigmaE",&m_jet_SigmaE);
	m_pTTree_3->Branch("lepton_SigmaTheta",&m_lepton_SigmaTheta);
	m_pTTree_3->Branch("lepton_SigmaPhi",&m_lepton_SigmaPhi);
	m_pTTree_3->Branch("lepton_SigmaInvpT",&m_lepton_SigmaInvpT);
	m_pTTree_4 = new TTree("Uncertainties","Uncertainties");
	m_pTTree_4->SetDirectory(m_pTFile);
	m_pTTree_4->Branch("probability_best",&m_probability_best,"probability_best/F") ;
	m_pTTree_4->Branch("nSLDecayTotal",&m_nSLDecayTotal,"nSLDecayTotal/I") ;
	m_pTTree_4->Branch("ISRE_total",&m_ISRE_total,"ISRE_total/F") ;
	m_pTTree_4->Branch("FSRE_total",&m_FSRE_total,"FSRE_total/F") ;
	m_pTTree_4->Branch("mcHmumuPx",&m_mcHmumuPx,"mcHmumuPx/F") ;
	m_pTTree_4->Branch("mcHmumuPy",&m_mcHmumuPy,"mcHmumuPy/F") ;
	m_pTTree_4->Branch("mcHmumuPz",&m_mcHmumuPz,"mcHmumuPz/F") ;
	m_pTTree_4->Branch("mcHmumuPt",&m_mcHmumuPt,"mcHmumuPt/F") ;
	m_pTTree_4->Branch("mcHmumuP",&m_mcHmumuP,"mcHmumuP/F") ;
	m_pTTree_4->Branch("mcHmumuE",&m_mcHmumuE,"mcHmumuE/F") ;
	m_pTTree_4->Branch("jet_startE_best",&m_jet_startE_best);
	m_pTTree_4->Branch("jet_startTheta_best",&m_jet_startTheta_best);
	m_pTTree_4->Branch("jet_startPhi_best",&m_jet_startPhi_best);
	m_pTTree_4->Branch("SigmaPx2",&m_SigmaPx2);
	m_pTTree_4->Branch("SigmaPxSigmaPy",&m_SigmaPxPy);
	m_pTTree_4->Branch("SigmaPxSigmaPz",&m_SigmaPxPz);
	m_pTTree_4->Branch("SigmaPxSigmaE",&m_SigmaPxE);
	m_pTTree_4->Branch("SigmaPy2",&m_SigmaPy2);
	m_pTTree_4->Branch("SigmaPySigmaPz",&m_SigmaPyPz);
	m_pTTree_4->Branch("SigmaPySigmaE",&m_SigmaPyE);
	m_pTTree_4->Branch("SigmaPz2",&m_SigmaPz2);
	m_pTTree_4->Branch("SigmaPzSigmaE",&m_SigmaPzE);
	m_pTTree_4->Branch("SigmaE2",&m_SigmaE2);
	m_pTTree_4->Branch("TotalSigmaPx2",&m_TotalSigmaPx2,"TotalSigmaPx2/F");
	m_pTTree_4->Branch("TotalSigmaPy2",&m_TotalSigmaPy2,"TotalSigmaPy2/F");
	m_pTTree_4->Branch("TotalSigmaPz2",&m_TotalSigmaPz2,"TotalSigmaPz2/F");
	m_pTTree_4->Branch("TotalSigmaE2",&m_TotalSigmaE2,"TotalSigmaE2/F");
	m_pTTree_4->Branch("TotalSigmaPxPy",&m_TotalSigmaPxPy,"TotalSigmaPxPy/F");
	m_pTTree_4->Branch("TotalSigmaPxPz",&m_TotalSigmaPxPz,"TotalSigmaPxPz/F");
	m_pTTree_4->Branch("TotalSigmaPxE",&m_TotalSigmaPxE,"TotalSigmaPxE/F");
	m_pTTree_4->Branch("TotalSigmaPyPz",&m_TotalSigmaPyPz,"TotalSigmaPyPz/F");
	m_pTTree_4->Branch("TotalSigmaPyE",&m_TotalSigmaPyE,"TotalSigmaPyE/F");
	m_pTTree_4->Branch("TotalSigmaPzE",&m_TotalSigmaPzE,"TotalSigmaPzE/F");
	m_pTTree_5 = new TTree("ISR_FSR","ISR_FSR");
	m_pTTree_5->SetDirectory(m_pTFile);
	m_pTTree_5->Branch("ISRPx",&m_ISRPx);
	m_pTTree_5->Branch("ISRPy",&m_ISRPy);
	m_pTTree_5->Branch("ISRPz",&m_ISRPz);
	m_pTTree_5->Branch("ISRPt",&m_ISRPt);
	m_pTTree_5->Branch("ISRE",&m_ISRE);
	m_pTTree_5->Branch("ISRPx_total",&m_ISRPx_total,"ISRPx_total/F") ;
	m_pTTree_5->Branch("ISRPy_total",&m_ISRPy_total,"ISRPy_total/F") ;
	m_pTTree_5->Branch("ISRPz_total",&m_ISRPz_total,"ISRPz_total/F") ;
	m_pTTree_5->Branch("ISRPt_total",&m_ISRPt_total,"ISRPt_total/F") ;
	m_pTTree_5->Branch("ISRE_total",&m_ISRE_total,"ISRE_total/F") ;
	m_pTTree_5->Branch("FSRPx",&m_FSRPx);
	m_pTTree_5->Branch("FSRPy",&m_FSRPy);
	m_pTTree_5->Branch("FSRPz",&m_FSRPz);
	m_pTTree_5->Branch("FSRPt",&m_FSRPt);
	m_pTTree_5->Branch("FSRE",&m_FSRE);
	m_pTTree_5->Branch("FSRPx_total",&m_FSRPx_total,"FSRPx_total/F") ;
	m_pTTree_5->Branch("FSRPy_total",&m_FSRPy_total,"FSRPy_total/F") ;
	m_pTTree_5->Branch("FSRPz_total",&m_FSRPz_total,"FSRPz_total/F") ;
	m_pTTree_5->Branch("FSRPt_total",&m_FSRPt_total,"FSRPt_total/F") ;
	m_pTTree_5->Branch("FSRE_total",&m_FSRE_total,"FSRE_total/F") ;
	m_pTTree_6 = new TTree("ErrorFlow","ErrorFlow");
	m_pTTree_6->SetDirectory(m_pTFile);
	m_pTTree_6->Branch("probability_best",&m_probability_best,"probability_best/F") ;
	m_pTTree_6->Branch("nSLDecayTotal",&m_nSLDecayTotal,"nSLDecayTotal/I") ;
	m_pTTree_6->Branch("ISRE_total",&m_ISRE_total,"ISRE_total/F") ;
	m_pTTree_6->Branch("FSRE_total",&m_FSRE_total,"FSRE_total/F") ;
	m_pTTree_6->Branch("mcHmumuPt",&m_mcHmumuPt,"mcHmumuPt/F") ;
	m_pTTree_6->Branch("mcHmumuP",&m_mcHmumuP,"mcHmumuP/F") ;
	m_pTTree_6->Branch("mcHmumuE",&m_mcHmumuE,"mcHmumuE/F") ;
	m_pTTree_6->Branch("Theta_mcQuark",&m_Theta_mcQuark);
	m_pTTree_6->Branch("Phi_mcQuark",&m_Phi_mcQuark);
	m_pTTree_6->Branch("Energy_mcQuark",&m_Energy_mcQuark);
	m_pTTree_6->Branch("Theta_mcJet_wInvisibles",&m_Theta_mcJet_wInvisibles);
	m_pTTree_6->Branch("Phi_mcJet_wInvisibles",&m_Phi_mcJet_wInvisibles);
	m_pTTree_6->Branch("Energy_mcJet_wInvisibles",&m_Energy_mcJet_wInvisibles);
	m_pTTree_6->Branch("Theta_mcJet_woInvisibles",&m_Theta_mcJet_woInvisibles);
	m_pTTree_6->Branch("Phi_mcJet_woInvisibles",&m_Phi_mcJet_woInvisibles);
	m_pTTree_6->Branch("Energy_mcJet_woInvisibles",&m_Energy_mcJet_woInvisibles);
	m_pTTree_6->Branch("Theta_recoJet",&m_Theta_recoJet);
	m_pTTree_6->Branch("Phi_recoJet",&m_Phi_recoJet);
	m_pTTree_6->Branch("Energy_recoJet",&m_Energy_recoJet);
	m_pTTree_6->Branch("nPFOs_NeutralHadrons",&m_nPFOs_NeutralHadrons);
	m_pTTree_6->Branch("nPFOs_Photons",&m_nPFOs_Photons);
	m_pTTree_6->Branch("nPFOs_ChargedPFOs",&m_nPFOs_ChargedPFOs);
	m_pTTree_6->Branch("nPFOs_Total",&m_nPFOs_Total);
	m_pTTree_6->Branch("PFOEnergy_NeutralHadrons",&m_PFOEnergy_NeutralHadrons);
	m_pTTree_6->Branch("PFOEnergy_Photons",&m_PFOEnergy_Photons);
	m_pTTree_6->Branch("PFOEnergy_ChargedPFOs",&m_PFOEnergy_ChargedPFOs);
	m_pTTree_6->Branch("PFOEnergy_Total",&m_PFOEnergy_Total);
	m_pTTree_6->Branch("EnergyFraction_NeutralHadrons",&m_EnergyFraction_NeutralHadrons);
	m_pTTree_6->Branch("EnergyFraction_Photons",&m_EnergyFraction_Photons);
	m_pTTree_6->Branch("EnergyFraction_ChargedPFOs",&m_EnergyFraction_ChargedPFOs);
	m_pTTree_6->Branch("ResidualTheta_mcQuark_mcJet",&m_ResidualTheta_mcQuark_mcJet);
	m_pTTree_6->Branch("ResidualPhi_mcQuark_mcJet",&m_ResidualPhi_mcQuark_mcJet);
	m_pTTree_6->Branch("ResidualEnergy_mcQuark_mcJet",&m_ResidualEnergy_mcQuark_mcJet);
	m_pTTree_6->Branch("ResidualTheta_mcJet_recoJet",&m_ResidualTheta_mcJet_recoJet);
	m_pTTree_6->Branch("ResidualPhi_mcJet_recoJet",&m_ResidualPhi_mcJet_recoJet);
	m_pTTree_6->Branch("ResidualEnergy_mcJet_recoJet",&m_ResidualEnergy_mcJet_recoJet);
	m_pTTree_6->Branch("ResidualTheta_mcQuark_recoJet",&m_ResidualTheta_mcQuark_recoJet);
	m_pTTree_6->Branch("ResidualPhi_mcQuark_recoJet",&m_ResidualPhi_mcQuark_recoJet);
	m_pTTree_6->Branch("ResidualEnergy_mcQuark_recoJet",&m_ResidualEnergy_mcQuark_recoJet);
	m_pTTree_6->Branch("NormalizedResidualTheta_mcQuark_mcJet",&m_NormalizedResidualTheta_mcQuark_mcJet);
	m_pTTree_6->Branch("NormalizedResidualPhi_mcQuark_mcJet",&m_NormalizedResidualPhi_mcQuark_mcJet);
	m_pTTree_6->Branch("NormalizedResidualEnergy_mcQuark_mcJet",&m_NormalizedResidualEnergy_mcQuark_mcJet);
	m_pTTree_6->Branch("NormalizedResidualTheta_mcJet_recoJet",&m_NormalizedResidualTheta_mcJet_recoJet);
	m_pTTree_6->Branch("NormalizedResidualPhi_mcJet_recoJet",&m_NormalizedResidualPhi_mcJet_recoJet);
	m_pTTree_6->Branch("NormalizedResidualEnergy_mcJet_recoJet",&m_NormalizedResidualEnergy_mcJet_recoJet);
	m_pTTree_6->Branch("NormalizedResidualTheta_mcQuark_recoJet",&m_NormalizedResidualTheta_mcQuark_recoJet);
	m_pTTree_6->Branch("NormalizedResidualPhi_mcQuark_recoJet",&m_NormalizedResidualPhi_mcQuark_recoJet);
	m_pTTree_6->Branch("NormalizedResidualEnergy_mcQuark_recoJet",&m_NormalizedResidualEnergy_mcQuark_recoJet);
	m_pTTree_6->Branch("JetErrorTheta_ErrorFlow",&m_JetErrorTheta_ErrorFlow);
	m_pTTree_6->Branch("JetErrorPhi_ErrorFlow",&m_JetErrorPhi_ErrorFlow);
	m_pTTree_6->Branch("JetErrorEnergy_ErrorFlow",&m_JetErrorEnergy_ErrorFlow);

	h_HDecayMode = new TH1I("h_HDecayMode", "Decay mode of Higgs boson", 4, 0, 4);
	h_HDecayMode->GetXaxis()->SetBinLabel(1,"H #rightarrow b#bar{b}");
	h_HDecayMode->GetXaxis()->SetBinLabel(2,"H #rightarrow c#bar{c}");
	h_HDecayMode->GetXaxis()->SetBinLabel(3,"H #rightarrow gg");
	h_HDecayMode->GetXaxis()->SetBinLabel(4,"H #rightarrow other");
	h_HDecayMode->SetDirectory(m_pTFile);
	h_Zmass_beforefit_woNu = new TH1F("h_Zmass_beforefit_woNu", "Z mass before fit without #nu correction", 400, 0., 200.);
	h_Zmass_beforefit_woNu->SetDirectory(m_pTFile);
	h_Hmass_beforefit_woNu = new TH1F("h_Hmass_beforefit_woNu", "H mass before fit without #nu correction", 400, 0., 200.);
	h_Hmass_beforefit_woNu->SetDirectory(m_pTFile);
	h_Zmass_beforefit_wNu = new TH1F("h_Zmass_beforefit_wNu", "Z mass before fit with #nu correction", 400, 0., 200.);
	h_Zmass_beforefit_wNu->SetDirectory(m_pTFile);
	h_Hmass_beforefit_wNu = new TH1F("h_Hmass_beforefit_wNu", "H mass before fit with #nu correction", 400, 0., 200.);
	h_Hmass_beforefit_wNu->SetDirectory(m_pTFile);
	h_Zmass_afterfit_woNu = new TH1F("h_Zmass_afterfit_woNu", "Z mass after fit without #nu correction", 400, 0., 200.);
	h_Zmass_afterfit_woNu->SetDirectory(m_pTFile);
	h_Hmass_afterfit_woNu = new TH1F("h_Hmass_afterfit_woNu", "H mass after fit without #nu correction", 400, 0., 200.);
	h_Hmass_afterfit_woNu->SetDirectory(m_pTFile);
	h_Zmass_afterfit_wNu = new TH1F("h_Zmass_afterfit_wNu", "Z mass after fit with #nu correction", 400, 0., 200.);
	h_Zmass_afterfit_wNu->SetDirectory(m_pTFile);
	h_Hmass_afterfit_wNu = new TH1F("h_Hmass_afterfit_wNu", "H mass after fit with #nu correction", 400, 0., 200.);
	h_Hmass_afterfit_wNu->SetDirectory(m_pTFile);
	h_fitError_wNu = new TH1I("h_fitError_wNu", "fit error with #nu correction", 5, -1.5, 3.5);
	h_fitError_wNu->SetDirectory(m_pTFile);
	h_fitError_woNu = new TH1I("h_fitError_woNu", "fit error without #nu correction", 5, -1.5, 3.5);
	h_fitError_woNu->SetDirectory(m_pTFile);
	h_ErrorCode_wNu_woNu = new TH2I("h_ErrorCode_wNu_woNu", "Error Code; Error Code_{with #nu}; Error Code_{without #nu}", 5, -1.5, 3.5, 5, -1.5, 3.5);
	h_ErrorCode_wNu_woNu->SetDirectory(m_pTFile);
	h_fitProbability_wNu = new TH1F("h_fitProbability_wNu", "fit probability with #nu correction", 100, 0., 1.);
	h_fitProbability_wNu->SetDirectory(m_pTFile);
	h_fitProbability_woNu = new TH1F("h_fitProbability_woNu", "fit probability without #nu correction", 100, 0., 1.);
	h_fitProbability_woNu->SetDirectory(m_pTFile);
	h_fitProbability_best = new TH1F("h_fitProbability_best", "best fit probability", 100, 0., 1.);
	h_fitProbability_best->SetDirectory(m_pTFile);
	h_nJets = new TH1I("h_nJets", "number of found jets per event; n_{Jets}; n_{events}", 8, -0.5, 7.5);
	h_nJets->SetDirectory(m_pTFile);
	h_nLeptons = new TH1I("h_nLeptons", "number of found isolated leptons per event; n_{ISOLeptons}; n_{events}", 8, -0.5, 7.5);
	h_nLeptons->SetDirectory(m_pTFile);
	h_nLeptons_nJets = new TH2I("h_nLeptons_nJets", "nISOleptons vs nJets per event; n_{ISOLeptons}; n_{Jets}", 8, -0.5, 7.5, 8, -0.5, 7.5);
	h_nLeptons_nJets->SetDirectory(m_pTFile);
	h_ISRE_1mcp_fit = new TH2F("h_ISRE_1mcp_fit", "ISR energy MCP(p_{z}^{max}) vs FIT; ISR_{MCP} E [GeV]; ISR_{FIT} E [GeV]", 100, 0., m_isrpzmax * 1.5, 100, 0., m_isrpzmax * 1.5);
	h_ISRE_1mcp_fit->SetDirectory(m_pTFile);
	h_ISRE_2mcp_fit = new TH2F("h_ISRE_2mcp_fit", "ISR energy MCP(#Sigma E_{ISR}) vs FIT; ISR_{MCP} E [GeV]; ISR_{FIT} E [GeV]", 100, 0., m_isrpzmax * 1.5, 100, 0., m_isrpzmax * 1.5);
	h_ISRE_2mcp_fit->SetDirectory(m_pTFile);
	h_ISRpz_1mcp_fit = new TH2F("h_ISRpz_1mcp_fit", "ISR p_{z} MCP(p_{z}^{max}) vs FIT; ISR_{MCP} p_{z} [GeV]; ISR_{FIT} p_{z} [GeV]", 100, -1.5 * m_isrpzmax , 1.5 * m_isrpzmax , 100, -1.5 * m_isrpzmax , 1.5 * m_isrpzmax);
	h_ISRpz_1mcp_fit->SetDirectory(m_pTFile);
	h_ISRpz_2mcp_fit = new TH2F("h_ISRpz_2mcp_fit", "ISR p_{z} MCP(#Sigma #vec{p}_{z}) vs FIT; ISR_{MCP} p_{z} [GeV]; ISR_{FIT} p_{z} [GeV]", 100, -1.5 * m_isrpzmax , 1.5 * m_isrpzmax , 100, -1.5 * m_isrpzmax , 1.5 * m_isrpzmax);
	h_ISRpz_2mcp_fit->SetDirectory(m_pTFile);
	h_ISR_pzc_wNu = new TH2F("h_ISR_pzc_wNu", "ISR momentum vs p_{z} constraint with #nu correction; #Sigma p_{z} [GeV]; p_{z}^{ISR} [GeV]", 80, -40., 40., 80, -40., 40.);
	h_ISR_pzc_wNu->SetDirectory(m_pTFile);
	h_ISR_pzc_woNu = new TH2F("h_ISR_pzc_woNu", "ISR momentum vs p_{z} constraint without #nu correction; #Sigma p_{z} [GeV]; p_{z}^{ISR} [GeV]", 80, -40., 40., 80, -40., 40.);
	h_ISR_pzc_woNu->SetDirectory(m_pTFile);
	h_ISR_pzc_best = new TH2F("h_ISR_pzc_best", "ISR momentum vs p_{z} constraint (best fit); #Sigma p_{z} [GeV]; p_{z}^{ISR} [GeV]", 80, -40., 40., 80, -40., 40.);
	h_ISR_pzc_best->SetDirectory(m_pTFile);
	h_mcpISRPx = new TH2F("h_mcpISRPx", "ISR Px; MCP ISR_{1} P_{x} [GeV]; MCP ISR_{2} P_{x} [GeV]", 1000, -m_isrpzmax * 1.5, m_isrpzmax * 1.5, 1000, -m_isrpzmax * 1.5, m_isrpzmax * 1.5);
	h_mcpISRPx->SetDirectory(m_pTFile);
	h_mcpISRPy = new TH2F("h_mcpISRPy", "ISR Py; MCP ISR_{1} P_{y} [GeV]; MCP ISR_{2} P_{y} [GeV]", 1000, -m_isrpzmax * 1.5, m_isrpzmax * 1.5, 1000, -m_isrpzmax * 1.5, m_isrpzmax * 1.5);
	h_mcpISRPy->SetDirectory(m_pTFile);
	h_mcpISRPz = new TH2F("h_mcpISRPz", "ISR Pz; MCP ISR_{1} P_{z} [GeV]; MCP ISR_{2} P_{z} [GeV]", 1000, -m_isrpzmax * 1.5, m_isrpzmax * 1.5, 1000, -m_isrpzmax * 1.5, m_isrpzmax * 1.5);
	h_mcpISRPz->SetDirectory(m_pTFile);
	h_mcpISRE = new TH2F("h_mcpISRE", "ISR Energy; MCP ISR_{1} Energy [GeV]; MCP ISR_{2} Energy [GeV]", 1000, 0., m_isrpzmax * 1.5, 1000, 0., m_isrpzmax * 1.5);
	h_mcpISRE->SetDirectory(m_pTFile);
	h_pull_jet_E = new TH1F("h_pull_jet_E", "pull of E for jets after successful fit; pull E_{jet} [GeV]; n_{jets}", 100, -10., 10.);
	h_pull_jet_E->SetDirectory(m_pTFile);
	h_pull_jet_theta = new TH1F("h_pull_jet_theta", "pull of #theta for jets after successful fit; pull #theta_{jet}; n_{jets}", 100, -10., 10.);
	h_pull_jet_theta->SetDirectory(m_pTFile);
	h_pull_jet_phi = new TH1F("h_pull_jet_phi", "pull of #phi for jets after successful fit; pull #phi_{jet}; n_{jets}", 100, -10., 10.);
	h_pull_jet_phi->SetDirectory(m_pTFile);
	h_pull_lepton_InvPt = new TH1F("h_pull_lepton_InvPt", "pull of #frac{1}{p_{T}} for leptons after successful fit; pull #frac{1}{p_{T}}_{lepton}; n_{leptons}", 100, -10., 10.);
	h_pull_lepton_InvPt->SetDirectory(m_pTFile);
	h_pull_lepton_theta = new TH1F("h_pull_lepton_theta", "pull of #theta for leptons after successful fit; pull #theta_{lepton}; n_{leptons}", 100, -10., 10.);
	h_pull_lepton_theta->SetDirectory(m_pTFile);
	h_pull_lepton_phi = new TH1F("h_pull_lepton_phi", "pull of #phi for leptons after successful fit; pull #phi_{lepton}; n_{leptons}", 100, -10., 10.);
	h_pull_lepton_phi->SetDirectory(m_pTFile);
	h_error_jet_E = new TH1F("h_error_jet_E", "; #sigma_{E_{jet}}; n_{jets}", 10000, 0., 100.);
	h_error_jet_E->SetDirectory(m_pTFile);
	h_error_jet_Theta = new TH1F("h_error_jet_Theta", "; #sigma_{#theta_{jet}}; n_{jet}", 10000, 0., 0.1);
	h_error_jet_Theta->SetDirectory(m_pTFile);
	h_error_jet_Phi = new TH1F("h_error_jet_Phi", "; #sigma_{#phi_{jet}}; n_{jets}", 10000, 0., 0.1);
	h_error_jet_Phi->SetDirectory(m_pTFile);
	h_error_lepton_InvpT = new TH1F("h_error_lepton_InvpT", "; #sigma_{#frac{1}{p_{T}}}; n_{lepton}", 10000, 0., 0.01);
	h_error_lepton_InvpT->SetDirectory(m_pTFile);
	h_error_lepton_Theta = new TH1F("h_error_lepton_Theta", "; #sigma_{#theta_{lepton}}; n_{lepton}", 10000, 0., 0.01);
	h_error_lepton_Theta->SetDirectory(m_pTFile);
	h_error_lepton_Phi = new TH1F("h_error_lepton_Phi", "; #sigma_{#phi_{lepton}}; n_{lepton}", 10000, 0., 0.01);
	h_error_lepton_Phi->SetDirectory(m_pTFile);
	h_SigmaPx2 = new TH1F("h_SigmaPx2", "; #sigma^{2}_{p_{x}}; n_{jet}", 400, 0., 0.01);
	h_SigmaPx2->SetDirectory(m_pTFile);
	h_SigmaPy2 = new TH1F("h_SigmaPy2", "; #sigma^{2}_{p_{y}}; n_{jet}", 400, 0., 0.01);
	h_SigmaPy2->SetDirectory(m_pTFile);
	h_SigmaPz2 = new TH1F("h_SigmaPz2", "; #sigma^{2}_{p_{z}}; n_{jet}", 400, 0., 0.01);
	h_SigmaPz2->SetDirectory(m_pTFile);
	h_SigmaE2 = new TH1F("h_SigmaE2", "; #sigma^{2}_{E}; n_{jet}", 400, 0., 50);
	h_SigmaE2->SetDirectory(m_pTFile);
	h_SigmaPxPy = new TH1F("h_SigmaPxPy", "; #sigma_{p_{x}p_{y}}; n_{jet}", 400, -0.01, 0.01);
	h_SigmaPxPy->SetDirectory(m_pTFile);
	h_SigmaPxPz = new TH1F("h_SigmaPxPz", "; #sigma_{p_{x}p_{z}}; n_{jet}", 400, -0.01, 0.01);
	h_SigmaPxPz->SetDirectory(m_pTFile);
	h_SigmaPxE = new TH1F("h_SigmaPxE", "; #sigma_{p_{x}E}; n_{jet}", 400, -0.01, 0.01);
	h_SigmaPxE->SetDirectory(m_pTFile);
	h_SigmaPyPz = new TH1F("h_SigmaPyPz", "; #sigma_{p_{y}p_{z}}; n_{jet}", 400, -0.01, 0.01);
	h_SigmaPyPz->SetDirectory(m_pTFile);
	h_SigmaPyE = new TH1F("h_SigmaPyE", "; #sigma_{p_{y}E}; n_{jet}", 400, -0.01, 0.01);
	h_SigmaPyE->SetDirectory(m_pTFile);
	h_SigmaPzE = new TH1F("h_SigmaPzE", "; #sigma(p_{z}E); n_{jet}", 400, -0.01, 0.01);
	h_SigmaPzE->SetDirectory(m_pTFile);
	h_fitProbability_diJetAngle = new TH2F("h_fitProbability_diJetAngle", "fit probability vs di-jet angle (best fit); di-jet angle [deg]; fit probability", 90, 0., 180., 200, 0., 1.0);
	h_fitProbability_diJetAngle->SetDirectory(m_pTFile);
	h_fitProbability_Ejet = new TH2F("h_fitProbability_Ejet", "fit probability vs jet energy (best fit); jet energy [GeV]; fit probability", 75, 0., 150., 200, 0., 1.0);
	h_fitProbability_Ejet->SetDirectory(m_pTFile);
	h_fitProbability_Thetajet = new TH2F("h_fitProbability_Thetajet", "fit probability vs jet theta (best fit); #theta_{jet} [deg]; fit probability", 90, 0., 180., 200, 0., 1.0);
	h_fitProbability_Thetajet->SetDirectory(m_pTFile);
	h_fitProbability_Phijet = new TH2F("h_fitProbability_Phijet", "fit probability vs jet phi (best fit); #phi_{jet} [deg]; fit probability", 90, -180., 180., 200, 0., 1.0);
	h_fitProbability_Phijet->SetDirectory(m_pTFile);
	h_fitProbability_SigmaEjet = new TH2F("h_fitProbability_SigmaEjet", "fit probability vs #sigma(E_{jet}) (best fit); #sigma(E_{jet}) [GeV]; fit probability", 50, 0., 50., 200, 0., 1.0);
	h_fitProbability_SigmaEjet->SetDirectory(m_pTFile);
	h_fitProbability_SigmaThetajet = new TH2F("h_fitProbability_SigmaThetajet", "fit probability vs #sigma(#theta_{jet}) (best fit); #sigma(#theta_{jet}) [radian]; fit probability", 1000, 0., m_JetResThetaCoeff * 0.2, 200, 0., 1.0);
	h_fitProbability_SigmaThetajet->SetDirectory(m_pTFile);
	h_fitProbability_SigmaPhijet = new TH2F("h_fitProbability_SigmaPhijet", "fit probability vs #sigma(#phi_{jet}) (best fit); #sigma(#phi_{jet}) [radian]; fit probability", 1000, 0., m_JetResPhiCoeff * 0.2, 200, 0., 1.0);
	h_fitProbability_SigmaPhijet->SetDirectory(m_pTFile);
	h_fitProbability_pullEjet = new TH2F("h_fitProbability_pullEjet", "fit probability vs pull E_{jet} (best fit); pull E_{jet}; fit probability", 100, -10., 10., 200, 0., 1.0);
	h_fitProbability_pullEjet->SetDirectory(m_pTFile);
	h_fitProbability_pullThetajet = new TH2F("h_fitProbability_pullThetajet", "fit probability vs pull #theta_{jet} (best fit); pull #theta_{jet}; fit probability", 100, -10., 10., 200, 0., 1.0);
	h_fitProbability_pullThetajet->SetDirectory(m_pTFile);
	h_fitProbability_pullPhijet = new TH2F("h_fitProbability_pullPhijet", "fit probability vs pull #phi_{jet} (best fit); pull #phi_{jet}; fit probability", 100, -10., 10., 200, 0., 1.0);
	h_fitProbability_pullPhijet->SetDirectory(m_pTFile);
	h_fitProbability_constraintPx = new TH2F("h_fitProbability_constraintPx", "fit probability vs #Sigma p_{x} (best fit); #Sigma p_{x} [GeV]; fit probability", 100, -50., 50., 200, 0., 1.0);
	h_fitProbability_constraintPx->SetDirectory(m_pTFile);
	h_fitProbability_constraintPy = new TH2F("h_fitProbability_constraintPy", "fit probability vs #Sigma p_{y} (best fit); #Sigma p_{y} [GeV]; fit probability", 100, -50., 50., 200, 0., 1.0);
	h_fitProbability_constraintPy->SetDirectory(m_pTFile);
	h_fitProbability_constraintPz = new TH2F("h_fitProbability_constraintPz", "fit probability vs #Sigma p_{z} (best fit); #Sigma p_{z} [GeV]; fit probability", 100, -50., 50., 200, 0., 1.0);
	h_fitProbability_constraintPz->SetDirectory(m_pTFile);
	h_fitProbability_constraintE = new TH2F("h_fitProbability_constraintE", "fit probability vs #Sigma E (best fit); #Sigma E [GeV]; fit probability", 100, -50., 50., 200, 0., 1.0);
	h_fitProbability_constraintE->SetDirectory(m_pTFile);
	h_fitProbability_sigmaPx = new TH2F("h_fitProbability_sigmaPx", ";total #sigma_{p_{x}} [GeV]; fit probability", 200, 0., 5., 200, 0., 1.0);
	h_fitProbability_sigmaPx->SetDirectory(m_pTFile);
	h_fitProbability_sigmaPy = new TH2F("h_fitProbability_sigmaPy", ";total #sigma_{p_{y}} [GeV]; fit probability", 200, 0., 5., 200, 0., 1.0);
	h_fitProbability_sigmaPy->SetDirectory(m_pTFile);
	h_fitProbability_sigmaPz = new TH2F("h_fitProbability_sigmaPz", ";total #sigma_{p_{z}} [GeV]; fit probability", 200, 0., 5., 200, 0., 1.0);
	h_fitProbability_sigmaPz->SetDirectory(m_pTFile);
	h_fitProbability_sigmaE = new TH2F("h_fitProbability_sigmaE", ";total #sigma_{E} [GeV]; fit probability", 200, 0., 10., 200, 0., 1.0);
	h_fitProbability_sigmaE->SetDirectory(m_pTFile);
	h_constraintPx_lowFitProb = new TH1F("h_constraintPx_lowFitProb", "fit Probability < 0.1; #Sigma p_{x} [GeV]; n_{events}", 400, -50.0, 50.0);
	h_constraintPx_lowFitProb->SetDirectory(m_pTFile);
	h_constraintPx_midFitProb = new TH1F("h_constraintPx_midFitProb", "0.1 < fit Probability < 0.9; #Sigma p_{x} [GeV]; n_{events}", 400, -50.0, 50.0);
	h_constraintPx_midFitProb->SetDirectory(m_pTFile);
	h_constraintPx_highFitProb = new TH1F("h_constraintPx_highFitProb", "0.9 < fit Probability; #Sigma p_{x} [GeV]; n_{events}", 400, -50.0, 50.0);
	h_constraintPx_highFitProb->SetDirectory(m_pTFile);
	h_constraintPy_lowFitProb = new TH1F("h_constraintPy_lowFitProb", "fit Probability < 0.1; #Sigma p_{y} [GeV]; n_{events}", 400, -50.0, 50.0);
	h_constraintPy_lowFitProb->SetDirectory(m_pTFile);
	h_constraintPy_midFitProb = new TH1F("h_constraintPy_midFitProb", "0.1 < fit Probability < 0.9; #Sigma p_{y} [GeV]; n_{events}", 400, -50.0, 50.0);
	h_constraintPy_midFitProb->SetDirectory(m_pTFile);
	h_constraintPy_highFitProb = new TH1F("h_constraintPy_highFitProb", "0.9 < fit Probability; #Sigma p_{y} [GeV]; n_{events}", 400, -50.0, 50.0);
	h_constraintPy_highFitProb->SetDirectory(m_pTFile);
	h_constraintPz_lowFitProb = new TH1F("h_constraintPz_lowFitProb", "fit Probability < 0.1; #Sigma p_{z} [GeV]; n_{events}", 400, -50.0, 50.0);
	h_constraintPz_lowFitProb->SetDirectory(m_pTFile);
	h_constraintPz_midFitProb = new TH1F("h_constraintPz_midFitProb", "0.1 < fit Probability < 0.9; #Sigma p_{z} [GeV]; n_{events}", 400, -50.0, 50.0);
	h_constraintPz_midFitProb->SetDirectory(m_pTFile);
	h_constraintPz_highFitProb = new TH1F("h_constraintPz_highFitProb", "0.9 < fit Probability; #Sigma p_{z} [GeV]; n_{events}", 400, -50.0, 50.0);
	h_constraintPz_highFitProb->SetDirectory(m_pTFile);
	h_constraintE_lowFitProb = new TH1F("h_constraintE_lowFitProb", "fit Probability < 0.1; #Sigma E [GeV]; n_{events}", 400, -50.0, 50.0);
	h_constraintE_lowFitProb->SetDirectory(m_pTFile);
	h_constraintE_midFitProb = new TH1F("h_constraintE_midFitProb", "0.1 < fit Probability < 0.9; #Sigma E [GeV]; n_{events}", 400, -50.0, 50.0);
	h_constraintE_midFitProb->SetDirectory(m_pTFile);
	h_constraintE_highFitProb = new TH1F("h_constraintE_highFitProb", "0.9 < fit Probability; #Sigma E [GeV]; n_{events}", 400, -50.0, 50.0);
	h_constraintE_highFitProb->SetDirectory(m_pTFile);
	h_constraintPx_uncertaintyPx = new TH2F("h_constraintPx_uncertaintyPx", "; #Sigma p_{x} [GeV]; #sigma_{p_{x}} [GeV]", 200, -10., 10., 200, 0., 5.);
	h_constraintPx_uncertaintyPx->SetDirectory(m_pTFile);
	h_constraintPy_uncertaintyPy = new TH2F("h_constraintPy_uncertaintyPy", "; #Sigma p_{y} [GeV]; #sigma_{p_{y}} [GeV]", 200, -10., 10., 200, 0., 5.);
	h_constraintPy_uncertaintyPy->SetDirectory(m_pTFile);
	h_constraintPz_uncertaintyPz = new TH2F("h_constraintPz_uncertaintyPz", "; #Sigma p_{z} [GeV]; #sigma_{p_{z}} [GeV]", 200, -10., 10., 200, 0., 5.);
	h_constraintPz_uncertaintyPz->SetDirectory(m_pTFile);
	h_constraintE_uncertaintyE = new TH2F("h_constraintE_uncertaintyE", "; #Sigma E [GeV]; #sigma_{E} [GeV]", 200, -50., 50., 200, 0., 10.);
	h_constraintE_uncertaintyE->SetDirectory(m_pTFile);
	h_constraintPx_uncertaintyPx_lowFitProb = new TH2F("h_constraintPx_uncertaintyPx_lowFitProb", "fit probability < 0.1; #Sigma p_{x} [GeV]; #sigma_{p_{x}} [GeV]", 200, -10., 10., 200, 0., 5.0);
	h_constraintPx_uncertaintyPx_lowFitProb->SetDirectory(m_pTFile);
	h_constraintPx_uncertaintyPx_midFitProb = new TH2F("h_constraintPx_uncertaintyPx_midFitProb", "0.1 #leq fit probability < 0.9; #Sigma p_{x} [GeV]; #sigma_{p_{x}} [GeV]", 200, -10., 10., 200, 0., 5.0);
	h_constraintPx_uncertaintyPx_midFitProb->SetDirectory(m_pTFile);
	h_constraintPx_uncertaintyPx_highFitProb = new TH2F("h_constraintPx_uncertaintyPx_highFitProb", "0.9 #leq fit probability; #Sigma p_{x} [GeV]; #sigma_{p_{x}} [GeV]", 200, -10., 10., 200, 0., 5.0);
	h_constraintPx_uncertaintyPx_highFitProb->SetDirectory(m_pTFile);
	h_constraintPy_uncertaintyPy_lowFitProb = new TH2F("h_constraintPy_uncertaintyPy_lowFitProb", "fit probability < 0.1; #Sigma p_{y} [GeV]; #sigma_{p_{y}} [GeV]", 200, -10., 10., 200, 0., 5.0);
	h_constraintPy_uncertaintyPy_lowFitProb->SetDirectory(m_pTFile);
	h_constraintPy_uncertaintyPy_midFitProb = new TH2F("h_constraintPy_uncertaintyPy_midFitProb", "0.1 #leq fit probability < 0.9; #Sigma p_{y} [GeV]; #sigma_{p_{y}} [GeV]", 200, -10., 10., 200, 0., 5.0);
	h_constraintPy_uncertaintyPy_midFitProb->SetDirectory(m_pTFile);
	h_constraintPy_uncertaintyPy_highFitProb = new TH2F("h_constraintPy_uncertaintyPy_highFitProb", "0.9 #leq fit probability; #Sigma p_{y} [GeV]; #sigma_{p_{y}} [GeV]", 200, -10., 10., 200, 0., 5.0);
	h_constraintPy_uncertaintyPy_highFitProb->SetDirectory(m_pTFile);
	h_constraintPz_uncertaintyPz_lowFitProb = new TH2F("h_constraintPz_uncertaintyPz_lowFitProb", "fit probability < 0.1; #Sigma p_{z} [GeV]; #sigma_{p_{z}} [GeV]", 200, -10., 10., 200, 0., 5.0);
	h_constraintPz_uncertaintyPz_lowFitProb->SetDirectory(m_pTFile);
	h_constraintPz_uncertaintyPz_midFitProb = new TH2F("h_constraintPz_uncertaintyPz_midFitProb", "0.1 #leq fit probability < 0.9; #Sigma p_{z} [GeV]; #sigma_{p_{z}} [GeV]", 200, -10., 10., 200, 0., 5.0);
	h_constraintPz_uncertaintyPz_midFitProb->SetDirectory(m_pTFile);
	h_constraintPz_uncertaintyPz_highFitProb = new TH2F("h_constraintPz_uncertaintyPz_highFitProb", "0.9 #leq fit probability; #Sigma p_{z} [GeV]; #sigma_{p_{z}} [GeV]", 200, -10., 10., 200, 0., 5.0);
	h_constraintPz_uncertaintyPz_highFitProb->SetDirectory(m_pTFile);
	h_constraintE_uncertaintyE_lowFitProb = new TH2F("h_constraintE_uncertaintyE_lowFitProb", "fit probability < 0.1; #Sigma E [GeV]; #sigma_{E} [GeV]", 200, -50., 50., 200, 0., 10.0);
	h_constraintE_uncertaintyE_lowFitProb->SetDirectory(m_pTFile);
	h_constraintE_uncertaintyE_midFitProb = new TH2F("h_constraintE_uncertaintyE_midFitProb", "0.1 #leq fit probability < 0.9; #Sigma E [GeV]; #sigma_{E} [GeV]", 200, -50., 50., 200, 0., 10.0);
	h_constraintE_uncertaintyE_midFitProb->SetDirectory(m_pTFile);
	h_constraintE_uncertaintyE_highFitProb = new TH2F("h_constraintE_uncertaintyE_highFitProb", "0.9 #leq fit probability; #Sigma E [GeV]; #sigma_{E} [GeV]", 200, -50., 50., 200, 0., 10.0);
	h_constraintE_uncertaintyE_highFitProb->SetDirectory(m_pTFile);
	h_pull_jet_theta_jet_phi_lowFitProb = new TH2F("h_pull_jet_theta_jet_phi_lowFitProb", "fit probability < 0.1; pull #theta_{jet}; pull #phi_{jet}", 200, -10.0, 10.0, 200, -10.0, 10.0);
	h_pull_jet_theta_jet_phi_lowFitProb->SetDirectory(m_pTFile);
	h_pull_jet_theta_jet_phi_midFitProb = new TH2F("h_pull_jet_theta_jet_phi_midFitProb", "0.1 #leq fit probability < 0.9; pull #theta_{jet}; pull #phi_{jet}", 200, -10.0, 10.0, 200, -10.0, 10.0);
	h_pull_jet_theta_jet_phi_midFitProb->SetDirectory(m_pTFile);
	h_pull_jet_theta_jet_phi_highFitProb = new TH2F("h_pull_jet_theta_jet_phi_highFitProb", "0.9 #leq fit probability; pull #theta_{jet}; pull #phi_{jet}", 200, -10.0, 10.0, 200, -10.0, 10.0);
	h_pull_jet_theta_jet_phi_highFitProb->SetDirectory(m_pTFile);

	h_NormalizedResidualTheta_mcQmcJ_ThetaMcJet = new TH2F("h_NormalizedResidualTheta_mcQmcJ_ThetaMcJet", "; cos#theta_{j}^{MC} [radian]; (#theta_{j}^{MC} - #theta_{q}^{MC}) / #sigma_{#theta_{j}}", 100, -1.0, 1.0, 200, -10.0, 10.0);
	h_NormalizedResidualTheta_mcQmcJ_ThetaMcJet->SetDirectory(m_pTFile);
	h_NormalizedResidualPhi_mcQmcJ_PhiMcJet = new TH2F("h_NormalizedResidualPhi_mcQmcJ_PhiMcJet", "; #phi_{j}^{MC} [radian]; (#phi_{j}^{MC} - #phi_{q}^{MC}) / #sigma_{#phi_{j}}", 100, -PI, PI, 200, -10.0, 10.0);
	h_NormalizedResidualPhi_mcQmcJ_PhiMcJet->SetDirectory(m_pTFile);
	h_NormalizedResidualEnergy_mcQmcJ_EnergyMcJet = new TH2F("h_NormalizedResidualEnergy_mcQmcJ_EnergyMcJet", "; E_{j}^{MC} [GeV]; (E_{j}^{MC} - E_{q}^{MC}) / #sigma_{E_{j}}", 150, 0.0, 150.0, 200, -10.0, 10.0);
	h_NormalizedResidualEnergy_mcQmcJ_EnergyMcJet->SetDirectory(m_pTFile);
	h_NormalizedResidualTheta_mcQmcJ_EnergyMcJet = new TH2F("h_NormalizedResidualTheta_mcQmcJ_EnergyMcJet", "; E_{j}^{MC} [GeV]; (#theta_{j}^{MC} - #theta_{q}^{MC}) / #sigma_{#theta_{j}}", 150, 0.0, 150.0, 200, -10.0, 10.0);
	h_NormalizedResidualTheta_mcQmcJ_EnergyMcJet->SetDirectory(m_pTFile);
	h_NormalizedResidualPhi_mcQmcJ_EnergyMcJet = new TH2F("h_NormalizedResidualPhi_mcQmcJ_EnergyMcJet", "; E_{j}^{MC} [GeV]; (#phi_{j}^{MC} - #phi_{q}^{MC}) / #sigma_{#phi_{j}}", 150, 0.0, 150.0, 200, -10.0, 10.0);
	h_NormalizedResidualPhi_mcQmcJ_EnergyMcJet->SetDirectory(m_pTFile);

	h_NormalizedResidualTheta_mcJrecJ_ThetaRecJet = new TH2F("h_NormalizedResidualTheta_mcJrecJ_ThetaRecJet", "; cos#theta_{j}^{REC} [radian]; (#theta_{j}^{REC} - #theta_{j}^{MC}) / #sigma_{#theta_{j}}", 100, -1.0, 1.0, 200, -10.0, 10.0);
	h_NormalizedResidualTheta_mcJrecJ_ThetaRecJet->SetDirectory(m_pTFile);
	h_NormalizedResidualPhi_mcJrecJ_PhiRecJet = new TH2F("h_NormalizedResidualPhi_mcJrecJ_PhiRecJet", "; #phi_{j}^{REC} [radian]; (#phi_{j}^{REC} - #phi_{j}^{MC}) / #sigma_{#phi_{j}}", 100, -PI, PI, 200, -10.0, 10.0);
	h_NormalizedResidualPhi_mcJrecJ_PhiRecJet->SetDirectory(m_pTFile);
	h_NormalizedResidualEnergy_mcJrecJ_EnergyRecJet = new TH2F("h_NormalizedResidualEnergy_mcJrecJ_EnergyRecJet", "; E_{j}^{REC} [GeV]; (E_{j}^{REC} - E_{j}^{MC}) / #sigma_{E_{j}}", 150, 0.0, 150.0, 200, -10.0, 10.0);
	h_NormalizedResidualEnergy_mcJrecJ_EnergyRecJet->SetDirectory(m_pTFile);
	h_NormalizedResidualTheta_mcJrecJ_EnergyRecJet = new TH2F("h_NormalizedResidualTheta_mcJrecJ_EnergyRecJet", "; E_{j}^{REC} [GeV]; (#theta_{j}^{REC} - #theta_{j}^{MC}) / #sigma_{#theta_{j}}", 150, 0.0, 150.0, 200, -10.0, 10.0);
	h_NormalizedResidualTheta_mcJrecJ_EnergyRecJet->SetDirectory(m_pTFile);
	h_NormalizedResidualPhi_mcJrecJ_EnergyRecJet = new TH2F("h_NormalizedResidualPhi_mcJrecJ_EnergyRecJet", "; E_{j}^{REC} [GeV]; (#phi_{j}^{REC} - #phi_{j}^{MC}) / #sigma_{#phi_{j}}", 150, 0.0, 150.0, 200, -10.0, 10.0);
	h_NormalizedResidualPhi_mcJrecJ_EnergyRecJet->SetDirectory(m_pTFile);

	h_NormalizedResidualTheta_mcQrecJ_ThetaRecJet = new TH2F("h_NormalizedResidualTheta_mcQrecJ_ThetaRecJet", "; cos#theta_{j}^{REC} [radian]; (#theta_{j}^{REC} - #theta_{q}^{MC}) / #sigma_{#theta_{j}}", 100, -1.0, 1.0, 200, -10.0, 10.0);
	h_NormalizedResidualTheta_mcQrecJ_ThetaRecJet->SetDirectory(m_pTFile);
	h_NormalizedResidualPhi_mcQrecJ_PhiRecJet = new TH2F("h_NormalizedResidualPhi_mcQrecJ_PhiRecJet", "; #phi_{j}^{REC} [radian]; (#phi_{j}^{REC} - #phi_{q}^{MC}) / #sigma_{#phi_{j}}", 100, -PI, PI, 200, -10.0, 10.0);
	h_NormalizedResidualPhi_mcQrecJ_PhiRecJet->SetDirectory(m_pTFile);
	h_NormalizedResidualEnergy_mcQrecJ_EnergyRecJet = new TH2F("h_NormalizedResidualEnergy_mcQrecJ_EnergyRecJet", "; E_{j}^{REC} [GeV]; (E_{j}^{REC} - E_{q}^{MC}) / #sigma_{E_{j}}", 150, 0.0, 150.0, 200, -10.0, 10.0);
	h_NormalizedResidualEnergy_mcQrecJ_EnergyRecJet->SetDirectory(m_pTFile);
	h_NormalizedResidualTheta_mcQrecJ_EnergyRecJet = new TH2F("h_NormalizedResidualTheta_mcQrecJ_EnergyRecJet", "; E_{j}^{REC} [GeV]; (#theta_{j}^{REC} - #theta_{q}^{MC}) / #sigma_{#theta_{j}}", 150, 0.0, 150.0, 200, -10.0, 10.0);
	h_NormalizedResidualTheta_mcQrecJ_EnergyRecJet->SetDirectory(m_pTFile);
	h_NormalizedResidualPhi_mcQrecJ_EnergyRecJet = new TH2F("h_NormalizedResidualPhi_mcQrecJ_EnergyRecJet", "; E_{j}^{REC} [GeV]; (#phi_{j}^{REC} - #phi_{q}^{MC}) / #sigma_{#phi_{j}}", 150, 0.0, 150.0, 200, -10.0, 10.0);
	h_NormalizedResidualPhi_mcQrecJ_EnergyRecJet->SetDirectory(m_pTFile);

	h_recE_fitE = new TH1F("h_recE_fitE", "; ( E_{j}^{fit} - E_{j}^{rec} ) / E_{j}^{rec}; n_{events}", 400, -1.0, 1.0);
	h_recE_fitE->SetDirectory(m_pTFile);
	h_recPhifitPhi_recThetafitTheta = new TH2F("h_recPhifitPhi_recThetafitTheta", "; #phi_{j}^{fit} - #phi_{j}^{rec} [radian]; #theta_{j}^{fit} - #theta_{j}^{rec} [radian]", 200, -0.10, 0.10, 200, -0.10, 0.10);
	h_recPhifitPhi_recThetafitTheta->SetDirectory(m_pTFile);
	h_recPhifitPhi_recThetafitTheta_lowfit = new TH2F("h_recPhifitPhi_recThetafitTheta_lowfit", "fit probability < 0.1; #phi_{j}^{fit} - #phi_{j}^{rec} [radian]; #theta_{j}^{fit} - #theta_{j}^{rec} [radian]", 200, -0.10, 0.10, 200, -0.10, 0.10);
	h_recPhifitPhi_recThetafitTheta_lowfit->SetDirectory(m_pTFile);
	h_recPhifitPhi_recThetafitTheta_midfit = new TH2F("h_recPhifitPhi_recThetafitTheta_midfit", "0.1 #leq fit probability < 0.9; #phi_{j}^{fit} - #phi_{j}^{rec} [radian]; #theta_{j}^{fit} - #theta_{j}^{rec} [radian]", 200, -0.10, 0.10, 200, -0.10, 0.10);
	h_recPhifitPhi_recThetafitTheta_midfit->SetDirectory(m_pTFile);
	h_recPhifitPhi_recThetafitTheta_highfit = new TH2F("h_recPhifitPhi_recThetafitTheta_highfit", "0.9 #leq fit probability; #phi_{j}^{fit} - #phi_{j}^{rec} [radian]; #theta_{j}^{fit} - #theta_{j}^{rec} [radian]", 200, -0.10, 0.10, 200, -0.10, 0.10);
	h_recPhifitPhi_recThetafitTheta_highfit->SetDirectory(m_pTFile);
	h_recPhifitPhi_recThetafitTheta_SigmaThetaPhi = new TH2F("h_recPhifitPhi_recThetafitTheta_SigmaThetaPhi", "; ( #phi_{j}^{fit} - #phi_{j}^{rec} ) / #sigma(#phi_{j}^{ErrFl}); (#theta_{j}^{fit} - #theta_{j}^{rec} ) / #sigma(#theta{j}^{ErrFl})", 100, -5.0, 5.0, 100, -5.0, 5.0);
	h_recPhifitPhi_recThetafitTheta_SigmaThetaPhi->SetDirectory(m_pTFile);
	h_recPhifitPhi_recThetafitTheta_lowfit_SigmaThetaPhi = new TH2F("h_recPhifitPhi_recThetafitTheta_lowfit_SigmaThetaPhi", "fit probability < 0.1; ( #phi_{j}^{fit} - #phi_{j}^{rec} ) / #sigma(#phi_{j}^{ErrFl}); (#theta_{j}^{fit} - #theta_{j}^{rec} ) / #sigma(#theta{j}^{ErrFl})", 100, -5.0, 5.0, 100, -5.0, 5.0);
	h_recPhifitPhi_recThetafitTheta_lowfit_SigmaThetaPhi->SetDirectory(m_pTFile);
	h_recPhifitPhi_recThetafitTheta_midfit_SigmaThetaPhi = new TH2F("h_recPhifitPhi_recThetafitTheta_midfit_SigmaThetaPhi", "0.1 #leq fit probability < 0.9; ( #phi_{j}^{fit} - #phi_{j}^{rec} ) / #sigma(#phi_{j}^{ErrFl}); (#theta_{j}^{fit} - #theta_{j}^{rec} ) / #sigma(#theta{j}^{ErrFl})", 100, -5.0, 5.0, 100, -5.0, 5.0);
	h_recPhifitPhi_recThetafitTheta_midfit_SigmaThetaPhi->SetDirectory(m_pTFile);
	h_recPhifitPhi_recThetafitTheta_highfit_SigmaThetaPhi = new TH2F("h_recPhifitPhi_recThetafitTheta_highfit_SigmaThetaPhi", "0.9 #leq fit probability; ( #phi_{j}^{fit} - #phi_{j}^{rec} ) / #sigma(#phi_{j}^{ErrFl}); (#theta_{j}^{fit} - #theta_{j}^{rec} ) / #sigma(#theta{j}^{ErrFl})", 100, -5.0, 5.0, 100, -5.0, 5.0);
	h_recPhifitPhi_recThetafitTheta_highfit_SigmaThetaPhi->SetDirectory(m_pTFile);
}

void ZHllqq5CFit::Clear()
{
	m_nSLDecayBHadron = 0;
	m_nSLDecayCHadron = 0;
	m_nSLDecayTotal = 0;
	m_nJets = 0;
	m_nLeptons = 0;
	m_HDecayMode = -1;
	bestprob = 0.;
	m_bestNuCombination.clear();
	m_iError_wNu.clear();
	m_probability_wNu.clear();
	m_chi2_wNu.clear();
	m_n_itter_wNu.clear();
	m_startmassZ_wNu.clear();
	m_startmassH_wNu.clear();
	m_beststartmassZ_wNu.clear();
	m_beststartmassH_wNu.clear();
	m_Zmass_after_fit_wNu.clear();
	m_Hmass_after_fit_wNu.clear();
	m_bestphotonenergy_wNu.clear();
	m_chi2startmassZ_wNu.clear();
	m_chi2startmassH_wNu.clear();
	m_jet_startPx_wNu.clear();
	m_jet_startPy_wNu.clear();
	m_jet_startPz_wNu.clear();
	m_jet_startE_wNu.clear();
	m_jet_startTheta_wNu.clear();
	m_jet_startPhi_wNu.clear();
	m_jet_fittedE_wNu.clear();
	m_jet_fittedTheta_wNu.clear();
	m_jet_fittedPhi_wNu.clear();
	m_dijet_angle_wNu.clear();
	m_jet_startPx_wNu_bestfit.clear();
	m_jet_startPy_wNu_bestfit.clear();
	m_jet_startPz_wNu_bestfit.clear();
	m_jet_startE_wNu_bestfit.clear();
	m_jet_startTheta_wNu_bestfit.clear();
	m_jet_startPhi_wNu_bestfit.clear();
	m_jet_fittedE_wNu_bestfit.clear();
	m_jet_fittedTheta_wNu_bestfit.clear();
	m_jet_fittedPhi_wNu_bestfit.clear();
	m_dijet_angle_wNu_bestfit = 0.;
	m_jet_startPx_woNu.clear();
	m_jet_startPy_woNu.clear();
	m_jet_startPz_woNu.clear();
	m_jet_startE_woNu.clear();
	m_jet_startTheta_woNu.clear();
	m_jet_startPhi_woNu.clear();
	m_jet_fittedE_woNu.clear();
	m_jet_fittedTheta_woNu.clear();
	m_jet_fittedPhi_woNu.clear();
	m_dijet_angle_woNu = 0.;
	m_jet_startPx_best.clear();
	m_jet_startPy_best.clear();
	m_jet_startPz_best.clear();
	m_jet_startE_best.clear();
	m_jet_startTheta_best.clear();
	m_jet_startPhi_best.clear();
	m_jet_fittedE_best.clear();
	m_jet_fittedTheta_best.clear();
	m_jet_fittedPhi_best.clear();
	m_dijet_angle_best = 0.;
	m_lepton_startPx_wNu.clear();
	m_lepton_startPy_wNu.clear();
	m_lepton_startPz_wNu.clear();
	m_lepton_startE_wNu.clear();
	m_lepton_startPx_wNu_bestfit.clear();
	m_lepton_startPy_wNu_bestfit.clear();
	m_lepton_startPz_wNu_bestfit.clear();
	m_lepton_startE_wNu_bestfit.clear();
	m_lepton_startPx_woNu.clear();
	m_lepton_startPy_woNu.clear();
	m_lepton_startPz_woNu.clear();
	m_lepton_startE_woNu.clear();
	m_lepton_startPx_best.clear();
	m_lepton_startPy_best.clear();
	m_lepton_startPz_best.clear();
	m_lepton_startE_best.clear();
	m_jet_SigmaTheta.clear();
	m_jet_SigmaPhi.clear();
	m_jet_SigmaE.clear();
	m_lepton_SigmaTheta.clear();
	m_lepton_SigmaPhi.clear();
	m_lepton_SigmaInvpT.clear();
	m_pull_jet_E_wNu.clear();
	m_pull_jet_th_wNu.clear();
	m_pull_jet_phi_wNu.clear();
	m_pull_lepton_InvpT_wNu.clear();
	m_pull_lepton_th_wNu.clear();
	m_pull_lepton_phi_wNu.clear();
	m_iError_wNu_bestfit = -5;
	m_probability_wNu_bestfit = 0.;
	m_chi2_wNu_bestfit = 0.;
	m_n_itter_wNu_bestfit = 0;
	m_startmassZ_wNu_bestfit = 0.;
	m_startmassH_wNu_bestfit = 0.;
	m_beststartmassZ_wNu_bestfit = 0.;
	m_beststartmassH_wNu_bestfit = 0.;
	m_Zmass_after_fit_wNu_bestfit = 0.;
	m_Hmass_after_fit_wNu_bestfit = 0.;
	m_bestphotonenergy_wNu_bestfit = 0.;
	m_chi2startmassZ_wNu_bestfit = 0.;
	m_chi2startmassH_wNu_bestfit = 0.;
	m_pull_jet_E_wNu_bestfit.clear();
	m_pull_jet_th_wNu_bestfit.clear();
	m_pull_jet_phi_wNu_bestfit.clear();
	m_pull_lepton_InvpT_wNu_bestfit.clear();
	m_pull_lepton_th_wNu_bestfit.clear();
	m_pull_lepton_phi_wNu_bestfit.clear();
	m_iError_woNu = -5;
	m_probability_woNu = 0.;
	m_chi2best_woNu = 0.;
	m_n_itter_woNu = 0;
	m_startmassZ_woNu = 0.;
	m_startmassH_woNu = 0.;
	m_beststartmassZ_woNu = 0.;
	m_beststartmassH_woNu = 0.;
	m_Zmass_after_fit_woNu = 0.;
	m_Hmass_after_fit_woNu = 0.;
	m_bestphotonenergy_woNu = 0.;
	m_chi2startmassZ_woNu = 0.;
	m_chi2startmassH_woNu = 0.;
	m_pull_jet_E_woNu.clear();
	m_pull_jet_th_woNu.clear();
	m_pull_jet_phi_woNu.clear();
	m_pull_lepton_InvpT_woNu.clear();
	m_pull_lepton_th_woNu.clear();
	m_pull_lepton_phi_woNu.clear();
	m_iError_best = -5;
	m_probability_best = 0.;
	m_chi2_best = 0.;
	m_n_itter_best = 0;
	m_startmassZ_best = 0.;
	m_startmassH_best = 0.;
	m_beststartmassZ_best = 0.;
	m_beststartmassH_best = 0.;
	m_Zmass_after_fit_best = 0.;
	m_Hmass_after_fit_best = 0.;
	m_bestphotonenergy_best = 0.;
	m_chi2startmassZ_best = 0.;
	m_chi2startmassH_best = 0.;
	m_pull_jet_E_best.clear();
	m_pull_jet_th_best.clear();
	m_pull_jet_phi_best.clear();
	m_pull_lepton_InvpT_best.clear();
	m_pull_lepton_th_best.clear();
	m_pull_lepton_phi_best.clear();
	m_pxc_before_ISR_wNu = 0.;
	m_pyc_before_ISR_wNu = 0.;
	m_pzc_before_ISR_wNu = 0.;
	m_ec_before_ISR_wNu = 0.;
	m_zc_before_ISR_wNu = 0.;
	m_pxc_before_fit_wNu = 0.;
	m_pyc_before_fit_wNu = 0.;
	m_pzc_before_fit_wNu = 0.;
	m_ec_before_fit_wNu = 0.;
	m_zc_before_fit_wNu = 0.;
	m_pxc_after_fit_wNu = 0.;
	m_pyc_after_fit_wNu = 0.;
	m_pzc_after_fit_wNu = 0.;
	m_ec_after_fit_wNu = 0.;
	m_zc_after_fit_wNu = 0.;
	m_pxc_before_ISR_woNu = 0.;
	m_pyc_before_ISR_woNu = 0.;
	m_pzc_before_ISR_woNu = 0.;
	m_ec_before_ISR_woNu = 0.;
	m_zc_before_ISR_woNu = 0.;
	m_pxc_before_fit_woNu = 0.;
	m_pyc_before_fit_woNu = 0.;
	m_pzc_before_fit_woNu = 0.;
	m_ec_before_fit_woNu = 0.;
	m_zc_before_fit_woNu = 0.;
	m_pxc_after_fit_woNu = 0.;
	m_pyc_after_fit_woNu = 0.;
	m_pzc_after_fit_woNu = 0.;
	m_ec_after_fit_woNu = 0.;
	m_zc_after_fit_woNu = 0.;
	m_pxc_before_ISR_best = 0.;
	m_pyc_before_ISR_best = 0.;
	m_pzc_before_ISR_best = 0.;
	m_ec_before_ISR_best = 0.;
	m_zc_before_ISR_best = 0.;
	m_pxc_before_fit_best = 0.;
	m_pyc_before_fit_best = 0.;
	m_pzc_before_fit_best = 0.;
	m_ec_before_fit_best = 0.;
	m_zc_before_fit_best = 0.;
	m_pxc_after_fit_best = 0.;
	m_pyc_after_fit_best = 0.;
	m_pzc_after_fit_best = 0.;
	m_ec_after_fit_best = 0.;
	m_zc_after_fit_best = 0.;
	m_mcHmumuPx = 0.;
	m_mcHmumuPy = 0.;
	m_mcHmumuPz = 0.;
	m_mcHmumuPt = 0.;
	m_mcHmumuP = 0.;
	m_mcHmumuE = 0.;
	m_mcISR1Px = 0.;
	m_mcISR1Py = 0.;
	m_mcISR1Pz = 0.;
	m_mcISR1E = 0.;
	m_mcISR2Px = 0.;
	m_mcISR2Py = 0.;
	m_mcISR2Pz = 0.;
	m_mcISR2E = 0.;
	m_mcISRPx = 0.;
	m_mcISRPy = 0.;
	m_mcISRPz = 0.;
	m_mcISRE = 0.;
	m_ISR_startPx_wNu = NAN;
	m_ISR_startPy_wNu = NAN;
	m_ISR_startPz_wNu = NAN;
	m_ISR_startPx_woNu = NAN;
	m_ISR_startPy_woNu = NAN;
	m_ISR_startPz_woNu = NAN;
	m_TotalSigmaPx2 = 0.;
	m_TotalSigmaPy2 = 0.;
	m_TotalSigmaPz2 = 0.;
	m_TotalSigmaE2 = 0.;
	m_TotalSigmaPxPy = 0.;
	m_TotalSigmaPxPz = 0.;
	m_TotalSigmaPxE = 0.;
	m_TotalSigmaPyPz = 0.;
	m_TotalSigmaPyE = 0.;
	m_TotalSigmaPzE = 0.;
	m_ISRPx.clear();
	m_ISRPy.clear();
	m_ISRPz.clear();
	m_ISRPt.clear();
	m_ISRE.clear();
	m_ISRPx_total = 0.;
	m_ISRPy_total = 0.;
	m_ISRPz_total = 0.;
	m_ISRPt_total = 0.;
	m_ISRE_total = 0.;
	m_FSRPx.clear();
	m_FSRPy.clear();
	m_FSRPz.clear();
	m_FSRPt.clear();
	m_FSRE.clear();
	m_FSRPx_total = 0.;
	m_FSRPy_total = 0.;
	m_FSRPz_total = 0.;
	m_FSRPt_total = 0.;
	m_FSRE_total = 0.;
	m_SigmaPx2.clear();
	m_SigmaPxPy.clear();
	m_SigmaPxPz.clear();
	m_SigmaPxE.clear();
	m_SigmaPy2.clear();
	m_SigmaPyPz.clear();
	m_SigmaPyE.clear();
	m_SigmaPz2.clear();
	m_SigmaPzE.clear();
	m_SigmaE2.clear();
	m_Theta_mcQuark.clear();
	m_Phi_mcQuark.clear();
	m_Energy_mcQuark.clear();
	m_Theta_mcJet_wInvisibles.clear();
	m_Phi_mcJet_wInvisibles.clear();
	m_Energy_mcJet_wInvisibles.clear();
	m_Theta_mcJet_woInvisibles.clear();
	m_Phi_mcJet_woInvisibles.clear();
	m_Energy_mcJet_woInvisibles.clear();
	m_Theta_recoJet.clear();
	m_Phi_recoJet.clear();
	m_Energy_recoJet.clear();
	m_nPFOs_NeutralHadrons.clear();
	m_nPFOs_Photons.clear();
	m_nPFOs_ChargedPFOs.clear();
	m_nPFOs_Total.clear();
	m_PFOEnergy_NeutralHadrons.clear();
	m_PFOEnergy_Photons.clear();
	m_PFOEnergy_ChargedPFOs.clear();
	m_PFOEnergy_Total.clear();
	m_EnergyFraction_NeutralHadrons.clear();
	m_EnergyFraction_Photons.clear();
	m_EnergyFraction_ChargedPFOs.clear();
	m_ResidualTheta_mcQuark_mcJet.clear();
	m_ResidualPhi_mcQuark_mcJet.clear();
	m_ResidualEnergy_mcQuark_mcJet.clear();
	m_ResidualTheta_mcJet_recoJet.clear();
	m_ResidualPhi_mcJet_recoJet.clear();
	m_ResidualEnergy_mcJet_recoJet.clear();
	m_ResidualTheta_mcQuark_recoJet.clear();
	m_ResidualPhi_mcQuark_recoJet.clear();
	m_ResidualEnergy_mcQuark_recoJet.clear();
	m_NormalizedResidualTheta_mcQuark_mcJet.clear();
	m_NormalizedResidualPhi_mcQuark_mcJet.clear();
	m_NormalizedResidualEnergy_mcQuark_mcJet.clear();
	m_NormalizedResidualTheta_mcJet_recoJet.clear();
	m_NormalizedResidualPhi_mcJet_recoJet.clear();
	m_NormalizedResidualEnergy_mcJet_recoJet.clear();
	m_NormalizedResidualTheta_mcQuark_recoJet.clear();
	m_NormalizedResidualPhi_mcQuark_recoJet.clear();
	m_NormalizedResidualEnergy_mcQuark_recoJet.clear();
	m_JetErrorEnergy_ErrorFlow.clear();
	m_JetErrorTheta_ErrorFlow.clear();
	m_JetErrorPhi_ErrorFlow.clear();
}

void ZHllqq5CFit::processRunHeader()
{
	_nRun++ ;
}

void ZHllqq5CFit::compcorrect() //finds the jet with cross checking
{
	if( std::abs(delta_theta[bestjet_phi]) < std::abs(delta_phi[bestjet_th]) )
	{
		bestjet=bestjet_phi;
	}
	if( std::abs(delta_theta[bestjet_phi]) > std::abs(delta_phi[bestjet_th]) )
	{
		bestjet=bestjet_th;
	}
}

void ZHllqq5CFit::Setvalues()
{

	Px=0.;
	Px2=0.;
	Py=0.;
	Py2=0.;
	Pz=0.;
	Pz2=0.;
	Pt=0.;
	Pt2=0.;
	P=0.;
	P2=0.;
	SigPx2=0.;
	SigPxSigPy=0.;
	SigPxSigPz=0.;
	SigPy2=0.;
	SigPySigPz=0.;
	SigPz2=0.;
	SigE2=0.;
	dth_dpx=0.;
	dth_dpy=0.;
	dth_dpz=0.;
	dphi_dpx=0.;
	dphi_dpy=0.;
	JetResE=0.;
	JetResTheta=0.;
	JetResPhi=0.;
	bestprob = 0.;
	bestnit = 0;
	beststartmassZ = 0., beststartmassH = 0.;
	startmassZ = 0., startmassH = 0.;
	bestphotonenergy = 0.;
 	besterr = 999;
 	bestzvalue = 10000.;
 	chi2startmassZ = 0.;
 	chi2startmassH = 0.;
	memset(Zmomentum, 0, sizeof(Zmomentum));
	memset(Hmomentum, 0, sizeof(Hmomentum));
	memset(ISRmomentum, 0, sizeof(ISRmomentum));
	Z_Energy=0.;
	H_Energy=0.;
	chi2best=0.;
	errorcode=0.;
	streamlog_out(DEBUG)  << "Values set to defaults" <<std::endl;

}

void ZHllqq5CFit::getHmumu4Momentum( EVENT::LCEvent *pLCEvent )
{
	LCCollection *MCPCol{};
	m_nEvt = pLCEvent->getEventNumber();
	try
	{
		MCPCol			= pLCEvent->getCollection( MCPCollection );
		MCParticle* ISR1 = dynamic_cast<MCParticle*>(MCPCol->getElementAt(6));
		MCParticle* ISR2 = dynamic_cast<MCParticle*>(MCPCol->getElementAt(7));
		MCParticle* Muon1 = dynamic_cast<MCParticle*>(MCPCol->getElementAt(8));
		MCParticle* Muon2 = dynamic_cast<MCParticle*>(MCPCol->getElementAt(9));
		MCParticle* Higgs = dynamic_cast<MCParticle*>(MCPCol->getElementAt(10));
		if ( ISR1->getPDG() == 22 )
		{
			m_mcISR1Px = ISR1->getMomentum()[0];
			m_mcISR1Py = ISR1->getMomentum()[1];
			m_mcISR1Pz = ISR1->getMomentum()[2];
			m_mcISR1E = ISR1->getEnergy();
		}
		if ( ISR2->getPDG() == 22 )
		{
			m_mcISR2Px = ISR2->getMomentum()[0];
			m_mcISR2Py = ISR2->getMomentum()[1];
			m_mcISR2Pz = ISR2->getMomentum()[2];
			m_mcISR2E = ISR2->getEnergy();
		}
		m_mcISRPx = ISR1->getMomentum()[0] + ISR2->getMomentum()[0];
		m_mcISRPy = ISR1->getMomentum()[1] + ISR2->getMomentum()[1];
		m_mcISRPz = ISR1->getMomentum()[2] + ISR2->getMomentum()[2];
		m_mcISRE = ISR1->getEnergy() + ISR2->getEnergy();
		h_mcpISRPx->Fill( m_mcISR1Px , m_mcISR2Px );
		h_mcpISRPy->Fill( m_mcISR1Py , m_mcISR2Py );
		h_mcpISRPz->Fill( m_mcISR1Pz , m_mcISR2Pz );
		h_mcpISRE->Fill( m_mcISR1E , m_mcISR2E );
		if ( abs(Muon1->getPDG()) == 13 && abs(Muon2->getPDG()) == 13 && abs(Higgs->getPDG()) == 25 )
		{
			m_mcHmumuPx = Muon1->getMomentum()[0] + Muon2->getMomentum()[0] + Higgs->getMomentum()[0];
			m_mcHmumuPy = Muon1->getMomentum()[1] + Muon2->getMomentum()[1] + Higgs->getMomentum()[1];
			m_mcHmumuPz = Muon1->getMomentum()[2] + Muon2->getMomentum()[2] + Higgs->getMomentum()[2];
			m_mcHmumuPt = std::sqrt( pow( m_mcHmumuPx , 2 ) + pow( m_mcHmumuPy , 2 ) );
			m_mcHmumuP = std::sqrt( pow( m_mcHmumuPx , 2 ) + pow( m_mcHmumuPy , 2 )  + pow( m_mcHmumuPz , 2 ) );
			m_mcHmumuE = Muon1->getEnergy() + Muon2->getEnergy() + Higgs->getEnergy();
		}
	}
	catch(DataNotAvailableException &e)
	{
		streamlog_out(MESSAGE) << "getHmumu4Momentum : Input collections not found in event " << m_nEvt << std::endl;
	}

}

void ZHllqq5CFit::getMCISRFSR( EVENT::LCEvent *pLCEvent )
{
	LCCollection *MCISRCol{};
	LCCollection *MCFSRCol{};
	try
	{
		MCISRCol		= pLCEvent->getCollection( mcISRcol );
		MCFSRCol		= pLCEvent->getCollection( mcFSRcol );
		int n_mcISR		= MCISRCol->getNumberOfElements();
		int n_mcFSR		= MCFSRCol->getNumberOfElements();
//		float ISRPx_total,ISRPy_total,ISRPz_total,ISRE_total;
//		float FSRPx_total,FSRPy_total,FSRPz_total,FSRE_total;
		for (int i_ISR = 0 ; i_ISR < n_mcISR ; ++i_ISR)
		{
			ReconstructedParticle *mcISR = dynamic_cast<ReconstructedParticle*>( MCISRCol->getElementAt( i_ISR ) ) ;
			m_ISRPx.push_back( mcISR->getMomentum()[0] );
			m_ISRPx_total += mcISR->getMomentum()[0];
			m_ISRPy.push_back( mcISR->getMomentum()[1] );
			m_ISRPy_total += mcISR->getMomentum()[1];
			m_ISRPz.push_back( mcISR->getMomentum()[2] );
			m_ISRPz_total += mcISR->getMomentum()[2];
			m_ISRPt.push_back( std::sqrt( pow( mcISR->getMomentum()[0] , 2 ) + pow( mcISR->getMomentum()[1] , 2 ) ) );
			m_ISRE.push_back( mcISR->getEnergy() );
			m_ISRE_total += mcISR->getEnergy();
		}
//		m_ISRPx_total.push_back( ISRPx_total );
//		m_ISRPy_total.push_back( ISRPy_total );
//		m_ISRPz_total.push_back( ISRPz_total );
		m_ISRPt_total = std::sqrt( pow(m_ISRPx_total , 2 ) + pow( m_ISRPy_total , 2 ) );
//		m_ISRE_total.push_back( ISRE_total );

		for (int i_FSR = 0 ; i_FSR < n_mcFSR ; ++i_FSR)
		{
			ReconstructedParticle *mcFSR = dynamic_cast<ReconstructedParticle*>( MCFSRCol->getElementAt( i_FSR ) ) ;
			m_FSRPx.push_back( mcFSR->getMomentum()[0] );
			m_FSRPx_total += mcFSR->getMomentum()[0];
			m_FSRPy.push_back( mcFSR->getMomentum()[1] );
			m_FSRPy_total += mcFSR->getMomentum()[1];
			m_FSRPz.push_back( mcFSR->getMomentum()[2] );
			m_FSRPz_total += mcFSR->getMomentum()[2];
			m_FSRPt.push_back( std::sqrt( pow( mcFSR->getMomentum()[0] , 2 ) + pow( mcFSR->getMomentum()[1] , 2 ) ) );
			m_FSRE.push_back( mcFSR->getEnergy() );
			m_FSRE_total += mcFSR->getEnergy();
		}
//		m_FSRPx_total.push_back( FSRPx_total );
//		m_FSRPy_total.push_back( FSRPy_total );
//		m_FSRPz_total.push_back( FSRPz_total );
		m_FSRPt_total = std::sqrt( pow(m_FSRPx_total , 2 ) + pow( m_FSRPy_total , 2 ) );
//		m_FSRE_total.push_back( FSRE_total );
		m_pTTree_5->Fill();
	}
	catch(DataNotAvailableException &e)
	{
		streamlog_out(MESSAGE) << "getMCISRFSR : Input collections not found in event " << m_nEvt << std::endl;
	}
}

void ZHllqq5CFit::getJetAngleResolution( EVENT::LCEvent *pLCEvent , int nSLDecayTotal )
{
	LCCollection *RecoJetCol{};
	LCCollection *MCPCol{};
	LCCollection *TrueJet_wNuCol{};
	LCCollection *TrueJet_woNuCol{};
	m_nEvt = pLCEvent->getEventNumber();
	try
	{
		RecoJetCol		= pLCEvent->getCollection( jetcollection );
		MCPCol			= pLCEvent->getCollection( MCPCollection );
		TrueJet_wNuCol		= pLCEvent->getCollection( mcjetcollection_wNu );
		TrueJet_woNuCol	= pLCEvent->getCollection( mcjetcollection_woNu );
		int n_MCPs		= MCPCol->getNumberOfElements();
		int n_RecoJets		= RecoJetCol->getNumberOfElements();
		int n_TrueJetswNu	= TrueJet_wNuCol->getNumberOfElements();
		int n_TrueJetswoNu	= TrueJet_woNuCol->getNumberOfElements();
		if ( !( n_RecoJets == 2 && n_TrueJetswoNu == 2 && n_TrueJetswNu == 2 ) )
		{
			streamlog_out(MESSAGE) << "Number of n_TrueJets mis-match n_RecoJets in event " << m_nEvt << std::endl;
			return;
		}
		int i_mcHiggs = 0;
		for (int i_mcp = 0 ; i_mcp < n_MCPs ; ++i_mcp )
		{
			MCParticle* mcp = dynamic_cast<MCParticle*>(MCPCol->getElementAt(i_mcp));
			if ( mcp->getPDG() == 25 ) i_mcHiggs = i_mcp;
		}
		streamlog_out(MESSAGE) << "Found Higgs boson at index " << i_mcHiggs << " in MCParticle Collection" << std::endl;
		MCParticle *mcHiggs = dynamic_cast<MCParticle*>(MCPCol->getElementAt(i_mcHiggs));
		std::vector<const EVENT::MCParticle*> mcQuark;
		std::vector<const EVENT::ReconstructedParticle *> mcJet_wNu;
		std::vector<const EVENT::ReconstructedParticle *> mcJet_woNu;
		std::vector<const EVENT::ReconstructedParticle *> recoJet;

		MCParticle *mcQ1 = mcHiggs->getDaughters()[0];
		MCParticle *mcQ2 = mcHiggs->getDaughters()[1];
		TVector3 QuarkMomentum1( mcQ1->getMomentum() ); QuarkMomentum1.SetMag(1.0);
		TVector3 QuarkMomentum2( mcQ2->getMomentum() ); QuarkMomentum2.SetMag(1.0);
		mcQuark.push_back( mcQ1 );
		mcQuark.push_back( mcQ2 );

		ReconstructedParticle *trueJet_wNu1 = dynamic_cast<ReconstructedParticle*>( TrueJet_wNuCol->getElementAt( 0 ) ) ;
		ReconstructedParticle *trueJet_wNu2 = dynamic_cast<ReconstructedParticle*>( TrueJet_wNuCol->getElementAt( 1 ) ) ;
		TVector3 trueJet_wNuMomentum1( trueJet_wNu1->getMomentum() ); trueJet_wNuMomentum1.SetMag(1.0);
		TVector3 trueJet_wNuMomentum2( trueJet_wNu2->getMomentum() ); trueJet_wNuMomentum2.SetMag(1.0);
		if ( trueJet_wNuMomentum1.Dot( QuarkMomentum1 ) >= trueJet_wNuMomentum2.Dot( QuarkMomentum1 ) && trueJet_wNuMomentum2.Dot( QuarkMomentum2 ) >= trueJet_wNuMomentum1.Dot( QuarkMomentum2 ) )
		{
			mcJet_wNu.push_back( trueJet_wNu1 );
			mcJet_wNu.push_back( trueJet_wNu2 );
			streamlog_out(MESSAGE) << "MCJet_wNu1 corresponds to MCQuark1 and MCJet_wNu2 corresponds to MCQuark2" << std::endl;
		}
		else
		{
			mcJet_wNu.push_back( trueJet_wNu2 );
			mcJet_wNu.push_back( trueJet_wNu1 );
			streamlog_out(MESSAGE) << "MCJet_wNu1 corresponds to MCQuark2 and MCJet_wNu2 corresponds to MCQuark1" << std::endl;
		}

		ReconstructedParticle *trueJet_woNu1 = dynamic_cast<ReconstructedParticle*>( TrueJet_woNuCol->getElementAt( 0 ) ) ;
		ReconstructedParticle *trueJet_woNu2 = dynamic_cast<ReconstructedParticle*>( TrueJet_woNuCol->getElementAt( 1 ) ) ;
		TVector3 trueJet_woNuMomentum1( trueJet_woNu1->getMomentum() ); trueJet_woNuMomentum1.SetMag(1.0);
		TVector3 trueJet_woNuMomentum2( trueJet_woNu2->getMomentum() ); trueJet_woNuMomentum2.SetMag(1.0);
		if ( trueJet_woNuMomentum1.Dot( QuarkMomentum1 ) >= trueJet_woNuMomentum2.Dot( QuarkMomentum1 ) && trueJet_woNuMomentum2.Dot( QuarkMomentum2 ) >= trueJet_woNuMomentum1.Dot( QuarkMomentum2 ) )
		{
			mcJet_woNu.push_back( trueJet_woNu1 );
			mcJet_woNu.push_back( trueJet_woNu2 );
			streamlog_out(MESSAGE) << "MCJet_woNu1 corresponds to MCQuark1 and MCJet_woNu2 corresponds to MCQuark2" << std::endl;
		}
		else
		{
			mcJet_woNu.push_back( trueJet_woNu2 );
			mcJet_woNu.push_back( trueJet_woNu1 );
			streamlog_out(MESSAGE) << "MCJet_woNu1 corresponds to MCQuark2 and MCJet_woNu2 corresponds to MCQuark1" << std::endl;
		}

		ReconstructedParticle *RecoJet1 = dynamic_cast<ReconstructedParticle*>( RecoJetCol->getElementAt( 0 ) ) ;
		ReconstructedParticle *RecoJet2 = dynamic_cast<ReconstructedParticle*>( RecoJetCol->getElementAt( 1 ) ) ;
		TVector3 RecoJetMomentum1( RecoJet1->getMomentum() ); RecoJetMomentum1.SetMag(1.0);
		TVector3 RecoJetMomentum2( RecoJet2->getMomentum() ); RecoJetMomentum2.SetMag(1.0);
		if ( RecoJetMomentum1.Dot( QuarkMomentum1 ) >= RecoJetMomentum2.Dot( QuarkMomentum1 ) && RecoJetMomentum2.Dot( QuarkMomentum2 ) >= RecoJetMomentum1.Dot( QuarkMomentum2 ) )
		{
			recoJet.push_back( RecoJet1 );
			recoJet.push_back( RecoJet2 );
			streamlog_out(MESSAGE) << "RecoJet1 corresponds to MCQuark1 and RecoJet2 corresponds to MCQuark2" << std::endl;
		}
		else
		{
			recoJet.push_back( RecoJet2 );
			recoJet.push_back( RecoJet1 );
			streamlog_out(MESSAGE) << "RecoJet1 corresponds to MCQuark2 and RecoJet2 corresponds to MCQuark1" << std::endl;
		}

		float px, px2, py, py2, pz, pz2, pt, pt2, p, p2;
		float Sigpx2, SigpxSigpy, Sigpy2, SigpxSigpz, SigpySigpz, Sigpz2, Sige2;
		float Dth_Dpx, Dth_Dpy, Dth_Dpz;
		float Dphi_Dpx, Dphi_Dpy;
		float JetEnergyErr, JetThetaErr, JetPhiErr;
		for ( int i_jet = 0 ; i_jet < 2 ; ++i_jet )
		{
			TVector3 mcQuarkP( ( mcQuark[ i_jet ] )->getMomentum() );
			TVector3 mcQuarkPt( mcQuarkP[0] , mcQuarkP[1] , 0.0 );
			mcQuarkP.SetMag(1.0); mcQuarkPt.SetMag(1.0);
			TVector3 mcJetwNuP( ( mcJet_wNu[ i_jet ] )->getMomentum() );
			TVector3 mcJetwNuPt( mcJetwNuP[0] , mcJetwNuP[1] , 0.0 );
			mcJetwNuP.SetMag(1.0); mcJetwNuPt.SetMag(1.0);
			TVector3 mcJetwoNuP( ( mcJet_woNu[ i_jet ] )->getMomentum() );
			TVector3 mcJetwoNuPt( mcJetwoNuP[0] , mcJetwoNuP[1] , 0.0 );
			mcJetwoNuP.SetMag(1.0); mcJetwoNuPt.SetMag(1.0);
			TVector3 recoJetP( ( recoJet[ i_jet ] )->getMomentum() );
			TVector3 recoJetPt( recoJetP[0] , recoJetP[1] , 0.0 );
			recoJetP.SetMag(1.0); recoJetPt.SetMag(1.0);

			float mcQuark_Theta = mcQuarkP.Theta();
			float mcQuark_Phi = mcQuarkP.Phi();
			m_Theta_mcQuark.push_back( mcQuark_Theta );
			m_Phi_mcQuark.push_back( mcQuark_Phi );
			m_Energy_mcQuark.push_back( ( mcQuark[ i_jet ] )->getEnergy() );
			float mcJetwNu_Theta = mcJetwNuP.Theta();
			float mcJetwNu_Phi = mcJetwNuP.Phi();
			m_Theta_mcJet_wInvisibles.push_back( mcJetwNu_Theta );
			m_Phi_mcJet_wInvisibles.push_back( mcJetwNu_Phi );
			m_Energy_mcJet_wInvisibles.push_back( ( mcJet_wNu[ i_jet ] )->getEnergy() );
			float mcJetwoNu_Theta = mcJetwoNuP.Theta();
			float mcJetwoNu_Phi = mcJetwoNuP.Phi();
			m_Theta_mcJet_woInvisibles.push_back( mcJetwoNu_Theta );
			m_Phi_mcJet_woInvisibles.push_back( mcJetwoNu_Phi );
			m_Energy_mcJet_woInvisibles.push_back( ( mcJet_woNu[ i_jet ] )->getEnergy() );
			float recoJet_Theta = recoJetP.Theta();
			float recoJet_Phi = recoJetP.Phi();
			m_Theta_recoJet.push_back( recoJet_Theta );
			m_Phi_recoJet.push_back( recoJet_Phi );
			m_Energy_recoJet.push_back( ( recoJet[ i_jet ] )->getEnergy() );

			px =	( recoJet[ i_jet ] )->getMomentum()[0];
			px2 =	std::pow( px , 2 );
			py =	( recoJet[ i_jet ] )->getMomentum()[1];
			py2 =	std::pow( py , 2 );
			pz =	( recoJet[ i_jet ] )->getMomentum()[2];
			pz2 =	std::pow( pz , 2 );
			pt2 =	px2 + py2;
			pt =	std::sqrt( pt2 );
			p =	std::sqrt( px2 + py2 + pz2 );
			p2 =	std::pow( p , 2 );

			Sigpx2 =		( recoJet[ i_jet ] )->getCovMatrix()[0];
			SigpxSigpy =	( recoJet[ i_jet ] )->getCovMatrix()[1];
			Sigpy2 =		( recoJet[ i_jet ] )->getCovMatrix()[2];
			SigpxSigpz =	( recoJet[ i_jet ] )->getCovMatrix()[3];
			SigpySigpz =	( recoJet[ i_jet ] )->getCovMatrix()[4];
			Sigpz2 =		( recoJet[ i_jet ] )->getCovMatrix()[5];
			Sige2 =		( recoJet[ i_jet ] )->getCovMatrix()[9];

			Dth_Dpx =	px * pz / ( p2 * pt );
			Dth_Dpy =	py * pz / ( p2 * pt );
			Dth_Dpz =	-pt / p2;

			Dphi_Dpx =	-py / pt2;
			Dphi_Dpy =	px / pt2;

			JetEnergyErr =	sigmaScaleFactor * std::sqrt( Sige2 );
			JetThetaErr =	m_JetResThetaCoeff * std::sqrt( std::fabs( Sigpx2 * std::pow( Dth_Dpx , 2 ) + Sigpy2 * std::pow( Dth_Dpy , 2 ) + Sigpz2 * std::pow( Dth_Dpz , 2 ) + 2 * ( SigpxSigpy * Dth_Dpx * Dth_Dpy ) + 2 * ( SigpySigpz * Dth_Dpy * Dth_Dpz ) + 2 * ( SigpxSigpz * Dth_Dpx * Dth_Dpz ) ) );
			JetPhiErr =	m_JetResPhiCoeff * std::sqrt( std::fabs( Sigpx2 * std::pow( Dphi_Dpx , 2 ) + Sigpy2 * std::pow( Dphi_Dpy , 2 ) + 2 * ( SigpxSigpy * Dphi_Dpx * Dphi_Dpy ) ) );

			m_ResidualTheta_mcQuark_mcJet.push_back( mcJetwNu_Theta - mcQuark_Theta );
			m_NormalizedResidualTheta_mcQuark_mcJet.push_back( ( mcJetwNu_Theta - mcQuark_Theta) / JetThetaErr );
//			m_ResidualPhi_mcQuark_mcJet.push_back( acos( mcQuarkPt.Dot(mcJetwNuPt) ) );
			streamlog_out(MESSAGE) << "mcJetwNu_Theta - mcQuark_Theta = " << mcJetwNu_Theta << " - " << mcQuark_Theta << std::endl;
			streamlog_out(MESSAGE) << "mcJetwNu_Phi - mcQuark_Phi = " << mcJetwNu_Phi << " - " << mcQuark_Phi << std::endl;
			streamlog_out(MESSAGE) << "mcJetwNu_Energy - mcQuark_Energy = " << ( mcJet_wNu[ i_jet ] )->getEnergy() << " - " << ( mcQuark[ i_jet ] )->getEnergy() << std::endl;
			if ( mcJetwNu_Phi > mcQuark_Phi )
			{
				m_ResidualPhi_mcQuark_mcJet.push_back( acos( mcQuarkPt.Dot(mcJetwNuPt) ) );
				m_NormalizedResidualPhi_mcQuark_mcJet.push_back( acos( mcQuarkPt.Dot(mcJetwNuPt) ) / JetPhiErr );
			}
			else
			{
				m_ResidualPhi_mcQuark_mcJet.push_back( -acos( mcQuarkPt.Dot(mcJetwNuPt) ) );
				m_NormalizedResidualPhi_mcQuark_mcJet.push_back( -acos( mcQuarkPt.Dot(mcJetwNuPt) ) / JetPhiErr );
			}
			m_ResidualEnergy_mcQuark_mcJet.push_back( ( mcJet_wNu[ i_jet ] )->getEnergy() - ( mcQuark[ i_jet ] )->getEnergy() );
			m_NormalizedResidualEnergy_mcQuark_mcJet.push_back( ( ( mcJet_wNu[ i_jet ] )->getEnergy() - ( mcQuark[ i_jet ] )->getEnergy() ) / JetEnergyErr );

			m_ResidualTheta_mcJet_recoJet.push_back( recoJet_Theta - mcJetwoNu_Theta );
			m_NormalizedResidualTheta_mcJet_recoJet.push_back( ( recoJet_Theta - mcJetwoNu_Theta ) / JetThetaErr );
//			m_ResidualPhi_mcJet_recoJet.push_back( acos( mcJetwoNuPt.Dot(recoJetPt) ) );
			streamlog_out(MESSAGE) << "recoJet_Theta - mcJetwoNu_Theta = " << recoJet_Theta << " - " << mcJetwoNu_Theta << std::endl;
			streamlog_out(MESSAGE) << "recoJet_Phi - mcJetwoNu_Phi = " << recoJet_Phi << " - " << mcJetwoNu_Phi << std::endl;
			streamlog_out(MESSAGE) << "recoJet_Energy - mcJetwoNu_Energy = " << ( recoJet[ i_jet ] )->getEnergy() << " - " << ( mcJet_woNu[ i_jet ] )->getEnergy() << std::endl;
			if ( recoJet_Phi > mcJetwoNu_Phi )
			{
				m_ResidualPhi_mcJet_recoJet.push_back( acos( mcJetwoNuPt.Dot(recoJetPt) ) );
				m_NormalizedResidualPhi_mcJet_recoJet.push_back( acos( mcJetwoNuPt.Dot(recoJetPt) ) / JetPhiErr );
			}
			else
			{
				m_ResidualPhi_mcJet_recoJet.push_back( -acos( mcJetwoNuPt.Dot(recoJetPt) ) );
				m_NormalizedResidualPhi_mcJet_recoJet.push_back( -acos( mcJetwoNuPt.Dot(recoJetPt) ) / JetPhiErr );
			}
			m_ResidualEnergy_mcJet_recoJet.push_back( ( recoJet[ i_jet ] )->getEnergy() - ( mcJet_woNu[ i_jet ] )->getEnergy() );
			m_NormalizedResidualEnergy_mcJet_recoJet.push_back( ( ( recoJet[ i_jet ] )->getEnergy() - ( mcJet_woNu[ i_jet ] )->getEnergy() ) / JetEnergyErr );

			m_ResidualTheta_mcQuark_recoJet.push_back( recoJet_Theta - mcQuark_Theta );
			m_NormalizedResidualTheta_mcQuark_recoJet.push_back( ( recoJet_Theta - mcQuark_Theta ) / JetThetaErr );
//			m_ResidualPhi_mcQuark_recoJet.push_back( acos( mcQuarkPt.Dot(recoJetPt) ) );
			streamlog_out(MESSAGE) << "recoJet_Theta - mcQuark_Theta = " << recoJet_Theta << " - " << mcQuark_Theta << std::endl;
			streamlog_out(MESSAGE) << "recoJet_Phi - mcQuark_Phi = " << recoJet_Phi << " - " << mcQuark_Phi << std::endl;
			streamlog_out(MESSAGE) << "recoJet_Energy - mcQuark_Energy = " << ( recoJet[ i_jet ] )->getEnergy() << " - " << ( mcQuark[ i_jet ] )->getEnergy() << std::endl;
			if ( recoJet_Phi > mcQuark_Phi )
			{
				m_ResidualPhi_mcQuark_recoJet.push_back( acos( mcQuarkPt.Dot(recoJetPt) ) );
				m_NormalizedResidualPhi_mcQuark_recoJet.push_back( acos( mcQuarkPt.Dot(recoJetPt) ) / JetPhiErr );
			}
			else
			{
				m_ResidualPhi_mcQuark_recoJet.push_back( -acos( mcQuarkPt.Dot(recoJetPt) ) );
				m_NormalizedResidualPhi_mcQuark_recoJet.push_back( -acos( mcQuarkPt.Dot(recoJetPt) ) / JetPhiErr );
			}
			m_ResidualEnergy_mcQuark_recoJet.push_back( ( recoJet[ i_jet ] )->getEnergy() - ( mcQuark[ i_jet ] )->getEnergy() );
			m_NormalizedResidualEnergy_mcQuark_recoJet.push_back( ( ( recoJet[ i_jet ] )->getEnergy() - ( mcQuark[ i_jet ] )->getEnergy() ) / JetEnergyErr );
			m_JetErrorEnergy_ErrorFlow.push_back( JetEnergyErr );
			m_JetErrorTheta_ErrorFlow.push_back( JetThetaErr );
			m_JetErrorPhi_ErrorFlow.push_back( JetPhiErr );

			if ( nSLDecayTotal == 0 )
			{
				h_NormalizedResidualTheta_mcQmcJ_ThetaMcJet->Fill( cos( mcJetwNu_Theta ) , m_ResidualTheta_mcQuark_mcJet[ i_jet ] / JetThetaErr );
				h_NormalizedResidualPhi_mcQmcJ_PhiMcJet->Fill( mcJetwNu_Phi , m_ResidualPhi_mcQuark_mcJet[ i_jet ] / JetPhiErr );
				h_NormalizedResidualEnergy_mcQmcJ_EnergyMcJet->Fill( ( mcJet_wNu[ i_jet ] )->getEnergy() , m_ResidualEnergy_mcQuark_mcJet[ i_jet ] / JetEnergyErr );
				h_NormalizedResidualTheta_mcQmcJ_EnergyMcJet->Fill( ( mcJet_wNu[ i_jet ] )->getEnergy() , m_ResidualTheta_mcQuark_mcJet[ i_jet ] / JetThetaErr );
				h_NormalizedResidualPhi_mcQmcJ_EnergyMcJet->Fill( ( mcJet_wNu[ i_jet ] )->getEnergy() , m_ResidualPhi_mcQuark_mcJet[ i_jet ] / JetPhiErr );

				h_NormalizedResidualTheta_mcJrecJ_ThetaRecJet->Fill( cos( recoJet_Theta ) , m_ResidualTheta_mcJet_recoJet[ i_jet ] / JetThetaErr );
				h_NormalizedResidualPhi_mcJrecJ_PhiRecJet->Fill( mcJetwNu_Phi , m_ResidualPhi_mcJet_recoJet[ i_jet ] / JetPhiErr );
				h_NormalizedResidualEnergy_mcJrecJ_EnergyRecJet->Fill( ( recoJet[ i_jet ] )->getEnergy() , m_ResidualEnergy_mcJet_recoJet[ i_jet ] / JetEnergyErr );
				h_NormalizedResidualTheta_mcJrecJ_EnergyRecJet->Fill( ( recoJet[ i_jet ] )->getEnergy() , m_ResidualTheta_mcJet_recoJet[ i_jet ] / JetThetaErr );
				h_NormalizedResidualPhi_mcJrecJ_EnergyRecJet->Fill( ( recoJet[ i_jet ] )->getEnergy() , m_ResidualPhi_mcJet_recoJet[ i_jet ] / JetPhiErr );

				h_NormalizedResidualTheta_mcQrecJ_ThetaRecJet->Fill( cos( recoJet_Theta ) , m_ResidualTheta_mcQuark_recoJet[ i_jet ] / JetThetaErr );
				h_NormalizedResidualPhi_mcQrecJ_PhiRecJet->Fill( mcJetwNu_Phi , m_ResidualPhi_mcQuark_recoJet[ i_jet ] / JetPhiErr );
				h_NormalizedResidualEnergy_mcQrecJ_EnergyRecJet->Fill( ( recoJet[ i_jet ] )->getEnergy() , m_ResidualEnergy_mcQuark_recoJet[ i_jet ] / JetEnergyErr );
				h_NormalizedResidualTheta_mcQrecJ_EnergyRecJet->Fill( ( recoJet[ i_jet ] )->getEnergy() , m_ResidualTheta_mcQuark_recoJet[ i_jet ] / JetThetaErr );
				h_NormalizedResidualPhi_mcQrecJ_EnergyRecJet->Fill( ( recoJet[ i_jet ] )->getEnergy() , m_ResidualPhi_mcQuark_recoJet[ i_jet ] / JetPhiErr );
			}
			ReconstructedParticleVec jetPFOs  = recoJet[ i_jet ]->getParticles();
			int nPFOs = jetPFOs.size();
			int n_Nh = 0;
			int n_Ph = 0;
			int n_Ch = 0;
			float E_Nh = 0;
			float E_Ph = 0;
			float E_Ch = 0;
			float E_tot = 0;
			for ( int i_pfo = 0 ; i_pfo < nPFOs ; ++i_pfo )
			{
				ReconstructedParticle *pfo = jetPFOs[ i_pfo ];
				E_tot += pfo->getEnergy();
				if ( 0 != pfo->getCharge() )
				{
					E_Ch += pfo->getEnergy();
					++n_Ch;
				}
				else if ( 22 == pfo->getType() )
				{
					E_Ph += pfo->getEnergy();
					++n_Ph;
				}
				else
				{
					E_Nh += pfo->getEnergy();
					++n_Nh;
				}
			}
			m_nPFOs_NeutralHadrons.push_back( n_Nh );
			m_nPFOs_Photons.push_back( n_Ph );
			m_nPFOs_ChargedPFOs.push_back( n_Ch );
			m_nPFOs_Total.push_back( nPFOs );
			m_PFOEnergy_NeutralHadrons.push_back( E_Nh );
			m_PFOEnergy_Photons.push_back( E_Ph );
			m_PFOEnergy_ChargedPFOs.push_back( E_Ch );
			m_PFOEnergy_Total.push_back( E_tot );
			m_EnergyFraction_NeutralHadrons.push_back( E_Nh / E_tot );
			m_EnergyFraction_Photons.push_back( E_Ph / E_tot );
			m_EnergyFraction_ChargedPFOs.push_back( E_Ch / E_tot );
		}
		streamlog_out(MESSAGE) << "Residuals of jet angles calculated in event " << m_nEvt << std::endl;

	}
	catch(DataNotAvailableException &e)
	{
		streamlog_out(MESSAGE) << "getJetAngleResolution : Input collections not found in event " << m_nEvt << std::endl;
	}

}

void ZHllqq5CFit::processEvent( EVENT::LCEvent *pLCEvent )
{


	this->Clear();
	m_nRun = pLCEvent->getRunNumber();
	m_nEvt = pLCEvent->getEventNumber();
	streamlog_out(MESSAGE) << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "///////////////////////////////////////////////// processing event " << m_nEvt << " in run " << m_nRun << " /////////////////////////////////////////////////////////////////////" << std::endl ;
	streamlog_out(MESSAGE) << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;

	LCCollection *inputJetCollection{};
	LCCollection *inputErrorFlowCollection{};
	LCCollection *HiggsDecayMode{};
	LCCollection *inputLeptonCollection{};
	LCCollection *inputSLDecayCollection{};
	LCCollection *inputNuEnergyB{};
	LCCollection *inputNuEnergyC{};
	std::vector<int> B_index{};
	std::vector<int> C_index{};
	std::vector<float> BHadENuPlus{};
	std::vector<float> BHadENuMinus{};
	std::vector<float> CHadENuPlus{};
	std::vector<float> CHadENuMinus{};
	int m_isDecayedTob = 0;
	int m_isDecayedToc = 0;
	int m_isDecayedTog = 0;
	int m_isDecayedToother = 0;
	float ISR1E_mcp = 0.;
	float ISR2E_mcp = 0.;
	float ISRE_mcp_max = 0.;
	float ISR1pz_mcp = 0.;
	float ISR2pz_mcp = 0.;
	float ISRpz_mcp_max = 0.;
	try
	{
		inputJetCollection = pLCEvent->getCollection( jetcollection );
		inputErrorFlowCollection = pLCEvent->getCollection( errorflowcollection );
		inputLeptonCollection = pLCEvent->getCollection( leptoncollection );
		HiggsDecayMode = pLCEvent->getCollection( hDecayMode );
		m_nJets = inputJetCollection->getNumberOfElements();
		m_nLeptons = inputLeptonCollection->getNumberOfElements();
		m_isDecayedTob = HiggsDecayMode->getParameters().getIntVal("isDecayedTob");
		m_isDecayedToc = HiggsDecayMode->getParameters().getIntVal("isDecayedToc");
		m_isDecayedTog = HiggsDecayMode->getParameters().getIntVal("isDecayedTog");
		m_isDecayedToother = HiggsDecayMode->getParameters().getIntVal("isDecayedToother");
		ISR1E_mcp = HiggsDecayMode->getParameters().getFloatVal("ISR1Energy");
		ISR2E_mcp = HiggsDecayMode->getParameters().getFloatVal("ISR2Energy");
		ISR1pz_mcp = HiggsDecayMode->getParameters().getFloatVal("ISR1Pz");
		ISR2pz_mcp = HiggsDecayMode->getParameters().getFloatVal("ISR2Pz");
		h_nJets->Fill(m_nJets);
		h_nLeptons->Fill(m_nLeptons);
		streamlog_out(DEBUG) << " found " << m_nJets << " jets in event: " << m_nEvt << " in run: " << m_nRun << std::endl;
		streamlog_out(DEBUG) << " found " << m_nLeptons << " isolated leptons in event: " << m_nEvt << " in run: " << m_nRun << std::endl;
		if ( m_isDecayedTob == 1 )
		{
			m_HDecayMode = 1;
		}
		if ( m_isDecayedToc == 1 )
		{
			m_HDecayMode = 2;
		}
		if ( m_isDecayedTog == 1 )
		{
			m_HDecayMode = 3;
		}
		if ( m_isDecayedToother == 1 )
		{
			m_HDecayMode = 4;
		}
		h_HDecayMode->Fill(m_HDecayMode-0.5);
		m_pTTree->Fill();

		if ( m_isDecayedTob == 1 )
		{
			if ( !m_includeHbb ) return;
		}
		if ( m_isDecayedToc == 1 )
		{
			if ( !m_includeHcc ) return;
		}
		if ( m_isDecayedTog == 1 )
		{
			if ( !m_includeHgg ) return;
		}
		if ( m_isDecayedToother == 1 )
		{
			if ( !m_includeHother ) return;
		}
		h_nLeptons_nJets->Fill( m_nLeptons , m_nJets );
		if (m_nJets != 2 || m_nLeptons != 2) return;

		this->getMCISRFSR( pLCEvent );
		this->getHmumu4Momentum( pLCEvent );

		if ( ISR1E_mcp > ISR2E_mcp)
		{
			ISRE_mcp_max = ISR1E_mcp;
			ISRpz_mcp_max = ISR1pz_mcp;
		}
		else
		{
			ISRE_mcp_max = ISR2E_mcp;
			ISRpz_mcp_max = ISR2pz_mcp;
		}
		if (m_fitNuE)
		{
			inputSLDecayCollection = pLCEvent->getCollection( SLDecayCollection );
			inputNuEnergyB = pLCEvent->getCollection( NuEnergyCollectionB );
			inputNuEnergyC = pLCEvent->getCollection( NuEnergyCollectionC );
			m_nSLDecayBHadron = inputSLDecayCollection->getParameters().getIntVal("nBSLD");
			m_nSLDecayCHadron = inputSLDecayCollection->getParameters().getIntVal("nCSLD");
			m_nSLDecayTotal = inputSLDecayCollection->getParameters().getIntVal("nSLD");
			B_index = inputSLDecayCollection->getParameters().getIntVals("BHadronIndex", B_index);
			C_index = inputSLDecayCollection->getParameters().getIntVals("CHadronIndex", C_index);
			BHadENuPlus = inputNuEnergyB->getParameters().getFloatVals("recEnergyENuPlusSLDB", BHadENuPlus);
			BHadENuMinus = inputNuEnergyB->getParameters().getFloatVals("recEnergyENuMinusSLDB", BHadENuMinus);
			CHadENuPlus = inputNuEnergyC->getParameters().getFloatVals("recEnergyENuPlusSLDC", CHadENuPlus);
			CHadENuMinus = inputNuEnergyC->getParameters().getFloatVals("recEnergyENuMinusSLDC", CHadENuMinus);
			streamlog_out(DEBUG)  << "Found " << m_nSLDecayBHadron << " semileptonic decay of B-hadron in event " << m_nEvt << std::endl;
			streamlog_out(DEBUG)  << "Found " << m_nSLDecayCHadron << " semileptonic decay of C-hadron in event " << m_nEvt << std::endl;
		}
		this->getJetAngleResolution( pLCEvent , m_nSLDecayTotal );

		float yminus = inputJetCollection->parameters().getFloatVal( "YMinus" );
		streamlog_out(DEBUG)  << " Yminus = " << yminus << std::endl;
		float yplus = inputJetCollection->parameters().getFloatVal( "YPlus" );
		streamlog_out(DEBUG)  << " Yplus = " << yplus << std::endl;

		for (int i_jet = 0; i_jet < m_nJets; i_jet++)
		{
			ReconstructedParticle* jet_EF = dynamic_cast<ReconstructedParticle*>( inputErrorFlowCollection->getElementAt( i_jet ) ) ;
			m_SigmaPx2.push_back( jet_EF->getCovMatrix()[0] );
			h_SigmaPx2->Fill( jet_EF->getCovMatrix()[0] );

			m_SigmaPxPy.push_back( jet_EF->getCovMatrix()[1] );
			h_SigmaPxPy->Fill( jet_EF->getCovMatrix()[1] );

			m_SigmaPxPz.push_back( jet_EF->getCovMatrix()[3] );
			h_SigmaPxPz->Fill( jet_EF->getCovMatrix()[3] );

			m_SigmaPxE.push_back( jet_EF->getCovMatrix()[6] );
			h_SigmaPxE->Fill( jet_EF->getCovMatrix()[6] );

			m_SigmaPy2.push_back( jet_EF->getCovMatrix()[2] );
			h_SigmaPy2->Fill( jet_EF->getCovMatrix()[2] );

			m_SigmaPyPz.push_back( jet_EF->getCovMatrix()[4] );
			h_SigmaPyPz->Fill( jet_EF->getCovMatrix()[4] );

			m_SigmaPyE.push_back( jet_EF->getCovMatrix()[7] );
			h_SigmaPyE->Fill( jet_EF->getCovMatrix()[7] );

			m_SigmaPz2.push_back( jet_EF->getCovMatrix()[5] );
			h_SigmaPz2->Fill( jet_EF->getCovMatrix()[5] );

			m_SigmaPzE.push_back( jet_EF->getCovMatrix()[8] );
			h_SigmaPzE->Fill( jet_EF->getCovMatrix()[8] );

			m_SigmaE2.push_back( jet_EF->getCovMatrix()[9] );
			h_SigmaE2->Fill( jet_EF->getCovMatrix()[9] );
		}
		for (int i_lep = 0; i_lep < m_nLeptons; i_lep++)
		{
			ReconstructedParticle* lep_EF = dynamic_cast<ReconstructedParticle*>( inputLeptonCollection->getElementAt( i_lep ) ) ;
			m_SigmaPx2.push_back( lep_EF->getCovMatrix()[0] );
			h_SigmaPx2->Fill( lep_EF->getCovMatrix()[0] );

			m_SigmaPxPy.push_back( lep_EF->getCovMatrix()[1] );
			h_SigmaPxPy->Fill( lep_EF->getCovMatrix()[1] );

			m_SigmaPxPz.push_back( lep_EF->getCovMatrix()[3] );
			h_SigmaPxPz->Fill( lep_EF->getCovMatrix()[3] );

			m_SigmaPxE.push_back( lep_EF->getCovMatrix()[6] );
			h_SigmaPxE->Fill( lep_EF->getCovMatrix()[6] );

			m_SigmaPy2.push_back( lep_EF->getCovMatrix()[2] );
			h_SigmaPy2->Fill( lep_EF->getCovMatrix()[2] );

			m_SigmaPyPz.push_back( lep_EF->getCovMatrix()[4] );
			h_SigmaPyPz->Fill( lep_EF->getCovMatrix()[4] );

			m_SigmaPyE.push_back( lep_EF->getCovMatrix()[7] );
			h_SigmaPyE->Fill( lep_EF->getCovMatrix()[7] );

			m_SigmaPz2.push_back( lep_EF->getCovMatrix()[5] );
			h_SigmaPz2->Fill( lep_EF->getCovMatrix()[5] );

			m_SigmaPzE.push_back( lep_EF->getCovMatrix()[8] );
			h_SigmaPzE->Fill( lep_EF->getCovMatrix()[8] );

			m_SigmaE2.push_back( lep_EF->getCovMatrix()[9] );
			h_SigmaE2->Fill( lep_EF->getCovMatrix()[9] );
		}
		streamlog_out(DEBUG) << " SigmaPx2[jet1]  = " << m_SigmaPx2[0] << std::endl ;
		streamlog_out(DEBUG) << " SigmaPx2[jet2]  = " << m_SigmaPx2[1] << std::endl ;
		streamlog_out(DEBUG) << " SigmaPx2[lep1]  = " << m_SigmaPx2[2] << std::endl ;
		streamlog_out(DEBUG) << " SigmaPx2[lep2]  = " << m_SigmaPx2[3] << std::endl ;



//		const int m_nSLDB = m_nSLDecayBHadron;
//		const int m_nSLDC = m_nSLDecayCHadron;
		std::vector<int> JetmatchSLDB{};//const int m_nSLDB = m_nSLDecayBHadron];
		std::vector<int> JetmatchSLDC{};//const int m_nSLDC = m_nSLDecayCHadron];
		for (unsigned int i_SLDB = 0; i_SLDB < B_index.size(); i_SLDB++) JetmatchSLDB.push_back(this->FindMatchingJettoSLD(pLCEvent, B_index[i_SLDB]));
		for (unsigned int i_SLDC = 0; i_SLDC < C_index.size(); i_SLDC++) JetmatchSLDC.push_back(this->FindMatchingJettoSLD(pLCEvent, C_index[i_SLDC]));

		const int n_B_perm = pow( 3 , m_nSLDecayBHadron );
//		std::vector<int> E_nu_B_permut[n_B_perm]{};
		std::vector<std::vector<int>> E_nu_B_permut;
		std::vector<int> rowB;
		for (int E_nu_B_permutation = 0; E_nu_B_permutation < n_B_perm; E_nu_B_permutation++)
		{
			int s = E_nu_B_permutation;
			int m = 0;
			int r = 0;
			for (int i = 0; i < m_nSLDecayBHadron; i++)
			{
				m = s / 3;
				r = s - 3 * m;
				s = m;
				rowB.push_back(r);
//				E_nu_B_permut[E_nu_B_permutation].push_back(r);
			}
			E_nu_B_permut.push_back(rowB);
			rowB.clear();
		}

		int n_C_perm = pow( 3 , m_nSLDecayCHadron );
//		std::vector<int> E_nu_C_permut[n_C_perm]{};
		std::vector<std::vector<int>> E_nu_C_permut;
		std::vector<int> rowC;
		for (int E_nu_C_permutation = 0; E_nu_C_permutation < n_C_perm; E_nu_C_permutation++)
		{
			int s = E_nu_C_permutation;
			int m = 0;
			int r = 0;
			for (int i = 0; i < m_nSLDecayCHadron; i++)
			{
				m = s / 3;
				r = s - 3 * m;
				s = m;
				rowC.push_back(r);
//				E_nu_C_permut[E_nu_C_permutation].push_back(r);
			}
			E_nu_C_permut.push_back(rowC);
			rowC.clear();
		}
		float ENeutrinoJet0_BSLD = 0.;
		float ENeutrinoJet1_BSLD = 0.;
		float ENeutrinoJet0_CSLD = 0.;
		float ENeutrinoJet1_CSLD = 0.;
//		float ENeutrinoJet0 = 0.;
//		float ENeutrinoJet1 = 0.;
		float permuted_ENu_B = 0;
		float permuted_ENu_C = 0;
		TLorentzVector BNutlv(0,0,0,0);
		TLorentzVector Jet0_BNutlv(0,0,0,0);
		TLorentzVector Jet1_BNutlv(0,0,0,0);
		TLorentzVector CNutlv(0,0,0,0);
		TLorentzVector Jet0_CNutlv(0,0,0,0);
		TLorentzVector Jet1_CNutlv(0,0,0,0);
		TLorentzVector Jet0_Nutlv(0,0,0,0);
		TLorentzVector Jet1_Nutlv(0,0,0,0);
		std::vector<std::vector<float>> FitResultwNu{};
		std::vector<std::vector<float>> FitResultwoNu{};

		std::vector<float> fitStartValueswNu;
		std::vector<float> fitOutputswNu;
		std::vector<float> fittedParticleswNu;
		std::vector<float> pullswNu;
		std::vector<float> constraintswNu;
		std::vector<float> uncertaintieswNu;
		std::vector<float> diJetSystemwNu;

		std::vector<float> fitStartValueswNu_bestfit;
		std::vector<float> fitOutputswNu_bestfit;
		std::vector<float> fittedParticleswNu_bestfit;
		std::vector<float> pullswNu_bestfit;
		std::vector<float> constraintswNu_bestfit;
		std::vector<float> uncertaintieswNu_bestfit;
		std::vector<float> diJetSystemwNu_bestfit;

		std::vector<float> fitStartValueswoNu;
		std::vector<float> fitOutputswoNu;
		std::vector<float> fittedParticleswoNu;
		std::vector<float> pullswoNu;
		std::vector<float> constraintswoNu;
		std::vector<float> uncertaintieswoNu;
		std::vector<float> diJetSystemwoNu;

		std::vector<float> fitStartValues_best;
		std::vector<float> fitOutputs_best;
		std::vector<float> fittedParticles_best;
		std::vector<float> pulls_best;
		std::vector<float> constraints_best;
		std::vector<float> uncertainties_best;
		std::vector<float> diJetSystem_best;

		int ierr = 0.;
		float bestfitprob_wNu = 0.;
		float bestfitprob_woNu = 0.;

		std::vector<int> B_permutation{};
		std::vector<int> C_permutation{};
		int best_B_Nu = 0;
		int best_C_Nu = 0;
		float m_pull_jet1_E_wNu_bestfit = 0.;
		float m_pull_jet2_E_wNu_bestfit = 0.;
		float m_pull_jet1_th_wNu_bestfit = 0.;
		float m_pull_jet2_th_wNu_bestfit = 0.;
		float m_pull_jet1_phi_wNu_bestfit = 0.;
		float m_pull_jet2_phi_wNu_bestfit = 0.;
		float m_pull_lepton1_InvpT_wNu_bestfit = 0.;
		float m_pull_lepton2_InvpT_wNu_bestfit = 0.;
		float m_pull_lepton1_th_wNu_bestfit = 0.;
		float m_pull_lepton2_th_wNu_bestfit = 0.;
		float m_pull_lepton1_phi_wNu_bestfit = 0.;
		float m_pull_lepton2_phi_wNu_bestfit = 0.;
		float ISRpx_wNu_bestfit = 0.;
		float ISRpy_wNu_bestfit = 0.;
		float ISRpz_wNu_bestfit = 0.;
		float Zpx_wNu_bestfit = 0.;
		float Zpy_wNu_bestfit = 0.;
		float Zpz_wNu_bestfit = 0.;
		float ZE_wNu_bestfit = 0.;
		float Hpx_wNu_bestfit = 0.;
		float Hpy_wNu_bestfit = 0.;
		float Hpz_wNu_bestfit = 0.;
		float HE_wNu_bestfit = 0.;
		float ISRpx_woNu = 0.;
		float ISRpy_woNu = 0.;
		float ISRpz_woNu = 0.;
		float Zpx_woNu = 0.;
		float Zpy_woNu = 0.;
		float Zpz_woNu = 0.;
		float ZE_woNu = 0.;
		float Hpx_woNu = 0.;
		float Hpy_woNu = 0.;
		float Hpz_woNu = 0.;
		float HE_woNu = 0.;
		float ISRpx_best = 0.;
		float ISRpy_best = 0.;
		float ISRpz_best = 0.;
		float Zpx_best = 0.;
		float Zpy_best = 0.;
		float Zpz_best = 0.;
		float ZE_best = 0.;
		float Hpx_best = 0.;
		float Hpy_best = 0.;
		float Hpz_best = 0.;
		float HE_best = 0.;

		Jet0_Nutlv = TLorentzVector(0.,0.,0.,0.);
		Jet1_Nutlv = TLorentzVector(0.,0.,0.,0.);

		streamlog_out(MESSAGE) << "***********************************************************************************************************************************************" << std::endl;
		streamlog_out(MESSAGE) << "********************************************* Perform fit without neutrino correction *********************************************************" << std::endl;
		streamlog_out(MESSAGE) << "***********************************************************************************************************************************************" << std::endl;
//		FitResultwoNu = this->performFIT( pLCEvent , 0. , 0. );
		FitResultwoNu = this->performFIT( pLCEvent , Jet0_Nutlv , Jet1_Nutlv );
		fitStartValueswoNu = FitResultwoNu[0];
		fitOutputswoNu = FitResultwoNu[1];
		fittedParticleswoNu = FitResultwoNu[2];
		pullswoNu = FitResultwoNu[3];
		constraintswoNu = FitResultwoNu[4];
		uncertaintieswoNu = FitResultwoNu[5];
		diJetSystemwoNu = FitResultwoNu[6];
		streamlog_out(DEBUG)  << "FitResult without neutrino correction received from fit successfully: " << std::endl ;

		m_iError_woNu = fitOutputswoNu[0];
		h_fitError_woNu->Fill(m_iError_woNu);
		if ( m_iError_woNu == 0 )
		{
			m_probability_woNu = fitOutputswoNu[1];
			bestfitprob_woNu = m_probability_woNu;
			m_n_itter_woNu = fitOutputswoNu[2];
			m_startmassZ_woNu = fitOutputswoNu[3];
			m_startmassH_woNu = fitOutputswoNu[4];
			m_beststartmassZ_woNu = fitOutputswoNu[5];
			m_beststartmassH_woNu = fitOutputswoNu[6];
			m_Zmass_after_fit_woNu = fitOutputswoNu[7];
			m_Hmass_after_fit_woNu = fitOutputswoNu[8];
			m_chi2startmassZ_woNu = fitOutputswoNu[9];
			m_chi2startmassH_woNu = fitOutputswoNu[10];
			m_chi2best_woNu = fitOutputswoNu[11];
			m_bestphotonenergy_woNu = fitOutputswoNu[12];
			m_pull_jet_E_woNu.push_back(pullswoNu[0]);
			m_pull_jet_E_woNu.push_back(pullswoNu[1]);
			m_pull_jet_th_woNu.push_back(pullswoNu[2]);
			m_pull_jet_th_woNu.push_back(pullswoNu[3]);
			m_pull_jet_phi_woNu.push_back(pullswoNu[4]);
			m_pull_jet_phi_woNu.push_back(pullswoNu[5]);
			m_pull_lepton_InvpT_woNu.push_back(pullswoNu[6]);
			m_pull_lepton_InvpT_woNu.push_back(pullswoNu[7]);
			m_pull_lepton_th_woNu.push_back(pullswoNu[8]);
			m_pull_lepton_th_woNu.push_back(pullswoNu[9]);
			m_pull_lepton_phi_woNu.push_back(pullswoNu[10]);
			m_pull_lepton_phi_woNu.push_back(pullswoNu[11]);
			ISRpx_woNu = fittedParticleswoNu[0];
			ISRpy_woNu = fittedParticleswoNu[1];
			ISRpz_woNu = fittedParticleswoNu[2];
			Zpx_woNu = fittedParticleswoNu[3];
			Zpy_woNu = fittedParticleswoNu[4];
			Zpz_woNu = fittedParticleswoNu[5];
			ZE_woNu = fittedParticleswoNu[6];
			Hpx_woNu = fittedParticleswoNu[7];
			Hpy_woNu = fittedParticleswoNu[8];
			Hpz_woNu = fittedParticleswoNu[9];
			HE_woNu = fittedParticleswoNu[10];
			m_pxc_before_ISR_woNu = constraintswoNu[0];
			m_pyc_before_ISR_woNu = constraintswoNu[1];
			m_pzc_before_ISR_woNu = constraintswoNu[2];
			m_ec_before_ISR_woNu = constraintswoNu[3];
			m_pxc_before_fit_woNu = constraintswoNu[4];
			m_pyc_before_fit_woNu = constraintswoNu[5];
			m_pzc_before_fit_woNu = constraintswoNu[6];
			m_ec_before_fit_woNu = constraintswoNu[7];
			m_pxc_after_fit_woNu = constraintswoNu[8];
			m_pyc_after_fit_woNu = constraintswoNu[9];
			m_pzc_after_fit_woNu = constraintswoNu[10];
			m_ec_after_fit_woNu = constraintswoNu[11];
			m_zc_before_ISR_woNu = constraintswoNu[12];
			m_zc_before_fit_woNu = constraintswoNu[13];
			m_zc_after_fit_woNu = constraintswoNu[14];
			m_jet_startPx_woNu.push_back(fitStartValueswoNu[0]);
			m_jet_startPx_woNu.push_back(fitStartValueswoNu[1]);
			m_jet_startPy_woNu.push_back(fitStartValueswoNu[2]);
			m_jet_startPy_woNu.push_back(fitStartValueswoNu[3]);
			m_jet_startPz_woNu.push_back(fitStartValueswoNu[4]);
			m_jet_startPz_woNu.push_back(fitStartValueswoNu[5]);
			m_jet_startE_woNu.push_back(fitStartValueswoNu[6]);
			m_jet_startE_woNu.push_back(fitStartValueswoNu[7]);
			m_lepton_startPx_woNu.push_back(fitStartValueswoNu[8]);
			m_lepton_startPx_woNu.push_back(fitStartValueswoNu[9]);
			m_lepton_startPy_woNu.push_back(fitStartValueswoNu[10]);
			m_lepton_startPy_woNu.push_back(fitStartValueswoNu[11]);
			m_lepton_startPz_woNu.push_back(fitStartValueswoNu[12]);
			m_lepton_startPz_woNu.push_back(fitStartValueswoNu[13]);
			m_lepton_startE_woNu.push_back(fitStartValueswoNu[14]);
			m_lepton_startE_woNu.push_back(fitStartValueswoNu[15]);
			m_jet_startTheta_woNu.push_back(diJetSystemwoNu[1]);
			m_jet_startTheta_woNu.push_back(diJetSystemwoNu[7]);
			m_jet_startPhi_woNu.push_back(diJetSystemwoNu[2]);
			m_jet_startPhi_woNu.push_back(diJetSystemwoNu[8]);
			m_jet_fittedE_woNu.push_back(fitOutputswoNu[13]);
			m_jet_fittedE_woNu.push_back(fitOutputswoNu[14]);
			m_jet_fittedTheta_woNu.push_back(fitOutputswoNu[15]);
			m_jet_fittedTheta_woNu.push_back(fitOutputswoNu[16]);
			m_jet_fittedPhi_woNu.push_back(fitOutputswoNu[17]);
			m_jet_fittedPhi_woNu.push_back(fitOutputswoNu[18]);
			m_dijet_angle_woNu = diJetSystemwoNu[12];
		}
		m_jet_SigmaTheta.push_back(uncertaintieswoNu[0]);
		m_jet_SigmaTheta.push_back(uncertaintieswoNu[1]);
		m_jet_SigmaPhi.push_back(uncertaintieswoNu[2]);
		m_jet_SigmaPhi.push_back(uncertaintieswoNu[3]);
		m_jet_SigmaE.push_back(uncertaintieswoNu[4]);
		m_jet_SigmaE.push_back(uncertaintieswoNu[5]);
		m_lepton_SigmaTheta.push_back(uncertaintieswoNu[6]);
		m_lepton_SigmaTheta.push_back(uncertaintieswoNu[7]);
		m_lepton_SigmaPhi.push_back(uncertaintieswoNu[8]);
		m_lepton_SigmaPhi.push_back(uncertaintieswoNu[9]);
		m_lepton_SigmaInvpT.push_back(uncertaintieswoNu[10]);
		m_lepton_SigmaInvpT.push_back(uncertaintieswoNu[11]);
		m_ISR_startPx_woNu = uncertaintieswoNu[12];
		m_ISR_startPy_woNu = uncertaintieswoNu[13];
		m_ISR_startPz_woNu = uncertaintieswoNu[14];
		h_error_jet_E->Fill(m_jet_SigmaE[0]);
		h_error_jet_E->Fill(m_jet_SigmaE[1]);
		h_error_jet_Theta->Fill(m_jet_SigmaTheta[0]);
		h_error_jet_Theta->Fill(m_jet_SigmaTheta[1]);
		h_error_jet_Phi->Fill(m_jet_SigmaPhi[0]);
		h_error_jet_Phi->Fill(m_jet_SigmaPhi[1]);
		h_error_lepton_InvpT->Fill(m_lepton_SigmaInvpT[0]);
		h_error_lepton_InvpT->Fill(m_lepton_SigmaInvpT[1]);
		h_error_lepton_Theta->Fill(m_lepton_SigmaTheta[0]);
		h_error_lepton_Theta->Fill(m_lepton_SigmaTheta[1]);
		h_error_lepton_Phi->Fill(m_lepton_SigmaPhi[0]);
		h_error_lepton_Phi->Fill(m_lepton_SigmaPhi[1]);

		streamlog_out(DEBUG) << "size of FitResult without neutrino = " << FitResultwoNu.size() << endl;
		streamlog_out(DEBUG) << "Fit probability without neutrino correction = " << FitResultwoNu[0][1] << endl;

		for (int i_perm_B = 0; i_perm_B < n_B_perm; i_perm_B++)
		{
			ENeutrinoJet0_BSLD = 0.;
			ENeutrinoJet1_BSLD = 0.;
			permuted_ENu_B = 0;
			BNutlv = TLorentzVector(0.,0.,0.,0.);
			for (unsigned int i_sldB = 0; i_sldB < E_nu_B_permut[i_perm_B].size(); i_sldB++)
			{
				EVENT::MCParticle *BNuPlus = dynamic_cast<EVENT::MCParticle*>(inputNuEnergyB->getElementAt( 2 * i_sldB ));
				EVENT::MCParticle *BNuMinus = dynamic_cast<EVENT::MCParticle*>(inputNuEnergyB->getElementAt( 2 * i_sldB + 1 ));
				if ( E_nu_B_permut[i_perm_B][i_sldB] == 0 )
				{
					permuted_ENu_B = BHadENuMinus[i_sldB];
					BNutlv = TLorentzVector(BNuMinus->getMomentum(),BNuMinus->getEnergy());
				}
				else if ( E_nu_B_permut[i_perm_B][i_sldB] == 1 )
				{
					permuted_ENu_B = BHadENuPlus[i_sldB];
					BNutlv = TLorentzVector(BNuPlus->getMomentum(),BNuPlus->getEnergy());
				}
				else
				{
					permuted_ENu_B = 0.;
					BNutlv = TLorentzVector(0.,0.,0.,0.);
				}
				if ( JetmatchSLDB[i_sldB] == 0 )
				{
					ENeutrinoJet0_BSLD += permuted_ENu_B;
					Jet0_BNutlv += BNutlv;
				}
				else if ( JetmatchSLDB[i_sldB] == 1 )
				{
					ENeutrinoJet1_BSLD += permuted_ENu_B;
					Jet1_BNutlv += BNutlv;
				}
			}
			for (int i_perm_C = 0; i_perm_C < n_C_perm; i_perm_C++)
			{
				ENeutrinoJet0_CSLD = 0;
				ENeutrinoJet1_CSLD = 0;
				permuted_ENu_C = 0;
				CNutlv = TLorentzVector(0.,0.,0.,0.);
				for (unsigned int i_sldC = 0; i_sldC < E_nu_C_permut[i_perm_C].size(); i_sldC++)
				{
					EVENT::MCParticle *CNuPlus = dynamic_cast<EVENT::MCParticle*>(inputNuEnergyC->getElementAt( 2 * i_sldC ));
					EVENT::MCParticle *CNuMinus = dynamic_cast<EVENT::MCParticle*>(inputNuEnergyC->getElementAt( 2 * i_sldC + 1 ));
					if ( E_nu_C_permut[i_perm_C][i_sldC] == 0 )
					{
						permuted_ENu_C = CHadENuMinus[i_sldC];
						CNutlv = TLorentzVector(CNuMinus->getMomentum(),CNuMinus->getEnergy());
					}
					else if ( E_nu_C_permut[i_perm_C][i_sldC] == 1 )
					{
						permuted_ENu_C = CHadENuPlus[i_sldC];
						CNutlv = TLorentzVector(CNuPlus->getMomentum(),CNuPlus->getEnergy());
					}
					else
					{
						permuted_ENu_C = 0.;
					}
					if ( JetmatchSLDC[i_sldC] == 0 )
					{
						ENeutrinoJet0_CSLD += permuted_ENu_C;
						Jet0_CNutlv += CNutlv;
					}
					else if ( JetmatchSLDC[i_sldC] == 1 )
					{
						ENeutrinoJet1_CSLD += permuted_ENu_C;
						Jet1_CNutlv += CNutlv;
					}
				}
//				ENeutrinoJet0 = ENeutrinoJet0_BSLD + ENeutrinoJet0_CSLD;
//				ENeutrinoJet1 = ENeutrinoJet1_BSLD + ENeutrinoJet1_CSLD;
				Jet0_Nutlv = Jet0_BNutlv + Jet0_CNutlv;
				Jet1_Nutlv = Jet1_BNutlv + Jet1_CNutlv;
//				if ( ENeutrinoJet0 == 0. && ENeutrinoJet1 == 0. ) continue;
				streamlog_out(MESSAGE) << "******************************************************************************************************************************************************************" << std::endl;
				streamlog_out(MESSAGE) << "************* Perform fit with Jet0_Nutlv = (" << Jet0_Nutlv.Px() << " , " << Jet0_Nutlv.Py() << " , " << Jet0_Nutlv.Pz() << " , " << Jet0_Nutlv.E() << ") and Jet1_Nutlv = (" << Jet1_Nutlv.Px() << " , " << Jet1_Nutlv.Py() << " , " << Jet1_Nutlv.Pz() << " , " << Jet1_Nutlv.E() << ") *************" << std::endl;
				streamlog_out(MESSAGE) << "******************************************************************************************************************************************************************" << std::endl;
//				FitResultwNu = this->performFIT( pLCEvent , ENeutrinoJet0 , ENeutrinoJet1 );
				FitResultwNu = this->performFIT( pLCEvent , Jet0_Nutlv , Jet1_Nutlv );
				fitStartValueswNu = FitResultwNu[0];
				fitOutputswNu = FitResultwNu[1];
				fittedParticleswNu = FitResultwNu[2];
				pullswNu = FitResultwNu[3];
				constraintswNu = FitResultwNu[4];
				uncertaintieswNu = FitResultwNu[5];
				diJetSystemwNu = FitResultwNu[6];
				streamlog_out(DEBUG)  << "FitResult with neutrino correction received from fit successfully: " << std::endl ;

				Jet0_BNutlv = TLorentzVector(0.,0.,0.,0.);
				Jet0_CNutlv = TLorentzVector(0.,0.,0.,0.);
				Jet1_BNutlv = TLorentzVector(0.,0.,0.,0.);
				Jet1_CNutlv = TLorentzVector(0.,0.,0.,0.);
				Jet0_Nutlv = TLorentzVector(0.,0.,0.,0.);
				Jet1_Nutlv = TLorentzVector(0.,0.,0.,0.);
				m_iError_wNu.push_back(fitOutputswNu[0]);
				ierr = fitOutputswNu[0];
				h_fitError_wNu->Fill(ierr);
				h_ErrorCode_wNu_woNu->Fill( ierr , m_iError_woNu );
				if ( ierr == 0 )//&& ( ENeutrinoJet0 != 0 || ENeutrinoJet1 !=0 ) )
				{
					m_probability_wNu.push_back(fitOutputswNu[1]);
					m_n_itter_wNu.push_back(fitOutputswNu[2]);
					m_startmassZ_wNu.push_back(fitOutputswNu[3]);
					m_startmassH_wNu.push_back(fitOutputswNu[4]);
					m_beststartmassZ_wNu.push_back(fitOutputswNu[5]);
					m_beststartmassH_wNu.push_back(fitOutputswNu[6]);
					m_Zmass_after_fit_wNu.push_back(fitOutputswNu[7]);
					m_Hmass_after_fit_wNu.push_back(fitOutputswNu[8]);
					m_chi2startmassZ_wNu.push_back(fitOutputswNu[9]);
					m_chi2startmassH_wNu.push_back(fitOutputswNu[10]);
					m_chi2_wNu.push_back(fitOutputswNu[11]);
					m_bestphotonenergy_wNu.push_back(fitOutputswNu[12]);
					m_pull_jet_E_wNu.push_back(pullswNu[0]);
					m_pull_jet_E_wNu.push_back(pullswNu[1]);
					m_pull_jet_th_wNu.push_back(pullswNu[2]);
					m_pull_jet_th_wNu.push_back(pullswNu[3]);
					m_pull_jet_phi_wNu.push_back(pullswNu[4]);
					m_pull_jet_phi_wNu.push_back(pullswNu[5]);
					m_pull_lepton_InvpT_wNu.push_back(pullswNu[6]);
					m_pull_lepton_InvpT_wNu.push_back(pullswNu[7]);
					m_pull_lepton_th_wNu.push_back(pullswNu[8]);
					m_pull_lepton_th_wNu.push_back(pullswNu[9]);
					m_pull_lepton_phi_wNu.push_back(pullswNu[10]);
					m_pull_lepton_phi_wNu.push_back(pullswNu[11]);
					m_jet_startPx_wNu.push_back(fitStartValueswNu[0]);
					m_jet_startPx_wNu.push_back(fitStartValueswNu[1]);
					m_jet_startPy_wNu.push_back(fitStartValueswNu[2]);
					m_jet_startPy_wNu.push_back(fitStartValueswNu[3]);
					m_jet_startPz_wNu.push_back(fitStartValueswNu[4]);
					m_jet_startPz_wNu.push_back(fitStartValueswNu[5]);
					m_jet_startE_wNu.push_back(fitStartValueswNu[6]);
					m_jet_startE_wNu.push_back(fitStartValueswNu[7]);
					m_jet_startTheta_wNu.push_back(diJetSystemwNu[1]);
					m_jet_startTheta_wNu.push_back(diJetSystemwNu[7]);
					m_jet_startPhi_wNu.push_back(diJetSystemwNu[2]);
					m_jet_startPhi_wNu.push_back(diJetSystemwNu[8]);
					m_jet_fittedE_wNu.push_back(fitOutputswNu[13]);
					m_jet_fittedE_wNu.push_back(fitOutputswNu[14]);
					m_jet_fittedTheta_wNu.push_back(fitOutputswNu[15]);
					m_jet_fittedTheta_wNu.push_back(fitOutputswNu[16]);
					m_jet_fittedPhi_wNu.push_back(fitOutputswNu[17]);
					m_jet_fittedPhi_wNu.push_back(fitOutputswNu[18]);
					m_dijet_angle_wNu.push_back(diJetSystemwNu[12]);
					m_lepton_startPx_wNu.push_back(fitStartValueswNu[8]);
					m_lepton_startPx_wNu.push_back(fitStartValueswNu[9]);
					m_lepton_startPy_wNu.push_back(fitStartValueswNu[10]);
					m_lepton_startPy_wNu.push_back(fitStartValueswNu[11]);
					m_lepton_startPz_wNu.push_back(fitStartValueswNu[12]);
					m_lepton_startPz_wNu.push_back(fitStartValueswNu[13]);
					m_lepton_startE_wNu.push_back(fitStartValueswNu[14]);
					m_lepton_startE_wNu.push_back(fitStartValueswNu[15]);
				}
				if ( fitOutputswNu[0] == 0 && bestfitprob_wNu <= fitOutputswNu[1] )
				{
					best_B_Nu = i_perm_B;
					best_C_Nu = i_perm_C;
					m_iError_wNu_bestfit = fitOutputswNu[0];
					m_probability_wNu_bestfit = fitOutputswNu[1];
					bestfitprob_wNu = m_probability_wNu_bestfit;
					m_n_itter_wNu_bestfit = fitOutputswNu[2];
					m_startmassZ_wNu_bestfit = fitOutputswNu[3];
					m_startmassH_wNu_bestfit = fitOutputswNu[4];
					m_beststartmassZ_wNu_bestfit = fitOutputswNu[5];
					m_beststartmassH_wNu_bestfit = fitOutputswNu[6];
					m_Zmass_after_fit_wNu_bestfit = fitOutputswNu[7];
					m_Hmass_after_fit_wNu_bestfit = fitOutputswNu[8];
					m_chi2startmassZ_wNu_bestfit = fitOutputswNu[9];
					m_chi2startmassH_wNu_bestfit = fitOutputswNu[10];
					m_chi2_wNu_bestfit = fitOutputswNu[11];
					m_bestphotonenergy_wNu_bestfit = fitOutputswNu[12];
					m_pull_jet1_E_wNu_bestfit = pullswNu[0];
					m_pull_jet2_E_wNu_bestfit = pullswNu[1];
					m_pull_jet1_th_wNu_bestfit = pullswNu[2];
					m_pull_jet2_th_wNu_bestfit = pullswNu[3];
					m_pull_jet1_phi_wNu_bestfit = pullswNu[4];
					m_pull_jet2_phi_wNu_bestfit = pullswNu[5];
					m_pull_lepton1_InvpT_wNu_bestfit = pullswNu[6];
					m_pull_lepton2_InvpT_wNu_bestfit = pullswNu[7];
					m_pull_lepton1_th_wNu_bestfit = pullswNu[8];
					m_pull_lepton2_th_wNu_bestfit = pullswNu[9];
					m_pull_lepton1_phi_wNu_bestfit = pullswNu[10];
					m_pull_lepton2_phi_wNu_bestfit = pullswNu[11];
					ISRpx_wNu_bestfit = fittedParticleswNu[0];
					ISRpy_wNu_bestfit = fittedParticleswNu[1];
					ISRpz_wNu_bestfit = fittedParticleswNu[2];
					Zpx_wNu_bestfit = fittedParticleswNu[3];
					Zpy_wNu_bestfit = fittedParticleswNu[4];
					Zpz_wNu_bestfit = fittedParticleswNu[5];
					ZE_wNu_bestfit = fittedParticleswNu[6];
					Hpx_wNu_bestfit = fittedParticleswNu[7];
					Hpy_wNu_bestfit = fittedParticleswNu[8];
					Hpz_wNu_bestfit = fittedParticleswNu[9];
					HE_wNu_bestfit = fittedParticleswNu[10];
					m_pxc_before_ISR_wNu = constraintswNu[0];
					m_pyc_before_ISR_wNu = constraintswNu[1];
					m_pzc_before_ISR_wNu = constraintswNu[2];
					m_ec_before_ISR_wNu = constraintswNu[3];
					m_pxc_before_fit_wNu = constraintswNu[4];
					m_pyc_before_fit_wNu = constraintswNu[5];
					m_pzc_before_fit_wNu = constraintswNu[6];
					m_ec_before_fit_wNu = constraintswNu[7];
					m_pxc_after_fit_wNu = constraintswNu[8];
					m_pyc_after_fit_wNu = constraintswNu[9];
					m_pzc_after_fit_wNu = constraintswNu[10];
					m_ec_after_fit_wNu = constraintswNu[11];
					m_zc_before_ISR_wNu = constraintswNu[12];
					m_zc_before_fit_wNu = constraintswNu[13];
					m_zc_after_fit_wNu = constraintswNu[14];
					m_jet_startPx_wNu_bestfit.push_back(fitStartValueswNu[0]);
					m_jet_startPx_wNu_bestfit.push_back(fitStartValueswNu[1]);
					m_jet_startPy_wNu_bestfit.push_back(fitStartValueswNu[2]);
					m_jet_startPy_wNu_bestfit.push_back(fitStartValueswNu[3]);
					m_jet_startPz_wNu_bestfit.push_back(fitStartValueswNu[4]);
					m_jet_startPz_wNu_bestfit.push_back(fitStartValueswNu[5]);
					m_jet_startE_wNu_bestfit.push_back(fitStartValueswNu[6]);
					m_jet_startE_wNu_bestfit.push_back(fitStartValueswNu[7]);
					m_jet_startTheta_wNu_bestfit.push_back(diJetSystemwNu[1]);
					m_jet_startTheta_wNu_bestfit.push_back(diJetSystemwNu[7]);
					m_jet_startPhi_wNu_bestfit.push_back(diJetSystemwNu[2]);
					m_jet_startPhi_wNu_bestfit.push_back(diJetSystemwNu[8]);
					m_jet_fittedE_wNu_bestfit.push_back(fitOutputswNu[13]);
					m_jet_fittedE_wNu_bestfit.push_back(fitOutputswNu[14]);
					m_jet_fittedTheta_wNu_bestfit.push_back(fitOutputswNu[15]);
					m_jet_fittedTheta_wNu_bestfit.push_back(fitOutputswNu[16]);
					m_jet_fittedPhi_wNu_bestfit.push_back(fitOutputswNu[17]);
					m_jet_fittedPhi_wNu_bestfit.push_back(fitOutputswNu[18]);
					m_dijet_angle_wNu_bestfit = diJetSystemwNu[12];
					m_lepton_startPx_wNu_bestfit.push_back(fitStartValueswNu[8]);
					m_lepton_startPx_wNu_bestfit.push_back(fitStartValueswNu[9]);
					m_lepton_startPy_wNu_bestfit.push_back(fitStartValueswNu[10]);
					m_lepton_startPy_wNu_bestfit.push_back(fitStartValueswNu[11]);
					m_lepton_startPz_wNu_bestfit.push_back(fitStartValueswNu[12]);
					m_lepton_startPz_wNu_bestfit.push_back(fitStartValueswNu[13]);
					m_lepton_startE_wNu_bestfit.push_back(fitStartValueswNu[14]);
					m_lepton_startE_wNu_bestfit.push_back(fitStartValueswNu[15]);
				}
				m_ISR_startPx_wNu = fitStartValueswNu[16];
				m_ISR_startPy_wNu = fitStartValueswNu[17];
				m_ISR_startPz_wNu = fitStartValueswNu[18];
				streamlog_out(DEBUG) << "size of FitResult with neutrino = " << FitResultwNu.size() << endl;
//				streamlog_out(DEBUG) << "Fit best probability with neutrino correction = " << FitResultwNu[1] << endl;
				FitResultwNu.clear();
			}
		}
		if ( m_iError_wNu_bestfit == 0 )
		{
			m_pull_jet_E_wNu_bestfit.push_back(m_pull_jet1_E_wNu_bestfit);
			m_pull_jet_E_wNu_bestfit.push_back(m_pull_jet2_E_wNu_bestfit);
			m_pull_jet_th_wNu_bestfit.push_back(m_pull_jet1_th_wNu_bestfit);
			m_pull_jet_th_wNu_bestfit.push_back(m_pull_jet2_th_wNu_bestfit);
			m_pull_jet_phi_wNu_bestfit.push_back(m_pull_jet1_phi_wNu_bestfit);
			m_pull_jet_phi_wNu_bestfit.push_back(m_pull_jet2_phi_wNu_bestfit);
			m_pull_lepton_InvpT_wNu_bestfit.push_back(m_pull_lepton1_InvpT_wNu_bestfit);
			m_pull_lepton_InvpT_wNu_bestfit.push_back(m_pull_lepton2_InvpT_wNu_bestfit);
			m_pull_lepton_th_wNu_bestfit.push_back(m_pull_lepton1_th_wNu_bestfit);
			m_pull_lepton_th_wNu_bestfit.push_back(m_pull_lepton2_th_wNu_bestfit);
			m_pull_lepton_phi_wNu_bestfit.push_back(m_pull_lepton1_phi_wNu_bestfit);
			m_pull_lepton_phi_wNu_bestfit.push_back(m_pull_lepton2_phi_wNu_bestfit);
		}
		std::vector<int> best_E_nu_permut{};
		int ss = best_C_Nu;
		int mm = 0;
		int rr = 0;
		for (int i = 0; i < m_nSLDecayCHadron; i++)
		{
			mm = ss / 3;
			rr = ss - 3 * mm;
			ss = mm;
			if (rr == 1)
			{
				best_E_nu_permut.push_back( i + 1 );
			}
			else if ( rr == 0 )
			{
				best_E_nu_permut.push_back( -1 - i );
			}
			else
			{
				best_E_nu_permut.push_back( 0 );
			}
		}
		ss = best_B_Nu;
		mm = 0;
		rr = 0;
		for (int i = 0; i < m_nSLDecayBHadron; i++)
		{
			mm = ss / 3;
			rr = ss - 3 * mm;
			ss = mm;
			if (rr == 1)
			{
				best_E_nu_permut.push_back( m_nSLDecayCHadron + i + 1 );
			}
			else if ( rr == 0 )
			{
				best_E_nu_permut.push_back( -1 - i - m_nSLDecayCHadron);
			}
			else
			{
				best_E_nu_permut.push_back( 0 );
			}
		}
		for (int i = 0; i < m_nSLDecayBHadron + m_nSLDecayCHadron ; i++)
		{
			m_bestNuCombination.push_back(best_E_nu_permut[m_nSLDecayBHadron + m_nSLDecayCHadron - i - 1]);
		}

//		m_iError_best = m_iError_wNu_bestfit;

		if ( m_iError_wNu_bestfit == 0 && bestfitprob_woNu < bestfitprob_wNu )
		{
			m_iError_best = m_iError_wNu_bestfit;
			m_probability_best = m_probability_wNu_bestfit;
			m_n_itter_best = m_n_itter_wNu_bestfit;
			m_startmassZ_best = m_startmassZ_wNu_bestfit;
			m_startmassH_best = m_startmassH_wNu_bestfit;
			m_beststartmassZ_best = m_beststartmassZ_wNu_bestfit;
			m_beststartmassH_best = m_beststartmassH_wNu_bestfit;
			m_Zmass_after_fit_best = m_Zmass_after_fit_wNu_bestfit;
			m_Hmass_after_fit_best = m_Hmass_after_fit_wNu_bestfit;
			m_chi2startmassZ_best = m_chi2startmassZ_wNu_bestfit;
			m_chi2startmassH_best = m_chi2startmassH_wNu_bestfit;
			m_chi2_best = m_chi2_wNu_bestfit;
			m_bestphotonenergy_best = m_bestphotonenergy_wNu_bestfit;
			m_pull_jet_E_best.push_back(m_pull_jet_E_wNu_bestfit[0]);
			m_pull_jet_E_best.push_back(m_pull_jet_E_wNu_bestfit[1]);
			h_pull_jet_E->Fill(m_pull_jet_E_wNu_bestfit[0]);
			h_pull_jet_E->Fill(m_pull_jet_E_wNu_bestfit[1]);

			m_pull_jet_th_best.push_back(m_pull_jet_th_wNu_bestfit[0]);
			m_pull_jet_th_best.push_back(m_pull_jet_th_wNu_bestfit[1]);
			h_pull_jet_theta->Fill(m_pull_jet_th_wNu_bestfit[0]);
			h_pull_jet_theta->Fill(m_pull_jet_th_wNu_bestfit[1]);

			m_pull_jet_phi_best.push_back(m_pull_jet_phi_wNu_bestfit[0]);
			m_pull_jet_phi_best.push_back(m_pull_jet_phi_wNu_bestfit[1]);
			h_pull_jet_phi->Fill(m_pull_jet_phi_wNu_bestfit[0]);
			h_pull_jet_phi->Fill(m_pull_jet_phi_wNu_bestfit[1]);

			m_pull_lepton_InvpT_best.push_back(m_pull_lepton_InvpT_wNu_bestfit[0]);
			m_pull_lepton_InvpT_best.push_back(m_pull_lepton_InvpT_wNu_bestfit[1]);
			h_pull_lepton_InvPt->Fill(m_pull_lepton_InvpT_wNu_bestfit[0]);
			h_pull_lepton_InvPt->Fill(m_pull_lepton_InvpT_wNu_bestfit[1]);

			m_pull_lepton_th_best.push_back(m_pull_lepton_th_wNu_bestfit[0]);
			m_pull_lepton_th_best.push_back(m_pull_lepton_th_wNu_bestfit[1]);
			h_pull_lepton_theta->Fill(m_pull_lepton_th_wNu_bestfit[0]);
			h_pull_lepton_theta->Fill(m_pull_lepton_th_wNu_bestfit[1]);

			m_pull_lepton_phi_best.push_back(m_pull_lepton_phi_wNu_bestfit[0]);
			m_pull_lepton_phi_best.push_back(m_pull_lepton_phi_wNu_bestfit[1]);
			h_pull_lepton_phi->Fill(m_pull_lepton_phi_wNu_bestfit[0]);
			h_pull_lepton_phi->Fill(m_pull_lepton_phi_wNu_bestfit[1]);
			m_jet_startPx_best.push_back(m_jet_startPx_wNu_bestfit[0]);
			m_jet_startPx_best.push_back(m_jet_startPx_wNu_bestfit[1]);
			m_jet_startPy_best.push_back(m_jet_startPy_wNu_bestfit[0]);
			m_jet_startPy_best.push_back(m_jet_startPy_wNu_bestfit[1]);
			m_jet_startPz_best.push_back(m_jet_startPz_wNu_bestfit[0]);
			m_jet_startPz_best.push_back(m_jet_startPz_wNu_bestfit[1]);
			m_jet_startE_best.push_back(m_jet_startE_wNu_bestfit[0]);
			m_jet_startE_best.push_back(m_jet_startE_wNu_bestfit[1]);
			m_jet_startTheta_best.push_back(m_jet_startTheta_wNu_bestfit[0]);
			m_jet_startTheta_best.push_back(m_jet_startTheta_wNu_bestfit[1]);
			m_jet_startPhi_best.push_back(m_jet_startPhi_wNu_bestfit[0]);
			m_jet_startPhi_best.push_back(m_jet_startPhi_wNu_bestfit[1]);
			m_jet_fittedE_best.push_back(m_jet_fittedE_wNu_bestfit[0]);
			m_jet_fittedE_best.push_back(m_jet_fittedE_wNu_bestfit[1]);
			m_jet_fittedTheta_best.push_back(m_jet_fittedTheta_wNu_bestfit[0]);
			m_jet_fittedTheta_best.push_back(m_jet_fittedTheta_wNu_bestfit[1]);
			m_jet_fittedPhi_best.push_back(m_jet_fittedPhi_wNu_bestfit[0]);
			m_jet_fittedPhi_best.push_back(m_jet_fittedPhi_wNu_bestfit[1]);
			m_dijet_angle_best = m_dijet_angle_wNu_bestfit;
			m_lepton_startPx_best.push_back(m_lepton_startPx_wNu_bestfit[0]);
			m_lepton_startPx_best.push_back(m_lepton_startPx_wNu_bestfit[1]);
			m_lepton_startPy_best.push_back(m_lepton_startPy_wNu_bestfit[0]);
			m_lepton_startPy_best.push_back(m_lepton_startPy_wNu_bestfit[1]);
			m_lepton_startPz_best.push_back(m_lepton_startPz_wNu_bestfit[0]);
			m_lepton_startPz_best.push_back(m_lepton_startPz_wNu_bestfit[1]);
			m_lepton_startE_best.push_back(m_lepton_startE_wNu_bestfit[0]);
			m_lepton_startE_best.push_back(m_lepton_startE_wNu_bestfit[1]);

			ISRpx_best = ISRpx_wNu_bestfit;
			ISRpy_best = ISRpy_wNu_bestfit;
			ISRpz_best = ISRpz_wNu_bestfit;
			Zpx_best = Zpx_wNu_bestfit;
			Zpy_best = Zpy_wNu_bestfit;
			Zpz_best = Zpz_wNu_bestfit;
			ZE_best = ZE_wNu_bestfit;
			Hpx_best = Hpx_wNu_bestfit;
			Hpy_best = Hpy_wNu_bestfit;
			Hpz_best = Hpz_wNu_bestfit;
			HE_best = HE_wNu_bestfit;
			m_pxc_before_ISR_best = m_pxc_before_ISR_wNu;
			m_pyc_before_ISR_best = m_pyc_before_ISR_wNu;
			m_pzc_before_ISR_best = m_pzc_before_ISR_wNu;
			m_ec_before_ISR_best = m_ec_before_ISR_wNu;
			m_zc_before_ISR_best = m_zc_before_ISR_wNu;
			m_pxc_before_fit_best = m_pxc_before_fit_wNu;
			m_pyc_before_fit_best = m_pyc_before_fit_wNu;
			m_pzc_before_fit_best = m_pzc_before_fit_wNu;
			m_ec_before_fit_best = m_ec_before_fit_wNu;
			m_zc_before_fit_best = m_zc_before_fit_wNu;
			m_pxc_after_fit_best = m_pxc_after_fit_wNu;
			m_pyc_after_fit_best = m_pyc_after_fit_wNu;
			m_pzc_after_fit_best = m_pzc_after_fit_wNu;
			m_ec_after_fit_best = m_ec_after_fit_wNu;
			m_zc_after_fit_best = m_zc_after_fit_wNu;
		}
		else if ( m_iError_woNu == 0 )
		{
			m_iError_best = m_iError_woNu;
			m_probability_best = m_probability_woNu;
			m_n_itter_best = m_n_itter_woNu;
			m_startmassZ_best = m_startmassZ_woNu;
			m_startmassH_best = m_startmassH_woNu;
			m_beststartmassZ_best = m_beststartmassZ_woNu;
			m_beststartmassH_best = m_beststartmassH_woNu;
			m_Zmass_after_fit_best = m_Zmass_after_fit_woNu;
			m_Hmass_after_fit_best = m_Hmass_after_fit_woNu;
			m_chi2startmassZ_best = m_chi2startmassZ_woNu;
			m_chi2startmassH_best = m_chi2startmassH_woNu;
			m_chi2_best = m_chi2best_woNu;
			m_bestphotonenergy_best = m_bestphotonenergy_woNu;
			m_pull_jet_E_best.push_back(m_pull_jet_E_woNu[0]);
			m_pull_jet_E_best.push_back(m_pull_jet_E_woNu[1]);
			h_pull_jet_E->Fill(m_pull_jet_E_woNu[0]);
			h_pull_jet_E->Fill(m_pull_jet_E_woNu[1]);

			m_pull_jet_th_best.push_back(m_pull_jet_th_woNu[0]);
			m_pull_jet_th_best.push_back(m_pull_jet_th_woNu[1]);
			h_pull_jet_theta->Fill(m_pull_jet_th_woNu[0]);
			h_pull_jet_theta->Fill(m_pull_jet_th_woNu[1]);

			m_pull_jet_phi_best.push_back(m_pull_jet_phi_woNu[0]);
			m_pull_jet_phi_best.push_back(m_pull_jet_phi_woNu[1]);
			h_pull_jet_phi->Fill(m_pull_jet_phi_woNu[0]);
			h_pull_jet_phi->Fill(m_pull_jet_phi_woNu[1]);

			m_pull_lepton_InvpT_best.push_back(m_pull_lepton_InvpT_woNu[0]);
			m_pull_lepton_InvpT_best.push_back(m_pull_lepton_InvpT_woNu[1]);
			h_pull_lepton_InvPt->Fill(m_pull_lepton_InvpT_woNu[0]);
			h_pull_lepton_InvPt->Fill(m_pull_lepton_InvpT_woNu[1]);

			m_pull_lepton_th_best.push_back(m_pull_lepton_th_woNu[0]);
			m_pull_lepton_th_best.push_back(m_pull_lepton_th_woNu[1]);
			h_pull_lepton_theta->Fill(m_pull_lepton_th_woNu[0]);
			h_pull_lepton_theta->Fill(m_pull_lepton_th_woNu[1]);

			m_pull_lepton_phi_best.push_back(m_pull_lepton_phi_woNu[0]);
			m_pull_lepton_phi_best.push_back(m_pull_lepton_phi_woNu[1]);
			h_pull_lepton_phi->Fill(m_pull_lepton_phi_woNu[0]);
			h_pull_lepton_phi->Fill(m_pull_lepton_phi_woNu[1]);
			m_jet_startPx_best.push_back(m_jet_startPx_woNu[0]);
			m_jet_startPx_best.push_back(m_jet_startPx_woNu[1]);
			m_jet_startPy_best.push_back(m_jet_startPy_woNu[0]);
			m_jet_startPy_best.push_back(m_jet_startPy_woNu[1]);
			m_jet_startPz_best.push_back(m_jet_startPz_woNu[0]);
			m_jet_startPz_best.push_back(m_jet_startPz_woNu[1]);
			m_jet_startE_best.push_back(m_jet_startE_woNu[0]);
			m_jet_startE_best.push_back(m_jet_startE_woNu[1]);
			m_jet_startTheta_best.push_back(m_jet_startTheta_woNu[0]);
			m_jet_startTheta_best.push_back(m_jet_startTheta_woNu[1]);
			m_jet_startPhi_best.push_back(m_jet_startPhi_woNu[0]);
			m_jet_startPhi_best.push_back(m_jet_startPhi_woNu[1]);
			m_jet_fittedE_best.push_back(m_jet_fittedE_woNu[0]);
			m_jet_fittedE_best.push_back(m_jet_fittedE_woNu[1]);
			m_jet_fittedTheta_best.push_back(m_jet_fittedTheta_woNu[0]);
			m_jet_fittedTheta_best.push_back(m_jet_fittedTheta_woNu[1]);
			m_jet_fittedPhi_best.push_back(m_jet_fittedPhi_woNu[0]);
			m_jet_fittedPhi_best.push_back(m_jet_fittedPhi_woNu[1]);
			m_dijet_angle_best = m_dijet_angle_woNu;
			m_lepton_startPx_best.push_back(m_lepton_startPx_woNu[0]);
			m_lepton_startPx_best.push_back(m_lepton_startPx_woNu[1]);
			m_lepton_startPy_best.push_back(m_lepton_startPy_woNu[0]);
			m_lepton_startPy_best.push_back(m_lepton_startPy_woNu[1]);
			m_lepton_startPz_best.push_back(m_lepton_startPz_woNu[0]);
			m_lepton_startPz_best.push_back(m_lepton_startPz_woNu[1]);
			m_lepton_startE_best.push_back(m_lepton_startE_woNu[0]);
			m_lepton_startE_best.push_back(m_lepton_startE_woNu[1]);

			m_bestNuCombination.push_back(0);
			ISRpx_best = ISRpx_woNu;
			ISRpy_best = ISRpy_woNu;
			ISRpz_best = ISRpz_woNu;
			Zpx_best = Zpx_woNu;
			Zpy_best = Zpy_woNu;
			Zpz_best = Zpz_woNu;
			ZE_best = ZE_woNu;
			Hpx_best = Hpx_woNu;
			Hpy_best = Hpy_woNu;
			Hpz_best = Hpz_woNu;
			HE_best = HE_woNu;
			m_pxc_before_ISR_best = m_pxc_before_ISR_woNu;
			m_pyc_before_ISR_best = m_pyc_before_ISR_woNu;
			m_pzc_before_ISR_best = m_pzc_before_ISR_woNu;
			m_ec_before_ISR_best = m_ec_before_ISR_woNu;
			m_zc_before_ISR_best = m_zc_before_ISR_woNu;
			m_pxc_before_fit_best = m_pxc_before_fit_woNu;
			m_pyc_before_fit_best = m_pyc_before_fit_woNu;
			m_pzc_before_fit_best = m_pzc_before_fit_woNu;
			m_ec_before_fit_best = m_ec_before_fit_woNu;
			m_zc_before_fit_best = m_zc_before_fit_woNu;
			m_pxc_after_fit_best = m_pxc_after_fit_woNu;
			m_pyc_after_fit_best = m_pyc_after_fit_woNu;
			m_pzc_after_fit_best = m_pzc_after_fit_woNu;
			m_ec_after_fit_best = m_ec_after_fit_woNu;
			m_zc_after_fit_best = m_zc_after_fit_woNu;
		}
//		if ( m_iError_woNu == 0 && m_iError_wNu_bestfit == 0 )
		if ( m_iError_woNu == 0 )
		{
			h_Zmass_beforefit_woNu->Fill(m_startmassZ_woNu);
			h_Hmass_beforefit_woNu->Fill(m_startmassH_woNu);
			h_Zmass_afterfit_woNu->Fill(m_Zmass_after_fit_woNu);
			h_Hmass_afterfit_woNu->Fill(m_Hmass_after_fit_woNu);
			h_fitProbability_woNu->Fill(m_probability_woNu);
			h_ISR_pzc_woNu->Fill( -1 * m_pzc_before_ISR_woNu , ISRpz_woNu );
		}
		if ( m_iError_wNu_bestfit == 0 )
		{
			h_Zmass_beforefit_wNu->Fill(m_startmassZ_wNu_bestfit);
			h_Hmass_beforefit_wNu->Fill(m_startmassH_wNu_bestfit);
			h_fitProbability_wNu->Fill(m_probability_wNu_bestfit);
			h_ISR_pzc_wNu->Fill( -1 * m_pzc_before_ISR_wNu , ISRpz_wNu_bestfit );
		}
		if ( m_iError_best == 0 )
		{
			streamlog_out(DEBUG) << "Fit best probability = " << FitResultwoNu[0][1] << endl;
			h_Zmass_afterfit_wNu->Fill(m_Zmass_after_fit_best);
			h_Hmass_afterfit_wNu->Fill(m_Hmass_after_fit_best);
			h_fitProbability_best->Fill(m_probability_best);
			h_ISR_pzc_best->Fill( -1 * m_pzc_before_ISR_best , ISRpz_best );
			h_ISRE_1mcp_fit->Fill( ISRE_mcp_max , m_bestphotonenergy_best );
			h_ISRpz_1mcp_fit->Fill( ISRpz_mcp_max , ISRpz_best );
			h_ISRE_2mcp_fit->Fill( ISR1E_mcp + ISR2E_mcp , m_bestphotonenergy_best );
			h_ISRpz_2mcp_fit->Fill( ISR1pz_mcp + ISR2pz_mcp , ISRpz_best );

			h_fitProbability_diJetAngle->Fill( m_dijet_angle_best , m_probability_best );
			h_fitProbability_Ejet->Fill( m_jet_startE_best[0] , m_probability_best );
			h_fitProbability_Ejet->Fill( m_jet_startE_best[1] , m_probability_best );
			h_fitProbability_Thetajet->Fill( m_jet_startTheta_best[0] * 45 / atan(1) , m_probability_best );
			h_fitProbability_Thetajet->Fill( m_jet_startTheta_best[1] * 45 / atan(1) , m_probability_best );
			h_fitProbability_Phijet->Fill( m_jet_startPhi_best[0] * 45 / atan(1) , m_probability_best );
			h_fitProbability_Phijet->Fill( m_jet_startPhi_best[1] * 45 / atan(1) , m_probability_best );
			h_fitProbability_SigmaEjet->Fill( m_jet_SigmaE[0] , m_probability_best );
			h_fitProbability_SigmaEjet->Fill( m_jet_SigmaE[1] , m_probability_best );
			h_fitProbability_SigmaThetajet->Fill( m_jet_SigmaTheta[0] * m_JetResThetaCoeff , m_probability_best );
			h_fitProbability_SigmaThetajet->Fill( m_jet_SigmaTheta[1] * m_JetResThetaCoeff , m_probability_best );
			h_fitProbability_SigmaPhijet->Fill( m_jet_SigmaPhi[0] * m_JetResPhiCoeff , m_probability_best );
			h_fitProbability_SigmaPhijet->Fill( m_jet_SigmaPhi[1] * m_JetResPhiCoeff , m_probability_best );
			h_fitProbability_pullEjet->Fill( m_pull_jet_E_best[0] , m_probability_best );
			h_fitProbability_pullEjet->Fill( m_pull_jet_E_best[1] , m_probability_best );
			h_fitProbability_pullThetajet->Fill( m_pull_jet_th_best[0] , m_probability_best );
			h_fitProbability_pullThetajet->Fill( m_pull_jet_th_best[1] , m_probability_best );
			h_fitProbability_pullPhijet->Fill( m_pull_jet_phi_best[0] , m_probability_best );
			h_fitProbability_pullPhijet->Fill( m_pull_jet_phi_best[1] , m_probability_best );
			h_fitProbability_constraintPx->Fill( m_pxc_before_ISR_wNu , m_probability_best );
			h_fitProbability_constraintPy->Fill( m_pyc_before_ISR_wNu , m_probability_best );
			h_fitProbability_constraintPz->Fill( m_pzc_before_ISR_wNu , m_probability_best );
			h_fitProbability_constraintE->Fill( m_ec_before_ISR_wNu , m_probability_best );

			m_TotalSigmaPx2 = m_SigmaPx2[0] + m_SigmaPx2[1] + m_SigmaPx2[2] + m_SigmaPx2[3];
			m_TotalSigmaPy2 = m_SigmaPy2[0] + m_SigmaPy2[1] + m_SigmaPy2[2] + m_SigmaPy2[3];
			m_TotalSigmaPz2 = m_SigmaPz2[0] + m_SigmaPz2[1] + m_SigmaPz2[2] + m_SigmaPz2[3];
			m_TotalSigmaE2 = m_SigmaE2[0] + m_SigmaE2[1] + m_SigmaE2[2] + m_SigmaE2[3];
			m_TotalSigmaPxPy = m_SigmaPxPy[0] + m_SigmaPxPy[1] + m_SigmaPxPy[2] + m_SigmaPxPy[3];
			m_TotalSigmaPxPz = m_SigmaPxPz[0] + m_SigmaPxPz[1] + m_SigmaPxPz[2] + m_SigmaPxPz[3];
			m_TotalSigmaPxE = m_SigmaPxE[0] + m_SigmaPxE[1] + m_SigmaPxE[2] + m_SigmaPxE[3];
			m_TotalSigmaPyPz = m_SigmaPyPz[0] + m_SigmaPyPz[1] + m_SigmaPyPz[2] + m_SigmaPyPz[3];
			m_TotalSigmaPyE = m_SigmaPyE[0] + m_SigmaPyE[1] + m_SigmaPyE[2] + m_SigmaPyE[3];
			m_TotalSigmaPzE = m_SigmaPzE[0] + m_SigmaPzE[1] + m_SigmaPzE[2] + m_SigmaPzE[3];

			streamlog_out(DEBUG) << " SigmaPx2[jet1]  = " << m_SigmaPx2[0] << std::endl ;
			streamlog_out(DEBUG) << " SigmaPx2[jet2]  = " << m_SigmaPx2[1] << std::endl ;
			streamlog_out(DEBUG) << " SigmaPx2[lep1]  = " << m_SigmaPx2[2] << std::endl ;
			streamlog_out(DEBUG) << " SigmaPx2[lep2]  = " << m_SigmaPx2[3] << std::endl ;
			streamlog_out(DEBUG) << " TotalSigmaPx2[lep1]  = " << m_TotalSigmaPx2 << std::endl ;

			h_fitProbability_sigmaPx->Fill( std::sqrt( m_TotalSigmaPx2 ) , m_probability_best );
			h_fitProbability_sigmaPy->Fill( std::sqrt( m_TotalSigmaPy2 ) , m_probability_best );
			h_fitProbability_sigmaPz->Fill( std::sqrt( m_TotalSigmaPz2 ) , m_probability_best );
			h_fitProbability_sigmaE->Fill( std::sqrt( m_TotalSigmaE2 ) , m_probability_best );

			h_constraintPx_uncertaintyPx->Fill( m_pxc_before_ISR_wNu , std::sqrt( m_TotalSigmaPx2 ) );
			h_constraintPy_uncertaintyPy->Fill( m_pyc_before_ISR_wNu , std::sqrt( m_TotalSigmaPy2 ) );
			h_constraintPz_uncertaintyPz->Fill( m_pzc_before_ISR_wNu , std::sqrt( m_TotalSigmaPz2 ) );
			h_constraintE_uncertaintyE->Fill( m_ec_before_ISR_wNu , std::sqrt( m_TotalSigmaE2 ) );

			h_recE_fitE->Fill( ( m_jet_fittedE_best[0] - m_jet_startE_best[0] ) / m_jet_startE_best[0] );
			h_recE_fitE->Fill( ( m_jet_fittedE_best[1] - m_jet_startE_best[1] ) / m_jet_startE_best[1] );
			h_recPhifitPhi_recThetafitTheta->Fill( m_jet_fittedPhi_best[0] - m_jet_startPhi_best[0] , m_jet_fittedTheta_best[0] - m_jet_startTheta_best[0] );
			h_recPhifitPhi_recThetafitTheta->Fill( m_jet_fittedPhi_best[1] - m_jet_startPhi_best[1] , m_jet_fittedTheta_best[1] - m_jet_startTheta_best[1] );
			h_recPhifitPhi_recThetafitTheta_SigmaThetaPhi->Fill( ( m_jet_fittedPhi_best[0] - m_jet_startPhi_best[0] ) / m_jet_SigmaPhi[0] , ( m_jet_fittedTheta_best[0] - m_jet_startTheta_best[0] ) / m_jet_SigmaTheta[0] );
			h_recPhifitPhi_recThetafitTheta_SigmaThetaPhi->Fill( ( m_jet_fittedPhi_best[1] - m_jet_startPhi_best[1] ) / m_jet_SigmaPhi[1] , ( m_jet_fittedTheta_best[1] - m_jet_startTheta_best[1] ) / m_jet_SigmaTheta[1] );

			if ( m_probability_best < 0.1 )
			{
				h_constraintPx_lowFitProb->Fill( m_pxc_before_ISR_wNu );
				h_constraintPy_lowFitProb->Fill( m_pyc_before_ISR_wNu );
				h_constraintPz_lowFitProb->Fill( m_pzc_before_ISR_wNu );
				h_constraintE_lowFitProb->Fill( m_ec_before_ISR_wNu );
				h_constraintPx_uncertaintyPx_lowFitProb->Fill( m_pxc_before_ISR_wNu , std::sqrt( m_TotalSigmaPx2 ) );
				h_constraintPy_uncertaintyPy_lowFitProb->Fill( m_pyc_before_ISR_wNu , std::sqrt( m_TotalSigmaPy2 ) );
				h_constraintPz_uncertaintyPz_lowFitProb->Fill( m_pzc_before_ISR_wNu , std::sqrt( m_TotalSigmaPz2 ) );
				h_constraintE_uncertaintyE_lowFitProb->Fill( m_ec_before_ISR_wNu , std::sqrt( m_TotalSigmaE2 ) );
				h_pull_jet_theta_jet_phi_lowFitProb->Fill( m_pull_jet_th_best[0] , m_pull_jet_phi_best[0] );
				h_pull_jet_theta_jet_phi_lowFitProb->Fill( m_pull_jet_th_best[1] , m_pull_jet_phi_best[1] );
				h_recPhifitPhi_recThetafitTheta_lowfit->Fill( m_jet_fittedPhi_best[0] - m_jet_startPhi_best[0] , m_jet_fittedTheta_best[0] - m_jet_startTheta_best[0] );
				h_recPhifitPhi_recThetafitTheta_lowfit->Fill( m_jet_fittedPhi_best[1] - m_jet_startPhi_best[1] , m_jet_fittedTheta_best[1] - m_jet_startTheta_best[1] );
				h_recPhifitPhi_recThetafitTheta_lowfit_SigmaThetaPhi->Fill( ( m_jet_fittedPhi_best[0] - m_jet_startPhi_best[0] ) / m_jet_SigmaPhi[0] , ( m_jet_fittedTheta_best[0] - m_jet_startTheta_best[0] ) / m_jet_SigmaTheta[0] );
				h_recPhifitPhi_recThetafitTheta_lowfit_SigmaThetaPhi->Fill( ( m_jet_fittedPhi_best[1] - m_jet_startPhi_best[1] ) / m_jet_SigmaPhi[1] , ( m_jet_fittedTheta_best[1] - m_jet_startTheta_best[1] ) / m_jet_SigmaTheta[1] );
			}
			else if ( 0.1 <= m_probability_best && m_probability_best < 0.9 )
			{
				h_constraintPx_midFitProb->Fill( m_pxc_before_ISR_wNu );
				h_constraintPy_midFitProb->Fill( m_pyc_before_ISR_wNu );
				h_constraintPz_midFitProb->Fill( m_pzc_before_ISR_wNu );
				h_constraintE_midFitProb->Fill( m_ec_before_ISR_wNu );
				h_constraintPx_uncertaintyPx_midFitProb->Fill( m_pxc_before_ISR_wNu , std::sqrt( m_TotalSigmaPx2 ) );
				h_constraintPy_uncertaintyPy_midFitProb->Fill( m_pyc_before_ISR_wNu , std::sqrt( m_TotalSigmaPy2 ) );
				h_constraintPz_uncertaintyPz_midFitProb->Fill( m_pzc_before_ISR_wNu , std::sqrt( m_TotalSigmaPz2 ) );
				h_constraintE_uncertaintyE_midFitProb->Fill( m_ec_before_ISR_wNu , std::sqrt( m_TotalSigmaE2 ) );
				h_pull_jet_theta_jet_phi_midFitProb->Fill( m_pull_jet_th_best[0] , m_pull_jet_phi_best[0] );
				h_pull_jet_theta_jet_phi_midFitProb->Fill( m_pull_jet_th_best[1] , m_pull_jet_phi_best[1] );
				h_recPhifitPhi_recThetafitTheta_midfit->Fill( m_jet_fittedPhi_best[0] - m_jet_startPhi_best[0] , m_jet_fittedTheta_best[0] - m_jet_startTheta_best[0] );
				h_recPhifitPhi_recThetafitTheta_midfit->Fill( m_jet_fittedPhi_best[1] - m_jet_startPhi_best[1] , m_jet_fittedTheta_best[1] - m_jet_startTheta_best[1] );
				h_recPhifitPhi_recThetafitTheta_midfit_SigmaThetaPhi->Fill( ( m_jet_fittedPhi_best[0] - m_jet_startPhi_best[0] ) / m_jet_SigmaPhi[0] , ( m_jet_fittedTheta_best[0] - m_jet_startTheta_best[0] ) / m_jet_SigmaTheta[0] );
				h_recPhifitPhi_recThetafitTheta_midfit_SigmaThetaPhi->Fill( ( m_jet_fittedPhi_best[1] - m_jet_startPhi_best[1] ) / m_jet_SigmaPhi[1] , ( m_jet_fittedTheta_best[1] - m_jet_startTheta_best[1] ) / m_jet_SigmaTheta[1] );
			}
			else
			{
				h_constraintPx_highFitProb->Fill( m_pxc_before_ISR_wNu );
				h_constraintPy_highFitProb->Fill( m_pyc_before_ISR_wNu );
				h_constraintPz_highFitProb->Fill( m_pzc_before_ISR_wNu );
				h_constraintE_highFitProb->Fill( m_ec_before_ISR_wNu );
				h_constraintPx_uncertaintyPx_highFitProb->Fill( m_pxc_before_ISR_wNu , std::sqrt( m_TotalSigmaPx2 ) );
				h_constraintPy_uncertaintyPy_highFitProb->Fill( m_pyc_before_ISR_wNu , std::sqrt( m_TotalSigmaPy2 ) );
				h_constraintPz_uncertaintyPz_highFitProb->Fill( m_pzc_before_ISR_wNu , std::sqrt( m_TotalSigmaPz2 ) );
				h_constraintE_uncertaintyE_highFitProb->Fill( m_ec_before_ISR_wNu , std::sqrt( m_TotalSigmaE2 ) );
				h_pull_jet_theta_jet_phi_highFitProb->Fill( m_pull_jet_th_best[0] , m_pull_jet_phi_best[0] );
				h_pull_jet_theta_jet_phi_highFitProb->Fill( m_pull_jet_th_best[1] , m_pull_jet_phi_best[1] );
				h_recPhifitPhi_recThetafitTheta_highfit->Fill( m_jet_fittedPhi_best[0] - m_jet_startPhi_best[0] , m_jet_fittedTheta_best[0] - m_jet_startTheta_best[0] );
				h_recPhifitPhi_recThetafitTheta_highfit->Fill( m_jet_fittedPhi_best[1] - m_jet_startPhi_best[1] , m_jet_fittedTheta_best[1] - m_jet_startTheta_best[1] );
				h_recPhifitPhi_recThetafitTheta_highfit_SigmaThetaPhi->Fill( ( m_jet_fittedPhi_best[0] - m_jet_startPhi_best[0] ) / m_jet_SigmaPhi[0] , ( m_jet_fittedTheta_best[0] - m_jet_startTheta_best[0] ) / m_jet_SigmaTheta[0] );
				h_recPhifitPhi_recThetafitTheta_highfit_SigmaThetaPhi->Fill( ( m_jet_fittedPhi_best[1] - m_jet_startPhi_best[1] ) / m_jet_SigmaPhi[1] , ( m_jet_fittedTheta_best[1] - m_jet_startTheta_best[1] ) / m_jet_SigmaTheta[1] );
			}

			LCCollectionVec *OutputCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
			ReconstructedParticleImpl* ISRfitrec = new ReconstructedParticleImpl;
			ReconstructedParticleImpl* Zfitrec = new ReconstructedParticleImpl;
			ReconstructedParticleImpl* Hfitrec = new ReconstructedParticleImpl;

			ISRmomentum[0] = ISRpx_best;
			ISRmomentum[1] = ISRpy_best;
			ISRmomentum[2] = ISRpz_best;
			ISRfitrec->setMomentum(ISRmomentum);
			ISRfitrec->setEnergy(m_bestphotonenergy_best);
			ISRfitrec->setType (22);
			streamlog_out(DEBUG) << " Energy ISR:   " << ISRfitrec->getEnergy() << std::endl ;
			streamlog_out(DEBUG) << " IS ISR:   " << ISRfitrec->getType() << std::endl ;
			OutputCol->addElement(ISRfitrec);

			Zmomentum[0] = Zpx_best;
			Zmomentum[1] = Zpy_best;
			Zmomentum[2] = Zpz_best;
			Zfitrec->setMomentum(Zmomentum);
			Zfitrec->setEnergy(ZE_best);
			Zfitrec->setMass(m_Zmass_after_fit_best);
			Zfitrec->setType (23);
			streamlog_out(DEBUG) << "  Zmomentum :   " << Zfitrec->getMomentum()[0] << "," << Zfitrec->getMomentum()[1]<<","<< Zfitrec->getMomentum()[2] << std::endl ;
			streamlog_out(DEBUG) << " Energy Z:   "	<< Zfitrec->getEnergy() << std::endl ;
			streamlog_out(DEBUG) << " Mass Z:   " << Zfitrec->getMass() << std::endl ;
			streamlog_out(DEBUG) << " IS Z :   " << Zfitrec->getType() << std::endl ;
			OutputCol->addElement(Zfitrec);

			Hmomentum[0] = Hpx_best;
			Hmomentum[1] = Hpy_best;
			Hmomentum[2] = Hpz_best;
			Hfitrec->setMomentum(Hmomentum);
			Hfitrec->setEnergy(HE_best);
			Hfitrec->setMass(m_Hmass_after_fit_best);
			Hfitrec->setType (25);
			streamlog_out(DEBUG) << " Energy H:   "	<< Hfitrec->getEnergy() << std::endl ;
			streamlog_out(DEBUG) << " Mass H:   " << Hfitrec->getMass() << std::endl ;
			streamlog_out(DEBUG) << " IS H:   " << Hfitrec->getType() << std::endl ;
			OutputCol->addElement(Hfitrec);

			pLCEvent->addCollection( OutputCol, outputFitcollection.c_str() );


			OutputCol->parameters().setValue("bestchisq", (float)m_chi2_best);
			streamlog_out(DEBUG) << " chi2:   " << m_chi2_best << std::endl ;
			OutputCol->parameters().setValue("best_prob", (float)m_probability_best);
			streamlog_out(DEBUG) << " prob:   " << m_probability_best << std::endl ;
			OutputCol->parameters().setValue("error_code", (int)m_iError_best);
			streamlog_out(DEBUG) << "Error Code:   " << m_iError_best << std::endl ;

			m_pTTree_0->Fill();
			m_pTTree_1->Fill();
			m_pTTree_2->Fill();
			m_pTTree_4->Fill();
			m_pTTree_6->Fill();
		}
		m_pTTree_3->Fill();
		FitResultwoNu.clear();
		m_nEvtSum++;
		m_nEvt++ ;
	}
	catch(DataNotAvailableException &e)
	{
		streamlog_out(MESSAGE) << "processEvent : Input collections not found in event " << m_nEvt << std::endl;
	}

}

int ZHllqq5CFit::FindMatchingJettoSLD(EVENT::LCEvent *pLCEvent, int had_index)
{
	LCCollection *inputJetCollection = pLCEvent->getCollection( jetcollection );
	LCCollection *inputMCPCollection = pLCEvent->getCollection( MCPCollection );
	int nJets = inputJetCollection->getNumberOfElements();
	MCParticle* mcpHadron = dynamic_cast<MCParticle*>( inputMCPCollection->getElementAt(had_index));
	MCParticleVec mcpDaughters = mcpHadron->getDaughters();
	HepLorentzVector lepton_3momentum;
	HepLorentzVector jetvec;
	float lepton_theta = 0.;
	float lepton_phi = 0.;
	float match_theta = 100.;
	float match_phi = 100.;

	int jetIndex_thetaMatch = -1;
	int jetIndex_phiMatch = -1;
	int jetIndex_bestMatch = -1;

	for(unsigned int i_MCPD = 0; i_MCPD < mcpDaughters.size(); i_MCPD++)
	{
		if((std::abs(mcpDaughters[i_MCPD]->getPDG()) == 11) || (std::abs(mcpDaughters[i_MCPD]->getPDG()) == 13) || (std::abs(mcpDaughters[i_MCPD]->getPDG()) == 15))
		{
			lepton_3momentum = HepLorentzVector( mcpDaughters[i_MCPD]->getMomentum()[0] , mcpDaughters[i_MCPD]->getMomentum()[1] , mcpDaughters[i_MCPD]->getMomentum()[2] , mcpDaughters[i_MCPD]->getEnergy() );
			lepton_theta = lepton_3momentum.theta();
			lepton_phi = lepton_3momentum.phi();
		}
	}
	for (int i_jet = 0; i_jet < nJets; i_jet++)
	{
		ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>( inputJetCollection->getElementAt( i_jet ) ) ;
		if (jet)
		{
			jetvec = HepLorentzVector ( (jet->getMomentum())[0] , (jet->getMomentum())[1] , (jet->getMomentum())[2] , jet->getEnergy() );
			delta_theta[i_jet] = jetvec.theta() - lepton_theta;
			delta_phi[i_jet] = jetvec.phi() - lepton_phi;
			if ( std::abs(delta_theta[i_jet]) < match_theta )
			{
				match_theta = std::abs(delta_theta[i_jet]);
				jetIndex_thetaMatch = i_jet;
			}
			if ( std::abs(delta_phi[i_jet]) < match_phi )
			{
				match_phi = std::abs(delta_phi[i_jet]);
				jetIndex_phiMatch = i_jet;
			}
		}

	}
	if ( jetIndex_thetaMatch != jetIndex_phiMatch )
	{
		if (abs(delta_theta[jetIndex_phiMatch]) < 0.5 && abs(delta_phi[jetIndex_phiMatch]) < 0.5)   // 1
		{
			if ( abs(delta_theta[jetIndex_thetaMatch]) > 0.5 && abs(delta_phi[jetIndex_thetaMatch]) > 0.5 ) //2
			{
				jetIndex_bestMatch = jetIndex_phiMatch;
			}
			else  //3
			{
				if( std::abs(delta_theta[jetIndex_phiMatch]) < std::abs(delta_phi[jetIndex_thetaMatch]) )
				{
					jetIndex_bestMatch = jetIndex_phiMatch;
				}
				if( std::abs(delta_theta[jetIndex_phiMatch]) > std::abs(delta_phi[jetIndex_thetaMatch]) )
				{
					jetIndex_bestMatch = jetIndex_thetaMatch;
				}
			}
		}
		else if (abs(delta_theta[jetIndex_thetaMatch]) < 0.5 && abs(delta_phi[jetIndex_thetaMatch]) <0.5 ) //4
		{
			jetIndex_bestMatch = jetIndex_thetaMatch;
		}
		else //5
		{
			if( std::abs(delta_theta[jetIndex_phiMatch]) < std::abs(delta_phi[jetIndex_thetaMatch]) )
			{
				jetIndex_bestMatch = jetIndex_phiMatch;
			}
			if( std::abs(delta_theta[jetIndex_phiMatch]) > std::abs(delta_phi[jetIndex_thetaMatch]) )
			{
				jetIndex_bestMatch = jetIndex_thetaMatch;
			}
		}
		streamlog_out(DEBUG)  << "BESTJET FINALLY IS   " << jetIndex_bestMatch <<std::endl;
	}
	else
	{
		jetIndex_bestMatch = jetIndex_thetaMatch;
		streamlog_out(DEBUG)  << "BESTJET FINALLY IS   " << jetIndex_bestMatch <<std::endl;
	}

	return jetIndex_bestMatch;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////					perform FIT
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<float>> ZHllqq5CFit::performFIT(EVENT::LCEvent *pLCEvent, TLorentzVector Jet0_Nutlv, TLorentzVector Jet1_Nutlv)
{
	Setvalues();
	std::vector<std::vector<float>> FitResult{};
	std::vector<float> fitStartValues;
	std::vector<float> fitOutputs;
	std::vector<float> fittedParticles;
	std::vector<float> pulls;
	std::vector<float> constraints;
	std::vector<float> uncertainties;
	std::vector<float> diJetSystem;
	LCCollection *inputJetCollection = pLCEvent->getCollection( jetcollection );
	LCCollection *inputErrorFlowCollection = pLCEvent->getCollection( errorflowcollection );
	LCCollection *inputLeptonCollection = pLCEvent->getCollection( leptoncollection );

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

	JetFitObject *jet[nJets];
	LeptonFitObject *lepton[nLeptons];
	int nSingularCovMatrix = 0;
	HepLorentzVector jetvec;
	HepLorentzVector Nuvec;
	HepLorentzVector leptonvec;
	Hep3Vector jet0_unit;
	Hep3Vector jet1_unit;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////					Set JetFitObjects
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for (int i_jet = 0; i_jet < nJets; i_jet++)
	{
		ReconstructedParticle *j;
		if ( m_useErrorFlow )
		{
			j = dynamic_cast<ReconstructedParticle*>( inputErrorFlowCollection->getElementAt( i_jet ) );
		}
		else
		{
			j = dynamic_cast<ReconstructedParticle*>( inputJetCollection->getElementAt( i_jet ) );
		}
		if ( i_jet == 0 ) Nuvec = HepLorentzVector( Jet0_Nutlv.Px() , Jet0_Nutlv.Py() , Jet0_Nutlv.Pz() , Jet0_Nutlv.E() );
		if ( i_jet == 1 ) Nuvec = HepLorentzVector( Jet1_Nutlv.Px() , Jet1_Nutlv.Py() , Jet1_Nutlv.Pz() , Jet1_Nutlv.E() );

		Px =	j->getMomentum()[0];
		Px2 =	std::pow( Px , 2 );
		Py =	j->getMomentum()[1];
		Py2 =	std::pow( Py , 2 );
		Pz =	j->getMomentum()[2];
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
			jet0_E = j->getEnergy();
			jet0_Theta = startJetMomentum.Theta();
			jet0_Phi = startJetMomentum.Phi();
		}
		else if ( i_jet == 1 )
		{
			jet1_Px = Px;
			jet1_Py = Py;
			jet1_Pz = Pz;
			jet1_E = j->getEnergy();
			jet1_Theta = startJetMomentum.Theta();
			jet1_Phi = startJetMomentum.Phi();
		}

		SigPx2 =	j->getCovMatrix()[0];
		SigPxSigPy =	j->getCovMatrix()[1];
		SigPy2 =	j->getCovMatrix()[2];
		SigPxSigPz =	j->getCovMatrix()[3];
		SigPySigPz =	j->getCovMatrix()[4];
		SigPz2 =	j->getCovMatrix()[5];
		SigE2 =		j->getCovMatrix()[9];

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
			JetResTheta = m_jetThetaError;
			JetResPhi = m_jetPhiError;
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

		if (m_useErrorFlow)
		{
			jet[i_jet] = new JetFitObject ( jetvec.e(), jetvec.theta() , jetvec.phi(),
					JetResE * sigmaScaleFactor , JetResTheta , JetResPhi , jetvec.m() );
			streamlog_out(DEBUG)  << " start four-vector of jet[" << i_jet << "]: " << *jet[i_jet]  << std::endl ;
			if ( i_jet == 0 )
			{
				jet[i_jet]->setName("Jet0");
			}
			else if (i_jet == 1)
			{
				jet[i_jet]->setName("Jet1");
			}
		}
		else
		{
			jet[i_jet] = new JetFitObject ( jetvec.e(), jetvec.theta() , jetvec.phi(),
					JetEnergyResolution( jetvec.e() ) , m_jetThetaError , m_jetPhiError , jetvec.m() );
			if ( i_jet == 0 )
			{
				jet[i_jet]->setName("Jet0");
			}
			else
			{
				jet[i_jet]->setName("Jet1");
			}
			streamlog_out(DEBUG)  << " start four-vector of jet[" << i_jet << "]: " << *jet[i_jet]  << std::endl ;
		}
		diJetSystem.push_back(jetvec.e());
		diJetSystem.push_back(jetvec.theta());
		diJetSystem.push_back(jetvec.phi());
		diJetSystem.push_back(JetResE * sigmaScaleFactor);
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
			theta_err = m_jetThetaError;
			phi = leptonvec.phi();
			phi_err = m_jetPhiError;
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
			<< leptonvec.theta() <<" +- " << m_jetThetaError << " , "
			<< leptonvec.phi() <<" +- " << m_jetPhiError << std::endl ;

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
	JetFitObject *jets[NJETS];
	for (int i = 0; i < NJETS; ++i)
	{
		jets[i] = &fitjets[i];
		streamlog_out(MESSAGE)  << "start four-vector of jet " << i << ": " << *(jets[i])  << std::endl ;
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
	for (int i = 0; i < NJETS; ++i) pxc.addToFOList(*(jets[i]));
	for (int i = 0; i < NLEPTONS; ++i) pxc.addToFOList(*(leptons[i]));

	MomentumConstraint pyc (0, 0, 1, 0, 0);
	pyc.setName("sum(p_y)");
	for (int i = 0; i < NJETS; ++i) pyc.addToFOList(*(jets[i]));
	for (int i = 0; i < NLEPTONS; ++i) pyc.addToFOList(*(leptons[i]));

	MomentumConstraint pzc (0, 0, 0, 1, 0);
	pzc.setName("sum(p_z)");
	for (int i = 0; i < NJETS; ++i) pzc.addToFOList(*(jets[i]));
	for (int i = 0; i < NLEPTONS; ++i) pzc.addToFOList(*(leptons[i]));

	E_lab= 2 * sqrt( std::pow( 0.548579909e-3 , 2 ) + std::pow( m_ECM / 2 , 2 ) + std::pow( target_p_due_crossing_angle , 2 ) + 0. + 0.);
	MomentumConstraint ec(1, 0, 0, 0, E_lab);
	ec.setName("sum(E)");
	for (int i = 0; i < NJETS; ++i) ec.addToFOList(*(jets[i]));
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
	h.addToFOList(*(jets[0]), 1);
	h.addToFOList(*(jets[1]), 1);

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
	for (int i = 0; i < NJETS; ++i) fitter.addFitObject(*(jets[i]));
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

	prob = fitter.fit();
	double chi2 = fitter.getChi2();
	nit = fitter.getIterations();

	streamlog_out(DEBUG) << "fit probability = " << prob << std::endl ;
	streamlog_out(DEBUG) << "fit chi2 = " << chi2  << std::endl ;
	streamlog_out(DEBUG) << "error code: " << fitter.getError() << std::endl ;

	for (int i = 0; i < NJETS; ++i)
	{
		streamlog_out(MESSAGE)  << "final four-vector of jet " << i << ": " << *(jets[i]) << std::endl ;
		streamlog_out(DEBUG)  << "final px of jet " << i << ": " << (jets[i]) << std::endl ;
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
					fitted = jets[ifo]->getParam(ipar);
					start = startjets[ifo].getParam(ipar);
					errfit = jets[ifo]->getError(ipar);
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
					fittedjet0_E = jets[ifo]->getParam(0);
					fittedjet0_Theta = jets[ifo]->getParam(1);
					fittedjet0_Phi = jets[ifo]->getParam(2);
				}
				else
				{
					fittedjet1_E = jets[ifo]->getParam(0);
					fittedjet1_Theta = jets[ifo]->getParam(1);
					fittedjet1_Phi = jets[ifo]->getParam(2);
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
			Hmomentum[0] = jets[0]->getPx() + jets[1]->getPx();
			Hmomentum[1] = jets[0]->getPy() + jets[1]->getPy();
			Hmomentum[2] = jets[0]->getPz() + jets[1]->getPz();
			Z_Energy = leptons[0]->getE() + leptons[1]->getE();
			H_Energy = jets[0]->getE() + jets[1]->getE();
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

	Zmass_before_fit=beststartmassZ;
	Zmass_after_fit=bestmassZ;
	Hmass_before_fit=beststartmassH;
	Hmass_after_fit=bestmassH;
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

void ZHllqq5CFit::check( LCEvent* )
{
//	nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ZHllqq5CFit::end()
{
	streamlog_out(MESSAGE) << "# of events: " << m_nEvt << std::endl;
//	streamlog_out(ERROR) << "# of nucorrection: " << correction<< std::endl;
//	streamlog_out(ERROR) << "# of Covariance failed: " << nCo<< std::endl;

	m_pTFile->cd();
	m_pTTree->Write();
	m_pTTree_0->Write();
	m_pTTree_1->Write();
	m_pTTree_2->Write();
	m_pTTree_3->Write();
	m_pTTree_4->Write();
	m_pTTree_5->Write();
	m_pTTree_6->Write();
	h_HDecayMode->Write();
	h_Zmass_beforefit_woNu->Write();
	h_Hmass_beforefit_woNu->Write();
	h_Zmass_beforefit_wNu->Write();
	h_Hmass_beforefit_wNu->Write();
	h_Zmass_afterfit_woNu->Write();
	h_Hmass_afterfit_woNu->Write();
	h_Zmass_afterfit_wNu->Write();
	h_Hmass_afterfit_wNu->Write();
	h_fitError_wNu->Write();
	h_fitError_woNu->Write();
	h_ErrorCode_wNu_woNu->Write();
	h_fitProbability_wNu->Write();
	h_fitProbability_woNu->Write();
	h_fitProbability_best->Write();
	h_nJets->Write();
	h_nLeptons->Write();
	h_nLeptons_nJets->Write();
	h_ISRE_1mcp_fit->Write();
	h_ISRE_2mcp_fit->Write();
	h_ISRpz_1mcp_fit->Write();
	h_ISRpz_2mcp_fit->Write();
	h_ISR_pzc_woNu->Write();
	h_ISR_pzc_wNu->Write();
	h_ISR_pzc_best->Write();
	h_mcpISRPx->Write();
	h_mcpISRPy->Write();
	h_mcpISRPz->Write();
	h_mcpISRE->Write();
	h_pull_jet_E->Write();
	h_pull_jet_theta->Write();
	h_pull_jet_phi->Write();
	h_pull_lepton_InvPt->Write();
	h_pull_lepton_theta->Write();
	h_pull_lepton_phi->Write();
	h_error_jet_E->Write();
	h_error_jet_Theta->Write();
	h_error_jet_Phi->Write();
	h_error_lepton_InvpT->Write();
	h_error_lepton_Theta->Write();
	h_error_lepton_Phi->Write();
	h_SigmaPx2->Write();
	h_SigmaPy2->Write();
	h_SigmaPz2->Write();
	h_SigmaE2->Write();
	h_SigmaPxPy->Write();
	h_SigmaPxPz->Write();
	h_SigmaPxE->Write();
	h_SigmaPyPz->Write();
	h_SigmaPyE->Write();
	h_SigmaPzE->Write();
	h_fitProbability_diJetAngle->Write();
	h_fitProbability_Ejet->Write();
	h_fitProbability_Thetajet->Write();
	h_fitProbability_Phijet->Write();
	h_fitProbability_SigmaEjet->Write();
	h_fitProbability_SigmaThetajet->Write();
	h_fitProbability_SigmaPhijet->Write();
	h_fitProbability_pullEjet->Write();
	h_fitProbability_pullThetajet->Write();
	h_fitProbability_pullPhijet->Write();
	h_fitProbability_constraintPx->Write();
	h_fitProbability_constraintPy->Write();
	h_fitProbability_constraintPz->Write();
	h_fitProbability_constraintE->Write();
	h_fitProbability_sigmaPx->Write();
	h_fitProbability_sigmaPy->Write();
	h_fitProbability_sigmaPz->Write();
	h_fitProbability_sigmaE->Write();
	h_constraintPx_lowFitProb->Write();
	h_constraintPx_midFitProb->Write();
	h_constraintPx_highFitProb->Write();
	h_constraintPy_lowFitProb->Write();
	h_constraintPy_midFitProb->Write();
	h_constraintPy_highFitProb->Write();
	h_constraintPz_lowFitProb->Write();
	h_constraintPz_midFitProb->Write();
	h_constraintPz_highFitProb->Write();
	h_constraintE_lowFitProb->Write();
	h_constraintE_midFitProb->Write();
	h_constraintE_highFitProb->Write();
	h_constraintPx_uncertaintyPx->Write();
	h_constraintPy_uncertaintyPy->Write();
	h_constraintPz_uncertaintyPz->Write();
	h_constraintE_uncertaintyE->Write();
	h_constraintPx_uncertaintyPx_lowFitProb->Write();
	h_constraintPx_uncertaintyPx_midFitProb->Write();
	h_constraintPx_uncertaintyPx_highFitProb->Write();
	h_constraintPy_uncertaintyPy_lowFitProb->Write();
	h_constraintPy_uncertaintyPy_midFitProb->Write();
	h_constraintPy_uncertaintyPy_highFitProb->Write();
	h_constraintPz_uncertaintyPz_lowFitProb->Write();
	h_constraintPz_uncertaintyPz_midFitProb->Write();
	h_constraintPz_uncertaintyPz_highFitProb->Write();
	h_constraintE_uncertaintyE_lowFitProb->Write();
	h_constraintE_uncertaintyE_midFitProb->Write();
	h_constraintE_uncertaintyE_highFitProb->Write();
	h_pull_jet_theta_jet_phi_lowFitProb->Write();
	h_pull_jet_theta_jet_phi_midFitProb->Write();
	h_pull_jet_theta_jet_phi_highFitProb->Write();
	h_NormalizedResidualTheta_mcQmcJ_ThetaMcJet->Write();
	h_NormalizedResidualPhi_mcQmcJ_PhiMcJet->Write();
	h_NormalizedResidualEnergy_mcQmcJ_EnergyMcJet->Write();
	h_NormalizedResidualTheta_mcQmcJ_EnergyMcJet->Write();
	h_NormalizedResidualPhi_mcQmcJ_EnergyMcJet->Write();
	h_NormalizedResidualTheta_mcJrecJ_ThetaRecJet->Write();
	h_NormalizedResidualPhi_mcJrecJ_PhiRecJet->Write();
	h_NormalizedResidualEnergy_mcJrecJ_EnergyRecJet->Write();
	h_NormalizedResidualTheta_mcJrecJ_EnergyRecJet->Write();
	h_NormalizedResidualPhi_mcJrecJ_EnergyRecJet->Write();
	h_NormalizedResidualTheta_mcQrecJ_ThetaRecJet->Write();
	h_NormalizedResidualPhi_mcQrecJ_PhiRecJet->Write();
	h_NormalizedResidualEnergy_mcQrecJ_EnergyRecJet->Write();
	h_NormalizedResidualTheta_mcQrecJ_EnergyRecJet->Write();
	h_NormalizedResidualPhi_mcQrecJ_EnergyRecJet->Write();
	h_recE_fitE->Write();
	h_recPhifitPhi_recThetafitTheta->Write();
	h_recPhifitPhi_recThetafitTheta_lowfit->Write();
	h_recPhifitPhi_recThetafitTheta_midfit->Write();
	h_recPhifitPhi_recThetafitTheta_highfit->Write();
	h_recPhifitPhi_recThetafitTheta_SigmaThetaPhi->Write();
	h_recPhifitPhi_recThetafitTheta_lowfit_SigmaThetaPhi->Write();
	h_recPhifitPhi_recThetafitTheta_midfit_SigmaThetaPhi->Write();
	h_recPhifitPhi_recThetafitTheta_highfit_SigmaThetaPhi->Write();
	m_pTFile->Close();
	delete m_pTFile;

}
