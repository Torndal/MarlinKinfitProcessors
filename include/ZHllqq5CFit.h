#ifndef ZHllqq5CFit_h
#define ZHllqq5CFit_h 1

#include <iostream>
#include <vector>
#include <string>

#include "marlin/Processor.h"
#include "lcio.h"
#include "TrueJet_Parser.h"
#include <EVENT/Vertex.h>
#include <EVENT/ReconstructedParticle.h>
#include "TLorentzVector.h"
#include "DDMarlinCED.h"

#include <GeometryUtil.h>
#include <CLHEP/Vector/LorentzVector.h>
#include "JetFitObject.h"
#include "LeptonFitObject.h"
#include "ISRPhotonFitObject.h"
#include "MomentumConstraint.h"
#include "OPALFitterGSL.h"
#include "NewFitterGSL.h"
#include "TextTracer.h"
#include "NewtonFitterGSL.h"
#include "MassConstraint.h"
#include "SoftGaussParticleConstraint.h"
#include "SoftGaussMassConstraint.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/LCCollectionVec.h"
#include <EVENT/LCCollection.h>
#include "TFile.h"
#include "TH1F.h"
#include "TH2I.h"
#include "TH2F.h"
#include "TTree.h"

using namespace lcio ;
using namespace marlin ;

class ZHllqq5CFit : public Processor , public TrueJet_Parser
{

	public:

		virtual Processor*  newProcessor()
		{
			return new ZHllqq5CFit;
		}
		ZHllqq5CFit() ;
		virtual ~ZHllqq5CFit() = default;
		ZHllqq5CFit(const ZHllqq5CFit&) = delete;
		ZHllqq5CFit& operator=(const ZHllqq5CFit&) = delete;
		typedef std::vector<EVENT::ReconstructedParticle*>	pfoVector;
		virtual void	init();
		virtual void	Clear();
		virtual void	processRunHeader();
		virtual void	processEvent( EVENT::LCEvent *pLCEvent );
		virtual void	getSLDsInJet( LCRelationNavigator JetSLDNav , EVENT::ReconstructedParticle* jet , pfoVector &Neutrinos );
		virtual void	getSLDCombination( std::vector<int> jet1nSLDSolutions , int iteration , std::vector<int> &jetSLDCombination );
		int		performFIT( 	TLorentzVector jet1FourMomentum , std::vector<float> jet1CovMat , TLorentzVector jet2FourMomentum , std::vector<float> jet2CovMat , pfoVector leptons ,
						float &fitProbability , float (&fitOutputs)[ 18 ] , std::vector< TLorentzVector > &fittedObjects , float (&pull)[ 12 ] , bool traceEvent );
//		std::vector<std::vector<float>>	performOldFIT(EVENT::LCEvent *pLCEvent, TLorentzVector Jet0_Nutlv, TLorentzVector Jet1_Nutlv , std::vector<float> nu1CovMat , std::vector<float> nu2CovMat);
		virtual void	getJetResolutions(	TLorentzVector jetFourMomentum , std::vector<float> jetCovMat , float &sigmaE , float &sigmaTheta , float &sigmaPhi );
		virtual void	getLeptonParameters( ReconstructedParticle* lepton , float (&parameters)[ 3 ] , float (&errors)[ 3 ] );
		virtual void	getNormalizedResiduals( EVENT::LCEvent *pLCEvent , ReconstructedParticleImpl* outJet1 , ReconstructedParticleImpl* outJet2 );
		virtual void	check( LCEvent * evt );
		virtual void	end();
		std::string get_recoMCTruthLink()
		{
			return _recoMCTruthLink;
		};

	private:

		std::string				m_inputIsolatedlaptonCollection{};
		std::string				m_inputJetCollection{};
		std::string				m_inputSLDVertexCollection{};
		std::string				m_inputJetSLDLink{};
		std::string				m_outputFitcollection{};
		std::string				m_outputJetCollection{};
		std::string				m_outputFile{};
		int					m_nAskedJets{};
		int					m_nAskedIsoLeps{};
		bool					m_fitISR = true;
		int					m_fitter{};
		bool					m_traceall{};
		int					m_ievttrace{};

		int					m_nJets{};
		int					m_nIsoLeps{};
		int					m_nCorrectedSLD{};
		int					m_nRun;
		int					m_nEvt;
		int					m_nRunSum;
		int					m_nEvtSum;
		float					m_Bfield;
		float					c;
		float					mm2m;
		float					eV2GeV;
		float					eB;
		float					m_ECM{};
		float					m_isrpzmax{};
		double					b{};
		double					ISRPzMaxB{};
		TFile					*m_pTFile{};
	        TTree					*m_pTTree{};
		int					m_nSLDecayBHadron{};
		int					m_nSLDecayCHadron{};
		int					m_nSLDecayTauLepton{};
		int					m_nSLDecayTotal{};
		int					m_FitErrorCode_woNu{};
		float					m_ZMassBeforeFit_woNu{};
		float					m_HMassBeforeFit_woNu{};
		float					m_ZMassAfterFit_woNu{};
		float					m_HMassAfterFit_woNu{};
		std::vector<float>			m_pullJetEnergy_woNu{};
		std::vector<float>			m_pullJetTheta_woNu{};
		std::vector<float>			m_pullJetPhi_woNu{};
		std::vector<float>			m_pullLeptonInvPt_woNu{};
		std::vector<float>			m_pullLeptonTheta_woNu{};
		std::vector<float>			m_pullLeptonPhi_woNu{};
		std::vector<float>			m_normalizedResidualJetEnergy_woNu{};
		std::vector<float>			m_normalizedResidualJetTheta_woNu{};
		std::vector<float>			m_normalizedResidualJetPhi_woNu{};
		std::vector<float>			m_normalizedResidualLeptonInvPt_woNu{};
		std::vector<float>			m_normalizedResidualLeptonTheta_woNu{};
		std::vector<float>			m_normalizedResidualLeptonPhi_woNu{};
		float					m_FitProbability_woNu{};
		int					m_FitErrorCode_wNu{};
		float					m_ZMassBeforeFit_wNu{};
		float					m_HMassBeforeFit_wNu{};
		float					m_ZMassAfterFit_wNu{};
		float					m_HMassAfterFit_wNu{};
		float					m_FitProbability_wNu{};
		std::vector<float>			m_pullJetEnergy_wNu{};
		std::vector<float>			m_pullJetTheta_wNu{};
		std::vector<float>			m_pullJetPhi_wNu{};
		std::vector<float>			m_pullLeptonInvPt_wNu{};
		std::vector<float>			m_pullLeptonTheta_wNu{};
		std::vector<float>			m_pullLeptonPhi_wNu{};
		std::vector<float>			m_normalizedResidualJetEnergy_wNu{};
		std::vector<float>			m_normalizedResidualJetTheta_wNu{};
		std::vector<float>			m_normalizedResidualJetPhi_wNu{};
		std::vector<float>			m_normalizedResidualLeptonInvPt_wNu{};
		std::vector<float>			m_normalizedResidualLeptonTheta_wNu{};
		std::vector<float>			m_normalizedResidualLeptonPhi_wNu{};
		int					m_FitErrorCode{};
		float					m_ZMassBeforeFit{};
		float					m_HMassBeforeFit{};
		float					m_ZMassAfterFit{};
		float					m_HMassAfterFit{};
		float					m_FitProbability{};
		std::vector<float>			m_pullJetEnergy{};
		std::vector<float>			m_pullJetTheta{};
		std::vector<float>			m_pullJetPhi{};
		std::vector<float>			m_pullLeptonInvPt{};
		std::vector<float>			m_pullLeptonTheta{};
		std::vector<float>			m_pullLeptonPhi{};
		std::vector<float>			m_normalizedResidualJetEnergy{};
		std::vector<float>			m_normalizedResidualJetTheta{};
		std::vector<float>			m_normalizedResidualJetPhi{};
		std::vector<float>			m_normalizedResidualLeptonInvPt{};
		std::vector<float>			m_normalizedResidualLeptonTheta{};
		std::vector<float>			m_normalizedResidualLeptonPhi{};

		double					ZEnergy{};
		double					Zmomentum[3]{0.0};
		double					HEnergy{};
		double					Hmomentum[3]{0.0};
		double					ISREnergy{};
		double					ISRmomentum[3]{0.0};
		float					Hmass_NoFit{};
		int					Error_code{};
		int					errorcode{};
		float					hpull_jet_E{};
		float					hpull_jet2_E{};
		float					hpull_jet_th{};
		float					hpull_jet2_th{};
		float					hpull_jet_phi{};
		float					hpull_jet2_phi{};
		float					hpull_lepton_InvpT{};
		float					hpull_lepton2_InvpT{};
		float					hpull_lepton_th{};
		float					hpull_lepton2_th{};
		float					hpull_lepton_phi{};
		float					hpull_lepton2_phi{};
};

#endif
