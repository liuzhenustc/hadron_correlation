#ifndef STMINITREEMAKER_HH
#define STMINITREEMAKER_HH

/***************************************************************************
 *
 * $Id: StMiniTreeMaker.h 2015/04/09  Exp $ 
 * StMiniTreeMaker - class to produce miniTree for mtd related analysis
 * Author: Shuai Yang
 *--------------------------------------------------------------------------
 *
 ***************************************************************************/

#include "StMaker.h"
#include "StTreeStructure.h"

#include "StThreeVectorF.hh"
#include "TLorentzVector.h"

#include "StMtdUtil/StMtdConstants.h"

#include <vector>
#include <map>
#ifndef ST_NO_NAMESPACES
using std::vector;
#endif

class TH1D;
class TH2D;
class TString;
class TTree;
class TFile;

class StMuDstMaker;
class StMuDst;
class StMuTrack;
class StMuMtdHit;
class StRefMultCorrVPDMBZDCNoVtx;

#if !defined(ST_NO_TEMPLATE_DEF_ARGS) || defined(__CINT__)
typedef vector<Bool_t> BoolVec;
typedef vector<Int_t> IntVec;
typedef vector<Double_t> DoubleVec;
typedef vector<TLorentzVector> LorentzVec;
typedef vector<StThreeVectorF> StThreeVecF;
#else
typedef vector<Bool_t, allocator<Bool_t>> BoolVec;
typedef vector<Int_t, allocator<Int_t>> IntVec;
typedef vector<Double_t, allocator<Double_t>> DoubleVec;
typedef vector<TLorentzVector, allocator<TLorentzVector>> LorentzVec;
typedef vector<StThreeVectorF, allocator<StThreeVectorF>> StThreeVecF;
#endif

const Float_t PxCutL[9]={-2.5,-4.5,-7.5,-11.5,-16.5,-20.5,-27.5,-32.5,-37.5};
const Float_t PxCutR[9]={-2.5,-4.5,-7.5,-11.5,-17.5,-21.5,-28.5,-33.5,-40.5};

class StMiniTreeMaker : public StMaker {
    public:
        StMiniTreeMaker(const Char_t *name = "StMiniTreeMaker");
        ~StMiniTreeMaker();

        Int_t    Init();
        Int_t    InitRun();
        Int_t    Make();
        Int_t    Finish();

        void     setMaxTpcVz(const Double_t max);
        void     setMaxDiffVz(const Double_t max);
        void     setMinTrackP(const Double_t min);
        void     setMaxTrackEta(const Double_t max);
        void     setMinNHitsFit(const Int_t min);
        void     setMinNHitsFitRatio(const Int_t min);
        void     setMinNHitsDedx(const Int_t min);
        void     setMaxnSigmaPi(const Double_t max);//
        void     setFillHisto(const Bool_t fill);
        void     setFillTree(const Bool_t fill);
        void     setOutFileName(const TString name);
        void     setPrintMemory(const Bool_t pMem);
        void     setPrintCpu(const Bool_t pCpu);
        void     setPrintConfig(const Bool_t print);

    protected:
        void     printConfig();
        void     bookTree();
        void     bookHistos();
        Bool_t   processEvent();
        Bool_t   IsValidTrack(StMuTrack *track);
        Bool_t   IsMtdTrack(StMuTrack *track);
        Bool_t   IsPion(StMuTrack *track);
        Bool_t   IsPionEvent();
        Bool_t   Is2PionEvent();
        Bool_t   PassPxCandidate(StMuTrack *track);
        Bool_t   PassPiCandidate(StMuTrack *track);
        void     ComparePt(StThreeVectorF* vec1 ,StThreeVectorF* vec2);
        void     GetTrigger();
        Float_t  PxCal(StThreeVecF trkVec);
        Bool_t   PassPxCut(Float_t PxL,Float_t PxR, Short_t cen);

    private:

        static Double_t constexpr mPiMass    = 0.13957018;
        static Double_t constexpr mPMass = 0.938272046;
        StMuDstMaker   *mMuDstMaker;          // Pointer to StMuDstMaker
        StMuDst        *mMuDst;              // Pointer to MuDst event

        StRefMultCorrVPDMBZDCNoVtx *refMultCorr; //decide centrality

        Bool_t         mFillHisto;           // Flag of fill the histogram 
        Bool_t         mFillTree;            // Flag of fill the event tree

        Bool_t         mPrintMemory;         // Flag to print out memory usage
        Bool_t         mPrintCpu;            // Flag to print out CPU usage
        Bool_t         mPrintConfig;         // Flag to print out task configuration
        Double_t       mMaxTpcVz;             // Maximum r
        Double_t       mMaxDiffVz;           // Maximum VpdVz-TpcVz 
        Double_t       mMinTrkP;            // Minimum track pt
        Double_t       mMaxTrkP;            // Minimum track pt
        Double_t       mMinTrkEta;           // Maximum track eta
        Double_t       mMaxTrkEta;           // Maximum track eta
        Double_t       mMinTrkPhi;           // Maximum track eta
        Double_t       mMaxTrkPhi;           // Maximum track eta
        Int_t          mMinNHitsFit;         // Minimum number of hits used for track fit
        Double_t       mMinNHitsFitRatio;    // Minimum ratio of hits used for track fit
        Int_t          mMinNHitsDedx;        // Minimum number of hits used for de/dx
        Double_t       mMaxDca;              // Maximum track dca
        Double_t       mMinNsigmaPi;         // Maximum nSigmaPi cut
        Double_t       mMaxNsigmaPi;         // Maximum nSigmaPi cut

        TFile          *fOutFile;            // Output file
        TString        mOutFileName;         // Name of the output file 
        StEvtData      mEvtData;
        TTree          *mEvtTree;            // Pointer to the event tree

        Int_t          mModuleToQT[gMtdNBacklegs][gMtdNModules];     // Map from module to QT board index
        Int_t          mModuleToQTPos[gMtdNBacklegs][gMtdNModules];  // Map from module to the position on QA board
        Int_t          mQTtoModule[4][8];                            // Map from QT board to module index
        Int_t          mQTSlewBinEdge[4][16][8];                     // Bin edges for online slewing correction for QT
        Int_t          mQTSlewCorr[4][16][8];                        // Slewing correction values for QT
        Int_t          mTrigQTpos[4][2];                             // corresponding QT position of fire trigger bit

        IntVec      mTriggerIDs;
        StThreeVecF PiCandidate;
        StThreeVecF TrkCandidateL;//-1<eta<0.5
        StThreeVecF TrkCandidateR;//0.5<eta<1

        //define histograms ongoing...
        TH1D           *hEvent;
        TH1D           *hDca;
      
        ClassDef(StMiniTreeMaker, 1)
};

inline void StMiniTreeMaker::setMaxTpcVz(const Double_t max) { mMaxTpcVz = max; }
inline void StMiniTreeMaker::setMaxDiffVz(const Double_t max) { mMaxDiffVz = max; }
inline void StMiniTreeMaker::setMinTrackP(const Double_t min){ mMinTrkP = min;}
inline void StMiniTreeMaker::setMaxTrackEta(const Double_t max){ mMaxTrkEta = max; }
inline void StMiniTreeMaker::setMinNHitsFit(const Int_t min) { mMinNHitsFit = min; }
inline void StMiniTreeMaker::setMinNHitsFitRatio(const Int_t min) { mMinNHitsFitRatio = min; }
inline void StMiniTreeMaker::setMinNHitsDedx(const Int_t min) { mMinNHitsDedx = min; }
inline void StMiniTreeMaker::setMaxnSigmaPi(const Double_t max) { mMaxNsigmaPi = max; }
inline void StMiniTreeMaker::setFillHisto(const Bool_t fill) { mFillHisto = fill; }
inline void StMiniTreeMaker::setFillTree(const Bool_t fill) { mFillTree = fill; }
inline void StMiniTreeMaker::setOutFileName(const TString name) { mOutFileName = name; }
inline void StMiniTreeMaker::setPrintMemory(const Bool_t pMem) { mPrintMemory = pMem; }
inline void StMiniTreeMaker::setPrintCpu(const Bool_t pCpu) { mPrintCpu = pCpu; }
inline void StMiniTreeMaker::setPrintConfig(const Bool_t print) { mPrintConfig = print; }
#endif
