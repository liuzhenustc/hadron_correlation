#include "headers.h"
#include "StMiniTreeMaker.h"

ClassImp(StMiniTreeMaker)

    //_____________________________________________________________________________
    StMiniTreeMaker::StMiniTreeMaker(const Char_t *name) :
        StMaker(name), 
        mFillTree(1), mFillHisto(1), 
        mPrintConfig(1), mPrintMemory(0), mPrintCpu(0),
        mMaxTpcVz(100.), mMaxDiffVz(3.), 
        mMinTrkP(0.15), mMaxTrkP(1e4),
        mMaxTrkEta(1.), mMinTrkEta(-1.),
        mMinTrkPhi(0.), mMaxTrkPhi(2*3.1415926),
        mMinNHitsFit(15), mMinNHitsFitRatio(0.52), mMinNHitsDedx(10.),
        mMaxDca(3.), 
        mMinNsigmaPi(0.), mMaxNsigmaPi(3.), 
        fOutFile(0), mOutFileName(""), mEvtTree(0)
{
    mTriggerIDs.clear();
    mTriggerIDs.push_back(450601); // dimuon
    mTriggerIDs.push_back(450611); // dimuon
    mTriggerIDs.push_back(450621); // dimuon
    mTriggerIDs.push_back(450631); // dimuon
    mTriggerIDs.push_back(450641); // dimuon
    mTriggerIDs.push_back(450604); // dimuon-30-hft (production_2014)
    mTriggerIDs.push_back(450605); // dimuon-5-hft (production_mid_2014; production_low_2014)
    mTriggerIDs.push_back(450606); // dimuon-5-hft (production_mid_2014)
    mTriggerIDs.push_back(450602); // e-mu (production_2014)
    mTriggerIDs.push_back(450612); // e-mu (production_2014; production_low_2014)
    mTriggerIDs.push_back(450622); // e-mu (production_2014; production_low_2014)
    mTriggerIDs.push_back(450632); // e-mu (production_mid_2014)
    mTriggerIDs.push_back(450642); // e-mu (production_2014; production_low_2014)
    mTriggerIDs.push_back(450600); // single-muon (production_2014)
    mTriggerIDs.push_back(450610); // single-muon (production_2014; production_low_2014)
    mTriggerIDs.push_back(450620); // single-muon (production_2014; production_low_2014)
    mTriggerIDs.push_back(450630); // single-muon (production_mid_2014)
    mTriggerIDs.push_back(450640); // single-muon (production_2014, production_low_2014)
}
//_____________________________________________________________________________
StMiniTreeMaker::~StMiniTreeMaker()
{
    // default destructor
}
//_____________________________________________________________________________
Int_t StMiniTreeMaker::Init()
{
    refMultCorr = new StRefMultCorrVPDMBZDCNoVtx("grefmult");
    if(!mOutFileName.Length())
    {
        LOG_ERROR << "StMiniTreeMaker:: no output file specified for tree and histograms." << endm;
        return kStERR;
    }
    fOutFile = new TFile(mOutFileName.Data(),"recreate");
    LOG_INFO << "StMiniTreeMaker:: create the output file to store the tree and histograms: " << mOutFileName.Data() << endm;

    if(mFillTree)    bookTree();
    if(mFillHisto)   bookHistos();

    return kStOK;
}
//_____________________________________________________________________________
Int_t StMiniTreeMaker::InitRun()
{
    //initial maps
    memset(mModuleToQT,-1,sizeof(mModuleToQT));
    memset(mModuleToQTPos,-1,sizeof(mModuleToQTPos));
    memset(mQTtoModule,-1,sizeof(mQTtoModule));

    // obtain maps from DB
    LOG_INFO << "Retrieving mtdModuleToQTmap table from database ..." << endm;
    TDataSet *dataset = GetDataBase("Geometry/mtd/mtdModuleToQTmap");
    //cout<<"output dataset pointer:"<<dataset<<endl;
    St_mtdModuleToQTmap *mtdModuleToQTmap = static_cast<St_mtdModuleToQTmap*>(dataset->Find("mtdModuleToQTmap"));
    //cout<<"output mtdModuleToQTmap pointer:"<<mtdModuleToQTmap<<endl;//usedfor debug 
    if(!mtdModuleToQTmap){
        LOG_ERROR << "No mtdModuleToQTmap table found in database" << endm;
        return kStErr;
    }
    mtdModuleToQTmap_st *mtdModuleToQTtable = static_cast<mtdModuleToQTmap_st*>(mtdModuleToQTmap->GetTable());


    for(Int_t i=0; i<gMtdNBacklegs; i++)
    {
        for(Int_t j=0; j<gMtdNModules; j++)
        {
            Int_t index = i*5 + j;
            Int_t qt = mtdModuleToQTtable->qtBoardId[index];
            Int_t channel = mtdModuleToQTtable->qtChannelId[index];
            mModuleToQT[i][j] = qt;
            if(channel<0)
            {
                mModuleToQTPos[i][j] = channel;
            }
            else
            {
                if(channel%8==1) mModuleToQTPos[i][j] = 1 + channel/8 * 2;
                else             mModuleToQTPos[i][j] = 2 + channel/8 * 2;
            }
            if(mModuleToQT[i][j]>0 && mModuleToQTPos[i][j]>0)
                mQTtoModule[mModuleToQT[i][j]-1][mModuleToQTPos[i][j]-1] = j + 1;
        }
    }

    if(Debug()){
        for(Int_t i=0;i<gMtdNBacklegs;i++){
            for(Int_t j=0;j<gMtdNModules;j++){
                LOG_INFO<<"BL:"<<i+1<<"\t Module:"<<j+1<<"\t QT:"<<mModuleToQT[i][j]<<"\t Pos:"<<mModuleToQTPos[i][j]<<endm;
            }
        }
    }

    // online slewing correction for QT board
    memset(mQTSlewBinEdge,-1,sizeof(mQTSlewBinEdge));
    memset(mQTSlewCorr,-1,sizeof(mQTSlewCorr));
    LOG_INFO << "Retrieving mtdQTSlewingCorr table from database ..." << endm;
    dataset = GetDataBase("Calibrations/mtd/mtdQTSlewingCorr");
    St_mtdQTSlewingCorr *mtdQTSlewingCorr = static_cast<St_mtdQTSlewingCorr*>(dataset->Find("mtdQTSlewingCorr"));
    if(!mtdQTSlewingCorr){
        LOG_ERROR << "No mtdQTSlewingCorr table found in database" << endm;
        return kStErr;
    }
    mtdQTSlewingCorr_st *mtdQTSlewingCorrtable = static_cast<mtdQTSlewingCorr_st*>(mtdQTSlewingCorr->GetTable());
    for(Int_t j=0; j<4; j++){
        for(Int_t i=0; i<16; i++){
            for(Int_t k=0; k<8; k++){
                Int_t index = j*16*8 + i*8 + k;
                mQTSlewBinEdge[j][i][k] = (Int_t) mtdQTSlewingCorrtable->slewingBinEdge[index];
                mQTSlewCorr[j][i][k] = (Int_t) mtdQTSlewingCorrtable->slewingCorr[index];
            }
        }
    }

    return kStOK;
}
//_____________________________________________________________________________
Int_t StMiniTreeMaker::Finish()
{
    if(fOutFile)
    {
        fOutFile->cd();
        fOutFile->Write();
        fOutFile->Close();
        LOG_INFO << "StMiniTreeMaker::Finish() -> write out tree in " << mOutFileName.Data() << endm;
    }

    if(mPrintConfig) printConfig();

    return kStOK;
}
//_____________________________________________________________________________
Int_t StMiniTreeMaker::Make()
{
    memset(&mEvtData, 0, sizeof(mEvtData)); //initial the mEvtData structure

    StTimer timer;
    if(mPrintMemory) StMemoryInfo::instance()->snapshot();
    if(mPrintCpu)    timer.start();

    mMuDstMaker = (StMuDstMaker *)GetMaker("MuDst");
    if(!mMuDstMaker){
        LOG_WARN<<"No StMuDstMaker !"<<endm;
        return kStOK;
    }

    mMuDst = mMuDstMaker->muDst();
    if(!mMuDst){
        LOG_WARN<<"No MuDst !"<<endm;
        return kStOK;
    }

    mMuEvent = mMuDst->event();
    if(!mMuEvent){
        LOG_WARN<<"No Muevent !"<<endm;
        return kFALSE;
    }

    PiCandidate.clear();
    TrkCandidateL.clear();
    TrkCandidateR.clear();

    if(!passEvent(mMuEvent)) return kStOK;
    if(!processEvent()) return kStOK;
    if(mFillTree) mEvtTree->Fill();
    if(mPrintMemory){
        StMemoryInfo::instance()->snapshot();
        StMemoryInfo::instance()->print();
    }
    if(mPrintCpu){
        timer.stop();
        LOG_INFO << "CPU time for StMiniTreeMaker::Make(): " 
            << timer.elapsedTime() << "sec " << endm;
    }

    return kStOK;
}
//_____________________________________________________________________________
Bool_t StMiniTreeMaker::passEvent(StMuEvent *ev)
{
    //if(Debug()) LOG_INFO<<"Already in, Great!"<<endm;
    //if(Debug()) LOG_INFO<<"nTriggerMtd: "<<nTriggerMtd<<endm;
    //if(Debug()) LOG_INFO<<"mTriggerIdMtd0: "<<Pico::mTriggerIdMtd[0]<<endm;
    if(mFillHisto) hEvent->Fill(0.5);
    StBTofHeader *mBTofHeader = mMuDst->btofHeader();
    //////////////////////////////////////
    //select the right vertex using VPD
    //////////////////////////////////////
    Float_t vzVpd = -999;
    if(mBTofHeader) vzVpd = mBTofHeader->vpdVz();
    for(unsigned int i=0;i<mMuDst->numberOfPrimaryVertices();i++) {
        StMuPrimaryVertex *vtx = mMuDst->primaryVertex(i);
        if(!vtx) continue;
        Float_t vz = vtx->position().z();
        if(fabs(vzVpd)<200 && fabs(vzVpd-vz)<3.) {
            mMuDst->setVertexIndex(i);
            break;
        }
    }

    if(!ev){
        LOG_INFO << "StPicoCut::passEvent  No StMuEvent" << endm;
        return kFALSE;
    }
    StThreeVectorF pVertex = ev->primaryVertexPosition();
    if(Debug()) LOG_INFO<<"vertex: "<<pVertex.x()<<endm;
    if(fabs(pVertex.x())<1.e-5 && fabs(pVertex.y())<1.e-5 && fabs(pVertex.z())<1.e-5){
        LOG_INFO << "StPicoCut::passEvent  bad vertices (x,y,z) = ("
            << pVertex.x() << ","
            << pVertex.y() << ","
            << pVertex.z() << ")"
            << endm;
        return kFALSE;
    }
    if(fabs(pVertex.z())>1e4){
        LOG_INFO << "StPicoCut::passEvent  z-vertex out of range, vz  (evtSum)=" << pVertex.z() << " (direct)=" << ev->primaryVertexPosition().z() << endm;
        return kFALSE;
    }
    const Float_t vx = pVertex.x() ;
    const Float_t vy = pVertex.y() ;
    if(sqrt(vx*vx+vy*vy)>1e4){
        LOG_INFO << "StPicoCut::passEvent  vr-vertex out of range, vr = " << sqrt(vx*vx+vy*vy)
            << ",  vx = " << vx
            << ",  vy = " << vy
            << endm;
        return kFALSE ;
    }

    bool isTrg = kFALSE;
    for(int i=0;i<nTrigger;i++) {
        if(ev->triggerIdCollection().nominal().isTrigger(Pico::mTriggerId[i])){
            isTrg = kTRUE;
            break;
        }
    }

    if(!isTrg){
        if(Debug()) LOG_INFO << "nTriggerMtd: " << nTriggerMtd <<endm;
        for(int i=0;i<nTriggerMtd;i++) {
            if(ev->triggerIdCollection().nominal().isTrigger(Pico::mTriggerIdMtd[i])){
                isTrg = kTRUE;
                break;
            }
        }
    }

    if(!isTrg) {
        LOG_INFO << "StPicoCut::passEvent trigger not fired " << endm;
        return kFALSE;
    }

    if(ev->refMult()<0) {
        LOG_INFO << "StPicoCut::passEvent refMult out of range, refMult = " << ev->refMult() << endm;
        return kFALSE;
    }

    if(mFillHisto) hEvent->Fill(1.5);

    return kTRUE;
}
//_____________________________________________________________________________
Bool_t StMiniTreeMaker::processEvent()
{
    //-----test
    if(Debug()) LOG_INFO << "test eventId: " << endm;
    Int_t testEventId = mMuEvent->eventNumber();
    if( testEventId != 1893484 ) return kFALSE;//test events;
    if(Debug()) LOG_INFO<<"test eventId: "<<testEventId<<endm;
    //-----test

    if(mFillHisto) hEvent->Fill(2.5);

    Bool_t validTrigger = kFALSE;
    Bool_t DiMuon       = kFALSE;
    Bool_t DiMuonHFT    = kFALSE;
    Bool_t SingleMuon   = kFALSE;
    Bool_t EMuon        = kFALSE;

    UInt_t triggerId    = 0;
    for(Int_t i=0;i<mTriggerIDs.size();i++){
        if(mMuEvent->triggerIdCollection().nominal().isTrigger(mTriggerIDs[i])){
            triggerId |= (1<<i);//didn't use right now
            if(i<5)       { DiMuon = kTRUE; }
            else if(i<8)  { DiMuonHFT = kTRUE; }
            else if(i<13) { EMuon = kTRUE; }
            else          { SingleMuon = kTRUE; }
            validTrigger = kTRUE;
        }
    }
    if(!validTrigger){
        LOG_WARN<<"No valid mtd related triggers !"<<endm;
        return kFALSE;
    }

    //----filter--Event---
    //StPicoDst* mPicoDst = new StPicoDst();
    //int nMthTrigHits = 0;
    //Int_t nMtdHits = mPicoDst->numberOfMtdHits();
    //Int_t nPidTraits = mPicoDst->numberOfMtdPidTraits();
    //for(Int_t i=0; i<nMtdHits; i++)
    //{
    //    StPicoMtdHit *hit = mPicoDst->mtdHit(i);
    //    if(!hit) continue;
    //    if(hit->triggerFlag()<1) continue;
    //    Int_t index = -1;
    //    for(Int_t j=0; j<nPidTraits; j++)
    //    {
    //        StPicoMtdPidTraits *mtdPid = mPicoDst->mtdPidTraits(j);
    //        if(mtdPid->backleg()==hit->backleg() &&
    //                mtdPid->module()==hit->module() &&
    //                mtdPid->cell()==hit->cell())
    //        {
    //            index = j;
    //            break;
    //        }
    //    }
    //    if(index<0) continue;
    //    nMthTrigHits ++;
    //}
    //cout << "filter Event: " << nMthTrigHits <<endl;
    //----------------------

    //reject dimuon event overlap with other triggers
    //but should have been rejected
    //if( DiMuon ){ 
    //    StMuMtdHeader *muMtdHeader = mMuDst->mtdHeader();
    //    if(muMtdHeader && muMtdHeader->shouldHaveRejectEvent()==1) DiMuon = kFALSE;
    //}
    //if( DiMuonHFT ){ 
    //    StMuMtdHeader *muMtdHeader = mMuDst->mtdHeader();
    //    if(muMtdHeader && muMtdHeader->shouldHaveRejectEvent()==1) DiMuonHFT = kFALSE;
    //}

    if(!DiMuon && !DiMuonHFT)  return kFALSE;//dimuon and dimuon hft event
    if(mFillHisto)  hEvent->Fill(3.5);
    if(Debug()) LOG_INFO<<"Pass Dimuon Trigger"<<endm;

    //select the right vertex using VPD
    Float_t vpdVz = -999; 
    StBTofHeader *mBTofHeader = mMuDst->btofHeader();
    if(mBTofHeader) vpdVz = mBTofHeader->vpdVz();//select vpdVz
    if(abs(vpdVz)>200)
    {
        LOG_DEBUG << "Bad vpd vertext: " << vpdVz << endm;
        return kStOK;
    }

    StMuPrimaryVertex* priVertex = mMuDst->primaryVertex();
    StThreeVectorF verPos = priVertex->position();
    Float_t tpcVz = verPos.z();
    if(Debug()) LOG_INFO << "primary vertex: " << tpcVz << endm;

    //---vertex--cut
    if(mMaxTpcVz<1e4 && abs(tpcVz)>mMaxTpcVz) return kStOK;
    if(mMaxDiffVz<1e4 && abs(tpcVz-vpdVz)>mMaxDiffVz) return kStOK;
    if(mFillHisto)  hEvent->Fill(4.5);

    //---pion--match-with-MTD 
    //if(!IsPionEvent()) return kStOK;
    //if(Debug()) LOG_INFO<<"Pass pion Trigger"<<endm;
    //if(!Is2PionEvent()) return kStOK;
    //if(Debug()) LOG_INFO<<"Pass 2pion Trigger"<<endm;
    //if(mFillHisto)  hEvent->Fill(5.5);

    //---get--centrality---
    Int_t runId = mMuEvent->runNumber();
    if(Debug()) cout<<"runId: "<<runId<<endl;
    Int_t eventId = mMuEvent->eventNumber();
    if(Debug()) cout<<"eventId: "<<eventId<<endl;
    StRunInfo runInfo = mMuEvent->runInfo();
    float zdcRate = runInfo.zdcCoincidenceRate();
    if(Debug()) cout<<"zdcRate: "<<zdcRate<<endl;
    Int_t refMult = mMuEvent->refMult();
    if(Debug()) cout<<"refMult: "<<refMult<<endl;
    Int_t gRefMult = mMuEvent->grefmult();
    if(Debug()) cout<<"gRefMult: "<<gRefMult<<endl;

    refMultCorr->init(runId);
    refMultCorr->initEvent(gRefMult, tpcVz, zdcRate);
    Float_t gRefMultCorr   = refMultCorr->getRefMultCorrVPDMBZDCNoVtx();
    Float_t evtWeight      = refMultCorr->getWeight();
    Short_t centrality     = refMultCorr->getCentralityBin16();
    if(Debug()) LOG_INFO<<"gRefMult: "<<gRefMult<<" \tgRefMultCorr: "<<gRefMultCorr<<" \tmCentrality: "<<centrality<<endm;

    //--for----track--check-------
    int nPrimary = mMuDst->numberOfPrimaryTracks(); 
    if(Debug()) cout<<"# of primary track: "<<nPrimary<<endl;
    int nTrack = 0;
    for(Int_t i=0;i<nPrimary;i++){
        StMuTrack* pTrack = mMuDst->primaryTracks(i);
        Int_t trkId = pTrack->id();
        if(trkId == 406 ){ 
            if(Debug()) cout<<"trkId: "<<trkId<<endl;
            if(Debug()) cout<<"nSigmaPi: "<<pTrack->nSigmaPion()<<endl;
            cout<<"trkFlag: "<<pTrack->globalTrack()->flag()<<endl;
            StThreeVectorF P = pTrack->momentum();
            if(Debug()) cout<<"P: "<< P.perp()<<endl;
            if(Debug()) cout<<"eta: "<< P.pseudoRapidity()<<endl;
            if(Debug()) cout<<"Phi: "<< P.phi()<<endl;
            cout<<"IsValidTrack: "<<IsValidTrack(pTrack)<<endl;
            cout<<"nHitsFit: "<<pTrack->nHitsFit(kTpcId)<<endl;
            cout<<"nHitsFit_global: "<<pTrack->globalTrack()->nHitsFit(kTpcId)<<endl;
            cout<<"kTpcId: "<<kTpcId<<endl;
            cout<<"nHitsDedx: "<<pTrack->nHitsDedx()<<endl;
            cout<<"nHitsRatio: "<<pTrack->nHitsFit(kTpcId)/(1.0*pTrack->nHitsPoss(kTpcId))<<endl;
            cout<<"IsMtdTrack: "<<IsMtdTrack(pTrack)<<endl;
            cout<<"iMtd: "<<pTrack->index2MtdHit()<<endl;
            cout<<"iMtd_global: "<<pTrack->globalTrack()->index2MtdHit()<<endl;
            cout<<"IsPion: "<<IsPion(pTrack)<<endl;
        }
        if(!PassPiCandidate(pTrack)) continue;
        mEvtData.mPiTrkId[nTrack] = pTrack->id();
        mEvtData.mnSigmaPi[nTrack] = pTrack->nSigmaPion();
        mEvtData.mgDca[nTrack] = pTrack->dcaGlobal().mag();
        nTrack++;
    }
    mEvtData.mNPiTrk  = nTrack;
    //for(Int_t i=0;i<nTrack;i++){
    //    mEvtData.mPiPt[i] = PiCandidate[i].perp();
    //    mEvtData.mPiEta[i] = PiCandidate[i].pseudoRapidity();
    //    mEvtData.mPiPhi[i] = PiCandidate[i].phi();
    //}
    //-----for---track---check----

    //--PxCut-centrality-dependent---
    //int nPrimary = mMuDst->numberOfPrimaryTracks();//Use primary tracks for Px candidate 
    //for(Int_t i=0;i<nPrimary;i++){
    //    StMuTrack* pTrack = mMuDst->primaryTracks(i);
    //    PassPxCandidate(pTrack);//--for TPC track
    //    PassPiCandidate(pTrack);//--for TPC track
    //}

    //GetTrigger();

    Float_t PxL = PxCal(TrkCandidateL);
    Float_t PxR = PxCal(TrkCandidateR);

    //if( !PassPxCut(PxL,PxR,centrality) )  return kFALSE;
    //if(mFillHisto)  hEvent->Fill(4.5);
    //if(Debug()) LOG_INFO<<"Pass Px Cut"<<endm;

    mEvtData.mRunId           =runId;
    mEvtData.mEventId         =eventId;
    mEvtData.mGRefMult        =gRefMult;
    mEvtData.mGRefMultCorr    =gRefMultCorr;
    mEvtData.mEvtWeight       =evtWeight;
    mEvtData.mCentrality      =centrality;
    mEvtData.mZDCRate         =zdcRate;
    mEvtData.mVpdVz           =vpdVz;
    mEvtData.mTpcVz           =tpcVz;

    mEvtData.mPxL   =PxL;
    mEvtData.mPxR   =PxR;

    //---track---level----
    //Int_t nPi = PiCandidate.size();
    //mEvtData.mNPiTrk  = nPi;
    //for(Int_t i=0;i<nPi;i++){
    //    mEvtData.mPiPt[i] = PiCandidate[i].perp();
    //    mEvtData.mPiEta[i] = PiCandidate[i].pseudoRapidity();
    //    mEvtData.mPiPhi[i] = PiCandidate[i].phi();
    //}

    return kTRUE;
}
//_________________________________________________________
Bool_t StMiniTreeMaker::IsValidTrack(StMuTrack *track)
{
    if(!track) return kFALSE;
    StThreeVectorF mom = track->momentum();
    Float_t p = mom.mag();
    Float_t eta = mom.pseudoRapidity();
    Float_t phi = mom.phi();
    if(phi<0) phi += 2*3.1415926;

    if(p < mMinTrkP   || p > mMaxTrkP)             return kFALSE;
    if(eta < mMinTrkEta || eta > mMaxTrkEta)           return kFALSE;
    if(phi < mMinTrkPhi || phi > mMaxTrkPhi)           return kFALSE;
    if(track->nHitsFit(kTpcId)<mMinNHitsFit)           return kFALSE;
    if(track->nHitsDedx()<mMinNHitsDedx)               return kFALSE;
    if(track->nHitsFit(kTpcId)/(1.0*track->nHitsPoss(kTpcId))<mMinNHitsFitRatio) return kFALSE;
    return kTRUE;
}
//_____________________________________________________________________________
Bool_t StMiniTreeMaker::IsMtdTrack(StMuTrack *track)
{
    Double_t gDca = track->dcaGlobal().mag();
    if( gDca > 3. ) {cout<<"IsMtdTrack_gDca: "<<gDca<<endl;return kFALSE;}
    int iMtd = track->index2MtdHit();
    if(iMtd<0) return kFALSE;
    return kTRUE;
}
//_____________________________________________________________________________
Bool_t StMiniTreeMaker::IsPion(StMuTrack *track)
{
    if(Debug()) cout<<"in IsPion "<<endl;
    Float_t nSgmPi = track->nSigmaPion();
    if(nSgmPi < mMinNsigmaPi || nSgmPi > mMaxNsigmaPi) return kFALSE;
    if(Debug())cout<<"IsPion: "<<nSgmPi<<endl;
    return kTRUE;
}
//_____________________________________________________________________________
Bool_t StMiniTreeMaker::IsPionEvent()
{
    Int_t nPi =0;
    int nPrimary = mMuDst->numberOfPrimaryTracks(); 
    for(Int_t i=0;i<nPrimary;i++){
        StMuTrack* pTrack = mMuDst->primaryTracks(i);
        if(!IsValidTrack(pTrack)) continue;
        if(!IsMtdTrack(pTrack)) continue;
        if(!IsPion(pTrack)) continue;
        nPi++;
    }
    if(nPi>0){
        if(Debug()) cout<<"IsPionEvent"<<endl;
        return kTRUE;
    } 
    else return kFALSE;
}
//_____________________________________________________________________________
Bool_t StMiniTreeMaker::Is2PionEvent()
{
    Int_t nPi =0;
    int nPrimary = mMuDst->numberOfPrimaryTracks(); 
    for(Int_t i=0;i<nPrimary;i++){
        StMuTrack* pTrack = mMuDst->primaryTracks(i);
        if(!IsValidTrack(pTrack)) continue;
        if(!IsMtdTrack(pTrack)) continue;
        if(!IsPion(pTrack)) continue;
        nPi++;
    }
    if(nPi>1){
        if(Debug()) cout<<"Is2PionEvent nPrimary: "<<endl;
        return kTRUE;
    }
    else return kFALSE;
}
//_____________________________________________________________________________
Bool_t StMiniTreeMaker::PassPxCandidate(StMuTrack *track)
{
    if(!IsValidTrack(track)) return kFALSE;
    StThreeVectorF mom = track->momentum();
    float eta = mom.pseudoRapidity();
    if(-1.<eta && eta<-0.5) {
        TrkCandidateL.push_back(mom);
    }
    if(0.5<eta && eta<1.){ 
        TrkCandidateR.push_back(mom);
    }
    return kTRUE;
}
//_____________________________________________________________________________
Bool_t StMiniTreeMaker::PassPiCandidate(StMuTrack *track)
{
    if(!IsValidTrack(track)) return kFALSE;
    if(!IsMtdTrack(track)) return kFALSE;
    if(!IsPion(track)) return kFALSE;
    //float eta = mom.pseudoRapidity();
    StThreeVectorF mom = track->momentum();
    PiCandidate.push_back(mom);
    return kTRUE;
}
//_____________________________________________________________________________
void StMiniTreeMaker::ComparePt(StThreeVectorF* vec1 ,StThreeVectorF* vec2)
{
    StThreeVectorF pi1=*vec1;
    StThreeVectorF pi2=*vec2;
    Float_t pt1 = pi1.perp();
    Float_t pt2 = pi2.perp();
    StThreeVectorF vec;
    if(pt1<pt2) {
        vec=*vec1;
        *vec1=*vec2;
        *vec2=vec;
    }   
}
//_____________________________________________________________________________
void StMiniTreeMaker::GetTrigger()
{
    Int_t nVec = PiCandidate.size();
    if(Debug()) cout<<"GetTrigger size: "<<nVec<<endl;
    for(Int_t i=0;i<nVec;i++){
        ComparePt(&PiCandidate[0],&PiCandidate[i]);
    }
}
//_____________________________________________________________________________
Float_t StMiniTreeMaker::PxCal(StThreeVecF trkVec)
{
    Int_t nVec = trkVec.size();
    Float_t Px = 0;
    Float_t PtTrig = PiCandidate[0].perp();
    Float_t PhiTrig = PiCandidate[0].phi();
    if(Debug()) cout<<"PxCal Trigger: "<<PtTrig<<endl;

    for(Int_t i=0;i<nVec;i++){
        Float_t PtCandidate = trkVec[i].perp();
        Float_t PhiCandidate = trkVec[i].phi();
        Float_t DeltaPhi = PhiCandidate-PhiTrig;//candidate-trigger
        if( fabs(DeltaPhi) < 0.5*3.1415926 ) continue;
        //if(Debug()) {cout<<"After DeltaPhi cut: "<<DeltaPhi<<endl;}
        Float_t par = PtCandidate*cos(DeltaPhi);
        Px = Px + par;
        //if(Debug()) cout<<"Px: "<<Px<<endl;
    }   
    return Px;
}
//_____________________________________________________________________________
Bool_t StMiniTreeMaker::PassPxCut(Float_t PxL,Float_t PxR, Short_t cen)
{
    //const Int_t Cen[9+1]={0,2,4,6,8,10,11,13,14,15};
    Int_t i=0;
    if(cen==0 || cen==1)i=0;
    else if(cen==2 || cen==3) i=1;
    else if(cen==4 || cen==5) i=2;
    else if(cen==6 || cen==7) i=3;
    else if(cen==8 || cen==9) i=4;
    else if(cen==10 || cen==11) i=5;
    else if(cen==12 || cen==13) i=6;
    else if(cen==14) i=7;
    else  i=8;

    //if( PxL<PxCutL[i] || PxR<PxCutR[i] ){ 
    //    if(Debug()) {
    //        cout<<"PxCut: "<<endl;
    //        cout<<"centrality: "<<cen<<"  PxL: "<<PxL<<"  PxR: "<<PxR<<endl;
    //        cout<<"PxCutL : "<<PxCutL[i]<<"  PxCutR: "<<PxCutR[i]<<endl;
    //    }
    //    return kTRUE;
    //}
}
//_____________________________________________________________________________
void StMiniTreeMaker::bookHistos()
{
    hEvent = new TH1D("hEvent","Event statistics",7,0,7);
    hEvent->GetXaxis()->SetBinLabel(1,"All events");
    hEvent->GetXaxis()->SetBinLabel(2,"dimuon events");
    hEvent->GetXaxis()->SetBinLabel(3,"After vertex cut");
    hEvent->GetXaxis()->SetBinLabel(4,"2 pion events");
    hEvent->GetXaxis()->SetBinLabel(5,"PxCut");

    hDca = new TH1D("hDca","Dca Distribution",1000,0,10);
}
//_____________________________________________________________________________
void StMiniTreeMaker::bookTree()
{
    LOG_INFO << "StMiniTreeMaker:: book the event tree to be filled."<< endm;
    mEvtTree = new TTree("miniDst","miniDst");
    mEvtTree->SetAutoSave(100000000); // 100 MB

    //  event information
    mEvtTree->Branch("mRunId", &mEvtData.mRunId, "mRunId/I");
    mEvtTree->Branch("mEventId", &mEvtData.mEventId, "mEventId/I");
    mEvtTree->Branch("mGRefMult", &mEvtData.mGRefMult, "mGRefMult/I");
    mEvtTree->Branch("mGRefMultCorr", &mEvtData.mGRefMultCorr, "mGRefMultCorr/F");
    mEvtTree->Branch("mEvtWeight", &mEvtData.mEvtWeight, "mEvtWeight/F");
    mEvtTree->Branch("mCentrality", &mEvtData.mCentrality, "mCentrality/S");
    mEvtTree->Branch("mZDCRate", &mEvtData.mZDCRate, "mZDCRate/F");
    mEvtTree->Branch("mVpdVz", &mEvtData.mVpdVz, "mVpdVz/F");
    mEvtTree->Branch("mTpcVz", &mEvtData.mTpcVz, "mTpcVz/F");
    // track information
    mEvtTree->Branch("mNPiTrk", &mEvtData.mNPiTrk, "mNPiTrk/S");
    mEvtTree->Branch("mPiTrkId", &mEvtData.mPiTrkId, "mPiTrkId[mNPiTrk]/S");
    mEvtTree->Branch("mnSigmaPi", &mEvtData.mnSigmaPi, "mnSigmaPi[mNPiTrk]/F");
    mEvtTree->Branch("mgDca", &mEvtData.mgDca, "mgDca[mNPiTrk]/F");
    mEvtTree->Branch("mPiPt", &mEvtData.mPiPt, "mPiPt[mNPiTrk]/F");
    mEvtTree->Branch("mPiEta", &mEvtData.mPiEta, "mPiEta[mNPiTrk]/F");
    mEvtTree->Branch("mPiPhi", &mEvtData.mPiPhi, "mPiPhi[mNPiTrk]/F");
    mEvtTree->Branch("mPxL", &mEvtData.mPxL, "mPxL/F");
    mEvtTree->Branch("mPxR", &mEvtData.mPxR, "mPxR/F");
}
//_____________________________________________________________________________
void StMiniTreeMaker::printConfig()
{
    const char *decision[2] = {"no","yes"};
    printf("=== Configuration for StMiniTreeMaker ===\n");
    printf("Fill the miniDst tree: %s\n",decision[mFillTree]);
    printf("Fill the QA histo: %s\n",decision[mFillHisto]);
    printf("Maximum |Vz|: %1.2f\n",mMaxTpcVz);
    printf("Maximum |VzDiff|: %1.2f\n",mMaxDiffVz);
    printf("Minimum Track p: %1.2f\n",mMinTrkP);
    printf("Maximum Track |eta| : %1.2f\n",mMaxTrkEta);
    printf("Minimum number of fit hits: %d\n",mMinNHitsFit);
    printf("=======================================\n");
}
