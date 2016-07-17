const Int_t mMax = 200;
struct StEvtData
{
	// event information
	Int_t    mRunId;
	Int_t    mEventId;
	Int_t    mGRefMult;
    Float_t    mGRefMultCorr;
	Float_t    mEvtWeight;
	Short_t    mCentrality;
	Float_t  mZDCRate;
	Float_t  mVpdVz;
    Float_t  mTpcVz;

	//track information
	Short_t    mNPiTrk;
	Short_t    mPiTrkId[mMax];
    Float_t  mPiPt[mMax];
    Float_t  mPiEta[mMax];
    Float_t  mPiPhi[mMax];
    Float_t  mPxL;
	Float_t  mPxR;
};
