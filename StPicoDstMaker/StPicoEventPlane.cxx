#include <cstdlib> 
#include <bitset>
#include "StPicoEventPlane.h"
#include "StPicoConstants.h"
#include "StEventTypes.h"
#include "StTree.h"
#include "StuRefMult.hh"
#include "TVector2.h"
#include "StMuDSTMaker/COMMON/StMuDst.h" 
#include "StMuDSTMaker/COMMON/StMuTrack.h" 
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StMuDSTMaker/COMMON/StMuMtdHeader.h"
#include "StBTofHeader.h"
#include "StMessMgr.h"
#include "StPicoUtilities.h"

ClassImp(StPicoEventPlane)

StPicoEventPlane::StPicoEventPlane()
{}

StPicoEventPlane::StPicoEventPlane(const StMuDst& muDst)
{
  mQx_eta_pos = 0.;
  mQy_eta_pos = 0.;
  mQx_eta_neg = 0.;
  mQy_eta_neg = 0.;
  mNtrk_eta_pos = 0;
  mNtrk_eta_neg = 0;
  mWeight_eta_pos = 0.;
  mWeight_eta_neg = 0.;

  mQx_chg_pos = 0.;
  mQy_chg_pos = 0.;
  mQx_chg_neg = 0.;
  mQy_chg_neg = 0.;
  mNtrk_chg_pos_eta_pos = 0;
  mNtrk_chg_pos_eta_neg = 0;
  mNtrk_chg_neg_eta_pos = 0;
  mNtrk_chg_neg_eta_neg = 0;
  mWeight_chg_pos = 0.;
  mWeight_chg_neg = 0.;

  mQx_ran_1 = 0.;
  mQy_ran_1 = 0.;
  mQx_ran_2 = 0.;
  mQy_ran_2 = 0.;
  mNtrk_ran_1_eta_pos = 0;
  mNtrk_ran_1_eta_neg = 0;
  mNtrk_ran_2_eta_pos = 0;
  mNtrk_ran_2_eta_neg = 0;
  mWeight_ran_1 = 0.;
  mWeight_ran_2 = 0.;

  Int_t nPrimary = muDst.numberOfPrimaryTracks();
  Int_t trkIndex[nPrimary];
  Int_t counter = 0;
  for(Int_t i=0; i<nPrimary; i++)
    {
      StMuTrack* pTrack = muDst.primaryTracks(i);
      if(!pTrack) continue;
      Float_t pt = pTrack->pt();
      Float_t eta = pTrack->eta();
      Float_t phi = pTrack->phi();

      if(pt<0.15 || pt>2.) continue;
      if(abs(eta)>1.) continue;
      if(pTrack->nHitsFit(kTpcId)<15) continue;
      if(pTrack->nHitsFit(kTpcId)/(1.0*pTrack->nHitsPoss(kTpcId))<0.52) continue;
      if(pTrack->dcaGlobal().mag()>1.) continue;
      trkIndex[counter] = i;
      counter++;

      Short_t charge = pTrack->charge();
      Float_t qx = pt * TMath::Cos(2*phi);
      Float_t qy = pt * TMath::Sin(2*phi);
      if(eta>0)
	{
	  mQx_eta_pos += qx;
	  mQy_eta_pos += qy;
	  mWeight_eta_pos += pt;
	  mNtrk_eta_pos++;
	  if(charge>0)
	    {
	      mQx_chg_pos += qx;
	      mQy_chg_pos += qy;
	      mWeight_chg_pos += pt;
	      mNtrk_chg_pos_eta_pos++;
	    }
	  else
	    {
	      mQx_chg_neg += qx;
	      mQy_chg_neg += qy; 
	      mWeight_chg_neg += pt;
	      mNtrk_chg_neg_eta_pos++;
	    }
	}
      else
	{
	  mQx_eta_neg += qx;
	  mQy_eta_neg += qy;
	  mWeight_eta_neg += pt;
	  mNtrk_eta_neg++;
	  if(charge>0)
	    {
	      mQx_chg_pos += qx;
	      mQy_chg_pos += qy;
	      mWeight_chg_pos += pt;
	      mNtrk_chg_pos_eta_neg++;
	    }
	  else
	    {
	      mQx_chg_neg += qx;
	      mQy_chg_neg += qy; 
	      mWeight_chg_neg += pt;
	      mNtrk_chg_neg_eta_neg++;
	    }
	}
    }

  unsigned int seed = (muDst.event()->eventInfo().runId()*muDst.event()->eventInfo().id()) & 0xffffffff;
  std::srand(seed);
  random_shuffle(trkIndex,trkIndex+counter);
  int nTrk = counter/2*2;
  for(int i=0; i<nTrk; i++)
    {
      StMuTrack* pTrack = muDst.primaryTracks(trkIndex[i]);
      Float_t pt = pTrack->pt();
      Float_t eta = pTrack->eta();
      Float_t phi = pTrack->phi();
      Float_t qx = pt * TMath::Cos(2*phi);
      Float_t qy = pt * TMath::Sin(2*phi);
      if(i<nTrk/2)
	{
	  mQx_ran_1 += qx;
	  mQy_ran_1 += qy;
	  mWeight_ran_1 += pt;
	  if(eta>0) mNtrk_ran_1_eta_pos ++;
	  else      mNtrk_ran_1_eta_neg ++;
	}
      else
	{
	  mQx_ran_2 += qx;
	  mQy_ran_2 += qy;
	  mWeight_ran_2 += pt;
	  if(eta>0) mNtrk_ran_2_eta_pos ++;
	  else      mNtrk_ran_2_eta_neg ++;
	}
    }

}

StPicoEventPlane::~StPicoEventPlane()
{ }
