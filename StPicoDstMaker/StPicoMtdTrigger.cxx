#include <iostream>
#include <bitset>

#include "StTriggerData.h"
#include "StMessMgr.h"
#include "StMuDSTMaker/COMMON/StMuDst.h" 
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuMtdHeader.h"
#include "StEnumerations.h"
#include "StPicoMtdTrigger.h"

ClassImp(StPicoMtdTrigger)

//----------------------------------------------------------------------------------
StPicoMtdTrigger::StPicoMtdTrigger()
{
  mVpdTacSum = 0;
  memset(mQTtacSum,0,sizeof(mQTtacSum));
  memset(mMT101Tac,0,sizeof(mMT101Tac));
  memset(mMT101Id,0,sizeof(mMT101Id));
  mTF201TriggerBit = 0;
  mShouldHaveRejectEvent = -1;
}

//----------------------------------------------------------------------------------
StPicoMtdTrigger::StPicoMtdTrigger(const StMuDst& muDst, const int QTtoModule[4][8],
				   const int QTSlewBinEdge[4][16][8], const int QTSlewCorr[4][16][8])
{
  StTriggerData *trigger = const_cast<StTriggerData*>(muDst.event()->triggerData());

  // VPD tac sum
  mVpdTacSum = trigger->vpdEarliestTDCHighThr(east) + trigger->vpdEarliestTDCHighThr(west);

  // QT information
  UShort_t mtdQTtac[4][16];
  UShort_t mtdQTadc[4][16];
  for(Int_t i=0; i<32; i++)
    {
      if((i/4)%2==1)
        {
          mtdQTtac[0][i-i/4*2-2] = trigger->mtdAtAddress(i,0);
          mtdQTtac[1][i-i/4*2-2] = trigger->mtdgemAtAddress(i,0);
          mtdQTtac[2][i-i/4*2-2] = trigger->mtd3AtAddress(i,0);
          mtdQTtac[3][i-i/4*2-2] = trigger->mtd4AtAddress(i,0);
        }
      else
	{
	  mtdQTadc[0][i-i/4*2] = trigger->mtdAtAddress(i,0);
	  mtdQTadc[1][i-i/4*2] = trigger->mtdgemAtAddress(i,0);
	  mtdQTadc[2][i-i/4*2] = trigger->mtd3AtAddress(i,0);
	  mtdQTadc[3][i-i/4*2] = trigger->mtd4AtAddress(i,0);
	}
    }

  UShort_t j[2], a[2];
  for(int im=0; im<4; im++)
    {
      for(int i=0; i<8; i++)
	{
	  mQTtacSum[im][i] = 0;
	  // slewing correction
	  for(int k=0; k<2; k++)
	    {
	      j[k] = mtdQTtac[im][i*2+k];
	      a[k] = mtdQTadc[im][i*2+k];

	      int slew_bin = -1;
	      if(a[k]>=0 && a[k]<=QTSlewBinEdge[im][i*2+k][0]) slew_bin = 0;
	      else
		{
		  for(int l=1; l<8; l++)
		    {
		      if(a[k]>QTSlewBinEdge[im][i*2+k][l-1] && a[k]<=QTSlewBinEdge[im][i*2+k][l])
			{
			  slew_bin = l;
			  break;
			}
		    }
		}
	      if(slew_bin>=0)
		j[k] += QTSlewCorr[im][i*2+k][slew_bin];
	    }
	  
          if(j[0]<mtd_qt_tac_min || j[0]>mtd_qt_tac_max || 
             j[1]<mtd_qt_tac_min || j[1]>mtd_qt_tac_max ||
             TMath::Abs(j[0]-j[1])>mtd_qt_tac_diff_range_abs) 
	    {
	      mQTtacSum[im][i] = 0;
	      continue;
	    }
	  
	  // position correction
	  int module = QTtoModule[im][i];
	  mQTtacSum[im][i] = UShort_t( j[0] + j[1] + abs(module-3)*1./8 * (j[0]-j[1]) );
	}
    }

  // MT101
  for(Int_t i = 0; i < 4; i++)
    {
      mMT101Tac[i][0] = (trigger->mtdDsmAtCh(3*i,0)) + ((trigger->mtdDsmAtCh(3*i+1,0)&0x3)<<8);
      mMT101Id[i][0]  = (trigger->mtdDsmAtCh(3*i+1,0)&0xc)>>2;
      mMT101Tac[i][1] = (trigger->mtdDsmAtCh(3*i+1,0)>>4) + ((trigger->mtdDsmAtCh(3*i+2,0)&0x3f)<<4);
      mMT101Id[i][1]  = (trigger->mtdDsmAtCh(3*i+2,0)&0xc0)>>6;
    }

  // TF201
  mTF201TriggerBit = 0;
  UShort_t  decision = trigger->dsmTF201Ch(0);
  for(Int_t i=0; i<8; i++)
    {
      mTF201TriggerBit |= ((decision>>(i+4))&0x1)<<i;
    }
  LOG_DEBUG << "input  = " << (std::bitset<16>) decision << endm;  
  LOG_DEBUG << "output = " << (std::bitset<8>)  mTF201TriggerBit << endm;

  // Header info
  StMuMtdHeader *muMtdHeader = muDst.mtdHeader();
  if(muMtdHeader)
    mShouldHaveRejectEvent = (Char_t)(muMtdHeader->shouldHaveRejectEvent());
  else
    mShouldHaveRejectEvent = -1;
}

//----------------------------------------------------------------------------------
StPicoMtdTrigger::~StPicoMtdTrigger()
{
}


//----------------------------------------------------------------------------------
void StPicoMtdTrigger::getMaximumQTtac(const Int_t qt, Int_t& pos1, Int_t& pos2)
{
  pos1 = 0;
  pos2 = 0;

  UShort_t max1 = 0, max2 = 0;

  for(Int_t i=0; i<8; i++)
    {
      if(max1 < mQTtacSum[qt-1][i])
	{
	  max2 = max1;
	  pos2 = pos1;
	  max1 = mQTtacSum[qt-1][i];
	  pos1 = i+1;
	}
      else if(max2 < mQTtacSum[qt-1][i])
	{
	  max2 = mQTtacSum[qt-1][i];
	  pos2 = i+1;
	}
    }
}
