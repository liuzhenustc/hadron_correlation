#ifndef StPicoMtdTrigger_hh
#define StPicoMtdTrigger_hh

class StMuDst;
class StTriggerData;
#include "TObject.h"

class StPicoMtdTrigger : public TObject {
 public:
  StPicoMtdTrigger();
  StPicoMtdTrigger(const StMuDst& muDst, const int QTtoModule[4][8],
		   const int QTSlewBinEdge[4][16][8], const int QTSlewCorr[4][16][8]);
  ~StPicoMtdTrigger();

  // qt: 1-4, pos: 1-8
  UShort_t   getVpdTacSum()                                 { return mVpdTacSum;             }
  UShort_t   getQTtacSum(const Int_t qt, const Int_t pos)   { return mQTtacSum[qt-1][pos-1]; }
  UShort_t   getMT101Tac(const Int_t qt, const Int_t index) { return mMT101Tac[qt-1][index]; }
  UShort_t   getMT101Id(const Int_t qt, const Int_t index)  { return mMT101Id[qt-1][index];  }
  UChar_t    getTF201TriggerBit()                           { return mTF201TriggerBit;       }
  Char_t     shouldHaveRejectEvent()                        { return mShouldHaveRejectEvent; }

  void       getMaximumQTtac(const Int_t qt, Int_t& pos1, Int_t& pos2);

 protected:
  static const UShort_t mtd_qt_tac_max = 4095;
  static const UShort_t mtd_qt_tac_min = 100; 
  static const UShort_t mtd_qt_tac_diff_range_abs = 600; 


 private:
  UShort_t      mVpdTacSum;      // VPD tac sum
  UShort_t      mQTtacSum[4][8]; // tacSum in 4 QT boards
  UShort_t      mMT101Tac[4][2]; // two largest tacSum from each QT board
  UChar_t       mMT101Id[4][2];  // id of largest tacSum -> position
  UChar_t       mTF201TriggerBit; // final trigger bit to TCU
  Char_t        mShouldHaveRejectEvent; // indication of event status in filtering
                                        // 0 - events not triggered di-muon
                                        // 1 - events should have been rejected 
                                        // if only triggered by di-muon
                                        // 2 - events pass filtering cuts

  ClassDef(StPicoMtdTrigger,3);
};
#endif
