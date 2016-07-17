#ifndef StPicoEventPlane_hh
#define StPicoEventPlane_hh

class StMuDst;

#include "TVector2.h"

class StPicoEventPlane : public TObject {
public:
  StPicoEventPlane();
  ~StPicoEventPlane();
  StPicoEventPlane(const StMuDst& muDst) ;
  void Clear(const Option_t*) {}
  
  TVector2 Q() const         { return TVector2(mQx_eta_pos+mQx_eta_neg,mQy_eta_pos+mQy_eta_neg); }

  TVector2 Q_eta_pos() const           { return TVector2(mQx_eta_pos,mQy_eta_pos); }
  TVector2 Q_eta_neg() const           { return TVector2(mQx_eta_neg,mQy_eta_neg); }
  Short_t nTrk_eta_pos() const         { return mNtrk_eta_pos;                     }
  Short_t nTrk_eta_neg() const         { return mNtrk_eta_neg;                     }
  Float_t weight_eta_pos() const       { return mWeight_eta_pos;                   }
  Float_t weight_eta_neg() const       { return mWeight_eta_neg;                   }

  TVector2 Q_chg_pos() const           { return TVector2(mQx_chg_pos,mQy_chg_pos); }
  TVector2 Q_chg_neg() const           { return TVector2(mQx_chg_neg,mQy_chg_neg); }
  Short_t nTrk_chg_pos_eta_pos() const { return mNtrk_chg_pos_eta_pos;             }
  Short_t nTrk_chg_pos_eta_neg() const { return mNtrk_chg_pos_eta_neg;             }
  Short_t nTrk_chg_neg_eta_pos() const { return mNtrk_chg_neg_eta_pos;             }
  Short_t nTrk_chg_neg_eta_neg() const { return mNtrk_chg_neg_eta_neg;             }
  Float_t weight_chg_pos() const       { return mWeight_chg_pos;                   }
  Float_t weight_chg_neg() const       { return mWeight_chg_neg;                   }

  TVector2 Q_ran_1() const             { return TVector2(mQx_ran_1,mQy_ran_1);     }
  TVector2 Q_ran_2() const             { return TVector2(mQx_ran_2,mQy_ran_2);     }
  Short_t nTrk_ran_1_eta_pos() const   { return mNtrk_ran_1_eta_pos;               }
  Short_t nTrk_ran_1_eta_neg() const   { return mNtrk_ran_1_eta_neg;               }
  Short_t nTrk_ran_2_eta_pos() const   { return mNtrk_ran_2_eta_pos;               }
  Short_t nTrk_ran_2_eta_neg() const   { return mNtrk_ran_2_eta_neg;               }
  Float_t weight_ran_1() const         { return mWeight_ran_1;                     }
  Float_t weight_ran_2() const         { return mWeight_ran_2;                     }
      
protected: //these are written out
  
  Float_t mQx_eta_pos;
  Float_t mQy_eta_pos;
  Float_t mQx_eta_neg;
  Float_t mQy_eta_neg;
  Short_t mNtrk_eta_pos;
  Short_t mNtrk_eta_neg;
  Float_t mWeight_eta_pos;
  Float_t mWeight_eta_neg;

  Float_t mQx_chg_pos;
  Float_t mQy_chg_pos;
  Float_t mQx_chg_neg;
  Float_t mQy_chg_neg;
  Short_t mNtrk_chg_pos_eta_pos;
  Short_t mNtrk_chg_pos_eta_neg;
  Short_t mNtrk_chg_neg_eta_pos;
  Short_t mNtrk_chg_neg_eta_neg;
  Float_t mWeight_chg_pos;
  Float_t mWeight_chg_neg;

  Float_t mQx_ran_1;
  Float_t mQy_ran_1;
  Float_t mQx_ran_2;
  Float_t mQy_ran_2;
  Short_t mNtrk_ran_1_eta_pos; 
  Short_t mNtrk_ran_1_eta_neg; 
  Short_t mNtrk_ran_2_eta_pos; 
  Short_t mNtrk_ran_2_eta_neg; 
  Float_t mWeight_ran_1;
  Float_t mWeight_ran_2;
    
  ClassDef(StPicoEventPlane,1)
};

#endif
