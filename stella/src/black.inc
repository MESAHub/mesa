  BLACK(Lbl,Tpbl)=(exp(-(FREQMN(Lbl)/Tpbl)))/
                    (1.d0-(exp(-(FREQMN(Lbl)/Tpbl))));
  BLACKD(Lbl,Tpbl)=(FREQMN(Lbl)/Tpbl)*(exp(-(FREQMN(Lbl)/Tpbl)))/
      (1.d0-(exp(-(FREQMN(Lbl)/Tpbl))))**2;
