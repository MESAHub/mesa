Real*8  WMOD(@zone, @Nitem), AZinv(Mzon);
Integer NZONE, Nitind, IndAbun(Natom);
Double precision Rhomod(@zone), Tpmod(@zone), Pmod(@zone), ammod(@zone);
Double precision Xnimod(@zone, 2); -- 1st col. Xni,
-- 2nd col. am (in solar units) from foreign model Woosley
dimension XABUN(Natom + 1, @zone); -- foreign model NTominaga
Logical Woomix;
Integer Ninm;
COMMON/nin/ Ninm;
COMMON/ProfMod/WMOD, AZinv, Tpmod, Rhomod, Pmod, ammod, XNimod, XABUN,
IndAbun,NZONE, Woomix;
-- IndAbun  --  Exist/NotExist flag for elements (see Fundrad.inc)
