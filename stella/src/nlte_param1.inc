--- F2T - Fortran to Trefor converter --
-- *** MODULE NLTE_PARAM1 ***
--...T anslated by Pacific-Sierra Research vf90 Personal 3.4N5 19:06:55 1/24/03
--...S itches:
--       _INCLUDE datadefs;
INTEGER KIS;
PARAMETER (KIS = 8);
INTEGER IIadi;
PARAMETER (IIadi = 149);
INTEGER LDR;
PARAMETER (LDR = 394);
INTEGER LDRL;
PARAMETER (LDRL = 50);
INTEGER LEX;
PARAMETER (LEX = 5);
INTEGER NFR;
PARAMETER (NFR = 27);
INTEGER IPHOTNM;
PARAMETER (IPHOTNM = 106);
--
--   NUMNT   - Maximum number of depth points
--   NTZO    - Maximum number of FM intervals
--   INMAX   - Maximum number of simultaneously calculated ionis.st.
--   KEL     - Maximum number of calculated elements
--   LPH     - Maximum number of levels per ionis.st.
--   LEVM    - Maximum number of levels per element
--   NCORE   - Number of core rays
INTEGER NUMNT;
PARAMETER (NUMNT = 41);
INTEGER NTZO;
PARAMETER (NTZO = 7);
INTEGER INMAX;
PARAMETER (INMAX = 75);
INTEGER KEL;
PARAMETER (KEL = 30);
INTEGER LPH;
PARAMETER (LPH = 50);
INTEGER LEVM;
PARAMETER (LEVM = 151);
INTEGER NCORE;
--
--   NFRT    - number of nue-points of the standard Kurucz grid
--   NMERGE  - step size of add. Kurucz grid points
--   NADD    - number of add. nue-points due to shocks
--   NFS     - number of add. nue-points due to line blocking - Sampling M.
--   IXBL    - number of add. nue-points due to line blocking - Formal Int.
--             if OPTADD = T (see MAININP) IXBL should be equal NFS
--   IX1     - max of total number of nue-points
INTEGER NFRT;
PARAMETER (NFRT = 342);
INTEGER NMERGE                             ;--MOD
PARAMETER (NMERGE = 1);
INTEGER NADD;
PARAMETER (NADD = 1200);
INTEGER NFS;
PARAMETER (NFS = 2000);
INTEGER IXBL;
PARAMETER (IXBL = 2000);
INTEGER IX1;
PARAMETER (IX1 = NFR * 6 * LEX + 21 + NFRT / NMERGE + 2 + NADD + NFS);
--                                       4375 !MOD
--
--   MODELN  - available number of Kurucz-models (EMTFLUX_NEW)
--   MODPL   - number of external  HOPF-PARAMETER-SETS
--   NFRE    - number of nue-points for exact photo-data
--   NQ      - max of main quantum number of approx. photo data
--   NFM     - number of grid-points for which the line-force is calc.
--   L3      - length of a record of the line-list
--   L3R     - Reduced L3 value used for temporary arrays
INTEGER MODELN;
PARAMETER (MODELN = 60);
INTEGER MODPL;
PARAMETER (MODPL = 75);
INTEGER NFRE;
PARAMETER (NFRE = 20);
INTEGER NQ;
PARAMETER (NQ = 10);
INTEGER NFM;
PARAMETER (NFM = NUMNT * 3 - 2);
INTEGER L3;
PARAMETER (L3 = 50000);
INTEGER L3R;
PARAMETER (L3R = L3);
;
;
--  parameters for storage of photoionization values
--   LP1     - Number of levels for 2. ioiz. - a change of  LP. requires
--   LP2     - Number of levels for 3. ioiz. - a modification of the
--   LP3     - Number of levels for 4. ioiz. - structure of LMERK, IGYM
--                                           - NOTE ADI IS WATCHING YOU!*?
--   LPN     - Max. number of photoioniz. per level
--   IIR     - Integr. parts of photoioniz. rates - changed from  2 to  3
--   IIG     - I. of ph. r. for special exponents - changed from 11 to 13
--   KAUM    - Number of elements for which Auger-ioniz. rates are calc.
--   IAUM    - Number of iois.st. for which Auger-ioniz. rates are calc.
--   LMERK   - Max. number of photoioniz. per iois.st.
INTEGER LP1                                ;--from 30-11.99
PARAMETER (LP1 = 34);
INTEGER LP2;
PARAMETER (LP2 = 22);
INTEGER LP3;
PARAMETER (LP3 = 8);
INTEGER LPN;
PARAMETER (LPN = 4);
INTEGER IIR;
PARAMETER (IIR = 3);
INTEGER IIG                                ;--MOD
PARAMETER (IIG = 13);
INTEGER KAUM;
PARAMETER (KAUM = 11);
INTEGER IAUM;
PARAMETER (IAUM = 3);
INTEGER LMERK;
PARAMETER (LMERK = LPH + LP1 + LP2 + LP3);
;
--  only for formal.f90
--   NUMWAVE - Max. number of nue-points representing the spectral resol.
--             at the moment NUMWAVE < IX1 (reason:XHCONT, XJCONT)
--   NDOUT   - Number of depth points to print in stdout
INTEGER NUMWAVE                      ;--<= IX1 from 601+IXBL 25.1.00
PARAMETER (NUMWAVE = IX1);
INTEGER NDOUT                              ;--9
PARAMETER (NDOUT = 9);
;
-- other parameters
--   NCALC   - number of most important elements
--   NUMEL1  - number of elements for NITC = 1
--   NUMEL   - number of elements for NITC = 2
--   IONMAX  - max number of ionis.st. for dep. coef. c.
--   LEVNUM  - max of total number of levels
--   LINM    - max number of lines per ionis.st.
--   LINMC   - max number of col. trans. per ionis.st.
--   LINMAX  - max number of transitions
--   LININC  - LPH*INMAX
--   NUMKEL  - NUMNT*KEL
--   KELKIS  - KEL*KIS
INTEGER NCALC;
PARAMETER (NCALC = 16);
INTEGER NUMEL1;
PARAMETER (NUMEL1 = 2);
INTEGER NUMEL;
PARAMETER (NUMEL = 9);
INTEGER IONMAX;
PARAMETER (IONMAX = 3);
INTEGER LEVNUM;
PARAMETER (LEVNUM = LPH * IIadi);
INTEGER LINM;
PARAMETER (LINM = LPH * (LPH - 1) / 4);
INTEGER LINMC;
PARAMETER (LINMC = LINM * 3 / 5);
INTEGER LINMAX;
PARAMETER (LINMAX = LPH * (LPH - 1) / 2 + LPH);
INTEGER LININC;
PARAMETER (LININC = LPH * INMAX);
INTEGER NUMKEL;
PARAMETER (NUMKEL = NUMNT * KEL);
INTEGER KELKIS;
PARAMETER (KELKIS = KEL * KIS);
;
-- the basic directory name for atomic data files is initialized
INTEGER LPNAM;
PARAMETER (LPNAM = 260);
CHARACTER ABDIRNAM*260;
COMMON /NLTE_PARAM1/NCORE;
COMMON /NLTE_PARAM1_C/ABDIRNAM;
-- *** END MODULE NLTE_PARAM1 ***
