_Define @Ntab 10 -- Peter Hoeflich (fict 12)
--        @Ntab 15 -- Ron table lower T & rho limit  -- use ronfitab15 for this!
-- for b30*
@Ntab 14 -- Ron, add fictitious points in table
--        @Ns    5  -- old tables
@Ns    6 -- new tables
--        @Nfreq 76
--        @Nfreq   96 -- main standard
--        @Nfreq   36 -- mesa standard
--        @Nfreq   196
--        @Nfreq   146
@Mfreq  130 --  Nfreq+14     -- this alias is needed for formats
@Npe        15         -- number of Envelope tangent rays
@Npc        15         -- number of core tangent rays
-- NYDIM use max(2(@Npe+@Npc), @Nfreq+4))
@Nitem 22    -- for Woosley/Nomoto models (was 20) must be = Natom+7
-- if radius is included
@1stAbun 6   -- it was 5
@Zone  1800    -- for foreign models
@term  *
@WRES  4     -- output channel for results
@wterm  "write(@term,*)" ;
--      @wterm  "print*,";
--NVARS - number of independent variables
include '../obj/nfreq.inc';
PARAMETER (NVARS = 3);                  -- 3 - NOCONV, 4 - CONV
PARAMETER (NFREQ = @Nfreq+4);
PARAMETER(Mzon = 400);
--  PARAMETER(Mzon = 90); -- h toy model
-- PARAMETER(Mzon = 171); -- zoning for W7
-- PARAMETER(Mzon = 340); -- zoning for W7+wind
--PARAMETER(Mzon = 86); -- zoning for W7fhh+wind
--PARAMETER(Mzon = 43); -- zoning for W7fhh
--   PARAMETER (Mzon = 521); -- zoning for ntomi Crab model
PARAMETER(NYDIM = (NVARS + 2 * NFREQ) * Mzon, MAXDER = 4);
-- PARAMETER(NYDIM = 2*NFREQ*Mzon, MAXDER =4);
-- Nydim must be  2*(@Npe+@Npc)*Mzon  to use in Feau
-- PARAMETER(NYDIM = 2*(@Npe+@Npc)*Mzon, MAXDER = 4); -- to use in Feau
Parameter (Is = 5); -- for test -- chego???
--   PARAMETER (NZ = 1200000); --  for Nfreq = 200, Mzon = 200
PARAMETER (NZ = 3000000); --  for Nfreq = 100, Mzon = 600
Parameter(Nstage = 28, Natom = 15);
PARAMETER(KOMAX = 80);  --MAX NUMBER OF STEPS IN PRESCRIBED MOMENTS
Logical LSYSTEM ;-- /.TRUE./; -- for IBM
Parameter(LSystem = .FALSE.);
