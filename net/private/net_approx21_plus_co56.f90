! ***********************************************************************
!
!   Copyright (C) 2014  Bill Paxton, Frank Timmes
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module net_approx21_plus_co56
      use const_def, only: dp, qp, avo, clight
      use utils_lib, only: is_bad


      implicit none

      

      logical, parameter :: plus_co56 = .true.      
      
      logical, parameter :: reduced_net_for_testing = .true.
      !logical, parameter :: reduced_net_for_testing = .false.




      integer, parameter :: species = 22
      integer :: iso_cid(species) ! these are corresponding chem ids for the isos
         ! e.g., iso_cid(ife52) is = the iso number for fe52 as defined in mesa/chem
      integer :: & ! these are indices in y vector
         ih1, &
         ihe3, &
         ihe4, &
         ic12, &
         in14, &
         io16, &
         ine20, &
         img24, &
         isi28, &
         is32, &
         iar36, &
         ica40, &
         iti44, &
         icr48, &
         icrx, &
         ife52, &
         ife54, &
         ife56, &
         ico56, &
         ini56, &
         ineut, &
         iprot
         
      ! not used, but must be defined for shared code
      integer :: ife53, ife55

      integer, parameter :: &
         approx21_plus_co56_num_mesa_reactions = 94, &
         approx21_plus_co56_nrat = 117
      
      integer, parameter :: num_mesa_reactions = approx21_plus_co56_num_mesa_reactions
      integer, parameter :: num_reactions = approx21_plus_co56_nrat
      
      integer :: rate_id(num_mesa_reactions) ! rate ids for the mesa reactions
         ! e.g., rate_id(ir3a) is reaction id for triple alpha as defined in mesa/rates
      integer :: & ! these are indices in rates arrays
         ir3a, &
         irg3a, &
         ircag, &
         ir1212, &
         ir1216, &
         ir1616, &
         iroga, &
         iroag, &
         irnega, &
         irneag, &
         irmgga, &
         irmgag, &
         irsiga, &
         irmgap, &
         iralpa, &
         iralpg, &
         irsigp, &
         irsiag, &
         irsga, &
         irsiap, &
         irppa, &
         irppg, &
         irsgp, &
         irsag, &
         irarga, &
         irsap, &
         irclpa, &
         irclpg, &
         irargp, &
         irarag, &
         ircaga, &
         irarap, &
         irkpa, &
         irkpg, &
         ircagp, &
         ircaag, &
         irtiga, &
         ircaap, &
         irscpa, &
         irscpg, &
         irtigp, &
         irtiag, &
         ircrga, &
         irtiap, &
         irvpa , &
         irvpg, &
         ircrgp, &
         ircrag, &
         irfega, &
         ircrap, &
         irmnpa, &
         irmnpg, &
         irfegp, &
         irfeag, &
         irniga, &
         irfeap, &
         ircopa, &
         ircopg, &
         irnigp, &

         ! for fe54 photodisintegration
         ir52ng, &
         ir53gn, &
         ir53ng, &
         ir54gn, &
         irfepg, &
         ircogp, &

         ! for he4 photodisintegration
         irheng, &
         irhegn, &
         irhng, &
         irdgn, &
         irdpg, &
         irhegp, &

         ! weak reactions
         irpen, &
         irnep, &
         irn56ec, &
         irco56ec, &

         ! ppchain
         irpp, &
         ir33, &
         irhe3ag, &
         ir_be7_wk_li7, &
         ir_be7_pg_b8, &

         ! cno cycles
         ircpg, &
         irnpg, &
         iropg, &
         irnag, &

         ! for reactions to fe56 
         ir54ng, &
         ir55gn, &
         ir55ng, &
         ir56gn, &
         irfe54ap, &
         irco57pa, &
         irfe56pg, &
         irco57gp, &

         ! for n15 branching
         irn15pa, &
         irn15pg, &

         ! the equilibrium links
         ifa, &
         ifg, &

         irr1, &
         irs1, &
         irt1, &
         iru1, &
         irv1, &
         irw1, &
         irx1, &

         ir1f54, &
         ir2f54, &
         ir3f54, &
         ir4f54, &
         ir5f54, &
         ir6f54, &
         ir7f54, &
         ir8f54, &

         iralf1, &
         iralf2, &

         irfe56_aux1, &
         irfe56_aux2, &
         irfe56_aux3, &
         irfe56_aux4

      ! names
      character (len=40) :: ratnam(num_reactions)


      contains

#define PLUS_CO56
#include "net_approx21_procs.inc"

      end module net_approx21_plus_co56
