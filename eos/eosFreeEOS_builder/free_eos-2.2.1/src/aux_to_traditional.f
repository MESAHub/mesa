C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: aux_to_traditional.f 842 2008-07-07 21:18:51Z airwin $
C
C       For the latest version of this source code, please contact
C       Alan W. Irwin
C       Department of Physics and Astronomy
C       University of Victoria,
C       Box 3055
C       Victoria, B.C., Canada
C       V8W 3P6
C       e-mail: irwin@beluga.phys.uvic.ca.
C
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C
C       End of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C*******************************************************************************
      subroutine aux_to_traditional(
     &  nions, max_index, naux, inv_aux,
     &  aux,
     &  h,
     &  hd,
     &  he,
     &  rl,
     &  h2plus,
     &  sum0,
     &  sum2,
     &  iextraoff, maxnextrasum, nxextrasum,
     &  extrasum,
     &  xextrasum)
      implicit none
      integer nions, max_index, naux, inv_aux(naux),
     &  iextraoff, maxnextrasum, nxextrasum
      double precision
     &  aux(naux), auxf(naux), auxt(naux), daux_dv(nions+2, naux),
     &  h, hf, ht, h_dv(nions+2),
     &  hd, hdf, hdt, hd_dv(nions+2),
     &  he, hef, het, he_dv(nions+2),
     &  rl, rf, rt, r_dv(nions+2),
     &  h2plus, h2plusf, h2plust, h2plus_dv(nions+2),
     &  sum0, sum0f, sum0t, sum0_dv(nions+2),
     &  sum2, sum2f, sum2t, sum2_dv(nions+2),
     &  extrasum(maxnextrasum), extrasumf(maxnextrasum),
     &  extrasumt(maxnextrasum), extrasum_dv(nions+2, maxnextrasum),
     &  xextrasum(nxextrasum), xextrasumf(nxextrasum),
     &  xextrasumt(nxextrasum), xextrasum_dv(nions+2, nxextrasum)
      logical ifcsum, ifextrasum, ifxextrasum, if_dv
      integer index, jndex, if_mc
C      Copy from aux to traditional form (exclude all derivatives
C      to reduce cross-talk as much as possible).
      if(inv_aux(1).ne.0) h = aux(1)
      if(inv_aux(2).ne.0) hd = aux(2)
      if(inv_aux(3).ne.0) he = aux(3)
      if(inv_aux(4).ne.0) rl = aux(4)
      if(inv_aux(5).ne.0) h2plus = aux(5)
      if(inv_aux(6).ne.0) sum0 = aux(6)
      if(inv_aux(7).ne.0) sum2 = aux(7)
C      hf = auxf(1)
C      hdf = auxf(2)
C      hef = auxf(3)
C      rf = auxf(4)
C      h2plusf = auxf(5)
C      sum0f = auxf(6)
C      sum2f = auxf(7)
C      ht = auxt(1)
C      hdt = auxt(2)
C      het = auxt(3)
C      rt = auxt(4)
C      h2plust = auxt(5)
C      sum0t = auxt(6)
C      sum2t = auxt(7)
      do jndex = 1, maxnextrasum
        if(inv_aux(iextraoff+jndex).ne.0)
     &    extrasum(jndex) = aux(iextraoff+jndex)
C        extrasumf(jndex) = auxf(iextraoff+jndex)
C        extrasumt(jndex) = auxt(iextraoff+jndex)
      enddo
      do jndex = 1, nxextrasum
        if(inv_aux(iextraoff+maxnextrasum+jndex).ne.0)
     &    xextrasum(jndex) = aux(iextraoff+maxnextrasum+jndex)
C        xextrasumf(jndex) = auxf(iextraoff+maxnextrasum+jndex)
C        xextrasumt(jndex) = auxt(iextraoff+maxnextrasum+jndex)
      enddo
C      n.b. daux_dv(index,jaux)  has compact index for
C      index and *uncompact* index for jaux
C      do index = 1, max_index
C        h_dv(index) = daux_dv(index,1)
C        hd_dv(index) = daux_dv(index,2)
C        he_dv(index) = daux_dv(index,3)
C        r_dv(index) = daux_dv(index,4)
C        h2plus_dv(index) = daux_dv(index,5)
C        sum0_dv(index) = daux_dv(index,6)
C        sum2_dv(index) = daux_dv(index,7)
C        do jndex = 1, maxnextrasum
C          extrasum_dv(index, jndex) = daux_dv(index,iextraoff+jndex)
C        enddo
C        do jndex = 1, nxextrasum
C          xextrasum_dv(index, jndex) =
C     &      daux_dv(index,iextraoff+maxnextrasum+jndex)
C        enddo
C      enddo
      return
      entry traditional_to_aux(
     &  if_mc, ifcsum, ifextrasum, ifxextrasum, if_dv,
     &  h, hf, ht, h_dv,
     &  hd, hdf, hdt, hd_dv,
     &  he, hef, het, he_dv,
     &  rl, rf, rt, r_dv,
     &  h2plus, h2plusf, h2plust, h2plus_dv,
     &  sum0, sum0f, sum0t, sum0_dv,
     &  sum2, sum2f, sum2t, sum2_dv,
     &  iextraoff, maxnextrasum, nxextrasum,
     &  extrasum, extrasumf, extrasumt, extrasum_dv,
     &  xextrasum, xextrasumf, xextrasumt, xextrasum_dv,
     &  nions, max_index, naux,
     &  aux, auxf, auxt, daux_dv)
C      Copy from traditional form of auxiliary variables *and their
C      derivatives* to aux (and its derivatives).
      aux(1) = h
      aux(2) = hd
      aux(3) = he
      aux(4) = rl
      aux(5) = h2plus
      auxf(1) = hf
      auxf(2) = hdf
      auxf(3) = hef
      auxf(4) = rf
      auxf(5) = h2plusf
      auxt(1) = ht
      auxt(2) = hdt
      auxt(3) = het
      auxt(4) = rt
      auxt(5) = h2plust
      if(if_mc.eq.1) then
C        if_mc.eq.1 is a special case where you need both sum0 and sum2
C        (which carry the completely ionized component of these sums
C        divided by (rho*avogadro) in this special case) copied back to aux.
        aux(6) = sum0
        aux(7) = sum2
      endif
      if(ifcsum) then
        aux(6) = sum0
        aux(7) = sum2
        auxf(6) = sum0f
        auxf(7) = sum2f
        auxt(6) = sum0t
        auxt(7) = sum2t
      endif
      if(ifextrasum) then
        do jndex = 1, maxnextrasum
          aux(iextraoff+jndex) = extrasum(jndex)
          auxf(iextraoff+jndex) = extrasumf(jndex)
          auxt(iextraoff+jndex) = extrasumt(jndex)
        enddo
      endif
      if(ifxextrasum) then
        do jndex = 1, nxextrasum
          aux(iextraoff+maxnextrasum+jndex) = xextrasum(jndex)
          auxf(iextraoff+maxnextrasum+jndex) = xextrasumf(jndex)
          auxt(iextraoff+maxnextrasum+jndex) = xextrasumt(jndex)
        enddo
      endif
C      n.b. daux_dv(index,jaux) has compact index for
C      index and *uncompact* index for jaux
      if(if_dv) then
        do index = 1, max_index
          daux_dv(index,1) = h_dv(index)
          daux_dv(index,2) = hd_dv(index)
          daux_dv(index,3) = he_dv(index)
          daux_dv(index,4) = r_dv(index)
          daux_dv(index,5) = h2plus_dv(index)
          if(ifcsum) then
            daux_dv(index,6) = sum0_dv(index)
            daux_dv(index,7) = sum2_dv(index)
          endif
          if(ifextrasum) then
            do jndex = 1, maxnextrasum
              daux_dv(index,iextraoff+jndex) =
     &          extrasum_dv(index, jndex)
            enddo
          endif
          if(ifxextrasum) then
            do jndex = 1, nxextrasum
              daux_dv(index,iextraoff+maxnextrasum+jndex) =
     &          xextrasum_dv(index, jndex)
            enddo
          endif
        enddo
      endif
      end
