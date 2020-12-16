      subroutine res_adipar(i_param,
     * el, els1, dels, dfsig1, dfsig2, sig1, sig2, dfsig,
     * eltrw1, eltrw2, sgtrw1, sgtrw2,
     * nsel, nsig1, nsig2, itrsig, nsig, istsig,
     * inomd1, iscan, irotkr)
c
c  set or reset parameters to be iterated in calls of adipls
c
c  For i_param = 1, store parameters in para_el, etc. from 
c  parameters in arguments. 
c
c  For i_param = 2, reset parameters in arguments from para_el, ...
c
c  Original version: 8/7/05
c
      implicit double precision(a-h, o-z)
c
c  commons for parameter transmission
c
      common/cadr_param/
     *  para_el, para_els1, para_dels, para_dfsig1, para_dfsig2,
     *  para_sig1, para_sig2, para_dfsig, para_eltrw1, para_eltrw2,
     *  para_sgtrw1, para_sgtrw2
      common/cadi_param/
     *  ipara_nsel, ipara_nsig1, ipara_nsig2, ipara_itrsig, ipara_nsig, 
     *  ipara_istsig, ipara_inomd1, ipara_iscan, ipara_irotkr
c
c  common defining standard input and output, standard error
c
      common/cstdio/ istdin, istdou, istdpr, istder
c 
      save
c
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Enter res_adipar with i_param =',i_param
c
c  test for storing original parameters in common /cadr_param/
c
      if(i_param.eq.1) then
        para_el = el
        para_els1 = els1
        para_dels = dels
        para_dfsig1 = dfsig1
        para_dfsig2 = dfsig2
        para_sig1 = sig1
        para_sig2 = sig2
        para_dfsig = dfsig
        para_eltrw1 = eltrw1
        para_eltrw2 = eltrw2
        para_sgtrw1 = sgtrw1
        para_sgtrw2 = sgtrw2
c
        ipara_nsel = nsel
        ipara_nsig1 = nsig1
        ipara_nsig2 = nsig2
        ipara_itrsig = itrsig
        ipara_nsig = nsig
        ipara_istsig = istsig
        ipara_inomd1 = inomd1
        ipara_iscan = iscan
        ipara_irotkr = irotkr
c
c  test for setting parameters from /cadr_param/ into internal 
c  parameters
c
      else if(i_param.eq.2) then
        el = para_el
        els1 = para_els1
        dels = para_dels
        dfsig1 = para_dfsig1
        dfsig2 = para_dfsig2
        sig1 = para_sig1
        sig2 = para_sig2
        dfsig = para_dfsig
        eltrw1 = para_eltrw1
        eltrw2 = para_eltrw2
        sgtrw1 = para_sgtrw1
        sgtrw2 = para_sgtrw2
c
        nsel = ipara_nsel
        nsig1 = ipara_nsig1
        nsig2 = ipara_nsig2
        itrsig = ipara_itrsig
        nsig = ipara_nsig
        istsig = ipara_istsig
        inomd1 = ipara_inomd1
        iscan = ipara_iscan
        irotkr = ipara_irotkr
      end if
      return 
      end
