Just a Module 
=============

This document describes how to use the MESA equation of state, opacity, and nuclear reaction networks outside of MESA.


| Download from Zenodo this set of directories and makefiles: `Just A Module <http://doi.org/10.5281/zenodo.4763740>`_.
| Unpack the zip file anywhere you prefer.


1. Equation of State
--------------------

| Before starting, explore :ref:`Overview of eos module` and :ref:`eos module controls`.
| From where you unpacked the zip file

::

   cd my_eos
   cp $MESA_DIR/eos/test/src/sample_eos.f90 .


Edit sample_eos.f90 and change the variable 


::

  my_mesa_dir = '../..'

to your $MESA_DIR, or use a blank string, in which case your $MESA_DIR is automagically used

::

  my_mesa_dir = ''

while you are editing sample_eos.f90, take some time to explore the source code. Save and exit sample_eos.f90.
Now take some time to explore the makefile, and use it to build the executable

:: 

     make

Run the executable, named sample_eos

.. code-block:: console

 ./sample_eos
 loading eos tables

 give the temperature, density, and mass fractions (h1, he4, c12, n14, o16, ne20, mg24) =>
 hit return for T = 1e9 K, Rho = 1e9 g/cc, x(c12) = 1 ; enter -1 to stop


 T     =  1.000000E+09 Rho   =  1.000000E+04 abar  =  1.200000E+01 zbar  =  6.000000E+00
 h1    =  0.000000E+00 he4   =  0.000000E+00 c12   =  1.000000E+00 n14   =  0.000000E+00
 o16   =  0.000000E+00 ne20  =  0.000000E+00 mg24  =  0.000000E+00
  
 quantity      value          d/d(Rho)       d/d(T)       d^2/d(Rho)^2   d^2/d(Rho)d(T) d^2/d(T)^2
 p tot   =   3.033664E+21   4.616282E+16   1.098950E+13   7.389591E+11  -2.632100E+07   3.148155E+04
 p gas   =   5.117530E+20   4.616282E+16   9.018528E+11   7.389591E+11  -2.632100E+07   1.218618E+03
 p rad   =   2.521911E+21   0.000000E+00   1.008764E+13   0.000000E+00   0.000000E+00   3.026293E+04

 e tot   =   8.606704E+17  -7.955833E+13   3.425527E+09   1.626134E+10  -3.604408E+05   1.309327E+01
 e gas   =   1.040971E+17  -3.900997E+12   3.992340E+08   1.129877E+09  -5.781150E+04   4.014391E+00
 e rad   =   7.565733E+17  -7.565733E+13   3.026293E+09   1.513147E+10  -3.026293E+05   9.078880E+00

 s tot   =   1.476688E+09  -1.098950E+05   3.425527E+00   2.186704E+01  -3.604408E-04   9.667743E-09
 s gas   =   4.679233E+08  -9.018528E+03   3.992340E-01   1.691755E+00  -5.781150E-05   3.615157E-09
 s rad   =   1.008764E+09  -1.008764E+05   3.026293E+00   2.017529E+01  -3.026293E-04   6.052587E-09

 n_ion   =   5.018451E+26   5.018451E+22   0.000000E+00
 n_ele   =   3.200158E+27   2.870797E+23   1.414171E+18
 eta_e   =  -4.182056E+00   9.451310E-05  -1.319764E-09
 cv      =   3.425527E+09  -3.604408E+05   1.309327E+01
 cp      =   2.958707E+10  -5.803942E+06   2.045230E+02
 gamma_1 =   1.314315E+00   5.192048E-06  -3.807384E-10
 gamma_2 =   1.322910E+00   9.783797E-07  -1.079536E-10
 gamma_3 =   1.320812E+00   2.002091E-06  -1.740076E-10
 grad_ad =   2.440906E-01   5.590452E-07  -6.168458E-11
 chi_t   =   3.622516E+00  -5.143301E-05   2.381262E-09
 chi_d   =   1.521685E-01   1.443976E-05  -5.143301E-10
 c_sound =   6.310334E+08  -2.640054E+04   1.141379E+00

 dsp   =  2.220446E-16 dpe   =  7.993606E-15 dsp   = -2.331468E-15

 give the temperature, density, and mass fractions (h1, he4, c12, n14, o16, ne20, mg24) =>
 hit return for T = 1e9 K, Rho = 1e9 g/cc, x(c12) = 1 ; enter -1 to stop

 -1
 STOP normal termination


For homework, edit sample_eos.f90 to write out :math:`\partial{T}/\partial{\rho}|_{S}`.
As mentioned in sample_eos.f90, it can be useful to look at the integer indices contained in $MESA_DIR/eos/public/eos_def.f90. 


 
2. Opacity 
----------

| Before starting, explore :ref:`Overview of kap module` and :ref:`kap module controls`.
| From where you unpacked the zip file

::

   cd my_kap
   cp $MESA_DIR/kap/test/src/sample_kap.f90 .
   cp $MESA_DIR/kap/test/sample_kap_agb.model .


Edit sample_kap.f90 and change

::

  my_mesa_dir = '../..'

to your $MESA_DIR or to a blank string (in which case your $MESA_DIR is automatically picked up)

::

  my_mesa_dir = ''

while you are editing sample_kap.f90, take some time to explore the source code. Save and exit sample_kap.f90.
Now take some time to explore the makefile, and use it to build the executable

:: 

     make

Run the executable, named sample_kap

.. code-block:: console

 ./sample_kap
  Npts        1331
 Nspec          31

 Z_init   1.0000000000000000E-002

 write kap_test.data

Exlore the output with, for example, 

.. code-block:: console

 head -4 kap_test.data
                        grid                     log_T                    log_Rho                      kappa                   kappa_CO                dlnK_dlnRho                  dlnK_dlnT
                           1   3.5585465937700458E+000   -8.4473997504616456E+000    1.7963661540128417E-003    1.7963661540128417E-003    5.9324713626960102E-001    7.7443291473465914E+000
                           2   3.5585885995787634E+000   -8.4471065383083204E+000    1.7984321714713182E-003    1.7984321714713182E-003    5.9331124051888307E-001    7.7411914141868570E+000
                           3   3.5586446634268447E+000   -8.4467158120772723E+000    1.8011910725944315E-003    1.8011910725944315E-003    5.9339708647788947E-001    7.7370154334784580E+000




3. Nuclear Reaction Networks
----------------------------

| Before starting, explore :ref:`Overview of net module` and :ref:`Reaction Networks`.
| From where you unpacked the zip file

::

   cd my_net
   cp $MESA_DIR/net/test/src/sample_net.f90 .


Edit sample_net.f90 and change

::

  my_mesa_dir = '../..'

to your $MESA_DIR or to a blank string (in which case your $MESA_DIR is automatically picked up)

::

  my_mesa_dir = ''

while you are editing sample_net.f90, take some time to explore the source code. Save and exit sample_net.f90.
Now take some time to explore the makefile, and use it to build the executable

:: 

     make

Run the executable, named sample 

.. code-block:: console

 ./sample_net 
 load basic.net
                                                   logT    8.0000000000000000D+00
                                                 logRho    6.0000000000000000D+00
                                                eps_nuc    7.0567990734355760D+08

 

For homework, edit sample_net.f90 to write out the initial composition and the net neutrino loss rate.
One should find


.. code-block:: console

 ./sample_net
 load basic.net
 initial 1H     7.587664E-01
 initial 4He    2.395223E-01
 initial 24Mg   1.711250E-03
                                                   logT    8.0000000000000000D+00
                                                 logRho    6.0000000000000000D+00
                                                eps_nuc    7.0567990734355760D+08
                                                eps_neu    1.7599406836404651D+08

