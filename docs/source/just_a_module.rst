Just a Module 
=============

This document describes how to use the MESA equation of state, opacity, and nuclear reaction networks outside of MESA.
This document is current as of MESA version r15140.

Download from Zenodo this set of directories and makefiles:
`Just A Module <http://doi.org/10.5281/zenodo.4760707>`_
Unpack the zip file anywhere you prefer.


1. Equation of State
--------------------

Before starting, explore :ref:`Overview of eos module` and :ref:`eos module controls`.
From where you unpacked the zip file

::

   cd my_eos
   cp $MESA_DIR/eos/test/src/sample_eos.f .


Edit sample_eos.f and change

::

  my_mesa_dir = '../..'

to your $MESA_DIR

::

  my_mesa_dir = '/Users/fxt/mesa_work/r12778'

while you are editing the file, explore the source code. Save and exit sample_eos.f.
Now explore the makefile, and use it to build the executable

:: 

     make

Run the executable, named sample 

.. code-block:: console

 ./sample
 loading eos tables
 call eosDT_get

         temperature     0.200000000000E+09
             density     0.100000000000E+03
                   Z     0.200000000000E-01
                   X     0.700000000000E+00
                abar     0.129663534103E+01
                zbar     0.110214003987E+01

             logPgas     0.184313660755E+02
             grad_ad     0.255248168879E+00
                 c_P     0.920285316872E+10


         temperature     0.200000000000E+09
             density     0.100000001167E+03
                   Z     0.200000000000E-01
                   X     0.700000000000E+00
                abar     0.129663534103E+01
                zbar     0.110214003987E+01

             logPgas     0.184313660806E+02
             grad_ad     0.255248168970E+00
                 c_P     0.920285300909E+10



For homework, edit sample_eos.f to write out the internal energy, specific entropy, partial derivative of 
the specific energy with respect tp temperature, and the free electron degeneracy parameter. For this assignment
it can be useful to look at the result arrays returned by subroutine eosDT_get in $MESA_DIR/eos/public/eos_lib.f90 
and the integer indices in contained $MESA_DIR/eos/public/eos_def.f90. One should find

.. code-block:: console

 ./sample
 loading eos tables
 call eosDT_get

         temperature     0.200000000000E+09
             density     0.100000000000E+03
                   Z     0.200000000000E-01
                   X     0.700000000000E+00
                abar     0.129663534103E+01
                zbar     0.110214003987E+01

             logPgas     0.184308197524E+02
             grad_ad     0.255246340077E+00
                 c_P     0.920820660614E+10
               log E     0.172104903449E+02
               log S     0.941743116140E+01
               dS/dT     0.131585701995E+02
              etaele    -0.564766316514E+01


 
2. Opacity 
----------

Before starting, explore :ref:`Overview of kap module` and :ref:`kap module controls`.
From where you unpacked the zip file

::

   cd my_kap
   cp $MESA_DIR/kap/test/src/sample_kap.f90 .
   cp $MESA_DIR/kap/test/sample_kap_agb.model .


Edit sample_kap.f90 and change

::

  my_mesa_dir = '../..'

to your $MESA_DIR

::

  my_mesa_dir = '/Users/fxt/mesa_work/r12778'

while you are editing the file, explore the source code. Save and exit sample_kap.f90.
Now explore the makefile, and use it to build the executable

:: 

     make

Run the executable, named sample 

.. code-block:: console

 ./sample
  Npts        1331
 Nspec          31

 Z_init   1.0000000000000000E-002

 write kap_test.data

Exlore the output with, for example, 

.. code-block:: console

 head -4 kap_test.data
                        grid                     log_T                    log_Rho                      kappa                   kappa_CO                dlnK_dlnRho                  dlnK_dlnT
                           1   3.5585465937700458E+000   -8.4473997504616456E+000    1.7963661540119481E-003    1.7963661540119481E-003    5.9324713626873704E-001    7.7443291473431390E+000
                           2   3.5585885995787634E+000   -8.4471065383083204E+000    1.7984321714704225E-003    1.7984321714704225E-003    5.9331124051801798E-001    7.7411914141834002E+000
                           3   3.5586446634268447E+000   -8.4467158120772723E+000    1.8011910725935314E-003    1.8011910725935314E-003    5.9339708647702316E-001    7.7370154334749985E+000




3. Nuclear Reaction Networks
----------------------------

Before starting, explore :ref:`Overview of net module` and :ref:`Reaction Networks`.
From where you unpacked the zip file

::

   cd my_net
   cp $MESA_DIR/net/test/src/sample_net.f .


Edit sample_net.f and change

::

  my_mesa_dir = '../..'

to your $MESA_DIR

::

  my_mesa_dir = '/Users/fxt/mesa_work/r12778'

while you are editing the file, explore the source code. Save and exit sample_net.f.
Now explore the makefile, and use it to build the executable

:: 

     make

Run the executable, named sample 

.. code-block:: console

 ./sample 
 load basic.net
                                                   logT    8.0000000000000000D+00
                                                 logRho    6.0000000000000000D+00
                                                eps_nuc    7.0567990734355760D+08

 

For homework, edit sample_net.f to add writing out the initial composition and the net neutrino loss rate.


.. code-block:: console

 ./sample 
 load basic.net
 initial 1H     7.587664E-01
 initial 4He    2.395223E-01
 initial 24Mg   1.711250E-03
                                                   logT    8.0000000000000000D+00
                                                 logRho    6.0000000000000000D+00
                                                eps_nuc    7.0567990734355760D+08
                                                eps_neu    1.7599406836404651D+08

