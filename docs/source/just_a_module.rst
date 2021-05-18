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


Edit sample_eos.f90 and change

::

  my_mesa_dir = '../..'

to your $MESA_DIR or to a blank string (in which case your $MESA_DIR is automatically picked up)

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
 call eosDT_get

         temperature     0.200000000000E+09
             density     0.100000000000E+03
                   Z     0.200000000000E-01
                   X     0.700000000000E+00
                abar     0.129663534103E+01
                zbar     0.110214003987E+01

             logPgas     0.184308195086E+02
             grad_ad     0.255246371899E+00
                 c_P     0.920814743376E+10


         temperature     0.200000000000E+09
             density     0.100000000856E+03
                   Z     0.200000000000E-01
                   X     0.700000000000E+00
                abar     0.129663534103E+01
                zbar     0.110214003987E+01

             logPgas     0.184308195123E+02
             grad_ad     0.255246371971E+00
                 c_P     0.920814731691E+10


For homework, edit sample_eos.f90 to write out the internal specific energy, specific entropy, partial derivative of 
the specific energy with respect tp temperature, and the free electron degeneracy parameter. 
It can be useful to look at the integer indices contained in $MESA_DIR/eos/public/eos_def.f90. One should find

.. code-block:: console

 ./sample_eos
 loading eos tables
 call eosDT_get


         temperature     0.200000000000E+09
             density     0.100000000000E+03
                   Z     0.200000000000E-01
                   X     0.700000000000E+00
                abar     0.129663534103E+01
                zbar     0.110214003987E+01

             logPgas     0.184308195086E+02
             grad_ad     0.255246371899E+00
                 c_P     0.920814743376E+10

               log E     0.172104888232E+02
               log S     0.941743030071E+01
               dS/dT     0.131585149852E+02
              etaele    -0.564766321638E+01


 
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

