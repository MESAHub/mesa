Just a Module 
=============

This docuemnt describes how to use the MESA equation of state, opacity, and nuclear reaction networks outside of MESA.


1. Equation of State
--------------------

Before starting, study :ref:`Overview of eos module`

::

   mkdir my_eos
   cd my_eos
   cp $MESA_DIR/eos/test/src/sample_eos.f .


Edit sample_eos.f and change

::

  my_mesa_dir = '../..'

to your $MESA_DIR

::

  my_mesa_dir = '/Users/fxt/mesa_work/r12778'

Now build the executable

:: 

     make

And finally run the executable, named sample 

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



For homework, edit sample_eos.f to add writing out the internal energy, specific entropy, partial derivative of 
the specific energy with respect tp temperature, and the free electron degeneracy parameter. For this assignment
it can be useful to look at the result arrays returned by subroutine eosDT_get in $MESA_DIR/eos/public/eos_lib.f90 
and the integer indices in contained $MESA_DIR/eos/public/eos_def.f90

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

3. Nuclear Reaction Networks
----------------------------

