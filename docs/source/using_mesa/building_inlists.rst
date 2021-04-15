****************
Building inlists
****************

Inlists for MESA are composed of five main sections labelled ``&star_job``, ``&controls``.
``&eos``,  ``&skap`` and ``&pgstar``. The ``&star_job`` section contains instructions about which MESA modules should be used, while the ``&controls`` section is where the star module options are specified. The ``&kap`` and ``&eos`` sections are where you specify controls for the opacity and the equation of state respectively.
The ``&pgstar`` section contains the commands for creating pgstar plots. 

&star_job
=========

The main modules of MESA (other than **star**) are the :ref:`eos`, the :ref:`kap`, the :ref:`atm`, the nuclear reactions.
In this section of the inlist, you'll have to make choices for which atmosphere and nuclear reactions network you want to use, as well as which nuclear reactions rates you want to use. 
You will also specify here some information about your starting model, and about the output of the evolution.
Here we describe only some of the most commonly used controls. For a complete list of available controls see :ref:`list-starjob`. 

starting model
--------------

To start an evolution you can either create a pre main sequence model and start from there (case 1) or you can start from a previously calculated model (case 2). In the latter case, it is highly recommended to start from a model which was calculated with the same MESA version as the one used for the subsequent evolution.

Case 1:

.. literalinclude:: inlist_example
   :start-after: ! Start Case 1
   :end-before: ! Start Case 2

Case 2:

.. literalinclude:: inlist_example
   :start-after: ! Start Case 2
   :end-before: ! Output history

output
------

There are a number of controls to specify what we want as MESA outputs.
The history file contains global information about the star at some timesteps. The profile files contain information about the interior profiles in the stellar model at a given time. The defaults for these files are located in ``$MESA_DIR/star/defaults``.
If you want output different from the defaults, copy those files in your working directory, rename them (here we use ``my_history_columns.list`` and ``my_profile_columns.list``) , and customize them to your liking. If you include the following lines in your inlist, MESA will use the files sitting in your working directory rather than the default files.

.. literalinclude:: inlist_example
   :start-after: ! Output history
   :end-before: ! Output models

You may want to save the final model of your evolution. In that case you have to tell MESA to do so (default is ``save_model_when_terminate=.false.``) and specify the name of the last model.

.. literalinclude:: inlist_example
   :start-after: ! Output models
   :end-before: ! Modifications to model

Initial composition
-------------------

There are several ways to specify the initial composition. 
(to be filled in).

You can also use pre-defined chemical compositions based on published data. These are set using the control ``ìnitial_zfracs``.
By default, the initial composition in MESA is ``initial_zfracs = 3`` which corresponds to the `GS98 <https://link.springer.com/article/10.1023%2FA%3A1005161325181>`__ metal fraction. There are 8 possible predefined choices for this control in MESA.

If you want for example to use the more recent available solar composition given in `AGSS09 <https://www.annualreviews.org/doi/pdf/10.1146/annurev.astro.46.060407.145222>`__ , you need to set ``initial_zfracs = 6``.
Since it is very important to use the opacity tables which are built using the solar composition used, we also have to set the ``kappa_file_prefix`` to the 2009 solar composition (the default table corresponds to the gs98 composition).

.. literalinclude:: inlist_example
   :start-after: ! Modifications to model
   :end-before: ! Nuclear reactions


nuclear reactions
-----------------

Choice of network of nuclear reactions. This network should be chosen according to the physics to be studied. Choosing a very comprehensive set of nuclear reactions when studying main sequence evolution is not necessary and will slow down the computation considerably. It would however be essential when studying advanced burning stages of evolution. The description of the available nuclear reactions networks in MESA is given in the README file in ``$MESA_DIR/data/net_data``. The default reactions network used by MESA is ``basic.net``.

For example when evolving a stellar model on the horizontal branch (helium burning) this net is insufficient. One possibility is to use the nuclear reactions network called ``pp_and_cno_extras.net``, which provides a more complete coverage for hydrogen and helium burning.

.. literalinclude:: inlist_example
   :start-after: ! Nuclear reactions
   :end-before: ! eos


&controls
=========

Energy equation
---------------

The energy equation can be written in the dLdm or the dedt form in MESA (see `MESAV <https://arxiv.org/pdf/1903.01426.pdf>`__). As explained in `MESAV <https://arxiv.org/pdf/1903.01426.pdf>`__, using the dedt form leads to much better energy conservation. 
The dLdm form is currently the default in MESA. If the dEdt form is preferred it has to be specified in the inlist. 

.. literalinclude:: inlist_example
   :start-after:  ! energy
   :end-before: ! mass and metallicity


Starting model
--------------

The main stellar parameters to specify are its initial mass M, metallicity Z, and helium fraction Y. If only M and Z are specified,the helium content is by default Y=0.24 + 2Z.

.. literalinclude:: inlist_example
   :start-after:  ! mass and metallicity
   :end-before: ! opacity controls



When to stop
------------

Output
------

Opacity controls
----------------

Type2 opacities should be used for extra C/O during and after He burning. To use Type2 opacities one needs to specify a base metallicity, Zbase, which gives the metal abundances previous to any CO enhancement. In regions where central hydrogen is above a given threshold, or the metallicity is not significantly higher than Zbase, Type1 tables are used instead, with blending regions to smoothly transition from one to the other. The reason is that Type1 tables cover a wider range of X and have a higher resolution in Z for each X. For more info, refer to the ``$MESA_DIR/star/defaults/controls.defaults`` file.

If the evolution includes helium burning, type2 opacities should be used.

.. literalinclude:: inlist_example
   :start-after:  ! opacity controls
   :end-before: ! convection 

Convection and Convective boundaries
------------------------------------

Convection in MESA is treated using the MLT theory of convection, and provides different formalisms. By default, MESA uses the `Cox&Giuli 1968 <https://ui.adsabs.harvard.edu/abs/1968pss..book.....C/abstract>`__ formalism.

If you want to use another formalism, for example the `Henyey theory of convection <http://articles.adsabs.harvard.edu/pdf/1965ApJ...142..841H>`__ it can be specified using the ``MLT_option`` control. Several parameters can be specified for this option. The main one is the mixing length parameter. Note that the default value for this parameter in MESA is ``mixing_length_alpha=2``. This value does not come from any calibration.

.. literalinclude:: inlist_example
   :start-after:  ! convection 
   :end-before: ! convective boundaries 

There are two possible criteria that can be used to determine the position of the convective boundaries: the Schwarzchild and Ledoux criteria. By default MESA uses the Schwarzchild criterion. If determined correctly, the position of the convective boundaries should not depend on which criterion is used. But using the Schwarzchild or the Ledoux criterion can lead to different abundance profiles outside the convective region. There are extensive discussions about this topic in the `MESAIV <https://arxiv.org/pdf/1710.08424.pdf>`__ and `MESAV <https://arxiv.org/pdf/1903.01426.pdf>`__ papers. 
Two new algorithms have been introduced in MESA, called **Predictive mixing** (described in `MESAIV <https://arxiv.org/pdf/1710.08424.pdf>`__ )and **Convective PreMixing**, described in `MESAV <https://arxiv.org/pdf/1903.01426.pdf>`__). By default, none of these are used in MESA, which can lead to very incorrectly determined convective boundaries, with important consequences on the evolution of the stellar model. It is therefore highly recommended to use one of these algorithms. 

If using **Convective premixing**, there is no additional parameter to specify.

.. literalinclude:: inlist_example
   :start-after:  ! convective boundaries
   :end-before: ! Predictive

If using **Predictive mixing**, there are additional controls. They are described in `MESAIV <https://arxiv.org/pdf/1710.08424.pdf>`__ .

.. literalinclude:: inlist_example
   :start-after: ! Predictive
   :end-before: ! temp 



Overshooting
------------

Timestep and grid controls
--------------------------

&kap
====

&eos
====
