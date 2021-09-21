Quickstart
==========

Many questions that appear on the mesa-users list have answers already available somewhere in the MESA world. Navigating through this world is not always straightforward. Here is a little guide on how to do that.

1. Installation
---------------

The first thing to do is to follow exactly the installation steps explained in :ref:`installation:Installing MESA`.
 
In case of problems installing MESA, look at the section on troubleshooting and consult the :ref:`MESA FAQ <FAQ>`.

If you cannot find the answer to your problems in these pages, search the `mesa-users list archive <https://lists.mesastar.org/pipermail/mesa-users/>`_ to see if someone has had a similar problem in the past.

If you still cannot find the answer to your installation problem, post your question to mesa-users@lists.mesastar.org describing the problem.
When you do that, **it is important to follow the clear procedure detailed in** :ref:`the troubleshooting guide <installation:Troubleshooting>`. This tells you which information you need to provide for people to be able to help you.

2. Using MESA
-------------

MESA is built to allow users to run experiments in stellar evolution. As such, a lot of freedom is given to the users in terms of choosing parameters and their values. This can be overwhelming when using MESA for the first time.

The first step is to go through the example detailed here: :ref:`Running MESA <running>`.

There are also more examples located in the subdirectory ``$MESA_DIR/star/test_suite``.
There you will find examples closer to the type of star you want to study.
Some of these test cases are documented in detail in the :ref:`Test suite documentation <test_suite:Test suite>`. 
You can use these cases as starting points for your inlists, and slowly modify them to fit your needs.
Please beware that the examples in the test suite are not necessarily at the correct resolution, since they need to run fast, and you may need to adjust the resolution to have results that are converged. Please note also, and this is very important, that MESA defaults or test_suite inlists will generally NOT be optimal or even acceptable for your particular science cases. It is your responsibility to ensure that the MESA options and controls you choose are appropriate for the physics you want to study. This will usually require appropriate testing and critical analysis of the models obtained.

There is a section in the documentation about `building inlists <using_mesa/building_inlists.html>`__.  

Another source of examples is the `MESA Marketplace <http://cococubed.asu.edu/mesa_market/>`__ and the `MESA Zenodo community <https://zenodo.org/communities/mesa/>`__.
There you will find inlists used in published papers and MESA Summer School material. You will also find MESA community written guidance on using MESA.

3. YouTube videos
-----------------

`Installing MESA on Linux  <https://youtu.be/NmaLHFxpALg>`_

`Installing MESA on osx  <https://youtu.be/mr_A0XrGqNA>`_

`Getting started with MESA  <https://youtu.be/b0bZ9FAgyrg>`_

`Using pgstar for your MESA project <https://youtu.be/JZFa4WURztI>`_


