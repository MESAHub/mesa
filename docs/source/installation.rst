.. highlight:: console

Installing MESA
===============

This page describes how to install MESA.

Prerequisites
-------------

Ensure your system meets the minimum hardware requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The minimum system requirements for MESA are:

-  Mac or Linux operating system
-  64-bit processor
-  8 GB RAM
-  20 GB free disk space
-  Windows users should :ref:`follow the instructions here <windows-install:Installing MESA on Windows>`.

Most laptop or desktop computers built in the last three years will 
satisfy these requirements.

Ensure you have Python (3.5 or newer) installed on your system
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. warning:: As of r24.08.1, building MESA now requires Python (3.5 or newer) to be installed.

Install the MESA SDK
^^^^^^^^^^^^^^^^^^^^

Before you install MESA, you need to get the prerequisites. The `MESA
SDK <http://user.astro.wisc.edu/~townsend/static.php?ref=mesasdk>`__
simplifies this process by providing a prebuilt set of compilers and
run-time libraries that should make your MESA install go
smoothly. Visit the `MESA SDK website
<http://user.astro.wisc.edu/~townsend/static.php?ref=mesasdk>`__ for
the details of setting it up.

If you would prefer to use ifort (the MESA SDK uses gfortran), that is
also an option, so long as you use ifort 14 or later. Even if you choose
to use ifort, you should still visit the MESA SDK website to get a feel
for the other MESA requirements.

Not using the MESA SDK means you'll need to replace the file
:file:`$MESA_DIR/utils/makefile_header` with a version customized to your
system. There's a template to get you started at
:file:`$MESA_DIR/utils/makefile_header_non_mesasdk`.

Regardless of whether you use the MESA SDK or ifort, and whether your
machine runs MacOS or linux, the output of MESA should be bit-for-bit
identical.  If it's not, this is considered to be a bug. (This has
been the case since Release 5819 in early January 2014.)

Installation
------------

Download MESA
^^^^^^^^^^^^^

The simplest way to get the MESA software is to download a zip file of
the `latest MESA release <https://doi.org/10.5281/zenodo.2602941>`__.

The compressed file is about 2GB, so don't worry if it takes a little
while to download.  

The unzipped and installed package will be large, so make sure you have
at least 20 GB free on your disk.

When you unzip the file, it will create a directory named
mesa-\ |version|. This will be your main MESA directory. You are
free to rename it, just make sure to set MESA_DIR accordingly (see the
next section).

You can also download zip files of `older MESA releases <https://doi.org/10.5281/zenodo.2602941>`__.
If you plan to do so, please read :ref:`this FAQ entry <faq:Installing Older Versions of MESA>`.

.. _environment:

Set your environment variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The easiest way to make sure that your system is always configured
appropriately is to define the necessary environment variables in
your `shell start-up file <https://kb.iu.edu/d/abdy>`__. The file that
you need to edit will depend on `which shell you're using
<http://askubuntu.com/questions/590899/how-to-check-which-shell-am-i-using>`__.
You can find out by running ``echo $0``. The default on most Linux
distros is bash, in which case you need to edit ``$HOME/.bashrc``. If
you don't set the environment variables in your shell start-up file,
you will need to re-define them each time you open a new shell.

The exact paths depend on where you installed MESA and which operating
system you are using. After you add these commands to your shell
startup file, don't forget to open a new shell (or ``source`` the
startup file in an existing one).

Here is an example from a machine that uses bash as its shell (and hence
uses export to set variables):

.. code-block:: bash

    # set MESA_DIR to be the directory to which you downloaded MESA
    # The directory shown is only an example and must be modified for your particular system.
    export MESA_DIR=/Users/jschwab/Software/mesa-r21.12.1

    # set OMP_NUM_THREADS to be the number of cores on your machine
    export OMP_NUM_THREADS=2

    # you should have done this when you set up the MESA SDK
    # The directory shown is only an example and must be modified for your particular system.
    export MESASDK_ROOT=/Applications/mesasdk
    source $MESASDK_ROOT/bin/mesasdk_init.sh

    # add shmesa (the MESA command line tool) to your PATH 
    export PATH=$PATH:$MESA_DIR/scripts/shmesa


If your machine uses csh as its shell, use ``setenv`` instead of ``export``.
    
One caveat is that if you initialize the MESA SDK in your shell
profile, you'll always be using the MESA SDK supplied version of gcc
which may be a compatibility issue if you work with other other codes.
Alternative (unsupported) initialization scripts are available `here
<https://github.com/jschwab/mesa-init>`__.

Compile MESA
^^^^^^^^^^^^

Now we are ready to compile the code. This will take a little while, so
do something else for a bit or get up and get a cup of coffee.

::

   cd $MESA_DIR
   ./install

.. warning::

   There is no reason to use ``sudo``. The MESA install does not
   require root privileges.


Once it is done, you should receive the message

::

   ************************************************
   ************************************************
   ************************************************

   MESA installation was successful

   ************************************************
   ************************************************
   ************************************************

If so, you can learn more about MESA by looking at other pages.

Read the linked page that summarizes some :ref:`best practices <using_mesa/best_practices:Best practices>`
to keep in mind throughout the lifecycle of your project.

Troubleshooting
---------------

First, confirm that you can reproduce the error. Do

::

   cd $MESA_DIR
   ./clean
   ./install

and see if you get the same error.

Check that your environment variables are set correctly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One of the most common issues is unset or incorrectly set environment
variables. In the same terminal window where you are trying to install
MESA, execute the command::

    echo $MESA_DIR


and if you're using the MESA SDK, execute the command::

    echo $MESASDK_ROOT

Confirm that
these showed the directories where you have installed MESA and the MESA
SDK. If they did not, please re-read the instructions on how to :ref:`environment`.

Confirm that you installed the MESA SDK correctly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Please check that you followed the MESA SDK `installation
instructions <http://user.astro.wisc.edu/~townsend/static.php?ref=mesasdk>`__.
Pay particular attention to the prerequisites for your system.

Consult the FAQ
^^^^^^^^^^^^^^^

Check to see if there is any information about your problem in the
:ref:`MESA FAQ <faq:FAQ>`.

If you are using the MESA SDK and are having a problem with
installation, you should also consult the `MESA SDK
FAQ <http://user.astro.wisc.edu/~townsend/static.php?ref=mesasdk#Frequently_Asked_Questions_.01FAQ.01>`__.

Search the mesa-users mailing list archive
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Search the `mailing list
archives <https://lists.mesastar.org/pipermail/mesa-users/>`__ to see if
someone has had a similar problem in the past.

Post a question to mesa-users
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the previous steps have not solved your problem, send an email
message to mesa-users@lists.mesastar.org describing the problem.

Please provide the following information:

-  What version of MESA are you trying to build?

-  Are you using the MESA SDK? If so, what version?

-  Describe your computer (machine type, operating system, operating
   system version).

-  What is the error message you received?

-  Attach the ``$MESA_DIR/build.log`` file.  This includes the output of the build process along with the output of each of the following commands ::

    uname -a
    gfortran -v
    echo $MESASDK_ROOT
    echo $PATH
    echo $MESA_DIR

