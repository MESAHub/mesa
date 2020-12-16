Installing MESA on Windows
==========================

This page describes two options for installing MESA on Windows.

Option 1: MESA-Docker
---------------------
MESA-Docker provides a pre built version of MESA inside a Docker container. This simplifies the setup as MESA and all 
its dependencies are already installed.

`MESA-Docker <https://github.com/evbauer/MESA-Docker>`__

Option 2: WSL2
--------------

WSL2 enables installing Linux inside your Windows operating system.

.. note::

    You must install version 2 of WSL and not version 1. WSL1 does not work with MESA.

`Follow the instructions from Microsoft for installing WSL2 <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`__

.. note::

    We recommend you choose to install Ubuntu as your Linux distribution.

Once this is set up, follow the instructions for :ref:`Installing MESA` on Linux.


Enabling pgplot
^^^^^^^^^^^^^^^

To enable pgplot to work you can install `VcXsrv <https://sourceforge.net/projects/vcxsrv>`__.
You must have ``VcXsrv`` running before you start a Linux terminal for it to work.

.. note::

    Sometimes VcXsrv is flaky and does not work well. This does not impact your science,
    only whether a pgplot window appears. If VcXsrv does not work you can always set your
    pgstar inlists to save plots instead of displaying them. Saving works whether you have a 
    working VcXsrv or not.

Once installed launch ``VcXsrv``. Accept the default choices except that you should select the option
``Disable access control``. The first time you run VcXsrv you may get a firewall prompt, in which 
case you should allow it access to public and private networks.







