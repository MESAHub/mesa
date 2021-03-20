****
FAQ
****

General Questions
=================

How do you say MESA?
--------------------

Say "MAY-sa".

Can I use the MESA logo in a talk or presentation?
--------------------------------------------------

Yes. It is available in a variety of sizes for your convenience.

-  :download:`30 pt </assets/logo/mesa_logo2_30pt.png>`
-  :download:`60 pt </assets/logo/mesa_logo2_60pt.png>`
-  :download:`100 pt </assets/logo/mesa_logo2_100pt.png>`
-  :download:`144 pt </assets/logo/mesa_logo2_144pt.png>`
-  :download:`200 pt </assets/logo/mesa_logo2_200pt.png>`
-  :download:`400 pt </assets/logo/mesa_logo2_400pt.png>`

What happens when an instrument paper is being written?
-------------------------------------------------------

The picture shows the intense authors hard at work making final edits.

.. figure:: /assets/images/paper_session.jpg
   :alt: Writing the mesa paper

A late night session working on |MESA I| around Bill's
dining room table with a projector, a makeshift screen, and ample red
wine (left to right, Frank Timmes, Aaron Dotter, Falk Herwig, Lars
Bildsten, and Bill Paxton). Photo taken by Bill's patient wife, Kathlyn,
who deserves a great deal of credit for the existence of MESA.

By |MESA V|, there's a bigger room and a better screen.
Left-to-right: Aaron Dotter, Jared Goldberg, Lars Bildsten, Radek
Smolec, Anne Thoul, Bill Paxton, Adam Jermyn, Josiah Schwab, Frank
Timmes, Croissants, Rich Townsend, Coffee, and Pablo Marchant. Photo
taken by Matteo Cantiello.

|MESAV-Photo|

.. |MESAV-Photo| image:: /assets/images/paper_v_session_small.jpg
   :target: /assets/images/paper_v_session.jpg


Why did my inlists and/or run_star_extras stop working?
-------------------------------------------------------

As MESA evolves, inlist options, variable names, and other aspects of
the code change in ways that are not backward-compatible. Here is a list
of some of the major renaming events; it is by no means exhaustive.
Consult the release notes for these versions for more information.

-  `r4411 <release/2012/08/25/r4411.html>`__: major renaming and
   reorganization of the code
-  `r8118 <release/2016/02/03/r8118.html>`__: renaming and restructuring
   of controls related to winds. ``run_star_extras`` routines changed to
   hooks.
-  `r8845 <release/2016/06/19/r8845.html>`__: most names using "cnt"
   were changed to use "interval" (e.g., ``terminal_cnt`` was renamed to
   ``terminal_interval``). the option ``photostep`` was renamed to
   ``photo_interval``. One exception is ``Abundance_legend_max_cnt``
   which remained unchanged since it is not an interval between events.

Prereqs & Installation
======================

Can MESA be run on 32-bit machines?
-----------------------------------

MESA is no longer supported on 32-bit systems. This isn't by design; we
simply haven't been able to get it to work! So if you REALLY want to use
an antique 32-bit system for MESA, you'll have to try to make it work
yourself.

Which SVN revisions were MESA release versions?
-----------------------------------------------

{% include releases.html %}

Why is my MESA zipfile download taking forever?
-----------------------------------------------

The MESA zipfiles are large (~1GB), so if you're on a slow connection,
there's nothing to be done. However, users have reported substantial
variations in download speed from the sourceforge site. If you are
experiencing slow download speeds, you may want to change which mirror
you are using. (After you select a zipfile to download, click the "try
another mirror" link near the top of the page.)

Why am I getting "wrong module version" errors?
-----------------------------------------------

Errors like::

   Fatal Error: Wrong module version '7' (expected '5')

typically arise because you have changed the gfortran compiler since you
last built MESA. (Sometimes a changed gfortran version is an indication
that you forgot to activate the MESA SDK.) To fix, run ./clean in the
MESA directory, and then try building again.

What do I need to do to run MESA on OS X?
-----------------------------------------

It is necessary to take the following actions on recent versions of OS X
(10.9 or later) , even if you had a working version of MESA on your
previous version of OS X.

-  Install or reinstall the current version of
   `Xquartz <http://xquartz.macosforge.org/landing/>`__.

-  Install the `command line
   tools <https://developer.apple.com/library/ios/technotes/tn2339/_index.html#//apple_ref/doc/uid/DTS40014588-CH1-WHAT_IS_THE_COMMAND_LINE_TOOLS_PACKAGE_>`__,
   using the command::

    xcode-select -install

-  If you are using macOS 10.14 (Mojave), there has been a change in how
   the system headers are provided (see section "Command Line Tools" in
   the `Xcode 10 release
   notes <https://developer.apple.com/documentation/xcode_release_notes/xcode_10_release_notes>`__).
   As a workaround, Apple provides an extra package that will install
   the headers to the base system. At present, one must do this:::

    open /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg

-  Install the latest version of the `MESA
   SDK <http://www.astro.wisc.edu/~townsend/static.php?ref=mesasdk#Mac_OS_X>`__
   appropriate for your version of OS X.

Why do I get the error "stdlib.h: No such file or directory"?
-------------------------------------------------------------

If yore using macOS 10.14 (Mojave), then this error indicates that you
did not install the system headers (see previous FAQ or `the SDK
website <http://www.astro.wisc.edu/~townsend/static.php?ref=mesasdk#Mac_OS_X>`__).

Why do I get the error "C compiler cannot create executables"?
--------------------------------------------------------------

If you receive the error

::

   configure: error: C compiler cannot create executables
   See `config.log' for more details.

during the MESA installation of crlibm, please check that you are not
using an unsupported 32-bit system. If you're using macOS, make sure
you've `installed Xcode and the command line tools <#osx>`__.

If you're using Linux and you see an error message in config.log like

::

   In file included from /opt/mesasdk/lib/gcc/x86_64-pc-linux-gnu/5.3.1/include-fixed/features.h:338:0,
                    from /usr/include/bits/libc-header-start.h:33,
                    from /usr/include/stdio.h:28,
                    from test.c:1:
   /usr/include/sys/cdefs.h:481:49: error: missing binary operator before token "("
    #if __GNUC_PREREQ (4,8) || __glibc_clang_prereq (3,5)
                                                    ^
   In file included from test.c:1:0:
   /usr/include/stdio.h:320:43: error: missing binary operator before token "("
    #if defined __USE_XOPEN2K8 || __GLIBC_USE (LIB_EXT2)
                                              ^
   /usr/include/stdio.h:399:17: error: missing binary operator before token "("
    #if __GLIBC_USE (LIB_EXT2)
                    ^
   /usr/include/stdio.h:657:43: error: missing binary operator before token "("
    #if defined __USE_XOPEN2K8 || __GLIBC_USE (LIB_EXT2)

you may need to `tell gcc to regenerate its fixed header
files <https://gcc.gnu.org/onlinedocs/gcc-7.2.0/gcc/Fixed-Headers.html#Fixed-Headers>`__.

You can do so with the commands:

::

   GCC_VERSION=`gcc --version | grep ^gcc | sed 's/^.* //g'`
   $MESASDK_ROOT/libexec/gcc/x86_64-pc-linux-gnu/$GCC_VERSION/install-tools/mkheaders $MESASDK_ROOT

If the error persists, please email mesa-users@lists.mesastar.org.
Follow `the instructions for posting a question to
mesa-users <prereqs.html#post-a-question-to-mesa-users>`__ and also
attach the file ``$MESA_DIR/crlibm/crlibm/config.log``.

Why do I get the error "Cannot utime: Function not implemented"?
----------------------------------------------------------------

This signifies that you are using WSL1 and not the newer WSL2. `Please update your WSL
to the newest version <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`__


Why do I get an error "0xc03a001a"?
-----------------------------------

You need to disable disk compression on the folder that stores the Linux distribution you
are using with WSL2.

Using file explorer navigate to the folder ``C:\Users\USERNAME\AppData\Local\Packages\``,
where ``USERNAME`` is your windows username. If you installed Ubuntu, then you need to find the folder
``CanonicalGroupLimited.Ubuntu``, other distributions will have different names. Once you find this folder
right click on its name to open it properties, go to ``Advanced`` and unclick ``Compress contents to save disk space`` 

`More information can be found here <https://github.com/microsoft/WSL/issues/4299>`__


Why do I get an error "#000: H5F.c line 444 in H5Fcreate(): unable to create file"?
-----------------------------------------------------------------------------------

You need to turn off HDF5 file locking:

::

   export HDF5_USE_FILE_LOCKING='FALSE'


Why do I get an error like "is smaller than expected for a file tracked by git LFS"?
------------------------------------------------------------------------------------

First check that you have git lfs installed and then try running:

::

   git lfs pull

If you are checking out MESA on a network file system you may find setting a short delay duirng the
install step to give time for the files transfered by git lfs to properly sync.

::

   export MESA_GIT_LFS_SLEEP=10



Installing Older Versions of MESA
=================================

Older versions of MESA may fail to compile with more recent versions of
the MESA SDK. There are a few possible solutions:

(1) Use the contemporary version of the SDK. There is a `list of old
    MESA
    SDKs <http://www.astro.wisc.edu/~townsend/static.php?ref=mesasdk-old>`__
    that indicates which version of MESA was current when each was
    released.

(2) Work in a Docker container configured to set up an environment
    appropriate for the desired MESA version. NuGrid's
    `NuDocker project <https://github.com/NuGrid/NuDocker>`__
    provides a convenient interface for MESA versions as far back
    as 4942.
    
(3) Fix the individual compilation errors (there are typically only a
    handful to work through). The basic workflow is

-  ./install until an error turns up,
-  cd to that module and edit and ./mk until it compiles, then
-  go back to ./install and repeat until everything works before
-  a final ./touch and ./install.

The following FAQs give examples of specific errors.

Why do I get an error like "Error: 'time0' may be used uninitialized in this function"?
---------------------------------------------------------------------------------------

This sort of error typically occurs when using an older MESA release
with a newer SDK. These warnings (which are treated as errors) can
safely be ignored. To do so, edit the file ``utils/makefile_header`` and
add ``-Wno-uninitialized`` to ``FCbasic``.

Why do I get an error like "Error: Blank required in STOP statement near (1)"?
------------------------------------------------------------------------------

This sort of error typically occurs when using an older MESA release
with a newer SDK. To work around this, simply insert the blank space as
requested. i.e. change ``stop'fixup'`` to ``stop 'fixup'``.

Why do I get an error like "Error: Array reference at (1) out of bounds (0 < 1) in loop beginning at (2) [-Werror=do-subscript]"?
---------------------------------------------------------------------------------------------------------------------------------

This sort of error typically occurs when using an older MESA release
with a newer SDK. These warnings (which are treated as errors) can
safely be ignored. To do so, edit the file ``utils/makefile_header`` and
delete ``-Werror`` from ``FCwarn``.

Why do I get an error like "Error: '__builtin_memset' specified size between 9223372036854775808 and 18446744073709551615 exceeds maximum object size 9223372036854775807 [-Werror=stringop-overflow=]"?
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

This sort of error typically occurs when using an older MESA release
with a newer SDK. These warnings (which are treated as errors) can
safely be ignored. To do so, edit the file ``utils/makefile_header`` and
delete ``-Werror`` from ``FCwarn``.

Why do I get a segfault in do_history_info?
-------------------------------------------

This sort of error typically occurs when using an older MESA release
(r10398 or before) with a newer SDK (that includes gfortran 7.3 or
later).

::

   Program received signal SIGSEGV: Segmentation fault - invalid memory reference. 

   Backtrace for this error: 
   #0  0x7f28c0a93a7f in ??? 
   #1  0x54313f in do_history_info 
          at ../private/history.f90:383 
   ...
      

To work around this, edit ``$MESA_DIR/star/private/history.f90`` and
replace the line

.. code-block:: fortran

  if (open_close_log) close(io)

.. code-block:: fortran

  if (open_close_log .and. write_flag) close(io)

Note that this only applies to this specific segfault and not segfaults
in general.

When I install MESA from a zipfile, why do I get svn errors?
------------------------------------------------------------

Older versions of MESA (r6794 or earlier) assumed that they would be
installed from the subversion repository. Therefore, if you install
these versions from a zipfile, you will see messages like

::

   svn: E155007: '/Users/fxt/mesa/mesa-r6794' is not a working copy

when you run MESA. These messages are safe to ignore.

PGPLOT and pgstar
=================

Why don't I see any PGPLOT output when I run MESA?
--------------------------------------------------

Make sure you have the option pgstar_flag = .TRUE. in the &star_job
section of your input file. Also, if you're not using the SDK, make sure
you enabled PGPLOT in your utils/makefile_header file (when using the
SDK, PGPLOT is enabled by default).

How can I make a movie from my pgstar output?
---------------------------------------------

As of the 20140713 release, the MESA SDK includes the ffmpeg encoder and
a simple script, images_to_movie.sh, which uses ffmpeg to create movies
from PNG files produced by MESA.

To make use of this capability, consult `Rich's
instructions <http://www.astro.wisc.edu/~townsend/static.php?ref=mesasdk#Making_Movies>`__.

Programming
===========

What programming language is MESA written in?
---------------------------------------------

Fortran. But MESA is written using advanced features of modern Fortran
which make it very different from Fortran 77 codes you might have
previously seen (or written!). If you are not already familiar with the
new and wonderful things that have been added, there are good resources
available on the web - here's one: `Fortran 95 language
features <http://en.wikipedia.org/wiki/Fortran_95_language_features>`__.

What does "thread-safe" mean?
-----------------------------

"Thread-safe" simply means that users can take advantage of multicore
processors.

For example, during stellar evolution, you need to evaluate the eos at
lots of points:

.. code-block:: fortran

  do k = 1, num_zones
    call eos(T(k), rho(k), ….)
  end do

Most fortran compilers (ifort and gfortran and others) support OpenMP,
so the loop can be made to run in parallel by adding 2 lines of magic:

.. code-block:: fortran

  !$OMP PARALLEL DO PRIVATE (k)
  do k = 1, num_zones
    call eos(T(k), rho(k), ....)
  end do
  !$OMP END PARALLEL DO

Now, if I have 4 cores, I'll have 4 threads of execution evaluating the
eos.

However, for it to work, the implementation needs to be thread-safe. In
practice this means, making shared data read-only after initialization.
Working memory must be allocated from the stack (as local variables of
routines) or allocated dynamically (using fortran95's allocate
primitive). So, basically it boils down to avoiding common blocks and
save's. It's easy to arrange for this in new code; it can be nasty
converting old code however. Stellar evolution is a good candidate for
making use of many cores. Just wrapping "parallel" directives around
some loops does it if the system is designed with that in mind - and
MESA is. But there is still much to be done to make good use of more
than 3 or 4 cores in MESA/star. It will be an ongoing effort to improve
the design of the code in that area; help with that is welcomed!

How can I use an individual MESA module?
----------------------------------------

The easiest way to get the idea is to look at a sample, and in the MESA
directory you will find a subdirectory called "sample". Make a copy of
the sample directory anywhere you'd like and give it whatever name you
want. Follow the instructions in the README file to make and test the
sample. You may also want to look at the following information from
previous MESA summer schools:

-  `Kevin Moore
   (2015) <https://mesa2015-external-code.readthedocs.io/en/latest/>`__
-  `Frank Timmes (2019) <https://zenodo.org/record/3372835>`__

Third-party Tools
=================

How can I read and/or plot MESA data using language X?
------------------------------------------------------

Users have posted numerous useful scripts in a variety of languages on
`the MESA
marketplace <http://cococubed.asu.edu/mesa_market/add-ons.html>`__.

How can I include the effects of chemical enhancements on the opacities?
------------------------------------------------------------------------

Ehsan Moravveji has developed a package available for `free access on
bitbucket <https://bitbucket.org/ehsan_moravveji/op_mono/overview>`__ to
recompute OP opacities for any desired mixture (that MESA supports) and
for any user-specified iron and nickel enhancement factors.

Other
=====

Why am I having difficulty receiving mesa-users messages?
---------------------------------------------------------

If you are using an address hosted at qq.com, this domain is blocked by
our mail host. Due to a deluge of spam, no mail will be accepted from or
delivered to qq.com.

