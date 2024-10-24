============
Build system
============

*This page documents the MESA build system. It is only relevant if you want to modify the build system, or use some of the MESA modules outside of MESA. For building MESA itself without modifications to the build system, see* :ref:`installation:Installing MESA`.

The MESA build system uses ``make``, with some additional (shell) scripts. All relevant files can be found in the ``make`` subfolder, except for the main ``Makefile`` in the repository root, and the module specific make files (see :ref:`developing/build-system:Module build configuration`).

Running ``make`` in the root folder will build all the modules, install relevant files, and run all the associated module tests (test suite tests will not be run). To build an individual module and all its dependencies, run ``make MODULE_NAME``. Append ``-install`` or ``-test`` to run those steps. If you don't want to try to build the dependencies, ``cd`` to the folder of that module and run ``make`` there.

Build options
=============

The build system exposes several options to control the build process. If you want to change these options, give them as arguments to make, e.g. ``make BUILD_DIR=my_build_dir``.

* **BUILD_DIR**: Where all compiled files will end up, default is ``build``. The output from each module's build will be in a subfolder with the module's name.
* **COMPILER**: Which compiler to use. Supported options are ``gfortran`` (default).
* **PROFILE**: Which set of compiler flags to use, mostly to control optimizations. Supported options are ``release`` (enable optimizations), ``release-with-dbg-info`` (enable optimizations, but keep debug info, default), and ``debug`` (disable most optimizations, keep debug info).
* **WITH_OPENMP**: Whether to use openmp. Defaults to ``yes`` (any other value will disable openmp).
* **WITH_CRLIBM**: Whether to use crlibm for certain math operations. Defaults to ``yes`` (any other value will disable crlibm).

Module build configuration
==========================

Each MESA module has its own configuration Makefile. This make file is located in the root directory of the module. For example, inn the case of the ``eos`` module, this make file is ``eos/Makefile``. This make file includes configurations for how to build, test, and install the module. This can be both in the form of make variables and rules.

The following properties control the main compilation of the module:

* **MODULE_NAME**: The name of the module. (e.g. ``eos``)
* **SRCS**: Source files (.f, .f90) that are part of this module. Include files (.inc, .dek) should not be included in this list.
* **INTERNAL_DEPENDS_ON**: Which MESA modules this module directly depends on. Transitive dependencies do not need to be filled in here.
* **EXTERNAL_DEPENDS_ON**: Which external packages (e.g. lapack) this module depends on.
* **BINTYPE**: What the eventual output of the build should be, can be ``static-lib``, ``dynamic-lib`` or ``executable``.
* **INCLUDE_DIRS**: In which directories are the source files found, prefixed by ``-I``. By default, this is ``-Ipublic -Iprivate``.

The following properties control building and running the tests of the module:

* **SRCS_CHECK**: Source files that are part of the test binary. Typically, all the .f or .f90 files under ``module-name/test/src``.
* **CHECK_RESULTS_GOLDEN**: File to compare the output from the test binary against. If this is left empty, no tests will be run.
* **CHECK_DIFF_PROG**: For setting a different a program to compare the output to the golden output (e.g. ndiff). By default this is `diff -b`.

Finally, we have the install settings, which control what files are made available to other modules:

* **MODULES**: Which mod files should be made available.
* **INSTALL_INCLUDES**: Which additional include files should be made available.

Conditional compilation
-----------------------

Sometimes you want to have certain files included in the compilation depending on other variables. This can be done with the usual make file syntax. For example, see the following snippet from ``math/Makefile`` which adds a source file and dependency based on whether the variable ``NOCRLIBM`` is set.

.. code-block:: make

    ifeq ($(WITH_CRLIBM),yes)
    SRCS += public/math_lib_crmath.f90
    EXTERNAL_DEPENDS_ON += crlibm-fortran
    else
    SRCS += public/math_lib_intrinsic.f90
    endif


All variables mentioned above can be modified in this way.

Preprocessing
-------------

Some of the modules require some preprocessing of the FORTRAN source code. This requires setting up a make target for these files, and add the generated files to the ``SRCS_GENERATED`` or ``SRCS_CHECK_GENERATED``. These source files will be automatically prefixed with the path to the build directory, so make sure that the output files from the preprocessing are written to the build directory. The following snippet shows an example from the ``mtx`` module:

.. code-block:: make

    SRCS_GENERATED := private/my_lapack95_dble.f90 \
                      private/my_lapack95_quad.f90
    SRCS_CHECK_GENERATED := \
            test/src/test_block_tridiagonal_dble.f90 \
            test/src/test_block_tridiagonal_quad.f90
    
    include ../make/Makefile
    
    # Custom build steps (this needs to come after the include statement,
    #  otherwise the BUILD_DIR_MODULE variable is not set)
    
    $(BUILD_DIR_MODULE)/private/my_lapack95_dble.f90: private/my_lapack95.F90 | $(BUILD_DIR_MODULE)/private/
        # Note: PREPROCESS just calls the C preprocessor, and is set in compile-settings-*.mk
    	$(PREPROCESS) -DDBLE $^ > $@
    
    $(BUILD_DIR_MODULE)/private/my_lapack95_quad.f90: private/my_lapack95.F90 | $(BUILD_DIR_MODULE)/private/
    	$(PREPROCESS) $^ > $@
    
    $(BUILD_DIR_MODULE)/test/src/test_block_tridiagonal_dble.f90: test/src/test_block_tridiagonal.f90 | $(BUILD_DIR_MODULE)/test/src/
    	$(PREPROCESS) -DDBLE $^ > $@
    
    $(BUILD_DIR_MODULE)/test/src/test_block_tridiagonal_quad.f90: test/src/test_block_tridiagonal.f90 | $(BUILD_DIR_MODULE)/test/src/
    	$(PREPROCESS) $^ > $@

Integrating modules in another piece of software
================================================

Since MESA is written in a modular way, it is possible to take some of the piece of the MESA source code and integrate them in your own project. You will need at least the ``make`` folder which contains all the build files, the main Makefile in the root of the repository, the source of the module(s) you want to include, and the sources of all the dependent modules. If these modules need data files (e.g. eos tables), you will need to ensure you have them available as well. It is recommended to keep the same source tree, otherwise certain paths in the various make files may no longer point to the right directories. One change that will be necessary, is to edit the list of all modules in ``make/subdirs.mk``. This is used by the make files to establish the right dependencies. However, if some of the folders in that list do not exist, you will get errors. In order to build and test the modules, run ``make -C MESA_SOURCE_TREE MODULE1 MODULE2 ...``, where ``MESA_SOURCE_TREE`` refers where in your project you have copied the MESA source code. If you want to consolidate build directories of your project and the MESA modules, add ``BUILD_DIR=LOCATION_OF_BUILD_DIR`` where ``LOCATION_OF_BUILD_DIR`` is relative to the MESA source directory.

The nitty-gritty details
========================

When running ``make`` on a clean MESA source tree, the following happens:

#. The main ``Makefile`` gets loaded, which on its turn will load all relevant dependent make files in ``make/``. These additional make files contain instructions for setting up the build dir, getting the version of the source code, and collecting dependencies between the modules.
#. For each module listed in ``make/subdirs.mk``, the main make process will call into the module specific make file (``MODULE/Makefile``) and request the dependent modules of that module. A make file containing dependencies between all the modules is then generated by the ``make/gen-folder-deps`` script.
#. The version number is written to ``data/version_number`` if it does not exist with ``git describe``. As source distributions will not contain any git information, this file should be precreated.
#. For each module that needs building, make will run their builds and tests in parallel, obeying the dependencies set in step 2.

Within each module's separate make process, the following steps will happen:

#. The main make file for each module will set the relevant settings and load ``make/Makefile``, which contains the core build logic.
#. Some initialization is repeated, such as the setup of the build directory to support building modules separately.
#. Load relevant other make files, such as those for setting the correct dependencies between source files, loading compiler settings, link commands, ...
#. ``PKG_CONFIG_PATH`` is set to all known module's pkg-config build folders (``build/module_name/lib/pkg-config``). This path is generated by the ``make/gen-pkgconfig-path`` script (since it is much cleaner to do this in bash than in the make files themselves).
#. ``make/gen-compile-tree`` is called with all FORTRAN source files for that module, and is given all the necessary include paths. It will parse (in a simplified way) the source files and generate a make file that contains the necessary dependencies between all the modules. Meaning that if a module ``abc`` contains ``use def``, the compilation of module ``abc`` will need to wait until the compilation of ``def`` is done. The resulting make file is then loaded by make.
#. The source files get compiled, and the final object file gets linked according to which type of object is selected in the module config.

After the module is built, all relevant files are installed and the module-specific tests will run:

#. Copy the public modules, include files, and compiled library to ``include`` and ``lib``
#. Compile the test binary from the test source files and link it with the object file generated previously.
#. Run the test binary and compare the output to a golden output specified by the module config.

Parallel compilation
--------------------

In order to improve parallelisation, the compilation of each source file is split into two steps: generating the ``.mod`` files, and actually compiling the source code. FORTRAN ``.mod`` files describe the interface of a certain module, but do not contain the actual compiled code. In order to build dependent modules, only this ``.mod`` file is necessary. Since this file is generated much faster than compiling the module, the build system generates them in a separate step, which allows the compilation of the dependent module to start almost immediately.
