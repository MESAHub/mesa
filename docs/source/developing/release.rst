Release checklist
=================

This is a guide to what needs to be done before a release can be made.

General steps
-------------

- Pick a release date 
- Pick a candidate release version
- Pick someone to be a release manager



Source code updates
-------------------

- Create :file:`data/version_number`

MESA release versions will have a human-readable version number
instead of the auto-generated git commit id.  Immediately before
release, this file should be created and commited.  Because this file
is normally ignored (via ``.gitignore``), it must be explicitly added
via ``git add -f``.

- Update version number in :file:`docs/source/conf.py`

The variables ``release`` and ``version`` in the ``Project
information`` section should be set.


Documentation
-------------

- The Changelog should be updated

.. note::
    At a minimum this should mention options that are removed/replaced and how to convert from a previous version to the newest version.

- A release notes document should be written


Testing
-------


- TestHub should report all tests pass for both Linux and macOS on multiple machines and with different OS versions
- The previous SDK version should be tested.

.. note::
    If the previous SDK does not pass we can decide whether to bump the minimum SDK version or fix the issues.

- A non-SDK machine should test the test_suite.
- At least one Windows machine should get tested.
- Recalibrate test suite cases (things like simplex_solar_calibration and example_astero)


Additional Testing
------------------

Additional checks that are not essential but should be done if there is time.

- Check the Brunt is smooth.
- Check test_memory runs and reports no memory leaks.
- Check MESA compiles with SHARED_LIBS=True
- Run with FPE checking on.
- Distribute the candidate version to several beta testers for science verification.



Release steps
-------------

- Upload release to Zenodo
- Send email to mesa-users




