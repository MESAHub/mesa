Documentation style
===================

This describes the style conventions used in the MESA documentation.

Here is a link to a `reStructured Text primer`_.

.. _reStructured Text primer: http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html

Click "View page source" in the upper right to see the rst file corresponding to this page.

The :ref:`Format for MESA defaults files` page has additional
information on ReST formatting from the perspective of writing
defaults files.

Names
-----

Write module names as literals.  My favorite module is ``eos``.

Option names also literals like ``my_option``.  Ideally someday these
would become hyperlinks to their reference entries.


References
----------

References should be written in the usual astrophysics style --
like `Dotter et al. (2009) <https://ui.adsabs.harvard.edu/abs/2009A%26A...507.1617D/abstract>`__
-- and hyperlinked to either their ADS entry (preferred) or a DOI.


Symbols
-------

Common astrophysics symbols have substitutions already defined in ``conf.py``.

.. list-table::
   :widths: 15 15
   :header-rows: 1

   * - Symbol
     - Substitution
   * - |Msun|
     - \|Msun\|
   * - |Lsun|
     - \|Lsun\|
   * - |Rsun|
     - \|Rsun\|
   * - |Teff|
     - \|Teff\|
   * - |chi^2|
     - \|chi^2\|


Instrument Papers
^^^^^^^^^^^^^^^^^

Instrument papers should be indicated using the abbreviated form |MESA I|, etc.   (This is achieved using the substitution \|MESA I\|.)

