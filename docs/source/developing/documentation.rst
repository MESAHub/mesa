=============
Documentation
=============

The main MESA documentation lives in ``$MESA_DIR/docs`` and is generated using `Sphinx <https://www.sphinx-doc.org/en/master/>`__.

Once you have `installed Sphinx <https://www.sphinx-doc.org/en/master/usage/installation.html>`__, 
you can generate the docs by doing

.. code-block:: console

    cd $MESA_DIR/docs
    make html

You can then view the docs in your browser by visiting the URL ``file://<MESA_DIR>/docs/build/html/index.html`` (replace <MESA_DIR> with the appropriate path).

File Locations
==============

The main body of the docs lives in ``$MESA_DIR/docs/source``.  However, some documentation lives outside of the ``docs`` subdirectory so that it can be close to the code that it documents.

That documentation can then be incorporated into the main docs either by symlinking the file into the docs tree in an appropriate location or by using ``..include:: <filename>`` statements.

Examples of documentation that lives elsewhere include:

* The defaults files live in ``<module>/defaults/`` are symlinked into ``docs`` and also pre-processed  (see :ref:`reference/format:Format for MESA defaults files`).
* The ``README.rst`` file in each test suite case (which is included in ``docs/test_suite``).


Style
=====

This describes the style conventions used in the MESA documentation.

Here is a link to a `reStructured Text primer`_.

.. _reStructured Text primer: http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html

Click "View page source" in the upper right to see the rst file corresponding to this page.

The :ref:`reference/format:Format for MESA defaults files` page has additional
information on ReST formatting from the perspective of writing
defaults files.


Cross-references
----------------

These docs use the `autosectionlabel
<https://www.sphinx-doc.org/en/master/usage/extensions/autosectionlabel.html>`__
Sphinx extension along with the settings (see ``conf.py``):

.. code:: python

    autosectionlabel_prefix_document = True
    autosectionlabel_maxdepth = 3

This means that you can refer to sections by their title.  However,
sections are disambiguated by their location in the docs.  For
example, :ref:`this section <developing/documentation:Cross-references>`
has the label ``developing/documentation:Cross-references``.


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
-----------------

Instrument papers should be indicated using the abbreviated form |MESA I|, etc.   (This is achieved using the substitution \|MESA I\|.)

