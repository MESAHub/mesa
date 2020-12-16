Overview of kap module
======================

The opacity values returned by the ``kap`` module combine opacity
data from many sources.  The most important opacity-related MESA
options select which opacity sources to use and control the location
of the blends between them.  The total opacity is the appropriately
combined radiative opacity :math:`\kappa_{\rm rad}` and the conductive
opacity :math:`\kappa_{\rm cond}`:

.. math::

    \frac{1}{\kappa} = \frac{1}{\kappa_{\rm rad}} + \frac{1}{\kappa_{\rm cond}} ~.


Radiative Opacities
-------------------

The radiative opacity is the Rosseland mean opacity.  The opacity
depends on the temperature, density, and composition.  A set of
opacity tables consists of a collection of individual opacity tables,
each at a different composition.  Each individual table is 
tabulated in :math:`\log T` and in :math:`\log R \equiv \log \rho - 3
\log T + 18` (cgs).

For tables that assume a fixed metal distribution, the composition can
be encoded via the variables :math:`(X, Z)`.  Tables assume a
particular metal abundance pattern (usually scaled solar).  Some opacity
table sets include additional parameters in order to allow the metal
abundance pattern (often the CNO abundances) to vary.  Such
variation naturally occurs in the stellar core during helium burning
(and beyond) and in the stellar envelope as a result of dredge-up
processes.


.. note::
   
   The value of the option :ref:`Zbase` provides the reference
   metallicity necessary to calculate element variations (e.g., carbon
   and oxygen enhancement) from the composition of a cell.  The
   default opacity configuration requires this value to be specified.
   Physically, this usually corresponds to the initial metallicity of
   the star.


In MESA, separate opacity table sets are used for high and low
temperature.  In the intermediate region, both opacities are evaluated
and blended.  The location of this blend is controlled with the
options :ref:`kap_blend_logT_upper_bdy` and
:ref:`kap_blend_logT_lower_bdy`.



------------------------------------------------
High temperature :math:`(T \gtrsim 10^4\,\rm K)`
------------------------------------------------

The OPAL tables (|OPAL|) with fixed metal distributions are called
Type 1 and cover the region :math:`0.0 \leq X \leq 1-Z` and
:math:`0.0\leq Z \leq 0.1`. Type 1 tables from the Opacity Project
(OP; |OP|) are also available.  The set of tables to be used are
selected by the option :ref:`kap_file_prefix`.

Additionally, there is support for the OPAL Type 2 tables that allow
for varying amounts of C and O beyond that accounted for by :math:`Z`;
these are needed during helium burning and beyond. These have a range
:math:`0.0 \leq X \leq 0.7`, :math:`0.0\leq Z\leq0.1`.  The set of
tables to be used are selected by the option :ref:`kap_CO_prefix`.

Type 2 tables on by default (see :ref:`use_Type2_opacities`) and The
blends between these table sets occur based on hydrogen fraction
(see :ref:`kap_Type2_full_off_X` and :ref:`kap_Type2_full_on_X` and )
and metal enhancement (controlled by :ref:`kap_Type2_full_off_dZ` and
:ref:`kap_Type2_full_on_dZ`).


-------------------------------------------------
Low temperature  :math:`(T \lesssim 10^4\,\rm K)`
-------------------------------------------------

Low temperature opacities are selected with the option
:ref:`kap_lowT_prefix`.

Tables based on the work of |Fergusson| include the effects of
molecules and grains and cover the range
:math:`2.7 \le \log T \le 4.5` and :math:`-8 \le \log R \le 1`.

Tables based on the work of |Freedman| include the effects of
molecules and cover the range :math:`1.88 \le \log T \le 4.5` and
:math:`-8 \le \log R \le 9`.  The table set was privately communicated
by R. S. Freedman in 2011.  Unlike other opacity sources, this is a 1D
sequence of tables in :math:`Z` as opposed to a 2D grid of
:math:`(X,Z)` values.  (The assumed H/He abundances scale with
:math:`Z`.)


Tables from ÆSOPUS (|AESOPUS|) include variation factors for the CNO
isotopes.  The opacity is evaluated using the global value of
:math:`Z_{\rm base}` and the local (cell) values of :math:`(X, X_{\rm
C}, X_{\rm N}, X_{\rm O})`.

The ÆSOPUS tables are provided at a set of reference metalicites.  In
order to interpolate to the provided :math:`Z_{\rm base}`, the opacity
is evaluated at an appropriate subset of these reference values (and
then interpolated).  For each such :math:`Z_{\rm ref}`, the ÆSOPUS
composition parameters

.. math::

   \begin{eqnarray*}
   f_{\rm CO} = \log(X_{\rm C}/X_{\rm O}) - f_{\rm CO, ref} \\
   f_{\rm C} = \log(X_{\rm C}/Z_{\rm ref}) - f_{\rm C, ref} \\
   f_{\rm N} = \log(X_{\rm N}/Z_{\rm ref}) - f_{\rm N, ref} \\
   \end{eqnarray*}

are calculated, the opacities evaluated the tables with bracketing
compositions, and the resulting opacities linearly interpolated.
(Note that this means that the interpolation in :math:`Z` occurs at
fixed :math:`X` and :math:`f_{\rm CO}`, but not at fixed :math:`f_{\rm
C}` or :math:`f_{\rm N}`.)
   
------------------
Compton Scattering
------------------

At sufficiently high temperature :math:`(T \gtrsim 10^8\,\rm K)`, the
opacity will be dominated by Compton scattering.  MESA calculates the
opacity of Compton scattering using the equations of |BY76|.  Near the
high-:math:`T` and low-:math:`R` edges of the high temperature opacity
tables, MESA smoothly blends the tabulated opacity values with the
Compton scattering values.  The location of these blends is not
user-controllable.

Conductive Opacities
--------------------

The conductive opacity :math:`(\kappa_{\rm cond})` is given by the
thermal conductivity :math:`(K)` appropriately recast such that the heat
transfer equation resembles the form of the equation used in radiative
diffusion (e.g., HKT Section 4.5).  This implies

.. math::

   \kappa_{\rm cond} = \frac{16 \sigma_{\rm SB} T^3}{\rho K} ~.

The thermal conductivities used in MESA are an extended version of the
results of |Cassisi| privately communicated by A.Y. Potekhin.  They
are tabulated for a set of :math:`1 \le \bar{Z} \le 60`.  Each table
spans :math:`-6 \le \log(\rho/\rm g\,cm^{-3}) \le 11.50` and :math:`3
\le \log(T/\rm K) \le 10`.


.. |BY76| replace:: `Buchler & Yueh (1976) <https://ui.adsabs.harvard.edu/abs/1976ApJ...210..440B/abstract>`__

.. |Fergusson| replace:: `Ferguson et al. (2005) <https://ui.adsabs.harvard.edu/abs/2005ApJ...623..585F/abstract>`__

.. |Freedman| replace:: `Freedman et al. (2008) <https://ui.adsabs.harvard.edu/abs/2008ApJS..174..504F/abstract>`__

.. |AESOPUS| replace:: `Marigo & Aringer 2009 <https://ui.adsabs.harvard.edu/abs/2009A%26A...508.1539M/abstract>`__

.. |OPAL| replace:: Iglesias & Rogers `1993 <https://ui.adsabs.harvard.edu/abs/1993ApJ...412..752I/abstract>`__, `1996 <https://ui.adsabs.harvard.edu/abs/1996ApJ...464..943I/abstract>`__


.. |OP| replace:: `Seaton 2005 <https://ui.adsabs.harvard.edu/abs/2005MNRAS.362L...1S/abstract>`__                    

.. |Cassisi| replace:: `Cassisi et al. (2007) <https://ui.adsabs.harvard.edu/abs/2007ApJ...661.1094C/abstract>`__

