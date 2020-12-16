Debugging
=========

Debugging MESA requires effort and tenacity.  There is no single
approach that is applicable in all circumstances.  This page collects
information on frequently useful tools and techniques.

General Techniques
------------------

Navigating through source code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Tools like ``grep`` --- or more advanced variants like `ack <https://beyondgrep.com/>`_, `ag <https://github.com/ggreer/the_silver_searcher>`_, or `rg <https://github.com/BurntSushi/ripgrep>`_ --- allow you to rapidly search through the code base for occurrences of a given string or pattern.  These tools are helpful for finding where a particular function/subroutine is defined or called or where a particular variable is used.

Your text editor (or an IDE) may also provide capabilities for semantic navigation (e.g., jump to definition).

There is auto-generated documentation from Doxygen available at `<http://mesa.sourceforge.net/dox/>`_.  This provides searchable lists of files, function/subroutines, and shows call graphs.


Using a debugger
~~~~~~~~~~~~~~~~

The MESA SDK contains the GNU debugger `gdb <https://www.gnu.org/software/gdb/>`_.  Detailed instructions for using ``gdb`` are beyond the scope of the documentation, but it allows you to interactively set breakpoints, step through code, and inspect values.

Using pgstar
~~~~~~~~~~~~

The built in pgstar plotting is itself an invaluable a debugging tool.  Something may jump out from a plot as strange behavior before it shows up as a problem in the run.  You may notice that an issue begins when the model reaches certain conditions and that can provide a useful guide to the problem.  After you identify an issue and configure some useful plots, consider setting ``pgstar_interval = 1`` and ``pause = .true.`` so that you can step through the model timestep-by-timestep.

Using a profiler
~~~~~~~~~~~~~~~~

Troubleshooting performance bugs or memory leaks benefits from
specialized tools (see :doc:`profiling`).


Diagnosing Solver Struggles
---------------------------

Identifying when the solver is having difficulties and understanding
the root cause is an important subset of MESA debugging.  When you see
a MESA model experiencing problems (e.g., taking unexpectedly small
timesteps due to many retries), you should at the least report the
issue to the developers' list and ideally investigate it a bit
yourself.  Paying attention to such problems has proved a valuable way
of identifying errors and inconsistencies within MESA.

Diagnosing solver struggles primarily relies on internal MESA tooling.
The most common cause of solver struggles is incorrect/inaccurate
derivatives.  This document will illustrate how to go from a failing
model down into the depths of the code in order to identify and fix
the problem.


Step -1: Introduce a bug
~~~~~~~~~~~~~~~~~~~~~~~~

For the purpose of this tutorial, we will introduce a bug in MESA.  In
more typical circumstances, you or someone else will have kindly already
taken care of this step.

Navigate to the subroutine ``Get_kap_Results`` in ``kap/private/kap_eval.f90``.
After the call to ``Get_kap_Results_blend_T``, add

.. code-block:: fortran

     dlnkap_rad_dlnRho = -dlnkap_rad_dlnRho
         
thereby making the derivatives of the radiative opacity with respect
to density incorrect.

After you've make this change, in ``kap`` do

.. code-block:: console

    ./mk
    ./export

to recompile and export the module.


Step 0: Notice your model has problems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are many ways in which a MESA model may have problems.  It might
fail to run with gold tolerances or it might take unexpectedly small
timesteps.  The timestep is often driven downwards by repeated
retries.  Information about the retry trigger is printed to the
terminal.  Or, if the solver requires more iterations than a specified
limit (``solver_iters_timestep_limit`` and related controls), it also
reduces the timestep and indicates the ``dt_limit`` as ``solver iters``.

.. note::

    The debugging steps that follow are not generally applicable when
    the solver is performing well.  For other kinds of bugs (e.g., an
    option doesn't produce the expected behavior) one often has
    clearer idea where the problem may be and so you can head straight
    to that part of the code.

We will use the ``1.3M_ms_high_Z`` test suite case to illustrate the
effect of the bug that we introduced.

In r14460, this case took 377 steps, required 0 retries, and there
were a total of 1850 iterations of the newton solver.  With the bug,
this test case stops when it reaches ``max_model_number`` after 500
steps, having done 12 retries and a total of 8145 newton iterations.

Navigate to this test (in ``star/test_suite``) and run it.  In the
terminal output, you can see that the timesteps are often set by
``solver iters``.


.. code-block:: console

           step    lg_Tmax     Teff     lg_LH      lg_Lnuc     Mass       H_rich     H_cntr     N_cntr     Y_surf   eta_cntr   zones  retry
       lg_dt_yr    lg_Tcntr    lg_R     lg_L3a     lg_Lneu     lg_Mdot    He_core    He_cntr    O_cntr     Z_surf   gam_cntr   iters  
         age_yr    lg_Dcntr    lg_L     lg_LZ      lg_Lphoto   lg_Dsurf   C_core     C_cntr     Ne_cntr    Z_cntr   v_div_cs       dt_limit
    __________________________________________________________________________________________________________________________________________________
    
             10   7.186162   6026.290   0.310381   0.310381   1.300000   1.300000   0.639691   0.003250   0.320000  -2.349937    805      0
       4.820702   7.186162   0.119685 -44.790170  -1.003419 -99.000000   0.000000   0.320132   0.018721   0.040000   0.071396     16
     8.4252E+05   1.739919   0.314263 -14.875409 -99.000000  -7.019694   0.000000   0.005826   0.004199   0.040176  0.000E+00  solver iters
                                   rel_E_err    5.8461078509513597D-13
                           log_rel_run_E_err      -11.0597023512447645
    
             13   7.185990   6017.566   0.300802   0.300802   1.300000   1.300000   0.639670   0.003358   0.320000  -2.355541    805      0
       4.873378   7.185990   0.113850 -44.803173  -1.013273 -99.000000   0.000000   0.320139   0.018721   0.040000   0.071285     15
     1.0438E+06   1.737310   0.300077 -14.750470 -99.000000  -7.010156   0.000000   0.005733   0.004199   0.040192  0.000E+00    varcontrol
                                   rel_E_err    7.3563886572871659D-13
                           log_rel_run_E_err      -10.9707092147706380
    
    save LOGS/profile2.data for model 13
             20   7.185235   5994.060   0.274258   0.274258   1.300000   1.300000   0.639610   0.003646   0.320000  -2.362919    806      0
       4.851161   7.185235   0.101141 -44.848473  -1.043016 -99.000000   0.000000   0.320157   0.018721   0.040000   0.071188     18
     1.5910E+06   1.733084   0.267860 -15.352530 -99.000000  -6.988274   0.000000   0.005486   0.004199   0.040233  0.000E+00  solver iters
                                   rel_E_err    5.1727204553814077D-13
                           log_rel_run_E_err      -10.1566890117973099
    

Step 1: Activate debugging options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are standard inlist commands available that can reveal more about
the struggle that the solver is going through.  They are stored in the
aptly-named ``debugging_stuff_for_inlists`` which lives in
``star/test_suite`` directory.  The full list of commands there is

.. literalinclude:: ../../../star/test_suite/debugging_stuff_for_inlists

Copy the contents of ``debugging_stuff_for_inlists`` to the controls
section of your inlist.  Note that many test suite inlists already
contain these options.

Begin by uncommenting the following lines to get info about the newton solver progress::

      report_solver_progress = .true. ! set true to see info about newton iterations
      report_ierr = .true. ! if true, produce terminal output when have some internal error



Step 2: Run the model and find the bad spot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Running the model again will now produce a lot more output, giving
information about each and every call to the solver.  That output will
look something like this:

.. code-block:: console

                                     solver_call_number          17
                                  gold tolerances level           2
                       correction tolerances: norm, max    3.0000000000000001D-05    3.0000000000000001D-03
                    tol1 residual tolerances: norm, max    9.9999999999999994D-12    1.0000000000000001D-09
    17    1  coeff  1.0000     avg resid  0.518E-05      max resid dv_dt     1  0.111E-01  mix type xx000     avg corr  0.474E-03    max corr lnd   556  0.710E-02  mix type 00000     avg+max corr+resid
    17    2  coeff  1.0000     avg resid  0.561E-06      max resid dv_dt     1  0.275E-02  mix type xx000     avg corr  0.140E-03    max corr lnd   493  0.140E-02  mix type 11100     avg corr, avg+max resid
    17    3  coeff  1.0000     avg resid  0.384E-06      max resid dv_dt     1  0.521E-02  mix type xx000     avg corr  0.262E-04    max corr lnd   494  0.253E-03  mix type 11000     avg+max resid
    17    4  coeff  1.0000     avg resid  0.681E-06      max resid dv_dt     1  0.819E-02  mix type xx000     avg corr  0.560E-04    max corr lnd   494 -0.555E-03  mix type 11000     avg corr, avg+max resid
                    tol2 residual tolerances: norm, max    1.0000000000000000D-08    1.0000000000000001D-05
    17    5  coeff  1.0000     avg resid  0.104E-05      max resid dv_dt     1  0.121E-01  mix type xx000     avg corr  0.780E-04    max corr lnd   494  0.758E-03  mix type 11000     avg corr, avg+max resid
    17    6  coeff  1.0000     avg resid  0.159E-05      max resid dv_dt     1  0.177E-01  mix type xx000     avg corr  0.113E-03    max corr lnd   494 -0.110E-02  mix type 11000     avg corr, avg+max resid
    17    7  coeff  1.0000     avg resid  0.251E-05      max resid dv_dt     1  0.263E-01  mix type xx000     avg corr  0.166E-03    max corr lnd   494  0.161E-02  mix type 11000     avg corr, avg+max resid
    17    8  coeff  1.0000     avg resid  0.403E-05      max resid dv_dt     1  0.383E-01  mix type xx000     avg corr  0.246E-03    max corr lnd   494 -0.238E-02  mix type 11000     avg corr, avg+max resid
    17    9  coeff  0.5000     avg resid  0.185E-05      max resid dv_dt     1  0.936E-02  mix type xx000     avg corr  0.359E-03    max corr lnd   494  0.348E-02  mix type 11000     avg+max corr+resid
    17   10  coeff  0.5000     avg resid  0.777E-06      max resid dv_dt     1  0.227E-02  mix type xx000     avg corr  0.857E-04    max corr lnd   494 -0.827E-03  mix type 11000     avg corr, avg+max resid
    17   11  coeff  1.0000     avg resid  0.262E-06      max resid dv_dt     1  0.329E-02  mix type xx000     avg corr  0.226E-04    max corr lnd   494  0.220E-03  mix type 11000     avg+max resid
    17   12  coeff  0.5000     avg resid  0.660E-07      max resid dv_dt     1  0.754E-03  mix type xx000     avg corr  0.306E-04    max corr lnd   494 -0.296E-03  mix type 11000     avg corr, avg+max resid
    17   13  coeff  0.5000     avg resid  0.190E-07      max resid dv_dt     1  0.185E-03  mix type xx000     avg corr  0.712E-05    max corr lnd   494  0.691E-04  mix type 11000     avg+max resid
    17   14  coeff  0.5000     avg resid  0.540E-08      max resid dv_dt     1  0.409E-04  mix type xx000     avg corr  0.172E-05    max corr lnd   494 -0.166E-04  mix type 11000     avg+max resid
    17   15  coeff  0.5000     avg resid  0.200E-08      max resid dv_dt     1  0.106E-04  mix type xx000     avg corr  0.382E-06    max corr lnd   494  0.369E-05  mix type 11000     max resid
    17   16  coeff  0.5000     avg resid  0.754E-09      max resid dv_dt     1  0.216E-05  mix type xx000     avg corr  0.103E-06    max corr lnd   494 -0.100E-05  mix type 11000     okay!


The first line gives ``solver_call_number``, which uniquely identifies
a particular solver call.  (It will differ from the model number after
retries have occurred.)  Then, there is an indication of what set of
solver tolerances are being used and what those tolerances are.  After
that, there is a series of lines that gives iteration-by-iteration
summaries of what the solver did.  The iteration number is the second
column.  After certain numbers of iterations (user-defined, see e.g.,
``gold_iter_for_resid_tol2``), the solver relaxes the residual
tolerances and the new tolerances are indicated.

Other columns show information about the residuals and corrections.
The "residual" is the left over difference between the left and right
hand sides of the equation we are trying to solve. We do iterations to
reduce that, but we are limited by the non-linearity of the problem
and the quality of the estimates for the derivatives.  The
"correction" is the change in the primary variable, calculated using
Newton's method, so the Jacobian and residuals give a correction that
would make the next residual vanish if the problem were linear and the
Jacobian were exact, neither of which are true. So the best we can
hope for is that the corrections will get smaller next iteration.
When the maximum residual/correction is shown, the equation/variable
and its cell index are indicated.  The mixing types for the cells
around these locations are also indicated.  The 1's are the convection
and the 0's are for no-mixing.  (The full list of mixing types is in
``const_def``.)


The last column is the simplest to understand.  It tells us that for
the early iterations, the solution was rejected because most
everything was unsatisfactory.  Then, as the iterations proceed and
the residuals and corrections shrink, eventually only the maximum
residual is unsatisfactory.  Finally, the solver decides the solution
is good enough.

.. note::

   Less artificial bugs may behave somewhat differently than this toy
   problem.  In this example, the solver is pretty much always
   struggling, so you could pick almost any step as the bad one.  In
   practice, models often only experience problems under specific
   conditions and so a much smaller subset of steps are bad.
   Generally, you want to pick a step just as things start to go
   wrong.  In this example, the solver also eventually succeeded.  In
   practice, one often see cases where the newton iterations will
   completely fail to drive the max residual below the tolerances and
   the solver will give up, forcing a retry and reducing the timestep.



Step 3: Check the partial derivatives
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

From the terminal output in the previous part, we can see that the max
residual of ``dv_dt`` equation is difficult to reduce and is located
in cell 1.

We want to find the source of this bad residual by checking
correctness of the partial derivatives of the corresponding equation.
MESA can compare the analytic derivative that it has been given with a
numerical estimate of the derivative based on finite differences and
Ridder's method of polynomial extrapolation (aka ``dfridr``).  To do
so, we must ask for the derivatives at a particular cell at a
particular iteration number of a particular solver call.  This can be
achieved with the inlist options::

      solver_test_partials_call_number = 17
      solver_save_photo_call_number = 17
      solver_test_partials_iter_number = 8
      solver_test_partials_k = 1
      solver_test_partials_equ_name = 'dv_dt'
      solver_test_partials_var_name = 'all'
      solver_test_partials_dx_0 = 1d-6


Setting ``solver_save_photo_call_number = 17`` produces a
photo at the preceding step (``x016``) and so re-running the bad step
becomes as easy as doing ``./re``. Debugging is itself an iterative
process, so you want to make sure you can rapidly and consistently
reproduce the problem.  A little time invested up front letting you
easily and rapidly trigger the problem pays dividends when you have
repeat the procedure a dozen times.


Now when you do ``./rn`` (or ``./re`` after the photo exists) and the
model will stop after the specified solver call and you'll see

.. code-block:: console

                                             equ name dv_dt
    *****                 log dfridr rel_diff partials wrt  lnd(k)   of dv_dt(k)      1     0.677   log uncertainty    -9.668   analytic   -2.9597633549276736E-01   numeric   -1.7033292388187307E+00
                          log dfridr rel_diff partials wrt  lnd(k+1) of dv_dt(k)      1   -99.000   log uncertainty   -99.000   analytic    0.0000000000000000E+00   numeric    0.0000000000000000E+00
    
                          log dfridr rel_diff partials wrt  lnT(k)   of dv_dt(k)      1   -10.058   log uncertainty   -10.644   analytic   -6.9708999322145324E+00   numeric   -6.9708999316051248E+00
                          log dfridr rel_diff partials wrt  lnT(k+1) of dv_dt(k)      1   -99.000   log uncertainty   -99.000   analytic    0.0000000000000000E+00   numeric    0.0000000000000000E+00
    
    ?????                   log dfridr rel_diff partials wrt  L(k)   of dv_dt(k)      1    -4.255   log uncertainty    -4.272   analytic    6.3131382934316702E-39   numeric    6.3127875526588140E-39
                            log dfridr rel_diff partials wrt  L(k+1) of dv_dt(k)      1   -99.000   log uncertainty   -99.000   analytic    0.0000000000000000E+00   numeric    0.0000000000000000E+00
    
                          log dfridr rel_diff partials wrt  lnR(k)   of dv_dt(k)      1   -11.110   log uncertainty   -10.064   analytic   -2.0080506394821409E+00   numeric   -2.0080506394665569E+00
                          log dfridr rel_diff partials wrt  lnR(k+1) of dv_dt(k)      1   -99.000   log uncertainty   -99.000   analytic    0.0000000000000000E+00   numeric    0.0000000000000000E+00


This output shows the derivatives in the cell for each equation with
respect to each variable.  The column after the cell number shows the
relative difference between the numerical (``dfridr``) and analytic
derivatives.  When things are working well, the relative differences
go down to < 1d-12, but don't panic if you see 1d-7.  That can still
work.  But 1d-4 or more is a bad sign.  Here, our bad derivative jumps
out with an order unity relative error and is highlighted by the
``****`` on the left side.  So we now know the bad derivative is the
one with respect to density.

The relative uncertainty in the numerical estimate of the derivative
is also indicated.  When this is comparable to the relative error
between the numerical and analytic derivatives, then you're probably
not learning anything about the quality of the analytic derivative.
(See for example the line marked with ``????`` above.)  Sometimes that
sort of result can be a false alarm caused by the particular choice of
starting difference for the dfridr search.  This is set by
``solver_test_partials_dx_0`` and is typically something around 1d-6.
Usually that works, but if you get a large err that makes no sense,
you should try both larger and smaller values of to see if the dfridr
result changes.

It is often the case that bad derivatives in an equation come from bad
derivatives in the microphysical inputs that enter into the equation.
Since this is common, MESA supports easily checking the basic ``eos``,
``kap``, and ``net`` quantities.  In particular, the value of
``solver_test_partials_equ_name`` can be set to ``'lnE'`` or ``'lnP'``
(testing ``eos``), to ``'eps_nuc'`` (testing ``net``), or
``'opacity'`` (testing ``kap``).

Apply your knowledge about which equation showed the issue.  If the
problem were in the energy equation, you would probably want to test
``lnE`` and ``eps_nuc``.  Since this issue appeared in the momentum
equation, we will want to test ``lnP`` and since it was in the surface
cell, we will want to test ``opacity`` as that enters into the
atmosphere boundary condition.

We already found the issue was in the density derivatives, so we can
now restrict our attention to that derivative by setting
``solver_test_partials_var_name = 'lnd'``.

Setting ``solver_test_partials_equ_name == 'lnP'`` and restarting
shows no problem

.. code-block:: console

                            log dfridr rel_diff partials wrt  lnd(k)   of lnP(k)      1    -9.791   log uncertainty   -12.139   analytic    9.9965278727210283E-01   numeric    9.9965278743368469E-01
                            log dfridr rel_diff partials wrt  lnd(k+1) of lnP(k)      1   -99.000   log uncertainty   -99.000   analytic    0.0000000000000000E+00   numeric    0.0000000000000000E+00


but ``solver_test_partials_equ_name == 'opacity'`` and restarting shows

.. code-block:: console

    *****               log dfridr rel_diff partials wrt  lnd(k)   of opacity(k)      1     0.301   log uncertainty    -9.108   analytic   -2.9110358120328123E-01   numeric    2.9110358102038758E-01
                        log dfridr rel_diff partials wrt  lnd(k+1) of opacity(k)      1   -99.000   log uncertainty   -99.000   analytic    0.0000000000000000E+00   numeric    0.0000000000000000E+00

Comparing the bad analytic derivative to the numerical one reveals the
sign error that we introduced.
                        
We've reached the end of the part of debugging that can be done with
inlists only.  We will need to edit the code in order to MESA to dive
deeper into a module or to check individual terms in the equation
derivatives.

As a test case maintainer, this may be the point at which you stop,
report the issue, and turn things over to others.  If the issue is in
the ``eos``, you might dig in a bit more to see which component of the
eos is active by setting ``solver_test_partials_write_eos_call_info =
.true.`` and include that info in your report.


Step 5: Descend into an individual module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the bug has already been found, this step isn't necessary, but
often one must continue deeper into a module.  For the modules
``atm``, ``eos``, ``kap``, ``net``, and ``rates``, it is possible to
continue to do some of this debugging from star.  This is extremely
helpful because it means that the precise properties of the
troublesome cell from star are passed into the module.

In the inlist set::

      solver_test_kap_partials = .true.

(Note that using this option automatically sets MESA to run in
single-thread mode.)

Since will specify which quantity to check in the code, we should stop
selecting the equation in the inlist, but can continue to specify
which derivative of that quantity to check::

      solver_test_partials_equ_name = ''
      solver_test_partials_var_name = 'lnd'
         
The variables ``kap_test_partials_val`` and
``kap_test_partials_dval_dx`` can now be set within ``kap`` and
thereby communicate their values into star and dfridr.

At the end of ``Get_kap_results`` in ``kap/private/kap_eval.f90`` ---
near where we introduced our bug --- add

.. code-block:: fortran

         kap_test_partials_val = log(kap_rad)
         kap_test_partials_dval_dx = dlnkap_rad_dlnRho

Again, the key point of this step is that we're passing up a value
(``kap_rad``) that we wouldn't normally have access to from star.
This procedure allows you to debug deep within a module from star.

We changed code in ``kap``, so we need to do ``./mk && ./export`` in
``star`` and then in the work directory also recompile before we
restart ``./mk && ./re``.

The output should clearly illustrate that this derivative has the
wrong sign.  Before making a change that corrects the value of
``dlnkap_rad_dlnRho``, consider putting what you think the right
answer is in ``kap_test_partials_dval_dx`` and re-running to check
that your new expression is OK.  (If you first fix the bad derivative
at the location of the bug, that change will mean that the newton
iterations won't be the same any more and you may not hit the same
conditions and or trigger ``dfridr``.)

After you've checked, make the bug fix, turn off the debugging
options, run your model, and hopefully see improved performance.


Advanced: Investigate the bad cell and equation in detail
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note:: 

    We will now pretend we didn't find the bug and re-use this setup
    to illustrate how to dig in at an even deeper level.

MESA equations and structure variables have integer indices.  The full
list is in ``star_data/public/star_data.inc``, but the most common
ones are

.. code-block:: fortran
  
      ! index definitions for the variables (= 0 if variable not in use)
         integer :: i_lnd ! ln(cell density average by mass)
         integer :: i_lnT ! ln cell temperature average by mass
         integer :: i_lnR ! ln radius at outer face of cell
         integer :: i_lum ! luminosity at outer face of cell
         integer :: i_v ! Lagrangian velocity at outer face of cell
         integer :: i_u ! Lagrangian velocity at center of cell
         integer :: i_lnPgas ! ln cell gas pressure average by mass

      ! index definitions for the equations (= 0 if equation not in use)
         integer :: i_dlnR_dt ! dlnR/dt = v/r
         integer :: i_dv_dt ! momentum conservation using v
         integer :: i_du_dt ! momentum conservation using u
         integer :: i_dlnd_dt ! mass conservation
         integer :: i_dlnE_dt ! energy conservation
         integer :: i_equL ! luminosity


In Step 3, we found that we want to investigate the ``lnd`` derivative
of the ``dv_dt`` equation in cell ``k = 1``.  If we hadn't been able
to so easily narrow things down to ``kap``, then we would need to make
a more detailed investigation of the equation in this cell.

To do so, go to ``star/private/hydro_eqns.f90`` and follow the flow
subroutine ``eval_equ_for_solver`` until you reach the equation that
you want investigate.

The surface cell (``k = 1``) is special and so we see that the
``dv_dt`` equation for this cell is provided by ``PT_eqns_surf``.  Go
to that subroutine.  Continue to follow the logic until you get the
place that defines the residual of the desired equation.  In this
case, that is yet one more level down in ``set_Psurf_BC``.

There you'll see the code that indicates that we've reached the
``dv_dt`` equation in cell 1.

.. code-block:: fortran

      if (s% u_flag) then
         i_eqn = i_du_dt
      else
         i_eqn = i_dv_dt
      end if
      
      equ(i_eqn,1) = lnP_bc - s% lnP(1)

This is an extremely simple equation.  By driving the residual to
zero, MESA requires that the surface cell pressure match the imposed
boundary condition.
    
Now, we want to tell MESA to test the partial derivatives.  Near the
start of the routine, set

.. code-block:: fortran
  
      test_partials = (s% solver_iter == s% solver_test_partials_iter_number)
      !test_partials = .false.

Almost all the equation routines will already have ``test_partials``
code like this.  And once you reach one routine that has
``test_partials``, most of the routines it calls are usually also set
up to check partials.  (But if you are in uncharted territory you
might have to add them yourself.)

.. note::
    
    Since ``PT_eqns_surf`` and ``set_Psurf_BC`` are only called for
    ``k=1``, we do not test for zone number.  In other equations that
    are evaluated for multiple ``k``, the analogous test will usually
    look like
    
    .. code-block:: fortran
      
          test_partials = (k == s% solve_test_partials_k .and. s% solver_iter == s% solver_test_partials_iter_number)
          !test_partials = .false.
    
Set ``solver_test_partials_val`` to a variable that is used in
calculating the equation residual.  Let's begin with the residual
itself, which should reproduce the bad entry from the previous step in
this tutorial.

.. code-block:: fortran

      if (test_partials) s% solver_test_partials_val = equ(i_eqn,1)
         
Find the call on ``e00`` for ``i_lnd`` to find the relevant partial with respect to density.

.. code-block:: fortran

      call e00(s, xscale, i_eqn, i_lnd, 1, nvar, dlnP_bc_dlnd - s% chiRho_for_partials(1))

.. note::

   Recall the MESA Jacobian has a block tridiagonal structure (Fig. 47
   in |MESA II|).  The ``e00`` idiom refers to forming the block
   corresponding to the derivatives of the equations in cell ``k``
   with respect to the variables in cell ``k``.  Similarly, ``em1``
   and ``ep1`` refer to the the derivatives of the equations in cell
   ``k`` with respect to the variables in cell ``k-1`` and in cell
   ``k+1`` respectively.  The 3rd and 4th arguments to these
   subroutines indicate that this call is setting the derivative of a
   specific equation (recall ``i_eqn = i_dv_dt`` earlier) with respect
   to a specific variable (here ``i_lnd``).
      
Set ``solver_test_partials_var`` to the thing you want to take a
derivative with respect to and ``solver_test_partials_dval_dx`` to the
analytic derivative.
         
.. code-block:: fortran

      if (test_partials) then
         s% solver_test_partials_var = i_lnd
         s% solver_test_partials_dval_dx = dlnP_bc_dlnd - s% chiRho_for_partials(1)
         write(*,*) 'set_Psurf_BC', s% solver_test_partials_var
      end if

      
We changed code in star, so we need to do ``./mk && ./export`` in
``star`` and then in the work directory also recompile before we
restart ``./mk && ./re``.  We will repeat this recompilation/restart
process many times, so do it in a way that makes it easy for you do
over and over.

Since we are now specifying which derivative to check in the code, we
must also stop selecting the equations/variables in the inlist.

.. code-block:: fortran

      solver_test_partials_equ_name = ''
      solver_test_partials_var_name = ''
      
      solver_test_kap_partials = .false.

After that edit, once you recompile and restart you should see

.. code-block:: console

        set_Psurf_BC           1
    
    
                      analytic and numeric partials wrt lnd   -2.9597633549276736D-01   -1.7033292388187307D+00
        log dfridr relative uncertainty for numeric partial   -9.6678068357602811D+00
                                                   rel_diff    4.7549507665296238D+00
    
    STOP done solver_test_partials


Again, this looks bad!  In the final line, we can see that the
analytic estimate of the derivative (what goes in the Jacobian) is
order-unity different from the numerically evaluated derivative.

Now that we've confirmed which derivative is the problem, we must look
for which term or terms are responsible.  There are only two in this
equation, ``lnP_bc - s% lnP(1)``, so you can test them one by one

.. code-block:: fortran

      s% solver_test_partials_val = lnP_bc

and

.. code-block:: fortran

      s% solver_test_partials_var = i_lnd
      s% solver_test_partials_dval_dx = dlnP_bc_dlnd

Recompile and restart to confirm that this is the bad term.

You can repeat this with the other term to confirm that ``lnP(1)`` and
its derivative ``s% chiRho_for_partials(1)`` are OK.  If they weren't,
then it would be off to the EOS to figure out why.  We'll discuss how
to dive deeply into a specific module in the next step.

Once you locate the bad term, find where that term is defined, so that
you can investigate what is causing it to be in error.  On our way
down, we saw that ``dlnP_bc`` and its derivatives are set in
``PT_eqns_surf``.  So now turn off the partial tests in
``set_Psurf_BC``, switch to ``PT_eqns_surf``, turn on the partial
tests there.

.. note::

    When you finish testing a partial in one place, don't forget to
    turn off the test in that place before moving on to the new place.
    Typically the code that reports results by setting
    ``solver_test_partials_val``, ``solver_test_partials_var``, and
    ``solver_test_partials_dval_dx``, will also have a write that says
    who it is.  This can be a lifesaver if/when you do forget to turn
    off one routine before turning on another.  The terminal output
    will show the problem and tell you where it is.

We can start by testing the partial of ``lnP_bc`` to confirm our
previous finding.  (Sample code omitted so you can practice setting
the variables yourself.)  That should give

.. code-block:: console

       set_Psurf_BC           1
    
                      analytic and numeric partials wrt lnd    7.0157923462638694D-01   -7.0157923451592796D-01
            dfridr relative uncertainty for numeric partial   -1.6481349835630704D-12
                                                   rel_diff    2.0000000001574434D+00
    
    STOP done solver_test_partials

and at this point you can clearly see the sign error we introduced.
With that confirmed, we want to start dissecting

.. code-block:: fortran

      dlnP_bc_dlnd = dlnP_bc_dlnPsurf*dlnPsurf_dlnkap*dlnkap_dlnd

For more complex expressions involving many terms, you might need to
test term by term.

.. note:: 

    The ``dfridr`` machinery in MESA tests the derivatives with
    respect to structure variables, so it can't magically check every
    derivative.  To check ``dlnP_bc_dlnPsurf`` you can look at the
    expression relating ``lnP_bc`` and ``lnP_surf`` in this routine
    and confirm it with pencil and paper.  To check
    ``dlnPsurf_dlnkap`` you would have to follow the code down into
    ``atm``.  Once there, you could check by hand or set up your own
    ``dfridr`` code that passes in varying values of ``kap`` and
    exercises the module independent of star.

Here, since there's only one ``dlnd`` derivative, ``dlnkap_dlnd`` is
the clear place to start (even if we didn't already know it was the
issue).  Set
    
.. code-block:: fortran

    s% solver_test_partials_val = log(s% opacity(1))
            
    s% solver_test_partials_var = s% i_lnd
    s% solver_test_partials_dval_dx = dlnkap_dlnd

and ``dfridr`` should confirm that the bad derivative is the density
derivative of the opacity.

The surface cell momentum equation is one of the simpler equations.
Two other types of derivatives often arise: derivatives with respect
to quantities in neighboring cells and derivatives with respect to
composition.  The energy conservation equation
(``get1_energy_equation`` in ``hydro_equ_l.f90``) provides a good
example of an equation with both these features.

When testing a derivative like ``d_esum_dlnTm1``, that is a partial of
a quantity ``esum(k)`` with respect to ``lnT(k-1)``, note that the
``dfridr`` code perturbs the structure of cell ``k =
solver_test_partials_k``.  That is, it can test the partials of
equations in cell ``k-1``, ``k``, and ``k+1`` with respect to the
structure variables in cell ``k``. (In practice if there is a problem
with partials of a quantity in cell ``k`` then there are also problems
in cells ``k+1`` and ``k-1``.  But if you really need to test one of
these derivatives in a specific cell, you can always adjust
``solver_test_partials_k`` in the inlist.)

Therefore, to check ``d_esum_dlnTm1`` one can do

.. code-block:: fortran
  
      test_partials = (k-1 == s% solver_test_partials_k .and. s% solver_iter == s% solver_test_partials_iter_number)
      !test_partials = .false.

and
      
.. code-block:: fortran
                
      s% solver_test_partials_val = esum
      s% solver_test_partials_var = s% i_lnT
      s% solver_test_partials_dval_dx = d_esum_dlnTm1

This will check the derivative of ``esum`` in cell
``solver_test_partials_k+1`` with respect to ``lnT`` in cell
``solver_test_partials_k``.  


Checking composition derivatives requires setting a sink isotope such
that the dfridr perturbations to the chosen isotope are offset by
changes in a "sink" isotope to keep the mass fraction sum = 1.  The
sink needs to be inert in the current situation and have enough
abundance to absorb the changes.

When debugging from via the inlist, this isotope can be specified
using the control ``solver_test_partials_sink_name``.  So for example,
addition to setting the usual k, call_number, iter_number, and dx_0,
to check partials of the energy equation with respect he4, you might
set the following

.. code-block:: console

      solver_test_partials_equ_name = 'dlnE_dt'
      solver_test_partials_var_name = 'he4'
      solver_test_partials_sink_name = 'fe56'

If ``solver_test_partials_dx_0`` is too large, you'll get an error if
the mg24 mass fraction is too small to deal with the desired increase
in the targeted isotope.

When you need to check composition derivatives of specific terms, you
must pick the various ``test_partials`` values in the source code.
This can be tricky in that in some places the abundance is given as an
index in xa and in others as a variable number for the solver.  This
offset is ``s% nvar_hydro``, that is the index of the first chem
equation is ``s% equchem1 == nvar_hydro + 1``.

If you set ``solver_test_partials_write_eos_call_info = .true.``, the
output includes also includes composition information like the
following so it is easy to find the needed numbers.

.. code-block:: console

     xa(j,k) neut           1           6        3312    1.0010960558047912D-99
       xa(j,k) h1           2           7        3312    0.0000000000000000D+00
     xa(j,k) prot           3           8        3312    1.0010960558047908D-99
      xa(j,k) he3           4           9        3312    0.0000000000000000D+00
      xa(j,k) he4           5          10        3312    6.7423592032960378D-01
      xa(j,k) c12           6          11        3312    1.0630746138489515D-01
      xa(j,k) n14           7          12        3312    6.1988280539418319D-03
      xa(j,k) o16           8          13        3312    2.0390257320339500D-01
     xa(j,k) ne20           9          14        3312    5.7660956389211995D-03
     ...
     xa(j,k) fe56          20          25        3312    1.4153093511237100D-03


So then the in the code you would do something like

.. code-block:: fortran

        if (test_partials) then   
           s% solver_test_partials_dx_sink = 20 ! fe56 (index in xa)
           s% solver_test_partials_var = 10 ! he4 (variable number) 
           s% solver_test_partials_dval_dx = s% d_epsnuc_dx(5,k) ! he4 (index in xa)
           write(*,*) 'get1_energy_eqn', s% solver_test_partials_var
        end if




Diagnosing Meshing Problems
---------------------------

.. note::

   This was the file DEBUG_mesh.notes.  It needs modernized.

to get info about why the mesh is doing what it is doing,

in controls section of inlist::
  
   show_mesh_changes = .true.
   
run to get current mesh call number, n

in controls section of inlist::
  
   mesh_dump_call_number = n
   
run to get plot data

mesh_plan.rb shows input data that was used to form the plan

mesh.rb shows new mesh (with some old values interpolated to the new mesh for comparison)
   

.. _fpe:

Floating point exceptions
-------------------------

Floating point exceptions (fpe) occur when the code attempts to perform an illegal math operation, for instance x/0, sqrt(-1), or log(0). In these cases any number that uses the result of one of these expressions is also undefined.

To test MESA for floating point exceptions we can turn on the compiler checks that will flag when an issue occurs:

If using MESA > 14503 then set the environment variable:

.. code-block:: sh

    MESA_FPE_CHECKS_ON=1

Then run ./clean and ./mk in MESA_DIR

If using MESA <= 14503

In utils/makefile_header we switch:

.. code-block:: makefile

    SKIP_TRAPS = YES

to

.. code-block:: makefile

    SKIP_TRAPS = NO

and insert the variable $(FCtrapNANs) into the FCbasic line

.. code-block:: makefile

    FCbasic = -Wno-uninitialized -fno-range-check -fmax-errors=7 $(SPECIAL_FC_FLAGS) $(FCbasic2) $(FCtrapNANs)

Then run clean and mk to rebuild MESA and fix any issues found in the compilation step. If they are all fixed then run a model or a test case.

You should also add to your controls inlist:

.. code-block:: fortran

      fill_arrays_with_nans = .true.

To make sure you catch any uninitialized array accesses.


Step -1: Introduce a bug (again)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lets add a fpe in a test case to see what happens:

Alter star/test_suite/1.3M_ms_high_Z/src/run_star_extras.f90 and in extras_finish_step add:

.. code-block:: fortran

    write(*,*) log(-1d0*s% model_number)

Step 0: Notice your code has a problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run ./mk and ./rn to run the test case and you'll get the following output:

.. code-block:: console

      Program received signal SIGFPE: Floating-point exception - erroneous arithmetic operation.

      Backtrace for this error:
      #0  0x7efca94096af in ???
      #1  0x7efca9cc6330 in ???
      #2  0x405f10 in __run_star_extras_MOD_extras_finish_step
            at ../src/run_star_extras.f90:155
      #3  0x425372 in __run_star_support_MOD_after_step_loop
            at ../job/run_star_support.f90:634
      #4  0x427d63 in __run_star_support_MOD_run1_star
            at ../job/run_star_support.f90:162
      #5  0x407076 in __run_star_MOD_do_run_star
            at ../../../../star/job/run_star.f90:26
      #6  0x4070eb in run
            at ../src/run.f90:13
      #7  0x407137 in main
            at ../src/run.f90:2
      ./rn1: line 6: 2166962 Floating point exception(core dumped) ./star

Step 1 Diagnose the issue
~~~~~~~~~~~~~~~~~~~~~~~~~

The important line is the first line that does not include ???, in this case it is

.. code-block:: console

      #2  0x405f10 in __run_star_extras_MOD_extras_finish_step
            at ../src/run_star_extras.f90:155

This is telling us that the problem occurs at line 155 in ../src/run_star_extras.f90

.. note:: 
      The ../src is because the file path is relative to the make/ folder in your work folder

At this point you should go to that line and work out why the number is undefined.


Uninitialised variables
-----------------------

Uninitialised variables occur when we use a variable before it has been defined, for example

.. code-block:: fortran

      real(dp) :: x,y
      y = 2 * x

If you are lucky then x will be zero, if you are unlucky it can be any random number that happened to be in the memory location that x occupies. This can lead to results that change either each time the model runs or gives different values between different machines.

By default MESA compiles with the -finit-real=snan option, this will set all scalar floating point variables to a NaN when they are declared.

.. note::
      This does not protect integers

As the variable will be a NaN we can use the floating point exceptions mechanism to track them down, so follow the instructions for debugging floating point exceptions.

Arrays of floating point numbers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These can also be protected by wrapping the array after its declared with:


      call fill_with_NaNs(array)

.. code-block:: fortran

      use utils_lib, only:: fill_with_NaNs
      real(dp),dimension(:),allocatable :: x

      allocate(x(1:10))
      call fill_with_NaNs(x)

.. note:: 
      Never try to fill the array yourself with a NaN, say by setting each element to 1/0. We want to track down the location where an array element is being **used** before being set, not the point where we set the element to a NaN. If you set the array yourself the fpe checking will warn you about that location, not where the array is being before being set.  fill_with_NaNs does some clever tricks to set each element to the bit pattern that corresponds to a NaN, which does not trigger the fpe checks.

Memory leaks
------------
