Using PGSTAR
============

Watching lots of numbers fly by in the terminal window can provide some
idea of how a run is going, but it is much better (and much more fun) to
have windows with plots that change at each step. PGstar is built on
PGPLOT for creating plots using XWindows for on-screen plots and PNG for
frames of movies to animate those plots.

Devote lots of time watching models evolve.  It's a fun way to learn!

If you're using the MESA SDK, you should have everything you need to use
PGSTAR. Otherwise, you need to have installed PGPLOT, XWindows, and a
few libraries that they need. See the makefile_header_non_mesasdk in
$MESA_DIR/utils for details. You'll also find pgplot.tar.bz2 in
$MESA_DIR/utils in case you need to install PGPLOT (thanks to Tim
Pearson of CalTech for creating PGPLOT and giving us permission to
distribute it with MESA).

Overview of PGSTAR
------------------

PGSTAR has an enormous number of features and controls. Here is a brief
overview.

Profile Plots and History Plots
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Profile plots show information about current model. You can display
anything that can be included in a profile_columns.list. You are not
limited to the items in your current profile_columns.list.

History plots show information specified by your current
history_columns.list. If it is in your history_columns.list, you can
show it in a history plot. The name for a history item must be the same
as one of the column headings in your current LOGS/history.data - these
are the same as the entries in the history_columns.list but with any
spaces replaced by "\_". For example, "center h1" in
history_columns.list changes to "center_h1" as a column heading in
LOGS/history.data, so you should use "center_h1" as the name in the
plot.

Single Panel, Multiple Panel, and Grid Plots
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Single panel plots contain one graph, optionally with a 2nd Y axis.
Multiple panel plots have a stack of several graphs using the same X
axis. Grid plots combine several other plots in a user-specified grid
layout.

General Controls
~~~~~~~~~~~~~~~~

You can change the width, aspect ratio, margins, text scale, title, and
axis limits of any plot by editing your pgstar_controls inlist file. The
file is reread at each timestep, so you can make changes while the
evolution is running. For user-specified plots, you can also change any
of the specs on-the-fly including what is plotted for the axes and the
number of panels or grid layout. You can change the "win_flag" from
.false. to .true. on-the-fly to create a new window. Changing the
win_flag from .true. to .false. removes the window, but because of
quirks of pgplot and/or Xwindows, it will pause and ask you to hit a
return.

You can have the code pause at each step by setting pause = .true. It
will pause every nth step if you set pause_interval = n.

You can slow it down by setting pgstar_sleep to the minimum number of
seconds to allow between plot updates.

It will update the windows every nth step if you set pgstar_interval = n.

To save the plot as a png file, set the "file_flag" for the window to
.true. You can specify a directory and prefix for the files, as well as
the frequency of writing (every nth step) and the size of the plot for
the file. You can even get files written when the star moves a given
distance on the HR diagram!

To switch from black background to white change the white_on_black_flag.
There are separate flags for windows and files.

The Inventory of Plots
----------------------

single panel profile plots
~~~~~~~~~~~~~~~~~~~~~~~~~~

-  TRho_Profile - current model in T-Rho plane
-  Summary_Profile - various profile properties
-  Summary_Burn - current model T, rho, eps burning, and eps neutrinos
-  Abundance - current model abundance profiles
-  Power - current model nuclear power profiles
-  Mixing - current model mixing diffusion coefficients
-  Dynamo - current model dynamo magnetic field info
-  Mode_Prop - asteroseimology mode propagation diagram

single panel history plots
~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Kipp - "Kippenhann" history of mixing, burning, and more
-  Summary_History - various history properties
-  TRho - history of central temperature vs. density
-  HR - history of log_L vs. log_Teff
-  logg_Teff - history of logg vs. Teff
-  logg_logT - history of logg vs. log_Teff
-  dPg_dnu - delta_Pg vs. delta_nu (for asteroseismology)

Keep in mind that if it isn't in your current history list, then it
won't be available for plotting in PGSTAR history plots. For example, if
your Kipp plots doesn't show burning and mixing regions, then add them
to your history_columns.list - e.g., add these:

::

   mixing_regions 40
   burning_regions 80

profile panel plots: Profile_Panels1 to Profile_Panels9
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You select axes from anything allowed in a profile_columns.list. In
addition, a panel can be any of these pre-defined plots: Abundance,
Power, Mixing, Dynamo, Mode_prop, or Summary_Profile.

-  Profile_Panels1 - default has X = mass, panels = (logT,entropy) and
   (logRho,logP).
-  Profile_Panels2 - default has X = logP, panels = Abundance and Power.
-  Profile_Panels3 - default has X = logP, panels = Abundance, Power,
   and Mixing.
-  Profile_Panels4 - default has X = logP, panels = Abundance, Power,
   Mixing, and Dynamo.
-  Profile_Panels5 - default has X = logP, panels = Summary_Profile,
   Mode_prop, and Mixing.

history track plots: History_Track1 to History_Track9
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Like TRho, HR, and other history plots, but you select the axes from
items in your current history_columns.list.

-  History_Track1 - default shows log_L vs. log_center_T

history panel plots: History_Panels1 to History_Panels9
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You select Y axes from the column heading in current LOGS/history.data.

-  History_Panels1 - default shows logT, logRho, log_L, log_Teff, mass,
   mdot

text summary plots: Text_Summary1 to Text_Summary9
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Text_Summary1 - default gives info similar to standard terminal
   output

grids: Grid1 to Grid9
~~~~~~~~~~~~~~~~~~~~~

-  Grid1 - default TRho_Profile, HR, TRho, Text_Summary1
-  Grid2 - default TRho_Profile, Kipp
-  Grid3 - default TRho_Profile, Summary_Profile, Kipp
-  Grid4 - default TRho_Profile, Summary_Profile, HR, TRho, Kipp,
   Text_Summary1
-  Grid5 - default TRho_Profile, Summary_Profile, HR, TRho, Kipp
-  Grid6 - default Summary_Burn, Abundance, HR, TRho, Text_Summary1
-  Grid7 - default Abundance, TRho, Text_Summary1
-  Grid8 - default Summary_Burn, Abundance, HR, TRho, TRho_Profile,
   Text_Summary1
-  Grid9 - default Abundance, HR, TRho, TRho_Profile, Text_Summary1

Hands-on Tutorial
-----------------

Do this tutorial while $MESA_DIR/star is running and pgstar_flag is set
.true. so you can watch as the plots change as you edit and save your
inlist.

I suggest using the test case 1M_pre_ms_to_wd. Go to the directory,
$MESA_DIR/star/test_suite/1M_pre_ms_to_wd, and open "inlist_1.0" in your
text editor. In the &star_job section, set pgstar_flag = .true. and
check the &pgstar section to make sure it starts empty. Then do the
usual ./mk and ./rn to start the test case.

First, open some plot windows and take a look at what's available. Just
edit your pgstar inlist to set the "win_flag" for the plot; the changes
will take effect on the next step after you save the file; you don't
need to stop and restart the run to change the plots.

Start by opening the pre-defined profile plots. Copy-and-paste the
following lines to your inlist, then do save. Take a quick look at each
to familiarize yourself with the options.::

  TRho_Profile_win_flag = .true.
  Summary_Profile_win_flag = .true.
  Abundance_win_flag = .true.
  Power_win_flag = .true.
  Mixing_win_flag = .true.
  Dynamo_win_flag = .true.
  Mode_Prop_win_flag = .true.


Next replace those lines by ones to open some of the history plots.
Copy-and-paste these lines to replace the previous ones. Then do save.
PGPLOT will ask you to hit return to close the windows that were
previously open.::

  Kipp_win_flag = .true.
  TRho_win_flag = .true.
  HR_win_flag = .true.
  Summary_History_win_flag = .true.


After looking at those, check out some of the defaults for a few of the
multi-panel plots. As before, cut-and-paste these lines and replace the
previous ones.::

  Profile_Panels1_win_flag = .true.
  Profile_Panels2_win_flag = .true.
  History_Track1_win_flag = .true.
  History_Panels1_win_flag = .true.
  Grid1_win_flag = .true.
  Grid2_win_flag = .true.

  
When you're finished admiring these beauties, delete the "win_flag"
lines from your inlist, and we'll take a look at how you can modify
plots.

Let's start by changing some of the Abundance plot controls; do the
changes one at a time so you can see each change. And make sure the plot
is visible when you do the save!::

  Abundance_win_flag = .true.
  Abundance_win_width = 6
  Abundance_win_aspect_ratio = 1.2
  Abundance_log_mass_frac_min = -4
  Abundance_xaxis_name = 'logP'
  Abundance_xmin = 12


You can change the set of lines on the Summary_Profile plot. There is a
control for the number of lines, and name and legend strings for each
line. Each value can be shown scaled or unscaled. The scaling is max to
1.0 and min to 0.0. Usually mass fractions are the only things shown
unscaled. Here are some examples from the default settings. Note the use
of the PGPLOT symbol for "omega".::

  Summary_Profile_win_flag = .true.
  Summary_Profile_num_lines = 11 ! <= 16

  Summary_Profile_name(1) = 'x' ! hydrogen mass fraction
  Summary_Profile_legend(1) = 'X'
  Summary_Profile_scaled_value(1) =.false.

  Summary_Profile_name(2) = 'y' ! helium mass fraction
  Summary_Profile_legend(2) = 'Y'
  Summary_Profile_scaled_value(2) = .false.

  Summary_Profile_name(3) = 'log_j_rot' ! specific angular momentum
  Summary_Profile_legend(3) = 'log j rel'

  Summary_Profile_name(4) = 'log_omega' ! angular velocity
  Summary_Profile_legend(4) = 'log (0650) rel'

Open the Summary_Profile window, then make some changes. Note that you
can remove a line just by setting the name to ''; you don't have to
renumber the other lines.::

  Summary_Profile_num_lines = 12 ! <= 16
  Summary_Profile_name(12) = 'zone'
  Summary_Profile_name(3) = ''
  Summary_Profile_name(4) = ''

You can change the set of lines on the Summary_History plot. There is a
control for the number of lines, and name and legend strings for each
line. Each value can be shown scaled or unscaled. The scaling is max to
1.0 and min to 0.0. Usually mass fractions are the only things shown
unscaled. Here are some examples from the default settings. Note the use
of PGPLOT text controls for the subscript "c" for center values.::

  Summary_History_num_lines = 7 ! <= 16

  Summary_History_name(1) = 'log_center_T'
  Summary_History_legend(1) = 'log T\dc\u rel'

  Summary_History_name(2) = 'log_center_Rho'
  Summary_History_legend(2) = 'log Rho\dc\u rel'

Open the Summary_History window, then make some changes.  Note that
you can remove a line just by setting the name to ''; you don't have
to renumber the other lines.::

  Summary_History_win_flag = .true.
  Summary_History_num_lines = 9
  Summary_History_name(8) = 'star_mass'
  Summary_History_name(9) = 'log_abs_mdot'
  Summary_History_name(6) = ''
  Summary_History_name(7) = ''

Next, turn on the History\_Track1 plot.  Then change what it is
plotting by editing the axis name and label. Here's an example.  After
that, try plotting some other combinations; just pick axis names from
the column headings in your current LOGS/history.data.::

  History_Track1_win_flag = .true.

  History_Track1_xname = 'log_center_P'
  History_Track1_xaxis_label = 'log P\dcenter'
  History_Track1_title = 'L vs Center P'

Turn on Profile\_Panels1 and History\_Panels1; they are set up with
defaults for the number of panels and axes.  Change the defaults to
show other things -- for the profiles you can select anything that can
be in a profile_columns.list; for the history you have to select one
of the column headings in your current LOGS/history.data.::

  History_Panels1_win_flag = .true.
  History_Panels1_other_yaxis_name(1) = 'log_center_P'

  Profile_Panels1_win_flag = .true.
  Profile_Panels1_xaxis_name = 'logP'
  Profile_Panels1_xaxis_reversed = .true.

Add another panel to the Profile plot.::

  Profile_Panels1_num_panels = 3
  Profile_Panels1_yaxis_name(3) = 'logtau'
  Profile_Panels1_other_yaxis_name(3) = 'log_opacity'

Increase the y margins.::

  Profile_Panels1_ymargin(1) = 0.2
  Profile_Panels1_other_ymargin(1) = 0.2
  Profile_Panels1_ymargin(2) = 0.2
  Profile_Panels1_other_ymargin(2) = 0.2
  Profile_Panels1_ymargin(3) = 0.2
  Profile_Panels1_other_ymargin(3) = 0.2

Change the aspect ratio, reduce the width, and fix the left and right margins.::

  Profile_Panels1_win_aspect_ratio = 1.0 ! aspect_ratio = height/width
  Profile_Panels1_win_width = 6
  Profile_Panels1_xleft = 0.18
  Profile_Panels1_xright = 0.82

Pick some other axis names and change what is shown in the panels.
You can use any name that is valid in a profile_columns.list; unlike
the history case, you don't have to limit yourself to the contents of
your current list.

Next, take a look at the following default definition for
Profile\_Panels3.::

  Profile_Panels3_xaxis_name = 'logP'
  Profile_Panels3_xaxis_reversed = .true.
  Profile_Panels3_num_panels = 3
  Profile_Panels3_yaxis_name(1) = 'Abundance'
  Profile_Panels3_yaxis_name(2) = 'Power'
  Profile_Panels3_yaxis_name(3) = 'Mixing'

Open the plot window and then change the number of panels and the
contents.  Revise the title and switch the xaxis to mass.::

  Profile_Panels3_win_flag = .true.
  Profile_Panels3_num_panels = 4
  Profile_Panels3_yaxis_name(4) = 'Dynamo'
  Profile_Panels3_title = 'Abundance-Power-Mixing-Dynamo'
  Profile_Panels3_xaxis_name = 'mass'
  Profile_Panels3_xaxis_reversed = .false.

Now, edit the definition of the Grid1 plot.  Replace the TRho\_Profile
plot by the Kipp plot and adjust the margins and text scale.::

  Grid1_win_flag = .true.
  Grid1_plot_name(1) = 'Kipp'
  Grid1_plot_pad_left(1) = 0.03 ! fraction of full window width for padding on left
  Grid1_plot_pad_right(1) = 0.03 ! fraction of full window width for padding on right
  Grid1_plot_pad_bot(1) = 0.12 ! fraction of full window height for padding at bottom
  Grid1_txt_scale_factor(1) = 0.6 ! multiply txt_scale for subplot by this

Move the text summary up to just below the Kipp plot, and increase the
number of rows to make the HR and TRho plots taller.  This will
temporarily mess us the spacing between the subplots, but we'll fix
that next.::

  Grid1_num_rows = 9 ! divide plotting region into this many equal height rows
  Grid1_plot_row(2) = 7 ! number from 1 at top
  Grid1_plot_rowspan(2) = 3 ! plot spans this number of rows
  Grid1_plot_row(3) = 7 ! number from 1 at top
  Grid1_plot_rowspan(3) = 3 ! plot spans this number of rows
  Grid1_plot_row(4) = 5 ! number from 1 at top
  Grid1_plot_rowspan(4) = 2 ! plot spans this number of rows

After that, fix the padding between the plots and adjust the text sizes.::

  Grid1_plot_pad_top(2) = 0.01 ! fraction of full window height for padding at top
  Grid1_plot_pad_bot(2) = 0.1 ! fraction of full window height for padding at bottom
  Grid1_plot_pad_left(2) = 0.05 ! fraction of full window width for padding on left
  Grid1_plot_pad_right(2) = 0.1 ! fraction of full window width for padding on right
  Grid1_txt_scale_factor(2) = 0.6 ! multiply txt_scale for subplot by this
  
  Grid1_plot_pad_top(3) = 0.01 ! fraction of full window height for padding at top
  Grid1_plot_pad_bot(3) = 0.1 ! fraction of full window height for padding at bottom
  Grid1_plot_pad_left(3) = 0.1 ! fraction of full window width for padding on left
  Grid1_plot_pad_right(3) = 0.05 ! fraction of full window width for padding on right
  Grid1_txt_scale_factor(3) = 0.6 ! multiply txt_scale for subplot by this
  
  Grid1_plot_pad_top(4) = 0.00 ! fraction of full window height for padding at top
  Grid1_plot_pad_bot(4) = 0.05 ! fraction of full window height for padding at bottom

Change the text summary to report 'log\_L\_div\_Ledd' instead of 'log\_LH'.::

 Text_Summary1_name(3,4) = 'log_L_div_Ledd'

Finally, take a quick look at the various multi-panel and grid
defaults (listed above).  You are not limited to those, but they might
give you ideas for your own personalized plots.  Here are some of my
favorites -- you might also find them useful.::

  Grid8_win_flag = .true. ! Summary_Burn, Abundance, HR, TRho, TRho_Profile, Text_Summary1
  Profile_Panels4_win_flag = .true. ! Abundance, Power, Mixing, and Dynamo
  History_Panels1_win_flag = .true. ! logT, logRho, log_L, log_Teff, mass, mdot

If you don't have rotation turned on, change from Profile\_Panels4 to
Profile\_Panels3, or edit your controls for Profile\_Panels4 to drop the
last panel.::

  Profile_Panels3_win_flag = .true. ! Abundance, Power, and Mixing

or::

  Profile_Panels4_num_panels = 3
  Profile_Panels4_title = 'Abundance-Power-Mixing'

PGSTAR has a number of options to control its file output.

The default output format is PNG::

  file_device = 'png'
  file_extension = 'png'

but you can use PostScript output by setting::

  file_device = 'vcps'
  file_extension = 'ps'

You can change the foreground/background color of your plots between
black/white and white/black::

  ! white_on_black flags -- true means white foreground color on black background
  file_white_on_black_flag = .true.

and control the number of digits that appear in the filenames::
  
  file_digits = 5 ! number of digits for model_number in filenames

In addition, each plot has controls for its own file output similar to
the following ones for the TRho_Profile:::

  TRho_Profile_file_dir = 'png'
  TRho_Profile_file_flag = .false.
  TRho_Profile_file_prefix = 'trho_profile_'
  TRho_Profile_file_interval = 5 ! output when mod(model_number,TRho_Profile_file_interval)==0
  TRho_Profile_file_width = -1 ! (inches) negative means use same value as for window
  TRho_Profile_file_aspect_ratio = -1 ! negative means use same value as for window


A directory with the name given by the value of
TRho\_Profile\_file\_dir (in this case, the default value, which is
"png") must exist in the work directory in order for the files to be
stored.  Otherwise MESA will run, but fail to write the png files.
Make sure to create such a directory first!  (If you're using
PostScript output, you probably want to create a directory named "ps"
and set TRho\_Profile\_file\_dir='ps'.)

Finally, there is an HR distance trigger for file output.::

  ! trigger file output by distance traveled on HR diagram
  delta_HR_limit_for_file_output = -1 ! negative means no limit

  ! HR distance since last file output = sum of dHR
  ! where per step dHR = same definition as used for timestep limits
  ! dHR = sqrt((delta_HR_ds_L*dlgL)**2 + (delta_HR_ds_Teff*dlgTeff)**2)
  ! dlgL = log10(L/L_prev)
  ! dlgTeff = log10(Teff/Teff_prev)


There are a variety of tools available for combining png files in
movies.

