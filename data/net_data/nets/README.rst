=================
Reaction Networks
=================

Summary of available nets
-------------------------

This is an incomplete summary of the nets included for use with
mesa/net.  See ``$MESA_DIR/data/net_data/nets`` for the full list.

This lists the isotopes that are included in each net, but not the
reactions.  To check on those, it is best to have the code list them .
In mesa/star you can do this by setting the star_job options
``show_net_reactions_info = .true.``; similarly, you can get a list of
the isotopes by setting ``show_net_species_info = .true.``.


"no frills" hydrogen and helium burning
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. describe:: basic

      h1, he3, he4, c12, n14, o16, ne20, mg24

      Assumes T is low enough so can ignore advanced burning and hot cno issues.


more complete coverage for hydrogen and helium burning
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. describe:: c13

   basic + c13

.. describe:: o18_and_ne22

   basic + o18 and ne22

.. describe:: o18_to_mg26

   o18_and_ne22 + mg26

.. describe:: hot_cno

   basic + c13, n13, n15, o15, o17, o18, f17, f18

   for high temperatures where cno becomes beta limited.

.. describe:: cno_extras

   hot_cno + o14, f19, ne18, ne19, mg22

   for high temperatures where start to breakout of hot cno.

.. describe::cno_extras_to_ni56

   cno_extras + s30, ni56

   extends rp breakout up to ni56

   see Wallace & Woosley, ApjS, 45:389-420, 1981, Appendix C

.. describe::_extras

   basic + h2, li7, be7, b8

.. describe:: pp_and_cno_extras

   pp_extras + cno_extras

.. describe:: pp_cno_extras_o18_ne22

   pp_and_cno_extras + o18_and_ne22



extensions of the basic net for C/O burning and alpha chains
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. describe:: co_burn

   basic + si28

.. describe:: alpha_s32

   co_burn + s32

.. describe:: alpha_ar36

   alpha_s32 + ar36

.. describe:: alpha_ca40

   alpha_ar36 + ca40

.. describe:: alpha_ti44

   alpha_ca40 + ti44

.. describe:: alpha_cr48

   alpha_ti44 + cr48

.. describe:: alpha_fe52

   alpha_cr48 + fe52

.. describe:: alpha_ni56

   alpha_fe52 + fe54, ni56, neut, and prot

.. describe:: approx19

   same isotopes as alpha_ni56; similar to Frank Timmes' APPROX19

.. describe:: approx20

   approx19 + fe56

.. describe:: approx21

   approx20 + cr56

.. describe:: ns_he

   approx21 + (a,p)*(p,g) intermediates for alpha links + "bypass" for c12


Creating a custom net
---------------------

If you don't find what you need here, it is fairly easy to create new nets.
First, study o18_and_ne22.net to see how it works.
The idea is to start from a net that is close to what you need, then add
isotopes and reactions as needed; you may also need to remove some reactions.
You can put your new net file in your work directory; mesa/net looks in the
current directory before it checks data/net_data/nets.  Before using
your new net, get a list of isotopes and reactions and carefully check them.

For reference, the current set of isotopes is in mesa/data/chem_data/isotopes.data.
The format of that file is described in a README file in the same directory.
The current set of reactions is in mesa/net_data/reactions.list.  The reaction names
in net files must be defined in the reactions.list file.  Those reactions are
implemented by the mesa/rates module.


Description of net format
-------------------------

A net is defined by a set of isotopes and a set of reactions.
These sets are specified in a net definition file.

A net definition file contains a sequence of commands.
It can have any number of the following commands in any order.

.. describe:: include 'net_filename'

   This lets you break up a definition into several files (e.g., one
   net can defined as an extension of another).  This first tries to
   open 'net_filename' in current directory.  If that fails, it tries
   in mesa/data/net_data/nets/.

.. describe:: add_isos(isos_list)

   For each iso in isos_list, add it to the current set of isos in the net.
   It really is a set, so multiple adds are okay.

.. describe:: remove_isos(isos_list)

   For each iso in isos_list, remove it from the current set of isos in the net.

.. describe:: add_reactions(reactions_list)

   For each reaction in reactions_list, add it to the current set of
   reactions in the net.  It really is a set, so multiple adds are
   okay.

.. describe:: remove_reactions(reactions_list)

   For each reaction in reactions_list, remove it from the current set
   of reactions in the net.

.. describe:: add_isos_and_reactions(isos_list)

   For each iso in isos_list, add the iso to the current set,
   and add reactions linking it to other isos in current set.

   The reactions (see below for descriptions of the reaction handle
   format) that automatically added by the add_isos_and_reactions
   command are:

    * p,a,n capture reactions and photodisintegrations (e.g., r_x_pg_y)
    * p,a,n exchanges (e.g., r_x_ap_y)
    * standard 1-to-1 weak reactions (e.g., r_x_wk_y)
    * other weak reactions (e.g., r_x_wk_h1_y)

   It also checks if there are special cases to be added. Here's an
   incomplete list.

       + r_h1_h1_wk_h2
       + r_h1_h1_ec_h2
       + r_h2_h2_to_h1_h3, r_h1_h3_to_h2_h2
       + r_he3_ec_h3
       + r_h2_h2_to_neut_he3, r_neut_he3_to_h2_h2
       + r_neut_neut_he4_to_h3_h3, r_h1_h1_he4_to_he3_he3
       + r_h3_he3_to_h2_he4, r_h2_he4_to_h3_he3
       + r_li6_to_neut_h1_he4, r_neut_h1_he4_to_li6
       + r_h1_li7_to_h2_li6, r_h2_li6_to_h1_li7
       + r_he4_he4_he4_to_c12, r_c12_to_he4_he4_he4
       + r_c12_c12_to_he4_ne20, r_he4_ne20_to_c12_c12
       + r_neut_mg23_to_c12_c12, r_c12_c12_to_neut_mg23
       + r_c12_ne20_to_h1_p31, r_h1_p31_to_c12_ne20
       + r_s27_wk_h1_h1_al25
       + r_he4_ca36_to_h1_h1_ca38, r_h1_h1_ca38_to_he4_ca36
       + r_pd89_to_h1_h1_ru87, r_h1_h1_ru87_to_pd89

The arguments to these commands have the following form:

.. describe:: isos_list

   A sequence of iso_spec values optionally separated by commas

.. describe:: iso_spec

   Either iso_name or element_name A_lo A_hi.

   iso_name is something like h1 or he4 defined in chem.

   element_name is something like h or he defined in chem.
   A_lo and A_hi define the range of A values to add for the element
   e.g., if iso_spec is he 3 4, then adds he3 and he4.

.. describe:: reactions_list

   A sequence of reaction_handle values optionally separated by commas.

.. describe:: reaction_handle

   A reaction name defined in reactions.list

   or

   A reaction name derived from the input and output isos as follows,
   where x and y are iso_names; e.g., r_c12_pg_n13.

   * p,a,n capture reactions and photodisintegrations

     + r_x_pg_y
     + r_y_gp_x
     + r_x_ag_y
     + r_y_ga_x
     + r_x_ng_y
     + r_y_gn_x

   * p,a,n exchanges

     + r_x_ap_y
     + r_y_pa_x
     + r_x_np_y
     + r_y_pn_x
     + r_x_na_y
     + r_y_an_x

   * standard 1-to-1 weak reactions

     + r_x_wk_y       (*positron emission or electron capture*)
     + r_x_wk-minus_y (*electron emission or positron capture*)

   * other weak reactions (not in weaklib)

     + r_x_wk_h1_y    (* beta decay with products h1 and y*)
     + r_x_wk_he4_y   (*beta decay with products he4 and y*)
     + r_h1_h1_ec_h2  (*electron capture*)
     + r_h1_h1_wk_h2  (*beta decay*)
     + r_he3_ec_h3
     + r_be7_ec_li7
     + r_h1_h1_wk_h2
     + r_h1_he3_wk_he4

   * other reactions

     + r_<inputs>_to_<outputs>

         where <inputs> and <outputs> are iso_name values separated by '_'
         If the same iso appears 2 or more times, repeat the name that many times.
         e.g., triple alpha is r_he4_he4_he4_to_c12.
         iso values are ordered by increasing Z and N. e.g., r_h3_be7_to_neut_h1_he4_he4.


Changing/selecting reaction rates
---------------------------------

.. warning::

   This section needs work.

change rates for existing reactions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   create a file of (T8, rate) pairs as in data/rates_data/rates
   add the file name to reactions.list
   (you can have a local copy of reactions.list as well as rates directory)


add a new reaction with rate give by table of (T8, rate) pairs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   get a name for the reaction by calling reaclib_create_handle
   then proceed same as for changing rate of an existing reaction
   use the reaction in a net in same way as for existing reactions.

selecting rates for reactions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

reactions in reactions.list have various rates available such as NACRE, REACLIB, or CF88.
several reactions have other options as well: see rates_def.f for listing of choices

other reactions use rates from weaklib or reaclib.

also can provide a table of values for any reaction.
arbitrarily many pairs of temperatures and rates; piecewise monotonic cubic interpolation.
specify files for rate tables in rate_list.txt


Adding new reactions
^^^^^^^^^^^^^^^^^^^^

advanced topics -- experts only

to add a reaction
   1) add to reactions.list
   2) add to rates_def
   3) add to rates_names
   4) add to raw_rates

For a weak rate, add 1/2 life and Qneu to weak_info.list in weaklib_data.
If possible, get 1/2 life = log(2)/exp(c1) with c1 from reaclib.

To get the reaction added automatically add it to init_special_case_reaction_info in net_initialize.
