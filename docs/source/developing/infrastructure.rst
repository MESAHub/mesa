==============
Infrastructure
==============

This summarizes some of the external MESA infrastructure and who
controls it.

GitHub
------

We have the `MESAHub <https://github.com/MESAHub>`__ GitHub
organization.  The members of the MTC/MAC all have owner-level
privileges.  MESAHub is classified as an educational organization and
this currently grants us a 100% discount, meaning this costs $0/yr.
Unsubsidized, it would cost about $1200/yr.

To classify MESAHub as an eduction account requires someone in the MESAHub
organization to sign up for a eduction plan. This requires that the user has
an University ID card with an expiry date (this is important). Steps to get
access to this:

- Add your university email address to Github.
- Follow instructions `Apply to Github Global campus <https://docs.github.com/en/education/explore-the-benefits-of-teaching-and-learning-with-github-education/github-global-campus-for-teachers/apply-to-github-global-campus-as-a-teacher>`__
- It is set up, if you have access to this https://github.com/community/Global-Campus-Teachers repository
- Github claims they will periodically require revalidation of your credentials, so this may need repeating in the future.

If MESAHub gets unsynced from this and is no long on the eduction plan:

- Goto `Global campus <https://education.github.com/globalcampus/teacher#>`__
- Find ``Upgrade your academic organizations``
- Click ``Upgrade to Github Team``
- Upgrade ``MESAHub``

This is currently tied to Rob's account.

mesastar.org
------------

This is a long-standing domain name used by MESA.

Domain Name
^^^^^^^^^^^

Matteo controls the domain name (using domain.com) and handles
renewing it, updating nameservers, etc.

DNS Servers
^^^^^^^^^^^

The mesastar.org domain is configured to use Digital Ocean nameservers
and the DNS records are managed by Josiah.

Mailing lists
^^^^^^^^^^^^^

The mesa-users@lists.mesastar.org and
mesa-developers@lists.mesastar.org lists are hosted by `MailmanLists
<https://www.mailmanlists.net/>`__.  The account is controlled by
Frank who pays the hosting fees (through grants at ASU).
This costs about $100/year.

New account signups are handled by Frank, Evan, and Rob.

Marketplace
^^^^^^^^^^^

The address mesastar.org redirects to the `MESA marketplace
<http://cococubed.com/mesa_market/>`__, which is controlled and
hosted by Frank at ASU.


ReadTheDocs
-----------

The Sphinx documentation is hosted by `ReadTheDocs
<https://readthedocs.org/>`__.  Rich controls the account (and Earl, Evan, and Joey have access to it).  This is
currently free for open source software.


Sourceforge
-----------

The original home of MESA development was sourceforge.  Most developers
involved c. 2015 or earlier have admin privileges.  This is a free
service.

Website
^^^^^^^

We still use the mesa.sourceforge.net website domian name, but only for a top-level
redirect to docs.mesastar.org.

The source code for the old sourceforge site lives at https://github.com/MESAHub/mesa-website.

Slack
-----

Rich controls the Slack and pays for it (through grants at UW
Madison).  This costs around $130 per year (about $10/user/yr).


TestHub
-------

The MESA TestHub runs on Heroku and the account is controlled by Bill
Wolf.  It is funded by Frank and Bill W.  The Heroku account costs
around $600/yr, though its cost has been higher during the GitHub
transition.

The testing `log archive <https://logs.mesastar.org/>`__ lives on a
server controlled by Josiah.  The marginal cost is $1/month, which he
covers.

Jenkins
^^^^^^^

The Flatiron `Scientific Computing Core <https://www.simonsfoundation.org/flatiron/scientific-computing-core>`__
provides continuous integration testing for MESA via jenkins.
This was set up in 2021 by Dylan Simon (Flatiron) and Josiah.  The configuration lives in the ``jenkins`` directory.
Matteo serves as the point of contact at Flatiron.

Zenodo
------

We have multiple resources archived on Zenodo.  This is a free service.

Community
^^^^^^^^^

The `MESA community <https://zenodo.org/communities/mesa/>`__ on Zenodo
is controlled by Pablo.

Records
^^^^^^^

* The `record for MESA releases <https://zenodo.org/record/4311514>`__ is controlled by Pablo.
* The `record for OP Mono data <https://zenodo.org/record/4390522>`__ is controlled by Josiah.
* The records for MESA SDK releases (`macOS <https://zenodo.org/record/4638654>`__, `linux <https://zenodo.org/record/4638535>`__) are controlled by Pablo.
