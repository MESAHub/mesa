==============
Infrastructure
==============

This summarizes some of the external MESA infrastructure and who
controls it.

GitHub
------

We have the `MESAHub <https://github.com/MESAHub>`__ GitHub
organization. The members of the MTC/MAC all have owner-level
privileges. MESAHub is classified as an educational organization and
this currently grants us a 100% discount, meaning this costs $0/yr.
Unsubsidized, it would cost about $1200/yr.

To classify MESAHub as an education account requires someone in the MESAHub
organization to sign up for an education plan. This requires that the user has
a university ID card with an expiry date (this is important). Steps to get
access to this:

- Add your university email address to GitHub.
- Follow instructions `Apply to GitHub Global campus <https://docs.github.com/en/education/explore-the-benefits-of-teaching-and-learning-with-github-education/github-global-campus-for-teachers/apply-to-github-global-campus-as-a-teacher>`__
- It is set up, if you have access to the `Global-Campus-Teachers repository <https://github.com/community/Global-Campus-Teachers>`__
- GitHub claims they will periodically require revalidation of your credentials, so this may need repeating in the future.

If MESAHub gets unsynced from this and is no longer on the education plan:

- Go to `Global Campus <https://education.github.com/globalcampus/teacher#>`__
- Find ``Upgrade your academic organizations``
- Click ``Upgrade to GitHub Team``
- Upgrade ``MESAHub``

This is currently tied to Rob's account.

mesastar.org
------------

This is a long-standing domain name used by MESA.
The `https://mesastar.org <https://mesastar.org>`__ redirects to a landing page hosted on `GitHub Pages <https://github.com/MESAHub/MESAHub.github.io>`__

Domain Name
^^^^^^^^^^^

Matteo controls the domain name (using `domain.com <https://domain.com>`__) and handles
renewing it, updating nameservers, etc.

DNS Servers
^^^^^^^^^^^

The mesastar.org domain is configured to use `Cloudflare <https://cloudflare.com>`__ nameservers
and the DNS records are managed by Philip.

Mailing lists
^^^^^^^^^^^^^

The mesa-users@lists.mesastar.org and
mesa-developers@lists.mesastar.org lists are hosted by `MailmanLists
<https://www.mailmanlists.net/>`__. The account is controlled by
Matteo who pays the hosting fees through Flatiron.
This costs about $100/year.

New account signups are handled by Earl, Eb, and Warrick.

Marketplace
^^^^^^^^^^^

The `MESA marketplace
<https://mesastar.org/marketplace/>`__ is an archive of shared resources prior to 2022,
hosted on `GitHub Pages <https://github.com/mesahub/marketplace>`__.
For more recent contributions, see the `MESA Zenodo Community <https://zenodo.org/communities/mesa/>`__.


ReadTheDocs
-----------

The Sphinx documentation is hosted by `ReadTheDocs <https://readthedocs.org/>`__.
Rich controls the account (and Earl, Joey, Meridith, and Philip have access to it). This is
currently free for open source software.


Sourceforge
-----------

The original home of MESA development was sourceforge. Most developers
involved c. 2015 or earlier have admin privileges. This is a free service.

Website
^^^^^^^

We still use the mesa.sourceforge.net website domain name, but only for a top-level
redirect to docs.mesastar.org.

The source code for the old sourceforge site lives at
`<https://github.com/MESAHub/mesa-website>`__.

Slack
-----

Rich controls the Slack and pays for it (through grants at UW
Madison). This costs around $130 per year (about $10/user/yr).


TestHub
-------

The MESA TestHub runs on Heroku and the account is controlled by Bill
Wolf. It is funded by Frank and Bill W. The Heroku account costs
around $600/yr, though its cost has been higher during the GitHub
transition.

The testing `log archive <https://mesa-logs.flatironinstitute.org/>`__ lives on a
server at the Flatiron Institute, with Philip Mocz as a point of contact.

Jenkins
^^^^^^^

The Flatiron `Scientific Computing Core <https://www.simonsfoundation.org/flatiron/scientific-computing-core>`__
provides continuous integration testing for MESA via Jenkins.
This was set up in 2021 by Dylan Simon (Flatiron) and Josiah.
The configuration lives in the ``jenkins`` directory.
Matteo and Philip serve as the point of contact at Flatiron.
View the `Jenkins test results <https://jenkins.flatironinstitute.org/job/mesa/job/main/>`__


Zenodo
------

We have multiple resources archived on Zenodo. This is a free service.

Community
^^^^^^^^^

The `MESA Zenodo Community <https://zenodo.org/communities/mesa/>`__ on Zenodo
is controlled by Philip.

Records
^^^^^^^

* The `record for MESA releases <https://zenodo.org/records/18023257>`__ is controlled by Philip.
* The records for MESA SDK releases (`macOS <https://zenodo.org/records/13768941>`__, `linux <https://zenodo.org/records/13768913>`__) are controlled by Philip.
* The `record for OP Mono data <https://zenodo.org/record/4390522>`__ is controlled by Josiah.
