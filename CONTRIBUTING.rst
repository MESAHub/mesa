============
Contributing
============

MESA welcomes contributions from the user community.  Please note that
contributing code to MESA (e.g., acceptance of a pull request) does
not entitle one to co-authorship of a MESA instrument paper.


Bug Reports
-----------

We welcome bug reports about MESA.  However, bug reports must describe
specific, actionable issues.  An undiagnosed difficulty with a MESA
model should not be treated as a bug.

If you are encountering problems with your MESA models, user support
occurs through the mesa-users@lists.mesastar.org mailing list.  Send a
message describing the problem and including enough detail (e.g., MESA
version, inlists, saved models) such that anyone can reproduce it.
Other users and the MESA developers can help you understand if your
problem represents a bug.


Code Changes / Pull Requests
----------------------------

Proposed changes will be carefully reviewed by the MESA developers.
Simple bugfixes or improvements to documentation are likely to be
readily accepted.  Modifications that can be easily realized via hooks
may be more appropriate for `mesa-contrib <https://github.com/MESAHub/mesa-contrib>`__
and the MESA developers may suggest adding your changes there.

User contributions are expected to follow the
`standard technical guidelines for MESA developers <https://docs.mesastar.org/en/latest/developing.html>`__.
We strongly recommend against undertaking any significant development
effort without first emailing mesa-developers@lists.mesastar.org with
an outline of your anticipated changes.

New contributions must be proposed through GitHub's *pull request*
(PR) system.  The process is roughly:

1. fork the ``mesa`` repo (click *Fork* in the top-right of the GitHub interface),
2. clone your fork to your computer,
3. create a new branch for your additions (e.g. ``git switch -c my-new-hook``),
4. make, commit and push your changes and
5. open a PR against the ``main`` branch.

Once the PR is opened, the developers will review your code and make a
decision about whether to merge it into the repo, perhaps after some
changes are made.

Becoming a developer
--------------------

If you simply wish to fix a bug or add new code to MESA then you do not need to become a developer.
However, if you wish to take a more active role in managing the project then please see the instructions :doc:`here <developing/new_developers>`.
