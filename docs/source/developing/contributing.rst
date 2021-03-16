=================
How to contribute
=================

.. warning::

   This document is written for MESA developers.  Users are not yet
   able to contribute via GitHub.

Obtaining MESA
==============

Join the MESA team on GitHub
----------------------------

You must create an account on `GitHub <https://github.com/>`_.  The
MESA organization is called `MESAHub <https://github.com/MESAHub/>`_.
After you are registered on GitHub, you can be invited to the
organization.


Set up Git
----------

Follow `this GitHub guide <https://help.github.com/en/github/getting-started-with-github/set-up-git>`_ to set up git.

When instructed to install git, it is likely simplest to do so using your system package manager on Linux (e.g., apt, yum) or a macOS package manager (e.g., Homebrew, MacPorts).

.. note::
   This documentation will mainly describe command line use of git.  However, there are many graphical git clients available.  For example, GitHub has its own `desktop app <https://desktop.github.com/>`_.  (See also :ref:`Graphical interfaces`).

   
Install Git LFS
---------------

MESA uses Large File Storage (LFS), a git extension that allows us to
easily store and version the large microphysics data files.

Obtain and install Git LFS from `its website <https://git-lfs.github.com/>`_.  You only need to do Step 1 in the Getting Started section.  Step 2 has already been done for MESA.


Download the MESA repository
----------------------------

Make a git clone of the MESA repository. (This is analogous to ``svn
checkout``.) This will automatically fetch the files from git LFS.  To
clone over HTTPS, do::

    git clone https://github.com/MESAHub/mesa.git mesa

This will create the MESA repository in a directory named ``mesa``.



Making changes to MESA
======================

.. note::

   Git is a powerful and complex tool.  This section is an extremely
   incomplete introduction to git on the command line.  There are many
   online resources at the tutorial level (e.g., `git
   <https://git-scm.com/docs/gittutorial>`_, from `GitHub
   <https://guides.github.com/introduction/git-handbook/>`_, from
   `Bitbucket <https://www.atlassian.com/git/tutorials>`_) and beyond.

The conventional name of the remote repository is ``origin``.  The
``main`` branch is the development version of MESA.

.. note::

   Prior to 2020, the default branch of a git repository was often
   canonically named ``master``.  You may still see this used in some
   older materials.  The MESA default branch is ``main``.

MESA development uses the following branching model:

* **main** --- Default branch.  This is equivalent of SVN trunk.  Small/innocent changes can directly be made in ``main``. The ``main`` branch should always be passing the test suite.

* **feature branches** --- Branched from ``main`` for large/significant changes. Feature branches should be pushed to GitHub and then merged via a pull request.  Confirm that the test suite passes before merging.  This also gives an opportunity for other developers to comment on the changes.

* **release branches** --- Branched from main at the time of the first release candidate.  Any final fixes can be backported (cherry-picked or merged) to ``main``.

The line between what should be a feature branch / pull request and what's committed straight onto the main branch can be drawn based on how disruptive a change is expected to be. If it has the potential to break test suite cases, don't commit it directly to main.  Use a branch.


Making a commit
---------------

When you want to add a new file or add changes to a file that is already tracked by git, do::

  git add <filename>

Note that this is different from ``svn add``, which is only used to
track new files.

Once you have added your desired changes::

  git commit -m <message>

If you want to commit all your changes to files already tracked by git, you can skip the git add and do::

  git commit -a -m <message>

Note that this is different than ``svn commit``, as it only alters
your local repository.  It does not transmit your changes to the
remote repository.

Sharing your changes
--------------------

Once you are ready to share your changes::

  git push origin main

The first argument is the remote repository (``origin`` = GitHub).
The second argument is the branch name.  If you are making changes in
the non-default branch (i.e., not ``main``), see :ref:`Branching`.
  
Fetching others changes
-----------------------

When you want to get others' changes::

  git pull origin main

The first argument is the remote repository (``origin`` = GitHub).
The second argument is the branch name.  If you want to get changes
from a non-default branch (i.e., not ``main``), see :ref:`Branching`.

If you have made changes to a branch and pull from a branch that
others have independently modified, git must decide how to reconcile
the divergent branches.  In order to avoid unnecessary merges (thereby
giving a simpler, more linear commit history), it is suggested that
you default to rebasing your changes.  Rebasing means that git will
pull others' changes and then replay your changes on top of them.  If
the changes conflict, you will have a chance to resolve the conflicts.
To make this the default behavior, issue the following command while
your current working directory is in the MESA git repository::

  git config pull.rebase true

You can also pass ``--rebase``, ``--no-rebase`` (combine changes with
a merge commit), or ``--ff-only`` (refuse to pull if there are other
changes) on the command line to override the configured default per
invocation.


If you want to get others' changes, but not immediately update your
local repository to match that content::

  git fetch origin


Checking out a revision
-----------------------

To checkout out a previous version of the repository you first need the hash of the commit.
This is a 40 character code looking like ``37cbee26a70574189d2e6169594af360215e18b6``, luckily though you do not need the full 40
characters, you only need the enough of the hash to be unique. This is usually only 6-8 characters long::

  git checkout 37cbee26

This will return your local repository to the state is was in the commit given by ``37cbee26``, but with your current uncommitted changes
on top.

  
Restoring a file
----------------

.. note::

  Recent versions of git include the new command ``git restore`` that
  splits out some of the functionality of ``git checkout``.  (If you
  already know how ``git checkout`` works, you can also use that
  command to accomplish the same goal.)

If you want to reset a file to its most recently committed state::

  git restore path/to/file

The ``path/to/file`` could also be something like the current
directory (``.``) or a list of files (``*.f90``).


``git restore`` can also be used to restore a file from another commit::

  git restore --source=37cbee26 path/to/file
  


Branching
---------

.. note::

  Recent versions of git include the new command ``git switch`` that
  splits out some of the functionality of ``git checkout``.  (If you
  already know how ``git checkout`` works, you can also use that
  command to accomplish the same goal.)


If you decided to make a new branch this can be done with::

  git branch mynewbranch
  git switch mynewbranch

or::

  git switch -b mynewbranch

Any changes you now make will not apply to ``main`` but instead to ``mynewbranch``.

To delete the branch::

  git branch -D mynewbranch

If you want to push that branch to GitHub to share it with others, do::

  git push --set-upstream origin mynewbranch

This will create a new branch on GitHub named ``mynewbranch`` and associate it with the local branch on your machine of the same name.

.. note::

  Give the branch a short, descriptive name.  To help others quickly
  see who a branch belongs to, you can prepend your initials (e.g.,
  ``jws/kap-compton`` or ``rf/rates-nullify``).

Once you have set the upstream branch, you may omit the branch name when you push additional changes to this branch::

  git push origin

or pull additional changes from others on this branch::

  git pull origin


When you are ready to merge the changes from ``mynewbranch`` into ``main`` then push ``mynewbranch`` to GitHub and make a pull request.


If someone else has created a new branch, you can switch to it with::
  
    git switch -c theirnewbranch origin/theirnewbranch

This will check out the branch ``theirnewbranch`` and associate it with the remote branch.


Stashing changes
----------------

Lets say you are working on the code and suddenly a bug report comes in and you decide to fix that code first before you finish your current work.
Because your initial work is still in progress you want to save it but do not want to commit it yet. This is where 
git stash comes in::

  git stash

This saves your current changes that have not been committed in a ``stash`` and resets your repository to the
current committed version. You can then make your changes to fix the new bug then re-apply the stash on top of the new
commit::

  git stash apply

This way your in progress changes do not get mixed in with unrelated changes. Note the ``stash`` still exists, so you need to drop
it once you no longer need it::

  git stash drop

You may have multiple stashes at once, in which case they are indexed by:: 

  git stash ${X}

where ``X`` is a number starting at 0 for the most recent ``stash``.


Graphical interfaces
--------------------

Not everything needs to be done by command line. There are at least two GUI tools that are usually shipped with git,
git gui and gitk. 
::

  git gui

This provides a convenient interface for making commits. You can select which files to commit, which lines of which
files, set the commit message, and make the commit. 
::

  gitk

This provides a convenient interface for viewing the history of the repository you can view the commits, files changed, and commit messages.
::

  gitk --all

By default ``gitk`` only shows the current branch ``--all`` shows all branches.
::

  gitk path/to/file

Will only show the commits that change the file.

Git testing tips
----------------

::

  git fetch --all

Fetches all commits over all branches

::

  $(git log -1) == *'[ci skip]'*

Tests to see if we should skip testing the test cases. Note we still want to compile test MESA even if we
skip the full test suite.


Pull requests
=============

Preparing to make a PR
----------------------

After you have made a branch and pushed it to GitHub (see
:ref:`Branching`), the test suite will automatically be run and the
results reported to MESA TestHub.  You should wait for the TestHub to
confirm that the test suite is passing before merging a PR.

If the changes in your branch conflict or interact with changes that
have occurred on ``main``, it is recommended that you merge ``main`` into your
branch (or rebase your branch to the tip of ``main``) before issuing the
PR. This allows you to handle conflicts in advance and ensure that the
test suite will re``main`` passing after you merge your branch back into
``main``.

     

Making a pull request
---------------------

After you have made a branch and pushed it to GitHub (see
:ref:`Branching`), you can issue a pull request for the code on your
branch to be merged into ``main``.

If you have recently pushed a branch, GitHub will offer you the option to make a PR on the `main page <https://github.com/MESAHub/mesa>`_.  Otherwise, the most general approach is to visit the `new pull request page <https://github.com/MESAHub/mesa/compare>`_, select the code you want to merge from the 'compare' dropdown, and then click the green 'Create pull request' button.  You will be asked to provide a title and description for the PR as well as other optional information (like selecting a reviewer).  Then click 'Create pull request'.  Once you have made the PR, it will show up in the `list of pull requests <https://github.com/MESAHub/mesa/pulls>`_.

A set of code reviewers is automatically selected for each PR based on the contents of the ``CODEOWNERS`` file.  For now, this request for review can be treated as a heads up that there are changes in a part of the code you may be interested in.  Reviewers are not required to complete requested reviews and reviews are not required before a PR is merged from a MESA developer.  However, please exercise good judgment and solicit feedback before merging, especially for significant changes or changes that you feel uncertain about.  You may want to ping relevant individuals or channels in Slack.


Merging a pull request
----------------------

Once the code is ready, it can be merged by visiting the page associated with the PR (e.g., `<https://github.com/MESAHub/mesa/pull/161>`_).

GitHub offers several strategies for merging pull requests.  Each one may be appropriate in different circumstances.  The merge strategy is selected by using the dropdown arrow on the big green button at the bottom of the PR.

* If the PR is a small set of simple, well-contained changes, the 'Rebase and merge' strategy is recommended.  This will take the commits and add them to the tip of ``main``, ensuring that the commit history of ``main`` remains linear.


* If the PR is a set of changes whose detailed history is not relevant, the 'Squash and merge' strategy is recommended.  This will take the commits, combine them into a single commit, and then add it to the tip of ``main``. This stragegy is most useful when the series of individual commits simply reflects the (possibly wandering) path to achieving the final cumulative change.


* If the PR is a set of changes where each commit is a meaningful, self-contained change, but the cumulative change is not simple enough for the 'Rebase and merge' strategy, then the 'Create a merge commit' strategy is appropriate.  This will preserve the full history of your branch when it is joined with ``main``.  If a change has this level of complexity, it is also recommended that its interaction with ``main`` should be tested by merging ``main`` into the branch.
