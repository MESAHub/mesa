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

Adding:
::

  fetch = +refs/pull/*/head:refs/remotes/origin/pr/*

To your ``.git/config`` in the ``[remote "origin"]`` section enables fetching pull requests

::

  $(git log -1) == *'[ci skip]'*

Tests to see if we should skip testing the test cases. Note we still want to compile test MESA even if we
skip the full test suite.
