=================
How to contribute
=================

.. warning::

   This document is written for MESA developers.  It currently targets
   the sandbox repository used for testing the migration to GitHub.
   Changes to this repository will not be included in MESA.

Obtaining MESA
==============

Join the MESA team on GitHub
----------------------------

You must create an account on `GitHub <https://github.com/>`_.  The MESA organization is
called `MESAHub <https://github.com/MESAHub/>`_.  After you are registered, you can be invited.


Set up Git
----------

Follow `this GitHub guide <https://help.github.com/en/github/getting-started-with-github/set-up-git>`_ to set up git.

When instructed to install git, it is likely simplest to do so using your system package manager on Linux (e.g., apt, yum) or a macOS package manager (e.g., Homebrew, MacPorts).

.. note::
   The documentation will describe command line use of git.  However, there are many graphical git clients available.  For example, GitHub has its own `desktop app <https://desktop.github.com/>`_.

   
Install Git LFS
---------------

MESA uses Large File Storage (LFS), a git extension that allows us to
easily store and version the large microphysics data files.

Obtain and install Git LFS from `its website <https://git-lfs.github.com/>`_.  

.. note::

   The Git LFS home page has additional information about configuring a repository.  This has already been done for MESA.  The only command you need to run is ``git lfs install``.

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
``master`` branch is the development version of MESA.



Making a commit
---------------

When you want to add a new file or add changes to an existing file, do::

  git add <filename>

Note that this is different from ``svn add``, which is only used to
track new files.

Once you have added your desired changes::

  git commit -m <message>

If you want to commit all your changes to existing files, you can skip the git add and do::

  git commit -a -m <message>

Note that this is different than ``svn commit``, as it only alters
your local repository.  It does not transmit your changes to the
remote repository.

Sharing your changes
--------------------

Once you are ready to share your changes::

  git push origin master

Fetching others changes
-----------------------

When you want to get others' changes::

  git pull origin master

Checking out a revision
-----------------------

To checkout out a previous version of the repository you first need the hash of the commit.
This is a 40 character code looking like ``37cbee26a70574189d2e6169594af360215e18b6``, luckily though you do not need the full 40
characters, you only need the enough of the hash to be unique. This is usually only 6-8 characters long::

  git checkout 37cbee26

This will return your local repository to the state is was in the commit given by ``37cbee26``, but with your current uncommitted changes
on top.

git checkout can also be used to checkout the code based on other identifiers::

  git checkout path/to/file

Resets the file to its previous committed state::

  git checkout mybranch

Will check out the branch ``mybranch``


Branching
---------

If you decided to make a new branch this can be done with::

  git branch mynewbranch
  git checkout mynewbranch

or::

  git checkout -b mynewbranch

Any changes you now make will not apply to ``master`` but instead to ``mynewbranch``

When you are ready to merge the changes from ``mynewbranch`` into ``master`` then you do

.. note::
  Decide on merging strategy


To delete the branch::

  git branch -D mynewbranch


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
