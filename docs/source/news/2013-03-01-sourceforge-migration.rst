=============================
SourceForge Repository Change
=============================

:Date:   2013-03-01

SourceForge.net switched to a new platform in late 2012 and required all
projects to migrate to it.

The most significant change for typical users will be that the MESA svn
repository has a new URL. The MESA project webpage “getting started”
section has been updated to reflect the change.

The simplest way to keep up to date with MESA is to check out a new copy
from the new repository; it’s as simple as that.

If you maintain a working copy of MESA under version control for
development purposes, you’ll need to do something akin to “svn relocate”
from the old repository to the new. If you have svn version 1.6 or
lower, then use “svn switch –relocate” instead. See ‘svn help relocate’
(v1.7 and up) or ‘svn help switch’ (v1.6 and lower).

Don’t just copy/paste the below commands. First, run ‘svn info’ in your
MESA directory to make sure that you have the correct URL for the old
repository.

(svn 1.7 and above)

::

   svn relocate http://mesa.svn.sourceforge.net/svnroot/mesa/trunk svn://svn.code.sf.net/p/mesa/code/trunk

(svn 1.6 and below)

::

   svn switch --relocate http://mesa.svn.sourceforge.net/svnroot/mesa/trunk svn://svn.code.sf.net/p/mesa/code/trunk

Finally, if you use a version control program other than svn, see the
`sourceforge advice
page <http://sourceforge.net/p/forge/community-docs/Repository%20Upgrade%20FAQ/>`__.

That should be all you need to know. If you have any questions about
this process, please send an email to mesa-users.
