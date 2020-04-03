=================
Change Management
=================

This documentation assumes the reader is already familiary with using ``git``
for managing file revisions.

.. contents::
   :local:

Getting started
===============

#.  Go to https://gitlab.com/gromacs/gromacs

.. todo:: Link to GitLab docs for setting up an account, setting up git, and cloning repositories.

Uploading a commit for review
-----------------------------

Make sure your HEAD is up to date (use ``git pull --rebase origin`` if
someone else has committed since you last pulled), check that your commit
message follows the :doc:`commitstyle`.

Uploading a Work-In-Progress (WIP) or Private commit for review
---------------------------------------------------------------

.. todo:: Note work-in-progress and privacy semantics.

Code Review
===========

Reviewing someone else's uploaded commit
----------------------------------------

The reviewing workflow is the following:

#. https://gitlab.com/gromacs/gromacs/-/issues shows all open changes
#. A change needs a +2 and usually +1 review, as well as a +2 verified
   to be allowed to be merged.
#. Usually a patch goes through several cycles of voting, commenting and
   updating before it becomes merged, with votes from the developers indicating
   if they think that change hat progressed enough to be included.
#. A change is submitted for merging and post-submit testing
   by clicking "Submit" by one of the main developers. This should be done by
   the reviewer after voting +2. After a patch is submitted it is
   replicated to the main git server.

Do not review your own code. The point of the policy is that at least
two non-authors have voted +1, and that the issues are resolved in the
opinion of the person who applies a +2 before a merge. If you have
uploaded a minor fix to someone else's patch, use your judgement in
whether to vote on the patch +1.

Guide for reviewing
-------------------

-  First and foremost, check correctness to the extent possible;
-  As portability and performance are the most important things (after
   correctness) do check for potential issues;
-  Check adherence to the :ref:`GROMACS coding
   standards <style-guidelines>`;
-  We should try to ensure that commits that implement bugfixes (as
   well as important features and tasks) get an `issue tracker`_ entry created
   and linked. The linking is done **automatically** through
   `special syntax <https://gitlab.com/help/user/markdown#special-gitlab-references>`__
-  If the commit is a **bugfix**\ :

   -  if present in the `issue tracker`_, it has to contain a valid reference to the
      issue;
   -  if it's a **major bug**, there has to be a bug report filed in the
      `issue tracker`_  (with urgent or
      immediate priority) and referenced appropriately.

-  If the commit is a **feature/task** implementation:

   -  if it's present in the `issue tracker`_ it
      has to contain a valid reference to the issue;
   -  If no current issue is currently present and the change
      would benefit of one for future explanation on why it was
      added, a new issue should be created.

More git tips
=============

.. rubric:: Q: Are there some other useful git configuration settings?

A: If you need to work with
branches that have large
differences (in particular, if a
lot of files have moved), it can
be helpful to set

::

    git config diff.renamelimit 5000

to increase the limit of inexact
renames that Git considers. The
default value is not sufficient,
for example, if you need to do a
merge or a cherry-pick from
a release branch to master.

.. rubric:: Q: How do I use git rebase (also ``git pull --rebase``)?

A: Assume you have a local
feature branch checked out, that
it is based on master, and master
has gotten new commits. You can
then do

::

    git rebase master

to move your commits on top of
the newest commit in master. This
will save each commit you did,
and replay them on top of master.
If any commit results in
conflicts, you need to resolve
them as usual (including marking
them as resolved using git add),
and then use

::

    git rebase --continue

Note that unless you are sure
about what you are doing, you
should not use any commands that
create or delete commits (git
commit, or git checkout or git
reset without paths). ``git rebase
--continue`` will create the commit
after conflicts have been
resolved, with the original
commit message (you will get a
chance to edit it).

If you realize that the conflicts
are too messy to resolve (or that
you made a mistake that resulted
in messy conflicts), you can use

::

    git rebase --abort

to get back into the state you
started from (before the
original git rebase master
invocation). If the rebase is
already finished, and you realize
you made a mistake, you can get
back where you started with
(use git
log <my-branch>@{1} and/or git
reflog <my-branch> to check that
this is where you want to go)

::

    git reset --hard <my-branch>@{1}

.. rubric:: Q: How do I prepare several commits at once?

A: Assume I have multiple independent changes in my working tree.
Use

::

    git add [-p] [file]

to add one independent change at
a time to the index. Use

::

    git diff --cached

to check that the index contains
the changes you want. You can
then commit this one change:

::

    git commit

 If you want to test that the
change works, use to temporarily
store away other changes, and do
your testing.

::

    git stash

If the testing fails, you can
amend your existing commit with
``git commit --amend``. After you are
satisfied, you can push the
commit for review. If
you stashed away your changes and
you want the next change to be
reviewed independently, do

::

    git reset --hard HEAD^
    git stash pop

(only do this if you pushed the
previous change upstream,
otherwise it is difficult to get
the old changes back!) and repeat
until each independent change is
in its own commit. If you skip
the ``git reset --hard`` step, you
can also prepare a local feature
branch from your changes.

.. rubric:: Q: How do I edit an earlier commit?

A: If you want to edit the latest
commit, you can simply do the
changes and use

::

    git commit --amend

If you want to edit some other
commit, and commits after that
have not changed the same lines,
you can do the changes as usual
and use

::

    git commit --fixup <commit>

or

::

    git commit --squash <commit>

where <commit> is the commit you
want to change (the difference is
that ``--fixup`` keeps the original
commit message, while ``--squash``
allows you to input additional
notes and then edit the original
commit message during ``git rebase
-i``). You can do multiple commits
in this way. You can also mix
``--fixup/--squash`` commits with
normal commits. When you are
done, use

::

    git rebase -i --autosquash <base-branch>

to merge the ``--fixup/--squash``
commits to the commits they
amend. See separate question on
``git rebase -i`` on how to choose
<base-branch>.

In this kind of workflow, you
should try to avoid to change the
same lines in multiple commits
(except in ``--fixup/--squash``
commits), but if you have already
changed some lines and want to
edit an earlier commit, you can
use

::

    git rebase -i <base-branch>

but you likely need to resolve
some conflicts later. See ``git
rebase -i`` question later.

.. rubric:: Q: How do I split a commit?

A: The instructions below apply
to splitting the HEAD commit; see
above how to use ``git rebase -i`` to
get an earlier commit as HEAD to
split it.

The simplest case is if you want
to split a commit A into a chain
A'-B-C, where A' is the first new
commit, and contains most of the
original commit, including the
commit message. Then you can do

::

    git reset -p HEAD^ [-- <paths>]
    git commit --amend

to selectively remove parts from
commit A, but leave them in your
working tree. Then you can create
one or more commits of the
remaining changes as described in
other tips.

If you want to split a commit A
into a chain where the original
commit message is reused for
something else than the first
commit (e.g., B-A'-C), then you
can do

::

    git reset HEAD^

to remove the HEAD commit, but
leave everything in your working
tree. Then you can create your
commits as described in other
tips. When you come to a point
where you want to reuse the
original commit message, you can
use

::

    git reflog

to find how to refer to your
original commit as ``HEAD@{n}``, and
then do

::

    git commit -c HEAD@{n}

.. rubric:: Q: How do I use git rebase -i to only edit local commits?

A: Assume that you have a local
feature branch checked out, this
branch has three commits, and
that it is based on master.
Further, assume that master has
gotten a few more commits after
you branched off. If you want to
use ``git rebase -i`` to edit your
feature branch (see above), you
probably want to do

::

    git rebase -i HEAD~3

followed by a separate

::

    git rebase master

The first command allows you to
edit your local branch without
getting conflicts from changes in
master. The latter allows you to
resolve those conflicts in a
separate rebase run. If you feel
brave enough, you can also do
both at the same time using

::

    git rebase -i master
