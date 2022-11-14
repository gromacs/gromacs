=================
Change Management
=================

This documentation assumes the reader is already familiar with using ``git``
for managing file revisions.

Getting started
===============

GROMACS development happens on gitlab at https://gitlab.com/gromacs/gromacs.
Create a user account at https://gitlab.com/users/sign_in#register-pane or use
an existing account at gitlab.com. For more information on how to use gitlab have
a look at their extensive user documentation at https://docs.gitlab.com/ee/user/index.html.
We follow the workflow described in https://docs.gitlab.com/ee/topics/gitlab_flow.html. 

If you do not already have a GROMACS repository set up, user 
``git clone git@gitlab.com:gromacs/gromacs.git`` to obtain the current GROMACS
repository from gitlab. Otherwise use 
``git remote add gitlab git@gitlab.com:gromacs/gromacs.git``. 

Using gitlab, new code enters GROMACS by merging git development branches into
the main branch. 

To automatically detect issues in new code, it is tested within continuous
integration (CI) with a large combination of settings.
See `gmx-codeformatting` for help meeting and testing the style guidelines.

Setting up login credentials with gitlab
----------------------------------------

You will need a public ssh key::

    ssh-keygen -t rsa -C "your.email@address.com"
    cat ~/.ssh/id_rsa.pub

Copy the output of the last command, got to gitlab.com, find you user in the
right top corner and select settings.

Chose SSH keys in the menu on the left and past your key in the text field.

Creating issues
---------------

The meta-level code design and discussions is organised in issues and visible at
https://gitlab.com/gromacs/gromacs/-/issues. Please check if if your issue or a
similar issue already exists before creating a new one.

Note that all Redmine issues have been transferred to gitlab with the same issue
numbers as used in gitlab. However, comments and discussion are now represented
by gitlab user @acmnpv - the original authors are found inline at the bottom of
the comments. 

Uploading code for review - creating a merge request
----------------------------------------------------

Issues are addressed with new code via "merge requests" (MR). Find the current
MRs at https://gitlab.com/gromacs/gromacs/-/merge_requests. 
There are two ways of creating a merge request - either via the gitlab graphical
user interface or via the command line. 

To use the GUI, find the relevant issue or open a new one, then find the 
"create merge request" button to create a merge request related to that issue in gitlab.
The default selection is to mark this a work in progress (WIP) merge-request.
We recommend keeping this setting until you are completely satisfied with the 
code yourself and all tests are passed.

Select milestone and assignees to make tracking of the progress easier. 
Keep the requirements for merging as they are set by default.

You can also use ``git push`` on the command line directly and create a merge request 
following the link that is output on the command line.

Your repository should be in sync with the GROMACS repository. To ensure this,
use ``git fetch`` to obtain the newest branches, then merge the main branch
into your branch with ``git merge main`` while on your branch.

Naming branches
---------------

Good names: documentation_UpdateDevelopersDocsTOGitLab, nbnxm_MakeNbnxmGPUIntoClass, pme_FEPPMEGPU. 
Bad names: branch1234, mybranch, test, etc

Documentation
-------------

Contributors and reviewers frequently overlook the effects of changes on the built documentation.
Contributors and reviewers should note that the build artifacts from the automated test jobs
are available for download through the GitLab CI web interface (``webpage:build`` job artifacts).
For earlier review or alternative preferences, consider building and sharing a Docker image
containing the built documentation. See
`docs/docs.dockerfile <https://gitlab.com/gromacs/gromacs/-/tree/main/docs/docs.dockerfile>`__
in the source tree.

Labels
======

`Labels <https://docs.gitlab.com/ee/user/project/labels.html>`__
help developers by allowing additional filtering of issues and merge requests.

The GROMACS project `defines many labels <https://gitlab.com/gromacs/gromacs/-/labels>`__.

.. Note: labeling guidelines TBD. See https://gitlab.com/gromacs/gromacs/-/issues/3949 and open new issues as appropriate.

To minimize duplicated documentation, refer to the
`GROMACS Labels <https://gitlab.com/gromacs/gromacs/-/labels>`__ web interface for label descriptions.

When creating a new label, please provide a short description
so that people can understand what the label is intended to convey,
and when they should apply it to their own issues or merge requests.

In general:

* Ongoing categorizations to help specify the GROMACS component or development area use the ``#7F8C8D`` color.
* Specific features or subproject areas targeting an upcoming release use the ``#8E44AD`` background color.
* Status labels use ``#428BCA``. Note that Status labels are also used for Issues,
  and are used according to
  :ref:`status label guidelines <status label guidelines>`

.. Best practices and labeling policies can be proposed as changes to this document. See https://gitlab.com/gromacs/gromacs/-/issues/3949

Code Review
===========

Reviewing someone else's uploaded code
--------------------------------------

The reviewing workflow is the following:

#. https://gitlab.com/gromacs/gromacs/-/merge_requests shows all open changes
#. A change needs two approvals to go in, of which one approval has to come from
   a member of either GMX Core or GMX Developers.
#. Usually a patch goes through several cycles of voting, commenting and
   updating before it becomes merged, with votes from the developers indicating
   if they think that change hat progressed enough to be included.
#. A change is submitted for merging and post-submit testing
   by clicking "Merge".

Do not review your own code. The point of the policy is that at least
two non-authors have approved, and that the issues are resolved in the
opinion of the person who applies an approval before a merge. If you have
uploaded a minor fix to someone else's patch, use your judgement in
whether to approve yourself.

Guide for reviewing
-------------------

-  First and foremost, check correctness to the extent possible;
-  As portability and performance are the next most important things do check 
   for potential issues;
-  Check adherence to the :ref:`GROMACS coding standards <style-guidelines>`;
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

.. _status label guidelines:

Update the Status label
"""""""""""""""""""""""

-  Please update the Status label :ref:`for the issue <issue workflow>` when a merge request is under review.
-  Please update the Status label :ref:`for the merge request <merge request status>` when it is closed.

.. _merge request status:

Closing Merge Requests
----------------------

A merge request that has had no updates for six months or more can acquire the status label "Status::Stale"
If the proposed change still seems important and the next steps are unclear,
contributors with stale issues *are encouraged...*

- to contact existing reviewers (or potential reviewers),
- to participate in the `developer discussion forum`_, and
- to attend the biweekly teleconference to coordinate.

If the future of the merge request has not become clear within a month
(especially if it has become stale multiple times),
developers may close the merge request with a label indicating why it has entered a "closed" state.
`"Status::MR::..." labels <https://gitlab.com/gromacs/gromacs/-/labels?subscribed=&search=status%3A%3Amr>`__
do not indicate that the merge request has been reviewed
unless it is explicitly rejected.

See :issue:`4126` for background discussion.

- `Status::MR::Inactive <https://gitlab.com/gromacs/gromacs/-/merge_requests?label_name%5B%5D=Status%3A%3AMR%3A%3AInactive>`__: No response from contributor or no reviewers available for over six months.
- `Status::MR::Superseded <https://gitlab.com/gromacs/gromacs/-/merge_requests?label_name%5B%5D=Status%3A%3AMR%3A%3ASuperseded>`__: This merge request is no longer necessary.
- `Status::MR::Rejected <https://gitlab.com/gromacs/gromacs/-/merge_requests?label_name%5B%5D=Status%3A%3AMR%3A%3ARejected>`__: The solution (or its associated issue) will not be accepted.
- `Status::MR::Needs discussion <https://gitlab.com/gromacs/gromacs/-/merge_requests?label_name%5B%5D=Status%3A%3AMR%3A%3ANeeds+discussion>`__: More discussion must take place at the tracked issue before a MR is opened.
- `Status::Stale <https://gitlab.com/gromacs/gromacs/-/labels?subscribed=&search=status%3A%3AStale>`__: No activity for over six months.

.. seealso:: :ref:`issue workflow` for use of Status labels in Issue management.

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
a release branch to main.

.. rubric:: Q: How do I use git rebase (also ``git pull --rebase``)?

A: Assume you have a local
feature branch checked out, that
it is based on main, and main
has gotten new commits. You can
then do

::

    git rebase main

to move your commits on top of
the newest commit in main. This
will save each commit you did,
and replay them on top of main.
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
original git rebase main
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
that it is based on main.
Further, assume that main has
gotten a few more commits after
you branched off. If you want to
use ``git rebase -i`` to edit your
feature branch (see above), you
probably want to do

::

    git rebase -i HEAD~3

followed by a separate

::

    git rebase main

The first command allows you to
edit your local branch without
getting conflicts from changes in
main. The latter allows you to
resolve those conflicts in a
separate rebase run. If you feel
brave enough, you can also do
both at the same time using

::

    git rebase -i main
