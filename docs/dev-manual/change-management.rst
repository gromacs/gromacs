.. _gmx-gerrit:

=========================
GROMACS change management
=========================

This documentation assumes the reader is already familiary with using ``git``
for managing file revisions.

.. contents::
   :local:

Getting started
===============

#.  Go to https://gerrit.gromacs.org
#.  Click Register (you can choose any OpenID provider including any
    existing Google/Yahoo account. If you manually enter the URL make sure
    to start with ``http(s)://``)
#.  Choose a username and add an ssh key

See `here <https://gerrit.gromacs.org/Documentation/intro-quick.html>`_ for
a quick intro into Gerrit.

Creating the SSH key for Gerrit
-------------------------------

In order to push your commits to gerrit server, you must have an SSH key
in your computer which matches with the one registered in your Gerrit
user account. To do so, you first need to create this unique SSH
key. You will be asked to enter a passphrase. *This is
optional with respect to Gerrit, but it is a good security practice to have
it.*

To proceed with the creation of the SSH key, type the following commands
from your terminal window:

::

    $ cd ~/.ssh

    $ ssh-keygen -t rsa -C "your.email@address.com"

Please substitute the email string in the command above
with the same email address which you used to register the account in
Gerrit.

Now you have created your public SSH key, which you need to copy/paste
into your Gerrit profile. First, open it with the following command:

::

    $ cat id_rsa.pub

Copy all the contents of the file id_rsa.pub in your clipboard, and
switch to your favorite web browser where you logged in to Gerrit
GROMACS page. Click on your username at the top right corner of the
Gerrit webpage and select "Settings". You should now be in your Gerrit
profile settings page, where you should see a vertical menu.

From this vertical menu, select "SSH Public Keys", then click the button
"Add Key ..." and an edit box will appear below the button. Here you
need to paste the contents of id_rsa.pub file, which you previously
copied to your clipboard.

Now you are ready to operate!

Setting up a local repository to work with gerrit
-------------------------------------------------

Either clone using::

    $ git clone ssh://USER@gerrit.gromacs.org/gromacs.git

(replace **USER** \ with your username)

or change the remote url using:

::

    $ git remote set-url origin ssh://USER@gerrit.gromacs.org/gromacs.git

(change **USER** with the username you've registered)

Or add a new remote url using:

::

    $ git remote add upload ssh://USER@gerrit.gromacs.org/gromacs.git

If you are working with a GROMACS repository other than the source code,
then you should substitute e.g. regressiontests.git or releng.git
instead of gromacs.git above.

Be sure to configure your user name and e-mail to match those registered to Gerrit::

       git config [--global] user.name "Your Name"
       git config [--global] user.email "your.name@domain.org"

It is optional if you want to set those settings for git on a global
level, or just for the current repository.

If necessary, register the e-mail address you want to use
with Gerrit.

Install the commit hook
-----------------------

Differently from a simple usage of git, with Gerrit a Change-ID is
needed at the end of each commit message. Gerrit uses Change-IDs to
understand whether your new commit is patching a previous commit or it
should be regarded as a separate, different patch, uncorrelated with
your previously pushed commits.

To allow git to append such Change-IDs automatically after each commit,
type the following command:

::

    $ scp -p USER@gerrit.gromacs.org:hooks/commit-msg .git/hooks/

(change **USER** with the username you've registered in Gerrit)

.. Note::

   This commit hook needs to be added to the repo where the
   commit will occur, not the repo where the push to upstream will occur
   (should they be different).

Uploading a commit for review
-----------------------------

Make sure your HEAD is up to date (use ``git pull --rebase origin`` if
someone else has committed since you last pulled), check that your commit
message follows the :doc:`commitstyle`, make your commit and then use

::

    $ git push origin HEAD:refs/for/BRANCH

Replace ``BRANCH`` with the branch it should be committed to.
Master has a number of sub branches that can be used to show
what the patch is relevant to such as OpenCL and tools-cleanup.
These can be pushed to by specifying them after the branch,
for example ``BRANCH/domdec-cleanup``.

When updating/replacing an existing change, make sure the commit message
has the same Change-ID. Please see the section `Ammending a change <gmx-ammend-change>`
below.

Uploading a Work-In-Progress (WIP) or Private commit for review
---------------------------------------------------------------

You can use the WIP or Private workflow on Gerrit to upload changes
that might not be ready yet for public review and merging.
Those changes will only be visible to people explicitly added as reviewers,
and will not automatically trigger Jenkins if the reviewer "Jenkins Buildbot"
is not added manually to them.

For uploading a new private change, push to refs/for/master%private
(substituting master with the branch you want to push to). To remove the private
flag when uploading a new patch set, use refs/for/master%remove-private.
To mark change as Work-In-Progress, push to refs/for/master%wip,
to unmark push to refs/for/master%ready.
You can also mark and unmark changes as Private or WIP in the Gerrit web-interface.

To manually trigger Jenkins on a WIP or Private change, you need to log in
to Jenkis after adding the "Jenkins Buildbot" reviewer. In Jenkins, navigate to
http://jenkins.gromacs.org/gerrit_manual_trigger/ and tell it to
search for the commit for which you want to trigger the build agents.
For example, https://gerrit.gromacs.org/#/c/1238/ is 1238 (but maybe
SHA or ChangeID will work, too).
Any change made to the commit after "Jenkins Buildbot" was added to the
list of reviewers will also trigger Jenkins.

After uploading a commit
------------------------

Use

::

    $ git reset --keep HEAD^

to reset your branch to the HEAD before the commit you just uploaded.
This allows you to keep your repo in sync with what every other repo
thinks is the HEAD. In particular, if you have another patch to upload
(or worse, have to pull in other people's patches, and then have a new
patch), you probably do not want to have the second patch depend on the
first one. If the first one is rejected, you have made extra work for
yourself sorting out the mess. Your repo still knows about the commit,
and you can cherry-pick it to somewhere if you want to use it.

Code Review
===========

Reviewing someone else's uploaded commit
----------------------------------------

The reviewing workflow is the following:

#. https://gerrit.gromacs.org/#q/status:open shows all open changes
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
   well as important features and tasks) get a `Redmine`_ entry created
   and linked. The linking is done **automatically** by
   `Redmine`_ **if the commit message contains** keyword
   "#issueID", the valid syntax is explained below.
-  If the commit is a **bugfix**\ :

   -  if present in Redmine it has to contain a valid reference to the
      issue;
   -  if it's a **major bug**, there has to be a bug report filed in
      `Redmine`_  (with urgent or
      immediate priority) and referenced appropriately.

-  If the commit is a **feature/task** implementation:

   -  if it's present in `Redmine`_ it
      has to contain a valid reference to the issue;
   -  If no current issue is currently present and the change
      would benefit of one for future explanation on why it was
      added, a new redmine issue should be created.

Use of Verify
-------------

Jenkins has been installed for automated build testing. So it isn't
required to vote "verify +2" anymore. As the testing is not always
perfect, and because test coverage can be spotty, developers can still
manually vote to indicate that a change performs as intended. Please note
that this should not be abused to bypass Jenkins testing. The vote from
the test suite should only be discarded if failures are caused by unrelated
issues.

Further information
-------------------

Currently it is possible to review your own code. It is undesirable to
review your own code, because that defeats the point. It will be
deactivated if it is being abused and those responsible may lose
their voting rights.

For further documentation:

-  |Gromacs| `specific manual <https://gerrit.gromacs.org/Documentation/index.html>`__
-  `General tutorials <https://gerrit-documentation.storage.googleapis.com/Documentation/2.15.3/index.html#_tutorials>`__

FAQs
====

How do I access gerrit behind a proxy?
--------------------------------------

If you are behind a firewall blocking port 22, you can use socat to
overcome this problem by adding the following block to your
``~/.ssh/config``

::

    Host gerrit.gromacs.org
           User USER
           Hostname gerrit.gromacs.org
           ProxyCommand socat - PROXY:YOURPROXY:gerrit.gromacs.org,proxyport=PORT

Replace ``YOURPROXY``, ``PORT`` and ``USER``, (but not ``PROXY``!) with your own
settings.

How do I link fixes with Redmine issues?
----------------------------------------

The linking of commits that relate to an existing issue is
done automatically by `Redmine`_ if
the git commit message contains a reference to the Redmine entry
through the issueID, the numeric ID of the respective issue (bug,
feature, task). The general syntax of a git comit reference is [keyword]
#issueID.

The following two types of refereces are possible:

-  For bugfix commits the issueID should be preceeded by
   the "Fixes" keyword;
-  For commits related to a general issue (e.g. partial implementation of
   feature or partial fix), the issueID should be preceeded by the "Refs" keyword;

An example commit message header::

    This commit refs #1, #2 and fixes #3

How can I submit conflicting changes?
-------------------------------------

When there are several, mutually conflicting changes in gerrit pending
for review, the submission of the 2nd and subsequent ones will fail.
Those need to be resolved locally and updated by

::

    $ git pull --rebase

Then fix the conflicts and use

::

    $ git push

Please add a comment (review without voting) saying that it was rebased
with/without conflicts, to help the reviewer.


.. _gmx-ammend-change:

How do I upload an update to a pending change?
----------------------------------------------

First, obtain the code you want to update. If you haven't changed your
local repository, then you already have it. Maybe you can check out the
branch again, or consult your git reflog. Otherwise, you should go to
gerrit, select the latest patch set (remembering that others may have
contributed to your work), and use the "Download" link to give you a
"Checkout" command that you can run, e.g.

::

    $ git fetch ssh://USER@gerrit.gromacs.org/gromacs refs/changes/?/?/? && git checkout FETCH_HEAD

Make your changes, then add them to the index, and use

::

    $ git commit --amend
    $ git push origin HEAD:refs/for/BRANCH

When amending the previous commit message, leave the "Change-Id" intact
so that gerrit can recognize this is an update and not open a new issue.

DO NOT rebase your patch set and update it in one step. If both are done
in one step, the diff between patch set versions has both kinds of
changes. This makes it difficult for the reviewer, because it is not
clear what parts have to be re-reviewed. If you need to update and
rebase your change please do it in two steps (order doesn't matter).
gerrit has a feature that allows you to rebase within gerrit, which
creates the desired independent patch for that rebase (if the rebase is
clean).

How do I get a copy of my commit for which someone else has uploaded a patch?
-----------------------------------------------------------------------------

Gerrit makes this easy. You can download the updated commit in various
ways, and even copy a magic git command to your clipboard to use in your
shell.

You can select the kind of git operation you want to do (cherry-pick is
best if you are currently in the commit that was the parent, checkout is
best if you just want to get the commit and not worry about the current
state of your checked out git branch) and how you want to get it. The
icon on the far right will paste the magic shell command into your
clipboard, for you to paste into a terminal to use.

How do I submit lots of independent commits (e.g. bug fixes)?
-------------------------------------------------------------

Simply pushing a whole commit tree of unrelated fixes creates
dependencies between them that make for trouble when one of them needs
to be changed. Instead, from an up-to-date repo, create and commit the
first change (or git cherry-pick it from an existing other branch).
Upload it to gerrit. Then do

::

    $ git reset --keep HEAD^

This will revert to the old HEAD, and allow you to work on a new commit
that will be independent of the one you've already uploaded. The one
you've uploaded won't appear in the commit history until it's been
reviewed and accepted on gerrit and you've pulled from the main repo,
however the version of it you uploaded still exists in your repo. You
can see it with git show or git checkout using its hash - which you can
get from the gerrit server or by digging in the internals of your repo.

How can I avoid needing to remember all these arcane git commands?
------------------------------------------------------------------

In your ``.gitconfig``, having set the git remote for the gerrit repo to
upload, use something like the following to make life easier:

::

    [alias]
            upload-r2018  = push origin HEAD:refs/for/release-2018
            upload-r2016  = push origin HEAD:refs/for/release-2016
            upload-master = push origin HEAD:refs/for/master
            upload-reset  = reset --keep HEAD^


How can I get my patch in gerrit to have a different parent?
------------------------------------------------------------

Sometimes, some other patch under review is a relevant point from which
to start work. For simple changes without conflicts to the previous
work, you can use the Gerrit web UI to either rebase or cherry-pick
the change you are working on.

If this is not possible, you can still use
the canned gerrit checkouts to (say) checkout out patch 2117 and start work:

::

    git fetch https://gerrit.gromacs.org/gromacs refs/changes/17/2117/2 && git checkout FETCH_HEAD

Other times you might have already uploaded a patch (e.g. patch 1 of
2145), but now see that some concurrent work makes more sense as a
parent commit (e.g. patch 2 of 2117), so check it out as above, and then
use the canned gerrit **cherry-pick**:

::

    git fetch https://gerrit.gromacs.org/gromacs refs/changes/45/2145/1 && git cherry-pick FETCH_HEAD

Resolve any merge commits, check things look OK, and then upload.
Because the ChangeId of 2145 hasn't changed, and nothing about 2117 has
changed, the second patch set of 2145 will reflect the state of 2145 now
having 2117 as a parent.

This can also be useful for constructing a short development branch
where the commits are somehow dependent, but should be separated for
review purposes. This technique is useful when constructing a series of
commits that will contribute to a release.

How can I revert a change back to an old patchset?
--------------------------------------------------

If a change accidentally gets updated or when a patchset is incorrect,
you might want to revert to an older patchset. This can be done by
fetching an old patchset, running git commit --amend to update the time
stamp in the commit and pushing the commit back up to gerrit. Note that
without the amending you will get an error from the remote telling you
that there are no new changes.

How do I handle common errors
-----------------------------

.. rubric:: error: server certificate verification failed. CAfile...

If you try to cherry-pick a change from the server, you'll probably get
the error:

::

    $ git fetch https://gerrit.gromacs.org/p/gromacs refs/changes/09/109/1 && git cherry-pick FETCH_HEAD
    error: server certificate verification failed.
    CAfile: /etc/ssl/certs/ca-certificates.crt
    CRLfile: none while accessing https://gerrit.gromacs.org/p/gromacs/info/refs

    fatal: HTTP request failed

As explained
`here <http://code.google.com/p/chromium-os/issues/detail?id=13402>`__,
the problem is with git not trusting the certificate and as a workaround
one can set globally

::

    $ git config --global --add http.sslVerify false

or prepend GIT_SSL_NO_VERIFY=1 to the command

::

    $ GIT_SSL_NO_VERIFY=1  git fetch https://gerrit.gromacs.org/p/gromacs refs/changes/09/109/1 \
     && git cherry-pick FETCH_HEAD

.. rubric:: Various error messages and their meanings

http://review.coreboot.org/Documentation/error-messages.html

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
commit into gerrit for review. If
you stashed away your changes and
you want the next change to be
reviewed independently, do

::

    git reset --hard HEAD^
    git stash pop

(only do this if you pushed the
previous change to gerrit,
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

.. rubric:: Interacting with Gerrit
   :name: interacting-with-gerrit
   :class: editable

This section is intended for
using git to interact with
gerrit; interacting with the web
UI may be better dealt with on a
separate page.

.. rubric:: Q: How do I move a change from a branch to another?

A: Moving one or a few changes is
most easily done using ``git
cherry-pick``. To move a single
change, first do

::

    git checkout <target-branch>

Then, open the change/patch set
in Gerrit that you want to move,
select "cherry-pick" in the
Download section for that patch
set, and copy/paste the given
command:

::

    git fetch ... refs/changes/... && git cherry-pick FETCH_HEAD

Resolve any conflicts and do

::

    git commit [-a]

You can also cherry-pick multiple
changes this way to move a small
topic branch. Before pushing the
change to Gerrit, remove the
lines about conflicts from the
commit message, as they don't
serve any useful purpose in the
history. You can type that
information into the change as a
Gerrit comment if it helps the
review process. Note that Gerrit
creates a new change for the
target branch, even if Change-Ids
are same in the commits. You need
to manually abandon the change in
the wrong branch.

