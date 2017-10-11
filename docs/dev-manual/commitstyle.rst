.. _code-commitstyle:

Guidelines for formatting of git commits
========================================

While there is no true correct way on how to submit new commits for 
code review for |Gromacs|, some rules should be observed.

General rules for newly submitted code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The new code should follow the other :ref:`style rules<style-guidelines>`
outlined above before submitting. This will make it less likely that your change
will be rejected due to that. 

Guidelines for git commit messages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Commit messages should contain a quick explanation in verb form on what has been
changed or what has been the purpose of the change. If available, the final
part of the message before the change id should be a short section like
**Fixes #redmine-id** to link the change to a possibly previously 
posted issue.

Concerning inline code comments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

New code should be sufficiently commented so that other people will be able to 
understand the purpose of the code, and less about the current operation.

For example, the following comment would be insufficient to explain the 
(made up example) of iteration over a list of interactions::

    /* Code takes each item and iterates over them in a loop
     * to store them.
     */

A much better example would be explaining why the iteration takes place::

    /*We iterate over the items in the list to get 
     *the specific interaction type for all of them
     *and store them in the new data type for future 
     *use in function foo
     */

From the second example, someone debugging might be able to deduce better
if an error observed in *foo* is actually caused by the previous assignment.
