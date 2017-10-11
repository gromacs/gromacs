.. _code-commitstyle:

Guidelines for formatting of git commits
========================================

While there is no true correct way on how to submit new commits for 
code review for |Gromacs|, following these guidelines will help the
review process go smoothly.

General rules for newly submitted code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

New code should follow the other :ref:`style rules<style-guidelines>`
outlined above before submitting. This will make it less likely that your change
will be rejected due to that. If your change modifies some existing
code that does not yet conform to the style, then a preliminary
patch that cleans up the surrounding area is a good idea. We like
to slowly improve the quality while we add or change functionality.

Guidelines for git commit messages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Commit messages should contain a quick explanation in verb form on what has been
changed or what has been the purpose of the change. If available, the final
part of the message before the ChangeId should be a short section like
**Fixes #redmine-id** to link the change to a possibly previously 
posted issue, or **Refs #redmine-id** if the present patch is somehow
related to that work without necessarily fixing the whole issue.

Concerning inline code comments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

New code should be sufficiently commented so that other people will be able to 
understand the purpose of the code, and less about the current operation.
Preferably the variable naming and code structure clarify the mechanics, and
comments should only refer to higher-level things, such as choice of algorithm,
or the desire to be consistent with some other part of the code.

For example, the following comment would be insufficient to explain the 
(made up example) of iteration over a list of interactions::

    /* Code takes each item and iterates over them in a loop
     * to store them.
     */

A much better example would be explaining why the iteration takes place::

    /* We iterate over the items in the list to get 
     * the specific interaction type for all of them
     * and store them in the new data type for future 
     * use in function foo
     */

From the second example, someone debugging might be able to deduce better
if an error observed in *foo* is actually caused by the previous assignment.
