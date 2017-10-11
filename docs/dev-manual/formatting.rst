.. _code-formatting:

Guidelines for code formatting
==============================

The following list provides the general formatting/indentation rules for
|Gromacs| code (C/C++):

* Basic indentation is four spaces.
* Keep lines at a reasonable length.  There is no hard limit, but use 80
  characters as a guideline.  If you end up indenting very deeply,
  consider splitting the code into functions.
* Do not use tabs, only spaces.  Most editors can be configured to generate
  spaces even when pressing tab.  Tabs (in particular when mixed with spaces)
  easily break indentation in contexts where settings are not exactly equal
  (e.g., in ``git diff`` output).
* No trailing whitespace.
* Use braces always for delimiting blocks, even when there is only a single
  statement in an ``if`` block or similar.
* Put braces on their own lines.  The only exception is short one-line inline
  functions in C++ classes, which can be put on a single line.
* Use spaces liberally.
* ``extern "C"`` and ``namespace`` blocks are not indented, but all others
  (including ``class`` and ``switch`` bodies) are.

Additionally:

* All source files and other non-trivial scripts should contain a copyright
  header with a predetermined format and license information (check existing
  files).  Copyright holder should be "the |Gromacs| development team" for the
  years where the code has been in the |Gromacs| source repository, but earlier
  years can hold other copyrights.
* Whenever you update a file, you should check that the current year is listed
  as a copyright year.

Most of the above guidelines are enforced using uncrustify, an automatic source
code formatting tool.  The copyright guidelines are enforced by a separate
Python script.  See :doc:`uncrustify` for details.  Note that due to the
nature of uncrustify (it only does all-or-nothing formatting), it enforces
several additional formatting rules in addition to those above.

Enforcing a consistent formatting has a few advantages:

* No one needs to manually review code for most of these formatting issues,
  and people can focus on content.
* A separate automatic script (see below) can be applied to re-establish the
  formatting after refactoring like renaming symbols or changing some
  parameters, without needing to manually do it all.

A number of user provided set-ups are available for the correct settings of your
favourite text editor. They are provided for convenience only, and may not
exactly conform to the expectations of uncrustify.

Emacs formatting set-up
-----------------------
Insert the following into your .emacs configuration file::

    (defun gromacs-c-mode-common-hook ()
    ;; GROMACS customizations for c-mode

    (c-set-offset 'substatement-open 0)
    (c-set-offset 'innamespace 0)
    ;; other customizations can go here

    (setq c++-tab-always-indent t)
    (setq c-basic-offset 4)                  ;; Default is 2
    (setq c-indent-level 4)                  ;; Default is 2
    (setq c-file-style "stroustrup")
    (setq tab-stop-list '(4 8 12 16 20 24 28 32 36 40 44 48 52 56 60))
    (setq tab-width 4)
    (setq indent-tabs-mode nil)  ; use tabs if t
    )
    (add-hook 'c-mode-common-hook 'gromacs-c-mode-common-hook)

    (defun gromacs-c++-mode-common-hook ()
    ;; GROMACS customizations for c++-moe

    (c++-set-offset 'substatement-open 0)
    (c++-set-offset 'innamespace 0)
    ;; other customizations can go here

    (setq c++-tab-always-indent t)
    (setq c++-basic-offset 4)                  ;; Default is 2
    (setq c++-indent-level 4)                  ;; Default is 2
    (setq c++-file-style "stroustrup")
    
    (setq tab-stop-list '(4 8 12 16 20 24 28 32 36 40 44 48 52 56 60))
    (setq tab-width 4)
    (setq indent-tabs-mode nil)  ; use tabs if t
    )
    
    (add-hook 'c++-mode-common-hook 'gromacs-c++-mode-common-hook)

This configuration is based on content from `stackoverflow`_.

.. _stackoverflow: http://stackoverflow.com/questions/663588/emacs-c-mode-incorrect-indentation

Eclipse/cdt formatting set-up
-----------------------------

For correct formatting, please use `this profile`_.

.. _this profile: https://gist.github.com/rolandschulz/74f4fae8985d65f33ff6
