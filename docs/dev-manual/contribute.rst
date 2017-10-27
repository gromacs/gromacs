Contribute to |Gromacs|
=======================

Please read this part first when you plan to contribute new functionality
to |Gromacs|. We like to encourage contributions, because new things can
lead to exciting science, but the subsequent maintenance is expensive,
and the success of the project is compromised if the overall quality
of code and methods is uneven.

If you have any questions regarding contribution or adding code to
the project, don't hesitate to ask on the `developer mailing list`_.

We also strongly encourage additions to documentation or the manual
for both the back-end and front-end of |Gromacs|. Those contributions
can be done using the same mechanism as the source code contributions,
and will be reviewed in similar ways.

Checklist
---------

When you plan to send in your code, please make sure that all the
points on this list have been fulfilled before you send in your code for review:

* Usefulness: Your code should have wide applicability within the scientific
  community. You are welcome to have smaller projects tracking our code,
  but we are not prepared to include and maintain code that will only have
  limited application. Evidence that people are already using your code or
  method is one good way to show that your code is useful.
  Scientific publications is another, but those publications should
  ideally come from several different research groups to show
  widespread adoption of the method.

* Advance discussion: Please communicate with the other developers, e.g.
  on the `gmx-developers <gmx-developers@gromacs.org>`_ mailing list,
  or `Redmine <https://redmine.gromacs.org>`_ to let them know
  of the general nature of your plans. This will prevent duplicate or
  wasted effort. It is also a good idea to search those resources as
  well as the literature and WWW for other projects that may be relevant.

* Verifiable: If you propose a new method that passes the first check,
  please make sure that we can easily verify that it will be correct
  from a physics point of view. That must include documentation (both 
  in the source code and as later additions to the reference manual) that
  a capable graduate student can read and understand well enough to use
  your method appropriately. The source code documentation will also
  help in maintenance and later development.

  This will be facilitated by the inclusions of unit tests for your code,
  as described in the section on how to write
  :ref:`new tests <gmx-make-new-tests>`.

* Coding style: Please make sure that your code follows all the
  :ref:`coding style guidelines <code-formatting>`. This will make
  the code review go more smoothly on both sides. There are a number of
  tools already included with |Gromacs| to facilitate this, please have
  a look at :ref:`the respective part of the documentation <gmx-uncrustify>`.
  This probably means high-level code documentation,
  as well as doxygen for all new functions and classes.

  In addition to coding style, please also follow the instructions given
  concerning the :ref:`commit style <code-commitstyle>`. This will also
  facilitate the code review process.


.. TODO add more points here to make things clear

If you have questions regarding those points, please feel free to contact
one of the |Gromacs| developers to clarify issues. Once we have decided that
your code would be of use for |Gromacs|, we will then assist you in the
process of including it in the main source tree. The best way to contact
the team is through the `developer mailing list`_.

.. _developer mailing list: gmx-developers@gromacs.org
