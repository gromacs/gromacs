Contribute to |Gromacs|
=======================

If you are planning to contribute new functionality or fixes to |Gromacs|,
we encourage you to get in contact with us, because new things can
lead to exciting science. However, as  the subsequent code maintenance is
time-consuming and requires long-term commitment, please read this part first.
We therefore strongly encourage you to communicate your plans with us and
read the information below already at an early stage of your project

We also strongly encourage additions to documentation or the manual
for both the back-end and front-end of |Gromacs|. Those contributions
can be done using the same mechanism as the source code contributions,
and will be reviewed in similar ways.

Checklist
---------

Before you send us your code for review and inclusion into |Gromacs|,
please make sure that you have checked all the points on this list:

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
  :doc:`coding style <style>` and :ref:`code formatting <code-formatting>`
  guidelines. This will make the code review go more smoothly on both sides. There are a number of
  tools already included with |Gromacs| to facilitate this, please have
  a look at :ref:`the respective part of the documentation <gmx-uncrustify>`.

* Code documentation: To ensure proper code documentation, please follow the 
  instructions provided for the use of :doc:`doxygen <doxygen>`. In addition to this,
  the new functionality should also be documented in the manual and possible the user guide .

* In addition to coding style, please also follow the instructions given
  concerning the :ref:`commit style <code-commitstyle>`. This will also
  facilitate the code review process.


.. TODO add more points here to make things clear

If you have questions regarding these points, or would like feedback on your ideas for contributing,
please feel free to contact us through the `developer mailing list`_.
If your code is of interest to the wider |Gromacs| community, we will be happy to assist you
in the process of including it in the main source tree.

.. _developer mailing list: gmx-developers@gromacs.org
