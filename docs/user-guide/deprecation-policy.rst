.. _deprecation-policy:

Policy for deprecating |Gromacs| functionality
==============================================

Occasionally functionality ceases being useful, is unable to be fixed
or maintained, or its user interface needs to be improved. The
development team does this sparingly. Broken functionality might be
removed without notice if nobody willing to fix it can be found.
Working functionality will be changed only after announcing in the
previous major release the intent to remove and/or change the form of
such functionality. Thus there is typically a year for users and
external tool providers to prepare for such changes, and contact the
|Gromacs| developers to see how they might be affected and how best to
adapt.

There is a current list of anticipated changes and deprecated functionality
in the "Major release" :ref:`notes <release-notes>`.

When environment variables are deprecated, it is up to the user to make
sure that their scripts are updated accordingly for the new release. In
cases where it is sensible, the development team should do the effort to
keep the old environment variables working for one extra release cycle,
before fully removing them. The user should be informed about this future
deprecation with a warning. If keeping the old environment variable is
not possible or highly problematic, setting the removed environment
variable should be triggering a warning during one release cycle.
