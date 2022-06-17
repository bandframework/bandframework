Release Process
===============

A release can be undertaken only by a project administrator.
A project administrator should have an administrator role on the bandframework
GitHub pages.

Before release
--------------

- A release branch should be taken off ``development`` (or ``development`` pulls
  controlled).

- Release notes for this version are added to the ``CHANGELOG.rst`` file.

- Version number is updated wherever it appears and ``+dev`` suffix is removed
  (in ``README.rst``).

- Check ``README.rst`` *Citing bandframework* for correctness.

- Tests are run with source to be released (this may iterate):

  - On-line CI (GitHub Actions) tests must pass.

  - Documentation must build and display correctly wherever hosted.

- Pull request from either the development or release branch to main requesting
  one or more reviewers (including at least one other administrator).

- Reviewer will check that all tests have passed and will then approve merge.

During release
--------------

An administrator will take the following steps.

- Merge the pull request into main.

- Once CI tests have passed on main:

  - A GitHub release will be taken from the main

- If the merge was made from a release branch (instead of development), merge this
  branch into development.

- Create a new commit on development that appends ``+dev`` to the version number
  (wherever it appears).

After release
-------------

- Ensure all relevant GitHub issues are closed.
