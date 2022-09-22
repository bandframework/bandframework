BAND Framework Release Process
===============

A release can be undertaken only by a project administrator.
A project administrator should have an administrator role on the `bandframework
GitHub <https://github.com/bandframework>`_.

Before release
--------------

- A release branch should be taken off ``develop`` (or ``develop`` pulls
  controlled).

- Release notes for this version are added to the `CHANGELOG.rst </CHANGELOG.rst>`_ file.

- Ensure that links and references have been updated (e.g., no occurrences of ``privateband``).

- Version number is updated wherever it appears and ``+dev`` suffix is removed
  (in `README.md </README.md>`_).

- Check `README.md </README.md>`_ *Citing bandframework* for correctness (e.g., ensure that author list matches `AUTHORS </AUTHORS>`_).

- Tests are run with source to be released (this may iterate):

  - Online CI (GitHub Actions) tests must pass.

  - Documentation must build and display correctly wherever hosted.

- Pull request from either the ``develop`` or release branch to ``main`` requesting
  one or more reviewers (including at least one other administrator).

- Reviewer will check that all tests have passed and will then approve merge.

- Send an email to bandframework@cels.anl.gov asking to confirm whether the list of email contacts in `CODE_OF_CONDUCT.md <CODE_OF_CONDUCT.md>`_ is up-to-date.

During release
--------------

An administrator will take the following steps.

- Merge the pull request into ``main``.

- Once CI tests have passed on ``main``:

  - A GitHub release will be taken from the ``main``

- If the merge was made from a release branch (instead of ``develop``), merge this
  branch into ``develop``.

- Create a new commit on ``develop`` that appends ``+dev`` to the version number
  (wherever it appears).

After release
-------------

- Ensure all relevant GitHub issues are closed.
