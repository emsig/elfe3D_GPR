reStructuredText Quick Reference
================================

This guide summarizes the most important reStructuredText syntax for the `docs/source` folder.

Headings
--------

Use underlines to define section titles. The character used indicates the heading level.

Basic example:

Heading Level 1
===============

Heading Level 2
---------------

Heading Level 3
^^^^^^^^^^^^^^^

Lists
-----

Unordered list:

- Item one
- Item two
- Item three

Ordered list:

1. First item
2. Second item
3. Third item

Nested list:

- Top level
  - Subitem
  - Subitem
- Another top level

Code blocks
-----------

Use `.. code-block::` for syntax-highlighted code.

.. code-block:: python

   from pathlib import Path
   print(Path('.') / 'docs')

For plain literal blocks, indent text by 3 spaces:

   This is a literal block.
   It keeps whitespace exactly.

Links
-----

External link:

`Read the Docs <https://readthedocs.org>`_

Internal page link:

:doc:`overview`

Cross-reference to a label:

.. _my-label:

Section title
-------------

See :ref:`my-label` for details.

Images
------

.. image:: ../_static/model-diagram.png
   :alt: Model diagram
   :align: center

Tables
------

Simple table using a list table directive:

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 80

   * - Field
     - Description
   * - ``air_eps_r``
     - Relative permittivity for the air layer
   * - ``layer_sigma``
     - Electrical conductivity for each layer

Citations and bibliography
--------------------------

Define citation references with ``.. [KEY]`` in a page such as `references.rst`:

.. [DOE2025] Doe, J. (2025). Example title. Publisher.

Use them in text with a standard footnote style if you have an extension that supports it.

Labels and references
---------------------

Define explicit labels for sections:

.. _installation-guide:

Installation
------------

Refer to the section with :ref:`installation-guide`.

To link to another page, use ``:doc:`` or normal inline links:

- :doc:`quickstart`
- `Quick Start <quickstart.html>`_

Directives and roles
--------------------

Common directives:

- ``.. toctree::`` for table of contents trees
- ``.. code-block::`` for code examples
- ``.. image::`` for images
- ``.. note::`` for notes
- ``.. warning::`` for warnings

Example note:

.. note::

   This is a helpful note for readers.

Inline markup
-------------

- ``*italic*`` for emphasis
- ``**bold**`` for strong emphasis
- ````inline code```` for code literals
- ``:math:`` if enabled in your Sphinx config

Sphinx-specific page structure
------------------------------

The landing page typically includes a ``.. toctree::`` block listing included pages.

Example:

.. toctree::
   :maxdepth: 2

   overview
   installation
   quickstart
   references

Tips
----

- Use consistent heading levels.
- Keep internal links with ``:doc:`` and ``:ref:`` for robustness.
- Use ``_`` after inline links to make them active.
- Prefer explicit labels on important sections when you need cross-page references.
