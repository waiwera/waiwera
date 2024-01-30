.. index:: Waiwera; input, input

.. _waiwera_input:

*************
Waiwera input
*************

Waiwera takes its main input from a file in `JSON <https://www.json.org/>`_ format. JSON is a lightweight data-interchange format, designed to be easily human-readable and writeable, as well as being easy for computers to generate and parse.

.. index:: JSON
.. index:: JSON; editing

Editing JSON files
==================

A file in JSON format is just a text file, so it can be read and edited using a text editor. Many text editors (especially those designed for programmers) have special modes for helping the user prepare and read JSON files. These will highlight special JSON characters and keywords in different colours, to make them easier to read.

JSON syntax makes a lot of use of brackets and braces (see below), so another useful text-editor feature is **bracket matching**. This allows the user to identify matching pairs of brackets within the file -- e.g. to find the closing bracket matching an opening bracket (and check that the closing bracket is not missing).

There is also dedicated software designed specifically for editing JSON files, which can take care of much of the syntax for you. Examples include `JSON editor online <https://jsoneditoronline.org/>`_.

.. index:: JSON; structure

JSON file structure
===================

JSON files can be structured in a **hierarchical** way, so that larger data structures in the file can contain smaller ones, which can in turn contain their own smaller data structures, and so on. The Waiwera input JSON files make use of this hierarchical structure to group related data together.

JSON allows some **flexibility** in the way data are specified. For example, data may sometimes be specified in different ways, depending on the problem. Also, some data may not be required to be present, and can be given default values if they are absent.

.. index:: JSON; data structures
.. _json_data_structures:

JSON data structures
====================

The data structures representable in JSON are defined in detail on the `JSON website <https://www.json.org/>`_. There are not many of them, but they can be combined together to create complex composite data structures.

The two main JSON data structures are the **array** and the **object**.

.. index:: JSON; arrays

JSON arrays
-----------

An **array** is an ordered collection of values, beginning and ending with square brackets (``[ ]``), and with the values delimited by commas. (If you are familiar with the Python scripting language, this structure is known in Python as a "list".)

.. index:: JSON; strings, JSON; numbers, JSON; boolean, JSON; null

The values in the array may be character **strings** (enclosed by double quotes), **numbers** (integer or floating point), **boolean** values (``true`` or ``false``) or **null** values (``null``). In addition (and importantly), the values may themselves be arrays or objects. The values in an array do not all have to be of the same type.

.. index:: JSON; objects
.. _json_objects:

JSON objects
------------

An **object** is an unordered collection of named values, i.e. pairs of names (or "keys") and their corresponding values. The object is enclosed by `braces` (curly brackets, ``{ }``). Each name (a character string) is followed by a colon (``:``) and then the value. The pairs of names and values are separated by commas. (This data structure is known in Python as a "dictionary".)

The values in an object can be of the same types as in an array (strings, numbers, etc.). Once again the values may themselves be arrays or objects.

Here is a simple example of a JSON object, in which the values are all numbers:

.. code-block:: json

  {"length": 20.5, "width": 14.8, "depth": 17.2}

Here is a more complex example of an object with a mix of different value types (an object, an array of numbers, a string and a boolean):

.. code-block:: json

  {"dimensions": {"length": 20.5, "width": 14.8, "depth": 17.2},
   "position": [161.2, -12.5, 405.1],
   "colour": "blue",
   "checked": true
  }

Given a JSON object which contains a sub-object, it is sometimes convenient to refer to values inside the sub-object directly. This is done by using a JSON **path**, consisting of the value names joined together with dots between them. So in the above example we might refer to the value "dimensions.width", which in this case would have the value 14.8.

.. index:: JSON; layout

JSON file layout
================

The exact positions of the data structures in a JSON file are not important. Arrays and object are delimited by brackets or braces as described above, and values are separated by commas, so the data can be located anywhere on each line. Usually they are arranged in a way that is easy for the user to read.

For example, the following two JSON files are completely equivalent:

.. code-block:: json

   {"name": "Ted Smith",
    "age": 76}

.. code-block:: json

   {"name": "Ted Smith", "age": 76}

.. index:: JSON; scripting, scripting; JSON
.. _json_script:

JSON files in scripts
=====================

Many scripting and programming languages are able to parse JSON files (either directly or via add-on libraries), so that the JSON file contents can then be manipulated using a script or program. Some are also able to write JSON data structures out to JSON files.

For example, the Python scripting language has a built-in ``json`` module for this. Python has data structures equivalent to those found in JSON files, although some have different names: JSON arrays are called "lists" in Python, and JSON objects are called "dictionaries".

In the example script below, the JSON object in the previous example (see :ref:`json_objects`) is read into a Python dictionary. The "dimensions.depth" value is then changed (note the syntax for accessing an object within an object), and a new value called "age" is added to the object. Then the edited object is written out to a new JSON file.

.. code-block:: python

   import json

   data = json.load(file('input.json', 'r'))

   data["dimensions"]["depth"] = 17.3
   data["age"] = 3.5

   json.dump(data, file('edited.json', 'w'))

.. index:: JSON; validation

JSON validation
===============

If you are editing a JSON file manually using a text editor, particularly if it contains complex nested data structures, it can be easy to make small syntactical mistakes, e.g. forgetting a closing bracket. JSON is not forgiving of such mistakes, so it is often wise to **validate** a JSON file before using it.

There are software tools available to check if a the contents of a file are valid JSON, and point out any errors. There are JSON validation tools available inside some text editors, or as stand-alone tools, or as online tools such as `this <https://jsonlint.com/>`_ one.

.. warning::
It is not a good idea to have duplicate keys in a JSON object (although the JSON specification does not explicitly disallow it). Different JSON software implementations handle duplicates differently, so you cannot rely on any particular duplicate (e.g. the last one in the file) being the one that will always be used. Some JSON validators (e.g. the one linked to above) will treat duplicate keys as invalid and produce warnings.

.. index:: JSON; schema

Validation against a JSON schema
================================

The exact expected structure of a JSON file for a particular application (e.g. Waiwera input) can be defined using a JSON `schema <http://json-schema.org/>`_, a file (actually a JSON file itself) containing the data structure specification. As well as providing a technical specification for the data structures, the schema can be used to check if a particular JSON file has the required structure. Software `tools <http://json-schema.org/implementations.html>`_ (e.g. libraries, text editor plugins or online validators like `this one <https://jsonschemalint.com>`_) are available for automatically validating a JSON file against the schema.

Schema validation is useful for avoiding input errors, e.g. missing data, or data with mis-spelled keys that would otherwise be ignored or given default values. (Basic JSON syntax validation, as described above, should be done first as well.)

The `utils/` directory of the Waiwera source code contains a schema file (`input_schema.json`) for Waiwera JSON input.

.. index:: Waiwera; JSON file structure

Waiwera JSON file structure
===========================

All input data for a Waiwera simulation are contained within one object in the input JSON file.

The various aspects of the simulation are specified by different named values within that object, as listed below. These values are of various types, and many are objects themselves, with their own internal structure.

.. note::

   **JSON object**: Waiwera JSON input

   +-----------------+-----------------+--------------------------+
   |**name**         |**type**         |**value**                 |
   +-----------------+-----------------+--------------------------+
   |"title"          |string           |simulation title          |
   +-----------------+-----------------+--------------------------+
   |"mesh"           |string | object  |simulation mesh           |
   |                 |                 |                          |
   +-----------------+-----------------+--------------------------+
   |"rock"           |object           |rock properties           |
   +-----------------+-----------------+--------------------------+
   |"boundaries"     |array            |boundary conditions       |
   +-----------------+-----------------+--------------------------+
   |"source"         |array            |source terms              |
   +-----------------+-----------------+--------------------------+
   |"initial"        |object           |initial conditions        |
   +-----------------+-----------------+--------------------------+
   |"gravity"        |number | array | |gravity                   |
   |                 |``null``         |                          |
   +-----------------+-----------------+--------------------------+
   |"thermodynamics" |string           |thermodynamic formulation |
   |                 |                 |                          |
   +-----------------+-----------------+--------------------------+
   |"eos"            |object           |equation of state         |
   +-----------------+-----------------+--------------------------+
   |"tracer"         |object           |tracers                   |
   +-----------------+-----------------+--------------------------+
   |"time"           |object           |time stepping             |
   +-----------------+-----------------+--------------------------+
   |"logfile"        |boolean | object |output log file           |
   +-----------------+-----------------+--------------------------+
   |"output"         |boolean | object |output results file       |
   +-----------------+-----------------+--------------------------+

.. index:: JSON; schema
