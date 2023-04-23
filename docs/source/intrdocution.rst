Intrdocution
==============

.. _installation:

Installation
------------

The latest stable release is available on PyPi, and you can install it by saying

.. code-block:: console

   (.venv) $ pip install bacteria

To build Bacteria from source, say ``python setup.py build``.
Then, to install Bacteria, say ``python setup.py install``.
If all went well, you should be able to execute the demo scripts (read docs)
(OS X users should follow the installation guide given below).

Alternatively, you can download or clone the repository and use ``pip`` to handle dependencies:

.. code-block:: console

   (.venv) $ unzip bacteria.zip
   (.venv) $ pip install -e bacteria

or

.. code-block:: console

   (.venv) $ git clone https://github.com/tuliofalmeida/bacteria
   (.venv) $ pip install -e bacteria

By calling ``pip list`` you should see ``bacteria`` now as an installed package:
``bacteria (1.x.x, /path/to/bacteria)``

.. _IO:

IO
---

The library's input is the SuperSegger data, the clist, and the organization used in 
pandas DataFrame. Its possible to load one clist and extract the 2D and 3D data and 
also concatenate clists from different fields of view into a single dataframe.

.. code-block:: console

   (.venv) $ pip install bacteria