Usage
=====

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

Alternatively, you can download or clone the repository and use `pip` to handle dependencies:

``unzip bacteria.zip``

``pip install -e bacteria``

or

``git clone https://github.com/tuliofalmeida/bacteria``

``pip install -e bacteria``

By calling ``pip list`` you should see ``bacteria`` now as an installed package:
``bacteria (1.x.x, /path/to/bacteria)``