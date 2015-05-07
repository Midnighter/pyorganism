=============================================
Organisational Principles in Living Organisms
=============================================


This Python package provides data structures to model and
simulate the organisation at (unicellular) system scale.

* Genetic regulation (``pyorganism.regulation``)
    * Analog and digital control (cf. [1_])
    * Transcriptional regulatory network (TRN): directed transcription factor-gene interaction network
    * Gene regulatory network (GRN): directed gene-gene interaction network
      (projection of the TRN)
    * Gene proximity network (GPN): undirected gene-gene proximity network
    * Couplons (intersection between sigma factor's and nucleoid associated
      proteins' targets)
    * Functional clusters (GO)
    * Metabolic control (cf. [2_])
    * Continuous control (cf. [3_])
* Metabolism (``pyorganism.metabolism``)
    * Metabolic systems
    * Network representations
    * Flux balance analysis (FBA)


Requirements
------------

* future_ for Python2/3 compatiblity
* Cython_
* networkx_
* numpy_
* pandas_
* reading/writing SBML files requires libsbml_ or the `PyPi package`_
* FBA requires either:
    * Gurobi_ (currently the only working interface),
    * GLPK_,
    * MOSEK_,
    * or any other linear programming solver that you are willing to provide the
      Python interface for (like gurobipy or cvxopt_).
* KEGG interface requires SOAPpy_
* (pytables_ for HDF5 storage can be replaced by pickling of numpy objects)

.. _future: http://python-future.org/
.. _Cython: http://cython.org/
.. _networkx: http://networkx.github.com/
.. _numpy: http://www.numpy.org/
.. _pandas: http://pandas.pydata.org/
.. _libsbml: http://sbml.org/Software/libSBML
.. _Gurobi: http://www.gurobi.com/
.. _GLPK: http://www.gnu.org/software/glpk/
.. _MOSEK: http://www.mosek.com/
.. _cvxopt: http://abel.ee.ucla.edu/cvxopt/
.. _SOAPpy: http://pywebsvcs.sourceforge.net/
.. _pytables: http://www.pytables.org/
.. _`PyPi package`: https://pypi.python.org/pypi/python-libsbml-experimental


Authors
-------


Former and current members of the `Computational Systems Biology workgroup`_ headed by Prof.  Marc-Thorsten Hütt at Jacobs University Bremen. In alphabetical order:

* Beber, Moritz Emanuel
* Grigore, Alexandra Mirela
* Kölling, Nils
* Sonnenschein, Nikolaus

.. _`Computational Systems Biology workgroup`: http://sysbio.jacobs-university.de/website/


Relevant References
-------------------


.. [1] Marr, C., Geertz, M., Hutt, M.-T. & Muskhelishvili, G. Dissecting the logical types of network control in gene expression profiles. *BMC Systems Biology* 2, 18 (2008).
.. [2] Sonnenschein, N., Geertz, M., Muskhelishvili, G. & Hütt, M.-T. Analog regulation of metabolic demand. *BMC Systems Biology* 5, 40 (2011).
.. [3] Beber, M. E., Sobetzko, P., Muskhelishvili, G. & Hütt, M.-T. Interplay of digital and analog control in time-resolved gene expression profiles. *BMC Systems Biology* submitted, (2015).
