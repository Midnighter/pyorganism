=============================================
Organisational Principles in Living Organisms
=============================================


This Python package consists of modules that represent ways to model and
simulate the organisation of whole organisms. The starting point are unicellular
organisms (or single cells) and the following modules and abstractions for
different levels of organisation.

The package depicts our current state of knowledge on and representation of
the regulation and its regulated functions in (unicellular) organisms.

* Genetic regulation
    * Transcriptional regulation (directed gene-gene interaction network - TRN)
    * Gene proximity network (undirected gene-gene proximity network - GPN)
    * Couplons (intersection between sigma factor's and nucleoid associated
      protein's targets)
    * Functional clusters (GO)
* Metabolism
    * Network representations
    * Flux balance analysis (FBA)
    * Metabolic control


Requirements
------------


* networkx_
* numpy_
* FBA requires either:
    * Gurobi_,
    * GLPK_,
    * MOSEK_,
    * or any other linear programming solver that you are willing to provide the
      Python interface for (like gurobipy or cvxopt_).
* (pytables_ for HDF5 storage can be replaced by pickling of numpy objects)

.. _networkx: http://networkx.github.com/
.. _numpy: http://www.numpy.org/
.. _Gurobi: http://www.gurobi.com/
.. _GLPK: http://www.gnu.org/software/glpk/
.. _MOSEK: http://www.mosek.com/
.. _cvxopt: http://abel.ee.ucla.edu/cvxopt/
.. _pytables: http://www.pytables.org/


Authors
-------


Former and current members of the `Computational Systems Biology workgroup`__ headed by Prof.  Marc-Thorsten Hütt at Jacobs University Bremen. In alphabetical order:

* Beber, Moritz Emanuel
* Grigore, Alexandra Mirela
* Kölling, Nils
* Sonnenschein, Nikolaus

.. __: http://sysbio.jacobs-university.de/website/


Relevant References
-------------------


.. [1] Marr, C., Geertz, M., Hütt, M.-T., Muskhelishvili, G., 2008. 'Dissecting the logical types of network control in gene expression profiles'.  BMC Systems Biology 2, 18.

