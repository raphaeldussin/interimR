Installation
------------

Before we start
^^^^^^^^^^^^^^^

**interimR** requires

+ a python3 environment::

    conda create -n forcing
    source activate forcing
    conda install numpy netcdf4 configparser

+ the ECMWF_ python API. You need to register and put obtained key in $HOME/.ecmwfapirc (see howto_access_dataset_).
  Then install API with pip (check for updated version)::

    pip install https://software.ecmwf.int/wiki/download/attachments/56664858/ecmwf-api-client-python.tgz

+ cdo_ (available with most package managers)

Install interimR
^^^^^^^^^^^^^^^^

Clone from github and run the install script::

    git clone https://github.com/raphaeldussin/interimR.git
    python setup.py install


.. _ECMWF: https://www.ecmwf.int
.. _howto_access_dataset: https://confluence.ecmwf.int/display/WEBAPI/Access+ECMWF+Public+Datasets
.. _cdo: https://code.mpimet.mpg.de/projects/cdo/
