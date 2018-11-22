Installation
------------

Before we start
^^^^^^^^^^^^^^^

**interimR** requires

+ a python3 environment::

    conda create -n forcing
    source activate forcing
    conda install numpy netcdf4 configparser

+ the ECMWF_ python API (needs registration, see howto_access_dataset_) and install with (check for updates)::

    pip install https://software.ecmwf.int/wiki/download/attachments/56664858/ecmwf-api-client-python.tgz


.. _ECMWF: https://www.ecmwf.int
.. _howto_access_dataset: https://confluence.ecmwf.int/display/WEBAPI/Access+ECMWF+Public+Datasets
