Cheat Sheet
============

Wopfile or Workflow definition file example
---------------------------------------------------

This is a full :download:`Wopfile example </data/Wopfile.yml>`. Here there are the first lines:

.. code-block:: bash

    rule SampleInformation:
        tool: wopmetabarcoding.wrapper.SampleInformation
        input:
            file:
                csv: data/1_sample_info.csv
    ...



