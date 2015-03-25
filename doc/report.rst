MIGEC analysis report
---------------------

You can run the following command to generate

.. code:: bash

    java -jar migec.jar Report [options] [output_path or .]
    
Command parameters set the checkout (``-c``), histogram (``-h``),
assemble (``-a``), cdrblast (``-b``) and final results (``-f``) folders. Missing
results will be excluded from the report.

.. warning:: Works for **batch** analysis only.
    
.. note::
    
    Running this command requires installing dependencies for R markdown compilation,
    see http://rmarkdown.rstudio.com/
    
Alternatively, you can install 
`Rstudio <http://www.rstudio.com/>`__ 
and use 
`this Rmd template <https://github.com/mikessh/migec/blob/master/util/migec_summary.Rmd>`__
to manually knit the report HTML.