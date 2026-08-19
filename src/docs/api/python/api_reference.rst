..
  Copyright Contributors to the rawtoaces Project.
  SPDX-License-Identifier: CC-BY-4.0

Python API Reference
====================

Classes
-------

.. autoclass:: rawtoaces.ImageConverter

   .. automethod:: rawtoaces.ImageConverter.process_image
   .. automethod:: rawtoaces.ImageConverter.configure
   .. automethod:: rawtoaces.ImageConverter.get_WB_multipliers
   .. automethod:: rawtoaces.ImageConverter.get_transform_matrix
   .. automethod:: rawtoaces.ImageConverter.get_IDT_matrix
   .. automethod:: rawtoaces.ImageConverter.get_CAT_matrix
   .. automethod:: rawtoaces.ImageConverter.get_supported_formats
   .. automethod:: rawtoaces.ImageConverter.get_supported_illuminants
   .. automethod:: rawtoaces.ImageConverter.get_supported_cameras
  
.. autoclass:: rawtoaces.ImageConverter.Settings
   :members:
   :member-order: groupwise

.. autoclass:: rawtoaces.Metadata
   :members:
   :member-order: groupwise
   
.. autoclass:: rawtoaces.MetadataSolver
   :members:
   :inherited-members:
   :exclude-members: calculate_CAT_matrix, calculate_IDT_matrix, calculate_transform
   :member-order: groupwise
   
   .. automethod:: rawtoaces.MetadataSolver.calculate_CAT_matrix
   .. automethod:: rawtoaces.MetadataSolver.calculate_IDT_matrix
   .. automethod:: rawtoaces.MetadataSolver.calculate_transform
    
.. autoclass:: rawtoaces.SpectralData
   :members:
   :member-order: groupwise
   :exclude-members: load
   
   .. automethod:: rawtoaces.SpectralData.load
   
.. autoclass:: rawtoaces.SpectralSolver
   :members:
   :inherited-members:
   :exclude-members: collect_data_files, find_camera, find_illuminant,
      calculate_WB, load_spectral_data, get_WB_multipliers,
      calculate_IDT_matrix, get_IDT_matrix, calculate_transform
   :member-order: groupwise
      
   .. automethod:: rawtoaces.SpectralSolver.collect_data_files
   .. automethod:: rawtoaces.SpectralSolver.find_camera
   .. automethod:: rawtoaces.SpectralSolver.find_illuminant
   .. automethod:: rawtoaces.SpectralSolver.calculate_WB
   .. automethod:: rawtoaces.SpectralSolver.load_spectral_data
   .. automethod:: rawtoaces.SpectralSolver.get_WB_multipliers
   .. automethod:: rawtoaces.SpectralSolver.calculate_IDT_matrix
   .. automethod:: rawtoaces.SpectralSolver.get_IDT_matrix
   .. automethod:: rawtoaces.SpectralSolver.calculate_transform

Utility Functions
-----------------

.. autofunction:: rawtoaces.collect_image_files
