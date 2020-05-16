To run the augmented reality viewer:
   -Simply open run the augmentedRealityViewer.m script in Matlab.

The files:
   - augmentedRealityViewer.m uses relative paths to access the
     folders COLMAP_Files and Scenes
   - if these files are moved, the paths need to be adjusted
     accordingly

The functions:
   - There is one function, CameraProjection3Dto2D, that our script
     calls. This function is located at the bottom of the script, it
     is NOT contained in a separate location