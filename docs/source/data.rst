Data Structure
===============

.. _data:

This section presents the dataframe columns of 2D and 3D data.

.. note::
   Most of the columns are direct from `SuperSegger`_ and some are 
   added during concatenation (bacteria.concatenate_clist()). 

2D data
--------

 * Cell ID - Unique cell identification
 * Cell birth time - The frame number in which the cell is born (min)
 * Cell age - The total time in which the cell "lived" or was identified
 * stat0 - Indicator of cell birth and division.
 * Long axis (L) birth - Length of the long axis of the cell at birth (px)
 * Long axis (L) death - Length of the long axis of the cell at division (px)
 * Short axis birth - Length of the short axis of the cell at birth (px)
 * Short axis death - Length of the short axis of the cell at division (px)
 * Area birth - Area of the cell at birth (px)
 * Area death - Area of the cell at division (px)
 * Fluor1 sum - Fluorescence intensity integrated over the cell area
 * Fluor1 mean - Average (by area) fluorescence intensity for the cell
 * Mother ID - "Cell ID" of this cell's mother cell
 * Daughter1 ID - "Cell ID" of one daughter cell
 * Daughter2 ID - "Cell ID" of second daughter cell
 * L death / L birth - Ratio of long axis lengths in "Cell Division Time" and "Cell Birth Time"
 * Fluor1 sum death - Fluorescence intensity integrated over the cell area at cell division
 * Fluor1 mean death - Mean (by area) fluorescence intensity for the cell at division
 * Long axis/Short axis birth - Ratio of long to short axis lengths at birth
 * Long axis/Short axis death - Ratio of long to short axis lengths at division
 * Growth Rate - Average cell growth rate
 * Volume birth - Cell volume at birth (:math:`um^3`)
 * Volume division - Cell volume at division (:math:`um^3`)
 * Vd-Vb - Volume at division minus Volume at birth (:math:`um^3`)
 * fov - The field of view that this cell belongs to
 * Time Division - "Cell birth time" + "Cell age" (min)

.. note::
   stat0 = 0 means no birth was observed, stat0 = 1 means no death/division was observed, 
   stat0 = 2 means both birth and death were observed.

3D data
--------

 * Cell ID - Unique cell identification
 * Long axis (L) - Length of the long axis at each time point (px)
 * Short axis - Length of the short axis at each time point (px)
 * Area - Area of the cell at each time point (px)
 * Fluor1 sum - Fluorescence intensity integrated over the cell area at each time point
 * Fluor1 mean - Average (by area) fluorescence intensity for the cell at each time point
 * Long axis length - Length of the long axis at each time point (px)
 * Time (Frames) - Frame reference for each feature (frames)
 * Age (Frames) - Cell age at each time point (frames)
 * Relative Age (Frames) - Cell age at each time point relative to the experiment (frames)
 * Time (fps) - Frame reference for each feature in minutes (min)
 * Age (fps) - Cell age at each time point in minutes (min) 
 * Volume - Cell volume at each time point (:math:`um^3`)
 * fov - The field of view that this cell belongs to
 * F/V - Ratio between Fluorescence and Volume at each time point
 * Cell Cycle - The cell cycle in bins between 0 and 1

.. _SuperSegger: https://github.com/tlo-bot/supersegger-omnipose