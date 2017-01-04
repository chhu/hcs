# The H-coordinate system
This include / template library is an extended, multi-dimensional recursive coordinate system, similar to the Z order space-filling curve or Morton-codes. It is more than a Z curve as it captures the recursive nature best visualized from the H fractal. Each recursive refinement is called level. The coordinate storage type canbe changed here and is intentionally not a template parameter. Also the part that stores the level can be configured (HCS_LEVEL_BIT) but you should know what you are doing.

The default coordinate type is unsigned 64bit, allowing a 3D recursion depth of 19 iterations (levels), resulting in a closest distance of 1/2^19 for a scale of one (1x1x1 box). The most significant bits determine the (iteration-) level, _the_ most significant bit marks a boundary coord. The lower bits describe sub-coordinates, every D (dimensions) bits. The sub-coordinate of the highest level is always at the least significant bits.

Because the H fractal best illustrates this, its called the H-coordinate-system (https://en.wikipedia.org/wiki/H_tree):


         6|        |7
          |        |
          |--------|
          |   /    |
      2| 4|  #  |3 |5
       |    /   |
       |--------|
       |        |
      0|        |1

   A single coordinate for a level (3-bit, 0-7 values for 3D) is as follows: LSB is X+/-, next bit Y+/-, next bit Z+/-
   Example:
   
   
  	    Level-1 coordinate:                                                  X        Y      Z
  	                    L1
  	     0b0000 01 ... 001            Would result in position (center) + (+scale, -scale, -scale) / 2
  	                   ZYX
  	              ^ HCS_LEVEL_BIT
  
  	    Level-2 coordinate:
  	                    L1  L2
  	     0b0000 10 ... 011 110        Would result in position (center) + (+scale, -scale, +scale) / 2 + (-scale, +scale, +scale) / 4
  	                   ZYX ZYX
  
Again: The LSBs (Least Significant Bit) are always the ones for the highest level (from level part). This guarantees shortest linear distance in memory for neighbors and is compatible with morton-codes.

The intention behind the HCS and the storage class Field is of numerical nature, ultimately to create a algebraic-mesh refinement multi-grid solver for PDEs. 

