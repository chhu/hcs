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

## Navigating through the HCS
### Scaled and unscaled positions
The HCS facility provides first the ability to encode / decode Cartesian positions into Morton-codes. The center and scales for each dimension can be changed, they default to a 1x1x1x... box originating at (0,0,0,...). The `createFromPosition(level, pos)` and `getPosition(coord)` do exactly that. Creating a coord from a position requires to specify which level the resulting coord should have. This does not mean the resulting coordinate will be at exactly the provided position, but instead at the one closest to the provided, which will be closer the higher the level. Unscaled means they operate on integers representing the whole level (`getUnscaled(coord)` and `createFromUnscaled(level, unscaled)`). A level-8 coordinate in 2D has 2^8 x 2^8 possible locations, so an unscaled level-8 coordinate would be in _unscaled_ Cartesian space from X= 0 -> 255 and Y=0 -> 255, while a level-9 has 2^9, so X= 0 -> 511 and Y= 0 -> 511. An unscaled coordinate has a completely different "true" location than the same coordinate at a different level! Nevertheless, they are fast and useful.
### Neighbors
An important task for numerical application is to find neighboring coordinates for stencils. For Cartesian coordinates this is straight forward, subtracting or adding one to the coordinates. For Morton-codes this is tricky. The fastest way seems to follow the [Moser de Bruijn sequence](https://en.wikipedia.org/wiki/Moser%E2%80%93de_Bruijn_sequence). To stay dimensionally independent, the forward or backward direction for each dimension are encoded like this: The least-significant bit of the direction parameter is 0 for forward / positive direction and 1 for negative direction along the dimension, which is encoded in the following bits. `getNeighbor(coord, direction)` is the method to use for that. It will return the neighboring coordinate in the desired direction for the same level as coord. 
### Boundaries
The neighbor algorithm takes care if you request a boundary. It encodes this in the resulting coord, marking it as boundary, which boundary was hit, using the same encoding as direction, so for example the `0` boundary would mean X+ (or right), `1` X-  (left), `2` Y+ (or upper) and so forth. It also leaves the coordinate that lead to that boundary intact. A `Field()` has an array of lambdas to provide values for the boundaries upon request. These lambdas are called with the boundary coordinate so they can reveal the position of the "hitting" coordinate to provide position-dependent boundary values. 
## Storage
The `Field()` class provides a sparse storage object for arbitrary data types. The following features are provided:

 - dedicated refinement / coarsening
 - lower-level coords always exist, but top-level get marked as such. (Top-Level-Coordinate = TLC)
 - iterator class that allows fast iteration over all top-level or all existing coords or existing coords of a specific level
 - bi-linear interpolation of non-existing coords, providing coefficients for TLC (preserving divergence of a vector field)
 - Supports arbitrary data types that need to support some basic arithmetic operations
 - Field objects supply basic arithmetic operators (!) This means you can add / multiply / ... two fields with different structures, the resulting Field will have the structure of the first operand, all missing coords from the second field will be interpolated 
 - A bracket operator for coordinates is implemented, with adjustable behavior for non-existing coords.
 - Convert and merge methods implemented with lambdas to change the data-type of a Field. 
 - The performance of exists() relies on STL's map::lower_bound O(log) complexity
 - The center coordinate (0) always exists
 - boundary conditions can be implemented as lambdas, defaulting to zero
_Look at the tests on how to use these features!_

## Goal of this project
The way the interpolation works enables finite-difference / -volume methods to use mesh-refinement, with the ultimate goal to build a CFD Multi-Grid solver with adaptive mesh refinement. Some day...

