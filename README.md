# MathLinear
Small size linear math library for 3D computer graphics related operations.

[Classes]
vec2d    : 2D vector data structure
vec3d    : 3D vector data structure
vec4d    : 4D vector data structure 
mat3x3   : 3x3 Matrix
mat4x4   : 4x4 Matrix
ray      : Parametric line segment
plane    : 3D plane structure
rectangle: 2D bounding rectangle 
AABB     : Axis-Aligned-Bounding-Box
frustum  : view frustum for visibility testing 
colorRGB : RGB triplet color structure

Note        : All matrix transformations are right handed and row major.
              In order to use them with OpenGL pipeline users need to use
              transpose() method of corresponding matrix class.
