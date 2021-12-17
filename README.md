# Chicken! Triangulate!

A very _very_ basic Delaunay triangulator in pure scheme
Follows the algoritm of ![paul burke](http://paulbourke.net/papers/triangulate/)

# Usage

(draw-points (gen-rand-points 10 500))

You can also read vertex data from an obj file as
(read-obj-file "path/to/file")

To just triangulate and get a list of the triangles use
(triangulate vertex-vector)

# Requirements

	srfi-133 (for vector ops)
	srfi-1   (for filter)
	sdl2     (for graphics)

# Triangles!

![some](https://github.com/ramrunner/chicken-triangulate/blob/main/img/10p.png?raw=true)
![triangles](https://github.com/ramrunner/chicken-triangulate/blob/main/img/50p.png?raw=true)
