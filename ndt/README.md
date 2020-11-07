# NDT scene

This directory contains a scene file for rendering mazes using the
[ndt](https://github.com/doing-science-to-stuff/ndt) n-dimensional ray-tracer.

### Rendering a maze using NDT
 1. Place the file `maze-cube.c` into the `scenes` directory of ndt.
 2. Rebuild ndt (`cmake .; make`).
 3. Run ndt (`./ndt -s scenes/maze-cube.so -u PathToMazeFile`)

Where `PathToMazeFile` is a file written using the `-o` flag of the puzzle maze
generator (pmg).

When rendering mazes with more than three spatial dimensions, the `-d` flag
must also be used to indicate the number of dimensions for the maze.
