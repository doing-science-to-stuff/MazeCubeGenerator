# NDT scene

This directory contains a scene plugin for rendering mazes using the
[ndt](https://github.com/doing-science-to-stuff/ndt) n-dimensional ray-tracer.

### Rendering a maze using NDT
 1. Place the file `maze-cube.c` into the `scenes` directory of ndt.
 2. Update path to `maze.c` near the top of `maze-cube.c`.
 3. Rebuild ndt (`cmake .; make`).
 4. Run ndt (`./ndt -s scenes/maze-cube.so -u PathToMazeFile`)

Where `PathToMazeFile` is a text file written using the `-o` flag of the maze
cube generator (mcg).  For best results, the maze file should have been
written with a solution using the `-s` flag.

When rendering mazes with more than three spatial dimensions, the `-d` flag
must also be used to indicate the number of dimensions for the maze.
