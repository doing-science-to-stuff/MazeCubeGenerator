# Maze Cube Generator

Maze Cube Generator creates random multi-dimensional mazes similar to the
[Oskar's Cube](https://oskarvandeventer.nl/bitsandpieces.html).

**Features:**
 * Custom maze sizes
 * Custom solution lengths
 * Hyper-dimensional mazes
 * Finds solutions
 * STL output (3D only)
 * Unfolded STL for easier 3D printing (3D Only)
 * Multiple metrics to guage maze difficulty

## Building
To build the maze cube generator on Linux or macOS you will first need to
install an appropriate compiler (e.g., [llvm](https://llvm.org/), [gcc](https://gcc.gnu.org/), or [Xcode](https://developer.apple.com/xcode/)) as well as [`cmake`](https://cmake.org/).

Once the build tools are installed, compile Maze Cube Generator with:
 1. `cmake .`
 2. `make`

## Running
Maze Cube Generator is a command-line program with options to control its
behavior.  Below is a list of command-line options and their purpose.

**Basic:**

`-d size` Generate a new maze with given size (`size` format: `l,w,h` or `l,w,h:options`, e.g., `11,11,11` or `11,11,11:r`).<br/>
`-h` Print usage information.<br/>
`-i filename.txt` Load an existing maze from `filename.txt` as produced by `-o`.<br/>
`-m filename.stl` Write maze as STL to `filename.stl`.<br/>
`-f filename.stl` Write flattened maze as STL to `filename.stl`.<br/>
`-u filename.stl` Write maze as more easily printable STL to `filename.stl`.<br/>
`-o filename.txt` Write maze as reloadable (-i) text to `filename.txt`.<br/>
`-p solution.stl` Write solution as STL to `solution.stl` (implies `-s`).<br/>
`-r num` Seed random number generator using `num`.<br/>
`-s` Find a solution to the maze.<br/>
`-x` Scale factor when exporting to STL.<br/>

**Advanced:**

`-l num` Generate mazes until a solution with length &geq;`num` is found.<br/>
`-k num` Generate mazes until at most `num` disconnected regions in maze.<br/>
`-g filename.gv` Writes an input file for [graphviz](https://graphviz.org/) to `filename.gv`.<br/>
`-e width` Specified minimum edge width for faces that are printed flat (-f &
-u) (default: 0.4).

## Examples
Basic random maze generation:<br/>
`$ ./mcg -r 12345 -d 11,11,11 -o output.txt`<br/>
Seeds the random number generator with `12345`, generates a random 11x11x11 maze and stores the results in `output.txt`.

Larger maze with long solution:<br/>
`$ ./mcg -r 1 -d 21,21,21:o -l 70 -k 1 -m maze21.stl -o maze21.txt`<br/>
Seeds the random number generate with `1`, generates random 21x21x21 mazes
until one with a solution of 70 or more is found and all maze positions belong
to the same region.  The resulting maze is then stored in `maze21.stl` and
`maze21.txt`.

Maze generation with optimized start & end locations:<br/>
`$ ./mcg -r 12345 -d 11,11,11:o -o optimal.txt`<br/>
Seeds the random number generator with `12345`, generates a random 11x11x11 maze, finds start & end points that maximize solution length, and stores the results in `optimal.txt`.

Solve an existing maze:<br/>
`$ ./mcg -i unsolved_maze.txt -s -o solved_maze.txt`<br/>
Loads a maze from `unsolved_maze.txt`, finds the solution, and writes the maze
along with the solution to `solved_maze.txt`.

Write 3D model of solution:<br/>
`$ ./mcg -i unsolved_maze.txt -s -p solution.stl`<br/>
Loads a maze from `unsolved_maze.txt`, finds the solution, and writes a 3D
model of the solution to `solution.stl`.

Write unfolded/flattened model:<br/>
`$ ./mcg -i maze.txt -f maze_flat.stl`<br/>
Loads a maze from `maze.txt` and writes a flat-pack style 3D model to
`maze_flat.stl` for fast 3D printing.

Write mostly assembled model:<br/>
`$ ./mcg -i maze.txt -u printable.stl -x 5.0`<br/>
Loads a maze from `maze.txt` and writes a mostly assembled 3D model to
`printable.stl` for easy to assemble 3D printing at 5x scale.

Hyper-dimensional maze:<br/>
`$ ./mcg -r 1 -d 11,11,11,11:o -s -o maze4d.txt`<br/>
Generates a 4-dimensional maze with size 11x11x11x11, and writes the maze
along with its solution to `maze4d.txt`.
