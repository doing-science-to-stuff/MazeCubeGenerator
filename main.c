#include <stdio.h>
#include "maze.h"

int main(int argc, char **argv) {
    maze_t maze;
    int sizes[3] = {50, 50, 50};
    maze_init(&maze, 3, sizes);

    maze_generate(&maze);

    maze_export_stl(&maze, "output.stl");

    return 0;
}
