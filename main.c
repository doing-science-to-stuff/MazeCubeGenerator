#include <stdio.h>
#include "maze.h"

int main(int argc, char **argv) {
    printf("%s starting.\n", argv[0]);
    maze_t maze;
    int sizes[3] = {50, 50, 50};
    maze_init(&maze, 3, sizes);

    maze_generate(&maze);

    maze_export_stl(&maze, "output.stl");
    printf("%s finished.\n", argv[0]);

    return 0;
}
