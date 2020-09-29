/*
 * main.c
 * PuzzleMaze: A 3D puzzle model generator
 *
 * Copyright (c) 2020 Bryan Franklin. All rights reserved.
 */
#include <stdio.h>
#include "maze.h"

int main(int argc, char **argv) {
    maze_t maze;
    int sizes[3] = {21, 21, 21};
    maze_init(&maze, 3, sizes);

    maze_generate(&maze);

    maze_export_stl(&maze, "output.stl");

    maze_free(&maze);

    return 0;
}
