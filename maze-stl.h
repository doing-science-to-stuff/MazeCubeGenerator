/*
 * maze-stl.h
 * MazeCubeGen: maze cube generator
 *
 * Copyright (c) 2020-2021 Bryan Franklin. All rights reserved.
 */
#ifndef MAZE_STL_H
#define MAZE_STL_H
#include "maze.h"

int maze_export_stl(maze_t *maze, char *filename);
int maze_export_stl_flat(maze_t *maze, char *filename);
int maze_export_stl_solution(maze_t *maze, char *filename);

#endif /* MAZE_H */
