/*
 * maze-export.h
 * MazeCubeGen: maze cube generator
 *
 * Copyright (c) 2020-2023 Bryan Franklin. All rights reserved.
 */
#ifndef MAZE_EXPORT_H
#define MAZE_EXPORT_H
#include "maze.h"

/* export maze shape to an STL file */
int maze_export_stl(maze_t *maze, char *filename);
int maze_export_stl_flat(maze_t *maze, char *filename);
int maze_export_stl_solution(maze_t *maze, char *filename);

/* export maze shape as series of triangles */
void *maze_export_trig_list(maze_t *maze);
void maze_export_trig_list_free(void *list);
int maze_export_num_dims(void *list);
int maze_export_num_trigs(void *list);
float maze_export_vertex_dim(void *list, int trig, int dim);

#endif /* MAZE_EXPORT_H */
