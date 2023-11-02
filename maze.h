/*
 * maze.h
 * MazeCubeGen: maze cube generator
 *
 * Copyright (c) 2020-2021 Bryan Franklin. All rights reserved.
 */
#ifndef MAZE_H
#define MAZE_H

typedef struct maze_face {
    int d1, d2;
    int rows, cols;
    char *cells;
} face_t;

int face_init(face_t *face, int *sizes, int d1, int d2);
int face_free(face_t *face);
int face_get_cell(face_t *face, int row, int col);

typedef int* position_t;
int position_init(position_t *pos, int dimensions);
int position_copy(position_t *dst, position_t *src, int dimensions);
int position_compare(position_t *pos1, position_t *pos2, int dimensions);
int position_free(position_t *pos);

typedef struct position_list {
    int numDimensions;
    int capacity;
    int num;
    position_t *positions;
} position_list_t;

typedef struct maze {
    int numFaces;
    int numDimensions;
    int *dimensions;
    int maxSegments;
    int minPathLength;
    position_t startPos;
    position_t endPos;
    position_list_t reachable;
    position_list_t solution;
    face_t *faces;
    char pathSelMode;
} maze_t;

int maze_init_str(maze_t *maze, char *dimStr);
int maze_init(maze_t *maze, int numDimensions, int *sizes, char *options);
int maze_free(maze_t *maze);
void maze_set_segments(maze_t *maze, int segs);
void maze_set_path_length(maze_t *maze, int len);
int maze_pick_goals(maze_t *maze);
int maze_solve(maze_t *maze);
int maze_generate(maze_t *maze);
int maze_write(maze_t *maze, char *filename);
int maze_load(maze_t *maze, char *filename);
int maze_export_gv(maze_t *maze, char *filename);
int maze_metrics(maze_t *maze);

#endif /* MAZE_H */
