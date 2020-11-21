/*
 * maze.h
 * MazeCubeGen: maze cube generator
 *
 * Copyright (c) 2020 Bryan Franklin. All rights reserved.
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

typedef int* position_t;
typedef struct position_list {
    int numDimensions;
    int posListCap;
    int posListNum;
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
} maze_t;

int maze_init(maze_t *maze, int numDimensions, int *sizes);
int maze_free(maze_t *maze);
void maze_set_segments(maze_t *maze, int segs);
void maze_set_path_length(maze_t *maze, int len);
int maze_pick_goals(maze_t *maze);
int maze_solve(maze_t *maze);
int maze_generate(maze_t *maze);
int maze_write(maze_t *maze, char *filename);
int maze_load(maze_t *maze, char *filename);
int maze_export_stl(maze_t *maze, char *filename);
int maze_export_stl_solution(maze_t *maze, char *filename);

#endif /* MAZE_H */
