#ifndef MAZE_H
#define MAZE_H

typedef struct maze {
    int numFaces;
    int *dimensions;
    int ***cells;
} maze_t;

int maze_init(maze_t *maze, int numDimensions, int *sizes);
int maze_generate(maze_t *maze);
int maze_export_stl(maze_t *maze, char *filename);


#endif // MAZE_H
