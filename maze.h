#ifndef MAZE_H
#define MAZE_H

typedef struct maze_face {
    int d1, d2;
    int rows, cols;
    char *cells;
} face_t;

int face_init(face_t *face, int *sizes, int d1, int d2);

typedef struct maze {
    int numFaces;
    int numDimensions;
    int *dimensions;
    face_t *faces;
} maze_t;

int maze_init(maze_t *maze, int numDimensions, int *sizes);
int maze_generate(maze_t *maze);
int maze_export_stl(maze_t *maze, char *filename);

#endif // MAZE_H
