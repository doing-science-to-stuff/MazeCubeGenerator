/*
 * main.c
 * MazeCubeGen: maze cube generator
 *
 * Copyright (c) 2020-2024 Bryan Franklin. All rights reserved.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "maze.h"
#include "maze-export.h"

static char *configStr = NULL;
static char *inputFile = NULL;
static char *outputFile = NULL;
static char *stlFile = NULL;
static char *stlPrintDir = NULL;
static char *stlFileFlat = NULL;
static char *stlSolFile = NULL;
static char *gvFile = NULL;

static void free_buffers() {
    if(configStr!=NULL) free(configStr);
    if(inputFile!=NULL) free(inputFile);
    if(outputFile!=NULL) free(outputFile);
    if(stlFile!=NULL) free(stlFile);
    if(stlPrintDir!=NULL) free(stlPrintDir);
    if(stlFileFlat!=NULL) free(stlFileFlat);
    if(stlSolFile!=NULL) free(stlSolFile);
    if(gvFile!=NULL) free(gvFile);
}


static void show_help(int argc, char **argv) {
    free_buffers();
    printf("%s {-d l,w,h[:olb] | -i file.txt} [-h] [-s] [-r seed] [-l length] {-m file.stl | -o file.txt} [-c] [-f flat.stl] [-u printable.stl] [-k maxSegments] [-p solution.stl] [-g graph.gv] [-e edgeWidth] [-x scale]\n\n",
            argv[0]);
    exit(0);
}


int main(int argc, char **argv) {
    maze_t maze;
    int genMaze = 0;
    int doSolve = 0;
    int minSolutionLen = -1;
    int seed = 0;
    int useSeed = 0;
    int maxSegments = -1;
    maze_output_opts_t maze_opts;

    /* set default values */
    maze_opts.scale = 1.0;
    maze_opts.edgeWidth = 0.4;
    maze_opts.separateMarkers = 1;

    /* set defaults */
    genMaze = 0;
    seed = 0;
    useSeed = 0;

    char ch='\0';
    /* remaining letters: abjnqtvwyz */
    while( (ch=getopt(argc, argv, "cd:e:f:g:hi:k:l:m:o:p:r:su:x:"))!=-1 ) {
        switch(ch) {
            case 'c':
                /* combine markers into face */
                maze_opts.separateMarkers = 0;
                break;
            case 'd':
                configStr = strdup(optarg);
                genMaze = 1;
                break;
            case 'e':
                maze_opts.edgeWidth = atof(optarg);
                break;
            case 'f':
                /* output flat STL model to file given as argument */
                stlFileFlat = strdup(optarg);
                break;
            case 'g':
                /* output graphviz gv to file given as argument */
                gvFile = strdup(optarg);
                break;
            case 'h':
            case '?':
                show_help(argc,argv);
                break;
            case 'i':
                /* load maze from file given as argument */
                inputFile = strdup(optarg);
                genMaze = 0;
                break;
            case 'l':
                /* minimum solution length */
                minSolutionLen = atoi(optarg);
                break;
            case 'k':
                maxSegments = atoi(optarg);
                break;
            case 'm':
                /* output STL model to file given as argument */
                stlFile = strdup(optarg);
                break;
            case 'u':
                /* output maze with the ends removed into directory */
                stlPrintDir = strdup(optarg);
                break;
            case 'o':
                /* output maze to file given as argument */
                outputFile = strdup(optarg);
                break;
            case 'p':
                /* output solution STL model to file given as argument */
                stlSolFile = strdup(optarg);
                break;
            case 'r':
                seed = atoi(optarg);
                useSeed = 1;
                break;
            case 's':
                /* generate a solution */
                doSolve = 1;
                break;
            case 'x':
                maze_opts.scale = atof(optarg);
                break;
        }
    }

    /* check that sufficient settings exist */
    if( inputFile==NULL && !genMaze ) {
        fprintf(stderr,"\n\nNo maze source given, use -i to specify an input filename or -g to generate a random maze.\n\n");
        show_help(argc,argv);
    }
    if( outputFile==NULL && stlFile==NULL && stlPrintDir==NULL && stlSolFile==NULL && stlFileFlat==NULL && gvFile==NULL ) {
        fprintf(stderr,"\n\nNo output specified, use -o, -m, and/or -u to specify an output filename.\n\n");
        show_help(argc,argv);
    }

    /* generate/load maze and solve if requested/needed */
    if( genMaze ) {
        /* create a new maze */
        if( configStr )
            maze_init_str(&maze, configStr);
        else
            maze_init_str(&maze, "11x11x11");
        if( maxSegments > 0 )
            maze_set_segments(&maze, maxSegments);
        if( minSolutionLen > 0 )
            maze_set_path_length(&maze, minSolutionLen);
        printf("Generating maze.\n");
        if( useSeed ) {
            printf("  using %i as the random seed.\n", seed);
            srand(seed);
        }
        int segments = maze_generate(&maze);
        if( maxSegments > 0 )
            printf("%i maze segment%s.\n", segments, (segments==1)?"":"s");
    }
    else if( inputFile ) {
        /* load a maze from an input file */
        printf("Loading maze from '%s'.\n", inputFile);
        if( maze_load(&maze,inputFile) < 0 ) {
            fprintf(stderr, "Failed to load maze from '%s'\n", inputFile);
            free_buffers();
            exit(1);
        }
    }

    /* solve the maze, if requested/needed */
    if( doSolve || minSolutionLen>0 ) {
        maze_solve(&maze);
        if( minSolutionLen>0 )
            printf("solution length is %i.\n", maze.solution.num);
    }

    /* write requested output */
    if( outputFile ) {
        printf("Writing %iD maze to `%s`.\n", maze.numDimensions, outputFile);
        maze_write(&maze, outputFile);
    }
    if( stlFile ) {
        printf("Exporting %iD maze to STL file `%s`.\n", maze.numDimensions, stlFile);
        maze_export_stl(&maze, stlFile, &maze_opts);
    }
    if( stlPrintDir ) {
        printf("Exporting %iD maze to easily to assemble STL files into `%s`.\n", maze.numDimensions, stlPrintDir);
        maze_export_stl_printable(&maze, stlPrintDir, &maze_opts);
    }
    if( stlSolFile ) {
        printf("Exporting %iD maze solution to STL file `%s`.\n", maze.numDimensions, stlFile);
        maze_export_stl_solution(&maze, stlSolFile, &maze_opts);
    }
    if( stlFileFlat ) {
        printf("Exporting %iD maze to flat-packed STL file `%s`.\n", maze.numDimensions, stlFileFlat);
        maze_export_stl_flat(&maze, stlFileFlat, &maze_opts);
    }
    if( gvFile ) {
        printf("Exporting %iD maze to graphviz gv file `%s`.\n", maze.numDimensions, gvFile);
        maze_export_gv(&maze, gvFile);
    }

    maze_metrics(&maze);

    /* free memory */
    maze_free(&maze);
    free_buffers();

    return 0;
}
