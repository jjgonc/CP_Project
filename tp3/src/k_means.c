#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <math.h>
#include<omp.h>
#include<mpi.h>
#include "../include/utils.h"

int N = 10000000;
int K = 4;
int New_N = 0;

/**
 * @brief Initializes cluster_point struct and initializes with a random value for its center
 * 
 * @param pointX X coordinate of the point
 * @param pointY Y coordinate of the point
 * @param centerX X coordinate of the center
 * @param centerY Y coordinate of the center
 * @param cluster_size Size of the cluster
*/
void inicializa(float *pointX, float *pointY, float *centerX, float *centerY, int *cluster_size, int nprocesses, int rank){
    srand(10);

    for(int i = 0; i < New_N; i++) {        //assigns random values to the coordinates
        pointX[i] = (float)rand()/ RAND_MAX;
        pointY[i] = (float)rand()/ RAND_MAX;
    }
    if(rank == 0){
        for(int i = 0; i < K; i++) {    // assigns the first K points as the initial centers
            centerX[i] = pointX[i];
            centerY[i] = pointY[i];
            cluster_size[i] = 0;
        }
    }
    //broadcast the values of the centers and the size of the clusters
    MPI_Bcast(&(*centerX), K, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(*centerY), K, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(*cluster_size), K, MPI_INT, 0, MPI_COMM_WORLD);
}


/**
 * @brief Reevaluate centroides. Calculates the mean of each cluster and assigns that value to the cluster center
 * 
 * @param centerX X coordinate of the center
 * @param centerY Y coordinate of the center
 * @param cluster_size Size of the cluster
 * @param sum_x Sum of the X coordinates of the points in the cluster
 * @param sum_y Sum of the Y coordinates of the points in the cluster
 */
void reevaluate_centers(float *centerX, float *centerY, int *cluster_size, float *sumX, float *sumY) {
    for(int i = 0; i < K; i++){
        centerX[i] = sumX[i] / cluster_size[i];
        centerY[i] = sumY[i] / cluster_size[i];
    }
}

/**
 * @brief Initializes the centers and the size of the clusters to 0
 * @param centerX array of the X coordinates of the centers
 * @param centerY array of the Y coordinates of the centers
 * @param cluster_size array of the size of the clusters
*/
void initialize_centers(float *centerX, float *centerY, int* size){
    for(int i = 0; i < K; i++) {  
        centerX[i] = 0;
        centerY[i] = 0;
        size[i] = 0;
    }
}
/**
 * @brief Calculates the distance between two points withouth using sqrt
 * 
 * @param pointX X coordinate of the point
 * @param pointY Y coordinate of the point
 * @param centerX X coordinate of the center
 * @param centerY Y coordinate of the center
 * @return float Distance between two points
 */
float distance(float pointX, float pointY, float centerX, float centerY) {
    return ((pointX - centerX) * (pointX - centerX)) + ((pointY - centerY) * (pointY - centerY));
}

/**
 * @brief Calculates new points and checks if the cluster center has changed
 * 
 * @param cluster Cluster to which the point belongs
 * @param pointX X coordinate of the point
 * @param pointY Y coordinate of the point
 * @param centerX X coordinate of the center
 * @param centerY Y coordinate of the center
 * @param cluster_size Size of the cluster
 * @return int  1 if the cluster center has changed, 0 otherwise
 */
int cluster_points(int* cluster, float *pointX, float *pointY, float *centerX, float *centerY, int *size, int nprocesses, int rank) {
    int changed = 0;
    float sum_x[K], sum_y[K];
    initialize_centers(sum_x, sum_y, size);
    for(int i = 0; i < New_N; i++){
        //calculates the distance to the center of each cluster and saves the cluster with the smallest distance
        int min_cluster = 0;
        float min_dist = distance(pointX[i], pointY[i], centerX[0], centerY[0]);
        for(int j = 1; j < K; j++) {
            if(distance(pointX[i], pointY[i], centerX[j], centerY[j]) < min_dist) {
                min_cluster = j;
                min_dist = distance(pointX[i], pointY[i], centerX[j], centerY[j]);
            }
        }
        //Verify if the cluster has changed and if so, change the cluster and the changed variable
        if(cluster[i] != min_cluster) {
            changed += 1;
            cluster[i] = min_cluster;
        }
        //sum the coordinates of the points in each cluster
        sum_x[min_cluster] += pointX[i];
        sum_y[min_cluster] += pointY[i];
        size[min_cluster]++;
    }

    float aux_sum_x[K], aux_sum_y[K];
    int aux_changed = 0, aux_size[K];

    //sum the values of the variables changed and sum_x and sum_y and the size of the clusters
    MPI_Reduce(&changed, &aux_changed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sum_x, &aux_sum_x, K, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sum_y, &aux_sum_y, K, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(*size), &aux_size, K, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    

    if (rank == 0){
        for(int i = 0; i < K; i++) {
            size[i] = aux_size[i];
        }
        reevaluate_centers(centerX, centerY, size, aux_sum_x, aux_sum_y);
    }
    //broadcast the new values of the centers and the changed variable
    MPI_Bcast(&(*centerX), K, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(*centerY), K, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&changed, 1, MPI_INT, 0, MPI_COMM_WORLD);

    return changed;
}


/**
 * @brief Prints the cluster centers
 * 
 * @param iter The iteration number
 * @param centerX X coordinate of the center
 * @param centerY Y coordinate of the center
 * @param size The size of the cluster
 */
void print_clusters(int iter, float *centerX, float *centerY, int *size) {
    for(int i = 0; i < K; i++) {
        printf("Center: (%.3f, %.3f)", centerX[i], centerY[i]);
        //printf("Center: (%.3f, %.3f) ", center[i].x, center[i].y);
        printf("Size : %d \n", size[i]);
    }
    printf("Iterations: %d \n", iter);
}

/**
 * @brief Function that joins all the steps of the k-means algorithm
 * 
 */
void k_means(int nprocesses, int rank) {
    //Initialize variables
    int size[K];
    float centerX[K], centerY[K];
    New_N = N / nprocesses;
    if (rank == nprocesses - 1) {
        New_N += N % nprocesses;
    }
    float *pointsX = (float *)malloc(New_N * sizeof(float));
    float *pointsY = (float *)malloc(New_N * sizeof(float));
    int *cluster = (int *)malloc((N / nprocesses) * sizeof(int));
    inicializa(pointsX, pointsY, centerX, centerY, size, nprocesses, rank);

    // start iteration counter
    int iter = 0;
    
    while (iter < 20 && cluster_points(cluster, pointsX, pointsY, centerX, centerY, size, nprocesses, rank) > 0) {
        iter++;
    }
    if (rank == 0) {
        print_clusters(iter, centerX, centerY, size);
    }

    // free memory
    free(pointsX);
    free(pointsY);
    free(cluster);
}

/**
 * @brief Main function
 * 
 * @return int Returns 0 if the program runs successfully
 */
int main(int argc, char *argv[]) {
    int nprocesses, myrank;
    if (argc >= 2) {
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocesses);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

        N = atoi(argv[1]);
        K = atoi(argv[2]);
    }else{
        printf("Usage: %s <Number of points> <Number of clusters> -np <Number of processes>\n", argv[0]);
        exit(1);
    }
    k_means(nprocesses, myrank);
    MPI_Finalize();
    return 0;
}
