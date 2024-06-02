#ifndef NEIGHBOURLIST_H
#define NEIGHBOURLIST_H

#include <stdlib.h>
#include <stdio.h>
typedef struct Node_list
{
    int id;
    Node_list *next;
} Node_list;

typedef struct neighbourList
{
    Node_list *head;
    Node_list *tail;
} neighbourList;


// Function to create a new Node_list
__device__ Node_list* createNode_list(int id) {
    Node_list* newNode_list = (Node_list*)malloc(sizeof(Node_list));
    if (newNode_list == NULL) {
        printf("Memory allocation createNode_list failed \n");
    }
    newNode_list->id = id;
    newNode_list->next = NULL;
    return newNode_list;
}

__device__ void initNode_list(neighbourList* newNode_list) {
    newNode_list = (neighbourList*)malloc(sizeof(neighbourList));
    if (newNode_list == NULL) {
        printf("Memory allocation initNode_list failed\n");
    }
    newNode_list->head = NULL;
}

// Function to append a Node_list to the end of the neighbourList
__device__ void appendNode_list(neighbourList* neighbourList, int id) {
    Node_list* newNode_list = createNode_list(id);
    if (neighbourList->head == NULL) {
        neighbourList->head = newNode_list;
        neighbourList->tail = newNode_list;
    } else {
        neighbourList->tail->next = newNode_list;
        neighbourList->tail = newNode_list;
        }
}

__device__ int isEmpty(neighbourList* head) {
    return head == NULL;
}

// Function to print the neighbourList
__device__ void printList(neighbourList* neighbourList) {
    Node_list* temp = neighbourList->head;
    while (temp != NULL) {
        printf("%d ", temp->id);
        temp = temp->next;
    }
    printf("\n");
}

// Function to free the neighbourList
__device__ void freeList(neighbourList* neighbourList) {
    Node_list* temp;
    while (neighbourList->head != NULL) {
        temp = neighbourList->head;
        neighbourList->head = neighbourList->head->next;
        free(temp);
    }
    neighbourList->tail = NULL;
}

// Function to initialize the array of lists
neighbourList* initializeArrayOfLists(int size) {
    neighbourList* array = (neighbourList*)malloc(size * sizeof(neighbourList));
    if (array == NULL) {
        printf("Memory allocation for array failed\n");
        exit(1);
    }
    for (int i = 0; i < size; i++) {
        array[i].head = NULL;
        array[i].tail = NULL;
    }
    return array;
}

// void printArrayOfLists(neighbourList* array, int size) {
//     for (int i = 0; i < size; i++) {
//         printf("neighbourList %d:\n", i);
//         printList(&array[i]);
//     }
// }

// Function to free the array of lists
__global__ void freeArrayOfLists(neighbourList* array, int num_molecules) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < num_molecules){
        freeList(&array[index]);
    }
    free(array);
}
#endif