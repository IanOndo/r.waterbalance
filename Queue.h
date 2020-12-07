/***************************************************************************************************************************************************************************************************************************
 *
 * MODULE:       r.flowrouting
 *
 * AUTHOR(S):    Ian Ondo
 *              
 * PURPOSE:      Ce programme propose une méthode permettant de modéliser la redistribution d'un flux d'eau le long d'un versant à partir de l'équation d'onde diffusive.
 *				 L'approche consiste à déterminer le temps de trajet d'un point de départ vers un point d'arrivée quelconque situé en aval en suivant un chemin d'écoulement.
 *               Une fonction de réponse basée sur la moyenne et la variance du temps d'écoulement, est modélisée par la fonction de densité du premier temps de passage.
 *               Elle permet de déterminer pour chaque point du paysage la quantité de ruissellement reçu à chaque instant t donné.
 *               Le module calcule pour un pas de temps donné la quantité d'eau drainant depuis chaque pixel vers chaque point situé en aval le long d'un chemin d'écoulement.				 
 *               La sortie du modèle est donc une carte raster représentant à un instant t la redistribution latérale d'un flux d'eau le long d'un versant.
 *
 ************************************************************************************************************************************************************************************************************************/

/***********************************************************************************************
 *
 *				Queue.h
 *				Ce fichier d'en-tête déclare les fonctions, variables et structures des données
 *				utilisées par la fonction principale du programme du module r.flowrouting
 *
 ***********************************************************************************************/
 
#include<stdio.h>
#include<stdlib.h>

#ifndef _QUEUE_H
#define _QUEUE_H

/*
 * Constants
 * ---------
 */

// ERROR_These signal error conditions in queue functions and are used as exit codes for the program.
#define ERROR_QUEUE   2
#define ERROR_MEMORY  3
 
// maxElements represents maximum number of items that can be stored into the queue.
#define maxElements   10000

//
#define priq_purge(q) (q)->n = 1
#define priq_size(q) ((q)->n - 1)


/*Queue has five properties. capacity stands for the maximum number of elements Queue can hold.
  Size stands for the current size of the Queue and elements is the array of elements. front is the
 index of first element (the index at which we remove the element) and rear is the index of last element
 (the index at which we insert the element) */

/*
 * Type: Queue
 * --------------
 * The actual implementation of a queue is completely
 * hidden.  Client will work with queueADT which is a
 * pointer to underlying queueCDT.
 */
typedef void *QueueElements; 
typedef struct Queue
{
        int capacity;
        int size;
        int front;
        int rear;
        QueueElements *elements;
}Queue;

/*
 * Function: CreateQueue
 * Usage: queue = CreateQueue();
 * -------------------------
 * A new empty queue is created and returned.
 */
Queue *CreateQueue(void);

/* Function: QueueDestroy
 * Usage: QueueDestroy(queue);
 * -----------------------
 * This function frees all memory associated with
 * the queue.  "queue" may not be used again unless
 * queue = QueueCreate() is called first.
 */
void DestroyQueue(Queue *Q);

/*
 * Functions: EnQueue, DeQueue
 * Usage: EnQueue(queue, element);
 *        DeQueuequeue);
 * --------------------------------------------
 * These are the fundamental queue operations that enter
 * elements in and delete elements from the queue.  A call
 * to DeQueue() on an empty queue or to EnQueue()
 * on a full queue is an error.  Make use of QueueIsFull()
 * and QueueIsEmpty() (see below) to avoid these errors.
 */
void EnQueue(Queue *Q, QueueElements elements);
void DeQueue(Queue *Q);

/*
 * Functions: Front()
 * Usage: Front(queue);
 * ---------------------
 * Return the element which is at the front */
QueueElements Front(Queue *Q);

/*
 * Functions: QueueIsEmpty, QueueIsFull
 * Usage: if (QueueIsEmpty(queue)) ...
 * -----------------------------------
 * These return a true/false value based on whether
 * the queue is empty or full, respectively.
 */
int QueueIsEmpty(Queue *Q);
int QueueIsFull(Queue *Q);

/* Fonction de réallocation sécurisée */
void *realloc_s(void **ptr, size_t size);


#endif  /* not defined _QUEUE_H */