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
 *				Queue.c
 *				Ce fichier d'en-tête déclare les fonctions, variables et structures des données
 *				utilisées par la fonction principale du programme du module r.flowrouting
 *
 ***********************************************************************************************/

/*The Queue has five properties - capacity stands for the maximum number of elements Queue can
hold, Size stands for the current size of the Queue, elements is the array of elements, front is the
index of first element (the index at which we remove the element) and rear is the index of last
element (the index at which we insert the element). Functions on Queue

1. createQueue function takes argument the maximum number of elements the Queue can
hold, creates a Queue according to it and returns a pointer to the Queue. It initializes Q-
>size to 0, Q->capacity to maxElements, Q->front to 0 and Q->rear to -1.

2. EnQueue function - This function takes the pointer to the top of the queue Q and the item
(element) to be inserted as arguments. Check for the emptiness of queue
a. If Q->size is equal to Q->capacity, we cannot push an element into Q as there is
no space for it.
b. Else, EnQueue an element at the end of Q, increase its size by one. Increase the
value of Q->rear to Q->rear + 1. As we fill the queue in circular fashion, if
Q->rear is equal to Q->capacity make Q->rear = 0. Now, Insert the element in its
rear side

Q->elements[Q->rear] = element

3. DeQueue function - This function takes the pointer to the top of the stack S as an
argument.
a. If Q->size is equal to zero, then it is empty. So, we cannot DeQueue.
b. Else, remove an element which is equivalent to incrementing index of front by
one. Decrease the size by 1. As we fill elements in circular fashion, if Q->front is
equal to Q->capacity make Q->front=0.

4. front function – This function takes the pointer to the top of the queue Q as an argument
and returns the front element of the queue Q. It first checks if the queue is empty
(Q->size is equal to zero). If it’s not it returns the element which is at the front of the
queue.

Q->elements[Q->front] */
 
#include <stdio.h>
#include <stdlib.h>
#include "Queue.h"

/* createQueue function takes argument the maximum number of elements the Queue can hold, creates
   a Queue according to it and returns a pointer to the Queue. */
   
Queue *CreateQueue(void)
{
        /* Create a Queue */
        Queue *Q;
        Q = (Queue *)malloc(sizeof(Queue));
		
		if (Q == NULL) {
		fprintf(stderr, "Insufficient Memory for new Queue.\n");
		exit(ERROR_MEMORY);  /* Exit program, returning error code. */
		}
				
        /* Initialise its properties */
        Q->elements = (QueueElements *)malloc(sizeof(QueueElements)*maxElements);
        Q->size     = 0;
        Q->capacity = maxElements;
        Q->front    = 0;
        Q->rear     = -1;
		
        /* Return the pointer */
        return Q;
};

void DestroyQueue(Queue *Q)
{
    if (Q == NULL) {
        printf("Queue is already NULL. Nothing to free.");
		exit(1);
    };
  free(Q->elements);
  free(Q);
};

void DeQueue(Queue *Q)
{
        /* If Queue size is zero then it is empty. So we cannot pop */
        if(QueueIsEmpty(Q))
        {
                printf("Queue is Empty\n");
                return;
        }
        /* Removing an element is equivalent to incrementing index of front by one */
        else
        {
				Q->size--;
                Q->front++;
				
                /* As we fill elements in circular fashion */
                if(Q->front==Q->capacity)
                {
                        Q->front=0;
                }
        }
        return;
};

QueueElements Front(Queue *Q)
{
        if(QueueIsEmpty(Q))
        {
                printf("Queue is Empty\n");
                exit(ERROR_QUEUE);
        }
        /* Return the element which is at the front*/
        return Q->elements[Q->front];
};

void EnQueue(Queue *Q, QueueElements element)
{
	int i;
        /* If the Queue is full, we cannot push an element into it as there is no space for it.*/
        if(QueueIsFull(Q))
        {			
                printf("Queue is Full\n");
				exit(ERROR_QUEUE);
		}
        else
        {
                Q->size++;
                Q->rear = Q->rear + 1;
                /* As we fill the queue in circular fashion */
                if(Q->rear == Q->capacity)
                {
                        Q->rear = 0;
                }
                /* Insert the element in its rear side */ 
                Q->elements[Q->rear] = element;
        }
        return;
};

int QueueIsEmpty(Queue *Q)
{
	return(Q->size == 0);
};

int QueueIsFull(Queue *Q)
{
	return(Q->size == Q->capacity);
};

void *realloc_s(void **ptr, size_t size)
{
	void *ptr_realloc = realloc(*ptr, size);
	
	if(ptr_realloc!=NULL)
	*ptr = ptr_realloc; /* Meme si ptr_realloc est nul on ne vide pas la mémoire. On laisse l'initiative au programme */
	
	return ptr_realloc;
};
