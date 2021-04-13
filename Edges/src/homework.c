#include "fem.h"


# ifndef NOEXPAND

void edgesExpand(femEdges *theEdges)
{
    femEdge* edge_list = theEdges->edges;
    int nLocal = theEdges->mesh->nLocalNode;
    for (size_t num_elem = 0; num_elem < theEdges->mesh->nElem; num_elem++){
        for (size_t num_node = 0; num_node < theEdges->mesh->nLocalNode; num_node++){
            if(theEdges->edges==NULL){printf("crash malloc in edges expand for i:%ld j:%ld\n",num_elem,num_node);}
            edge_list->elem[0] = num_elem;
            edge_list->elem[1] = -1;
            //printf("elem cava for i:%ld j:%ld\n",num_elem,num_node);
            edge_list->node[0] = theEdges->mesh->elem[num_elem * nLocal + num_node];
            edge_list->node[1] = theEdges->mesh->elem[num_elem * nLocal + ((num_node+1)%nLocal)];
            edge_list++;
        }
    }
}

# endif
# ifndef NOSORT

void edgesSort(femEdges *theEdges)
{
	qsort(theEdges->edges, theEdges->nEdge, sizeof(femEdge), edgesCompare);
}

# endif
# ifndef NOCOMPARE

int edgesCompare(const void* e0, const void *e1)
{   
    femEdge* edge0 = (femEdge*) e0;
    femEdge* edge1 = (femEdge*) e1;
    int min_edge0 = edge0->node[0] < edge0->node[1] ? edge0->node[0] : edge0->node[1];
    int min_edge1 = edge1->node[0] < edge1->node[1] ? edge1->node[0] : edge1->node[1];
    //printf("e0 n0:%d n1:%d index_min_edge0:%d // e1 n0:%d n1:%d index_min_edge1:%d\n",edge0->node[0],edge0->node[1],min_edge0,edge1->node[0],edge1->node[1],min_edge1);


    if(min_edge0 < min_edge1){ 
        return -1;
    }else if (min_edge0 > min_edge1){ 
        return 1;
    }else{
        int second_e0 = edge0->node[0] == min_edge0 ? edge0->node[1] : edge0->node[0];
        int second_e1 = edge1->node[0] == min_edge1 ? edge1->node[1] : edge1->node[0];
        if(second_e0 < second_e1){
            return -1;
        }else if(second_e0 > second_e1){
            return 1;
        }
    }
    return 0;
}

# endif
# ifndef NOSHRINK

void edgesShrink(femEdges *theEdges)
{
    int n = theEdges->nEdge;
    int nBoundary = theEdges->nEdge;
    for (size_t i = 0; i < theEdges->nEdge -1; i++)
    {
        if((theEdges->edges[i].node[0] == theEdges->edges[i+1].node[1] && theEdges->edges[i].node[1] == theEdges->edges[i+1].node[0])
            || (theEdges->edges[i].node[0] == theEdges->edges[i+1].node[0] && theEdges->edges[i].node[1] == theEdges->edges[i+1].node[1]))
        {
                theEdges->edges[i].elem[1] = theEdges->edges[i+1].elem[0];
                theEdges->edges[i+1].elem[0] = -1;
                n--;
                nBoundary-=2;
        }
    }
    
    femEdge* new_edges = theEdges->edges;
    int index_edges = 0;
    for (size_t i = 0; i < theEdges->nEdge; i++)
    {
        if(theEdges->edges[i].elem[0]!=-1){
            new_edges[index_edges] = theEdges->edges[i];
            index_edges++;
        }
    }
    
    // Reallocation du tableau des edges
    //theEdges->edges = realloc(theEdges->edges, n * sizeof(femEdge));
    theEdges->edges = realloc(new_edges,sizeof(femEdge)*n);
    theEdges->nEdge = n;
    theEdges->nBoundary = nBoundary;
    femEdgesPrint(theEdges);
}

# endif
# ifndef NOBOUNDARYLENGTH

double edgesBoundaryLength(femEdges *theEdges)
{
    double L = 0;
    for (size_t i_edge = 0; i_edge < theEdges->nEdge; i_edge++)
    {
        if(theEdges->edges[i_edge].elem[1]==-1){
            int n0 = theEdges->edges[i_edge].node[0];
            int n1 = theEdges->edges[i_edge].node[1];
            L+=sqrtf(
                powf(theEdges->mesh->X[n0] - theEdges->mesh->X[n1],2) +
                powf(theEdges->mesh->Y[n0] - theEdges->mesh->Y[n1],2)
            );
        }
    }
    
    return L;
}

# endif
