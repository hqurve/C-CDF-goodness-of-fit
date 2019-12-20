#ifndef __RESOURCE_STACK__
#define __RESOURCE_STACK__

#include <stdio.h>
#define MAX_RESOURCE_COUNT 20

#define INTIALIZE_RESOURCE_MANAGER struct resource_stack resource_stack; resource_stack.size = 0;
#define R_EXIT_FAILURE freeAllResources(&resource_stack); return 1;
#define R_EXIT_SUCCESS freeAllResources(&resource_stack); return 0;

#define FREE_RESOURCE(resource) freeResource(&resource_stack, resource)
#define PUSH_RESOURCE(resource, funcFree) if(!pushResource(&resource_stack, resource, funcFree)){printf("Resource stack overflow"); R_EXIT_FAILURE}
#define PUSH_RESOURCE_MESSAGE(resource, funcFree, message) if(resource == NULL){printf(message); R_EXIT_FAILURE;}else{PUSH_RESOURCE(resource, funcFree);}


struct resource{
    void* item;
    void (*funcFree)(void*);
};

struct resource_stack{
    int size;
    struct resource stack[MAX_RESOURCE_COUNT];
};

int pushResource(struct resource_stack* resources, void* resource, void (*funcFree)(void*)){
    if (resources->size >= MAX_RESOURCE_COUNT){
        return 0;
    }
    resources->stack[resources->size].item = resource;
    resources->stack[resources->size].funcFree = funcFree;
    return 1;
}
int freeResource(struct resource_stack* resources, void* resource){
    for (int i = 0; i<resources->size; i++){
        if (resources->stack[i].item == resource){
            resources->stack[i].funcFree(resources->stack[i].item);
            resources->size--;
            for (int j = i; j <resources->size; j++){
                resources->stack[j] = resources->stack[j+1];
            }
            return 1;
        }
    }
    return 0;
}
void freeAllResources(struct resource_stack* resources){
    for (int i = 0; i < resources->size; i++){
        resources->stack[i].funcFree(resources->stack[i].item);
    }
}

#endif