// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "cluster_internal.h"

#include <ctype.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "util/alea_mcnp_input.h"

alea_cluster_status_t alea_cluster_read_mcnp_input(
        alea_cluster_t* cluster,const char* path,char** data,size_t* length){
    if(!cluster)return ALEA_CLUSTER_INVALID_ARGUMENT;
    char* expanded=NULL;size_t size=0;
    alea_cluster_status_t local=ALEA_CLUSTER_OK;
    if(cluster->rank==0){
        int rc=path?alea_read_mcnp_input(path,NULL,&expanded,&size,NULL):-1;
        if(rc==-1)local=ALEA_CLUSTER_IO_ERROR;
        else if(rc==-2)local=ALEA_CLUSTER_OUT_OF_MEMORY;
        else if(rc)local=ALEA_CLUSTER_COMPUTE_ERROR;
    }
    return alea_cluster_broadcast_owned_bytes(cluster,expanded,size,local,data,length);
}
