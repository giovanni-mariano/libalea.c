// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_mcnp.h"
#include "alea_openmc.h"
#include "alea_xml.h"
#include "alea_fingerprint.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
static alea_fingerprint_set_t* load(const char* path,const alea_fingerprint_options_t* options) {
    alea_fingerprint_set_t* f=NULL;
    size_t n=strlen(path);
    if(n>=9 && !strcmp(path+n-9,".alea.xml")) {
        alea_model_t* m=alea_xml_load(path);
        if(m) {
            f=alea_fingerprint_set_create(alea_model_system(m),options);
            alea_model_destroy(m);
        }
    }
    else if(n>=4 && !strcmp(path+n-4,".xml")) {
        openmc_model_t* m=openmc_load(path);
        if(m) {
            f=alea_fingerprint_set_create(m->sys,options);
            openmc_model_destroy(m);
        }
    }
    else {
        mcnp_model_t* m=mcnp_load(path);
        if(m) {
            f=alea_fingerprint_set_create(m->sys,options);
            mcnp_model_destroy(m);
        }
    }
    return f;
}
int main(int argc,char** argv) {
    int json=0;
    alea_fingerprint_options_t options={sizeof(options),ALEA_FINGERPRINT_LEGACY};
    const char* paths[2]= {
        0
    };
    int count=0;
    for(int i=1;i<argc;i++) {
        if(!strcmp(argv[i],"--json"))json=1;
        else if(!strcmp(argv[i],"--profile")) {
            if(++i==argc)return 2;
            if(!strcmp(argv[i],"semantic"))options.profile=ALEA_FINGERPRINT_SEMANTIC;
            else if(!strcmp(argv[i],"legacy"))options.profile=ALEA_FINGERPRINT_LEGACY;
            else {fprintf(stderr,"Unknown comparison profile: %s\n",argv[i]);return 2;}
        }
        else if(argv[i][0]=='-')return 2;
        else if(count<2)paths[count++]=argv[i];
        else return 2;
    }
    if(!count) {
        fprintf(stderr,"Usage: mc_compare [--json] [--profile legacy|semantic] model [other-model]\n");
        return 2;
    }
    alea_fingerprint_set_t* sa=load(paths[0],&options);
    alea_fingerprint_set_t* sb=count==2?load(paths[1],&options):NULL;
    const alea_fingerprints_t* a=alea_fingerprint_set_values(sa);
    const alea_fingerprints_t* b=alea_fingerprint_set_values(sb);
    int rc=2;
    if(!a || (count==2 && !b)) {
        if(json)puts("{\"error\":\"cannot fingerprint input\"}");
        else fprintf(stderr,"%s\n",alea_error());
        goto done;
    }
    if(!b) {
        if(json)printf("{\"version\":%u,\"profile\":\"%s\",\"fingerprints\":[",a->version,options.profile==ALEA_FINGERPRINT_LEGACY?"legacy":"semantic");
        size_t total=a->cell_count+a->universe_count;
        for(size_t i=0;i<total;i++) {
            int u=i>=a->cell_count;
            const alea_fingerprint_t* f=u?&a->universes[i-a->cell_count]:&a->cells[i];
            char g[65],h[65];
            alea_digest_hex(&f->geometry,g);
            alea_digest_hex(&f->content,h);
            if(json)printf("%s{\"kind\":\"%s\",\"id\":%d,\"geometry\":\"%s\",\"content\":\"%s\"}",i?",":"",u?"universe":"cell",f->id,g,h);
            else printf("%s %d geometry=%s content=%s\n",u?"universe":"cell",f->id,g,h);
        }
        if(json)puts("]}");
        rc=0;
    }
    else {
        alea_model_diff_t* diff=alea_fingerprint_set_diff(sa,sb);
        if(!diff)goto done;
        rc=0;
        static const char* names[]= {
            "unchanged","geometry","content","renumbered","removed","added","ambiguous"
        };
        if(json)printf("{\"version\":%u,\"profile\":\"%s\",\"changes\":[",a->version,options.profile==ALEA_FINGERPRINT_LEGACY?"legacy":"semantic");
        for(size_t i=0;i<diff->count;i++) {
            alea_model_change_t* c=&diff->changes[i];
            if(c->kind!=ALEA_DIFF_UNCHANGED)rc=1;
            if(json)printf("%s{\"kind\":\"%s\",\"old_id\":%d,\"new_id\":%d,\"change\":\"%s\"}",i?",":"",c->universe?"universe":"cell",c->old_id,c->new_id,names[c->kind]);
            else printf("%s %d -> %d: %s\n",c->universe?"universe":"cell",c->old_id,c->new_id,names[c->kind]);
        }
        if(json)puts("]}");
        alea_model_diff_destroy(diff);
    }
    done:alea_fingerprint_set_destroy(sa);
    alea_fingerprint_set_destroy(sb);
    return rc;
}
