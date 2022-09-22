/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2015-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "atom.h"
#include "mpi.h"
#include "math.h"
#include "modify.h"
#include "properties.h"
#include "error.h"
#include "memory.h"
//#include "fix_multisphere.h"
//#include "multisphere.h"
#include "fix_property_particle.h"
#include "fix_property_global.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   find a property that was requested
   called e.g. from CFD data exchange model
   property may come from atom class, from a fix property, or fix ms
   last 2 args are the data length and are used for all data
------------------------------------------------------------------------- */

void* Properties::find_property_pu(const char *name, const char *type, int &len1, int &len2)
{
        
    void *ptr = NULL;

    // possiblility 1
    // may be atom property - look up in atom class

    ptr = atom->extract(name,len2);
    // if nlocal == 0 and property found atom->extract returns NULL and len2 >= 0
    if(ptr || len2 >= 0)
    {
        len1 = atom->nlocal; //atom->tag_max();
        // check if length correct
        if(((strcmp(type,"scalar-atom") == 0) && (len2 != 1)) || ((strcmp(type,"vector-atom") == 0) && (len2 != 3)))
            return NULL;
        return ptr;
    }

    // possiblility 2
    // may come from a fix multisphere
    // also handles scalar-multisphere and vector-multisphere

    /*ms_ = static_cast<FixMultisphere*>(modify->find_fix_style_strict("multisphere",0));
    if(ms_) ms_data_ = &ms_->data();

    if(ms_)
    {
        ptr = ms_->extract(name,len1,len2);
        if(((strcmp(type,"scalar-multisphere") == 0) && (len2 != 1)) || ((strcmp(type,"vector-multisphere") == 0) && (len2 != 3)))
            return NULL;

        if(ptr || ((len1 >= 0) && (len2 >= 0)))
        return ptr;
    }*/

    // possiblility 3
    // may be fix property per atom - look up in modify class

    Fix *fix = NULL;
    FixPropertyParticle* fix2;
    
    if(strcmp(type,"scalar-atom") == 0)
    {
       fix = modify->find_fix_property(name,"property/particle","scalar",0,0,"cfd coupling",false);
       if(fix)
       {
           len1 = atom->nlocal; //atom->tag_max();
           len2 = 1;
           return (void*) static_cast<FixPropertyParticle*>(fix)->vector_atom;
       }
    }
    else if(strcmp(type,"vector-atom") == 0)
    {
       fix2 = static_cast<FixPropertyParticle*>(modify->find_fix_property(name,"property/particle","vector",0,0,"cfd coupling",false));
       if(fix2)
       {
           len1 = atom->nlocal; //atom->tag_max();
           //len2 = 3;
	   len2 = fix2->nvalues;
           return (void*) fix2->array_atom;
       }
    }
    else if(strcmp(type,"scalar-global") == 0)
    {
       fix = modify->find_fix_property(name,"property/global","scalar",0,0,"cfd coupling",false);
       len1 = len2 = 1;
       if(fix) return (void*) static_cast<FixPropertyGlobal*>(fix)->values;
    }
    else if(strcmp(type,"vector-global") == 0)
    {
       fix = modify->find_fix_property(name,"property/global","vector",0,0,"cfd coupling",false);
       if(fix)
       {
           len1 = static_cast<FixPropertyGlobal*>(fix)->nvalues;
           len2 = 1;
           return (void*) static_cast<FixPropertyGlobal*>(fix)->values;
       }
    }
    else if(strcmp(type,"matrix-global") == 0)
    {
       fix = modify->find_fix_property(name,"property/global","matrix",0,0,"cfd coupling",false);
       if(fix)
       {
           len1  = static_cast<FixPropertyGlobal*>(fix)->size_array_rows;
           len2  = static_cast<FixPropertyGlobal*>(fix)->size_array_cols;
           return (void*) static_cast<FixPropertyGlobal*>(fix)->array;
       }
    }
    else if(strcmp(name,"ex") == 0) 
    {
        // possiblility 4A - Dipole is specified as atom property (requires DIPOLE package)
        ptr = atom->extract("mu",len2);
        if(ptr)
        {
            len1 = atom->nlocal; //atom->tag_max();
            // check if length correct
            if( (strcmp(type,"vector-atom") == 0) && (len2 != 3) )
                return NULL;
            return ptr;
        }

        // possiblility 4B - Quaternion is specified as atom property (requires ASPHERE package)
        // requires a fix that computes orientation data from quaternion,
        // or another fix that provides orientationEx (e.g., from POEMS)
        //TODO: Write fix that computes orientation from quaternion
        //TODO: Check if correct data is drawn in case a FixPOEMS is used
        fix = modify->find_fix_property("orientationEx","property/particle","vector",0,0,"cfd coupling",false);
        if(fix)
        {
               len1 = atom->nlocal; //atom->tag_max();
               len2 = 3;
               return (void*) static_cast<FixPropertyParticle*>(fix)->array_atom;
        }
        else printf("WARNING: Fix with name 'orientationEx' not found that stores orientation information. \n");
    }
    return NULL;
}
