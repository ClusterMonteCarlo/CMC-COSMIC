/*
* Copyright 1993-2006 NVIDIA Corporation.  All rights reserved.
*
* NOTICE TO USER:   
*
* This source code is subject to NVIDIA ownership rights under U.S. and 
* international Copyright laws.  
*
* NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE 
* CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR 
* IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH 
* REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF 
* MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.   
* IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL, 
* OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS 
* OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE 
* OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE 
* OR PERFORMANCE OF THIS SOURCE CODE.  
*
* U.S. Government End Users.  This source code is a "commercial item" as 
* that term is defined at 48 C.F.R. 2.101 (OCT 1995), consisting  of 
* "commercial computer software" and "commercial computer software 
* documentation" as such terms are used in 48 C.F.R. 12.212 (SEPT 1995) 
* and is provided to the U.S. Government only as a commercial end item.  
* Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through 
* 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the 
* source code with only those rights set forth herein.
*/


/* CUda UTility Library :: additional functionality for graphics
 *                         interoperability */

// includes, file
#include <cutil_interop.h>

// includes, system
#include <stdio.h>
#include <stdio.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

// namespace unnamed (internal)
namespace 
{  
#ifndef _WIN32
    // definitions for lib nvidia-cfg
    typedef enum {
        NVCFG_TRUE = 1,
        NVCFG_FALSE = 0,
    } NvCfgBool;

    typedef struct {
        int bus;
        int slot;
    } NvCfgDevice;
#endif

} // end, namespace unnamed
//////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
//! Tests if graphics interop is supported on the current system
//! @return  CUTTrue if graphics interop is supported
//////////////////////////////////////////////////////////////////////////////
CUTBoolean CUTIL_API
isInteropSupported() {
#ifdef _WIN32
    int num_devices = 0;
    int dev = 0;

    do {
        DISPLAY_DEVICE dispDevice;
        memset(&dispDevice, 0, sizeof(dispDevice));
        dispDevice.cb = sizeof(dispDevice);
        if (!EnumDisplayDevices(NULL, dev, &dispDevice, 0))
            break;
        if ((dispDevice.StateFlags & DISPLAY_DEVICE_ATTACHED_TO_DESKTOP) ||
            (!(dispDevice.StateFlags & DISPLAY_DEVICE_MIRRORING_DRIVER) &&
            !(dispDevice.StateFlags & DISPLAY_DEVICE_MODESPRUNED)))
            ++num_devices;
        ++dev;
    } while (dev);

#else // Linux

    void* handle = 0;

    // NvCfgBool (*ptr_nvCfgGetDevices)( int*, NvCfgDevice**);
    typedef NvCfgBool (*nvCfgGetDevicesType)( int*, NvCfgDevice**);
    nvCfgGetDevicesType ptr_nvCfgGetDevices;

    // acquire handle to shared library
    handle = dlopen( "libnvidia-cfg.so", RTLD_LAZY);
    if( ! handle) {
        // no NVIDIA driver installed
        fprintf( stderr, "Cannot find NVIDIA driver.\n" );
        fprintf( stderr, "Graphics interoperability not supported.\n" );
        return CUTFalse;
    }

    ptr_nvCfgGetDevices = (nvCfgGetDevicesType) 
                                        dlsym( handle, "nvCfgGetDevices");
    char* error = 0;
    if( 0 != (error = dlerror())) {
        fprintf( stderr, "Cannot query number of devices in the system.\n" );
        fprintf( stderr, "Graphics interoperability not supported.\n" );
        return CUTFalse;
     }

     // call function
     int num_devices;
     NvCfgDevice* devs;
   
     (*ptr_nvCfgGetDevices)( &num_devices, &devs);

#endif // _WIN32
    if( 1 != num_devices) {
        fprintf( stderr, "Graphics interoperability on multi GPU systems currently not supported.\n" );
        return CUTFalse;

    }
    return CUTTrue;
}


