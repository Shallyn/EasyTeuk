

#include "myFileIO.h"
#include <stdio.h>
#include <unistd.h>
#include <string.h>

INT get_REAL8TimeSeries(REAL8TimeSeries ** Series, FILE * fp)
{
    const size_t block = 1024;
    REAL8 start;
    REAL8 end;
    REAL8 dt;
    REAL8 *data = NULL;
    size_t bufsz = 0;
    size_t n, l;
    CHAR line[LINE_MAX];
    CHAR t0[LINE_MAX];
    CHAR t1[LINE_MAX];
    
    for (l = 0, n = 0; fgets(line, sizeof(line), fp); ++l) {
        int c;
        if (*line == '#')
            continue;
        if (n == bufsz) {       /* allocate more memory */
            bufsz += block;
            data = realloc(data, bufsz * sizeof(*data));
        }
        c = sscanf(line, "%le %le", n ? &end : &start, data + n);
        if (c != 2) {
            print_err("error: format error on line %zd: %s\n", l, line);
            exit(1);
        }
        ++n;
    }

    data = realloc(data, n * sizeof(*data));
    dt = (end - start) / (n - 1);
    *Series = CreateREAL8TimeSeries (start, (REAL8)dt, n);
    memcpy((*Series)->data->data, data, n * sizeof(*data));
    
    free(data);
    return CEV_SUCCESS;
}



INT get_REAL8TimeSeries_waveform(REAL8TimeSeries ** hplus, REAL8TimeSeries ** hcross, FILE * fp)
{
    const size_t block = 1024;
    REAL8 start;
    REAL8 end;
    REAL8 dt;
    REAL8 *hp = NULL;
    REAL8 *hc = NULL;
    size_t bufsz = 0;
    size_t n, l;
    CHAR line[LINE_MAX];
    
    for (l = 0, n = 0; fgets(line, sizeof(line), fp); ++l) {
        int c;
        if (*line == '#')
            continue;
        if (n == bufsz) {       /* allocate more memory */
            bufsz += block;
            hp = realloc(hp, bufsz * sizeof(*hp));
            hc = realloc(hc, bufsz * sizeof(*hc));
        }
        c = sscanf(line, "%le %le %le", n ? &start : &end, hp + n, hc + n);
        if (c != 3) {
            print_err("error: format error on line %zd: %s\n", l, line);
            exit(1);
        }
        ++n;
    }
    hp = realloc(hp, n * sizeof(*hp));
    hc = realloc(hc, n * sizeof(*hp));
    dt = (end - start) / (n - 1);
    *hplus = CreateREAL8TimeSeries (start, dt, n);
    *hcross = CreateREAL8TimeSeries (start, dt, n);
    memcpy((*hplus)->data->data, hp, n * sizeof(*hp));
    memcpy((*hcross)->data->data, hc, n * sizeof(*hc));
    
    free(hp);
    free(hc);
    return CEV_SUCCESS;
}

INT get_COMPLEX16TimeSeries_waveform(COMPLEX16TimeSeries ** hout, FILE * fp)
{
    const size_t block = 1024;
    REAL8 start;
    REAL8 end;
    REAL8 dt;
    REAL8 *hp = NULL;
    REAL8 *hc = NULL;
    size_t bufsz = 0;
    size_t n, l;
    CHAR line[LINE_MAX];
    COMPLEX16TimeSeries *hLM;
    
    for (l = 0, n = 0; fgets(line, sizeof(line), fp); ++l) {
        int c;
        if (*line == '#')
            continue;
        if (n == bufsz) {       /* allocate more memory */
            bufsz += block;
            hp = realloc(hp, bufsz * sizeof(*hp));
            hc = realloc(hc, bufsz * sizeof(*hc));
        }
        c = sscanf(line, "%le %le %le", n ? &start : &end, hp + n, hc + n);
        if (c != 3) {
            print_err("error: format error on line %zd: %s\n", l, line);
            exit(1);
        }
        ++n;
    }
    hp = realloc(hp, n * sizeof(*hp));
    hc = realloc(hc, n * sizeof(*hp));
    dt = (end - start) / (n - 1);
    hLM = CreateCOMPLEX16TimeSeries(0, dt, n);
    for (l=0; l < n; l++)
    {
        hLM->data->data[l] = hp[l] + I*hc[l];
    }    
    free(hp);
    free(hc);
    *hout = hLM;
    return CEV_SUCCESS;
}

// void CreateFolder(CHAR *folderName)
// {
//     // CHAR folderName[] = "RunData";

//     if (access(folderName, 0) == -1)
//     {
//         _mkdir(folderName);
//     }
// }

INT cmd_mkdir(CHAR *folderName)
{
    CHAR cmd[256] = "mkdir ";
    strcat(cmd, folderName);
    if (access(folderName, 0) == -1)
    {
        system(cmd);
    }
    return CEV_SUCCESS;
}

INT read_waveform(REAL8Vector **time, 
                  REAL8Vector **hreal, 
                  REAL8Vector **himag,
                  FILE *file)
{
    const size_t block = 1024;
    REAL8 start;
    REAL8 end;
    double *t = NULL;
    double *hp = NULL;
    double *hc = NULL;
    size_t bufsz = 0;
    size_t n, l;
    char line[LINE_MAX];

    for (l = 0, n = 0; fgets(line, sizeof(line), file); ++l) 
    {
        int c;
        if (*line == '#')
            continue;
        if (n == bufsz) 
        {       
            /* allocate more memory */
            bufsz += block;
            t = realloc(t, bufsz * sizeof(*t));
            hp = realloc(hp, bufsz * sizeof(*hp));
            hc = realloc(hc, bufsz * sizeof(*hc));
        }
        c = sscanf(line, "%le %le %le", t + n, hp + n, hc + n);
        if (c != 3) {
            fprintf(stderr, "error: format error on line %zd: %s\n", l, line);
            free(t);
            free(hp);
            free(hc);
            return CEV_FAILURE;
        }
        ++n;
    }
    t = realloc(t, n * sizeof(*t));
    hp = realloc(hp, n * sizeof(*hp));
    hc = realloc(hc, n * sizeof(*hp));
    REAL8Vector *tVec = CreateREAL8Vector(n);
    REAL8Vector *hrVec = CreateREAL8Vector(n);
    REAL8Vector *hiVec = CreateREAL8Vector(n);
    memcpy(tVec->data, t, n*sizeof(*t));
    memcpy(hrVec->data, hp, n * sizeof(*hp));
    memcpy(hiVec->data, hc, n * sizeof(*hc));
    free(t);
    free(hp);
    free(hc);
    *time = tVec;
    *hreal = hrVec;
    *himag = hiVec;
    return CEV_SUCCESS;
}


/**
* Create a hdf5 file
*/
INT DumpREAL8ArrayTohdf5(CHAR *fname, CHAR *dname, REAL8Array *array, INT is_delete)
{
    hid_t       file, group, subgroup, gcpl;        /* Handles */
    hid_t       datatype, dataspace, dataset;
    herr_t      status;
    H5G_info_t  ginfo;
    
    if (!array)
        return CEV_FAILURE;
    // for(i=0; i<3; i++)
    //     for(j=0; j<4; j++)
    //         for(k=0; k<5; k++)
    //             data[i][j][k] = i+j+k;
    /*
     * Create a new file using the default properties.
     */
    if (access(fname, 0) == -1)
    {
        file = H5Fcreate (fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    } else if (is_delete) {
        remove(fname);
        file = H5Fcreate (fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    } else {
        file = H5Fopen (fname, H5F_ACC_RDWR, H5P_DEFAULT);
    }

    /*
     * Describe the size of the array and create the data space.
     */
    // dims[0] = 3;
    // dims[1] = 4;
    // dims[2] = 5;
    // int         dimLen;
    // dimLen = array->dimLength->length;
    hsize_t     dims[array->dimLength->length];       /* Index */
    int i;
    for (i=0; i<array->dimLength->length; i++)
        dims[i] = array->dimLength->data[i];
    dataspace = H5Screate_simple(array->dimLength->length, dims, NULL);

    /*
     * Define datatype for the data in the file
     */
    datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    status = H5Tset_order(datatype, H5T_ORDER_LE);

    /*
     * Create a new dataset within the file using defined dataspace
     */
    dataset = H5Dcreate2(file, dname, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Write the data into the dataset using default transfer properties
     */
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, array->data);

#if 0
    /*
     * Create group creation property list and enable link creation
     * order tracking.  Attempting to track by creation order in a
     * group that does not have this property set will result in an
     * error.
     */
    gcpl = H5Pcreate (H5P_GROUP_CREATE);
    status = H5Pset_link_creation_order( gcpl, H5P_CRT_ORDER_TRACKED |
                H5P_CRT_ORDER_INDEXED );

    /*
     * Create primary group using the property list.
     */
    group = H5Gcreate (file, gname, H5P_DEFAULT, gcpl, H5P_DEFAULT);
	dataset = H5Dopen(group, "Signal", H5P_DEFAULT);

    /*
     * Create subgroups in the primary group.  These will be tracked
     * by creation order.  Note that these groups do not have to have
     * the creation order tracking property set.
     */
    subgroup = H5Gcreate (group, "H", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Gclose (subgroup);

    /*
     * Get group info.
     */
    status = H5Gget_info (group, &ginfo);

    /*
     * Traverse links in the primary group using alphabetical indices
     * (H5_INDEX_NAME).
     */
    printf("Traversing group using alphabetical indices:\n\n");
    for (i=0; i<ginfo.nlinks; i++) {

        /*
         * Get size of name, add 1 for null terminator.
         */
        size = 1 + H5Lget_name_by_idx (group, ".", H5_INDEX_NAME, H5_ITER_INC,
                    i, NULL, 0, H5P_DEFAULT);

        /*
         * Allocate storage for name.
         */
        name = (char *) malloc (size);

        /*
         * Retrieve name, print it, and free the previously allocated
         * space.
         */
        size = H5Lget_name_by_idx (group, ".", H5_INDEX_NAME, H5_ITER_INC, i, name,
                    (size_t) size, H5P_DEFAULT);
        printf ("Index %d: %s\n", (int) i, name);
        free (name);
    }

    /*
     * Traverse links in the primary group by creation order
     * (H5_INDEX_CRT_ORDER).
     */
    printf("\nTraversing group using creation order indices:\n\n");
    for (i=0; i<ginfo.nlinks; i++) {

        /*
         * Get size of name, add 1 for null terminator.
         */
        size = 1 + H5Lget_name_by_idx (group, ".", H5_INDEX_CRT_ORDER,
                    H5_ITER_INC, i, NULL, 0, H5P_DEFAULT);

        /*
         * Allocate storage for name.
         */
        name = (char *) malloc (size);

        /*
         * Retrieve name, print it, and free the previously allocated
         * space.
         */
        size = H5Lget_name_by_idx (group, ".", H5_INDEX_CRT_ORDER, H5_ITER_INC, i,
                    name, (size_t) size, H5P_DEFAULT);
        printf ("Index %d: %s\n", (int) i, name);
        free (name);
    }
    // status = H5Pclose (gcpl);
    // status = H5Gclose (group);
#endif
    /*
     * Close and release resources.
     */
    status = H5Sclose(dataspace);
    status = H5Tclose(datatype);
    status = H5Dclose(dataset);
    status = H5Fclose (file);

    return CEV_SUCCESS;
}
