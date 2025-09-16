#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include "hdf5.h"

/*
 * This is a simple C program to convert the HR file, which stores the data of the Wannier tight-binding models, to
 * HDF5 Format.
 *
 * How to use
 * =============================
 *  This program takes two arguments:
 *  1) the name of the HR file
 *  2) the name of the HDF5 file output
 *
 *
 * Author : Nabil
 */

#define MAX_BUFF_SIZE 1024

void strtod_arr(char *s, double *vals) {
    int counter = 0;
    for (const char *tok = strtok(s, " ");
                tok != nullptr; tok = strtok(nullptr, " \t\n\r")) {
        vals[counter++] = strtod(tok, nullptr);
                }
}


int read_int_from_nextline_fstream(FILE* f) {
    // This just reads a line containing an integer.
    // Need to do some error checking. Otherwise, it will unsafe to use this function, as it will
    // return 0 when a read error is encountered


    char buf[MAX_BUFF_SIZE] = "";

    if (fgets(buf, MAX_BUFF_SIZE, f) == NULL) {
        fprintf(stderr, "Error reading an integer from the data file. Please double check the hr file. \n");
        return -1;
    }

    errno = 0; char *endptr;
    const long numwann_tol = strtol(buf, &endptr, 10);

    if (errno != 0) {
        fprintf(stderr, "Error encountered in strtol.. \n");
        return -1;
    }

    if (endptr == buf) {
        fprintf(stderr, "Invalid input at line: %s \n", buf);
    }

    // Invalid input
    return (int) numwann_tol;

}

int main(const int argc, char *argv[]) {


    // Check if the number of arguments is exactly two
    if (argc != 3) {
        printf("Usage: <Wannier hr file path>  <HDF5 file path>\n");
        return 1;
    }
    // first argument
    const char *hrfile   = argv[1];
    fprintf(stdout, "--- Name of the file containing the TB model data: %s \n", hrfile);
    // second argument
    const char *hdf5file = argv[2];
    fprintf(stdout, "--- Name of the HDF5 file                        : %s \n", hdf5file);
    fprintf(stdout, "\n");


    // Attempt to open the file
    FILE *hrfp = fopen(hrfile, "r");

    clock_t before = clock();

    // NOTE: If the fopen fails to open the file (it does not exist), then errno will store the error number
    //          and hrfp will be NULL
    if (hrfp == NULL) {
        fprintf(stderr, "Could not open file: %s\n", hrfile);
        return 1;
    }
    fprintf(stdout, "--- successfully opened the hr data file\n");
    fprintf(stdout, "\n");


    // Skip the first line, which is a comment line.
    int c = 0;
    do { c = fgetc(hrfp); } while (c != '\n' && c != EOF);

    // The next line is the number of Wannier orbitals (number of bands in the TB model).
    const int num_wann = read_int_from_nextline_fstream(hrfp);
    fprintf(stdout, "--- Number of wannier orbitals: %i \n", num_wann);

    // Number of R Vectors
    const int nrvecs = read_int_from_nextline_fstream(hrfp);
    fprintf(stdout, "--- Number of R vectors       : %i \n", nrvecs);

    // Next lines record the degeneracy factor for each R vectors.
    /*
     * This is organized as several lines where each line contains no more than 15 integers. The total number of
     * lines that I need to read is nrvecs / 15 + 1. We read these integers into a dynamically allocated large array
     */
    const int num_lines_degeneracy_facts    = (nrvecs / 15) + 1;
    int* degeneracy_facts                   = malloc(nrvecs * sizeof(int));
    int  current_line                       = 0;
    char buffer[MAX_BUFF_SIZE]; int counter = 0;

    while (fgets(buffer, MAX_BUFF_SIZE, hrfp)) {

        // Make sure we are not going beyond the last line containing degeneracy information
        if (current_line >= num_lines_degeneracy_facts) {
            break;
        }

        // Now, we need to get the integers written in the current buffer.
        char *endptr = nullptr;
        for (const char *tok = strtok(buffer, " ");
                    tok != nullptr; tok = strtok(nullptr, " \t\n\r")) {
                    degeneracy_facts[counter++] = (int)strtol(tok, &endptr, 10);
                    }
        current_line++;
    }

    // Now, we can read the hopping amplitudes //
    const unsigned long n_doubles_ham = num_wann * num_wann * nrvecs;
    const unsigned long n_bytes_ham   = n_doubles_ham * sizeof(double);

    fprintf(stdout, "--- Size of the model data in memory: %lu bytes \n", 2 * n_bytes_ham);

    // Allocate Arrays //
    double* rvecs = malloc(nrvecs * 3 * sizeof(double));     // Row major: Address using rvecs[index * 3 + direction]

    double* reH = malloc( n_bytes_ham );    // Real part of the Hopping amplitudes
    double* imH = malloc( n_bytes_ham );    // Imaginary part of the Hopping amplitudes


    if (!reH || !imH || !rvecs) {
        fprintf(stderr, "Memory Allocation Failure: Cannot allocate arrays for the real space model.\n");
        free(degeneracy_facts);
        free(rvecs);
        free(reH);
        free(imH);
        fclose(hrfp);
        return 1;
    }

    // first data entry //
    double vals[7];
    strtod_arr(buffer, vals);

    rvecs[0]   = vals[0]; rvecs[1]   = vals[1]; rvecs[2]   = vals[2];   // R Vectors

    // First Hopping amplitude [orbitals 1, 1]
    reH[0] = vals[5] / degeneracy_facts[0]; imH[0] = vals[6] / degeneracy_facts[0];

    unsigned int rvecs_counter3 = 0;
    unsigned int rvecs_counter = 0;
    while (fgets(buffer, MAX_BUFF_SIZE, hrfp)) {

        strtod_arr(buffer, vals);

        /*
         * Data stored in each row:
         * 0 => n1  [R = n1 * A1 + n2 * A2 + n3 * A3]
         * 1 => n2
         * 2 => n3
         * 3 => Alpha
         * 4 => Beta
         * 5 => Re H
         * 6 => Im H
         */

        // Update the R vectors //
        if (rvecs[rvecs_counter3] != vals[0]
            || rvecs[rvecs_counter3 + 1] != vals[1]
            || rvecs[rvecs_counter3 + 2] != vals[2]) {

            rvecs_counter3 += 3;
            rvecs[rvecs_counter3]     = vals[0];
            rvecs[rvecs_counter3 + 1] = vals[1];
            rvecs[rvecs_counter3 + 2] = vals[2];

            rvecs_counter++;

        }
        const int alpha = (int)vals[3] - 1;
        const int beta  = (int)vals[4] - 1;

        reH[rvecs_counter + num_wann * nrvecs * alpha + nrvecs * beta] = vals[5] / degeneracy_facts[rvecs_counter];
        imH[rvecs_counter + num_wann * nrvecs * alpha + nrvecs * beta] = vals[6] / degeneracy_facts[rvecs_counter];

    }


    const clock_t elapsed_time_reading_input = clock() - before;
    fprintf(stdout, "--- Time taken to read the Hr file  : %f seconds\n",
        (double)elapsed_time_reading_input / CLOCKS_PER_SEC);

    // Put Data in HDF5 Format //
    fprintf(stdout, "\n");
    fprintf(stdout, "--- Creating HDF5 data file ---\n");

    const hid_t fi = H5Fcreate(hdf5file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (fi < 0) {
        fprintf(stderr, "Could not open HDF5 file for writing: %s\n", hdf5file);
    }

    const hsize_t dims[1]       = { n_doubles_ham };
    const hsize_t dims_rvecs[1] = { nrvecs * 3 };

    // Re H
    const hid_t dspace_id_reH = H5Screate_simple(1, dims, nullptr);
    const hid_t dset_id_reH   = H5Dcreate2(fi, "/reH", H5T_NATIVE_DOUBLE, dspace_id_reH,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const hid_t status_reH    = H5Dwrite(dset_id_reH, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, reH);

    // Im H
    const hid_t dspace_id_imH = H5Screate_simple(1, dims, nullptr);
    const hid_t dset_id_imH   = H5Dcreate2(fi, "/imH", H5T_NATIVE_DOUBLE, dspace_id_imH,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const hid_t status_imH    = H5Dwrite(dset_id_imH, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, imH);

    // R vectors
    const hid_t dspace_rvecs = H5Screate_simple(1, dims_rvecs, nullptr);
    const hid_t dset_id_rvecs= H5Dcreate2(fi, "/rvecs", H5T_NATIVE_DOUBLE, dspace_rvecs,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const hid_t status_rvecs = H5Dwrite(dset_id_rvecs, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, rvecs);

    // Number of Wannier functions
    const hid_t dspace_numwann = H5Screate(H5S_SCALAR);
    const hid_t dset_numwann   = H5Dcreate2(fi, "/nw", H5T_NATIVE_INT,
                        dspace_numwann, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const hid_t status_nwann   = H5Dwrite(dset_numwann, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                  H5P_DEFAULT, &num_wann);

    const hid_t dspace_nr = H5Screate(H5S_SCALAR);
    const hid_t dset_nr   = H5Dcreate2(fi, "/nr", H5T_NATIVE_INT,
                        dspace_nr, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const hid_t status_nr = H5Dwrite(dset_nr, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                  H5P_DEFAULT, &nrvecs);




    if (status_imH < 0 || status_reH < 0 || status_rvecs < 0 || status_nwann < 0 || status_nr < 0) {
        fprintf(stderr, "Could not write the model array or metadata to HDF5 file: %s\n", hdf5file);
    }
    else {
        fprintf(stdout, "--- Done --- \n");
    }




    // Free allocated memory and close the input data file //
    free(degeneracy_facts);
    free(reH);
    free(imH);
    free(rvecs);
    fclose(hrfp);
    H5Fclose(fi);
    return 0;
}