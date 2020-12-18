## Requirements

- Eigen
- FFTW

Define environment variables EIGEN_INC, FFTW_INC and FFTW_LIB (see Makefile)

## Example

In folder test/
    
    ../ffth --input IM64.raw  --nx 64 --ny 64 --nz 64 --loading load.inp --material mat.inp  --phase phases.inp

## Input files

- `IM64.raw`

    Unsigned char (uint8) : phases from 0 to n (max. 255)

- `load.inp`

    ```
    E11 E22 E33 E12 E23 E13
    ```

-  `mat.inp`

    ```
    nmat
    
    C^0_ij [6x6]
    a^0_i [1x6]

    C^1_ij [6x6]
    a^1_i [1x6]

    ...
    
    ```

    Indices are contracted using a modified Voigt notation 11 --> 1, 22 --> 2, 33 --> 3, 12 --> 4, 23 --> 5, 13 --> 6

- `phases.inp`

    ```
    nphase
    matnum theta_x theta_y theta_z
    ...

    ```

    where theta_x, theta_y ,theta_z are rotation angles around x, y, z axis. (cf. `rotate_matrix` in `src/utils.h`)


