#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <petsc.h>

PetscErrorCode WriteMat(Mat mat,char const *name){

    PetscViewer viewer;
    PetscErrorCode ierr;

    char filename[64] = "Mat_"; char pfix[12] = ".m";
    strcat(filename,name); strcat(filename,pfix);
    ierr = PetscObjectSetName((PetscObject)mat,name); CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);
    ierr = MatView(mat,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}



#endif // FUNCTIONS_H
