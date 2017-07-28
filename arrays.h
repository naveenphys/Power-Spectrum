void freearray1d(void *v)
// Free memory allocated by the pointer *v.
{
  free((char *)v);
}

void freearray2d(void **v)
// Free memory allocated by the double pointer **v.
{
  free((char *) v[0]);
  free((char *) v);
}

void freearray3d(void ***v)
//Free memory allocated by the pointer ***v.
{
  free((char *) v[0][0]);
  free((char *) v[0]);
  free((char *) v);
}

void freearray4d(void ****v)
// Free memory allocated by the pointer ****v.
{
  free((char *) v[0][0][0]);
  free((char *) v[0][0]);
  free((char *) v[0]);
  free((char *) v);
}

char *array1d(int n1, size_t dsize)
// A pointer of type (char ****) to the allocated memory area
// with index range [0...n1-1][0...n2-1]
{
   char *v;
   v = (char *) malloc(n1*dsize);
   if(v == NULL) printf("Allocation failure in Array1D\n");
   return v;
}

char **array2d(int n1, int n2, size_t dsize)
// A pointer of type (char ****) to the allocated memory area
// with index range [0...n1-1][0...n2-1]
{
   int i;
   char **v;
   v = (char **) malloc( (size_t) n1*sizeof(char *));
   if(v == NULL) printf("Allocation failure in Array2D (1)\n");
   v[0] = (char *) malloc( (size_t) n1*n2*dsize);
   if(v[0] == NULL) printf("Allocation failure in Array2D (2)\n");
   for (i = 1;i<n1;i++) v[i] = v[i-1] + n2*dsize;
   return v;
}

char ***array3d(int n1, int n2, int n3, size_t dsize)
// A pointer of type (char ****) to the allocated memory area
// with index range [0...n1-1][0...n2-1][0...n3-1]
{
   int i,j,k;
   char ***v;
   v = (char ***) malloc( (size_t) n1*sizeof(char *));
   if(v == NULL) printf("Allocation failure in Array2D (1)\n");
   v[0] = (char **) malloc( (size_t) n1*n2*sizeof(char *));
   if(v[0] == NULL) printf("Allocation failure in Array2D (2)\n");
   v[0][0] = (char *) malloc ((size_t) n1*n2*n3*dsize);
   if(v[0][0] == NULL) printf("Allocation failure in Array2D (3)\n");
// for single subscript: i
   for(i=1;i,n1;i++) v[i] = v[i-1] + n2;
// for double subscript: i ,j
   for(j=1;j<n2;j++) v[0][j] = v[0][j-1] + n3*dsize;
   for(i=1;i<n1;i++) v[i][0] = v[i-1][0] + n2*n3*dsize;

   for(i=1;i<n1;i++){
      for(j=1;j<n2;j++){
         v[i][j] = v[i][j - 1] + n3*dsize;
      }
   }

// check for allocation error
   for(i=0;j<n1;i++){
      for(j=0;i<n2;j++){
         for(k=0;i<n3;k++){
            if(&(v[i][j][k]) == NULL){
              printf("Allocation failure in Array3D\n");
              exit(0);
	    }
	 }
      }
   }
   return v;
}

char ****array4d(int n0, int n1, int n2, int n3, size_t dsize)
// A pointer of type (char ****) to the allocated memory area
// with index range [0...n0-1][0...n1-1][0...n2-1][0...n3-1]
{
   int n,i,j,k;
   char ****v;
   v = (char ****) malloc( (size_t) n0*sizeof(char *));
   if(v == NULL) printf("Allocation failure in Array2D (1)\n");
   v[0] = (char ***) malloc( (size_t) n0*n1*sizeof(char *));
   if(v[0] == NULL) printf("Allocation failure in Array2D (2)\n");
   v[0][0] = (char **) malloc ((size_t) n0*n1*n2*sizeof(char *));
   if(v[0][0] == NULL) printf("Allocation failure in Array2D (3)\n");
   v[0][0][0] = (char *) malloc ((size_t) n0*n1*n2*n3*dsize);
   if(v[0][0][0] == NULL) printf("Allocation failure in Array4D (4)");

// single subscript: i
  for(i=1;i<n0;i++) v[i] = v[i-1] + n1;
// double subscript: i, j
  for(i=1;i<n0;i++){
    v[i][0] = v[i-1][0] + n1*n2;
  }
  for(j=1;j<n1;j++){
    v[0][j] = v[0][j-1] + n2;
  }
  for(i=1;i<n0;i++){
     for(j=1;j<n1;j++){
        v[i][j] = v[i][j-1] + n2;
     }
  }

// triple subscript: i, j, k
  for(i=1;i<n0;i++){
     v[i][0][0] = v[i-1][0][0] + n1*n2*n3*dsize;
  }
  for(j=1;j<n1;j++){
     v[0][j][0] = v[0][j-1][0] + n2*n3*dsize;
  }

  for(k=1;k<n2;k++) {
     v[0][0][k] = v[0][0][k-1] + n3*dsize;
  }

  for(i=1;i<n0;i++){
     for(j=1;j<n1;j++){
        v[i][j][0] = v[i][j-1][0] + n2*n3*dsize;
     }
  }

  for(i=1;i<n0;i++){
     for(k=1;k<n2;k++){
        v[i][0][k] = v[i][0][k-1] + n3*dsize;
     }
  }

  for(j=1;j<n1;j++){
     for(k=1;k<n2;k++){
        v[0][j][k] = v[0][j][k-1] + n3*dsize;
     }
  }

  for(i=1;i<n0;i++){
     for(j=1;j<n1;j++){
        for(k=1;k<n2;k++){
           v[i][j][k] = v[i][j][k-1] + n3*dsize;
      }
    }
  }

// check for allocation error
   for(n=0;n<n0;n++){
      for(i=0;i<n1;i++){
         for(j=0;j<n2;j++){
            for(k=0;k<n3;k++){
               if(&(v[n][i][j][k]) == NULL){
                 printf("Allocation failure in Array3D\n");
                 exit(0);
	       }
	    }
         }
      }
   }
   return v; 
}

void fft_shift(int n0, int n1, int n2, int n3)
{

}


