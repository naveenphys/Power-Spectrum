void read_dbl(int f, double ****vel)
{
   FILE *fp;
   double d;
   long int offset;
   char filenumb[5];
   char filename[100];
   int i,j,k,dir;

   sprintf(filenumb,"%04d",f);
   strcat(filename,datdir);
   strcat(filename,prefix);
   strcat(filename,filenumb);
   strcat(filename,suffix);

   printf("%s\n",filename);
   fp = fopen(filename,"rb");
   if(fp==NULL){
     printf("fopen error (1)");
     exit(0);
   }

   for(i=0;i<nx;i++){
      for(j=0;j<ny;j++){
         for(k=0;k<nz;k++){
            for(dir=0;dir<nv;dir++){
	       offset = i*ny*nz+j*ny+k + (dir+1)*nx*ny*nz;
	       fseek(fp,sizeof(double)*offset,SEEK_SET);
	       fread(&d,sizeof(double),1,fp);
	       vel[dir][i][j][k] = d;
	    }
         }
      }
   }
   fclose(fp);
   return;
}


void printarray4d(int n0, int n1, int n2, int n3, double ****vel)
{
   int n,i,j,k;
   for(n=0;n<n0;n++){
       for(i=0;i<n1;i++){
          for(j=0;j<n2;j++){
             for(k=0;k<n3;k++){
                printf("%d,%2d,%2d,%2d,%16.12f\n",n,i,j,k,vel[n][i][j][k]);
             }
          }
       }
    }
   return;
}

void write_velk(int dir, double ****velk, fftw_complex *out)
{
   int i,j,k;
   double real, imag;
   for(i=0;i<nx;i++){
      for(j=0;j<ny;j++){
         for(k=0;k<nz/2+1;k++){
            real = creal(out[i*nx*ny+j*ny+k]);
	    imag = cimag(out[i*nx*ny+j*ny+k]);
            velk[dir][i][j][k]=sqrt(real*real+imag*imag);
         }
      }
   }
   return;
}

void write_output(int type, int f, int n0, int n1, int n2, int n3, double ****velk)
{
   FILE *fp;
   double d;
   long int offset;
   char filenumb[5];
   char filename[100];

   sprintf(filenumb,"%04d",f);
   strcat(filename,datdir);
   strcat(filename,prefix);
   strcat(filename,filenumb);
   
   if(type == 0){
     strcat(filename,".dbl");
   } 
   if(type == 1){
     strcat(filename,".vtk");
   }

   printf("%s\n",filename);

   int n,i,j,k;
   if(type == 0){
     fp = fopen(filename,"wb");
     for(n=0;n<n0;n++){
	for(i=0;i<n1;i++){
           for(j=0;j<n2;j++){
	      for(k=0;k<n3;k++){
		 fwrite(velk[n][i][j][k],sizeof(double),1,fp);
	      }
	   }
	}
     }
     fclose(fp);
   }


   if(type == 1){
     fp = fopen(filename,"r");


     fclose(fp);
     for(n=0;n<n0;n++){
	fp = fopen(filename,"r");

	for(i=0;i<n1;i++){
           for(j=0;j<n2;j++){
	      for(k=0;k<n3;k++){
		 fwrite(velk[n][i][j][k],sizeof(double),1,fp);
	      }
	   }
	}
     }
     fclose(fp);
   }
}










}
