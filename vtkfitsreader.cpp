#include "vtkfitsreader.h"
#include "vtkCommand.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkFloatArray.h"
#include <cmath>
#include "vtkPointData.h"
#include "vtkfitsreader.h"

#include <stdlib.h>     /* atof */
#include <string>
#include <vector>
#include <sstream>

#include <boost/algorithm/string.hpp>

#include <sys/time.h>  // Per utilizzare la funzione "gettimeofday" per prendere il tempo.

#include <omp.h>

//vtkCxxRevisionMacro(vtkFitsReader, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkFitsReader);

//----------------------------------------------------------------------------
vtkFitsReader::vtkFitsReader()
{
    this->filename[0]='\0';
    this->xStr[0]='\0';
    this->yStr[0]='\0';
    this->zStr[0]='\0';
    this->title[0]='\0';
    this->SetNumberOfInputPorts( 0 );
    this->SetNumberOfOutputPorts( 1 );

    for (int i=0; i<3; i++)
    {
        crval[i]=0;
        cpix[i]=0;
        cdelt[i]=0;
        naxes[i]= 10;
    }
    
    this->is3D=false;
    this->isMoment3D=false;

}

//----------------------------------------------------------------------------
vtkFitsReader::~vtkFitsReader()
{
}

void vtkFitsReader::SetFileName(std::string name) {


    if (name.empty()) {
        vtkErrorMacro(<<"Null Datafile!");
        return;
    }

    filename = name;
    this->Modified();
}


//----------------------------------------------------------------------------
vtkStructuredPoints* vtkFitsReader::GetOutput()
{
    return this->GetOutput(0);
}

//----------------------------------------------------------------------------
vtkStructuredPoints* vtkFitsReader::GetOutput(int port)
{
    return vtkStructuredPoints::SafeDownCast(this->GetOutputDataObject(port));
}



//----------------------------------------------------------------------------
int vtkFitsReader::FillOutputPortInformation(
        int vtkNotUsed(port), vtkInformation* info)
{
    // now add our info
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkStructuredPoints");
    return 1;
}




void vtkFitsReader::ReadHeader() {



    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */

    int status, nkeys, keypos, hdutype, ii, jj;
    char card[FLEN_CARD];   /* standard string lengths defined in fitsioc.h */
    
    
    char crval1[80];
    char crval2[80];
    char crval3[80];
    char crpix1[80];
    char crpix2[80];
    char crpix3[80];
    char cdelt1[80];
    char cdelt2[80];
    char cdelt3[80];
    char naxis1[80];
    char naxis2[80];
    char naxis3[80];
    
    
    crval1[0] ='\0';
    crval2[0] ='\0';
    crval3[0] ='\0';
    crpix1[0] ='\0';
    crpix2[0] ='\0';
    crpix3[0] ='\0';
    cdelt1[0] ='\0';
    cdelt2[0] ='\0';
    cdelt3[0] ='\0';
    
    std::string val1, val2, val3, pix1,pix2, pix3, delt1, delt2, delt3, nax1, nax2, nax3;

    status = 0;


    char *fn=new char[filename.length() + 1];;
    strcpy(fn, filename.c_str());

    if ( fits_open_file(&fptr, fn, READONLY, &status) )
        printerror( status );
    delete []fn;

    /* attempt to move to next HDU, until we get an EOF error */
    for (ii = 1; !(fits_movabs_hdu(fptr, ii, &hdutype, &status) ); ii++)
    {

        /* get no. of keywords */
        if (fits_get_hdrpos(fptr, &nkeys, &keypos, &status) )
            printerror( status );

        for (jj = 1; jj <= nkeys; jj++)  {

            if ( fits_read_record(fptr, jj, card, &status) )
                printerror( status );

            if (!strncmp(card, "CTYPE", 5)) {
                // cerr << card << endl;
                char *first = strchr(card, '\'');
                char *last = strrchr(card, '\'');

                *last = '\0';
                if (card[5] == '1')
                    strcpy(xStr, first+1);
                if (card[5] == '2')
                    strcpy(yStr, first+1);
                if (card[5] == '3')
                    strcpy(zStr, first+1);

            }

            if (!strncmp(card, "OBJECT", 6)) {
                cerr << card << endl;
                char *first = strchr(card, '\'');
                char *last = strrchr(card, '\'');
                *last = '\0';
                strcpy(title, first+1);
            }

            if (!strncmp(card, "CRVAL", 5)) {
                char *first = strchr(card, '=');
                char *last = strrchr(card, '=');
                *last = '\0';

                // char *last = strrchr(card, '/');
                //*last = '\0';

                if (card[5] == '1')
                {
                    strcpy(crval1, first+1);
                    char *pch = strtok (crval1," ,");
                    strcpy(crval1, pch);
                    
                }
                
                if (card[5] == '2')
                {
                    strcpy(crval2, first+1);
                    char *pch = strtok (crval2," ,");
                    strcpy(crval2, pch);

                }
                
                if (card[5] == '3')
                {
                    strcpy(crval3, first+1);
                    char *pch = strtok (crval3," ,");
                    strcpy(crval3, pch);

                }
            }

            if (!strncmp(card, "CRPIX", 5)) {
                char *first = strchr(card, '=');
                char *last = strrchr(card, '=');
                *last = '\0';
                
                
                if (card[5] == '1')
                {
                    strcpy(crpix1, first+1);

                    char *pch = strtok (crpix1," ,");
                    strcpy(crpix1, pch);
                }
                
                if (card[5] == '2')
                {
                    strcpy(crpix2, first+1);

                    char *pch = strtok (crpix2," ,");
                    strcpy(crpix2, pch);
                }
                if (card[5] == '3')
                {
                    strcpy(crpix3, first+1);

                    char *pch = strtok (crpix3," ,");
                    strcpy(crpix3, pch);
                }
            }

            if (!strncmp(card, "CDELT", 5)) {
                char *first = strchr(card, '=');
                char *last = strrchr(card, '=');
                *last = '\0';
                
                if (card[5] == '1')
                {
                    strcpy(cdelt1, first+1);
                    char *pch = strtok (cdelt1," ,");
                    strcpy(cdelt1, pch);
                    
                }
                
                if (card[5] == '2')
                {
                    strcpy(cdelt2, first+1);
                    char *pch = strtok (cdelt2," ,");
                    strcpy(cdelt2, pch);
                }
                
                if (card[5] == '3')
                {
                    strcpy(cdelt3, first+1);
                    char *pch = strtok (cdelt3," ,");
                    strcpy(cdelt3, pch);
                }
            }
            
            

        }
    }


    val1=crval1;
    val2=crval2;
    val3=crval3;
    pix1=crpix1;
    pix2=crpix2;
    pix3=crpix3;
    delt1=cdelt1;
    delt2=cdelt2;
    delt3=cdelt3;


    
    crval[0]= atof(val1.c_str());
    crval[1]= atof(val2.c_str());
    crval[2]= atof(val3.c_str());
    cpix[0]= atof(pix1.c_str());
    cpix[1]= atof(pix2.c_str());
    cpix[2]=atof(pix3.c_str());
    cdelt[0]= atof(delt1.c_str());
    cdelt[1]= atof(delt2.c_str());
    cdelt[2]= atof(delt3.c_str());
    initSlice=crval[2]-(cdelt[2]*(cpix[2]-1));


    
}

double vtkFitsReader::GetRMS() {
    return rms;
}

/* CUDA kernels */

//__global__ void min_max_Kernel ()


void vtkFitsReader::ReadDataAndCalculateRMS() {
    
    struct timeval start, start1, end, end1;
    double time_taken, time_taken1;
    
    gettimeofday(&start, NULL);
    
    ReadHeader();
    
    vtkStructuredPoints *output = (vtkStructuredPoints *) this->GetOutput();
    fitsfile *fptr;
    int status = 0, nfound = 0, anynull = 0;
    long fpixel, nbuffer, npixels, npixelstot, n=0; // npixelstot -> Modifica per la parallelizzazione
//    double meansquare=0;
    /* buffsize: Dimensione di una slice del cubo fits dal leggere. La lettura del cubo fits viene fatta slice per slice. */
//    const int buffsize = 500;
//    const int buffsize = 1000;
//    const int buffsize = 2000;
//    const int buffsize = 4000;
//    const int buffsize = 10000;
//    const int buffsize = 20000;
    const int buffsize = 40000; // (1)
//    const int buffsize = 100000;
//    const int buffsize = 200000;
//    const int buffsize = 400000;
//    const int buffsize = 1000000;
//    const int buffsize = 2000000;
//    const int buffsize = 4000000;
//    const int buffsize = 10420224; /* npixels in Small FITS */
//    const int buffsize = 10000000;
//    const int buffsize = 20000000;
//    const int buffsize = 40000000;
//    const int buffsize = 99081304; // (2)
//    const int buffsize = 100000000;
//    const int buffsize = 200000000;
//    const int buffsize = 396325216; //const int buffsize = naxes[0] * naxes[1] * naxes[2];
    /* npixels in Big FITS */
//    const int buffsize = 400000000;
//    const int buffsize = 792650432;
//    const int buffsize = 1000000000;
//    const int buffsize = 1981626080;
    
//    float nullval, buffer[buffsize];
    float nullval;
    float *buffer;
    buffer = new float[buffsize];
    char *fn=new char[filename.length() + 1];
    strcpy(fn, filename.c_str());
    
    if ( fits_open_file(&fptr, fn, READONLY, &status) )
        printerror( status );
    
    delete []fn;
    vtkFloatArray *scalars = vtkFloatArray::New();
    vtkFloatArray *scalars_test_1 = vtkFloatArray::New(); // VC
    vtkFloatArray *scalars_test_2 = vtkFloatArray::New(); // VC
    if ( fits_read_keys_lng(fptr, "NAXIS", 1, 3, naxes, &nfound, &status) )
        printerror( status );
    
    npixels  = naxes[0] * naxes[1] * naxes[2]; /* Numero totale di pixels del cubo fits dato dal prodotto del numero di pixels lungo ogni asse (naxes[0], naxes[1] e naxes[2]). */
    npixelstot = npixels;  // Modifica per la parallelizzazione
    cout<<"Occupied memory in bytes = "<<npixels*sizeof(float)<<endl;
    cout<<"Occupied memory in GB = "<<npixels*sizeof(float)/(1000*1000*1000)<<endl;
//    cout<<"Occupied memory in bytes = "<<10*npixels*sizeof(float)<<endl;
//    cout<<"Occupied memory in GB = "<<10*npixels*sizeof(float)/(1000*1000*1000)<<endl;
    n=npixels;
    
    fpixel   = 0; //fpixel   = 1; // Modifica per la parallelizzazione.
    nullval  = 0;
    datamin  = 1.0E30;
    datamax  = -1.0E30;

    output->SetDimensions(naxes[0], naxes[1], naxes[2]);
    output->SetOrigin(0.0, 0.0, 0.0);
    
    scalars->Allocate(npixels);
    
    int slice;
    int num=0;

    minmaxslice=new float*[naxes[2]];
    for(int i=0;i< naxes[2];i++)
    {
        minmaxslice[i] = new float[2];
        minmaxslice[i][0]= 1.0E30;
        minmaxslice[i][1]= -1.0E30;
    }

    int itncounter = 0;
    
    gettimeofday(&end, NULL);
        time_taken = (end.tv_sec - start.tv_sec) * 1e6;
        time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
    
    cout<<"Time taken before the while loop' = "<<time_taken<<" s"<<endl;
    
    gettimeofday(&start1, NULL);
    
    //For every pixel
    cout<<"npixels = "<<npixels<<endl;
//    cout<<"npixels%buffsize = "<<npixels%buffsize<<endl;
//    cout<<"itntotal_before = "<<npixels/buffsize + 1<<endl;
//    cout<<"naxes[0]*naxes[1] = "<<naxes[0]*naxes[1]<<endl;
    
    long itntotal = 0;
    
    if (npixels%buffsize != 0) {
        itntotal = npixels/buffsize + 1;
    } else {
        itntotal = npixels/buffsize;
    }
    
//    cout<<"itntotal_after = "<<itntotal<<endl;
//    cout<<"itntotal/2 = "<<itntotal/2<<endl;
    
    double *bad_local;   // Modifica per la parallelizzazione
    bad_local = new double[itntotal];   // Modifica per la parallelizzazione
    int bad=0;
    
    double *meansquare_local;   // Modifica per la parallelizzazione
    meansquare_local = new double[itntotal];   // Modifica per la parallelizzazione
    double meansquare=0;
    

    cout<<"The code runs on "<<omp_get_num_threads( )<<" OpenMP threads before the parallel region!"<<endl;
    
#pragma omp parallel private(npixels, nbuffer, fpixel, buffer, num, slice, datamin, datamax) shared(bad_local, meansquare_local)
{
    cout<<"The code runs on "<<omp_get_num_threads()<<" OpenMP threads inside the parallel region!"<<endl;
    
    /* Inizializzazione delle variabili "private" su ogni thread di OpenMP */
    npixels = npixelstot;  // In ogni thread di OpenMP, npixels viene inizializzato al numero totale di pixels nel cubo FITS, npixelstot.
    nbuffer = 0;
    fpixel = 0;
//    for (int i = 0; i < buffsize; i++) buffer[i] = 0.0;
    num = 0;
    slice = 0;
    datamin  = 1.0E30;
    datamax  = -1.0E30;
    
    #pragma omp for reduction(+: bad) reduction(+: meansquare)
    for (long itncounter = 0; itncounter < itntotal; itncounter++) {
        
        gettimeofday(&start, NULL);

        nbuffer = npixels;
        if (npixels > buffsize)
            nbuffer = buffsize;
        
        fpixel = itncounter*nbuffer + 1;
        
        if ( fits_read_img(fptr, TFLOAT, fpixel, nbuffer, &nullval,
                           buffer, &anynull, &status) )
            printerror( status );
        
        num = itncounter*nbuffer;
        for (long ii = 0; ii < nbuffer; ii++)  {
            slice = ((num + ii)/ (naxes[0]*naxes[1]) );

            if (std::isnan(buffer[ii]))
                buffer[ii] = -1000000.0;
            
            scalars->InsertValue(itncounter*buffsize + ii,buffer[ii]);

            if ( buffer[ii]!=-1000000.0)
            {
                if ( buffer[ii] < datamin )
                    datamin = buffer[ii];
                if ( buffer[ii] > datamax   )
                    datamax = buffer[ii];

                if ( buffer[ii] < minmaxslice[slice][0] )
                    minmaxslice[slice][0] = buffer[ii];
                if ( buffer[ii] > minmaxslice[slice][1]   )
                    minmaxslice[slice][1] = buffer[ii];

                // meansquare+=buffer[ii]*buffer[ii];
                // media+=buffer[ii];
                meansquare_local[itncounter]+=buffer[ii]*buffer[ii];

            }
            else
                bad_local[itncounter] = bad_local[itncounter] + 1;   // Modifica per la parallelizzazione
        }
        
        bad = bad + bad_local[itncounter];
        meansquare = meansquare + meansquare_local[itncounter];
        
        npixels = npixelstot - (itncounter + 1)*nbuffer;
//        npixels -= nbuffer;
//        fpixel  += nbuffer;
        
        gettimeofday(&end, NULL);
            time_taken = (end.tv_sec - start.tv_sec) * 1e6;
            time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;

        cout<<"Time taken by iteration "<<itncounter<<" of the while loop' = "<<time_taken<<" s"<<endl;
    }
}
    gettimeofday(&end1, NULL);
        time_taken1 = (end1.tv_sec - start1.tv_sec) * 1e6;
        time_taken1 = (time_taken1 + (end1.tv_usec - start1.tv_usec)) * 1e-6;
    
    cout<<"Time taken by the while loop' = "<<time_taken1<<" s"<<endl;

    gettimeofday(&start, NULL);
    
    n=n-bad;
    double means=meansquare/n;
    rms=sqrt(means);
    // sigma=qSqrt(sigma/n);

    if ( fits_close_file(fptr, &status) )
        printerror( status );
    
    output->GetPointData()->SetScalars(scalars);
    
    gettimeofday(&end, NULL);
        time_taken = (end.tv_sec - start.tv_sec) * 1e6;
        time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
    
    cout<<"Time taken after the while loop' = "<<time_taken<<" s"<<endl;

    return;
}

// Note: from cookbook.c in fitsio distribution.
void vtkFitsReader::printerror(int status) {

    cerr << "vtkFitsReader ERROR.";
    if (status) {
        fits_report_error(stderr, status); /* print error report */
        exit( status );    /* terminate the program, returning error status */
    }
    return;
}


int vtkFitsReader::GetNaxes(int i)
{

    return naxes[i];

}

