// .NAME vtkFitsReader - read structured points from FITS file.
// .SECTION Description
// vtkFitsReader is a source object that reads FITS data files
// .SECTION Caveats
// Uses CFITSIO v2.0 (http://heasarc.gsfc.nasa.gov/docs/software/fitsio)

#ifndef __vtkFitsReader_h
#define __vtkFitsReader_h


#include "vtkAlgorithm.h"
#include "vtkStructuredPoints.h"
#include "vtkFloatArray.h"

extern "C" {
#include "fitsio.h"
}


//class VTK_EXPORT vtkFitsReader : public vtkStructuredPointsSource
class VTK_EXPORT vtkFitsReader : public vtkAlgorithm
{
public:
    //    vtkFitsReader();
    //  static vtkFitsReader *New() {return new vtkFitsReader;}

    static vtkFitsReader *New();
    vtkFitsReader();
    ~vtkFitsReader();

    void SetFileName(std::string name);
    std::string GetFileName(){return filename;}

    vtkTypeMacro(vtkFitsReader,vtkAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);
    void PrintHeader(ostream& os, vtkIndent indent);
    void PrintTrailer(std::ostream& os , vtkIndent indent);

    
    bool is3D;
    bool isMoment3D;
    void PrintDetails();

    float GetMin(){return datamin;}
    float GetMax(){return datamax;}
    
    float GetCrval(int i){
        if(i<3)
            return crval[i];
        else return -1;
    }
    float GetCpix(int i){
        if(i<3) return cpix[i];
        else return -1;}
    float GetCdelt(int i){if(i<3) return cdelt[i];
        else return -1;}

    int GetNaxes(int i);
    

    // Description:
    // Get the output data object for a port on this algorithm.
    vtkStructuredPoints* GetOutput();
    vtkStructuredPoints* GetOutput(int);
    virtual void SetOutput(vtkDataObject* d);
    // Description:
    // see vtkAlgorithm for details
    virtual int ProcessRequest(vtkInformation*,
                               vtkInformationVector**,
                               vtkInformationVector*);

    int getMomentOrder();

    void setMomentOrder(int order);

protected:
   
    float datamin;
    float datamax;
    float epoch;   // Part of FITS Header file
    int *dimensions;   // [x,y,z]
    int naxis;
    int points;
    double crval[3];
    double cpix[3];
    double cdelt[3];
    double sigma=0;
    double media=0;
    double rms=0;
    long npix=0;
    long naxes[3];

    float **minmaxslice;

    // Description:
    // This is called by the superclass.
    // This is the method you should override.
    virtual int RequestDataObject(
            vtkInformation* request,
            vtkInformationVector** inputVector,
            vtkInformationVector* outputVector );

    // convenience method
    virtual int RequestInformation(
            vtkInformation* request,
            vtkInformationVector** inputVector,
            vtkInformationVector* outputVector );

    // Description:
    // This is called by the superclass.
    // This is the method you should override.
    virtual int RequestData(
            vtkInformation* request,
            vtkInformationVector** inputVector,
            vtkInformationVector* outputVector );

    // Description:
    // This is called by the superclass.
    // This is the method you should override.
    virtual int RequestUpdateExtent(
            vtkInformation*,
            vtkInformationVector**,
            vtkInformationVector* );

    virtual int FillOutputPortInformation( int port, vtkInformation* info );

private:
    //  vtkFitsReader( const vtkFitsReader& ); // Not implemented.
    //  void operator = ( const vtkFitsReader& );  // Not implemented.
    std::string filename;
    char title[80];
    char xStr[80];
    char yStr[80];
    char zStr[80];
//    char cunit3[1][200];
    double initSlice;

    void ReadHeader();
    void printerror(int status); // from fitsio distribution
};



#endif
