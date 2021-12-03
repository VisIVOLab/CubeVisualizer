// Derived from VTK/Examples/Cxx/Medical4.cxx
// This example reads a volume dataset and displays it via volume rendering.
//

#include <vtkCamera.h>
#include <vtkColorTransferFunction.h>
#include <vtkFixedPointVolumeRayCastMapper.h>
#include <vtkMetaImageReader.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPiecewiseFunction.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkVolume.h>
#include <vtkVolumeProperty.h>
#include "vtkfitsreader.h"
#include <vtkOutlineFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkMarchingCubes.h>
#include <vtkProperty.h>
#include <vtkAxesActor.h>
#include <vtkTextProperty.h>
#include <vtkCaptionActor2D.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkTextActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkLegendScaleActor.h>
#include <vtkConeSource.h>
#include  <vtkFlyingEdges3D.h>
#include <vtkStripper.h>
#include <vtkPlanes.h>
#include <array>
#include <vtkFrustumSource.h>
#include <vtkAlgorithmOutput.h>

//SetRequestedRenderModeToGPU
int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    cout << "Usage: " << argv[0] << "cube.fits" << endl;
    return EXIT_FAILURE;
  }

  vtkNew<vtkNamedColors> colors;

  std::array<unsigned char, 4> bkg{{51, 77, 102, 255}};
  colors->SetColor("BkgColor", bkg.data());

  // Create the renderer, the render window, and the interactor. The renderer
  // draws into the render window, the interactor enables mouse- and
  // keyboard-based interaction with the scene.
  vtkNew<vtkRenderer> ren;
  vtkNew<vtkRenderWindow> renWin;
  renWin->AddRenderer(ren);
  vtkNew<vtkRenderWindowInteractor> iren;
  iren->SetRenderWindow(renWin);
    

  auto reader = vtkSmartPointer<vtkFitsReader>::New();
  reader->SetFileName(argv[1]);
  reader->ReadDataAndCalculateRMS();

  vtkStructuredPoints *out=reader->GetOutput();

    // outline
    vtkNew<vtkOutlineFilter> outlineF;
    outlineF->SetInputData(out);
    vtkNew<vtkPolyDataMapper> outlineM;
    outlineM->SetInputConnection(outlineF->GetOutputPort());
    outlineM->ScalarVisibilityOff();
    vtkNew<vtkActor> outlineA;
    outlineA->SetMapper(outlineM);
            
    // isosurface
    auto shellE = vtkFlyingEdges3D::New();
    shellE->SetInputData(reader->GetOutput());
    shellE->ComputeNormalsOn();
    shellE->SetValue(0, 3*reader->GetRMS());
    vtkPolyDataMapper *shellM = vtkPolyDataMapper::New();
    shellM->SetInputConnection(shellE->GetOutputPort());
    shellM->ScalarVisibilityOff();
    vtkActor *shellA = vtkActor::New();
    shellA->SetMapper(shellM);
    shellA->GetProperty()->SetColor(1.0, 0.5, 1.0);

    ren->AddActor(shellA);
    ren->AddActor(outlineA);
    
    // Set a background color for the renderer
    ren->SetBackground(0.21,0.23,0.25);


  // Increase the size of the render window
  renWin->SetSize(1280 , 960);
  renWin->SetWindowName("FitsCubeDemo");

  // Interact with the data.
  renWin->Render();
  iren->Initialize();
    

  iren->Start();

  return EXIT_SUCCESS;
}

