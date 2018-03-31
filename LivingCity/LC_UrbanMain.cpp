//---------------------------------------------------------------------------------------------------------------------
// Copyright 2017, 2018 Purdue University, Ignacio Garcia Dorado, Daniel Aliaga
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
// following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the
// following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the
// following disclaimer in the documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote
// products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS AS IS AND ANY EXPRESS OR IMPLIED WARRANTIES,
// INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
// USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//---------------------------------------------------------------------------------------------------------------------

#define BOOST_TT_HAS_OPERATOR_HPP_INCLUDED

#include "LC_UrbanMain.h"
#include "LC_GLWidget3D.h"
#include "qthread.h"
#include "qfiledialog.h"

#include "roadGraphDynameqLoader.h"

namespace LC {

LCUrbanMain::LCUrbanMain(QWidget *parent,
                         Qt::WindowFlags flags) : QMainWindow(parent, flags) {
  ui.setupUi(this);

  glWidget3D = new LCGLWidget3D(this);
  this->setCentralWidget(glWidget3D);
  initModeActive = true; //for proper initialization
  init();
  initModeActive = false;
}//

LCUrbanMain::~LCUrbanMain() {
}//


/**
      * Initialize geometry and generate the entire geometry objects.
      */
void LCUrbanMain::init() {
  G::global()["figureMode"] = 0;
  G::global()["mouseState"] = M_NORMAL;
  G::global()["selectedBoundingPolygonVertexIdx"] = 0;

  G::global()["render_lighting"] = true;
  G::global()["view_place_types_centept_render"] = true;
  G::global()["view_buildings_render"] = true;
  G::global()["view_arterials_edges_render"] = true;
  G::global()["view_arterials_vertices_render"] = true;
  G::global()["view_arterials_arrows_render"] = true;
  G::global()["view_trees_render"] = true;

  G::global()["roadLaneWidth"] = 4.5f;
  G::global()["render_lighting"] = true;

  //------ Procedural modeling

  connect(ui.proceduralModelingCheckBox, SIGNAL(stateChanged(int)),
          this, SLOT(onProceduralModeling(int)));
  //------ Place types ------
  connect(ui.numPlaceTypesSpinBox, SIGNAL(valueChanged(int)),
          this, SLOT(numPlaceTypesSlot(int)));

  //------ Roads ------
  //connect(ui.ArterialEdgesWidthSlider, SIGNAL(valueChanged(int)),
  //	this, SLOT(arterialEdgesWidthSlot(int)));
  connect(ui.cudaLanesNumArterialSpinBox, SIGNAL(valueChanged(int)),
          this, SLOT(cudaLanesNumArterialSlot(int)));
  connect(ui.cudaLanesNumRoadSpinBox, SIGNAL(valueChanged(int)),
          this, SLOT(cudaLanesNumRoadSlot(int)));


  connect(ui.arterialEdgesIrregularitySlider, SIGNAL(valueChanged(int)),
          this, SLOT(arterialEdgesIrregularitySlot(int)));
  connect(ui.arterialEdgesCurvatureSlider, SIGNAL(valueChanged(int)),
          this, SLOT(arterialEdgesCurvatureSlot(int)));
  connect(ui.arterialNumDerpartingDial, SIGNAL(valueChanged(int)),
          this, SLOT(arterialNumDepartingSlot(int)));

  connect(ui.localEdgesLengthUSlider, SIGNAL(valueChanged(int)),
          this, SLOT(localEdgesLengthUSlot(int)));
  connect(ui.localEdgesLengthVSlider, SIGNAL(valueChanged(int)),
          this, SLOT(localEdgesLengthVSlot(int)));
  connect(ui.localEdgesIrregularitySlider, SIGNAL(valueChanged(int)),
          this, SLOT(localEdgesIrregularitySlot(int)));
  connect(ui.localEdgesCurvatureSlider, SIGNAL(valueChanged(int)),
          this, SLOT(localEdgesCurvatureSlot(int)));

  // ------ Parcels ------
  connect(ui.parcelAreaMeanSlider, SIGNAL(valueChanged(int)),
          this, SLOT(parcelAreaMeanSlot(int)));
  connect(ui.parcelAreaDeviationSlider, SIGNAL(valueChanged(int)),
          this, SLOT(parcelAreaDeviationSlot(int)));
  connect(ui.parcelSplitDeviationSlider, SIGNAL(valueChanged(int)),
          this, SLOT(parcelSplitDeviationSlot(int)));
  connect(ui.parcelSetbackFrontSlider, SIGNAL(valueChanged(int)),
          this, SLOT(parcelSetbackFrontSlot(int)));
  connect(ui.parcelSetbackSidesSlider, SIGNAL(valueChanged(int)),
          this, SLOT(parcelSetbackSidesSlot(int)));


  // ------ Buildings ------
  connect(ui.buildingHeightMeanSlider, SIGNAL(valueChanged(int)),
          this, SLOT(buildingHeightMeanSlot(int)));
  connect(ui.buildingHeightDeviationSlider, SIGNAL(valueChanged(int)),
          this, SLOT(buildingHeightDeviationSlot(int)));

  //--LU
  connect(ui.landUseParkPercentSlider, SIGNAL(valueChanged(int)),
          this, SLOT(landUseParkPercentageSlot(int)));

  // GUI
  //connect(ui.action2D, SIGNAL(triggered()), this, SLOT(onRender2DSlot()));
  //connect(ui.action3D, SIGNAL(triggered()), this, SLOT(onRender3DSlot()));

  // TERRAIN
  connect(ui.terrainPaint_sizeSlider, SIGNAL(valueChanged(int)), this,
          SLOT(updateTerrainLabels(int)));
  connect(ui.terrainPaint_changeSlider, SIGNAL(valueChanged(int)), this,
          SLOT(updateTerrainLabels(int)));

  connect(ui.actionNewTerrain, SIGNAL(triggered()), this, SLOT(onNewTerrain()));
  connect(ui.actionLoadTerrain, SIGNAL(triggered()), this, SLOT(onLoadTerrain()));
  connect(ui.actionSaveTerrain, SIGNAL(triggered()), this, SLOT(onSaveTerrain()));

  updateTerrainLabels(-1);

  // LAYERS
  connect(ui.layer_enableCheckBox, SIGNAL(stateChanged(int)), this,
          SLOT(onLayerEnable(int)));
  connect(ui.layer_numberSpinBox, SIGNAL(valueChanged(int)), this,
          SLOT(onLayerEnable(int)));

  connect(ui.layer_loadLayersButton, SIGNAL(clicked(bool)), this,
          SLOT(onLoadLayers(bool)));
  connect(ui.layer_saveLayersButton, SIGNAL(clicked(bool)), this,
          SLOT(onSaveLayers(bool)));
  connect(ui.layer_clearLayersButton, SIGNAL(clicked(bool)), this,
          SLOT(onClearLayers(bool)));

  // TRAFFIC
  connect(ui.cudaCPUSimulateButton, SIGNAL(clicked(bool)), this,
          SLOT(onSimulatePressed(bool)));
  connect(ui.ArterialEdgesSpeedSlider, SIGNAL(valueChanged(int)), this,
          SLOT(onArterialEdgesSpeedSlot(int)));
  ////////////////////////////////////////////////////////
  // TOOL BAR
  ui.topDock->setTitleBarWidget(new QWidget());

  connect(ui.TB_but2DButton, SIGNAL(clicked(bool)), this,
          SLOT(onRender2DSlot(bool)));
  connect(ui.TB_but3DButton, SIGNAL(clicked(bool)), this,
          SLOT(onRender3DSlot(bool)));
  connect(ui.TB_roadLabelButton, SIGNAL(clicked(bool)), this,
          SLOT(onLabelPressed(bool)));

  // B TRAFFIC
  connect(ui.bCPUSimulateButton, SIGNAL(clicked(bool)), this,
          SLOT(onBSimulateCPUPressed(bool)));
  connect(ui.bDynameqDemandCheckBox, SIGNAL(stateChanged(int)), this,
          SLOT(onDynameqDemand(int)));

  connect(ui.bCreatePeopleButton, SIGNAL(clicked(bool)), this,
          SLOT(onCreatePeople(bool)));
  connect(ui.bLoadPeopleButton, SIGNAL(clicked(bool)), this,
          SLOT(onLoadPeople(bool)));
  connect(ui.bSavePeopleButton, SIGNAL(clicked(bool)), this,
          SLOT(onSavePeople(bool)));
  //GPU
  connect(ui.bGPUSimulateButton, SIGNAL(clicked(bool)), this,
          SLOT(onBSimulateGPUPressed(bool)));
  ui.progressBar->hide();

  // Hide menu
  ui.placeTypeEditWidget->hide();
}//

/////////////////////////////////
// PM
/////////////////////////////////

void LCUrbanMain::numPlaceTypesSlot(int val) {


  // Make sure that we have enought place types, otherwise append copied ones
  if (!initModeActive) {
    int num_place_types = val;
    G::global()["num_place_types"] = num_place_types;

    while (glWidget3D->cg.geoPlaceTypes.myPlaceTypes.size() <= num_place_types) {
      glWidget3D->cg.geoPlaceTypes.myPlaceTypes.push_back(
        glWidget3D->cg.geoPlaceTypes.myPlaceTypes.back());
    }

    //glWidget3D->generateGeometry(ClientGeometry::kStartFromRoads);
    glWidget3D->generateGeometry(ClientGeometry::kJustRoads);
    glWidget3D->updateGL();
  }
}//

void LCUrbanMain::arterialEdgesIrregularitySlot(int val) {

  if (!initModeActive) {
    int selected = glWidget3D->cg.selectedPlaceTypeIdx;
    glWidget3D->cg.geoPlaceTypes[selected]["pt_edges_irregularity"] = 0.01f *
        (float)val;;
    //glWidget3D->generateGeometry(ClientGeometry::kStartFromRoads);
    glWidget3D->generateGeometry(ClientGeometry::kJustRoads);
    glWidget3D->updateGL();
  }
}//

void LCUrbanMain::arterialEdgesCurvatureSlot(int val) {
  if (!initModeActive) {
    int selected = glWidget3D->cg.selectedPlaceTypeIdx;
    glWidget3D->cg.geoPlaceTypes[selected]["pt_edges_curvature"] = 0.01f *
        (float)val;;
    //glWidget3D->generateGeometry(ClientGeometry::kStartFromRoads);
    glWidget3D->generateGeometry(ClientGeometry::kJustRoads);
    glWidget3D->updateGL();
  }
}//

void LCUrbanMain::arterialNumDepartingSlot(int val) {

  ui.arterialNumDerpartingGroupBox->setTitle("#Rad: " + QString::number(val));

  if (!initModeActive) {
    int selected = glWidget3D->cg.selectedPlaceTypeIdx;
    glWidget3D->cg.geoPlaceTypes[selected]["pt_num_departing"] = val;
    //glWidget3D->generateGeometry(ClientGeometry::kStartFromRoads);
    glWidget3D->generateGeometry(ClientGeometry::kJustRoads);
    glWidget3D->updateGL();
  }
}

void LCUrbanMain::localEdgesLengthUSlot(int val) {
  if (!initModeActive) {

    int selected = glWidget3D->cg.selectedPlaceTypeIdx;
    glWidget3D->cg.geoPlaceTypes[selected]["pt_loc_edges_lengthU"] = 0.01f *
        (float)val;
    //glWidget3D->generateGeometry(ClientGeometry::kStartFromRoads);
    glWidget3D->generateGeometry(ClientGeometry::kJustRoads);
    glWidget3D->updateGL();

  }
}//

void LCUrbanMain::localEdgesLengthVSlot(int val) {

  if (!initModeActive) {
    int selected = glWidget3D->cg.selectedPlaceTypeIdx;
    glWidget3D->cg.geoPlaceTypes[selected]["pt_loc_edges_lengthV"] = 0.01f *
        (float)val;
    //glWidget3D->generateGeometry(ClientGeometry::kStartFromRoads);
    glWidget3D->generateGeometry(ClientGeometry::kJustRoads);
    glWidget3D->updateGL();
  }
}//

void LCUrbanMain::localEdgesIrregularitySlot(int val) {
  if (!initModeActive) {
    int selected = glWidget3D->cg.selectedPlaceTypeIdx;
    glWidget3D->cg.geoPlaceTypes[selected]["pt_loc_edges_irregularity"] = 0.01f *
        (float)val;
    //glWidget3D->generateGeometry(ClientGeometry::kStartFromRoads);
    glWidget3D->generateGeometry(ClientGeometry::kJustRoads);
    glWidget3D->updateGL();
  }
}

void LCUrbanMain::localEdgesCurvatureSlot(int val) {
  if (!initModeActive) {
    int selected = glWidget3D->cg.selectedPlaceTypeIdx;
    glWidget3D->cg.geoPlaceTypes[selected]["pt_loc_edges_curvature"] = 0.01f *
        (float)val;
    //glWidget3D->generateGeometry(ClientGeometry::kStartFromRoads);
    glWidget3D->generateGeometry(ClientGeometry::kJustRoads);
    glWidget3D->updateGL();
  }
}//



/*void LCUrbanMain::arterialEdgesWidthSlot(int val)
      {
      LC::misctools::Global::global()->arterial_edges_width =
      (float)val*0.5f;
      ui.ArterialEdgesWidthGroupBox->setTitle(
      "Width: "
      + QString::number(LC::misctools::Global::global()->arterial_edges_width)
      + " m"
      );

      if(!initModeActive){
      cg.generateGeometry(ClientGeometry::kStartFromRoads,clientMain->mGLWidget_3D);
      mGLWidget_3D->updateGL();
      }
      }*/

void LCUrbanMain::cudaLanesNumArterialSlot(int val) {
  G::global()["cuda_arterial_numLanes"] = val;
  G::global()["arterial_edges_width"] = (2 * val + 1) *
                                        G::global().getFloat("roadLaneWidth");// REMV

  if (!initModeActive) {
    //glWidget3D->generateGeometry(ClientGeometry::kStartFromRoads);
    glWidget3D->generateGeometry(ClientGeometry::kJustRoads);
    glWidget3D->updateGL();
  }
}//

void LCUrbanMain::cudaLanesNumRoadSlot(int val) {
  G::global()["cuda_road_numLanes"] = val;

  if (!initModeActive) {
    //glWidget3D->generateGeometry(ClientGeometry::kStartFromRoads);
    glWidget3D->generateGeometry(ClientGeometry::kJustRoads);
    glWidget3D->updateGL();
  }
}//

void LCUrbanMain::parcelAreaMeanSlot(int val) {
  float parcel_area_mean = ui.parcelAreaMeanSlider->value() *
                           ui.parcelAreaMeanSlider->value() * 5.0f;
  float parcel_area_deviation = ui.parcelAreaDeviationSlider->value();
  ui.parcelAreaGroupBox->setTitle(
    "Area: "
    + QString::number(parcel_area_mean)
    + " m2 +- "
    + QString::number(parcel_area_deviation)
    + " %"
  );

  if (!initModeActive) {
    int selected = glWidget3D->cg.selectedPlaceTypeIdx;
    glWidget3D->cg.geoPlaceTypes[selected]["pt_parcel_area_mean"] =
      parcel_area_mean;
    //glWidget3D->generateGeometry(ClientGeometry::kStartFromRoads);
    glWidget3D->generateGeometry(ClientGeometry::kJustRoads);
    glWidget3D->updateGL();
  }
}//

void LCUrbanMain::parcelAreaDeviationSlot(int val) {
  float parcel_area_mean = ui.parcelAreaMeanSlider->value() *
                           ui.parcelAreaMeanSlider->value() * 5.0f;
  float parcel_area_deviation = ui.parcelAreaDeviationSlider->value();
  ui.parcelAreaGroupBox->setTitle(
    "Area: "
    + QString::number(parcel_area_mean)
    + " m2 +- "
    + QString::number(parcel_area_deviation)
    + " %"
  );

  if (!initModeActive) {
    int selected = glWidget3D->cg.selectedPlaceTypeIdx;
    glWidget3D->cg.geoPlaceTypes[selected]["pt_parcel_area_deviation"] =
      parcel_area_deviation;
    //glWidget3D->generateGeometry(ClientGeometry::kStartFromRoads);
    glWidget3D->generateGeometry(ClientGeometry::kJustRoads);
    glWidget3D->updateGL();
  }
}//

void LCUrbanMain::parcelSplitDeviationSlot(int val) {
  if (!initModeActive) {
    int selected = glWidget3D->cg.selectedPlaceTypeIdx;
    glWidget3D->cg.geoPlaceTypes[selected]["pt_parcel_split_deviation"] = val *
        0.01f;
    //glWidget3D->generateGeometry(ClientGeometry::kStartFromRoads);
    glWidget3D->generateGeometry(ClientGeometry::kJustRoads);
    glWidget3D->updateGL();
  }
}//

void LCUrbanMain::parcelSetbackFrontSlot(int val) {

  float parcel_setback_front = 0.1f * (float)val;
  ui.ParcelSetbackFrontGroupBox->setTitle(
    "Front: "
    + QString::number(parcel_setback_front)
    + " m"
  );

  if (!initModeActive) {
    int selected = glWidget3D->cg.selectedPlaceTypeIdx;
    glWidget3D->cg.geoPlaceTypes[selected]["pt_parcel_setback_front"] =
      parcel_setback_front;
    //glWidget3D->generateGeometry(ClientGeometry::kStartFromRoads);
    glWidget3D->generateGeometry(ClientGeometry::kJustRoads);
    glWidget3D->updateGL();
  }
}

void LCUrbanMain::parcelSetbackSidesSlot(int val) {
  float parcel_setback_sides = 0.1f * (float)val;
  ui.ParcelSetbackSidesGroupBox->setTitle(
    "Sides: "
    + QString::number(parcel_setback_sides)
    + " m"
  );

  if (!initModeActive) {
    int selected = glWidget3D->cg.selectedPlaceTypeIdx;
    glWidget3D->cg.geoPlaceTypes[selected]["pt_parcel_setback_sides"] =
      parcel_setback_sides;
    //glWidget3D->generateGeometry(ClientGeometry::kStartFromRoads);
    glWidget3D->generateGeometry(ClientGeometry::kJustRoads);
    glWidget3D->updateGL();
  }
}//

void LCUrbanMain::landUseParkPercentageSlot(int val) {
  float land_use_park_percentage = 0.01f * (float)val;
  ui.LandUseParkPercentGroupBox->setTitle(
    "% Park: "
    + QString::number(val)
    + " %"
  );

  if (!initModeActive) {
    int selected = glWidget3D->cg.selectedPlaceTypeIdx;
    glWidget3D->cg.geoPlaceTypes[selected]["pt_park_percentage"] =
      land_use_park_percentage;
    //glWidget3D->generateGeometry(ClientGeometry::kStartFromRoads);
    glWidget3D->generateGeometry(ClientGeometry::kJustRoads);
    glWidget3D->updateGL();
  }
}//


void LCUrbanMain::buildingHeightMeanSlot(int val) {
  float building_height_mean = (int)(0.01f * (float)(
                                       ui.buildingHeightMeanSlider->value() * ui.buildingHeightMeanSlider->value()));
  float building_height_deviation = 1.0f * (float)
                                    ui.buildingHeightDeviationSlider->value();
  ui.BuildingHeightGroupBox->setTitle(
    "Num Floors: "
    + QString::number(building_height_mean)
    + " +- "
    + QString::number(building_height_deviation)
    + " %"
  );

  if (!initModeActive) {
    int selected = glWidget3D->cg.selectedPlaceTypeIdx;
    glWidget3D->cg.geoPlaceTypes[selected]["pt_building_height_mean"] =
      building_height_mean;
    glWidget3D->generateGeometry(ClientGeometry::kJustRoads);
    //glWidget3D->generateGeometry(ClientGeometry::kStartFromRoads);
    glWidget3D->updateGL();
  }
}//

void LCUrbanMain::buildingHeightDeviationSlot(int val) {
  float building_height_mean = (int)(0.01f * (float)(
                                       ui.buildingHeightMeanSlider->value() * ui.buildingHeightMeanSlider->value()));
  float building_height_deviation = 1.0f * (float)
                                    ui.buildingHeightDeviationSlider->value();

  ui.BuildingHeightGroupBox->setTitle(
    "Num Floors: "
    + QString::number(building_height_mean)
    + " � "
    + QString::number(building_height_deviation)
    + " %"
  );

  if (!initModeActive) {
    int selected = glWidget3D->cg.selectedPlaceTypeIdx;
    glWidget3D->cg.geoPlaceTypes[selected]["pt_building_height_deviation"] =
      building_height_deviation;
    //glWidget3D->generateGeometry(ClientGeometry::kStartFromRoads);
    glWidget3D->generateGeometry(ClientGeometry::kJustRoads);
    glWidget3D->updateGL();
  }
}//

void LCUrbanMain::onRender2DSlot(bool) {
  printf("onRender2DSlot\n");
  ui.action3D->setChecked(false);
  G::global()["render_mode"] = 1; //2D
  glWidget3D->changeRender2D3D();
  glWidget3D->updateGL();
}//

void LCUrbanMain::onRender3DSlot(bool) {
  printf("onRender3DSlot\n");
  ui.action2D->setChecked(false);
  G::global()["render_mode"] = 0; //3D
  glWidget3D->changeRender2D3D();
  glWidget3D->updateGL();
}//

// TERRAIN

void LCUrbanMain::updateTerrainLabels(int newValue) {
  int size = ui.terrainPaint_sizeSlider->value();
  ui.terrainPaint_sizeLabel->setText("Size: " + QString::number(size) + "%");
  G::global()["2DterrainEditSize"] = size / 100.0f;

  float change = ui.terrainPaint_changeSlider->value() * 500 /
                 100.0f; //*1785/100.0f;
  ui.terrainPaint_changeLabel->setText("Ch: " + QString::number(change, 'f',
                                       0) + "m");
  G::global()["2DterrainEditChange"] = change;
}//

void LCUrbanMain::onNewTerrain() {
  printf("No implemented yet\n");
  /*TerrainSizeInputDialog dlg(this);
              if (dlg.exec() == QDialog::Accepted) {
              glWidget->vboRenderManager.changeTerrainDimensions(dlg.side,dlg.cellResolution);
              glWidget->updateGL();
              }*/
}

void LCUrbanMain::onLoadTerrain() {
  QString filename = QFileDialog::getOpenFileName(this,
                     tr("Load Terrain file..."), "", tr("Terrain Files (*.png *.jpg)"));

  if (filename.isEmpty()) {
    return;
  }

  glWidget3D->vboRenderManager.vboTerrain.loadTerrain(filename);
  glWidget3D->updateGL();
}

void LCUrbanMain::onSaveTerrain() {
  QString filename = QFileDialog::getSaveFileName(this,
                     tr("Save Terrain file..."), "", tr("Terrain Files (*.png)"));

  if (filename.isEmpty()) {
    return;
  }

  glWidget3D->vboRenderManager.vboTerrain.saveTerrain(filename);
}

////////////////////////////////////
// LAYER
//////////////////////////////////////

void LCUrbanMain::onLoadLayers(bool) {
  printf("onLoadLayers\n");
  glWidget3D->vboRenderManager.layers.loadFromFile();
  int actiTex = ui.layer_numberSpinBox->value() % 2; //just two textures right now
  glWidget3D->vboRenderManager.layers.layers[actiTex].updateImage();
  //refresh
  glWidget3D->updateGL();
}//
void LCUrbanMain::onSaveLayers(bool) {
  printf("onSaveLayers\n");
  glWidget3D->vboRenderManager.layers.saveToFile();
}//

void LCUrbanMain::onClearLayers(bool) {
  printf("onClearLayers\n");
  glWidget3D->vboRenderManager.layers.clearLayers();
  int actiTex = ui.layer_numberSpinBox->value() % 2; //just two textures right now
  glWidget3D->vboRenderManager.layers.layers[actiTex].updateImage();
  //refresh
  glWidget3D->updateGL();
}//

void LCUrbanMain::onLayerEnable(int) {
  printf(">>onLayerEnable\n");

  if (ui.layer_enableCheckBox->isChecked()) {
    //activate layer mode
    printf("onLayerEnable Active\n");
    ui.action2D->setChecked(true);
    ui.action3D->setChecked(false);
    G::global()["render_mode"] = 2; //layer
    glWidget3D->changeRender2D3D();
    //update texture
    int actiTex = ui.layer_numberSpinBox->value() % 2; //just two textures right now
    glWidget3D->vboRenderManager.layers.layers[actiTex].updateImage();
    //refresh
    glWidget3D->updateGL();
  } else {
    //activate normal mode
    printf("onLayerEnable DESActive\n");

    if (ui.action2D->isChecked()) {
      G::global()["render_mode"] = 1; //2D
    } else {
      G::global()["render_mode"] = 0; //3D
    }

    glWidget3D->changeRender2D3D();
    //refresh
    glWidget3D->updateGL();
  }
}//


////////////////////////////////////
// TRAFFIC
//////////////////////////////////////

void LCUrbanMain::onArterialEdgesSpeedSlot(int val) {
  G::global()["cuda_arterial_edges_speed_ms"] = 1000.0f * val /
      3600.0f; //km/h--> m/s
  ui.ArterialEdgesSpeedGroupBox->setTitle(
    "Arterial Speed: "
    + QString::number(val)
    + " km/h"
  );

  if (!initModeActive) {
    //glWidget3D->generateGeometry(ClientGeometry::kStartFromRoads);
    glWidget3D->generateGeometry(ClientGeometry::kJustRoads);
    glWidget3D->updateGL();
  }
}//

void LCUrbanMain::onSimulatePressed(bool) {
  /////////////////////////
  printf(">>People and layers\n");
  int numPeople = ui.cudaNumPeopleSpinBox->value();

  //if(glWidget3D->infoLayers.layers.size()<=0)
  //	glWidget3D->infoLayers.initInfoLayer();
  printf(">>Init simulator\n");

  float deltaTime =
    ui.cudaDeltaTimeSpinBox->value(); //LC::misctools::Global::global()->cuda_delta_time;
  float cellSize = ui.cudaCellSizeSpinBox->value(); //meter

  printf(">>onSimulatePressed deltaTime %f cellSize %f numPeople %d\n", deltaTime,
         cellSize, numPeople);
  //exit(0);

  glWidget3D->cudaTrafficSimulator.initSimulator(deltaTime, cellSize,
      &glWidget3D->cg.roadGraph, numPeople, glWidget3D->vboRenderManager.layers,
      this);

  //////////////////////////

  //if(ui.cudaCPUSimulateButton==(QPushButton*)QObject::sender()){
  printf(">>CUDACPUSimulatePressed\n");
  int numPasses = ui.cudaShortestPathNumPassesSpinBox->value();
  glWidget3D->cudaTrafficSimulator.simulateInCPU_MultiPass(numPasses);
  /*}else{
                      printf(">>CUDAGPUSimulatePressed\n");
                      printf(">>Make calls\n");
                      if(mGLWidget_3D->cudaTrafficSimulator.trafficPersonVec.size()!=numPeople)
                              mGLWidget_3D->cudaTrafficSimulator.generateTrafficPersonJob();
                      else
                              mGLWidget_3D->cudaTrafficSimulator.resetPeopleJobANDintersections();

                      mGLWidget_3D->cudaTrafficSimulator.createLaneMap();
                      mGLWidget_3D->cudaTrafficSimulator.generateCarPaths();
                      mGLWidget_3D->cudaTrafficSimulator.simulateInGPU();
                      //mGLWidget_3D->cudaTrafficSimulator.calculateAndDisplayTrafficDensity();
              }*/

  //calculate density automatically
  //CUDATimeSliderSlot(ui.cudaTimeSlider->value());//update timeSlider
  printf("<<CUDACPUSimulatePressed\n");
}//



void LCUrbanMain::onLabelPressed(bool) {
  if (ui.TB_roadLabelButton->isChecked() == true) {
    //generate labels
    VBORoadLabels::updateRoadLabels(glWidget3D->vboRenderManager,
                                    glWidget3D->cg.roadGraph);
  }
}//

void LCUrbanMain::onDynameqDemand(int) {
  // when dynameq demand, we do not choose the dumber of people
  if (ui.bDynameqDemandCheckBox->isChecked()) {
    ui.bNumPeopleSpinBox->setEnabled(false);
  } else {
    ui.bNumPeopleSpinBox->setEnabled(true);

  }
}//

void LCUrbanMain::onBSimulateCPUPressed(bool) {
  // create/load things if neccesary
  onCreatePeople(true);

  ///////////////////////////////////////
  // SIMULATION
  // time
  QTime sTime = ui.stratTimeEdit->time();
  QTime eTime = ui.endTimeEdit->time();
  float startTime = sTime.hour() + sTime.minute() / 60.0f;
  float endTime = eTime.hour() + eTime.minute() / 60.0f;
  //simulate
  glWidget3D->bTrafficSimulator.simulateCPU(startTime, endTime);
}//

void LCUrbanMain::onBSimulateGPUPressed(bool) {
  // create/load things if neccesary
  onCreatePeople(true);

  ///////////////////////////////////////
  // SIMULATION
  // time
  QTime sTime = ui.stratTimeEdit->time();
  QTime eTime = ui.endTimeEdit->time();
  float startTime = sTime.hour() + sTime.minute() / 60.0f;
  float endTime = eTime.hour() + eTime.minute() / 60.0f;
  //simulate
  glWidget3D->bTrafficSimulator.simulateGPU(startTime, endTime);
}//

void LCUrbanMain::onCreatePeople(bool) {
  QTime timer;
  timer.start();

  /////////////////////////////////////////////
  // INIT
  if (glWidget3D->bTrafficSimulator.initialized == false) {
    //float deltaTime = ui.bDeltaTimeSpinBox->value();
    //float cellSize = ui.bCellSizeSpinBox->value();
    printf(">>1. Init simulator\n");
    glWidget3D->bTrafficSimulator.init(&glWidget3D->cg.roadGraph,
                                       this); //deltaTime, cellSize,
    printf("<<1. Init simulator %d ms\n", timer.elapsed());
  } else {
    printf("<<1. Already initiated\n");
  }

  ////////////////////////////////////////////
  // CREATE PEOPLE
  if (glWidget3D->bTrafficSimulator.people.numPeople <= 0) {
    printf(">>2. create people\n");
    timer.restart();
    QTime sDemandTime = ui.demandStratTimeEdit->time();
    QTime eDemadTime = ui.demandEndTimeEdit->time();
    float startDemandTime = sDemandTime.hour() + sDemandTime.minute() / 60.0f;
    float endDemandTime = eDemadTime.hour() + eDemadTime.minute() / 60.0f;

    if (ui.bDynameqDemandCheckBox->isChecked() == false) {
      printf(">>2.1 create people\n");
      int numPeople = ui.bNumPeopleSpinBox->value();

      if (glWidget3D->vboRenderManager.layers.layersEmpty() == true) {
        printf(">>2.Layers empty--> Try to load from file\n");
        onLoadLayers(true);
      } else {
        printf(">>2.Layers not empty\n");
      }

      printf(">>2. createPeople %d\n", numPeople);
      glWidget3D->bTrafficSimulator.createPeople(startDemandTime, endDemandTime,
          numPeople, glWidget3D->vboRenderManager.layers);
    } else {
      printf(">>2. createPeopleDynameq\n");
      glWidget3D->bTrafficSimulator.createPeopleDynameq(startDemandTime,
          endDemandTime);
    }

    printf("<<2. Create people %d ms\n", timer.elapsed());
  } else {
    printf("<<2. Already people\n");
  }

  printf("<<2. Already routes %d\n",
         glWidget3D->bTrafficSimulator.people.nextEdge.size());

  ////////////////////////////////////////////
  // CREATE MAP
  if (glWidget3D->bTrafficSimulator.laneMapL[0].size() <= 0) {
    printf(">>3. Create lane map\n");
    timer.restart();
    glWidget3D->bTrafficSimulator.createLaneMap();
    printf("<<3. Create lane map %d ms\n", timer.elapsed());
  } else {
    printf("<<3. Already lane map\n");
  }

  ////////////////////////////////////////////
  // CREATE INTERSECTIONS
  if (glWidget3D->bTrafficSimulator.intersec.numIn.size() <= 0) {
    printf(">>3. Create Intersections\n");
    timer.restart();
    glWidget3D->bTrafficSimulator.createIntersections();
    printf("<<3. Create Intersections %d ms\n", timer.elapsed());
  } else {
    printf("<<3. Already Intersections\n");
  }

  ////////////////////////////////////////////
  // CREATE ROUTES
  glWidget3D->bTrafficSimulator.routesChanged = true;

  if (glWidget3D->bTrafficSimulator.people.nextEdge.size() <= 0) {
    printf(">>4. Create routes\n");
    timer.restart();
    bool dijkstraRouting = ui.bRoutingDijkstraCheckBox->isChecked();

    if (dijkstraRouting == true) {
      printf("dijkstraRouting\n");
      glWidget3D->bTrafficSimulator.generateCarDijstra();
    } else {
      printf("johnsonRouting\n");
      glWidget3D->bTrafficSimulator.generateCarJohnson();
    }

    printf("<<4. Create routes D %d ms\n", timer.elapsed());
  } else {
    printf("<<4. Already routes %d\n",
           glWidget3D->bTrafficSimulator.people.nextEdge.size());
  }

  printf("<< onCreatePeople\n");
}//

void LCUrbanMain::onLoadPeople(bool) {
  glWidget3D->bTrafficSimulator.people.loadFromFile();
  printf("onLoadPeople: Loaded %d\n",
         glWidget3D->bTrafficSimulator.people.numPeople);
}//

void LCUrbanMain::onSavePeople(bool) {
  glWidget3D->bTrafficSimulator.people.saveToFile();
  printf("onSavePeople: Saved %d\n",
         glWidget3D->bTrafficSimulator.people.numPeople);
}//

void LCUrbanMain::onProceduralModeling(int) {
  glWidget3D->bTrafficSimulator.clearSimulator();//clear people

  if (ui.proceduralModelingCheckBox->isChecked() == true) {
    glWidget3D->generateGeometry(ClientGeometry::kStartFromRoads);
  } else {
    RoadGraphDynameq::loadDynameqRoadGraph(glWidget3D->cg.roadGraph, glWidget3D);
    glWidget3D->changeRender2D3D();//force generate VBOs
  }

  printf("<< onProceduralModeling\n");
  glWidget3D->updateGL();
}//

/*void LCUrbanMain::updateViewCheckboxesSlot(bool){
      //--interactivity
      LC::misctools::Global::global()->render_lighting =
      ui.RenderLightingCheckbox->isChecked();
      //LC::misctools::Global::global()->selectedParametricBuildings=LC::misctools::Global::global()->render_lighting;

      //--place types
      LC::misctools::Global::global()->view_place_types_centept_render =
      ui.RenderPlaceTypesCheckbox->isChecked();

      //--arterials

      LC::misctools::Global::global()->view_arterials_vertices_render =
      ui.RenderArterialsVerticesCheckbox->isChecked();

      LC::misctools::Global::global()->view_arterials_edges_render =
      ui.RenderArterialsEdgesCheckbox->isChecked();//

      LC::misctools::Global::global()->view_arterials_arrows_render =
      ui.RenderRoadArrowsCheckbox->isChecked();//

      //--buildings
      LC::misctools::Global::global()->view_buildings_render =
      ui.RenderBuildingsCheckbox->isChecked();

      //--parcels
      LC::misctools::Global::global()->render_parcelTextures =
      ui.RenderParcelTexturesCheckbox->isChecked();

      LC::misctools::Global::global()->figureMode=
      ui.cudaFigureModeCheckbox->isChecked();

      if(!initModeActive){
      mGLWidget_3D->updateGL();
      }
      }//*/


} // namespace LC
