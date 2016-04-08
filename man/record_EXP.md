# RECORD_EXP

The class `record_EXP()` is used to process image data of single records generated at Technische Universität Berlin.

Created by Gregor Kendzierski, December 2015

1. [general remarks](#general-remarks)
 1.1 [example](#example)
 1.2 [units](#units)
 1.3 [coordinate system](#coordinate-system)
2. [class description](#class-description)
 2.1. [properties](#properties)
 2.2. [public methods](#public-methods)
 2.3. [private methods](#private-methods)

# general remarks

## example

```matlab
%% process one record
config;
path = 'full_path_to_the_record';
rec = record_EXP(path);
rec.run_all();

%% access to processed data of the record
rec.drops(2).centers
rec.Events.bottomDetachmentTime
rec.Events.bottomDetachmentTime

%% access to saved processed data of the record
load([path,'\exp_data.mat']);
data.drops(2).centers
data.Events.bottomDetachmentTime
data.Events.bottomDetachmentTime
```
## units

| physical quantity | unit   |
|--- | --- |
| length | px   |
| time | frame  |
| velocitiy | px/frames|
| framerate | s/frame |
| px2mm | mm/px |

## coordinate system
The point of origin lays in the upper left corner of the rotated image:
- when the `x`-coordinate rises the direction is downwards (**v**)
- when the `y`-coordinate rises the direction is rightwards (**>**)

Therefore, the `x`-velocity is negative when the lower drop rises.

Points are defined as complex numbers:  `P(x,y) = x + y i`

# class description

## properties
| name             | default | type | discription|
| ---              | ---     | ---  | ---        |
|reduce_usageOfRam | true    | bool | reduce the usage of RAM by deleting big variables|
| path | | char | directory in which the experiment data  with the parameter.csv/mat can be found|
| parameters | | struct array| struct array of parameters of the test facility
| parameters_csv_file | 'parameters.csv'| char| CSV-Name of the parameters file|
parameters_mat_file | 'parameters.mat' |char | MAT-Name of the parameters file|
imageFileTypes | {'bmp';'jpg'}| char-array | all supported image file type|
|STDfactor| 2| double(1,1) | factor of standard deviation, needed for the estimation of the threshold, see [get\_threshold()](#get95thresholdobji) |
|min_EventDistance | 1| double(1,1) [px] |  minimal distance between objects for event detection, see [get\_bottomDetachmentTime()](#get95bottomDetachmentTime) and [get\_collisionTime()](#get95collisionTime)|
|smoothing\_range | 10 | double(1,1) [frames] |  PostProcessing smoothing range for moving average, see [smooth_dropsDynamic()](#smooth95dropsdynamicobj)|
| log | empty array | char-array | logs methods which are done |
|display_logs| false| bool| display msg which are added|

### obj.parameters
 struct array of the record's parameter, imported by parameter-file

 used parameters:

| name             | type | discription|
| ---              |  --- | ---        |
| Framerate | double(1,1) | framerate of the record |
| d_disp_top_target | double(1,1) | diameter of the top drop|
| d_disp_bot_target | double(1,1) | diameter of the bottom drop |
| hashID | char(32) | ID of the record |
| units  | sturct array | stores the units of all given parameters |

### obj.cam_system
 cell-array table which stores the properties of different camera systems

| name             | type | discription|
| ---              |  --- | ---        |
| obj.cam_system(:,1)| double(1,1) | camera system ID |
| obj.cam_system(:,1)| double(1,1) | maxRelStdNoises of the camera |

### obj.frames
 struct array of image data

| name       | type | discription|
| ---        | ---  | ---        |
| filepaths  | cell-array of char | single image paths |
| imgOriginal| cell-array of double | original images  |
| imgFiltered| cell-array of double | result of the filtered original images |
| imgBinary  | cell-array of logical| matrix which masks pixels of the objects, result after applied threshold|
| numberFrames|double(1,1) | number of image frames of the experiment|
| size| double(1,2)| dimension scale of the image frame matrix|
| px2mm| double(1,1)| mm/px-relation, scale to SI-Units |
| maxRelStdNoises| double(1,1) | maximum of the relative standard deviation of filtered pixel values, default value is set to 0.025 |

### obj.cannulas
 struct array of cannulas data

| name     | type   | discription|
| ---      | ---    | ---        |
| imgKum   | double | cumulated image of all grames |
| imgBinary| logical / cell array of logical | binarise images of `imgKum` / masks of cannuals pixel|
| tops_m| complex double(2,1)| coordinates of the midway cannula's top |
| tops| complex double(2,2)| coordinates of the edges of the cannual's top |
| roots| complex double(2,2) | coordinates of the edges of the cannual's root |
| type | double(1,1) | type of the cannual , (1 -> top, 2 -> bottom) |
| diameters |double(2,1) | diameter of the cannulas |
| length | double(2,1) | length of the cannuals |
| diameter_mm |double(1,1) | diameter of cannulas according to producer, default value is set to  0.8 |

### obj.drops
 struct array of drops for each drop

| name     | type   | discription|
| ---      | ---    | ---        |
| centers  | complex double(1,n) | coordinates of drop's center |
| edgePoints |cell array of complex double | coordinates of drop's edge points |
|mean_radii | complex double(1,n) | mean drop radius |
|min_radii  | complex double(1,n) | minimal drop radius |
|max_radii  | complex double(1,n) | maximal drop radius |
|velocities | complex double(1,n) | velocity of the drops center|
|validates  | logical(1,n)  | fulfill drop's criteria |

```matlab
    % get the center of the bottom drop at frame number 42
    rec.drops(2).center(42)
```

### obj.Events
 struct array of events data

Time intervall of each event is set to [0,0] as default.

| name     | type   | discription|
| ---      | ---    | ---        |
| bottomDetachmentTime| double(1,2) | Time when the bottom drop detachs from the lower cannula |
| collisionTime| double(1,2)| Time when both drop's surfaces touch each other|
| coalescenceTime|double(1,2)| Time when both drop's coalescence |
| repulsionTime | double(1,2)| Time when the drop contact breaks |

### obj.manualValues
 manual evalutated event data and preprocessing parameters

Time intervall of each event is set to [0,0] as default.

| name     | type   | discription|
| ---      | ---    | ---        |
| bottomDetachmentTime| double(1,2) | Time when the bottom drop detachs from the lower cannula |
| collisionTime| double(1,2)| Time when both drop's surfaces touch each other|
| coalescenceTime|double(1,2)| Time when both drop's coalescence |
| repulsionTime | double(1,2)| Time when the drop contact breaks |
| crop_coordinates| complex double(1,2)| diagonal rectangle edges to crope all frames in preprocessing |
| man_threshold | double(1,1)| set threshold by hand for binary frames after filtering|

## public methods

### obj = record\_EXP ( path, CAM\_ID)
  contructor method of the record\_EXP class

| name  | datatype    | description                |
| ---            | ---         | ---                        |
| path           | char        | path to record's directory |
| CAM\_ID        | double(1,1) | id of the camera system (optional)|
| obj            | record\_EXP |                            |

### out = record\_EXP.create\_ParameterMAT(path,framerate)
  static method to create parameterMAT-file from directory name with a minimum of information

Use this method if the parameter.csv-file is not available for an experiment. Take care to define the correct framerate!

| name      | datatype    | description |
| ---                | ---         | ---         |
| path               | char        | path to record's directory |
| framerate          | double(1,1) | framerate of the record    |
| out                | bool        | true on success |

### out = record\_EXP.create\_manualValueFile()
  static method to create the empty manualValue file

### run\_all(obj)
  execute all processing methods in the right order

### out =  detect\_imageCorruption(obj)
  method to detect image corruption caused by recording with maladjusted ratio of framerate and image size

| name      | datatype    | description        |
| ---                | ---         | ---                |
| out                | bool        | true on success |

### out = load\_frames(obj, i)
  method to load the images data of the experiment to obj.frames.imgOriginal as cellarray if i is not given

| name      | datatype    | description        |
| ---                | ---         | ---                |
| i                  | double(1,1) | timestep to get the i-th image of the frame|
| out                | bool/double        | true on success/image if `i` is given

### out = filter\_frames(obj,i)
  method to filter the images to reduce large range gradients along image length

| name      | datatype    | description        |
| ---                | ---         | ---                |
| i                  | double(1,1) | timestep to filter the i-th image of the frame|
| out                | bool/double        | true on success/image if `i` is given |

### binary\_frames(obj)
  method to generate image masks to select object pixels of every frame via threshold estimated for each frame

### detect\_cannuals(obj)
  method to find cannulas through selecting object pixel via specific criterion for cannulas

### rotate2top\_frames(obj)
  method to rotate the images, so that the interested drop is rising.

### reset\_DropProperties(obj)
  method to reset and initialize values for drop detection

### get\_dropCentersAfterCollsion(obj)
  method to  seperate drop pixel up into their associated drop types after collision via moving dividing line constructed through the collsion point

### get\_bottomDetachmentTime(obj)
  method to restrict the time interval of the bottomDetachmentTime

### get\_collisionTime(obj)
  method to restrict the interval of the collisionTime

### get\_coalescenceTime(obj)
  method to restrict the interval of the coalescenceTime

### get\_repulsionTime(obj)
  method to restrict the interval of the repulsionTime

### get\_contactTime(obj)
  method to restrict the interval of the contactTime

### get\_drainageTime(obj)
  method to restrict the interval of the drainageTime

### out = crop\_frames(obj,coords)
  method to crop frame by given two coordinates

| name      | datatype    | description        |
| ---                | ---         | ---                |
| coords             | complex double(2,1) | two points limits the croping rectangle|
| out                | bool        | is true on success |

### save\_manualValuesFiles(obj)
  method to save values set by hand

### disp\_fullLog(obj)
  method to diplay a list of messages created by method used with unusual adjustments for replicability

## private methods

### get\_threshold(obj,I)
  method to estimate threshold for given image based on properties of the camera system which is used

| name | datatype    | description        |
| ---           | ---         | ---                |
| I             | double(n,n) | is a 2-dim double array as image |

### set\_cannulasImgKum(obj)
  method to generate average image by arithmetic averaging pixel values with same image position

*input:*  `obj.frames.imgFiltered`

*output:* `obj.cannuals.imgKum`

### set\_cannulasImgBinary(obj)
  method to generate image masks to select object pixels of `obj.cannuals.imgKum` via threshold

*output:* `obj.cannuals.imgBinary`

### get\_cannulasAreas(obj)
  method to select and label cannulas pixel

*input:*  `obj.cannuals.imgBinary` double(n,n)

*output:* `obj.cannuals.imgBinary` cell array of the masks

### improve\_cannulasShape(obj)
  method to eliminate defects on the cannula masks

*input:*  `obj.cannuals.imgBinary` cell array of the masks

*output:* `obj.cannuals.imgBinary` cell array of the masks

### set\_cannulasExtrema(obj)
  method to set extrem quadrangle edge points of each cannula

*input:*  `obj.cannuals.imgBinary` cell array of the masks

*output:* `obj.cannuals.tops` , `obj.cannuals.roots`

### set\_cannulasDiameter(obj)
  method to estimate cannulas' diameter

*input:*  `obj.cannuals.imgBinary` cell array of the masks

*output:* `obj.cannuals.diameters`

### set\_cannulasLength(obj)
  method to estimate cannulas' length

*input:*  `obj.cannuals.imgBinary` cell array of the masks

*output:* `obj.cannuals.length`

### set\_cannulasTopsM(obj)
  method to estimate cannula mean top

*input:*  `obj.cannuals(:).tops`

*output:* `obj.cannuals(:).tops_m`

### set\_px2mm(obj)
  method to estimate the ratio of millimeter per pixel length

*input:*  `obj.cannuals.diameters`,`obj.cannulas.diameter_mm`

*output:* `obj.frames.px2mm`

### get\_dropCenters(obj)
  method to estimate and adjust drop center through drop edge

*output:* `obj.drops`

### get\_dropsEdgePoints(obj,i)
  method to get the the edge drop points which have the maximal distant to the drop centers

*output:* `obj.drops`

  TODO: split get\_DropCenters up to get\_dropsEdgePoints(obj,i)

### out =  get\_firstDropFrame(obj)
  method to get the frame number when the lower drop has been generated

| name | datatype    | description        |
| ---           | ---         | ---                |
| out           | double(1,1) | timestep of the event |

### out = cleared\_cannulasImgBinary(obj,i)
  method to eliminate cannulas and disturbing artefacts from binary image

| name | datatype    | description        |
| ---           | ---         | ---                |
| i             | double(1,1) | timestep of the image which should be cleared |

### smooth\_dropsDynamic(obj)
  method to smooth drop trajectory by moving average, see smoothing\_range property for adjustments
### get\_dropsVelocities(obj)
  method to estimate drop velocities via central difference quotient
### get\_dropsRadii(obj)
  method to estimate drop radii from given drop centers and drop edge

### add\_log(obj,msg)
  method to add message to log for replicability when methods are used with unusual adjustments
### load\_cameraSystem(obj, CAM\_ID)
  method to load properties of the used camera system. Properties can be added as

| name | datatype    | description        |
| ---           | ---         | ---                |
| CAM\_ID       | double(1,1) | id of the camera system, see `obj.cam_syst`|

### load\_parameters(obj)
  method to load parameter values
### load\_parametersCSV(obj)
  method to load parameter values from given CSV-file
### load\_parametersMAT(obj)
  method to load paramter values from converted MAT-file
### set\_hashID(obj)
  method to generate hashID from the first image of the record as unique identifier
### load\_manualValues(obj)
  method to load the manual values for this record based on its hashID
