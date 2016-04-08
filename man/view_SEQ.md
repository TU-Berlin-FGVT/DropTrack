# VIEW\_SEQ

The class `view_SEQ()` is used to generate averaged plots of records with same parameters of rising drops after evaluation with [`record_EXP()`](record_EXP.md).

Created by Gregor Kendzierski,
December 2015

# general remarks

## example

### 1. Loading

```matlab
config;
path = 'full_path_to_parameter_directory';
seq = view_SEQ(path);
seq.load_data();
```

### 2.Select record

```matlab
%% with same parameterSet
parameterSets = seq.label_parameterSet();

% select first parameterSet
selected_parameterSet = parameterSet == 1;
```

### 3. Generate averaged plot of selected records

```matlab
fg = seq.subplot_averaged(selected_parameterSet);
```

# Units

| physical quantity | Image-Units | SI-Units |
|--- | --- | --- |
| length | px   | mm |
| time | frame  | s |
| velocitiy | px/frames| mm/s |

# class description

## properties

| name | default | type | discription|
| ---  | --- | ---  | ---  |
|path  | - | char| path in which the records with same parameter can be found|
|DATA | - |cell array of struct | cell array of single records' data |
|DATA_path | - | cell array of char |cell array of paths for each record |
|dataFileName | 'exp_data.mat'| char | name of the records' datafile |
|SIunits      | false | bool(1,1)| Output Units: </br>true - SI-Units, </br>false - Image-Units |

## public methods

### obj = view\_SEQ(path)
 constructor of the view\_SEQ class

| name  | datatype  | description  |
| ---  | --- | ---  |
| path | char  | path to parameterSet's directory |
| obj  | view\_SEQ() |  |

### out = load\_data()
 method to load all data and reduce to required only data

### [out,parameterSetValues] = label\_same\_parameterSet()
 method to return labeling array to mask records with same parameterSets

| name  | datatype  | description  |
| ---  | --- | ---  |
| out   | double(1,n)| out = 0 , not selected record</br> out > 0 , same paramterSets have the same numbered label |
|parameterSetValues| struct | contains labels and parameterValues of each Set of selected record |

### out = label\_same\_temperature(arraylabel)
 method to return labeling array to mask records with same temperature

| name  | datatype  | description  |
| ---  | --- | ---  |
| arraylabel | bool(n,1)| selected records |
| out   | double(n,1)| out =-1 , noTemperature is given </br>out = 0 , not selected record <br>out > 0 , same temperatures of selected records have the same numbered label |

### out = label\_same\_gitVersion(arraylabel)
 method to return labeling array to mask records with same gitVersions

| name  | datatype  | description  |
| ---  | --- | ---  |
| arraylabel | bool(n,1)| selected records |
| out   | double(1,n)| out =-1 , no gitVersion is given </br>out = 0 , not selected record </br>out > 0 , same gitVersion of selected records have the same numbered label |

### set_SIunits(value)
 method to set output units to SI-Units or Image-Units

| name  | datatype  | description  |
| ---   | --- | ---  |
| value | bool(n,1)| true -> SI-Units,</br> false -> Image-Units |

### [px2mm, fps] = get_mean_SIunitsFactors(obj,labelarray)
 method to estimate the converting factors to SI-Units

| name  | datatype  | description  |
| ---  | --- | ---  |
| arraylabel | bool(n,1)| selected records |
| px2mm | double(1,1)| [mm/px] factor to convert to SI-Units |
| fps | double(1,1)| [frames/s] factor to convert to SI-Units |

### fh = subplot\_averaged(arraylabel)
 method to plot the averaged data of selected records

**See**:
  - [plot\_averaged\_data()](#ah-plot95averaged95dataah-data-varargin),
  - [plot\_averaged\_height\_over\_time()](#ah-plot95averaged95height95over95timearraylabel-ah-varargin),
  - [plot\_averaged\_risingVelocity\_over\_time](#ah-plot95averaged95risingvelocity95over95timearraylabel-ah-varargin),
  - [plot\_averaged\_risingVelocity\_over\_height](#ah-plot95averaged95risingvelocity95over95heightarraylabel-ah-varargin)

| name  | datatype  | description  |
| ---  | --- | ---  |
| arraylabel | bool(n,1)| selected records |
| fh | figure handle | handle of the current figure |

### ah = plot\_averaged\_data(ah, data, varargin)
 method to plot averaged data generally

| name  | datatype  | description  |
| ---  | --- | ---  |
| ah | axes handle | handle of the current axes |
| data | struct | contains the data arrays that will be plotted
| varargin | cell array | contains switch off keywords to adjust the plot |

**See**:
  - [plot\_averaged\_height\_over\_time()](#ah-plot95averaged95height95over95timearraylabel-ah-varargin),
  - [plot\_averaged\_risingVelocity\_over\_time](#ah-plot95averaged95risingvelocity95over95timearraylabel-ah-varargin),
  - [plot\_averaged\_risingVelocity\_over\_height](#ah-plot95averaged95risingvelocity95over95heightarraylabel-ah-varargin)

#### data struct

|varname        | datatype | description |
| ---           | ---| --- |
|X              | doubel(n,1) | scale on the x-axes |
|Y_mean         | doubel(n,1) | mean data over x-axes |
|Y_std          | doubel(n,1) | mean +/-std data over x-axes|
|Y_max          | doubel(n,1) | max(mean) data over x-axes|
|Y_min          | doubel(n,1) | min(mean) data over x-axes|
|bottomDetachmentEvent|doubel(1,1) |  vertical line|
|collisionEvent | doubel(1,1) | vertical line |
|numberOfRecords| doubel(n,1) | graph  as vertical line over x-axes |
|arraylabel     | bool(n,1) | selected records |

#### plot varargin keywords

|keyword        | description |
| ---           | ---|
|`OverwriteOFF` | switch off to overwrite the plot, just adding to the plot |
|`LegendOFF`    | switch off the the legend |
|`NumberOfRecordOFF` |switch off the number of record which are used for averaging of each timestep |

See also [warning varargin keywords](#warning-varargin-keywords) to surpress warnings.

### ah = plot\_averaged\_height\_over\_time(arraylabel, ah, varargin)
 method to plot the averaged height over time

| name  | datatype  | description  |
| ---  | --- | ---  |
| ah | axes handle | handle of the current axes |
| arraylabel | bool(n,1)| selected records |
| varargin | cell array | contains switch off keywords to adjust the plot, see [plot varargin keywords](#plot-varargin-keywords) |

**See**:
  - [plot\_averaged\_data()](#ah-plot95averaged95dataah-data-varargin)
  - [warning varargin keywords](#warning-varargin-keywords) to surpress warnings

### ah = plot\_averaged\_risingVelocity\_over\_time(arraylabel, ah, varargin)
 method to plot the averaged risingVelocity over time

| name  | datatype  | description  |
| ---  | --- | ---  |
| ah | axes handle | handle of the current axes |
| arraylabel | bool(n,1)| selected records |
| varargin | cell array | contains switch off keywords to adjust the plot, see [plot varargin keywords](#plot-varargin-keywords) |

**See**:
  - [plot\_averaged\_data()](#ah-plot95averaged95dataah-data-varargin)
  - [warning varargin keywords](#warning-varargin-keywords) to surpress warnings

### ah = plot\_averaged\_risingVelocity\_over\_height(arraylabel, ah, varargin)
 method to plot the averaged risingVelocity over height

| name  | datatype  | description  |
| ---  | --- | ---  |
| ah | axes handle | handle of the current axes |
| arraylabel | bool(n,1)| selected records |
| varargin | cell array | contains switch off keywords to adjust the plot, see [plot varargin keywords](#plot-varargin-keywords) |

**See**:
  - [plot\_averaged\_data()](#ah-plot95averaged95dataah-data-varargin)
  - [warning varargin keywords](#warning-varargin-keywords) to surpress warnings

### out = get\_bottomDetachmentTime(arraylabel)
 method to get the mean of the bottomDetachmentTime intervall for each selected record

| name  | datatype  | description  |
| ---  | --- | ---  |
| arraylabel | bool(n,1) | selected records |
| out | cell array | contains the mean bottomDetachmentTime of each seleced record |

### [timeshift, arraylabel] = get\_timeshift\_bottomDetachmentEvent(arraylabel)
 method to get the alginement in time with the bottomDetachmentEvent for each selected record, records with no detected event will be ignored

| name  | datatype  | description  |
| ---  | --- | ---  |
| arraylabel | bool(n,1)| selected records, </br> unselect record which have no detected bottomDetachmentEvent |
| timeshift | cell array | contains the required alignement for each record |

### out = get\_averaged\_bottomDetachmentHeight(arraylabel, heightshift)
 method to get the averaged height of the bottomDetachmentEvent of selected records

| name  | datatype  | description  |
| ---  | --- | ---  |
| arraylabel | bool(n,1)|  selected records |
| heightshift| double(1,1) | contains the required alignement for each record |
| out | double(1,1)| averaged bottomDetachmentHeight of selected records |

### out = get\_collisionTime(arraylabel)
 method to get the mean of the collisionTime intervall for each selected record

| name  | datatype  | description  |
| ---  | --- | ---  |
| arraylabel | bool(n,1)|  selected records |
| out | double(1,1)| mean collisonTime of each record |

### out = get\_averaged\_collisionTime(arraylabel, timeshift)
 method to get the averaged time of the collisionEvent of selected records

| name  | datatype  | description  |
| ---  | --- | ---  |
| arraylabel | bool(n,1)|  selected records |
| timeshift | cell array | contains the required alignement for each record |
| out | double(1,1)| averaged collisonTime of selected records |

### out = get\_averaged\_collisionHeight(arraylabel, heightshift)
 method to get the averaged collisionHeight of selected records

| name  | datatype  | description  |
| ---  | --- | ---  |
| arraylabel | bool(n,1)|  selected records |
| heightshift| double(1,1) | contains the required alignement for each record |
| out | double(1,1)| averaged collisonHeight of selected records |

### out = get\_averaged\_height\_over\_time(arraylabel)
 method to average height of seleced records for each timestep

| name  | datatype  | description  |
| ---  | --- | ---  |
| arraylabel | bool(n,1)|  selected records |
| out | struct | contains the data arrays that will be plotted, see [data struct](#data-struct)|

### out = get\_averaged\_risingVelocity\_over\_time(arraylabel)
 method to average risingVelocity of seleced records for each timestep

| name  | datatype  | description  |
| ---  | --- | ---  |
| arraylabel | bool(n,1)|  selected records |
| out | struct | contains the data arrays that will be plotted, see [data struct](#data-struct)|

### [height, risingVelocity,arraylabel] = get\_strict\_monotone\_height(arraylabel)
 method to force strict monotone decreasing height by croping intervalls

| name  | datatype  | description  |
| ---  | --- | ---  |
| arraylabel | bool(n,1)| *IN*: selected records records</br> *OUT*: records which can be enforced to strict monotonie |
| height | cellarray of double | contains the strict monotone decreasing height|
| risingVelocity | cellarray of double | contains the associated risingVelocity to the height |

### out = get\_averaged\_risingVelocity\_over\_height(arraylabel)
 method to average risingVelocity of seleced records for each height which is forced to stric monotonie

| name  | datatype  | description  |
| ---  | --- | ---  |
| arraylabel | bool(n,1)|  selected records |
| out | struct | contains the data arrays that will be plotted, see [data struct](#data-struct)|

## private methods

### out = is\_cleanVersions(arraylabel)
 method to check if selected sequences created by commit clean gitVersion

| name  | datatype  | description  |
| ---  | --- | ---  |
| arraylabel | bool(n,1)|  selected records |
| out | bool(1,1) | wheter all selected records created by clean gitVersion |

### out = is\_same\_gitVersion(arraylabel)
 method to check if selected sequences have the same commitID

| name  | datatype  | description  |
| ---  | --- | ---  |
| arraylabel | bool(n,1)|  selected records |
| out | bool(1,1) | wheter all selected records have the same commitID |

### plot\_warnings( arraylabel, varargin)
 method to promp up warnings when plotting records with different properties

| name  | datatype  | description  |
| ---  | --- | ---  |
| arraylabel | bool(n,1)|  selected records |

#### warning varargin keywords

|keyword        | description |
| ---           | ---|
|`warningOFF`| switch off all warnings|
|`checkVersionsOFF` | switch off gitVersion checking |

See also [plot varargin keywords](#plot-varargin-keywords)
