# DropTrack

DropTrack implements a pipeline to detect and track drops in high-speed image
sequences of a test cell developed at the Chair of Chemical and Process Engineering
of Technische Universität Berlin. This automated image analysis tool was developed
in cooperation with the Computer Vision and Remote Sensing laboratory of Technische
Universität Berlin to analyse the huge amount of recorded image sequences with
varying resolutions and qualities. It is able to determine the trajectories of two
colliding drops as well as the important events of drop detachment from cannulas
and their collision. With this information the drop velocity in each sequence is
calculated and mean values of multiple drop collisions are determined for serial
examinations of single drop collisions.

For further details and downloads please visit http://www.rhaensch.de/droptrack.html
For questions or remarks don't hesitate to contact: johannes.kamp@tu-berlin.de

If you use this code in your work, please refer to the following paper:

``` bibtex
@ARTICLE{Kamp2016,
    author = {Kamp, Johannes and H\"ansch, Ronny and Kendzierski, Gregor and Kraume, Matthias and Hellwich, Olaf},
    title  = {Automated image analysis for trajectory determination of single drop collisions},
    journal= {Computers \& Chemical Engineering},
    year   = {2016},
    doi    = {10.1016/j.compchemeng.2016.03.033},
    url    = {http://dx.doi.org/10.1016/j.compchemeng.2016.03.033}
}
```

DropTrack (c) 2015, Johannes Kamp, Ronny Hänsch, Gregor Kendzierski and contributors.

DropTrack is free software: you can redistribute it and/or modify it under the terms
of the BSD 2-clause License.
DropTrack is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the license file provided with the code for the license terms and more details.

Created by Gregor Kendzierski, December 2015

# Requirements and Dependency

- MatLab 2014b with Image Processing Toolbox

# Installation
- download zip-File and unzip
- copy `config_example.m` and rename it to `config.m`
- set `ROOT` variable as global path to DropTrack files in `config.m`

# Usage and Manuals
- automated image analysis: [record_EXP](man/record_EXP.md) (Please take notice of the selected [coordinate system](man/record_EXP.md#coordinate-system).)
- post-processing: [view_SEQ](man/view_SEQ.md)

# Examples
- to evaluate a single record, see [examples/run_record_EXP.m](examples/run_record_EXP.m)
- to post-process parameterSets, see [examples/run_view_SEQ.m](examples/run_view_SEQ.m)

# Troubleshooting

### strjoin() is undefined

```matlab
??? Undefined function or method 'strjoin' for input arguments of type 'cell'.
```
If your Matlab Version < 2013a then you have to import the `strjoin()` function from the [Matlab FileExchange](http://www.mathworks.com/matlabcentral/fileexchange/31862-strjoin)
