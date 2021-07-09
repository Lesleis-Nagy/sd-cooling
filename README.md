# Introduction
The `sd-cooling` utility calculates TRM acquisition as a function of
temperature, grain size & shape, field strength and direction. It uses simple 
single domain Stoner-Wohlfarth grains.

# Compilation
This code can be compiled on linux Ubuntu (version 20.04.2 LTS) using the CMake
([https://cmake.org](https://cmake.org)) crossplatform make utility (minimum
version 3.17).  Additionally this code requires the Boost
([https://boost.org](https://boost.org)) library with `program-options` to be
compiled and available on the build system. The current version of the code has
been compiled with Boost version 1.76.0.

Once the code has been downloaded from the repository, a build directory
should be created somewhere on the system, we subsequently refer to this 
directory as `<build-dir>`. Then in order to compile 
```bash
    $> cd <build-dir>
    $> cmake -DCMAKE_BUILD_TYPE=Release <source-dir>
```
where `<source-dir>` is the directory which was acquired from the repository.
The `-DCMAKE_BUILD_TYPE=Release` flag is optional, however it produces a code
that should run quicker and if it is not present debugging information will be
sent to `stdout`.

It is possible that CMake encounters a problem finding the Boost libraries, 
in this case it may be useful to set the `-DBOOST_ROOT=<boost-dir>` flag
to point to the root directory of the Boost installation for example
```bash
    $> cmake -DCMAKE_BUILD_TYPE=Release -DBOOST_ROOT=<boost-dir> <source-dir>
```

Once CMake has completed successfully, you can build the executable using
make
```bash
    $> make
```
This command should result in an executable called `cooling` in the
`<build-dir>/src_cmdline/cooling` directory. 

# Usage
The `cooling` utility can be called from the command line
```bash
    $> cooling model.json --real-run
```
where `model.json` is a file in JSON format (described below) that is the input
for a cooling model and `--real-run` indicates that we actually want to run the
model as opposed to displaying diagnostic information.

## JSON model file format
An input file for the `cooling` utility is a JSON file with the following 
format
```
{
    "material": <choice of "magnetite", "iron">,
    "sizes": <JSON fragment (see below)>,
    "elongations": <JSON fragment (see below)>,
    "directions": <JSON fragment (see below)>,
    "applied_field": {
        "strength": <integer>,
        "direction": [<float>, <float>, <float>],
        "unit": <choice of: "T", "mT", "uT">
    },
    "cooling_regime": {
        "ambient_temperature": <float>,
        "initial_temperature": <float>,
        "reference_time": <float>,
        "temperature_reference_time": <float>,
        "allowable_percentage_drop": <float>,
        "stopping_temperature": 50
    },
    "outputs": {
        "model": <string>,
        "sizes": <optional string>,
        "elongations": <optional string>,
        "directions": <optional string>,
        "cooling_regime": <optional string>
    },
    "tau0": <optional float>,
    "epsilon": <optional float>,
    "n_polish": <optional float>
}
```

### JSON sizes fragment
The sizes fragment should be either a list or a distribution of individual 
sizes. The following fragment specifies a list where each size and the
proportion of grains that make up an assembly of that size are specified
```
{
    "list": [
        {"value": <float>, "fraction": <float>},
        ...
        {"value": <float>, "fraction": <float>}
    ],
    "unit": <choice of "m", "cm", "mm", "um", "nm">
}
```
where `list` denotes a list of sizes followed by a list of one or more
name-value pairs where `value` is the grain size and `fraction` is the
proportion of grains at that size. Internally the values given for `fraction`
are normalized so that they sum to unity.

Alternatively we may specify a lognormal distribution for sizes
```
{
    "distribution": {
        "type": "lognormal",
        "nbins": <integer>,
        "shape": <float>,
        "location": <float>,
        "scale": <float>,
        "start": <float>,
        "end": <float>
    },
    "unit": <choice of "m", "cm", "mm", "um", "nm">
}
```
where `nbins` is the number of grain sizes, `start` is the size range start
value and `end` is the size range value. The other parameters: `shape`,
`location` and `scale`, govern the log normal distribution from which values 
are drawn as shown in the equation

<img src="doc/images/lognormal/lognormal.png" alt="Lognormal" width="70%"/>

### JSON elongations fragment

The elongations fragment looks exactly like the sizes fragment described above
however the `unit` field is no longer needed.

### JSON directions fragment
The direction fragment should either be a list or a selection of directions
drawn from a function. The list fragment is very similar to the list fragments
for sizes and elongations, with the exception that the `value` field is now
a triple of three floating point values as opposed to a scalar and the `unit`
field is not required
```
{
    "list": [
        {"value": [<float>, <float>, <float>], "fraction": <float>},
        ...
        {"value": [<float>, <float>, <float>], "fraction": <float>}
    ]
}
```

Only one direction distribution function is currently supported called
*fibonacci*. This function will uniformly sample the sphere in `nbins`
directions and assign a fraction taken from a uniform distribution to each
grain orientation. The JSON fragment should look like
```
{
    "distribution": {
        "type": "fibonacci",
        "nbins": <int>
    }
}
```

### JSON cooling\_regime fragment
The `cooling_regime` fragment tells us how quickly a sample cools. It follows 
the Newtonian cooling rate equation

<img src="doc/images/newton-cooling/newton-cooling.png" alt="Newton cooling" width="70%"/>

* `ambient_temperature` is identified with T<sub>amb</sub>,
* `initial_temperature` is identified with T<sub>0</sub>,
* `reference_time` is identified with t<sub>1</sub>,
* `temperature_at_reference_time` is identified with T<sub>1</sub>,
* `allowable_percentage_drop` is the maximum allowed drop in temperature
  (as a fraction) for each time step - this quantity governs the time step
  at which we calculate energy barriers.
* `stopping_temperature` is an optional value which is a hard stop to the 
  calculation


### JSON outputs fragment
The `outputs` fragments is a list of files that are written when the 
model has completed execution.

* `model` is the name of the final model file in comma separated file format,
  the fields are: 1) time in seconds, 2) temperature 3)
  cooled normalized remanence 4) equilibrium remanence;
* `sizes` is the name of a file that contains a histogram of the sizes used in
  the model, the fields are: 1) size, 2) fraction;
* `elongations` is the name of a file that contains a histogram of the
  elongations used in the model, the fields are: 1) elongation, 2) fraction;
* `directions` is a file containing each direction and the fraction of grains
  occupying that direction, the fields are: 1) direction x component, 2)
  direction y component, 3) direction z component, 4) fraction;
* `cooling_regime` is a file containing the cooling curve, the fields are 1)
  time in seconds, 2) temperature.

### JSON other fields
The only other fields in a model JSON input are 

* `tau0` which is the magnetic reordering time usually taken to be 1E-9 or 
  1E-10 seconds,
* `epsilon` a value used to check whether floating point quantities in the code
  are sufficiently close to zero (this does not affect multiprecision routines).
* `n_polish` the number of Newton-Raphson polishing steps used when finding 
  roots to the Stoner-Wohlfarth equation (usually zero).

## JSON model file examples

The following two examples show possible uses of the JSON model file format.
The first example illustrates how the user can explicitly supply `sizes`,
`elongations` and `directions` fields with list fragments. The second example
illustrates the use of distributions for `sizes`, `elongations` and
`directions`. Of course it is possible to mix lists and distributions, for
example you can have a distribution of `direcions` and `sizes` and examine
specific elongations using a list.

### Example JSON - explicit lists
```json
{
  "material": "magnetite",

  "sizes": {
    "list": [{"value": 30, "fraction": 1.0}],
    "unit": "nm"
  },

  "elongations": {
    "list": [{"value": 60, "fraction": 1.0}]
  },

  "directions": {
    "list": [{"value": [1,0,0], "fraction": 0.2}]
  },

  "applied_field": {
    "strength": 30,
    "direction": [1,0,0],
    "unit": "uT"
  },

  "cooling_regime": {
    "ambient_temperature": 15.00,
    "initial_temperature": 579.999,
    "reference_time": 6E15,
    "temperature_at_reference_time": 15.15,
    "allowable_percentage_drop": 0.1,
    "stopping_temperature": 50
  },

  "outputs": {
    "model": "example_cooling2.csv",
    "sizes": "example_cooling_sizes2.csv",
    "elongations": "example_cooling_elongs2.csv",
    "directions": "example_directions2.csv",
    "cooling_regime": "example_cooling_regime2.csv"
  },

  "tau0": 1E-10,
  "epsilon": 1E-15,
  "n_polish": 0

}
```

### Example JSON - distributions

```json
{
  "material": "magnetite",
  "sizes": {
    "distribution": {
      "type": "lognormal",
      "nbins": 30,
      "shape": 1.0,
      "location": 20.0,
      "scale": 10.0,
      "start": 20.0,
      "end": 70
    },
    "unit": "nm"
  },

  "elongations": {
    "distribution": {
      "type": "lognormal",
      "nbins": 30,
      "shape": 1.0,
      "location": 5,
      "scale": 10.0,
      "start": 5,
      "end": 100
    }
  },

  "directions": {
    "distribution": {
      "type": "fibonacci",
      "nbins": 30
    }
  },

  "applied_field": {
    "strength": 30,
    "direction": [1,0,0],
    "unit": "uT"
  },

  "cooling_regime": {
    "ambient_temperature": 15.00,
    "initial_temperature": 579.999,
    "reference_time": 6E1,
    "temperature_at_reference_time": 15.15,
    "allowable_percentage_drop": 0.1,
    "stopping_temperature": 50
  },

  "outputs": {
    "model": "example_cooling2.csv",
    "sizes": "cooling_sizes2.csv",
    "elongations": "example_cooling_elongs2.csv",
    "directions": "example_directions2.csv",
    "cooling_regime": "example_cooling_regime2.csv"
  },

  "tau0": 1E-10,
  "epsilon": 1E-15,
  "n_polish": 0

}
```

# Credits
This software is written and developed by

* Les Nagy     (l1nagy@ucsd.edu)        
    - Scripps Institution of Oceanography, La Jolla CA, USA
* Wyn Williams (wyn.williams@ed.ac.uk)  
    - School of GeoSciences, University of Edinburgh, Edinburgh, Scotland
* Lisa Tauxe   (ltauxe@ucsd.edu)       
    - Scripps Institution of Oceanography, La Jolla CA, USA

## Third party software
The software was made possible because of the contribution of several open
source software projects. The following projects are included in the
`third_party` directory.

### Eigen
The Eigen library is maintained on the web page [https://eigen.tuxfamily.org](https://eigen.tuxfamily.org).
A complete set of contributors to Eigen can be found on the [credits
page](https://eigen.tuxfamily.org/index.php?title=Main_Page#Credits).

### JSON
The JSON library is developed by Niels Lohmann. It can be found on
([https://github.com/nlohmann/json](https://github.com/nlohmann/json)).
        
