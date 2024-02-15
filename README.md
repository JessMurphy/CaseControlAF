# CaseControlAF
Case Control AF Reconstruction R Package

This repository contains the source code for the CaseControlAF R package which can be used to reconstruct the AF for cases and controls separately given commonly available summary statistics. 

The package contains two functions:

1) CaseControl_AF
2) CaseControl_SE

## CaseControl_AF
Use this function when you have the following statistics (for each variant)
* Number of cases
* Number of controls
* OR or beta
* **AF** for the whole sample (cases and controls combined)

### Usage
**N_case**: an integer for the number of case samples
**N_control**: an integer for the number of control samples
**OR**: a vector of doubles with the OR for each variant
**AF_population**: a vector of doubles with the AF for each variant

**Return**: Returns a dataframe with two columns with names: AF_case and AF_control. The number of rows is equal to the number of variants

```
nCase = 1000
nControl = 5000
af = c(.82, .43, .22, .40, .04)
ors = c(1.1, 1.22, 1.01, 1.05, 0.95)

CaseControl_AF(N_case = nCase, N_control = nControl, AF_population = af, OR = ors)
```

## CaseControl_SE
Use this function when you have the following statistics (for each variant)
* Number of cases
* Number of controls
* OR or beta
* **SE** of the log(OR) for each variant

*Code adapted from ReACt GroupFreq function available here: [(https://github.com/Paschou-Lab/ReAct/blob/main/GrpPRS_src/CountConstruct.c)]*

### Usage
**N_case**: an integer for the number of case samples
**N_control**: an integer for the number of control samples
**OR**: a vector of doubles with the OR for each variant
**SE**: a vector of doubles with the SE(log(OR)) for each variant

**Return**: Returns a dataframe with three columns with names: pCase, pControl and pPop. The number of rows is equal to the number of variants

```
nCase = 1000
nControl = 5000
af = c(.82, .43, .22, .40, .04)
se = c(.0025, .004, .0054, .002, .051)

CaseControl_AF(N_case = nCase, N_control = nControl, AF_population = af, SE = se)
```
