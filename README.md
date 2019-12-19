# SovDebtWelfare

December 2019 

This repository contains the code that generates the results used for the tables
and figures of the simulations in the paper:

_"On the Welfare Losses from External Sovereign Borrowing"_, by Mark Aguiar,
Manuel Amador and Stelios Fourakis 

## The code 

The code contained in this repository does the following: 

1. Solves the model of Chatterjee and Eyigungor (2012).

2. Generates the ergodic joint distribution of debt, income, and default in that
   model.

3. Generates the mean policy functions and mean paths of equilibrium variables
   starting at 0 debt and the ergodic distribution of income.

4. Computes the implied welfare variation (relative to autarky) for a consumer
   who has different preferences (discounting and/or risk aversion) from the
   government and/or experiences the output losses from default differently from
   the government.

5. Decomposes the results of 4. into a term measuring variation due purely to
   default costs, a term measuring variation due to changes in the volatility of
   consumption, and a term measuring variation due purely to the government's
   choice of time path for consumption.

6. Decomposes the term in 5. which measures variation due purely to default
   costs into a term measuring variation related to the consumer and the
   government valuing losses after a default differently and a term measuring
   variation related to the consumer valuing, at time 0, the losses entailed by
   a default that the government experiences differently.

7. Generates the figures. 

The code should run without error in any version Julia newer than v1.0.
Results used in the paper were computed in Julia v1.2.0.

## Description of the files

1. *LongTermDebt_Methods.jl*: contains model structures, setup, and solution
   methods.

2. *LongTermDebt_PoliciesMethods.jl*: contains structures, setup, and
   computation methods for 3. above.

3. *LongTermDebt_WelfareMethods.jl*: contains structures, setup, and computation
   methods for 2., 4., 5., and 6. above.

4. *LongTermDebtRunAll.jl*: solves the models, produces related objects, and
   produces and writes to disk (OUTPUT/CSV) the welfare decompositions used in
   the paper.

5. *generate_figures.ipynb*: Jupyter notebook that uses the generated csv files
   to make the figures. (uses Python)

## Outputs

File 4 generates the following output files in the OUTPUT/CSV directory:

1. CEBenchmark12_welfareGainDecomposition.csv
2. ArBenchmark08_welfareGainDecomposition.csv
3. AGTransitory06_welfareGainDecomposition.csv
4. CEBenchmark12_lambdaDDecomposition.csv
5. asset_over_Y_limit_ND.csv

The columns of the first 3 files contain:

1. beta=the consumer's discount factor;
2. (1+lambda)=the total consumption equivalent welfare variation of openness
   relative to autarky;
3. (1+lambda_D)=the part of (1+lambda) due purely to default costs;
4. (1+lambda_V)=the part of (1+lambda) due purely to changes in the volatility
   of consumption;
5. (1+lambda_T)=the part of (1+lambda) due purely to the government's choice of
   time path for consumption.

The columns of the fourth file contains:

1. beta=the consumer's discount factor;
2. (1+lambdaH_D)=(E\[V_H(s,0)]/E\[V_H^ND(s,0)])^(1/(1-gamma)) (=(1_lambda_D)
   above)
3. (1+lambdaHHat_D)=(E\[VHat_H(s,0)]/E\[VHat_H^ND(s,0)])^(1/(1-gamma))
4. (1+lambdaG_D)=(E\[V_G(s,0)]/E\[V_G^ND(s,0)])^(1/(1-gamma)) (invariant to
   consumer beta)
5. (1+lambdaH_D)/(1+lambdaHHat_D)
6. (1+lambdaHHat_D)/(1+lambdaG_D)

where V_G (V_G^ND) denotes the value function computed using the government's
discount factor with default costs (with default costs); where V_H (V_H^ND)
denotes the value function computed using the consumer's discount factor with
default costs (with default costs); and VHat_H (VHat_H^ND) denotes the value
function computed using the consumer's discount factor until the first period of
default and the government's discount factor in that period and thereafter with
default costs (with default costs).

The rows of the last file contains the asset limits conditional on no default 
for the CE2012, Ar2008, and AG2006 models respectively. This is used to generate
the benchmark comparison in the figures. 

Comments in files 2-4 are sparse, may be inaccurate/obsolete, and are subject to
change.
