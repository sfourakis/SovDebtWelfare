# SovDebtWelfare

This repository contains codes for:

1. solving the model of Chatterjee and Eyigungor (2012);
2. generating the ergodic joint distribution of debt, income, and default in that model;
3. generating the mean policy functions and mean paths of equilibrium variables starting at 0 debt and the ergodic distribution of income;
4. computing the implied welfare variation (relative to autarky) for a consumer who has different preferences (discounting and/or risk aversion) from the government and/or experiences the output losses from default differently from the government.
5. decomposing the results of 4. into a term measuring variation due purely to default costs, a term measuring variation due changes in the volatility of consumption, and a term measuring variation due purely to the government's choice of time path for consumption.
6. decomposing the term in 5. which measures variation due purely to default costs into a term measuring variation related to the consumer and the government valuing losses after a default differently and a term measuring variation related to the consumer valuing, at time 0, the losses entailed by a default that the government experiences differently.

The files are:
1. LongTermDebt_Methods.jl: contains model structures, setup, and solution methods;
2. LongTermDebt_PoliciesMethods.jl: contains structures, setup, and computation methods for 3. above;
3. LongTermDebt_WelfareMethods.jl: contains structures, setup, and computation methods for 2., 4., 5., and 6. above;
4. LongTermDebtRunAll.jl: solves the models, produces related objects, and produces and writes to disk the welfare decompositions used in the paper.

Comments in files 2-4 are sparse, may be inaccurate/obsolete, and are subject to change.

These codes should run without error in any version Julia newer than v1.0. Results used in the paper were computed in Julia v1.2.0.
