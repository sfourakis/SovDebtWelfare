# Julia 1.2.0 Code

include("LongTermDebt_Methods.jl")
include("LongTermDebt_PoliciesMethods.jl")
include("LongTermDebt_WelfareMethods.jl")

################################################################################
# Initializing values that are constant across all simulations

# asset grid
aBounds = [-1.0, 0.0]
aPoints = 200  # 350

# output process
yPoints = 100 # 200
rho = 0.948503
eta2 = 0.027092^2
mu = 0.0
stdSpan = 3.0
inflateEndpoints = false

# m shock parameters
mPoints = 12
epsilon2 = 0.003^2
mMu = 0.0
mStdSpan = 2.0

# output and m structs
yParams = ar1Params(yPoints, rho, eta2, mu, stdSpan, inflateEndpoints)
mParams = iidParams(mPoints, epsilon2, mMu, mStdSpan)


################################################################################
# Chatterjee and Eyigungor (2012)
# The following are exactly the parameters of Chatterjee and Eyigungor (2012)

println("")
println("Computing CE 2012")
println("")

# model parameters
betaCE = 0.9540232420
thetaCE = 0.0385
gammaCE = 2.0
hPen0CE = -0.1881927550
hPen1CE = 0.2455843389
RCE = 1.01
lambdaCE = 0.05
coupCE = 0.03

simplePenCE = false
mixFacQCE = 0.5
mixFacVCE = 0.5

betaCBoundsCE = [0.9, 0.999]
betaCPointsCE = 20   # 199
gammaBoundsCE = [2.0, 2.0]
gammaPointsCE = 1
penMultCBoundsCE = [0.0, 1.0]
penMultCPointsCE = 2

# computations
LTBSpecCE = longTermBondSpec(
    betaCE,
    thetaCE,
    gammaCE,
    hPen0CE,
    hPen1CE,
    RCE,
    lambdaCE,
    coupCE,
    aBounds,
    aPoints,
    yParams,
    mParams,
    simplePenCE,
    mixFacQCE,
    mixFacVCE,
)

LTBEvalCE = makeEval(LTBSpecCE)
GC.gc()
@time vfiGOneStep!(LTBSpecCE, LTBEvalCE, 1e-10, 2000)

@time LTBPolCE = makeMeanPolicies(LTBSpecCE, LTBEvalCE)
@time LTBPathsCE, APCheckCE, APDivYCheckCE =
    simulatePaths(LTBSpecCE, LTBEvalCE, LTBPolCE, 0.0, 1000, 1000)

welfareRangeParamsCE = welfareRangeParams(
    betaCBoundsCE,
    betaCPointsCE,
    gammaBoundsCE,
    gammaPointsCE,
    penMultCBoundsCE,
    penMultCPointsCE,
    true,
    true,
    true,
    true,
)

@time welfareResAllCE =
    getWelfareRange(LTBSpecCE, LTBEvalCE, welfareRangeParamsCE, true, 1e-10, 100000)

welfareResPenUseGCE, welfareResPenUseCCE, welfareResNoPenUseGCE, welfareResNoPenUseCCE =
    splitWelfareResults(welfareResAllCE)

decompResultsCE = decomposeWelfareBeta(
    LTBSpecCE,
    LTBEvalCE,
    welfareResPenUseCCE,
    welfareResNoPenUseCCE,
    LTBPathsCE,
    0.0,
    1000,
)

wGNoPenSpecCE = debtWelfareSpec(betaCE, gammaCE, 0.0, false, false, betaCE)

wGNoPenEvalCE = makeWelfareEval(LTBSpecCE, wGNoPenSpecCE, LTBEvalCE)

estimateWelfare!(LTBSpecCE, wGNoPenSpecCE, LTBEvalCE, wGNoPenEvalCE, 1e-10, 5000)

lambdaDDecompResultsCE = decomposeLambdaD(
    LTBSpecCE,
    LTBEvalCE,
    welfareResPenUseCCE,
    welfareResNoPenUseCCE,
    welfareResPenUseGCE,
    welfareResNoPenUseGCE,
    wGNoPenEvalCE,
    0.0,
    1000,
)


################################################################################
# Arellano (2008)
#
# The following are exactly the parameters of Arellano (2008)
# except for the m shock, which is set to be identical to the one in Chatterjee
# and Eyigungor (2012)

println("")
println("Computing AR 2008")
println("")

# model parameters
betaAr = 0.95282
thetaAr = 0.282
gammaAr = 2.0
hPen0Ar = 0.969
hPen1Ar = 0.0
RAr = 1.017
lambdaAr = 1.0
coupAr = 0.03

simplePenAr = true
mixFacQAr = 0.5
mixFacVAr = 0.5

betaCBoundsAr = [0.9, 0.999]
betaCPointsAr = 20 # 199
gammaBoundsAr = [2.0, 2.0]
gammaPointsAr = 1
penMultCBoundsAr = [0.0, 1.0]
penMultCPointsAr = 2

LTBSpecAr = longTermBondSpec(
    betaAr,
    thetaAr,
    gammaAr,
    hPen0Ar,
    hPen1Ar,
    RAr,
    lambdaAr,
    coupAr,
    aBounds,
    aPoints,
    yParams,
    mParams,
    simplePenAr,
    mixFacQAr,
    mixFacVAr,
)

# computations
LTBEvalAr = makeEval(LTBSpecAr)
GC.gc()
@time vfiGOneStep!(LTBSpecAr, LTBEvalAr, 1e-10, 2000)

@time LTBPolAr = makeMeanPolicies(LTBSpecAr, LTBEvalAr)
@time LTBPathsAr, APCheckAr, APDivYCheckAr =
    simulatePaths(LTBSpecAr, LTBEvalAr, LTBPolAr, 0.0, 1000, 1000)

welfareRangeParamsAr = welfareRangeParams(
    betaCBoundsAr,
    betaCPointsAr,
    gammaBoundsAr,
    gammaPointsAr,
    penMultCBoundsAr,
    penMultCPointsAr,
    true,
    true,
    false,
    true,
)

@time welfareResAllAr =
    getWelfareRange(LTBSpecAr, LTBEvalAr, welfareRangeParamsAr, true, 1e-10, 100000)

welfareResPenUseCAr, welfareResNoPenUseCAr = splitWelfareResults(welfareResAllAr)

decompResultsAr = decomposeWelfareBeta(
    LTBSpecAr,
    LTBEvalAr,
    welfareResPenUseCAr,
    welfareResNoPenUseCAr,
    LTBPathsAr,
    0.0,
    1000,
)


################################################################################
# Aguiar-Gopinath (2006)
#
# The following are exactly the parameters of Aguiar-Gopinath (2008)
# except for the m shock, which is set to be identical to the one in Chatterjee
# and Eyigungor (2012)

println("")
println("Computing AG 2006")
println("")

# model parameters
betaAG = 0.8
thetaAG = 0.1
gammaAG = 2.0
hPen0AG = 0.02
hPen1AG = 0.0
RAG = 1.01
lambdaAG = 1.0
coupAG = 0.03

simplePenAG = false
mixFacQAG = 0.5
mixFacVAG = 0.5

betaCBoundsAG = [0.7, 0.999]
betaCPointsAG = 20 # 300
gammaBoundsAG = [2.0, 2.0]
gammaPointsAG = 1
penMultCBoundsAG = [0.0, 1.0]
penMultCPointsAG = 2

# computations
LTBSpecAG = longTermBondSpec(
    betaAG,
    thetaAG,
    gammaAG,
    hPen0AG,
    hPen1AG,
    RAG,
    lambdaAG,
    coupAG,
    aBounds,
    aPoints,
    yParams,
    mParams,
    simplePenAG,
    mixFacQAG,
    mixFacVAG,
)

LTBEvalAG = makeEval(LTBSpecAG)
GC.gc()
@time vfiGOneStep!(LTBSpecAG, LTBEvalAG, 1e-10, 2000)

@time LTBPolAG = makeMeanPolicies(LTBSpecAG, LTBEvalAG)
@time LTBPathsAG, APCheckAG, APDivYCheckAG =
    simulatePaths(LTBSpecAG, LTBEvalAG, LTBPolAG, 0.0, 1000, 1000)

welfareRangeParamsAG = welfareRangeParams(
    betaCBoundsAG,
    betaCPointsAG,
    gammaBoundsAG,
    gammaPointsAG,
    penMultCBoundsAG,
    penMultCPointsAG,
    true,
    true,
    false,
    true,
)

@time welfareResAllAG =
    getWelfareRange(LTBSpecAG, LTBEvalAG, welfareRangeParamsAG, true, 1e-10, 100000)

welfareResPenUseCAG, welfareResNoPenUseCAG = splitWelfareResults(welfareResAllAG)

decompResultsAG = decomposeWelfareBeta(
    LTBSpecAG,
    LTBEvalAG,
    welfareResPenUseCAG,
    welfareResNoPenUseCAG,
    LTBPathsAG,
    0.0,
    1000,
)

println("")
println("Savings the output..")

writecsv(
    joinpath("OUTPUT", "CSV", "CEBenchmark12_welfareGainDecomposition.csv"),
    decompResultsCE,
)
writecsv(
    joinpath("OUTPUT", "CSV", "ArBenchmark08_welfareGainDecomposition.csv"),
    decompResultsAr,
)
writecsv(
    joinpath("OUTPUT", "CSV", "AGTransitory06_welfareGainDecomposition.csv"),
    decompResultsAG,
)
writecsv(
    joinpath("OUTPUT", "CSV", "CEBenchmark12_lambdaDDecomposition.csv"),
    lambdaDDecompResultsCE,
)
writecsv(
    joinpath("OUTPUT", "CSV", "apDivYNDLimAll.csv"),
    [LTBPathsCE.assets.apDivY.varNDPath[end],LTBPathsAr.assets.apDivY.varNDPath[end],LTBPathsAG.assets.apDivY.varNDPath[end]],
)

println("Done.")
