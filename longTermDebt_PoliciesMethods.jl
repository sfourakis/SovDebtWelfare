
struct consumptionPolicy{F<:Real}
    consRGrid::Array{F,2}
    consGrid::Array{F,2}
    consNoPen::Array{F,2}
    consNoPenD::Array{F,1}
    consRVar::Array{F,2}
    consVar::Array{F,2}
    consDVar::Array{F,1}
end

struct consumptionDivYPolicy{F<:Real}
    consRDivYGrid::Array{F,2}
    consDivYGrid::Array{F,2}
    consDivYDInitGrid::Array{F,1}
    consDivYDFutGrid::Array{F,1}
    consRDivYVar::Array{F,2}
    consDivYVar::Array{F,2}
    consDivYDInitVar::Array{F,1}
    consDivYDFutVar::Array{F,1}
end

struct utilityPolicy{F<:Real}
    uGrid::Array{F,2}
    uRGrid::Array{F,2}
    uDInitGrid::Array{F,1}
    uDFutGrid::Array{F,1}
    uNoPen::Array{F,2}
    uNoPenD::Array{F,1}
    utYObs::Array{F,2}
    utYObsD::Array{F,1}
end

struct outputPolicy{F<:Real}
    yObs::Array{F,2}
    yObsD::Array{F,1}
end

struct aPrimeDefPolicy{F<:Real,S<:Integer}
    apRGrid::Array{F,2}
    apRDivYGrid::Array{F,2}
    defaultProb::Array{F,2}
    repayProb::Array{F,2}
    apProbability::Array{SparseMatrixCSC{F,S},1}
    lastAlwaysDefInd::Array{S,1}
end


struct meanPolicies{F<:Real,S<:Integer}
    cPol::consumptionPolicy{F}
    cDivYPol::consumptionDivYPolicy{F}
    apDefPol::aPrimeDefPolicy{F,S}
    uPol::utilityPolicy{F}
    yPol::outputPolicy{F}
end


struct singlePath{F<:Real}
    varPath::Array{F,1}
    varDistPath::Array{F,2}
    varNDPath::Array{F,1}
    varNDDistPath::Array{F,2}
end

struct consumptionPaths{F<:Real}
    cY::singlePath{F}
    cRY::singlePath{F}
    cA::singlePath{F}
    cRA::singlePath{F}
    cCert::singlePath{F}
    cYNoPen::singlePath{F}
    cYNoPenCert::singlePath{F}
    cYVar::singlePath{F}
end

struct consumptionDivYPaths{F<:Real}
    cDivY::singlePath{F}
    cRDivY::singlePath{F}
    cDivYVar::singlePath{F}
end

struct assetPaths{F<:Real}
    a::singlePath{F}
    aR::singlePath{F}
    ap::singlePath{F}
    apDivY::singlePath{F}
end

struct outputPaths{F<:Real}
    y::singlePath{F}
    yR::singlePath{F}
    outputWithPen::singlePath{F}
    outputCert::singlePath{F}
end

struct utilityPaths{F<:Real}
    ut::singlePath{F}
    utNoPen::singlePath{F}
    utOutput::singlePath{F}
end

struct meanPaths{F<:Real}
    consVal::consumptionPaths{F}
    consDivY::consumptionDivYPaths{F}
    assets::assetPaths{F}
    yVal::outputPaths{F}
    uVal::utilityPaths{F}
    RMassPath::Array{F,1}
    DMassPath::Array{F,1}
    NDMassPath::Array{F,1}
end

function initializePath(m::longTermBondSpec{F,S},distDim::S,bigT::S) where{F<:Real,S<:Integer}


    varPath=zeros(F,bigT)
    varDistPath=zeros(F,distDim,bigT)
    varNDPath=zeros(F,bigT)
    varNDDistPath=zeros(F,distDim,bigT)

    return singlePath(varPath,varDistPath,varNDPath,varNDDistPath)
end

function makeConsPaths(m::longTermBondSpec{F,S},bigT::S) where{F<:Real,S<:Integer}
    aDimPathDummy=initializePath(m,m.aPoints,bigT)
    yDimPathDummy=initializePath(m,m.yParams.nPoints,bigT)
    cY=deepcopy(yDimPathDummy)
    cRY=deepcopy(yDimPathDummy)
    cA=deepcopy(aDimPathDummy)
    cRA=deepcopy(aDimPathDummy)
    cCert=deepcopy(yDimPathDummy)
    cYNoPen=deepcopy(yDimPathDummy)
    cYNoPenCert=deepcopy(yDimPathDummy)
    cYVar=deepcopy(yDimPathDummy)
    return consumptionPaths(cY,cRY,cA,cRA,cCert,cYNoPen,cYNoPenCert,cYVar)
end


function makeConsDivYPaths(m::longTermBondSpec{F,S},bigT::S) where{F<:Real,S<:Integer}
    yDimPathDummy=initializePath(m,m.yParams.nPoints,bigT)

    cDivY=deepcopy(yDimPathDummy)
    cRDivY=deepcopy(yDimPathDummy)
    cDivYVar=deepcopy(yDimPathDummy)


    return consumptionDivYPaths(cDivY,cRDivY,cDivYVar)
end


function makeAssetPaths(m::longTermBondSpec{F,S},bigT::S) where{F<:Real,S<:Integer}
    yDimPathDummy=initializePath(m,m.yParams.nPoints,bigT)

    a=deepcopy(yDimPathDummy)
    aR=deepcopy(yDimPathDummy)
    ap=deepcopy(yDimPathDummy)
    apDivY=deepcopy(yDimPathDummy)

    return assetPaths(a,aR,ap,apDivY)

end

function makeOutputPaths(m::longTermBondSpec{F,S},bigT::S) where{F<:Real,S<:Integer}
    aDimPathDummy=initializePath(m,m.aPoints,bigT)
    yDimPathDummy=initializePath(m,m.yParams.nPoints,bigT)
    y=deepcopy(aDimPathDummy)
    yR=deepcopy(aDimPathDummy)
    outputWithPen=deepcopy(yDimPathDummy)
    outputCert=deepcopy(yDimPathDummy)

    return outputPaths(y,yR,outputWithPen,outputCert)
end

function makeUtilityPaths(m::longTermBondSpec{F,S},bigT::S) where{F<:Real,S<:Integer}
    aDimPathDummy=initializePath(m,m.aPoints,bigT)
    yDimPathDummy=initializePath(m,m.yParams.nPoints,bigT)
    ut=deepcopy(yDimPathDummy)
    utNoPen=deepcopy(yDimPathDummy)
    utOutput=deepcopy(yDimPathDummy)

    return utilityPaths(ut,utNoPen,utOutput)

end


function makeMeanPaths(m::longTermBondSpec{F,S},bigT::S) where{F<:Real,S<:Integer}
    consVal=makeConsPaths(m,bigT)
    consDivY=makeConsDivYPaths(m,bigT)
    assets=makeAssetPaths(m,bigT)
    yVal=makeOutputPaths(m,bigT)
    uVal=makeUtilityPaths(m,bigT)
    RMassPath=zeros(F,bigT)
    DMassPath=zeros(F,bigT)
    NDMassPath=zeros(F,bigT)


    return meanPaths(consVal,consDivY,assets,yVal,uVal,RMassPath,DMassPath,NDMassPath)
end

function makeConsPolicy(m::longTermBondSpec{F,S}) where{F<:Real,S<:Integer}
    consRGrid=zeros(F,m.aPoints,m.yParams.nPoints)
    consGrid=zeros(F,m.aPoints,m.yParams.nPoints)
    consNoPen=zeros(F,m.aPoints,m.yParams.nPoints)
    consNoPenD=zeros(F,m.yParams.nPoints)
    consRVar=zeros(F,m.aPoints,m.yParams.nPoints)
    consVar=zeros(F,m.aPoints,m.yParams.nPoints)
    consDVar=zeros(F,m.yParams.nPoints)

    return consumptionPolicy(consRGrid,consGrid,consNoPen,consNoPenD,consRVar,consVar,consDVar)

end


function makeConsDivYPolicy(m::longTermBondSpec{F,S}) where{F<:Real,S<:Integer}
    consDivYGrid=zeros(F,m.aPoints,m.yParams.nPoints)
    consRDivYGrid=zeros(F,m.aPoints,m.yParams.nPoints)
    consDivYDInitGrid=zeros(F,m.yParams.nPoints)
    consDivYDFutGrid=zeros(F,m.yParams.nPoints)
    consRDivYVar=zeros(F,m.aPoints,m.yParams.nPoints)
    consDivYVar=zeros(F,m.aPoints,m.yParams.nPoints)
    consDivYDFutVar=zeros(F,m.yParams.nPoints)
    consDivYDInitVar=zeros(F,m.yParams.nPoints)

    return consumptionDivYPolicy(consDivYGrid,consRDivYGrid,consDivYDInitGrid,consDivYDFutGrid,consRDivYVar,consDivYVar,consDivYDFutVar,consDivYDInitVar)

end
function makeUtilityPolicy(m::longTermBondSpec{F,S}) where{F<:Real,S<:Integer}
    uGrid=zeros(F,m.aPoints,m.yParams.nPoints)
    uRGrid=zeros(F,m.aPoints,m.yParams.nPoints)
    uDInitGrid=zeros(F,m.yParams.nPoints)
    uDFutGrid=zeros(F,m.yParams.nPoints)
    uNoPen=zeros(F,m.aPoints,m.yParams.nPoints)
    uNoPenD=zeros(F,m.yParams.nPoints)
    utYObs=zeros(F,m.aPoints,m.yParams.nPoints)
    utYObsD=zeros(F,m.yParams.nPoints)

    return utilityPolicy(uGrid,uRGrid,uDInitGrid,uDFutGrid,uNoPen,uNoPenD,utYObs,utYObsD)
end

function makeOutputPolicy(m::longTermBondSpec{F,S}) where{F<:Real,S<:Integer}
    yObs=zeros(F,m.aPoints,m.yParams.nPoints)
    yObsD=zeros(F,m.yParams.nPoints)

    return outputPolicy(yObs,yObsD)
end

function makeAPDefPolicy(m::longTermBondSpec{F,S}) where{F<:Real,S<:Integer}
    apRGrid=zeros(F,m.aPoints,m.yParams.nPoints)
    apRDivYGrid=zeros(F,m.aPoints,m.yParams.nPoints)
    defaultProb=zeros(F,m.aPoints,m.yParams.nPoints)
    repayProb=zeros(F,m.aPoints,m.yParams.nPoints)
    apProbability=Array{SparseMatrixCSC{F,S},1}(undef,m.yParams.nPoints)
    lastAlwaysDefInd=zeros(S,m.yParams.nPoints)

    SDummyF=sparse(zeros(F,m.aPoints,m.aPoints))
    for i in 1:(m.yParams.nPoints)
        apProbability[i]=deepcopy(SDummyF)
    end

    return aPrimeDefPolicy(apRGrid,apRDivYGrid,defaultProb,repayProb,apProbability,lastAlwaysDefInd)

end


function makeMeanPolicies(m::longTermBondSpec{F,S},s::longTermBondEval{F,S}) where{F<:Real,S<:Integer}


    cPol=makeConsPolicy(m)
    cDivYPol=makeConsDivYPolicy(m)
    apDefPol=makeAPDefPolicy(m)
    uPol=makeUtilityPolicy(m)
    yPol=makeOutputPolicy(m)

    copy!(uPol.uDInitGrid,s.VF.vDInitFlow)
    copy!(uPol.uDFutGrid,s.VF.vDFutFlow)

    copy!(cPol.consNoPenD,s.income.yGrid)
    copy!(yPol.yObsD,s.income.yDefGrid)
    copy!(uPol.utYObsD,s.VF.vDFutFlow)


    for i in 1:m.yParams.nPoints
        tempADInd=0
        while tempADInd<m.aPoints
            if s.pol.alwaysDefault[tempADInd+1,i]==false
                break
            else
                tempADInd+=1
            end
        end
        apDefPol.lastAlwaysDefInd[i]=tempADInd
    end



    for i in 1:m.yParams.nPoints
        tempPSum=0.0
        for mInd in 1:(m.mParams.nPoints-1)
            cDivYPol.consDivYDInitGrid[i]+=s.income.mProb[mInd]*(s.income.yDefGrid[i]+s.income.mBounds[1])/(s.income.yGrid[i]+s.income.mMidPoints[mInd])
            cDivYPol.consDivYDFutGrid[i]+=s.income.mProb[mInd]*(s.income.yDefGrid[i]+s.income.mMidPoints[mInd])/(s.income.yGrid[i]+s.income.mMidPoints[mInd])
            uPol.uNoPenD[i]+=s.income.mProb[mInd]*u(m,s.income.yGrid[i]+s.income.mMidPoints[mInd])
            tempPSum+=s.income.mProb[mInd]
        end
        tempPSumInv=tempPSum^(-1)
        cDivYPol.consDivYDInitGrid[i]*=tempPSumInv
        cDivYPol.consDivYDFutGrid[i]*=tempPSumInv
        uPol.uNoPenD[i]*=tempPSumInv
    end



    for i in 1:m.yParams.nPoints
        tempPSum=0.0
        for mInd in 1:(m.mParams.nPoints-1)
            cPol.consDVar[i]+=s.income.mProb[mInd]*(s.income.yDefGrid[i]+s.income.mMidPoints[mInd]-s.income.yDefGrid[i])^2
            cDivYPol.consDivYDFutVar[i]+=s.income.mProb[mInd]*((s.income.yDefGrid[i]+s.income.mMidPoints[mInd])/(s.income.yGrid[i]+s.income.mMidPoints[mInd])-cDivYPol.consDivYDFutGrid[i])^2
            cDivYPol.consDivYDInitVar[i]+=s.income.mProb[mInd]*((s.income.yDefGrid[i]+s.income.mBounds[1])/(s.income.yGrid[i]+s.income.mMidPoints[mInd])-cDivYPol.consDivYDInitGrid[i])^2
            tempPSum+=s.income.mProb[mInd]
        end
        tempPSumInv=tempPSum^(-1)
        cPol.consDVar[i]*=tempPSumInv
        cDivYPol.consDivYDFutVar[i]*=tempPSumInv
        cDivYPol.consDivYDInitVar[i]*=tempPSumInv

    end




    for i in 1:(m.yParams.nPoints)
        for j in 1:(m.aPoints)
            tempPSum=0.0
            if s.pol.neverDefault[j,i]==true
                mUBInd=2
                aPolicyInd=1
                lastMVal=s.income.mGrid[1]
                while (mUBInd<=(m.mParams.nPoints))&&(aPolicyInd<=(s.pol.mListLength[j,i]))
                    if (s.pol.mAPThreshold[i][aPolicyInd,j])<(s.income.mGrid[mUBInd])
                        apDefPol.apProbability[i][aPolicyInd,j]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes
                        cPol.consRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])
                        apDefPol.apRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*s.aGrid[s.pol.apPolicy[i][aPolicyInd,j]]
                        apDefPol.apRDivYGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*s.aGrid[s.pol.apPolicy[i][aPolicyInd,j]]/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        cDivYPol.consRDivYGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        uPol.uRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*u(m,c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1]))

                        tempPSum+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes
                        lastMVal=s.pol.mAPThreshold[i][aPolicyInd,j]
                        aPolicyInd+=1
                    else
                        apDefPol.apProbability[i][aPolicyInd,j]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes
                        cPol.consRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])
                        apDefPol.apRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*s.aGrid[s.pol.apPolicy[i][aPolicyInd,j]]
                        apDefPol.apRDivYGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*s.aGrid[s.pol.apPolicy[i][aPolicyInd,j]]/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        cDivYPol.consRDivYGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        uPol.uRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*u(m,c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1]))

                        tempPSum+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes
                        lastMVal=s.income.mGrid[mUBInd]
                        mUBInd+=1
                    end
                end
                tempPSumInv=tempPSum^(-1)
                cPol.consRGrid[j,i]*=tempPSumInv
                apDefPol.apRGrid[j,i]*=tempPSumInv
                apDefPol.apRDivYGrid[j,i]*=tempPSumInv
                cDivYPol.consRDivYGrid[j,i]*=tempPSumInv
                uPol.uRGrid[j,i]*=tempPSumInv
                for apInd in 1:(s.pol.mListLength[j,i])
                    apDefPol.apProbability[i][apInd,j]*=tempPSumInv
                end
                cPol.consGrid[j,i]=cPol.consRGrid[j,i]
                cDivYPol.consDivYGrid[j,i]=cDivYPol.consRDivYGrid[j,i]
                uPol.uGrid[j,i]=uPol.uRGrid[j,i]

                uPol.utYObs[j,i]=uPol.uNoPenD[i]
                cPol.consNoPen[j,i]=cPol.consGrid[j,i]
                uPol.uNoPen[j,i]=uPol.uGrid[j,i]
                yPol.yObs[j,i]=s.income.yGrid[i]
            elseif s.pol.alwaysDefault[j,i]==false
                #If the government does default only sometimes, find the position in the m grid of the upper bound of the interval in which the threshold level of m for default falls
                #and set the last value of m observed to the value directly preceding that upper bound
                mUBInd=searchsortedfirst(s.income.mGrid,s.pol.defThreshold[j,i])
                lastMVal=s.income.mGrid[mUBInd-1]
                #If there are any entire intervals in which default occurs, add their contribution to the output value
                if mUBInd>2
                    for k in 3:mUBInd
                        apDefPol.defaultProb[j,i]+=s.income.mProb[k-2]
                        tempPSum+=s.income.mProb[k-2]
                        cDivYPol.consDivYGrid[j,i]+=s.income.mProb[k-2]*(s.income.yDefGrid[i]+s.income.mBounds[1])/(s.income.yGrid[i]+s.income.mMidPoints[k-2])

                        uPol.utYObs[j,i]+=s.income.mProb[k-2]*u(m,s.income.yDefGrid[i]+s.income.mBounds[1])
                        yPol.yObs[j,i]+=s.income.mProb[k-2]*(s.income.yDefGrid[i]+s.income.mBounds[1])
                        uPol.uNoPen[j,i]+=s.income.mProb[k-2]*u(m,s.income.yGrid[i]+s.income.mMidPoints[k-2])
                        cPol.consNoPen[j,i]+=s.income.mProb[k-2]*(s.income.yGrid[i]+s.income.mMidPoints[k-2])
                    end
                end
                #Set the index of the first relevant entry of the borowing policy function
                aPolicyInd=s.pol.firstRepayInd[j,i]

                #Add the contribution of default in the interval in which the threshold lies to the output value
                apDefPol.defaultProb[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes
                cDivYPol.consDivYGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes*(s.income.yDefGrid[i]+s.income.mBounds[1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])


                uPol.utYObs[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes*u(m,s.income.yDefGrid[i]+s.income.mBounds[1])
                yPol.yObs[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes*(s.income.yDefGrid[i]+s.income.mBounds[1])
                uPol.uNoPen[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes*u(m,s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                cPol.consNoPen[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes*(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])

                tempPSum+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes

                #Update the last value of m observed to the threshold level of m at which default occurs
                lastMVal=s.pol.defThreshold[j,i]

                while (mUBInd<=(m.mParams.nPoints))&&(aPolicyInd<=(s.pol.mListLength[j,i]))
                    if (s.pol.mAPThreshold[i][aPolicyInd,j])<(s.income.mGrid[mUBInd])
                        apDefPol.apProbability[i][aPolicyInd,j]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes
                        cPol.consRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])
                        apDefPol.apRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*s.aGrid[s.pol.apPolicy[i][aPolicyInd,j]]
                        apDefPol.apRDivYGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*s.aGrid[s.pol.apPolicy[i][aPolicyInd,j]]/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        cDivYPol.consRDivYGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        uPol.uRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*u(m,c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1]))

                        uPol.utYObs[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*u(m,s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        yPol.yObs[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        uPol.uNoPen[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*u(m,c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1]))
                        cPol.consNoPen[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])

                        tempPSum+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes
                        lastMVal=s.pol.mAPThreshold[i][aPolicyInd,j]
                        aPolicyInd+=1
                    else
                        apDefPol.apProbability[i][aPolicyInd,j]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes
                        cPol.consRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])
                        apDefPol.apRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*s.aGrid[s.pol.apPolicy[i][aPolicyInd,j]]
                        apDefPol.apRDivYGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*s.aGrid[s.pol.apPolicy[i][aPolicyInd,j]]/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        cDivYPol.consRDivYGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        uPol.uRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*u(m,c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1]))

                        uPol.utYObs[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*u(m,s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        yPol.yObs[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        uPol.uNoPen[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*u(m,c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1]))
                        cPol.consNoPen[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])


                        tempPSum+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes
                        lastMVal=s.income.mGrid[mUBInd]
                        mUBInd+=1
                    end
                end
                cDivYPol.consDivYGrid[j,i]+=cDivYPol.consRDivYGrid[j,i]

                tempPSumInv=tempPSum^(-1)
                apDefPol.defaultProb[j,i]*=tempPSumInv
                cDivYPol.consDivYGrid[j,i]*=tempPSumInv

                uPol.utYObs[j,i]*=tempPSumInv
                yPol.yObs[j,i]*=tempPSumInv
                uPol.uNoPen[j,i]*=tempPSumInv
                cPol.consNoPen[j,i]*=tempPSumInv

                cPol.consRGrid[j,i]*=(1.0-apDefPol.defaultProb[j,i])^(-1)
                apDefPol.apRGrid[j,i]*=(1.0-apDefPol.defaultProb[j,i])^(-1)
                apDefPol.apRDivYGrid[j,i]*=(1.0-apDefPol.defaultProb[j,i])^(-1)
                cDivYPol.consRDivYGrid[j,i]*=(1.0-apDefPol.defaultProb[j,i])^(-1)
                uPol.uRGrid[j,i]*=(1.0-apDefPol.defaultProb[j,i])^(-1)

                cPol.consGrid[j,i]=(1.0-apDefPol.defaultProb[j,i])*cPol.consRGrid[j,i]+apDefPol.defaultProb[j,i]*(s.income.yDefGrid[i]+s.income.mBounds[1])
                uPol.uGrid[j,i]=(1.0-apDefPol.defaultProb[j,i])*uPol.uRGrid[j,i]+apDefPol.defaultProb[j,i]*uPol.uDInitGrid[i]

                for apInd in (s.pol.firstRepayInd[j,i]):(s.pol.mListLength[j,i])
                    apDefPol.apProbability[i][apInd,j]*=tempPSumInv
                end
            else
                apDefPol.defaultProb[j,i]=1.0
                cPol.consGrid[j,i]=(s.income.yDefGrid[i]+s.income.mBounds[1])
                cDivYPol.consDivYGrid[j,i]=cDivYPol.consDivYDInitGrid[i]
                uPol.uGrid[j,i]=uPol.uDInitGrid[i]
                uPol.utYObs[j,i]=u(m,s.income.yDefGrid[i]+s.income.mBounds[1])
                yPol.yObs[j,i]=s.income.yDefGrid[i]+s.income.mBounds[1]
                uPol.uNoPen[j,i]=uPol.uNoPenD[i]
                cPol.consNoPen[j,i]=cPol.consNoPenD[i]
            end

        end
    end


    for i in 1:(m.yParams.nPoints)
        for j in 1:(m.aPoints)
            tempPSum=0.0
            if s.pol.neverDefault[j,i]==true
                mUBInd=2
                aPolicyInd=1
                lastMVal=s.income.mGrid[1]
                while (mUBInd<=(m.mParams.nPoints))&&(aPolicyInd<=(s.pol.mListLength[j,i]))
                    if (s.pol.mAPThreshold[i][aPolicyInd,j])<(s.income.mGrid[mUBInd])
                        cPol.consRVar[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*(c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])-cPol.consRGrid[j,i])^2
                        cDivYPol.consRDivYVar[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*(c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])-cDivYPol.consRDivYGrid[j,i])^2

                        tempPSum+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes
                        lastMVal=s.pol.mAPThreshold[i][aPolicyInd,j]
                        aPolicyInd+=1
                    else
                        cPol.consRVar[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*(c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])-cPol.consRGrid[j,i])^2
                        cDivYPol.consRDivYVar[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*(c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])-cDivYPol.consRDivYGrid[j,i])^2

                        tempPSum+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes
                        lastMVal=s.income.mGrid[mUBInd]
                        mUBInd+=1
                    end
                end
                tempPSumInv=tempPSum^(-1)
                cPol.consRVar[j,i]*=tempPSumInv
                cDivYPol.consRDivYVar[j,i]*=tempPSumInv

                cPol.consVar[j,i]=cPol.consRVar[j,i]
                cDivYPol.consDivYVar[j,i]=cDivYPol.consRDivYVar[j,i]

            elseif s.pol.alwaysDefault[j,i]==false
                #If the government does default only sometimes, find the position in the m grid of the upper bound of the interval in which the threshold level of m for default falls
                #and set the last value of m observed to the value directly preceding that upper bound
                mUBInd=searchsortedfirst(s.income.mGrid,s.pol.defThreshold[j,i])
                lastMVal=s.income.mGrid[mUBInd-1]
                #If there are any entire intervals in which default occurs, add their contribution to the output value
                if mUBInd>2
                    for k in 3:mUBInd
                        tempPSum+=s.income.mProb[k-2]
                        cPol.consVar[j,i]+=s.income.mProb[k-2]*(s.income.yDefGrid[i]+s.income.mBounds[1]-cPol.consGrid[j,i])^2
                        cDivYPol.consDivYVar[j,i]+=s.income.mProb[k-2]*((s.income.yDefGrid[i]+s.income.mBounds[1])/(s.income.yGrid[i]+s.income.mMidPoints[k-2])-cDivYPol.consDivYGrid[j,i])^2
                    end
                end
                #Set the index of the first relevant entry of the borowing policy function
                aPolicyInd=s.pol.firstRepayInd[j,i]

                #Add the contribution of default in the interval in which the threshold lies to the output value
                cPol.consVar[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes*(s.income.yDefGrid[i]+s.income.mBounds[1]-cPol.consGrid[j,i])^2
                cDivYPol.consDivYVar[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes*((s.income.yDefGrid[i]+s.income.mBounds[1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])-cDivYPol.consDivYGrid[j,i])^2

                tempPSum+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes

                #Update the last value of m observed to the threshold level of m at which default occurs
                lastMVal=s.pol.defThreshold[j,i]

                while (mUBInd<=(m.mParams.nPoints))&&(aPolicyInd<=(s.pol.mListLength[j,i]))
                    if (s.pol.mAPThreshold[i][aPolicyInd,j])<(s.income.mGrid[mUBInd])
                        cPol.consRVar[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*(c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])-cPol.consRGrid[j,i])^2
                        cDivYPol.consRDivYVar[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*(c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])-cDivYPol.consRDivYGrid[j,i])^2

                        cPol.consVar[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*(c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])-cPol.consGrid[j,i])^2
                        cDivYPol.consDivYVar[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*(c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])-cDivYPol.consDivYGrid[j,i])^2


                        tempPSum+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes
                        lastMVal=s.pol.mAPThreshold[i][aPolicyInd,j]
                        aPolicyInd+=1
                    else
                        cPol.consRVar[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*(c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])-cPol.consRGrid[j,i])^2
                        cDivYPol.consRDivYVar[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*(c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])-cDivYPol.consRDivYGrid[j,i])^2

                        cPol.consVar[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*(c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])-cPol.consGrid[j,i])^2
                        cDivYPol.consDivYVar[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*(c(m,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])-cDivYPol.consDivYGrid[j,i])^2

                        tempPSum+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes
                        lastMVal=s.income.mGrid[mUBInd]
                        mUBInd+=1
                    end
                end
                tempPSumInv=tempPSum^(-1)
                cPol.consVar[j,i]*=tempPSumInv
                cDivYPol.consDivYVar[j,i]*=tempPSumInv

                cPol.consRVar[j,i]*=(1.0-apDefPol.defaultProb[j,i])^(-1)
                cDivYPol.consRDivYVar[j,i]*=(1.0-apDefPol.defaultProb[j,i])^(-1)

            else
                cPol.consVar[j,i]=0.0
                cDivYPol.consDivYVar[j,i]=cDivYPol.consDivYDInitVar[i]
            end

        end
    end




    apDefPol.repayProb.=1.0.-apDefPol.defaultProb


    return meanPolicies(cPol,cDivYPol,apDefPol,uPol,yPol)
end



function iterateDistributionCell!(m::longTermBondSpec{F,S},s::longTermBondEval{F,S},mPol::meanPolicies{F,S},i::S,j::S,RDistEA::Array{F,2},RDistEP::Array{F,2},RNDDistEA::Array{F,2},RNDDistEP::Array{F,2},DDistEP::Array{F,1}) where{F<:Real,S<:Integer}

    if s.pol.alwaysDefault[j,i]==true
        DDistEP[i]+=RDistEA[j,i]
    elseif s.pol.neverDefault[j,i]==true
        for k in 1:(s.pol.mListLength[j,i])
            RDistEP[s.pol.apPolicy[i][k,j],i]+=mPol.apDefPol.apProbability[i][k,j]*RDistEA[j,i]
            RNDDistEP[s.pol.apPolicy[i][k,j],i]+=mPol.apDefPol.apProbability[i][k,j]*RNDDistEA[j,i]
        end


    else
        DDistEP[i]+=mPol.apDefPol.defaultProb[j,i]*RDistEA[j,i]
        for k in (s.pol.firstRepayInd[j,i]):(s.pol.mListLength[j,i])
            RDistEP[s.pol.apPolicy[i][k,j],i]+=mPol.apDefPol.apProbability[i][k,j]*RDistEA[j,i]
            RNDDistEP[s.pol.apPolicy[i][k,j],i]+=mPol.apDefPol.apProbability[i][k,j]*RNDDistEA[j,i]
        end
    end

    return RDistEP,RNDDistEP,DDistEP
end
#simulateDistribution generates the stationary, ergodic joint distribution of default state, income, and borrowing implied by an income process and government policy functions. Two versions are supplied. The first sets the initial distribution to have equal mass in every state with zero borrowing and not in default and all other states to have zero mass. The second accepts as an initial distribution an object of type debtDist. The arguments of the first version are:
#1. m: a model specification
#2. s: a collection of objects used in solving the model
#3. tol: the tolerance for convergence
#4. maxIter: the maximum number of iterations should convergence not be achieved
#5. scaleFactor: a number by which to multiply the joint distribution during the calculation in order to insure that small masses which appear early in the process do not disappear before convergence due to issues of the numerical precision available to the computer for the specified type
#Its output is an object of type debtDist

function simulatePaths(m::longTermBondSpec{F,S},s::longTermBondEval{F,S},mPol::meanPolicies{F,S},tol::F,maxIter::S,bigT::S) where{F<:Real,S<:Integer}

    #initialize variables containing the dimensions of the income grid and borrowing grid
    ydim=m.yParams.nPoints
    adim=m.aPoints

    #Initialize new and old copies of:
    #1. the noninflated joint distribution of income and borrowing at the beginning of the period, conditional on having access to financial markets, once uncertainty regarding regaining access to financial markets and income transitions in the current period has been resolved
    #2. the noninflated joint distribution of income and borrowing at the end of the period, conditional on having access to financial markets, before uncertainty regarding regaining access to financial markets and income transitions in the next period has been resolved
    #3. the noninflated distribution of income at the beginning of each period, conditional on being in default, once uncertainty regarding regaining access to financial markets and income transitions in the current period has been resolved
    #4. the noninflated distribution of income at the end of each period, conditional on being in default, before uncertainty regarding regaining access to financial markets and income transitions in the next period has been resolved

    #The naming convention is that myXDistYZ is the noninflated distribution for access state X (R or D, with R indicating access to financial markets and D indicating no access), Y timing (EA for beginning of period, EP for end of period), and Z update marker (New for newest iteration or Old for previous iteration)
    myRDistEAOld=deepcopy(s.VF.vGrid)
    myRDistEANew=deepcopy(s.VF.vGrid)

    myRDistEPOld=deepcopy(s.VF.vGrid)
    myRDistEPNew=deepcopy(s.VF.vGrid)

    myRNDDistEAOld=deepcopy(s.VF.vGrid)
    myRNDDistEANew=deepcopy(s.VF.vGrid)

    myRNDDistEPOld=deepcopy(s.VF.vGrid)
    myRNDDistEPNew=deepcopy(s.VF.vGrid)


    myDDistEAOld=deepcopy(s.VF.vDFutGrid)
    myDDistEANew=deepcopy(s.VF.vDFutGrid)

    myDDistEPOld=deepcopy(s.VF.vDFutGrid)
    myDDistEPNew=deepcopy(s.VF.vDFutGrid)

    stIncDist=genStIncDist(m,s,tol,maxIter)
    pathsOut=makeMeanPaths(m,bigT)

    #Initialize an iteration counter and a measure of distance between successive iterations
    dREPDist=tol+one(F)
    dDEPDist=tol+one(F)
    dREADist=tol+one(F)
    dDEADist=tol+one(F)
    dDist=tol+one(F)

    #Put full mass, equally divided, on states in which the government has access to financial markets and 0.0 incoming borrowing.
    for i in 1:ydim
        for j in 1:adim
            myRDistEPOld[j,i]=0.0
            myRNDDistEPOld[j,i]=0.0
        end
        myDDistEPOld[i]=0.0
    end


    for i in 1:ydim
        myRDistEAOld[:,i].=0.0
        myRNDDistEAOld[:,i].=0.0

        myDDistEAOld[i]=0.0


        myRDistEAOld[s.a0Ind,i]=stIncDist[i]
        myRNDDistEAOld[s.a0Ind,i]=stIncDist[i]




        for j in 1:adim

            iterateDistributionCell!(m,s,mPol,i,j,myRDistEAOld,myRDistEPOld,myRNDDistEAOld,myRNDDistEPOld,myDDistEPOld)

        end
    end
    APCheck=initializePath(m,m.yParams.nPoints,bigT)
    APDivYCheck=initializePath(m,m.yParams.nPoints,bigT)

    #Incoming assets with the distribution variables containing group mean values by persistent income level
    pathsOut.assets.a.varDistPath[:,1].=0.0
    pathsOut.assets.a.varPath[1]=0.0
    pathsOut.assets.a.varNDDistPath[:,1].=0.0
    pathsOut.assets.a.varNDPath[1]=0.0

    #Incoming assets with the distribution variables containing group mean values by persistent income level conditional on repayment
    pathsOut.assets.aR.varDistPath[:,1].=0.0
    pathsOut.assets.aR.varPath[1]=0.0
    pathsOut.assets.aR.varNDDistPath[:,1].=0.0
    pathsOut.assets.aR.varNDPath[1]=0.0

    #Outgoing assets with the distribution variables containing group mean values by persistent income level conditional on repayment
    pathsOut.assets.ap.varDistPath[:,1].=mPol.apDefPol.apRGrid[s.a0Ind,:]
    pathsOut.assets.ap.varPath[1]=dot(stIncDist,pathsOut.assets.ap.varDistPath[:,1])
    pathsOut.assets.ap.varNDDistPath[:,1].=pathsOut.assets.ap.varDistPath[:,1]
    pathsOut.assets.ap.varNDPath[1]=pathsOut.assets.ap.varPath[1]

    APCheck.varDistPath[:,1].=mPol.apDefPol.apRGrid[s.a0Ind,:]
    APCheck.varPath[1]=dot(stIncDist,pathsOut.assets.ap.varDistPath[:,1])
    APCheck.varNDDistPath[:,1].=pathsOut.assets.ap.varDistPath[:,1]
    APCheck.varNDPath[1]=pathsOut.assets.ap.varPath[1]

    #Outgoing assets/current GDP with the distribution variables containing group mean values by persistent income level conditional on repayment
    pathsOut.assets.apDivY.varDistPath[:,1].=mPol.apDefPol.apRDivYGrid[s.a0Ind,:]
    pathsOut.assets.apDivY.varPath[1]=dot(stIncDist,pathsOut.assets.apDivY.varDistPath[:,1])
    pathsOut.assets.apDivY.varNDDistPath[:,1].=pathsOut.assets.apDivY.varDistPath[:,1]
    pathsOut.assets.apDivY.varNDPath[1]=pathsOut.assets.apDivY.varPath[1]

    APDivYCheck.varDistPath[:,1].=mPol.apDefPol.apRDivYGrid[s.a0Ind,:]
    APDivYCheck.varPath[1]=dot(stIncDist,pathsOut.assets.ap.varDistPath[:,1])
    APDivYCheck.varNDDistPath[:,1].=pathsOut.assets.ap.varDistPath[:,1]
    APDivYCheck.varNDPath[1]=pathsOut.assets.ap.varPath[1]

    #Consumption with the distribution variables containing group mean values by persistent income value
    pathsOut.consVal.cY.varDistPath[:,1].=mPol.cPol.consGrid[s.a0Ind,:]
    pathsOut.consVal.cY.varPath[1]=dot(stIncDist,pathsOut.consVal.cY.varDistPath[:,1])
    pathsOut.consVal.cY.varNDDistPath[:,1].=pathsOut.consVal.cY.varDistPath[:,1]
    pathsOut.consVal.cY.varNDPath[1]=pathsOut.consVal.cY.varPath[1]

    pathsOut.consDivY.cDivY.varDistPath[:,1].=mPol.cDivYPol.consDivYGrid[s.a0Ind,:]
    pathsOut.consDivY.cDivY.varPath[1]=dot(stIncDist,pathsOut.consDivY.cDivY.varDistPath[:,1])
    pathsOut.consDivY.cDivY.varNDDistPath[:,1].=pathsOut.consDivY.cDivY.varDistPath[:,1]
    pathsOut.consDivY.cDivY.varNDPath[1]=pathsOut.consDivY.cDivY.varPath[1]

    #Consumption with the distribution variables containing group mean values by persistent income value conditional on repayment
    pathsOut.consVal.cRY.varDistPath[:,1].=mPol.cPol.consGrid[s.a0Ind,:]
    pathsOut.consVal.cRY.varPath[1]=dot(stIncDist,pathsOut.consVal.cRY.varDistPath[:,1])
    pathsOut.consVal.cRY.varNDDistPath[:,1].=pathsOut.consVal.cRY.varDistPath[:,1]
    pathsOut.consVal.cRY.varNDPath[1]=pathsOut.consVal.cRY.varPath[1]

    pathsOut.consDivY.cRDivY.varDistPath[:,1].=mPol.cDivYPol.consRDivYGrid[s.a0Ind,:]
    pathsOut.consDivY.cRDivY.varPath[1]=dot(stIncDist,pathsOut.consDivY.cRDivY.varDistPath[:,1])
    pathsOut.consDivY.cRDivY.varNDDistPath[:,1].=pathsOut.consDivY.cRDivY.varDistPath[:,1]
    pathsOut.consDivY.cRDivY.varNDPath[1]=pathsOut.consDivY.cRDivY.varPath[1]

    #Consumption with the distribution variables containing group mean values by asset level
    pathsOut.consVal.cA.varDistPath[s.a0Ind,1]=dot(stIncDist,mPol.cPol.consGrid[s.a0Ind,:])
    pathsOut.consVal.cA.varPath[1]=dot(stIncDist,mPol.cPol.consGrid[s.a0Ind,:])
    pathsOut.consVal.cA.varNDDistPath[:,1].=pathsOut.consVal.cA.varDistPath[:,1]
    pathsOut.consVal.cA.varNDPath[1]=pathsOut.consVal.cA.varPath[1]

    #Consumption with the distribution variables containing group mean values by asset level conditional on repayment
    pathsOut.consVal.cRA.varDistPath[s.a0Ind,1]=dot(stIncDist,mPol.cPol.consGrid[s.a0Ind,:])
    pathsOut.consVal.cRA.varPath[1]=dot(stIncDist,mPol.cPol.consGrid[s.a0Ind,:])
    pathsOut.consVal.cRA.varNDDistPath[:,1].=pathsOut.consVal.cRA.varDistPath[:,1]
    pathsOut.consVal.cRA.varNDPath[1]=pathsOut.consVal.cRA.varPath[1]

    #Persistent income level with the distribution variables containing group mean values by asset level
    pathsOut.yVal.y.varDistPath[s.a0Ind,1]=dot(stIncDist,s.income.yGrid)
    pathsOut.yVal.y.varPath[1]=dot(stIncDist,s.income.yGrid)
    pathsOut.yVal.y.varNDDistPath[:,1].=pathsOut.yVal.y.varDistPath[:,1]
    pathsOut.yVal.y.varNDPath[1]=pathsOut.yVal.y.varPath[1]

    #Persistent income level with the distribution variables containing group mean values by incoming asset level, conditional on exiting the period in good standing
    pathsOut.yVal.yR.varDistPath[s.a0Ind,1]=dot(stIncDist,s.income.yGrid)
    pathsOut.yVal.yR.varPath[1]=dot(stIncDist,s.income.yGrid)
    pathsOut.yVal.yR.varNDDistPath[:,1].=pathsOut.yVal.yR.varDistPath[:,1]
    pathsOut.yVal.yR.varNDPath[1]=pathsOut.yVal.yR.varPath[1]


    pathsOut.uVal.ut.varDistPath[:,1].=mPol.uPol.uGrid[s.a0Ind,:]
    pathsOut.uVal.ut.varPath[1]=dot(stIncDist,pathsOut.uVal.ut.varDistPath[:,1])
    pathsOut.uVal.ut.varNDDistPath[:,1].=pathsOut.uVal.ut.varDistPath[:,1]
    pathsOut.uVal.ut.varNDPath[1]=pathsOut.uVal.ut.varPath[1]


    pathsOut.consVal.cYVar.varDistPath[:,1].=mPol.cPol.consVar[s.a0Ind,:]
    pathsOut.consVal.cYVar.varPath[1]=dot(stIncDist,pathsOut.consVal.cYVar.varDistPath[:,1])+dot(stIncDist,(pathsOut.consVal.cY.varDistPath[:,1].-pathsOut.consVal.cY.varPath[1]).^2)
    pathsOut.consVal.cYVar.varNDDistPath[:,1].=pathsOut.consVal.cYVar.varDistPath[:,1]
    pathsOut.consVal.cYVar.varNDPath[1]=pathsOut.consVal.cYVar.varPath[1]

    pathsOut.consDivY.cDivYVar.varDistPath[:,1].=mPol.cDivYPol.consDivYVar[s.a0Ind,:]
    pathsOut.consDivY.cDivYVar.varPath[1]=dot(stIncDist,pathsOut.consDivY.cDivYVar.varDistPath[:,1])+dot(stIncDist,(pathsOut.consDivY.cDivY.varDistPath[:,1].-pathsOut.consDivY.cDivY.varPath[1]).^2)
    pathsOut.consDivY.cDivYVar.varNDDistPath[:,1].=pathsOut.consDivY.cDivYVar.varDistPath[:,1]
    pathsOut.consDivY.cDivYVar.varNDPath[1]=pathsOut.consDivY.cDivYVar.varPath[1]


    #Consumption with the distribution variables containing group mean values by persistent income value
    pathsOut.consVal.cYNoPen.varDistPath[:,1].=mPol.cPol.consNoPen[s.a0Ind,:]
    pathsOut.consVal.cYNoPen.varPath[1]=dot(stIncDist,pathsOut.consVal.cYNoPen.varDistPath[:,1])
    pathsOut.consVal.cYNoPen.varNDDistPath[:,1].=pathsOut.consVal.cYNoPen.varDistPath[:,1]
    pathsOut.consVal.cYNoPen.varNDPath[1]=pathsOut.consVal.cYNoPen.varPath[1]

    #Consumption with the distribution variables containing group mean values by persistent income value
    pathsOut.yVal.outputWithPen.varDistPath[:,1].=mPol.yPol.yObs[s.a0Ind,:]
    pathsOut.yVal.outputWithPen.varPath[1]=dot(stIncDist,pathsOut.yVal.outputWithPen.varDistPath[:,1])
    pathsOut.yVal.outputWithPen.varNDDistPath[:,1].=pathsOut.yVal.outputWithPen.varDistPath[:,1]
    pathsOut.yVal.outputWithPen.varNDPath[1]=pathsOut.yVal.outputWithPen.varPath[1]

    #Consumption with the distribution variables containing group mean values by persistent income value
    pathsOut.uVal.utNoPen.varDistPath[:,1].=mPol.uPol.uNoPen[s.a0Ind,:]
    pathsOut.uVal.utNoPen.varPath[1]=dot(stIncDist,pathsOut.uVal.utNoPen.varDistPath[:,1])
    pathsOut.uVal.utNoPen.varNDDistPath[:,1].=pathsOut.uVal.utNoPen.varDistPath[:,1]
    pathsOut.uVal.utNoPen.varNDPath[1]=pathsOut.uVal.utNoPen.varPath[1]

    #Consumption with the distribution variables containing group mean values by persistent income value
    pathsOut.uVal.utOutput.varDistPath[:,1].=mPol.uPol.utYObs[s.a0Ind,:]
    pathsOut.uVal.utOutput.varPath[1]=dot(stIncDist,pathsOut.uVal.utOutput.varDistPath[:,1])
    pathsOut.uVal.utOutput.varNDDistPath[:,1].=pathsOut.uVal.utOutput.varDistPath[:,1]
    pathsOut.uVal.utOutput.varNDPath[1]=pathsOut.uVal.utOutput.varPath[1]



    #Total mass in default at the end of the period
    pathsOut.DMassPath[1]=0.0
    #Total mass with access to financial markets at the end of the period
    pathsOut.RMassPath[1]=1.0
    #Total mass which has never lost access to financial markets
    pathsOut.NDMassPath[1]=1.0

    littleT=2
    #Iterate until the maximum number of iterations has been reached or convergence (as defined by the tolerance) has been achieved
    while (littleT<=bigT)&(dDist>=tol)

        #Set the new ex ante distributions according to the values implied by the income state transition matrix only

        mul!(myRDistEANew,myRDistEPOld,s.income.yTMat')
        mul!(myDDistEANew,s.income.yTMat,myDDistEPOld)
        mul!(myRNDDistEANew,myRNDDistEPOld,s.income.yTMat')



        #Iterate over income states, adding proportion theta of each mass in default in that state to the mass in that state with access to financial markets and 0 debt while subtracting it from the mass in default
        for i in 1:ydim
            myRDistEANew[s.a0Ind,i]+=m.theta*myDDistEANew[i]
            myDDistEANew[i]=(1.0-m.theta)*myDDistEANew[i]
        end






        #Iterate over income states
        for i in 1:ydim
            #Reset the current column of the ex post distribution for repayment to 0.0
            myRDistEPNew[:,i].=0.0
            myRNDDistEPNew[:,i].=0.0

            #Set the ex post mass in default equal to the ex ante mass in default
            myDDistEPNew[i]=myDDistEANew[i]

            #Iterate over incoming borrowing states
            for j in 1:adim
                iterateDistributionCell!(m,s,mPol,i,j,myRDistEANew,myRDistEPNew,myRNDDistEANew,myRNDDistEPNew,myDDistEPNew)
            end
        end

        pathsOut.DMassPath[littleT]=sum(myDDistEPNew)
        pathsOut.RMassPath[littleT]=sum(myRDistEPNew)
        pathsOut.NDMassPath[littleT]=sum(myRNDDistEPNew)

        #Incoming assets with the distribution variables containing group mean values by persistent income level
        pathsOut.assets.a.varDistPath[:,littleT]=(s.aGrid'*myRDistEANew)./max.(sum(myRDistEANew,dims=1),eps(zero(F)))
        pathsOut.assets.a.varPath[littleT]=dot(sum(myRDistEANew,dims=1),pathsOut.assets.a.varDistPath[:,littleT])/max(sum(myRDistEANew),eps(zero(F)))
        pathsOut.assets.a.varNDDistPath[:,littleT]=(s.aGrid'*myRNDDistEANew)./max.(sum(myRNDDistEANew,dims=1),eps(zero(F)))
        pathsOut.assets.a.varNDPath[littleT]=dot(sum(myRNDDistEANew,dims=1),pathsOut.assets.a.varNDDistPath[:,littleT])/max(sum(myRNDDistEANew),eps(zero(F)))

        #Incoming assets with the distribution variables containing group mean values by persistent income level conditional on repayment
        pathsOut.assets.aR.varDistPath[:,littleT]=(s.aGrid'*(mPol.apDefPol.repayProb.*myRDistEANew))./max.(sum(mPol.apDefPol.repayProb.*myRDistEANew,dims=1),eps(zero(F)))
        pathsOut.assets.aR.varPath[littleT]=dot(sum(mPol.apDefPol.repayProb.*myRDistEANew,dims=1),pathsOut.assets.aR.varDistPath[:,littleT])/max(pathsOut.RMassPath[littleT],eps(zero(F)))
        pathsOut.assets.aR.varNDDistPath[:,littleT]=(s.aGrid'*(mPol.apDefPol.repayProb.*myRNDDistEANew))./max.(sum(mPol.apDefPol.repayProb.*myRNDDistEANew,dims=1),eps(zero(F)))
        pathsOut.assets.aR.varNDPath[littleT]=dot(sum(mPol.apDefPol.repayProb.*myRNDDistEANew,dims=1),pathsOut.assets.a.varNDDistPath[:,littleT])/max(pathsOut.NDMassPath[littleT],eps(zero(F)))

        #Outgoing assets with the distribution variables containing group mean values by persistent income level conditional on repayment
        pathsOut.assets.ap.varDistPath[:,littleT]=(s.aGrid'*myRDistEPNew)./max.(sum(myRDistEPNew,dims=1),eps(zero(F)))
        pathsOut.assets.ap.varPath[littleT]=dot(sum(myRDistEPNew,dims=2),s.aGrid)/max(pathsOut.RMassPath[littleT],eps(zero(F)))
        pathsOut.assets.ap.varNDDistPath[:,littleT]=(s.aGrid'*myRNDDistEPNew)./max.(sum(myRNDDistEPNew,dims=1),eps(zero(F)))
        pathsOut.assets.ap.varNDPath[littleT]=dot(sum(myRNDDistEPNew,dims=2),s.aGrid)/max(pathsOut.NDMassPath[littleT],eps(zero(F)))

        APCheck.varDistPath[:,littleT]=sum((mPol.apDefPol.repayProb.*myRDistEANew.*mPol.apDefPol.apRGrid),dims=1)./max.(sum(mPol.apDefPol.repayProb.*myRDistEANew,dims=1),eps(zero(F)))
        APCheck.varPath[littleT]=dot(sum(mPol.apDefPol.repayProb.*myRDistEANew,dims=1),APCheck.varDistPath[:,littleT])/max(sum(mPol.apDefPol.repayProb.*myRDistEANew),eps(zero(F)))
        APCheck.varNDDistPath[:,littleT]=sum((mPol.apDefPol.repayProb.*myRNDDistEANew.*mPol.apDefPol.apRGrid),dims=1)./max.(sum(mPol.apDefPol.repayProb.*myRNDDistEANew,dims=1),eps(zero(F)))
        APCheck.varNDPath[littleT]=dot(sum(mPol.apDefPol.repayProb.*myRNDDistEANew,dims=1),APCheck.varNDDistPath[:,littleT])/max(sum(mPol.apDefPol.repayProb.*myRNDDistEANew),eps(zero(F)))


        #Outgoing assets with the distribution variables containing group mean values by persistent income level conditional on repayment
        APDivYCheck.varDistPath[:,littleT]=(s.aGrid'*myRDistEPNew)./max.(sum(myRDistEPNew,dims=1),eps(zero(F)))./(s.income.yGrid')
        APDivYCheck.varPath[littleT]=dot(sum(myRDistEPNew,dims=1),APDivYCheck.varDistPath[:,littleT])/max(pathsOut.RMassPath[littleT],eps(zero(F)))
        APDivYCheck.varNDDistPath[:,littleT]=(s.aGrid'*myRNDDistEPNew)./max.(sum(myRNDDistEPNew,dims=1),eps(zero(F)))./(s.income.yGrid')
        APDivYCheck.varNDPath[littleT]=dot(sum(myRNDDistEPNew,dims=1),APDivYCheck.varNDDistPath[:,littleT])/max(pathsOut.NDMassPath[littleT],eps(zero(F)))

        pathsOut.assets.apDivY.varDistPath[:,littleT]=sum((mPol.apDefPol.repayProb.*myRDistEANew.*mPol.apDefPol.apRDivYGrid),dims=1)./max.(sum(mPol.apDefPol.repayProb.*myRDistEANew,dims=1),eps(zero(F)))
        pathsOut.assets.apDivY.varPath[littleT]=dot(sum(mPol.apDefPol.repayProb.*myRDistEANew,dims=1),pathsOut.assets.apDivY.varDistPath[:,littleT])/max(sum(mPol.apDefPol.repayProb.*myRDistEANew),eps(zero(F)))
        pathsOut.assets.apDivY.varNDDistPath[:,littleT]=sum((mPol.apDefPol.repayProb.*myRNDDistEANew.*mPol.apDefPol.apRDivYGrid),dims=1)./max.(sum(mPol.apDefPol.repayProb.*myRNDDistEANew,dims=1),eps(zero(F)))
        pathsOut.assets.apDivY.varNDPath[littleT]=dot(sum(mPol.apDefPol.repayProb.*myRNDDistEANew,dims=1),pathsOut.assets.apDivY.varNDDistPath[:,littleT])/max(sum(mPol.apDefPol.repayProb.*myRNDDistEANew),eps(zero(F)))


        #Persistent income level with the distribution variables containing group mean values by asset level
        pathsOut.yVal.y.varDistPath[:,littleT]=(myRDistEANew*s.income.yGrid)./max.(sum(myRDistEANew,dims=2),eps(zero(F)))
        pathsOut.yVal.y.varPath[littleT]=dot(sum(myRDistEANew,dims=1),s.income.yGrid)/max(sum(myRDistEANew),eps(zero(F)))
        pathsOut.yVal.y.varNDDistPath[:,littleT]=(myRNDDistEANew*s.income.yGrid)./max.(sum(myRNDDistEANew,dims=2),eps(zero(F)))
        pathsOut.yVal.y.varNDPath[littleT]=dot(sum(myRNDDistEANew,dims=1),s.income.yGrid)/max(sum(myRNDDistEANew),eps(zero(F)))

        #Persistent income level with the distribution variables containing group mean values by incoming asset level conditional on exiting the period in good standing
        pathsOut.yVal.yR.varDistPath[:,littleT]=((mPol.apDefPol.repayProb.*myRDistEANew)*s.income.yGrid)./max.(sum(mPol.apDefPol.repayProb.*myRDistEANew,dims=2),eps(zero(F)))
        pathsOut.yVal.yR.varPath[littleT]=dot(sum(mPol.apDefPol.repayProb.*myRDistEANew,dims=1),s.income.yGrid)/max(sum(mPol.apDefPol.repayProb.*myRDistEANew),eps(zero(F)))
        pathsOut.yVal.yR.varNDDistPath[:,littleT]=((mPol.apDefPol.repayProb.*myRNDDistEANew)*s.income.yGrid)./max.(sum(mPol.apDefPol.repayProb.*myRNDDistEANew,dims=2),eps(zero(F)))
        pathsOut.yVal.yR.varNDPath[littleT]=dot(sum(mPol.apDefPol.repayProb.*myRNDDistEANew,dims=1),s.income.yGrid)/max(sum(mPol.apDefPol.repayProb.*myRNDDistEANew),eps(zero(F)))

        #Consumption with the distribution variables containing group mean values by persistent income value
        pathsOut.consVal.cY.varDistPath[:,littleT]=(sum(myRDistEANew.*mPol.cPol.consGrid,dims=1).+myDDistEANew'.*s.income.yDefGrid')./max.(sum(myRDistEANew,dims=1).+myDDistEANew',eps(zero(F)))
        pathsOut.consVal.cY.varPath[littleT]=dot((sum(myRDistEANew,dims=1).+myDDistEANew'),pathsOut.consVal.cY.varDistPath[:,littleT])/(sum(myRDistEANew)+sum(myDDistEANew))
        pathsOut.consVal.cY.varNDDistPath[:,littleT]=sum(myRNDDistEANew.*mPol.cPol.consGrid,dims=1)./max.(sum(myRNDDistEANew,dims=1),eps(zero(F)))
        pathsOut.consVal.cY.varNDPath[littleT]=sum(myRNDDistEANew.*mPol.cPol.consGrid)/max(sum(myRNDDistEANew),eps(zero(F)))

        pathsOut.consDivY.cDivY.varDistPath[:,littleT]=(sum(myRDistEANew.*mPol.cDivYPol.consDivYGrid,dims=1).+myDDistEANew'.*mPol.cDivYPol.consDivYDFutGrid')./max.(sum(myRDistEANew,dims=1).+myDDistEANew',eps(zero(F)))
        pathsOut.consDivY.cDivY.varPath[littleT]=dot((sum(myRDistEANew,dims=1).+myDDistEANew'),pathsOut.consDivY.cDivY.varDistPath[:,littleT])/(sum(myRDistEANew)+sum(myDDistEANew))
        pathsOut.consDivY.cDivY.varNDDistPath[:,littleT]=sum(myRNDDistEANew.*mPol.cDivYPol.consDivYGrid,dims=1)./max.(sum(myRNDDistEANew,dims=1),eps(zero(F)))
        pathsOut.consDivY.cDivY.varNDPath[littleT]=sum(myRNDDistEANew.*mPol.cDivYPol.consDivYGrid)/max(sum(myRNDDistEANew),eps(zero(F)))

        #Consumption with the distribution variables containing group mean values by persistent income value conditional on repayment
        pathsOut.consVal.cRY.varDistPath[:,littleT]=sum(mPol.apDefPol.repayProb.*myRDistEANew.*mPol.cPol.consRGrid,dims=1)./max.(sum(mPol.apDefPol.repayProb.*myRDistEANew,dims=1),eps(zero(F)))
        pathsOut.consVal.cRY.varPath[littleT]=dot(sum(mPol.apDefPol.repayProb.*myRDistEANew,dims=1),pathsOut.consVal.cRY.varDistPath[:,littleT])/max(sum(mPol.apDefPol.repayProb.*myRDistEANew),eps(zero(F)))
        pathsOut.consVal.cRY.varNDDistPath[:,littleT]=sum(mPol.apDefPol.repayProb.*myRNDDistEANew.*mPol.cPol.consRGrid,dims=1)./max.(sum(mPol.apDefPol.repayProb.*myRNDDistEANew,dims=1),eps(zero(F)))
        pathsOut.consVal.cRY.varNDPath[littleT]=sum(mPol.apDefPol.repayProb.*myRNDDistEANew.*mPol.cPol.consRGrid)/max(sum(mPol.apDefPol.repayProb.*myRNDDistEANew),eps(zero(F)))

        pathsOut.consDivY.cRDivY.varDistPath[:,littleT]=sum(mPol.apDefPol.repayProb.*myRDistEANew.*mPol.cDivYPol.consRDivYGrid,dims=1)./max.(sum(mPol.apDefPol.repayProb.*myRDistEANew,dims=1),eps(zero(F)))
        pathsOut.consDivY.cRDivY.varPath[littleT]=dot(sum(mPol.apDefPol.repayProb.*myRDistEANew,dims=1),pathsOut.consDivY.cRDivY.varDistPath[:,littleT])/max(sum(mPol.apDefPol.repayProb.*myRDistEANew),eps(zero(F)))
        pathsOut.consDivY.cRDivY.varNDDistPath[:,littleT]=sum(mPol.apDefPol.repayProb.*myRNDDistEANew.*mPol.cDivYPol.consRDivYGrid,dims=1)./max.(sum(mPol.apDefPol.repayProb.*myRNDDistEANew,dims=1),eps(zero(F)))
        pathsOut.consDivY.cRDivY.varNDPath[littleT]=sum(mPol.apDefPol.repayProb.*myRNDDistEANew.*mPol.cDivYPol.consRDivYGrid)/max(sum(mPol.apDefPol.repayProb.*myRNDDistEANew),eps(zero(F)))

        #Consumption with the distribution variables containing group mean values by asset level
        pathsOut.consVal.cA.varDistPath[:,littleT]=sum(myRDistEANew.*mPol.cPol.consGrid,dims=2)./max.(sum(myRDistEANew,dims=2),eps(zero(F)))
        pathsOut.consVal.cA.varPath[littleT]=dot(sum(myRDistEANew,dims=2),pathsOut.consVal.cA.varDistPath[:,littleT])/max(sum(myRDistEANew),eps(zero(F)))
        pathsOut.consVal.cA.varNDDistPath[:,littleT]=sum(myRNDDistEANew.*mPol.cPol.consGrid,dims=2)./max.(sum(myRNDDistEANew,dims=2),eps(zero(F)))
        pathsOut.consVal.cA.varNDPath[littleT]=sum(myRNDDistEANew.*mPol.cPol.consGrid)/max(sum(myRNDDistEANew),eps(zero(F)))

        #Consumption with the distribution variables containing group mean values by asset level conditional on repayment
        pathsOut.consVal.cRA.varDistPath[:,littleT]=(sum(mPol.apDefPol.repayProb.*myRDistEANew.*mPol.cPol.consRGrid,dims=2))./max.(sum(mPol.apDefPol.repayProb.*myRDistEANew,dims=2),eps(zero(F)))
        pathsOut.consVal.cRA.varPath[littleT]=dot(sum(mPol.apDefPol.repayProb.*myRDistEANew,dims=2),pathsOut.consVal.cRA.varDistPath[:,littleT])/max(sum(mPol.apDefPol.repayProb.*myRDistEANew),eps(zero(F)))
        pathsOut.consVal.cRA.varNDDistPath[:,littleT]=sum(mPol.apDefPol.repayProb.*myRNDDistEANew.*mPol.cPol.consRGrid,dims=2)./max.(sum(mPol.apDefPol.repayProb.*myRNDDistEANew,dims=2),eps(zero(F)))
        pathsOut.consVal.cRA.varNDPath[littleT]=sum(mPol.apDefPol.repayProb.*myRNDDistEANew.*mPol.cPol.consRGrid)/max(sum(mPol.apDefPol.repayProb.*myRNDDistEANew),eps(zero(F)))

        pathsOut.uVal.ut.varDistPath[:,littleT]=(sum((myRDistEANew.*mPol.uPol.uGrid),dims=1)+myDDistEANew'.*mPol.uPol.uDFutGrid')./max.(sum(myRDistEANew,dims=1).+myDDistEANew',eps(zero(F)))
        pathsOut.uVal.ut.varPath[littleT]=dot(sum(myRDistEANew,dims=1).+myDDistEANew',pathsOut.uVal.ut.varDistPath[:,littleT])/max(sum(myRDistEANew)+sum(myDDistEANew),eps(zero(F)))
        pathsOut.uVal.ut.varNDDistPath[:,littleT]=sum((myRNDDistEANew.*mPol.uPol.uGrid),dims=1)./max.(sum(myRNDDistEANew,dims=1),eps(zero(F)))
        pathsOut.uVal.ut.varNDPath[littleT]=dot(sum(myRNDDistEANew,dims=1),pathsOut.uVal.ut.varNDDistPath[:,littleT])/max(sum(myRNDDistEANew),eps(zero(F)))

        #Consumption with the distribution variables containing group mean values by persistent income value
        pathsOut.consVal.cYNoPen.varDistPath[:,littleT]=(sum(myRDistEANew.*mPol.cPol.consNoPen,dims=1).+myDDistEANew'.*mPol.cPol.consNoPenD')./max.(sum(myRDistEANew,dims=1).+myDDistEANew',eps(zero(F)))
        pathsOut.consVal.cYNoPen.varPath[littleT]=dot((sum(myRDistEANew,dims=1).+myDDistEANew'),pathsOut.consVal.cYNoPen.varDistPath[:,littleT])/(sum(myRDistEANew)+sum(myDDistEANew))
        pathsOut.consVal.cYNoPen.varNDDistPath[:,littleT]=sum(myRNDDistEANew.*mPol.cPol.consNoPen,dims=1)./max.(sum(myRNDDistEANew,dims=1),eps(zero(F)))
        pathsOut.consVal.cYNoPen.varNDPath[littleT]=sum(myRNDDistEANew.*mPol.cPol.consNoPen)/max(sum(myRNDDistEANew),eps(zero(F)))

        #Consumption with the distribution variables containing group mean values by persistent income value
        pathsOut.yVal.outputWithPen.varDistPath[:,littleT]=(sum(myRDistEANew.*mPol.yPol.yObs,dims=1).+myDDistEANew'.*mPol.yPol.yObsD')./max.(sum(myRDistEANew,dims=1).+myDDistEANew',eps(zero(F)))
        pathsOut.yVal.outputWithPen.varPath[littleT]=dot((sum(myRDistEANew,dims=1).+myDDistEANew'),pathsOut.yVal.outputWithPen.varDistPath[:,littleT])/(sum(myRDistEANew)+sum(myDDistEANew))
        pathsOut.yVal.outputWithPen.varNDDistPath[:,littleT]=sum(myRNDDistEANew.*mPol.yPol.yObs,dims=1)./max.(sum(myRNDDistEANew,dims=1),eps(zero(F)))
        pathsOut.yVal.outputWithPen.varNDPath[littleT]=sum(myRNDDistEANew.*mPol.yPol.yObs)/max(sum(myRNDDistEANew),eps(zero(F)))

        #Consumption with the distribution variables containing group mean values by persistent income value
        pathsOut.uVal.utNoPen.varDistPath[:,littleT]=(sum(myRDistEANew.*mPol.uPol.uNoPen,dims=1).+myDDistEANew'.*mPol.uPol.uNoPenD')./max.(sum(myRDistEANew,dims=1).+myDDistEANew',eps(zero(F)))
        pathsOut.uVal.utNoPen.varPath[littleT]=dot((sum(myRDistEANew,dims=1).+myDDistEANew'),pathsOut.uVal.utNoPen.varDistPath[:,littleT])/(sum(myRDistEANew)+sum(myDDistEANew))
        pathsOut.uVal.utNoPen.varNDDistPath[:,littleT]=sum(myRNDDistEANew.*mPol.uPol.uNoPen,dims=1)./max.(sum(myRNDDistEANew,dims=1),eps(zero(F)))
        pathsOut.uVal.utNoPen.varNDPath[littleT]=sum(myRNDDistEANew.*mPol.uPol.uNoPen)/max(sum(myRNDDistEANew),eps(zero(F)))

        #Consumption with the distribution variables containing group mean values by persistent income value
        pathsOut.uVal.utOutput.varDistPath[:,littleT]=(sum(myRDistEANew.*mPol.uPol.utYObs,dims=1).+myDDistEANew'.*mPol.uPol.utYObsD')./max.(sum(myRDistEANew,dims=1).+myDDistEANew',eps(zero(F)))
        pathsOut.uVal.utOutput.varPath[littleT]=dot((sum(myRDistEANew,dims=1).+myDDistEANew'),pathsOut.uVal.utOutput.varDistPath[:,littleT])/(sum(myRDistEANew)+sum(myDDistEANew))
        pathsOut.uVal.utOutput.varNDDistPath[:,littleT]=sum(myRNDDistEANew.*mPol.uPol.utYObs,dims=1)./max.(sum(myRNDDistEANew,dims=1),eps(zero(F)))
        pathsOut.uVal.utOutput.varNDPath[littleT]=sum(myRNDDistEANew.*mPol.uPol.utYObs)/max(sum(myRNDDistEANew),eps(zero(F)))






        for i in 1:m.yParams.nPoints
            #get weighted average of within group variation E[var(y_ijk|i)]
            pathsOut.consVal.cYVar.varDistPath[i,littleT]=dot(myRDistEANew[:,i],mPol.cPol.consVar[:,i])+myDDistEANew[i]*mPol.cPol.consDVar[i]
            #add between group variation var(E[y_ijk|i])
            pathsOut.consVal.cYVar.varDistPath[i,littleT]+=dot(myRDistEANew[:,i],(mPol.cPol.consGrid[:,i].-pathsOut.consVal.cY.varDistPath[i,littleT]).^2)+myDDistEANew[i]*(s.income.yDefGrid[i]-pathsOut.consVal.cY.varDistPath[i,littleT])^2
            #inflate using total probability of state i
            pathsOut.consVal.cYVar.varDistPath[i,littleT]*=(max(eps(zero(F)),sum(myRDistEANew[:,i])+myDDistEANew[i]))^(-1)

            pathsOut.consVal.cYVar.varNDDistPath[i,littleT]=dot(myRNDDistEANew[:,i],mPol.cPol.consVar[:,i])
            pathsOut.consVal.cYVar.varNDDistPath[i,littleT]+=dot(myRNDDistEANew[:,i],(mPol.cPol.consGrid[:,i].-pathsOut.consVal.cY.varNDDistPath[i,littleT]).^2)
            pathsOut.consVal.cYVar.varNDDistPath[i,littleT]*=(max(eps(zero(F)),sum(myRNDDistEANew[:,i])))^(-1)
        end
        #get weighted average of within group variation E[var(y_ijk|i,j)]
        pathsOut.consVal.cYVar.varPath[littleT]=sum(myRDistEANew.*mPol.cPol.consVar)+sum(myDDistEANew.*mPol.cPol.consDVar)
        #add between group variation var(E[y_ijk|i,j])
        pathsOut.consVal.cYVar.varPath[littleT]+=sum(myRDistEANew.*(mPol.cPol.consGrid.-pathsOut.consVal.cY.varPath[littleT]).^2)+sum(myDDistEANew.*(s.income.yDefGrid.-pathsOut.consVal.cY.varPath[littleT]).^2)
        #inflate using total probability
        pathsOut.consVal.cYVar.varPath[littleT]*=(max(eps(zero(F)),sum(myRDistEANew)+sum(myDDistEANew)))^(-1)

        #get weighted average of within group variation E[var(y_ijk|i,j)]
        pathsOut.consVal.cYVar.varNDPath[littleT]=sum(myRNDDistEANew.*mPol.cPol.consVar)
        #add between group variation var(E[y_ijk|i,j])
        pathsOut.consVal.cYVar.varNDPath[littleT]+=sum(myRNDDistEANew.*(mPol.cPol.consGrid.-pathsOut.consVal.cY.varNDPath[littleT]).^2)
        #inflate using total probability
        pathsOut.consVal.cYVar.varNDPath[littleT]*=(max(eps(zero(F)),sum(myRNDDistEANew)))^(-1)

        for i in 1:m.yParams.nPoints
            #get weighted average of within group variation E[var(y_ijk|i)]
            pathsOut.consDivY.cDivYVar.varDistPath[i,littleT]=dot(myRDistEANew[:,i],mPol.cDivYPol.consDivYVar[:,i])+myDDistEANew[i]*mPol.cDivYPol.consDivYDFutVar[i]
            #add between group variation var(E[y_ijk|i])
            pathsOut.consDivY.cDivYVar.varDistPath[i,littleT]+=dot(myRDistEANew[:,i],(mPol.cDivYPol.consDivYGrid[:,i].-pathsOut.consDivY.cDivY.varDistPath[i,littleT]).^2)+myDDistEANew[i]*(mPol.cDivYPol.consDivYDFutGrid[i]-pathsOut.consVal.cY.varDistPath[i,littleT])^2
            #inflate using total probability of state i
            pathsOut.consDivY.cDivYVar.varDistPath[i,littleT]*=(max(eps(zero(F)),sum(myRDistEANew[:,i])+myDDistEANew[i]))^(-1)

            pathsOut.consDivY.cDivYVar.varNDDistPath[i,littleT]=dot(myRNDDistEANew[:,i],mPol.cDivYPol.consDivYVar[:,i])
            pathsOut.consDivY.cDivYVar.varNDDistPath[i,littleT]+=dot(myRNDDistEANew[:,i],(mPol.cDivYPol.consDivYGrid[:,i].-pathsOut.consDivY.cDivY.varNDDistPath[i,littleT]).^2)
            pathsOut.consDivY.cDivYVar.varNDDistPath[i,littleT]*=(max(eps(zero(F)),sum(myRNDDistEANew[:,i])))^(-1)
        end
        #get weighted average of within group variation E[var(y_ijk|i,j)]
        pathsOut.consDivY.cDivYVar.varPath[littleT]=sum(myRDistEANew.*mPol.cDivYPol.consDivYVar)+sum(myDDistEANew.*mPol.cDivYPol.consDivYDFutVar)
        #add between group variation var(E[y_ijk|i,j])
        pathsOut.consDivY.cDivYVar.varPath[littleT]+=sum(myRDistEANew.*(mPol.cDivYPol.consDivYGrid.-pathsOut.consDivY.cDivY.varPath[littleT]).^2)+sum(myDDistEANew.*(mPol.cDivYPol.consDivYDFutGrid.-pathsOut.consDivY.cDivY.varPath[littleT]).^2)
        #inflate using total probability
        pathsOut.consDivY.cDivYVar.varPath[littleT]*=(max(eps(zero(F)),sum(myRDistEANew)+sum(myDDistEANew)))^(-1)

        #get weighted average of within group variation E[var(y_ijk|i,j)]
        pathsOut.consDivY.cDivYVar.varNDPath[littleT]=sum(myRNDDistEANew.*mPol.cDivYPol.consDivYVar)
        #add between group variation var(E[y_ijk|i,j])
        pathsOut.consDivY.cDivYVar.varNDPath[littleT]+=sum(myRNDDistEANew.*(mPol.cDivYPol.consDivYGrid.-pathsOut.consDivY.cDivY.varNDPath[littleT]).^2)
        #inflate using total probability
        pathsOut.consDivY.cDivYVar.varNDPath[littleT]*=(max(eps(zero(F)),sum(myRNDDistEANew)))^(-1)





        #Calculate the various sup norm measures of distance
        dREPDist=maximum(abs.(myRDistEPNew-myRDistEPOld))
        dDEPDist=maximum(abs.(myDDistEPNew-myDDistEPOld))
        dREADist=maximum(abs.(myRDistEANew-myRDistEAOld))
        dDEADist=maximum(abs.(myDDistEANew-myDDistEAOld))
        dDist=max(dREPDist,dDEPDist,dREADist,dDEADist)

        #Update the old distributions
        copy!(myRDistEPOld,myRDistEPNew)
        copy!(myRNDDistEPOld,myRNDDistEPNew)
        copy!(myDDistEPOld,myDDistEPNew)
        copy!(myRDistEAOld,myRDistEANew)
        copy!(myRNDDistEAOld,myRNDDistEANew)
        copy!(myDDistEAOld,myDDistEANew)

        #At each iteration which is a multiple of 10% of the maximum number of iterations, print the current iteration number and variables tracking the convergence criteria
        if div(littleT,max(div(bigT,10),1))==(littleT/max(div(bigT,10),1))
            println([littleT,dREPDist,dDEPDist,dREADist,dDEADist])
        end
        #Increment the iteration counter
        littleT+=1
    end
    for littleT in 1:bigT
        if pathsOut.uVal.ut.varPath[littleT]!=0.0
            pathsOut.consVal.cCert.varPath[littleT]=uInverse(m,pathsOut.uVal.ut.varPath[littleT])
        end
        if pathsOut.uVal.ut.varNDPath[littleT]!=0.0
            pathsOut.consVal.cCert.varNDPath[littleT]=uInverse(m,pathsOut.uVal.ut.varNDPath[littleT])
        end
        for i in 1:(m.yParams.nPoints)
            if pathsOut.uVal.ut.varDistPath[i,littleT]!=0.0
                pathsOut.consVal.cCert.varDistPath[i,littleT]=uInverse(m,pathsOut.uVal.ut.varDistPath[i,littleT])
            end
            if pathsOut.uVal.ut.varNDDistPath[i,littleT]!=0.0
                pathsOut.consVal.cCert.varNDDistPath[i,littleT]=uInverse(m,pathsOut.uVal.ut.varNDDistPath[i,littleT])
            end
        end
        if pathsOut.uVal.utNoPen.varPath[littleT]!=0.0
            pathsOut.consVal.cYNoPenCert.varPath[littleT]=uInverse(m,pathsOut.uVal.utNoPen.varPath[littleT])
        end
        if pathsOut.uVal.utNoPen.varNDPath[littleT]!=0.0
            pathsOut.consVal.cYNoPenCert.varNDPath[littleT]=uInverse(m,pathsOut.uVal.utNoPen.varNDPath[littleT])
        end
        for i in 1:(m.yParams.nPoints)
            if pathsOut.uVal.utNoPen.varDistPath[i,littleT]!=0.0
                pathsOut.consVal.cYNoPenCert.varDistPath[i,littleT]=uInverse(m,pathsOut.uVal.utNoPen.varDistPath[i,littleT])
            end
            if pathsOut.uVal.utNoPen.varNDDistPath[i,littleT]!=0.0
                pathsOut.consVal.cYNoPenCert.varNDDistPath[i,littleT]=uInverse(m,pathsOut.uVal.utNoPen.varNDDistPath[i,littleT])
            end
        end
        if pathsOut.uVal.utOutput.varPath[littleT]!=0.0
            pathsOut.yVal.outputCert.varPath[littleT]=uInverse(m,pathsOut.uVal.utOutput.varPath[littleT])
        end
        if pathsOut.uVal.utOutput.varNDPath[littleT]!=0.0
            pathsOut.yVal.outputCert.varNDPath[littleT]=uInverse(m,pathsOut.uVal.utOutput.varNDPath[littleT])
        end
        for i in 1:(m.yParams.nPoints)
            if pathsOut.uVal.utOutput.varDistPath[i,littleT]!=0.0
                pathsOut.yVal.outputCert.varDistPath[i,littleT]=uInverse(m,pathsOut.uVal.utOutput.varDistPath[i,littleT])
            end
            if pathsOut.uVal.utOutput.varNDDistPath[i,littleT]!=0.0
                pathsOut.yVal.outputCert.varNDDistPath[i,littleT]=uInverse(m,pathsOut.uVal.utOutput.varNDDistPath[i,littleT])
            end
        end
    end



    println([littleT,dREPDist,dDEPDist,dREADist,dDEADist])


    return pathsOut,APCheck,APDivYCheck

end


#Generates the stationary distribution of the persistent component of the income process
function genStIncDist(m::longTermBondSpec{F,S},s::longTermBondEval{F,S},tol::F,maxIter::S) where{F<:Real,S<:Integer}

    oldDist=ones(m.yParams.nPoints)/m.yParams.nPoints
    newDist=zeros(m.yParams.nPoints)

    iCount=1
    dDist=1.0+tol

    while (iCount<=maxIter)&&(dDist>tol)
        mul!(newDist,s.income.yTMat,oldDist)

        dDist=maximum(abs.(newDist-oldDist))



        if div(iCount,max(div(maxIter,10),1))==(iCount/max(div(maxIter,10),1))
            println([iCount,dDist])
        end
        copy!(oldDist,newDist)
        iCount+=1
    end
    deflFac=sum(newDist)^(-1)
    return newDist*deflFac
end



function writePolicies(s::longTermBondEval{F,S},mPol::meanPolicies{F,S},filePrefix::String,fileDir::String) where{F<:Real,S<:Integer}
    writecsv(fileDir*filePrefix*"_apDefPol.apRGrid.csv",mPol.apDefPol.apRGrid)
    writecsv(fileDir*filePrefix*"_apDefPol.apRDivYGrid.csv",mPol.apDefPol.apRDivYGrid)
    writecsv(fileDir*filePrefix*"_cPol.consRGrid.csv",mPol.cPol.consRGrid)
    writecsv(fileDir*filePrefix*"_cDivYPol.consRDivYGrid.csv",mPol.cDivYPol.consRDivYGrid)
    writecsv(fileDir*filePrefix*"_cPol.consGrid.csv",mPol.cPol.consGrid)
    writecsv(fileDir*filePrefix*"_cDivYPol.consDivYGrid.csv",mPol.cDivYPol.consDivYGrid)
    writecsv(fileDir*filePrefix*"_cDivYPol.consDivYDInitGrid.csv",mPol.cDivYPol.consDivYDInitGrid)
    writecsv(fileDir*filePrefix*"_cDivYPol.consDivYDFutGrid.csv",mPol.cDivYPol.consDivYDFutGrid)
    writecsv(fileDir*filePrefix*"_uPol.uGrid.csv",mPol.uPol.uGrid)
    writecsv(fileDir*filePrefix*"_uPol.uDInitGrid.csv",mPol.uPol.uDInitGrid)
    writecsv(fileDir*filePrefix*"_uPol.uDFutGrid.csv",mPol.uPol.uDFutGrid)
    writecsv(fileDir*filePrefix*"_apDefPol.defaultProb.csv",mPol.apDefPol.defaultProb)
    writecsv(fileDir*filePrefix*"_apDefPol.repayProb.csv",mPol.apDefPol.repayProb)
    writecsv(fileDir*filePrefix*"_apDefPol.lastAlwaysDefInd.csv",mPol.apDefPol.lastAlwaysDefInd)
    writecsv(fileDir*filePrefix*"_yGrid.csv",s.income.yGrid)
    writecsv(fileDir*filePrefix*"_yDefGrid.csv",s.income.yDefGrid)
    writecsv(fileDir*filePrefix*"_aGrid.csv",s.aGrid)
end

function writeSinglePath(currentPath::singlePath{F},varName::String,filePrefix::String,fileDir::String) where{F<:Real}
    writecsv(fileDir*filePrefix*"_"*varName*"Path.csv",currentPath.varPath)
    writecsv(fileDir*filePrefix*"_"*varName*"DistPath.csv",currentPath.varDistPath)
    writecsv(fileDir*filePrefix*"_"*varName*"NDPath.csv",currentPath.varNDPath)
    writecsv(fileDir*filePrefix*"_"*varName*"NDDistPath.csv",currentPath.varNDDistPath)
end


function writeAllPaths(allPaths::meanPaths{F},filePrefix::String,fileDir::String) where{F<:Real}
    writeSinglePath(allPaths.a,"a",filePrefix,fileDir)
    writeSinglePath(allPaths.aR,"aR",filePrefix,fileDir)
    writeSinglePath(allPaths.ap,"ap",filePrefix,fileDir)
    writeSinglePath(allPaths.apDivY,"apDivY",filePrefix,fileDir)
    writeSinglePath(allPaths.cY,"cY",filePrefix,fileDir)
    writeSinglePath(allPaths.cDivY,"cDivY",filePrefix,fileDir)
    writeSinglePath(allPaths.cRY,"cRY",filePrefix,fileDir)
    writeSinglePath(allPaths.cRDivY,"cRDivY",filePrefix,fileDir)
    writeSinglePath(allPaths.cA,"cA",filePrefix,fileDir)
    writeSinglePath(allPaths.cRA,"cRA",filePrefix,fileDir)
    writeSinglePath(allPaths.y,"y",filePrefix,fileDir)
    writeSinglePath(allPaths.yR,"yR",filePrefix,fileDir)
    writeSinglePath(allPaths.ut,"ut",filePrefix,fileDir)
    writeSinglePath(allPaths.cCert,"cCert",filePrefix,fileDir)
    writecsv(fileDir*filePrefix*"_DMassPath.csv",allPaths.DMassPath)
    writecsv(fileDir*filePrefix*"_RMassPath.csv",allPaths.RMassPath)
    writecsv(fileDir*filePrefix*"_NDMassPath.csv",allPaths.NDMassPath)
    writeSinglePath(allPaths.cYNoPen,"cYNoPen",filePrefix,fileDir)
    writeSinglePath(allPaths.cYNoPenCert,"cYNoPenCert",filePrefix,fileDir)
    writeSinglePath(allPaths.utNoPen,"utNoPen",filePrefix,fileDir)
    writeSinglePath(allPaths.outputWithPen,"outputWithPen",filePrefix,fileDir)
    writeSinglePath(allPaths.outputCert,"outputCert",filePrefix,fileDir)
    writeSinglePath(allPaths.utOutput,"utOutput",filePrefix,fileDir)
end
