function changeBeta(m::longTermBondSpec{F,S},newBeta::F) where{F<:Real,S<:Integer}

    return longTermBondSpec(newBeta,m.theta,m.gamma,m.hpen0,m.hpen1,m.R,m.lambda,m.coup,m.aBounds,m.aPoints,m.yParams,m.mParams,m.simplePen,m.mixFacQ,m.mixFacV)
end


function getNextX!(xGrid::Array{F,1},fGrid::Array{F,1}) where{F<:Real}
    newBound=xGrid[2]
    tempSlope=0.0
    if fGrid[2]<0.0
        tempSlope=(fGrid[3]-fGrid[2])/(xGrid[3]-xGrid[2])

        #f=(f3-f2)/(x3-x2)*x+b
        #f3=(f3-f2)/(x3-x2)*x3+b
        #f3-(f3-f2)/(x3-x2)*x3=b
        tempIntercept=fGrid[3]-tempSlope*xGrid[3]

        #0=m*x+b
        #x=-b/m
        xGrid[2]=-tempIntercept/tempSlope
        xGrid[1]=newBound
        fGrid[1]=fGrid[2]

    else
        tempSlope=(fGrid[2]-fGrid[1])/(xGrid[2]-xGrid[1])

        #f=(f2-f1)/(x2-x1)*x+b
        #f2=(f2-f1)/(x2-x1)*x2+b
        #f2-(f2-f1)/(x2-x1)*x2=b
        tempIntercept=fGrid[2]-tempSlope*xGrid[2]

        #0=m*x+b
        #x=-b/m
        xGrid[2]=-tempIntercept/tempSlope
        xGrid[3]=newBound
        fGrid[3]=fGrid[2]

    end

    if (xGrid[2]>xGrid[3])||(xGrid[2]<xGrid[1])||(isinf(xGrid[2])==true)||(isnan(xGrid[2])==true)
        xGrid[2]=0.5*(xGrid[3]+xGrid[1])
    end
    return nothing
end


function findBetaIndiffRange(mOrig::longTermBondSpec{F,S},betaGBounds::Array{F,1},betaGPoints::S,tolG::F,maxIterG::S,tolC::F,maxIterC::S,xtol::F,ftol::F,maxIterSearch::S) where{F<:Real,S<:Integer}
    s=makeEval(mOrig)
    return findBetaIndiffRange(mOrig,s,betaGBounds,betaGPoints,tolG,maxIterG,tolC,maxIterC,xtol,ftol,maxIterSearch)
end

function findBetaIndiffRange(mOrig::longTermBondSpec{F,S},s::longTermBondEval{F,S},betaGBounds::Array{F,1},betaGPoints::S,tolG::F,maxIterG::S,tolC::F,maxIterC::S,xtol::F,ftol::F,maxIterSearch::S) where{F<:Real,S<:Integer}

    @assert betaGBounds[1]<betaGBounds[2]
    @assert betaGBounds[2]<1

    betaGRange=collect(LinRange(betaGBounds[1],betaGBounds[2],betaGPoints))
    indiffVAut=zeros(F,betaGPoints)
    indiffVOpenA0=zeros(F,betaGPoints)
    indiffCBeta=zeros(F,betaGPoints)
    ergodicGainUnexpected=zeros(F,betaGPoints)
    ergodicDMass=zeros(F,betaGPoints)
    ergodicAPDivY=zeros(F,betaGPoints)

    stIncDist=genStIncDist(mOrig,s,tolG,maxIterG)

    workingX=zeros(F,3)
    workingF=zeros(F,3)



    autEVFlow=zeros(mOrig.yParams.nPoints)

    for i in 1:(mOrig.yParams.nPoints)
        autEVFlow[i]=0.0
        for mInd in 1:(mOrig.mParams.nPoints-1)
            autEVFlow[i]+=s.income.mProb[mInd]*u(mOrig,s.income.yGrid[i]+s.income.mMidPoints[mInd])
        end
    end


    autVGrid=zeros(F,mOrig.yParams.nPoints)
    autVSolnMat=zeros(F,mOrig.yParams.nPoints,mOrig.yParams.nPoints)


    for i in 1:betaGPoints
        workingX[1]=betaGRange[i]
        workingX[3]=0.8+0.2*betaGRange[i]

        m=changeBeta(mOrig,betaGRange[i])

        vfiGOneStep!(m,s,tolG,maxIterG)

        autVSolnMat.=inv(Matrix(I,mOrig.yParams.nPoints,mOrig.yParams.nPoints)-betaGRange[i].*s.income.yTMat')
        mul!(autVGrid,autVSolnMat,autEVFlow)

        workingF[1]=(dot(stIncDist,autVGrid)/dot(stIncDist,s.VF.vGrid[s.a0Ind,:]))^(1.0/(1.0-m.gamma))-1.0

        println([0.0,workingX[1],workingF[1]])
        workingF[3]=1.0

        iCount=1

        workingX[3]=0.8+0.2*betaGRange[i]

        wSpec=debtWelfareSpec(workingX[3],m.gamma,1.0,true,false,workingX[3])

        w=makeWelfareEval(m,wSpec,s)

        estimateWelfare!(m,wSpec,s,w,tolC,maxIterC)

        autVSolnMat.=inv(Matrix(I,mOrig.yParams.nPoints,mOrig.yParams.nPoints)-workingX[3].*s.income.yTMat')
        mul!(autVGrid,autVSolnMat,autEVFlow)

        workingF[3]=(dot(stIncDist,autVGrid)/dot(stIncDist,w.VF.vGrid[s.a0Ind,:]))^(1.0/(1.0-m.gamma))-1.0
        println([0.0,workingX[3],workingF[3]])

        while workingF[3]<0.0
            workingX[3]=0.8+0.2*workingX[3]

            wSpec=debtWelfareSpec(workingX[3],m.gamma,1.0,true,false,workingX[3])

            estimateWelfare!(m,wSpec,s,w,tolC,maxIterC)

            autVSolnMat.=inv(Matrix(I,mOrig.yParams.nPoints,mOrig.yParams.nPoints)-workingX[3].*s.income.yTMat')
            mul!(autVGrid,autVSolnMat,autEVFlow)

            workingF[3]=(dot(stIncDist,autVGrid)/dot(stIncDist,w.VF.vGrid[s.a0Ind,:]))^(1.0/(1.0-m.gamma))-1.0
            println([0.0,workingX[3],workingF[3]])

        end

        workingX[2]=0.5*(workingX[1]+workingX[3])


        xOld=workingX[3]

        xDist=xtol+1.0
        fDist=ftol+1.0
        f0Val=0.0


        while ((xDist>xtol)||(fDist>ftol))&&(iCount<=maxIterSearch)

            wSpec=debtWelfareSpec(workingX[2],m.gamma,1.0,true,false,workingX[2])

            estimateWelfare!(m,wSpec,s,w,tolC,maxIterC)


            autVSolnMat.=inv(Matrix(I,mOrig.yParams.nPoints,mOrig.yParams.nPoints)-workingX[2].*s.income.yTMat')

            mul!(autVGrid,autVSolnMat,autEVFlow)


            workingF[2]=(dot(stIncDist,autVGrid)/dot(stIncDist,w.VF.vGrid[s.a0Ind,:]))^(1.0/(1.0-m.gamma))-1.0


            fDist=minimum(abs.(workingF))

            xDist=abs(workingX[2]-xOld)

            xOld=workingX[2]

            getNextX!(workingX,workingF)


            println([betaGRange[i],workingX[2],iCount,fDist,xDist])

            iCount+=1
        end

        indiffVAut[i]=dot(stIncDist,autVGrid)
        indiffVOpenA0[i]=dot(stIncDist,w.VF.vGrid[s.a0Ind,:])
        ergodicGainUnexpected[i]=(dot(stIncDist,autVGrid)/dot(stIncDist,w.VF.vGrid[s.a0Ind,:]))^(1.0/(1.0-m.gamma))
        indiffCBeta[i]=xOld

        dEA,dEP=simulateDistribution(m,s,tolG,maxIterG,true)

        pol=makeMeanPolicies(m,s)

        ergodicDMass[i]=sum(dEP.defaultDist)

        ergodicAPDivY[i]=sum(dEA.repayDist.*pol.apDefPol.repayProb.*pol.apDefPol.apRDivYGrid)/sum(dEA.repayDist.*pol.apDefPol.repayProb)

        println("Government Beta=")
        println(betaGRange[i])
        println("Consumer Beta=")
        println(indiffCBeta[i])
        println("Ergodic Default Mass=")
        println(ergodicDMass[i])
        println("ergodic AP/Y=")
        println(ergodicAPDivY[i])

    end


    return hcat(betaGRange,indiffCBeta,indiffVAut,indiffVOpenA0,ergodicGainUnexpected,ergodicDMass,ergodicAPDivY)

end
