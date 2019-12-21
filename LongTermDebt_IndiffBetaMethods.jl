function changeBeta(m::longTermBondSpec{F,S},newBeta::F) where{F<:Real,S<:Integer}

    return longTermBondSpec(newBeta,m.theta,m.gamma,m.hpen0,m.hpen1,m.R,m.lambda,m.coup,m.aBounds,m.aPoints,m.yParams,m.mParams,m.simplePen,m.mixFacQ,m.mixFacV)
end


function findBetaIndiffRange(mOrig::longTermBondSpec{F,S},betaGBounds::Array{F,1},betaGPoints::S,betaCMax::F,tolG::F,maxIterG::S,tolC::F,maxIterC::S,xtol::F,ftol::F,maxIterSearch::S) where{F<:Real,S<:Integer}
    @assert betaCMax<1
    @assert betaGBounds[1]<betaGBounds[2]
    @assert betaGBounds[2]<1

    betaGRange=collect(linspace(betaGBounds[1],betaGBounds[2],betaGPoints))
    indiffVAut=zeros(F,betaGPoints)
    indiffVOpenA0=zeros(F,betaGPoints)
    indiffCBeta=zeros(F,betaGPoints)
    ergodicGainUnexpected=zeros(F,betaGPoints)
    ergodicDMass=zeros(F,betaGPoints)
    ergodicAPDivY=zeros(F,betaGPoints)

    s=getEval(mOrig)


    stIncDist=genStIncDist(mOrig,s,tolG,maxIterG)

    workingX=zeros(F,3)



    autEVFlow=zeros(mOrig.yParams.yPoints)

    for i in 1:(mOrig.yParams.yPoints)
        autEVFlow[i]=0.0
        for mInd in 1:(mOrig.mParams.mPoints-1)
            autEVFlow[i]+=s.mProb[mInd]*u(mOrig,s.yGrid[i]+s.mMidPoints[mInd])
        end
    end


    autVGrid=zeros(mOrig.yParams.yPoints)
    changeGrid=zeros(mOrig.yParams.yPoints)

    for i in 1:betaGPoints
        workingX[1]=betaGRange[i]
        workingX[3]=betaCMax
        workingX[2]=0.5*(workingX[1]+workingX[3])

        m=changeBeta(mOrig,betaGRange[i])

        vfiGOneStep!(m,s,tolG,maxIterG)

        xDist=xtol+1.0
        fDist=ftol+1.0
        f0Val=0.0

        xOld=betaGRange[i]


        iCount=1

        w0=getWelfareEval(m,s,workingX[2])
        w1=getWelfareEval(w0,workingX[2])


        while ((xDist>xtol)||(fDist>ftol))&&(iCount<=maxIterSearch)


            w1=getWelfareEval(w0,workingX[2])

            estimateWelfare!(m,s,w1,tolC,maxIterC)


            autSolnMat=inv(eye(m.yParams.yPoints)-workingX[2]*s.yTMat)

            mul!(autVGrid,autSolnMat,autEVFlow)

            #changeGrid.=(autVGrid./w1.vGrid[:,s.a0Ind]).^(1.0/(1.0-m.gamma))

            #f0Val=(dot(stIncDist,autVGrid)/dot(stIncDist,w1.vGrid[:,s.a0Ind]))^(1.0/(1.0-m.gamma))
            changeGrid.=autVGrid-w1.vGrid[:,s.a0Ind]


            f0Val=dot(stIncDist,changeGrid)

            fDist=abs(f0Val)

            xDist=abs(workingX[2]-xOld)

            if f0Val==0.0
                xOld=workingX[2]
                break
            elseif f0Val<0.0
                workingX[1]=workingX[2]
            else
                workingX[3]=workingX[2]
            end
            println([betaGRange[i],workingX[2],iCount,f0Val,xDist])
            xOld=workingX[2]
            workingX[2]=0.5*(workingX[1]+workingX[3])
            iCount+=1
            w0=deepcopy(w1)
        end
        indiffVAut[i]=dot(stIncDist,autVGrid)
        indiffVOpenA0[i]=dot(stIncDist,w1.vGrid[:,s.a0Ind])
        ergodicGainUnexpected[i]=(dot(stIncDist,autVGrid)/dot(stIncDist,w1.vGrid[:,s.a0Ind]))^(1.0/(1.0-m.gamma))
        indiffCBeta[i]=xOld

        dEA,dEP=simulateDistribution(m,s,tolG,maxIterG,1000.0,true)

        pol=makeMeanPolicies(m,s)

        ergodicDMass[i]=sum(dEP.defaultDist)

        ergodicAPDivY[i]=sum(dEA.repayDist.*pol.repayProb.*pol.apRDivYGrid)/sum(dEA.repayDist.*pol.repayProb)

        println("Government Beta=")
        println(betaGRange[i])
        println("Consumer Beta=")
        println(indiffCBeta[i])
        println("Erogodic Default Mass=")
        println(ergodicDMass[i])
        println("ergodic AP/Y=")
        println(ergodicAPDivY[i])

    end


    return hcat(betaGRange,indiffCBeta,indiffVAut,indiffVOpenA0,ergodicGainUnexpected,ergodicDMass,ergodicAPDivY)

end



function findBetaIndiffRange(mOrig::longTermBondSpec{F,S},s::longTermBondEval{F,S},betaGBounds::Array{F,1},betaGPoints::S,betaCMax::F,tolG::F,maxIterG::S,tolC::F,maxIterC::S,xtol::F,ftol::F,maxIterSearch::S) where{F<:Real,S<:Integer}
    @assert betaCMax<1
    @assert betaGBounds[1]<betaGBounds[2]
    @assert betaGBounds[2]<1

    betaGRange=collect(linspace(betaGBounds[1],betaGBounds[2],betaGPoints))
    indiffVAut=zeros(F,betaGPoints)
    indiffVOpenA0=zeros(F,betaGPoints)
    indiffCBeta=zeros(F,betaGPoints)
    ergodicGainUnexpected=zeros(F,betaGPoints)
    ergodicDMass=zeros(F,betaGPoints)
    ergodicAPDivY=zeros(F,betaGPoints)

    stIncDist=genStIncDist(mOrig,s,tolG,maxIterG)

    workingX=zeros(F,3)



    autEVFlow=zeros(mOrig.yParams.yPoints)

    for i in 1:(mOrig.yParams.yPoints)
        autEVFlow[i]=0.0
        for mInd in 1:(mOrig.mParams.mPoints-1)
            autEVFlow[i]+=s.mProb[mInd]*u(mOrig,s.yGrid[i]+s.mMidPoints[mInd])
        end
    end


    autVGrid=zeros(mOrig.yParams.yPoints)
    changeGrid=zeros(mOrig.yParams.yPoints)

    for i in 1:betaGPoints
        workingX[1]=betaGRange[i]
        workingX[3]=betaCMax
        workingX[2]=0.5*(workingX[1]+workingX[3])

        m=changeBeta(mOrig,betaGRange[i])

        vfiGOneStep!(m,s,tolG,maxIterG)

        xDist=xtol+1.0
        fDist=ftol+1.0
        f0Val=0.0

        xOld=betaGRange[i]


        iCount=1

        w0=getWelfareEval(m,s,workingX[2])
        w1=getWelfareEval(w0,workingX[2])


        while ((xDist>xtol)||(fDist>ftol))&&(iCount<=maxIterSearch)


            w1=getWelfareEval(w0,workingX[2])

            estimateWelfare!(m,s,w1,tolC,maxIterC)


            autSolnMat=inv(eye(m.yParams.yPoints)-workingX[2]*s.yTMat)

            mul!(autVGrid,autSolnMat,autEVFlow)

            #changeGrid.=(autVGrid./w1.vGrid[:,s.a0Ind]).^(1.0/(1.0-m.gamma))

            #f0Val=(dot(stIncDist,autVGrid)/dot(stIncDist,w1.vGrid[:,s.a0Ind]))^(1.0/(1.0-m.gamma))
            changeGrid.=autVGrid-w1.vGrid[:,s.a0Ind]


            f0Val=dot(stIncDist,changeGrid)

            fDist=abs(f0Val)

            xDist=abs(workingX[2]-xOld)

            if f0Val==0.0
                xOld=workingX[2]
                break
            elseif f0Val<0.0
                workingX[1]=workingX[2]
            else
                workingX[3]=workingX[2]
            end
            println([betaGRange[i],workingX[2],iCount,f0Val,xDist])
            xOld=workingX[2]
            workingX[2]=0.5*(workingX[1]+workingX[3])
            iCount+=1
            w0=deepcopy(w1)
        end
        indiffVAut[i]=dot(stIncDist,autVGrid)
        indiffVOpenA0[i]=dot(stIncDist,w1.vGrid[:,s.a0Ind])
        ergodicGainUnexpected[i]=(dot(stIncDist,autVGrid)/dot(stIncDist,w1.vGrid[:,s.a0Ind]))^(1.0/(1.0-m.gamma))
        indiffCBeta[i]=xOld

        dEA,dEP=simulateDistribution(m,s,tolG,maxIterG,1000.0,true)

        pol=makeMeanPolicies(m,s)

        ergodicDMass[i]=sum(dEP.defaultDist)

        ergodicAPDivY[i]=sum(dEA.repayDist.*pol.repayProb.*pol.apRDivYGrid)/sum(dEA.repayDist.*pol.repayProb)

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
