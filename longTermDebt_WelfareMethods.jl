#Will be fully commented shortly



#Define some types
struct debtDist{F<:Real}
    repayDist::Array{F,2}
    defaultDist::Array{F,1}
end

struct debtWelfareSpec{F<:Real}
    betaC::F
    gammaC::F
    penMultC::F
    penM::Bool
    useVGPostDef::Bool
    betaCPostDef::F
end

struct debtWelfareEval{F<:Real}
    yGridC::Array{F,1}
    yDefGridC::Array{F,1}
    netRevM0A0GridC::Array{F,2}
    EVFlow::Array{F,2}
    EVA0::Array{F,1}
    VF::vFunc{F}
end

struct welfareResults{F<:Real}
    summaryVals::Array{F,2}
    autVGridAll::Array{F,2}
    openA0VGrid::Array{F,2}
    gainA0Grid::Array{F,2}
end

struct welfareRangeParams{F<:Real,S<:Integer}
    betaCBounds::Array{F,1}
    betaCPoints::S
    gammaCBounds::Array{F,1}
    gammaCPoints::S
    penMultCBounds::Array{F,1}
    penMultCPoints::S
    produceWithMPen::Bool
    produceWithoutMPen::Bool
    produceWithBetaGPostDef::Bool
    produceWithoutBetaGPostDef::Bool
end



#c is simply a convenience function for retrieving (the argument x denotes the value of the m shock)
#when dMark is false (so the government is not in default):
#1. s.netRevM0A0Grid[j,i]-s.qGrid[apInd,i]*s.aGridIncr[apInd,j]+x
#and when dMark is true (so the government is in default):
#2. s.income.yDefGrid[i]+x


function c(w::debtWelfareEval{F},s::longTermBondEval{F,S},i::S,j::S,apInd::S,x::U,dMark::Bool) where{F<:Real,U<:Real,S<:Integer}
    if dMark==false
        return w.netRevM0A0GridC[j,i]-s.qGrid[apInd,i]*s.aGridIncr[apInd,j]+x
    else
        return w.yDefGridC[i]+x
    end
end


#u is the CRRA utility function. Its arguments are just:
#1. m: the model specification
#2. x: the value of consumption
function u(wSpec::debtWelfareSpec{F},x::U) where{F<:Real,U<:Real}
    #If x is positive return the value of CRRA utility
    if x>0.0
        if wSpec.gammaC!=one(F)
            return x^(1.0-wSpec.gammaC)/(1.0-wSpec.gammaC)
        else
            return log(x)
        end
    #Otherwise return either a very negative value, scaled by how negative x is when utility is negative near x=0
    #or a large negative value which gets smaller as x rises to 0.
    #This is done to ensure that any optimization algorithm knows to try to make consumption positive at essentially any cost.
    else
        if wSpec.gammaC>=one(F)
            return u(m,1e-10)*(1.0+abs(x))
        else
            return -1e10*abs(x)
        end
    end
end


#Function to create welfare object.
function makeWelfareEval(m::longTermBondSpec{F,S},wSpec::debtWelfareSpec{F},s::longTermBondEval{F,S}) where{F<:Real,S<:Integer}
    time1=time()
    VF=deepcopy(s.VF)

    yGridC=deepcopy(s.income.yGrid)

    yDefPenC=wSpec.penMultC*(s.income.yGrid.-s.income.yDefGrid)
    yDefGridC=yGridC.-yDefPenC
    @assert minimum(yDefGridC)>-s.income.mBounds[1]

    if wSpec.useVGPostDef==false
        for i in 1:(m.yParams.nPoints)
            VF.vDFutFlow[i]=0.0
            VF.vDInitFlow[i]=u(wSpec,yDefGridC[i]+s.income.mBounds[1])
            for k in 1:(m.mParams.nPoints-1)
                VF.vDFutFlow[i]+=s.income.mProb[k]*u(wSpec,yDefGridC[i]+s.income.mMidPoints[k])
            end
        end

    else
        for i in 1:(m.yParams.nPoints)
            VF.vDFutFlow[i]=0.0
            VF.vDInitFlow[i]=u(m,yDefGridC[i]+s.income.mBounds[1])
            for k in 1:(m.mParams.nPoints-1)
                VF.vDFutFlow[i]+=s.income.mProb[k]*u(m,yDefGridC[i]+s.income.mMidPoints[k])
            end
        end

    end
    if wSpec.penM==false
        VF.vDInitFlow.=VF.vDFutFlow
    end
    time2=time()
    netRevM0A0GridC=deepcopy(s.netRevM0A0Grid)
    time3=time()
    EVFlow=zeros(m.aPoints,m.yParams.nPoints)
    EVA0=zeros(m.yParams.nPoints)
    EVA0.=VF.EVGrid[s.a0Ind,:]



    w=debtWelfareEval(yGridC,yDefGridC,netRevM0A0GridC,EVFlow,EVA0,VF)
    time4=time()
    for i in 1:(m.yParams.nPoints)
        for j in 1:(m.aPoints)
            w.EVFlow[j,i]=0.0
            if s.pol.alwaysDefault[j,i]==true
                w.EVFlow[j,i]+=w.VF.vDInitFlow[i]
            elseif s.pol.neverDefault[j,i]==true
                #If the government never defaults, then we begin the integration over the result under repayment at the first interval.
                #We set the location of the current upper bound in the m grid to the second point, the index of the relevant borrowing policy entry to the first,
                #and the most recent value of m to its lower bound
                mUBInd=2
                aPolicyInd=1
                lastMVal=s.income.mGrid[1]
                #Loop until we exit the upper boundary of the m space (the second check is to ensure that we never exit the range of the current borrowing policy function;
                #it should never matter, and if it does, the process should in general not lead to convergence).
                while (mUBInd<=(m.mParams.nPoints))&&(aPolicyInd<=(s.pol.mListLength[j,i]))
                    #If the upper bound for the range in which the current entry of the borrowing policy function is valid is strictly less than the upper bound of the current
                    #m interval, add the relevant contribution of that entry, update the last value of m reached, and increment only the variable marking our position in the borrowing
                    #policy function
                    if s.pol.mAPThreshold[i][aPolicyInd,j]<(s.income.mGrid[mUBInd])
                        w.EVFlow[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*u(wSpec,c(w,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1],false))
                        lastMVal=s.pol.mAPThreshold[i][aPolicyInd,j]
                        aPolicyInd+=1
                    else
                        #Otherwise, add the relevant contribution of the current entry which lies in this m interval, update the last value of m reached to be its upper bound,
                        #and then increment only the variable marking our position in the m grid
                        w.EVFlow[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*u(wSpec,c(w,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1],false))
                        lastMVal=s.income.mGrid[mUBInd]
                        mUBInd+=1
                    end
                end
            else
                #If the government does default only sometimes, find the position in the m grid of the upper bound of the interval in which the threshold level of m for default falls
                #and set the last value of m observed to the value directly preceding that upper bound
                mUBInd=searchsortedfirst(s.income.mGrid,s.pol.defThreshold[j,i])
                lastMVal=s.income.mGrid[mUBInd-1]
                #If there are any entire intervals in which default occurs, add their contribution to the output value
                if mUBInd>2
                    if wSpec.penM==true
                        for k in 3:mUBInd
                            w.EVFlow[j,i]+=s.income.mProb[k-2]*w.VF.vDInitFlow[i]
                        end
                    elseif wSpec.useVGPostDef==true
                        for k in 3:mUBInd
                            w.EVFlow[j,i]+=s.income.mProb[k-2]*u(m,w.yDefGridC[i]+s.income.mMidPoints[k-2])
                        end
                    else
                        for k in 3:mUBInd
                            w.EVFlow[j,i]+=s.income.mProb[k-2]*u(wSpec,w.yDefGridC[i]+s.income.mMidPoints[k-2])
                        end
                    end
                end
                #Set the index of the first relevant entry of the borowing policy function
                aPolicyInd=s.pol.firstRepayInd[j,i]
                #Add the contribution of default in the interval in which the threshold lies to the output value
                if wSpec.penM==true
                    w.EVFlow[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes*w.VF.vDInitFlow[i]
                elseif wSpec.useVGPostDef==true
                    w.EVFlow[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes*u(m,w.yDefGridC[i]+s.income.mMidPoints[mUBInd-1])
                else
                    w.EVFlow[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes*u(wSpec,w.yDefGridC[i]+s.income.mMidPoints[mUBInd-1])
                end
                #Update the last value of m observed to the threshold level of m at which default occurs
                lastMVal=s.pol.defThreshold[j,i]
                #Loop until we exit the upper boundary of the m space (the second check is to ensure that we never exit the range of the current borrowing policy function;
                #it should never matter, and if it does, the process should in general not lead to convergence).
                while (mUBInd<=(m.mParams.nPoints))&&(aPolicyInd<=(s.pol.mListLength[j,i]))
                    #If the upper bound for the range in which the current entry of the borrowing policy function is valid is strictly less than the upper bound of the current
                    #m interval, add the relevant contribution of that entry, update the last value of m reached, and increment only the variable marking our position in the borrowing
                    #policy function
                    if s.pol.mAPThreshold[i][aPolicyInd,j]<(s.income.mGrid[mUBInd])
                        w.EVFlow[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*u(wSpec,c(w,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1],false))
                        lastMVal=s.pol.mAPThreshold[i][aPolicyInd,j]
                        aPolicyInd+=1
                    else
                        #Otherwise, add the relevant contribution of the current entry which lies in this m interval, update the last value of m reached to be its upper bound,
                        #and then increment only the variable marking our position in the m grid
                        w.EVFlow[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*u(wSpec,c(w,s,i,j,s.pol.apPolicy[i][aPolicyInd,j],s.income.mMidPoints[mUBInd-1],false))
                        lastMVal=s.income.mGrid[mUBInd]
                        mUBInd+=1
                    end
                end
            end
        end
    end

    time5=time()
    println([time2-time1,time3-time2,time4-time3,time5-time4])

    return w
end


#Function to modify existing welfare objects
function makeWelfareEval(m::longTermBondSpec{F,S},s::longTermBondEval{F,S},w::debtWelfareEval{F},wSpecNew::debtWelfareSpec{F},wSpecOld::debtWelfareSpec{F}) where{F<:Real,S<:Integer}

    if (wSpecOld.gammaC!=wSpecNew.gammaC)||(wSpecOld.penMultC!=wSpecNew.penMultC)||(wSpecOld.penM!=wSpecNew.penM)||(wSpecOld.useVGPostDef!=wSpecNew.useVGPostDef)

        wOut=makeWelfareEval(m,wSpecNew,s)

        copyto!(wOut.VF.EVGrid,w.VF.EVGrid)
        copyto!(wOut.VF.vGrid,w.VF.vGrid)
        #if wSpecNew.useVGPostDef==false
            copyto!(wOut.VF.EVDGrid,w.VF.EVDGrid)
            copyto!(wOut.VF.vDInitGrid,w.VF.vDInitGrid)
            copyto!(wOut.VF.vDFutGrid,w.VF.vDFutGrid)
            copyto!(wOut.EVA0,w.EVA0)
        #end
        return wOut
    else
        wOut=deepcopy(w)

        inflFac=(1.0-wSpecOld.betaC)/(1.0-wSpecNew.betaC)

        wOut.VF.EVGrid.*=inflFac

        wOut.VF.vGrid.*=inflFac
        if wSpecNew.useVGPostDef==false
            wOut.VF.EVDGrid.*=inflFac
            wOut.VF.vDInitGrid.*=inflFac
            wOut.VF.vDFutGrid.*=inflFac
            wOut.EVA0.*=inflFac
        end


        return wOut
    end


end




#integrateMVApprox performs the integration step in V(y,b)=E_m[W(y,m,b)] for a specific value of the persistent component of income and incoming debt, given government policies.
#This function performs the integration exactly as described in Chatterjee and Eyigungor (2012). Its arguments are:
#1. m: a model specification
#2. wSpec: a consumer welfare specification
#3. s: a collection of objects used to solve the model
#4. w: a consumer welfare object
#5. i: the index of the persistent component of income
#6. j: the index of the incoming level of borrowing
#It returns a single number which is the value of the integral.

#The structure and internal logic of this function are essentially identical to those of the one directly above it.
function integrateMVApprox(m::longTermBondSpec{F,S},wSpec::debtWelfareSpec{F},s::longTermBondEval{F,S},w::debtWelfareEval{F},i::S,j::S) where{F<:Real,S<:Integer}
    outValV=0.0
    if s.pol.alwaysDefault[j,i]==true
        outValV+=w.VF.vDInitGrid[i]
    elseif s.pol.neverDefault[j,i]==true
        outValV+=w.EVFlow[j,i]
        #If the government never defaults, then we begin the integration over the result under repayment at the first interval.
        #We set the location of the current upper bound in the m grid to the second point, the index of the relevant borrowing policy entry to the first,
        #and the most recent value of m to its lower bound
        mUBInd=2
        aPolicyInd=1
        lastMVal=s.income.mGrid[1]
        #Loop until we exit the upper boundary of the m space (the second check is to ensure that we never exit the range of the current borrowing policy function;
        #it should never matter, and if it does, the process should in general not lead to convergence).
        while (mUBInd<=(m.mParams.nPoints))&&(aPolicyInd<=(s.pol.mListLength[j,i]))
            #If the upper bound for the range in which the current entry of the borrowing policy function is valid is strictly less than the upper bound of the current
            #m interval, add the relevant contribution of that entry, update the last value of m reached, and increment only the variable marking our position in the borrowing
            #policy function
            if s.pol.mAPThreshold[i][aPolicyInd,j]<(s.income.mGrid[mUBInd])
                outValV+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*wSpec.betaC*w.VF.EVGrid[s.pol.apPolicy[i][aPolicyInd,j],i]
                lastMVal=s.pol.mAPThreshold[i][aPolicyInd,j]
                aPolicyInd+=1
            else
                #Otherwise, add the relevant contribution of the current entry which lies in this m interval, update the last value of m reached to be its upper bound,
                #and then increment only the variable marking our position in the m grid
                outValV+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*wSpec.betaC*w.VF.EVGrid[s.pol.apPolicy[i][aPolicyInd,j],i]
                lastMVal=s.income.mGrid[mUBInd]

                mUBInd+=1
            end
        end

    else
        outValV+=w.EVFlow[j,i]
        #If the government does default only sometimes, find the position in the m grid of the upper bound of the interval in which the threshold level of m for default falls
        #and set the last value of m observed to the value directly preceding that upper bound
        mUBInd=searchsortedfirst(s.income.mGrid,s.pol.defThreshold[j,i])
        lastMVal=s.income.mGrid[mUBInd-1]
        #If there are any entire intervals in which default occurs, add their contribution to the output value
        if mUBInd>2
            for k in 3:mUBInd
                outValV+=s.income.mProb[k-2]*wSpec.betaCPostDef*w.VF.EVDGrid[i]

            end
        end
        #Set the index of the first relevant entry of the borowing policy function
        aPolicyInd=s.pol.firstRepayInd[j,i]

        #Add the contribution of default in the interval in which the threshold lies to the output value
        outValV+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes*wSpec.betaCPostDef*w.VF.EVDGrid[i]

        #Update the last value of m observed to the threshold level of m at which default occurs
        lastMVal=s.pol.defThreshold[j,i]
        #Loop until we exit the upper boundary of the m space (the second check is to ensure that we never exit the range of the current borrowing policy function;
        #it should never matter, and if it does, the process should in general not lead to convergence).
        while (mUBInd<=(m.mParams.nPoints))&&(aPolicyInd<=(s.pol.mListLength[j,i]))
            #If the upper bound for the range in which the current entry of the borrowing policy function is valid is strictly less than the upper bound of the current
            #m interval, add the relevant contribution of that entry, update the last value of m reached, and increment only the variable marking our position in the borrowing
            #policy function
            if s.pol.mAPThreshold[i][aPolicyInd,j]<(s.income.mGrid[mUBInd])
                outValV+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*wSpec.betaC*w.VF.EVGrid[s.pol.apPolicy[i][aPolicyInd,j],i]
                lastMVal=s.pol.mAPThreshold[i][aPolicyInd,j]
                aPolicyInd+=1
            else
                #Otherwise, add the relevant contribution of the current entry which lies in this m interval, update the last value of m reached to be its upper bound,
                #and then increment only the variable marking our position in the m grid
                outValV+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*wSpec.betaC*w.VF.EVGrid[s.pol.apPolicy[i][aPolicyInd,j],i]
                lastMVal=s.income.mGrid[mUBInd]
                mUBInd+=1
            end
        end
    end

    return outValV
end

#updateVD! is simple function which just updates the all the default value functions. Its arguments are:
#The arguments of the second are:
#1. m: a model specification
#2. wSpec: a consumer welfare specification
#3. s: a collection of objects used to solve the model
#4. w: a consumer welfare object
#5. VFNew: a value function object in which the result is to be stored
function updateVD!(m::longTermBondSpec{F,S},wSpec::debtWelfareSpec{F},s::longTermBondEval{F,S},w::debtWelfareEval{F},VFNew::vFunc{F}) where{F<:Real,S<:Integer}

    VFNew.vDFutGrid.=w.VF.vDFutFlow.+wSpec.betaCPostDef.*m.theta.*w.EVA0.+wSpec.betaCPostDef.*(1.0-m.theta).*s.income.yTMat'*w.VF.vDFutGrid
    VFNew.EVDGrid.=m.theta.*w.EVA0.+(1.0-m.theta).*s.income.yTMat'*VFNew.vDFutGrid
    VFNew.vDInitGrid.=w.VF.vDInitFlow.+wSpec.betaCPostDef.*VFNew.EVDGrid

    return VFNew
end


#updateV! performs the full integration across states to calculate V(y,b)=E_m[W(y,m,b)]. Its arguments are:
#1. m: a model specification
#2. wSpec: a consumer welfare specification
#3. s: a collection of objects used to solve the model
#4. w: a consumer welfare object
#5. VFNew: a value function object in which the result is to be stored

#Its output is simply a modified version of VFNew
function updateV!(m::longTermBondSpec{F,S},wSpec::debtWelfareSpec{F},s::longTermBondEval{F,S},w::debtWelfareEval{F},VFNew::vFunc{F})  where{F<:Real,S<:Integer}
    Threads.@threads for j in 1:(m.aPoints)
    #for j in 1:(m.aPoints)
        for i in 1:(m.yParams.nPoints)
            VFNew.vGrid[j,i]=integrateMVApprox(m,wSpec,s,w,i,j)
        end
    end

    return VFNew
end


#Updates the continuation value function under repayment
function updateEV!(m::longTermBondSpec{F,S},s::longTermBondEval{F,S},VFNew::vFunc{F})  where{F<:Real,S<:Integer}
    mul!(VFNew.EVGrid,VFNew.vGrid,s.income.yTMat)
    return VFNew
end

#Iterates the contraction mapping for the consumer value functions until convergence
function estimateWelfare!(m::longTermBondSpec{F,S},wSpec::debtWelfareSpec{F},s::longTermBondEval{F,S},w::debtWelfareEval{F},tol::F,maxIter::S) where{F<:Real,S<:Integer}

    VFNew=deepcopy(w.VF)


    #Initialize  a counter for the number of iterations and measures of sup norm distance between various objects
    iCount=1
    vFDist=tol+1.0
    EVDist=tol+1.0
    EVDDist=tol+1.0
    vDFutDist=tol+1.0
    vDInitDist=tol+1.0
    maxDist=tol+1.0

    while (iCount<=maxIter)&(maxDist>=tol)
        #Update the default value functions and calculate the distance between the update and the previous version
        updateVD!(m,wSpec,s,w,VFNew)
        vDFutDist=maximum(abs.(VFNew.vDFutGrid.-w.VF.vDFutGrid))
        EVDDist=maximum(abs.(VFNew.EVDGrid.-w.VF.EVDGrid))
        vDInitDist=maximum(abs.(VFNew.vDInitGrid.-w.VF.vDInitGrid))

        w.VF.vDFutGrid.=VFNew.vDFutGrid
        w.VF.vDInitGrid.=VFNew.vDInitGrid


        #Perform the integration step to obtain new versions of V and Z
        updateV!(m,wSpec,s,w,VFNew)
        updateEV!(m,s,VFNew)

        #Calculate the distance between successive iterations for V and Z
        vFDist=maximum(abs.(VFNew.vGrid.-w.VF.vGrid))
        EVDist=maximum(abs.(VFNew.EVGrid.-w.VF.EVGrid))

        #Update the guess for Z, Z^D, and V
        w.VF.EVDGrid.=VFNew.EVDGrid
        w.VF.EVGrid.=VFNew.EVGrid
        w.VF.vGrid.=VFNew.vGrid

        if wSpec.useVGPostDef==false
            w.EVA0.=w.VF.EVGrid[s.a0Ind,:]
        end

        maxDist=max(vFDist,EVDist,EVDDist,vDFutDist,vDInitDist)

        if div(iCount,max(div(maxIter,20),1))==(iCount/max(div(maxIter,20),1))
            println([iCount,EVDist,EVDDist,vFDist,vDFutDist,vDInitDist])
        end
        #Increment the iteration counter
        iCount+=1


    end
    println([iCount-1,EVDist,EVDDist,vFDist,vDFutDist,vDInitDist])

    return w
end

#Calculates the value of autarky
function estimateAutarkyFlow(m::longTermBondSpec{F,S},wSpec::debtWelfareSpec{F},s::longTermBondEval{F,S}) where{F<:Real,S<:Integer}

    autEVFlow=zeros(m.yParams.nPoints)
    for i in 1:(m.yParams.nPoints)
        for mInd in 1:(m.mParams.nPoints-1)
            autEVFlow[i]+=s.income.mProb[mInd]*u(wSpec,s.income.yGrid[i]+s.income.mMidPoints[mInd])
        end
    end

    return autEVFlow
end

function makeWelfareSpecList(m::longTermBondSpec{F,S},wParams::welfareRangeParams{F,S}) where{F<:Real,S<:Integer}
    betaCRange=ifelse(wParams.betaCPoints!=1,collect(LinRange(wParams.betaCBounds[1],wParams.betaCBounds[2],wParams.betaCPoints)),[wParams.betaCBounds[1]])
    gammaCRange=ifelse(wParams.gammaCPoints!=1,collect(LinRange(wParams.gammaCBounds[1],wParams.gammaCBounds[2],wParams.gammaCPoints)),[wParams.gammaCBounds[1]])
    penMultCRange=ifelse(wParams.penMultCPoints!=1,collect(LinRange(wParams.penMultCBounds[1],wParams.penMultCBounds[2],wParams.penMultCPoints)),[wParams.penMultCBounds[1]])

    paramPoints=wParams.betaCPoints*wParams.gammaCPoints*wParams.penMultCPoints
    dupMPenMult=(wParams.produceWithMPen+wParams.produceWithoutMPen)
    dupPostDBetaMult=(wParams.produceWithBetaGPostDef+wParams.produceWithoutBetaGPostDef)
    dupParamMult=dupMPenMult*dupPostDBetaMult
    paramPointsTot=paramPoints*dupParamMult
    paramGrid=zeros(F,paramPointsTot,6)

    if (wParams.produceWithMPen==true)&&(wParams.produceWithoutMPen==true)
        paramGrid[1:dupPostDBetaMult*paramPoints,4].=1.0
    elseif (wParams.produceWithMPen==true)
        paramGrid[:,4].=1.0
    end


    if (wParams.produceWithBetaGPostDef==true)&&(wParams.produceWithoutBetaGPostDef==true)
        paramGrid[1:paramPoints,5].=1.0
        if (wParams.produceWithMPen==true)&&(wParams.produceWithoutMPen==true)
            paramGrid[2*paramPoints+1:3*paramPoints,5].=1.0
        end
    elseif (wParams.produceWithBetaGPostDef==true)
        paramGrid[:,5].=1.0
    end

    tempParamInd=1
    for gammaInd in 1:wParams.gammaCPoints
        for penMultInd in 1:wParams.penMultCPoints
            for betaInd in 1:wParams.betaCPoints
                paramGrid[tempParamInd,1]=betaCRange[betaInd]
                paramGrid[tempParamInd,2]=gammaCRange[gammaInd]
                paramGrid[tempParamInd,3]=penMultCRange[penMultInd]
                tempParamInd+=1
            end
        end
    end
    if dupParamMult>1
        for dupNum in 2:dupParamMult
            paramGrid[(dupNum-1)*paramPoints+1:dupNum*paramPoints,1:3].=paramGrid[1:paramPoints,1:3]
        end
    end
    for paramInd in 1:paramPointsTot
        if paramGrid[paramInd,5]==1.0
            paramGrid[paramInd,6]=m.beta
        else
            paramGrid[paramInd,6]=paramGrid[paramInd,1]
        end
    end
    vDParamPointsTot=dupMPenMult*wParams.penMultCPoints
    vDParamGrid=zeros(vDParamPointsTot,2)
    if (wParams.produceWithMPen==true)
        vDParamGrid[1:wParams.penMultCPoints,2].=1.0
    end
    vDParamGrid[1:wParams.penMultCPoints,1].=penMultCRange
    if dupMPenMult>1
        vDParamGrid[wParams.penMultCPoints+1:2*wParams.penMultCPoints,1].=penMultCRange
    end
    vDSpecMap=zeros(S,paramPointsTot)

    for i in 1:paramPointsTot
        vDSpecMap[i]=findfirst((paramGrid[i,3].==vDParamGrid[:,1]).*(paramGrid[i,4].==vDParamGrid[:,2]))
    end

    welfareSpecList=[debtWelfareSpec(paramGrid[specInd,1],paramGrid[specInd,2],paramGrid[specInd,3],Bool(paramGrid[specInd,4]),Bool(paramGrid[specInd,5]),paramGrid[specInd,6]) for specInd in 1:paramPointsTot]
    vDSpecList=[debtWelfareSpec(m.beta,m.gamma,vDParamGrid[specInd,1],Bool(vDParamGrid[specInd,2]),false,m.beta) for specInd in 1:vDParamPointsTot]

    return paramGrid,welfareSpecList,paramPointsTot,vDSpecMap,vDSpecList,vDParamPointsTot
end

#Calculates consumer welfare across a range of parameters
function getWelfareRange(m::longTermBondSpec{F,S},s::longTermBondEval{F,S},wParams::welfareRangeParams{F,S},allOrNoPen::Bool,tol::F,maxIter::S) where{F<:Real,S<:Integer}
    @assert wParams.betaCBounds[1]<=wParams.betaCBounds[2]
    @assert wParams.gammaCBounds[1]<=wParams.gammaCBounds[2]
    @assert wParams.penMultCBounds[1]<=wParams.penMultCBounds[2]
    @assert wParams.betaCBounds[1]>0.0
    @assert wParams.betaCBounds[2]<1.0
    @assert wParams.gammaCBounds[1]>=0.0
    @assert wParams.penMultCBounds[1]>=0.0
    @assert wParams.betaCPoints>=1
    @assert wParams.gammaCPoints>=1
    @assert wParams.penMultCPoints>=1
    @assert (wParams.produceWithMPen==true)||(wParams.produceWithoutMPen==true)
    @assert (wParams.produceWithBetaGPostDef==true)||(wParams.produceWithoutBetaGPostDef==true)
    if wParams.betaCPoints==1
        @assert wParams.betaCBounds[1]==wParams.betaCBounds[2]
    end
    if wParams.gammaCPoints==1
        @assert wParams.gammaCBounds[1]==wParams.gammaCBounds[2]
    end
    if wParams.penMultCPoints==1
        @assert wParams.penMultCBounds[1]==wParams.penMultCBounds[2]
    end
    if allOrNoPen==true
        @assert (wParams.penMultCPoints==1)||(wParams.penMultCPoints==2)
        if (wParams.penMultCPoints==2)
            @assert (wParams.penMultCBounds[1]==0.0)&&(wParams.penMultCBounds[2]==1.0)
        else
            @assert ((wParams.penMultCBounds[1]==0.0)&&(wParams.produceWithoutMPen==true))||((wParams.penMultCBounds[1]==1.0)&&(wParams.produceWithMPen==true))
        end
    end
    paramGrid,welfareSpecList,paramPointsTot,vDSpecMap,vDSpecList,vDParamPointsTot=makeWelfareSpecList(m,wParams)

    if allOrNoPen==true
        keepLocs=paramGrid[:,3].==paramGrid[:,4]
        paramGrid=paramGrid[keepLocs,:]
        welfareSpecList=welfareSpecList[keepLocs]
        paramPointsTot=sum(keepLocs)
        vDSpecMap=vDSpecMap[keepLocs]
    end

    vDFutGridG=zeros(m.yParams.nPoints,vDParamPointsTot)
    vDInitGridG=zeros(m.yParams.nPoints,vDParamPointsTot)
    EVDGridG=zeros(m.yParams.nPoints,vDParamPointsTot)
    EVA0G=zeros(m.yParams.nPoints,vDParamPointsTot)

    wVD=makeWelfareEval(m,vDSpecList[1],s)
    estimateWelfare!(m,vDSpecList[1],s,wVD,tol^2,2*maxIter)

    EVDGridG[:,1].=wVD.VF.EVDGrid
    vDFutGridG[:,1].=wVD.VF.vDFutGrid
    vDInitGridG[:,1].=wVD.VF.vDInitGrid
    EVA0G[:,1].=wVD.VF.EVGrid[s.a0Ind,:]

    for k in 2:vDParamPointsTot
        wVD=makeWelfareEval(m,vDSpecList[k],s)
        estimateWelfare!(m,vDSpecList[k],s,wVD,tol^2,2*maxIter)

        EVDGridG[:,k].=wVD.VF.EVDGrid
        vDFutGridG[:,k].=wVD.VF.vDFutGrid
        vDInitGridG[:,k].=wVD.VF.vDInitGrid
        EVA0G[:,k].=wVD.VF.EVGrid[s.a0Ind,:]
    end

    changeGrid=zeros(m.yParams.nPoints)

    stJointDist=simulateDistribution(m,s,tol,maxIter,false)
    stIncDist=genStIncDist(m,s,tol,maxIter)



    meanIncInd=Int64(floor(m.yParams.nPoints/2))
    minIncInd=1
    maxIncInd=m.yParams.nPoints

    meanIncGain=zeros(paramPointsTot)
    minIncGain=zeros(paramPointsTot)
    maxIncGain=zeros(paramPointsTot)
    meanGain=zeros(paramPointsTot)
    minGain=zeros(paramPointsTot)
    maxGain=zeros(paramPointsTot)

    ergodicGainActual=zeros(paramPointsTot)
    ergodicGainUnexpected=zeros(paramPointsTot)

    autVGridAll=zeros(m.yParams.nPoints,paramPointsTot)
    openA0VGrid=zeros(m.yParams.nPoints,paramPointsTot)
    gainA0Grid=zeros(m.yParams.nPoints,paramPointsTot)



    w0=makeWelfareEval(m,welfareSpecList[1],s)
    if welfareSpecList[1].useVGPostDef==true
        w0.VF.EVDGrid.=EVDGridG[:,vDSpecMap[1]]
        w0.VF.vDFutGrid.=vDFutGridG[:,vDSpecMap[1]]
        w0.VF.vDInitGrid.=vDInitGridG[:,vDSpecMap[1]]
        w0.EVA0.=EVA0G[:,vDSpecMap[1]]
    end

    estimateWelfare!(m,welfareSpecList[1],s,w0,tol,maxIter)

    autVGrid=zeros(m.yParams.nPoints)
    autEVFlow=estimateAutarkyFlow(m,welfareSpecList[1],s)
    autSolnMat=inv(Matrix(I,m.yParams.nPoints,m.yParams.nPoints)-welfareSpecList[1].betaC*s.income.yTMat')

    mul!(autVGrid,autSolnMat,autEVFlow)


    if welfareSpecList[1].gammaC!=1.0
        changeGrid.=(autVGrid./w0.VF.vGrid[s.a0Ind,:]).^(1.0/(1.0-welfareSpecList[1].gammaC))
    else
        changeGrid.=exp.((1.0-welfareSpecList[1].betaC)*(autVGrid.-w0.VF.vGrid[s.a0Ind,:]))
    end
    autVGridAll[:,1].=autVGrid
    openA0VGrid[:,1].=w0.VF.vGrid[s.a0Ind,:]

    meanIncGain[1]=changeGrid[meanIncInd]
    minIncGain[1]=changeGrid[minIncInd]
    maxIncGain[1]=changeGrid[maxIncInd]

    meanGain[1]=mean(changeGrid)
    minGain[1]=minimum(changeGrid)
    maxGain[1]=maximum(changeGrid)

    if welfareSpecList[1].gammaC!=1.0
        ergodicGainActual[1]=(dot(stIncDist,autVGrid)/(sum(stJointDist.repayDist.*w0.VF.vGrid)+dot(stJointDist.defaultDist,w0.VF.vDFutGrid)))^(1.0/(1.0-welfareSpecList[1].gammaC))
        ergodicGainUnexpected[1]=(dot(stIncDist,autVGridAll[:,1])/dot(stIncDist,openA0VGrid[:,1]))^(1.0/(1.0-welfareSpecList[1].gammaC))
    else
        ergodicGainActual[1]=exp((1.0-welfareSpecList[1].betaC)*(dot(stIncDist,autVGrid)-(sum(stJointDist.repayDist.*w0.VF.vGrid)+dot(stJointDist.defaultDist,w0.VF.vDFutGrid))))
        ergodicGainUnexpected[1]=exp((1.0-welfareSpecList[1].betaC)*(dot(stIncDist,autVGridAll[:,1])-dot(stIncDist,openA0VGrid[:,1])))
    end

    println([1,paramGrid[1,1],paramGrid[1,2],paramGrid[1,3],paramGrid[1,4],paramGrid[1,5],meanIncGain[1],minIncGain[1],maxIncGain[1],meanGain[1],minGain[1],maxGain[1],ergodicGainActual[1],ergodicGainUnexpected[1]])
    w1=deepcopy(w0)

    for k in 2:paramPointsTot
        w1=makeWelfareEval(m,s,w0,welfareSpecList[k],welfareSpecList[k-1])

        if welfareSpecList[k].useVGPostDef==true
            w1.VF.EVDGrid.=EVDGridG[:,vDSpecMap[k]]
            w1.VF.vDFutGrid.=vDFutGridG[:,vDSpecMap[k]]
            w1.VF.vDInitGrid.=vDInitGridG[:,vDSpecMap[k]]
            w1.EVA0.=EVA0G[:,vDSpecMap[k]]
        end

        #if k==16
        #    return welfareSpecList[k],w1
        #end

        estimateWelfare!(m,welfareSpecList[k],s,w1,tol,maxIter)
        if welfareSpecList[k].gammaC!=welfareSpecList[k-1].gammaC
            autEVFlow.=estimateAutarkyFlow(m,welfareSpecList[k],s)
        end
        if welfareSpecList[k].betaC!=welfareSpecList[k-1].betaC
            autSolnMat.=inv(Matrix(I,m.yParams.nPoints,m.yParams.nPoints).-welfareSpecList[k].betaC.*s.income.yTMat')
        end

        mul!(autVGrid,autSolnMat,autEVFlow)

        autVGridAll[:,k].=autVGrid
        openA0VGrid[:,k].=w1.VF.vGrid[s.a0Ind,:]

        if welfareSpecList[k].gammaC!=1.0
            changeGrid.=(autVGrid./w1.VF.vGrid[s.a0Ind,:]).^(1.0/(1.0-welfareSpecList[k].gammaC))
        else
            changeGrid.=exp.((1.0-welfareSpecList[k].betaC)*(autVGrid.-w1.VF.vGrid[s.a0Ind,:]))
        end

        gainA0Grid[:,k].=changeGrid

        meanIncGain[k]=changeGrid[meanIncInd]
        minIncGain[k]=changeGrid[minIncInd]
        maxIncGain[k]=changeGrid[maxIncInd]

        meanGain[k]=mean(changeGrid)
        minGain[k]=minimum(changeGrid)
        maxGain[k]=maximum(changeGrid)

        if welfareSpecList[k].gammaC!=1.0
            ergodicGainActual[k]=(dot(stIncDist,autVGrid)/(sum(stJointDist.repayDist.*w1.VF.vGrid)+dot(stJointDist.defaultDist,w1.VF.vDFutGrid)))^(1.0/(1.0-welfareSpecList[k].gammaC))
            ergodicGainUnexpected[k]=(dot(stIncDist,autVGridAll[:,k])/dot(stIncDist,openA0VGrid[:,k]))^(1.0/(1.0-welfareSpecList[k].gammaC))
        else
            ergodicGainActual[k]=exp((1.0-welfareSpecList[k].betaC)*(dot(stIncDist,autVGrid)-(sum(stJointDist.repayDist.*w1.VF.vGrid)+dot(stJointDist.defaultDist,w1.VF.vDFutGrid))))
            ergodicGainUnexpected[k]=exp((1.0-welfareSpecList[k].betaC)*(dot(stIncDist,autVGridAll[:,k])-dot(stIncDist,openA0VGrid[:,k])))
        end
        w0=deepcopy(w1)
        println([k,paramGrid[k,1],paramGrid[k,2],paramGrid[k,3],paramGrid[k,4],paramGrid[k,5],meanIncGain[k],minIncGain[k],maxIncGain[k],meanGain[k],minGain[k],maxGain[k],ergodicGainActual[k],ergodicGainUnexpected[k]])
    end


    summaryVals=hcat(paramGrid,meanIncGain,minIncGain,maxIncGain,meanGain,minGain,maxGain,ergodicGainActual,ergodicGainUnexpected)


    return welfareResults(summaryVals,autVGridAll,openA0VGrid,gainA0Grid)

end




#Applies the transition operator implied by government policy functions to a single cell of the ex ante distribution in good standing using approximate integration
function iterateDistributionCellApprox!(m::longTermBondSpec{F,S},s::longTermBondEval{F,S},i::S,j::S,RDistEA::Array{F,2},RDistEP::Array{F,2},DDistEP::Array{F,1}) where{F<:Real,S<:Integer}

    if s.pol.alwaysDefault[j,i]==true
        DDistEP[i]+=RDistEA[j,i]
    elseif s.pol.neverDefault[j,i]==true
        #If the government never defaults, then we begin the integration over the result under repayment at the first interval.
        #We set the location of the current upper bound in the m grid to the second point, the index of the relevant borrowing policy entry to the first,
        #and the most recent value of m to its lower bound
        mUBInd=2
        aPolicyInd=1
        lastMVal=s.income.mGrid[1]
        #Loop until we exit the upper boundary of the m space (the second check is to ensure that we never exit the range of the current borrowing policy function;
        #it should never matter, and if it does, the process should in general not lead to convergence).
        while (mUBInd<=(m.mParams.nPoints))&&(aPolicyInd<=(s.pol.mListLength[j,i]))
            #If the upper bound for the range in which the current entry of the borrowing policy function is valid is strictly less than the upper bound of the current
            #m interval, add the relevant contribution of that entry, update the last value of m reached, and increment only the variable marking our position in the borrowing
            #policy function
            if s.pol.mAPThreshold[i][aPolicyInd,j]<(s.income.mGrid[mUBInd])
                RDistEP[s.pol.apPolicy[i][aPolicyInd,j],i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*RDistEA[j,i]
                lastMVal=s.pol.mAPThreshold[i][aPolicyInd,j]
                aPolicyInd+=1
            else
                #Otherwise, add the relevant contribution of the current entry which lies in this m interval, update the last value of m reached to be its upper bound,
                #and then increment only the variable marking our position in the m grid
                RDistEP[s.pol.apPolicy[i][aPolicyInd,j],i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*RDistEA[j,i]
                lastMVal=s.income.mGrid[mUBInd]

                mUBInd+=1
            end
        end

    else
        #If the government does default only sometimes, find the position in the m grid of the upper bound of the interval in which the threshold level of m for default falls
        #and set the last value of m observed to the value directly preceding that upper bound
        mUBInd=searchsortedfirst(s.income.mGrid,s.pol.defThreshold[j,i])
        lastMVal=s.income.mGrid[mUBInd-1]
        #If there are any entire intervals in which default occurs, add their contribution to the output value
        if mUBInd>2
            for k in 3:mUBInd
                DDistEP[i]+=s.income.mProb[k-2]*RDistEA[j,i]

            end
        end
        #Set the index of the first relevant entry of the borowing policy function
        aPolicyInd=s.pol.firstRepayInd[j,i]

        #Add the contribution of default in the interval in which the threshold lies to the output value
        DDistEP[i]+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes*RDistEA[j,i]

        #Update the last value of m observed to the threshold level of m at which default occurs
        lastMVal=s.pol.defThreshold[j,i]
        #Loop until we exit the upper boundary of the m space (the second check is to ensure that we never exit the range of the current borrowing policy function;
        #it should never matter, and if it does, the process should in general not lead to convergence).
        while (mUBInd<=(m.mParams.nPoints))&&(aPolicyInd<=(s.pol.mListLength[j,i]))
            #If the upper bound for the range in which the current entry of the borrowing policy function is valid is strictly less than the upper bound of the current
            #m interval, add the relevant contribution of that entry, update the last value of m reached, and increment only the variable marking our position in the borrowing
            #policy function
            if s.pol.mAPThreshold[i][aPolicyInd,j]<(s.income.mGrid[mUBInd])
                RDistEP[s.pol.apPolicy[i][aPolicyInd,j],i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*RDistEA[j,i]
                lastMVal=s.pol.mAPThreshold[i][aPolicyInd,j]
                aPolicyInd+=1
            else
                #Otherwise, add the relevant contribution of the current entry which lies in this m interval, update the last value of m reached to be its upper bound,
                #and then increment only the variable marking our position in the m grid
                RDistEP[s.pol.apPolicy[i][aPolicyInd,j],i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*RDistEA[j,i]
                lastMVal=s.income.mGrid[mUBInd]
                mUBInd+=1
            end
        end
    end

    return RDistEP,DDistEP
end



#simulateDistribution generates the stationary, ergodic joint distribution of default state, income, and borrowing implied by an income process and government policy functions. Two versions are supplied. The first sets the initial distribution to have equal mass in every state with zero borrowing and not in default and all other states to have zero mass. The second accepts as an initial distribution an object of type debtDist. The arguments of the first version are:
#1. m: a model specification
#2. s: a collection of objects used in solving the model
#3. tol: the tolerance for convergence
#4. maxIter: the maximum number of iterations should convergence not be achieved
#Its output is an object of type debtDist

function simulateDistribution(m::longTermBondSpec{F,S},s::longTermBondEval{F,S},tol::F,maxIter::S,bothEAEP::Bool) where{F<:AbstractFloat,S<:Integer}

    #initialize variables containing the dimensions of the income grid and borrowing grid
    ydim=m.yParams.nPoints
    adim=m.aPoints

    #Initialize new and old copies of:
    #1. the noninflated joint distribution of income and borrowing at the beginning of the period, conditional on having access to financial markets, once uncertainty regarding regaining access to financial markets and income transitions in the current period has been resolved
    #2. the noninflated joint distribution of income and borrowing at the end of the period, conditional on having access to financial markets, before uncertainty regarding regaining access to financial markets and income transitions in the next period has been resolved
    #3. the noninflated distribution of income at the beginning of each period, conditional on being in default, once uncertainty regarding regaining access to financial markets and income transitions in the current period has been resolved
    #4. the noninflated distribution of income at the end of each period, conditional on being in default, before uncertainty regarding regaining access to financial markets and income transitions in the next period has been resolved

    #The naming convention is that myXDistYZ is the noninflated distribution for access state X (R or D, with R indicating access to financial markets and D indicating no access), Y timing (EA for beginning of period, EP for end of period), and Z update marker (New for newest iteration or Old for previous iteration)
    myRDistEAOld=zeros(adim,ydim)
    myRDistEANew=zeros(adim,ydim)

    myRDistEPOld=zeros(adim,ydim)
    myRDistEPNew=zeros(adim,ydim)


    myDDistEAOld=zeros(ydim)
    myDDistEANew=zeros(ydim)

    myDDistEPOld=zeros(ydim)
    myDDistEPNew=zeros(ydim)


    #Initialize an iteration counter and a measure of distance between successive iterations
    dREPDist=tol+one(F)
    dDEPDist=tol+one(F)
    dREADist=tol+one(F)
    dDEADist=tol+one(F)
    dDist=tol+one(F)
    iCount=1

    #Put full mass, equally divided, on states in which the government has access to financial markets and 0.0 incoming borrowing.
    for i in 1:ydim
        for j in 1:adim
            myRDistEPOld[j,i]=0.0
        end
        myDDistEPOld[i]=0.0
    end
    myRDistEAOld.=rand(adim,ydim)
    myDDistEAOld.=rand(ydim)
    inflFac=1/(sum(myRDistEAOld)+sum(myDDistEAOld))
    myRDistEAOld.*=inflFac
    myDDistEAOld.*=inflFac
    for i in 1:ydim

        for j in 1:adim

            iterateDistributionCellApprox!(m,s,i,j,myRDistEAOld,myRDistEPOld,myDDistEPOld)

        end
    end

    #Iterate until the maximum number of iterations has been reached or convergence (as defined by the tolerance) has been achieved
    while (iCount<=maxIter)&(dDist>=tol)

        #Set the new ex ante distributions according to the values implied by the income state transition matrix only

        mul!(myRDistEANew,myRDistEPOld,s.income.yTMat')
        mul!(myDDistEANew,s.income.yTMat,myDDistEPOld)

        #Iterate over income states, adding proportion theta of each mass in default in that state to the mass in that state with access to financial markets and 0 debt while subtracting it from the mass in default
        for i in 1:ydim
            myRDistEANew[s.a0Ind,i]+=m.theta*myDDistEANew[i]
            myDDistEANew[i]=(1.0-m.theta)*myDDistEANew[i]
        end

        #Iterate over income states
        for i in 1:ydim
            #Reset the current column of the ex post distribution for repayment to 0.0
            myRDistEPNew[:,i].=0.0
            #Set the ex post mass in default equal to the ex ante mass in default
            myDDistEPNew[i]=myDDistEANew[i]

            #Iterate over incoming borrowing states
            for j in 1:adim
                iterateDistributionCellApprox!(m,s,i,j,myRDistEANew,myRDistEPNew,myDDistEPNew)
            end
        end

        #Calculate the various sup norm measures of distance
        dREPDist=maximum(abs.(myRDistEPNew-myRDistEPOld))
        dDEPDist=maximum(abs.(myDDistEPNew-myDDistEPOld))
        dREADist=maximum(abs.(myRDistEANew-myRDistEAOld))
        dDEADist=maximum(abs.(myDDistEANew-myDDistEAOld))
        dDist=max(dREPDist,dDEPDist,dREADist,dDEADist)

        #Update the old distributions
        copyto!(myRDistEPOld,myRDistEPNew)
        copyto!(myDDistEPOld,myDDistEPNew)
        copyto!(myRDistEAOld,myRDistEANew)
        copyto!(myDDistEAOld,myDDistEANew)

        #At each iteration which is a multiple of 10% of the maximum number of iterations, print the current iteration number and variables tracking the convergence criteria
        if div(iCount,max(div(maxIter,10),1))==(iCount/max(div(maxIter,10),1))
            println([iCount,dREPDist,dDEPDist,dREADist,dDEADist])
        end
        #Increment the iteration counter
        iCount+=1
    end
    #println(sum(myDDistEPNew-myDDistEANew)/sum(myRDistEANew))
    #println(sum(myDDistEPNew-myDDistEANew)/sum(myRDistEPOld))

    println([iCount,dREPDist,dDEPDist,dREADist,dDEADist])
    #Create and return the joint distribution, remembering to uninflate both pieces by the scale factor used in the beginning
    if bothEAEP==false
        deflFac=(sum(myRDistEANew)+sum(myDDistEANew))^(-1)

        return debtDist(myRDistEANew*deflFac,myDDistEANew*deflFac)
    else
        deflFacEA=(sum(myRDistEANew)+sum(myDDistEANew))^(-1)
        deflFacEP=(sum(myRDistEPNew)+sum(myDDistEPNew))^(-1)

        return debtDist(myRDistEANew*deflFacEA,myDDistEANew*deflFacEA),debtDist(myRDistEPNew*deflFacEP,myDDistEPNew*deflFacEP)
    end
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
        copyto!(oldDist,newDist)
        iCount+=1
    end
    deflFac=sum(newDist)^(-1)
    return newDist*deflFac
end


#Writes a distribution object to a .csv file
function writeDistribution(d::debtDist{F},filePrefix::String,fileDir::String) where{F<:Real}
    writecsv(fileDir*filePrefix*"_repayDist.csv",d.repayDist)
    writecsv(fileDir*filePrefix*"_defaultDist.csv",d.defaultDist)
end

#Writes a welfare results object to a .csv file
function writeWelfareResults(w::welfareResults{F},filePrefix::String,fileSuffix::String,fileDir::String) where{F<:Real}
    writecsv(fileDir*filePrefix*"_summaryVals_"*fileSuffix*".csv",w.summaryVals)
    writecsv(fileDir*filePrefix*"_autVGridAll_"*fileSuffix*".csv",w.autVGridAll)
    writecsv(fileDir*filePrefix*"_openA0VGrid_"*fileSuffix*".csv",w.openA0VGrid)
    writecsv(fileDir*filePrefix*"_gainA0Grid_"*fileSuffix*".csv",w.gainA0Grid)
end


#Reads a welfare results object from a set of .csv files
function readWelfareResults(filePrefix::String,fileSuffix::String,fileDir::String)
    summaryVals=readcsv(fileDir*filePrefix*"_summaryVals_"*fileSuffix*".csv")
    autVGridAll=readcsv(fileDir*filePrefix*"_autVGridAll_"*fileSuffix*".csv")
    openA0VGrid=readcsv(fileDir*filePrefix*"_openA0VGrid_"*fileSuffix*".csv")
    gainA0Grid=readcsv(fileDir*filePrefix*"_gainA0Grid_"*fileSuffix*".csv")
    if eltype(summaryVals)!=Float64
        return welfareResults(Float64.(summaryVals[2:end,:]),autVGridAll,openA0VGrid,gainA0Grid)
    else
        return welfareResults(summaryVals,autVGridAll,openA0VGrid,gainA0Grid)
    end
end

function splitWelfareResults(w::welfareResults{F}) where{F<:Real}
    splitParamsUnique=unique(w.summaryVals[:,3:5],dims=1)
    nUnique=size(splitParamsUnique)[1]
    return [splitWelfareResults(w,splitParamsUnique[i,:]...) for i in 1:nUnique]
end

function splitWelfareResults(w::welfareResults{F},penMultC::F,penM::F,useVGPostDef::F) where{F<:Real}
    nTot=size(w.summaryVals)[1]
    specMarks=zeros(Bool,nTot)
    for i in 1:nTot
        if (w.summaryVals[i,3]==penMultC)&&(w.summaryVals[i,4]==penM)&&(w.summaryVals[i,5]==useVGPostDef)
            specMarks[i]=true
        end
    end
    summaryVals=w.summaryVals[specMarks,:]

    autVGridAll=w.autVGridAll[:,specMarks]
    openA0VGrid=w.openA0VGrid[:,specMarks]
    gainA0Grid=w.gainA0Grid[:,specMarks]

    return welfareResults(summaryVals,autVGridAll,openA0VGrid,gainA0Grid)
end




function decomposeWelfareBeta(m::longTermBondSpec{F,S},s::longTermBondEval{F,S},wPen::welfareResults{F},wNoPen::welfareResults{F},allPaths::meanPaths{F},tol::F,maxIter::S) where{F<:Real,S<:Integer}
    @assert m.gamma!=1
    @assert minimum(wPen.summaryVals[:,2])==m.gamma
    @assert maximum(wPen.summaryVals[:,2])==m.gamma
    bigT=length(allPaths.RMassPath)
    ydim=m.yParams.nPoints
    adim=m.aPoints

    wdim=size(wPen.summaryVals)[1]

    betaGrid=wPen.summaryVals[:,1]

    stIncDist=genStIncDist(m,s,tol,maxIter)


    vAlongMeanNoPen=zeros(1,wdim)
    vAlongMeanAut=zeros(1,wdim)
    yBar=dot(stIncDist,s.income.yGrid)
    for wInd in 1:wdim
        for littleT in 1:bigT
            vAlongMeanNoPen[1,wInd]+=betaGrid[wInd]^(littleT-1)*u(m,allPaths.consVal.cYNoPen.varPath[littleT])
            vAlongMeanAut[1,wInd]+=betaGrid[wInd]^(littleT-1)*u(m,yBar)
        end
        vAlongMeanNoPen[1,wInd]+=betaGrid[wInd]^(bigT)/(1-betaGrid[wInd])*u(m,allPaths.consVal.cYNoPen.varPath[bigT])
        vAlongMeanAut[1,wInd]+=betaGrid[wInd]^(bigT)/(1-betaGrid[wInd])*u(m,yBar)
    end


    onePlusLambdaD=((stIncDist'*wPen.openA0VGrid)./(stIncDist'*wNoPen.openA0VGrid)).^(1/(1.0-m.gamma))
    onePlusLambdaV=((stIncDist'*wNoPen.openA0VGrid)./vAlongMeanNoPen).^(1/(1.0-m.gamma)).*(vAlongMeanAut./(stIncDist'*wPen.autVGridAll)).^(1/(1.0-m.gamma))
    onePlusLambdaT=(vAlongMeanNoPen./vAlongMeanAut).^(1/(1.0-m.gamma))


    onePlusTrueLambda=((stIncDist'*wPen.openA0VGrid)./(stIncDist'*wPen.autVGridAll)).^(1/(1.0-m.gamma))


    return hcat(betaGrid,onePlusTrueLambda',onePlusLambdaD',onePlusLambdaV',onePlusLambdaT')

end



function decomposeLambdaD(m::longTermBondSpec{F,S},s::longTermBondEval{F,S},wPenUseC::welfareResults{F},wNoPenUseC::welfareResults{F},wPenUseG::welfareResults{F},wNoPenUseG::welfareResults{F},wGNoPen::debtWelfareEval{F},tol::F,maxIter::S) where{F<:Real,S<:Integer}
    @assert m.gamma!=1
    @assert minimum(wPenUseC.summaryVals[:,2])==m.gamma
    @assert maximum(wPenUseC.summaryVals[:,2])==m.gamma
    betaGrid=wPenUseC.summaryVals[:,1]

    stIncDist=genStIncDist(m,s,tol,maxIter)

    onePlusLambdaHD=((stIncDist'*wPenUseC.openA0VGrid)./(stIncDist'*wNoPenUseC.openA0VGrid)).^(1/(1.0-m.gamma))
    onePlusLambdaHHatD=((stIncDist'*wPenUseG.openA0VGrid)./(stIncDist'*wNoPenUseG.openA0VGrid)).^(1/(1.0-m.gamma))

    onePlusLambdaGD=ones(F,1,length(onePlusLambdaHD))*(dot(stIncDist,s.VF.vGrid[s.a0Ind,:])/dot(stIncDist,wGNoPen.VF.vGrid[s.a0Ind,:]))^(1/(1.0-m.gamma))



    return hcat(betaGrid,onePlusLambdaHD',onePlusLambdaHHatD',onePlusLambdaGD',(onePlusLambdaHD./onePlusLambdaHHatD)',(onePlusLambdaHHatD./onePlusLambdaGD)')

end
