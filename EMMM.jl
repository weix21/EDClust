using Distributions
using SpecialFunctions
using StatsBase

#EM method
function EMPolya(alpha0::AbstractArray,delta::AbstractArray,alpha::AbstractArray,Y::AbstractArray,TS::AbstractArray,p::AbstractArray,L::Int64,I::AbstractArray,J::Int64,K::Int64,EMNum::Int64,MMNum::Int64,stopc::Float64=1e-4)
    mu=Array{Array{Float64}}(undef,L)
    MTS=zeros(Int64,L)
    MY=zeros(Int64,J,L)
    dLike=0
    dLikeNew=0
    MLike=zeros(Float64,K)
    lp=log.(p)

    #Average the initial delta
    if L>1
        delta[:,:,2:L].=mean(delta[:,:,2:L],dims=3)
        for l in 2:L
            alpha[:,:,l]=alpha0+delta[:,:,l]
        end
    end

        @inbounds    for l in 1:L
            mu[l]=zeros(Float64,K,I[l])
            MTS[l]=maximum(TS[l])
            @inbounds    for i in 1:I[l]
                @inbounds for k in 1:K
                    MLike[k]=lp[k,l]+logpdf(DirichletMultinomial(TS[l][i],alpha[:,k,l]),view(Y[l],:,i))
                end
                m=maximum(MLike)
                @inbounds for k in 1:K
                    mu[l][k,i]=MLike[k]-m-log(sum(exp.(MLike[:].-m)))
                end
                dLike+=m+log(sum(exp.(MLike[:].-m)))

        end
        lp[:,l]=maximum(mu[l],dims=2).+log.(sum(exp.(mu[l].-maximum(mu[l],dims=2)),dims=2)).-log(I[l])
        @inbounds for j in 1:J
            MY[j,l]=maximum(Y[l][j,:])
        end
    end
    Newmu=copy(mu)
    Newlp=copy(lp)

    Flag=1
    @inbounds for t1 in 1:EMNum

        #For the first ten EM iterations, Recalculate and average the delta
        if 2<=t1<=10
            if L>1
                delta[:,:,2:L].=mean(delta[:,:,2:L],dims=3)
                for l in 2:L
                    alpha[:,:,l]=alpha0+delta[:,:,l]
                end
            end

            @inbounds    for l in 1:L
                @inbounds    for i in 1:I[l]
                    @inbounds for k in 1:K
                        MLike[k]=lp[k,l]+logpdf(DirichletMultinomial(TS[l][i],alpha[:,k,l]),view(Y[l],:,i))
                    end
                    m=maximum(MLike)
                    @inbounds for k in 1:K
                        Newmu[l][k,i]=MLike[k]-m-log(sum(exp.(MLike[:].-m)))
                    end
                    dLike+=m+log(sum(exp.(MLike[:].-m)))

            end

            end

        end

        #MM step
        alpha0, delta, alpha=MMPolya(alpha0,delta,alpha,Y,TS,MY,MTS,Newmu,L,I,J,K,MMNum)

        #Compute the updated likelihood
        dLikeNew=0
        @inbounds for l in 1:L
            @inbounds for i in 1:I[l]
                s=0
                @inbounds for k in 1:K
                    MLike[k]=lp[k,l]+logpdf(DirichletMultinomial(TS[l][i],alpha[:,k,l]),view(Y[l],:,i))
                end
                m=maximum(MLike)
                @inbounds for k in 1:K
                    Newmu[l][k,i]=MLike[k]-m-log(sum(exp.(MLike[:].-m)))
                end
                dLikeNew+=m+log(sum(exp.(MLike[:].-m)))
            end
            Newlp[:,l]=maximum(mu[l],dims=2).+log.(sum(exp.(mu[l].-maximum(mu[l],dims=2)),dims=2)).-log(I[l])
        end

        #Stop criteria
        if dLikeNew>dLike
            if abs((dLikeNew-dLike)/dLike)> stopc
                dLike=dLikeNew
                mu=copy(Newmu)
                lp=copy(Newlp)
            else
                Flag=0
                println(dLike)
                println(dLikeNew)
            end
        else
            Flag=0
            println("decreasing likelihood!")
            println(dLike)
            println(dLikeNew)
        end

        if(Flag==0)
            println("Niter=",t1)
            break
        end

    end

    #Determine cell label based on mu
    groups=mapslices(argmax, mu[1], dims=1)[:]
    if L>1
        @inbounds for l in 2:L
            groups=vcat(groups,mapslices(argmax, mu[l], dims=1)[:])
        end
    end

    #println("Niter=",t1)
    #for l in 1:L
        #println(freqtable(W[l],groups[l][1,:]))
    #end

    return groups, dLikeNew, alpha0, delta, alpha, lp
end






#MM method
function MMPolya(alpha0::AbstractArray,delta::AbstractArray,alpha::AbstractArray,Y::AbstractArray,TS::AbstractArray,MY::AbstractArray,MTS::AbstractArray,mu::AbstractArray,L::Int64,I::AbstractArray,J::Int64,K::Int64,MMNum::Int64)
    Newdelta=similar(delta)
    Newalpha0=similar(alpha0)
    for t2 in 1:MMNum
        DS=zeros(J,K,L)
        DS2=zeros(J,K,L)
        DS1=zeros(K,L)
    @inbounds  for l in 1:L
                    @inbounds for k in 1:K
                        salpha=sum(view(alpha,:,k,l))

                        n2=maximum(view(mu[l],k,:)[TS[l].>=1])
                        @inbounds for c1 in 0:MTS[l]-1
                            DS1[k,l]+=sum(exp.(view(mu[l],k,:)[TS[l].>=c1+1].-n2))/(salpha+c1)
                        end
                        DS1[k,l]=n2+log(DS1[k,l])

                        @inbounds   for j in 1:J
                            if length(view(mu[l],k,:)[Y[l][j,:].>=1])==0
                                DS2[j,k,l]=-10000
                                DS[j,k,l]=0
                            else
                                n3=maximum(view(mu[l],k,:)[Y[l][j,:].>=1])
                                @inbounds   for c2 in 0:MY[j,l]-1
                                    DS2[j,k,l]+=sum(exp.(view(mu[l],k,:)[Y[l][j,:].>=c2+1].-n3))/(view(alpha,j,k,l).+c2)
                                end
                                DS2[j,k,l]=n3+log(DS2[j,k,l])
                                DS[j,k,l]=exp(DS2[j,k,l]-DS1[k,l])
                             end

                        end
                    end
                end
        Newdelta=delta.*DS
        Newalpha0=alpha0.*(sum(exp.(DS2),dims=3)./sum(exp.(DS1),dims=2)')
        delta=copy(Newdelta)
        alpha0=copy(Newalpha0)
        alpha0[alpha0.<1e-100].=1e-100
        for l in 1:L
            alpha[:,:,l]=alpha0+delta[:,:,l]
        end
    end
    alpha0, delta, alpha
end


#EMMM method
function fitPolya(Count::AbstractArray,Sample::AbstractArray,InitVal::AbstractArray,EMNum::Int64=100,MMNum::Int64=5,stopc::Float64=1e-4)
    L=length(unique(Sample))
    I=counts(Sample)
    J=size(InitVal)[1]
    K=size(InitVal)[2]

    Y=Array{Array{Int64}}(undef,L)
    for l in 1:L
        Y[l]=Count[:,(Sample.==l)[:]]
    end

    TS=Array{Array{Int64}}(undef,L)
    for l in 1:L
        TS[l]=sum(Y[l],dims=1)[:]
    end

    Salpha0=InitVal[:,:]
    Salpha0[Salpha0.<1e-10].=1e-10
    Sdelta=zeros(Float64,J,K,L)
    Salpha=zeros(Float64,J,K,L)

    for l in 1:L
        if l!=1
            Sdelta[:,:,l].=1e-5
        end
        Salpha[:,:,l]=Salpha0+Sdelta[:,:,l]
    end

    Sp=1/K.*ones(K,L)
    return EMPolya(Salpha0,Sdelta,Salpha,Y,TS,Sp,L,I,J,K,EMNum,MMNum,stopc)
end
