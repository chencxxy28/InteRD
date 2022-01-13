#'@title Variable selection in generalized linear models (GLM)
#'@description The function \code{\link{ELCIC.glm.single}} provides values of several model selection criteria including AIC, BIC, GIC and ELCIC, given a candidate mean model.
#'@usage ELCIC.glm.single(x, y, index.var=NULL, name.var = NULL, dist)
#'@param x A matrix containing covariates. The first column should contain all ones corresponding to the intercept if an intercept is considered in a candidate mean model.
#'@param y A vector containing outcomes.
#'@param index.var A vector containing index corresponding to covariates (including the intercept) in the candidate mean model. See more in Details section.
#'@param name.var A vector containing names of covariates in the candidate mean model. The names should be a subset of the column names of the x matrix. See more in Details section.
#'@param dist A specified distribution. It can be "gaussian", "poisson", and "binomial".
#'
#'@details "x" and "y" should be all observed. The corresponding individual data will be omitted in analysis if any missingness is detected.
#'
#'Either arguments "index.var" or "name.var" is used to identify the candidate mean model. If both arguments are provided, only the argument "name.var" will be used.
#'
#'@return A vector containing values of information criteria including AIC, BIC, GIC and ELCIC.
#'
#'@examples
#'## tests
#'# load data
#'data(glmsimdata)
#'x<-glmsimdata$x
#'y<-glmsimdata$y
#'#candidate model index
#'name.var<-c("intercept","x1","x2")
#'index.var<-c(1,2,3)
#'criteria<-ELCIC.glm.single(x, y, index.var =index.var, name.var = NULL, dist="poisson")
#'criteria
#'
#'@export
#'@import MASS
#'@importFrom stats as.formula glm rbinom rnbinom rnorm rpois

ELCIC.glm.single<-function (x,y,index.var=NULL,name.var=NULL,dist)
{
    if(!is.matrix(x))
    {stop("x should be in a matrix format")}else if(!is.vector(y)){stop("y should be in a vector format")}else
        if (!dist %in% c("gaussian","binomial","poisson")){stop("Invalid type of dist. It should be one of gaussian,binomial,poisson")}

    #delete individual data with missingness
    data.comb<-cbind(x,y)
    na.index<-which(is.na(data.comb),arr.ind=TRUE)
    if(nrow(na.index)>0)
    {
    y<-y[-unique(na.index[,1])]
    x<-x[-unique(na.index[,1]),]
    }

    if(!is.null(name.var))
    {
        if(length(name.var)>ncol(x)|length(unique(name.var))<length(name.var)){stop("Invalid candidate model provided")}
        index.var<-match(name.var,colnames(x))
    }
        if(length(index.var)>ncol(x)|max(index.var)>ncol(x)|length(unique(index.var))<length(index.var)){stop("Invalid candidate model provided")}
        samplesize<-nrow(x)
        beta<-rep(0,ncol(x))
        x_candidate<-x[,index.var]
        reg<-glm(y~x_candidate-1, family = dist)
        beta[index.var]<-reg$coefficients
        p<-length(reg$coefficients)


    lambda<-lambda.find.glm(x=x,y=y,beta=beta, dist=dist)
    Z<-ee.glm(x,y,beta,dist)

    likelihood<-apply(Z,2,function(x) {
    2*log(1+t(lambda)%*%x)})%*%rep(1,samplesize)
    ELCIC<-likelihood+p*log(samplesize)

    #get aic, bic, and gic
    AIC<-reg$aic
    tr_matrix<-II(x=x_candidate,y=y,beta=reg$coefficients,size=length(index.var),samplesize=samplesize, dist=dist)%*%ginv(JJ(x=x_candidate,y=y,beta=reg$coefficients,size=length(index.var),samplesize=samplesize,dist=dist))
    GIC<-AIC-2*p+2*sum(diag(tr_matrix))
    BIC<-AIC-2*p+p*log(samplesize)
    return(c(ELCIC=ELCIC, AIC=AIC, BIC=BIC, GIC=GIC))
}





#'@title Calculate ELCIC for a given candidate model under the GEE framework with longitudinal data that is complete or missing completely at random
#'@description The function \code{\link{ELCIC.gee.single}} calculates ELCIC for a given marginal candidate mean model with a specified "working" correlation structure. It is able to simultaneously evaluate the marginal mean and "working" correlation structures.
#'@usage ELCIC.gee.single(x, y, r, id, time, index.var=NULL, name.var = NULL,
#'                      dist, corstr, joints=TRUE)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept. If y and x are missing completely at random, use NA to indicate missingness and specify them in the argument "r".
#'@param y A vector containing outcomes. If y and x are missing completely at random, use NA to indicate missingness and specify them in the argument "r".
#'@param r A vector indicating the observation of data: 1 for an observed record (both outcome and covariates are observed for a given subject), and 0 for an unobserved record. The default setup is that all elements of r are ones, indicating all data are observed.
#'@param id A vector indicating subject id.
#'@param time The number of observations in total for each subject
#'@param index.var A vector containing index corresponding to covariates in the candidate model. See more in Details section.
#'@param name.var A vector containing names of covariates in the candidate model. The names should be a subset of the column names of the x matrix. See more in Details section.
#'@param dist A specified distribution. It can be "gaussian", "poisson", and "binomial".
#'@param corstr A condidate correlation structure. It can be either of "independence", "exchangeable", and "ar1".
#'@param joints A logic value for joint selection of marginal mean and "working" correlation structures. The default is TRUE. See more in Details section.
#'
#'@return A ELCIC value for a given candidate model.
#'
#'@details Either of the arguments, "index.var" and "name.var", is used to identify the candidate mean model. If both arguments are provided, only the argument "name.var" will be used.
#'
#'When the argument "joints" is TRUE, \code{\link{ELCIC.gee.single}} will calculate ELCIC value based on the function \code{\link{lambda.find.gee}} and \code{\link{ee.gee}}, which involve estimating equations for both marginal mean and correlation coefficients. When the argument "joints" is FALSE, \code{\link{ELCIC.gee.single}} will calculate ELCIC value based on the function \code{\link{lambda.find.gee.mean}} and \code{\link{ee.gee.mean}}, which only involve estimating equations for marginal mean.
#'
#'@examples
#'## tests
#'# load data
#'data(geesimdata)
#'x<-geesimdata$x
#'y<-geesimdata$y
#'id<-geesimdata$id
#'r<-rep(1,nrow(x))
#'time<-3
#'corstr<-"exchangeable"
#'dist<-"poisson"
#'criteria<-ELCIC.gee.single(x=x,y=y,r=r,id=id,time=time,index.var=c(1,2,3),
#'            name.var=NULL,dist=dist,corstr=corstr)
#'criteria
#'
#'@export
#'@importFrom geepack geeglm

ELCIC.gee.single<-function(x,y,r,id,time,index.var=NULL,name.var=NULL,dist,corstr,joints=TRUE)
{
    if(!is.matrix(x))
    {stop("x should be in a matrix format")}else if(!is.vector(y)){stop("y should be in a vector format")}else
        if (!dist %in% c("gaussian","binomial","poisson")){stop("Invalid type of dist. It should be one of gaussian,binomial,poisson")}
    if (!corstr %in% c("ar1","exchangeable","independence")){stop("Invalid type of correlation structure for outcomes. It should be one of ar1,exchangeable,independence")}


    if(!is.null(name.var))
    {
        if(length(name.var)>ncol(x)|length(unique(name.var))<length(name.var)){stop("Invalid candidate model provided")}
        index.var<-match(name.var,colnames(x))
    }
    if(length(index.var)>ncol(x)|max(index.var)>ncol(x)|length(unique(index.var))<length(index.var)){stop("Invalid candidate model provided")}

    if(joints)
    {
    samplesize<-length(unique(id))
    beta<-rep(0,ncol(x))

    #delete individual data with missingness, in order to run GEE
    data.comb<-cbind(x,y,id)
    na.index<-which(is.na(data.comb),arr.ind=TRUE)

    y.gee<-y
    id.gee<-id
    x.gee<-x

    if(nrow(na.index)>0)
    {
        y.gee<-y[-unique(na.index[,1])]
        id.gee<-id[-unique(na.index[,1])]
        x.gee<-x[-unique(na.index[,1]),]
        y[unique(na.index[,1])]<-0  #only use observed values
        x[unique(na.index[,1]),]<-0  #only use observed values
    }

    x_candidate<-(x.gee[,index.var])
    colnames(x_candidate)<-seq_len(ncol(x_candidate))
    y.gee<-as.matrix(y.gee,col=1)
    data<-data.frame(y.gee,x_candidate)

        fit<-geeglm(y.gee~x_candidate-1,data=data,family =dist,id=id.gee,corstr = corstr)
        beta[index.var]<-fit$coefficients
        pbeta<-length(fit$coefficients)

        if(corstr=="independence")
        {
            rho<-0
            p<-pbeta
        }else{
            rho<-unlist(summary(fit)$corr[1]) #correlation coefficients
            p<-pbeta+1
        }
        phihat<-unlist(summary(fit)$dispersion[1])
        Z<-as.matrix(ee.gee(y=y,x=x,r=r,id=id,beta=beta,rho=rho,phi=phihat,dist=dist,corstr=corstr))
        # epi<-1/samplesize
        # model<-function(lambda)
        # {
        #     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
        #     {x/(1+t(lambda)%*%x)}
        #     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,samplesize)
        # }
        lambda<-lambda.find.gee(x=x,y=y,id=id,beta=beta,r=r,dist=dist,rho=rho,phi=phihat,corstr=corstr)
        # lambda<-multiroot(f = model, start = rep(0,length(beta)+time-1))$root
        likelihood<-apply(Z,2,function(x) {
            2*log(1+t(lambda)%*%x)})%*%rep(1,samplesize)
        ELCIC<-likelihood+p*log(samplesize)
    return(ELCIC=ELCIC)
    }else{
        samplesize<-length(unique(id))
        beta<-rep(0,ncol(x))

        #delete individual data with missingness, in order to run GEE
        data.comb<-cbind(x,y,id)
        na.index<-which(is.na(data.comb),arr.ind=TRUE)

        y.gee<-y
        id.gee<-id
        x.gee<-x

        if(nrow(na.index)>0)
        {
            y.gee<-y[-unique(na.index[,1])]
            id.gee<-id[-unique(na.index[,1])]
            x.gee<-x[-unique(na.index[,1]),]
            y[unique(na.index[,1])]<-0
            x[unique(na.index[,1]),]<-0
        }

        x_candidate<-(x.gee[,index.var])
        colnames(x_candidate)<-seq_len(ncol(x_candidate))
        y.gee<-as.matrix(y.gee,col=1)
        data<-data.frame(y.gee,x_candidate)

        fit<-geeglm(y.gee~x_candidate-1,data=data,family =dist,id=id.gee,corstr = corstr)
        beta[index.var]<-fit$coefficients
        p<-length(fit$coefficients)

        if(corstr=="independence")
        {
            rho<-0
            #p<-pbeta
        }else{
            rho<-unlist(summary(fit)$corr[1]) #correlation coefficients
            #p<-pbeta+1
        }
        phihat<-unlist(summary(fit)$dispersion[1])
        Z<-as.matrix(ee.gee.mean(y=y,x=x,r=r,id=id,beta=beta,rho=rho,phi=phihat,dist=dist,corstr=corstr))
        # epi<-1/samplesize
        # model<-function(lambda)
        # {
        #     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
        #     {x/(1+t(lambda)%*%x)}
        #     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,samplesize)
        # }
        lambda<-lambda.find.gee.mean(x=x,y=y,id=id,beta=beta,r=r,dist=dist,rho=rho,phi=phihat,corstr=corstr)
        # lambda<-multiroot(f = model, start = rep(0,length(beta)+time-1))$root
        likelihood<-apply(Z,2,function(x) {
            2*log(1+t(lambda)%*%x)})%*%rep(1,samplesize)
        ELCIC<-likelihood+p*log(samplesize)
        return(ELCIC=ELCIC)
    }
}




#'@title Calculate ELCIC for a given candidate model under the WGEE framework for longitudinal data with drop-out and the assumption of missing at random
#'@description The function \code{\link{ELCIC.wgee.single}} to calculate ELCIC for a given candidate mean model with a specified "working" correlation structure. It is able to simultaneously evaluate the marginal mean and "working" correlation structures. The longutindal data has dropout missingness with the assumption of missing at random.
#'@usage ELCIC.wgee.single(x,y,x_mis,r,id,time,index.var=NULL,
#'       name.var=NULL,dist,corstr,joints=TRUE,lag=1)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept if the intercept is considered in the candidate model. The covariate matrix should be complete. See more in Details section.
#'@param y A vector containing outcomes. Use NA to indicate missing outcomes.
#'@param x_mis A matrix containing covariates for the missing data model. The first column should be all ones corresponding to the intercept. This covariate matrix should be complete and all observed. See more in Details section.
#'@param r A vector indicating the observation of outcomes: 1 for an observed record, and 0 for an unobserved record.
#'@param id A vector indicating subject id.
#'@param time The number of observations in total for each subject
#'@param index.var A vector containing index corresponding to the covariates in the candidate model. See more in Details section.
#'@param name.var A vector containing names of covariates in the candidate model. The names should be a subset of the column names of the x matrix. See more in Details section.
#'@param dist A specified distribution. It can be "gaussian", "poisson", and "binomial".
#'@param corstr A condidate correlation structure. It can be "independence","exchangeable", and "ar1".
#'@param joints A logic value for joint selection of marginal mean and "working" correlation structures. The default is TRUE. See more in Details section.
#'@param lag A numeric value indicating lag-response involved in the missing data model. It can be one of 0, 1, and 2. The default is 1.

#'@return A matrix containing values of calculated estimating equations.
#'
#'@details The covariate matrix "x" should be complete. If missing data are present in "x", the elements in covariate vector will be replaced by zeros for individuals who have missing covariates.
#'
#'The argument "x_mis" includes all covariates to fit the missing data model. It does not contains a lag variable based on the outcome y. The argument "lag" in this function will automatically add lag-response variables to indicate the data with missing at random????.
#'
#'Either of the arguments, "index.var" and "name.var", is used to identify the candidate mean model. If both arguments are provided, only the argument "name.var" will be used.
#'
#'When the argument "joints" is TRUE, \code{\link{ELCIC.wgee.single}} will calculate ELCIC based on the function \code{\link{lambda.find.wgee}} and \code{\link{ee.wgee}}, which involve estimating equations for both marginal mean and correlation coefficients. When the argument "joints" is FALSE, \code{\link{ELCIC.wgee.single}} will calculate ELCIC value based on the function \code{\link{lambda.find.wgee.mean}} and \code{\link{ee.wgee.mean}}, which only involve estimating equations for marginal mean.
#'
#'@examples
#'## tests
#'# load data
#'data(wgeesimdata)
#'corstr<-"exchangeable"
#'dist<-"binomial"
#'x<-wgeesimdata$x
#'y<-wgeesimdata$y
#'x_mis<-wgeesimdata$x_mis
#'r<-wgeesimdata$obs_ind
#'id<-wgeesimdata$id
#'time<-3
#'index.var<-c(1,2,3)
#'ELCIC_value<-ELCIC.wgee.single(x,y,x_mis,r,id,time,index.var,name.var=NULL,
#'                      dist,corstr,joints=TRUE)
#'ELCIC_value

#'@export
#'@importFrom wgeesel wgee ylag
#'
ELCIC.wgee.single<-function(x,y,x_mis,r,id,time,index.var=NULL,name.var=NULL,dist,corstr,joints=TRUE,lag=1)
{
    if(!is.matrix(x))
    {stop("x should be in a matrix format")}else if(!is.vector(y)){stop("y should be in a vector format")}else if(!is.matrix(x_mis)){stop("x_mis should be in a matrix format")}else
        if (!dist %in% c("gaussian","binomial","poisson")){stop("Invalid type of dist. It should be one of gaussian,binomial,poisson")}
    if (!corstr %in% c("ar1","exchangeable","independence")){stop("Invalid type of correlation structure for outcomes. It should be one of ar1,exchangeable,independence")}

    #replace missingness with 0 in x matrix
    na.index<-which(is.na(x),arr.ind=TRUE)
    r.x<-rep(1,nrow(x))
    r.x[unique(na.index[,1])]<-0
    if(nrow(na.index)>0)
    {
        if (sum(r.x==r)!=nrow(x))
        {
            stop("missingness pattern in x should be the same to the missing pattern in y.")
        }
        x[unique(na.index[,1]),]<-0
        warning("Covariate matrix x should be fully observed. The elements in covariate vector are replaced by zeros for individuals who have missing covariates.")
    }

    if(!is.null(name.var))
    {
        if(length(name.var)>ncol(x)|length(unique(name.var))<length(name.var)){stop("Invalid candidate model provided")}
        index.var<-match(name.var,colnames(x))
    }
    if(length(index.var)>ncol(x)|max(index.var)>ncol(x)|length(unique(index.var))<length(index.var)){stop("Invalid candidate model provided")}

    if(lag>2 | lag>=time)
    {
    stop("Invalid type of lag. It should be less than 3 and time")
    }else if(lag==0)
    {
        warning("No lag of response added may indicate missing completely at random. GEE may be used.")
    }

    #generate ylag1
    if(lag!=0)
    {
    samplesize<-length(unique(id))
    ylag.cov<-ylag(y=y,id=id,lag=1)
    ylag.cov[is.na(ylag.cov)]<-0
    ylag.cov[rep(seq_len(time),time=samplesize)==1]<-NA
    x_mis<-cbind(x_mis,ylag1.cov=ylag.cov)
    }
    if(lag==2){
        #generate ylag2
        ylag.cov<-ylag(y=y,id=id,lag=lag)
        x_mis<-cbind(x_mis,ylag2.cov=ylag.cov)
    }

    if(joints)
    {
        samplesize<-length(unique(id))
        beta<-rep(0,ncol(x))
        x_candidate<-(x[,index.var])
        colnames(x_candidate)<-seq_len(ncol(x_candidate))
        y<-as.matrix(y,col=1)
        #data_wgee<-data.frame(y,x_candidate,r,x_mis)
        data_wgee<-data.frame(r,x_mis)

        mismodel_formula<-paste("r~",colnames(data_wgee)[3])
        if(ncol(data_wgee)>3)
        {
        for(jj in 4:ncol(data_wgee))
        {
        mismodel_formula<-paste(mismodel_formula,"+",colnames(data_wgee)[jj])
        }
        }
        mismodel_formula<-as.formula(mismodel_formula)
        fit<-wgee(y~x_candidate-1,data_wgee,id,family=dist,corstr =corstr,scale = NULL,mismodel =mismodel_formula)

        beta[index.var]<-as.vector(summary(fit)$beta)
        pbeta<-length(summary(fit)$beta)

        if(corstr=="independence")
        {
            rho<-0
            p<-pbeta
        }else{
            rho<-summary(fit)$corr
            p<-pbeta+1
        }
        phi<-summary(fit)$phi
        gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1])
        pi<-cond.prob(x_mis,gamma,id,time)
        Z<-as.matrix(ee.wgee(y,x,r, pi,id,time,beta,rho,phi,dist,corstr))
        # epi<-1/samplesize
        # model<-function(lambda)
        # {
        #     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
        #     {x/(1+t(lambda)%*%x)}
        #     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,samplesize)
        # }
        lambda<-lambda.find.wgee(y,x,r, pi,id,time,beta,rho,phi,dist,corstr)
        # lambda<-multiroot(f = model, start = rep(0,length(beta)+time-1))$root
        likelihood<-apply(Z,2,function(x) {
            2*log(1+t(lambda)%*%x)})%*%rep(1,samplesize)
        ELCIC<-likelihood+p*log(samplesize)

    }else{
        samplesize<-length(unique(id))
        beta<-rep(0,ncol(x))
        x_candidate<-(x[,index.var])
        colnames(x_candidate)<-seq_len(ncol(x_candidate))
        y<-as.matrix(y,col=1)
        #data_wgee<-data.frame(y,x_candidate,r,x_mis)
        data_wgee<-data.frame(r,x_mis)

        mismodel_formula<-paste("r~",colnames(data_wgee)[3])
        if(ncol(data_wgee)>3)
        {
            for(jj in 4:ncol(data_wgee))
            {
                mismodel_formula<-paste(mismodel_formula,"+",colnames(data_wgee)[jj])
            }
        }
        mismodel_formula<-as.formula(mismodel_formula)
        fit<-wgee(y~x_candidate-1,data_wgee,id,family=dist,corstr =corstr,scale = NULL,mismodel =mismodel_formula)

        beta[index.var]<-as.vector(summary(fit)$beta)
        p<-length(summary(fit)$beta)

        if(corstr=="independence")
        {
            rho<-0
            #p<-pbeta
        }else{
            rho<-summary(fit)$corr
            #p<-pbeta+1
        }
        phi<-summary(fit)$phi
        gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1])
        pi<-cond.prob(x_mis,gamma,id,time)
        Z<-as.matrix(ee.wgee.mean(y,x,r, pi,id,time,beta,rho,phi,dist,corstr))
        # epi<-1/samplesize
        # model<-function(lambda)
        # {
        #     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
        #     {x/(1+t(lambda)%*%x)}
        #     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,samplesize)
        # }
        lambda<-lambda.find.wgee.mean(y,x,r, pi,id,time,beta,rho,phi,dist,corstr)
        # lambda<-multiroot(f = model, start = rep(0,length(beta)+time-1))$root
        likelihood<-apply(Z,2,function(x) {
            2*log(1+t(lambda)%*%x)})%*%rep(1,samplesize)
        ELCIC<-likelihood+p*log(samplesize)

    }
    return(ELCIC=ELCIC)
}






#'@title The whole variable selection procedure in GLM
#'@description The function \code{\link{ELCIC.glm}} provides the overall procedure for variable selection in GLM.
#'@usage ELCIC.glm(x,y,candidate.sets,name.var.sets=NULL,dist)
#'@param x A matrix containing covariates. The first column should contain all ones corresponding to the intercept if the intercept is expected in the mean structure.
#'@param y A vector containing outcomes.
#'@param candidate.sets A list containing index corresponding to covariates in each candidate model. See more in details section.
#'@param name.var.sets A list containing names of covariates coresponding to each candidate model. The names should be a subset of the column names of the x matrix. See more in details section.
#'@param dist A specified distribution. It can be "gaussian", "poisson", and "binomial".
#'
#'@return A matrix with each element containing ELCIC value for each candidate model (in columns) and (in rows)????
#'
#'@details "x" and "y" should be all observed. The corresponding individual data will be deleted if any missingness is detected.
#'
#'Either of the arguments, "candidate.sets" and "name.var.sets", is used to identify the set of candidate mean models. If both arguments are provided, only the argument "name.var.sets" will be used.
#'
#'
#'@examples
#'## tests
#'# load data
#'data(glmsimdata)
#'x<-glmsimdata$x
#'y<-glmsimdata$y
#'#candidate model index
#'candidate.sets<-list(c(1,2),c(1,2,3),c(1,2,3,4))
#'criteria<-ELCIC.glm(x, y, candidate.sets, name.var.sets = NULL, dist="poisson")
#'criteria
#'
#'@export

ELCIC.glm<-function(x,y,candidate.sets,name.var.sets=NULL,dist)
{
    if(!is.null(name.var.sets))
    {
    criterion.all<-rep()
    for (i in seq_len(length(name.var.sets)))
    {
        criterion<-ELCIC.glm.single(x=x,y=y,index.var=NULL,name.var=name.var.sets[[i]],dist=dist)
        criterion.all<-cbind(criterion.all,criterion)
    }
    colnames(criterion.all)<-name.var.sets
    criterion.all
    }else{
        criterion.all<-rep()
        for (i in seq_len(length(candidate.sets)))
        {
            criterion<-ELCIC.glm.single(x=x,y=y,index.var=candidate.sets[[i]],name.var=NULL,dist=dist)
            criterion.all<-cbind(criterion.all,criterion)
        }
        colnames(criterion.all)<-candidate.sets
        criterion.all
    }
}





#'@title The whole procedure for joint selection of marginal mean and correlation structures in longitudinal data that is complete or missing completely at random
#'@description The function \code{\link{ELCIC.gee}} provides the overall procedure for joint selection of marginal mean and correlation structures in longitudinal data that is complete or missing completely at random.
#'@usage ELCIC.gee(x,y,r,id,time,candidate.sets=NULL,name.var.sets=NULL,dist,
#'       candidate.cor.sets=c("independence","exchangeable", "ar1"), joints=TRUE)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept if the intercept is expected in the marginal mean. The covariate matrix should be complete. NA values will be replaced by 0 if missingness is detected in x.
#'@param y A vector containing outcomes. If y is missing completely at random, use NA to indicate missing outcomes and specify them in the argument "r".
#'@param r A vector indicating the observation of data: 1 for an observed record (both outcome and covariates are observed for a given subject), and 0 for an unobserved record. The default setup is that all elements of r are ones, indicating all data are observed.
#'@param id A vector indicating subject id.
#'@param time The number of observations in total for each subject
#'@param candidate.sets A list containing index corresponding to covariates in the candidate model. See more in Details section.
#'@param name.var.sets A list containing names of ccovariates in the candidate model. The names should be a subset of the column names of the x matrix. See more in Details section.
#'@param dist A specified distribution. It can be "gaussian", "poisson", and "binomial".
#'@param candidate.cor.sets A vector containing candidate correlation structures. When joints=TRUE, it can be any subset of c("independence", "exchangeable", "ar1"). The default is c("independence", "exchangeable", "ar1"). When joints=FALSE, it should be either of "independence", "exchangeable", and "ar1". See more in Details section.
#'@param joints A logic value for joint selection of marginal mean and "working" correlation structures. The default is TRUE. See more in Details section.
#'
#'@return A matrix with each element containing ELCIC for each candidate model.
#'
#'@details Either arguments "candidate.sets" or "name.var.sets" is used to identify the set of candidate mean models. If both arguments are provided, only the argument "name.var.sets" will be used.
#'
#'When joints=TRUE, the argument "candidate.cor.sets" can contain multiple correlation structures; however, when joints=FALSE, it should contain either of "independence","exchangeable", "ar1". If multiple correlation structures are provided, only the first one will be used.
#'
#'@examples
#'## tests
#'# load data
#'data(geesimdata)
#'x<-geesimdata$x
#'y<-geesimdata$y
#'id<-geesimdata$id
#'r<-rep(1,nrow(x))
#'time<-3
#'candidate.sets<-list(c(1,2),c(1,2,3))
#'candidate.cor.sets<-c("exchangeable")
#'dist<-"poisson"
#'criterion.elcic<-ELCIC.gee(x=x,y=y,r=r,id=id,time=time,candidate.sets=candidate.sets,
#'                 name.var.sets=NULL,dist=dist,candidate.cor.sets=candidate.cor.sets)
#'criterion.elcic
#'
#'@export

ELCIC.gee<-function(x,y,r,id,time,candidate.sets=NULL,name.var.sets=NULL,dist,candidate.cor.sets=c("independence",
                                                                                                   "exchangeable", "ar1"), joints=TRUE)
{
    if(joints)
    {
        if(!is.null(name.var.sets))
        {
            criterion.elcic<-rep()
            for(j in seq_len(length(candidate.cor.sets)))
            {
                corstr<-candidate.cor.sets[j]
                criterion.all<-rep()
                for (i in seq_len(length(name.var.sets)))
                {
                    criterion<-ELCIC.gee.single(x=x,y=y,r=r,id=id,time=time,index.var=NULL,name.var=name.var.sets[[i]],dist=dist,corstr=corstr)
                    criterion.all<-c(criterion.all,criterion)
                }
                criterion.elcic<-rbind(criterion.elcic,criterion.all)
                #print(j)
            }
            colnames(criterion.elcic)<-name.var.sets
            rownames(criterion.elcic)<-candidate.cor.sets
            criterion.elcic
        }else{
    criterion.elcic<-rep()
    for(j in seq_len(length(candidate.cor.sets)))
    {
        corstr<-candidate.cor.sets[j]
        criterion.all<-rep()
        for (i in seq_len(length(candidate.sets)))
        {
            criterion<-ELCIC.gee.single(x=x,y=y,r=r,id=id,time=time,index.var=candidate.sets[[i]],name.var=NULL,dist=dist,corstr=corstr)
            criterion.all<-c(criterion.all,criterion)
        }
        criterion.elcic<-rbind(criterion.elcic,criterion.all)
        #print(j)
    }
    colnames(criterion.elcic)<-candidate.sets
    rownames(criterion.elcic)<-candidate.cor.sets
    criterion.elcic
        }
    }else{
        if(!is.null(name.var.sets))
        {
                criterion.elcic<-rep()
                corstr<-candidate.cor.sets[1]
                for (i in seq_len(length(name.var.sets)))
                {
                    criterion<-ELCIC.gee.single(x=x,y=y,r=r,id=id,time=time,index.var=NULL,name.var=name.var.sets[[i]],dist=dist,corstr=corstr)
                    criterion.elcic<-c(criterion.elcic,criterion)
                }

            names(criterion.elcic)<-name.var.sets
            criterion.elcic
        }else{
                criterion.elcic<-rep()
                corstr<-candidate.cor.sets[1]
                for (i in seq_len(length(candidate.sets)))
                {
                    criterion<-ELCIC.gee.single(x=x,y=y,r=r,id=id,time=time,index.var=candidate.sets[[i]],name.var=NULL,dist=dist,corstr=corstr)
                    criterion.elcic<-c(criterion.elcic,criterion)
                }

            names(criterion.elcic)<-candidate.sets
            criterion.elcic
        }
    }
}




#'@title The whole procedure for joint selection of marginal mean and correlation structures for longitudinal data with dropout missingness under the mechanism of missing at random
#'@description The function \code{\link{ELCIC.wgee}} provides the overall procedure for joint selection of marginal mean and correlation structures for longitudinal data with dropout missingness under the mechanism of missing at random. It is also able to implement the marginal mean structure selection given a prespecified working correlation structure. The data is dropout missing and missing at random.
#'@usage ELCIC.wgee(x,y,x_mis,r,id,time,candidate.sets=NULL,name.var.sets=NULL,
#'      dist,candidate.cor.sets=c("independence","exchangeable", "ar1"), joints=TRUE,lag=1)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept if the intercept is considered in the marginal mean. The covariate matrix should be complete.
#'@param y A vector containing outcomes. Use NA to indicate missing outcomes.
#'@param x_mis A matrix containing covariates for the missing data model. The first column should be all ones corresponding to the intercept. This covariate matrix should be complete and all observed. See more in Details section.
#'@param r A vector indicating the observation of outcomes: 1 for an observed record, and 0 for an unobserved record.
#'@param id A vector indicating subject id.
#'@param time The number of observations in total for each subject
#'@param candidate.sets A list containing index corresponding to covariates in the candidate models. See more in Details section.
#'@param name.var.sets A list containing names of in the candidate models. The names should be a subset of the column names of the x matrix. See more in Details section.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param candidate.cor.sets A vector containing candidate correlation structures. When joints=TRUE, it can be any subset of c("independence", "exchangeable", "ar1"). The default is c("independence", "exchangeable", "ar1"). When joints=FALSE, it should be either of "independence","exchangeable", "ar1". See more in Details section.
#'@param joints A logic value for joint selection of marginal mean and working correlation structures. The default is TRUE. See more in Details section.
#'@param lag A numeric value indicating lag-response involved in the missing data model. It can be either 1 or 2. The default is 1.

#'@return A matrix with each element containing ELCIC value for each candidate model.
#'
#'@details Covariate matrix "x" should be complete. If missing data are present in "x", the elements in covariate vector will be replaced by zeros for individuals who have missing covariates.
#'
#'The argument "x_mis" includes all covariates to fit the missing data model. It does not contains a lag variable based on the outcome y. The argument "lag" in this function will automatically add lag-response variables to indicate the data missing at random.
#'
#'Either arguments "candidate.sets" or "name.var.sets" is used to identify the set of candidate mean models. If both arguments are provided, only the argument "name.var.sets" will be used.
#'
#'When joints=TRUE, the argument "candidate.cor.sets" can contain multiple correlation structures; however, when joints=FALSE, it should contain either of "independence", "exchangeable", and "ar1". If multiple correlation structures are provided, only the first one will be used.

#'@examples
#'## tests
#'# load data
#'data(wgeesimdata)
#'dist<-"binomial"
#'x<-wgeesimdata$x
#'y<-wgeesimdata$y
#'x_mis<-wgeesimdata$x_mis
#'r<-wgeesimdata$obs_ind
#'id<-wgeesimdata$id
#'time<-3
#'candidate.sets<-list(c(1,2,3))
#'candidate.cor.sets<-c("exchangeable")
#'criterion.elcic<-ELCIC.wgee(x,y,x_mis,r,id,time,candidate.sets,name.var.sets=NULL,
#'                                     dist,candidate.cor.sets,joints=TRUE)
#'criterion.elcic
#'
#'@export

ELCIC.wgee<-function(x,y,x_mis,r,id,time,candidate.sets=NULL,name.var.sets=NULL,dist,candidate.cor.sets=c("independence","exchangeable", "ar1"),joints=TRUE,lag=1)
{
    if(joints)
    {
        if(!is.null(name.var.sets))
        {
            criterion.elcic<-rep()
            for(j in seq_len(length(candidate.cor.sets)))
            {
                corstr<-candidate.cor.sets[j]
                criterion.all<-rep()
                for (i in seq_len(length(name.var.sets)))
                {
                    criterion<-ELCIC.wgee.single(x,y,x_mis,r,id,time,index.var=NULL,name.var=name.var.sets[[i]],dist,corstr,joints=TRUE,lag=lag)
                    criterion.all<-c(criterion.all,criterion)
                }
                criterion.elcic<-rbind(criterion.elcic,criterion.all)
                #print(j)
            }
            colnames(criterion.elcic)<-name.var.sets
            rownames(criterion.elcic)<-candidate.cor.sets
            criterion.elcic
        }else{
            criterion.elcic<-rep()
            for(j in seq_len(length(candidate.cor.sets)))
            {
                corstr<-candidate.cor.sets[j]
                criterion.all<-rep()
                for (i in seq_len(length(candidate.sets)))
                {
                    criterion<-ELCIC.wgee.single(x,y,x_mis,r,id,time,index.var=candidate.sets[[i]],name.var=NULL,dist,corstr,joints=TRUE,lag=lag)
                    criterion.all<-c(criterion.all,criterion)
                }
                criterion.elcic<-rbind(criterion.elcic,criterion.all)
                #print(j)
            }
            colnames(criterion.elcic)<-candidate.sets
            rownames(criterion.elcic)<-candidate.cor.sets
            criterion.elcic
        }
    }else{
        if(!is.null(name.var.sets))
        {
                criterion.elcic<-rep()
                corstr<-candidate.cor.sets[1]
                for (i in seq_len(length(name.var.sets)))
                {
                    criterion<-ELCIC.wgee.single(x,y,x_mis,r,id,time,index.var=NULL,name.var=name.var.sets[[i]],dist,corstr,joints=FALSE,lag=lag)
                    criterion.elcic<-c(criterion.elcic,criterion)
                }

            names(criterion.elcic)<-name.var.sets
            criterion.elcic
        }else{
            criterion.elcic<-rep()
                corstr<-candidate.cor.sets[1]
                for (i in seq_len(length(candidate.sets)))
                {
                    criterion<-ELCIC.wgee.single(x,y,x_mis,r,id,time,index.var=candidate.sets[[i]],name.var=NULL,dist,corstr,joints=FALSE,lag=lag)
                    criterion.elcic<-c(criterion.elcic,criterion)
                }

            names(criterion.elcic)<-candidate.sets
            criterion.elcic
        }
    }
}



#'@title the joint selection procedure of marginal mean and correlation structures in longitudinal data based on QIC
#'@description This function provides the joint selection of marginal mean and correlation structures in longitudinal data based on QIC.
#'@usage QICc.gee(x,y,id,dist,candidate.sets=NULL, name.var.sets=NULL,
#'    candidate.cor.sets=c("independence","exchangeable", "ar1"), joints=TRUE)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept. The covariate matrix should be complete.
#'@param y A vector containing outcomes.
#'@param id A vector indicating subject id.
#'@param candidate.sets A list containing index corresponding to covariates in the candidate models.
#'@param name.var.sets A list containing names of covariates in the candidate models. The names should be a subset of the column names of the x matrix.
#'@param dist A specified distribution. It can be "gaussian", "poisson", and "binomial".
#'@param candidate.cor.sets A vector containing condidate correlation structures. When joints=TRUE, it can be any subset of c("independence", "exchangeable", "ar1"). The default is c("independence", "exchangeable", "ar1"). When joints=FALSE, it should be either of "independence","exchangeable", "ar1". See more in Details section.
#'@param joints A logic value for joint selection of marginal mean and working correlation structures. The default is TRUE.
#'
#'@return A vector with each element containing QIC value for each candidate model. The row name of this vector is the selected correlation structure.
#'
#'@details Either arguments "index.var" or "name.var" is used to identify the candidate mean models. If both arguments are provided, only the argument "name.var" will be used.
#'
#'When joints=TRUE, the argument "candidate.cor.sets" can contain multiple correlation structures; however, when joints=FALSE, it should contain either of "independence", "exchangeable", and "ar1". If multiple correlation structures are provided, only the first one will be used.
#'
#'@examples
#'## tests
#'# load data
#'data(geesimdata)
#'x<-geesimdata$x
#'y<-geesimdata$y
#'id<-geesimdata$id
#'r<-rep(1,nrow(x))
#'time<-3
#'candidate.sets<-list(c(1,2),c(1,2,3))
#'candidate.cor.sets<-c("exchangeable")
#'dist="poisson"
#'criterion.qic<-QICc.gee(x=x,y=y,id=id,dist=dist,candidate.sets=candidate.sets,
#'                     name.var.sets=NULL,candidate.cor.sets=candidate.cor.sets)
#'criterion.qic
#'
#'@export
#'@importFrom geepack QIC

QICc.gee<-function (x,y,id,dist,candidate.sets=NULL,name.var.sets=NULL,candidate.cor.sets=c("independence","exchangeable", "ar1"),joints=TRUE)
{
    if(joints)
    {
        if(!is.null(name.var.sets))
        {
    data<-data.frame(y,x)
    modelQ<-rep()
    for(i in seq_len(length(candidate.cor.sets)))
    {
        fit<-geeglm(y~x-1,data=data,family =dist,id=id,corstr = candidate.cor.sets[i])
        modelQ[i]<-as.numeric(QIC(fit)[4])##look at QIC and QICu
    }
    QICcorstr<-candidate.cor.sets[which.min(modelQ)]

    criterion.mean<-rep()
    for (i in seq_len(length(name.var.sets)))
    {
        x.candidate<-x[,name.var.sets[[i]]]
        data<-data.frame(y,x.candidate)
        fit<-geeglm(y~x.candidate-1,data=data,family =dist,id=id,corstr = QICcorstr)
        criterion.mean<-c(criterion.mean,QIC(fit)[1])
    }
    criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
    rownames(criterion.mean)<-QICcorstr
    colnames(criterion.mean)<-name.var.sets
    criterion.mean
        }else{
            data<-data.frame(y,x)
            modelQ<-rep()
            for(i in seq_len(length(candidate.cor.sets)))
            {
                fit<-geeglm(y~x-1,data=data,family =dist,id=id,corstr = candidate.cor.sets[i])
                modelQ[i]<-as.numeric(QIC(fit)[4])##look at QIC and QICu
            }
            QICcorstr<-candidate.cor.sets[which.min(modelQ)]

            criterion.mean<-rep()
            for (i in seq_len(length(candidate.sets)))
            {
                x.candidate<-x[,candidate.sets[[i]]]
                data<-data.frame(y,x.candidate)
                fit<-geeglm(y~x.candidate-1,data=data,family =dist,id=id,corstr = QICcorstr)
                criterion.mean<-c(criterion.mean,QIC(fit)[1])
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            rownames(criterion.mean)<-QICcorstr
            colnames(criterion.mean)<-candidate.sets
            criterion.mean
        }
    }else{
        if(!is.null(name.var.sets))
        {
            data<-data.frame(y,x)
            QICcorstr<-candidate.cor.sets[1]

            criterion.mean<-rep()
            for (i in seq_len(length(name.var.sets)))
            {
                x.candidate<-x[,name.var.sets[[i]]]
                data<-data.frame(y,x.candidate)
                fit<-geeglm(y~x.candidate-1,data=data,family =dist,id=id,corstr = QICcorstr)
                criterion.mean<-c(criterion.mean,QIC(fit)[1])
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            rownames(criterion.mean)<-QICcorstr
            colnames(criterion.mean)<-name.var.sets
            criterion.mean
        }else{
            data<-data.frame(y,x)
            QICcorstr<-candidate.cor.sets[1]

            criterion.mean<-rep()
            for (i in seq_len(length(candidate.sets)))
            {
                x.candidate<-x[,candidate.sets[[i]]]
                data<-data.frame(y,x.candidate)
                fit<-geeglm(y~x.candidate-1,data=data,family =dist,id=id,corstr = QICcorstr)
                criterion.mean<-c(criterion.mean,QIC(fit)[1])
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            rownames(criterion.mean)<-QICcorstr
            colnames(criterion.mean)<-candidate.sets
            criterion.mean
        }
    }
}




#'@title The whole MLIC procedure for joint selection of marginal mean and correlation structures for longitudinal data with dropout missingness under the mechanism of missing at random
#'@description This function provides the overall MLIC procedure for joint selection of marginal mean and correlation structures for longitudinal data with dropout missingness under the mechanism of missing at random. It is also able to implement the marginal mean structure selection given a prespecified "working" correlation structure. The data is dropout missing and missing at random.
#'@usage MLIC.wgee(x,y,x_mis,r,id,time,candidate.sets=NULL, name.var.sets=NULL,dist,
#'       candidate.cor.sets=c("independence","exchangeable", "ar1"), joints=TRUE,lag=1)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept if the intercept is expected in the marginal mean model. The covariate matrix should be complete.
#'@param y A vector containing outcomes. Use NA to indicate missing outcomes.
#'@param x_mis A matrix containing covariates for the missing data model. The first column should be all ones corresponding to the intercept. This covariate matrix should be complete and all observed. See more in Details section.
#'@param r A vector indicating the observation of outcomes: 1 for an observed record, and 0 for an unobserved record.
#'@param id A vector indicating subject id.
#'@param time The number of observations in total for each subject.
#'@param candidate.sets A list containing index corresponding to covariates in the candidate models. See more in Details section.
#'@param name.var.sets A list containing names of covariates in the candidate models. The names should be a subset of the column names of x matrix. See more in Details section.
#'@param dist A specified distribution. It can be "gaussian", "poisson", and "binomial".
#'@param candidate.cor.sets A vector containing condidate correlation structures. When joints=TRUE, it is c("independence", "exchangeable", "ar1") as default. When joints=FALSE, it should be either of "independence", "exchangeable", "ar1". See more in Details section.
#'@param joints A logic value for joint selection of marginal mean and working correlation structures. The default is TRUE. See more in Details section.
#'@param lag A numeric value indicating lag-response involved in the missing data model. It can be either 1 or 2. The default is 1.

#'@return A vector with each element containing MLIC value for each candidate model. The row name of this vector is the selected correlation structure.
#'
#'@details The covariate matrix "x" should be complete. If missing data are present in "x", the elements in covariate vector will be replaced by zeros for individuals who have missing covariates.
#'
#'The argument "x_mis" includes all covariates to fit the missing data model. It does not contains a lag variable based on the outcome y. The argument "lag" in this function will automatically add lag-response variables to indicate the data missing at random.
#'
#'Either arguments "candidate.sets" or "name.var.sets" is used to identify the set of candidate mean models. If both arguments are provided, only the argument "name.var.sets" will be used.
#'
#'When joints=TRUE, the argument "candidate.cor.sets" can contain multiple correlation structures; however, when joints=FALSE, it should contain either of "independence","exchangeable", and "ar1". If multiple correlation structures are provided, only the first one will be used.
#'
#'@examples
#'## tests
#'# load data
#'data(wgeesimdata)
#'dist<-"binomial"
#'x<-wgeesimdata$x
#'y<-wgeesimdata$y
#'x_mis<-wgeesimdata$x_mis
#'r<-wgeesimdata$obs_ind
#'id<-wgeesimdata$id
#'time=3
#'candidate.sets<-list(c(1,2,3))
#'candidate.cor.sets<-c("exchangeable")
#'criterion.mlic<-MLIC.wgee(x,y,x_mis,r,id,time,candidate.sets,
#'             name.var.sets=NULL,dist,candidate.cor.sets,joints=FALSE)
#'criterion.mlic

#'@export
#'@importFrom wgeesel wgee data_sim MLIC.gee QICW.gee

MLIC.wgee<-function(x,y,x_mis,r,id,time,candidate.sets=NULL,name.var.sets=NULL,dist,candidate.cor.sets=c("independence","exchangeable", "ar1"),joints=TRUE,lag=1)
{
    #replace missingness with 0 in x matrix
    na.index<-which(is.na(x),arr.ind=TRUE)
    r.x<-rep(1,nrow(x))
    r.x[unique(na.index[,1])]<-0
    if(nrow(na.index)>0)
    {
        if (sum(r.x==r)!=nrow(x))
        {
            stop("missingness pattern in x should be the same to the missing pattern in y.")
        }
        x[unique(na.index[,1]),]<-0
        warning("Covariate matrix x should be fully observed. The elements in covariate vector are replaced by zeros for individuals who have missing covariates.")
    }

    if(lag>2 | lag>=time)
    {
        stop("Invalid type of lag. It should be less than 3 and time")
    }else if(lag==0)
    {
        warning("No lag of response added may indicate missing completely at random. GEE may be used.")
    }

    #generate ylag1
    if(lag!=0)
    {
        samplesize<-length(unique(id))
        ylag.cov<-ylag(y=y,id=id,lag=1)
        ylag.cov[is.na(ylag.cov)]<-0
        ylag.cov[rep(seq_len(time),time=samplesize)==1]<-NA
        x_mis<-cbind(x_mis,ylag1.cov=ylag.cov)
    }
    if(lag==2){
        #generate ylag2
        ylag.cov<-ylag(y=y,id=id,lag=lag)
        x_mis<-cbind(x_mis,ylag2.cov=ylag.cov)
    }

    if(joints)
    {
        if(!is.null(name.var.sets))
        {
            y<-as.matrix(y,col=1)
            #data_wgee<-data.frame(y,x,r,x_mis)
            data_wgee<-data.frame(r,x_mis)

            mismodel_formula<-paste("r~",colnames(data_wgee)[3])
            if(ncol(data_wgee)>3)
            {
                for(jj in 4:ncol(data_wgee))
                {
                    mismodel_formula<-paste(mismodel_formula,"+",colnames(data_wgee)[jj])
                }
            }
            mismodel_formula<-as.formula(mismodel_formula)
            modelQ<-rep()
            for(i in seq_len(length(candidate.cor.sets)))
            {
                fit<-wgee(y~x-1,data_wgee,id,family=dist,corstr =candidate.cor.sets[[i]],scale = NULL,mismodel =mismodel_formula)
                modelQ[i]<-MLIC.gee(fit,fit)$MLICc##look at MLICc
            }
            MLICcorstr<-candidate.cor.sets[which.min(modelQ)]

            #fit a full model
            fitf<-wgee(y~x-1,data_wgee,id,family=dist,corstr = MLICcorstr,scale = NULL,mismodel =mismodel_formula)

            criterion.mean<-rep()
            for (i in seq_len(length(name.var.sets)))
            {
                x.candidate<-x[,name.var.sets[[i]]]
                fitm<-wgee(y~x.candidate-1,data_wgee,id,family=dist,corstr =MLICcorstr,scale = NULL,mismodel =mismodel_formula)
                criterion.mean<-c(criterion.mean,MLIC.gee(fitm,fitf)$MLIC)
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            rownames(criterion.mean)<-MLICcorstr
            colnames(criterion.mean)<-name.var.sets
            criterion.mean
        }else{
            y<-as.matrix(y,col=1)
            #data_wgee<-data.frame(y,x,r,x_mis)
            data_wgee<-data.frame(r,x_mis)

            mismodel_formula<-paste("r~",colnames(data_wgee)[3])
            if(ncol(data_wgee)>3)
            {
                for(jj in 4:ncol(data_wgee))
                {
                    mismodel_formula<-paste(mismodel_formula,"+",colnames(data_wgee)[jj])
                }
            }
            mismodel_formula<-as.formula(mismodel_formula)
            modelQ<-rep()
            for(i in seq_len(length(candidate.cor.sets)))
            {
                fit<-wgee(y~x-1,data_wgee,id,family=dist,corstr =candidate.cor.sets[[i]],scale = NULL,mismodel =mismodel_formula)
                modelQ[i]<-MLIC.gee(fit,fit)$MLICc##look at MLICc
            }
            MLICcorstr<-candidate.cor.sets[which.min(modelQ)]

            #fit a full model
            fitf<-wgee(y~x-1,data_wgee,id,family=dist,corstr = MLICcorstr,scale = NULL,mismodel =mismodel_formula)

            criterion.mean<-rep()
            for (i in seq_len(length(candidate.sets)))
            {
                x.candidate<-x[,candidate.sets[[i]]]
                fitm<-wgee(y~x.candidate-1,data_wgee,id,family=dist,corstr =MLICcorstr,scale = NULL,mismodel =mismodel_formula)
                criterion.mean<-c(criterion.mean,MLIC.gee(fitm,fitf)$MLIC)
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            rownames(criterion.mean)<-MLICcorstr
            colnames(criterion.mean)<-candidate.sets
            criterion.mean
        }
    }else{
        if(!is.null(name.var.sets))
        {
            y<-as.matrix(y,col=1)
            #data_wgee<-data.frame(y,x,r,x_mis)
            data_wgee<-data.frame(r,x_mis)

            mismodel_formula<-paste("r~",colnames(data_wgee)[3])
            if(ncol(data_wgee)>3)
            {
                for(jj in 4:ncol(data_wgee))
                {
                    mismodel_formula<-paste(mismodel_formula,"+",colnames(data_wgee)[jj])
                }
            }
            mismodel_formula<-as.formula(mismodel_formula)

            MLICcorstr<-candidate.cor.sets[1]

            #fit a full model
            fitf<-wgee(y~x-1,data_wgee,id,family=dist,corstr = MLICcorstr,scale = NULL,mismodel =mismodel_formula)

            criterion.mean<-rep()
            for (i in seq_len(length(name.var.sets)))
            {
                x.candidate<-x[,name.var.sets[[i]]]
                fitm<-wgee(y~x.candidate-1,data_wgee,id,family=dist,corstr =MLICcorstr,scale = NULL,mismodel =mismodel_formula)
                criterion.mean<-c(criterion.mean,MLIC.gee(fitm,fitf)$MLIC)
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            colnames(criterion.mean)<-name.var.sets
            rownames(criterion.mean)<-MLICcorstr
            criterion.mean
        }else{
            y<-as.matrix(y,col=1)
            #data_wgee<-data.frame(y,x,r,x_mis)
            data_wgee<-data.frame(r,x_mis)

            mismodel_formula<-paste("r~",colnames(data_wgee)[3])
            if(ncol(data_wgee)>3)
            {
                for(jj in 4:ncol(data_wgee))
                {
                    mismodel_formula<-paste(mismodel_formula,"+",colnames(data_wgee)[jj])
                }
            }
            mismodel_formula<-as.formula(mismodel_formula)

            MLICcorstr<-candidate.cor.sets[1]

            #fit a full model
            fitf<-wgee(y~x-1,data_wgee,id,family=dist,corstr = MLICcorstr,scale = NULL,mismodel =mismodel_formula)

            criterion.mean<-rep()
            for (i in seq_len(length(candidate.sets)))
            {
                x.candidate<-x[,candidate.sets[[i]]]
                fitm<-wgee(y~x.candidate-1,data_wgee,id,family=dist,corstr =MLICcorstr,scale = NULL,mismodel =mismodel_formula)
                criterion.mean<-c(criterion.mean,MLIC.gee(fitm,fitf)$MLIC)
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            colnames(criterion.mean)<-candidate.sets
            rownames(criterion.mean)<-MLICcorstr
            criterion.mean
        }
    }
}




#'@title The whole QICW procedure for joint selection of marginal mean and correlation structures for longitudinal data with dropout missingness under the mechanism of missing at random
#'@description This function provides the overall QICW procedure for joint selection of marginal mean and correlation structures for longitudinal data with dropout missingness under the mechanism of missing at random. It is also able to implement the marginal mean structure selection given a prespecified "working" correlation structure. The data is dropout missing and missing at random.
#'@usage QICW.wgee(x,y,x_mis,r,id,time,candidate.sets,name.var.sets=NULL,
#'      dist,candidate.cor.sets=c("independence","exchangeable", "ar1"), joints=TRUE,lag=1)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept if the intercept is expected in the marginal mean model. The covariate matrix should be complete.
#'@param y A vector containing outcomes. Use NA to indicate missing outcomes.
#'@param x_mis A matrix containing covariates for the missing data model. The first column should be all ones corresponding to the intercept. This covariate matrix should be complete and all observed. See more in Details section.
#'@param r A vector indicating the observation of outcomes: 1 for an observed record, and 0 for an unobserved record.
#'@param id A vector indicating subject id.
#'@param time The number of observations in total for each subject.
#'@param candidate.sets A list containing index corresponding to covariates in the candidate models. See more in details section.
#'@param name.var.sets A list containing names of covariates in the candidate models. The names should be a subset of the column names of the x matrix. See more in Details section.
#'@param dist A specified distribution. It can be "gaussian", "poisson", and "binomial".
#'@param candidate.cor.sets A vector containing candidate correlation structures. When joints=TRUE, it is c("independence", "exchangeable", "ar1") as default. When joints=FALSE, it should be either of "independence", "exchangeable", "ar1". See more in Details section.
#'@param joints A logic value for joint selection of marginal mean and working correlation structures. The default is TRUE. See more in Details section.
#'@param lag A numeric value indicating lag-response involved in the missing data model. It can be either 1 or 2. The default is 1.

#'@return A vector with each element containing QICW value for each candidate model. The row name of this vector is the selected correlation structure.
#'
#'@details The covariate matrix "x" should be complete. If missing data are present in "x", the elements in covariate vector will be replaced by zeros for individuals who have missing covariates.
#'
#'The argument "x_mis" includes all covariates to fit the missing data model. It does not contains a lag variable based on the outcome y. The argument "lag" in this function will automatically add lag-response variables to indicate the data missing at random.
#'
#'Either arguments "candidate.sets" or "name.var.sets" is used to identify the set of candidate mean models. If both arguments are provided, only the argument "name.var.sets" will be used.
#'
#'When joints=TRUE, the argument "candidate.cor.sets" can contain multiple correlation structures; however, when joints=FALSE, it should contain either of "independence","exchangeable", and "ar1". If multiple correlation structures are provided, only the first one will be used.
#'
#'@examples
#'## tests
#'# load data
#'data(wgeesimdata)
#'dist="binomial"
#'x<-wgeesimdata$x
#'y<-wgeesimdata$y
#'x_mis<-wgeesimdata$x_mis
#'r<-wgeesimdata$obs_ind
#'id<-wgeesimdata$id
#'time=3
#'candidate.sets<-list(c(1,2,3))
#'candidate.cor.sets<-c("exchangeable")
#'criterion.qicw<-QICW.wgee(x,y,x_mis,r,id,time,candidate.sets,
#'           name.var.sets=NULL,dist,candidate.cor.sets,joints=FALSE)
#'criterion.qicw
#'
#'@export
#'@importFrom wgeesel wgee data_sim MLIC.gee QICW.gee

QICW.wgee<-function(x,y,x_mis,r,id,time,candidate.sets=NULL,name.var.sets=NULL,dist,candidate.cor.sets=c("independence","exchangeable", "ar1"),joints=TRUE,lag=1)
{
    #replace missingness with 0 in x matrix
    na.index<-which(is.na(x),arr.ind=TRUE)
    r.x<-rep(1,nrow(x))
    r.x[unique(na.index[,1])]<-0
    if(nrow(na.index)>0)
    {
        if (sum(r.x==r)!=nrow(x))
        {
            stop("missingness pattern in x should be the same to the missing pattern in y.")
        }
        x[unique(na.index[,1]),]<-0
        warning("Covariate matrix x should be fully observed. The elements in covariate vector are replaced by zeros for individuals who have missing covariates.")
    }

    if(lag>2 | lag>=time)
    {
        stop("Invalid type of lag. It should be less than 3 and time")
    }else if(lag==0)
    {
        warning("No lag of response added may indicate missing completely at random. GEE may be used.")
    }

    #generate ylag1
    if(lag!=0)
    {
        samplesize<-length(unique(id))
        ylag.cov<-ylag(y=y,id=id,lag=1)
        ylag.cov[is.na(ylag.cov)]<-0
        ylag.cov[rep(seq_len(time),time=samplesize)==1]<-NA
        x_mis<-cbind(x_mis,ylag1.cov=ylag.cov)
    }
    if(lag==2){
        #generate ylag2
        ylag.cov<-ylag(y=y,id=id,lag=lag)
        x_mis<-cbind(x_mis,ylag2.cov=ylag.cov)
    }

    if(joints)
    {
        if(!is.null(name.var.sets))
        {
            y<-as.matrix(y,col=1)
            #data_wgee<-data.frame(y,x,r,x_mis)
            data_wgee<-data.frame(r,x_mis)

            mismodel_formula<-paste("r~",colnames(data_wgee)[3])
            if(ncol(data_wgee)>3)
            {
                for(jj in 4:ncol(data_wgee))
                {
                    mismodel_formula<-paste(mismodel_formula,"+",colnames(data_wgee)[jj])
                }
            }
            mismodel_formula<-as.formula(mismodel_formula)
            modelQ<-rep()
            for(i in seq_len(length(candidate.cor.sets)))
            {
                fit<-wgee(y~x-1,data_wgee,id,family=dist,corstr =candidate.cor.sets[[i]],scale = NULL,mismodel =mismodel_formula)
                modelQ[i]<-QICW.gee(fit)$QICWr##look at QICWr
            }
            QICWcorstr<-candidate.cor.sets[which.min(modelQ)]

            criterion.mean<-rep()
            for (i in seq_len(length(name.var.sets)))
            {
                x.candidate<-x[,name.var.sets[[i]]]
                fitm<-wgee(y~x.candidate-1,data_wgee,id,family=dist,corstr =QICWcorstr,scale = NULL,mismodel =mismodel_formula)
                criterion.mean<-c(criterion.mean,QICW.gee(fitm)$QICWr)
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            rownames(criterion.mean)<-QICWcorstr
            colnames(criterion.mean)<-name.var.sets
            criterion.mean
        }else{
            y<-as.matrix(y,col=1)
            #data_wgee<-data.frame(y,x,r,x_mis)
            data_wgee<-data.frame(r,x_mis)

            mismodel_formula<-paste("r~",colnames(data_wgee)[3])
            if(ncol(data_wgee)>3)
            {
                for(jj in 4:ncol(data_wgee))
                {
                    mismodel_formula<-paste(mismodel_formula,"+",colnames(data_wgee)[jj])
                }
            }
            mismodel_formula<-as.formula(mismodel_formula)
            modelQ<-rep()
            for(i in seq_len(length(candidate.cor.sets)))
            {
                fit<-wgee(y~x-1,data_wgee,id,family=dist,corstr =candidate.cor.sets[[i]],scale = NULL,mismodel =mismodel_formula)
                modelQ[i]<-QICW.gee(fit)$QICWr##look at QICWr
            }
            QICWcorstr<-candidate.cor.sets[which.min(modelQ)]

            criterion.mean<-rep()
            for (i in seq_len(length(candidate.sets)))
            {
                x.candidate<-x[,candidate.sets[[i]]]
                fitm<-wgee(y~x.candidate-1,data_wgee,id,family=dist,corstr =QICWcorstr,scale = NULL,mismodel =mismodel_formula)
                criterion.mean<-c(criterion.mean,QICW.gee(fitm)$QICWr)
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            rownames(criterion.mean)<-QICWcorstr
            colnames(criterion.mean)<-candidate.sets
            criterion.mean
        }
    }else{
        if(!is.null(name.var.sets))
        {
            y<-as.matrix(y,col=1)
            #data_wgee<-data.frame(y,x,r,x_mis)
            data_wgee<-data.frame(r,x_mis)

            mismodel_formula<-paste("r~",colnames(data_wgee)[3])
            if(ncol(data_wgee)>3)
            {
                for(jj in 4:ncol(data_wgee))
                {
                    mismodel_formula<-paste(mismodel_formula,"+",colnames(data_wgee)[jj])
                }
            }
            mismodel_formula<-as.formula(mismodel_formula)
            QICWcorstr<-candidate.cor.sets[1]

            criterion.mean<-rep()
            for (i in seq_len(length(name.var.sets)))
            {
                x.candidate<-x[,name.var.sets[[i]]]
                fitm<-wgee(y~x.candidate-1,data_wgee,id,family=dist,corstr =QICWcorstr,scale = NULL,mismodel =mismodel_formula)
                criterion.mean<-c(criterion.mean,QICW.gee(fitm)$QICWr)
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            rownames(criterion.mean)<-QICWcorstr
            colnames(criterion.mean)<-name.var.sets
            criterion.mean
        }else{
            y<-as.matrix(y,col=1)
            #data_wgee<-data.frame(y,x,r,x_mis)
            data_wgee<-data.frame(r,x_mis)

            mismodel_formula<-paste("r~",colnames(data_wgee)[3])
            if(ncol(data_wgee)>3)
            {
                for(jj in 4:ncol(data_wgee))
                {
                    mismodel_formula<-paste(mismodel_formula,"+",colnames(data_wgee)[jj])
                }
            }
            mismodel_formula<-as.formula(mismodel_formula)
            QICWcorstr<-candidate.cor.sets[1]

            criterion.mean<-rep()
            for (i in seq_len(length(candidate.sets)))
            {
                x.candidate<-x[,candidate.sets[[i]]]
                fitm<-wgee(y~x.candidate-1,data_wgee,id,family=dist,corstr =QICWcorstr,scale = NULL,mismodel =mismodel_formula)
                criterion.mean<-c(criterion.mean,QICW.gee(fitm)$QICWr)
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            rownames(criterion.mean)<-QICWcorstr
            colnames(criterion.mean)<-candidate.sets
            criterion.mean
        }
    }
}


